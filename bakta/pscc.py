import logging
from concurrent.futures import ThreadPoolExecutor
from typing import Sequence, Tuple

import bakta.config as cfg
import bakta.constants as bc
import bakta.features.orf as orf
import bakta.utils as bu


############################################################################
# PSCC DB columns
############################################################################
DB_PSCC_COL_UNIREF50 = 'uniref50_id'
DB_PSCC_COL_PRODUCT = 'product'


log = logging.getLogger('PSCC')


def search(cdss: Sequence[dict]) -> Tuple[Sequence[dict], Sequence[dict], Sequence[dict]]:
    """Conduct homology search of CDSs against PSCC db."""
    cds_aa_path = cfg.tmp_path.joinpath('cds.pscc.faa')
    orf.write_internal_faa(cdss, cds_aa_path)
    diamond_output_path = cfg.tmp_path.joinpath('diamond.pscc.tsv')
    diamond_db_path = cfg.db_path.joinpath('pscc.dmnd')
    if(cfg.use_mmseqs):
        cmd = bu.build_mmseqs_easy_search_command(
            database=cfg.db_path.joinpath('pscc.fasta'),
            query=cds_aa_path,
            output=diamond_output_path,
            output_fields=('query', 'target', 'qlen', 'tlen', 'alnlen', 'pident', 'evalue', 'bits'),
            identity=bc.MIN_PSCC_IDENTITY,
            query_cover=bc.MIN_PSC_COVERAGE,
            subject_cover=bc.MIN_PSC_COVERAGE,
        )
    else:
        cmd = bu.build_diamond_blastp_command(
            database=diamond_db_path,
            query=cds_aa_path,
            output=diamond_output_path,
            identity=bc.MIN_PSCC_IDENTITY,
            query_cover=bc.MIN_PSC_COVERAGE,
            subject_cover=bc.MIN_PSC_COVERAGE,
        )
    log.debug('cmd=%s', cmd)
    proc = bu.run_tool(
        cmd,
        cwd=str(cfg.tmp_path),
        env=cfg.env,
    )
    if(proc.returncode != 0):
        log.debug('stdout=\'%s\', stderr=\'%s\'', proc.stdout, proc.stderr)
        tool = 'mmseqs' if cfg.use_mmseqs else 'diamond'
        log.warning('PSCC failed! %s-error-code=%d', tool, proc.returncode)
        raise Exception(f'{tool} error! error code: {proc.returncode}')

    cds_by_hexdigest = orf.get_orf_dictionary(cdss)
    with diamond_output_path.open() as fh:
        for line in fh:
            (aa_identifier, cluster_id, query_length, subject_length, alignment_length, identity, evalue, bitscore) = line.split('\t')
            cds = cds_by_hexdigest[aa_identifier]
            query_cov = int(alignment_length) / len(cds['aa'])
            subject_cov = int(alignment_length) / int(subject_length)
            identity = float(identity) / 100
            bitscore = float(bitscore)
            evalue = float(evalue)
            if(query_cov >= bc.MIN_PSC_COVERAGE and subject_cov >= bc.MIN_PSC_COVERAGE and identity >= bc.MIN_PSCC_IDENTITY):
                cds['pscc'] = {
                    DB_PSCC_COL_UNIREF50: cluster_id,
                    'query_cov': query_cov,
                    'subject_cov': subject_cov,
                    'identity': identity,
                    'score': bitscore,
                    'evalue': evalue
                }
                log.debug(
                    'homology: seq=%s, start=%i, stop=%i, strand=%s, aa-length=%i, query-cov=%0.3f, subject-cov=%0.3f, identity=%0.3f, score=%0.1f, evalue=%1.1e, UniRef50=%s',
                    cds['sequence'], cds['start'], cds['stop'], cds['strand'], len(cds['aa']), query_cov, subject_cov, identity, bitscore, evalue, cluster_id
                )

    psccs_found = []
    cds_not_found = []
    for cds in cdss:
        if('pscc' in cds):
            psccs_found.append(cds)
        else:
            cds_not_found.append(cds)
    log.info('found: PSCC=%i', len(psccs_found))
    return psccs_found, cds_not_found


def lookup(features: Sequence[dict], pseudo: bool = False):
    """Lookup PSCC information"""
    no_pscc_lookups = 0
    try:
        rec_futures = []
        with ThreadPoolExecutor(max_workers=max(10, cfg.threads)) as tpe:  # use min 10 threads for IO bound non-CPU lookups
            for feature in features:
                uniref50_id = None
                if(pseudo):  # if pseudogene use pseudogene info
                    if('psc' in feature[bc.PSEUDOGENE]):
                        uniref50_id = feature[bc.PSEUDOGENE]['psc'].get(DB_PSCC_COL_UNIREF50, None)
                else:
                    if('psc' in feature):
                        uniref50_id = feature['psc'].get(DB_PSCC_COL_UNIREF50, None)
                    elif('pscc' in feature):
                        uniref50_id = feature['pscc'].get(DB_PSCC_COL_UNIREF50, None)
                if(uniref50_id is not None):
                    if(bc.DB_PREFIX_UNIREF_50 in uniref50_id):
                        uniref50_id = uniref50_id[9:]  # remove 'UniRef50_' prefix
                    future = tpe.submit(fetch_db_pscc_result, uniref50_id)
                    rec_futures.append((feature, future))

        for (feature, future) in rec_futures:
            rec = future.result()
            if(rec is not None):
                pscc = parse_annotation(rec)
                if(pseudo):
                    feature[bc.PSEUDOGENE]['pscc'] = pscc
                else:
                    if('pscc' in feature):
                        feature['pscc'] = {**feature['pscc'], **pscc}  # merge dicts, add PSCC annotation info to PSCC alignment info
                    else:
                        feature['pscc'] = pscc  # add PSCC annotation info
                no_pscc_lookups += 1
                log.debug(
                    'lookup: seq=%s, start=%i, stop=%i, strand=%s, UniRef50=%s, product=%s',
                    feature['sequence'], feature['start'], feature['stop'], feature['strand'], pscc.get(DB_PSCC_COL_UNIREF50, ''), pscc.get(DB_PSCC_COL_PRODUCT, '')
                )
            else:
                log.debug('lookup: ID not found! uniref50_id=%s', uniref50_id)
    except Exception as ex:
        log.exception('Could not read PSCCs from db!', ex)
        raise Exception('SQL error!', ex)
    log.info('looked-up=%i', no_pscc_lookups)


def fetch_db_pscc_result(uniref50_id: str):
    with bu.get_db_connection() as conn:
        c = conn.cursor()
        c.execute('select * from pscc where uniref50_id=?', (uniref50_id,))
        rec = c.fetchone()
        c.close()
        return rec


def parse_annotation(rec) -> dict:
    uniref_full_id = bc.DB_PREFIX_UNIREF_50 + rec[DB_PSCC_COL_UNIREF50]
    pscc = {
        DB_PSCC_COL_UNIREF50: uniref_full_id,  # must not be NULL/None
        'db_xrefs': [
            'SO:0001217',
            f'{bc.DB_XREF_UNIREF}:{uniref_full_id}'
        ]
    }
    # add non-empty PSCC annotations and attach database prefixes to identifiers
    if(rec[DB_PSCC_COL_PRODUCT]):
        pscc[DB_PSCC_COL_PRODUCT] = rec[DB_PSCC_COL_PRODUCT]
    return pscc
