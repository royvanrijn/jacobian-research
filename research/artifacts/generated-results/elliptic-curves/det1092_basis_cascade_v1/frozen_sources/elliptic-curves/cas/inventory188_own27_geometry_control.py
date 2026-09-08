#!/usr/bin/env python3
"""Geometry-only control using the prepublication local27 seed and an existing policy."""
import argparse
from pathlib import Path
import certify_compact_r17_candidates as cert
from memory_rank_certificate import checked_rank
from research_runtime.store import checkpoint,digest

ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas'
ART=ROOT/'artifacts/generated-results/elliptic-curves'
D=ROOT/'artifacts/local/elliptic-curves/inventory188-own27-geometry-control-v1'
SEED=D/'seed.json';ROW={'initial_rank':27}
SOURCE=ART/'full11952_late64_r17_results_v1.json'


def sources():
    names=['inventory188_own27_geometry_control.py','prepare_inventory188_own27_geometry.sage',
           'prepare_fresh_r17_pari_batch.sage','prospective_half_lattice_v2.sage',
           'half_lattice_pointed_sieve.py','memory_rank_certificate.py',
           'search_observability.py','research_runtime/store.py']
    return {str(p.relative_to(ROOT)):cert.hashed(p) for p in [SOURCE,*(CAS/n for n in names)]}


def prepare():
    if (D/'protocol.json').exists() or SEED.exists():raise FileExistsError('preserve known-only geometry control')
    row=next(r for r in cert.read(SOURCE)['curves'] if r['id']=='11952-0959582')
    if row['parameter']!='110314/102227' or row['rank_lower_bound']!=27 or len(row['generic_points'])!=17 or row['points'][:17]!=row['generic_points']:
        raise ArithmeticError('old local27 and generic17 prefix required')
    model=tuple(map(cert.F,row['curve']));points=[tuple(map(cert.F,P)) for P in row['points']];proof=row['rank_certificate']
    if digest(checked_rank(model,points,[r['prime'] for r in proof['signatures']],proof['no_rational_2_torsion_prime']))!=digest(proof):
        raise ArithmeticError('old27 basis proof differs')
    seed={k:row[k] for k in ('family','parameter','curve','generic_points','points','rank_certificate')};checkpoint(SEED,seed)
    p={'schema':'elliptic-curves.inventory188-own27-geometry.v1','sources':sources(),'seed_sha256':cert.hashed(SEED),
       'sample_size':2048,'sample_domain':'full11952-specialized-followup-v1','initial_rank':27,'generic_prefix':17,
       'charts':49,'height_for_later_visibility_audit':125000,'geometry_seconds':300,'rss_bytes':2147483648,
       'scope':'Geometry only: reuse the existing specialized-parity policy with its original fixed sample domain. Exactly2048 distinct27-bit masks with nonzero quotient above the17 generic prefix,384-bit canonical heights rounded at10^6, unimodular LLL and CVP, then49 largest computed norms. Inputs are only the earlier independently certified local27 seed and generic17 prefix. No public points, external witness coordinates, new parameter, score or point enumeration. Freeze all49 maps before a separate retrospective witness audit. No adaptive mask expansion or automatic point-search stage.'}
    digest(p);checkpoint(D/'protocol.json',p)


def protocol():
    p=cert.read(D/'protocol.json')
    if p['sources']!=sources() or p['seed_sha256']!=cert.hashed(SEED):raise ArithmeticError('geometry protocol changed')
    return p


def masks(p):
    result=[];i=0
    while len(result)<p['sample_size']:
        m=int(digest([p['sample_domain'],i]),16)%(1<<27);i+=1
        if m>>17 and m not in result:result.append(m)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare']);p.parse_args();prepare()
