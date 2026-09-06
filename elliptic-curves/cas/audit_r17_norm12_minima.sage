#!/usr/bin/env sage-python
"""Exact minimum audit on the98 computed-norm12 classes in two R17 censuses."""
import argparse
from pathlib import Path
import sys
from sage.all import matrix,ZZ
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';sys.path.insert(0,str(CAS))
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint,digest
from exact_parity_ellipsoid import enumerate_coset
ATLAS=ROOT/'artifacts/generated-results/elliptic-curves/compact_six_r17_atlas_v1.json'
PARENT=ROOT/'artifacts/local/elliptic-curves/compact-six-r17-h4096-v1'
D=ROOT/'artifacts/local/elliptic-curves/r17-norm12-exact-minima-v1'


def sources():
    return {str(p.relative_to(ROOT)):cert.hashed(p) for p in (Path(__file__).resolve(),
        CAS/'exact_parity_ellipsoid.py',ROOT/'elliptic-curves/tests/test_exact_parity_ellipsoid.py',ATLAS)}


def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve exact-minimum protocol')
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.r17-norm12-minima.v1','sources':sources(),
        'families':['11952','08f72'],'census_hashes':{f:cert.hashed(PARENT/f/'generic-census.json') for f in ('11952','08f72')},
        'queries_per_family':49,'worker_wall_seconds':300,'worker_rss_bytes':1073741824,'maximum_workers':2,
        'gate':'Fresh numerical CVP censuses record49 norm12 representatives in each of these two families. Check whether these queried representatives are actually minimal, including the six omitted by each fixed43 cap, before using the numerical geometry to schedule further work.',
        'method':'Integer unimodular LLL is only a change of coordinates. Check its integral inverse and Gram identity exactly. Exhaust all vectors in each transformed parity coset of norm<=10 using rational LDL intervals. Evenness of the original integral Gram and a norm12 witness then certify the exact queried minimum.',
        'claim_boundary':'Only98 queried parity minima, not all2^17 minima, the full covering radius, or any specialized rank. Censored queries remain UNKNOWN.'})


def run(family):
    protocol=cert.read(D/'protocol.json')
    if protocol['sources']!=sources():raise ArithmeticError('exact-minimum sources changed')
    path=PARENT/family/'generic-census.json'
    if cert.hashed(path)!=protocol['census_hashes'][family]:raise ArithmeticError('census changed')
    census=cert.read(path);f=next(r for r in cert.read(ATLAS)['families'] if r['family']==family)
    G=matrix(ZZ,f['generic_height_gram'])
    if [list(map(int,row)) for row in G]!=f['generic_height_gram'] or any(G[i,i]%2 for i in range(17)):
        raise ArithmeticError('integral even Gram required')
    U=G.LLL_gram();V=U.inverse()
    if any(q.denominator()!=1 for q in V.list()):raise ArithmeticError('basis transform not unimodular')
    V=matrix(ZZ,V);H=U.transpose()*G*U
    if U*V!=matrix.identity(ZZ,17) or V*U!=matrix.identity(ZZ,17):raise ArithmeticError('unimodular identities failed')
    convert=lambda M:[list(map(int,r)) for r in M]
    chosen=[r for r in census['records'] if cert.F(r['norm'])==12]
    if len(chosen)!=49:raise ArithmeticError('query roster changed')
    original_used={r['mask'] for r in census['selected']}
    output=D/family/'result.json'
    if output.exists():raise FileExistsError('preserve minimum audit')
    data={'schema':'elliptic-curves.r17-exact-queried-parity-minima.v1','status':'RUNNING','protocol_hash':digest(protocol),
          'family':family,'gram':convert(G),'transform':convert(U),'inverse_transform':convert(V),'reduced_gram':convert(H),'rows':[]}
    checkpoint(output,data)
    for row in chosen:
        mask=row['mask'];parity=[(mask>>j)&1 for j in range(17)]
        residue=[sum(int(V[i,j])*parity[j] for j in range(17))%2 for i in range(17)]
        enumeration=enumerate_coset(convert(H),residue,10)
        minimum=enumeration['minimum_within_bound']
        if minimum is None:minimum=12;witness=row['representative']
        else:witness=[sum(int(U[i,j])*enumeration['minimum_witness'][j] for j in range(17)) for i in range(17)]
        if any((witness[i]-parity[i])%2 for i in range(17)) or sum(witness[i]*int(G[i,j])*witness[j] for i in range(17) for j in range(17))!=minimum:
            raise ArithmeticError('exact minimum witness failed')
        data['rows'].append({'mask':mask,'used_in_initial43':mask in original_used,'reported_census_norm':row['norm'],'exact_minimum':minimum,'minimum_witness':witness,'transformed_parity':residue,'enumeration':enumeration})
        checkpoint(output,data);print('EXACT PARITY MINIMUM',family,mask,minimum,flush=True)
    data['status']='COMPLETE_DECLARED_EXACT_AUDIT';checkpoint(output,data)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('stage',choices=['prepare','run']);p.add_argument('--family');a=p.parse_args();prepare() if a.stage=='prepare' else run(a.family)
