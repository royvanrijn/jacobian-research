#!/usr/bin/env sage-python
"""Exact scaled minimum audit on42 computed maximum classes in five MW16 censuses."""
import argparse
from pathlib import Path
import sys
from sage.all import matrix,ZZ
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';sys.path.insert(0,str(CAS))
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint,digest
from exact_parity_ellipsoid import enumerate_coset
ATLAS=ROOT/'artifacts/generated-results/elliptic-curves/compact_five_mw16_atlas_v1.json'
PARENT=ROOT/'artifacts/local/elliptic-curves/prospective-mw16-h4096-v1'
D=ROOT/'artifacts/local/elliptic-curves/mw16-exact-maximum-parities-v1'


def sources():
    return {str(p.relative_to(ROOT)):cert.hashed(p) for p in (Path(__file__).resolve(),
        CAS/'exact_parity_ellipsoid.py',ROOT/'elliptic-curves/tests/test_exact_parity_ellipsoid.py',ATLAS)}


def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve exact MW16 parity audit')
    families=[f'a1-fibration-{i:02}' for i in range(1,6)]
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.mw16-exact-maximum-parities.v1','sources':sources(),'families':families,'census_hashes':{f:cert.hashed(PARENT/f/'generic-census.json') for f in families},'query_counts':[12,4,8,10,8],'scale':2,'scaled_upper_bound':23,'worker_wall_seconds':300,'worker_rss_bytes':1073741824,'maximum_workers':1,'gate':'The saved complete MW16 censuses report a maximum witness norm23/2 and42 classes at that value across five families. The fixed43 chart cap is not an exact maximum-class theorem. Audit these42 potential maximum classes, and independently check every65536-class upper witness in each family, before claiming whether any maximum class was omitted. No specialized point, candidate rank or public target enters the audit.','method':'Scale the rational Gram by2 to an integral matrix. Check exact unimodular LLL coordinate changes. For each of42 proposed norm23 classes exhaust every vector of scaled norm<=22 by rational LDL enumeration. An empty ellipsoid plus the explicit norm23 witness certifies its minimum. Other parities have directly checked witnesses of scaled norm<=23.','claim_boundary':'At most42 queried exact minima and327680 parity upper witnesses. No continuous covering radius, saturation, specialization rank or automatic point campaign. Censored queries remain UNKNOWN.'})


def run(family):
    protocol=cert.read(D/'protocol.json')
    if protocol['sources']!=sources():raise ArithmeticError('exact-minimum sources changed')
    path=PARENT/family/'generic-census.json'
    if cert.hashed(path)!=protocol['census_hashes'][family]:raise ArithmeticError('census changed')
    census=cert.read(path);f=next(r for r in cert.read(ATLAS)['families'] if r['fibration_id']==family)
    G=matrix(ZZ,[[2*cert.F(v) for v in r] for r in f['generic_height_gram']])
    if [[str(cert.F(int(v),2)) for v in row] for row in G]!=[[str(cert.F(v)) for v in row] for row in f['generic_height_gram']]:
        raise ArithmeticError('scaled integral Gram required')
    U=G.LLL_gram();V=U.inverse()
    if any(q.denominator()!=1 for q in V.list()):raise ArithmeticError('basis transform not unimodular')
    V=matrix(ZZ,V);H=U.transpose()*G*U
    if U*V!=matrix.identity(ZZ,16) or V*U!=matrix.identity(ZZ,16):raise ArithmeticError('unimodular identities failed')
    convert=lambda M:[list(map(int,r)) for r in M]
    chosen=[r for r in census['records'] if 2*cert.F(r['norm'])==23]
    if len(chosen)!=protocol['query_counts'][protocol['families'].index(family)]:raise ArithmeticError('query roster changed')
    original_used={r['mask'] for r in census['selected']}
    output=D/family/'result.json'
    if output.exists():raise FileExistsError('preserve minimum audit')
    data={'schema':'elliptic-curves.mw16-exact-queried-parity-minima.v1','status':'RUNNING','protocol_hash':digest(protocol),
          'family':family,'gram':convert(G),'transform':convert(U),'inverse_transform':convert(V),'reduced_gram':convert(H),'rows':[]}
    checkpoint(output,data)
    for row in chosen:
        mask=row['mask'];parity=[(mask>>j)&1 for j in range(16)]
        residue=[sum(int(V[i,j])*parity[j] for j in range(16))%2 for i in range(16)]
        enumeration=enumerate_coset(convert(H),residue,22)
        minimum=enumeration['minimum_within_bound']
        if minimum is None:minimum=23;witness=row['representative']
        else:witness=[sum(int(U[i,j])*enumeration['minimum_witness'][j] for j in range(16)) for i in range(16)]
        if any((witness[i]-parity[i])%2 for i in range(16)) or sum(witness[i]*int(G[i,j])*witness[j] for i in range(16) for j in range(16))!=minimum:
            raise ArithmeticError('exact minimum witness failed')
        data['rows'].append({'mask':mask,'used_in_initial43':mask in original_used,'reported_census_norm':row['norm'],'exact_minimum':minimum,'minimum_witness':witness,'transformed_parity':residue,'enumeration':enumeration})
        checkpoint(output,data);print('EXACT PARITY MINIMUM',family,mask,minimum,flush=True)
    data['status']='COMPLETE_DECLARED_EXACT_AUDIT';checkpoint(output,data)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('stage',choices=['prepare','run']);p.add_argument('--family');a=p.parse_args();prepare() if a.stage=='prepare' else run(a.family)
