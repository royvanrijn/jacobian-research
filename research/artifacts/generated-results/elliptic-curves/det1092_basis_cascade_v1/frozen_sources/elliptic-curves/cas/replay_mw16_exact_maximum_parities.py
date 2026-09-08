#!/usr/bin/env python3
"""Replay queried minima and every generic parity upper bound without Sage."""
import argparse
from pathlib import Path
import certify_compact_r17_candidates as cert
import compact_mw16_specialization as spec
from exact_parity_ellipsoid import enumerate_coset,ldl
from research_runtime.store import checkpoint,digest
ROOT=Path(__file__).resolve().parents[2]
PARENT=ROOT/'artifacts/local/elliptic-curves/prospective-mw16-h4096-v1'


def multiply(A,B):
    return [[sum(A[i][k]*B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]


def replay(path,output):
    if output.exists():raise FileExistsError('preserve exact-minimum replay')
    data=cert.read(path);protocol_path=path.parents[1]/'protocol.json';protocol=cert.read(protocol_path)
    for name,h in protocol['sources'].items():
        if cert.hashed(ROOT/name)!=h:raise ArithmeticError('frozen enumeration source differs')
    if data['status']!='COMPLETE_DECLARED_EXACT_AUDIT' or data['protocol_hash']!=digest(protocol):raise ArithmeticError('incomplete enumeration')
    family=data['family'];f=next(r for r in cert.read(spec.ATLAS)['families'] if r['fibration_id']==family)
    G,U,V,H=(data[k] for k in ('gram','transform','inverse_transform','reduced_gram'))
    for M in (G,U,V,H):
        if len(M)!=16 or any(len(r)!=16 or any(type(v) is not int for v in r) for r in M):raise ArithmeticError('integer matrix dimensions differ')
    if G!=[[2*cert.F(v) for v in r] for r in f['generic_height_gram']]:raise ArithmeticError('scaled integral generic Gram differs')
    I=[[int(i==j) for j in range(16)] for i in range(16)]
    if multiply(U,V)!=I or multiply(V,U)!=I or multiply(list(zip(*U)),multiply(G,U))!=H:raise ArithmeticError('unimodular Gram identities failed')
    ldl(H)
    census_path=PARENT/family/'generic-census.json';census=cert.read(census_path)
    if cert.hashed(census_path)!=protocol['census_hashes'][family] or [[2*cert.F(v) for v in r] for r in census['gram']]!=G or len(census['records'])!=65536:raise ArithmeticError('complete census binding failed')
    queried=[]
    for mask,row in enumerate(census['records']):
        z=row['representative']
        if row['mask']!=mask or len(z)!=16 or any(type(v) is not int for v in z) or any((z[j]-(mask>>j))%2 for j in range(16)):raise ArithmeticError('census parity witness differs')
        nonzero=[(i,v) for i,v in enumerate(z) if v]
        norm=sum(v*G[i][j]*w for i,v in nonzero for j,w in nonzero)
        if norm!=2*cert.F(row['norm']) or norm<0 or norm>23:raise ArithmeticError('generic parity upper bound failed')
        if norm==23:queried.append(mask)
    if queried!=[r['mask'] for r in data['rows']]:raise ArithmeticError('maximum-candidate roster differs')
    selected={r['mask'] for r in census['selected']}
    for row in data['rows']:
        mask=row['mask'];p=[(mask>>i)&1 for i in range(16)]
        residue=[sum(V[i][j]*p[j] for j in range(16))%2 for i in range(16)]
        if residue!=row['transformed_parity'] or row['used_in_initial43']!=(mask in selected):raise ArithmeticError('parity transport differs')
        enumeration=enumerate_coset(H,residue,22)
        if enumeration!=row['enumeration']:raise ArithmeticError('complete rational ellipsoid enumeration differs')
        minimum=enumeration['minimum_within_bound']
        if minimum is None:minimum=23
        witness=row['minimum_witness']
        if len(witness)!=16 or any(type(v) is not int for v in witness) or any((witness[j]-p[j])%2 for j in range(16)):raise ArithmeticError('minimum witness parity failed')
        if minimum!=row['exact_minimum'] or sum(witness[i]*G[i][j]*witness[j] for i in range(16) for j in range(16))!=minimum:raise ArithmeticError('minimum witness norm failed')
    count=sum(r['exact_minimum']==23 for r in data['rows'])
    checkpoint(output,{'schema':'elliptic-curves.mw16-exact-parity-radius-replay.v1','status':'PASS','family':family,
        'sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in (Path(__file__).resolve(),ROOT/'elliptic-curves/cas/exact_parity_ellipsoid.py',path,protocol_path,census_path,spec.ATLAS)},
        'parity_upper_bounds_checked':65536,'queried_minima_checked':len(data['rows']),
        'exact_maximum_parity_coset_minimum':'23/2' if count else None,'classes_with_exact_minimum23_over2':count,
        'omitted_classes_with_exact_minimum23_over2':sum(r['exact_minimum']==23 and not r['used_in_initial43'] for r in data['rows']),
        'claim_boundary':'Every parity class has a checked norm<=23/2 witness. Every candidate maximum class is checked by exhaustive rational LDL enumeration of scaled integral norm<=22. A surviving norm23 witness proves discrete maximum23/2 on Z^16/2Z^16. No continuous covering radius, saturation or specialized rank.'})
    print('REPLAYED EXACT DISCRETE PARITY RADIUS',family,'norm23/2 classes',count,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--input',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();replay(a.input.resolve(),a.output.resolve())
