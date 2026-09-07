#!/usr/bin/env sage-python
"""Eighty-eight fixed MW14 triangles anchored at P6; exact frames and inverses.

Prepare freezes 88 explicit words and their rational pencils before target
tests. Probe uses exact generic cubic conversion at eight frozen primes.
One worker, 600 seconds per command; one checkpoint per completed model.
A residue survivor remains UNKNOWN.
The roster is a finite construction, not a complete fibration classification.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import runpy
import signal
from sage.all import QQ, ZZ, GF, EllipticCurve, matrix, vector, block_diagonal_matrix, pari
from sage.env import SAGE_VERSION

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT/'artifacts/generated-results'
SOURCE = ART/'elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json'
COMPILER = ROOT/'elkies-k3/scripts/search_curve302_shortword_triangles.sage'
RANK = ROOT/'elkies-k3/scripts/certify_curve302_triangle_mw14.sage'
CONVERT = ROOT/'elkies-k3/scripts/certify_curve302_triangle_generic_conversion.sage'
PUBLIC = ROOT/'elliptic-curves/cas/icarm_curve302.py'
PROTOCOL = ART/'elkies-k3-curve302-anchor6-mw14-protocol-v1.json'
OUT = ART/'elkies-k3-curve302-anchor6-mw14-v1.json'
ROSTER = ART/'elkies-k3-curve302-anchor6-mw14-roster-v1.json'
LOCAL = ROOT/'artifacts/local/elkies-k3/curve302-anchor6-mw14'
PRIMES = [1013,1021,1009,1031,1033,1039,1049,1051]


def digest(p): return sha256(p.read_bytes()).hexdigest()
def rows(m): return [list(map(str,r)) for r in m.rows()]

def prepare():
    source = json.loads(SOURCE.read_text())
    G = matrix(ZZ,source['sections']['height_gram'])
    N = block_diagonal_matrix(matrix(ZZ,[[-2,1],[1,0]]),-G)
    e = matrix.identity(ZZ,17); P = e[5]
    def sec(w): return vector(ZZ,[1,(w*G*w)//2]+list(w))
    O = sec(vector(ZZ,17)); compiler = runpy.run_path(str(COMPILER))
    exact_roots = runpy.run_path(str(RANK))['exact_roots']; models = []
    roster = json.loads(ROSTER.read_text())['records']
    assert len(roster)==88
    for record in roster:
        index=record['index']; word=record['Q_word']
        assert record['P_word']==list(P)
        Q = vector(ZZ,word); z = vector(ZZ,record['Z_word'])
        D = O+sec(P)+sec(Q); Z = sec(z)
        assert P*G*P == Q*G*Q == 6 and P*G*Q == 3
        assert D*N*D == 0 and D*N*Z == 1 and Z*N*Z == -2
        C = matrix(ZZ,[D*N,Z*N]).right_kernel_matrix()
        U = matrix(ZZ,pari(-C*N*C.transpose()).qflllgram())
        assert abs(U.det()) == 1
        C = U.transpose()*C; H = -C*N*C.transpose()
        assert H.det() == 948 and H.is_positive_definite()
        data = pari(H).qfminim(2,10000,2); V = matrix(ZZ,data[2]).transpose()
        assert int(data[0]) == 2*V.nrows() == 8 and V.rank() == 3
        exact,nodes = exact_roots(H)
        assert exact == {tuple(v) for r in V.rows() for v in (r,-r)}
        B = V.row_module().basis_matrix(); RG = B*H*B.transpose()
        assert RG.det() == 6
        smith,left,right = B.smith_form()
        assert list(smith.diagonal()) == [1,1,1]
        complement = right.inverse()[3:,:]
        projection = complement-complement*H*B.transpose()*RG.inverse()*B
        MW = projection*H*projection.transpose()
        assert MW.nrows() == 14 and MW.det() == 158
        m = dict(index=index,P_word=list(map(int,P)),Q_word=word,Z_word=list(map(int,z)))
        K,f,c,zmap,meet = compiler['make_pencil'](QQ,m); rat = compiler['rat']
        m.update(frame_basis=rows(C),frame_gram=rows(H),root_vectors=rows(V),
                 fibre_class=list(map(str,D)),zero_class=list(map(str,Z)),
                 generic_rank=14,root_type='A2+A1',torsion_order=1,
                 abstract_MW_gram=rows(MW),abstract_MW_determinant=158,
                 exact_root_nodes=nodes,
                 raw_cubic_coefficients=[[rat(a) for a in f[i].list()] for i in range(4)],
                 pencil_coefficients=list(map(rat,c)),zero_parameter=rat(zmap),
                 intersection_parameters=list(map(str,meet)))
        models.append(m); print('PREPARED_MW14',index,flush=True)
    return dict(schema='curve302.anchor6-mw14.protocol.v1',
                status='88_EXACT_MW14_PENCILS_PREPARED',sage_version=SAGE_VERSION,
                input_sha256={str(p.relative_to(ROOT)):digest(p) for p in [Path(__file__),ROSTER,SOURCE,COMPILER,RANK,CONVERT,PUBLIC]},
                models=models,primes=PRIMES,
                limits=dict(workers=1,seconds_per_command=600),
                boundary='Eighty-eight fixed triangle configurations on X948. No completeness or pairwise inequivalence claim. Full abstract generic MW lattices, but no Weierstrass-coordinate MW bases or302 specialization yet.')

def probe():
    protocol = json.loads(PROTOCOL.read_text())
    for p,h in protocol['input_sha256'].items(): assert digest(ROOT/p)==h,p
    convert = runpy.run_path(str(CONVERT))['convert']
    target = EllipticCurve(QQ,runpy.run_path(str(PUBLIC))['GENERAL_WEIERSTRASS_COEFFICIENTS'])
    records = []
    LOCAL.mkdir(parents=True,exist_ok=True)
    for model in protocol['models']:
        results = []
        for p in PRIMES:
            field = GF(p)
            try:
                assert field(target.discriminant())
                data = convert(field,model); j = data['j']; n,d = j.numerator(),j.denominator()
                f = n-field(target.j_invariant())*d
                roots = [int(a) for a in field if f(a)==0]; infinity = f[24]==0
                result = dict(prime=p,j_numerator=list(map(int,n.list())),j_denominator=list(map(int,d.list())),
                              comparison=list(map(int,f.list())),finite_roots=roots,infinity_root=bool(infinity),
                              status='EXCLUDED_ALL_RATIONAL_PARAMETERS' if not roots and not infinity else 'MODULAR_SURVIVOR')
            except (AssertionError,ArithmeticError,ValueError,NotImplementedError) as error:
                result = dict(prime=p,status='UNKNOWN',reason=type(error).__name__+': '+str(error))
            results.append(result)
            print('INVERSE',model['index'],p,result['status'],result.get('finite_roots'),flush=True)
            if result['status']=='EXCLUDED_ALL_RATIONAL_PARAMETERS': break
        status = 'EXCLUDED_ALL_RATIONAL_PARAMETERS' if any(r['status']=='EXCLUDED_ALL_RATIONAL_PARAMETERS' for r in results) else 'UNKNOWN'
        record=dict(index=model['index'],status=status,results=results)
        records.append(record)
        checkpoint=LOCAL/('row'+str(model['index'])+'.json')
        if checkpoint.exists(): assert record==json.loads(checkpoint.read_text())
        else: checkpoint.write_text(json.dumps(record,indent=2)+'\n')
    return dict(schema='curve302.anchor6-mw14.inverse.v1',protocol_sha256=digest(PROTOCOL),records=records,
                excluded_count=sum(r['status']=='EXCLUDED_ALL_RATIONAL_PARAMETERS' for r in records),
                boundary='Exact all-rational exclusions apply only to certified degree-preserving reductions of these 88 pencils. Survivors are UNKNOWN; no rational parameter is inferred from residues.')

if __name__=='__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('command',choices=['prepare','probe']); parser.add_argument('--check',action='store_true'); args=parser.parse_args()
    signal.alarm(600)
    result = prepare() if args.command=='prepare' else probe()
    path = PROTOCOL if args.command=='prepare' else OUT
    if args.check: assert result==json.loads(path.read_text())
    else:
        assert not path.exists(),'Preserve existing evidence'
        path.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    print('PASS',args.command,flush=True)
