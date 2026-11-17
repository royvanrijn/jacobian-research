"""PF2RULED31 discovery: exact finite reduction and rational Dixon inversion.
Run with sage -python from the repository. Temporary dense matrices are deleted;
retained polynomial certificates are checked by the separate FLINT replay.
"""
from pathlib import Path
from time import monotonic
import gzip,json,itertools,subprocess,tempfile
from sage.all import QQ,NumberField,PolynomialRing,matrix,prod
from sage.libs.singular.function_factory import singular_function
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]

def main():
    start=monotonic();primary=json.loads(gzip.decompress((HERE/'degree8_quartic_descent_certificates.json.gz').read_bytes()))['sectors'];sectors=[]
    for old in primary:
        k=old['multiplicities'][0];T=PolynomialRing(QQ,'t');P4=T(old['P4']);K=NumberField(P4,'t');t=K.gen()
        A5=PolynomialRing(K,names=old['core_variables']);V5=PolynomialRing(A5,'w');A=PolynomialRing(K,names=old['core_variables'][:4],order='degrevlex');V=PolynomialRing(A,'w');inject=A5.hom([*A.gens(),0],A)
        transfer=lambda text:V([inject(c) for c in V5(text)])
        W,M,L=[transfer(old[key]) for key in ['normalized_W','normalized_M','normalized_L']];s=K(old['normalized_cusp_parameter']);inputs=[M[j] for j in range(2,6)];I=A.ideal(inputs)
        G=list(I.groebner_basis(algorithm='libsingular:slimgb'));assert len(G)==35
        lift=singular_function('lift');witness=lift(I,A.ideal(G));assert matrix(A,1,4,inputs)*witness==matrix(A,1,35,G)
        Q=A.quotient(I,names=list(A.variable_names()));VQ=PolynomialRing(Q,'w');w=VQ.gen();qpoly=lambda f:VQ([Q(c) for c in f]);WW,MM,LL=map(qpoly,[W,M,L]);wh,rem=WW.quo_rem(w-Q(s));assert rem==0
        e5=A(((MM*wh).quo_rem(LL)[1][5]).lift())
        exponents=sorted(e for e in itertools.product(range(4),repeat=4) if sum(e)<=3);basis=[prod(q**a for q,a in zip(A.gens(),e)) for e in exponents];columns=[]
        for i,b in enumerate(basis):
            nf=A(Q(e5*b).lift());columns.append([nf.monomial_coefficient(bb) for bb in basis])
            if i%10==0:print(f'k{k} exact multiplication column{i}/35; elapsed{monotonic()-start:.1f}s',flush=True)
        rational=matrix(QQ,140,140)
        for i in range(35):
            for j in range(35):
                for pw in range(4):
                    c=columns[j][i]*t**pw
                    for pp in range(4):rational[4*i+pp,4*j+pw]=c[pp]
        rhs=['0']*140;rhs[4*exponents.index((0,0,0,0))]='1'
        with tempfile.TemporaryDirectory() as temp:
            mp=Path(temp)/'matrix.json.gz';sp=Path(temp)/'solution.json'
            mp.write_bytes(gzip.compress(json.dumps(dict(rows=[[str(c) for c in row] for row in rational.rows()],rhs=rhs)).encode(),mtime=0))
            subprocess.run([str(ROOT/'.venv/bin/python'),str(HERE/'compute_degree8_linear_remainder_dixon.py'),str(mp),str(sp)],check=True,timeout=30)
            solution=json.loads(sp.read_text())
        def encode(f):return [[list(e),[str(c[j]) for j in range(4)]] for e,c in A(f).dict().items()]
        sectors.append(dict(sector=k,P4=str(P4),P4_coefficients=[str(P4[j]) for j in range(5)],variables=list(A.variable_names()),
            generators=[encode(g) for g in inputs],rules=[encode(g) for g in G],rules_in_inputs=[[encode(witness[i,j]) for i in range(4)] for j in range(35)],
            W=[encode(c) for c in W],M=[encode(c) for c in M],L=[encode(c) for c in L],s=[str(s[j]) for j in range(4)],reduced_E5=encode(e5),basis_exponents=exponents,
            inverse=[[list(ex),solution[4*i:4*i+4]] for i,ex in enumerate(exponents)]))
        print(f'k{k} complete exact polynomial checkpoint; elapsed{monotonic()-start:.1f}s',flush=True)
    out=HERE/'degree8_linear_remainder_certificates.json.gz'
    out.write_bytes(gzip.compress((json.dumps(dict(as_of='2026-11-14',sectors=sectors),indent=2)+'\n').encode(),mtime=0))
    subprocess.run([str(ROOT/'.venv/bin/python'),str(HERE/'verify_degree8_linear_remainder.py')],check=True,timeout=60)
if __name__=='__main__':main()
