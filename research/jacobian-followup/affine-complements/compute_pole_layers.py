"""Exact bounded linear algebra for pole orders (2,2) and (2,3).

S=Bl_[0:0:1] P2, B is the strict transform of XZ-Y^2=0.
No Plucker/bilinear feasibility or finite-map existence is asserted.
"""
from pathlib import Path
from itertools import combinations
from time import monotonic
import hashlib,json,platform
import sympy as S

HERE=Path(__file__).resolve().parent
x,y,z=S.symbols("x y z")
conic=x-y*y

def basis(m):
    # O(mB): degree 2m homogeneous forms with multiplicity at least m at p.
    monomials=[x**i*y**j for degree in range(m,2*m+1)
               for i in range(degree+1) for j in [degree-i]]
    # b^m has coefficient 1 at x^m. Replace that monomial by b^m (the constant).
    monomials.remove(x**m)
    return monomials

def wedge_numerator(A,C,m,n):
    return S.Poly(S.expand(
        conic*(S.diff(A,x)*S.diff(C,y)-S.diff(A,y)*S.diff(C,x))
        -n*C*(S.diff(A,x)*S.diff(conic,y)-S.diff(A,y)*S.diff(conic,x))
        -m*A*(S.diff(conic,x)*S.diff(C,y)-S.diff(conic,y)*S.diff(C,x))),x,y)

def layer(m,n):
    start=monotonic()
    bm,bn=basis(m),basis(n)
    pairs=list(combinations(range(len(bm)),2)) if m==n else [
        (i,j) for i in range(len(bm)) for j in range(len(bn))]
    cols=[]
    for i,j in pairs:
        W=wedge_numerator(bm[i],bn[j],m,n).as_expr()
        # conic=z and x=z+y^2: pole <=3 iff z^(m+n-2) divides W.
        shifted=S.Poly(S.expand(W.subs(x,z+y*y)),z,y)
        cols.append({power:coeff for power,coeff in shifted.terms()
                     if power[0] < m+n-2})
    rows=sorted({power for col in cols for power in col})
    M=S.polys.matrices.DomainMatrix.from_Matrix(
        S.MutableSparseMatrix(len(rows),len(cols),{
            (i,j):cols[j].get(power,0)
            for j in range(len(cols)) for i,power in enumerate(rows)
            if cols[j].get(power,0)})).convert_to(S.QQ)
    rref,pivots=M.rref()
    kernel=M.nullspace().to_Matrix()
    dense=M.to_Matrix()
    assert dense*kernel.T==S.zeros(len(rows),kernel.rows)
    assert len(pivots)+kernel.rows==len(cols)
    # Rank also certified by a full-column-rank matrix made from pivot columns.
    assert dense[:,list(pivots)].rank()==len(pivots)
    return dict(pole_orders=[m,n],basis_m=list(map(str,bm)),basis_n=list(map(str,bn)),
                section_dimensions=[len(bm)+1,len(bn)+1],
                columns=pairs,row_monomials=[list(r) for r in rows],
                sparse_matrix=[[i,j,str(value)] for (i,j),value in
                               sorted(dense.todok().items())],
                constraint_rank=len(pivots),kernel_dimension=kernel.rows,
                pivot_columns=list(pivots),
                kernel_basis=[[str(v) for v in kernel.row(i)] for i in range(kernel.rows)],
                elapsed_seconds=round(monotonic()-start,6))

def main():
    layers=[layer(2,2),layer(2,3)]
    out=dict(as_of="2026-11-14",outcome="PASS",layers=layers,
             scope="Exact necessary differential-pole linear constraints; no decomposability, degree-four, ramification-support or finiteness solution",
             checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
             python=platform.python_version(),sympy=S.__version__)
    (HERE/"pole_layers.json").write_text(json.dumps(out,indent=2)+"\n")
    for r in layers:
        print("pole orders",r["pole_orders"],"sections",r["section_dimensions"],
              "columns",len(r["columns"]),"rank",r["constraint_rank"],
              "kernel",r["kernel_dimension"],"seconds",r["elapsed_seconds"])

if __name__=="__main__":
    main()

