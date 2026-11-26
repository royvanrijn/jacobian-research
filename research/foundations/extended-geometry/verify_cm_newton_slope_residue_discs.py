"""Check the characteristic-zero and residue orders of two slope minors."""

from collections import defaultdict
from math import comb, lcm

from sage.all import GF, QQ, ZZ, PolynomialRing, TermOrder, matrix

from verify_cm_newton_slope_eliminant import (
    homrows,
    reduced_residual_rows,
    weighted_exponents,
)


def compressed_setup(rows):
    fs=homrows(rows)
    R=PolynomialRing(QQ,names=("r","h"),order="degrevlex")
    leadmap=lambda f:R({tuple(e[:2]): QQ(c) for e,c in f.dict().items() if not any(e[2:])})
    leading=R.ideal([leadmap(fs[i]) for i in (0,2)])
    assert leading.vector_space_dimension()==12
    G=leading.groebner_basis()
    L=PolynomialRing(QQ,names=("r","h","w","x","t","s"),order=TermOrder("degrevlex",2)+TermOrder("degrevlex",4))
    phi=lambda f:L({tuple(e)+(j,): c for e,p in f.dict().items() for j,c in p.dict().items()})
    psi=R.hom(L.gens()[:2],L)
    fs=[phi(f) for f in fs]
    lifted=[]
    for g in G:
        u,v=g.lift(leading)
        assert not (g-u*leading.gen(0)-v*leading.gen(1))
        # retain homogeneous coefficient degrees in r,h.
        u=R({e:c for e,c in u.dict().items() if sum(e)==g.total_degree()-3})
        v=R({e:c for e,c in v.dict().items() if sum(e)==g.total_degree()-4})
        assert not (g-u*leading.gen(0)-v*leading.gen(1))
        w=psi(u)*fs[0]+psi(v)*fs[2]
        assert w.lm()==psi(g.lm()) and w.lc()==1
        lifted.append(w)
    assert all(c.denominator()%331 for f in lifted for c in f.coefficients())
    s0=PolynomialRing(QQ,"s"); s=s0.gen()
    lms=[tuple(g.lm().exponents()[0]) for g in lifted]
    divis=lambda e:any(all(i>=j for i,j in zip(e,g[:5])) for g in lms)
    mon9=[e for e in weighted_exponents((1,)*5,9) if sum(e)==9 and not divis(e)]
    assert len(mon9)==394
    pos={e:i for i,e in enumerate(mon9)}
    coeffmat=[]
    for i,f in enumerate(fs):
        if i in (0,2): continue
        d=(i+6)//2
        cache={}
        for a in weighted_exponents((1,)*5,9-d):
            if sum(a)!=9-d or divis(a): continue
            rh=tuple(a[:2]); rest=tuple(a[2:])
            if rh not in cache:
                cache[rh]=(f*L({rh+(0,)*4: 1})).reduce(lifted)
            g=cache[rh]
            if not g: continue
            terms=defaultdict(dict)
            for (r,h,w,x,t,z),v in g.dict().items():
                j=pos[r,h,w+rest[0],x+rest[1],t+rest[2]]
                terms[j][z]=v
            coeffmat.append({j:s0(k) for j,k in terms.items()})
    return coeffmat,L,lifted

def series_mats(coef,alpha,maxk,field):
    out=[matrix(field,len(coef),394,sparse=False) for k in range(maxk+1)]
    for i,row in enumerate(coef):
        for j,c in row.items():
            tc=tuple(sum(field(a)*comb(k,l)*field(alpha)**(k-l) for k,a in c.dict().items() if k>=l) for l in range(maxk+1))
            for k,x in enumerate(tc):
                if x:out[k][i,j]=x
    return out

def block(m,I,J):
    return m.matrix_from_rows_and_columns(I,J)

def local_pivots(mats,r_generic):
    ranks=[]; IP=[]; JP=[]; ir=list(range(mats[0].nrows())); jr=list(range(mats[0].ncols()))
    while mats:
        I=mats[0].transpose().pivots(); J=mats[0].pivots(); rank=len(I); assert rank==len(J)
        ranks.append(rank)
        if rank==0:
            mats.pop(0); continue
        IP.extend(ir[i] for i in I); JP.extend(jr[j] for j in J)
        ni=[i for i in range(mats[0].nrows()) if i not in I];nj=[j for j in range(mats[0].ncols()) if j not in J]
        ir=[ir[i] for i in ni]; jr=[jr[j] for j in nj]
        if sum(ranks)==r_generic: break
        Us=[block(m,I,J) for m in mats]; Vs=[block(m,I,nj) for m in mats]; Ws=[block(m,ni,J) for m in mats]; Zs=[block(m,ni,nj) for m in mats]
        U0inv=Us[0].inverse(); Ys=[]
        for k in range(len(mats)):
            Ys.append(U0inv*(Vs[k]-sum(Us[j]*Ys[k-j] for j in range(1,k+1))))
        Bs=[Zs[k]-sum(Ws[j]*Ys[k-j] for j in range(k+1)) for k in range(len(mats))]
        assert not Bs[0]
        mats=Bs[1:]
    return ranks,IP,JP


def solve_series(pivot, right, precision):
    """Solve a matrix equation over the Taylor ring to the requested order."""
    positive = [coefficient.change_ring(QQ) for coefficient in pivot[1:precision+1]]
    solution = []
    for k in range(precision + 1):
        residual = right[k] - sum(positive[j-1]*solution[k-j] for j in range(1, k+1))
        denominator = ZZ(residual.denominator())
        integral = (denominator*residual).change_ring(ZZ)
        answer = pivot[0]._solve_right_nonsingular_square(
            integral, check_rank=False, algorithm="flint",
        ) / denominator
        assert pivot[0]*answer == residual
        solution.append(answer)
    return solution


def exact_schur(coefficients, rows, columns, rank, alpha, order):
    assert len(rows) == len(columns) == 364
    selected = {j: k for k, j in enumerate(columns)}
    matrices = [matrix(ZZ, len(rows), len(columns)) for _ in range(order + 1)]
    for i, index in enumerate(rows):
        series = {
            selected[j]: [
                QQ(sum(a*comb(k, q)*alpha**(k-q) for k, a in c.dict().items() if k >= q))
                for q in range(order + 1)
            ]
            for j, c in coefficients[index].items() if j in selected
        }
        denominator = lcm(*(int(value.denominator()) for terms in series.values() for value in terms))
        assert denominator % 331
        for j, terms in series.items():
            for q, value in enumerate(terms):
                if value:
                    matrices[q][i, j] = ZZ(denominator*value)

    extra = len(rows) - rank
    A = [M.submatrix(0, 0, rank, rank) for M in matrices]
    B = [M.submatrix(0, rank, rank, extra) for M in matrices]
    C = [M.submatrix(rank, 0, extra, rank) for M in matrices]
    D = [M.submatrix(rank, rank, extra, extra) for M in matrices]
    assert A[0].change_ring(GF(331)).det()
    # C A^(-1) B agrees with C X + Y(B-A X) to twice the
    # precision of the two simpler linear solves.
    steps = order//2
    X = solve_series(A, B, steps)
    print(f"Slope {alpha}: right jet of order {steps} checked", flush=True)
    Y = [Z.transpose() for Z in solve_series(
        [M.transpose() for M in A[:steps+1]],
        [M.transpose() for M in C[:steps+1]], steps,
    )]
    residuals = [
        B[k] - sum(A[k-j]*X[j] for j in range(min(k, steps)+1))
        for k in range(steps+1, order+1)
    ]
    return [
        D[k] - sum(C[k-j]*X[j] for j in range(min(k, steps)+1))
        - sum(Y[j]*residuals[k-j-steps-1] for j in range(max(0, k-order), min(steps, k-steps-1)+1))
        for k in range(order+1)
    ]


def verify_residue_discs(rows):
    coefficients, _, _ = compressed_setup(rows)
    assert len(coefficients) == 492
    cases = [
        (QQ(17)/2, [359, 2, 0, 0, 0, 2, 0, 1], 19),
        (QQ(56)/5, [361, 1, 0, 1, 0, 0, 1], 10),
    ]
    for alpha, expected, vanishing_order in cases:
        order = len(expected)-1
        modular = series_mats(coefficients, alpha, order, GF(331))
        levels, row_indices, column_indices = local_pivots(modular, 364)
        assert levels == expected
        assert sum(i*count for i, count in enumerate(levels)) == vanishing_order
        exact = exact_schur(coefficients, row_indices, column_indices, levels[0], alpha, order)
        exact_levels, _, _ = local_pivots(exact, 364-levels[0])
        assert exact_levels == [0, *expected[1:]]
        print(f"Slope {alpha}: exact and residue orders both {vanishing_order}", flush=True)


if __name__ == "__main__":
    verify_residue_discs(reduced_residual_rows(5, 16, even_linear=True, primitive=False))
