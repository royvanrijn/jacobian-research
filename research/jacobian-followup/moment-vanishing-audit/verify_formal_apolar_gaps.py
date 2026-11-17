#!/usr/bin/env python3
"""Bounded exact controls for written all-formal/ring reductions; no gap solver."""
from pathlib import Path
from itertools import product
import hashlib,json
import sympy as S

HERE=Path(__file__).resolve().parent
x,y,z,X,Y,Z,e=S.symbols('x y z X Y Z e')
xs=(x,y,z);zs=(X,Y,Z)


def action(symbol,f):
    assert not S.sympify(symbol).has(S.Float) and not S.sympify(f).has(S.Float)
    f=S.expand(f);assert S.Poly(f,*xs).total_degree()<=32
    result=0
    for powers,coefficient in S.Poly(S.expand(symbol),*zs).terms():
        term=f
        for variable,k in zip(xs,powers):term=S.diff(term,variable,k)
        result+=coefficient*term
    return S.expand(result)


def iterate(symbol,f,k):
    assert isinstance(k,int) and 0<=k<=16
    for _ in range(k):f=action(symbol,f)
    return f


def ad_value(symbol,Q,f,j):
    return S.expand(sum((-1)**i*S.binomial(j,i)*iterate(symbol,Q*iterate(symbol,f,i),j-i)
                        for i in range(j+1)))


def height(symbol,phase,m):
    current=S.expand(phase**m);assert current!=0
    answer=0
    for k in range(1,17):
        current=action(symbol,current)
        if current==0:return answer
        answer=k
    raise AssertionError('nilpotence exceeds bounded control')


def apolar_dimension(f):
    degree=S.Poly(f,x,y).total_degree()
    mons=[x**i*y**j for i in range(degree+1) for j in range(degree+1-i)]
    derivatives=[S.diff(f,x,i,y,j) for i in range(degree+1) for j in range(degree+1-i)]
    matrix=S.Matrix([[S.expand(d).coeff(x,i).coeff(y,j)
                       for d in derivatives]
                      for i in range(degree+1) for j in range(degree+1-i)])
    return matrix.rank()


def main():
    comm=[]
    models=[(X-Y**2-Y**6-Y**24,x**2*y,x**3+y**2+z),
            (X**2+Y*Z,x*y+z**2,x**2*y*z+y**4)]
    for symbol,Q,f in models:
        delta=S.Poly(Q,*xs).total_degree()
        first=action(symbol,Q*f)-Q*action(symbol,f)
        taylor=0
        for gamma in product(range(delta+1),repeat=3):
            if not 0<sum(gamma)<=delta:continue
            q=Q;lam=symbol
            for v,w,k in zip(xs,zs,gamma):
                q=S.diff(q,v,k);lam=S.diff(lam,w,k)/S.factorial(k)
            taylor+=q*action(lam,f)
        assert S.expand(first-taylor)==0
        assert ad_value(symbol,Q,f,delta+1)==0
        for m in (delta,delta+1):
            lhs=iterate(symbol,Q*f,m)
            rhs=sum(S.binomial(m,j)*ad_value(symbol,Q,iterate(symbol,f,m-j),j)
                    for j in range(delta+1))
            assert S.expand(lhs-rhs)==0
            comm.append(dict(symbol=str(symbol),marker=str(Q),m=m,
                             taylor_equal=True,binomial_equal=True,ad_nilpotent=True))
    lemma=[]
    for symbol,f,k in [(X**2,x**3*y**2,2),(X-Y**2-Y**6-Y**24,z**3,1)]:
        assert iterate(symbol,f,k)==0
        for Q in (x,y**2,x*y*z):
            delta=S.Poly(Q,*xs).total_degree();m=k+delta
            assert iterate(symbol,Q*f,m)==0
            lemma.append(dict(marker=str(Q),nilpotent_power=k,m=m,zero=True))
    heights=[]
    for m in range(1,7):
        for tag,symbol,phase,expected in [
            ('unbalanced_safe',X**2,x*y**3,m//2),
            ('unit_multiple',X**2*(1+Y+Y**2),x*y**3,m//2),
            ('nonalgebraic_graph_jet',X-Y**2-Y**6-Y**24,z**3,0),
            ('purity_required_boundary',X,x,m),
        ]:
            h=height(symbol,phase,m);assert h==expected
            heights.append(dict(model=tag,m=m,height=h,deficit=m-h,
                                pure_at_index=h<m))
    reduce_dual=lambda f:S.Poly(S.expand(f),e).nth(0)+e*S.Poly(S.expand(f),e).nth(1)
    assert reduce_dual(iterate(e+X,y**2,2))==0
    assert reduce_dual(iterate(e+X,x**2*y**2,4))==0
    assert reduce_dual(iterate(e+X,y**2,1))==e*y**2
    dims={'x2':apolar_dimension(x**2),'x2_plus_y':apolar_dimension(x**2+y),
          'x4':apolar_dimension(x**4),'x2_plus_y_squared':apolar_dimension((x**2+y)**2)}
    assert dims==dict(x2=3,x2_plus_y=3,x4=5,x2_plus_y_squared=6)
    assert action(Y-X**2/2,x**2+y)==0 and action(X**3,x**2+y)==0
    inequalities=[]
    for N,delta,m in [(0,1,2),(1,2,5),(2,3,10),(3,2,9)]:
        assert m>(N+1)*delta
        assert m-delta>N*delta
        inequalities.append(dict(N=N,delta=delta,m=m,k=m-delta))
    guards=0
    for symbol,f in [(X,S.Float(1)),(X,x**33)]:
        try:action(symbol,f)
        except AssertionError:guards+=1
        else:raise AssertionError('guard failed')
    try:iterate(X,x,17)
    except AssertionError:guards+=1
    else:raise AssertionError('guard failed')
    deps=['verify_formal_symbol_source_rigidity.py','verify_filtered_principal_transfer.py']
    report=dict(completed=True,checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        supporting_checker_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        sympy_version=S.__version__,formal_jet_order=24,commutator_controls=comm,
        annihilation_lemma_controls=lemma,nilpotence_height_controls=heights,
        square_zero_constant_symbol_control=True,apolar_dimensions=dims,
        first_apolar_algebra_relation_control=True,strict_index_controls=inequalities,
        guard_refusals=guards,exact_bounded_calculator_implemented=True,
        universal_reductions_proved_by_finite_controls=False,
        general_eventuality_to_linearity_decided=False,
        unbounded_deficit_implies_positive_gap_proved=False,
        arbitrary_deficit_sequence_realized=False,
        symmetric_power_apolar_identification_without_hypothesis=False,
        ordinary_cutoff_ideals_claimed_equal=False,
        general_ideal_or_gap_solver_implemented=False,
        full_exceptional_MZ_decided=False,
        independent_replay=False,external_review=False,formal_verification=False,novelty_claimed=False)
    target=HERE/'formal_apolar_gaps_verification.json';assert not target.exists()
    target.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(completed=True,commutators=len(comm),lemma=len(lemma),
                          heights=len(heights),apolar_dimensions=dims,guards=guards)))


if __name__=='__main__':main()
