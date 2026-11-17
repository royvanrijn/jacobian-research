#!/usr/bin/env python3
"""Small exact controls for the written all-order height/gap construction."""
from pathlib import Path
import hashlib,json
import sympy as S

HERE=Path(__file__).resolve().parent


def action(symbol,poly,variables,dual):
    assert not S.sympify(symbol).has(S.Float) and not S.sympify(poly).has(S.Float)
    poly=S.expand(poly);assert S.Poly(poly,*variables).total_degree()<=64
    result=0
    for powers,coefficient in S.Poly(S.expand(symbol),*dual).terms():
        term=poly
        for variable,k in zip(variables,powers):term=S.diff(term,variable,k)
        result+=coefficient*term
    return S.expand(result)


def iterate(symbol,poly,k,variables,dual):
    assert isinstance(k,int) and 0<=k<=12
    for _ in range(k):poly=action(symbol,poly,variables,dual)
    return poly


def height(symbol,phase,m,variables,dual):
    if m==0:return 0
    poly=S.expand(phase**m);assert poly!=0
    for k in range(1,13):
        poly=action(symbol,poly,variables,dual)
        if poly==0:return k-1
    raise AssertionError('bounded nilpotence height exceeded')


def main():
    x,y,X,Y=S.symbols('x y X Y');variables=(x,y);dual=(X,Y)
    convolution=[]
    for m in (1,2):
        h1=height(X,x**2,m,variables,dual)
        h2=height(Y**2,y**3,m,variables,dual)
        hs=height(X+Y**2,x**2*y**3,m,variables,dual)
        hp=height(X*Y**2,x**2*y**3,m,variables,dual)
        hq=height((X+Y**2)**2,x**2*y**3,m,variables,dual)
        assert hs==h1+h2 and hp==min(h1,h2) and hq==hs//2
        convolution.append(dict(m=m,heights=[h1,h2],sum=hs,product=hp,power2=hq))
    safe=[]
    for m in range(1,7):
        h=height(X**3+Y**3,x*y,m,variables,dual)
        assert h==2*(m//3) and h<m
        safe.append(dict(m=m,height=h,nilpotence_index=h+1))
    assert safe[2]['nilpotence_index']>safe[0]['nilpotence_index']+safe[1]['nilpotence_index']
    sharp=iterate(X**3+Y**3,x**3*(x*y)**3,3,variables,dual)
    assert sharp==12960
    a,b,c,d,t,A,B,C,D,T=S.symbols('a b c d t A B C D T')
    seed=(t+c)*(a*d+b*t);seed_symbol=T*(A*D-B*C)
    sv=(a,b,c,d,t);sd=(A,B,C,D,T)
    seed_controls=[]
    for m in (1,2,3):
        assert iterate(seed_symbol,seed**m,m,sv,sd)==0
        assert iterate(seed_symbol,seed**m,m-1,sv,sd)!=0
        seed_controls.append(dict(m=m,height=m-1,pure=True,penultimate_nonzero=True))
    u,v,w,U,V,W=S.symbols('u v w U V W');av=(u,v,w);ad=(U,V,W)
    auxiliary=[]
    for r,indices in [(5,[1,2,3,4,5,6]),(7,[1,4,7,8])]:
        phase=u*v*w**(r-2);symbol=U**r+V**r+W**r
        for m in indices:
            actual=height(symbol,phase,m,av,ad)
            expected=2*(m//r)+((r-2)*m)//r
            assert actual==expected
            s=m%r
            assert expected==m-(2*s+r-1)//r
            auxiliary.append(dict(r=r,m=m,height=actual))
    # The family is verified by exact factor heights, not an eight-variable sweep.
    periods=[]
    for M in (1,2,3):
        r=2*M+3
        h=lambda m:(2*m-1-(2*(m%r)+r-1)//r)//2
        for m in range(1,2*r+1):
            expected=m-1 if m%r<=M+1 else m-2
            assert h(m)==expected and h(m)<m
        assert all(h(k+1)<k for k in range(M+1,2*M+1))
        assert all(h(k+2)>=k for k in range(1,2*r+1))
        assert h(r)==r-1
        periods.append(dict(M=M,r=r,pure_indices_checked=2*r,
                            zero_first_shift_window=[M+1,2*M],
                            square_marker_nonzero_indices_checked=2*r))
    # Actual M=1 window at k=2: every tensor allocation is zero directly.
    phase=u*v*w**3;symbol=U**5+V**5+W**5
    for j in range(5):
        left=iterate(seed_symbol,seed**3,j,sv,sd)
        right=iterate(symbol,phase**3,4-j,av,ad)
        assert left==0 or right==0
    defect=[]
    for n in (2,3,4):
        r=n+1
        ind=lambda m:n*(m//r)+1
        value=ind(r)-ind(1)-ind(r-1);assert value==n-1
        defect.append(dict(n=n,r=r,defect=value))
    guards=0
    for symbol,poly in [(X,S.Float(1)),(X,x**65)]:
        try:action(symbol,poly,variables,dual)
        except AssertionError:guards+=1
        else:raise AssertionError('guard failed')
    try:iterate(X,x,13,variables,dual)
    except AssertionError:guards+=1
    else:raise AssertionError('guard failed')
    deps=['verify_formal_apolar_gaps.py','verify_gvc_paired_recomposition.py']
    report=dict(completed=True,checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        supporting_checker_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        sympy_version=S.__version__,independent_block_controls=convolution,
        safe_subadditivity_controls=safe,sharp_boundary_value=str(sharp),
        retained_seed_height_controls=seed_controls,auxiliary_direct_height_controls=auxiliary,
        periodic_gap_integer_controls=periods,actual_M1_tensor_window_zero=True,
        dimension_free_defect_controls=defect,guard_refusals=guards,
        exact_bounded_calculator_implemented=True,
        universal_height_calculus_proved_by_samples=False,
        all_M_gap_theorem_inferred_from_fitting=False,
        general_formal_gap_decided=False,
        new_frontier_endpoint_claimed=False,
        genuine_symmetric_gap_theorem_refuted=False,
        dimension_dependent_quasisubadditivity_disproved=False,
        all_ring_height_calculus_claimed=False,
        full_exceptional_MZ_decided=False,
        independent_replay=False,external_review=False,formal_verification=False,novelty_claimed=False)
    target=HERE/'gvc_false_gaps_verification.json';assert not target.exists()
    target.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(completed=True,convolution=len(convolution),seed=len(seed_controls),
                          auxiliary=len(auxiliary),periods=len(periods),guards=guards)))


if __name__=='__main__':main()
