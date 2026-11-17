#!/usr/bin/env python3
"""Bounded exact transcription controls for the written quotient theorem."""
from pathlib import Path
import hashlib,json
import sympy as S

HERE=Path(__file__).resolve().parent

def action(symbol,poly,variables,dual):
    assert not S.sympify(symbol).has(S.Float) and not S.sympify(poly).has(S.Float)
    assert S.Poly(poly,*variables).total_degree()<=32
    result=0
    for powers,coefficient in S.Poly(S.expand(symbol),*dual).terms():
        term=poly
        for x,k in zip(variables,powers):term=S.diff(term,x,k)
        result+=coefficient*term
    return S.expand(result)

def iterate(symbol,poly,k,variables,dual):
    assert isinstance(k,int) and 0<=k<=8
    for _ in range(k):poly=action(symbol,poly,variables,dual)
    return poly

def main():
    t,y,T,Y,e=S.symbols('t y T Y e');v=(t,y);z=(T,Y)
    phase=y**2;mu=Y**3;tail=T*(1+Y+T*Y**2);lam=mu+tail
    shifted=[]
    for j,k in [(0,0),(0,1),(1,1),(2,1),(2,2),(3,2),(3,3)]:
        a=iterate(lam,phase**j,k,v,z);b=iterate(mu,phase**j,k,v,z)
        assert a==b
        shifted.append(dict(j=j,k=k,equal=True,nonzero=a!=0))
    # Nilpotent coefficient survives in the marked difference, while pure
    # output identities hold before any reduction of coefficients.
    nil_lam=mu+e*T
    for j,k in [(1,1),(2,1),(3,2)]:
        assert iterate(nil_lam,phase**j,k,v,z)==iterate(mu,phase**j,k,v,z)
    difference=iterate(nil_lam,t*phase,1,v,z)-iterate(mu,t*phase,1,v,z)
    reduced=S.rem(difference,e**2,e)
    assert reduced==e*y**2 and reduced!=0
    sharp=iterate(Y**3+T,t*y**6,3,v,z)
    assert sharp==2160
    assert iterate(Y**3+T,t*y**8,4,v,z)==0
    assert iterate(Y**3,y*y**2,1,v,z)==6
    a,b,c,s,A,B,C,SS=S.symbols('a b c s A B C SS')
    cv=(a,b,c,s);cz=(A,B,C,SS)
    cp=a+2*s*b+s**2*c;ann=A*C-B**2/4
    concise=[]
    for j in (0,1,2,3,4):
        assert action(ann,cp**j,cv,cz)==0
        # The same congruence works with an active operator, not only zero.
        assert iterate(ann+SS,cp**j,2,cv,cz)==iterate(SS,cp**j,2,cv,cz)
        concise.append(j)
    va,vb,vc,vs=S.symbols('va vb vc vs')
    derivative=va*S.diff(cp,a)+vb*S.diff(cp,b)+vc*S.diff(cp,c)+vs*S.diff(cp,s)
    equations=S.Poly(derivative,*cv).coeffs()
    assert S.linsolve(equations,(va,vb,vc,vs))==S.FiniteSet((0,0,0,0))
    # A finite low jet of the lacunary example; no transcendence claim is
    # deduced from this finite polynomial calculation.
    x,w,X,W=S.symbols('x w X W');lv=(x,w,y);lz=(X,W,Y)
    phi=W**2+W**6;ll=Y**3+X-phi
    lacunary=[]
    for j,k in [(1,1),(2,1),(3,2)]:
        assert iterate(ll,y**(2*j),k,lv,lz)==iterate(Y**3,y**(2*j),k,lv,lz)
        lacunary.append([j,k])
    assert iterate(ll,x*y**6,3,lv,lz)==2160
    assert action(Y**2,y,v,z)==0 and action(Y**2,y**2,v,z)==2
    assert action(T,t*y,v,z)==y and action(0,t*y,v,z)==0
    guards=0
    for sym,poly in [(Y,S.Float(1)),(Y,y**33)]:
        try:action(sym,poly,v,z)
        except AssertionError:guards+=1
        else:raise AssertionError('guard failed')
    try:iterate(Y,y,9,v,z)
    except AssertionError:guards+=1
    else:raise AssertionError('guard failed')
    deps=['verify_formal_apolar_gaps.py','verify_formal_binary_gvc.py',
          'verify_algebraic_formal_symbols.py','verify_formal_symbol_source_rigidity.py']
    report=dict(completed=True,checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        supporting_checker_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        shifted_output_controls=shifted,dual_number_marked_difference=str(reduced),
        sharp_boundary_value=str(sharp),beyond_sharp_boundary_zero=True,
        concise_annihilator_power_controls=concise,concise_linear_kernel_zero=True,
        lacunary_low_jet_controls=lacunary,guard_refusals=guards,
        exact_bounded_calculator_implemented=True,
        universal_quotient_theorem_proved_by_samples=False,
        arbitrary_annihilator_computation_implemented=False,
        general_algebraic_representative_recognition_implemented=False,
        least_cutoffs_preserved=False,full_raw_source_made_holonomic=False,
        general_formal_gap_decided=False,full_exceptional_MZ_decided=False,
        new_frontier_endpoint_claimed=False,independent_replay=False,
        external_review=False,formal_verification=False,novelty_claimed=False)
    (HERE/'phase_power_quotients_verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ['completed','checker_sha256','sharp_boundary_value','guard_refusals']}))

if __name__=='__main__':main()
