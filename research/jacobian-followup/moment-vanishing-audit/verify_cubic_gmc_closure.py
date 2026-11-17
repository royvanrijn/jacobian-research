#!/usr/bin/env python3
"""Exact reductions/sharpness regressions for the written cubic GMC closure."""
from pathlib import Path
import hashlib,json,math
import sympy as S


def wick(poly,variables,pair_count=0,real_indices=()):
    result=S.Integer(0)
    for powers,coefficient in S.Poly(S.expand(poly),*variables).terms():
        factor=S.Integer(1)
        for j in range(pair_count):
            x,y=powers[2*j:2*j+2]
            if x!=y:factor=0;break
            factor*=math.factorial(x)
        if not factor:continue
        for j in real_indices:
            n=powers[j]
            if n%2:factor=0;break
            factor*=S.factorial2(n-1) if n else 1
        result+=coefficient*factor
    return S.expand(result)


def main():
    directory=Path(__file__).resolve().parent
    t,u=S.symbols('t u')
    aa,ab,bb,p1,p2=S.symbols('aa ab bb p1 p2')
    A=S.Matrix([[aa,ab],[ab,bb]]);p=S.Matrix([p1,p2]);I=S.eye(2)
    determinant=S.factor((I-2*t*A).det())
    C=(I-2*t*A).inv()
    quartic=S.factor(t*t*(p.T*C*p)[0]/2)
    pp=(p.T*p)[0];pAp=(p.T*A*p)[0]
    expected=t*t*(pp-2*t*(S.trace(A)*pp-pAp))/(2*determinant)
    assert S.factor(quartic-expected)==0
    coefficients=S.Poly(S.cancel(2*determinant*quartic/t**2),t).all_coeffs()
    assert all(S.expand(x-y)==0 for x,y in zip(coefficients,[-2*S.trace(A)*pp+2*pAp,pp]))
    # A null p with pAp=0 is exactly the needed transverse frame direction.
    n=S.Matrix([1,S.I])
    assert (n.T*n)[0]==0
    assert S.factor((n.T*A*n)[0]-(aa-bb+2*S.I*ab))==0
    nilpotent=n*n.T
    assert nilpotent*nilpotent==S.zeros(2)
    assert S.expand((p.T*nilpotent*p)[0]-(p1+S.I*p2)**2)==0
    # Exact infinity alternatives for the only quartic degree drop.
    assert S.factor(S.limit(quartic/t,t,S.oo)+(p.T*A.inv()*p)[0]/4)==0
    a,d,pn=S.symbols('a d pn',nonzero=True)
    # Resolvent numerator linear in t: if its t coefficient vanishes,
    # the nonzero remaining quartic tends to pn/(8 det A).
    assert S.limit(t*t*pn/(2*(1-2*a*t+4*d*t*t)),t,S.oo)==pn/(8*d)

    # Isotropic repeated factor: generic-t quadratic elimination, not a
    # conditional t-adic Gaussian average in the null W coordinate.
    Z,T,W=S.symbols('Z T W')
    lam,c,At,qt,b0,b1=S.symbols('lam c At qt b0 b1',nonzero=True)
    H=Z*Z+b1*Z+c*T+b0
    G2=S.Symbol('g20')*Z*Z+S.Symbol('g11')*Z*T+qt*T*T+S.Symbol('g10')*Z+S.Symbol('g01')*T+S.Symbol('g00')
    cubic_a=S.Symbol('cubic_a')
    original=t*(lam*W*W+W*H+cubic_a*Z**3+G2)-Z*W-T*T/2
    shift=(Z-t*H)/(2*lam*t)
    effective=t*(cubic_a*Z**3+G2-H*H/(4*lam))+Z*H/(2*lam)-T*T/2-Z*Z/(4*lam*t)
    assert S.expand(original.subs(W,W+shift)-(lam*t*W*W+effective))==0
    qT=t*(qt-c*c/(4*lam))-S.Rational(1,2)
    qL=S.expand(effective).coeff(T,1)
    qK=S.expand(effective).coeff(T,0)
    univariate=S.cancel(qK-qL*qL/(4*qT))
    c4=S.factor(S.expand(qK).coeff(Z,4)-S.expand(qL*qL).coeff(Z,4)/(4*qT))
    assert S.factor(c4-(-t/(4*lam)-c*c*t*t/(16*lam*lam*qT)))==0
    assert S.factor(c4.subs(qt,0)+t/(2*(c*c*t+2*lam)))==0
    assert S.limit(c4.subs(qt,0),t,S.oo)==-1/(2*c*c)
    assert S.factor(S.limit(c4/t,t,S.oo)+qt/(4*lam*(qt-c*c/(4*lam))))==0
    # For leading Z^2*T, integrating T gives a nonzero quartic in either
    # the growing or the flattened quadratic branch.
    other=t*(Z*Z*T+G2-(b1*Z+c*T+b0)**2/(4*lam))+Z*(b1*Z+c*T+b0)/(2*lam)-T*T/2-Z*Z/(4*lam*t)
    cT=S.expand(other).coeff(T,2);lT=S.expand(other).coeff(T,1)
    other_c4=S.factor(-S.expand(lT*lT).coeff(Z,4)/(4*cT))
    assert S.factor(other_c4+t*t/(4*cT))==0
    assert S.factor(other_c4.subs(qt,c*c/(4*lam))-t*t/2)==0

    # All distinct degree-drop balances used below have shrinking widths
    # or a nonzero escaping action; count bounded+escaping dimensions.
    balances=[]
    for big_degree,small_degree,big_power,small_power in [(4,1,1,0),(4,2,1,0),(4,3,1,0),(4,2,2,1),(3,2,2,1),(3,2,3,1),(3,1,2,1)]:
        scale=S.Rational(big_power-small_power,big_degree-small_degree)
        action_power=big_power+small_degree*scale
        assert scale>0 and action_power>big_power
        assert (big_degree-small_degree)+(small_degree-1)==big_degree-1
        balances.append({'degree':big_degree,'leading_degree':small_degree,'powers':[big_power,small_power],
                         'escape_scale':str(scale),'escape_action_power':str(action_power),
                         'escape_count':big_degree-small_degree,'bounded_rank':small_degree-1})
    # Sharp mixed cubic shear and the known minimal degree/dimension examples.
    shear=W+Z*T-Z**3/2
    shear_rows=[]
    for n in range(1,9):
        assert wick(shear**n,(Z,W,T),1,(2,))==0
    for degree in range(1,4):
        value=wick(W**degree*shear**(3*degree),(Z,W,T),1,(2,))
        expected=S.factorial(3*degree)*(-S.Rational(1,2))**degree
        assert value==expected
        shear_rows.append({'multiplier_degree':degree,'power':3*degree,'mixed':str(value)})
    Z1,W1,Z2,W2=S.symbols('Z1 W1 Z2 W2')
    four=(1+Z2)*(W1*(1-Z1)+W2)
    three=(1+Z)*(W-(2+Z)*T*T/2)
    for n in range(1,6):
        assert wick(four**n,(Z1,W1,Z2,W2),2)==0
        assert wick(Z2*four**n,(Z1,W1,Z2,W2),2)==math.factorial(n)
        assert wick(three**n,(Z,W,T),1,(2,))==0
        assert wick(Z*three**n,(Z,W,T),1,(2,))==math.factorial(n)
    prior=json.loads((directory/'affine_line_conic_profiles_verification.json').read_text())
    assert prior['completed']
    result={'completed':True,'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'prior_squarefree_checker_sha256':prior['checker_sha256'],
            'quartic_identity':'t^2*(p.p-2t*(tr(A)*p.p-p.A.p))/(2 det(I-2tA))',
            'quartic_zero_frame':'For p!=0, identically zero iff p.p=p.A.p=0, giving an isotropic singular direction with F2=0',
            'degree_drop_balances':balances,'sharp_shear_regressions':shear_rows,
            'Long_regressions':{'orders':5,'three_variable_degree':4,'four_variable_degree':3,'mixed':'m!'},
            'resolved':'Every complex polynomial of total degree at most3 in three Gaussian coordinates satisfies GMC. Pure nullity forces a Gaussian W-linear frame; sufficient mixed cutoff m>3*deg Q is sharp.',
            'minimal_failure_consequence':'Degree4 is minimal for a three-variable GMC failure; four variables are minimal for a cubic GMC failure, using Long\u2019s explicit examples and the quadratic theorem.',
            'written_proof_certified_by_checker':False,
            'proof_boundary':'Universal resolvent/quartic and isotropic elimination identities, model degree/action/count arithmetic and stated bounded Wick regressions. The one-variable period exclusions, exhaustive normal-form coverage and all-order closure/corollaries are written proofs, not machine verified.',
            'remaining':'Independent specialist replay and external review; explicit global pure prefixes and global term-count optimality; coefficient schemes and lower-power GVC.'}
    (directory/'cubic_gmc_closure_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({key:result[key] for key in ['completed','resolved','minimal_failure_consequence','written_proof_certified_by_checker']}))


if __name__=='__main__':main()
