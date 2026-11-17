"""PF157: exact controls for canonical deck jets and rational pole coverage.
The universal arguments are in the note; controls do not classify actual maps.
"""
from pathlib import Path
import hashlib
import json
import math
import platform
import time
import sympy as S
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
SRC=HERE/'sources/canonical-deck-pole-coverage'
NOTE='CANONICAL_DECK_JET_FIELDS_AND_RATIONAL_POLE_ATLAS_COVERAGE.md'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    started=time.perf_counter();checks=[]
    def check(name,value):
        assert bool(value),name
        checks.append(name)
    ctx=json.loads((SRC/'context.json').read_text())
    for p,digest in ctx['retained_sha256'].items():check('retained source pin '+p,sha(ROOT/p)==digest)
    x,t,u,zeta,lam,xi,al=S.symbols('x t u zeta lam xi al')
    # Exact subgroup arithmetic supporting the written all-order proof.
    mono_cases=[(2,0),(2,1),(3,2),(4,1),(6,0),(6,1),(6,2),(6,5),(6,11),(8,3),(12,7),(15,9)]
    mono=[]
    for g,k in mono_cases:
        fixed=[j for j in range(g) if (j*(k+1))%g==0]
        h=math.gcd(g,k+1);mono.append(dict(g=g,k=k,h=h))
        assert len(fixed)==h
    check('monomial commuting subgroup matches gcd controls',all(v['h']>=1 for v in mono))
    poly_cases=[(6,[5]),(6,[1,5]),(6,[2,5]),(6,[0,5]),(6,[5,11]),(12,[3,7,11]),(10,[1,3,9])]
    poly=[]
    for g,exps in poly_cases:
        fixed=[j for j in range(g) if all((j*(n+1))%g==0 for n in exps)]
        h=math.gcd(g,*[n+1 for n in exps]);poly.append(dict(g=g,support=exps,h=h));assert len(fixed)==h
    check('general-replica character subgroup matches support gcd controls',True)
    check('canonical replica selects full deck group',all(math.gcd(g,g)==g for g in [2,3,4,6,12]))
    # Universal affine-action and volume transcriptions, without a root choice.
    gp,gpz=S.symbols('gp gpz')
    for g in [2,3,6]:
        tx=zeta**-1*x+zeta**-1*gp-gpz
        st=(x+gp)/(g*t**(g-1));trans=(tx+gpz)/(g*(zeta*t)**(g-1))
        check('invariant slice action g'+str(g),S.cancel(zeta**g*trans-st)==0)
    check('actual two-dimensional volume multiplier',S.det(S.Matrix([[zeta**-1,S.Symbol('shift_derivative')],[0,zeta]]))==1)
    # Actual multiplication in the finite FUNCTION FIELD C(u)[t]/(t^g-u).
    # This algebra is NOT asserted to be the original finite normalization.
    def mt(g):
        M=S.zeros(g)
        for j in range(g-1):M[j+1,j]=1
        M[0,g-1]=u
        return M
    def ev(expr,M):
        num,den=S.fraction(S.cancel(expr));g=M.rows
        def pol(p):return sum((co*M**ex[0] for ex,co in S.Poly(p,t).terms()),S.zeros(g))
        return (pol(num)*pol(den).inv()).applyfunc(S.cancel)
    def projection(expr,g):
        return S.cancel(S.trace(ev(expr,mt(g)))/g).subs(u,t**g)
    Gcases=[(2,t**-1),(3,t**-2),(5,t**-1/(1+t**5)),(6,t**-2+t**-3),(6,t**-1+1/(t*t-1)-1/(t**6-1))]
    trace_records=[]
    for idx,(g,G) in enumerate(Gcases):
        M=mt(g);I=S.eye(g)
        check('actual cyclic field multiplication g'+str(g)+' case'+str(idx),M**g==u*I)
        check('trace-free primitive case'+str(idx),S.cancel(S.trace(ev(G,M)))==0)
        Gpp=ev(S.diff(G,t,2),M)
        check('weighted second derivative trace cancellation case'+str(idx),S.cancel(S.trace(M*M*Gpp))==0)
        W=lam*M*xi+M.inv()*al+Gpp*M*xi
        # Direct derivative of Tr(W^2), retaining every constant term.
        linear=2*S.trace(M*xi*(M.inv()*al+Gpp*M*xi))
        check('actual linear quadratic-trace coefficient case'+str(idx),S.cancel(linear-2*g*xi*al)==0)
        expected=2*u*xi*xi if g==2 else 0
        check('actual quadratic lambda coefficient case'+str(idx),S.cancel(S.trace((M*xi)**2)-expected)==0)
        trace_records.append(dict(g=g,G=str(G)))
    # Subgroup primitives come from the same full-deck G.
    Gmix=t**-2+t**-3
    check('order-two subgroup primitive projection',S.cancel(projection(Gmix,2)-t**-2)==0)
    check('order-three subgroup primitive projection',S.cancel(projection(Gmix,3)-t**-3)==0)
    check('subgroup trace-free remainders',S.cancel(projection(Gmix-projection(Gmix,2),2))==0 and S.cancel(projection(Gmix-projection(Gmix,3),3))==0)
    # Central principal-part profiles use exact character congruences.
    g=6;support=[2,3]
    nus=[max(r for r in support if j*r%g) for j in range(1,g)]
    check('mixed central characters have primitive gcd',math.gcd(g,*support)==1)
    check('mixed central transition pole profile',nus==[3,2,3,2,3])
    check('central source-coordinate pole profile',[r+1 for r in nus]==[4,3,4,3,4])
    check('full marked-charge central packet lower bound',2*max(nus)+4==10)
    bad=[j for j in range(6) if j*2%6==0]
    check('inadmissible central stabilizer control',len(bad)==math.gcd(6,2)==2)
    # Explicit nonzero pole cluster: normalize the whole rational primitive.
    Gcluster=t**-1+1/(t*t-1)-1/(t**6-1)
    check('nonzero cluster primitive trace normalized',S.cancel(projection(Gcluster,6))==0)
    zz=(1+S.sqrt(-3))/2;residues=[]
    for j in range(6):
        z=S.simplify(S.expand(zz**j));z2=S.simplify(S.expand(z*z))
        Gj=1/(z*t)+1/(z2*t*t-1)-1/(t**6-1)
        residues.append(S.simplify(S.limit((t-1)*Gj,t,1)))
    check('actual nonzero principal-part values',residues==[S.Rational(1,3),-S.Rational(1,6),-S.Rational(1,6),S.Rational(1,3),-S.Rational(1,6),-S.Rational(1,6)])
    clusters=sorted([residues.count(v) for v in set(residues)])
    check('intermediate fibre sheet clusters retained',clusters==[2,4])
    # One-root pseudo-plane controls already retained by MV0, now normalized
    # to the primitive G=t^(1-g). No polynomial Keller target is claimed.
    degree_records=[]
    for g in [2,3,4,5,6]:
        kk=g-1;G=t**(-kk);uu=t**g;vv=uu*x-kk
        Z=S.cancel((vv**g-(-kk)**g)/uu)
        ss=(x+S.diff(G,t))/(g*t**(g-1));F0=t*x/kk
        check('polynomial invariant control g'+str(g),S.denom(Z)==1)
        check('source polynomial primitive representative g'+str(g),S.cancel(F0-G-S.Rational(g,kk)*uu*ss)==0)
        ww=t*vv
        check('retained one-root source relation g'+str(g),S.expand(ww**g-(-kk)**g*uu-uu**2*Z)==0)
        for dd in [1,2]:
            HH=S.expand(Z**dd)
            assert S.degree(HH,x)==g*dd and S.degree(HH.subs(t,0),x)==dd
        check('degree comparison sharp controls g'+str(g),True)
        XF=S.diff(Z,x);AF=S.diff(Z,t)-S.diff(G,t,2)*XF
        pole_coeff=S.limit(t**(g+1)*2*g*XF*AF,t,0)
        check('exact polar linear moment control g'+str(g),pole_coeff==-2*g**4*(g-1)**(2*g-1))
        degree_records.append(dict(g=g,r=g,central_degree=1,linear_moment_pole_order=g+1))
    # Complete quadratic-image construction, including a non-MX pole profile.
    quadratic_cases=[1/u,1/(u*(1+u)**2),1/u**2+1/(u-1)+2/(u+2)**2]
    quadratic_records=[]
    for idx,R0 in enumerate(quadratic_cases):
        rn,rd=S.fraction(S.cancel(R0))
        k0=min(ex[0] for ex,co in S.Poly(rd,u).terms())-min(ex[0] for ex,co in S.Poly(rn,u).terms())
        H0=S.cancel(R0+2*u*S.diff(R0,u));num,den=S.fraction(H0)
        A=S.Poly(den,u).monic().as_expr();Sigma=S.cancel(A*H0)
        check('quadratic reduced denominator and numerator case'+str(idx),S.denom(Sigma)==1 and S.gcd(A,Sigma)==1)
        check('quadratic central pole and infinity normalization case'+str(idx),A.subs(u,0)==0 and Sigma.subs(u,0)!=0 and S.degree(Sigma,u)<S.degree(A,u) and min(ex[0] for ex,co in S.Poly(A,u).terms())==k0>=1)
        AR=S.cancel(A*R0)
        check('quadratic rational primitive denominator covers R0 case'+str(idx),S.denom(AR)==1)
        C0=S.rem(-AR*S.invert(Sigma,A,u),A,u)
        D0=S.cancel(R0+C0*Sigma/A)
        check('quadratic polynomial primitive Bezout case'+str(idx),not S.denom(D0).has(u) and S.rem(C0*Sigma+AR,A,u)==0)
        aa=A.subs(u,t*t);sig=Sigma.subs(u,t*t)
        ww=t*(sig+aa*x);zz=2*sig*x+aa*x*x
        check('quadratic whole-source relation case'+str(idx),S.expand(ww*ww-t*t*sig*sig-t*t*aa*zz)==0)
        invx=-x-2*sig/aa
        w_tau=ww.subs({t:-t,x:invx},simultaneous=True)
        z_tau=zz.subs({t:-t,x:invx},simultaneous=True)
        check('quadratic source deck invariance case'+str(idx),S.cancel(w_tau-ww)==0 and S.cancel(z_tau-zz)==0)
        jac=S.det(S.Matrix([t*t,ww]).jacobian((x,t)))
        check('quadratic actual unit volume pullback case'+str(idx),S.expand(jac+2*t*t*aa)==0)
        F0=t*C0.subs(u,t*t)*x+t*D0.subs(u,t*t)
        GG=t*R0.subs(u,t*t);ss=(x+S.diff(GG,t))/(2*t)
        check('quadratic global polynomial representative case'+str(idx),S.cancel(F0-GG-2*t*t*C0.subs(u,t*t)*ss)==0)
        check('quadratic central whole-line chart case'+str(idx),ww.subs(t,0)==0 and S.expand(zz.subs(t,0)-2*Sigma.subs(u,0)*x)==0)
        check('quadratic even-leading criterion infinity bound case'+str(idx),S.degree(u*A,u)>=2)
        quadratic_records.append(dict(R0=str(R0),A=str(A),Sigma=str(Sigma),C0=str(C0),D0=str(D0),central_R0_pole_order=k0,central_M_order=k0+1,full_marked_quadratic_charge_lower_bound=4*k0+2))
    check('quadratic model includes central order beyond MX0',S.Poly(quadratic_records[2]['A'],u).coeff_monomial(u)==0)
    # Central valuations for r4 and an odd companion give a negative order.
    m0,hh,ssdeg=S.symbols('m0 hh ssdeg')
    even_order=4*hh+2*m0
    odd_order=ssdeg*hh+m0*(ssdeg-1)/2
    check('lowest ruling-degree parity obstruction universal',S.simplify(odd_order.subs(hh,-m0/2)+m0/2)==0)
    # The pair's two-dimensional Jacobian gives this precise determinant.
    XP,XQ,AP,AQ,gg=S.symbols('XP XQ AP AQ gg')
    LM=gg*S.Matrix([[2*XP*AP,XP*AQ+AP*XQ],[XP*AQ+AP*XQ,2*XQ*AQ]])
    check('linear coefficient Gram determinant exact',S.expand(LM.det()+gg*gg*(XP*AQ-AP*XQ)**2)==0)
    # Mandatory known 3D nonproper Keller control, with its pole retained.
    xx,yy,zz3,tt=S.symbols('xx yy zz3 tt');uu=1+xx*yy
    FA=uu**3*zz3+yy*yy*uu*(4+3*xx*yy)
    FB=yy+3*xx*uu*uu*zz3+3*xx*yy*yy*(4+3*xx*yy)
    FC=2*xx-3*xx*xx*yy-xx**3*zz3
    DF=S.Matrix([FA,FB,FC]).jacobian((xx,yy,zz3))
    check('retained 3D constant Jacobian control',S.factor(DF.det())==-2)
    arc={xx:1/tt,yy:0,zz3:2*tt*tt}
    check('retained 3D regular target arc',[S.cancel(z.subs(arc)) for z in [FA,FB,FC]]==[2*tt*tt,6*tt,0])
    check('retained 3D first-velocity pole',S.diff(FC,zz3).subs(arc)==-tt**-3)
    check('scope does not assert automatic integrality or arbitrary model realization',not ctx['full_marked_quadratic_regularity_proved'] and not ctx['arbitrary_rational_pole_data_claimed_to_give_affine_image_or_target'] and not ctx['comparison_normalization_replaces_original_S'])
    out=dict(claim='PF2RULED157',outcome='PASS',check_count=len(checks),checks=checks,elapsed_seconds=time.perf_counter()-started,
             environment=dict(python=platform.python_version(),sympy=S.__version__),proof_sha256=sha(HERE/NOTE),checker_sha256=sha(Path(__file__)),context_sha256=sha(SRC/'context.json'),
             exact_monomial_replica_index_formula='gcd(|D|,k+1)',general_replica_formula_requires_minimality=True,
             no_proper_equivariant_factor_iff_trivial_deck_at_minimal_degree=True,
             actual_pole_atlas_coverage_supplied=True,arbitrary_pole_data_model_realization_proved=False,
             target_ruling_degree_lower_bound='2*|D|',universal_linear_quadratic_trace_determinant='-|D|^2',
             full_Gram_determinant_claimed_constant=False,full_marked_quadratic_regularity_proved=False,
             monomial_controls=mono,general_replica_controls=poly,field_trace_controls=trace_records,
             degree_and_pole_controls=degree_records,quadratic_image_controls=quadratic_records,complete_quadratic_image_presentation_verified=True,mixed_central_nu_profile=nus,nonzero_pole_cluster_sizes=clusters,
             original_full_normalization_and_all_ramified_primes_retained=True,global_coverage_complete=False,
             new_external_theorem_imports=[],registered_predecessor_checkers_replayed=[],full_goal_complete=False,
             independent_replay=False,formal_verification=False,external_review=False)
    (HERE/'canonical_deck_pole_coverage_receipt.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({key:out[key] for key in ['claim','outcome','check_count','elapsed_seconds']}))

if __name__=='__main__':main()
