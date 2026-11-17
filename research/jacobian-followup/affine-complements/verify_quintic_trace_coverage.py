"""PF155: complete quintic trace coverage and exact chart certificates.
Universal field proof plus actual polynomial identities. No parameter sampling.
"""
from pathlib import Path
import hashlib
import importlib.util
import json
import platform
import time
import sympy as S
from sympy.polys.rings import ring
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
SRC=HERE/'sources/second-jet-coverage'
NOTE='QUINTIC_GRAPH_COVERAGE_AND_TWO_CONJUGATE_TRACE_CLASSES.md'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    start=time.perf_counter();checks=[]
    def check(name,value):
        assert bool(value),name
        checks.append(name)
    ctx=json.loads((SRC/'context.json').read_text())
    for p,v in ctx['retained_sha256'].items():check('retained source pin '+p,sha(ROOT/p)==v)
    for p,v in ctx['certificate_sha256'].items():check('certificate pin '+p,sha(SRC/p)==v)
    spec=importlib.util.spec_from_file_location('hc4_trace_witness_derivation',SRC/'derive_trace_witnesses.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    result=module.derive_witnesses()
    x,y,p,q=result['source_symbols']
    kap,A,D,h,i,e,rho,a,b,c,d,g,j,k,l,m,n,C=result['params']
    loc={str(z):z for z in result['params']}
    # Independently derive the potential jet from the actual parent, so the
    # witness extraction is tied to that potential rather than a copied jet.
    tr,ss=S.symbols('tr ss');uu=1+x*y
    AA=uu**3*p+3*x*uu**2*q-x**3*(tr+x*ss)
    BB=y*y*uu*(4+3*x*y)*p+(y+3*x*y*y*(4+3*x*y))*q+(2*x-3*x*x*y)*(tr+x*ss)
    parent_jet=S.Poly(S.expand(AA*AA+13*AA+2*BB),x)
    check('actual parent potential jet through x3 independently derived',
          all(S.expand(parent_jet.nth(jj).subs({tr:result['trace'],ss:result['normal']})-result['potential_coefficients'][jj])==0 for jj in range(4)))
    data=json.loads((SRC/'witnesses.json').read_text())
    ws=result['witnesses']
    for name,w in zip(('y7','y5p','y4q'),ws):
        check('actual full trace Hessian-jet witness '+name,S.expand(w-S.sympify(data['witnesses'][name],locals=loc))==0)
    # Cross-control from PF154's exact quadratic graph, independently of
    # the determinant-witness extraction and without a full normal prefix.
    sample={kap:-S.Rational(27,2),A:0,D:1,h:0,i:0,e:-S.Rational(3,2),rho:-6,
            a:-S.Rational(3,2),b:0,c:0,d:0,g:-S.Rational(39,4),j:0,
            k:S.Rational(1,4),l:-S.Rational(39,4),m:0,n:0,C:4}
    A0=uu**3*p+3*x*uu**2*q
    B0=y*y*uu*(4+3*x*y)*p+(y+3*x*y*y*(4+3*x*y))*q
    bb=18*y**3+1
    Fctrl=p*p+(13+8*y*y)*p+2*y*q+x*(bb*p-54*y**5)+x*x*(bb*bb-1)/4
    cc=S.Poly(S.expand(2*(2-3*x*y)-13*x*x-2*x*x*A0),x)
    dd=S.Poly(S.cancel((A0*A0+13*A0+2*B0-Fctrl)/x),x)
    tt=S.expand(-dd.nth(0)/4)
    sn=S.expand(-(cc.nth(1)*tt+dd.nth(1))/4)
    a2=S.Poly(S.expand(-(cc.nth(1)*sn+cc.nth(2)*tt+dd.nth(2))/4),y,p,q)
    check('three witnesses match independent PF154 exact quadratic control',
          all(S.expand(w.subs(sample)-192*a2.coeff_monomial(slot))==0 for w,slot in zip(ws,(y**7,y**5*p,y**4*q))))
    f,V=S.symbols('f V');loc.update(f=f,V=V)
    T=result['trace']+f*p**3
    Ky=3*(2*p+13)*p+2*(21*y*y*p+24*y*q+2*S.diff(T,y))
    Kp=12*y*p+6*q+39*y+14*y**3+4*S.diff(T,p)
    Kq=6*p+39+24*y*y+4*S.diff(T,q)
    HH=30*y*y*p*p+60*y*p*q+18*q*q+78*y*y*p+156*y*q+12*y**4*p+36*y**3*q-12*y*T
    plane=S.Poly(S.expand(-8*HH+8*Ky*Kq+4*Kp*Kp-64*y*Kp*Kq-(32*p-256*y*y)*Kq*Kq),y,p,q)
    W=2*a+3;RR=20*rho+123;P=8*rho*rho+99*rho+279;Q=160*rho*rho+1968*rho+6021
    DS=(16*rho+89-V)/2
    equations=[V*V+2*kap*RR,
               32*V*a+40*a*kap-8*V*e-12*V*rho-39*V+60*kap+Q,
               c*(4*V+5*kap),f*V,
               2*V*i-(16*rho+99)*A-(16*V+20*kap)*b]
    slots=[y**6,y**4*p,y**4*q,y**3*p*p,y**5];scales=[16,16,256,-192,-32]
    for idx,(slot,sc,eq) in enumerate(zip(slots,scales,equations)):
        check('complete field normal-form plane equation '+str(idx),S.expand(plane.coeff_monomial(slot).subs(D,DS)-sc*eq)==0)
    S4=S.Poly(S.expand(sum(co*y**ex[0]*p**ex[1]*q**ex[2] for ex,co in plane.terms() if sum(ex)==4)/64),y,p,q)
    check('full trace top quartic forces f on nonzero p-kernel',S4.coeff_monomial(p**4)==9*f*f and S4.coeff_monomial(p**3*q)==0)
    check('full trace q-square kernel coefficient',S4.coeff_monomial(y*y*q*q)==256*c*c)
    Z=6*D*a+9*D+64*a*a-32*a*e-16*a*rho+69*a+4*e*e+2*e*rho-S.Rational(45,2)*e-21*rho-S.Rational(69,4)
    KK=12*D+256*a-64*e-32*rho+138
    G=(P-W*(16*rho-2*D+89))/2+16*A*c+256*b*c-32*c*i
    check('actual full quartic quadratic coefficient Z',S.expand(S4.coeff_monomial(y*y*p*p).subs(f,0)-Z)==0)
    check('actual full quartic mixed coefficient cK',S.expand(S4.coeff_monomial(y*y*p*q).subs(f,0)-c*KK)==0)
    check('actual full quartic linear q coefficient G',S.expand(S4.coeff_monomial(y**3*q)-G)==0)
    cancelled={A,i,b,h,d,g,j,k,l,m,n,C}
    c0=[S.expand(w.subs(c,0)) for w in ws]
    for name,w in zip(('y7','y5p','y4q'),c0):check('all omitted trace parameters cancel on c0 '+name,not (w.free_symbols&cancelled))
    ke=-V*V/(2*RR)
    ee=(V*(32*a-12*rho-39)+20*ke*W+Q)/(8*V)
    generic=[S.factor(S.cancel(w.subs({D:DS,kap:ke,e:ee}))) for w in c0]
    # Exact V membership, with inputs reidentified from actual witnesses.
    numerators=[S.fraction(S.factor(generic[0]/(24*V)))[0],S.fraction(generic[1])[0],S.fraction(generic[2])[0],S.fraction(S.factor(S.cancel(Z.subs({D:DS,e:ee}))))[0]]
    cert=json.loads((SRC/'V_membership_certificate.json').read_text())
    va,vt,vr=S.symbols('a V rho');vloc={'a':va,'V':vt,'rho':vr}
    fi=[S.sympify(z,locals=vloc) for z in cert['input']]
    for idx,(got,expected) in enumerate(zip(numerators,fi)):
        ratio=S.cancel(got/expected)
        check('V certificate input from actual witness '+str(idx),ratio.is_Rational and ratio!=0)
    Ring,ra,rv,rr=ring('a,V,rho',S.QQ,order='grevlex')
    fp=[Ring.from_expr(z) for z in fi];hp=[Ring.from_expr(S.sympify(z,locals=vloc)) for z in cert['multipliers']]
    check('explicit polynomial membership sum equals V',sum((u*v for u,v in zip(fp,hp)),Ring.zero)==rv)
    # q-kernel generic normalization: no tacit loss at W=0 or P=0.
    qgen=[S.factor(S.cancel(z.subs(V,P/W))) for z in generic]
    E=[S.factor(S.cancel(qgen[0]*W*W*RR/(24*P))),S.factor(S.cancel(-qgen[1]*W*W*RR*P/24)),S.factor(-qgen[2]/48)]
    qdata=json.loads((SRC/'q_kernel_equations.json').read_text());qloc={'a':a,'rho':rho}
    for idx,(got,ex) in enumerate(zip(E,qdata['equations'])):check('full-trace q-kernel certificate input '+str(idx),S.expand(got-S.sympify(ex,locals=qloc))==0)
    qc=json.loads((SRC/'q_kernel_unit_certificate.json').read_text())
    Rq,qa,qr=ring('a,rho',S.QQ,order='grevlex')
    qp=[Rq.from_expr(S.sympify(z,locals=qloc)) for z in qc['input']]
    qh=[Rq.from_expr(S.sympify(z,locals=qloc)) for z in qc['multipliers']]
    check('q certificate inputs identical',all(S.expand(S.sympify(z,locals=qloc)-w)==0 for z,w in zip(qc['input'],E)))
    check('explicit q-kernel polynomial identity equals one',sum((u*v for u,v in zip(qp,qh)),Rq.zero)==Rq.one)
    # Missing q-kernel denominator chart W=0.
    boundary_w5=S.factor(generic[2].subs(a,-S.Rational(3,2)))
    check('q-kernel W0 exact witness independent V',not boundary_w5.has(V))
    check('q-kernel W0 witness remainder',S.expand(S.rem(boundary_w5,P,rho)-540*(25*rho+141))==0)
    check('q-kernel W0 resultant',S.resultant(P,25*rho+141,rho)==-15552)
    check('q-kernel V0 incompatible with nonzero trace',P.subs(rho,-S.Rational(123,20))==-S.Rational(2727,100))
    # The other kernel chart's V=0 boundary has an actual unit witness.
    b7=S.factor(c0[0].subs({rho:-S.Rational(123,20),D:-S.Rational(47,10),kap:S.Rational(153,100)/W}))
    check('non-q kernel V0 exact nonzero witness',S.cancel(b7-S.Rational(66096,125)/W)==0)
    # Nonzero c chart: retain complex parameters until a Bezout certificate.
    vv=8*RR/5;ks=-32*RR/25;es=-(1120*rho*rho+15888*rho+55503)/(64*RR)
    sub={D:(16*rho+89-vv)/2,kap:ks,e:es,i:(16*rho+99)*A/(2*vv)}
    poly_rank=4000*rho*rho+49520*rho+157969
    poly_trans=160*rho*rho+1968*rho+5841
    check('c nonzero rank identity',S.cancel((1024*Z-KK*KK).subs(sub)+S.Rational(24,25)*poly_rank)==0)
    check('c nonzero complete-trace witness cancels all extra parameters',S.cancel(ws[0].subs(sub)-S.Rational(192,25)*RR*poly_trans)==0)
    check('c nonzero exact resultant',S.resultant(poly_rank,poly_trans,rho)==2544276275200)
    bez=(S.Rational(91,1129384)-5*rho/1552903)*poly_rank+(125*rho/1552903-S.Rational(24945,12423224))*poly_trans
    check('c nonzero explicit Bezout identity',S.expand(bez)==1)
    # kappa=0 branch is retained, including both conjugates and arbitrary f.
    check('residual A coefficient unit modulo Q',S.resultant(Q,16*rho+99,rho)==-7776)
    for sign in (1,-1):
        rs=-S.Rational(123,20)+sign*3*S.sqrt(34)/40
        ds=-S.Rational(47,10)+sign*3*S.sqrt(34)/5
        check('retained conjugate residual rho '+str(sign),S.simplify(Q.subs(rho,rs))==0)
        check('retained conjugate residual D '+str(sign),S.simplify((16*rs+89)/2-ds)==0)
    check('residual quadratic discriminant requires sqrt34',S.discriminant(Q,rho)==576*34 and S.integer_nthroot(34,2)[1] is False)
    check('scheme and global coverage limitations retained',not ctx['nonreduced_kappa_nilpotence_certificate_proved'] and not ctx['A2_alone_claimed_empty'] and not ctx['global_minimal_counterexample_coverage_proved'])
    out=dict(claim='PF2RULED155',outcome='PASS',check_count=len(checks),checks=checks,elapsed_seconds=time.perf_counter()-start,
             environment=dict(python=platform.python_version(),sympy=S.__version__),proof_sha256=sha(HERE/NOTE),checker_sha256=sha(Path(__file__)),context_sha256=sha(SRC/'context.json'),
             field_scope=ctx['field_scope'],complete_18_parameter_trace_coverage=True,
             all_nonzero_quintic_trace_graphs_excluded=True,rational_degree_at_most_five_graphs_excluded=True,
             conjugate_residual_trace_classes=2,nonreduced_kappa_nilpotence_proved=False,A2_alone_empty_proved=False,
             explicit_V_membership_identity_verified=True,explicit_q_kernel_unit_identity_verified=True,
             complex_c_nonzero_branch_excluded_by_Bezout=True,global_coverage_proved=False,
             full_normal_increasing_power_containment_proved=False,new_external_theorem_imports=[],registered_predecessor_checkers_replayed=[],
             full_goal_complete=False,independent_replay=False,formal_verification=False,external_review=False)
    (HERE/'quintic_trace_coverage_receipt.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k] for k in ['claim','outcome','check_count','elapsed_seconds']}))

if __name__=='__main__':main()
