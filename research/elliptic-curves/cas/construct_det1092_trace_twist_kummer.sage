#!/usr/bin/env sage-python
"""Exact trace/anti-trace Kummer identity; one old conic and nine old fibres.

No rational-point search, exceptional-point input, new address, number field,
Selmer dimension or class group. Freeze inputs first; 25-second total cap.
"""
import hashlib,json,signal
from pathlib import Path
from sage.all import QQ,PolynomialRing,EllipticCurve

ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
DIR=ART/'det1092_trace_twist_kummer_v1'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def retain(p,d):
    payload=json.dumps(d,indent=2,sort_keys=True)+'\n'
    if p.exists():assert p.read_text()==payload
    else:p.write_text(payload)
def rec(f):return {'numerator':list(map(str,f.numerator().list())),
                  'denominator':list(map(str,f.denominator().list()))}
def construct():
    paths=[ART/'curve302_recovered_mw17_parent_v1.json',
           ART/'det1092_orbit8044_rank18_base_change_v2.json',
           ART/'det1092_rr_generic_point_controls_v2/protocol.json']
    parent,cover,roster=map(read,paths)
    protocol={'classification':'new explicit obstruction and verified application',
      'rule':'Use the unchanged generic conic8044 and its inherited trace. Construct the rational anti-trace point on its quadratic twist, and compare its norm-square cubic Kummer class with the trace class by an explicit square identity. Evaluate only the nine old original-parameter addresses, with no new point or parameter selection.',
      'limits':{'wall_seconds':25,'generic_symbolic_identities':1,'old_conics':1,
                'old_addresses':9,'new_addresses':0,'point_searches':0,
                'exceptional_point_inputs':0,'V3_inputs':0,'pilot_changes':0,
                'class_groups':0,'Selmer_dimensions':0,'number_fields':0},
      'inputs':{str(p.relative_to(ROOT)):sha(p) for p in paths+[Path(__file__)]}}
    retain(DIR/'protocol.json',protocol)
    for name,digest in cover['inputs'].items():assert sha(ROOT/name)==digest

    # Universal chord calculation, with denominators cleared before expansion.
    B=PolynomialRing(QQ,names=('a','m','u','v','d','n'))
    a,m,u,v,d,n=B.gens();V=PolynomialRing(B,'X');X=V.gen()
    c=m*m-a-2*u;L=m*X+n;g=(X-u)**2-v*v*d
    f=L*L+(X-c)*g;w=m*u+n
    rnum=w*w-v*v*d*(a+2*u)
    Nnum=w*(X-u)+d*m*v*v
    assert f[3]==1 and f[2]==a
    assert d*v*v*f-Nnum*Nnum==g*(d*v*v*X-rnum)
    assert ((rnum-v*v*d*X)*L*L-(c-X)*Nnum*Nnum)%f==0
    universal={'classification':'new explicit algebraic deduction',
      'variables':['a','m','u','v','d','n'],
      'definitions':{'c':'m^2-a-2u','L':'mX+n','g':'(X-u)^2-v^2*d',
        'f':'L^2+(X-c)*g','w':'mu+n','rnum':'w^2-v^2*d*(a+2u)',
        'Nnum':'w*(X-u)+d*m*v^2'},
      'polynomial_identity':'d*v^2*f-Nnum^2=g*(d*v^2*X-rnum)',
      'kummer_square_identity':'d*(r-theta)=(c-theta)*(Nnum(theta)/(v*L(theta)))^2',
      'r':'rnum/(d*v^2)',
      'patch':'char(K)!=2, d*v!=0, smooth f, and (c-theta),L(theta),Nnum(theta) units',
      'interpretation':'For Q=(u+v*sqrt(d),w+m*v*sqrt(d)), trace Z=Q+sigma(Q) and R=Q-sigma(Q), delta_Ed(R)=delta_E(Z) in the naturally identified cubic norm kernels.',
      'coefficient_residuals_zero':True}
    retain(DIR/'universal.json',universal)

    T=PolynomialRing(QQ,'t');K=T.fraction_field()
    def dec(item):return K(T(item['numerator']))/T(item['denominator'])
    ai=list(map(dec,parent['a_invariants']));E=EllipticCurve(K,ai)
    rr_c,rr_b,rr_a=map(dec,cover['lift']['residual_coefficients'])
    h=K(T(cover['splitting']['discriminant_square_factor']))
    d=K(T(cover['curve_over_Q']['q_coefficients']))
    f0,f1,f2=[K(T(c)) for c in cover['lift']['line_coefficients']]
    x0=-rr_b/(2*rr_a);x1=h/(2*rr_a)
    y0=-(f0+f1*x0)/f2;y1=-f1*x1/f2
    u,v=4*x0,4*x1
    w,z=8*y0+4*ai[0]*x0+4*ai[2],8*y1+4*ai[0]*x1
    m=z/v;n=w-m*u;lam=w/v
    W=PolynomialRing(K,'X');X=W.gen()
    cubic=[E.b2(),8*E.b4(),16*E.b6()]
    a,b,c0=cubic;f=X**3+a*X*X+b*X+c0;L=m*X+n;g=(X-u)**2-v*v*d
    quotient,remainder=(f-L*L).quo_rem(g)
    assert not remainder and quotient.degree()==1 and quotient[1]==1
    cx=-quotient[0];cy=-L(cx)
    basis=[E(list(map(dec,P))) for P in parent['basis_weierstrass_coordinates']]
    trace=sum((int(n)*P for n,P in zip(cover['lift']['trace_word'],basis)),E(0))
    assert cx==4*trace[0] and cy==8*trace[1]+4*ai[0]*trace[0]+4*ai[2]
    r=lam*lam/d-a-2*u;eta=lam*(u-r)/d-z
    N=lam*(X-u)+d*z
    assert d*eta*eta==f(r)
    assert d*f-N*N==g*d*(X-r)
    assert (d*(r-X)*L*L-(cx-X)*N*N)%f==0
    # Unit claims are verified separately in every rational replay below.
    payload={'classification':'equation-only generic construction',
      'cubic_coefficients':list(map(rec,cubic)),'d':rec(d),
      'quadratic_point':{'X':[rec(u),rec(v)],'Y':[rec(w),rec(z)]},
      'trace':[rec(cx),rec(cy)],'trace_word':cover['lift']['trace_word'],
      'anti_twist_point':[rec(r),rec(eta)],'line':[rec(n),rec(m)],
      'square_numerator':[rec(N[0]),rec(N[1])],
      'square_denominator':[rec(L[0]),rec(L[1])],
      'twist_model':'d*eta^2=f(r); monic model coordinates (d*r,d^2*eta)',
      'kummer_representatives':'alpha_twist=d*(r-theta); alpha_trace=cx-theta',
      'inputs':{str(p.relative_to(ROOT)):sha(p) for p in paths+[DIR/'protocol.json']}}
    retain(DIR/'construction.json',payload)
    print('CHECKPOINT_UNIVERSAL_AND_GENERIC_SQUARE_IDENTITIES',flush=True)
    results=[]
    for i,case in enumerate(roster['cases']):
        tau=QQ(case['parameter']);dv=d(tau)
        aa,bb,cc=[v(tau) for v in cubic]
        Z=PolynomialRing(QQ,'X');XX=Z.gen();ff=XX**3+aa*XX**2+bb*XX+cc
        assert ff.discriminant() and dv
        rr,ee=r(tau),eta(tau);cxv,cyv=cx(tau),cy(tau)
        LL=Z([n(tau),m(tau)]);NN=Z([N[0](tau),N[1](tau)])
        assert LL.gcd(ff)==1 and NN.gcd(ff)==1 and (cxv-XX).gcd(ff)==1
        gamma=(NN*LL.inverse_mod(ff))%ff
        assert (dv*(rr-XX)-(cxv-XX)*gamma*gamma)%ff==0
        assert dv*ee*ee==ff(rr) and cyv*cyv==ff(cxv)
        assert dv>0
        dn,dd=dv.numerator(),dv.denominator();sn,sd=dn.isqrt(),dd.isqrt()
        assert not (sn*sn==dn and sd*sd==dd)
        row={'classification':'verified application candidate; independent replay required',
          'index':i,'label':case['label'],'parameter':str(tau),'d':str(dv),
          'cubic_coefficients':list(map(str,[aa,bb,cc])),
          'trace':list(map(str,[cxv,cyv])),
          'anti_twist_point':list(map(str,[rr,ee])),
          'square_witness':list(map(str,gamma.list())),
          'nonsquare_bounds':list(map(str,[sn,sd])),
          'relative_elliptic_2_kummer_class':'ZERO_INHERITED_TRACE',
          'elliptic_seed_existence':'UNKNOWN_NOT_EXCLUDED'}
        retain(DIR/f'case-{i:02d}.json',row);results.append(row)
    retain(DIR/'panel.json',{'status':'CANDIDATE_TRACE_TWIST_KUMMER_COLLAPSE',
      'classification':'new explicit obstruction and verified application',
      'cases':results,'limits':protocol['limits'],
      'inputs':{str(p.relative_to(ROOT)):sha(p) for p in
                [DIR/'protocol.json',DIR/'universal.json',DIR/'construction.json',Path(__file__)]}})
    print('CANDIDATE_TRACE_TWIST_KUMMER_COLLAPSE','all nine classes inherited',flush=True)

if __name__=='__main__':
    signal.alarm(25);DIR.mkdir(parents=True,exist_ok=True);construct()
