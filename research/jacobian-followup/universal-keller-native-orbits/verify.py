#!/usr/bin/env python3
"""Selected exact native-orbit identities and fixed cubic descent certificates."""
from pathlib import Path
import ast,hashlib,json,platform,time
import sympy as s
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];checks=[]
def ck(label,v):
    assert v,label
    checks.append(label)
def zero(f):return s.cancel(f)==0
def main():
    begun=time.monotonic();x,y,z,lam=s.symbols('x y z lambda');xyz=[x,y,z];constructor=ROOT/'jacobian-followup/universal-keller-compression/verify.py';nodes=[n for n in ast.parse(constructor.read_text()).body if isinstance(n,ast.FunctionDef) and n.name=='source_map'];ck('reuse pure source constructor only',len(nodes)==1)
    env={'x':x,'y':y,'z':z};exec(compile(ast.Module(body=nodes,type_ignores=[]),str(constructor),'exec'),env);compiler=env['source_map'];K=s.Matrix(compiler(3,{}));t=1+x*y;q=t*t*z+y*y*(1+3*t);sig=t*q;X=x*q;LL=5-3*t-x*x*z;C=K[2]
    ck('actual base cubic Jacobian one',zero(K.jacobian(xyz).det()-1))
    ck('actual C has native x factor and order one',zero(C-x*LL) and LL.subs(x,0)==2)
    ck('native image of x-zero source plane is full target c-zero plane',K.subs(x,0)==s.Matrix([z+4*y*y,-y/2,0]))
    ck('q-zero component has nonconstant C residue',q.subs({y:0,z:0})==0 and C.subs({y:0,z:0})==2*x)
    # The whole nonlinear-r native family is a literal target shear.
    I=lambda a:a*(3*a-1)*(4*a+1);factor=lam*I(sig)**2;A=-3*factor*t**3*x*q*LL**3;B=2*factor*t**3*LL**3;v4=s.Matrix([-2*sig**2*X**2,-2*X**4]);v5=s.Matrix([-s.Rational(5,2)*sig**2*X**3,-3*X**5]);cor=A*v4+B*v5
    ck('actual nonlinear native controls equal the exact target shear',all(zero(v) for v in cor-s.Matrix([factor*sig**5*C**3,0])))
    actual=s.Matrix([sig,K[1]+factor*sig**5*C**3,C]);point={x:1,y:0,z:1}
    ck('literal nonlinear-r native Keller family Jacobian point one',actual.jacobian(xyz).subs(point).det()==1)
    ck('literal nonlinear controls have polynomial landing',all(s.denom(s.cancel(v))==1 for v in [A,B]))
    a,b,c,S,r=s.symbols('a b c S r');f=lam*a**5*c**3*I(a)**2;E=s.Matrix([a,b+f,c]);invE=s.Matrix([a,b-f,c])
    ck('target shear actual polynomial inverse exact',E.subs({a:invE[0],b:invE[1],c:invE[2]},simultaneous=True)==s.Matrix([a,b,c]))
    ck('target shear determinant one globally',E.jacobian([a,b,c]).det()==1)
    Cfield=(r+1)*S-a*S**3
    for m in [3,5]:
        af=-3*lam*a*Cfield**m/S**2;bf=2*lam*Cfield**m/S**3
        ck('native orbit effective r degree '+str(m),s.degree(af,r)==m and s.degree(bf,r)==m)
    # Native module membership is stronger than the source-jet condition.
    bad=2*sig**3*C**2/X**3;bad_slice=s.cancel(bad.subs({y:0,z:1}))
    ck('quadratic target shear fails polynomial native control',zero(bad_slice-2*(2-x*x)**2/x))
    ck('quadratic target shear native control has nonzero source residue',s.limit(x*bad_slice,x,0)==8)
    cshear=(sig**5*K[1]**3).subs(x,0)
    ck('a nonzero C target shear violates even native x contact',cshear!=0)
    # Exact marked-contact ideal, including both branches of max(5,ell).
    for ell in [1,2,6]:
        expected=a**max(5,ell)*(3*a-1)**ell*(4*a+1)**ell*c**3;got=s.lcm(a**5*c**3,I(a)**ell)
        ck('marked target orbit lcm formula order '+str(ell),s.Poly(got,a,c).monic()==s.Poly(expected,a,c).monic())
    ck('actual complete marked restrictions share first-output slices',all(I(v)==0 for v in [s.Integer(0),s.Rational(1,3),-s.Rational(1,4)]))
    ck('whole nonlinear correction divisible by the marked thickening',zero((actual[1]-K[1])/I(sig)**2-lam*sig**5*C**3))
    # All-degree leading-coefficient reductions and a sharp noncertificate.
    pp=s.Function('p')(S);kap,N=s.symbols('kappa N');vn=kap*pp**N
    ck('generic leading cancellation coefficient is constant over signal field',s.simplify(N*s.diff(pp,S)*vn-pp*s.diff(vn,S))==0)
    Bfield=(r-1-3*a*S*S)/(2*S);mixed=Bfield+lam*a**5*Cfield**3
    ck('nonlinear orbit has polynomial-r degree pair three and one',s.degree(mixed,r)==3 and s.degree(Cfield,r)==1)
    ck('elementary degree reduction removes the full nonlinear artifact',s.expand(mixed-lam*a**5*Cfield**3-Bfield)==0)
    toy=s.Matrix([S**2*r**2,S**3*r**3])
    ck('nondividing degree pair alone is not a Keller certificate',toy.jacobian([S,r]).det()==0)
    # Restricted stable profile moduli have the exact weight action.
    eta,theta,u=s.symbols('eta theta u',nonzero=True);profile=1+u*a+a*a
    action=lambda h,e:s.expand(e**-5*h.subs(a,e**-2*a))
    ck('profile weight action composes exactly',zero(action(action(profile,eta),theta)-action(profile,eta*theta)))
    ck('quartic signal coefficients have weights five seven nine',zero(action(profile,eta)-(eta**-5+u*a*eta**-7+a*a*eta**-9)))
    scaled=K.subs({x:eta*x,y:y/eta,z:z/eta**2},simultaneous=True)
    ck('stable source signal normalization is old-output recoverable',all(zero(v) for v in scaled-s.diag(eta**-2,eta**-1,eta)*K))
    # General quasi-finite saturation has a necessary hypothesis.
    ck('non-quasi-finite map permits forbidden rational descent',s.cancel((x*y)/x)==y)
    ck('its exceptional fibre is positive dimensional',x.subs(x,0)==0 and (x*y).subs(x,0)==0)
    # One fixed cubic field detects every polynomial old-output descent.
    FF=s.QQ.frac_field(a,b,c,lam);EP=s.Poly(a*S**3+b*S**2+S-c/2,S,domain=FF);rp=1+2*b*S+3*a*S*S
    def nf(expr):
        num,den=s.fraction(s.cancel(expr));pn=s.Poly(num,S,domain=FF);pd=s.Poly(den,S,domain=FF);ii=s.invert(pd,EP);return (pn*ii).rem(EP).as_expr()
    sinv=2*(1+b*S+a*S*S)/c;rinv=s.invert(s.Poly(rp,S,domain=FF),EP).as_expr()
    ck('generic cubic is degree three',EP.degree()==3)
    ck('generic S inverse exact',zero(nf(S*sinv)-1))
    ck('generic derivative inverse exact',zero(nf(rp*rinv)-1))
    ck('all derivative inverses reduce to fixed cubic dimension',s.degree(rinv,S)<=2)
    ck('old B reconstructs as constant target coefficient',zero(nf((rp-1-3*a*S*S)/(2*S))-b))
    C0F=(rp+1)*S-a*S**3;ck('old C reconstructs as constant target coefficient',zero(nf(C0F)-c))
    ck('a source coordinate is genuinely non-descended',zero(nf((1-rp)/S)+2*b+3*a*S))
    RR=nf(b+lam*a**5*C0F**3)
    ck('cubic reduction certifies nonlinear target-orbit descent',zero(RR-b-lam*a**5*c**3))
    quartic=nf(b-2*a**4*S*S)
    ck('cubic reduction rejects the genuinely new quartic covering degree',s.Poly(quartic,S).coeff_monomial(S**2)==-2*a**4)
    RF=a**4*S**3*(rp+s.Rational(1,3)-s.Rational(3,5)*a*S*S);dB=nf(-2*a**4*S*S*RF);dC=nf(-2*a**4*S**4*RF)
    ck('scalar source-dependent tangent fails exact old-field descent',any(s.Poly(v,S).coeff_monomial(S**j)!=0 for v in [dB,dC] for j in [1,2]))
    ck('non-descended tangent remains nonzero at generic parameter',any(s.diff(nf(lam*v),lam)!=0 for v in [dB,dC]))
    # The actual collision witness supports the field test via polynomial saturation.
    Rsrc=X**3*(q+sig/3)-s.Rational(3,5)*X**5;delta=s.Matrix([0,-2*sig**2*X**2*Rsrc,-2*X**4*Rsrc]);v0={x:0,y:0,z:-1};v1={x:-s.Rational(1,2),y:3,z:26}
    ck('actual same-old-output witness',list(K.subs(v0))==list(K.subs(v1))==[-1,0,0])
    ck('actual tangent separates the old-output fibre',list(delta.subs(v0))==[0,0,0] and list(delta.subs(v1))==[0,s.Rational(32,15),s.Rational(32,15)])
    Rind=a**4*S**3*(r+s.Rational(1,3)-s.Rational(3,5)*a*S*S)
    ck('formal tangent changes the C r coefficient excluded by profile orbits',s.diff(-2*a**4*S**4*Rind,r)!=0)
    ck('all literal maps use exact rational coefficients',not any(v.has(s.Float) for v in list(actual)+list(delta)))
    manifest=json.loads((HERE/'sources.json').read_text())
    for src in manifest['local_sources']:ck('local source pin '+src['path'],hashlib.sha256((ROOT/src['path']).read_bytes()).hexdigest()==src['sha256'])
    for src in manifest['external_sources']:
        blob=(HERE/src['path']).read_bytes();ck('primary source pin '+src['path'],len(blob)==src['bytes'] and hashlib.sha256(blob).hexdigest()==src['sha256'])
    receipt={'schema_version':1,'result':'pass','check_count':len(checks),'checks':checks,'python':platform.python_version(),'sympy':s.__version__,'runtime_seconds':round(time.monotonic()-begun,3),'old_factor_field_dimension':3,'nonlinear_orbit_degree_notice':'The selected nonlinear-r family has generic covering degree three because it is an explicit target automorphism of the cubic. r degree is not generic covering or total coordinate degree.','assurance':'Selected exact native target shears, whole-source control formulas, source-pole and marked-contact boundaries, all-degree leading reduction identities and exact cubic field certificates accompany written general proofs. Not independent whole-proof replay, formal verification, external review or publication-priority evidence. No canonical atomic receiver or all-factor/ unrestricted stable-equivalence classification.','sha256':{f:hashlib.sha256((HERE/f).read_bytes()).hexdigest() for f in ['PROOF.md','LITERATURE_AND_GAPS.md','verify.py','sources.json']}}
    (HERE/'receipt.json').write_text(json.dumps(receipt,indent=2,ensure_ascii=False)+'\n');print(json.dumps({'result':'pass','checks':len(checks),'seconds':receipt['runtime_seconds']}))
if __name__=='__main__':main()
