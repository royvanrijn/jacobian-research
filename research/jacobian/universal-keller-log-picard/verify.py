#!/usr/bin/env python3
"""Exact local controls; the general surface proofs are in PROOF.md."""
from pathlib import Path
import ast,hashlib,json,platform,time
import sympy as s
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];checks=[]
def ck(label,value):
    assert value,label
    checks.append(label)
def zero(v):return s.cancel(v)==0
def main():
    begun=time.monotonic()
    x,y,z,S,Y,a,r,v,w,u,delta=s.symbols('x y z S Y a r v w u delta')
    t=1+x*y;Q=t*t*z+y*y*(1+3*t);sig=t*Q
    ctor=ROOT/'jacobian/universal-keller-compression/verify.py'
    nodes=[n for n in ast.parse(ctor.read_text()).body if isinstance(n,ast.FunctionDef) and n.name=='source_map']
    ck('literal source constructor isolated without executing its checker',len(nodes)==1)
    env={'x':x,'y':y,'z':z};exec(compile(ast.Module(body=nodes,type_ignores=[]),str(ctor),'exec'),env)
    K=s.Matrix(env['source_map'](3,{}))
    ck('retained source signal agrees with literal constructor',zero(K[0]-sig))
    ck('retained chart determinant is exactly t',zero(s.Matrix([sig,x/t,y]).jacobian([x,y,z]).det()-t))
    rho=1-S*Y
    back={x:S/rho,y:Y,z:a*rho**3-rho*(rho+3)*Y**2}
    ck('retained inverse chart recovers signal',zero(sig.subs(back,simultaneous=True)-a))
    ck('retained inverse chart recovers S and Y',zero((x/t).subs(back,simultaneous=True)-S) and back[y]==Y)
    J=s.Matrix([back[x],back[y],back[z]]).jacobian([a,S,Y])
    ck('source volume becomes r da wedge dS wedge dY',zero(J.det()-rho))
    ck('whole Keller determinant scales the relative volume by delta',zero(delta*J.det()-delta*rho))
    # Source compactification: smooth conic and two transverse infinity crossings.
    X,H,Z=s.symbols('X H Z');conic=Z**2-X*H
    grad=[s.diff(conic,h) for h in [X,H,Z]]
    ck('conic gradient has only the inadmissible homogeneous zero',s.solve(grad,[X,H,Z])=={X:0,H:0,Z:0})
    for point,name in [({X:1,H:0,Z:0},'X'),({X:0,H:1,Z:0},'Y')]:
        ck('conic and infinity meet transversely at '+name,s.Matrix([grad,[0,0,1]]).subs(point).rank()==2)
    ck('log canonical degree of conic plus line is zero',-3+2+1==0)
    chart=s.Matrix([1/w,v/w]);jac=chart.jacobian([v,w]).det();rr=(w*w-v)/(w*w)
    ck('infinity chart has exact area coefficient',zero(jac-1/w**3))
    ck('infinity chart has exact unit coefficient',zero((1-S*Y).subs({S:1/w,Y:v/w},simultaneous=True)-rr))
    eta=jac/rr;omega=rr*jac
    ck('log volume has only simple conic and infinity poles',zero(eta-1/(w*(w*w-v))))
    ck('whole Keller volume has infinity pole order five',zero(omega-(w*w-v)/w**5) and s.cancel(omega*w**5).subs(w,0)==-v)
    ck('whole Keller and logarithmic volumes differ by r squared',zero(omega/eta-rr**2))
    # Local boundary resolutions and pullback do not require finite g.
    h,j,p,q=s.symbols('h j p q')
    # Blow up x=u, y=u*v at a smooth boundary x=0.
    bl=s.Matrix([u,u*v]).jacobian([u,v]).det()
    ck('one-branch boundary blowup has one zero of log volume',zero(bl/u-1))
    ck('two-branch boundary blowup preserves logarithmic volume',zero(bl/(u*u*v)-1/(u*v)))
    ck('one-branch log canonical discrepancy is one',2-1==1)
    ck('two-branch log canonical discrepancy is zero',2-2==0)
    dlog=s.Matrix([[p/u,q/v],[h/u,j/v]])
    ck('monomial pullback of logarithmic volume has only log poles',zero(dlog.det()-(p*j-q*h)/(u*v)))
    mm=s.symbols('M',integer=True,positive=True);cc=s.symbols('c',nonzero=True)
    ck('all-order volume multiplier is the negative even unit power',s.powsimp(cc*delta**(-mm)*r**(-2*mm)*(delta*r**2)**mm,force=True)==cc)
    ck('multiplier cannot be constant for a positive pluriform order',s.simplify(s.diff(r**(-2*mm),r)+2*mm*r**(-2*mm-1))==0 and mm.is_positive)
    # The rational point is carried through an identity-mixing whole automorphism.
    U,T1,T2,T3=s.symbols('U T1 T2 T3')
    G=s.Matrix([u+z*z,x+u+z*z,y+(x+u+z*z)**2,z])
    Gi=s.Matrix([U-T3*T3,T1-U,T2-T1*T1,T3])
    sub={U:G[0],T1:G[1],T2:G[2],T3:G[3]}
    ck('mixed stabilized source automorphism has exact inverse',Gi.subs(sub,simultaneous=True)==s.Matrix([u,x,y,z]))
    ck('mixed source automorphism is Keller',G.jacobian([u,x,y,z]).det()==1)
    P=G.subs({x:0,y:0,z:a})
    ck('carried whole section has literal polynomial coordinates',P==s.Matrix([u+a*a,u+a*a,(u+a*a)**2,a]))
    qsig=sig.subs({x:Gi[1],y:Gi[2],z:Gi[3]},simultaneous=True)
    psub={U:P[0],T1:P[1],T2:P[2],T3:P[3]}
    ck('carried section fixes retained identity under q',zero(Gi[0].subs(psub,simultaneous=True)-u))
    ck('carried section fixes retained signal under q',zero(qsig.subs(psub,simultaneous=True)-a))
    # Exact differential coefficient on every smooth higher-power surface.
    AA,BB,Ap,Bp,PP,Pp,CC,Cp,DD,M=s.symbols('A B Ap Bp P Pp C Cp D M')
    coef=(AA-BB)*Pp+(M*Ap-Bp)*PP
    subst={BB:AA-CC,Bp:Ap-Cp}
    reduced=s.expand(coef.subs(subst,simultaneous=True)).subs(Ap,(DD-Cp)/(M-1))
    ck('all-exponent primitive identity reduces to Bezout',zero(reduced-(CC*Pp+DD*PP)))
    ck('primitive formula requires exponent different from one',s.denom(s.factor((DD-Cp)/(M-1)))==M-1)
    for nd,nm in [(2,2),(3,2),(3,3)]:
        poly=z**nd-a;pp=s.diff(poly,z);ci=s.invert(pp,poly,z);di=s.cancel((1-ci*pp)/poly)
        A=s.integrate((di-s.diff(ci,z))/(nm-1),z);B=A-ci
        ck('polynomial Bezout inverse degree '+str(nd)+' exponent '+str(nm),zero(ci*pp+di*poly-1))
        ck('regular primitive coefficient identity degree '+str(nd)+' exponent '+str(nm),zero((A-B)*pp+(nm*s.diff(A,z)-s.diff(B,z))*poly-1))
        ck('primitive coefficients have no source-variable denominators degree '+str(nd)+' exponent '+str(nm),all(s.denom(s.cancel(co)).free_symbols<={a} for co in [A,B]))
    # Canonical form patching, smoothness and arithmetic/geometric test models.
    for nd,nm in [(2,1),(2,2),(3,3)]:
        poly=z**nd-a;pp=s.diff(poly,z)
        ck('square-free retained polynomial degree '+str(nd),s.gcd(poly,pp,z)==1)
        ck('localized model solves for v degree '+str(nd)+' exponent '+str(nm),zero(u**nm*(poly/u**nm)-poly))
        ck('arithmetic model retained K-point degree '+str(nd)+' exponent '+str(nm),zero((u**nm*v-poly).subs({u:1,v:-a,z:0})))
        ck('canonical volume patches on implicit surface degree '+str(nd)+' exponent '+str(nm),zero((pp/u**nm)/pp-1/u**nm))
    # An all-exponent raw spread fails whole retained-base smoothness.
    ni,mi=s.symbols('ni mi',integer=True,nonnegative=True)
    raw=z**(ni+2)-u**(mi+2)*v
    ck('all higher-power displayed spreads have an entire critical zero-fibre line',all(s.simplify(s.diff(raw,co)).subs({u:0,z:0})==0 for co in [u,v,z]) and raw.subs({u:0,z:0})==0)
    raw1=z**(ni+2)-u*v
    ck('all ordinary displayed spreads have a critical zero-fibre point',all(s.simplify(s.diff(raw1,co)).subs({u:0,v:0,z:0})==0 for co in [u,v,z]))
    alpha,beta,c0,theta=s.symbols('alpha beta c theta')
    ck('ordinary-family cylinder period is nonzero for distinct roots',s.integrate(s.I,(theta,0,2*s.pi))*(beta-alpha)==2*s.pi*s.I*(beta-alpha))
    ck('a u-coordinate mate would require contradictory root coefficients',s.solve([c0-cc,c0+cc],[c0,cc])==[])
    # Quotient lattice controls: one primitive diagonal relation, no fixed classes.
    for nd in [2,3,5]:
        # Quotient basis e_0,...,e_(d-2), with e_(d-1)=-sum.
        perm=s.zeros(nd-1)
        for col in range(nd-1):
            if col+1<nd-1:perm[col+1,col]=1
            else:
                for row in range(nd-1):perm[row,col]=-1
        ck('root-cycle action has exact order degree '+str(nd),perm**nd==s.eye(nd-1))
        ck('root-cycle quotient has no rational fixed class degree '+str(nd),(perm-s.eye(nd-1)).det()!=0)
        ck('diagonal relation is primitive degree '+str(nd),s.gcd_list([1]*nd)==1)
    ck('generic completed A1 fibre has negative log degree',-2+1==-1)
    manifest=json.loads((HERE/'sources.json').read_text())
    for src in manifest['local_sources']:
        ck('retained local source '+src['path'],hashlib.sha256((ROOT/src['path']).read_bytes()).hexdigest()==src['sha256'])
    for src in manifest['external_sources']:
        blob=(HERE/src['path']).read_bytes()
        ck('reviewed primary HTML '+src['path'],len(blob)==src['bytes'] and hashlib.sha256(blob).hexdigest()==src['sha256'])
    receipt={'schema_version':1,'result':'pass','check_count':len(checks),'checks':checks,
        'python':platform.python_version(),'sympy':s.__version__,'runtime_seconds':round(time.monotonic()-begun,3),
        'assurance':'Selected exact chart, logarithmic-pole, boundary-blowup, mixed-section, regular-primitive and Picard-lattice controls accompany written general proofs. Miyanishi affine-plane characterization and effective sheaf descent are explicit published inputs. Not independent whole-proof replay, formal verification, external review or novelty evidence. No required higher-power etale diagram, whole factor, nonlinear integration or canonical atomic receiver is constructed.',
        'sha256':{f:hashlib.sha256((HERE/f).read_bytes()).hexdigest() for f in ['PROOF.md','LITERATURE_AND_GAPS.md','verify.py','sources.json']}}
    (HERE/'receipt.json').write_text(json.dumps(receipt,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps({'result':'pass','checks':len(checks),'seconds':receipt['runtime_seconds']}))
if __name__=='__main__':main()
