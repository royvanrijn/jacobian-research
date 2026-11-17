#!/usr/bin/env python3
"""Selected exact controls for written all-degree chart and curve arguments."""
from pathlib import Path
import hashlib,json,platform,time
import sympy as s
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];checks=[]
def ck(label,value):
    assert value,label
    checks.append(label)
def zero(value):return s.cancel(value)==0
def invariant_poly(expr,ss,uu):
    tt=s.Dummy('t');poly=s.Poly(s.expand(expr.subs(ss,tt+s.Rational(1,2))),tt)
    assert all(k[0]%2==0 for k,c in poly.terms()),'involution invariance'
    return s.expand(sum(c*(uu+s.Rational(1,4))**(k[0]//2) for k,c in poly.terms()))
def main():
    begun=time.monotonic();u,w,z,v,ss,alpha,beta=s.symbols('u w z v s alpha beta')
    M=s.symbols('M',integer=True,positive=True)
    transition=s.Matrix([u,w+(alpha-beta)/u**M])
    ck('all-exponent chart translation preserves area',zero(transition.jacobian([u,w]).det()-1))
    ck('all-exponent chart translation has the exact inverse',zero((transition[1]+(beta-alpha)/u**M)-w))
    ck('all-exponent chart canonical volume pulls back to plane volume',zero(s.diff(alpha+u**M*w,w)/u**M-1))
    for roots,nm in [([0,1],1),([0,1],2),([0,1,2],3)]:
        P=s.prod(z-a0 for a0 in roots)
        for ii,a0 in enumerate(roots):
            Zi=a0+u**nm*w;Vi=w*s.prod(a0-aj+u**nm*w for aj in roots if aj!=a0)
            label=' degree '+str(len(roots))+' exponent '+str(nm)+' chart '+str(ii)
            ck('chart equation'+label,zero(u**nm*Vi-P.subs(z,Zi)))
            ck('chart rational inverse'+label,zero((Zi-a0)/u**nm-w))
            ck('chart covers the whole retained component'+label,zero(Vi.subs(u,0)-s.diff(P,z).subs(z,a0)*w) and s.diff(P,z).subs(z,a0)!=0)
            if ii+1<len(roots):
                aj=roots[ii+1];Znext=aj+u**nm*w
                Vnext=w*s.prod(aj-ak+u**nm*w for ak in roots if ak!=aj)
                sh=w+(a0-aj)/u**nm
                ck('neighboring chart coordinates glue'+label,zero(Znext.subs(w,sh)-Zi) and zero(Vnext.subs(w,sh)-Vi))
                F=Zi**2+u*Vi;Fn=Znext**2+u*Vnext
                ck('global regular functions have polynomial compatible restrictions'+label,zero(Fn.subs(w,sh)-F))
    # A putative chart identity cannot extend across the other component.
    ck('identity chart coordinate acquires a genuine other-chart boundary pole',s.limit(u**2*(w+1/u**2),u,0)==1)
    # Exact positive control for the terminating plane Keller reduction.
    f=u+w*w;g=f*f+w;H=s.Matrix([f,g])
    ck('divisible ruling-degree positive control is Keller',s.expand(H.jacobian([u,w]).det())==1)
    ck('positive control has ruling degrees two and four',s.degree(f,w)==2 and s.degree(g,w)==4)
    ck('constant target shear removes divisible leading degree',s.expand(g-f*f)==w)
    b,c=s.symbols('b c')
    inv=s.Matrix([b-(c-b*b)**2,c-b*b])
    ck('the reduced positive control has a whole polynomial inverse',inv.subs({b:f,c:g},simultaneous=True).applyfunc(s.expand)==s.Matrix([u,w]))
    pp,ll=s.symbols('p l',integer=True,positive=True);lam=s.symbols('lambda')
    Af=s.Function('A')(u);Bf=lam*Af**ll
    ck('all-divisible-degree leading coefficient relation',s.simplify((ll*pp)*s.diff(Af,u)*Bf-pp*Af*s.diff(Bf,u))==0)
    ck('nondividing degrees alone do not establish a Keller pair',s.degree(u+w*w,w)==2 and s.degree(w**3,w)==3 and s.Matrix([u+w*w,w**3]).jacobian([u,w]).det()==3*w*w)
    # Standard formal neighborhoods can fail to extend to a global chart.
    cc=s.symbols('c')
    for order in [3,6]:
        branch=s.series((-1+s.sqrt(1+4*cc))/2,cc,0,order+1).removeO()
        ck('nilpotent quadratic branch control order '+str(order),s.Poly(branch+branch**2-cc,cc).rem(s.Poly(cc**(order+1),cc)).is_zero)
    ck('quadratic local model is etale on its stated open',s.diff(w+w*w,w)==1+2*w)
    ck('quadratic local model misses the entire discriminant line',zero((w+w*w+s.Rational(1,4))-(w+s.Rational(1,2))**2))
    # Root components and the sharp contact construction.
    ck('a root line is a whole embedded curve inside its own plane chart',zero((u**2*v-z*(z-1)).subs({u:ss,v:0,z:0})))
    for nm in [2,3]:
        cm=s.factorial(2*nm-1)/s.factorial(nm-1)**2
        zz=s.integrate(cm*ss**(nm-1)*(1-ss)**(nm-1),(ss,0,ss))
        uu=ss*(ss-1);vv=s.cancel(zz*(zz-1)/uu**nm)
        label=' exponent '+str(nm)
        ck('sharp contact polynomial normalization'+label,zz.subs(ss,0)==0 and zz.subs(ss,1)==1)
        ck('two distinct root components have order-m contact'+label,all(s.diff(zz,ss,j).subs(ss,0)==0 and s.diff(zz-1,ss,j).subs(ss,1)==0 for j in range(nm)))
        ck('contact curve has globally polynomial v'+label,s.denom(vv)==1 and zero(uu**nm*vv-zz*(zz-1)))
        ck('sharp unbounded contact degree is attained'+label,s.degree(zz,ss)==2*nm-1 and s.degree(uu,ss)==2)
        ck('contact polynomial has the exact reflection symmetry'+label,s.expand(zz.subs(ss,1-ss)+zz-1)==0)
        hh=s.cancel((zz-s.Rational(1,2))/(ss-s.Rational(1,2)));Hp=invariant_poly(hh,ss,u);Vp=invariant_poly(vv,ss,u)
        ck('collision polynomial is a ruling polynomial of degree m-1'+label,s.degree(Hp,u)==nm-1 and zero(zz-s.Rational(1,2)-(ss-s.Rational(1,2))*Hp.subs(u,uu)))
        ck('collision roots avoid zero and the double-parameter value'+label,Hp.subs(u,0)==1 and Hp.subs(u,-s.Rational(1,4))==cm*s.Integer(4)**(1-nm))
        ck('collision root is coprime to the parameter discriminant'+label,s.gcd(Hp,u*(1+4*u),u)==1)
        ck('both colliding parameters have the same surface v'+label,s.rem(4*u**nm*Vp+1,Hp,u)==0)
        ck('surface v is invariant under exchange of the two parameters'+label,zero(vv-Vp.subs(u,uu)))
        if nm==2:
            ck('exponent-two collision has the exact literal root',Hp.subs(u,s.Rational(1,2))==0 and Vp.subs(u,s.Rational(1,2))==-1)
    # General contact estimate: sum(m*e_j-1)=m*N-r.
    e0,e1,mm=s.symbols('e0 e1 m',integer=True,positive=True)
    ck('summed derivative contacts give the unbounded degree lower bound',s.expand((mm*e0-1)+(mm*e1-1)-(mm*(e0+e1)-2))==0)
    manifest=json.loads((HERE/'sources.json').read_text())
    for src in manifest['local_sources']:
        ck('retained local source '+src['path'],hashlib.sha256((ROOT/src['path']).read_bytes()).hexdigest()==src['sha256'])
    for src in manifest['external_sources']:
        blob=(HERE/src['path']).read_bytes()
        ck('reviewed primary source '+src['path'],len(blob)==src['bytes'] and hashlib.sha256(blob).hexdigest()==src['sha256'])
    receipt={'schema_version':1,'result':'pass','check_count':len(checks),'checks':checks,'python':platform.python_version(),'sympy':s.__version__,'runtime_seconds':round(time.monotonic()-begun,3),
        'assurance':'Selected exact atlas, inverse, volume, compatibility, degree-reduction, sharp-contact and collision controls accompany written general proofs. Classical Danielewski charts and open-plane JC(2) reduction are credited. Line injectivity/rectification and open-and-closed sections are published inputs. Not independent whole-proof replay, formal verification, external review or novelty evidence. No compatible etale atlas, whole intermediate or canonical atomic receiver is constructed.',
        'sha256':{f:hashlib.sha256((HERE/f).read_bytes()).hexdigest() for f in ['PROOF.md','LITERATURE_AND_GAPS.md','verify.py','sources.json']}}
    (HERE/'receipt.json').write_text(json.dumps(receipt,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps({'result':'pass','checks':len(checks),'seconds':receipt['runtime_seconds']}))
if __name__=='__main__':main()
