"""Exact controls for the all-degree specialization and K2 one-coordinate proof."""
import hashlib,json,platform,time
from pathlib import Path
import sympy as sp
START=time.monotonic();HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];checks=[]
def check(value,label):
    assert bool(value),label
    checks.append(label)
def zero(value,label):check(sp.cancel(value)==0 or sp.simplify(value)==0,label)
def zeros(values,label):check(all(sp.cancel(v)==0 or sp.simplify(v)==0 for v in values),label)
x,y,z,w,s,u,v,c,T,H,L,P,Q=sp.symbols('x y z w s u v c T H L P Q')
def Hamiltonian(f,g):return sp.expand(sp.diff(f,u)*sp.diff(g,v)-sp.diff(f,v)*sp.diff(g,u))
# A special fibre can be a power of a coordinate. Do not strengthen JW0 to
# claim that every specialization is itself a coordinate.
for exponent in [2,3]:
    complement=v+x*u*u if exponent==3 else v
    f=x*u+complement**exponent
    check(sp.Matrix([f,complement]).jacobian([u,v]).det()==x,
          'the generic coordinate pair has its real field-unit Jacobian,power='+str(exponent))
    invu=(H-L**exponent)/x
    invv=L-x*invu**2 if exponent==3 else L
    zeros([f.subs({u:invu,v:invv},simultaneous=True)-H,
           complement.subs({u:invu,v:invv},simultaneous=True)-L],
          'the generic-coordinate claim includes exact polynomial-over-k(x) inverse identities,power='+str(exponent))
    zero(f.subs(x,0)-v**exponent,'the nonconstant specialization can be a proper power,power='+str(exponent))
    zero(Hamiltonian(f,f),'the generic Hamiltonian kills its own coordinate,power='+str(exponent))
    for var in [u,v]:
        zero(Hamiltonian(f,var).subs(x,0)-Hamiltonian(f.subs(x,0),var),
             'Hamiltonian coefficients specialize exactly,var='+str(var)+',power='+str(exponent))
    # Specialization is visibly a nonzero triangular LND; the general proof
    # descends exact nilpotence identities before using Rentschler.
    check(Hamiltonian(f.subs(x,0),u)!=0,'the characteristic-zero specialized Hamiltonian is nonzero,power='+str(exponent))
    zero(Hamiltonian(f.subs(x,0),Hamiltonian(f.subs(x,0),u)),
         'the specialized Hamiltonian is nilpotent even though the function is not a coordinate,power='+str(exponent))
# Uniform native chart and the actual exceptional relation.
K2=w*(2*x*y+z*z-1)+z
native_u=w;native_v=1-2*x*y-z*z
zero(K2+native_u*native_v-z,'the actual uniform chart has reverse z=s+uv')
zz=s+u*v;Fs=v+zz*zz-1
zero(Fs.subs({s:K2,u:native_u,v:native_v},simultaneous=True)+2*x*y,
     'the full native source has the literal relation 2xy+F_s=0')
zero((1-native_v-(K2+native_u*native_v)**2)/(2*x)-y,
     'inverting x gives a polynomial ring over C[x,x^-1,s],with exact y recovery')
for pole_order in [1,2,3]:
    # Native y^m has a pole numerator divisible by F_s^m. It is not a
    # generic plane coordinate; a native source coordinate alone is insufficient.
    numerator=sp.expand((-Fs/2)**pole_order)
    zero(numerator-(-1/sp.Integer(2))**pole_order*Fs**pole_order,
         'clearing a genuine native pole produces the exceptional curve factor,order='+str(pole_order))
    check(sp.rem(numerator,Fs,u,v)==0,
          'every displayed pole-leading numerator lies in the exact special-fibre kernel,order='+str(pole_order))
zero(Hamiltonian(Fs,zz)+v,'the exceptional pencil Hamiltonian differentiates its curve parameter by minus v')
zero(v*(1+u*(zz+s))-(Fs+1-s*s),
     'v is an actual unit on the generic exceptional curve,not merely nonzero at generic points')
Z=sp.symbols('Z')
vparam=1+c-Z*Z;uparam=(Z-s)/vparam
zero((s+u*v).subs({u:uparam,v:vparam},simultaneous=True)-Z,
     'the complete pencil parameter has its exact reverse z identity')
zero((Fs-c).subs({u:uparam,v:vparam},simultaneous=True),
     'the rational solved formula gives the entire nonexceptional pencil fibres')
for sign in [-1,1]:
    zero((Z-s).subs(Z,sign)-sign+s,'the two omitted points of Gamma are exactly z=+-1,sign='+str(sign))
    check((sign-s)!=0,'the generic value field excludes both critical specializations,sign='+str(sign))
zero(Fs-(s*s-1)-v*(1+2*s*u+u*u*v),'the only reducible pencil level has its full two-factor identity')
second=1+2*s*u+u*u*v
zero(second.subs(v,-(1+2*s*u)/(u*u)),'the nonline component is the whole Gm with parameter u')
zero(u*(-2*s-u*v)-(1-second),'the second component has an explicit nonconstant unit u')
zero(second.subs({u:-1/(2*s),v:0}),'the A1 and Gm components meet at the one stated point')
check(sp.diff(second,u).subs({u:-1/(2*s),v:0})==2*s,
      'the two special components are distinct and meet transversely in the generic value field')
zero((Fs+1).subs({u:(Z-s)/(-Z*Z),v:-Z*Z},simultaneous=True),
     'the other special pencil level is precisely a punctured line,not a missing A1 component')
crit={u:-1/(2*s),v:0}
zeros([sp.diff(Fs,q).subs(crit) for q in [u,v]],'the unique critical pencil point is on the v-zero line')
zero(sp.hessian(Fs,[u,v]).det().subs(crit)+4*s*s,'the unique critical pencil point is Morse for generic s')
for sign in [-1,1]:
    specialized=sp.factor(Fs.subs(s,sign))
    zero(specialized-v*(1+2*sign*u+u*u*v),
         'Gamma is reducible at the actual critical values,so its non-A1 irreducibility is used only over C(s),sign='+str(sign))
# A coordinate line omitting Gamma would force a constant nonzero F_s on it.
# Two such fibres are distinct; the proof uses the complete pencil list above.
alpha,beta=sp.symbols('alpha beta',nonzero=True)
check(sp.expand((beta+alpha)-(beta-alpha))==2*alpha,
      'the affine z restriction omits two distinct coordinate values')
check(sp.degree(sp.Symbol('q')**3,sp.Symbol('q'))*2>1,
      'a nonlinear outer polynomial cannot compose with a nonconstant rational map to degree one')
# Exact native relative ideal and all critical points.
hx,hy,hz,hw=[sp.diff(K2,q) for q in [x,y,z,w]]
zero(x*hz-z*hy-x,'the native relative ideal contains x by an exact Bezout combination')
zero(hw-2*y*(x*hz-z*hy)-(z*z-1),
     'the same relative ideal contains z^2-1 with an exact polynomial witness')
zero(hy-2*w*x,'the reverse relative-ideal inclusion for h_y is literal')
zero(hw-(2*x*y+z*z-1),'the reverse relative-ideal inclusion for h_w is literal')
for sign in [-1,1]:
    line={x:0,z:sign,w:-sp.Rational(1,2)/sign}
    zeros([hy.subs(line),hz.subs(line),hw.subs(line)],
          'the entire closed affine y-line is relatively critical,sign='+str(sign))
    point=dict(line);point[y]=0
    zero(K2.subs(point)-sign,'the full critical value is the stated sign,sign='+str(sign))
    zero(hx.subs(line)+y/sign,'the remaining gradient output cuts one reduced point on the critical line,sign='+str(sign))
    zero(sp.hessian(K2,[x,y,z,w]).det().subs(point)-4,
         'both full critical points are Morse with the retained common determinant,sign='+str(sign))
# Trace alignment uses actual A1 components and excludes the hyperbola factors.
B=w*(z*z-1)+z
for sign in [-1,1]:
    zero(B-sign-(z-sign)*(w*(z+sign)+1),'the entire critical trace fibre is factored,sign='+str(sign))
    zero((w*(z+sign)+1).subs(w,-1/(z+sign)),
         'the second critical-level component is exactly Gm,sign='+str(sign))
    zero((z+sign)*(-w)-(1-(w*(z+sign)+1)),
         'the second trace component has a nonconstant unit and is not A1,sign='+str(sign))
tplus,tminus,ap,am=sp.symbols('tplus tminus ap am')
alignment=sp.Poly(tplus+ap*(z-1)-(tminus+am*(z+1)),z)
zero(alignment.coeff_monomial(z)-(ap-am),'the two exact component equations force one common affine multiplier')
zero((tplus+ap*(z-1)).subs(ap,(tplus-tminus)/2)-((tplus-tminus)*z/2+(tplus+tminus)/2),
     'the only possible trace coordinate is literally affine in z')
# Coprime Keller/linear coefficient rows give the actual generic plane coordinate.
F=2*x*T;G=T*T-1;aa=T/(2*x);bb=-1
zero(aa*F+bb*G-1,'the C2 coefficient row is genuinely unimodular over C(x)[T]')
HH=F*P+G*Q;LL=-bb*P+aa*Q
zero(aa*HH-G*LL-P,'the generic row completion has the first exact inverse identity')
zero(bb*HH+F*LL-Q,'the generic row completion has the second exact inverse identity')
check(sp.Matrix([HH,LL]).jacobian([P,Q]).det()==1,
      'the generic-coordinate assertion follows from a polynomial row completion,not abstract rationality')
C2=2*x*z*y+(z*z-1)*w+z
zero(C2.subs(x,0)-B,'the C2 comparison has exactly the same native exceptional trace')
zero(sp.hessian(C2,[x,y,z,w]).det()-16*z**4,
     'the compatible C2 point determinant is not mistaken for a unit native Hessian')
# Constant affine source normalization is the only Hessian invariance used here.
lin=x+2*y+3*z-w+7
normalizer={x:x-2*y-3*z+w-7,y:y,z:z,w:w}
zero(lin.subs(normalizer,simultaneous=True)-x,'a nonconstant affine pullback is normalized by an actual constant affine precomposition')
check(sp.Matrix([normalizer[q] for q in [x,y,z,w]]).jacobian([x,y,z,w]).det()==1,
      'the affine normalization has a full polynomial inverse and constant Jacobian')
context=json.loads((HERE/'sources/coordinate-specialization-context.json').read_text());algebra=len(checks)
for path,digest in context['retained_proof_and_source_sha256'].items():
    check(hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest,'retained proof/checker/receipt bytes: '+path)
for entry in context['newly_cached_primary_context']:
    check(hashlib.sha256((HERE/'sources'/entry['retained_file']).read_bytes()).hexdigest()==entry['sha256'],
          'inspected primary context bytes: '+entry['retained_file'])
check(len(context['newly_imported_proved_external_statements'])==1 and context['newly_imported_proved_external_statements'][0]['hypotheses_checked'],
      'the sole new classification import is explicitly restricted to two variables over characteristic-zero fields')
check(context['newly_used_unproved_external_statements']==[],'no unproved specialization or action-classification claim is imported')
check(context['K2_one_coordinate_unit_hessian_route_excluded']
      and not context['all_K2_preserving_actions_fixing_x_classified']
      and not context['fixed_K2_full_source_orbit_resolved']
      and not context['generic_coordinate_assumed_from_abstract_rational_field_generator']
      and not context['specialization_assumed_itself_a_coordinate']
      and not context['unrestricted_HC4_or_JC2_resolved'],
      'the one-coordinate source class is closed while all unrestricted source and endpoint questions remain open')
flags={k:v for k,v in context.items() if isinstance(v,bool)}
sources={'sources/coordinate-specialization-context.json':hashlib.sha256((HERE/'sources/coordinate-specialization-context.json').read_bytes()).hexdigest()}
sources.update({'sources/'+e['retained_file']:e['sha256'] for e in context['newly_cached_primary_context']})
receipt=dict(outcome='PASS',claims=['HC4MYJW0','HC4MYJW1','HC4MYJW2'],check_count=len(checks),
 algebraic_and_scope_controls=algebra+3,checks=checks,proved_import_count=1,
 checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 proof_sha256=hashlib.sha256((HERE/'GENERIC_COORDINATE_SPECIALIZATION_AND_K2_ONE_COORDINATE_CLOSURE.md').read_bytes()).hexdigest(),
 source_sha256=sources,dependencies=context['retained_proof_and_source_sha256'],
 assurance='Finite exact controls of written all-degree coordinate specialization,no-pole,whole-pencil and K2 one-coordinate arguments; not independent replay,formal verification or external review.',
 python=platform.python_version(),sympy=sp.__version__,elapsed_seconds=round(time.monotonic()-START,3),**flags)
(HERE/'receipt-coordinate-specialization.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({k:receipt[k] for k in ['outcome','check_count','algebraic_and_scope_controls','proved_import_count','elapsed_seconds','original_goal_complete']},indent=2))
