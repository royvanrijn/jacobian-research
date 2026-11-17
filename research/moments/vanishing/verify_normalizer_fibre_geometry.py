#!/usr/bin/env python3
"""Exact normalizer fibre maps, j-divisors, and derived-series controls.

Universal birational descent statements are written proofs, not a finite search.
Requires Sympy; imports existing exact binary-polyhedral group arithmetic.
"""
from pathlib import Path
from collections import deque,Counter
import hashlib,json
import sympy as S
from verify_degree12_four_weights import groups

ROOT=Path(__file__).resolve().parent

def derived_series(name):
 field,gens,whole=groups(name);add,mul,scale=field.add,field.mul,field.scale
 identity=(field.one,field.zero,field.zero,field.one)
 minus=tuple(scale(a,-1) for a in identity)
 def mm(a,b):
  return (add(mul(a[0],b[0]),mul(a[1],b[2])),add(mul(a[0],b[1]),mul(a[1],b[3])),
          add(mul(a[2],b[0]),mul(a[3],b[2])),add(mul(a[2],b[1]),mul(a[3],b[3])))
 def inv(a):return a[3],scale(a[1],-1),scale(a[2],-1),a[0]
 def comm(a,b):return mm(mm(mm(a,b),inv(a)),inv(b))
 def closure(generating):
  out={identity};queue=deque([identity])
  while queue:
   a=queue.popleft()
   for g in generating:
    b=mm(a,g)
    if b not in out:out.add(b);queue.append(b)
  return out
 def derived(generating,group):
  dgens=set(comm(a,b) for a in generating for b in generating)-{identity}
  while True:
   out=closure(dgens)
   extra={mm(mm(g,a),inv(g)) for a in dgens for g in generating}-out
   if not extra:break
   dgens|=extra
  assert out<=group
  assert all(mm(mm(g,a),inv(g)) in out for a in dgens for g in generating)
  assert all(comm(a,b) in out for a in generating for b in generating)
  return list(dgens),out
 series=[len(whole)];current=whole;currentgens=gens;levels=[whole]
 while len(current)>1:
  newgens,new=derived(currentgens,current)
  series.append(len(new));levels.append(new)
  if new==current:break
  currentgens,current=newgens,new
 assert series==([48,24,8,2,1] if name=='2O' else [120,120])
 second=levels[2] if name=='2O' else levels[1]
 order4=[a for a in second if mm(a,a)==minus]
 pair=next((a,b) for a in order4 for b in order4 if comm(a,b)==minus)
 q8=closure(pair);assert len(q8)==8 and q8<=second
 orders=Counter()
 for a in q8:
  power=identity
  for k in range(1,9):
   power=mm(power,a)
   if power==identity:orders[k]+=1;break
 assert dict(orders)=={1:1,2:1,4:6}
 return {'derived_series_orders':series,'second_derived_contains_Q8':True,
         'Q8_element_orders':dict(sorted(orders.items()))}

def main():
 u,v,a,t,b,k,x=S.symbols('u v a t b k x')
 U,V,T=S.symbols('U V T')
 circuit=(1+U)*(V-T*(U+2)*(2-T*(U**2+2*U+2)))
 compact=(1+U)*(1-(1+U)**2*T)**2/U
 assert S.cancel(circuit.subs(V,(1-T)**2/U)-compact)==0
 # On the right-T quotient u*v+a^2=1.
 phase=S.cancel(compact.subs({U:u**2,T:a**2}))
 beta=(1+u**2)*a;rho=1-beta**2;kappa=rho/u
 assert S.cancel(kappa**2-phase+(1-beta**2)**2)==0
 curve=k**2+(1-b**2)**2-t
 def zero_on_curve(expr):
  numerator,denominator=S.fraction(S.cancel(expr))
  domain=S.QQ.frac_field(b,t)
  n=S.Poly(numerator,k,domain=domain).rem(S.Poly(curve,k,domain=domain))
  d=S.Poly(denominator,k,domain=domain).rem(S.Poly(curve,k,domain=domain))
  assert not n.as_expr() and d.as_expr()!=0
 inverse_u=(1-b**2)/k;inverse_a=b*k**2/t
 zero_on_curve(phase.subs({u:inverse_u,a:inverse_a},simultaneous=True)-t)
 zero_on_curve(beta.subs({u:inverse_u,a:inverse_a},simultaneous=True)-b)
 zero_on_curve(kappa.subs({u:inverse_u,a:inverse_a},simultaneous=True)-k)
 # Forward round-trip uses the actual phase value, without specializations.
 assert S.cancel(inverse_u.subs({b:beta,k:kappa},simultaneous=True)-u)==0
 assert S.cancel(inverse_a.subs({b:beta,k:kappa,t:phase},simultaneous=True)-a)==0
 # Recover the full SL2 field by a free right-torus coordinate x.
 y=-u/(2*x);z=x*(1-a)/u;w=(1+a)/(2*x)
 assert S.cancel(x*w-y*z-1)==0
 assert S.cancel(-2*x*y-u)==0 and S.cancel(2*z*w-(1-a**2)/u)==0
 assert S.cancel(x*w+y*z-a)==0
 # The original, non-Weyl Hopf phase has rational rather than elliptic fibres.
 hopf=(1+u)*(1-(1+u)**2*a**2)/u
 hopf_beta=(1+u)*a
 hu=(1-b**2)/(t-1+b**2);ha=b*(t-1+b**2)/t
 assert S.cancel(hopf.subs({u:hu,a:ha},simultaneous=True)-t)==0
 assert S.cancel(hopf_beta.subs({u:hu,a:ha},simultaneous=True)-b)==0
 assert S.cancel(hu.subs({b:hopf_beta,t:hopf},simultaneous=True)-u)==0
 assert S.cancel(ha.subs({b:hopf_beta,t:hopf},simultaneous=True)-a)==0
 # A generic right translate of x/z is a rational left-Ga slice.
 xx,yy,zz,ww,aa,cc=S.symbols('xx yy zz ww aa cc')
 translated=(xx*aa+yy*cc)/(zz*aa+ww*cc)
 assert S.cancel(zz*S.diff(translated,xx)+ww*S.diff(translated,yy)-1)==0
 identity=translated-xx/zz+cc/(zz*(zz*aa+ww*cc))
 num=S.fraction(S.cancel(identity))[0]
 assert S.rem(num,S.Poly(xx*ww-yy*zz-1,xx),xx)==0
 quartic=(1-b**2)**2-t
 discriminant=S.factor(S.discriminant(quartic,b))
 I=16-12*t;J=144*t-128
 assert S.expand(discriminant-256*t**2*(1-t))==0
 assert S.expand(4*I**3-J**2-27*discriminant)==0
 A=-27*I;B=-27*J
 j=S.factor((-48*A)**3/(-16*(4*A**3+27*B**2)))
 assert S.cancel(j-64*(4-3*t)**3/(t**2*(1-t)))==0
 numerator,denominator=S.fraction(j)
 assert S.factor(numerator)==64*(3*t-4)**3
 assert S.factor(denominator)==t**2*(t-1)
 assert S.limit(j,t,S.oo)==1728
 # These three unique marked points force every Mobius deck map of j to be1.
 divisor={'poles':{'0':2,'1':1},'zeros':{'4/3':3},'j_at_infinity':1728}
 assert S.cancel(j-j.subs(t,-t)-128*(81*t**2-80)/(t*(1-t**2)))==0
 record={'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
  'sympy_version':S.__version__,'generic_curve':'k^2=t-(1-b^2)^2',
  'inverse_coordinates':{'u':'(1-b^2)/k','tau':'b*k^2/t'},
  'quartic_polynomial_discriminant':str(discriminant),'I':str(I),'J':str(J),
  'j':str(j),'j_divisors':divisor,
  'generic_birational_maps_and_matrix_chart_verified':True,
  'hopf_rational_coordinate_maps_verified':True,
  'generic_translated_slice_derivative_and_average_formula_verified':True,
  'group_controls':{name:derived_series(name) for name in ['2O','2I']},
  'imported_implementations_sha256':{n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in [
   'verify_degree12_four_weights.py','verify_haar_shear_witnesses.py','verify_octahedral_minimal_type.py']},
  'completed':True,'universal_birational_group_proof_machine_verified':False,
  'all_exceptional_MZ_decisions':'OPEN','independent_replay':False,'external_review':False,'formal_verification':False}
 (ROOT/'normalizer_fibre_verification.json').write_text(json.dumps(record,indent=2)+'\n')
 print(json.dumps(record,indent=2))

if __name__=='__main__':main()
