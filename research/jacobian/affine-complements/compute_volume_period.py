"""Universal two-chart period identities for the forced F1 volume form.

The global exactness obstruction and family coverage are written in
VOLUME_PERIOD_OBSTRUCTION.md; the remaining exact form is not a finite map.
"""
from pathlib import Path
import sympy as S,json,hashlib
HERE=Path(__file__).resolve().parent
x,y,u,v,a,b=S.symbols('x y u v a b')
rows=[]
for name,r,f0,f1,F1 in [
 ('contact_at_marked',1,b+a*x+x**3-x**4*y,-b*u-a+v,(b*u+a)*v-v*v/2),
 ('contact_away_marked',2,b+a*x-x*x+x**4*y,-b*u-a-v,(b*u+a)*v+v*v/2)]:
    transition={x:1/u,y:u**3*v+u**r}
    assert S.expand(-u*f0.subs(transition,simultaneous=True)-f1)==0
    F0=S.integrate(f0,y)
    difference=S.expand(F0.subs(transition,simultaneous=True)/u**2-F1)
    period=S.expand(difference*u).coeff(u,0)
    assert period==(b if r==1 else a)
    rows.append(dict(name=name,cocycle_exponent=r,
      transition={str(k):str(e) for k,e in transition.items()},
      local_volume_f0=str(f0),local_volume_f1=str(f1),local_antiderivative_F0=str(F0),
      local_primitive_on_u_chart=str(F1),primitive_difference=str(difference),period=str(period)))
f=b-x*x+x**4*y
eta0=-(b*y-x*x*y+x**4*y*y/2+S.Rational(1,2))
eta1=b*u*v+b+v*v/2
assert S.diff(eta0,y)==-f
assert S.cancel(-eta0.subs({x:1/u,y:u**3*v+u*u},simultaneous=True)/u**2-eta1)==0
payload=dict(as_of='2026-11-14',outcome='EXACT_PERIOD_IDENTITIES_COMPUTED',cases=rows,
    remaining_family=dict(a=0,b='nonzero',global_primitive_x_coefficient=str(eta0),
                         global_primitive_u_coefficient=str(eta1)),
    assurance=dict(all_height_obstruction='written exactness argument',
                   finite_map_on_remaining_family_constructed=False,remaining_family_excluded=False),
    source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
(HERE/'volume_period_certificates.json').write_text(json.dumps(payload,indent=2)+'\n')
for row in rows:print(row['name'],row['primitive_difference'],'period',row['period'])
# On the unique surviving model, identify the entire affine ring inside
# the actual polynomial source ring, not merely the abstract surface.
X,Y,F0,F1,F2,F3,A,B,D=S.symbols('X Y F0 F1 F2 F3 A B D')
al=1+X*X*Y;be=-al*al-Y;tt=1+X*X*al
actual=[be*tt*tt,be*X*tt,al*tt,al*X]
C=1+B+D*D
forward=[-Y-C*C-2*A*D-B*D*D,-A-D-2*B*D-D**3,C,D]
sub={A:X*Y,B:X*X*Y,D:X**3*Y+X}
for value,expected in zip(actual,forward):assert S.expand(value-expected.subs(sub))==0
BB=F2-F3*F3-1
AA=-F1-F3-2*BB*F3-F3**3
YY=-F0-F2*F2-2*AA*F3-BB*F3*F3
reverse=[YY,AA,BB,F3]
source_sub=dict(zip((F0,F1,F2,F3),actual))
for value,expected in zip(reverse,(Y,X*Y,X*X*Y,X**3*Y+X)):
    assert S.expand(value.subs(source_sub,simultaneous=True)-expected)==0
# Check that the three scroll equations generate the whole source kernel.
relations=[A*A-Y*B,A*(B+1)-Y*D,B*(B+1)-A*D]
full=S.groebner([A-X*Y,B-X*X*Y,D-X**3*Y-X],X,D,B,A,Y,order='lex')
kernel=[e.as_expr() for e in full.polys if not e.as_expr().has(X)]
gb=S.groebner(relations,D,B,A,Y,order='lex')
assert list(S.groebner(kernel,D,B,A,Y,order='lex'))==list(gb)
marking=dict(as_of='2026-11-14',outcome='EXACT_MARKED_RING_IDENTITIES_COMPUTED',
  normalized_boundary=dict(B=['t','s**2'],E=['t**2-s**2','s**2*t']),
  old_generators=list(map(str,actual)),old_generators_in_wright=list(map(str,forward)),
  wright_generators_in_old=list(map(str,reverse)),
  kernel_generators=list(map(str,relations)),kernel_groebner_basis=list(map(str,gb)),
  omitted_prime_generators=['Y','A','B+1'],
  wright_overlap=dict(X='1/u',Y='u**3*v-u**2',A='u**2*v-u',B='u*v-1',D='v'),
  omitted_valuations=dict(X=-1,Y=2,A=1,B=0,D=0),
  assurance=dict(ring_identities_preserve_source_inclusion=True,
    source_boundary_reduced=True,finite_map_constructed=False,remaining_wright_problem_excluded=False),
  source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
(HERE/'wright_marking_certificates.json').write_text(json.dumps(marking,indent=2)+'\n')
print('Exact marked ring equality: C[F0,F1,F2,F3]=C[Y,XY,X^2Y,X^3Y+X] inside C[X,Y]')
print('Remaining zero-period model is not excluded; no finite map constructed')
