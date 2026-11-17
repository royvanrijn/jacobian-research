#!/usr/bin/env python3
"""Exact structural controls for the written degree6 cubic-component proof.
No coefficient sampling or independent universal/formal proof claim.
"""
from pathlib import Path
from math import gcd
import hashlib,json,time
import sympy as S
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 start=time.monotonic();ledger={e['id']:e for e in json.loads((ROOT/'MATH_STATUS.json').read_text())['entries']};premises={}
 for key in ['GLDC7','GLDC8','KCUB23','KPAIR1']:
  e=ledger[key];assert e['state']=='proved' and e['artifact_hash']=='sha256:'+sha(ROOT/e['checker']);premises[key]=dict(checker_sha256=sha(ROOT/e['checker']),proof_sha256=sha(ROOT/e['canonical_source']))
 u,v,z,R,x,y=S.symbols('u v z R x y');k,H,lam,A,B,C=S.symbols('k H lambda A B C');r=x+x*x*y;D=lambda f:S.expand(-x*x*S.diff(f,x)+(1+2*x*y)*S.diff(f,y));assert D(r)==0 and D(1/x)==1;assert S.gcd(S.diff(r,x),S.diff(r,y))==1;assert S.gcd(x*x,x-S.Symbol('c'))==1
 # General monic(3,2) coefficient equations, not an imposed candidate.
 a,b,c,g=[S.Function(t)(u) for t in ['a','b','c','g']];P=v**3+a*v*v+b*v+c;Q=v*v+g;br=S.expand(S.diff(P,u)*S.diff(Q,v)-S.diff(P,v)*S.diff(Q,u));poly=S.Poly(br,v);assert poly.nth(3)==2*S.diff(a,u);assert S.expand(poly.nth(2)-(2*S.diff(b,u)-3*S.diff(g,u)))==0;assert S.expand(poly.nth(1)-(2*S.diff(c,u)-2*a*S.diff(g,u)))==0;assert S.expand(poly.nth(0)+b*S.diff(g,u))==0
 b0,a0,c0=S.symbols('b0 a0 c0');gg=S.Rational(2,3)*(b-b0);PP=v**3+a0*v*v+b*v+a0*gg+c0;QQ=v*v+gg;expr=S.expand(S.diff(PP,u)*S.diff(QQ,v)-S.diff(PP,v)*S.diff(QQ,u));assert S.simplify(expr+S.Rational(2,3)*b*S.diff(b,u))==0;assert S.discriminant(-S.Rational(3,2)*k*u*u+H,u)==6*k*H
 # Chain rule for v=x*z=z/u, including differentiation at fixed z.
 PF=S.Function('P')(u,v);QF=S.Function('Q')(u,v);du_fixed_z=lambda f:S.diff(f,u)-v/u*S.diff(f,v);dz=lambda f:S.diff(f,v)/u;chain=S.expand(du_fixed_z(PF)*dz(QF)-dz(PF)*du_fixed_z(QF));assert S.simplify(chain-(S.diff(PF,u)*S.diff(QF,v)-S.diff(PF,v)*S.diff(QF,u))/u)==0
 # The Laurent residue differential identity before specializing a constant bracket.
 PF=S.Function('P')(u,z);QF=S.Function('Q')(u,z);SS=PF**S.Rational(-1,3);br=S.diff(PF,u)*S.diff(QF,z)-S.diff(PF,z)*S.diff(QF,u);identity=S.diff(QF*S.diff(PF,z)*SS,u)-S.diff(QF*S.diff(PF,u)*SS,z);assert S.simplify(identity+br*SS)==0
 # Rational primitive obstruction for every nonconstant affine-linear ell.
 denominator=B*R*u**3-B*u*u+C*u+A;disc=S.discriminant(denominator,u);assert S.Poly(disc,R).nth(2)==-27*A*A*B*B;assert S.discriminant(B*R*u*u-B*u+C,u)==B*B-4*B*C*R;assert S.simplify((u/(A+C*u)*(u+A/C)).cancel().subs(u,-A/C)+A/C**2)==0;assert S.diff(u*u/(2*A),u)==u/A
 # Invariant proportional cubic leading: the z^2 equation is a derivative.
 qq,lambda0,hh=S.symbols('q lambda0 h');bb,ee,ff=[S.Function(t)(u) for t in ['b','e','f']];pp2=3*lambda0*bb+hh;P=lambda0*qq*z**3+pp2*z*z+ee*z+ff;Q=qq*z+bb;br=S.Poly(S.expand(S.diff(P,u)*S.diff(Q,z)-S.diff(P,z)*S.diff(Q,u)),z);assert br.nth(2)==0;assert S.expand(br.nth(1)-(qq*S.diff(ee,u)-2*pp2*S.diff(bb,u)))==0
 # Exhaust every positive z-degree pair through5 with the stated total cap6.
 cases=[]
 for m in range(1,6):
  for n in range(1,m+1):
   if m==n:route='affine-pair' if m==1 else ('common-linear-gcd-and-invariant-target-shear' if m==2 else 'constant-linear-degree-reduction')
   elif n>=3:assert m//gcd(m,n)>6-m and 6-n<=3;route='fibre-or-UFD-degree-contradiction'
   elif n==2:assert m in [3,4,5];route='quadratic-leading-square-class-and-rational-cusp'
   else:route='affine-leading-fibre-division-or-proportional-invariant-obstruction'
   cases.append(dict(m=m,n=n,leading_p_degree_at_most=6-m,leading_q_degree_at_most=6-n,route=route))
 assert len(cases)==15
 parity=[]
 for i in range(5):
  for j in range(3):
   if i+2*j<=4 and i%2==j%2:
    odd=i%2;deg_s=(i-odd)//2+2*((j-odd)//2);assert 3*odd+2*deg_s<=4;assert (odd==0 and deg_s<=2) or(odd==1 and deg_s==0);parity.append(dict(x_exponent=i,one_plus_xy_exponent=j,complete_odd_fibre=bool(odd),square_root_degree=deg_s))
 prior=json.loads((HERE/'validation.json').read_text());assert prior['status']=='pass' and prior['remaining_count']==20;closed=[p for p in prior['remaining_profiles'] if p[0]<=6 and p[2]<=3];assert closed==[[6,3,3],[6,4,3],[6,5,3],[6,6,3]];remaining=[p for p in prior['remaining_profiles'] if p not in closed];assert len(remaining)==16;control=ROOT/'jacobian-followup/equal-negative-completion/source-shear-validation.json';cv=json.loads(control.read_text());assert cv['status']=='pass' and cv['helper_sha256']==sha(ROOT/'jacobian-followup/equal-negative-completion/verify.py') and cv['component_degrees']==[3,7,9] and cv['determinant_one_checked'];assert 6*5==30<=100
 out=dict(status='pass',source_sha256=sha(Path(__file__)),sympy_version=S.__version__,proof_sha256=sha(HERE/'DEGREE6_COMPONENT_EXCLUSION.md'),premises=premises,monic_cubic_quadratic_coefficients_checked=True,chain_rule_v_equals_xz_checked=True,formal_residue_differential_identity_checked=True,all_affine_leading_primitive_strata_checked=True,discriminant_polynomial_cases_checked=True,proportional_invariant_cubic_identity_checked=True,z_degree_cases=cases,special_fibre_square_parity_cases=parity,universal_UFD_fibre_field_proofs='Written, not formalized or independently replayed in full; exact controls are not substituted for the proof.',additional_profiles=closed,total_profiles_excluded=58,remaining_profiles=remaining,remaining_count=16,three_variable_ordinary_degree_interval=[5,7],cubic_component_degree_interval=[7,9],known_cubic_component_control_sha256=sha(control),plane_fibre_bound=30,degree7_extension_proved=False,global_positive_witness_bounds_changed=False,global_simplicity_minima_closed=False,seconds=round(time.monotonic()-start,3));(HERE/'degree6-validation.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print('PASS degree6 structural identities and15 cases; 58excluded16remaining; cubic-component degree[7,9]',flush=True)
if __name__=='__main__':main()
