#!/usr/bin/env python3
"""Derive the exact PF160 modular branch certificate; no numerical roots."""
from pathlib import Path
import json
import sympy as S

BASE=Path(__file__).resolve().parent
v,w=S.symbols('v w'); u=1+v
prime,p0,tau0=101,3,1
F=(-u**4*v**30*w**5+(3*u**4*v**20+tau0*v**41)*w**4
   -3*u**4*v**10*w**3+(u**4*(1-v**12*(1-2*v))-2*tau0*u*v**41)*w**2
   +u**4*(v**2*(1-2*v)+v**11)*w+u**4*(v-p0)+tau0*u**2*v**41)
D=S.Poly(S.discriminant(F,w,domain=S.GF(prime).poly_ring(v)),v,modulus=prime)
base=S.Poly(v**96*(v+1)**10,v,modulus=prime)
R,rem=D.div(base)
assert rem.is_zero and D.degree()==321 and R.degree()==215
U,V,G=S.gcdex(R,R.diff())
assert G.is_one and (U*R+V*R.diff()).is_one
exception_value=int(R.resultant(R.diff()))%prime
assert exception_value!=0
coeff=lambda poly:[int(poly.nth(i))%prime for i in range(poly.degree()+1)]
cert={'claim':'PF2RULED160','prime':prime,'p0':p0,'tau0':tau0,
 'discriminant_degree':D.degree(),'base_factor_powers':{'v':96,'v_plus_1':10},
 'reduced_degree':R.degree(),'reduced_coefficients_ascending_mod_prime':coeff(R),
 'bezout_U_ascending_mod_prime':coeff(U),'bezout_V_ascending_mod_prime':coeff(V),
 'reduced_at_zero':int(R.eval(0))%prime,'reduced_at_minus_one':int(R.eval(-1))%prime,
 'reduced_leading_coefficient':int(R.LC())%prime,
 'exception_polynomial_at_one_mod_prime':exception_value,
 'simple_generic_branch_points_at_least':215,
 'geometric_genus_at_least':104,
 'primitive_degree_from_three_double_poles':3,
 'birational_bidegree_genus_at_most_if_mate':8,
 'universal_exception_polynomial_definition':'B(tau)=Res_v(R(v,3,tau),partial_v R(v,3,tau)); R=Disc_w(F)/(v^96*(v+1)^10)',
 'universal_exception_polynomial_expanded':False,
 'universal_discriminant_expanded_over_Z_p_tau':False,
 'independent_replay':False,'formal_verification':False,'external_review':False}
(BASE/'affine_discriminant_certificate.json').write_text(json.dumps(cert,indent=2)+'\n')
print(json.dumps({k:cert[k] for k in ['claim','prime','reduced_degree','reduced_at_zero','reduced_at_minus_one','reduced_leading_coefficient','exception_polynomial_at_one_mod_prime','geometric_genus_at_least']},indent=2))
