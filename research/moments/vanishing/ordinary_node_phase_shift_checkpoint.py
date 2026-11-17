#!/usr/bin/env python3
"""Universal exact normalization for the next ordinary-node prefix gate.

This does not construct a twisted quotient basis or prove a finite prefix.
"""
from pathlib import Path
import hashlib,json
import sympy as S


def main():
    p=Path(__file__).resolve().parent
    u,v,t=S.symbols('u v t')
    alpha,m,a,d,b,c,e,f,h,k=S.symbols('alpha m a d b c e f h k')
    A,B,C,ju,jv,j0,lu,lv,q0=S.symbols('A B C ju jv j0 lu lv q0')
    D=m-2*alpha*t
    G=a*u**3+d*v**3+b*u*u+c*v*v+e*u+f*v+k
    q=(A*u*u+B*v*v)/2+C*u*v+lu*u+lv*v+q0
    L=t*(u*v+h)-ju*u-jv*v-j0
    numerator=D*(t*G-q)+L*L/2
    shift={u:u+jv/t,v:v+ju/t}
    shifted_L=S.expand(L.subs(shift,simultaneous=True))
    assert S.expand(shifted_L-(t*u*v+t*h-j0-ju*jv/t))==0
    V=S.expand(t*t*numerator.subs(shift,simultaneous=True))
    delta=t*t*D
    assert S.Poly(V,t).degree()<=4
    pv=S.Poly(V,u,v)
    assert pv.coeff_monomial(u*u*v)==pv.coeff_monomial(u*v*v)==0
    coefficients={name:S.factor(pv.coeff_monomial(term)) for name,term in [('quartic',u*u*v*v),('u_cubic',u**3),('v_cubic',v**3),('u_quadratic',u*u),('v_quadratic',v*v)]}
    expected={'quartic':t**4/2,'u_cubic':a*t**3*D,'v_cubic':d*t**3*D,
              'u_quadratic':t*t*D*(t*b+3*a*jv-A/2),
              'v_quadratic':t*t*D*(t*c+3*d*ju-B/2)}
    assert all(S.expand(coefficients[name]-value)==0 for name,value in expected.items())
    source=S.expand(S.diff(V,t)*delta-V*S.diff(delta,t)-t**4*D*S.diff(D,t)/2)
    assert S.Poly(source,t).degree()<=6
    result={'completed':True,'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'normalization':'P=alpha*T^2+T*(uv+h)+a*u^3+d*v^3+b*u^2+c*v^2+eu+fv+k; Q=m*T^2/2+T*(ju*u+jv*v+j0)+q2',
            'D':str(D),'translation':{'old_u':'u+jv/t','old_v':'v+ju/t'},'Jacobian':'1',
            'normalized_phase':'V/(t^2*D)','determinant_prefactor':'constant*D^(-1/2)',
            'V':str(V),'delta':str(delta),'normalized_connection_source_numerator':str(source),
            'mixed_cubic_coefficients':'identically zero','parameter_degree_bounds':{'V':4,'delta':3,'connection_source_numerator':6},
            'selected_coefficients':{name:str(value) for name,value in coefficients.items()},
            'required_nonzero_principal_coefficients':{
                'irreducible_nodal':'quartic,u_cubic,v_cubic; a*d!=0,D!=0',
                'transverse_line_conic':'quartic,u_cubic,v_quadratic; a!=0,d=0 and(c!=0 or B!=0)',
                'independent_triangle':'quartic,u_quadratic,v_quadratic; a=d=0 and(b!=0 or A!=0),(c!=0 or B!=0)'},
            'next_obligation':'Construct full twisted quotient divisions using the retained derivative/divergence corrections, with bounded t denominators and all coefficient strata. A commutative critical Groebner basis alone does not certify the connection or prefix.',
            'finite_prefix_proved':False,'twisted_basis_constructed':False,'new_GMC_exclusion':False}
    (p/'ordinary_node_phase_shift_checkpoint.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({key:result[key] for key in ['completed','mixed_cubic_coefficients','parameter_degree_bounds','finite_prefix_proved','twisted_basis_constructed','new_GMC_exclusion']}))


if __name__=='__main__':main()
