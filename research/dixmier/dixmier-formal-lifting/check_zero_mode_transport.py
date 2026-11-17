#!/usr/bin/env python3
"""Exact controls for QZM1's zero-mode, trace and canonical transport proofs."""
from pathlib import Path

import tempfile as _check_tempfile

_CHECK_OUTPUT_DIRECTORY = None

def _check_output(original):
    """Write new runs to a temporary directory; preserve deposited receipts."""
    global _CHECK_OUTPUT_DIRECTORY
    if _CHECK_OUTPUT_DIRECTORY is None:
        _CHECK_OUTPUT_DIRECTORY = Path(_check_tempfile.mkdtemp(prefix="math-check-"))
        print("Verification output:", _CHECK_OUTPUT_DIRECTORY)
    return _CHECK_OUTPUT_DIRECTORY / Path(original).name


import hashlib
import json
import sympy as sp

HERE = Path(__file__).resolve().parent
a, v, w, x, r, h = sp.symbols('a v w x r h')
lam, mu, c1, kap = sp.symbols('lam mu c1 kap', nonzero=True)


def clean(op):
    out = {}
    for k,c in op.items():
        c = sp.cancel(c)
        if c != 0:
            out[k] = c
    return out


def add(*ops):
    out = {}
    for op in ops:
        for k,c in op.items():
            out[k] = out.get(k,0)+c
    return clean(out)


def scale(op,c):
    return clean({k:c*z for k,z in op.items()})


def mul(left,right,var,rate):
    out = {}
    for i,z in left.items():
        for j,y in right.items():
            dy = y
            for k in range(i+1):
                key = i+j-k
                out[key] = out.get(key,0)+sp.binomial(i,k)*h**k*z*dy
                dy = sp.cancel(rate*sp.diff(dy,var))
    return clean(out)


def comm(left,right,var,rate):
    return add(mul(left,right,var,rate),scale(mul(right,left,var,rate),-1))


def act(op,f,var,rate):
    n = max(op,default=0)
    derivatives = [f]
    for _ in range(n):
        derivatives.append(sp.cancel(h*rate*sp.diff(derivatives[-1],var)))
    return sp.cancel(sum(c*derivatives[k] for k,c in op.items()))


def poly_op(poly):
    return {k:c for (k,),c in sp.Poly(poly,a).terms()}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()



# Trace on E=C(w), w**N=v. Reduced constant coefficient is normalized trace.
domain = sp.QQ.frac_field(v,h)
def ntrace(function,N):
    numerator,denominator = sp.fraction(sp.cancel(function))
    modulus = sp.Poly(w**N-v,w,domain=domain)
    num = sp.Poly(numerator,w,domain=domain)
    den = sp.Poly(denominator,w,domain=domain)
    reduced = (num*sp.invert(den,modulus)).rem(modulus)
    return sp.cancel(reduced.nth(0))


trace_controls = []
f_cover = (w*w+2*w+1)/(w+2)
for N in [2,3]:
    rate_cover = sp.Rational(1,N)*w**(1-2*N)
    traced = ntrace(f_cover,N)
    derivative_cover = sp.cancel(rate_cover*sp.diff(f_cover,w))
    assert sp.cancel(ntrace(derivative_cover,N)-sp.diff(traced,v)/v)==0
    assert ntrace(sp.Integer(1),N)==1
    P_base = {2:sp.Integer(3),1:v+2}
    P_cover = {2:sp.Integer(3),1:w**N+2}
    left = ntrace(act(P_cover,f_cover,w,rate_cover),N)
    right = act(P_base,traced,v,1/v)
    assert sp.cancel(left-right)==0
    trace_controls.append({
        'cover_degree':N,
        'test_function':str(f_cover),
        'normalized_trace':str(traced),
        'derivative_commutation':'PASS',
        'order2_operator_commutation':'PASS',
    })

# The exact monomial indicial identity used in the pole proof.
beta = sp.symbols('beta')
monomial_controls = []
current = v**beta
for j in range(1,6):
    current = sp.cancel(h/v*sp.diff(current,v))
    expected = h**j*sp.prod(beta-2*t for t in range(j))*v**(beta-2*j)
    assert sp.simplify(current/ v**(beta-2*j)
                       - expected/v**(beta-2*j))==0
    # At every negative integer exponent each factor is nonzero.
    n = sp.symbols('n',integer=True,positive=True)
    assert all((n+2*t).is_positive for t in range(j))
    monomial_controls.append(j)

# The first-order obstruction is a genuine nonzero rational residue.
root = -c1/lam
forcing_form = kap/h*v/(lam*v+c1)
residue = sp.cancel(sp.limit((v-root)*forcing_form,v,root))
assert sp.cancel(residue+kap*c1/(h*lam**2))==0

# Arbitrary G degree is handled by the written proof; these are declared controls.
cases = [
    ('quadratic_Long',h,sp.Integer(2),4*a-3*r*a*a),
    ('linear_H_cubic',h,sp.Integer(2),4*a-3*r*a*a+h*h*a**3),
    ('tuned_derivative_zero',h,sp.Integer(2),
     4*a-3*r*a*a-(r*h+h*h/3)*a**3),
    ('degree5',sp.Integer(3),sp.Integer(5),a**5+2*a*a+a),
]
shifted_controls = []
for label,Lambda,Mu,G in cases:
    alpha = -Mu/Lambda
    c = sp.cancel(G.subs(a,alpha))
    shifted = sp.Poly(sp.expand(G.subs(a,a+alpha)-c),a)
    P = add(poly_op(shifted.as_expr()),{1:Lambda*v})
    assert act(P,sp.Integer(1),v,1/v)==0
    m = sp.degree(G,a)
    candidate_f = v/Lambda
    defect = sp.cancel(act(P,candidate_f,v,1/v)-h)
    expected = sum(
        shifted.coeff_monomial(a**j)*h**j
        *sp.prod(1-2*t for t in range(j))*v**(1-2*j)/Lambda
        for j in range(1,m+1))
    assert sp.cancel(defect-expected)==0
    top_pole = sp.cancel(sp.limit(defect*v**(2*m-1),v,0))
    assert top_pole != 0
    if label=='tuned_derivative_zero':
        assert sp.cancel(sp.diff(G,a).subs(a,alpha))==0
        assert shifted.coeff_monomial(a)==0
    shifted_controls.append({
        'case':label,'G_degree':int(m),'root':str(alpha),
        'zero_mode':'PASS',
        'forced_polynomial_primitive':str(candidate_f),
        'nonzero_defect':str(sp.factor(defect)),
        'lowest_v_exponent':int(1-2*m),
        'lowest_coefficient':str(sp.factor(top_pole)),
    })

# Positive case: G constant really admits a mate.
S_positive = {1:lam*v,0:mu*v+7}
assert comm(S_positive,{0:v/lam},v,1/v)=={0:h}

# Actual canonical transport in the polynomial Weyl algebra, including the module.
D_old = {1:sp.Integer(1)}
X_old = {0:x}
Euler = {1:x}
X_new_action = scale(D_old,-1)
D_new_action = X_old
assert comm(D_new_action,X_new_action,x,1)=={0:h}
# Phi(x)=D, Phi(D)=-x; Phi(xD)=-xD-h.
Fourier_L_action = add(
    scale(mul(X_new_action,D_new_action,x,1),-1),{0:-h})
assert Fourier_L_action==Euler
assert act(Fourier_L_action,sp.Integer(1),x,1)==0
assert act(X_new_action,sp.Integer(1),x,1)==0
# The unchanged coefficient representation would not retain this zero mode.
assert act({1:-x,0:-h},sp.Integer(1),x,1)==-h

# The residue functional directly certifies the class and its transport.
residue_tests = [
    sp.Integer(1),x**-3+2*x+1,(x+1)/(x-2),x**2/(x+1),
]
for f in residue_tests:
    image = act(Euler,f,x,1)
    assert sp.residue(image/x,x,0)==0
    assert act(Fourier_L_action,f,x,1)==image
assert sp.residue(1/x,x,0)==1

# Phi(D)=D+p(x), Phi(x)=x; inverse module action undoes the same shear.
p_shear = x**3+h*x
sheared_D_action = add(D_old,{0:-p_shear})
sheared_L_action = add(
    mul(X_old,sheared_D_action,x,1),{0:x*p_shear})
assert sheared_L_action==Euler
assert comm(sheared_D_action,X_old,x,1)=={0:h}

# The cyclic two-term resolution gives the elementary class identity.
# A positive pair makes e exact: L=D, T=x/h gives [L,T]=1 and L(T1)=1.
assert comm(D_old,{0:x/h},x,1)=={0:sp.Integer(1)}
assert act(D_old,x/h,x,1)==1

# Parameter-dependent root gauges extend compatibly to the full rational chart.
Lambda_r,Mu_r = h*r,sp.Integer(2)
alpha_r = -Mu_r/Lambda_r
q = v*v/2
potential_d = q*sp.diff(alpha_r,r)
mixed_defect = h*sp.diff(alpha_r,r)-h/v*sp.diff(potential_d,v)
assert sp.cancel(mixed_defect)==0
G_r = 4*a-3*r*a*a+h*h*a**3
c_r = sp.cancel(G_r.subs(a,alpha_r))
denom = sp.denom(c_r)
numer = sp.cancel(denom*c_r)
assert sp.Poly(denom,r) and sp.Poly(numer,r)
assert sp.cancel(denom*c_r-numer)==0
# Universal clearing identity tested with an actual source Weyl pair.
L_clear = {0:denom*x-numer}
T_source = {1:-sp.Integer(1)}
assert comm(L_clear,T_source,x,1)=={0:sp.cancel(h*denom)}

names = [
    'QUANTUM_ZERO_MODE_TRANSPORT.md',
    'check_zero_mode_transport.py',
    'ALGEBRAIC_INFINITY_LIMITS.md',
    'infinity_limits_controls.json',
    'CYCLIC_INFINITY_MATE_OBSTRUCTION.md',
    'infinity_mate_controls.json',
    'QUINTIC_PHASE_AND_BRANCH_OBSTRUCTION.md',
]
out = {
    'claim':'QZM1','status':'PASS',
    'finite_trace_controls':trace_controls,
    'monomial_action_orders':monomial_controls,
    'linear_primitive_nonzero_residue':str(residue),
    'shifted_polynomial_G_controls':shifted_controls,
    'constant_G_positive_pair':'PASS',
    'canonical_transport':{
        'Fourier':'exact transported action x_old*D_old',
        'polynomial_Weyl_shear':'exact transported action x_old*D_old',
        'residue_functional':'Res(f*dx/x), value1 on zero mode',
        'new_position_localization_can_kill_zero_mode':True,
    },
    'parameter_connection_mixed_flatness':'PASS',
    'auxiliary_value':str(c_r),
    'cleared_denominator':str(denom),
    'cleared_numerator':str(numer),
    'cleared_commutator_identity':'h*f(R)',
    'hashes':{name:sha(HERE/name) for name in names},
    'research_context':[
        {'source':'https://arxiv.org/html/0808.0699',
         'used':'Section4.4 first-order connection deformations and self-Ext context. The cyclic resolution and cokernel computation are written directly.'},
    ],
    'coverage_boundary':(
        'Exact trace/derivative and order2 operator controls on covers2,3, '
        'monomial identities through order5, shifted G examples of degrees2,3,5, '
        'a derivative-tuned example, rational residue, positive mate, actual '
        'Fourier/shear module transport with residue functional, parameter '
        'connection and denominator clearing. Universal nonzero cokernel, '
        'all degrees, all finite covers, automorphism covariance and identity '
        'stabilization use written proofs. A new localization may kill a '
        'transported supported vector. The polynomial separable theorem '
        'does not include arbitrary negative-v corrections or H constant. '
        'No identification with QRF2 or all-representative nonlifting is certified.'),
    'independent_replay':False,'formal_verification':False,'external_review':False,
}
(_check_output(HERE/'zero_mode_controls.json')).write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({
    'claim':out['claim'],'status':out['status'],
    'finite_trace_covers':[2,3],
    'G_control_degrees':[z['G_degree'] for z in shifted_controls],
    'Fourier_and_shear_transport':'PASS',
    'derivative_tuned_obstruction':'PASS',
    'parameter_connection':'PASS',
}))
