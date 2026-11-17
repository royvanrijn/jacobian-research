#!/usr/bin/env python3
"""Exact controls for QVL1. Universal infinity statements use written proofs."""
from pathlib import Path
import hashlib
import json
import sympy as sp

HERE = Path(__file__).resolve().parent
a, v, w, u, z, h, r, s = sp.symbols('a v w u z h r s')
half = sp.Rational(1, 2)


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
    return clean({k:c*x for k,x in op.items()})


def mul(left,right,var,rate):
    out = {}
    for i,x in left.items():
        for j,y in right.items():
            dy = y
            for k in range(i+1):
                key = i+j-k
                out[key] = out.get(key,0)+sp.binomial(i,k)*h**k*x*dy
                dy = sp.cancel(rate*sp.diff(dy,var))
    return clean(out)


def comm(left,right,var,rate):
    return add(mul(left,right,var,rate),scale(mul(right,left,var,rate),-1))


def power(op,n,var,rate):
    out = {0:sp.Integer(1)}
    for _ in range(n):
        out = mul(out,op,var,rate)
    return out


def poly_at(poly,op,var,rate):
    out = {}
    for (n,),c in sp.Poly(poly,a).terms():
        out = add(out,scale(power(op,n,var,rate),c))
    return out


def op_expr(op):
    return sp.expand(sum(c*a**k for k,c in op.items()))


def degree(coefficient,var):
    numerator,denominator = sp.fraction(sp.cancel(coefficient))
    return sp.degree(numerator,var)-sp.degree(denominator,var)


def max_degree(op,var):
    return max(degree(c,var) for c in op.values())


def limit_op(op,var):
    return clean({k:sp.limit(c,var,sp.oo) for k,c in op.items()})


def coefficient_at(op,var,exponent):
    out = sp.Integer(0)
    for k,c in op.items():
        for term in sp.Add.make_args(sp.expand(c)):
            coefficient,degree_term = term.as_coeff_exponent(var)
            assert not coefficient.has(var),(term,var)
            if degree_term == exponent:
                out += coefficient*a**k
    return sp.expand(out)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


prior = json.loads((HERE/'infinity_mate_controls.json').read_text())
assert prior['status']=='PASS'
for name,digest in prior['hashes'].items():
    assert sha(HERE/name)==digest,name

# The ordinary limit can have smaller order than the actual operator.
F = a*a+3*a+1
ordinary = poly_at(F,{1:sp.Integer(1)},z,1/(2*z))
ordinary = add(ordinary,{5:z**-2,3:2*z**-3})
ordinary_controls = []
for M in [half,sp.Integer(1),sp.Integer(3)]:
    K = a**3+2
    Q = scale(poly_at(K,{1:sp.Integer(1)},z,1/(2*z)),z**(2*M))
    Q = add(Q,{1:z**(2*M-1)})
    result = comm(ordinary,Q,z,1/(2*z))
    expected = sp.expand(h*M*sp.diff(F,a)*K)
    assert sp.expand(coefficient_at(result,z,2*(M-1))-expected)==0
    assert max_degree(result,z)<=2*(M-1)
    ordinary_controls.append({
        'mate_u_exponent':str(M),
        'commutator_u_exponent':str(M-1),
        'leading_identity':'PASS',
    })
assert max(ordinary)==5 and max(limit_op(ordinary,z))==2

# Euler leading coefficients use finite differences, including M=0.
Fe = a*a+h*a+1
Pe = add(poly_at(Fe,{1:sp.Integer(1)},z,z/2),{4:z**-2})
euler_controls = []
for M in [half,sp.Integer(2),sp.Integer(0),sp.Integer(-1)]:
    K = a+2
    Q = add(scale(poly_at(K,{1:sp.Integer(1)},z,z/2),z**(2*M)),
            {0:z**(2*M-1)})
    result = comm(Pe,Q,z,z/2)
    expected = sp.expand((Fe.subs(a,a+h*M)-Fe)*K)
    assert sp.expand(coefficient_at(result,z,2*M)-expected)==0
    if M==0:
        assert max_degree(result,z)<0
    else:
        assert max_degree(result,z)==2*M
    euler_controls.append({
        'mate_v_exponent':str(M),
        'finite_difference_identity':'PASS',
    })

# Repeated-root transport for an entire declared proper-power family.
G = 4*a-3*r*a*a
power_controls = []
for d,k in [(2,2),(2,3),(2,5),(3,2)]:
    alpha = -1/h
    rate = sp.Rational(1,k)*w**(2-k*d)
    U = sp.Rational(k,k*d-1)*w**(k*d-1)
    assert sp.cancel(rate*sp.diff(U,w)-1)==0
    old_a = {1:1/w,0:alpha}
    old_v = {0:w**k}
    assert comm(old_a,old_v,w,rate)=={0:h*w**(k*(1-d))}
    H = 2*(1+h*a)**k
    S = add(scale(poly_at(H,old_a,w,rate),w**k),
            poly_at(G,old_a,w,rate))
    assert max_degree(S,w)<=0
    expected_limit = {k:2*h**k,0:G.subs(a,alpha)}
    assert add(limit_op(S,w),scale(expected_limit,-1))=={}
    assert sp.expand(H.subs(h,0)-2)==0
    power_controls.append({
        'cyclic_degree':d,'root_multiplicity':k,
        'ordinary_coordinate':str(U),
        'generic_order':max(S),
        'limit':str(op_expr(limit_op(S,w))),
    })
    # Higher differential order with decaying coefficients is also allowed.
    if d==2 and k==3:
        perturbed = add(S,scale(poly_at(a**7,old_a,w,rate),h**2/w**k))
        assert max(perturbed)==7
        assert max_degree(perturbed,w)<=0
        assert add(limit_op(perturbed,w),scale(expected_limit,-1))=={}

# The previous fractional-cancellation control has two repeated roots.
alpha = sp.I
rate = sp.Rational(1,2)*w**-2
old_a = {1:1/w,0:alpha}
repeated = add(scale(poly_at((a*a+1)**2,old_a,w,rate),w*w),
               poly_at(G,old_a,w,rate))
assert max_degree(repeated,w)<=0
expected = {2:-4,0:G.subs(a,sp.I)}
assert add(limit_op(repeated,w),scale(expected,-1))=={}

# A Newton certificate can include positive intermediate v powers.
rho_one_a = {1:1/v,0:sp.Integer(1)}
newton = add(
    scale(poly_at((a-1)**2,rho_one_a,v,1),v**2),
    scale(poly_at((a-1)**3,rho_one_a,v,1),v),
    poly_at(1+a*a,rho_one_a,v,1),
    scale(poly_at((a-1)**5,rho_one_a,v,1),v**-2))
assert max_degree(newton,v)<=0
assert limit_op(newton,v)=={2:sp.Integer(1),0:sp.Integer(2)}

# Euler root normalization applies even when the limit has degree one.
alpha = -2/h
euler_a = {1:v**-2,0:alpha}
euler_root = add(
    scale(poly_at(2+h*a,euler_a,v,v),v*v),
    poly_at(G,euler_a,v,v))
assert max_degree(euler_root,v)<=0
assert add(limit_op(euler_root,v),scale({1:h,0:G.subs(a,alpha)},-1))=={}
assert comm(euler_a,{0:v},v,v)=={0:h/v}

# At the Euler boundary internal ordering survives in the limit.
euler_a_one = {1:v**-2,0:sp.Integer(1)}
euler_square = add(
    scale(poly_at((a-1)**2,euler_a_one,v,v),v**4),
    scale(poly_at(a**3,euler_a_one,v,v),1/v),
    poly_at(G,euler_a_one,v,v))
assert max_degree(euler_square,v)<=0
expected_euler = {2:sp.Integer(1),1:-2*h,0:G.subs(a,1)}
assert add(limit_op(euler_square,v),scale(expected_euler,-1))=={}

# The proper-power cubic still passes the earlier minimum rank calculation.
H3 = sp.Poly(2*(1+h*a)**3,a)
leading = H3.coeff_monomial(a**3)
coefficient2 = v*H3.coeff_monomial(a*a)-3*r
coefficient1 = v*H3.coeff_monomial(a)+4
coefficient0 = v*H3.coeff_monomial(1)
C = sp.Matrix([
    [0,0,(s-coefficient0)/(v*leading)],
    [1,0,-coefficient1/(v*leading)],
    [0,1,-coefficient2/(v*leading)],
])
A = (v/h*C).applyfunc(sp.cancel)
assert all(sp.degree(x,v)<=1 for x in A if x!=0)

# A surviving leading-linear cubic has an additional infinity contribution.
linear_cubic_symbol = v*(2+h*a)+4*a-3*r*a*a+h*h*a**3-s
growth = sp.expand(linear_cubic_symbol.subs({v:w*w,a:w*a})/w**3)
growth_limit = sp.limit(growth,w,sp.oo)
assert sp.expand(growth_limit-(h*a+h*h*a**3))==0
assert sp.degree(growth_limit/a,a)==2
assert sp.expand((2+h*a).subs(a,-2/h))==0
linear_cubic_A = sp.Matrix([
    [0,0,v*(s-2*v)/h**3],
    [v/h,0,-v*(4+h*v)/h**3],
    [0,v/h,3*r*v/h**3],
])
assert all(sp.degree(x,v)<=2 for x in linear_cubic_A if x!=0)

# Sharp positive finite pairs.
lam,mu,c = sp.symbols('lam mu c',nonzero=True)
positive_linear = comm({1:lam*v,0:mu*v},{0:v/lam},v,1/v)
assert positive_linear=={0:h}
positive_constant = comm({4:h,2:-3*r,1:4,0:c*v*v},
                         {1:-1/(2*c)},v,1/v)
assert positive_constant=={0:h}

# A different proper-power face has a genuine Weyl mate.
E = {1:v**3}
U = {0:-1/v}
assert comm(E,U,v,1/v)=={0:h}
P = add(power(E,2,v,1/v),U)
assert P=={2:v**6,1:3*h*v**4,0:-1/v}
assert comm(P,scale(E,-1),v,1/v)=={0:h}

# A rational commuting Darboux algebra has nonlinear constant limits.
D = {1:sp.Integer(1)}
A_darboux = {1:sp.Integer(1),0:-h/u}
L2 = {2:sp.Integer(1),0:-2*h*h/u**2}
L3 = {3:sp.Integer(1),1:-3*h*h/u**2,0:3*h**3/u**3}
assert mul(L2,A_darboux,u,1)==mul(A_darboux,power(D,2,u,1),u,1)
assert mul(L3,A_darboux,u,1)==mul(A_darboux,power(D,3,u,1),u,1)
assert comm(L2,L3,u,1)=={}
assert add(power(L3,2,u,1),scale(power(L2,3,u,1),-1))=={}
darboux_samples = []
for label,operator in [
    ('L2',L2),('L3',L3),
    ('L2+L3',add(L2,L3)),
    ('L2**2+L2*L3',add(power(L2,2,u,1),mul(L2,L3,u,1)))]:
    assert max_degree(operator,u)<=0
    limiting = limit_op(operator,u)
    assert max(limiting)>=2
    darboux_samples.append({'operator':label,'limit':str(op_expr(limiting))})

# A finite P with a nonterminating formal mate, consistent with the obstruction.
formal_P = {1:sp.Integer(1),2:h}
formal_controls = []
for N in [1,4]:
    Qjet = {j:u*(-2*h)**j for j in range(N+1)}
    defect = add(comm(formal_P,Qjet,u,1),{0:-h})
    assert defect=={N+1:-h*(-2*h)**(N+1)}
    formal_controls.append({'jet_index':N,'defect_h_order':N+2,
                            'mate_momentum_order':N})
# This identity is the exact geometric-series telescoping mechanism.
x = sp.symbols('x')
assert sp.expand((1+2*h*x)*sum((-2*h*x)**j for j in range(5))
                 -(1-(-2*h*x)**5))==0

names = [
    'ALGEBRAIC_INFINITY_LIMITS.md',
    'check_infinity_limits.py',
    'CYCLIC_INFINITY_MATE_OBSTRUCTION.md',
    'infinity_mate_controls.json',
    'FOURIER_CYCLIC_RANK_BUDGET.md',
    'fourier_rank_budget_controls.json',
]
out = {
    'claim':'QVL1','status':'PASS',
    'ordinary_limit_controls':ordinary_controls,
    'generic_order_exceeds_limit_order':{'generic':5,'limit':2},
    'Euler_finite_difference_controls':euler_controls,
    'proper_power_transport_controls':power_controls,
    'decaying_order7_correction_same_limit':'PASS',
    'two_repeated_root_polynomial_limit':str(op_expr(limit_op(repeated,w))),
    'ordinary_Newton_intermediate_positive_power':'PASS',
    'Euler_linear_limit':'PASS',
    'Euler_ordering_correction':'E*(E-2*h)',
    'proper_power_cubic_companion_degree':1,
    'proper_power_cubic_rank_uses_written_QKR1_framework':{'pushforward':6,'inverse_Fourier':3},
    'leading_linear_cubic_Newton_control':{
        'growing_branch_polynomial':str(sp.factor(growth_limit)),
        'growing_momentum_branches':2,
        'ordinary_companion_degree':2,
        'pushforward_rank':6,
        'inverse_Fourier_rank':4,
        'rank_uses_written_Newton_and_QKR1_proof':True,
    },
    'positive_finite_controls':['v_linear_H','v_square_constant_H','different_proper_power_face'],
    'commuting_Darboux_controls':darboux_samples,
    'formal_mate_controls':formal_controls,
    'hashes':{name:sha(HERE/name) for name in names},
    'research_context':[
        {'source':'https://arxiv.org/html/math-ph/0206029',
         'used':'Section2 bounded-coefficient bispectral context only; no classification theorem imported.'},
        {'source':'https://arxiv.org/html/2210.00257',
         'used':'Theorem2.8 leading-bracket framework; context only, no classification theorem imported.'},
        {'source':'https://arxiv.org/html/1111.6100',
         'used':'Remark1.13 and Section2 weighted/fractional Weyl leading calculus; the infinity-limit and whole-curve proof is written directly.'},
    ],
    'coverage_boundary':(
        'Exact finite scalar Ore identities, full repeated-root transforms, '
        'ordinary/Euler limits, finite Newton examples, a proper-power cubic '
        'companion, positive finite pairs and finite formal-mate controls. '
        'Universal no-mate, all degrees and algebraic covers use written '
        'Puiseux/transport proofs and QVM1. The complete v-linear class '
        'requires bounded lower coefficients, unlike the primitive QVM1 '
        'extension to arbitrary fractional growth below one. Other leading '
        'faces, H of degree<=1, unrestricted deformed first images and '
        'non-scalar carriers remain outside the exclusion.'),
    'independent_replay':False,'formal_verification':False,'external_review':False,
}
(HERE/'infinity_limits_controls.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({
    'claim':out['claim'],'status':out['status'],
    'proper_power_transports':len(power_controls),
    'ordinary_and_Euler_limits':'PASS',
    'positive_finite_controls':len(out['positive_finite_controls']),
    'formal_mate_boundary':'PASS',
}))
