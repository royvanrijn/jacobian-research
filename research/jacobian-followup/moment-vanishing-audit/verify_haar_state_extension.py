"""Exact source/state/extension controls, not universal geometric verification."""
from pathlib import Path
from hashlib import sha256
import json
import sympy as s

HERE=Path(__file__).resolve().parent
t,u,x,c=s.symbols('t u x c')
P=x**3+x**2
a=s.Function('a')(x,t)
Lx=lambda q:s.diff(q,x)+t*s.diff(P,x)*q
Lt=lambda q:s.diff(q,t)+P*q
assert s.expand(Lt(Lx(a))-Lx(Lt(a)))==0
S=lambda q:P*q-s.diff(q,u)
U=lambda q:u*q
q=x**2*u**3+u+1
assert s.expand(U(S(q))-S(U(q))-q)==0

# Taylor projector on a point-supported Weyl module: T=-d/du.
C=s.symbols('c0:6')
v=sum(C[j]*u**j for j in range(6))
T=lambda q:-s.diff(q,u)


def iterate(fun,q,j):
    for _ in range(j):q=fun(q)
    return q


def projector(q):
    return s.expand(sum(u**j*iterate(T,q,j)/s.factorial(j)
                        for j in range(s.degree(q,u)+1))) if q!=0 else s.Integer(0)


assert projector(v)==C[0]
restored=sum((-1)**j*u**j*projector(iterate(T,v,j))/s.factorial(j) for j in range(6))
assert s.expand(restored-v)==0
assert T(projector(v))==0
for j in range(1,6):assert iterate(T,u**j,j)==(-1)**j*s.factorial(j)

# Formula(4) for the rank-one exponential connection, at nontrivial vectors.
vec=1+t+2*t**3
for m in range(7):
    D=lambda q:s.diff(q,t)+c*q
    expected=iterate(D,vec,m).subs(t,0)/s.factorial(m)
    actual=s.expand((vec*s.exp(c*t)).series(t,0,7).removeO()).coeff(t,m)
    assert s.expand(actual-expected)==0

x,y,z,w=s.symbols('x y z w')
fplus=x**8+14*x**4*y**4+y**8
fminus=z**8+14*z**4*w**4+w**8
E=lambda q:x*s.diff(q,z)+y*s.diff(q,w)
F=lambda q:z*s.diff(q,x)+w*s.diff(q,y)


def haar(poly):
    ans=0
    for (aa,bb,cc,dd),coeff in s.Poly(s.expand(poly),x,y,z,w).terms():
        if cc==bb and dd==aa:
            ans+=coeff*(-1)**bb*s.factorial(aa)*s.factorial(bb)/s.factorial(aa+bb+1)
    return s.expand(ans)


Pactual=fplus+fminus;amplitude=F(fplus)
stokes=[]
for m in (1,3):
    first=haar(E(amplitude)*Pactual**m)
    second=m*haar(amplitude*E(Pactual)*Pactual**(m-1))
    assert first!=0 and s.expand(first+second)==0
    stokes.append({'m':m,'derivative_amplitude_term':str(first),
                   'phase_derivative_term':str(second)})

# A nonsplit ENTIRE lattice whose localized connection is semisimple.
f=(s.exp(t)-1)/t
assert s.simplify(t*s.diff(f,t)-(t-1)*f-1)==0
jet=s.series(f,t,0,9).removeO()
assert all(jet.coeff(t,j)==1/s.factorial(j+1) for j in range(9))
recurrences=[]
for j in range(6):
    if j==0:
        relation=t*s.diff(f,t)+(1-t)*f-1
    else:
        relation=t*s.diff(f,t,j+1)+(j+1-t)*s.diff(f,t,j)-j*s.diff(f,t,j-1)
    assert s.simplify(relation)==0
    recurrences.append({'derivative_order':j,'formal_constant':str(s.Rational(1,j+1))})
for j in range(1,6):
    assert s.diff(t**-j,t)==-j*t**(-j-1)
assert s.simplify(1+t*f-s.exp(t))==0
assert s.simplify(-1/t+s.exp(t)/t-f)==0

# Nonzero state ideal is not normalization over a ring.
ee,bb=s.symbols('ee bb')
ring=s.groebner([ee*bb],ee,bb,domain=s.QQ)
for j in (1,2,5):
    assert ring.reduce(bb*ee**j)[1]==0
    assert ring.reduce(ee**j)[1]!=0
assert s.groebner([bb],ee,bb,domain=s.QQ).reduce(s.Integer(1))[1]==1
fe=sum(ee**m*t**m/s.factorial(m+1) for m in range(8))
res=s.expand(t*s.diff(fe,t)-(ee*t-1)*fe-1)
assert all(res.coeff(t,m)==0 for m in range(8))

result={
 'completed':True,'sympy_version':s.__version__,
 'checker_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
 'twisted_source_connection_commutator_exact':True,
 'ordinary_graph_Weyl_commutator_exact':True,
 'point_supported_Taylor_projector_and_reconstruction_exact':True,
 'formal_state_Taylor_formula_controls_exact':True,
 'actual_exceptional_W8_nonzero_Stokes_cancellations':stokes,
 'entire_nonsplit_lattice_relation_and_localized_decomposition_exact':True,
 'entire_lattice_jet_and_fibre_recurrence_controls':recurrences,
 'nonzero_nonunit_state_ideal_ring_boundary_exact':True,
 'intrinsic_formal_state_scope':'ANY commutativeC-algebra;H/tH=R;formal statesR*lambda;all-original-marker safetyiff normalized polynomial state,not merely nonzero ideal',
 'actual_failure_scope':'indecomposable actual period module;nonzero C[t]-torsionfree quotient byC[t],t-bijective,and nonsplit extension even asC[t]-modules',
 'semisimple_sufficiency_scope':'UNLOCALIZED ordinary topGauss-Manin semisimple implies pure phase all-original-marker safety',
 'generic_fibre_criterion':'H^(n-1)generic fibre0 and H^n(X)=C imply topGauss-Manin delta_c;eventual purity forcesc0 and all-original-marker safety',
 'localized_semisimplicity_silently_promoted_to_unlocalized':False,
 'generic_rank_of_full_source_identified_with_actual_period_rank':False,
 'generic_fibre_acyclicity_verified_for_all_remaining_exceptional_phases':False,
 'unlocalized_semisimplicity_verified_for_all_remaining_exceptional_phases':False,
 'all_pure_exceptional_unit_extensions_proved_split':False,
 'ordinary_coefficient_rings_reduced_or_nonzero_Hom_treated_as_normalized':False,
 'universal_topology_direct_image_state_and_extension_proofs_machine_verified':False,
 'new_degree_support_density_or_prefix_campaign':False,
 'external_review_next_task':False,
 'full_exceptional_MZ_decisions':'OPEN on BOTH models;exclude or realize the actual nonsplit unit extension for all remaining pure phases',
 'independent_replay':False,'formal_verification':False,'external_review':False}
(HERE/'haar_state_extension_verification.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
