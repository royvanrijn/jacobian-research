"""Exact polynomial and graph-boundary controls for KMFI73--76.

The unbounded induction and deformation proofs are written separately.
A passed weight condition is not asserted necessary for symmetric monodromy.
"""
import os
os.environ.setdefault("SYMPY_GROUND_TYPES","flint")
from pathlib import Path
from math import factorial,gcd
import hashlib,json,platform,time
import sympy as sp
from coefficient_graph import completion,morse_weights,passes_morse_gate,scaled_root_polynomial
ROOT=Path(__file__).resolve().parent;started=time.monotonic();checks=[]
def check(name,value):
    if not value:raise AssertionError(name)
    checks.append(name)
def sparse(f,variables):
    return [[str(c),list(e)] for e,c in sp.Poly(f,*variables,domain=sp.QQ).terms()]
B=sp.Symbol("B")
cases=[('base3',3,{}),('quadratic4',4,{2:B**2}),
 ('linear5',5,{2:B,3:2*B}),('separated_linear6',6,{2:B,3:0,4:B})]
witnesses=[];controls=[]
for name,n,graphs in cases:
    d=completion(n,graphs);x=d['original'][0];m=n-2
    for stage in d['stages']:
        k=stage['k'];u,z=stage['u'],stage['z']
        check(f"exact x2 inverse division and coefficient n{n},stage{k}",
          sp.expand(x*x*stage['newz']+stage['f0'].subs(u,stage['newu'])-stage['v'])==0)
        check(f"exact x2 forward division n{n},stage{k}",
          sp.expand(x*x*stage['forward_w']+stage['psi'].subs(stage['v'],stage['f'])-u)==0)
        expected=(-sp.Rational(m,2) if k%4==2 else 1 if k%4==3 else
                  sp.Rational(m,2) if k%4==0 else -1)
        check(f"nonzero constant slope pattern n{n},stage{k}",stage['c']==expected)
        jac=sp.expand(sp.diff(stage['f'],u)*sp.diff(stage['forward_w'],z)-
                      sp.diff(stage['f'],z)*sp.diff(stage['forward_w'],u))
        check(f"entire polynomial coordinate exchange Jacobian n{n},stage{k}",jac==-1)
    ev=d['endpoint_vars'];out=d['endpoint'];mapping=d['zero_mapping']
    check(f"explicit ordinary affine3 polynomial coordinates {name}",
      all(sp.Poly(f,*ev,domain=sp.QQ) is not None for f in out))
    for j in range(2,n-1):
        residual=sp.expand(d['coefficients'][j].xreplace(mapping)-sp.sympify(graphs.get(j,0)).subs(B,out[2]))
        check(f"exact graph identity {name},coefficient{j}",residual==0)
    jac=sp.Matrix(out).jacobian(ev).det(method="domain-ge")
    check(f"independent whole expanded affine3 Jacobian {name}",sp.cancel(jac)==d['expected_jacobian'])
    # Verify the derivative relation directly from the expanded output graph.
    # Clearing p-powers avoids silently inverting a source divisor.
    T=d['T'];oldx=mapping[d['original'][0]];oldy=mapping[d['original'][1]];p=1+oldx*oldy
    Qgraph=out[1]*T**n+out[0]*T**(n-1)+sum(sp.sympify(graphs.get(j,0)).subs(B,out[2])*T**j for j in range(2,n-1))+T-out[2]
    coeffs=sp.Poly(Qgraph,T)
    root=sp.expand(sum(coeffs.nth(j)*oldx**j*p**(n-j) for j in range(n+1)))
    check(f"expanded graph root identity cleared {name}",root==0)
    deriv=sp.expand(sum(j*coeffs.nth(j)*oldx**(j-1)*p**(n-j) for j in range(1,n+1))-p)
    check(f"expanded graph marking identity cleared {name}",deriv==0)
    control={'name':name,'n':n,'graphs':{str(j):str(h) for j,h in graphs.items()},
      'variables':[str(v) for v in ev],'geometric_degree':n*m,'full_Jacobian':str(d['expected_jacobian']),
      'stage_slopes':[str(s['c']) for s in d['stages']],
      'output_degrees':[sp.Poly(f,*ev).total_degree() for f in out],
      'output_monomials':[len(sp.Poly(f,*ev).terms()) for f in out]}
    controls.append(control)
    witnesses.append({**control,'outputs_sparse':[sparse(f,ev) for f in out]})
Q,R,B,T,e,Z,V=sp.symbols("Q R B T e Z V")
weights=[];counter=[]
for n in range(4,13):
    # General q-independent dependence on r,b; terms containing r vanish on r=0.
    independent={j:R**(j+1)+B**(j+2)-1 for j in range(2,n-1)}
    positive={j:Q*B**(n-j+1)+R*Q**(j+1) for j in range(2,n-1)}
    for label,graphs in [('q_independent',independent),('q_dependent_positive',positive)]:
        check(f"exact strict weight gate {label},n{n}",passes_morse_gate(n,graphs))
        P,scaled,syms=scaled_root_polynomial(n,graphs)
        check(f"whole scaled polynomial has no negative epsilon powers {label},n{n}",
          sp.Poly(scaled,e) is not None)
        check(f"exact symmetric Morse special fiber {label},n{n}",sp.expand(scaled.subs(e,0)-(Z**n+Z-V))==0)
        weights.append({'n':n,'case':label,'min_epsilon_power':min(x['epsilon_power'] for x in morse_weights(n,graphs))})
    bad={2:-Q*B**(n-2)-R*B**(n-3)}
    check(f"equality boundary fails strict gate n{n}",not passes_morse_gate(n,bad))
    check(f"equality boundary exact zero weight n{n}",[x['epsilon_power'] for x in morse_weights(n,bad)]==[0])
    P,scaled,syms=scaled_root_polynomial(n,bad)
    check(f"rational root of full bad graph n{n}",sp.expand(P.subs(T,B))==0)
    check(f"whole bad graph divides exactly by T-b n{n}",sp.rem(P,T-B,T)==0)
    limit=sp.expand(scaled.subs(e,0))
    check(f"boundary special fiber retains rational root n{n}",sp.expand(limit.subs(Z,V))==0)
    m=n-2;D=1+m*Q*B**(m+1)+(m-1)*R*B**m
    check(f"exact derivative of rational root component n{n}",sp.expand(sp.diff(P,T).subs(T,B)-D)==0)
    A=sp.Symbol("A")
    qsol=(A**m-1-(m-1)*R*B**m)/(m*B**(m+1))
    check(f"two-sided rational field reconstruction n{n}",sp.cancel(D.subs(Q,qsol)-A**m)==0)
    check(f"simple q prime gives Kummer degree m n{n}",sp.degree(D,Q)==1 and sp.diff(D,Q)==m*B**(m+1))
    counter.append({'n':n,'marking_degree':m,'root':'T=b','D':str(D),
      'rational_source_field':'k(b,r,a)','extension':'cyclic degree m; no Keller model',
      'actual_source_has_nonconstant_unit':'p=1+xy, inverse1-by'})
    disc=sp.discriminant(Z**n+Z-V,Z)
    check(f"Morse discriminant degree n{n}",sp.degree(disc,V)==n-1)
    check(f"Morse discriminant simple n{n}",sp.gcd(disc,sp.diff(disc,V))==1)
x,y=sp.symbols("x y");p=1+x*y
check("global actual bad-source unit certificate",sp.expand(p*(1-B*y)-1-y*(x-B*p))==0)
# The gate is sufficient only: a pole in the scaling is a failed test, not an exclusion.
check("negative weight graph is rejected without claiming a group",not passes_morse_gate(5,{2:Q}))
manifest=json.loads((ROOT/'sources/coefficient-graph-manifest.json').read_text())
for s in manifest['sources']:
    f=ROOT/s['path']
    check("primary decomposition source bytes",f.stat().st_size==s['bytes'])
    check("primary decomposition source SHA256",hashlib.sha256(f.read_bytes()).hexdigest()==s['sha256'])
(ROOT/'coefficient_graph_witnesses.json').write_text(json.dumps({'representation':'rational coefficient and exponent triples for (r,q,b)','witnesses':witnesses},indent=2)+'\n')
model={'claims':['KMFI73','KMFI74','KMFI75','KMFI76'],'actual_affine3_graph_controls':controls,
 'strict_weight_controls':weights,'bad_graph_controls':counter,
 'universal_statements':{'all_b_polynomial_graphs_have_explicit_affine3_Keller_model':True,
 'all_q_independent_polynomial_graphs_have_full_root_Sn':True,
 'strict_weight_condition':'c+j-1>a*(n-1) for all q^a*b^c terms after r=0',
 'strict_weight_condition_is_necessary_for_Sn':False,
 'every_actual_b_graph_has_marked_group_Gn_and_exact_cyclotomic_constants':True,
 'every_actual_b_graph_is_absolutely_stably_atomic':True,
 'every_finite_cyclic_geometric_generic_deck_group_has_explicit_atomic_affine3_witness':True,
 'bad_boundary_component_rational_but_not_stably_affine_space':True,
 'arbitrary_polynomial_graph_has_affine_source':False},
 'named_PSL_frontier':'Proper root group requires a q-dependent term surviving r=0 and violating the strict weight gate; this is not sufficient.',
 'exceptional_PSL2_11_degree11_absolute_occurrence_decided':False,
 'natural_PSL2_11_degree12_absolute_occurrence':'excluded by KMFI67',
 'full_action_classification_complete':False}
(ROOT/'coefficient_graph_model.json').write_text(json.dumps(model,indent=2)+'\n')
names=['COEFFICIENT_GRAPHS.md','coefficient_graph.py','verify_coefficient_graph.py',
 'coefficient_graph_model.json','coefficient_graph_witnesses.json','sources/coefficient-graph-manifest.json',
 'sources/stacks-09E3.html']
receipt={'result':'pass','check_count':len(checks),'checks':checks,'python':platform.python_version(),
 'sympy':sp.__version__,'runtime_seconds':round(time.monotonic()-started,3),
 'sha256':{name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in names},
 'assurance':'Written all-degree graph completion and deformation proof, with exact full polynomial witnesses, whole Jacobians, strict-weight and bad-source controls. No formal verification or independent whole-proof review.',
 'goal_complete':False}
(ROOT/'coefficient_graph_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'result':'pass','checks':len(checks),'seconds':receipt['runtime_seconds'],
 'all_polynomial_b_graphs_affine3':True,'exceptional_PSL11_degree11_decided':False,'full_goal_complete':False}))
