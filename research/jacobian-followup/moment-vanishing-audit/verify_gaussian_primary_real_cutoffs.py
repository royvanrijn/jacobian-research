"""Focused exact controls for the general PBW-induced coefficient backend.

Run before registration. Writes this packet's own examples/receipt only.
It does not rerun or overwrite any historical sealed producer or UNKNOWN job.
"""
from pathlib import Path
from fractions import Fraction
import copy, hashlib, json, subprocess, sys, time
import gaussian_primary_real_cutoff_backend as B

p = Path(__file__).resolve().parent
sha = lambda f:hashlib.sha256(Path(f).read_bytes()).hexdigest()
dump = lambda x:json.dumps(x,indent=2)+'\n'
pair = [['0','1'],['1','0']]
cases = [
    ('free_linear',dict(dimension=2,parameter_count=1,phase_terms=[[[1,0,1],'1']],parameter_ideal=[],precision_matrix=pair),[['1','0']]),
    ('dual_linear',dict(dimension=2,parameter_count=1,phase_terms=[[[1,0,0],'1'],[[0,1,1],'1']],parameter_ideal=[[[[2],'1']]],precision_matrix=pair),[['1','2']]),
    ('embedded_linear',dict(dimension=2,parameter_count=2,phase_terms=[[[1,0,0,0],'1'],[[0,1,0,1],'1']],parameter_ideal=[[[[0,2],'1']],[[[1,1],'1']]],precision_matrix=pair),[['1','2']]),
    ('algebraic_residue',dict(dimension=2,parameter_count=1,phase_terms=[[[1,0,1],'1']],parameter_ideal=[[[[2],'1'],[[0],'-2']]],precision_matrix=pair),[['1','0']]),
    ('nilpotent_linear',dict(dimension=2,parameter_count=1,phase_terms=[[[1,0,1],'1']],parameter_ideal=[[[[2],'1']]],precision_matrix=pair),[['0','1'],['1','0']]),
    ('free_quadratic',dict(dimension=2,parameter_count=1,phase_terms=[[[2,0,1],'1']],parameter_ideal=[],precision_matrix=pair),[['1/2','0']]),
]
root = p/'gaussian_primary_real_cutoff_examples'; root.mkdir(exist_ok=True)
controls = []; identity_count = 0; heat_count = 0; results = {}; positive_envelopes = 0; negative_moments = 0
for name,data,expected in cases:
    out = root/name; out.mkdir(exist_ok=True)
    (out/'input.json').write_text(dump(data))
    result = B.solve(data,out/'work',timeout=25)
    assert result['status'] == 'DECIDED', (name,result)
    assert result['region']['vertices'] == expected
    assert result['actual_image_kernel_computed'] is False
    assert result['general_presented_ring_backend_implemented'] is True
    for item in result['component_results']:
        component = item['finite_algebra']['component']
        A = B.O.FiniteAlgebra(data['parameter_count'],component,data['dimension'])
        assert A.manifest() == item['finite_algebra']
        family = dict(data,parameter_ideal=component['primary_ideal'])
        order = lambda m:B.G.heat_order(family,m)
        if 'certificate' in item:
            source = B.source_input(data,A)
            assert source == item['source_input']
            assert B.verify(source,item['certificate']); identity_count += 2
            law = B.order_law(source,item['certificate'],order)
            assert law == item['order_law']
            for step in item['region']['trace']:
                C,E = Fraction(step['C']),Fraction(step['E']); answer = step['result']
                assert answer['status'] == 'DECIDED'
                if answer['zero']:
                    assert answer == B.line_envelope(law,order,C,E)
                    positive_envelopes += 1
                else:
                    w = answer['witness']; m = w['power']; alpha = w['marker']
                    cs = B.s.symbols('c:'+str(data['parameter_count']))
                    gb = B.s.groebner([B.O.polynomial(g,cs) for g in component['primary_ideal']],*cs,domain=B.s.QQ)
                    original = B.O.original_moment(family,alpha,m)
                    recorded = B.O.polynomial(w['polynomial'],cs)
                    assert gb.reduce(recorded-original)[1] == 0 and gb.reduce(recorded)[1] != 0
                    assert m > C*sum(alpha)+E
                    negative_moments += 1
        else:
            assert not B.O.Moments(data,A).power(A.rank)
            assert item['original_orders'] == [order(m) for m in range(1,A.rank+1)]
    for m in range(1,9):
        expected_order = (m if name in ('free_linear','algebraic_residue') else
                          2*m if name == 'free_quadratic' else
                          (1 if m == 1 else None) if name == 'nilpotent_linear' else
                          (1 if m == 1 else m-2))
        assert B.G.heat_order(data,m) == expected_order,(name,m)
        heat_count += 1
    points = [(x['power'],x['least_marker_degree']) for x in result['region']['supporting_original_points']]
    assert all(B.G.heat_order(data,m) == q for m,q in points)
    assert [[str(C),str(E)] for C,E in B.G.vertices(Fraction(result['region']['asymptotic_coefficient']),points)] == expected
    (out/'result.json').write_text(dump(result)); results[name] = result
    controls.append(dict(name=name,status=result['status'],vertices=expected,
        pure_linear_minimum=result['region']['pure_linear_minimum'],
        component_lengths=[z['finite_algebra']['rank'] for z in result['component_results']],
        result_sha256=sha(out/'result.json')))

embedded = results['embedded_linear']
assert len(embedded['primary_components']) == 2
assert sorted(z['dimension'] for z in embedded['primary_components']) == [0,1]
assert results['dual_linear']['region']['pure_linear_minimum'] is None
assert results['nilpotent_linear']['region']['pure_linear_minimum'] == '1'
assert results['algebraic_residue']['component_results'][0]['finite_algebra']['rank'] == 2
assert results['free_quadratic']['component_results'][0]['certificate']['ordinary']['method'] == 'exhaustive QQ(U) graph search'

# Faithful derivative-left model: Y=M_Y-d/dzeta, D_Y=zeta.
item = results['dual_linear']['component_results'][0]
A = B.O.FiniteAlgebra(1,item['finite_algebra']['component'],2)
zero = (A.zero,)*A.rank
def add(*values):
    out = {}
    for v in values:
        for e,c in v.items(): out[e] = A.add(out.get(e,zero),c)
    return {e:c for e,c in out.items() if any(c)}
def scale(v,c): return {e:A.scale(z,c) for e,z in v.items() if any(A.scale(z,c))}
def derivative(v): return {e+1:c for e,c in v.items()}
def coordinate(v):
    out = {}; y = A.vector(A.dependent[0])
    for e,c in v.items():
        out[e] = A.add(out.get(e,zero),A.mul(y,c))
        if e: out[e-1] = A.add(out.get(e-1,zero),A.scale(c,-e))
    return {e:c for e,c in out.items() if any(c)}
induced_controls = 0
for e in range(4):
    for j in range(A.rank):
        v = {e:A.basis_vector(j)}
        assert add(derivative(coordinate(v)),scale(coordinate(derivative(v)),-1)) == v
        induced_controls += 1
seed = {0:A.unit}; nil = coordinate(seed)
assert nil and not coordinate(nil) and derivative(nil)
assert coordinate(coordinate(derivative(seed))) != {}  # I kills the seed, not the entire module.

# QQ(U) scalar projection and fractional late-line congruences.
projection,projection_audit = B.rational_component([[[1,0],'c0/(c0+1)'],[[0,1],'-c0/(c0+1)']],['c0'])
assert projection.as_expr() == B.S.X-B.S.Y
assert projection_audit['cleared_denominator'] == 'c0 + 1'
assert projection_audit['selected_parameter_monomial'] == [1]
law = dict(marker_slope='1/2',exact_line_threshold=4,late_offsets=['-1/2','0'])
order = lambda m:m//2
assert B.line_envelope(law,order,2,1)['late_line_bounds'][0]['first_admissible_power'] == 5
assert B.line_envelope(law,order,2,0) is None

# General graph search also detects unsafe real linear phases.
boundary_controls = []
for name,data,status in [
    ('zero_ring',dict(cases[0][1],parameter_ideal=[[[[0],'1']]]),'DECIDED'),
    ('zero_phase',dict(cases[0][1],phase_terms=[]),'DECIDED'),
    ('unsafe_real',dict(dimension=1,parameter_count=1,phase_terms=[[[1,0],'1']],parameter_ideal=[[[[1],'1']]],precision_matrix=[['1']]),'UNSAFE')]:
    out = root/name; out.mkdir(exist_ok=True);(out/'input.json').write_text(dump(data))
    result = B.solve(data,out/'work',timeout=25)
    assert result['status'] == status,(name,result)
    if status == 'UNSAFE':
        assert result['component_results'][0]['order_law']['marker_slope'] == '0'
        assert B.verify(result['component_results'][0]['source_input'],result['component_results'][0]['certificate'])
        identity_count += 2
    (out/'result.json').write_text(dump(result))
    boundary_controls.append(dict(name=name,status=status,result_sha256=sha(out/'result.json')))

# Corruptions must fail against ORIGINAL source generators.
item = results['dual_linear']['component_results'][0]; rejections = []
def reject(name,mutate):
    data,cert = copy.deepcopy(item['source_input']),copy.deepcopy(item['certificate']);mutate(data,cert)
    try: B.verify(data,cert)
    except (ValueError,TypeError): rejections.append(name)
    else: raise AssertionError(name)
reject('wrong specialization start',lambda d,c:c.__setitem__('specialization_start',1))
reject('zero ordinary scalar',lambda d,c:c['ordinary'].__setitem__('scalar_terms',[]))
reject('wrong phase',lambda d,c:d['phase_terms'][0].__setitem__(1,'2'))
reject('changed coefficient relation',lambda d,c:d['coefficient_relations'][0][0].__setitem__(1,'2'))
def change_witness(d,c):
    w = next(v for v in c['ordinary']['left_witnesses'] if v);w[0][1] = str(B.s.sympify(w[0][1])+1)
reject('left witness changed',change_witness)
reject('nonraising Euler remainder',lambda d,c:c['Euler']['operator_terms'].append([[0]*7,'1']))
budget = B.O.Budget(0);budget.end = time.monotonic()-1
try: B.search_identity(item['source_input'],'ordinary',budget)
except B.O.BudgetExpired: rejections.append('expired graph budget')
else: raise AssertionError('expired graph budget')

# The CLI uses the same complete presented-ring route, with independent output.
out = root/'cli';out.mkdir(exist_ok=True)
proc = subprocess.run([sys.executable,str(p/'gaussian_primary_real_cutoff_backend.py'),str(root/'dual_linear/input.json'),str(out/'result.json'),'--timeout','25'],capture_output=True,text=True,timeout=35)
assert proc.returncode == 0,proc.stderr
assert json.loads((out/'result.json').read_text())['region']['vertices'] == [['1','2']]
runtime = dict(python=sys.version,sympy=B.s.__version__,Macaulay2_version=subprocess.run(['M2','--version'],capture_output=True,text=True,timeout=10,check=True).stdout.strip(),
    imported_source_sha256={f:sha(p/f) for f in ['gaussian_real_cutoff_producer.py','gaussian_ordinary_cutoff_producer.py','gaussian_family_producer.py','gaussian_generic_cutoff_producer.py','algebraic_bernstein_generator.py','bernstein_slope_producer.py']},
    primary_decomposition='Trusted CAS,with exact intersection and prime-radical checks in retained transcript',
    transitive_environment_fully_pinned=False)
(p/'gaussian_primary_real_cutoff_runtime.json').write_text(dump(runtime))
files = {str(f.relative_to(p)):sha(f) for f in sorted(root.rglob('*')) if f.is_file()}
packet = dict(schema_version=1,universal_claims_use_written_proofs=True,
    examples=[dict(name=name,input_file=str((root/name/'input.json').relative_to(p)),
                   result_file=str((root/name/'result.json').relative_to(p)),
                   result_sha256=sha(root/name/'result.json')) for name,data,vertices in cases],
    boundary_controls=boundary_controls,primary_decomposition_trusted_CAS=True,
    source_contract='D_Y/(D_Y I)=direct_sum D_Y^alpha tensor T;T injects at derivative degree zero.',
    specialization_contract='Both exact unlocalized graph identities specialize from m=0;step start excludes ordinary rational-component roots.',
    scalar_projection_control=projection_audit,line_envelope_requires_certified_law=True,
    exact_graph_identity_count=identity_count,general_presented_ring_backend_implemented=True,
    retained_example_files_sha256=files,independent_replay=False,formal_verification=False,external_review=False)
(p/'gaussian_primary_real_cutoff_certificate.json').write_text(dump(packet))
receipt = dict(completed=True,kind='finite exact controls plus written universal PBW/source/termination proof',
    canonical_proof_sha256=sha(p/'PBW_INDUCED_GAUSSIAN_CERTIFICATES_OVER_PRESENTED_RINGS.md'),
    producer_sha256=sha(p/'gaussian_primary_real_cutoff_backend.py'),checker_sha256=sha(Path(__file__)),
    runtime_sha256=sha(p/'gaussian_primary_real_cutoff_runtime.json'),
    certificate_sha256=sha(p/'gaussian_primary_real_cutoff_certificate.json'),
    end_to_end_controls=controls,boundary_controls=boundary_controls,
    exact_graph_identity_count=identity_count,original_heat_order_controls=heat_count,
    positive_line_envelope_controls=positive_envelopes,original_negative_moment_controls=negative_moments,
    induced_Weyl_commutator_controls=induced_controls,ordinary_nilpotent_seed_retained=True,
    nonsplit_residue_field_retained=True,embedded_primary_component_retained=True,
    scalar_projection_control=projection_audit,fractional_line_congruence_control=True,
    rejection_controls=rejections,CLI_control=True,retained_example_files_sha256=files,
    general_presented_ring_backend_implemented=True,primary_decomposition_trusted_CAS=True,
    arbitrary_coefficient_kernel_oracle_supported=False,actual_image_kernel_computed=False,
    old_terminal_UNKNOWN_jobs_restarted=False,polynomial_convergent_novelty_claimed=False,
    independent_replay=False,formal_verification=False,external_review=False)
(p/'gaussian_primary_real_cutoff_verification.json').write_text(dump(receipt))
print(json.dumps({k:receipt[k] for k in ['completed','exact_graph_identity_count','original_heat_order_controls','positive_line_envelope_controls','original_negative_moment_controls','induced_Weyl_commutator_controls']}))
