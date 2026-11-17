"""Small exact graph-certificate controls,not an all-input computational proof."""
from pathlib import Path
import copy
import hashlib
import json
import subprocess
import sys
import algebraic_bernstein_generator as G

p=Path(__file__).resolve().parent
proof=p/'OAKU_GRAPH_CERTIFICATES_AND_ACTUAL_POLYNOMIAL_DECISION.md'
sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
dump=lambda x:json.dumps(x,indent=2)+'\n'
root=p/'algebraic_bernstein_examples';root.mkdir(exist_ok=True)
cases=[
    ('unary_homogeneous',dict(dimension=1,phase_terms=[[[2],'1']],symbol_terms=[[[3],'1']]),True,'3/2'),
    ('unary_inhomogeneous',dict(dimension=1,phase_terms=[[[2],'1'],[[1],'1']],symbol_terms=[[[3],'1']]),True,'3/2'),
    ('inactive_direction',dict(dimension=2,phase_terms=[[[1,0],'1']],symbol_terms=[[[0,1],'1']]),True,None),
    ('coupled_phase_symbol',dict(dimension=2,phase_terms=[[[1,1],'1']],symbol_terms=[[[1,1],'1']]),False,'1')]
controls=[];files={}
for name,input_data,pure,rho in cases:
    out=root/name;out.mkdir(exist_ok=True)
    (out/'input.json').write_text(dump(input_data))
    pipeline=G.decide(input_data,out/'cas',20,all_positive_pure=pure)
    certificate=pipeline['certificate'];decision=pipeline['decision']
    assert certificate['exact_rational_Weyl_witness_verification'] is True
    assert G.verify_certificate(certificate)
    if rho is None:
        assert decision['alternative']=='eventually infinite first-survival orders'
        assert certificate['eta']['torsion_polynomial_terms']==[[[1],'1']]
    else:assert decision['selected_rational_slope']==rho
    assert decision['all_marker_eventual_safety']==(pure is True)
    if pure:assert type(decision['pure_linear_cutoff'])==int
    (out/'certificate.json').write_text(dump(certificate));(out/'pipeline.json').write_text(dump(pipeline))
    for f in out.rglob('*'):
        if f.is_file():files[str(f.relative_to(p))]=sha(f)
    controls.append(dict(name=name,slope=rho,decision=decision,specialization_thresholds=certificate['specialization_thresholds'],
        scalar_torsion_factor=certificate['eta']['torsion_polynomial_terms'],
        certificate_file=str((out/'certificate.json').relative_to(p)),certificate_sha256=sha(out/'certificate.json')))

elementary=[]
for name,phase,symbol,safe in [
    ('zero phase',[],[[[0],'1']],True),('zero operator',[[[2],'1']],[],True),
    ('invertible operator',[[[2],'1']],[[[0],'1'],[[1],'1']],False),
    ('constant phase',[[[0],'2']],[[[1],'1']],True)]:
    data=dict(dimension=1,phase_terms=phase,symbol_terms=symbol)
    result=G.decide(data,root/'unused_elementary',20)
    assert result['decision']['all_marker_eventual_safety']==safe
    elementary.append(dict(name=name,result=result['decision']))

# The direct verifier uses rational PBW multiplication,not Singular normal forms.
commutators=[]
A=G.Weyl(2);v=[A.var(i) for i in range(4)]
for i in range(2):
    for j in range(2):
        assert A.sub(A.mul(v[2+i],v[j]),A.mul(v[j],v[2+i]))==A.scalar(int(i==j))
        commutators.append([i,j])
f=A.add(A.power(v[0],2),v[1]);g=A.add(v[2],A.power(v[3],2))
assert A.mul(A.mul(f,g),f)==A.mul(f,A.mul(g,f))

original=json.loads((root/'unary_homogeneous/certificate.json').read_text());rejections=[]
def reject(name,mutate):
    c=copy.deepcopy(original);mutate(c)
    try:G.verify_certificate(c)
    except ValueError:rejections.append(name)
    else:raise AssertionError(name)
reject('scalar polynomial changed',lambda c:c['ordinary']['B_terms'].append([[0,0,0,0,0],'1']))
reject('original phase changed',lambda c:c['input']['phase_terms'].__setitem__(0,[[2],'2']))
reject('order constraint changed',lambda c:c['certificate_summary']['polynomial_terms'].append([[0,0],'1']))
reject('specialization threshold changed',lambda c:c['specialization_thresholds'].__setitem__('ordinary_graph_relation_start',1))
def change_witness(c):
    w=next(w for w in c['ordinary']['graph_witness'] if w)
    w[0][1]=str(G.Fraction(w[0][1])+1)
reject('left graph witness changed',change_witness)

# Run every documented CLI mode on the same original finite input.
tmp=p.parents[2]/'research/.scratch/algebraic_bernstein_cli_20261116';tmp.mkdir(parents=True,exist_ok=True)
cli=[]
for mode,flags,source in [('generate',[],root/'unary_homogeneous/input.json'),
                          ('decide',['--decide','--pure'],root/'unary_homogeneous/input.json'),
                          ('verify',['--verify'],root/'unary_homogeneous/certificate.json')]:
    output=tmp/(mode+'.json')
    command=[sys.executable,str(p/'algebraic_bernstein_generator.py'),str(source),str(output)]+flags
    result=subprocess.run(command,capture_output=True,text=True,timeout=30)
    assert result.returncode==0,result.stderr
    value=json.loads(output.read_text());assert value['completed']
    if mode=='generate':assert G.verify_certificate(value)
    elif mode=='decide':assert value['decision']['selected_rational_slope']=='3/2' and value['decision']['pure_linear_cutoff']>=1
    else:assert value['exact_rational_Weyl_witness_verification']
    cli.append(dict(operation=mode,output_sha256=sha(output)))

runtime=G.runtime_manifest();(p/'algebraic_bernstein_runtime.json').write_text(dump(runtime))
receipt=dict(completed=True,kind='exact finite graph certificates and end-to-end controls,not universal formal verification',
    canonical_proof_sha256=sha(proof),generator_sha256=sha(p/'algebraic_bernstein_generator.py'),checker_sha256=sha(Path(__file__)),
    credit_note_sha256=sha(p/'BERNSTEIN_CLASSICAL_CREDIT_AND_FORMAL_EXTENSIONS.md'),runtime_sha256=sha(p/'algebraic_bernstein_runtime.json'),
    graph_certificate_controls=controls,retained_example_files_sha256=files,
    elementary_boundary_controls=elementary,
    rational_Weyl_commutator_controls=commutators,associativity_control=True,rejection_controls=rejections,CLI_controls=cli,
    general_finite_QQ_certificate_generator_implemented=True,actual_image_kernel_computed=False,
    arbitrary_nonconvergent_formal_certificate_generator_implemented=False,all_ring_rational_generator_claimed=False,
    polynomial_convergent_novelty_claimed=False,classical_credit=['Sabbah','Maisonobe','Oaku','Mebkhout-Narvaez-Macarro'],
    universal_termination_and_specialization_verified_by='written proof with precise classical imports;finite identities checked separately',
    independent_replay=False,formal_verification=False,external_review=False)
(p/'algebraic_bernstein_generator_verification.json').write_text(dump(receipt))
print(json.dumps(dict(completed=True,graph_inputs=len(controls),Weyl_commutators=len(commutators),rejected_tampers=len(rejections),CLI=len(cli),example_files=len(files),actual_image_kernel=False)))
