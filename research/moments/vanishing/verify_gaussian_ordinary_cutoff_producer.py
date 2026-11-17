"""Focused ordinary-ideal, primary, rational IBP and least-integer controls."""
from pathlib import Path
import argparse,ast,copy,hashlib,json,shutil,subprocess,sys
import sympy as sp
import gaussian_ordinary_cutoff_producer as O

p=Path(__file__).resolve().parent;root=p/'gaussian_ordinary_cutoff_examples'
prototype=p.parents[2]/'research/.scratch/gaussian_ordinary_20261116'
sha=lambda f:hashlib.sha256(Path(f).read_bytes()).hexdigest()
linear=dict(dimension=1,parameter_count=2,phase_terms=[[[0,1,0],'1'],[[1,0,1],'1']])
constant=dict(dimension=1,parameter_count=1,phase_terms=[[[0,1],'1']])
embedded=dict(dimension=1,parameter_count=2,phase_terms=[[[0,1,0],'1']],parameter_ideal=[[[[2,0],'1']],[[[1,1],'1']]])
noncoordinate=dict(dimension=1,parameter_count=2,phase_terms=[[[1,2,0],'1'],[[1,0,2],'1']],parameter_ideal=[[[[4,0],'1'],[[2,2],'2'],[[0,4],'1']]])
split=dict(dimension=2,parameter_count=1,phase_terms=[[[3,0,1],'1']],precision_matrix=[[0,1],[1,0]])
cases=[('linear_pure',linear,0,0,True,['c0','c1**2'],2),
       ('linear_cutoff',linear,1,0,False,['c0','c1**2'],2),
       ('linear_affine',linear,1,1,False,['c0**2+c1**2','c0*c1**2','c1**4'],3),
       ('linear_zero_slope',linear,0,1,False,['c0**2','c0*c1','c1**2'],4),
       ('constant_thickness',constant,3,3,False,['c0**4'],4),
       ('embedded_constant',embedded,1,0,True,['c0'],1),
       ('noncoordinate_nilpotent',noncoordinate,1,0,True,['(c0**2+c1**2)**2'],0),
       ('split_cubic',split,1,0,False,[],0),
       ('zero_algebra',dict(dimension=1,parameter_count=1,phase_terms=[[[1,1],'1']],parameter_ideal=[[[[0],'1']]]),1,0,False,['1'],0)]
least_cases=[('nilpotent_nonpure_safe',False,False,'DECIDED',2),('nilpotent_nonpure_safe',True,False,'NO_LINEAR_BOUND',None),
             ('split_safe_cubic',True,False,'DECIDED',1),('split_safe_cubic',True,True,'DECIDED',1),
             ('universal_affine',False,True,'DECIDED',0),('universal_affine',False,False,'NO_FINITE_BOUND',None)]


def ideal_equal(encoded,expected,r):
    cs=sp.symbols('c:'+str(r));a=[O.polynomial(f,cs) for f in encoded]
    b=[sp.sympify(f,locals={str(c):c for c in cs}) for f in expected]
    return sp.groebner(a,*cs,domain=sp.QQ)==sp.groebner(b,*cs,domain=sp.QQ)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--generate',action='store_true');args=parser.parse_args()
    root.mkdir(exist_ok=True);results={};controls=[]
    for name,data,N,E,pure,expected,prefix in cases:
        out=root/name
        if args.generate:
            out.mkdir(exist_ok=True);(out/'input.json').write_text(O.dump(data))
            result=O.run(data,out/'cas',N,E,pure,20)
            if result['status']=='DECIDED':result['least_generating_prefix']=O.least_generating_prefix(result)
            (out/'result.json').write_text(O.dump(result))
        result=json.loads((out/'result.json').read_text());results[name]=result
        assert result['status']=='DECIDED' and result['input']==data and O.verify_result(result)
        assert ideal_equal(result['ideal_lift'],expected,data['parameter_count'])
        assert result['least_generating_prefix']==O.least_generating_prefix(result)
        assert result['least_generating_prefix']['least_total_index_prefix']==prefix
        for i,step in enumerate(result.get('stopping_trace',[])):
            oracle=step['oracle']
            if 'primary_components' in oracle:
                stdout=(out/'cas'/('step_'+str(i))/'primary/primary.stdout').read_text()
                components=json.loads(next(l.split('|',1)[1] for l in stdout.splitlines() if l.startswith('PRIMARY_JSON|')))
                assert components==oracle['primary_components']
        controls.append(dict(name=name,ideal=result['ideal_lift'],least_total_index_prefix=prefix))
    least=[]
    for name,pure,fields,status,N in least_cases:
        tag=name+('_pure' if pure else '_affine')+('_fields' if fields else '_ring');out=root/('least_'+tag)
        if args.generate:
            out.mkdir(exist_ok=True)
            numerical=json.loads((p/'gaussian_generic_cutoff_examples'/('family_'+name)/'result.json').read_text())
            result=O.least_bound(numerical['input'],out/'cas',pure,fields,20,numerical)
            (out/'result.json').write_text(O.dump(result))
        result=json.loads((out/'result.json').read_text());results['least_'+tag]=result
        assert result['status']==status and result.get('least_integer_coefficient')==N and O.verify_least_result(result)
        least.append(dict(name=tag,status=status,least_integer_coefficient=N))
    geometry=root/'embedded_primary'
    if args.generate:
        geometry.mkdir(exist_ok=True)
        components=O.primary_components(embedded['parameter_ideal'],2,geometry/'cas',O.Budget(10))
        tested=[O.component_test(embedded,c,1,0,True,O.Budget(10)) for c in components]
        w=dict(marker=[0],power=1,polynomial=[[[1,0],'1']])
        oracle=dict(status='DECIDED',zero=False,primary_components=components,components=tested,witness=w)
        (geometry/'result.json').write_text(O.dump(oracle))
    embedded_oracle=json.loads((geometry/'result.json').read_text())
    assert [c['zero'] for c in embedded_oracle['components']]==[True,False]
    assert O.verify_zero_test(embedded,embedded['parameter_ideal'],embedded_oracle,1,0,True)
    ibp=root/'original_IBP';ibp.mkdir(exist_ok=True)
    if args.generate:
        for name in ['pure_quadratic_IBP.json','direct_zero_cone_IBP.json','zero_cone_IBP.json']:
            shutil.copyfile(prototype/name,ibp/name)
        (ibp/'cap_manifest.json').write_text(O.dump(dict(status='UNKNOWN',timeout_seconds=20,terminal=True,restarted=False,
            historical_source_snapshot_available=False,assurance='Cap is retained as performance boundary, not as a successful computation; no historical full source hash was captured.')))
    quadratic=json.loads((ibp/'pure_quadratic_IBP.json').read_text());quadratic_data=dict(dimension=1,parameter_count=1,phase_terms=[[[2,0],'1']])
    assert O.verify_component(quadratic_data,quadratic,0,0,True) and not quadratic['zero']
    cone=json.loads((ibp/'direct_zero_cone_IBP.json').read_text());A=O.FiniteAlgebra(1,cone['component'],1)
    assert O.verify_ibp(cone['input'],A,cone['certificate'],0,0,False)
    operator={tuple(e):A.parse_coefficient(c) for e,c in cone['certificate']['operator']}
    assert O.indicial(operator,A)==cone['origin'] and cone['origin']['bound']==1
    assert len(cone['certificate']['identity_terms'])==16
    cap=json.loads((ibp/'cap_manifest.json').read_text())
    assert json.loads((ibp/'zero_cone_IBP.json').read_text())['status']=='UNKNOWN' and cap['terminal'] and not cap['restarted']
    rejected=[]
    def reject(name,fn):
        try:fn()
        except (ValueError,TypeError,KeyError,IndexError,ZeroDivisionError):rejected.append(name)
        else:raise AssertionError('accepted false certificate: '+name)
    damaged=copy.deepcopy(results['linear_pure']);damaged['ideal_lift']=[[[[1,0],'1']],[[[0,1],'1']]]
    reject('ordinary ideal replaced by radical',lambda:O.verify_result(damaged))
    damaged=copy.deepcopy(results['linear_zero_slope']);damaged['pairing_certificate']['coefficient_from_moment_matrix'][0][0]='0'
    reject('wrong Wick pairing inverse',lambda:O.verify_result(damaged))
    damaged=copy.deepcopy(embedded_oracle);damaged['primary_components'].pop();damaged['components'].pop()
    reject('lost embedded primary component',lambda:O.verify_zero_test(embedded,embedded['parameter_ideal'],damaged,1,0,True))
    damaged=copy.deepcopy(quadratic);damaged['zero']=True;damaged['witness']=None
    reject('nonzero actual scalar period declared zero',lambda:O.verify_component(quadratic_data,damaged,0,0,True))
    damaged=copy.deepcopy(quadratic);damaged['origin_indicial']['bound']=0
    reject('wrong scalar origin threshold',lambda:O.verify_component(quadratic_data,damaged,0,0,True))
    damaged=copy.deepcopy(quadratic);damaged['zero_series_operator']=damaged['IBP_certificate']['operator']
    reject('failure to subtract the pure unit',lambda:O.verify_component(quadratic_data,damaged,0,0,True))
    damaged=copy.deepcopy(cone['certificate']);damaged['identity_terms'].pop(0)
    reject('missing original rational primitive',lambda:O.verify_ibp(cone['input'],A,damaged))
    damaged=copy.deepcopy(cone['certificate']);damaged['denominator_polynomial']=[]
    reject('wrong rational localization polynomial',lambda:O.verify_ibp(cone['input'],A,damaged))
    damaged=copy.deepcopy(cone['certificate']);damaged['operator']=[]
    reject('zero scalar annihilator',lambda:O.verify_ibp(cone['input'],A,damaged))
    damaged=copy.deepcopy(results['least_split_safe_cubic_pure_ring']);damaged['ordered_cutoff_tests'].pop(0)
    reject('missing smaller integer obstruction',lambda:O.verify_least_result(damaged))
    damaged=copy.deepcopy(results['least_split_safe_cubic_pure_fields']);damaged['target_ideal']=[[[[0],'1']]]
    reject('wrong field-safe target',lambda:O.verify_least_result(damaged))
    damaged=copy.deepcopy(results['least_split_safe_cubic_pure_ring']);damaged['numerical_result']['all_safe_field_fibers_cutoff']['slope']='0'
    reject('understated generic numerical bound',lambda:O.verify_least_result(damaged))
    for value in [-1,sp.Rational(1,2),True]:
        reject('invalid cutoff '+str(value),lambda:O.normalized_input(linear,value,0,False))
    if args.generate:
        cli=root/'CLI';cli.mkdir(exist_ok=True)
        proc=subprocess.run([sys.executable,str(p/'gaussian_ordinary_cutoff_producer.py'),str(root/'constant_thickness/input.json'),
            str(cli/'ordinary.json'),'--slope','3','--intercept','3','--prefix','--workdir',str(cli/'ordinary')],capture_output=True,text=True,timeout=10)
        assert proc.returncode==0,proc.stdout+proc.stderr
        (cli/'split_input.json').write_text(O.dump(split))
        proc=subprocess.run([sys.executable,str(p/'gaussian_ordinary_cutoff_producer.py'),str(cli/'split_input.json'),str(cli/'least.json'),
            '--least','--pure','--numerical-result',str(p/'gaussian_generic_cutoff_examples/family_split_safe_cubic/result.json'),
            '--workdir',str(cli/'least')],capture_output=True,text=True,timeout=10)
        assert proc.returncode==0,proc.stdout+proc.stderr
    assert O.verify_result(json.loads((root/'CLI/ordinary.json').read_text()))
    assert O.verify_least_result(json.loads((root/'CLI/least.json').read_text()))
    for name in ['gaussian_ordinary_cutoff_producer.py','verify_gaussian_ordinary_cutoff_producer.py']:ast.parse((p/name).read_text())
    files={str(f.relative_to(p)):sha(f) for f in root.rglob('*') if f.is_file()}
    receipt=dict(completed=True,canonical_proof_sha256=sha(p/'ORDINARY_GAUSSIAN_CUTOFF_IDEALS_AND_LEAST_INTEGER_BOUNDS.md'),
        generator_sha256=sha(p/'gaussian_ordinary_cutoff_producer.py'),core_sha256=sha(O.G.CORE),checker_sha256=sha(Path(__file__)),runtime_sha256=sha(p/'gaussian_ordinary_runtime.json'),
        ordinary_ideal_controls=controls,least_integer_controls=least,embedded_primary_control=True,scalar_IBP_search_control=True,
        supplied_rational_cone_IBP_terms=16,origin_bound=1,rejection_controls=rejected,CLI_routes=2,
        cone_enumerator_cap=cap,historical_cusp_restarted=False,retained_example_files_sha256=files,
        ordinary_full_cutoff_ideal_producer_implemented=True,least_integer_and_generating_prefix_producers_implemented=True,
        least_real_slopes_claimed=False,independent_replay=False,formal_verification=False,external_review=False,
        assurance_boundary='Finite original identities, ordinary intersections and actual jets replay; primary status/independence and earlier complete family decomposition remain trusted CAS. Unlimited termination is the written holonomic/ACC proof, not inferred from these examples.')
    (p/'gaussian_ordinary_cutoff_producer_verification.json').write_text(O.dump(receipt))
    print(json.dumps(dict(completed=True,ideals=len(controls),least_integer_controls=len(least),rejections=len(rejected),retained_files=len(files),cone_search='UNKNOWN; not restarted')))


if __name__=='__main__':main()
