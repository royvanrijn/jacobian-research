"""Focused family/normalization controls; --generate executes retained CAS jobs."""
from pathlib import Path
import argparse
import ast
import hashlib
import json
import shutil
import subprocess
import sys

import sympy as sp
import gaussian_family_producer as G

p=Path(__file__).resolve().parent
root=p/'gaussian_family_examples'
sha=lambda f:hashlib.sha256(Path(f).read_bytes()).hexdigest()
dump=lambda x:json.dumps(x,indent=2)+'\n'
cases=[
    ('universal_affine',G.universal_family(1,1),['c0','c1'],[],False,[0,1,2]),
    ('singular_safe_cusp',dict(dimension=1,parameter_count=2,phase_terms=[[[0,3,0],'-1'],[[0,0,2],'1']]),['c0^3-c1^2'],[],False,[0,1]),
    ('zero_divisor_normalization',dict(dimension=1,parameter_count=2,phase_terms=[[[1,1,0],'1']],parameter_ideal=[[[[1,1],'1']]]),['c0'],['c1'],False,[1,1,2]),
    ('nilpotent_nonpure_safe',dict(dimension=1,parameter_count=1,phase_terms=[[[1,1],'1']],parameter_ideal=[[[[3],'1']]]),['c0'],['1'],True,[1]),
    ('split_safe_cubic',dict(dimension=2,parameter_count=1,phase_terms=[[[3,0,1],'1']],precision_matrix=[[0,1],[1,0]]),[],['1'],True,[0]),
    ('ordinary_thick_torsion',dict(dimension=1,parameter_count=2,phase_terms=[[[1,1,0],'1']],parameter_ideal=[[[[2,2],'1']]]),['c0'],['c1^2'],False,[1,1,2]),
    ('algebraic_coefficient',dict(dimension=2,parameter_count=1,phase_terms=[[[1,0,0],'1'],[[0,1,1],'1']],parameter_ideal=[[[[2],'1'],[[0],'1']]]),['c0^2+1'],['1'],True,[1]),
    ('zero_algebra',dict(dimension=1,parameter_count=1,phase_terms=[[[1,1],'1']],parameter_ideal=[[[[0],'1']]]),['1'],['1'],True,[]),
    ('zero_phase',dict(dimension=1,parameter_count=1,phase_terms=[]),[],['1'],True,[0])]


def ideal_basis(values,r,encoded=True):
    cs=sp.symbols('c:'+str(r))
    if encoded:
        values=[sum(sp.Rational(c)*sp.prod(x**e for x,e in zip(cs,es)) for es,c in f) for f in values]
    else:
        values=[sp.sympify(x.replace('^','**'),locals={str(c):c for c in cs}) for x in values]
    if not any(values):return ()
    return tuple(sp.groebner(values,*cs,domain=sp.QQ).polys)


def helper_script():
    return f'''load "{G.CORE.as_posix()}";
W=QQ[x,y,dx,dy,WeylAlgebra=>{{x=>dx,y=>dy}}];
C=W^1/ideal(y-x^2,dx+2*x*dy);
O=W^1/ideal(dx,dy); B=QQ[x,y]; S=ideal(y-x^2);
r1=gGenericRank(C,S);r2=gGenericRank(C++C,S);r3=gGenericRank(C++O,S);
assert(r1#0==1 and r2#0==2 and r3#0==1);
W=QQ[t,c0,dt,dc0,WeylAlgebra=>{{t=>dt,c0=>dc0}}];
O=W^1/ideal(dt,dc0); Delta=W^1/ideal(dt,c0);
L1=localCohom({{0,1}},ideal(c0),O,Strategy=>OaTa);
L2=localCohom({{0,1}},ideal(c0),Delta,Strategy=>OaTa);
assert(gIsZero(L1#0) and not gIsZero(L1#1));
assert(not gIsZero(L2#0) and gIsZero(L2#1));
B=QQ[c0];
Graph=W^1/ideal(W_1-W_0^2,W_2+2*W_0*W_3);
bad=gSourceBad(Graph,ideal(0_B));
assert(#bad==1 and first bad==ideal(c0));
assert(not isHolonomic(W^1/ideal(W_0*W_2-W_1)));
W=QQ[x,dx,WeylAlgebra=>{{x=>dx}}];
K=Drestriction(ideal(x^2*dx+1),{{1}});
assert(gIsZero(K#0) and gIsZero(K#1));
H=new HashTable from {{
  "nonlinear_conormal_ranks"=>{{r1#0,r2#0,r3#0}},
  "ordinary_restriction_uses_codimension_cohomology"=>true,
  "punctual_source_degree_zero_would_be_wrong"=>true,
  "vertical_graph_critical_parameter"=>gSerializeIdeal first bad,
  "central_Euler_is_not_jointly_holonomic"=>true,
  "irregular_origin_restriction_is_zero"=>true}};
stdio << "HELPERS_COMPLETE|" << toJSON H << endl;
exit 0;
'''


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--generate',action='store_true')
    args=parser.parse_args()
    root.mkdir(exist_ok=True)
    controls=[];files={}
    for name,data,safe,normal,is_safe,codims in cases:
        out=root/name
        if args.generate:
            out.mkdir(exist_ok=True);(out/'input.json').write_text(dump(data))
            reusable=False
            if (out/'result.json').exists() and (out/'cas/family.m2').exists():
                prior=json.loads((out/'result.json').read_text())
                reusable=(prior.get('core_sha256')==sha(G.CORE) and prior.get('input')==data
                          and (out/'cas/family.m2').read_text()==G.make_script(data))
            if not reusable:
                result=G.run(data,out/'cas',20);(out/'result.json').write_text(dump(result))
        result=json.loads((out/'result.json').read_text())
        assert json.loads((out/'input.json').read_text())==data
        assert result['status']=='DECIDED'
        parsed=G.parse_output((out/'cas/family.stdout').read_text())
        assert all(result[k]==v for k,v in parsed.items())
        for f,h in result['cas_files_sha256'].items():assert sha(out/'cas'/f)==h
        assert result['core_sha256']==sha(G.CORE)
        r=data['parameter_count']
        assert ideal_basis(result['safe_radical_ideal'],r)==ideal_basis(safe,r,False)
        assert ideal_basis(result['polynomial_Hom_normalization_ideal_lift'],r)==ideal_basis(normal,r,False)
        assert result['coefficient_ring_phase_safe'] is is_safe
        assert sorted(n['codimension'] for n in result['nodes'])==codims
        assert result['normalization_colon_stable'] and result['normalization_torsion_exponent']>=1
        for node in result['nodes']:
            assert node['generic_polynomial_Hom_dimension'] in [0,1]
            assert node['whole_closure_safe']==(node['generic_polynomial_Hom_dimension']==1)
            if node['whole_closure_safe']:assert not node['children']
            for child in node['children']:
                matches=[n for n in result['nodes'] if ideal_basis(n['stratum_closure_ideal'],r)==ideal_basis(child,r)]
                assert len(matches)==1 and matches[0]['dimension']<node['dimension']
        if name=='ordinary_thick_torsion':assert result['normalization_torsion_exponent']==2
        if name=='nilpotent_nonpure_safe':assert result['normalization_torsion_exponent']==3
        controls.append(dict(name=name,safe_ideal=result['safe_radical_ideal'],normalization_ideal=result['polynomial_Hom_normalization_ideal_lift'],
            coefficient_ring_safe=is_safe,normalization_torsion_exponent=result['normalization_torsion_exponent'],
            node_dimensions=[n['dimension'] for n in result['nodes']]))

    helpers=root/'geometry_controls';helpers.mkdir(exist_ok=True)
    if args.generate:
        (helpers/'controls.m2').write_text(helper_script())
        run=subprocess.run(['M2','--no-readline','--silent','--stop',str(helpers/'controls.m2')],capture_output=True,text=True,timeout=20)
        (helpers/'controls.stdout').write_text(run.stdout+run.stderr)
        assert run.returncode==0
    hs=(helpers/'controls.stdout').read_text()
    line=next(l for l in hs.splitlines() if l.startswith('HELPERS_COMPLETE|'))
    helper=json.loads(line.split('|',1)[1]);assert helper['nonlinear_conormal_ranks']==[1,2,1]

    # Preserve the already terminal prototype cap; do not restart it.
    cap=root/'capped_cusp_degeneration'
    if args.generate:
        original=p.parents[2]/'research/.scratch/gaussian_family_20261116/unsafe_cusp_stratum'
        cap.mkdir(exist_ok=True)
        for name in ['family.m2','family.stdout','core_at_timeout.m2','timeout_manifest.json']:
            shutil.copyfile(original/name,cap/name)
        result=json.loads((original.parent/'unsafe_cusp_stratum.json').read_text())
        (cap/'result.json').write_text(dump(result))
    capped=json.loads((cap/'result.json').read_text());manifest=json.loads((cap/'timeout_manifest.json').read_text())
    assert capped['status']=='UNKNOWN' and manifest['process_terminal'] and manifest['timeout_seconds']==30
    assert sha(cap/'core_at_timeout.m2')==manifest['core_snapshot_sha256']

    rejections=[]
    for data in [dict(dimension=0,parameter_count=1,phase_terms=[]),
                 dict(dimension=1,parameter_count=0,phase_terms=[]),
                 dict(dimension=1,parameter_count=1,phase_terms=[[[1],'1']]),
                 dict(dimension=1,parameter_count=1,phase_terms=[],precision_matrix=[[0]]),
                 dict(dimension=2,parameter_count=1,phase_terms=[],precision_matrix=[[1,1],[0,1]]),
                 dict(dimension=1,parameter_count=1,phase_terms=[[[1,1],0.5]])]:
        try:G.normalize(data)
        except ValueError:rejections.append(data)
        else:raise AssertionError('accepted invalid family input')
    for output in ['', 'GAUSSIAN_FAMILY_COMPLETE', 'RESULT_JSON|{}\nGAUSSIAN_FAMILY_COMPLETE']:
        try:G.parse_output(output)
        except ValueError:rejections.append('incomplete CAS output')
        else:raise AssertionError('accepted incomplete CAS output')
    universal_checks=[]
    import math
    for n in range(1,4):
        for d in range(4):
            data=G.universal_family(n,d)
            assert data['parameter_count']==math.comb(n+d,d)
            assert len({tuple(e[:n]) for e,c in data['phase_terms']})==len(data['phase_terms'])
            G.normalize(data);universal_checks.append([n,d,data['parameter_count']])
    if args.generate:
        out=p.parents[2]/'research/.scratch/gaussian_family_cli_20261116';out.mkdir(exist_ok=True)
        proc=subprocess.run([sys.executable,str(p/'gaussian_family_producer.py'),str(root/'split_safe_cubic/input.json'),
                             str(out/'result.json'),'--workdir',str(out/'cas')],capture_output=True,text=True,timeout=20)
        assert proc.returncode==0 and json.loads((out/'result.json').read_text())['coefficient_ring_phase_safe']
    for f in root.rglob('*'):
        if f.is_file():files[str(f.relative_to(p))]=sha(f)
    ast.parse((p/'gaussian_family_producer.py').read_text())
    receipt=dict(completed=True,canonical_proof_sha256=sha(p/'GENERAL_GAUSSIAN_FAMILY_LOCUS_AND_NORMALIZATION.md'),
        generator_sha256=sha(p/'gaussian_family_producer.py'),core_sha256=sha(G.CORE),checker_sha256=sha(Path(__file__)),
        runtime_sha256=sha(p/'gaussian_family_runtime.json'),family_controls=controls,geometry_controls=helper,
        capped_control=dict(status='UNKNOWN',timeout_seconds=30,restarted=False,historical_core_revision=True),
        universal_full_space_controls=universal_checks,rejection_controls=rejections,CLI_executed_in_initial_generation=True,
        retained_example_files_sha256=files,assurance='small CAS executions and trace/ideal controls, not universal independent proof assurance',
        numerical_universal_cutoff_computed=False,ordinary_full_cutoff_ideal_producer_implemented=False)
    (p/'gaussian_family_producer_verification.json').write_text(dump(receipt))
    print(json.dumps(dict(completed=True,families=len(controls),geometry_controls=helper,universal_constructors=len(universal_checks),rejections=len(rejections),capped='UNKNOWN;not restarted',retained_files=len(files))))


if __name__=='__main__':main()
