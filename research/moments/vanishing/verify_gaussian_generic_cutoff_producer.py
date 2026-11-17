"""Focused generic-field original certificates; --generate prepares CAS artifacts."""
from pathlib import Path
import argparse
import ast
import copy
import hashlib
import json
import shutil
import subprocess
import sys

import sympy as sp
import gaussian_generic_cutoff_producer as G

p=Path(__file__).resolve().parent
root=p/'gaussian_generic_cutoff_examples'
prototype=p.parents[2]/'research/.scratch/gaussian_generic_20261116'
sha=lambda f:hashlib.sha256(Path(f).read_bytes()).hexdigest()
dump=lambda d:json.dumps(d,indent=2)+'\n'
cases=[('isotropic_linear',True,'1'),('generic_real_line',False,None),
       ('generic_split_cubic',True,'1/3'),('elliptic_split_phase',True,'1/2'),
       ('quartic_two_markers',False,None),('zero_generic_point',True,'0'),
       ('algebraic_isotropic',True,'1')]
family_cases=[('nilpotent_nonpure_safe','0','2'),('ordinary_thick_torsion','0','3'),
              ('split_safe_cubic','1/3','0'),('universal_affine','0','0')]


def import_generic(name,out):
    prior=json.loads((prototype/(name+'.json')).read_text());data=prior['input']
    out.mkdir(parents=True,exist_ok=True);cas=out/'cas';cas.mkdir(exist_ok=True)
    for f in ['generic.m2','generic.stdout']:shutil.copyfile(prototype/name/f,cas/f)
    assert (cas/'generic.m2').read_text()==G.make_script(data)
    source=json.loads(next(l.split('|',1)[1] for l in (cas/'generic.stdout').read_text().splitlines() if l.startswith('RESULT_JSON|')))
    family,n,r,phase,precision,field,generic=G.normalize(data)
    verified=G.verify_source(source,n,generic,precision,field)
    selected=G.actual_periods(source,verified,n,generic,precision,field)
    result=dict(prior,source=source,verified_original_certificate=verified,actual_period_selection=selected,
                core_sha256=sha(G.CORE),cas_files_sha256={f:sha(cas/f) for f in ['generic.m2','generic.stdout']})
    (out/'input.json').write_text(dump(data));(out/'result.json').write_text(dump(result))


def root_controls():
    field=G.GenericField(2,[[[[2,0],'1'],[[0,2],'1']]])
    x,y=field.parameters;z=sp.Symbol('z')
    assert G.Fraction(field,y*y,x*x)==-1
    controls=[]
    for expression,wanted in [(z-2,[2]),((z-2)*(z+3),[-3,2]),
                              ((3*z-2)*(z*z-2),[]),(z, [0])]:
        scale=G.Fraction(field,x+y,x)
        coefficients={int(e[0]):scale*c for e,c in sp.Poly(expression,z).terms()}
        roots,cert=field.integer_roots(coefficients)
        assert roots==wanted
        controls.append(dict(polynomial=str(expression),roots=roots,certificate=cert))
    roots,cert=field.integer_roots({1:field.one,0:G.Fraction(field,-x)})
    assert roots==[];controls.append(dict(polynomial='z-c0',roots=roots,certificate=cert))
    try:G.Fraction(field,1,x*x+y*y)
    except ZeroDivisionError:pass
    else:raise AssertionError('zero denominator modulo prime accepted')
    return controls


def rejection_controls(results,families):
    rejected=[]
    def reject(name,fn):
        try:fn()
        except (ValueError,ZeroDivisionError,KeyError,IndexError):rejected.append(name)
        else:raise AssertionError('accepted damaged certificate: '+name)
    base=results['quartic_two_markers']
    for name,mutate in [
        ('missing coordinate',lambda d:d['source']['coordinate_annihilators'].pop()),
        ('missing coordinate witness',lambda d:d['source']['coordinate_annihilator_left_witnesses'].pop()),
        ('missing relation row',lambda d:d['source']['relations'].pop()),
        ('wrong integration threshold',lambda d:d['source'].__setitem__('integration_cutoff',0)),
        ('lost original recurrence witness',lambda d:d['source']['original_recurrence_left_witnesses'].pop()),
        ('wrong generic prime contract',lambda d:d['source'].__setitem__('generic_prime_checked',False)),
        ('wrong actual safety decision',lambda d:d['actual_period_selection'].__setitem__('safe',True))]:
        damaged=copy.deepcopy(base);mutate(damaged)
        reject(name,lambda:G.verify_result(damaged))
    damaged=copy.deepcopy(results['isotropic_linear'])
    damaged['verified_original_certificate']['nonconstant_coefficient_denominators']=[]
    reject('lost specialization denominators',lambda:G.verify_result(damaged))
    damaged=copy.deepcopy(results['generic_split_cubic'])
    damaged['actual_period_selection']['cutoff']['slope']='0'
    reject('wrong numerical cutoff',lambda:G.verify_result(damaged))
    for name,d in families.items():
        if name not in ['nilpotent_nonpure_safe','ordinary_thick_torsion']:continue
        damaged=copy.deepcopy(d['family_result']);damaged['normalization_torsion_exponent']=1
        reject(name+' false ordinary torsion exponent',lambda:G.verify_family_torsion(d['input'],damaged))
    d=families['nilpotent_nonpure_safe'];damaged=copy.deepcopy(d['family_result'])
    damaged['coefficient_ring_phase_safe']=False
    reject('false normalization-one decision',lambda:G.verify_family_torsion(d['input'],damaged))
    reject('generic prime excludes parameter ideal',lambda:G.normalize(dict(
        dimension=1,parameter_count=1,phase_terms=[],parameter_ideal=[[[[1],'1']]],generic_prime=[])))
    reject('unit generic ideal',lambda:G.GenericField(1,[[[[0],'1']]]))
    reject('zero generic root polynomial',lambda:G.GenericField(1,[]).integer_roots({}))
    return rejected


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--generate',action='store_true')
    args=parser.parse_args();root.mkdir(exist_ok=True)
    controls=[];results={};families={}
    for name,safe,slope in cases:
        out=root/name
        if args.generate:import_generic(name,out)
        result=json.loads((out/'result.json').read_text());results[name]=result
        assert json.loads((out/'input.json').read_text())==result['input']
        assert G.verify_result(result) and result['core_sha256']==sha(G.CORE)
        for f,h in result['cas_files_sha256'].items():assert sha(out/'cas'/f)==h
        assert (out/'cas/generic.m2').read_text()==G.make_script(result['input'])
        source=json.loads(next(l.split('|',1)[1] for l in (out/'cas/generic.stdout').read_text().splitlines() if l.startswith('RESULT_JSON|')))
        assert source==result['source']
        selection=result['actual_period_selection'];assert selection['safe']==safe
        if safe:
            assert selection['cutoff']['slope']==slope and selection['cutoff']['intercept']=='0'
            assert selection['cutoff']['extends_to_ALL_specializations_of_component']
        controls.append(dict(name=name,safe=safe,slope=slope,rank=source['rank'],
            denominator_count=len(result['verified_original_certificate']['nonconstant_coefficient_denominators'])))
    assert results['quartic_two_markers']['source']['rank']==2
    assert results['isotropic_linear']['verified_original_certificate']['nonconstant_coefficient_denominators']
    assert results['elliptic_split_phase']['verified_original_certificate']['nonconstant_coefficient_denominators']
    for name,C,E in family_cases:
        out=root/('family_'+name)
        if args.generate:
            out.mkdir(exist_ok=True)
            prior=json.loads((p/'gaussian_family_examples'/name/'result.json').read_text());data=prior['input']
            result=G.family_bounds(data,out/'cas',30,prior)
            (out/'input.json').write_text(dump(data));(out/'result.json').write_text(dump(result))
        result=json.loads((out/'result.json').read_text());families[name]=result
        assert result['status']=='DECIDED'
        assert json.loads((out/'input.json').read_text())==result['input']
        assert result['family_result']==json.loads((p/'gaussian_family_examples'/name/'result.json').read_text())
        assert G.verify_family_torsion(result['input'],result['family_result'])==result['verified_ordinary_torsion']
        for component in result['component_certificates']:assert G.verify_result(component)
        bound=result['presented_ring_normalization_cutoff']
        assert bound['slope']==C and bound['intercept']==E
        for i,component in enumerate(result['component_certificates']):
            for f,h in component['cas_files_sha256'].items():assert sha(out/'cas'/('component_'+str(i))/f)==h
    roots=root_controls();rejects=rejection_controls(results,families)
    failures=root/'backend_boundaries'
    if args.generate:
        failures.mkdir(exist_ok=True)
        for f in ['field_test.m2','field_test.stdout','field_test.sing','field_test_singular.stdout',
                  'root_patch_test.m2','root_patch_test.stdout','generic_roots.m2','orientation.m2']:
            shutil.copyfile(prototype/f,failures/f)
        # This is a NEW rejected input, not a restart of any capped computation.
        bad=dict(dimension=1,parameter_count=1,phase_terms=[[[1,1],'1']],generic_prime=[[[[2],'1'],[[1],'-1']]])
        (failures/'nonprime.m2').write_text(G.make_script(bad))
        proc=subprocess.run(['M2','--no-readline','--silent','--stop',str(failures/'nonprime.m2')],capture_output=True,text=True,timeout=10)
        (failures/'nonprime.stdout').write_text(proc.stdout+proc.stderr)
        assert proc.returncode!=0 and 'must be PRIME' in proc.stdout+proc.stderr
        # Exercise both CLI routes, retaining all outputs.
        cli=root/'cli';cli.mkdir(exist_ok=True)
        for name,extra in [('generic_split_cubic',[]),('family_split_safe_cubic',
              ['--family','--family-result',str(p/'gaussian_family_examples/split_safe_cubic/result.json')])]:
            proc=subprocess.run([sys.executable,str(p/'gaussian_generic_cutoff_producer.py'),str(root/name/'input.json'),
                str(cli/(name+'.json')),'--workdir',str(cli/name)]+extra,capture_output=True,text=True,timeout=20)
            assert proc.returncode==0,proc.stdout+proc.stderr
    assert 'expected polynomial over QQ' in (failures/'field_test.stdout').read_text()
    assert 'only univarite minpoly allowed' in (failures/'field_test_singular.stdout').read_text()
    assert 'Drestriction.m2:874' in (failures/'root_patch_test.stdout').read_text()
    assert 'must be PRIME' in (failures/'nonprime.stdout').read_text()
    cli=root/'cli'
    assert G.verify_result(json.loads((cli/'generic_split_cubic.json').read_text()))
    assert json.loads((cli/'family_split_safe_cubic.json').read_text())['presented_ring_normalization_cutoff']['slope']=='1/3'
    for name in ['gaussian_generic_cutoff_producer.py','verify_gaussian_generic_cutoff_producer.py']:ast.parse((p/name).read_text())
    files={str(f.relative_to(p)):sha(f) for f in root.rglob('*') if f.is_file()}
    receipt=dict(completed=True,canonical_proof_sha256=sha(p/'GENERIC_FIELD_CERTIFICATES_AND_NUMERICAL_GAUSSIAN_FAMILY_BOUNDS.md'),
        generator_sha256=sha(p/'gaussian_generic_cutoff_producer.py'),core_sha256=sha(G.CORE),checker_sha256=sha(Path(__file__)),
        runtime_sha256=sha(p/'gaussian_generic_runtime.json'),generic_controls=controls,
        family_controls=[dict(name=n,cutoff=d['presented_ring_normalization_cutoff'],ordinary_torsion=d['verified_ordinary_torsion']) for n,d in families.items()],
        rational_content_controls=roots,rejection_controls=rejects,nonprime_CAS_rejection_retained=True,CLI_routes=2,
        retained_example_files_sha256=files,capped_cusp_restarted=False,
        numerical_cutoffs_sufficient_not_least=True,ordinary_full_cutoff_ideal_producer_implemented=False,
        independent_replay=False,formal_verification=False,external_review=False,
        assurance_boundary='Finite original identities and selected jets replay over supplied prime domains; primality and full joint family decomposition/normalization retain trusted CAS boundaries.')
    (p/'gaussian_generic_cutoff_producer_verification.json').write_text(dump(receipt))
    print(json.dumps(dict(completed=True,generic=len(controls),families=len(families),roots=len(roots),rejections=len(rejects),retained_files=len(files))))


if __name__=='__main__':main()
