#!/usr/bin/env python3
"""Replay or audit exact Singular geometry, with a Hilbert completeness check."""
from pathlib import Path
from math import comb
import argparse,hashlib,json,re,subprocess,sys,time
import sympy as S

P=Path(__file__).resolve().parent
EXPECTED_WICK='7254eff7dad198b2e15fa727da0703bab5f8e3b8704697b27a0ed212527821d6'
assert hashlib.sha256((P/'verify_nodal_product_reduction.py').read_bytes()).hexdigest()==EXPECTED_WICK
from verify_nodal_product_reduction import moments


def parse_polynomial(raw,variables):
    raw=raw.strip().strip(',').replace(' ','')
    if raw=='0':return S.Integer(0)
    answer=0
    for token in re.findall(r'[+-]?[^+-]+',raw):
        sign=-1 if token.startswith('-') else 1;token=token.lstrip('+-')
        number=re.match(r'\d+(?:/\d+)?',token)
        coefficient=S.Rational(number[0]) if number else S.Integer(1)
        rest=token[len(number[0]):] if number else token;monomial=S.Integer(1)
        while rest:
            part=re.match(r'([a-f])(\d*)',rest);assert part,(raw,rest)
            monomial*=variables[part[1]]**(int(part[2]) if part[2] else 1)
            rest=rest[len(part[0]):]
        answer+=sign*coefficient*monomial
    return S.expand(answer)


def main():
    start=time.monotonic();parser=argparse.ArgumentParser();parser.add_argument('--audit-existing',action='store_true');args=parser.parse_args()
    if not args.audit_existing:
        run=subprocess.run([sys.executable,str(P/'run_three_factor_colon_geometry_v3.py')],text=True,capture_output=True)
        assert run.returncode==0,(run.stdout,run.stderr)
    rec=json.loads((P/'three_factor_colon_geometry_v3_run.json').read_text())
    assert rec['terminal'] and rec['completed'] and rec['exit_code']==0 and not rec['timed_out']
    input_path=P/'three_factor_colon_geometry_v3.sing';log_path=P/'three_factor_colon_geometry_v3_stdout.txt'
    assert rec['input_sha256']==hashlib.sha256(input_path.read_bytes()).hexdigest()
    assert rec['log_sha256']==hashlib.sha256(log_path.read_bytes()).hexdigest()
    source=input_path.read_text();log=log_path.read_text();assert '?' not in log and log.rstrip().endswith('END')
    a,b,c,d,e,f=S.symbols('a b c d e f');coords=[a,b,c,d,e,f];variables=dict(zip('abcdef',coords))
    mm=moments([[a,d,e],[d,b,f],[e,f,c]],10)
    generated='ideal I='+','.join(str(mm[k-1]).replace('**','^') for k in [2,4,6,8])+';'
    assert generated in source and 'poly P='+str(mm[9]).replace('**','^')+';' in source
    assert 'ideal Iall=I,P;' in source and 'ideal G10=std(Iall);' in source
    assert 'K=std(K*J)' in source and 'list Lsat=sat_with_exp(Q,ideal(a+b+c));' in source
    markers=list(re.finditer(r'(?m)^([A-Z][A-Z0-9_]*)$',log));blocks={}
    for i,m in enumerate(markers):
        body=log[m.end():markers[i+1].start() if i+1<len(markers) else len(log)].strip()
        blocks.setdefault(m[1],[]).append(body)
    def one(name):assert len(blocks[name])==1;return blocks[name][0]
    def integers(name):return [int(x) for x in re.findall(r'-?\d+',one(name))]
    def vector(name):
        v=integers(name)
        while v and v[-1]==0:v.pop()
        return v
    t=S.symbols('t')
    def polynomial(v):return S.Add(*(c*t**i for i,c in enumerate(v)))
    H8=polynomial(vector('HILBERT_FIRST_I8'));H10=polynomial(vector('HILBERT_FIRST_I10'));HQ=polynomial(vector('HILBERT_FIRST_COLON'))
    assert S.expand(H8-H10-t**15*HQ)==0
    assert integers('COLON_FORWARD_MULTIPLIER')==[0]
    hq=vector('COLON_HILBERT');hs=vector('SATURATION_HILBERT')
    assert hq==[1,5,15,32,52,63,47,4,-24,-8,-1]
    assert hs==[1,5,15,32,52,57,23,1]
    assert S.expand(HQ-polynomial(hq)*(1-t)**5)==0
    h8=S.cancel(H8/(1-t)**3);h10=S.cancel(H10/(1-t)**3)
    assert S.denom(h8)==S.denom(h10)==1
    assert h8.subs(t,1)==h10.subs(t,1)==90
    torsion=S.cancel((polynomial(hq)-polynomial(hs))/(1-t));assert S.denom(torsion)==1
    th=[int(S.expand(torsion).coeff(t,k)) for k in range(S.degree(torsion,t)+1)]
    assert th==[0,0,0,0,0,6,30,33,9,1] and sum(th)==79
    assert sum(hq)==sum(hs)==186
    for label in ['COLON_RADICAL_FORWARD','COLON_RADICAL_REVERSE','AXIS_RADICAL_FORWARD','AXIS_RADICAL_REVERSE','AXIS_PARAMETER_REGULAR','SATURATION_EXPONENT_CHECK']:
        assert integers(label)==[0],label
    assert integers('SATURATION_INTERSECTION_CHECK')==[0,0]
    assert integers('PRIMARY_DECOMPOSITION_CHECK')==[0,0]
    assert integers('SATURATION_POWER')==[4] and integers('AXIS_GENERIC_LENGTH')==[62]
    assert integers('ORIGIN_COMPONENT_LENGTH')==[823] and 4*186+79==823
    qgens=[parse_polynomial(v,variables) for v in one('COLON_GENERATORS').split(',') if v.strip()]
    qagens=[parse_polynomial(v,variables) for v in one('AXIS_PRIMARY_BASIS').split(',') if v.strip()]
    assert len(qgens)==55 and len(qagens)==38
    assert all(S.Poly(v,*coords).is_homogeneous for v in qgens+qagens)
    assert min(S.Poly(v,*coords).total_degree() for v in qgens)==3
    cubics=[v for v in qgens if S.Poly(v,*coords).total_degree()==3]
    assert len(cubics)==3 and set(map(S.expand,cubics))=={a*b*c,d*e*f,a*f*f+b*e*e+c*d*d}
    # Exact low-degree annihilator boundary from the complete colon Hilbert function.
    cumulative=[sum(hq[:k+1]) for k in range(4)]
    assert cumulative==[1,6,21,53]
    assert [comb(k+5,5)-cumulative[k] for k in range(4)]==[0,0,0,3]
    powers=[int(v.splitlines()[0]) for v in blocks['NIL_POWER']];assert [int(v.splitlines()[1]) for v in blocks['NIL_POWER']]==[819,1015,1240];r8=[int(v) for v in blocks['I8_REMAINDER_SIZE']];r10=[int(v) for v in blocks['I10_REMAINDER_SIZE']]
    assert powers==[12,13,14] and r8==[33,6,0] and r10==[6,0,0]
    lower8=[parse_polynomial(v,variables) for v in one('I8_LOWER_WITNESS').splitlines()]
    lower10=[parse_polynomial(v,variables) for v in one('I10_LOWER_WITNESS').splitlines()]
    assert lower8==[(a*f)**12*(b*e)]*2 and lower10==[(a*f)**11*(b*e)]*2
    assert integers('GLOBAL_BASES')==[815,644]
    certificate={'field':'QQ, characteristic zero','algorithm':'Exact standard bases, verified colon inclusion plus full Hilbert-series equality, saturation and monomial radical-power reductions.',
                 'colon_generators':list(map(str,qgens)),'axis_a_primary_generators':list(map(str,qagens)),
                 'other_axis_components':'Factor-label permutations of the a-axis component',
                 'origin_component':'Q+((a+b+c)^4)',
                 'colon_radical':'(d,e,f,a*b,a*c,b*c), union of the three diagonal axes',
                 'first_Hilbert_numerators':{'I8':vector('HILBERT_FIRST_I8'),'I10':vector('HILBERT_FIRST_I10'),'Q':vector('HILBERT_FIRST_COLON')},
                 'Hilbert_completeness_identity':'H_I8-H_I10=t^15*H_Q, all over denominator(1-t)^6; inclusion Q*M10 in I8 is also checked.',
                 'reduced_Hilbert_numerators':{'Q':hq,'saturated_Q':hs},'origin_torsion_Hilbert_function':th,
                 'axis_generic_lengths':[62,62,62],'origin_torsion_length':79,'chosen_origin_primary_component_length':823,
                 'low_degree_annihilator_dimensions':[0,0,0,3],
                 'exact_nilradical_indices':{'I8':14,'I10_equals_I_infinity':13},
                 'nonmembership_witnesses':{'I8_J_power13':str(lower8[0]),'I10_J_power12':str(lower10[0])},
                 'provenance':{'input_sha256':rec['input_sha256'],'log_sha256':rec['log_sha256'],
                               'runner_sha256':hashlib.sha256((P/'run_three_factor_colon_geometry_v3.py').read_bytes()).hexdigest(),
                               'wick_input_sha256':EXPECTED_WICK},
                 'boundary':'This is the full annihilator of I10/I8 and its primary decomposition, not a full primary decomposition of I10. The raw syzygy route timed out; the Hilbert identity supplies completeness instead.'}
    artifact=P/'three_factor_colon_geometry_certificate.json';artifact.write_text(json.dumps(certificate,indent=2)+'\n')
    singular_version=subprocess.run(['/usr/bin/Singular','--version'],text=True,capture_output=True).stdout.splitlines()[0]
    receipt={'claim':'G3PRODUCTCOLON20261114','audit_existing_only':args.audit_existing,
             'exact_nilradical_indices':certificate['exact_nilradical_indices'],'colon_standard_basis_generators':55,
             'primary_components':4,'associated_primes_of_I10_over_I8':['(b,c,d,e,f)','(a,c,d,e,f)','(a,b,d,e,f)','(a,b,c,d,e,f)'],
             'axis_generic_length':62,'intrinsic_origin_torsion_length':79,'chosen_origin_primary_length':823,
             'colon_completeness':'Forward ideal membership plus exact full Hilbert-series identity, not merely candidate generators.',
             'certificate_sha256':hashlib.sha256(artifact.read_bytes()).hexdigest(),
             'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
             'input_hashes':certificate['provenance'],'singular_version':singular_version,'sympy_version':S.__version__,
             'elapsed_seconds':time.monotonic()-start,'assurance':'Internal exact symbolic/CAS proof, no independent specialist replay, formal verification or external review. Invalid v1/v2 runs and terminal raw-syzygy timeout are retained.'}
    (P/'three_factor_colon_geometry_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))


if __name__=='__main__':main()
