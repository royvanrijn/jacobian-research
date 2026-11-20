#!/usr/bin/env python3
"""Complete native finite principal-atom bridge and exact analytic constants.

The native loop checks every required prime/target below2^30. This is not
an extrapolation from the small-anchor control range. Requires g++/C++17.
"""
if not __debug__:
    raise SystemExit('Run without -O/-OO: exact assertions are required.')
import argparse
from fractions import Fraction as F
import hashlib
import itertools
import json
from pathlib import Path
import shutil
import subprocess
import tempfile

HERE=Path(__file__).resolve().parent
parser=argparse.ArgumentParser()
parser.add_argument('--repo-root',type=Path,default=None)
parser.add_argument('--output',type=Path,default=HERE/'minus23_principal_atoms_and_core_five_verification.json')
args=parser.parse_args()
if args.repo_root is not None:
    root=args.repo_root.resolve()
else:
    candidates=[HERE,*HERE.parents,Path.cwd(),*Path.cwd().parents]
    root=next((v for v in candidates if(v/'research/MATH_STATUS.json').is_file()),None)
    if root is None:
        parser.error('Could not locate the repository; pass --repo-root explicitly.')
proof=HERE/'MINUS23_PRINCIPAL_ATOMS_AND_INERT_CORE_FIVE.md'
cpp=HERE/'verify_minus23_principal_atoms_finite.cpp'

def digest(p):
    h=hashlib.sha256()
    with p.open('rb')as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()

ledger_path=HERE/'minus23_principal_atoms_and_core_five_source_ledger.json'
ledger=json.loads(ledger_path.read_text())
for item in ledger['files']:
    target=(root if item['location']=='repository'else HERE)/item['path']
    assert digest(target)==item['sha256'],('stale packet or primary source',str(target))
mandatory={('packet',proof.name),('packet',Path(__file__).name),('packet',cpp.name),
           ('repository','research/arithmetic/sources/rosser_schoenfeld_1962.pdf')}
assert mandatory<={(v['location'],v['path'])for v in ledger['files']}

# The actual selected other-block ideal count, with orientations capped at2,3.
split_small=[2,3,13,29,31,41,47,59]
for p in split_small:
    assert p==2 or pow((-23)%p,(p-1)//2,p)==1
C6=F(2**48,2*3*13*29*31*41*47*59)
assert C6==F(140737488355328,3986190273)<6**6
assert F(3,2)**6<13
# Exact logarithm and onset comparisons, without numerical real arithmetic.
assert F(87,32)**73<1536**10
assert F(3,2)-F(1,1024)>=F(149,100)
margin=F(68,125)*32-(F(64)+F(8000,149))/F(73,10)-F(1152,1024)
assert margin==F(1750191,10877000)>0
assert F(8,3)-1==F(5,3)

# Exact finite neutral phase law, including all outside-bank active counts.
missing=[]
def allowed(e,R,shift):return R+sum((a+b)%3!=0 for a,b in zip(e,shift))!=1
seed_vectors={1:(0,0),2:(1,0),3:(0,1),4:(2,0)}
for e in itertools.product(range(3),repeat=2):
    for R in (0,1,2):
        covered=any(all(allowed(e,R,seed_vectors[v])for v in pair)for pair in [(1,4),(2,3)])
        if not covered:missing.append([*e,R])
assert missing==[[0,0,0],[0,1,0],[1,0,0],[1,2,0]]
def Q(x,y):return x*x+x*y+6*y*y
assert Q(2,1)+Q(8,1)==90
smallP={Q(x,y)for x in range(-20,21)for y in range(-3,4)if __import__('math').gcd(x,y)==1 and 1<=Q(x,y)<=15}
assert smallP=={1,6,8,12}
assert all(not any(a+b==n for a in smallP for b in smallP)for n in [5,10,15])

compiler=shutil.which('g++')
assert compiler,'g++ is needed for the exact native finite bridge'
with tempfile.TemporaryDirectory(prefix='minus23-principal-exact-')as td:
    work=Path(td);exe=work/('checker.exe'if __import__('os').name=='nt'else'checker')
    subprocess.run([compiler,'-O3','-std=c++17','-Wall','-Wextra',str(cpp),'-o',str(exe)],check=True)
    prefix=work/'full'
    subprocess.run([str(exe),str(2**30-1),str(prefix),'5000'],check=True)
    count_path=Path(str(prefix)+'_counts.json')
    counts=json.loads(count_path.read_text());counts.pop('elapsed_seconds')
    assert counts=={
        'inclusive_norm_bound':2**30-1,'prime_cap':214748364,
        'all_primes_through_cap':11849083,'principal_split_primes_through_cap':1974507,
        'primitive_norm_value_count':132325868,'positive_y_ellipse_candidates':703340510,
        'primitive_positive_y_points':427599363,'anchor_limit':5000,'anchor_count':834,
        'coefficient_domain':[5,10,15],'covered_by_coefficient':[1974507,1026431,700565],
        'maximum_anchor_trials':[76,63,67],'independent_square_regression_bound':100000,'failures':[]}
    flag_path=Path(str(prefix)+'_primitive_flags_le64.bin')
    assignment_path=Path(str(prefix)+'_assignments.bin')
    assert flag_path.stat().st_size==2**27
    assert assignment_path.stat().st_size==12*sum(counts['covered_by_coefficient'])
    flag_sha=digest(flag_path);assignment_sha=digest(assignment_path)
    assert flag_sha=='52539112b7b448b497c968da7e79910a6204d370aad1acd3de1f79c41ba73d07'
    assert assignment_sha=='9345ce971a4aad2dcad25e12ad23f2d7242454d48fe5bf8deeb427d5ec22e731'
    # Actual comparison to the old complete coordinate-vs-genus authenticated prefix.
    old=json.loads((root/'research/arithmetic/primitive_norm_minus23_sumset_checks.json').read_text())
    n=old['bound'];assert n==2**24
    with flag_path.open('rb')as f:raw=f.read(n//8+1)
    expansion=[bytes((a>>j)&1 for j in range(8))for a in range(256)]
    byteflags=b''.join(expansion[a]for a in raw)[:n+1]
    assert hashlib.sha256(byteflags).hexdigest()==old['primitive_flags_sha256']==old['ideal_flags_sha256']
    assert sum(byteflags)==old['primitive_value_count']
    # Negative control: keeping only anchor1 really loses required principal cases.
    negative=work/'missing-anchor'
    result=subprocess.run([str(exe),'1000',str(negative),'1'],capture_output=True,text=True)
    assert result.returncode==1
    assert json.loads(Path(str(negative)+'_counts.json').read_text())['failures']

sources=[
 'research/arithmetic/MINUS23_NEUTRAL_CUBE_ATOM_CLOSURE.md',
 'research/arithmetic/PRIMITIVE_NORM_CLASS_NUMBER_THREE.md',
 'research/arithmetic/PRIMITIVE_NORM_MINUS23_EFFECTIVE.md',
 'research/arithmetic/primitive_norm_minus23_sumset_checks.json',
 'research/arithmetic/sources/rosser_schoenfeld_1962.pdf',
]
payload={'status':'PASS','scope':'Actual complete native principal-atom finite bridge below2^30 plus exact improved analytic constants and complete neutral phase reduction. Combined with the written all-prime sieve and retained cube theorem, this proves the entire inert-core5 statement; other cores/global40 remain open.',
 'complete_native_counts':counts,'principal_target_onset_power_of_two':30,'scaled_positive_surplus':str(margin),
 'selected_other_binary_representation_constant':12,'selected_ideal_constant_sixth_power':str(C6),
 'norm_flags_encoding':'little-endian uint64 bitset, bitn means n is a primitiveQ value',
 'assignment_encoding':'little-endian uint32 triples(p,c,anchor), in increasingp and coefficient5/10/15 order',
 'full_primitive_flags_sha256':flag_sha,'full_assignment_sha256':assignment_sha,
 'old_complete_coordinate_vs_genus_prefix_bound':2**24,'old_prefix_primitive_flags_sha256':old['primitive_flags_sha256'],
 'neutral_phase_missing_cells':missing,'inert_core_five_exact_exceptions':[5,10,15,115,230,345],
 'negative_controls':['ordinary_norm4_is_not_primitive','independent_signed_square_through100000','missing_anchor_family_has_actual_failures'],
 'source_sha256':{'proof':digest(proof),'python_checker':digest(Path(__file__).resolve()),'native_checker':digest(cpp),'source_ledger':digest(ledger_path)},
 'retained_source_sha256':{s:digest(root/s)for s in sources},
 'compile_command':['g++','-O3','-std=c++17','-Wall','-Wextra'],'native_runtime_is_not_a_proof_premise':True,
 'external_review_or_formal_verification':False}
args.output.parent.mkdir(parents=True,exist_ok=True)
args.output.write_text(json.dumps(payload,indent=2)+'\n')
print(json.dumps(payload,indent=2))
