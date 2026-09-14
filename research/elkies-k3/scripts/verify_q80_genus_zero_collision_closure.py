#!/usr/bin/env python3
"""Independent finite point census and Hensel matrices for Q80 genus0 k0."""
import argparse
from collections import defaultdict
from fractions import Fraction
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import resource
import time

ROOT=Path(__file__).resolve().parents[2]
DEFAULT=ROOT/'artifacts/generated-results/elkies-k3-q80-genus-zero-integral-reduction-v1'
SOURCE='artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json'
CENSUS='artifacts/generated-results/elkies-k3-q80-alternate-rootless-bisection-orbits.json'
ALT='artifacts/generated-results/q80-alternate-fifth-q6-rootless-transport.json'
FROB='artifacts/generated-results/elliptic-curves/rank_jump_residual_quartic_at_131_v1.json'
FROB_REPLAY='artifacts/generated-results/elliptic-curves/rank_jump_residual_quartic_at_131_verification_v1.json'
PRIOR='artifacts/generated-results/elkies-k3-q80-degree-two-reciprocity-v1'
PRODUCER='elkies-k3/scripts/certify_q80_genus_zero_collision_closure.sage'
HELPERS=['elkies-k3/scripts/verify_q80_degree_two_reciprocity.py',
 'elkies-k3/scripts/verify_q80_single_branch_reciprocity.py',
 'elkies-k3/scripts/verify_q80_complete_genus_one_pairs.py',
 'elkies-k3/scripts/verify_r17_mestre_shared_twist.py',
 'elkies-k3/scripts/verify_r17_correlated_genus_one_pencils.py']
def load_module(name,rel):
    spec=importlib.util.spec_from_file_location(name,ROOT/rel)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod
V=load_module('q80_collision_polynomials',HELPERS[2])
P=131
require=V.require
mul=lambda a,b:V.F.ff_mul(a,b,P)
sub=lambda a,b:V.F.ff_sub(a,b,P)
scale=lambda a,n:V.F.trim([v*n%P for v in a])
add=lambda a,b:sub(a,scale(b,-1))
ev=lambda a,t:V.ev(a,t,P)
def read(path):return json.loads(path.read_text())
def digest(path):return sha256(path.read_bytes()).hexdigest()
def derivative(a):return V.F.trim([i*a[i]%P for i in range(1,len(a))])
def det_mod(rows,p=P):
    a=[[c%p for c in r] for r in rows];n=len(a);answer=1
    require(all(len(r)==n for r in a),'square matrix required')
    for j in range(n):
        pivot=next((i for i in range(j,n) if a[i][j]),None)
        if pivot is None:return 0
        if pivot!=j:a[j],a[pivot]=a[pivot],a[j];answer=-answer
        value=a[j][j];answer=answer*value%p;inv=pow(value,-1,p)
        for i in range(j+1,n):
            factor=a[i][j]*inv%p
            for k in range(j,n):a[i][k]=(a[i][k]-factor*a[j][k])%p
    return answer%p
def det_integer(rows):
    a=[[Fraction(c) for c in r] for r in rows];out=Fraction(1);n=len(a)
    for j in range(n):
        pivot=next((i for i in range(j,n) if a[i][j]),None)
        if pivot is None:return 0
        if pivot!=j:a[j],a[pivot]=a[pivot],a[j];out=-out
        out*=a[j][j]
        for i in range(j+1,n):
            factor=a[i][j]/a[j][j]
            for k in range(j,n):a[i][k]-=factor*a[j][k]
    require(out.denominator==1,'integral matrix determinant');return out.numerator
def congruence(C,H):
    n=len(H)
    return [[sum(C[i][a]*H[a][b]*C[j][b] for a in range(n) for b in range(n)) for j in range(n)] for i in range(n)]
def bit_rank(codes):
    pivots={}
    for c in codes:
        while c:
            i=c.bit_length()-1
            if i in pivots:c^=pivots[i]
            else:pivots[i]=c;break
    return len(pivots)
def hensel_determinant(a,y,b):
    u=[-b%P,1];s=V.divexact(y,u,P)
    left=scale(mul(mul(u,u),s),-2);right=scale(mul(s,s),-1)
    columns=[[0]*j+a for j in range(5)]+[[0]*j+left for j in range(6)]+[[0]*j+right for j in range(2)]
    require(all(len(c)<=13 for c in columns),'degree12 linearization')
    return det_mod([[c[i] if i<len(c) else 0 for c in columns] for i in range(13)])
def point_identity(x,y,A,B):
    require(len(x)<=5 and 0<len(y)<=7,'polynomial height4 degree bounds')
    require(mul(y,y)==add(add(mul(mul(x,x),x),mul(A,x)),B),'finite section equation')
def repeated_labels(rows):
    fibres=defaultdict(set)
    for row in rows:fibres[row['b']].add(row['root'])
    require(all(len(v)==1 for v in fibres.values()),'two repeated contacts select distinct torsion roots')
    return {str(b):sorted(v) for b,v in sorted(fibres.items())}
def verify(out):
    started=time.monotonic();packet=read(out/'input.json');result=read(out/'result.json')
    expected={SOURCE,CENSUS,ALT,FROB,FROB_REPLAY,PRODUCER,str(Path(__file__).relative_to(ROOT)),*HELPERS}
    expected.update(PRIOR+'/'+s for s in ['input.json','result.json','rational-fibres.json','independent-replay.json'])
    require(set(packet['bindings'])==expected,'exact source scope')
    for s,h in packet['bindings'].items():require(digest(ROOT/s)==h,'source binding '+s)
    require(packet['preflight_sha256']==digest(out/'preflight-norm4.json'),'preserved discovery')
    require(packet['prime']==P and packet['unoriented_norm4_count']==1313 and packet['cpu_seconds']==40 and packet['memory_bytes']==4*1024**3,'frozen scope')
    require(result['input_sha256']==digest(out/'input.json') and result['status']=='PASS','result binding')
    require(set(result['records'])=={'norm4-sections.json','finite-lattice.json'},'complete records')
    for s,h in result['records'].items():require(digest(out/s)==h,'record binding')
    D=load_module('q80_collision_degree_two',HELPERS[0]);prior=D.verify(ROOT/PRIOR)
    require(packet['collision_sites']==prior['allowed_double_branch_residues'],'prior branch-collision theorem')
    source=read(ROOT/SOURCE);G=source['sections']['height_gram'];C=source['sections']['coordinate_matrix_in_compiled_frame']
    H=source['frame_certificate']['frame_gram'];U=source['frame_certificate']['integral_isometry_to_alternate_Q80']
    alt=read(ROOT/ALT);census=read(ROOT/CENSUS)
    require(abs(det_integer(C))==abs(det_integer(U))==1 and congruence(C,H)==G and congruence(U,alt['rootless_frame'])==H,'unimodular census transport')
    require(det_integer(G)==948,'Q80 Gram determinant')
    require(census['input']['rootless_frame_sha256']==alt['rootless_frame_sha256'] and census['streaming_enumeration']['signed_shell_counts']['4']==2626,'retained complete norm4 census')
    words=packet['words'];require(len(words)==len({tuple(w) for w in words})==1313,'complete distinct words')
    for w in words:
        require(len(w)==17 and all(type(c) is int for c in w),'integral word')
        require(next(c for c in w if c)>0 and sum(w[i]*G[i][j]*w[j] for i in range(17) for j in range(17))==4,'unoriented norm4 word')
    context_packet={'A':source['weierstrass_model']['A_coefficients_low_to_high'],
        'B':source['weierstrass_model']['B_coefficients_low_to_high'],
        'basis':[{'X':r['X'],'Y':r['Y']} for r in source['sections']['records']],
        'gram':G,'words':words}
    ctx=V.context(context_packet,P) # Independent rational-function basis Gram.
    A,B=ctx['A'],ctx['B']
    codes=[v for r in read(ROOT/PRIOR/'rational-fibres.json')['rows'] for v in r['rational_codes']]
    require(bit_rank(codes)==17,'finite-basis two-saturation')
    old,replay=read(ROOT/FROB),read(ROOT/FROB_REPLAY)
    for data in [old,replay]:
        require(data['status']=='PASS','retained Frobenius status')
        for s,h in data['bindings'].items():require(digest(ROOT/s)==h,'retained Frobenius binding')
    require(old['Frobenius_traces']==[1884,319520],'inherited Frobenius moments')
    t1,t2=1884-17*P,319520-17*P*P
    require((t1,t2)==(-343,27783) and (t1+P)**2==t2+P*P and t1+P==-212,'residual orthogonal factor')
    residual=[P**3,P*P+212*P,212+P,1]
    require(residual==[2248091,44933,343,1] and 2*P*P-212*P!=0,'simple negative-p eigenspace')
    lattice=read(out/'finite-lattice.json')
    require(lattice=={'height_gram':G,'determinant':948,'kummer_character_rank':17,
        'inherited_norm4_signed_count':2626,'residual_H2_eigenvalue_polynomial_coefficients':residual,
        'constant_nonsquare_twist_rank_upper_bound':1},'finite lattice certificate')
    table=read(out/'norm4-sections.json');rows=table['records']
    require(table['prime']==P and [r['index'] for r in rows]==list(range(1313)),'section roster coverage')
    require(len({tuple(r['x']) for r in rows})==1313,'duplicate section up to sign')
    repeated=[];incidences=0;simple=0;matches=0
    for r in rows:
        i=r['index'];require(r['word']==words[i],'word attachment index');x,y=r['x'],r['y']
        require(all(type(c) is int and 0<=c<P for c in x+y),'finite coefficient encoding')
        require(x==V.F.trim(x) and y==V.F.trim(y),'canonical polynomial encoding')
        point_identity(x,y,A,B);a=add(scale(mul(x,x),3),A)
        require(len(V.F.ff_gcd(a,y,P))==1 and (len(y)==7 or len(a)==9),'homogeneous derivative/ordinate coprimality')
        # Two height4 points can meet at no more than6 fibres unless equal.
        # Seven actual fibrewise group sums therefore certify each word.
        for t,sum_word in zip(ctx['sites'][:7],ctx['word_sums'][:7]):
            require((ev(x,t),ev(y,t))==sum_word(i),'seven-fibre norm4 word identification')
            matches+=1
        roots=[]
        for b in packet['collision_sites']:
            if ev(y,b):continue
            order=1;dy=derivative(y)
            while not ev(dy,b):order+=1;dy=derivative(dy);require(dy,'zero ordinate polynomial')
            det=hensel_determinant(a,y,b)
            require(bool(det)==(order==1),'Hensel rank and contact order')
            row={'b':b,'root':ev(x,b),'order':order,'jacobian_determinant':det}
            roots.append(row);incidences+=1
            if order>1:repeated.append({'index':i,**row})
            else:simple+=1
        require(roots==r['collision_roots'],'complete contact roster')
    require((incidences,simple,len(repeated),matches)==(443,439,4,9191),'complete incidence counts')
    require(result['repeated_root_records']==repeated and result['repeated_roots_by_fibre']==repeated_labels(repeated),'repeated-root attachment')
    for key,value in {'unoriented_norm4_sections':1313,'collision_branch_incidences':443,
        'simple_root_hensel_exclusions':439,'repeated_root_pairs_selecting_distinct_torsion':0}.items():require(result[key]==value,'result count '+key)
    require(result['genus0_k0_all_five_allocations']=='EXCLUDED over Q_131 by the written proof' and result['genus1_allocations']=='UNKNOWN','theorem scope')
    require(result['positive_mw17_target_complete'] is False and result['old_norm8_singular_replay_upgraded'] is False,'unproved endpoint or assurance')
    return {'status':'PASS','input_sha256':digest(out/'input.json'),'result_sha256':digest(out/'result.json'),
        'basis_gram_replayed':True,'finite_kummer_character_rank':17,'norm4_sections':1313,
        'exact_fibrewise_word_identifications':9191,'contact_matrices':443,'nonsingular_contact_matrices':439,
        'repeated_contacts':4,'repeated_roots_by_fibre':repeated_labels(repeated),
        'independent_finite_arithmetic_replay':True,'written_integrality_and_hensel_proof_formally_verified':False,
        'positive_mw17_target_complete':False,'old_norm8_singular_replay_upgraded':False,
        'elapsed_seconds':round(time.monotonic()-started,6)}
if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--input',type=Path,default=DEFAULT)
    parser.add_argument('--record',type=Path);args=parser.parse_args()
    resource.setrlimit(resource.RLIMIT_CPU,(40,45));resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3))
    result=verify(args.input)
    if args.record:
        with args.record.open('x') as f:json.dump(result,f,indent=2,sort_keys=True);f.write('\n')
    print(json.dumps(result,sort_keys=True),flush=True)
