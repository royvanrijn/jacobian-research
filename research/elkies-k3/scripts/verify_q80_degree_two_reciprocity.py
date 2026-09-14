#!/usr/bin/env python3
"""Independent F_131[i] arithmetic for the complete degree-two reduction gate."""
import argparse
from collections import Counter, defaultdict
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import resource
import time

ROOT=Path(__file__).resolve().parents[2]
DEFAULT=ROOT/'artifacts/generated-results/elkies-k3-q80-degree-two-reciprocity-v1'
SOURCE='artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json'
PRODUCER='elkies-k3/scripts/certify_q80_degree_two_reciprocity.sage'
OLD_CHECKER='elkies-k3/scripts/verify_q80_single_branch_reciprocity.py'
OLD_INPUT='artifacts/generated-results/elkies-k3-q80-single-branch-reciprocity-v1/input.json'
P=131
ORDER=P*P
def require(v,message):
    if not v:raise ValueError(message)
def read(path):return json.loads(path.read_text())
def digest(path):return sha256(path.read_bytes()).hexdigest()
def add(a,b):return (a%P+b%P)%P+P*((a//P+b//P)%P)
def sub(a,b):return (a%P-b%P)%P+P*((a//P-b//P)%P)
def mul(a,b):
    ar,ai=a%P,a//P;br,bi=b%P,b//P
    return (ar*br-ai*bi)%P+P*((ar*bi+ai*br)%P)
def inv(a):
    ar,ai=a%P,a//P
    n=(ar*ar+ai*ai)%P;require(n,'zero field denominator')
    d=pow(n,-1,P)
    return ar*d%P+P*((-ai*d)%P)
def square(a):
    ar,ai=a%P,a//P
    return pow((ar*ar+ai*ai)%P,(P-1)//2,P)==1
def trim(a):
    while a and not a[-1]:a.pop()
    return a
def poly_mul(a,b):
    c=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):c[i+j]=add(c[i+j],mul(x,y))
    return trim(c)
def poly_sub(a,b):
    return trim([sub(a[i] if i<len(a) else 0,b[i] if i<len(b) else 0) for i in range(max(len(a),len(b)))])
def poly_rem(a,b):
    a=trim(a[:]);require(b,'zero polynomial modulus');lead=inv(b[-1])
    while len(a)>=len(b):
        v=mul(a[-1],lead);shift=len(a)-len(b)
        for i,c in enumerate(b):a[shift+i]=sub(a[shift+i],mul(v,c))
        trim(a)
    return a
def poly_pow(a,n,f):
    out=[1]
    while n:
        if n&1:out=poly_rem(poly_mul(out,a),f)
        n>>=1
        if n:a=poly_rem(poly_mul(a,a),f)
    return out
def poly_gcd(a,b):
    while b:a,b=b,poly_rem(a,b)
    return a
def ev(a,t):
    value=0
    for c in reversed(a):value=add(mul(value,t),c)
    return value
def check_roster(a,b,roots):
    require(roots==sorted(roots) and all(type(e) is int and 0<=e<ORDER for e in roots),'field root encoding')
    f=[b,a,0,1]
    if len(roots)==3:
        product=[1]
        for e in roots:product=poly_mul(product,[sub(0,e),1])
        require(product==f,'complete three-root factorization')
    elif len(roots)==1:
        e=roots[0];require(ev(f,e)==0,'claimed cubic root')
        disc=sub(0,add(mul(3,mul(e,e)),mul(4,a)))
        require(disc and not square(disc),'omitted quadratic roots')
    elif not roots:
        fixed=poly_sub(poly_pow([0,1],ORDER,f),[0,1])
        require(len(poly_gcd(f,fixed))==1,'omitted cubic roots')
    else:raise ValueError('incomplete multiplicity roster')
def derivative(a):return trim([mul(i%P,a[i]) for i in range(1,len(a))])
def leading_at(a,t):
    a=trim(a[:]);order=0
    while a:
        value=ev(a,t)
        if value:return order,value
        order+=1;a=derivative(a)
    return float('inf'),0
def rational_at(n,d,t,weight):
    n,d=trim(n[:]),trim(d[:]);require(d,'zero reduced denominator')
    if t is None:
        gap=len(d)-1+weight-(len(n)-1)
        if gap<0:return None
        return 0 if gap else mul(n[-1],inv(d[-1]))
    vn,an=leading_at(n,t);vd,ad=leading_at(d,t)
    if vn<vd:return None
    if vn>vd:return 0
    return mul(an,inv(ad))
def character_value(a,e,x):
    value=1 if x is None else (add(mul(3,mul(e,e)),a) if x==e else sub(x,e))
    require(value,'section meets singular node')
    return value
def characters(a,b,roots,points,degree):
    for x,y in points:
        require((x is None)==(y is None),'coordinate pole mismatch')
        require(x is None or mul(y,y)==add(add(mul(x,mul(x,x)),mul(a,x)),b),'specialized section identity')
    codes=[]
    for e in roots:
        code=0
        for i,(x,y) in enumerate(points):
            value=character_value(a,e,x)
            if degree==1:
                require(value<P,'nonrational value in rational character')
                is_square=pow(value,(P-1)//2,P)==1
            else:is_square=square(value)
            if not is_square:code|=1<<i
        codes.append(code)
    return codes
def verify(out):
    start=time.monotonic();packet=read(out/'input.json');result=read(out/'result.json')
    expected={SOURCE,PRODUCER,OLD_CHECKER,OLD_INPUT,str(Path(__file__).relative_to(ROOT))}
    require(set(packet['bindings'])==expected,'exact source scope')
    for s,h in packet['bindings'].items():require(digest(ROOT/s)==h,'frozen source binding')
    require(packet['source']==SOURCE and packet['prime']==P and packet['field_modulus_low_to_high']==[1,0,1],'fixed parent and field')
    require(packet['section_indices']==list(range(17)) and packet['cpu_seconds']==60 and packet['memory_bytes']==4*1024**3,'frozen scope and limits')
    require(pow(P-1,(P-1)//2,P)==P-1,'quadratic field modulus')
    require(result['input_sha256']==digest(out/'input.json') and result['status']=='PASS','result binding')
    require(set(result['records'])=={'rational-fibres.json','quadratic-fibres.json'},'record scope')
    for s,h in result['records'].items():require(digest(out/s)==h,'record binding')
    spec=importlib.util.spec_from_file_location('q80_reciprocity_prior',ROOT/OLD_CHECKER)
    old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
    prior=old.verify((ROOT/OLD_INPUT).parent)
    source=read(ROOT/SOURCE);model=source['weierstrass_model']
    A0=list(map(old.F,model['A_coefficients_low_to_high']));B0=list(map(old.F,model['B_coefficients_low_to_high']))
    A,B=[[old.mod(c,P) for c in values] for values in [A0,B0]]
    sections=[old.section_data(s,A0,B0,P) for s in source['sections']['records']]
    rational=read(out/'rational-fibres.json')['rows'];quadratic=read(out/'quadratic-fibres.json')['rows']
    require([r['t'] for r in rational]==list(range(P))+[None],'complete projective rational base')
    require([r['t'] for r in quadratic]==[i+P*j for i in range(P) for j in range(1,66)],'complete quadratic base orbit representatives')
    for degree,rows in [(1,rational),(2,quadratic)]:
        for r in rows:
            t=r['t'];a,b=(A[8],B[12]) if t is None else (ev(A,t),ev(B,t))
            delta=add(mul(4,mul(a,mul(a,a))),mul(27,mul(b,b)))
            require(r['smooth'] is bool(delta),'fibre smoothness attachment')
            check_roster(a,b,r['roots'])
            points=[(rational_at(xn,xd,t,4),rational_at(yn,yd,t,6)) for xn,xd,yn,yd in sections]
            codes=characters(a,b,r['roots'],points,2)
            require(codes==r['norm_codes'],'all seventeen norm characters')
            if len(codes)==3:require(codes[0]^codes[1]^codes[2]==0,'fibre norm relation')
            if degree==1:
                require(delta,'rational singular fibre')
                rr=[e for e in r['roots'] if e<P]
                require(rr==r['rational_roots'],'rational root descent')
                require(characters(a,b,rr,points,1)==r['rational_codes'],'rational characters')
                require(all(r['rational_codes']),'single-branch obstruction')
    groups=defaultdict(list);split=[]
    for r in rational:
        if len(r['rational_roots'])==3:
            split.append(r['t'])
            for i in range(3):
                for j in range(3):
                    if i!=j:groups[(r['rational_codes'][i],r['rational_codes'][j])].append([r['t'],r['roots'][i],r['roots'][j]])
    require(len(groups)==96 and all(len(v)==1 for v in groups.values()),'rational ordered-pair code injectivity')
    require(all(r['norm_codes'].count(0)<2 for r in quadratic if len(r['roots'])==3),'quadratic zero-code pair')
    require(all(r['norm_codes'].count(0)<2 for r in rational if len(r['roots'])==3 and len(r['rational_roots'])==1),'quadratic roots over rational base')
    hist=Counter((len(r['roots']),r['smooth']) for r in quadratic)
    require(result['quadratic_histogram']==[[n,s,v] for (n,s),v in sorted(hist.items())],'quadratic census')
    require(result['bad_quadratic_fibres']==[r for r in quadratic if not r['smooth']],'nodal reduction coverage')
    require(result['rational_full_splitting_values']==split and len(split)==16,'allowed collision fibres')
    for key,value in {'rational_base_values':132,'quadratic_base_orbits':8515,'rational_ordered_pair_codes':96,
                      'rational_ordered_pair_collisions':0,'quadratic_two_zero_code_triples':0,
                      'rational_nonsplit_two_zero_code_triples':0}.items():require(result[key]==value,'result count '+key)
    require(result['genus0_branch_reduction']=='double finite point in the sixteen full-splitting residues','theorem scope')
    require(all(result[k] is False for k in ['collision_strata_excluded','coefficient_lift_complete','positive_mw17_target_complete']),'unproved endpoint')
    return {'status':'PASS','input_sha256':digest(out/'input.json'),'result_sha256':digest(out/'result.json'),
        'projective_rational_base_values':132,'quadratic_base_orbits':8515,'inherited_section_identities':17,
        'rational_ordered_pair_codes':96,'allowed_double_branch_residues':split,
        'prior_independent_replay_status':prior['status'],'independent_extension_field_replay':True,
        'written_reciprocity_and_specialization_formally_verified':False,
        'collision_strata_excluded':False,'positive_mw17_target_complete':False,
        'elapsed_seconds':round(time.monotonic()-start,6)}
if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--input',type=Path,default=DEFAULT)
    parser.add_argument('--record',type=Path);args=parser.parse_args()
    resource.setrlimit(resource.RLIMIT_CPU,(60,65));resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3))
    result=verify(args.input)
    if args.record:
        with args.record.open('x') as f:json.dump(result,f,indent=2,sort_keys=True);f.write('\n')
    print(json.dumps(result,sort_keys=True),flush=True)
