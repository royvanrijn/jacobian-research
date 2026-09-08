#!/usr/bin/env sage-python
"""Exact local Galois-set descent through pair and theta resolvents.

Small integer proxy polynomials are proved 2-adically root-equivalent to
the input sextics. No global field or class group is constructed.
"""
import argparse,hashlib,json,signal
from itertools import combinations
from pathlib import Path
from sage.all import QQ,ZZ,PolynomialRing,matrix,identity_matrix,floor
ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'det1092_rr_dyadic_theta_v1'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def retain(p,d):
    payload=json.dumps(d,indent=2,sort_keys=True)+'\n'
    if p.exists():assert p.read_text()==payload
    else:p.write_text(payload)
def exterior_sum(A,k):
    basis=list(combinations(range(6),k));lookup={b:i for i,b in enumerate(basis)}
    B=matrix(ZZ,len(basis),len(basis))
    for j,b in enumerate(basis):
        for slot,old in enumerate(b):
            for new in range(6):
                if not A[new,old]:continue
                word=list(b);word[slot]=new
                if len(set(word))<k:continue
                sign=(-1)**sum(word[u]>word[v] for u in range(k) for v in range(u+1,k))
                B[lookup[tuple(sorted(word))],j]+=sign*A[new,old]
    return B
def root_tree(F,max_nodes=4096,max_depth=64):
    R=F.parent();x=R.gen();nodes=[];roots=[];pending=[(ZZ(0),0,F)]
    while pending:
        assert len(nodes)<max_nodes,'NODE_CAP'
        a,depth,P=pending.pop();assert depth<=max_depth,'DEPTH_CAP'
        content=min(c.valuation(2) for c in P if c);P=R(P/2**content)
        assert all(c.denominator()==1 for c in P)
        record={'a':str(a),'depth':depth,'content_v2':int(content),
                'normalized_polynomial_sha256':hashlib.sha256(str(P).encode()).hexdigest(),'children':[]}
        for r in [0,1]:
            value=P(r);derivative=P.derivative()(r)
            if value%2:status='NO_ROOT_MOD_TWO'
            elif derivative%2:
                status='UNIQUE_HENSEL_ROOT';roots.append({'a':str(a+2**depth*r),'depth':depth+1})
            else:
                status='SUBDIVIDE';pending.append((a+2**depth*r,depth+1,R(P(r+2*x))))
            record['children'].append({'residue':r,'status':status})
        nodes.append(record)
    return {'status':'COMPLETE_EXACT_Q2_ROOT_TREE','nodes':nodes,'root_balls':roots,'root_count':len(roots),
            'node_count':len(nodes),'maximum_depth':max(n['depth'] for n in nodes)}
def construct(i):
    paths=[ART/'det1092_rr_good_local_images_v1'/('case-%02d.json'%j) for j in range(10)]
    protocol={'classification':'frozen local arithmetic construction',
              'inputs':{str(p.relative_to(ROOT)):sha(p) for p in paths},
              'limits':{'cases':10,'seconds_per_case':25,'maximum_proxy_precision':512,
                        'root_tree_nodes_per_polynomial':4096,'root_tree_depth':64,
                        'resolvent_degrees':[6,15,10],'point_searches':0,'class_groups':0,
                        'parameter_searches':0,'remote_arithmetic':0,'pilot_changes':0},
              'proxy_rule':'Largest integral monic homothety; coefficientwise least residues modulo2^N with N=max(64,2*v2(disc)+16).',
              'script_sha256':sha(Path(__file__))}
    OUT.mkdir(parents=True,exist_ok=True);retain(OUT/'protocol.json',protocol)
    source=json.loads(paths[i].read_text());R=PolynomialRing(QQ,'x');x=R.gen();q=R(source['q']).monic()
    scale=min(floor(a.valuation(2)/(6-j)) for j,a in enumerate(q.list()[:-1]) if a)
    f=R(q(2**scale*x)/2**(6*scale));assert f.is_monic() and all(a.valuation(2)>=0 for a in f if a)
    vd=int(f.discriminant().valuation(2));N=max(64,2*vd+16);assert N<=512
    modulus=ZZ(2)**N
    lift=lambda a: ZZ(a.numerator()*a.denominator().inverse_mod(modulus)%modulus)
    g=R([lift(a) for a in f]);assert g.is_monic()
    assert all((a-b).valuation(2)>=N for a,b in zip(f,g) if a!=b)
    assert g.discriminant().valuation(2)==vd and N>2*vd
    # Multiplication by x, then its additive action on exterior powers.
    A=matrix(ZZ,6,6)
    for j in range(5):A[j+1,j]=1
    for j in range(6):A[j,5]=-ZZ(g[j])
    pairs=R(exterior_sum(A,2).charpoly().list())
    triples=exterior_sum(A,3);shifted=2*triples-A.trace()*identity_matrix(ZZ,20)
    triple_char=R(shifted.charpoly().list());assert all(triple_char[j]==0 for j in range(1,20,2))
    theta=R([triple_char[2*j] for j in range(11)])
    assert pairs.degree()==15 and theta.degree()==10
    assert pairs.gcd(pairs.derivative())==theta.gcd(theta.derivative())==1
    prefix={'classification':'exact local Galois resolvent construction; independent replay required',
            'case_index':i,'source':str(paths[i].relative_to(ROOT)),'source_sha256':sha(paths[i]),
            'T_equals_two_power_times_x':int(scale),'integral_monic_model':list(map(str,f.list())),
            'proxy_polynomial':list(map(str,g.list())),'coefficient_precision':N,'discriminant_v2':vd,
            'pair_sum_resolvent':list(map(str,pairs.list())),
            'theta_squared_difference_resolvent':list(map(str,theta.list())),
            'limits':protocol['limits']}
    retain(OUT/('case-%02d-resolvents.json'%i),prefix)
    trees={}
    for name,F in [('branch',g),('pair',pairs),('theta',theta)]:
        tree=root_tree(F);retain(OUT/('case-%02d-%s-tree.json'%(i,name)),tree);trees[name]=tree
    torsion_order=1+trees['pair']['root_count'];dim=torsion_order.bit_length()-1
    assert 2**dim==torsion_order
    rational_theta=bool(trees['branch']['root_count'] or trees['theta']['root_count'])
    kernel=int(not rational_theta)
    result={**prefix,'status':'PASS_EXACT_DYADIC_TORSION_AND_INHERITED_THETA_CLASS',
            'local_rational_2_torsion_dimension':dim,
            'local_true_Kummer_dimension':dim+2,'local_fake_Kummer_dimension':dim+2-kernel,
            'inherited_canonical_minus_two_basepoint_is_locally_divisible_by_two':rational_theta,
            'local_fake_kernel_dimension':kernel,
            'rational_branch_points':trees['branch']['root_count'],
            'rational_even_theta_characteristics':trees['theta']['root_count'],
            'tree_summaries':{key:{k:tree[k] for k in ['root_count','node_count','maximum_depth']} for key,tree in trees.items()},
            'full_local_Kummer_image':'NOT_COMPUTED; dimensions and inherited kernel class only',
            'full_global_Selmer_group':'NOT_COMPUTED',
            'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [paths[i],OUT/'protocol.json',Path(__file__)]}}
    retain(OUT/('case-%02d.json'%i),result)
    print('case',i,'true/fake dimensions',dim+2,dim+2-kernel,'D0 local2-divisible',rational_theta,
          'trees',result['tree_summaries'],flush=True)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--case',type=int,required=True,choices=range(10));args=ap.parse_args()
    signal.alarm(25);construct(args.case)
