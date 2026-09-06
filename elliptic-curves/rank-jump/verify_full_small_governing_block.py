#!/usr/bin/env python3
"""Exact norm, order, joint-group, and full-matrix verification; no curve search."""
import argparse
from collections import Counter
from fractions import Fraction as Q
from itertools import product
from math import prod
from pathlib import Path
import retrospective as r
import full_small_governing_block as source
import verify_unpointed_governing_norm as norm
import verify_explicit_governing_octic as octic

OUTPUT = r.OUT/'rank_jump_full_small_governing_block_v1.json'
SMALL = r.OUT/'rank_jump_small_quotient_block_v1.json'
DYADIC = r.OUT/'rank_jump_strict_cochain_dyadic_switch_v1.json'
BASES = [
    [['1'], ['0','1'], ['0','0','1'], ['0','0','0','1'], ['-1/2','0','0','0','1/2'],
     ['1/4','-1/4','-1/2','-1/2','-1/4','1/4'],
     ['39/100','0','1/20','0','3/20','0','1/100'],
     ['0','1089/2900','0','221/580','0','-67/580','0','1/2900']],
    [['1'], ['0','1'], ['0','0','1'], ['0','0','0','1'], ['-1/2','0','0','0','1/2'],
     ['1/4','3/20','-1/2','-1/2','-1/4','1/20'],
     ['-1/4','0','-7/100','0','-1/4','0','1/100'],
     ['1/4','1/10','-1/2','143/500','-1/4','0','0','1/500']]
]


def inverse_and_det(matrix):
    n = len(matrix); a = [[Q(v) for v in row]+[Q(i == j) for j in range(n)] for i, row in enumerate(matrix)]
    det = Q(1)
    for i in range(n):
        j = next(j for j in range(i, n) if a[j][i])
        if j != i: a[i], a[j] = a[j], a[i]; det = -det
        pivot = a[i][i]; det *= pivot; a[i] = [v/pivot for v in a[i]]
        for j in range(n):
            if j != i:
                pivot = a[j][i]; a[j] = [x-pivot*y for x, y in zip(a[j], a[i])]
    return [row[n:] for row in a], det


def arithmetic_pair(row, basis):
    K = norm.Algebra([-1,-14,-11,1])
    a, b, x, y = [K.elt(row[k]) for k in ('alpha_ascending','beta_ascending','X_ascending','Y_ascending')]
    assert K.add(K.mul(x,x), K.neg(K.mul(a,K.mul(y,y)))) == b
    assert K.norm(a) == 625 and K.norm(b) == 1
    delta = K.add(K.mul(x,x), K.neg(b)); D = K.norm(delta); P = K.norm(x)
    ratio = K.mul(K.add(K.mul(x,x),b), K.inverse(delta)); matrix = K.matrix(ratio)
    C = 2*D*sum(matrix[i][i] for i in range(3))
    h = [D*D,0,-4*D*(P-1),0,C,0,-4*(P+1),0,1]
    assert list(map(Q,row['rational_octic_ascending'])) == h
    g = list(map(int,row['reduced_octic_ascending'])); A = norm.Algebra(g)
    change = A.elt(row['old_root_in_reduced_algebra']); assert A.evaluate(h,change) == A.elt([])
    inverse_and_det(list(zip(*(A.power(change,i) for i in range(8)))))
    discriminant = octic.discriminant(g); assert str(discriminant) == row['reduced_octic_discriminant']
    assert prod(int(p)**e for p,e in row['discriminant_factorization']) == discriminant
    vectors = [A.elt(v) for v in basis]
    inverse, determinant = inverse_and_det(list(zip(*vectors)))
    for u,v in product(vectors,repeat=2):
        w = A.mul(u,v)
        coordinates = [sum(a*b for a,b in zip(line,w)) for line in inverse]
        assert all(c.denominator == 1 for c in coordinates)
    assert vectors[0] == A.elt([1])
    order_disc = determinant**2*discriminant
    assert order_disc == 2**18*163**4
    return {'class_indices':row['class_indices'],'norm_identity':True,'degree8_reduction_verified':True,
            'overorder_basis':basis,'integral_multiplication_products':64,
            'order_discriminant':str(order_disc),'finite_ramification_support_contained_in':[2,163]}


def joint_group():
    V=(0,3,5,6); G=((0,1,2),(1,2,0),(2,0,1)); roots=list(product(V,range(2)))
    def act(g,v): return sum(((v>>i)&1)<<g[i] for i in range(3))
    def lam(v): return act(G[2],v)
    def dot(a,b): return (a&b).bit_count()%2
    def compose(g,h): return tuple(g[h[i]] for i in range(3))
    def multiply(x,y):
        e,c,g,z=x; ee,cc,h,w=y
        omega=dot(e,act(g,lam(ee))) | (dot(e,act(g,cc))<<1) | (dot(lam(e),act(g,cc))<<2)
        return e^act(g,ee),c^act(g,cc),compose(g,h),z^w^omega
    group=list(product(V,V,G,range(8))); unit=(0,0,G[0],0)
    lookup={v:i for i,v in enumerate(group)}
    table=[[lookup[multiply(x,y)] for y in group] for x in group]
    identity=lookup[unit]; inverses=[next(j for j in range(384) if table[i][j]==identity and table[j][i]==identity) for i in range(384)]
    def action(x):
        e,c,g,z=x; answer=[]
        for index,(a,b) in enumerate(((e,lam(e)),(e,c),(lam(e),c))):
            answer += [8*index+roots.index((act(g,v)^b,s^((z>>index)&1)^dot(a,act(g,v)))) for v,s in roots]
        return tuple(answer)
    actions=[action(x) for x in group]; assert len(set(actions))==384
    for i,j in product(range(384),repeat=2):
        assert actions[table[i][j]] == tuple(actions[i][actions[j][k]] for k in range(24))
    commutators={table[table[table[i][j]][inverses[i]]][inverses[j]] for i,j in product(range(384),repeat=2)}
    derived={identity}; pending=[identity]
    while pending:
        i=pending.pop()
        for j in commutators:
            k=table[i][j]
            if k not in derived: derived.add(k);pending.append(k)
    assert len(derived)==128
    translations=[i for i,x in enumerate(group) if x[2]==G[0]]
    central_commutators={table[table[table[i][j]][inverses[i]]][inverses[j]] for i,j in product(translations,repeat=2)}
    assert {group[i] for i in central_commutators} == {(0,0,G[0],z) for z in range(8)}
    counts=Counter(); records=[]
    for e,c,g,z in group:
        if g==G[0]:continue
        value=[]
        for index,(a,b) in enumerate(((e,lam(e)),(e,c),(lam(e),c))):
            u=next(v for v in V if act(g,v)^v==a)
            value.append(dot(u,b)^((z>>index)&1))
        counts[''.join(map(str,value))]+=1
    assert dict(counts)=={''.join(map(str,v)):32 for v in product(range(2),repeat=3)}
    return {'joint_class_field_degree':48,'governing_compositum_degree':384,
            'faithful_joint_action_degree':24,'action_composition_checks':384**2,
            'central_commutator_subgroup_order':8,'commutator_subgroup_order':128,'abelianization':'C3',
            'inert_governing_vector_counts':dict(sorted(counts.items())),
            'conditional_density_per_full_CT_matrix':'1/8'}


def factor_degrees(coeff,p):
    trim,sub,mul,power,gcd=octic.finite_polynomials(p);h=trim(coeff)
    assert len(gcd(h,trim([i*h[i] for i in range(1,len(h))])))==1
    xpk=[0,1];counts={};degrees=[]
    for k in range(1,len(h)):
        xpk=power(xpk,p,h);total=len(gcd(sub(xpk,[0,1]),h))-1
        prior=sum(d*n for d,n in counts.items() if k%d==0)
        assert (total-prior)%k==0;counts[k]=(total-prior)//k
        assert counts[k]>=0;degrees += [k]*counts[k]
    assert sum(degrees)==len(h)-1
    return degrees


def radical_value(row,p):
    trim,sub,mul,power,gcd=octic.finite_polynomials(p);f=trim([-1,-14,-11,1])
    b=trim(row['beta_ascending']);N=p**3;odd=N-1;s=0
    while odd%2==0:odd//=2;s+=1
    nonsquare=next(v for v in range(2,p) if pow(v,(p-1)//2,p)==p-1)
    x=power(b,(odd+1)//2,f);t=power(b,odd,f);c=power([nonsquare],odd,f);m=s
    while t!=[1]:
        probe=t;i=0
        while probe!=[1]:probe=mul(probe,probe,f);i+=1
        assert i<m
        z=power(c,2**(m-i-1),f);x=mul(x,z,f);c=mul(z,z,f);t=mul(t,c,f);m=i
    assert mul(x,x,f)==b
    norm_root=25 if row['beta_ascending']==[12,-13,1] else 1
    exponent=p*p+p+1
    if power(x,exponent,f)!=[norm_root%p]:x=trim([-v for v in x])
    assert power(x,exponent,f)==[norm_root%p]
    xv=trim([r.mod(v,p) for v in row['X_ascending']])
    value=power(sub(xv,[-v for v in x]),exponent,f)
    if not value:return None
    assert len(value)==1
    return {'norm_mod_prime':value[0],'psi':int(pow(value[0],(p-1)//2,p)==p-1)}


def class_certificate():
    K=norm.Algebra([-1,-14,-11,1]);tau=K.elt([-2,-12,1]);eta=K.elt([-1,-1])
    eta1=K.evaluate(eta,tau);eta2=K.evaluate(eta1,tau)
    assert K.mul(K.mul(eta,eta1),eta2)==K.elt([1])
    intervals=[(Q(-11,10),Q(-101,100)),(Q(-1,13),Q(-1,14)),(Q(121,10),Q(13))]
    def evaluate(coeff,t):return sum(Q(c)*t**i for i,c in enumerate(coeff))
    for lo,hi in intervals:assert evaluate(K.f,lo)*evaluate(K.f,hi)<0
    signs=[]
    for coeff in (eta,eta1):
        row=[]
        for lo,hi in intervals:
            values=[evaluate(coeff,t) for t in (lo,hi)]
            if coeff[2] and lo < -coeff[1]/(2*coeff[2]) < hi:values.append(evaluate(coeff,-coeff[1]/(2*coeff[2])))
            assert min(values)>0 or max(values)<0
            row.append(int(max(values)<0))
        signs.append(row)
    assert r.rank(map(r.pack,signs))==2
    return {'eta_ascending':list(map(str,eta)),'tau_eta_ascending':list(map(str,eta1)),
            'eta_conjugate_product':1,'real_sign_rows':signs,'unit_orbit_dimension':2,
            'strict_orbit_dimension':2,'combined_dimension_over_cubic_field':4,
            'independence_reason':'The retained strict pair is independent and totally positive; the two unit conjugates have independent real signs.'}


def compute():
    inp=r.read(source.OUTPUT);prior=r.read(source.PRIOR);prior_verify=r.read(norm.OUTPUT)
    small=r.read(SMALL);dyadic=r.read(DYADIC);spec=r.read(source.PROTOCOL)
    for data in (inp,prior,prior_verify,dyadic):
        assert data['status']=='PASS'
        for path,sha in data['bindings'].items():assert r.digest((r.ROOT/path).read_bytes())==sha
    assert small['full_2_Selmer_dimensions']==[3,3]
    assert dyadic['CT_difference_matrix']==[[0,1],[1,0]]
    # eta is the Kummer class of the fixed original point (-1,1).
    assert (-1)**3-11*(-1)**2-14*(-1)-1==1
    orders=[arithmetic_pair(row,basis) for row,basis in zip(inp['rows'],BASES)]
    polynomials=[prior_verify['reduced_octic_ascending']]+[list(map(int,row['reduced_octic_ascending'])) for row in inp['rows']]
    norm_rows=[prior]+inp['rows'];table=[];excluded=[]
    for p in r.primes(spec['limits']['prime_bound']):
        if any(octic.discriminant(g)%p==0 for g in polynomials):
            excluded.append({'prime':p,'reason':'non-squarefree polynomial reduction; field ramification may be smaller'});continue
        f=[-1,-14,-11,1]
        if any(sum(c*x**i for i,c in enumerate(f))%p==0 for x in range(p)):continue
        bits=[];details=[]
        for row,g in zip(norm_rows,polynomials):
            degrees=factor_degrees(g,p);assert degrees in ([1,1,3,3],[2,6])
            bit=int(degrees==[2,6]);bits.append(bit)
            try:independent=radical_value(row,p)
            except ValueError:independent=None
            if independent is not None:assert bit==independent['psi']
            details.append({'factor_degrees':degrees,'independent_radical':independent})
        eligible=p%8==1 and pow(p%163,81,163)==1
        entry={'prime':p,'psi_vector':bits,'pair_replays':details,'twist_locally_square_at_2_163_infinity':eligible}
        if eligible:
            a,b,c=bits[0]^1,bits[1],bits[2]
            M=[[0,a,b],[a,0,c],[b,c,0]];rank=r.rank(map(r.pack,M))
            entry.update({'full_Selmer_dimension':3,'full_CT_matrix':M,'full_CT_rank':rank,
                          'full_CT_radical_dimension':3-rank,'MW_rank_upper_bound':3-rank,
                          'MW_rank_exact':'UNKNOWN','rational_solubility_of_radical':'UNKNOWN'})
        table.append(entry)
    eligible=[x for x in table if x['twist_locally_square_at_2_163_infinity']]
    return {'schema':'rank-jump.full-small-governing-block.v1','status':'PASS',
            'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in
                        (Path(__file__),source.OUTPUT,source.PRIOR,source.PROTOCOL,norm.OUTPUT,SMALL,DYADIC,Path(norm.__file__),Path(octic.__file__),Path(r.__file__))},
            'order_certificates':orders,'class_certificate':class_certificate(),'joint_group':joint_group(),
            'original_full_CT_matrix':[[0,1,0],[1,0,0],[0,0,0]],'joint_prime_table':table,'excluded_prime_reductions':excluded,
            'summary':{'joint_inert_primes':len(table),'independent_radical_replays':sum(d['independent_radical'] is not None for x in table for d in x['pair_replays']),
                       'eligible_prime_twists':[x['prime'] for x in eligible],
                       'eligible_zero_full_CT':[x['prime'] for x in eligible if x['full_CT_rank']==0],
                       'eligible_rank_at_most_one':[x['prime'] for x in eligible if x['full_CT_rank']==2]},
            'boundary':'Full Selmer/CT conclusions use the stated twist comparison theorem. No elliptic point search or twist descent is performed; zero full CT is not rational rank three.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);args=p.parse_args()
    result=compute()
    if args.mode=='build':r.write_new(OUTPUT,result)
    else:assert r.read(OUTPUT)==result
    print('PASS full governing block',result['summary'])
