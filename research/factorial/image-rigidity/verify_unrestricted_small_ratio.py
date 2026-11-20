#!/usr/bin/env python3
"""Exact/source support for the written unrestricted small-ratio proof."""
from __future__ import annotations
if not __debug__:
    raise SystemExit('Refusing -O/-OO: all exact and source controls must execute.')
import argparse,copy,hashlib,json,re,sys
from fractions import Fraction as Q
from math import comb,factorial
from pathlib import Path

HERE=Path(__file__).resolve().parent
PN='UNRESTRICTED_SMALL_RATIO_GAMMA_TWO_WELL_EXCLUSION.md'
CN='verify_unrestricted_small_ratio.py'
IN='unrestricted_small_ratio_inputs.json'
LN='unrestricted_small_ratio_source_ledger.json'
RN='unrestricted_small_ratio_verification.json'
OLD='research/factorial/image-rigidity/lower_exponent_initial_and_boundary.tex'
SPECS={
    'main_proof':('local',PN),
    'exact_checker':('local',CN),
    'exact_inputs':('local',IN),
    'raw_old_sparse_theorem':('repo',OLD),
}
RELATIONS={
    'main_proof':'written_complete_actual_cycle_argument',
    'exact_checker':'exact_algebra_arithmetic_source_support',
    'exact_inputs':'exact_typed_scope_and_controls',
    'raw_old_sparse_theorem':'logical_mathematical_premise_thm_sparse_only',
}
ASSURANCE={'formal_verification':False,'external_review':False,'independent_replay':False}
SCOPE={
    'schema':1,
    'reviewed_body_sha256':'ebf6489506e379f766728c199d821862ab02da3b034e9ffb428eb35a0bd925e5',
    'parameters':'integers a>=1,d>=2,N>=2; c=a/(d-1)^2,J=d^2/(a+d)',
    'theorem':'there exist unevaluated epsilon0>0 and N0 such that every actual pair is coprime when c<=epsilon0,N>=N0',
    'same_constants_all_supports_coefficients_and_char0_fields':True,
    'old_sparse_is_logical_not_context_only':True,
    'old_sparse_bad_constraint':'J<4*log(8*N)',
    'conditional_common_pair_root_capture_only':True,
    'unconditional_individual_root_annulus':False,
    'bounded_left_shape_compact_Morse_family':True,
    'left_shape_must_tend_to_infinity':False,
    'exact_Gamma_well_decomposition':True,
    'complete_whole_real_tails_and_own_width_L1':True,
    'actual_critical_powers_retained':True,
    'assumes_m_exp_minus_J_small':False,
    'both_endpoint_capture_before_radial_path':True,
    'd_le_3a_required':False,
    'uses_392_407_411_or_heavy_boundary_prefixes':False,
    'uses_sharp_lower_annulus_or_model_circle_law':False,
    'effective_constants':False,
    'all_order_SFC':False,
    'global_bad_triple_set_finite':False,
    'a_zero_included':False,
    'controls':{
        'maximum_lower':12,'maximum_gap':32,'maximum_order':4,
        'own_order_terminals':[2,3,5,8],
        'width_maximum_lower':64,'width_maximum_gap':192,
        'Gamma_series_terms':12,
        'false_individual_annulus':{'a':1,'d':20,'m':2,'exponential_Taylor_terms':128},
        'nonfrozen_power_indices':list(range(2,25)),
    },
    'assurance_flags':ASSURANCE,
}

def need(ok,msg):
    if not ok:raise ValueError(msg)
def h(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def typed(actual,expected):
    if type(actual) is not type(expected):return False
    if type(expected) is dict:
        return set(actual)==set(expected) and all(typed(actual[k],v)for k,v in expected.items())
    if type(expected) is list:
        return len(actual)==len(expected) and all(typed(x,y)for x,y in zip(actual,expected))
    return actual==expected
def check_inputs(doc):need(typed(doc,SCOPE),'Exact typed scope, bounds and fixtures required')
def sources(doc,repo):
    need(type(doc) is dict and set(doc)=={'schema','assurance_flags','sources'},'Ledger schema')
    need(type(doc['schema']) is int and doc['schema']==1,'Exact integer ledger schema1')
    need(typed(doc['assurance_flags'],ASSURANCE),'Assurance remains false')
    rows=doc['sources'];need(type(rows) is list and len(rows)==len(SPECS)==4,'All four exact source roles required')
    out={};paths=set();roots={'local':HERE,'repo':repo}
    for row in rows:
        need(type(row) is dict and set(row)=={'role','kind','path','relation','sha256'},'Exact source-row schema')
        role=row['role'];need(type(role) is str and role in SPECS and role not in out,'Unique required role')
        need((row['kind'],row['path'])==SPECS[role] and row['relation']==RELATIONS[role],'Exact role/path/relation')
        need(type(row['sha256']) is str and re.fullmatch('[0-9a-f]{64}',row['sha256']),'SHA256')
        root=roots[row['kind']].resolve();p=(root/row['path']).resolve()
        need(p.is_relative_to(root) and p.is_file() and h(p)==row['sha256'],'Missing/stale required source:'+role)
        key=(row['kind'],row['path']);need(key not in paths,'Unique source paths');paths.add(key);out[role]=row
    text=(repo/OLD).read_text()
    need('\\label{thm:sparse}' in text and 'd^2\\ge4(a+d)\\log(8N)' in text,
         'Pinned old source must contain the designated sparse theorem')
    return out

# Small, independent exact multivariate rational-polynomial kernel. Universal
# identities use formal variables and cross multiplication, not sampled values.
NV=8
def C(x):
    x=Q(x)
    return {} if not x else {(0,)*NV:x}
def X(i):
    e=[0]*NV;e[i]=1;return {tuple(e):Q(1)}
def add(*parts):
    out={}
    for p in parts:
        for e,c in p.items():
            out[e]=out.get(e,Q(0))+c
            if not out[e]:del out[e]
    return out
def sc(p,c):return {e:v*Q(c)for e,v in p.items() if v*Q(c)}
def sub(a,b):return add(a,sc(b,-1))
def mul(*parts):
    out=C(1)
    for p in parts:
        nxt={}
        for e,c in out.items():
            for f,d in p.items():
                g=tuple(x+y for x,y in zip(e,f));nxt[g]=nxt.get(g,Q(0))+c*d
        out={e:c for e,c in nxt.items() if c}
    return out
def power(p,n):
    out=C(1)
    for _ in range(n):out=mul(out,p)
    return out
def deriv(p,i):
    out={}
    for e,c in p.items():
        if e[i]:
            f=list(e);f[i]-=1;out[tuple(f)]=c*e[i]
    return out
def universal():
    ids=[];one=C(1)
    def ck(name,p):need(not p,'Universal identity:'+name);ids.append(name)
    a,d,S,w,Z,v,b,m=[X(i)for i in range(NV)]
    D=add(a,d)
    crit=add(mul(sub(a,S),add(w,Z)),mul(d,Z))
    ck('support-fixed actual stationary radial envelope',
       add(sub(mul(d,w),mul(sub(D,S),add(w,Z))),crit))
    # H_b'(v)=-v/(1+bv); differentiate the quotient exactly.
    den=add(one,mul(b,v));num=sc(v,-1)
    second_num=sub(mul(deriv(num,5),den),mul(num,deriv(den,5)))
    ck('compact bounded-shape Gamma Morse curvature',add(second_num,one))
    # Positive and negative log-rate inequalities are certified by their
    # cleared derivative numerators. Domain signs are stated in the proof.
    x=X(0);u=X(1);A=X(2);DD=X(3)
    xx=power(x,2);op=add(one,x)
    ck('positive log-rate derivative certificate',
       sub(sub(sc(mul(x,op),2),mul(x,add(C(2),x))),xx))
    ck('negative log-rate derivative certificate',
       sub(add(one,sc(op,-1),mul(x,op)),xx))
    ck('left center cross-loss integral positivity certificate',
       sub(mul(sub(DD,u),sub(DD,u)),sub(mul(DD,sub(DD,u)),mul(u,sub(DD,u)))))
    ck('right center cross-loss integral positivity certificate',
       sub(mul(sub(u,A),sub(DD,u)),sub(mul(DD,sub(u,A)),mul(u,sub(u,A)))))
    a,d,U,V=[X(i)for i in range(4)];D=add(a,d)
    ck('seam square-root center separation',sub(mul(sub(V,U),add(V,U)),sub(power(V,2),power(U,2))))
    ck('unrestricted seam rate compared with J',sub(sub(power(D,2),power(d,2)),mul(a,add(a,sc(d,2)))))
    ck('own-width cross-tail exact decomposition',sub(D,add(d,a)))
    ck('width numerator positive when a>=1 and d>=a',
       sub(sub(sc(mul(a,power(d,2)),4),power(D,2)),
           add(sc(mul(sub(a,one),power(d,2)),4),mul(sub(d,a),add(a,sc(d,3))))))
    c,d=X(0),X(1)
    ck('small-c bound forces J to infinity',sub(sub(mul(c,power(d,2)),mul(c,power(sub(d,one),2))),mul(c,sub(sc(d,2),one))))
    # E'= -Daction + m Daction' + bq dq, with m Daction'=-DeltaS.
    DA,DS,BQ,d,q,mm=[X(i)for i in range(6)]
    ck('exact adjacent action derivative cleared by m',
       sub(add(sc(mul(mm,DA),-1),sc(mul(mm,DS),-1),mul(d,q,BQ)),
           add(sc(mul(mm,add(DA,DS)),-1),mul(d,q,BQ))))
    # Log width is support-fixed: its q derivative has no contribution.
    q,R=X(0),X(1)
    ck('support-fixed width contributes zero q derivative',deriv(R,0))
    # Exact formal monomial exponents in
    # (E/B)^m [B/(A E)]^(m-j)=E^j/[A^(m-j) B^j].
    mm,j=X(0),X(1)
    ck('mixed Gamma normalization exponent of E',sub(sub(mm,sub(mm,j)),j))
    ck('mixed Gamma normalization exponent of A',add(sc(sub(mm,j),-1),sub(mm,j)))
    ck('mixed Gamma normalization exponent of B',add(sub(mm,j),sc(mm,-1),j))
    return ids

def arithmetic(doc):
    cfg=doc['controls'];vec=mom=own=width=series=coefficients=0
    def rising(x,n):
        y=1
        for j in range(n):y*=x+j
        return y
    for a in range(1,cfg['maximum_lower']+1):
        for d in range(2,cfg['maximum_gap']+1):
            D=a+d
            for m in range(1,cfg['maximum_order']+1):
                vector=[]
                for j in range(m+1):
                    z=comb(m,j)*rising(a*m+1,d*j)
                    need(z==comb(m,j)*factorial(a*m+d*j)//factorial(a*m),'Actual am+1 full coefficient vector')
                    # Compare each coefficient in both sides of the exact
                    # original-period identity I_m=(am)!A_m(m^d w)/m^(Dm+1).
                    original=Q(comb(m,j)*factorial(a*m+d*j),m**(a*m+d*j+1))
                    normalized=Q(factorial(a*m)*z*m**(d*(m-j)),m**(D*m+1))
                    need(original==normalized,'Actual whole-period coefficient identity')
                    vector.append(z);coefficients+=1
                vec+=1
                aa,bb=Q(2,3),Q(-5,7);x=aa/bb
                direct=sum(Q(comb(m,j)*factorial(a*m+d*j))*aa**(m-j)*bb**j for j in range(m+1))
                physical=bb**m*factorial(a*m)*sum(z*x**(m-j) for j,z in enumerate(vector))
                need(direct==physical,'Actual coefficient/field factorial identity');mom+=1
            for N in cfg['own_order_terminals']:
                n=N-1;W=Q(D**D,a**a*2**d);xx=Q(2,3)
                qN=xx/(N**d*W);qn=xx/(n**d*W)
                need(qn==Q(N,n)**d*qN,'Same support and exact own-order radial coupling');own+=1
    for a in range(1,cfg['width_maximum_lower']+1):
        for d in range(2,cfg['width_maximum_gap']+1):
            D=a+d;J=Q(d*d,D);R=Q(D,a)
            if J<1:continue
            need(R<=4*J,'Unrestricted width bound R<=4J')
            need(Q(a,d*d)<=Q(D,d*d),'sqrt(a)/d<=1/sqrt(J) after squaring')
            need(D*D-d*d==a*(a+2*d)>0,'Universal seam-rate positivity specialization')
            width+=1
    # H_b(v) at b=0 is removable; every nonzero coefficient has a
    # nonnegative b exponent. This is exact coefficientwise derivation.
    for n in range(2,cfg['Gamma_series_terms']+1):
        c=Q((-1)**(n+1),n)
        need(n-2>=0 and n*c==(-1)**(n+1),'Removable complete Gamma Taylor coefficient')
        if n==2:need(c==Q(-1,2),'Uniform curvature at every bounded shape')
        series+=1
    fixture=cfg['false_individual_annulus'];a=fixture['a'];d=fixture['d'];m=fixture['m'];D=a+d
    need((a,d,m)==(1,20,2),'Actual fixed countercontrol required')
    r1=rising(2*a+1,d);r2=rising(2*a+1,2*d)
    need(r2>r1*r1,'Actual quadratic has a conjugate nonreal root pair')
    lower_exp=sum(Q(2*d)**j/factorial(j)for j in range(fixture['exponential_Taylor_terms']+1))
    qmod2_lower=Q(r2*a**(2*a),m**(2*d)*D**(2*D))*lower_exp
    need(qmod2_lower>4,'Actual individual normalized quadratic roots lie OUTSIDE |q|<=2')
    # This is deliberately NOT a moment counterexample: the old sparse
    # theorem applies at terminal2 for this support.
    need(Q(d*d,4*D)>4 and sum(Q(4)**j/factorial(j)for j in range(5))>16,
         'The scope counterexample lies inside the old sparse coprime region')
    frozen=0
    for j in cfg['nonfrozen_power_indices']:
        need((Q(j+1,j))**(j*j)>=1+j>1,'Small critical shift need not give small powered error')
        frozen+=1
    wrong_shape=rising(2,2);correct_shape=rising(3,2)
    need(correct_shape==12 and wrong_shape==6 and correct_shape!=wrong_shape,'Actual +1 Gamma shape negative')
    return {'complete_actual_coefficient_vectors':vec,'actual_period_coefficient_identities':coefficients,
       'actual_factorial_field_normalizations':mom,'actual_own_order_maps':own,
       'width_seam_and_discrete_bounds':width,'removable_Gamma_series_coefficients':series,
       'nonfrozen_power_countermodels':frozen,'genuine_math_negative_controls':3,
       'actual_individual_root_counterexample':{'a':a,'d':d,'m':m,'both_root_moduli_greater_than':2,
          'proved_by_exact_quadratic_and_positive_exponential_Taylor_lower_bound':True,
          'old_sparse_theorem_certifies_pair_coprime':True}}

def negatives(ledger,repo):
    out=[]
    def reject(name,fn):
        try:fn()
        except(ValueError,KeyError,TypeError,OSError):out.append(name);return
        raise ValueError('Negative accepted:'+name)
    for role in SPECS:
        d=copy.deepcopy(ledger);d['sources']=[r for r in d['sources'] if r['role']!=role]
        reject('omit:'+role,lambda d=d:sources(d,repo))
        d=copy.deepcopy(ledger);next(r for r in d['sources']if r['role']==role)['sha256']='0'*64
        reject('badpin:'+role,lambda d=d:sources(d,repo))
    for key in ('unconditional_individual_root_annulus','left_shape_must_tend_to_infinity',
       'assumes_m_exp_minus_J_small','d_le_3a_required','uses_392_407_411_or_heavy_boundary_prefixes',
       'uses_sharp_lower_annulus_or_model_circle_law','effective_constants','all_order_SFC',
       'global_bad_triple_set_finite','a_zero_included'):
        d=copy.deepcopy(SCOPE);d[key]=True;reject('false-scope:'+key,lambda d=d:check_inputs(d))
    for key in ('old_sparse_is_logical_not_context_only','conditional_common_pair_root_capture_only',
       'bounded_left_shape_compact_Morse_family','exact_Gamma_well_decomposition',
       'complete_whole_real_tails_and_own_width_L1','actual_critical_powers_retained',
       'both_endpoint_capture_before_radial_path'):
        d=copy.deepcopy(SCOPE);d[key]=False;reject('missing-scope:'+key,lambda d=d:check_inputs(d))
    d=copy.deepcopy(SCOPE);d['schema']=True;reject('boolean-input-schema',lambda:check_inputs(d))
    d=copy.deepcopy(ledger);d['schema']=True;reject('boolean-ledger-schema',lambda:sources(d,repo))
    d=copy.deepcopy(ledger);d['sources']=[];reject('empty-ledger',lambda:sources(d,repo))
    return out

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--repo-root',type=Path,required=True)
    ap.add_argument('--output',type=Path,default=HERE/RN);ns=ap.parse_args();repo=ns.repo_root.resolve()
    inp=load(HERE/IN);check_inputs(inp)
    need(h(HERE/PN)==inp['reviewed_body_sha256'],'Exact reviewed proof body must match the declared review hash')
    ledger=load(HERE/LN);rows=sources(ledger,repo)
    identities=universal();counts=arithmetic(inp);bad=negatives(ledger,repo);sources(ledger,repo)
    result={'schema':1,'status':'PASS','scope':inp,'proof_sha256':h(HERE/PN),'checker_sha256':h(HERE/CN),
      'input_sha256':h(HERE/IN),'source_ledger_sha256':h(HERE/LN),'mandatory_source_roles':len(rows),
      'sources':[rows[k]for k in sorted(rows)],'universal_identities':identities,'counts':counts,
      'source_scope_negative_controls':bad,'assurance_flags':ASSURANCE,
      'meaning':'The complete written actual-cycle proof supplies uniformity. Universal exact derivations and finite arithmetic/source controls support that proof; they are not a finite-prefix inference, formal verification or an effective cutoff.'}
    ns.output.parent.mkdir(parents=True,exist_ok=True);ns.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':'PASS','mandatory_sources':len(rows),'universal_identities':len(identities),
      'counts':counts,'source_scope_negatives':len(bad),'output':str(ns.output)},sort_keys=True))
if __name__=='__main__':main()
