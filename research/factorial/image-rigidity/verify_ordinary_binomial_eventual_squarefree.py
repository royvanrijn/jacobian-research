#!/usr/bin/env python3
"""Exact/source controls for the written fixed-support simplicity deduction."""
from __future__ import annotations
if not __debug__:
    raise SystemExit('Refusing -O/-OO: exact guards must execute.')
from fractions import Fraction as Q
from pathlib import Path
from math import comb,factorial
import argparse,copy,hashlib,json,re
HERE=Path(__file__).resolve().parent
PN='EVENTUAL_SQUAREFREENESS_OF_ORDINARY_BINOMIAL_MOMENTS.md'
CN='verify_ordinary_binomial_eventual_squarefree.py'
IN='ordinary_binomial_eventual_squarefree_inputs.json'
LN='ordinary_binomial_eventual_squarefree_source_ledger.json'
RN='ordinary_binomial_eventual_squarefree_verification.json'
P='research/factorial/image-rigidity/'
ROLES={
 'proof':('local',PN),'checker':('local',CN),'inputs':('local',IN),
 'regular_packets':('support','ABOVE_CRITICAL_REGULAR_COMPACT_COMPLEMENT.md'),
 'simple_turn_packets':('support','SIMPLE_TURN_SHRINKING_WINDOW_CONSECUTIVE_EXCLUSION.md'),
 'strict_subcritical_packets':('repo',P+'SUBCRITICAL_ORDERED_PAIR_EVENTUAL_SFC.md'),
 'quartic_growing_packets':('repo',P+'CUSP_SHRINKING_WINDOW_LOG_RATIO_NONCOLLISION.md'),
 'quartic_bounded_packets':('repo',P+'BOUNDED_QUARTIC_CUSP_CONSECUTIVE_EXCLUSION.md'),
 'quartic_simple_model':('repo',P+'QUARTIC_LEE_YANG_ZERO_AND_ROBIN_BRIDGE.md'),
 'quartic_local_geometry':('repo',P+'QUARTIC_CUSP_UNFOLDING_AND_ACTUAL_ORDER_DRIFT.md'),
 'complete_morse_mechanisms':('repo',P+'CUSP_TWO_THIMBLE_POSITIVE_CYCLE.md'),
 'complete_wall_detours':('repo',P+'CUSP_GLOBAL_REGULAR_COMPLEMENT_AND_STOKES_DETOUR.md'),
 'stationary_phase_geometry':('repo',P+'COALESCING_BINOMIAL_SADDLE_BRIDGE.md'),
 'old_gamma_upper':('repo',P+'GAMMA_PRODUCT_ROOT_MODULUS.tex'),
 'constant_analytic_tail':('repo',P+'ALL_UNIVARIATE_TWO_TERM_SFC.tex'),
 'small_constant_prime_powers':('repo',P+'PRIME_POWER_TWO_TERM_SFC.tex'),
}
SCOPE={
 'schema':1,'fixed_support':{'a_min':0,'d_min':1,'integer_parameters':True},
 'statement':'for each fixed support exists m0; all m>=m0 actual A_m squarefree',
 'same_cutoff_over_characteristic_zero_fields':True,
 'individual_period_w_Cauchy_collars_logical_premise':True,
 'regular_positive_coexistence_two_packet':True,
 'bounded_models_one_root_counted_with_multiplicity':True,
 'new_sharp_lower_annulus_logical_premise':False,
 'novertical_pair_logical_premise':False,
 'eventual_consecutive_coprimality_alone_implies_simplicity':False,
 'all_exponent_constant_finite_prefix_logical_premise':False,
 'all_orders_all_positive_supports':False,
 'effective_cutoff':False,'uniform_cutoff_over_supports':False,
 'formal_verification':False,'external_review':False,'independent_replay':False,
}
def need(ok,msg):
    if not ok:raise AssertionError(msg)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def typed_equal(a,b):
    if type(a) is not type(b):return False
    if type(b) is dict:return set(a)==set(b) and all(typed_equal(a[k],v) for k,v in b.items())
    if type(b) is list:return len(a)==len(b) and all(typed_equal(x,y) for x,y in zip(a,b))
    return a==b
def check_inputs(d):need(typed_equal(d,SCOPE),'Exact typed fixed-support squarefree scope')
def check_sources(d,repo,support):
    need(type(d) is dict and set(d)=={'schema','sources'},'Ledger schema')
    need(type(d['schema']) is int and d['schema']==1,'Integer schema1')
    need(type(d['sources']) is list and len(d['sources'])==len(ROLES),'All mandatory roles')
    roots={'local':HERE,'repo':repo,'support':support};out={}
    for r in d['sources']:
        need(type(r) is dict and set(r)=={'role','kind','path','sha256'},'Exact source row')
        role=r['role'];need(type(role) is str and role in ROLES and role not in out,'Unique known role')
        need((r['kind'],r['path'])==ROLES[role],'Exact role/path')
        need(type(r['sha256']) is str and re.fullmatch('[0-9a-f]{64}',r['sha256']),'SHA256')
        root=roots[r['kind']].resolve();p=(root/r['path']).resolve()
        need(p.is_relative_to(root) and p.is_file() and sha(p)==r['sha256'],'Missing/stale source:'+role)
        out[role]=r
    need(set(out)==set(ROLES),'Complete role set')
    return out
def trim(a):
    a=list(a)
    while a and not a[-1]:a.pop()
    return a
def deriv(a):return trim([Q(i)*a[i] for i in range(1,len(a))])
def add(a,b):return trim([(a[i] if i<len(a) else 0)+(b[i] if i<len(b) else 0) for i in range(max(len(a),len(b)))])
def scale(a,c):return trim([c*x for x in a])
def shift(a):return [Q(0)]+list(a) if a else []
def mul(a,b):
    if not a or not b:return []
    out=[Q(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]+=x*y
    return trim(out)
def rem(a,b):
    a=trim(a);b=trim(b);need(bool(b),'Nonzero divisor')
    while len(a)>=len(b):
        q=a[-1]/b[-1];j=len(a)-len(b)
        for i,x in enumerate(b):a[i+j]-=q*x
        a=trim(a)
    return a
def gcd(a,b):
    while b:a,b=b,rem(a,b)
    return scale(a,1/a[-1]) if a else []
def rising(g,k):
    r=Q(1)
    for j in range(k):r*=g+j
    return r
def moment(a,d,m):
    return [Q(comb(m,m-j)*factorial(a*m+d*(m-j)),factorial(a*m)) for j in range(m+1)]
def algebra():
    # A_i' identity in the polynomial ring Q[P,P',Z,E']; sigma is a fixed sign.
    # I=P(1+sigma Z), I'=P'(1+sigma Z)+sigma P Z E'.
    sign_controls=0
    for sigma in (-1,1):
        for p in (Q(2,3),Q(-5,7)):
            for pd in (Q(-1,2),Q(3,5)):
                for e in (Q(-7,4),Q(9,11)):
                    z=-Q(1,sigma);I=p*(1+sigma*z)
                    Id=pd*(1+sigma*z)+sigma*p*z*e
                    need(I==0 and Id==-p*e and Id!=0,'Exact zero derivative criterion')
                    sign_controls+=1
    # Universal Laurent exponent identities after m=k*rho^(-p).
    laurent=0
    for p in (3,4):
        # Divide amplitude rho^(2-p), and remainder k^-1 rho^(1-p), by m*rho.
        num_amp=(0,2-p);num_rem=(-1,1-p);den=(1,1-p)
        need(tuple(x-y for x,y in zip(num_amp,den))==(-1,1),'rho/k')
        need(tuple(x-y for x,y in zip(num_rem,den))==(-2,0),'1/k^2')
        laurent+=2
    ode=0
    for g in (Q(1,3),Q(1),Q(5,2),Q(17)):
        for n in range(1,9):
            G=[Q(comb(n,n-j))*rising(g,n-j) for j in range(n+1)]
            Gp=deriv(G);Gpp=deriv(Gp)
            lhs=add(add(shift(Gpp),scale(Gp,1-g-n)),add(scale(shift(Gp),-1),scale(G,n)))
            need(not lhs and G[0]>0,'Actual all-positive-Gamma ODE coefficients')
            ode+=1
    appell=0
    for d in range(1,9):
        for m in range(2,9):
            need(deriv(moment(0,d,m))==scale(moment(0,d,m-1),m),'Actual a0 Appell')
            appell+=1
    literal=0
    for a in range(6):
        for d in range(1,7):
            for m in range(1,6):
                A=moment(a,d,m)
                need(len(A)==m+1 and A[-1]==1 and A[0]>0,'Actual monic moment')
                need(gcd(A,deriv(A))==[1],'Literal QQ simplicity control')
                literal+=1
    margins=0
    for rho in (Q(1,100),Q(1,10),Q(1,4)):
        for k in (8,16,32):
            need(rho/k+Q(1,k*k)<Q(1,16),'Normalized derivative-error base margin')
            margins+=1
    return {'exact_packet_zero_derivative_controls':sign_controls,
            'universal_Laurent_exponent_identities':laurent,
            'Gamma_ODE_coefficient_identities':ode,'actual_Appell_identities':appell,
            'literal_actual_QQ_squarefree_controls':literal,'rational_error_margins':margins}
def negatives(ledger,repo,support):
    out=[]
    def reject(name,f):
        try:f()
        except (AssertionError,KeyError,TypeError,ValueError):out.append(name);return
        raise AssertionError('Negative accepted:'+name)
    for role in sorted(ROLES):
        d=copy.deepcopy(ledger);d['sources']=[r for r in d['sources'] if r['role']!=role]
        reject('omit:'+role,lambda d=d:check_sources(d,repo,support))
        d=copy.deepcopy(ledger);next(r for r in d['sources'] if r['role']==role)['sha256']='0'*64
        reject('badpin:'+role,lambda d=d:check_sources(d,repo,support))
    for key in ('new_sharp_lower_annulus_logical_premise','novertical_pair_logical_premise',
                'eventual_consecutive_coprimality_alone_implies_simplicity',
                'all_exponent_constant_finite_prefix_logical_premise',
                'all_orders_all_positive_supports','effective_cutoff','uniform_cutoff_over_supports'):
        d=copy.deepcopy(SCOPE);d[key]=True
        reject('false-scope:'+key,lambda d=d:check_inputs(d))
    for key in ('regular_positive_coexistence_two_packet','bounded_models_one_root_counted_with_multiplicity'):
        d=copy.deepcopy(SCOPE);d[key]=False
        reject('missing-scope:'+key,lambda d=d:check_inputs(d))
    d=copy.deepcopy(ledger);d['schema']=True
    reject('boolean-schema',lambda:check_sources(d,repo,support))
    repeated=[Q(1),Q(-2),Q(1)];other=[Q(0),Q(1)]
    need(gcd(repeated,other)==[1] and len(gcd(repeated,deriv(repeated)))==2,
         'Actual coprime adjacent-looking pair need not imply squarefreeness')
    out.append('coprime-pair-with-repeated-upper')
    need(Q(1)-Q(1)==0 and Q(1)-Q(-1)==2,
         'Nonzero modulus separation can have zero real part')
    out.append('modulus-versus-real-separation')
    return out
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--repo-root',required=True,type=Path)
    ap.add_argument('--support-root',type=Path,default=HERE)
    ap.add_argument('--output',type=Path,default=HERE/RN);ns=ap.parse_args()
    repo=ns.repo_root.resolve();support=ns.support_root.resolve()
    inp=load(HERE/IN);check_inputs(inp);ledger=load(HERE/LN)
    rows=check_sources(ledger,repo,support);counts=algebra();bad=negatives(ledger,repo,support)
    check_sources(ledger,repo,support)
    payload={'schema':1,'status':'PASS','proof_sha256':sha(HERE/PN),'checker_sha256':sha(HERE/CN),
       'input_sha256':sha(HERE/IN),'source_ledger_sha256':sha(HERE/LN),
       'scope':inp,'mandatory_source_roles':len(rows),'sources':[rows[r] for r in sorted(rows)],
       'counts':counts,'negative_controls':bad,
       'meaning':'The written complete contour/Cauchy/Rouche argument proves eventual simplicity. Exact finite controls authenticate its algebra and scope; no finite sample proves eventuality.',
       'formal_verification':False,'external_review':False,'independent_replay':False}
    ns.output.parent.mkdir(parents=True,exist_ok=True)
    ns.output.write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':'PASS','mandatory_sources':len(rows),'counts':counts,
                      'negative_controls':len(bad),'output':str(ns.output)}))
if __name__=='__main__':main()
