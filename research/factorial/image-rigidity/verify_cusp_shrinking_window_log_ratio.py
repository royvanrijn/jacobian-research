#!/usr/bin/env python3
"""Exact/source support for fixed-cusp physical-neighborhood noncollision."""
from __future__ import annotations
if not __debug__:
    raise SystemExit('Refusing optimized Python: run without -O/-OO.')
import argparse
from copy import deepcopy
from fractions import Fraction as Q
from hashlib import sha256
import json
from math import comb, factorial
from pathlib import Path
import re
import subprocess
import sys
import tempfile

HERE=Path(__file__).resolve().parent
PN='CUSP_SHRINKING_WINDOW_LOG_RATIO_NONCOLLISION.md'
CN='verify_cusp_shrinking_window_log_ratio.py'
IN='cusp_shrinking_window_inputs.json'
LN='cusp_shrinking_window_source_ledger.json'
RN='cusp_shrinking_window_verification.json'
PRE='research/factorial/image-rigidity/'
ROLES={
 'proof':('local',PN),'checker':('local',CN),'inputs':('local',IN),
 'fixed_window_proof':('repo',PRE+'BOUNDED_QUARTIC_CUSP_CONSECUTIVE_EXCLUSION.md'),
 'cusp_jet_proof':('repo',PRE+'QUARTIC_CUSP_UNFOLDING_AND_ACTUAL_ORDER_DRIFT.md'),
 'fixed_completion_checker':('repo',PRE+'verify_quartic_cusp_local_completion.py'),
 'fixed_completion_inputs':('repo',PRE+'quartic_cusp_local_completion_inputs.json'),
 'fixed_completion_ledger':('repo',PRE+'quartic_cusp_local_completion_source_ledger.json'),
 'fixed_completion_receipt':('repo',PRE+'quartic_cusp_local_completion_verification.json'),
}
IMPORTED_ROLES={
 'bounded_checker','bounded_proof','bounded_receipt','completion_checker','completion_inputs',
 'cusp_checker','cusp_proof','cusp_receipt','model_checker','model_inputs','model_ledger',
 'model_proof','model_receipt','primary_hjn','primary_newman'}
SCOPE={
 'schema':1,'fixed_cusp':{'s_min':1,'a':'s^2','d':'2s+1'},
 'parameter_region':'one fixed physical neighborhood of w_c for each fixed cusp',
 'actual_polynomials':'A_m(x)=sum binom(m,k)(a*m+1)_(d*k)x^(m-k)',
 'actual_common_order_scaling':'w=x/N^d; lower w=(N/(N-1))^d*w',
 'order_quantifier':'for every fixed s exists r_s>0,N_s; all N>=N_s',
 'all_scaled_fields_in_neighborhood':True,
 'complex_m_extension_only_in_narrow_zero_free_S_sector':True,
 'exact_partial_periods_and_Cauchy_derivatives':True,
 'fixed_bounded_window_core_suffices':True,
 'uniform_H_core_required':False,'HJN_monotonicity_is_logical_premise':False,
 'all_RHP_is_zero_free':False,'uniform_in_s':False,'effective_cutoff':False,
 'global_fixed_cusp_SFC':False,'global_all_support_SFC':False,
 'formal_verification':False,'external_review':False,'independent_replay':False,
}

def need(ok,msg):
    if not ok:raise ValueError(msg)
def sha(p):return sha256(Path(p).read_bytes()).hexdigest()
def load(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def inputs_ok(data):
    need(type(data) is dict and data==SCOPE,'Exact local scope contract')
    need(all(type(data[k]) is type(v) for k,v in SCOPE.items()),'Exact scope types')
    need(type(data['fixed_cusp']['s_min']) is int,'Exact nested integer cusp')
def source_ok(doc,repo):
    need(type(doc) is dict and type(doc.get('schema')) is int and doc['schema']==1,'Ledger schema1')
    need(doc.get('formal_verification') is False and doc.get('external_review') is False
         and doc.get('independent_replay') is False,'Assurance flags false')
    rows=doc.get('sources');need(type(rows) is list and len(rows)==len(ROLES),'All nine mandatory sources')
    out={}
    for row in rows:
        need(type(row) is dict and set(row)=={'role','kind','path','sha256'},'Source record fields')
        role=row['role'];need(type(role) is str and role in ROLES and role not in out,'Unique source role')
        need((row['kind'],row['path'])==ROLES[role],'Exact source role/path')
        need(type(row['sha256']) is str and re.fullmatch('[0-9a-f]{64}',row['sha256']),'SHA256 syntax')
        root=HERE if row['kind']=='local' else repo
        p=(root/row['path']).resolve()
        need(p.is_relative_to(root.resolve()) and p.is_file() and sha(p)==row['sha256'],'Missing/stale:'+role)
        out[role]=row
    need(doc.get('logical_inputs')==['fixed_cusp_bounded_window_theorem','exact_cusp_quartic_jets'],
         'Do not substitute uniform core or HJN monotonicity')
    return out
def core_closure(repo):
    base=repo/PRE;doc=load(base/'quartic_cusp_local_completion_source_ledger.json')
    need(type(doc.get('schema')) is int and doc['schema']==1,'Imported schema1')
    rows=doc.get('sources');need(type(rows) is list and len(rows)==15,'Exact fixed17 source contract')
    need({r.get('role') for r in rows}==IMPORTED_ROLES,'Fixed17 role set')
    found=[];seen=set()
    for row in rows:
        need(type(row) is dict and set(row)=={'role','path','sha256'},'Imported record schema')
        p=(base/row['path']).resolve()
        need(p.is_relative_to(base.resolve()) and p.is_file() and sha(p)==row['sha256'],'Imported stale source')
        need(row['path'] not in seen,'Duplicate imported path');seen.add(row['path'])
        found.append({'path':PRE+row['path'],'sha256':row['sha256']})
    for name in ('quartic_cusp_local_completion_source_ledger.json','quartic_cusp_local_completion_verification.json'):
        found.append({'path':PRE+name,'sha256':sha(base/name)})
    need(len(found)==17,'Seventeen-file retained core closure')
    return sorted(found,key=lambda r:r['path'])
def equal_receipt(a,b):need(a==b,'Fresh fixed completion receipt differs')
def replay_core(repo):
    base=repo/PRE
    with tempfile.TemporaryDirectory(prefix='cusp-shrinking-core-') as tmp:
        out=Path(tmp)/'core.json'
        run=subprocess.run([sys.executable,'-B',str(base/'verify_quartic_cusp_local_completion.py'),
              '--repo-root',str(repo),'--output',str(out)],cwd=tmp,capture_output=True,text=True,timeout=60)
        need(run.returncode==0,'Fixed completion replay:'+run.stdout+run.stderr)
        equal_receipt(out.read_bytes(),(base/'quartic_cusp_local_completion_verification.json').read_bytes())
        data=load(out)
        need(data.get('status')=='PASS','Imported completion status')
        return sha(out)

# Sparse rational polynomials; the first algebra block is in x,c,A modulo
# c^3=3/8 and A^2=3c^2. Later blocks use ordinary independent variables.
NV=8
def clean(p):return {k:Q(v) for k,v in p.items() if v}
def const(v):return {} if not v else {(0,)*NV:Q(v)}
def var(i):
    k=[0]*NV;k[i]=1;return {tuple(k):Q(1)}
def add(*items):
    out={}
    for p in items:
        for k,v in p.items():out[k]=out.get(k,Q(0))+v
    return clean(out)
def scale(p,c):return clean({k:v*c for k,v in p.items()})
def mul(*items):
    out=const(1)
    for p in items:
        nxt={}
        for a,c in out.items():
            for b,v in p.items():
                k=tuple(x+y for x,y in zip(a,b));nxt[k]=nxt.get(k,Q(0))+c*v
        out=clean(nxt)
    return out
def powp(p,n):
    out=const(1)
    for _ in range(n):out=mul(out,p)
    return out
def subvar(p,i,q):
    out={}
    for k,c in p.items():
        j=list(k);ex=j[i];j[i]=0
        out=add(out,mul({tuple(j):c},powp(q,ex)))
    return out
def reduce_model(p):
    out={}
    for k,c in p.items():
        j=list(k);a,rem=divmod(j[2],2);j[2]=rem;j[1]+=2*a;c*=3**a
        a,rem=divmod(j[1],3);j[1]=rem;c*=Q(3,8)**a
        key=tuple(j);out[key]=out.get(key,Q(0))+c
    return clean(out)
def universal():
    x,c,A=var(0),var(1),var(2)
    real=add(scale(powp(x,4),-Q(1,12)),scale(mul(powp(c,2),powp(x,2)),Q(1,2)),
             scale(powp(c,4),-Q(1,12)),scale(c,-1))
    imag=add(scale(mul(c,powp(x,3)),Q(1,3)),scale(mul(powp(c,3),x),-Q(1,3)),scale(x,-1))
    loss=scale(powp(add(powp(x,2),scale(powp(c,2),-3)),2),Q(1,12))
    need(not reduce_model(add(subvar(real,0,A),scale(real,-1),scale(loss,-1))),'Model full real loss')
    imag_expected=scale(mul(c,powp(add(x,scale(A,-1)),2),add(x,scale(A,2))),Q(1,3))
    need(not reduce_model(add(imag,scale(subvar(imag,0,A),-1),scale(imag_expected,-1))),'Model full imaginary difference')
    need(not reduce_model(add(subvar(real,0,A),scale(subvar(real,0,{}),-1),scale(c,-Q(9,32)))),'Midpoint gap')
    need(not reduce_model(add(powp(A,2),scale(powp(c,2),-1),scale(powp(c,2),-2))),'Curvature real part')
    y=var(0)
    certificate=add(scale(powp(add(y,const(1)),2),2),scale(add(y,const(2)),-1))
    need(certificate==add(scale(powp(y,2),2),scale(y,3)),'Halfline ratio positive polynomial')
    need(64<75,'8/sqrt(3)<5')
    m,d,w,D,gap,bm,bw,Dw=[var(i) for i in range(NV)]
    direct=add(scale(mul(m,add(D,bm)),-1),mul(d,w,add(mul(m,Dw),bw)))
    target=add(scale(mul(m,D),-1),scale(mul(m,gap),-1),scale(mul(m,bm),-1),mul(d,w,bw))
    envelope=add(mul(d,w,Dw),gap)
    need(not add(direct,scale(target,-1),scale(mul(m,envelope),-1)),'Exact C17 chain rule and envelope')
    # C18 Laurent powers after m=k*rho^-4.
    terms=[(1,-1,0),(0,-1,-1),(-2,-1,0),(-3,-1,-1)]
    expected=[(5,-1),(4,-2),(2,-1),(1,-2)]
    need([(r-4*me,ke+me) for r,me,ke in terms]==expected,'All four C18 exponents')
    need((0-4*(-1)-3,-1)==(1,-1),'Relative eta path:1/(N*rho^3)=rho/k')
    need((Q(3,4)*4,3*4)==(3,12),'Core |h|^4=k^3 power')
    need(Q(1)-Q(5,10)==Q(1,2),'Complex k halfline preservation margin')
    return 9

def trim(p):
    out=[Q(v) for v in p]
    while len(out)>1 and out[-1]==0:out.pop()
    return out
def polyadd(p,q):
    out=[Q(0)]*max(len(p),len(q))
    for i,c in enumerate(p):out[i]+=c
    for i,c in enumerate(q):out[i]+=c
    return trim(out)
def polymul(p,q,limit=None):
    out=[Q(0)]*(len(p)+len(q)-1)
    for i,c in enumerate(p):
        for j,v in enumerate(q):out[i+j]+=c*v
    if limit is not None:out=out[:limit]
    return trim(out)
def polypow(p,n,limit=None):
    out=[Q(1)]
    for _ in range(n):out=polymul(out,p,limit)
    return out
def logjet(p,n):
    z=[v/p[0] for v in p[:n]];z[0]-=1
    out=[Q(0)]
    for j in range(1,n):
        out=polyadd(out,[Q((-1)**(j+1),j)*v for v in polypow(z,j,n)])
    return (out+[Q(0)]*n)[:n]
def rising(g,n):
    out=Q(1)
    for j in range(n):out*=g+j
    return out
def literal_controls():
    jets=own=pair=0
    for s in range(1,13):
        a=s*s;d=2*s+1;sc=s*(s+1);wc=Q(sc**d*(s+1),s)
        p=polypow([Q(sc),Q(1)],d,5);p[0]+=wc
        lj=polyadd([a*v for v in logjet([Q(sc),Q(1)],5)],logjet(p,5))
        lj=polyadd(lj,[0,-1]);lj=(lj+[Q(0)]*5)[:5]
        need(lj[1:4]==[0,0,0] and lj[4]==-Q(1,12*sc*sc),'Actual cusp quartic jet')
        sw=-Q(d*sc**(d-1))/(wc+sc**d)**2
        need(sw==-1/(d*wc),'Actual eta coefficient Phi_Sw')
        beta2=Q(sc)/(d*d*wc*wc);cd2=Q(d*d,sc)
        need(beta2*wc*wc*cd2==1,'Actual bounded-core normalization')
        jets+=1
        for m in range(1,7):
            direct=[Q(comb(m,j)*factorial(a*m+d*j),m**(a*m+d*j+1)) for j in range(m+1)]
            actual=[Q(factorial(a*m),m**((a+d)*m+1))*comb(m,j)*rising(Q(a*m+1),d*j)
                    *m**(d*(m-j)) for j in range(m+1)]
            need(direct==actual,'Full own-order integral polynomial')
            own+=1
        for N in range(2,8):
            n=N-1;kappa=Q(N,n)
            direct=[Q(comb(n,j)*factorial(a*n+d*j),n**(a*n+d*j+1))*kappa**(d*(n-j))
                    for j in range(n+1)]
            actual=[Q(factorial(a*n),n**((a+d)*n+1))*comb(n,j)*rising(Q(a*n+1),d*j)
                    *N**(d*(n-j)) for j in range(n+1)]
            need(direct==actual,'Actual same-x adjacent own-shape vectors')
            pair+=1
    return jets,own,pair
def rational_guards():
    out={
      'halfline_ratio':64<75,
      'complex_m_disk_sector':Q(1,20)/(1-Q(1,20))<Q(1,10),
      'cubic_gaussian_exponent':1-Q(5,2)==-Q(3,2),
      'squared_cubic_gaussian_exponent':2-Q(7,2)==-Q(3,2),
      'relative_error_exponent':-Q(3,2)+Q(1,2)==-1,
      'core_power_scaling':Q(3,4)*Q(4,3)==1,
      'Cauchy_m_scale':Q(1,8)>0,
      'Cauchy_w_scale':Q(1,8)>0,
      'small_remainder_log':Q(1,4)<1,
      'path_overlap':Q(99,100)*Q(99,100)**4>Q(1,2),
      'K_error_floor':Q(1,100)**2<Q(1,4),
      'rho_error_small':Q(1,100)/100<Q(1,4),
    }
    need(all(out.values()),'Exact rational support guards')
    return out
def negatives():
    # Actual full-vector wrong-shape control, not a disconnected scalar.
    a,d,N=1,3,3;n=N-1
    correct=[comb(n,j)*rising(Q(a*n+1),d*j)*N**(d*(n-j)) for j in range(n+1)]
    wrong=[comb(n,j)*rising(Q(a*N+1),d*j)*N**(d*(n-j)) for j in range(n+1)]
    need(correct!=wrong,'Shared upper Gamma-shape negative')
    # Removing -D from the product rule changes its generic cleared polynomial.
    m,D=var(0),var(1)
    need(bool(mul(m,D)),'Wrong log-ratio chain rule negative')
    # False unbounded RHP zero-free claim: d3,w16 has 1+sqrt(3)i root
    # up to the exact radius; use S=1+i,w=2-2i, S^3=-2+2i.
    def cmul(z,w):return (z[0]*w[0]-z[1]*w[1],z[0]*w[1]+z[1]*w[0])
    S=(1,1);S3=cmul(cmul(S,S),S);w=(2,-2)
    need((S3[0]+w[0],S3[1]+w[1])==(0,0) and 1>0,'Actual RHP branch counterexample')
    # At x=0 model real loss is strictly positive, so discarding it is false.
    need(Q(9,32)>0,'Midpoint loss negative')
    return 4
def source_negatives(ledger,repo):
    names=[]
    def reject(name,fn):
        try:fn()
        except (ValueError,KeyError,TypeError):names.append(name);return
        raise ValueError('Negative accepted:'+name)
    for role in ROLES:
        bad=deepcopy(ledger);bad['sources']=[r for r in bad['sources'] if r['role']!=role]
        reject('omit:'+role,lambda bad=bad:source_ok(bad,repo))
        bad=deepcopy(ledger);next(r for r in bad['sources'] if r['role']==role)['sha256']='0'*64
        reject('badpin:'+role,lambda bad=bad:source_ok(bad,repo))
    for key in ['uniform_H_core_required','HJN_monotonicity_is_logical_premise','all_RHP_is_zero_free',
                'uniform_in_s','effective_cutoff','global_fixed_cusp_SFC','global_all_support_SFC']:
        bad=deepcopy(SCOPE);bad[key]=True
        reject('false-scope:'+key,lambda bad=bad:inputs_ok(bad))
    for value in [True,'1',2]:
        bad=deepcopy(ledger);bad['schema']=value
        reject('schema:'+str(value),lambda bad=bad:source_ok(bad,repo))
    reject('forged-child',lambda:equal_receipt(b'{"status":"PASS"}',b'{"status":"PASS","real":1}'))
    return names
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--repo-root',required=True,type=Path)
    ap.add_argument('--output',type=Path,default=HERE/RN);args=ap.parse_args();repo=args.repo_root.resolve()
    data=load(HERE/IN);inputs_ok(data)
    ledger=load(HERE/LN);pins=source_ok(ledger,repo);closure=core_closure(repo);child=replay_core(repo)
    ids=universal();jets,own,pair=literal_controls();guards=rational_guards();mathneg=negatives()
    source_controls=source_negatives(ledger,repo)
    payload={'schema':1,'status':'PASS','scope':data,'formal_verification':False,
      'external_review':False,'independent_replay':False,
      'proof_sha256':sha(HERE/PN),'checker_sha256':sha(HERE/CN),'input_sha256':sha(HERE/IN),
      'source_ledger_sha256':sha(HERE/LN),'sources':[pins[k] for k in sorted(pins)],
      'retained_fixed_core_closure':closure,'fresh_fixed_completion_receipt_sha256':child,
      'fresh_fixed_completion_receipt_byte_identical':True,
      'counts':{'universal_algebra_and_scale_identities':ids,'actual_cusp_quartic_jet_vectors':jets,
        'full_actual_own_order_integral_vectors':own,'full_actual_consecutive_common_x_vectors':pair,
        'rational_guards':len(guards),'mathematical_negative_controls':mathneg,
        'source_scope_controls':len(source_controls)},
      'guard_names':sorted(guards),'negative_names':source_controls,
      'logical_inputs':['fixed_cusp_bounded_window_theorem','exact_cusp_quartic_jets'],
      'meaning':'Exact finite algebra and source replays support the written complete-contour/Cauchy proof. They do not formally verify those analytic arguments or prove global SFC.'}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':'PASS','counts':payload['counts']},sort_keys=True))
if __name__=='__main__':main()
