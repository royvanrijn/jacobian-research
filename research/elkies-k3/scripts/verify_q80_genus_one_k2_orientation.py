#!/usr/bin/env python3
"""Independent finite replay; the unbounded nodal orientation proof is written mathematics."""
import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import resource
import time

ROOT=Path(__file__).resolve().parents[2]
DEFAULT=ROOT/'artifacts/generated-results/elkies-k3-q80-genus-one-k2-orientation-v1'
SOURCE='artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json'
BOUNDARY='artifacts/generated-results/elkies-k3-q80-genus-one-two-disagreement-v1'
DEGREE='artifacts/generated-results/elkies-k3-q80-degree-two-reciprocity-v1'
PRODUCER='elkies-k3/scripts/certify_q80_genus_one_k2_orientation.sage'
TEST='tests/test_q80_genus_one_k2_orientation.py'
P=131
def require(condition,message):
    if not condition:raise ValueError(message)
def read(path):return json.loads(path.read_text())
def digest(path):return sha256(path.read_bytes()).hexdigest()
def trim(v):
    v=list(v)
    while v and not v[-1]:v.pop()
    return v
def coeff(values):
    return trim([Fraction(v).numerator*pow(Fraction(v).denominator,-1,P)%P for v in values])
def fa(a,b):return ((a[0]+b[0])%P,(a[1]+b[1])%P)
def fs(a,c):return (a[0]*c%P,a[1]*c%P)
def fm(a,b):
    return ((a[0]*b[0]-88*a[1]*b[1])%P,
            (a[0]*b[1]+a[1]*b[0]-62*a[1]*b[1])%P)
def fn(a):return (a[0]*a[0]-62*a[0]*a[1]+88*a[1]*a[1])%P
def fi(a):
    require(a!=(0,0),'nonzero field element')
    return fs(((a[0]-62*a[1])%P,-a[1]%P),pow(fn(a),-1,P))
def fp(a,n):
    r=(1,0)
    while n:
        if n&1:r=fm(r,a)
        a=fm(a,a);n>>=1
    return r
def fe(a,t=(0,1)):
    r=(0,0)
    for v in reversed(a):r=fa(fm(r,t),(v,0))
    return r
def deriv(a):return [i*a[i]%P for i in range(1,len(a))]
def check_nodal(record,source):
    require(all((t*t+62*t+88)%P for t in range(P)),'irreducible nodal quadratic')
    model=source['weierstrass_model']
    A=coeff(model['A_coefficients_low_to_high']);B=coeff(model['B_coefficients_low_to_high'])
    aa,bb=fe(A),fe(B);ap,bp=fe(deriv(A)),fe(deriv(B))
    node=(29,128);rho=(73,6);derivative=fa(fs(fm(rho,rho),3),aa)
    rho_prime=fs(fm(fa(fm(ap,rho),bp),fi(derivative)),-1)
    vprime=fa(fs(ap,-1),fs(fm(rho,rho_prime),-3*pow(2,-1,P)))
    delta=fa(fs(fp(aa,3),4),fs(fp(bb,2),27))
    delta_prime=fa(fs(fm(fp(aa,2),ap),12),fs(fm(bb,bp),54))
    require(rho==fs(node,-2) and rho!=(0,0) and derivative!=(0,0),'separated simple and nonzero double cubic roots')
    require(fa(fa(fp(rho,3),fm(aa,rho)),bb)==(0,0),'simple cubic root')
    require(fa(fs(fp(node,2),3),aa)==(0,0),'double cubic root')
    require(delta==(0,0) and delta_prime!=(0,0),'simple discriminant root')
    require(fa(fs(aa,-1),fs(fp(rho,2),-3*pow(4,-1,P)))==(0,0),'V vanishes')
    require(vprime!=(0,0) and fs(delta_prime,-1)==fs(fm(vprime,fp(derivative,2)),4),'actual implicit derivative and simple V zero')
    codes=[]
    for root in [node,rho]:
        code=0
        for j,section in enumerate(source['sections']['records']):
            x=section['X'];num=fe(coeff(x['numerator_coefficients_low_to_high']))
            den=fe(coeff(x['denominator_coefficients_low_to_high']))
            require(den!=(0,0),'literal basis regular at nodal root')
            value=fa(fm(num,fi(den)),fs(root,-1))
            require(value!=(0,0),'no unhandled regularized root value')
            norm=fn(value);require(norm!=0,'nonzero evaluation norm')
            code |= int(pow(norm,(P-1)//2,P)==P-1)<<j
        codes.append(code)
    require(len(source['sections']['records'])==17 and codes==[13412,0],'all seventeen inherited nodal norm characters')
    expected={'q':[88,62,1],'A_mod_p':A,'B_mod_p':B,'node':list(node),'rho':list(rho),
              'rho_prime':list(rho_prime),'simple_root_derivative':list(derivative),
              'V_prime':list(vprime),'delta_prime':list(delta_prime),'norm_codes_node_simple':codes}
    require(record==expected,'exact nodal first-derivative/character record')
    return expected

# Sparse integer polynomial arithmetic: exponent tuples -> coefficients.
def sa(a,b):
    c=dict(a)
    for e,v in b.items():c[e]=c.get(e,0)+v
    return {e:v for e,v in c.items() if v}
def ss(a,k):return {e:v*k for e,v in a.items() if v*k}
def sm(a,b):
    c={}
    for e,v in a.items():
        for f,w in b.items():
            g=tuple(x+y for x,y in zip(e,f));c[g]=c.get(g,0)+v*w
    return {e:v for e,v in c.items() if v}
def sp(a,n,dimension):
    r={(0,)*dimension:1}
    while n:
        if n&1:r=sm(r,a)
        a=sm(a,a);n>>=1
    return r
def sc(c,n):return {(0,)*n:c} if c else {}
def sv(k,n):return {tuple(int(i==k) for i in range(n)):1}
def check_symbolic(record):
    r,a,x=sv(0,3),sv(1,3),sv(2,3);power=lambda u,n:sp(u,n,3)
    b=ss(sa(power(r,3),sm(a,r)),-1)
    cubic=sa(sa(power(x,3),sm(a,x)),b)
    factor=sm(sa(x,ss(r,-1)),sa(sa(power(sa(ss(x,2),r),2),ss(a,4)),ss(power(r,2),3)))
    require(ss(cubic,4)==factor,'symbolic cubic factorization')
    disc=sa(ss(power(a,3),-4),ss(power(b,2),-27))
    factordisc=sm(sa(ss(a,-4),ss(power(r,2),-3)),power(sa(ss(power(r,2),3),a),2))
    require(disc==factordisc,'symbolic nodal discriminant identity')
    t,z=sv(0,2),sv(1,2);power=lambda u,n:sp(u,n,2)
    c=sa(power(z,2),sc(-3,2));s=ss(sm(z,c),4);beta=power(c,2)
    gamma=ss(sm(c,sa(power(z,2),sc(1,2))),-3)
    tb=sa(t,ss(beta,-1));tg=sa(t,ss(gamma,-1))
    Q=sa(sa(power(tb,2),ss(sm(sm(s,z),tb),4)),ss(sm(power(s,2),sa(power(z,2),sc(-1,2))),6))
    H=sa(power(sa(power(tg,2),ss(power(s,2),-3)),2),ss(sm(power(s,4),t),-1))
    require(H==sm(power(tb,2),Q),'symbolic local section identity')
    terms=[]
    for k in range(3):
        rows={e[1]:v for e,v in Q.items() if e[0]==k}
        terms.append(trim([rows.get(j,0) for j in range(max(rows,default=-1)+1)]))
    q0={e:v for e,v in Q.items() if e[0]==0}
    q1={(0,e[1]):v for e,v in Q.items() if e[0]==1}
    require(sa(power(q1,2),ss(q0,-4))==ss(sm(power(s,2),c),-8),'symbolic toy discriminant')
    expected={'cubic_factorization':True,'nodal_discriminant_identity':True,'toy_factorization':True,
              'toy_discriminant_identity':True,'checked_as_polynomial_identities':True,
              'variable_order':['t','z'],'toy_Q_coefficients':[[str(v) for v in row] for row in terms]}
    require(record==expected,'symbolic coefficient record')
    return expected

def v(a):
    a=Fraction(a)
    if not a:return float('inf')
    def vp(n):
        n=abs(n);count=0
        while n%P==0:n//=P;count+=1
        return count
    return vp(a.numerator)-vp(a.denominator)
def peval(coeffs,t):
    y=Fraction(0)
    for a in reversed(coeffs):y=y*t+a
    return y
def sqrt_hensel(unit,digits,seed):
    require(digits>=1 and 0<=seed<P and (seed*seed-unit)%P==0 and seed!=0,'simple square-root seed')
    root=seed;mod=P
    for _ in range(1,digits):
        require((unit-root*root)%mod==0,'exact Hensel lift divisibility')
        root+=((unit-root*root)//mod*pow(2*root,-1,P)%P)*mod
        mod*=P
    require((root*root-unit)%mod==0,'exact lifted square congruence')
    return root
def check_controls(record,packet):
    expected=[]
    for z in packet['control_z']:
        c=z*z-3;s=4*z*c;beta=c*c;gamma=-3*c*(z*z+1)
        q=[beta*beta-4*s*z*beta+6*s*s*(z*z-1),4*s*z-2*beta,1]
        disc=q[1]*q[1]-4*q[0];require(disc==-8*s*s*c,'exact specialized discriminant')
        vd=v(disc);unit=disc//P**vd;legendre=pow(unit%P,(P-1)//2,P)
        square=legendre==1
        row={'z':z,'contact_valuation':v(c),'Q':[str(vv) for vv in q],'disc_valuation':vd,
             'disc_unit_mod_p':unit%P,'disc_unit_square_mod_p':square,
             'branch_field_over_Qp':'ramified_quadratic' if vd%2 else ('split' if square else 'unramified_quadratic')}
        if not vd%2 and square:
            seed=min(i for i in range(1,P) if (i*i-unit)%P==0)
            root=sqrt_hensel(unit,packet['hensel_unit_digits'],seed)
            centers=[Fraction(-q[1]+sign*P**(vd//2)*root,2) for sign in [1,-1]]
            radii=[]
            for b in centers:
                vf=v(peval(q,b));vdq=v(2*b+q[1])
                require(vf>2*vdq,'strict Hensel root ball')
                radii.append(vf-vdq)
            radius=min(radii);require(radius>=19,'retained branch precision')
            x=lambda b:Fraction(-2)+(b-gamma)**2/Fraction(s*s)
            xv=[x(b) for b in centers]
            xradius=min(min(radius+v(2*(b-gamma)/Fraction(s*s)),2*radius-v(s*s)) for b in centers)
            values={'beta':v(beta),'gamma':v(gamma),'branch1':v(centers[0]),'branch2':v(centers[1]),
                    'branch_separation':v(centers[0]-centers[1]),'x_difference':v(xv[0]-xv[1]),
                    'node_offsets_sum':v(xv[0]+xv[1]-2),'node_offset1':v(xv[0]-1)}
            require(values=={'beta':4,'gamma':2,'branch1':2,'branch2':2,'branch_separation':3,
                             'x_difference':1,'node_offsets_sum':2,'node_offset1':1},'opposite node-root orientation and separation')
            require(radius>values['branch_separation'] and xradius>max(values['x_difference'],values['node_offsets_sum']),'valuations certified beyond center precision')
            row.update({'sqrt_unit_mod_p_digits':str(root),'sqrt_unit_seed':seed,'branch_centers':[str(b) for b in centers],
                        'branch_error_valuation_lower_bound':radius,'x_error_valuation_lower_bound':xradius,'valuations':values})
        expected.append(row)
    require([row['branch_field_over_Qp'] for row in expected]==['ramified_quadratic','unramified_quadratic','split'],'three distinct arithmetic control boundaries')
    require(record=={'rows':expected},'exact control and Hensel-ball record')
    return expected
def verify(out):
    start=time.process_time();wall=time.monotonic();packet=read(out/'input.json');result=read(out/'result.json')
    expected={SOURCE,PRODUCER,TEST,str(Path(__file__).relative_to(ROOT)),
              'elkies-k3/scripts/verify_q80_genus_one_k2_nodal_boundary.py',
              *[BOUNDARY+'/'+n for n in ['input.json','result.json','character-gate.json','independent-replay.json']],
              *[DEGREE+'/'+n for n in ['input.json','result.json','quadratic-fibres.json','independent-replay.json']]}
    require(set(packet['bindings'])==expected,'exact source binding set')
    for s,h in {**packet['bindings'],**packet['preserved_preflight']}.items():require(digest(ROOT/s)==h,'source/preflight binding '+s)
    require((packet['schema'],packet['prime'],packet['cpu_seconds'],packet['memory_bytes'])==(1,P,20,2*1024**3),'arithmetic and resource scope')
    require(packet['control_z']==[38,13793,30954] and packet['hensel_unit_digits']==16,'declared finite controls and precision')
    require(result['input_sha256']==digest(out/'input.json') and result['status']=='PASS','bound result')
    require(set(result['records'])=={'nodal-hypotheses.json','symbolic-identities.json','local-controls.json'},'all finite result records')
    for s,h in result['records'].items():require(digest(out/s)==h,'result record '+s)
    old=read(ROOT/BOUNDARY/'result.json');gate=read(ROOT/BOUNDARY/'character-gate.json')
    require(old['status']=='PASS' and old['input_sha256']==digest(ROOT/BOUNDARY/'input.json'),'retained boundary result binding')
    require(old['records']['character-gate.json']==digest(ROOT/BOUNDARY/'character-gate.json'),'retained norm-character evidence')
    require(old['remaining_nodal_denominator_boundary']=='UNKNOWN' and old['both_abscissas_require_negative_gauss_valuation'],'earlier result scope preserved')
    require(old['necessary_D_reduction']==[15,39,90,124,1],'retained q0-squared boundary')
    require(gate['quadratic_zero_rows'][0]['norm_codes']==[13412,13412,0],'retained unused-root character')
    nodal=check_nodal(read(out/'nodal-hypotheses.json'),read(ROOT/SOURCE))
    check_symbolic(read(out/'symbolic-identities.json'))
    controls=check_controls(read(out/'local-controls.json'),packet)
    require(result['genus1_allocations_remaining']==9 and not result['positive_mw17_target_complete']
            and not result['formal_verification'] and not result['old_norm8_replay_upgraded'],'unclosed positive endpoint and proof assurances')
    return {'status':'PASS','input_sha256':digest(out/'input.json'),'result_sha256':digest(out/'result.json'),
            'checker_sha256':digest(Path(__file__)),'finite_checks':{'nodal_norm_codes':nodal['norm_codes_node_simple'],
            'symbolic_identities':4,'controls':len(controls),'branch_fields':[r['branch_field_over_Qp'] for r in controls]},
            'written_analytic_proof_required':True,'earlier_complete_censuses_rerun':False,
            'cpu_seconds':time.process_time()-start,'wall_seconds':time.monotonic()-wall,
            'cpu_limit':20,'memory_limit_bytes':2*1024**3}
if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--directory',type=Path,default=DEFAULT)
    parser.add_argument('--receipt',type=Path);args=parser.parse_args()
    resource.setrlimit(resource.RLIMIT_CPU,(20,25));resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    receipt=verify(args.directory)
    if args.receipt:
        with args.receipt.open('x') as f:json.dump(receipt,f,indent=2,sort_keys=True);f.write('\n')
    print(json.dumps(receipt,sort_keys=True),flush=True)
