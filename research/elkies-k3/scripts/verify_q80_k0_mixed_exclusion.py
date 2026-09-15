#!/usr/bin/env python3
"""Independent integer-polynomial and unramified-ring replay of Q80 mixed-case premises."""
from fractions import Fraction
from hashlib import sha256
import argparse,importlib.util,json
from pathlib import Path
import resource,time
ROOT=Path(__file__).resolve().parents[2]
DEFAULT=ROOT/'artifacts/generated-results/elkies-k3-q80-genus-one-k0-orthogonal-contact-v1'
SOURCE='artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json'
NORM4='artifacts/generated-results/elkies-k3-q80-genus-zero-integral-reduction-v1'
NODE='artifacts/generated-results/elkies-k3-q80-genus-one-k2-orientation-v1'
POLE='artifacts/generated-results/elkies-k3-q80-genus-one-k0-pole-boundary-v1'
PRODUCER='elkies-k3/scripts/certify_q80_k0_mixed_exclusion.sage'
HELPER='elkies-k3/scripts/verify_q80_genus_one_k2_nodal_boundary.py'
P=131
s=importlib.util.spec_from_file_location('q80_mixed_polynomial_primitives',ROOT/HELPER);h=importlib.util.module_from_spec(s);s.loader.exec_module(h)
read=h.read;digest=h.digest;require=h.require;add=h.add;sub=h.sub;mul=h.mul;scale=h.scale;exact=h.exact;divrem=h.divrem;rank=h.rank;ev=h.ev
def det(a):
    a=[[v%P for v in r] for r in a];d=1
    for j in range(len(a)):
        k=next((i for i in range(j,len(a)) if a[i][j]),None)
        if k is None:return 0
        if j!=k:a[j],a[k]=a[k],a[j];d=-d
        d=d*a[j][j]%P;u=pow(a[j][j],-1,P)
        for i in range(j+1,len(a)):
            c=a[i][j]*u%P
            for k in range(j,len(a)):a[i][k]=(a[i][k]-c*a[j][k])%P
    return d%P
class Unramified:
    """(Z/modulus)[theta]/(theta^2+62theta+88), exact unit inverses only."""
    def __init__(self,modulus):self.m=modulus
    def add(self,a,b):return ((a[0]+b[0])%self.m,(a[1]+b[1])%self.m)
    def scale(self,a,c):return (a[0]*c%self.m,a[1]*c%self.m)
    def sub(self,a,b):return self.add(a,self.scale(b,-1))
    def mul(self,a,b):return ((a[0]*b[0]-88*a[1]*b[1])%self.m,(a[0]*b[1]+a[1]*b[0]-62*a[1]*b[1])%self.m)
    def inverse(self,a):
        n=(a[0]*a[0]-62*a[0]*a[1]+88*a[1]*a[1])%self.m
        require(n%P!=0,'unit inverse in unramified ring')
        return self.scale(((a[0]-62*a[1])%self.m,-a[1]%self.m),pow(n,-1,self.m))
    def div(self,a,b):return self.mul(a,self.inverse(b))
    def power(self,a,n):
        r=(1,0)
        while n:
            if n&1:r=self.mul(r,a)
            a=self.mul(a,a);n>>=1
        return r
    def eval(self,a,t):
        r=(0,0)
        for c in reversed(a):r=self.add(self.mul(r,t),(c,0))
        return r
def derivative(a,m):return [i*a[i]%m for i in range(1,len(a))]
def check_constant(record,A,B,point):
    q=[88,62,1];D=mul(q,q);x0=record['x0']
    require(len(x0)<=4 and divrem(x0,q)[1]==[73,6],'unique simple-root residue lift chart')
    require(not divrem(add(add(mul(mul(x0,x0),x0),mul(A,x0)),B),D)[1],'lift to q0 squared')
    require(len(h.gcd(add(scale(mul(x0,x0),3),A),q))==1,'simple-root Hensel uniqueness')
    accepted=[]
    for c in range(P):
        x=add(x0,scale(D,c));f=exact(add(add(mul(mul(x,x),x),mul(A,x)),B),D)
        if h.square_up_to_scalar(f):accepted.append((c,x,f))
    require(record['tried']==131 and record['q0']==q and len(accepted)==len(record['rows'])==1,'complete constant-twist contact census')
    row=record['rows'][0];c,x,f=accepted[0]
    require(c==row['c']==82 and x==row['x']==point['x'],'unique integral abscissa')
    require(f==scale(mul(row['r_monic'],row['r_monic']),row['scalar']) and row['r_monic'][-1]==1,'exact scalar-square quotient')
    require(row['scalar']==9 and row['scalar_square'] and pow(row['scalar'],65,P)==1,'no nonsquare twist contact')
def check_deformation(record,A,B,point,source):
    x,y=point['x'],point['y'];q=[88,62,1];D=mul(q,q);r=exact(y,q);fp=add(scale(mul(x,x),3),A)
    require(add(add(mul(mul(x,x),x),mul(A,x)),B)==mul(y,y),'exact reduced seed section')
    require((record['point_index'],record['word'],record['x'],record['y'],record['r'])==(732,point['word'],x,y,r),'exact deformation seed')
    w=point['word'];G=source['sections']['height_gram'];height=sum(Fraction(G[i][j])*w[i]*w[j] for i in range(17) for j in range(17))
    require(height==record['word_height']==4 and record['word_identity_over_F131_t'],'norm-four inherited word premise')
    columns=[[0]*k+fp for k in range(5)]+[[0]*k+scale(mul(D,r),-2) for k in range(5)]+[[0]*k+scale(mul(r,r),-1) for k in range(4)]
    J=[[columns[j][i] if i<len(columns[j]) else 0 for j in range(14)] for i in range(13)]
    minor=[[v for j,v in enumerate(row) if j!=4] for row in J];d=det(minor)
    require(rank(J)==record['jacobian_rank']==13 and d==record['fixed_x4_minor_determinant'] and d!=0,'unit formal implicit-function minor')
    tangent=record['normalized_tangent'];require(tangent['dx']==r and r[4]==3,'normalized one-parameter tangent')
    require(sub(sub(fp,scale(mul(D,tangent['dr']),2)),mul(r,tangent['dD']))==[] and len(tangent['dD'])<=4,'exact tangent divided identity')
    for key,f in [('r_mod_q0',r),('fprime_mod_q0',fp),('dD_mod_q0',tangent['dD'])]:
        require(record[key]==divrem(f,q)[1] and record[key],'unit branch-discriminant derivative '+key)
def check_nodal(record,source,point):
    R=Unramified(P*P);F=Unramified(P);model=source['weierstrass_model']
    A=h.fraction_poly(model['A_coefficients_low_to_high'],P*P);B=h.fraction_poly(model['B_coefficients_low_to_high'],P*P)
    theta=(0,1);aa=R.eval(A,theta);bb=R.eval(B,theta)
    delta=R.add(R.scale(R.power(aa,3),4),R.scale(R.power(bb,2),27))
    dp=R.add(R.scale(R.mul(R.power(aa,2),R.eval(derivative(A,P*P),theta)),12),R.scale(R.mul(bb,R.eval(derivative(B,P*P),theta)),54))
    alpha=R.sub(theta,R.div(delta,dp));aa=R.eval(A,alpha);bb=R.eval(B,alpha);e=R.div(R.scale(bb,-3),R.scale(aa,2))
    require(aa==R.scale(R.power(e,2),-3) and bb==R.scale(R.power(e,3),2),'exact nodal fibre modulo131 squared')
    require(list(alpha)==record['alpha_mod_p2'] and list(e)==record['node_mod_p2'],'independent discriminant Hensel lift')
    accum_a=(1,0);accum_b=(0,0);rows=[]
    for i,w in enumerate(point['word']):
        if not w:continue
        require(w in [-1,1],'literal short seed word');s=source['sections']['records'][i];values={}
        for key in ['X','Y']:
            n=h.fraction_poly(s[key]['numerator_coefficients_low_to_high'],P*P);d=h.fraction_poly(s[key]['denominator_coefficients_low_to_high'],P*P)
            values[key]=R.div(R.eval(n,alpha),R.eval(d,alpha))
        x,y=values['X'],values['Y'];require(R.power(y,2)==R.add(R.add(R.power(x,3),R.mul(aa,x)),bb),'literal basis point at the nodal fibre')
        a=R.scale(y,w);b=R.sub(x,e)
        # Projective coth addition; no square root of3e and no p-adic floating precision.
        na=R.add(R.mul(accum_a,a),R.scale(R.mul(e,R.mul(accum_b,b)),3))
        nb=R.add(R.mul(accum_a,b),R.mul(accum_b,a));accum_a,accum_b=na,nb
        rows.append({'index':i,'word':w,'x_mod_p2':list(x),'y_mod_p2':list(y)})
    u=R.div(accum_a,accum_b);yt=R.mul(u,R.sub(R.power(u,2),R.scale(e,3)))
    require(rows==record['basis_rows'] and list(u)==record['T_u_mod_p2'] and list(yt)==record['T_y_mod_p2'],'independent projective nodal group word')
    require(all(c%P==0 for c in u+yt) and any(c for c in u) and any(c for c in yt),'exact first-order contact, not an approximate zero')
    yover=tuple(c//P for c in yt);yp=F.eval(derivative(point['y'],P),theta);g=F.scale(F.div(yover,yp),-1);qp=(62,2)
    required=F.div(F.scale((29,128),12),F.mul(F.power(qp,2),F.power(g,2)))
    require(list(yover)==record['T_y_over_p_mod_q0']==[63,6] and list(yp)==record['T_yprime_mod_q0']==[21,32],'exact contact numerator and derivative')
    require(list(g)==record['contact_shift_over_p']==[91,113],'exact contact displacement')
    require(list(required)==record['required_primitive_scalar']==[69,40] and F.power(required,P)!=required,'necessary scalar lies outside F131')
    require(not record['required_scalar_in_F131'] and record['T_y_valuation']==1 and record['T_x_minus_simple_valuation']==2 and record['necessary_negative_abscissa_valuation']==-2,'stated valuation controls')
    return {'u_over_p':[c//P for c in u],'y_over_p':list(yover),'contact_shift':list(g),'required_scalar':list(required)}
def verify(out):
    started=time.process_time();packet=read(out/'input.json');result=read(out/'result.json')
    expected={SOURCE,PRODUCER,HELPER,str(Path(__file__).relative_to(ROOT)),
       *[NORM4+'/'+n for n in ['input.json','result.json','norm4-sections.json','complete-replay-input.json','complete-independent-replay.json']],
       *[NODE+'/'+n for n in ['input.json','result.json','nodal-hypotheses.json','independent-replay.json']],
       *[POLE+'/'+n for n in ['input.json','receipt.json']],
       'elkies-k3/Q80_GENUS_ONE_K0_POLE_BOUNDARY_2026-09-14.md'}
    require(set(packet['bindings'])==expected and packet['prime']==131 and packet['producer_padic_precision']==20,'exact finite input scope')
    for path,hsh in {**packet['bindings'],**packet['preserved_preflight']}.items():require(digest(ROOT/path)==hsh,'source/preflight binding '+path)
    require(result['status']=='PASS' and result['input_sha256']==digest(out/'input.json'),'bound result')
    require(set(result['records'])=={'constant-contacts.json','formal-deformation.json','nodal-contact.json'},'complete finite records')
    for path,hsh in result['records'].items():require(digest(out/path)==hsh,'exact finite record '+path)
    source=read(ROOT/SOURCE);point=read(ROOT/NORM4/'norm4-sections.json')['records'][732];model=source['weierstrass_model']
    A=h.fraction_poly(model['A_coefficients_low_to_high'],P);B=h.fraction_poly(model['B_coefficients_low_to_high'],P)
    check_constant(read(out/'constant-contacts.json'),A,B,point);check_deformation(read(out/'formal-deformation.json'),A,B,point,source)
    finite=check_nodal(read(out/'nodal-contact.json'),source,point)
    require(result['constant_twist_contacts']==1 and result['nonsquare_contacts']==0 and result['mixed_valuation_excluded_by_written_proof'],'finite versus written conclusion')
    require(result['remaining_genus1_k0_j']==[2] and result['remaining_valuation_pattern']=='both negative' and result['remaining_boundary']=='UNKNOWN','exact surviving case')
    require(result['written_deformation_and_valuation_argument_required'] and not result['new_analytic_argument_independently_verified'] and not result['formal_verification'] and not result['external_review'] and not result['positive_mw17_target_complete'],'analytic assurance and positive endpoint limits')
    return {'status':'PASS','input_sha256':digest(out/'input.json'),'result_sha256':digest(out/'result.json'),'checker_sha256':digest(Path(__file__)),
            'constant_contact_coefficients_checked':131,'constant_contacts':1,'nonsquare_contacts':0,'formal_jacobian_rank':13,
            'exact_unramified_modulus':P*P,'independent_nodal_word':finite,
            'independent_finite_replay':True,'analytic_proof_independently_verified':False,'positive_mw17_target_complete':False,'cpu_seconds':time.process_time()-started}
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--directory',type=Path,default=DEFAULT);p.add_argument('--receipt',type=Path);args=p.parse_args()
    resource.setrlimit(resource.RLIMIT_CPU,(30,35));resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    record=verify(args.directory)
    if args.receipt:
        with args.receipt.open('x') as f:json.dump(record,f,indent=2,sort_keys=True);f.write('\n')
    print(json.dumps(record),flush=True)
