#!/usr/bin/env sage-python
"""Three frozen generic rows: RR maps and constructive odd-degree descent.

One worker per selected orbit, capped at25s. No rational point search or
integer factorization. The old8044 equations are reused as a regression.
"""
import argparse,hashlib,json,runpy,signal
from pathlib import Path
from sage.all import QQ,ZZ,PolynomialRing,matrix,vector,gcd
signal.alarm(25)
ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves';CAS=ROOT/'elliptic-curves/cas'
OUT=ART/'det1092_rational_bisection_index_v1'
PARENT=ART/'curve302_recovered_mw17_parent_v1.json'
SELECT=OUT/'index-and-selection.json'
RR=CAS/'construct_curve302_parent_cheapest_lattice_bisection.sage'
LOADER=CAS/'load_curve302_recovered_parent.sage'
OLD=ART/'det1092_orbit8044_rank18_base_change_v2.json'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(name,data):
    p=OUT/name
    if p.exists():assert read(p)==data
    else:
        with p.open('x') as stream:json.dump(data,stream,indent=2,sort_keys=True);stream.write('\n')
args=argparse.ArgumentParser();args.add_argument('orbit',type=int);args=args.parse_args()
selection=read(SELECT);row=next(r for r in selection['selected'] if r['orbit']==args.orbit)
for name,digest in selection['inputs'].items():assert sha(ROOT/name)==digest
save('orbit-%d-protocol.json'%args.orbit,{'classification':'generic-only bounded multisection construction',
    'orbit':args.orbit,'selection_sha256':sha(SELECT),
    'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [PARENT,RR,LOADER,OLD,Path(__file__)]},
    'limits':{'seconds':25,'new_RR_systems':int(args.orbit!=8044),'point_searches':0,'integer_factorizations':0},
    'rational_point_rule':'First odd-intersection generic section; convert its odd divisor to a conic point by polynomial norm descent. No conic qfsolve or exceptional point.'})
E,basis,_=runpy.run_path(str(LOADER))['load_curve302_recovered_parent'](PARENT)
R=E.base_ring().ring();F=E.base_ring();t=R.gen()
rr=runpy.run_path(str(RR));w=vector(ZZ,row['word'])
G=matrix(ZZ,read(PARENT)['generic_height_gram']);assert w*G*w==10
negative_trace=-sum((v*Q for v,Q in zip(w,basis)),E(0))
if args.orbit==8044:
    old=read(OLD)
    def decode(v):return F(R(v['numerator']))/R(v['denominator'])
    c,b,a=[decode(v) for v in old['lift']['residual_coefficients']]
    line=[R(v) for v in old['lift']['line_coefficients']]
    q=R(old['curve_over_Q']['q_coefficients']);h=R(old['splitting']['discriminant_square_factor'])
    matrix_rank='REUSED_OLD_RANK19_CERTIFICATE'
else:
    relation=rr['primitive_kernel_relation'](negative_trace,R)
    residual=rr['residual_quadratic'](E,negative_trace,relation)
    c,b,a=residual.list();line=[relation[k] for k in ['f0','f1','f2']]
    disc=R(b*b-4*a*c);factorization=disc.factor();q=R(factorization.unit());h=R.one()
    for v,e in factorization:q*=v**(e%2);h*=v**(e//2)
    matrix_rank=int(relation['rank'])
assert b*b-4*a*c==h*h*q and q.degree()==2 and q.discriminant()!=0
f0,f1,f2=line
x0=-b/(2*a);x1=h/(2*a)
y0=-(f0+f1*x0)/f2;y1=-f1*x1/f2
coordinates=list(map(R,[x0,x1,y0,y1]))
x0,x1,y0,y1=coordinates
assert gcd(x1,y1).degree()==0
d,A,B=x1.xgcd(y1);A/=d;B/=d
assert A*x1+B*y1==1
index=row['first_odd_section'];section=basis[index]
intersection=row['section_intersection_degrees'][index]
assert intersection>0 and intersection%2
steps=[]
if QQ(q[2]).is_square():
    sqrt=QQ(q[2]).sqrt()
    point={'type':'rational_infinity','t_over_w':str(1/sqrt),'leading_sqrt':str(sqrt)}
    U=PolynomialRing(QQ,'u');K=U.fraction_field();u=U.gen()
    T=K(q[0]-u*u)/(2*sqrt*u-q[1]);W=sqrt*T+u
else:
    sr=F(A)*(section[0]-x0)+F(B)*(section[1]-y0)
    errors=[section[0]-x0-x1*sr,section[1]-y0-y1*sr,sr*sr-q]
    divisor=gcd([R(v.numerator()) for v in errors if v]).monic()
    forbidden=R(sr.denominator())*R(section[0].denominator())*R(section[1].denominator())
    while True:
        common=gcd(divisor,forbidden)
        if common.degree()==0:break
        divisor=(divisor//common).monic()
    assert divisor.degree()==intersection and divisor.degree()%2
    ordinate=(R(sr.numerator())*R(sr.denominator()).inverse_mod(divisor))%divisor
    assert (ordinate*ordinate-q)%divisor==0
    starting_divisor=divisor;starting_ordinate=ordinate
    while divisor.degree()>1:
        quotient,remainder=(ordinate*ordinate-q).quo_rem(divisor)
        assert remainder==0 and 0<quotient.degree()<divisor.degree() and quotient.degree()%2
        steps.append({'divisor':list(map(str,divisor.list())),'ordinate':list(map(str,ordinate.list())),
                      'quotient':list(map(str,quotient.list()))})
        divisor=quotient.monic();ordinate%=divisor
    t0=-divisor[0]/divisor[1];s0=ordinate(t0)
    assert s0*s0==q(t0)
    point={'type':'finite','t':str(t0),'s':str(s0),
           'starting_divisor':list(map(str,starting_divisor.list())),
           'starting_ordinate':list(map(str,starting_ordinate.list()))}
    U=PolynomialRing(QQ,'u');K=U.fraction_field();u=U.gen()
    T=K(t0)+(q.derivative()(t0)-2*s0*u)/(u*u-q[2]);W=K(s0)+u*(T-t0)
assert W*W==q(T) and max(T.numerator().degree(),T.denominator().degree())==2
def rec(v):return {'numerator':list(map(str,v.numerator().list())),
                   'denominator':list(map(str,v.denominator().list()))}
output={'status':'PASS_Q_RATIONAL_BISECTION_BY_ODD_DIVISOR_DESCENT',
    'classification':'generic-only exact construction; independent replay pending',
    'orbit':args.orbit,'word':row['word'],'RR_rank':matrix_rank,
    'line_coefficients':[list(map(str,v.list())) for v in line],
    'negative_trace':[rec(v) for v in negative_trace.xy()],
    'residual_coefficients':[rec(v) for v in [c,b,a]],
    'q':list(map(str,q.list())),'h':list(map(str,h.list())),
    'elliptic_quadratic_maps':[list(map(str,v.list())) for v in coordinates],
    'odd_section_index':index,'intersection_degree':intersection,
    'rational_point':point,'norm_descent_steps':steps,
    'base_map':rec(T),'conic_ordinate':rec(W),
    'function_field_rank_lower_bound':18,
    'independence':'The completed all-prime theorem makes every genuine geometrically irreducible bisection independent over the inherited generic span.',
    'boundary':'No specialization outcome was used to construct or parametrize this curve.',
    'inputs':read(OUT/('orbit-%d-protocol.json'%args.orbit))['inputs'],
    'protocol_sha256':sha(OUT/('orbit-%d-protocol.json'%args.orbit))}
save('orbit-%d.json'%args.orbit,output)
print(output['status'],args.orbit,'odd degree',intersection,'descent steps',len(steps),flush=True)
