#!/usr/bin/env sage-python
"""Transport the exact14 historical seed directions through the class1 pencil.

Fixed input representatives; no translations, new points or parameter search.
Per-case transport and good-prime certificates are checkpointed. A finite miss
leaves membership in the specialized generic MW17 span UNKNOWN.
"""
import argparse,hashlib,json,sys
from pathlib import Path
from sage.all import QQ,ZZ,GF,EllipticCurve,PolynomialRing,matrix,vector,pari,prime_range
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
WORK=ROOT/'artifacts/local/elliptic-curves/curve302-class1-point-transport-v1'
OUTPUT=ART/'curve302_class1_point_transport_v1.json'
PROTOCOL=ART/'curve302_class1_point_transport_protocol_v1.json'
sys.path.insert(0,str(Path(__file__).parent))
import icarm_curve302 as curve
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,d):
    if p.exists():assert read(p)==d,p
    else:
        p.parent.mkdir(parents=True,exist_ok=True)
        with p.open('x')as out:json.dump(d,out,indent=2,sort_keys=True);out.write('\n')
def polynomial(coeffs,t):
    ans=QQ(0)
    for x in reversed(coeffs):ans=ans*t+QQ(x)
    return ans
def dec(record,t):return polynomial(record['numerator'],t)/polynomial(record['denominator'],t)
def jpoint(P):return list(map(str,P.xy()))if P else None
def freeze():
    names=['curve302_recovered_mw17_parent_v1.json','curve302_residual_visibility_geometry_v1.json']
    names += ['x1092_class1_realization_'+n+'_v1.json'for n in ['marking','trace','rr','equation','sections','parent','compact_parent']]
    paths=[Path(__file__),Path(curve.__file__),*[ART/n for n in names]]
    parent=read(ART/names[0]);visibility=read(ART/names[1]);M=matrix(ZZ,parent['basis_embedding_in_public_D'])
    rows=[{'id':r['id'],'public_word':r['public_word'],'strict_mod_2':r['strict_mod_2']}for r in visibility['directions']]
    C=matrix(ZZ,[r['public_word']for r in rows]).transpose();det=M.augment(C).det()
    assert len(rows)==14 and abs(det)==1
    policy={'schema':'curve302.class1-point-transport.protocol.v1',
            'bindings':{str(p.relative_to(ROOT)):sha(p)for p in paths},'rows':rows,
            'integral_basis_determinant':int(det),'source_parameter':'0',
            'limits':{'workers':1,'wall_seconds_per_case':60,'good_prime_bound':2000,'points':14,
                      'new_points':0,'new_parameter_searches':0,'translations_tried':0},
            'scope':'The same fourteen representatives used in the seed-universality panel, forming an integral quotient basis. Exact surface transport, then finite tests against all17 specialized generic B-sections.',
            'boundary':'The surface map is not a homomorphism on the original302 fibre. Each B-base parameter is retained; a finite membership miss is UNKNOWN.'}
    save(PROTOCOL,policy);print('FROZEN14 quotient determinant',det,flush=True)
def build_case(index):
    policy=read(PROTOCOL)
    for name,digest in policy['bindings'].items():assert sha(ROOT/name)==digest,name
    row=policy['rows'][index];case=WORK/row['id'];case.mkdir(parents=True,exist_ok=True)
    parent=read(ART/'curve302_recovered_mw17_parent_v1.json')
    prefix='x1092_class1_realization_'
    tr,rr,eq,normalized,compact=[read(ART/(prefix+n+'_v1.json'))for n in ['trace','rr','equation','parent','compact_parent']]
    E=EllipticCurve(QQ,list(map(QQ,curve.GENERAL_WEIERSTRASS_COEFFICIENTS)))
    public=[E(list(map(QQ,p)))for p in curve.POINTS];t=QQ(0)
    assert [dec(a,t)for a in parent['a_invariants']]==list(E.a_invariants())
    M=matrix(ZZ,parent['basis_embedding_in_public_D'])
    genericA=[E([dec(c,t)for c in p])for p in parent['basis_weierstrass_coordinates']]
    for i,P in enumerate(genericA):
        assert P==sum((int(v)*Q for v,Q in zip(M.column(i),public)),E(0))
    P=sum((int(v)*Q for v,Q in zip(row['public_word'],public)),E(0))
    assert P
    # Exact completion of square and cube, without integral rescaling.
    x,y=P.xy();xs=x+E.b2()/12;ys=y+(E.a1()*x+E.a3())/2
    A,B=[polynomial(tr[k],t)for k in ['A','B']]
    assert A==-E.c4()/48 and B==-E.c6()/864 and ys**2==xs**3+A*xs+B
    h,nx,ny=[polynomial(tr[k],t)for k in ['h','Nx','Ny']]
    a0,b0,a1,b1=[polynomial(rr[k],t)for k in ['a0','b0','a1','b1']]
    numer=a1*(xs*h*h-nx)+b1*(ys*h**3+ny)
    denom=a0*(xs*h*h-nx)+b0*(ys*h**3+ny)
    assert denom!=0 and h!=0
    u=numer/denom;m=(a1-u*a0)/((u*b0-b1)*h);xp,yp=nx/h**2,ny/h**3
    assert ys==m*(xs-xp)-yp
    qc=[dec(c,u)for c in eq['quartic_coefficients']]
    t0,q0=[dec(c,u)for c in eq['quartic_zero']]
    sf=polynomial([dec(c,u)for c in eq['radical_square_factor']],t)
    ds=polynomial([dec(c,u)for c in eq['radical_denominator_sqrt']],t)
    v=(2*xs-(m*m-xp))*ds/sf;z=t-t0
    assert v*v==polynomial(qc,t) and q0*q0==polynomial(qc,t0)
    R=PolynomialRing(QQ,'z');shift=R(qc)(R.gen()+t0);e,d,c,b,a=[shift[i]for i in range(5)]
    assert z and q0
    X=(2*q0*(v+q0)+d*z)/z**2;Y=((X*X-4*e*a)*z-d*X-2*e*b)/(2*q0)
    EB=EllipticCurve(QQ,[dec(c,u)for c in eq['a_invariants']]);PB=EB(X,Y)
    zi=(2*q0*Y+d*X+2*e*b)/(X*X-4*e*a)
    vi=(X*zi*zi-d*zi-2*e)/(2*q0)
    assert zi+t0==t and vi==v
    assert (m*m-xp+vi*sf/ds)/2==xs
    un,ud=compact['u_of_s']['numerator'],compact['u_of_s']['denominator']
    assert len(un)==len(ud)==2
    s=(QQ(un[0])-u*QQ(ud[0]))/(u*QQ(ud[1])-QQ(un[1]))
    assert dec(compact['u_of_s'],s)==u
    gauge=dec(normalized['gauge_from_pointed_cubic'],u)
    scale=QQ(compact['scale']);base_den=polynomial(ud,s)
    mult=gauge*base_den**2/scale;assert mult
    XC=(X+EB.b2()/12)*mult**2;YC=Y*mult**3
    EC=EllipticCurve(QQ,[dec(c,s)for c in compact['a_invariants']])
    assert EC.a_invariants()[:3]==(0,0,0)
    PC=EC(XC,YC)
    assert EC.a4()==(-EB.c4()/48)*mult**4 and EC.a6()==(-EB.c6()/864)*mult**6
    genericB=[EC([dec(c,s)for c in p])for p in compact['basis_weierstrass_coordinates']]
    transported={'index':index,'id':row['id'],'public_word':row['public_word'],'source_parameter':str(t),
                 'source_point':jpoint(P),'source_short_point':list(map(str,[xs,ys])),
                 'u':str(u),'s':str(s),'quartic_point':list(map(str,[t,v])),
                 'pointed_cubic_point':jpoint(PB),'compact_scale':str(mult),
                 'compact_model':list(map(str,EC.a_invariants())),'compact_point':jpoint(PC),
                 'generic_points':[jpoint(P)for P in genericB],'maps_round_trip':True}
    save(case/'transport.json',transported)
    print('TRANSPORT',row['id'],'s bits',s.numerator().nbits()+s.denominator().nbits(),flush=True)
    # Cubic-descent characters at fixed ascending good primes. Discard a prime
    # if any displayed point hits a zero or pole in these simple coordinates.
    points=genericB+[PC];signatures=[[]for _ in points];blocks=[];torsion=None
    for prime in prime_range(3,2001):
        p=int(prime);Fp=GF(p)
        if any(c.denominator()%p==0 for c in EC.a_invariants()):continue
        ar,br=Fp(EC.a4()),Fp(EC.a6())
        if 4*ar**3+27*br**2==0:continue
        pol=R([EC.a6(),EC.a4(),0,1]).change_ring(Fp);roots=sorted(map(int,pol.roots(multiplicities=False)))
        if not roots:
            if torsion is None:torsion=p
            continue
        if any(P[0].denominator()%p==0 for P in points):continue
        values=[[Fp(P[0])-r for r in roots]for P in points]
        if any(not v for row0 in values for v in row0):continue
        bits=[[int(not v.is_square())for v in row0]for row0 in values]
        for target,part in zip(signatures,bits):target.extend(part)
        blocks.append({'prime':p,'roots':roots,'bits':bits})
        rank=int(matrix(GF(2),signatures).rank())
        if rank==18 and torsion is not None:break
    rank=int(matrix(GF(2),signatures).rank());grank=int(matrix(GF(2),signatures[:17]).rank())
    assert rank<=18 and grank<=17
    conclusion='OUTSIDE_GENERIC_RATIONAL_SPAN'if rank==18 and torsion else 'UNKNOWN_GENERIC_MEMBERSHIP'
    result={**transported,'finite_blocks':blocks,'no_rational_2_torsion_prime':torsion,
            'generic_mod2_rank':grank,'augmented_mod2_rank':rank,'conclusion':conclusion,
            'protocol_sha256':sha(PROTOCOL),'status':'PASS_EXACT_TRANSPORT_AND_FINITE_TEST'}
    save(case/'result.json',result)
    print('FINITE',row['id'],grank,rank,conclusion,flush=True)
def collect():
    policy=read(PROTOCOL);rows=[]
    for r in policy['rows']:
        p=WORK/r['id']/'result.json'
        rows.append(read(p)if p.exists()else {'id':r['id'],'status':'UNKNOWN_INCOMPLETE_CASE'})
    complete=all(r['status']=='PASS_EXACT_TRANSPORT_AND_FINITE_TEST'for r in rows)
    result={'schema':'curve302.class1-point-transport.v1','protocol_sha256':sha(PROTOCOL),
            'cases':rows,'complete':complete,'distinct_B_parameters':len({r['s']for r in rows if 's'in r}),
            'outside_generic_rational_span_count':sum(r.get('conclusion')=='OUTSIDE_GENERIC_RATIONAL_SPAN'for r in rows)}
    save(OUTPUT,result);print('COLLECT',complete,result['outside_generic_rational_span_count'],flush=True)
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['freeze','case','collect'])
    parser.add_argument('--index',type=int);args=parser.parse_args()
    if args.stage=='freeze':freeze()
    elif args.stage=='case':build_case(args.index)
    else:collect()
