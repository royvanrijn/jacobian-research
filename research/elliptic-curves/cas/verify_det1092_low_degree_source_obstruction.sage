#!/usr/bin/env sage-python
"""Independent affine-ellipsoid, modular-map and exact survivor replay.

The producer enumerates a parity coset of radius16; this verifier enumerates
an affine integer lattice of radius4. Finite point reduction uses primitive
projective polynomials, and words use manual reverse group law, not the
producer's affine Gauss valuation routine or Sage point additions. <=25s.
"""
import argparse,hashlib,json,signal
from fractions import Fraction as Rat
from math import isqrt
from pathlib import Path
from sage.all import QQ,ZZ,GF,PolynomialRing,EllipticCurve,matrix,vector,gcd,lcm
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'det1092_low_degree_source_obstruction_v1'
PARENT=ART/'curve302_recovered_mw17_parent_v1.json';PENCIL=ART/'det1092_norm8_seed_cover_v2/generic.json'
PICARD=ART/'det1092_genus1_picard_image_v1/replay.json'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def rat(x):return Rat(int(QQ(x).numerator()),int(QQ(x).denominator()))
def save(name,data):
    path=OUT/name
    if path.exists():assert read(path)==data,'immutable replay changed'
    else:
        with path.open('x') as f:json.dump(data,f,indent=2,sort_keys=True);f.write('\n')
def check_bound_inputs(row):
    for name,h in row['inputs'].items():assert sha(ROOT/name)==h,'input changed: '+name
def add(P,Q,A):
    if P is None:return Q
    if Q is None:return P
    x,y=P;u,v=Q
    if x==u and y==-v:return None
    slope=(v-y)/(u-x) if x!=u else (3*x*x+A)/(2*y)
    z=slope*slope-x-u
    return z,slope*(x-z)-y
def main():
    geometry=read(OUT/'geometry-protocol.json');protocol=read(OUT/'incidence-protocol-v2.json')
    lifts=read(OUT/'lift-protocol.json');lift_summary=read(OUT/'lift-summary.json')
    for record in [geometry,protocol,lifts,lift_summary]:check_bound_inputs(record)
    source=read(OUT/'sources.json');parent=read(PARENT);pencil=read(PENCIL)
    assert source['geometry_protocol_sha256']==sha(OUT/'geometry-protocol.json')
    G=matrix(ZZ,parent['generic_height_gram']);w=vector(ZZ,source['word']);U=matrix(ZZ,source['LLL_columns'])
    assert w*G*w==8 and abs(U.det())==1
    frame=read(ART/'det1092_genus1_picard_image_v1/frame.json');assert vector(QQ,frame['D'])==vector(QQ,[2,4,*list(w)])
    n=17;H=U.transpose()*G*U;center=[rat(a/2) for a in U.inverse()*w]
    L=[[Rat(i==j) for j in range(n)] for i in range(n)];D=[]
    for j in range(n):
        d=rat(H[j,j])-sum(L[j][k]**2*D[k] for k in range(j));assert d>0;D.append(d)
        for i in range(j+1,n):L[i][j]=(rat(H[i,j])-sum(L[i][k]*L[j][k]*D[k] for k in range(j)))/d
    assert all(rat(H[i,j])==sum(L[i][k]*D[k]*L[j][k] for k in range(n)) for i in range(n) for j in range(n))
    z=[0]*n;found=[];nodes=0
    def visit(j,used):
        nonlocal nodes
        nodes+=1;assert nodes<=100000,'UNKNOWN_INDEPENDENT_NODE_CAP'
        if j<0:
            v=U*vector(ZZ,z);degree=v*G*v-w*G*v
            assert rat(degree+2)==used and degree in [0,1,2]
            found.append((int(degree),tuple(map(int,v))));return
        shift=-center[j]+sum((L[k][j]*(Rat(z[k])-center[k]) for k in range(j+1,n)),Rat(0))
        budget=(Rat(4)-used)/D[j]
        if budget<0:return
        b,a=shift.denominator,shift.numerator
        radius=isqrt((budget*b*b).numerator//(budget*b*b).denominator)
        for value in range(-((radius+a)//b),(radius-a)//b+1):
            z[j]=value;visit(j-1,used+D[j]*(Rat(value)+shift)**2)
    visit(16,Rat(0))
    expected=[(r['degree'],tuple(r['word'])) for r in source['sources']]
    assert sorted(found)==expected and len(set(found))==len(found)==1308
    counts={str(d):sum(v[0]==d for v in found) for d in range(3)}
    assert counts==source['counts']=={'0':10,'1':130,'2':1168}
    assert all(r['centered_norm']==4*r['degree']+8 for r in source['sources'])
    R=PolynomialRing(QQ,'t');K=R.fraction_field()
    def dec(v):return K(R(v['numerator']))/R(v['denominator'])
    original=EllipticCurve(K,[dec(v) for v in parent['a_invariants']]);Aq,Bq=-original.c4()/48,-original.c6()/864
    basis=[]
    for row in parent['basis_weierstrass_coordinates']:
        x,y=map(dec,row);X,Y=x+original.b2()/12,y+(original.a1()*x+original.a3())/2
        assert Y*Y==X**3+Aq*X+Bq;basis.append((X,Y))
    hq=R(pencil['pole_h']);shiftq=R(pencil['shift'])
    cxq,cyq=K(R(pencil['nx']))/hq**2,K(R(pencil['ny']))/hq**3
    zstarq=QQ(protocol['carrier_label'])
    # Clear polynomial coefficient denominators, then integer content.
    # This is projective reduction over the Gauss DVR, with no valuation code.
    def primitive(polys):
        den=lcm([v.denominator() for p in polys for v in p.list()])
        integers=[R(p*den) for p in polys];content=gcd([ZZ(v) for p in integers for v in p.list()])
        assert content
        return [R(p/content) for p in integers]
    def pair(f):
        f=K(f);return primitive([f.numerator(),f.denominator()])
    triples=[]
    for X,Y in basis:
        den=lcm(X.denominator(),Y.denominator());triples.append(primitive([R(X*den),R(Y*den),R(den)]))
    pairs=[pair(v) for v in [Aq,Bq,hq,shiftq,cxq,cyq,zstarq]]
    active=[i for i,r in enumerate(source['sources']) if r['degree']==2];closed=[];prime_records=[]
    for p in protocol['primes']:
        saved=read(OUT/('prime-%d.json'%p));assert saved['incoming']==active
        assert saved['protocol_sha256']==sha(OUT/'incidence-protocol-v2.json')
        assert ZZ(p).is_prime(proof=True) and p>2
        fp=GF(p);Rp=PolynomialRing(fp,'t');Fp=Rp.fraction_field()
        reductions=[]
        for N,D0 in pairs:
            nn,dd=Rp(N),Rp(D0);assert nn or dd
            reductions.append(None if not dd else Fp(nn)/dd)
        A,B,h,shift,cx,cy,zstar=reductions
        if A is None or B is None or not (4*A**3+27*B**2):
            assert saved['status']=='UNKNOWN_BAD_GENERIC_REDUCTION' and saved['survivors']==active and not saved['trials']
            prime_records.append(dict(prime=p,status=saved['status'],remaining=len(active)));continue
        if any(v is None for v in [h,shift,cx,cy,zstar]) or not h:
            assert saved['status']=='UNKNOWN_PENCIL_GAUSS_CHART' and saved['survivors']==active and not saved['trials']
            prime_records.append(dict(prime=p,status=saved['status'],remaining=len(active)));continue
        reduced_basis=[]
        for triple in triples:
            X,Y,Z=map(Rp,triple);assert X or Y or Z
            if not Z:assert not X and Y;reduced_basis.append(None)
            else:
                xx,yy=Fp(X)/Z,Fp(Y)/Z;assert yy*yy==xx**3+A*xx+B;reduced_basis.append((xx,yy))
        cache={(0,)*17:None}
        def word_mod(word):
            key=tuple(word)
            if key not in cache:
                i=max(j for j,a in enumerate(key) if a);sgn=1 if key[i]>0 else -1
                smaller=list(key);smaller[i]-=sgn;P=reduced_basis[i]
                signed=None if P is None else (P[0],sgn*P[1])
                cache[key]=add(word_mod(smaller),signed,A)
            return cache[key]
        assert [r['source_index'] for r in saved['trials']]==active
        remaining=[];nexcluded=0
        for row in saved['trials']:
            i=row['source_index'];P=word_mod(source['sources'][i]['word'])
            if P is None:assert row['status']=='UNKNOWN_REDUCED_SECTION_IS_O'
            elif P[0]==cx:assert row['status']=='UNKNOWN_PENCIL_DENOMINATOR'
            else:
                zmap=(h*(P[1]+cy)/(P[0]-cx)+shift)/(h*h)
                num,den=zmap.numerator(),zmap.denominator();degree=max(num.degree(),den.degree())
                assert num==Rp(row['map_numerator']) and den==Rp(row['map_denominator']) and degree==row['degree']
                if degree!=2:assert row['status']=='UNKNOWN_MAP_DEGREE_DROP'
                else:
                    f=Rp(num-zstar*den);delta=f[1]**2-4*f[0]*f[2]
                    assert [int(f[j]) for j in range(3)]==row['incidence_coefficients'] and int(delta)==row['discriminant']
                    excluded=bool(f[2] and delta and pow(int(delta),(p-1)//2,p)==p-1)
                    assert (row['status']=='EXCLUDED_NO_PROJECTIVE_ROOT')==excluded
                    if not excluded:assert row['status']=='UNKNOWN_LOCAL_ROOT'
            if row['status']=='EXCLUDED_NO_PROJECTIVE_ROOT':closed.append(i);nexcluded+=1
            else:remaining.append(i)
        assert remaining==saved['survivors'];active=remaining
        prime_records.append(dict(prime=p,status=saved['status'],new_exclusions=nexcluded,remaining=len(active)))
    assert active==lifts['sources']==[179,1268] and len(closed)==len(set(closed))==1166
    def word_q(word):
        total=None
        for n,P in reversed(list(zip(word,basis))):
            for unused in range(abs(n)):total=add(total,(P[0],P[1] if n>0 else -P[1]),Aq)
        return total
    maps=[]
    for i in active:
        data=read(OUT/('lift-%d.json'%i));word=source['sources'][i]['word'];P=word_q(word)
        assert data['word']==word and tuple(dec(v) for v in data['section'])==P
        z=(hq*(P[1]+cyq)/(P[0]-cxq)+shiftq)/(hq*hq);assert z==dec(data['z_map'])
        assert max(z.numerator().degree(),z.denominator().degree())==2;maps.append(z)
        f=R(z.numerator()-zstarq*z.denominator());den=lcm([a.denominator() for a in f.list()])
        f=R(f*den);f=R(f/gcd([ZZ(a) for a in f.list()]))
        if f.leading_coefficient()<0:f=-f
        assert f==R(data['primitive_incidence_polynomial']) and f.degree()==2
        assert not data['infinity_preimage'] and not data['rational_preimages']
        cert=data['square_certificate'];delta=f[1]**2-4*f[0]*f[2];assert delta==QQ(cert['discriminant'])
        if delta<0:assert cert['negative']
        else:
            n,d=delta.numerator(),delta.denominator();a,b=ZZ(cert['numerator_floor_sqrt']),ZZ(cert['denominator_floor_sqrt'])
            assert n==ZZ(cert['numerator']) and d==ZZ(cert['denominator'])
            assert a*a<=n<(a+1)**2 and b*b<=d<(b+1)**2 and (a*a!=n or b*b!=d)
        assert data['status']=='EXACT_NO_RATIONAL_SOURCE_POINT'
    assert maps[0]==maps[1]
    assert vector(ZZ,source['sources'][active[0]]['word'])+vector(ZZ,source['sources'][active[1]]['word'])==w
    picard=read(PICARD);check_bound_inputs(picard)
    assert picard['status']=='PASS_INDEPENDENT_FULL_PICARD_RANK12_AND_FIRST_CLASS_RANK13'
    assert picard['first_augmented_rank']==13
    return dict(status='PASS_INDEPENDENT_INTRINSIC_LOW_DEGREE_SOURCE_ORBIT_OBSTRUCTION',
        classification='new intrinsic-source obstruction, retrospectively evaluated on first302 carrier',
        degree_counts=counts,independent_affine_nodes=nodes,modular_exclusions=1166,exact_quadratic_exclusions=2,
        prime_records=prime_records,conjugate_survivors=active,
        deduction='Neither first seed nor carrier conjugate lies on the orbit, under arbitrary generic alternate translations and inversion, of any original section of alternate degree<=2.',
        boundary='Not all original sections, alternating-fibration automorphisms, other302 seeds, or all representatives of this original quotient direction. No prospective positive302 seed.',
        checker_sha256=sha(Path(__file__)),
        inputs={str(p.relative_to(ROOT)):sha(p) for p in [OUT/'geometry-protocol.json',OUT/'sources.json',OUT/'incidence-protocol-v2.json',
            OUT/'lift-protocol.json',OUT/'lift-summary.json',PICARD,*[OUT/('prime-%d.json'%p) for p in protocol['primes']],
            *[OUT/('lift-%d.json'%i) for i in active]]})
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args()
    signal.alarm(25);result=main()
    if a.write:save('independent-replay.json',result)
    else:assert read(OUT/'independent-replay.json')==result
    print(result['status'],result['degree_counts'],flush=True)
