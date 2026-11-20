#!/usr/bin/env python3
"""Exact finite controls for the global sphere purity recognition proof.

The universal algorithm is proved in the note. This small replay does not
implement a general multivariable creative-telescoping engine.
"""
import argparse
from collections import defaultdict
from fractions import Fraction as Q
from hashlib import sha256
from math import comb, factorial
from pathlib import Path
import json
import subprocess
import sys
import tempfile

if not __debug__:
    raise SystemExit('Assertions must be enabled; do not run with -O.')
HERE=Path(__file__).resolve().parent

def clean(p): return {e:Q(c) for e,c in p.items() if c}
def add(*ps):
    out=defaultdict(Q)
    for p in ps:
        for e,c in p.items(): out[e]+=c
    return clean(out)
def scale(p,a): return clean({e:c*a for e,c in p.items()})
def mul(p,q):
    out=defaultdict(Q)
    for a,c in p.items():
        for b,d in q.items(): out[tuple(x+y for x,y in zip(a,b))]+=c*d
    return clean(out)
def power(p,n):
    one={(0,)*len(next(iter(p))):Q(1)}
    while n:
        if n&1: one=mul(one,p)
        p=mul(p,p); n//=2
    return one
def deriv(p,i):
    out={}
    for e,c in p.items():
        if e[i]:
            ee=list(e); ee[i]-=1; out[tuple(ee)]=c*e[i]
    return out
def decode(rows):
    out={}
    for num,den,s,t in rows:
        assert den>0 and s>=0 and t>=0 and (s,t) not in out
        out[s,t]=Q(num,den)
    assert clean(out)==out
    return out
def rising_half(k):
    a=Q(1)
    for i in range(k): a*=Q(2*i+1,2)
    return a

def moment_sphere(p):
    out=Q(0)
    for (i,j,k),c in p.items():
        if i==j and k%2==0:
            r=k//2
            out+=c*Q(2**i*factorial(i)*factorial(2*r)*factorial(i+r),factorial(r)*factorial(2*(i+r)+1))
    return out

def contract(p,n):
    out=Q(0)
    for e,c in p.items():
        if e[:n]==e[n:]:
            for v in e[:n]: c*=factorial(v)
            out+=c
    return out

def operator_series(cs,seq,k):
    ans=Q(0)
    for j,p in enumerate(cs):
        for (v,b),c in p.items():
            assert b==0
            ell=k-v+j
            if ell>=j and ell<len(seq): ans+=c*seq[ell]*Q(factorial(ell),factorial(ell-j))
    return ans

def check(packet):
    A={(0,2):Q(1),(0,6):Q(-1)}
    Y={(0,0):Q(1),(1,2):Q(-1),(1,6):Q(1)}
    assert decode(packet['A'])==A and decode(packet['Y'])==Y
    r=packet['order']; assert r==4
    cs=[decode(v) for v in packet['operator']]
    R=decode(packet['certificate_R'])
    lhs=add(*[mul(cs[k],scale(mul(power(A,k),power(Y,r-k)),rising_half(k))) for k in range(r+1)])
    rhs=add(mul(deriv(R,1),Y),scale(mul(R,deriv(Y,1)),-Q(2*r-1,2)))
    assert lhs==rhs
    for endpoint in (0,1):
        boundary=defaultdict(Q)
        for (s,t),c in R.items(): boundary[s]+=c*endpoint**t
        assert not clean(boundary)
    corrupt=dict(R); corrupt[(0,0)]=corrupt.get((0,0),Q(0))+1
    badrhs=add(mul(deriv(corrupt,1),Y),scale(mul(corrupt,deriv(Y,1)),-Q(2*r-1,2)))
    assert badrhs!=lhs
    # Direct polynomial moments, independent of the algebraic integral derivation.
    p={(1,0,1):Q(1),(0,1,1):Q(2),(0,1,3):Q(2)}
    seq=[]
    for k in range(29):
        integral=sum((Q((-1)**j*comb(k,j),2*k+4*j+1) for j in range(k+1)),Q(0))
        seq.append(Q(comb(2*k,k),4**k)*integral)
    for m in range(17):
        raw=moment_sphere(power(p,m))
        expected=0 if m%2 else 4**(m//2)*seq[m//2]
        assert raw==expected
    assert seq[0]==1 and seq[1]==Q(2,21)
    assert moment_sphere(power(p,2))==Q(8,21)
    for k in range(23): assert operator_series(cs,seq,k)==0
    # Homogenize J-1: q=-L(1), M=(qD-q')L. Verify indicial polynomial literally.
    q=scale(cs[0],-1); qp=deriv(q,0)
    ms=[]
    for j in range(6):
        a=cs[j] if j<len(cs) else {}
        prev=cs[j-1] if j else {}
        ms.append(add(mul(q,prev),mul(q,deriv(a,0)),scale(mul(qp,a),-1)))
    valuations=[min(e[0] for e in a) if a else None for a in ms]
    r0=min(v-j for j,v in enumerate(valuations) if v is not None)
    assert r0==-2
    ind={}
    for j,v in enumerate(valuations):
        if v is None or v-j!=r0: continue
        fall={(0,0):Q(1)}
        for ell in range(j): fall=mul(fall,{(0,1):Q(1),(0,0):Q(-ell)})
        ind=add(ind,scale(fall,ms[j][(v,0)]))
    expected={(0,0):Q(-10584)}
    for a,b in [(1,0),(1,-1),(2,-1),(6,-7),(6,1)]:
        expected=mul(expected,{(0,1):Q(a),(0,0):Q(b)})
    assert ind==expected
    h=list(seq);h[0]-=1
    for k in range(22): assert operator_series(ms,h,k)==0
    # Literal bihomogeneous lift F=T Z R^2+2T W R^2+2T^3 W, degree(4,4).
    vv=[]
    for i in range(4):
        e=[0]*4;e[i]=1;vv.append({tuple(e):Q(1)})
    xi1,xi2,z1,z2=vv
    RR=add(mul(xi1,z1),mul(xi2,z2))
    TT=add(mul(xi1,z1),scale(mul(xi2,z2),-1))
    ZZ=mul(xi1,z2);WW=scale(mul(xi2,z1),2)
    assert add(power(TT,2),scale(mul(ZZ,WW),2))==power(RR,2)
    FF=add(mul(mul(TT,ZZ),power(RR,2)),scale(mul(mul(TT,WW),power(RR,2)),2),scale(mul(power(TT,3),WW),2))
    for e in FF: assert sum(e[:2])==4 and sum(e[2:])==4
    for m in range(7):
        actual=contract(power(FF,m),2)/factorial(4*m+1)
        assert actual==moment_sphere(power(p,m))
    # A direct regular rational certificate proving a genuinely three-height
    # positive-weight phase pure: z*d_z F=d_x(xF)-d_y(yF).
    q4={(0,0,0,2):Q(1),(0,1,1,0):Q(2),(0,0,0,0):Q(-1)}
    p4={(0,1,0,0):Q(1),(0,1,0,2):Q(1),(0,1,0,4):Q(1)}
    zz={(1,0,0,0):Q(1)};xx={(0,1,0,0):Q(1)};yy={(0,0,1,0):Q(1)}
    DD=mul(q4,add({(0,0,0,0):Q(1)},scale(mul(zz,p4),-1)))
    assert mul(zz,deriv(DD,0))==add(mul(xx,deriv(DD,1)),scale(mul(yy,deriv(DD,2)),-1))
    positive={(1,0,0):Q(1),(1,0,2):Q(1),(1,0,4):Q(1)}
    for m in range(1,11): assert moment_sphere(power(positive,m))==0
    # Universal centered-linear coefficient family: exact symbolic moments.
    linear={(1,0,0,1,0,0):Q(1),(0,1,0,0,1,0):Q(1),(0,0,1,0,0,1):Q(1)}
    quad={(1,1,0):Q(2),(0,0,2):Q(1)}
    for m in range(13):
        out=defaultdict(Q)
        for e,c in power(linear,m).items():
            out[e[:3]]+=c*moment_sphere({e[3:]:Q(1)})
        want={} if m%2 else scale(power(quad,m//2),Q(1,m+1))
        assert clean(out)==want
    # General markers: retain cancellations between bidegrees and output degrees.
    one4={(0,0,0,0):Q(1)}
    marker=add(one4,RR,z1,scale(mul(RR,z1),2),xi1,
               mul(xi2,power(z1,2)),power(z1,2),mul(RR,power(z2,2)))
    parts={}
    for ee,c in marker.items():
        rr,ss=sum(ee[:2]),sum(ee[2:])
        if ss>=rr: parts.setdefault((ss-rr,ss),{})[ee]=c
    def image(poly):
        out=defaultdict(Q)
        for ee,c in poly.items():
            aa,bb=ee[:2],ee[2:]
            if all(a<=b for a,b in zip(aa,bb)):
                out[tuple(b-a for a,b in zip(aa,bb))]+=c*Q(factorial(bb[0])*factorial(bb[1]),factorial(bb[0]-aa[0])*factorial(bb[1]-aa[1]))
        return clean(out)
    marker_tests=0
    for m in range(6):
        fm=power(TT,m);im=image(mul(marker,fm))
        for ell in sorted({a for a,b in parts}):
            emin=min(b for a,b in parts if a==ell)
            for g0 in range(ell+1):
                gamma=(g0,ell-g0)
                xig=mul(power(xi1,g0),power(xi2,ell-g0))
                agg=Q(0)
                for (ll,ss),part in parts.items():
                    if ll!=ell: continue
                    AA=contract(mul(mul(xig,part),fm),2)/factorial(m+ss+1)
                    ratio=Q(1)
                    for jj in range(emin+2,ss+2): ratio*=m+jj
                    agg+=ratio*AA
                expected=Q(factorial(g0)*factorial(ell-g0))*im.get(gamma,Q(0))/factorial(m+emin+1)
                assert agg==expected
                marker_tests+=1
    # General n-pair factorial normalization on the radial powers.
    radial_tests=0
    for n in range(1,5):
        rad={}
        for i in range(n):
            e=[0]*(2*n);e[i]=e[n+i]=1;rad[tuple(e)]=Q(1)
        for k in range(13):
            assert Q(factorial(n-1))*contract(power(rad,k),n)/factorial(k+n-1)==1
            radial_tests+=1
    return {'literal_certificate_terms':sum(map(len,cs))+len(R),'raw_sphere_moments':17,'ode_series_coefficients':23,'homogenized_series_coefficients':22,'original_balanced_lift_moments':7,'radial_n_pair_fixtures':radial_tests,'three_height_first_nonzero_order':2,'three_height_mu2':'8/21','homogenized_indicial':'-10584*N*(N-1)*(2*N-1)*(6*N-7)*(6*N+1)','J_minus_one_cutoff':1,'sphere_G_minus_one_cutoff':2,'corrupt_certificate_rejected':True,'pure_three_height_regular_certificate':True,'pure_three_height_raw_moments':10,'symbolic_linear_family_moments':13,'arbitrary_marker_euler_aggregation_fixtures':marker_tests}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--input',type=Path,default=HERE/'global_sphere_three_height_certificates.json');ap.add_argument('--output',type=Path,default=HERE/'global_sphere_three_height_verification.json');args=ap.parse_args()
    packet=json.loads(args.input.read_text(encoding='utf-8-sig'))
    counts=check(packet)
    ledgerpath=HERE/'global_sphere_purity_source_ledger.json'
    sourcechecks=[]
    ledger=json.loads(ledgerpath.read_text(encoding='utf-8-sig'))
    for record in ledger['files']:
        path=HERE/record['path']
        digest=sha256(path.read_bytes()).hexdigest()
        assert digest==record['sha256'], 'Source pin mismatch: '+record['path']
        sourcechecks.append({'path':record['path'],'sha256':digest})
    with tempfile.TemporaryDirectory(prefix='sphere-indicial-') as tmp:
        outpath=Path(tmp)/'receipt.json'
        run=subprocess.run([sys.executable,str(HERE/'verify_formal_indicial_purity_certificates.py'),'--output',str(outpath)],capture_output=True,text=True)
        assert run.returncode==0,run.stdout+run.stderr
        independent=json.loads(outpath.read_text())
        assert independent['passed'] is True
    receipt={'status':'PASS','scope':'Finite exact controls for the proof; not a general creative-telescoping implementation.','checker_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),'input_sha256':sha256(args.input.read_bytes()).hexdigest(),'checks':counts,'verified_source_pins':sourcechecks,'independent_indicial_controls':independent,'ledger_sha256':sha256(ledgerpath.read_bytes()).hexdigest()}
    args.output.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8');print(json.dumps(receipt,indent=2))
if __name__=='__main__':main()

