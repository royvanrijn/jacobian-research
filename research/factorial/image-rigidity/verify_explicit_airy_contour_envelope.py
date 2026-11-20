#!/usr/bin/env python3
"""Exact support for the written global joint-Airy contour bound, not an SFC cutoff."""
if not __debug__:raise SystemExit('Refusing -O/-OO: exact assertions are required.')
from fractions import Fraction as F
from pathlib import Path
import argparse,hashlib,json
HERE=Path(__file__).resolve().parent
ROWS=[]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def ck(name,left,right,relation='<'):
    left,right=F(left),F(right)
    assert (left<right if relation=='<'else left==right),(name,left,relation,right)
    ROWS.append({'name':name,'left':str(left),'relation':relation,'right':str(right)})
def exp_upper(x):
    x=F(x);assert 0<=x<66
    term=total=F(1)
    for k in range(1,65):term*=x/k;total+=term
    return total+term*x/65/(1-x/66)
def add(*polys):
    out={}
    for p in polys:
        for m,c in p.items():out[m]=out.get(m,F(0))+c
    return {m:c for m,c in out.items()if c}
def sc(p,c):return {m:F(c)*v for m,v in p.items()if c*v}
def mul(p,q):
    out={}
    for u,c in p.items():
        for v,d in q.items():
            m=tuple(i+j for i,j in zip(u,v));out[m]=out.get(m,F(0))+c*d
    return {m:c for m,c in out.items()if c}
def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--repo-root',type=Path,required=True)
    ap.add_argument('--output',type=Path,default=HERE/'explicit_airy_contour_envelope_verification.json');args=ap.parse_args()
    repo=args.repo_root.resolve();lp=HERE/'explicit_airy_contour_envelope_source_ledger.json';ledger=json.loads(lp.read_text())
    expected={('local','EXPLICIT_AIRY_CONTOUR_ENVELOPE.md'),('local',Path(__file__).name),
       ('local','sources/dlmf-9-2-endpoint.html'),('local','sources/dlmf-9-7-endpoint.html'),
       ('repo','research/factorial/image-rigidity/GAP_TWO_LARGE_LOWER_EXPONENT_ENDPOINT_THEOREM.md')}
    assert ledger['schema_version']==1 and len(ledger['files'])==len(expected)
    assert {(e['root'],e['path'])for e in ledger['files']}==expected
    for e in ledger['files']:
        p=(HERE if e['root']=='local'else repo)/e['path'];assert p.is_file()and sha(p)==e['sha256'],e['path']
    ck('Gamma2/3 integral lower',1/(F(2,3)*F(5,3)),F(9,10),'=')
    ck('Gamma1/3 integral lower',1/(F(1,3)*F(4,3)),F(9,4),'=')
    ck('Bi initial root upper',3,F(4,3)**6)
    ck('Bi initial majorant',F(10,9)+F(8,27),F(3,2))
    ck('exp8 exact series',exp_upper(8),3000)
    ck('Bi disk majorant',F(3,2)*3000,4500,'=')
    ck('Bi derivative disk majorant',2*4500,9000,'=')
    ck('compact exponential',exp_upper(F(16,3)),225)
    ck('compact bend bracket at6',5-F(8,9)*36+6,0)
    ck('compact bend bracket derivative',1-F(16,9)*6,0)
    ck('compact arc factor squared',F(4,3),F(49,36))
    ck('sqrtpi compact bound squared',F(22,7),F(81,25))
    ck('compatible compact integral',F(7,6)*F(9,5)*225,500)
    ck('compact joint ratio',500*36000,2*10**7)
    ck('chi Ai parameter',F(7,6),2)
    ck('chi derivative parameter',1,2)
    ck('last Ai multiplier square-root bound',F(7,3)*F(22,7),F(49,4))
    ck('last Ai multiplier upper',F(7,2)*2+3,10,'=')
    ck('last derivative multiplier upper',2*2+3,7,'=')
    ck('primary action modulus at r4',F(2,3)*8,F(16,3),'=')
    errA=10*F(5,72)*F(3,16);errD=7*F(7,72)*F(3,16)
    ck('Ai leading relative error',errA,F(1,4));ck('derivative leading relative error',errD,F(1,4))
    ck('positive contour ratio',6*4/F(3,4),32,'=')
    rotA=3*F(5,72)*F(3,16);rotD=3*F(7,72)*F(3,16)
    ck('rotated Ai relative error',rotA,F(5,128),'=')
    ck('rotated derivative relative error',rotD,F(7,128),'=')
    ck('rotated joint lower',1-rotA-rotD,F(29,32),'=')
    ck('rotated lower margin',F(3,4),F(29,32))
    ck('negative contour ratio',50*4/F(3,4),300)
    ck('conditional normalized local CR',2300*F(2*10**7,6),8*10**9)
    # Universal compact phase and ODE majorant identities.
    def var(j):return {tuple(int(i==j)for i in range(4)):F(1)}
    X,Y,U,V=[var(j)for j in range(4)];ONE={(0,0,0,0):F(1)}
    def cm(p,q):return(add(mul(p[0],q[0]),sc(mul(p[1],q[1]),-1)),add(mul(p[0],q[1]),mul(p[1],q[0])))
    def phase_re(real,imag):
        pair=(real,imag);cube=cm(cm(pair,pair),pair);linear=cm((U,V),pair)
        return add(sc(cube[0],F(1,3)),sc(linear[0],-1))
    difference=add(phase_re(add(ONE,X),Y),sc(phase_re(ONE,Y),-1))
    target=add(mul(X,add(ONE,sc(U,-1),sc(mul(Y,Y),-1))),mul(X,X),sc(mul(mul(X,X),X),F(1,3)))
    assert difference==target
    exponent=add(sc(ONE,F(1,3)),sc(U,-1),sc(add(sc(ONE,16),sc(mul(U,U),-1)),F(1,4)))
    complete_square=add(sc(ONE,F(16,3)),sc(mul(add(U,sc(ONE,2)),add(U,sc(ONE,2))),-F(1,4)))
    assert exponent==complete_square
    # For B=Bi and C=Bi', 2N-N'=(2-x/2)B.
    ode_difference=add(sc(U,2),V,sc(V,-1),sc(mul(X,U),-F(1,2)))
    assert ode_difference==mul(add(sc(ONE,2),sc(X,-F(1,2))),U)
    # Connection coefficients as powers of exp(i*pi/12), sign encoded by adding12.
    def phase(k,sign=1):return(k+(12 if sign==-1 else 0))%24
    assert phase(16-1,-1)==phase(3)
    assert phase(8+1,-1)==phase(-3)
    assert phase(8+1)==phase(6+3)
    assert phase(16-1)==phase(6-3,-1)
    # The coefficient matrix [[1,1],[i,-i]] has determinant-2i.
    def ga(p,q):return(p[0]+q[0],p[1]+q[1])
    def gm(p,q):return(p[0]*q[0]-p[1]*q[1],p[0]*q[1]+p[1]*q[0])
    ii=(F(0),F(1));mi=(F(0),F(-1));oo=(F(1),F(0))
    assert ga(gm(oo,mi),gm((-oo[0],-oo[1]),ii))==(0,-2)
    # Actual sign controls: a left bend grows cubically; Ai alone can cancel.
    assert F(-1,24)+F(1,2)==F(11,24)>0
    u=oo;v=(-F(1),F(0))
    assert ga(u,v)==(0,0)and ga(gm(ii,u),gm(mi,v))==(0,2)
    result={'schema_version':1,'status':'PASS','scope':'Global absolute Ai-contour bound C_A=2e7, with compatible bent/rotated cycles. No actual CFU remainder constant or numerical SFC cutoff supplied.',
       'rational_guards':ROWS,'universal_identity_count':8,'negative_controls':2,
       'primary_error_scope':'Actual DLMF9.7.17 n1, sigma1/6Ai and0Aiprime, exact error sectors; no asymptotic inference or zero sampling.',
       'source_ledger_sha256':sha(lp),'proof_sha256':sha(HERE/'EXPLICIT_AIRY_CONTOUR_ENVELOPE.md')}
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':'PASS','rational_guards':len(ROWS),'universal_identities':8,'negative_controls':2}))
if __name__=='__main__':main()
