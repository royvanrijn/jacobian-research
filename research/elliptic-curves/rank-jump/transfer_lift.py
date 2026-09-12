"""Fresh-cover reduction, bounded blind recovery, and exact independent maps."""
import argparse
from fractions import Fraction as F
from pathlib import Path
import re
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

from transfer_common import *
from blind_constructed_cover_check import PAIRS, seq, det, mul, norm, evaluate, substitute, primitive
from blind_constructed_cover_quartic import add, scale, conv, ternary, binary, evalb, square_root


def quadrics(beta,f):
    qs=[[F(0)]*10 for _ in range(3)]
    for k,(i,j) in enumerate(PAIRS):
        if j==3:continue
        a=[F(i==n) for n in range(3)];b=[F(j==n) for n in range(3)]
        v=mul(beta,mul(a,b,f),f)
        for n in range(3):qs[n][k]=v[n]*(1 if i==j else 2)
    source=[qs[2].copy(),qs[1].copy()];source[1][-1]+=1
    return qs,source


def polynomial(q,names):
    pairs=[(i,j) for i in range(len(names)) for j in range(i,len(names))]
    return ' + '.join('(%s)*%s*%s'%(c,names[i],names[j]) for c,(i,j) in zip(q,pairs) if c) or '0'


def remote(folder,name,code):
    header='SetColumns(0); SetSeed(20260912);\nassert GetMemoryLimit() gt 0 and GetMemoryLimit() le 2147483648;\n'
    code=header+code
    assert len(code.encode())<=50000,'public service input limit; no split-job workaround'
    path=folder/(name+'.m');path.write_text(code)
    record=dict(request_sha256=sha(path),url='https://magma.maths.usyd.edu.au/xml/calculator.xml',service_cpu_limit_seconds=60,service_memory_limit_bytes=2147483648,remote_hardware='UNAVAILABLE',started_unix=time.time())
    req=urllib.request.Request(record['url'],data=urllib.parse.urlencode({'input':code}).encode(),headers={'Content-Type':'application/x-www-form-urlencoded','Referer':'https://magma.maths.usyd.edu.au/calc/'})
    tick=time.monotonic()
    try:
        with urllib.request.urlopen(req,timeout=65) as response:raw=response.read()
        tree=ET.fromstring(raw)
        result='\n'.join(''.join(line.itertext()) for line in tree.findall('.//results/line'))
        record['headers']={e.tag:e.text for e in tree.findall('.//headers/*')}
        record['reported_remote_cpu_seconds']=float(record['headers']['time']) if record['headers'].get('time') else None
    except Exception as exc:
        result=type(exc).__name__+': '+str(exc);raw=result.encode();record['error']=result
    record['request_wall_seconds']=time.monotonic()-tick
    path.with_suffix('.response.xml').write_bytes(raw)
    path.with_suffix('.result.txt').write_text(result)
    write(path.with_suffix('.meter.json'),record)
    if 'VERIFIED_TRANSFORMATION' not in result or 'DONE' not in result:
        raise ArithmeticError('bounded remote reduction unresolved; exact request and response retained')
    return result


def reduce_cover(folder,number):
    d=read(folder/'field.json');c=read(folder/f'cover-{number}.json')
    assert c['max_coefficient_bits']<=read(folder/'policy.json')['lifting']['maximum_cover_coefficient_bits']
    f=list(map(F,d['cubic_ascending']));beta=list(map(F,c['beta_ascending']))
    assert norm(beta,f)==F(c['norm'])==F(c['positive_norm_square_root'])**2
    qs,source=quadrics(beta,f)
    code='P<u,v,w,s>:=PolynomialRing(Rationals(),4);\nqs:=['+','.join(polynomial(q,'uvws') for q in source)+'];\n'
    code+='''M:=GenusOneModel(qs); m,g,levels:=Minimise(M); r,h:=Reduce(m);
print "REDUCED_MODEL"; print Eltseq(r);
print "COMPOSED_EQUATION_TRANS"; print Eltseq(Tuple(h*g)[1]);
print "COMPOSED_COORDINATE_TRANS"; print Eltseq(Tuple(h*g)[2]);
print "MINIMUM_LEVELS"; print levels;
assert (h*g)*M eq r; print "VERIFIED_TRANSFORMATION"; print "DONE";
'''
    text=remote(folder,f'minred-{number}',code)
    model=seq(text,'REDUCED_MODEL');A=seq(text,'COMPOSED_EQUATION_TRANS');S=seq(text,'COMPOSED_COORDINATE_TRANS')
    assert len(model)==20 and len(A)==4 and len(S)==16
    model=[model[:10],model[10:]];A=[A[:2],A[2:]];S=[S[4*i:4*i+4] for i in range(4)]
    assert det(A) and det(S)
    sub=[substitute(q,S) for q in source]
    assert model==[[sum(A[i][k]*sub[k][j] for k in range(2)) for j in range(10)] for i in range(2)]
    write(folder/f'reduction-{number}.json',dict(source=source,reduced=model,A=A,S=S))
    # Select any rational singular pencil member and move its vertex to the
    # first coordinate. This avoids assuming Magma chose the old conic chart.
    from sage.all import QQ, PolynomialRing, matrix, vector
    def symmetric(q):
        M=matrix(QQ,4)
        for c,(i,j) in zip(q,PAIRS):M[i,j]=M[j,i]=QQ(c)/(1 if i==j else 2)
        return M
    m0,m1=map(symmetric,model)
    T=PolynomialRing(QQ,'lambda');la=T.gen()
    candidates=[(F(1),F(0)),(F(0),F(1))]+[(F(str(r)),F(1)) for r in (la*m0.change_ring(T)+m1.change_ring(T)).det().roots(multiplicities=False)]
    transform=None
    for a,b in candidates:
        first=a*m0+b*m1
        if first.rank()!=3:continue
        vertex=first.right_kernel().basis()[0]
        second=m1 if a else m0
        if (vertex*second*vertex.column())[0]==0:continue
        rows=[vertex]
        for j in range(4):
            row=vector(QQ,[int(i==j) for i in range(4)])
            if matrix(QQ,rows+[row]).rank()>len(rows):rows.append(row)
            if len(rows)==4:break
        U=[list(map(lambda x:F(str(x)),row)) for row in rows]
        new0=substitute([a*x+b*y for x,y in zip(*model)],U)
        new1=substitute(model[1] if a else model[0],U)
        assert new0[:4]==[0]*4 and new1[0]
        transform=dict(U=U,pencil=[a,b],second=1 if a else 0,quadrics=[new0,new1])
        break
    if transform is None:raise ArithmeticError('no rational rank-three conic chart within fixed pencil extraction')
    write(folder/f'conic-chart-{number}.json',transform)
    q0,q1=transform['quadrics'];lead=q1[0]
    code='P<x,y,z>:=PolynomialRing(Rationals(),3);\nf:='+polynomial(q0[4:],'xyz')+';\nr:='+polynomial(q1[4:],'xyz')+';\nlin:='+' + '.join('(%s)*%s'%(v,n) for v,n in zip(q1[1:4],'xyz'))+';\n'
    code+='''C:=Conic(ProjectiveSpace(Rationals(),2),f); cm:=ParametrizationMatrix(C);
B<a,b>:=PolynomialRing(Rationals(),2); basis:=[a^2,a*b,b^2]; found:=false; pm:=cm; par:=basis;
for candidate in [cm,Transpose(cm),cm^(-1),Transpose(cm^(-1))] do
 par:=[ &+[basis[i]*candidate[i,j]:i in [1..3]] :j in [1..3]];
 if Evaluate(f,par) eq 0 then found:=true; pm:=candidate; break; end if;
end for;
assert found; print "CONIC_PARAMETRIZATION_MATRIX"; print Eltseq(pm);
'''
    code+='disc:=Evaluate(lin,par)^2 - (4*('+str(lead)+'))*Evaluate(r,par);\n'
    code+='''print "RAW_QUARTIC"; print [MonomialCoefficient(disc,a^(4-i)*b^i):i in [0..4]];
M:=GenusOneModel(disc); mm,g,levels:=Minimise(M); rr,h:=Reduce(mm);
print "REDUCED_QUARTIC"; print Eltseq(rr);
print "COMPOSED_SCALAR"; print Tuple(h*g)[1];
print "COMPOSED_BINARY_MATRIX"; print Eltseq(Tuple(h*g)[2]);
print "MINIMUM_LEVELS"; print levels;
assert (h*g)*M eq rr; print "VERIFIED_TRANSFORMATION"; print "DONE";
'''
    txt=remote(folder,f'quartic-minred-{number}',code)
    pm=seq(txt,'CONIC_PARAMETRIZATION_MATRIX');assert len(pm)==9
    pm=[pm[3*i:3*i+3] for i in range(3)];par=list(map(list,zip(*pm)))
    assert det(pm) and ternary(q0[4:],par)==[0]*5
    raw=seq(txt,'RAW_QUARTIC');reduced=seq(txt,'REDUCED_QUARTIC')
    scalar=F(re.search(r'COMPOSED_SCALAR\s*([^\s]+)',txt).group(1))
    T=seq(txt,'COMPOSED_BINARY_MATRIX');assert len(T)==4
    T=[T[:2],T[2:]]
    lin=[sum(q1[i+1]*par[i][j] for i in range(3)) for j in range(3)]
    assert len(raw)==len(reduced)==5 and det(T) and scalar
    assert raw==add(conv(lin,lin),scale(ternary(q1[4:],par),-4*lead))
    assert scale(binary(raw,T),scalar*scalar)==reduced
    write(folder/f'quartic-{number}.json',dict(pm=pm,par=par,raw=raw,reduced=reduced,scalar=scalar,T=T,lead=lead,linear=lin))


def search_cover(folder,number):
    from sage.all import pari
    q=read(folder/f'quartic-{number}.json');reduced=list(map(F,q['reduced']))
    y=square_root(reduced[0]);point=[F(1),y,F(0)] if y is not None else None
    attempts=[]
    for bound in read(folder/'policy.json')['lifting']['quartic_bounds']:
        if point is not None:break
        tick=time.monotonic();p=pari('Polrev(['+','.join(map(str,reversed(reduced)))+'])')
        pts=pari.hyperellratpoints(p,bound,1)
        attempts.append(dict(bound=bound,wall_seconds=time.monotonic()-tick,returned=len(pts)))
        if len(pts):point=[F(str(pts[0][0])),F(str(pts[0][1])),F(1)]
        write(folder/f'lifting-{number}-progress.json',dict(attempts=attempts,status='POINT_FOUND' if point else 'UNKNOWN_BOUNDED_MISS'))
    write(folder/f'lift-{number}.json',dict(status='CANDIDATE_PENDING_REPLAY' if point else 'UNKNOWN_BOUNDED_MISS',point=point,attempts=attempts,quartic_sha256=sha(folder/f'quartic-{number}.json')))


def verify_lift(folder,number):
    """Fraction map replay; no search or access to the baseline output."""
    lift=read(folder/f'lift-{number}.json')
    if lift['point'] is None:return
    assert lift['quartic_sha256']==sha(folder/f'quartic-{number}.json')
    d=read(folder/'field.json');c=read(folder/f'cover-{number}.json')
    q=read(folder/f'quartic-{number}.json');red=read(folder/f'reduction-{number}.json');conic=read(folder/f'conic-chart-{number}.json')
    scalar=F(q['scalar']);T=[list(map(F,row)) for row in q['T']];pm=[list(map(F,row)) for row in q['pm']]
    raw=list(map(F,q['raw']));reduced=list(map(F,q['reduced']));par=list(map(list,zip(*pm)))
    model=[list(map(F,row)) for row in red['reduced']];S=[list(map(F,row)) for row in red['S']];A=[list(map(F,row)) for row in red['A']]
    U=[list(map(F,row)) for row in conic['U']];a,b=map(F,conic['pencil']);qs=[list(map(F,row)) for row in conic['quadrics']]
    assert det(U) and qs[0]==substitute([a*x+b*y for x,y in zip(*model)],U) and qs[1]==substitute(model[conic['second']],U)
    assert ternary(qs[0][4:],par)==[0]*5
    lin=[sum(qs[1][i+1]*par[i][j] for i in range(3)) for j in range(3)]
    assert raw==add(conv(lin,lin),scale(ternary(qs[1][4:],par),-4*qs[1][0]))
    assert det(T) and scale(binary(raw,T),scalar*scalar)==reduced
    r,Y,t=map(F,lift['point']);assert Y*Y==evalb(reduced,r,t)
    aa=r*T[0][0]+t*T[1][0];bb=r*T[0][1]+t*T[1][1]
    oldY=Y/scalar;assert oldY*oldY==evalb(raw,aa,bb)
    tail=[evalb(p,aa,bb) for p in par]
    point=[(oldY-sum(qs[1][i+1]*tail[i] for i in range(3)))/(2*qs[1][0])]+tail
    assert all(evaluate(q,point)==0 for q in qs)
    point=[sum(point[i]*U[i][j] for i in range(4)) for j in range(4)]
    assert all(evaluate(q,point)==0 for q in model)
    point=[sum(point[i]*S[i][j] for i in range(4)) for j in range(4)]
    f=list(map(F,d['cubic_ascending']));beta=list(map(F,c['beta_ascending']));_,source=quadrics(beta,f)
    sub=[substitute(q,S) for q in source]
    assert det(A) and det(S) and model==[[sum(A[i][k]*sub[k][j] for k in range(2)) for j in range(10)] for i in range(2)]
    assert all(evaluate(q,point)==0 for q in source)
    u,v,w,s=point;assert s
    xi=[u/s,v/s,w/s];value=mul(beta,mul(xi,xi,f),f)
    assert value[2]==0 and value[1]==-1
    Z=value[0];W=F(c['positive_norm_square_root'])*norm(xi,f)
    assert norm(beta,f)==F(c['positive_norm_square_root'])**2
    assert W*W==sum(x*Z**i for i,x in enumerate(f))
    a1,a2,a3,a4,a6=map(F,d['minimal_model'])
    x=Z/4;y=(W-4*a1*x-4*a3)/8
    assert y*y+a1*x*y+a3*y==x**3+a2*x*x+a4*x+a6
    tu,tr,ts,tt=map(F,d['input_to_minimal'])
    X=tu*tu*x+tr;Y=tu**3*y+ts*tu*tu*x+tt
    model=list(map(F,read(folder/'seed.json')['curve']))
    assert Y*Y==X**3+model[3]*X+model[4]
    write(folder/f'point-{number}-verified.json',dict(status='PASS_EXACT_COVER_AND_MODEL_TRANSPORT',point=[X,Y],xi=xi,cover_sha256=sha(folder/f'cover-{number}.json'),point_on_compact_class=True))


def verify_compact(folder,number):
    from sage.all import pari, QQ
    from transfer_constructor import setup
    d=read(folder/'field.json');c=read(folder/f'cover-{number}.json');circuit=read(folder/f'compact-{number}-circuit.json');classes=read(folder/'classes.json')
    assert c['circuit_sha256']==sha(folder/f'compact-{number}-circuit.json')
    R,f,nf=setup(d)
    initial=[]
    for kind,i in classes['classes'][number]['factor_labels']:
        if kind=='generic':value=pari.Mod(pari(R(d['generic_classes'][i]['beta_ascending'])),pari(f))
        else:
            atom=classes['atoms'][i];value=int(atom['norm'])*pari.Mod(pari(R(atom['alpha_ascending'])),pari(f))
        initial.append(str(pari.nfalgtobasis(nf,value)))
    assert circuit['bases'][:len(initial)]==initial and circuit['exponents']==[1]*len(initial)+[-2]*len(circuit['steps'])
    assert circuit['bases'][len(initial):]==[step['multiplier'] for step in circuit['steps']]
    ideal=pari.idealhnf(nf,1)
    for step in circuit['steps']:
        assert pari(step['ideal_before'])==ideal
        assert pari.idealmul(nf,ideal,pari(step['factor_ideal']))==pari.idealmul(nf,pari(step['ideal_after']),pari(step['multiplier']))
        ideal=pari(step['ideal_after'])
    # Independently form the factored identity in the cubic algebra. A product
    # tree differs from producer nffactorback and retains exact denominators.
    terms=[pari.nfbasistoalg(nf,pari(b))**e for b,e in zip(circuit['bases'],circuit['exponents'])]
    while len(terms)>1:terms=[terms[i]*terms[i+1] if i+1<len(terms) else terms[i] for i in range(0,len(terms),2)]
    beta=pari.Mod(pari(R(c['beta_ascending'])),pari(f))
    assert terms[0]==beta and pari.idealhnf(nf,beta)==pari.idealpow(nf,ideal,2)
    assert QQ(pari.nfeltnorm(nf,beta))==QQ(c['norm'])==QQ(c['positive_norm_square_root'])**2
    write(folder/f'compact-{number}-verified.json',dict(status='PASS_EXACT_SQUARE_EQUIVALENCE_PRODUCT_TREE',cover_sha256=sha(folder/f'cover-{number}.json'),classes_sha256=sha(folder/'classes.json')))


def verify_rank(folder):
    initial=read(folder/'seed.json');model=list(map(F,initial['curve']));points=[list(map(F,p)) for p in initial['points']]
    count=read(folder/'classes-verified.json')['additional_classes'];recovered=[]
    for i in range(count):
        p=folder/f'point-{i}-verified.json'
        if not p.exists():continue
        comp=read(folder/f'compact-{i}-verified.json');rec=read(p)
        assert comp['cover_sha256']==rec['cover_sha256']==sha(folder/f'cover-{i}.json')
        assert comp['classes_sha256']==sha(folder/'classes.json')
        point=list(map(F,rec['point']));recovered.append(point)
    _,proof=admission(model,points+recovered)
    write(folder/'verified.json',dict(status='PASS_FRESH_CLASS_TO_POINT' if recovered else 'CLASSES_CERTIFIED_LIFTING_UNRESOLVED',transfer_success=bool(recovered),initial_rank=16,rank_lower_bound=len(points)+len(recovered),additional_classes=count,points=points+recovered,curve=model,proof=proof,ordinary_ideal_class_increment='UNKNOWN',solubility_boundary='A failed bounded lift is not a Sha certificate.'))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('mode',choices=['reduce_cover','search_cover','verify_lift','verify_compact','verify_rank']);ap.add_argument('--folder',type=Path,required=True);ap.add_argument('--number',type=int,default=0)
    args=ap.parse_args();folder=args.folder.resolve();guard(folder)
    from sage.all import pari
    pari.allocatemem(64000000,268435456,silent=True)
    if args.mode=='verify_rank':verify_rank(folder)
    else:globals()[args.mode](folder,args.number)
