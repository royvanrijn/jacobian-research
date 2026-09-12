#!/usr/bin/env sage-python
"""Compile the frozen 256-vector prefix using trace-section Riemann--Roch.

Adapts the existing H0(3O+kF) residual-chord method, adding the physical I2
component condition. Each candidate gets one child, 30 CPU seconds, no retry.
The geometry process reads no control classes or exceptional points.
"""
import hashlib
import json
import os
from pathlib import Path
import resource
import sys
import time
import traceback
from sage.all import QQ, ZZ, PolynomialRing, EllipticCurve, matrix, vector, lcm, gcd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/'artifacts/local/elliptic-curves/marked-two-class-v1'


def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def write(path, data):
    with path.with_suffix(path.suffix+'.tmp').open('w') as stream:
        json.dump(data, stream, indent=2, sort_keys=True)
        stream.write('\n')
    path.with_suffix(path.suffix+'.tmp').replace(path)


def ptext(p):
    return list(map(str, p.list())) or ['0']


def ratrec(p):
    return {'numerator': ptext(p.numerator()), 'denominator': ptext(p.denominator())}


def setup():
    g = json.loads((OUT/'geometry.json').read_text())
    R = PolynomialRing(QQ, 't'); K = R.fraction_field(); t = R.gen()
    f = g['family']
    A, B = [R(f[k]) for k in ['A_coefficients_low_to_high','B_coefficients_low_to_high']]
    E = EllipticCurve(K, [A, B])
    def dec(rec):
        return K(R(rec['numerator_coefficients_low_to_high']))/R(rec['denominator_coefficients_low_to_high'])
    points = [E(dec(s['X']), dec(s['Y'])) for s in f['sections']]
    Q = matrix(ZZ, g['source_ns_gram'])
    F, O, root = [vector(ZZ,g[k]) for k in ['fibre','zero','nonidentity_I2_component']]
    S = [vector(ZZ,s) for s in g['generic_divisors']]
    integral_basis = matrix(ZZ,[F,O,root,*S])
    assert abs(integral_basis.det()) == 1
    inverse = integral_basis.inverse()
    # The I2 fibre is lambda=infinity in the old pencil, with the exact
    # compact base matrix already verified by the inherited export.
    a,b,c,d = map(QQ,f['base_matrix_a_b_c_d'])
    if c:
        bad = -d/c
        aa,bb = A(bad),B(bad)
        location = {'chart':'finite','t':str(bad)}
    else:
        aa,bb = A[8],B[12]
        location = {'chart':'infinity'}
    assert aa != 0 and 4*aa**3+27*bb**2 == 0
    node_x = -3*bb/(2*aa)
    assert 3*node_x**2+aa == 0 and node_x**3+aa*node_x+bb == 0
    location['node_x'] = str(node_x); location['node_y'] = '0'
    return g,R,K,A,B,E,points,Q,F,O,root,inverse,location


def compile_one(index, row, ctx):
    start = time.process_time()
    g,R,K,A,B,E,points,Q,F,O,root,inverse,location = ctx
    D = vector(ZZ,row['native_ns_coordinates']); h=(row['norm']-10)//4
    assert D*Q*D == -2 and D*Q*F == 2 and D*Q*O == h
    m = ZZ(D*Q*root)
    assert 0 <= m <= 2 and D*Q*(F-root) == 2-m
    coefficients = D*inverse
    assert all(q in ZZ for q in coefficients)
    word = list(map(ZZ,coefficients[3:])); eps=m%2
    trace_height = QQ(row['norm'])-QQ(m*m)/2
    c = ZZ((trace_height-4+QQ(eps)/2)/2)
    k = h+c+6; vertical = (m+eps)//2
    T = sum((a*P for a,P in zip(word,points) if a),E(0))
    result = {'index':index,'divisor':row,'trace_word':list(map(int,word)),
        'I2_intersection':int(m),'trace_height':str(trace_height),'trace_zero_intersection':int(c),
        'RR_vertical_degree':int(k),'subtract_I2_multiplicity':int(vertical),
        'I2_location':location}
    if not T:
        return {**result,'status':'UNSUPPORTED_ZERO_TRACE','cpu_seconds':time.process_time()-start}
    Tminus = -T
    X,Y = Tminus[0],Tminus[1]
    den = R(X.denominator())
    assert den.is_square()
    hh=R(den.sqrt()).monic(); nx=R(X*hh**2); ny=R(Y*hh**3)
    assert Y.denominator() == hh**3
    c_actual=max(hh.degree(),(nx.degree()-4+1)//2,(ny.degree()-6+2)//3)
    assert c_actual == c
    # This is an identity of marked divisor classes, independently of RR.
    phi = sum((a*(S-O-(2+S*Q*O)*F+QQ(S*Q*root)/2*root)
               for a,S in zip(word,map(lambda s:vector(ZZ,s),g['generic_divisors']))),vector(QQ,[0]*19))
    minus_class = O+(2+c)*F-phi-QQ(eps)/2*root
    assert D+minus_class == 3*O+k*F-vertical*root
    bounds=[int(k),int(k-4),int(k-6)]
    polys=[]
    for bound,term in zip(bounds,[hh**3,nx*hh,ny]):
        polys += [term*R.gen()**j for j in range(bound+1)]
    degree=max(p.degree() for p in polys)
    equations=[[p[j] for p in polys] for j in range(degree+1)]
    if vertical:
        assert vertical == 1
        node_x=QQ(location['node_x'])
        node_row=[]
        for bound,channel in zip(bounds,[QQ(1),node_x,QQ(0)]):
            if location['chart']=='finite':
                node_row += [channel*QQ(location['t'])**j for j in range(bound+1)]
            else:
                node_row += [channel if j==bound else QQ(0) for j in range(bound+1)]
        equations.append(node_row)
    rr=matrix(QQ,equations)
    kernel=rr.right_kernel()
    result['RR_dimensions']={'rows':rr.nrows(),'columns':rr.ncols(),'kernel':kernel.dimension()}
    if kernel.dimension() != 1:
        return {**result,'status':'RR_KERNEL_NOT_ONE','cpu_seconds':time.process_time()-start}
    v=kernel.basis()[0];v*=lcm(q.denominator() for q in v);v=vector(ZZ,v);v/=gcd(list(v))
    if next(a for a in v if a)<0:v=-v
    relations=[];offset=0
    for bound in bounds:
        relations.append(R(list(v[offset:offset+bound+1])));offset+=bound+1
    f0,f1,f2=relations
    assert f0*hh**3+f1*nx*hh+f2*ny == 0
    result['RR_relation']=[ptext(f) for f in relations]
    if not f2:
        return {**result,'status':'NO_DEGREE_TWO_RESIDUAL','cpu_seconds':time.process_time()-start}
    XX=PolynomialRing(K,'x');x=XX.gen()
    elim=(K(f0)+K(f1)*x)**2-K(f2)**2*(x**3+K(A)*x+K(B))
    res,rem=elim.quo_rem(x-X)
    assert not rem and res.degree()==2
    res=res.monic(); discriminant=K(res[1]**2-4*res[0])
    if not discriminant:
        return {**result,'status':'NONREDUCED_RESIDUAL','cpu_seconds':time.process_time()-start}
    square=K(1);branch=K(1)
    for polynomial,sgn in [(discriminant.numerator(),1),(discriminant.denominator(),-1)]:
        fac=polynomial.factor();branch*=K(fac.unit())**sgn
        for factor,exponent in fac:
            exponent*=sgn
            square*=K(factor)**(exponent//2)
            if exponent%2:branch*=factor
    assert branch.denominator()==1
    branch=R(branch)
    assert discriminant == square**2*branch
    if branch.degree() not in (1,2):
        return {**result,'status':'RESIDUAL_NOT_GEOMETRIC_GENUS_ZERO','branch':ptext(branch),'cpu_seconds':time.process_time()-start}
    assert branch.gcd(branch.derivative()).degree()==0
    x0=-K(res[1])/2;x1=square/2
    y0=-(K(f0)+K(f1)*x0)/f2;y1=-K(f1)*x1/f2
    assert y0*y0+y1*y1*branch == x0**3+3*x0*x1*x1*branch+A*x0+B
    assert 2*y0*y1 == 3*x0*x0*x1+x1**3*branch+A*x1
    result.update({'status':'EXACT_GEOMETRIC_GENUS_ZERO_BISECTION',
        'branch':ptext(branch),'point_map':{key:ratrec(val) for key,val in [('x0',x0),('x1',x1),('y0',y0),('y1',y1)]},
        'trace_minus_map':{'x':ratrec(X),'y':ratrec(Y)},
        'residual_quadratic':[ratrec(K(res[i])) for i in range(3)],
        'discriminant_square_factor':ratrec(square),
        'trace_identity':'P_plus + P_minus = T over Q(t)',
        'irreducibility':'The nonconstant squarefree branch has degree1 or2, hence is nonsquare over Qbar(t). The RR divisor contains the certified trace and specified I2 component. Any residual vertical part V has D.V>=0 and C.V>=0, forcing V=0 from D^2=-2 and adjunction. Thus the actual bisection has class D, is smooth and geometrically rational.',
        'cpu_seconds':time.process_time()-start})
    return result


def main():
    directory=OUT/'compiled';directory.mkdir(exist_ok=False)
    start=time.process_time()
    protocol=json.loads((OUT/'protocol.json').read_text())
    assert digest(OUT/'geometry.json')==protocol['geometry_sha256']
    candidates=json.loads((OUT/'frozen-candidates.json').read_text())
    assert candidates['geometry_sha256']==protocol['geometry_sha256']
    write(OUT/'compilation-protocol.json',{'sources':{str(Path(__file__).relative_to(ROOT)):digest(Path(__file__))},
        'candidate_sha256':digest(OUT/'frozen-candidates.json'),'geometry_sha256':digest(OUT/'geometry.json'),
        'maximum_workers':1,'per_candidate_cpu_seconds':30,'retry':False,
        'method':'Marked trace section, polynomial H0(3O+kF) interpolation and resolved A1 first-vanishing constraint; exact residual quadratic and coefficient identities.'})
    snapshot=OUT/'source-snapshots'/Path(__file__).relative_to(ROOT)
    snapshot.write_bytes(Path(__file__).read_bytes())
    ctx=setup();receipts=[]
    for index,row in enumerate(candidates['rows']):
        path=directory/f'{index:03d}.json'
        wall=time.monotonic();pid=os.fork()
        if pid==0:
            resource.setrlimit(resource.RLIMIT_CPU,(30,31))
            try:
                result=compile_one(index,row,ctx)
                write(path,result)
                os._exit(0)
            except BaseException as exc:
                write(path,{'index':index,'divisor':row,'status':'COMPILATION_ERROR',
                    'exception':repr(exc),'traceback':traceback.format_exc(),'cpu_seconds':time.process_time()})
                os._exit(1)
        _,status,usage=os.wait4(pid,0)
        if not path.exists():
            write(path,{'index':index,'divisor':row,'status':'CPU_TIMEOUT_OR_WORKER_EXIT','wait_status':status})
        result=json.loads(path.read_text())
        receipt={'index':index,'status':result['status'],'wait_status':status,'cpu_seconds':usage.ru_utime+usage.ru_stime,
                 'wall_seconds':time.monotonic()-wall,'max_rss_kib':usage.ru_maxrss,'sha256':digest(path)}
        receipts.append(receipt)
        write(OUT/'compilation-progress.json',{'completed_candidates':len(receipts),'receipts':receipts})
        print(json.dumps(receipt),flush=True)
    from collections import Counter
    write(OUT/'compilation.json',{'status':'FROZEN_PREFIX_ATTEMPTED','candidates':len(receipts),
        'histogram':dict(Counter(r['status'] for r in receipts)), 'receipts':receipts,
        'parent_cpu_seconds':time.process_time()-start,'child_cpu_seconds':sum(r['cpu_seconds'] for r in receipts),
        'shells_complete':False,'boundary':'Only the frozen retained prefixes were compiled. No control splitting or labels influenced selection; no timeouts were retried.'})


if __name__=='__main__':
    main()
