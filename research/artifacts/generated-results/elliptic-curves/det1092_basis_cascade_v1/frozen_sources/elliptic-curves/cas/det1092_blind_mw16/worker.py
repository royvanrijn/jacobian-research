#!/usr/bin/env python3
"""Frozen blind RR worker. Reads exactly /input/fixture.json; no repo imports."""
import hashlib, json, os, resource, sys, time, traceback
from pathlib import Path
from sage.all import QQ, ZZ, PolynomialRing, EllipticCurve, matrix, vector
import sage.version

START = time.monotonic()
OUT = Path(os.environ['BLIND_OUTPUT'])
FIXTURE = Path(os.environ['BLIND_FIXTURE'])
def assert_redaction():
    try:
        open('/home/royvanrijn/src/jacobian-research/research/artifacts/generated-results/elliptic-curves/curve302_recovered_mw17_parent_v1.json').read(1)
    except PermissionError:
        return
    raise AssertionError('Full parent is accessible')

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def encode(v):
    return {'numerator':list(map(str,v.numerator().list())) or ['0'],
            'denominator':list(map(str,v.denominator().list()))}
def height(P):
    if P.is_zero(): return ZZ(0)
    x=P[0]
    return ZZ(max(x.numerator().degree(),x.denominator().degree()+4))
def gram(points):
    hs=list(map(height,points))
    return matrix(QQ,len(points),lambda i,j:(hs[i]+hs[j]-height(points[i]-points[j]))/2)
def usage():
    r=resource.getrusage(resource.RUSAGE_SELF)
    return {'wall_seconds':time.monotonic()-START,'user_seconds':r.ru_utime,
            'system_seconds':r.ru_stime,'max_rss_KiB':r.ru_maxrss}

def choose_centres(G):
    U=G.LLL_gram()
    assert abs(U.det())==1
    H=U.transpose()*G*U
    g=[[int(x) for x in row] for row in H.rows()]
    dirs=[(i,-1,0) for i in range(16)]+[(i,j,s) for i in range(16) for j in range(i+1,16) for s in [-1,1]]
    records=[]; pool={}
    for sample in range(512):
        mask=int.from_bytes(hashlib.sha256(('det1092-blind-mw16-v1:'+str(sample)).encode()).digest()[:2],'big')
        w=[(mask>>i)&1 for i in range(16)]
        z=[sum(g[i][j]*w[j] for j in range(16)) for i in range(16)]
        steps=0; capped=True
        for step in range(128):
            best=None
            for i,j,s in dirs:
                dot=z[i]+(s*z[j] if j>=0 else 0)
                dnorm=g[i][i]+(g[j][j]+2*s*g[i][j] if j>=0 else 0)
                delta=4*(dnorm-abs(dot))
                if delta<0 and (best is None or delta<best[0]): best=(delta,i,j,s,-1 if dot>0 else 1)
            if best is None: capped=False; break
            _,i,j,s,sign=best
            w[i]+=2*sign
            if j>=0: w[j]+=2*sign*s
            for k in range(16):z[k]+=2*sign*(g[k][i]+(s*g[k][j] if j>=0 else 0))
            steps+=1
        norm=sum(a*b for a,b in zip(w,z))
        original=tuple(map(int,U*vector(ZZ,w)))
        if next((a for a in original if a),1)<0:original=tuple(-a for a in original)
        records.append({'sample':sample,'mask_in_reduced_basis':mask,'norm':norm,'steps':steps,'step_cap_reached':capped})
        if 8<=norm<=14:
            key=tuple(a%2 for a in original)
            rankkey=(-norm,sum(map(abs,original)),original)
            if key not in pool or rankkey<pool[key][0]:pool[key]=(rankkey,original,norm)
    selected=sorted(pool.values())[:32]
    return {'LLL_columns':[list(map(int,row)) for row in U.rows()], 'samples':records,
            'eligible_distinct_parities':len(pool),
            'centres':[{'word':list(w),'height':n} for _,w,n in selected]}

def rr_lines(T,n,R,K):
    t=R.gen(); bounds=[n,n-4,n-6]
    funcs=[K(t**j)*v for v,b in zip([K(1),T[0],T[1]],bounds) for j in range(b+1)]
    den=R(1)
    for f in funcs:den=den.lcm(f.denominator())
    polys=[R(f*den) for f in funcs]
    M=matrix(QQ,max(p.degree() for p in polys)+1,len(polys),lambda i,j:polys[j][i])
    ker=M.right_kernel().basis_matrix();lines=[]
    for row in ker.rows():
        pos=0;fs=[]
        for b in bounds:fs.append(R(list(row[pos:pos+b+1])));pos+=b+1
        assert fs[0]+fs[1]*T[0]+fs[2]*T[1]==0
        lines.append(fs)
    return M,lines

def split_line(E,T,fs):
    K=E.base_ring(); X=PolynomialRing(K,'x');x=X.gen()
    f0,f1,f2=fs
    if not f2:return {'status':'VERTICAL_LINE_INHERITED','candidates':[]}
    a1,a2,a3,a4,a6=E.a_invariants();L=f0+f1*x
    residual,rem=(L**2-(a1*x+a3)*L*f2-f2**2*(x**3+a2*x*x+a4*x+a6)).quo_rem(x-T[0])
    assert not rem and residual.degree()==2
    c,b,a=residual.list();D=b*b-4*a*c
    row={'residual_coefficients':list(map(encode,[c,b,a])), 'discriminant':encode(D),'candidates':[]}
    if not D.is_square():row['status']='NONSPLIT_OVER_Q_T';return row
    s=D.sqrt();row['square_root']=encode(s);row['status']='SPLIT_OVER_Q_T'
    for sign in [1,-1]:
        xx=(-b+sign*s)/(2*a);P=E([xx,-(f0+f1*xx)/f2])
        if not any(P==Q for Q in row['candidates']):row['candidates'].append(P)
    assert len(row['candidates'])==1 or sum(row['candidates'],E(0))==-T
    return row

def selftest():
    R=PolynomialRing(QQ,'t');K=R.fraction_field()
    E=EllipticCurve(K,[0,0,0,0,1]);T=E([0,1])
    d=split_line(E,T,[R(-1),R(-1),R(1)])
    assert d['status']=='SPLIT_OVER_Q_T' and len(d['candidates'])==2
    assert sum(d['candidates'],E(0))==-T
    assert_redaction()
    assert not FIXTURE.exists()
    print('PASS synthetic residual and isolated runtime',sage.version.version,flush=True)

if '--selftest' in sys.argv:
    selftest();sys.exit(0)

resource.setrlimit(resource.RLIMIT_CPU,(180,185))
resource.setrlimit(resource.RLIMIT_AS,(8*1024**3,8*1024**3))
resource.setrlimit(resource.RLIMIT_FSIZE,(128*1024**2,128*1024**2))
D=json.loads(FIXTURE.read_text())
BASE={'starting_commit':D['starting_commit'],'arm_id':D['arm_id'],'fixture_sha256':sha(FIXTURE),'worker_sha256':sha(__file__)}
def event(kind,**kw):
    row=dict(BASE,event=kind,resource=usage(),**kw)
    with (OUT/'events.jsonl').open('a') as f:f.write(json.dumps(row,sort_keys=True)+'\n');f.flush();os.fsync(f.fileno())
def save(name,obj):
    with (OUT/name).open('x') as f:json.dump(dict(BASE,**obj),f,indent=2,sort_keys=True);f.write('\n')

result={'terminal_status':'EXPERIMENT_FAILURE','candidates':[],'stages':[],'software':sage.version.version}
try:
    assert set(D)=={'starting_commit','arm_id','a_invariants','sections','gram','policy'} and len(D['sections'])==16
    assert_redaction()
    event('START',visible_inputs=[str(FIXTURE),str(Path(__file__).resolve())],network_namespace='isolated',repository_visible=False)
    R=PolynomialRing(QQ,'t');K=R.fraction_field()
    def dec(v):return K(R(v['numerator']))/R(v['denominator'])
    E=EllipticCurve(K,list(map(dec,D['a_invariants'])))
    points=[E(list(map(dec,row))) for row in D['sections']]
    G=gram(points);assert G==matrix(QQ,D['gram']) and G.is_positive_definite() and G.rank()==16
    disc=R(E.discriminant());assert disc.degree()==24 and disc.gcd(disc.derivative()).degree()==0
    assert all(a.denominator().degree()==0 and a.numerator().degree()<=2*w for a,w in zip(E.a_invariants(),[1,2,3,4,6]))
    result['input_rank_certificate']={'gram':[list(map(str,row)) for row in G.rows()], 'determinant':str(G.det()),'rank':16,'positive_definite':True,'discriminant_degree':24,'discriminant_squarefree':True,'height_method':'24-I1 K3: h(P)=max(deg numerator x, deg denominator x+4)'}
    event('INPUT_RANK16_VERIFIED',determinant=str(G.det()))
    selection=choose_centres(G);save('selection.json',selection)
    event('CENTRES_FROZEN',count=len(selection['centres']))
    success=False
    for ci,centre in enumerate(selection['centres']):
        event('CENTRE_START',centre=ci,word=centre['word'],height=centre['height'])
        T=sum((c*P for c,P in zip(centre['word'],points)),E(0));assert height(T)==centre['height']
        h=centre['height'];first=(3*h+3)//4+1
        for n in [first,first+1]:
            stage={'centre':ci,'n':n,'status':'STARTED','members':[]};result['stages'].append(stage)
            event('RR_START',centre=ci,n=n)
            try:
                M,lines=rr_lines(T,n,R,K)
                stage.update(matrix_shape=list(M.dimensions()),matrix_rank=int(M.rank()),kernel_dimension=len(lines))
                for li,fs in enumerate(lines):
                    member={'index':li,'line_coefficients':[list(map(str,f.list())) for f in fs]}
                    split=split_line(E,T,fs);candidates=split.pop('candidates');member.update(split)
                    member['candidate_ids']=[];stage['members'].append(member)
                    for Q in candidates:
                        xy=list(map(encode,Q.xy()));old=next((j for j,c in enumerate(result['candidates']) if c['coordinates']==xy),None)
                        if old is not None:member['candidate_ids'].append(old);continue
                        cross=vector(QQ,[(height(Q)+height(P)-height(Q-P))/2 for P in points])
                        schur=height(Q)-cross*G.inverse()*cross
                        assert schur>=0
                        H=G.augment(matrix(QQ,16,1,list(cross))).stack(matrix(QQ,1,17,list(cross)+[height(Q)]))
                        rank=int(H.rank());assert rank in [16,17]
                        candidate={'coordinates':xy,'height':int(height(Q)), 'pairings':list(map(str,cross)), 'schur_complement':str(schur),'rank':rank,'gram_determinant':str(H.det()),'producer':[ci,n,li]}
                        idx=len(result['candidates']);result['candidates'].append(candidate);member['candidate_ids'].append(idx)
                        event('CANDIDATE_CERTIFIED',candidate_id=idx,candidate=candidate)
                        if rank==17:success=True
                    if success:break
                stage['status']='COMPLETE_THROUGH_FIRST_SUCCESS' if success else 'COMPLETE'
                event('RR_COMPLETE',stage=stage)
            except Exception as exc:
                stage['status']='FAILED';stage['error']=repr(exc);stage['traceback']=traceback.format_exc();event('RR_FAILED',stage=stage)
            if success:break
        if success:break
    result['terminal_status']='GENERIC_QUOTIENT_RECOVERY' if success else ('EXPERIMENT_FAILURE' if any(s['status']=='FAILED' for s in result['stages']) else 'BOUNDED_NO_RECOVERY')
    result['stopping_policy']='Stop at the first member with exact rank17, retaining all its rational roots; unexecuted selected centres, bounds and kernel rows are censored after success.'
    result['selected_centres']=len(selection['centres'])
    result['censored_after_success']=success
except BaseException as exc:
    result['error']=repr(exc);result['traceback']=traceback.format_exc();event('WORKER_FAILED',error=repr(exc))
finally:
    result['resource']=usage();save('result.json',result);event('TERMINAL',status=result['terminal_status']);print(D['starting_commit'],D['arm_id'],result['terminal_status'],usage(),flush=True)
