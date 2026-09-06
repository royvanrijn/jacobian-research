#!/usr/bin/env python3
"""One frozen native triple intersection; no parameter or point search."""
import argparse
from fractions import Fraction as F
from pathlib import Path
import subprocess
import retrospective as r

HERE=Path(__file__).resolve().parent
PROTOCOL=HERE/'NATIVE_TRIPLE_INTERSECTION_PROTOCOL.json'
BASE=r.OUT/'rank_jump_native_pair_collapse_locus_inputs_v1.json'
ATLAS=r.ROOT/'artifacts/generated-results/elkies-2026-equation-bisections-full.json'
RELATIONS=r.OUT/'rank_jump_paired_quartet_relations_verification_v1.json'
TRANSPORT=r.OUT/'compact_published_r17_generic_transport_v1.json'
LATTICE=r.OUT/'rank_jump_norm_six_carrier_solubility_inputs_v1.json'
INPUT=r.OUT/'rank_jump_native_triple_intersection_inputs_v1.json'
OUTPUT=r.OUT/'rank_jump_native_triple_intersection_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-native-triple-intersection-v1'


def capture():
    base=r.read(BASE);old=next(x for x in r.read(RELATIONS)['rows'] if x['id']=='08234-003')
    assert old['kernel_integer_vectors']==[[1,-1,0,1]]
    mat=r.read(TRANSPORT)['compact_sections_in_published_basis']
    s=[sum(F(v)*row[i] for v,row in zip(old['kernel_generic_coordinates'][0],mat,strict=True)) for i in range(17)]
    atlas=r.read(ATLAS)
    covers=[next(c for c in atlas['bisections'] if c['label']==label) for label in ('orbit-01333','orbit-0b2d0','orbit-19e45')]
    r.write_new(INPUT,{'schema':'rank-jump.native-triple-intersection-inputs.v1',
        **{k:base[k] for k in ('A','B','sections','parameters')},'covers':covers,
        'generic_word':list(map(str,s)),'gram':r.read(LATTICE)['gram'],
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (PROTOCOL,BASE,ATLAS,RELATIONS,TRANSPORT,LATTICE)},
        'boundary':r.read(PROTOCOL)['boundary']})


def compute():
    from sage.all import QQ,ZZ,PolynomialRing,EllipticCurve,matrix,vector
    inp=r.read(INPUT)
    for path,sha in inp['bindings'].items():assert r.digest((r.ROOT/path).read_bytes())==sha
    G=matrix(ZZ,inp['gram']);sv=vector(ZZ,[ZZ(x) for x in inp['generic_word']])
    traces=[vector(ZZ,c['published_basis_w']) for c in inp['covers']]
    z=2*sv+traces[1]-traces[2];rv=z-traces[0]
    assert all(w*G*w==10 for w in traces) and sv*G*sv==10 and rv*G*rv==10
    geometry={'translate_height':10,'residual_trace_vector':list(map(int,rv)),
              'residual_trace_height':10,'pair_image_trace':list(map(int,2*z)),
              'pair_image_fibre_degree':4,'pair_image_zero_intersection':int((z*G*z-4)/2),
              'pair_image_square':16,'pair_image_arithmetic_genus':9,'pair_carrier_genus':1,
              'predicted_proper_intersection_degree':12}
    print('Geometry gate',geometry,flush=True)
    R=PolynomialRing(QQ,'t');K=R.fraction_field();A=R(inp['A']);B=R(inp['B']);E=EllipticCurve(K,[A,B]);points=[]
    for c in inp['sections']:
        x=R(c['x_coefficients_low_to_high'])
        if 'y_coefficients_low_to_high' in c:y=R(c['y_coefficients_low_to_high'])
        else:
            ch=c['chord'];ref=points[ch['reference_basis_index']]
            y=ref[1]+R(ch['slope_coefficients_low_to_high'])*(x-ref[0])
        points.append(E(x,y))
    # A deterministic small-height ordering avoids large intermediate sections.
    terms=[(int(n/abs(n))*vector(ZZ,[int(j==i) for j in range(17)]),int(n/abs(n))*points[i]) for i,n in enumerate(sv) if n for _ in range(abs(int(n)))]
    acc=vector(ZZ,[0]*17);S=E(0)
    while terms:
        j=min(range(len(terms)),key=lambda j:((acc+terms[j][0])*G*(acc+terms[j][0]),j))
        v,P=terms.pop(j);acc+=v;S+=P
    assert acc==sv
    def enc(value):
        value=K(value);return {'numerator':list(map(str,value.numerator().list())),'denominator':list(map(str,value.denominator().list()))}
    a,b,c=inp['covers'];qs=[R(x['residual_chord']['q_coefficients']) for x in inp['covers']]
    assert all(q.is_squarefree() and q.degree()==2 for q in qs)
    assert all(qs[i].gcd(qs[j])==1 for i in range(3) for j in range(i))
    Z=PolynomialRing(K,'z');L=K.extension(Z.gen()**2-qs[1],'u');u=L.gen()
    V=PolynomialRing(L,'w');M=L.extension(V.gen()**2-L(qs[2]),'v');v=M.gen();EM=E.base_extend(M)
    def lift(c):return [R(c['lifted_section'][key+'_coefficients']) for key in ('x0','x1','y0','y1')]
    bx0,bx1,by0,by1=lift(b);cx0,cx1,cy0,cy1=lift(c)
    PB=EM(bx0+M(u)*bx1,by0+M(u)*by1);PC=EM(cx0+v*cx1,cy0+v*cy1)
    Q=EM(S)+PB-PC
    print('Biquadratic sum formed',flush=True)
    h=R(a['trace_section']['h_coefficients']);slope=K(R(a['residual_chord']['M_coefficients']))/h
    tx=K(R(a['trace_section']['Nx_coefficients']))/h**2;ty=K(R(a['trace_section']['Ny_coefficients']))/h**3
    intercept=-ty-slope*tx;residual=Q[1]-slope*Q[0]-intercept
    aa,bb=[L(residual[i]) for i in range(2)];assert bb
    vn=-aa/bb;first=aa*aa-L(qs[2])*bb*bb
    f0,f1=[K(first[i]) for i in range(2)];assert f1
    un=-f0/f1
    def evl(value):return K(L(value)[0])+un*K(L(value)[1])
    vr=evl(vn)
    def evm(value):return evl(M(value)[0])+vr*evl(M(value)[1])
    qx,qy=map(evm,(Q[0],Q[1]));ax0,ax1,ay0,ay1=lift(a);ar=(qx-ax0)/ax1
    f=R(a['quadratic_cover']['leading_coefficients'])*qx**2+R(a['quadratic_cover']['linear_coefficients'])*qx+R(a['quadratic_cover']['constant_coefficients'])
    norm=un*un-qs[1];common=norm.numerator().gcd(f.numerator()).monic()
    print('Norm degree',norm.numerator().degree(),'common degree',common.degree(),flush=True)
    functions=[un,vr,ar,qx,qy,slope,intercept,K(S[0]),K(S[1]),f0,f1,evl(bb)]
    for value in (Q[0],Q[1],residual):
        for j in range(2):
            for i in range(2):functions.append(K(L(M(value)[j])[i]))
    excluded=-16*(4*A**3+27*B**2)*h*ax1*f1.numerator()*evl(bb).numerator()
    for q in qs:excluded*=q
    for value in functions:excluded*=value.denominator()
    good=common
    while good.gcd(excluded).degree()>0:good=good//good.gcd(excluded)
    good=good.monic();assert good.degree()>0 and good.is_squarefree()
    def zero(value):
        value=K(value);assert value.denominator().gcd(good)==1 and value.numerator()%good==0
    zero(un*un-qs[1]);zero(vr*vr-qs[2]);zero(ar*ar-qs[0]);zero(qy-ay0-ar*ay1)
    zero(qy*qy-qx**3-A*qx-B)
    roots=[ar,un,vr];maps=[]
    for cover,root in zip(inp['covers'],roots,strict=True):
        x0,x1,y0,y1=lift(cover);maps.append((x0+root*x1,y0+root*y1))
    # Verify the relation independently in the finite parameter algebra.
    T=R.quotient(good,'tt')
    def red(f):f=K(f);return T(f.numerator())/T(f.denominator())
    def add(P,Q):
        x,y=P;xx,yy=Q;m=(yy-y)/(xx-x);xn=m*m-x-xx;return xn,m*(x-xn)-y
    pp=[tuple(map(red,P)) for P in maps]
    added=add(add(pp[0],(pp[1][0],-pp[1][1])),pp[2]);assert added==tuple(map(red,(S[0],S[1])))
    print('Finite algebra exact relation verified; factoring',flush=True)
    factors=[{'coefficients':list(map(str,f.list())),'multiplicity':int(e)} for f,e in good.factor()]
    checks=[]
    for row in inp['parameters']:
        t0=QQ(row['published_parameter']);op=bool(excluded(t0));hit=op and good(t0)==0
        check=row|{'on_certified_open':op,'on_relation_locus':bool(hit)}
        if hit:check['cover_roots']=list(map(str,[root(t0) for root in roots]))
        checks.append(check)
    return {'schema':'rank-jump.native-triple-intersection.v1','status':'PASS','layer':'solubility',
        'geometry':geometry,'generic_translate':{'x':enc(S[0]),'y':enc(S[1]),'word':inp['generic_word']},
        'norm_numerator_degree':int(norm.numerator().degree()),'unsaturated_common_degree':int(common.degree()),
        'intersection_degree':int(good.degree()),'intersection_polynomial':list(map(str,good.list())),
        'factorization':factors,'excluded_polynomial':list(map(str,excluded.squarefree_part().monic().list())),
        'rational_roots':[enc(root) for root in roots],
        'point_maps':[{'x':enc(x),'y':enc(y)} for x,y in maps],
        'parameter_checks':checks,
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (INPUT,PROTOCOL,Path(__file__),HERE/'retrospective.py')},
        'boundary':r.read(PROTOCOL)['boundary']}


def run():
    WORK.mkdir(parents=True,exist_ok=True);log=WORK/'worker.log';execution=WORK/'execution.json'
    if not log.exists():
        with log.open('x') as out:
            try:
                p=subprocess.run(['/home/royvanrijn/.local/bin/sage','-python',str(Path(__file__).resolve()),'worker'],stdout=out,stderr=out,timeout=60)
                result={'status':'COMPLETE' if p.returncode==0 else 'FAILED','returncode':p.returncode}
            except subprocess.TimeoutExpired:result={'status':'TIMEOUT'}
        r.write_new(execution,result)
    print(r.read(execution));print(log.read_text())


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['capture','run','worker','check']);a=p.parse_args()
    if a.mode=='capture':capture()
    elif a.mode=='run':run()
    else:
        data=compute()
        if a.mode=='worker':r.write_new(OUTPUT,data)
        else:assert data==r.read(OUTPUT)
        print('PASS intersection degree',data['intersection_degree'],'factor degrees',[len(x['coefficients'])-1 for x in data['factorization']],flush=True)
        print('Frozen hits',[x['source_id'] for x in data['parameter_checks'] if x['on_relation_locus']],flush=True)
