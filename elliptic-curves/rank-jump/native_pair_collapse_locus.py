#!/usr/bin/env python3
"""A fixed translated-bisection intersection, with no parameter search."""
import argparse
from pathlib import Path
import subprocess
import retrospective as r

HERE=Path(__file__).resolve().parent
PROTOCOL=HERE/'NATIVE_PAIR_COLLAPSE_LOCUS_PROTOCOL.json'
ATLAS=r.ROOT/'artifacts/generated-results/elkies-2026-equation-bisections-full.json'
MODEL=r.ROOT/'elkies-k3/data/fibrations/elkies_2026_published_r17_model.json'
SECTIONS=r.ROOT/'elkies-k3/data/fibrations/elkies_2026_published_r17_sections.json'
TRANSPORT=r.OUT/'compact_published_r17_generic_transport_v1.json'
RELATIONS=r.OUT/'rank_jump_paired_quartet_relations_verification_v1.json'
PARAMETERS=r.OUT/'rank_jump_solubility_first_inputs_v1.json'
INPUT=r.OUT/'rank_jump_native_pair_collapse_locus_inputs_v1.json'
OUTPUT=r.OUT/'rank_jump_native_pair_collapse_locus_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-native-pair-collapse-locus-v1'


def capture():
    from fractions import Fraction as F
    atlas=r.read(ATLAS);relations=r.read(RELATIONS);transport=r.read(TRANSPORT)
    old=next(row for row in relations['rows'] if row['id']=='08234-009')
    assert old['kernel_integer_vectors']==[[1,0,1,0]]
    word=[sum(F(v)*row[i] for v,row in zip(old['kernel_generic_coordinates'][0],transport['compact_sections_in_published_basis'],strict=True)) for i in range(17)]
    assert word==[0,0,0,0,0,0,0,0,0,0,0,1,0,-1,1,-1,0]
    covers=[next(c for c in atlas['bisections'] if c['label']==label) for label in ('orbit-0911e','orbit-1795d')]
    model=r.read(MODEL)
    r.write_new(INPUT,{'schema':'rank-jump.native-pair-collapse-locus-inputs.v1',
        'covers':covers,'generic_word':list(map(str,word)),'sections':r.read(SECTIONS)['sections'],
        'A':model['A_coefficients_low_to_high'],'B':model['B_coefficients_low_to_high'],
        'parameters':r.read(PARAMETERS)['parameters'],
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (PROTOCOL,ATLAS,MODEL,SECTIONS,TRANSPORT,RELATIONS,PARAMETERS)},
        'boundary':r.read(PROTOCOL)['boundary']})


def compute():
    from sage.all import QQ,PolynomialRing,EllipticCurve
    inp=r.read(INPUT)
    for path,sha in inp['bindings'].items():assert r.digest((r.ROOT/path).read_bytes())==sha
    R=PolynomialRing(QQ,'t');K=R.fraction_field();A=R(inp['A']);B=R(inp['B']);E=EllipticCurve(K,[A,B])
    points=[]
    for c in inp['sections']:
        x=R(c['x_coefficients_low_to_high'])
        if 'y_coefficients_low_to_high' in c:y=R(c['y_coefficients_low_to_high'])
        else:
            chord=c['chord'];ref=points[chord['reference_basis_index']]
            y=ref[1]+R(chord['slope_coefficients_low_to_high'])*(x-ref[0])
        points.append(E(x,y))
    S=sum((int(QQ(n))*P for n,P in zip(inp['generic_word'],points,strict=True)),E(0));assert not S.is_zero()
    c,d=inp['covers'];q=R(c['residual_chord']['q_coefficients']);q2=R(d['residual_chord']['q_coefficients'])
    Z=PolynomialRing(K,'z');L=K.extension(Z.gen()**2-q,'u');u=L.gen()
    def lift(c):
        return [R(c['lifted_section'][key+'_coefficients']) for key in ('x0','x1','y0','y1')]
    x0,x1,y0,y1=lift(c);EL=E.base_extend(L);P=EL(x0+u*x1,y0+u*y1);Q=EL(S)-P
    h=R(d['trace_section']['h_coefficients']);M=R(d['residual_chord']['M_coefficients'])
    tx=K(R(d['trace_section']['Nx_coefficients']))/h**2
    ty=K(R(d['trace_section']['Ny_coefficients']))/h**3
    m=K(M)/h;b=-ty-m*tx
    assert E(tx,ty)
    residual=Q[1]-m*Q[0]-b
    aa,bb=[K(residual[i]) for i in range(2)];assert bb
    norm=aa*aa-q*bb*bb;assert norm
    ur=-aa/bb
    xx=K(Q[0][0])+ur*K(Q[0][1]);yy=K(Q[1][0])+ur*K(Q[1][1])
    dx0,dx1,dy0,dy1=lift(d);vr=(xx-dx0)/dx1
    # The residual quadratic rejects the third point of the chord.
    f=R(d['quadratic_cover']['leading_coefficients'])*xx**2+R(d['quadratic_cover']['linear_coefficients'])*xx+R(d['quadratic_cover']['constant_coefficients'])
    common=norm.numerator().gcd(f.numerator()).monic()
    functions=[aa,bb,ur,xx,yy,vr,m,b,K(S[0]),K(S[1]),K(Q[0][0]),K(Q[0][1]),K(Q[1][0]),K(Q[1][1])]
    delta=-16*(4*A**3+27*B**2)
    excluded=delta*q*q2*h*dx1*bb.numerator()
    for value in functions:excluded*=value.denominator()
    removed=common.gcd(excluded)
    good=common
    while good.gcd(excluded).degree()>0:good=good//good.gcd(excluded)
    good=good.monic();assert good.degree()>0 and good.is_squarefree()
    # All identities in the reduced finite parameter algebra are exact.
    def zero(value):
        value=K(value);assert value.denominator().gcd(good)==1
        assert value.numerator()%good==0
    zero(ur**2-q);zero(vr**2-q2)
    zero(yy-dy0-vr*dy1);zero(yy**2-xx**3-A*xx-B)
    rx=x0+ur*x1;ry=y0+ur*y1
    zero(ry**2-rx**3-A*rx-B)
    # P+Q=S, without trusting extension-field group arithmetic.
    slope=(yy-ry)/(xx-rx)
    zero(slope*slope-rx-xx-S[0]);zero(slope*(rx-S[0])-ry-S[1])
    excluded*= (xx-rx).numerator()
    assert excluded.gcd(good)==1
    def enc(value):
        value=K(value);return {'numerator':list(map(str,value.numerator().list())),
                              'denominator':list(map(str,value.denominator().list()))}
    factors=[{'polynomial':list(map(str,f.list())),'multiplicity':int(e)} for f,e in good.factor()]
    rows=[]
    for par in inp['parameters']:
        t0=QQ(par['published_parameter']);on_open=excluded(t0)!=0
        hit=on_open and good(t0)==0
        row=par|{'on_certified_open':bool(on_open),'on_translated_intersection':bool(hit)}
        if hit:
            roots=[ur(t0),vr(t0)];P0=EllipticCurve(QQ,[A(t0),B(t0)])(rx(t0),ry(t0))
            Q0=P0.curve()(xx(t0),yy(t0));S0=P0.curve()(S[0](t0),S[1](t0))
            assert P0+Q0==S0
            row['cover_roots']=list(map(str,roots));row['exact_sum_equals_generic_translate']=True
        rows.append(row)
    return {'schema':'rank-jump.native-pair-collapse-locus.v1','status':'PASS','layer':'solubility',
        'generic_translate':{'word':inp['generic_word'],'x':enc(S[0]),'y':enc(S[1])},
        'chord_residual':{'constant':enc(aa),'root_coefficient':enc(bb)},
        'norm_numerator_degree':int(norm.numerator().degree()),'unsaturated_common_degree':int(common.degree()),
        'removed_degree':int(removed.degree()),'intersection_polynomial':list(map(str,good.list())),
        'intersection_degree':int(good.degree()),'factorization':factors,
        'excluded_polynomial':list(map(str,excluded.squarefree_part().monic().list())),
        'rational_maps':{key:enc(value) for key,value in [('u',ur),('v',vr),('px',rx),('py',ry),('qx',xx),('qy',yy)]},
        'parameter_checks':rows,
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (INPUT,PROTOCOL,Path(__file__),HERE/'retrospective.py')},
        'boundary':r.read(PROTOCOL)['boundary']}


def run():
    WORK.mkdir(parents=True,exist_ok=True);log=WORK/'worker.log';execution=WORK/'execution.json'
    if not log.exists():
        with log.open('x') as out:
            try:
                result=subprocess.run(['/home/royvanrijn/.local/bin/sage','-python',str(Path(__file__).resolve()),'worker'],stdout=out,stderr=out,timeout=60)
                status={'status':'COMPLETE' if result.returncode==0 else 'FAILED','returncode':result.returncode}
            except subprocess.TimeoutExpired:status={'status':'TIMEOUT'}
        r.write_new(execution,status)
    print(r.read(execution));print(log.read_text())


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['capture','run','worker','check']);a=p.parse_args()
    if a.mode=='capture':capture()
    elif a.mode=='run':run()
    else:
        data=compute()
        if a.mode=='worker':r.write_new(OUTPUT,data)
        else:assert r.read(OUTPUT)==data
        print('PASS intersection degree',data['intersection_degree'],'factor degrees',[len(f['polynomial'])-1 for f in data['factorization']],flush=True)
        print('Frozen hits',[p['source_id'] for p in data['parameter_checks'] if p['on_translated_intersection']],flush=True)
