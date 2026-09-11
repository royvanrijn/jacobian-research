#!/usr/bin/env sage-python
"""Eight old-R17 fibres as genus-one carriers over the compact11952 base.

Matched alternate-fibration supplement, using only the retained direct degree2
map. No historical degree11511 transport, new fibration, or point search.
"""
import importlib.machinery,importlib.util
from pathlib import Path
from sage.all import QQ,PolynomialRing,EllipticCurve
path=Path(__file__).with_name('rank_triangle_geometry.sage')
loader=importlib.machinery.SourceFileLoader('triangle_geometry',str(path));spec=importlib.util.spec_from_loader(loader.name,loader);g=importlib.util.module_from_spec(spec);loader.exec_module(g)
R,F,t=g.R,g.F,g.t;ART,ROOT,OUT=g.ART,g.ROOT,g.OUT
def main():
    source=ROOT/'artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json';atlas=ART/'compact_six_r17_atlas_v1.json'
    literal=g.read(source);compact=next(r for r in g.read(atlas)['families'] if r['family']=='11952')
    assert g.sha(source)==compact['source_sha256']
    g.save(OUT/'native-carrier-protocol.json',{'sources':{str(p):g.sha(p) for p in [Path(__file__),path,source,atlas,OUT/'geometry.json']},
      'limits':{'seconds':180,'carriers':8,'workers':1},'purpose':'Balance the302 alternate-fibration atlas with the existing direct degree2 old-R17/11952 fibration map. Preserve the original matched chord atlas.'})
    a,b,c,d=map(QQ,compact['base_matrix_a_b_c_d']);scale=QQ(compact['total_scale_from_literal_source']);u=F((a*t+b)/(c*t+d))
    substitute=lambda rec:g.dec(rec)(u)
    gm=literal['genus_one_model'];coeff=[substitute(r) for r in gm['q_coefficients_in_t_low_to_high']]
    tzero=substitute(gm['distinguished_point_from_old_zero']['t0']);vzero=substitute(gm['distinguished_point_from_old_zero']['W0'])
    gauge=substitute(literal['weierstrass_model']['gauge'])
    Z=PolynomialRing(F,'z');z=Z.gen();qpoly=Z(coeff);shift=qpoly(z+tzero);ee,dd,cc,bb,aa=[shift[i] for i in range(5)]
    assert ee==vzero*vzero
    a1=dd/vzero;a2=cc-dd*dd/(4*vzero*vzero);a3=2*vzero*bb;a4=-4*vzero*vzero*aa;a6=a2*a4;b2=a1*a1+4*a2
    A,B,sections,tstar=g.parent('11952')
    factor=F(3*gauge*(c*t+d)**2/scale)
    # This global model identity pins every coordinate scaling before targets.
    Egen=EllipticCurve(F,[a1,a2,a3,a4,a6])
    assert -Egen.c4()/48*factor**4==A and -Egen.c6()/864*factor**6==B
    rows=[]
    for original in [r for r in g.read(OUT/'geometry.json')['rows'] if r['fibre']=='11952']:
        px,py=map(QQ,original['parent_point']);at=lambda f:g.at(f,tstar)
        xg=px/at(factor)**2-at(b2)/12;yg=py/at(factor)**3-(at(a1)*xg+at(a3))/2
        old_t=at(tzero)+(2*at(vzero)*(xg+at(cc))-at(dd)**2/(2*at(vzero)))/yg
        zs=old_t-at(tzero);W=(xg*zs**2-at(dd)*zs-2*at(vzero)**2)/(2*at(vzero))
        raw=F(qpoly(old_t));assert at(raw)==W*W
        zc=F(old_t-tzero)
        x0=(2*vzero*vzero+dd*zc)/zc**2;x1=2*vzero/zc**2
        y0=(4*vzero**3+2*vzero*dd*zc+(2*vzero*cc-dd*dd/(2*vzero))*zc**2)/zc**3;y1=4*vzero*vzero/zc**3
        y0=(y0+(a1*x0+a3)/2)*factor**3;y1=(y1+a1*x1/2)*factor**3
        x0=(x0+b2/12)*factor**2;x1=x1*factor**2
        assert y0*y0+y1*y1*raw==x0**3+3*x0*x1*x1*raw+A*x0+B
        assert 2*y0*y1==3*x0*x0*x1+x1**3*raw+A*x1
        assert at(x0)+at(x1)*W==px and at(y0)+at(y1)*W==py
        twist,q,s=g.squarefree_cover(raw);v=W/at(s);assert v*v==twist*q(tstar)
        rec={'target':original['target'],'fibre':'11952','old_R17_fibre_parameter':str(old_t),
             'tag':{'family':'old_R17_fibre'},'degree_over_base':2,'normalization_genus':int((q.degree()-1)//2),
             'branch_monic_polynomial':list(map(str,q.list())),'branch_at_infinity':bool(q.degree()%2),
             'constant_twist_representative':str(twist),'square_multiplier':g.enc(s),'raw_radical':g.enc(raw),
             'lift':[str(tstar),str(v)],'maps':{k:g.enc(f) for k,f in zip(['x0','x1','y0','y1'],[x0,x1*s,y0,y1*s])}}
        assert rec['normalization_genus']==1
        g.save(OUT/'native-carriers'/f'E{original["target"]}.json',rec);rows.append(rec)
        print('NATIVE GENUS1',original['target'],flush=True)
    assert len({r['old_R17_fibre_parameter'] for r in rows})==8
    assert len({tuple(r['branch_monic_polynomial']) for r in rows})==8
    g.save(OUT/'native-carriers.json',{'status':'PASS_EIGHT_DISTINCT_GENUS_ONE_CARRIERS','rows':rows,
      'boundary':'Eight distinct covers in one alternate-fibration pencil, not eight sections on one common quadratic base change. Global genus0 ancestry remains UNKNOWN.'})
if __name__=='__main__':main()
