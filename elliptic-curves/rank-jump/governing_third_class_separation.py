#!/usr/bin/env python3
"""A finite-prime separation witness on a matched pair and historic control."""
import argparse
from pathlib import Path
import retrospective as r
import fresh_governing_panel as panel
import fresh_governing_octics as octics

PROTOCOL=Path(__file__).with_name('GOVERNING_THIRD_CLASS_SEPARATION_PROTOCOL.json')
OUTPUT=r.OUT/'rank_jump_governing_third_class_separation_v1.json'


def compute():
    from sage.all import QQ,GF,PolynomialRing
    data={x['token']:x for x in r.read(octics.OUTPUT)['rows']};R=PolynomialRing(QQ,'z');z=R.gen();rows=[]
    for token in ('case-02','case-03','case-13'):
        f,pts,scale=panel.model_data(token);x,y=pts[2];A=f[1]
        assert y*y==f(x) and len(pts)>=16
        q=z**4-3*x*z*z/2-y*z-3*x*x/16-A/4
        h=R(data[token]['integral_octic_ascending']);assert h.degree()==8
        witness=None;tested=0
        for p in r.primes(10009):
            try:
                P=PolynomialRing(GF(p),'z');hp=P(h);qp=P(q);fp=P(f.list());xx=GF(p)(x)
            except (ValueError,ZeroDivisionError):continue
            if not hp.discriminant() or not qp.discriminant() or not fp.discriminant():continue
            tested+=1
            hroots=hp.roots(multiplicities=False)
            if len(hroots)!=8:continue
            degrees=sorted(int(a.degree()) for a,e in qp.factor() for _ in range(e))
            if degrees!=[2,2]:continue
            roots=sorted(map(int,fp.roots(multiplicities=False)));assert len(roots)==3
            signs=[int(not (xx-root).is_square()) for root in roots]
            assert sum(signs)==2
            witness={'prime':p,'octic_roots':sorted(map(int,hroots)),
                'quartic_factor_degrees':degrees,'cubic_roots':roots,'third_kummer_signature':signs};break
        rows.append({'token':token,'status':'PASS' if witness else 'UNKNOWN','generic_index':2,
            'third_generic_point':[str(x),str(y)],'cubic_ascending':list(map(str,f.list())),
            'third_class_quartic_ascending':list(map(str,q.list())),
            'admissible_primes_tested':tested,'witness':witness,
            'boundary':'Third inherited class only. This proves field noncontainment; it does not discover an exceptional class.'})
    paths=(Path(__file__),PROTOCOL,panel.INPUT,octics.OUTPUT,Path(panel.__file__),Path(r.__file__))
    return {'schema':'rank-jump.governing-third-class-separation.v1','rows':rows,
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in paths},
        'scope':'Point-independent calibration using only previously marked generic sections on fixed retrospective fibres.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);args=p.parse_args();result=compute()
    if args.mode=='build':r.write_new(OUTPUT,result)
    else:assert result==r.read(OUTPUT)
    print([(x['token'],x['status'],x['witness']['prime'] if x['witness'] else None) for x in result['rows']])
