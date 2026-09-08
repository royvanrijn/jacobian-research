#!/usr/bin/env sage-python
"""Oracle-only generic translations and exact span audit of retained MW18 boxes.

Floating canonical heights propose one nearest generic translate for each sign,
point and retained centre. Every proposed group word, curve point and coordinate
is checked exactly. No optimality or complete-coset claim is made.
"""
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import sys
import zipfile
from collections import Counter
from math import log2

from sage.all import EllipticCurve, QQ, ZZ, matrix, vector, RealField, pari
from fpylll import GSO, IntegerMatrix, Enumeration

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'elliptic-curves/cas'))
from audit_mw18_retained_visibility import ART,BUNDLE,PUBLIC,transport
from elkies_rank28 import GENERAL_WEIERSTRASS_COEFFICIENTS,POINTS
from search_observability import point_visibility
from research_runtime.store import checkpoint,digest

LOCAL=ROOT/'artifacts/local/mw18-translated-visibility-v1'


def affine(p):
    return None if p.is_zero() else dict(x=str(p[0]),y=str(p[1]))


def combination(E,points,word):
    return sum((int(k)*p for k,p in zip(word,points) if k),E(0))


def run():
    z=zipfile.ZipFile(BUNDLE); protocol=json.loads(z.read('protocol.json'))
    rawaudit=json.loads((ART/'elliptic-curves/mw18_retained_visibility_v1.json').read_text())
    public={r['id']:r for r in json.loads(PUBLIC.read_text())['snapshot']['records']}
    result=dict(schema='elliptic-curves.mw18-translated-visibility.v1',cases={},retrospective_only=True,
        inputs={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in
            (BUNDLE,PUBLIC,ART/'elliptic-curves/mw18_retained_visibility_v1.json',Path(__file__))},
        translation_policy='One floating-CVP proposal for each signed representative and retained centre, '
        'using the specialized canonical Gram rounded after multiplication by 2^20. '
        'All group identities and chart visibility are exact; no CVP optimality claim.',
        claim_boundary='Pointwise witness minima over the declared proposals only. '
        'Exact rational-span membership in the displayed public subgroup, not saturation or full Mordell-Weil groups.')
    for case_id,case in sorted(protocol['cases'].items()):
        cp=LOCAL/digest(result['inputs'])/f'{case_id}.json'
        if cp.exists():
            result['cases'][case_id]=json.loads(cp.read_text());continue
        if case_id.startswith('historical'):model,points=GENERAL_WEIERSTRASS_COEFFICIENTS,POINTS
        else:
            row=public[int(case_id.split('-')[1])];model,points=row['ainvs'],row['points']
        E=EllipticCurve(QQ,list(map(QQ,model))); pub=[E(list(map(QQ,p))) for p in points]
        ER=EllipticCurve(QQ,case['curve'])
        # Choose the isomorphism orientation by an exact public-point match to
        # the Sage-free audit, so signs/centre words agree across implementations.
        expected=transport(model,points,case['curve'])[0]
        phi=next(f for f in E.isomorphisms(ER) if affine(f(pub[0]))==expected)
        inv=~phi
        generic=[inv(ER(QQ(p['x']),QQ(p['y']))) for p in case['points']]
        h=E.height_pairing_matrix(pub+generic,precision=192)
        hp=h[:28,:28]; proposed=hp.inverse()*h[:28,28:]
        words=matrix(ZZ,[[ZZ(x.round()) for x in r] for r in proposed.rows()])
        if any(combination(E,pub,words.column(j))!=p for j,p in enumerate(generic)):
            raise ArithmeticError('numerically proposed public/generic relation failed exact group law')
        assert words.rank()==18
        gm=words.transpose()*hp*words
        scaled=matrix(ZZ,[[ZZ((x*2**20).round()) for x in r] for r in gm.rows()])
        u=matrix(ZZ,pari(scaled).qflllgram()).transpose(); assert abs(u.det())==1
        reduced=u*scaled*u.transpose(); gso=GSO.Mat(IntegerMatrix.from_matrix([list(map(int,r)) for r in reduced.rows()]),gram=True,float_type='dd',update=True)
        mu=matrix(RealField(53),[[gso.get_mu(i,j) if i>j else int(i==j) for j in range(18)] for i in range(18)])
        projections=gm.inverse()*words.transpose()*hp
        inverse_u=u.inverse(); labels=list(rawaudit['cases'][case_id]['oracle'])
        memo={}; public_words={}; policies={}
        for policy in protocol['policies']:
            cell=json.loads(z.read(f'anchor-trial/cells/{case_id}--{policy}.json'))
            recovered={tuple(p.values()):p for c in cell['charts'] for p in c['search']['finite_curve_points']}
            rw=[]; replayed=0
            for i,ch in enumerate(cell['charts']):
                for p in ch['search']['finite_curve_points']:
                    audit=point_visibility(ch['search'],p)
                    if audit['status'] not in ('VISIBLE_AND_RECORDED','KNOWN_POINTED_ENDPOINT'):
                        raise ArithmeticError('retained discovery fails independent coverage replay')
                    replayed+=1
            for key,p in recovered.items():
                if key not in public_words:
                    point=inv(ER(QQ(p['x']),QQ(p['y'])))
                    hnew=E.height_pairing_matrix(pub+[point],precision=192)
                    v=(hp.inverse()*hnew[:28,28:29]).column(0)
                    w=vector(ZZ,[ZZ(x.round()) for x in v])
                    if combination(E,pub,w)!=point:raise ArithmeticError('recovered/public group relation failed')
                    public_words[key]=list(map(int,w))
                rw.append(public_words[key])
            span=words.augment(matrix(QQ,rw).transpose()) if rw else words
            gain=int(span.rank()-18)
            if gain!=cell['certified_gain']:raise ArithmeticError('rational-span gain differs from retained independence')
            minima={}; counts=Counter(); discrepancies=[]
            for label in labels:
                index=int(label[1:])-1; unit=vector(ZZ,[int(i==index) for i in range(28)])
                known=span.augment(matrix(QQ,28,1,list(unit))).rank()==span.rank()
                best=None
                for i,ch in enumerate(cell['charts']):
                    centre=vector(QQ,ch['centre']['coefficients'])
                    for sign in (1,-1):
                        target=(centre/2-sign*projections.column(index))*inverse_u
                        gs_target=tuple(map(float,target*mu))
                        # A rounded vector supplies a finite, valid initial radius.
                        trial=vector(ZZ,[round(float(x)) for x in target]); delta=trial-target
                        radius=float(delta*reduced*delta)+1
                        _,cvp=Enumeration(gso).enumerate(0,18,radius,0,target=gs_target)[0]
                        word=vector(ZZ,[round(x) for x in cvp])*u
                        full=words*word+sign*unit; key=tuple(full)
                        if key not in memo:memo[key]=affine(phi(combination(E,pub,full)))
                        p=memo[key]
                        v=point_visibility(ch['search'],p);counts[v['status']]+=1
                        if v['status']=='VISIBLE_NOT_RECORDED':discrepancies.append(dict(point=label,chart=i,sign=sign,**v))
                        height=v.get('minimum_affine_height')
                        if height is not None and (best is None or height<best['minimum_affine_height']):
                            best=dict(chart=i,sign=sign,translation_in_generic_basis=list(map(int,word)),
                                translated_point=p,**v)
                if best is None:
                    if not known:raise ArithmeticError('unexplained endpoint-only direction')
                    best=dict(status='ONLY_KNOWN_POINTED_ENDPOINTS',minimum_affine_height=None)
                best['already_in_recovered_rational_span']=bool(known)
                best['height_ratio']=str(Fraction(best['minimum_affine_height'],100000)) if best['minimum_affine_height'] else None
                best['log2_height']=log2(best['minimum_affine_height']) if best['minimum_affine_height'] else None
                minima[label]=best
            witnessed_gains={}
            for bound in (100000,200000,1000000,10000000):
                eligible=[int(k[1:])-1 for k,v in minima.items()
                    if v['minimum_affine_height'] is not None and v['minimum_affine_height']<=bound]
                extension=span
                for index in eligible:
                    extension=extension.augment(matrix(QQ,28,1,[int(i==index) for i in range(28)]))
                witnessed_gains[str(bound)]=int(extension.rank()-span.rank())
            policies[policy]=dict(recovered_public_basis_words=rw,exact_rational_span_gain=gain,
                recovered_points_coverage_replays=replayed,minima=minima,status_counts=dict(counts),
                discrepancies=discrepancies,witnessed_additional_rational_span_at_height=witnessed_gains,
                missing_basis_representatives=[k for k,v in minima.items() if not v['already_in_recovered_rational_span']])
            missing=[v['minimum_affine_height'] for v in minima.values() if not v['already_in_recovered_rational_span']]
            print('TRANSLATED',case_id,policy,'gain',gain,'statuses',dict(counts),'missing_heights',missing,flush=True)
        cr=dict(generic_words_in_public_basis=[list(map(int,r)) for r in words.rows()],policies=policies)
        checkpoint(cp,cr);result['cases'][case_id]=cr
    result['status']='PASS_NO_COVERAGE_DISCREPANCY' if not any(p['discrepancies'] for c in result['cases'].values() for p in c['policies'].values()) else 'COVERAGE_DISCREPANCY'
    checkpoint(ART/'elliptic-curves/mw18_translated_visibility_v1.json',result)


if __name__=='__main__':run()
