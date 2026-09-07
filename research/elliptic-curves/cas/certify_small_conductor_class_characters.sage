#!/usr/bin/env sage-python
"""Certify16 everywhere-unramified quadratic characters and independent ideal anchors."""
import argparse
from pathlib import Path
import runpy
from sage.all import AA,QQ,GF,PolynomialRing,matrix,pari
import extend_small_conductor_norm_batch as batch
from research_runtime.store import checkpoint

ROOT,ART,cert=batch.ROOT,batch.ART,batch.cert
SOURCE=Path(__file__).resolve()
LOWER_SOURCE=ROOT/'elliptic-curves/cas/certify_small_conductor_class_lower16.sage'
LOWER=ART/'small_conductor_class_lower16_v1.json'
BASE=ART/'small_conductor_norm_batch_relations_v1.json'
OUT=ART/'small_conductor_class_characters_v1.json'


def expected():
    low=runpy.run_path(str(LOWER_SOURCE))
    lower=low['expected']()
    if lower!=cert.read(LOWER):raise ArithmeticError('lower-bound replay differs')
    old=runpy.run_path(str(low['OLD_SOURCE']))
    data,F,nf,theta,w=old['setup']();proof=cert.read(low['PROOF'])
    points=[(QQ(x),QQ(y)) for x,y in proof['integral_points']]
    bs=[pari(4*x)-theta for x,y in points]
    K=matrix(GF(2),lower['even_valuation_kernel_basis'])
    local=[]; local_records=[]
    for P in pari.idealprimedec(nf,2):
        e,f=int(P[2]),int(P[3])
        if f!=1:raise ArithmeticError('this unit-filtration certificate needs residue field F2')
        pi=pari(2) if e==1 else pari.nfbasistoalg(nf,P[1])
        if pari.idealval(nf,pi,P)!=1:raise ArithmeticError('not a uniformizer')
        for k in range(1,2*e+1):
            unit=1+pi**k
            if pari.idealval(nf,unit,P)!=0:raise ArithmeticError('not a local unit')
            raw=[int(pari.nfhilbert(nf,b,unit,P)==-1) for b in bs]
            v=matrix(GF(2),[raw])*K.transpose()
            row=[int(v[0,j]) for j in range(18)]+[int(pari.nfhilbert(nf,-1,unit,P)==-1)]
            local.append(row)
            local_records.append({'p':2,'hnf':str(pari.idealhnf(nf,P)),'e':e,'f':f,'k':k,
                                  'uniformizer_power_basis':[str(pari.polcoef(pari.lift(pi),i)) for i in range(3)],
                                  'hilbert_bits_on_beta':raw,'constraint_on_19_classes':row})
    R=PolynomialRing(QQ,'t');roots=R([int(pari.polcoef(F,i)) for i in range(4)]).roots(AA,multiplicities=False)
    signs=[]
    for root in roots:
        raw=[int(4*x-root<0) for x,y in points]
        v=matrix(GF(2),[raw])*K.transpose()
        row=[int(v[0,j]) for j in range(18)]+[1]
        signs.append({'beta_sign_bits':raw,'constraint_on_19_classes':row});local.append(row)
    L=matrix(GF(2),local);H=L.right_kernel().basis_matrix()
    if L.rank()!=3 or H.nrows()!=16 or H.rank()!=16 or L*H.transpose()!=0:
        raise ArithmeticError('unramified character kernel differs')
    G=(H.matrix_from_columns(list(range(18)))*K).augment(H.matrix_from_columns([18]))
    if G.rank()!=16:raise ArithmeticError('unramified generators dependent')
    # beta_1,...,beta_22,-1 are independent, by the prior residue matrix and norm sign.
    base=cert.read(BASE);canonical={max(r['parity_columns']) for r in base['canonical_rational_relations']}
    anchors=[];labels=[]
    a=int(data['integral_norm_generator']['fixed_a'])
    for i,col in enumerate(base['columns']):
        q=col['p']
        if i in canonical or not 100<=q<=37638 or col['e']!=1 or col['f']!=1:
            continue
        if int(nf.disc())%q==0 or int(data['defining_order_index'])%q==0 or a%q==0:
            continue
        if not pari.isprime(q):raise ArithmeticError('anchor rational prime not proved')
        root=next(int(pari.lift(r)) for r in pari.polrootsmod(F,q)
                  if str(pari.idealhnf(nf,q,theta-int(pari.lift(r))))==col['hnf'])
        residues=[]
        for x,y in points:
            z=4*x
            if int(z.denominator())%q==0:break
            residues.append((int(z.numerator())*pow(int(z.denominator()),-1,q)-root)%q)
        if len(residues)!=22 or 0 in residues:continue
        raw=[int(pow(v,(q-1)//2,q)==q-1) for v in residues]+[int(q%4==3)]
        label=[int(z) for z in (matrix(GF(2),[raw])*G.transpose()).row(0)]
        if matrix(GF(2),labels+[label]).rank()>len(labels):
            labels.append(label)
            anchors.append({'column':i,'p':q,'hnf':col['hnf'],'theta_root':root,
                            'beta_residues':residues,'frobenius_character_bits':label})
        if len(labels)==16:break
    if len(anchors)!=16:raise ArithmeticError('independent ideal anchors not found')
    return {'schema':'elliptic-curves.small-conductor-class-characters.v1','status':'PASS',
        'sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in [SOURCE,LOWER_SOURCE,LOWER,BASE,batch.forms.OUT]},
        'dyadic_unit_constraints':local_records,'real_sign_constraints':signs,'local_constraint_rank':3,
        'unramified_basis_in_19_classes':[[int(z) for z in row] for row in H.rows()],
        'unramified_basis_in_beta_and_minus_one':[[int(z) for z in row] for row in G.rows()],
        'independent_unramified_characters':16,'anchors':anchors,'anchor_character_rank':16,
        'unit_generation_argument':'For each dyadic completion, residue field is F2 and v(2)=e. The units1+pi^k for1<=k<=2e generate O_v^*/O_v^{*2}: each successive U^k/U^{k+1} is F2, and U^{2e+1} consists of squares by Hensel. Vanishing Hilbert symbols with all these units makes the quadratic character unramified at that dyadic prime.',
        'global_argument':'The prior certificate supplies19 independent field Selmer classes, hence even valuations everywhere. Impose dyadic unit Hilbert symbols and positive signs at every real embedding; the resulting16-dimensional kernel consists of everywhere-unramified quadratic extensions. By global class field theory these give characters of the ordinary ideal class group, trivial on principal ideals. Their residue characters on the16 listed ideal classes form an invertible matrix. Those classes are independent in Cl(K)/2 unconditionally.',
        'claim_boundary':'Unconditional16 independent ideal classes and quadratic characters. This does not provide a class-group upper bound. They can be protected during row elimination without sacrificing any valid principal relation.'}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--check',action='store_true');a=p.parse_args()
    result=expected()
    if a.check:
        if cert.read(OUT)!=result:raise ArithmeticError('character certificate differs')
    else:
        if OUT.exists():raise FileExistsError('preserve certificate')
        checkpoint(OUT,result)
    print('16 INDEPENDENT UNRAMIFIED CHARACTERS AND IDEAL ANCHORS: PASS',flush=True)
