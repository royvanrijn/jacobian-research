#!/usr/bin/env sage-python
"""Oracle-conditioned nearest words over every extension of V2's 32 anchors."""
from pathlib import Path
from fractions import Fraction as F
from importlib.machinery import SourceFileLoader
from sage.all import matrix, vector, ZZ, QQ, pari
from research_runtime.store import checkpoint
from terminal_affine_cvp import AffineCVP
from visibility_lattice_v2 import ExactParity

CAS=Path(__file__).resolve().parent
l=SourceFileLoader('bank_landscape',str(CAS/'v2_terminal_landscape.sage')).load_module()
b=l.b;D=l.D
fixed=b.read(D/'fixed-M30.json');data=b.read(D/'coordinates.json')
points=tuple(tuple(map(F,p)) for p in fixed['points']);target=tuple(map(F,data['missing_point']));model=tuple(map(F,fixed['curve']))
g,u,inv,t,reduced,constant=l.setup(points,target,model)
diagonal=matrix.diagonal(ZZ,[2]*17+[1]*13)
bank_g=diagonal*g*diagonal
v=matrix(ZZ,pari(bank_g).qflllgram()).transpose()
bank_reduced=v*bank_g*v.transpose()
oracle=AffineCVP(bank_reduced.rows());centre_oracle=ExactParity(reduced.rows())
selection=b.read(b.V2/'calibration302/epoch-13/selection.json')
out=D/'retained-bank';out.mkdir(exist_ok=True)
for index,anchor in enumerate(selection['anchors']):
    path=out/f'anchor-{index:02d}.json'
    if path.exists():continue
    a=vector(ZZ,[int(x)%2 for x in anchor['anchor']['representative'][:17]]+[0]*13)
    projection=(t*u-a)*diagonal.inverse()*v.inverse()
    nearest=oracle.nearest(projection,count=16,node_limit=30000000)
    candidates=[];seen=set()
    for row in nearest['rows']:
        q=vector(ZZ,row['word'])*v*diagonal+a
        p=tuple(int(x)%2 for x in q)
        if p in seen:continue
        seen.add(p)
        rp=[int(x)%2 for x in vector(ZZ,p)*inv]
        seed,_=centre_oracle.babai([rp]);proof=centre_oracle.solve(rp,seed[0],2000000)
        choices={}
        for minimum in proof['minima']:
            word=list(map(int,vector(ZZ,minimum)*u));point=l.group.linear_combination(model,points,word)
            if point[1]<0:word=[-x for x in word];point=(point[0],-point[1])
            choices[point]=word
        witnesses=[l.evaluate(model,points,target,word,list(map(int,q))) for point,word in sorted(choices.items())]
        delta=vector(QQ,q)*inv-t
        candidates.append({'parity':list(p),'extension':sum(p[17+i]<<i for i in range(13)),
            'half_distance_numerator':str(delta*reduced*delta+constant),
            'centre_CVP':proof,'witnesses':witnesses,'best':min(witnesses,key=lambda r:int(r['height']))})
    result={'anchor_index':index,'anchor':anchor['anchor'],'extension_classes_covered_by_domain':8192,
        'nearest':nearest,'candidates':candidates,'best':min((r['best'] for r in candidates),key=lambda r:int(r['height']))}
    checkpoint(path,result)
    print('BANK',index,'best',result['best']['height'],result['best']['coordinate'],flush=True)
