#!/usr/bin/env sage-python
"""Global nonzero-translation degree minima from the norm-four shell.

PARI is a candidate producer only. Exact shell completeness is replayed
separately. No parameter or exceptional point input. External 25-second cap.
"""
import hashlib,json
from pathlib import Path
from sage.all import QQ,ZZ,matrix,vector,pari
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'det1092_two_fibration_action_v1';PARENT=ART/'curve302_recovered_mw17_parent_v1.json'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(name,data):
    path=OUT/name
    if path.exists():assert read(path)==data
    else:
        with path.open('x') as s:json.dump(data,s,indent=2,sort_keys=True);s.write('\n')
paths=[PARENT,OUT/'degree-minima-v2-minus.json',OUT/'degree-minima-v2-plus.json',Path(__file__)]
save('nontrivial-gap-protocol.json',{'classification':'generic-only exact degree obstruction',
    'rule':'Use the complete original norm-four shell to minimize each quadratic degree formula. Prove by an exact Cauchy bound that norms at least six cannot improve the resulting witness. No change of initial source or final alternate translating section.',
    'limits':{'seconds':25,'shell_norm':4,'maximum_stored_sign_pairs':5000,'parameter_evaluations':0,'point_searches':0},
    'inputs':{str(p.relative_to(ROOT)):sha(p) for p in paths}})
G=matrix(ZZ,read(PARENT)['generic_height_gram'])
raw=pari(G).qfminim(4,5000,2);V=matrix(ZZ,raw[2]);assert int(raw[0])==2*V.ncols()
words=sorted(set(tuple(s*x for x in v) for v in V.columns() for s in [-1,1]))
assert len(words)==int(raw[0]) and all(vector(ZZ,q)*G*vector(ZZ,q)==4 for q in words)
save('norm-four-shell.json',{'status':'EXACT_WORDS_COMPLETENESS_PENDING_REPLAY','count':len(words),'words':[list(map(int,w)) for w in words]})
reports=[]
for sign,path in zip([-1,1],paths[1:3]):
    d=read(path);center=vector(QQ,d['center']);c=QQ(d['quadratic_coefficient']);b=QQ(d['constant']);R=center*G*center
    values=[b+c*((vector(QQ,w)-center)*G*(vector(QQ,w)-center)) for w in words]
    best=min(values);winners=[list(map(int,w)) for w,val in zip(words,values) if val==best]
    # For n >= 6 > R, (sqrt(n)-sqrt(R))^2 increases. At n=6,
    # b+c*(6+R-2sqrt(6R)) > best iff the following rational gate holds.
    gap=6+R-(best-b)/c
    assert gap>0 and gap*gap>24*R
    reports.append({'sign':sign,'minimum_nonzero_translation_degree':int(best),
        'number_of_minimizers':len(winners),'minimizers':winners,
        'center_norm':str(R),'cauchy_gap':str(gap),'cauchy_squared_margin':str(gap*gap-24*R),
        'zero_translation_degree':int(b+c*R)})
save('nontrivial-gap.json',{'status':'EXACT_GLOBAL_DEGREE_MINIMA_PENDING_SHELL_REPLAY',
    'norm_four_count':len(words),'branches':reports,
    'proof':'The rootless elliptic K3 has even generic heights with every nonzero height at least four. Exhaust norm four. Cauchy excludes every height at least six using the displayed strictly positive rational margin.',
    'boundary':'Only the two indicated switching families; no high-degree incidence or new section rank assertion.'})
for row in reports:print(row['sign'],'nonzero minimum',row['minimum_nonzero_translation_degree'],'multiplicity',row['number_of_minimizers'])
