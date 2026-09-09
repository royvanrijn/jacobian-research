#!/usr/bin/env sage-python
"""Uniform conic-index gate and the three minimum stored-coordinate rows.

No new CVPs, orbit census or point search. Read the immutable quotient.
"""
import csv,hashlib,json,signal
from pathlib import Path
from sage.all import GF,ZZ,matrix,vector
signal.alarm(25)
ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'det1092_rational_bisection_index_v1'
PARENT=ART/'curve302_recovered_mw17_parent_v1.json'
LATTICE=ART/'curve302_parent_degree2_multisection_lattice_v1.json'
ORBITS=ART/'curve302_parent_degree2_multisection_orbits_v1.tsv'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(name,data):
    OUT.mkdir(exist_ok=True);p=OUT/name
    if p.exists():assert read(p)==data
    else:
        with p.open('x') as stream:json.dump(data,stream,indent=2,sort_keys=True);stream.write('\n')
save('protocol.json',{'classification':'generic-only conic-index deduction and bounded construction selection',
    'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [PARENT,LATTICE,ORBITS]},
    'script_sha256':sha(Path(__file__)),
    'limits':{'seconds':25,'new_CVPs':0,'new_orbit_enumerations':0,'point_searches':0,'selected_orbits_cap':3},
    'selection':'All rational rows tied for the smallest l1 of their stored parent-basis representative; fail if there are more than3. This is coordinate-dependent and does not assert minimum l1 across all representatives in an orbit. No specialized point or parameter enters selection.'})
p=read(PARENT);lattice=read(LATTICE)
assert lattice['status']=='PASS_COMPLETE_DEGREE2_TRANSLATION_QUOTIENT'
assert sha(ORBITS)==lattice['orbits_tsv_sha256']
G=matrix(ZZ,p['generic_height_gram']);assert G.det()==1092 and G.nrows()==17
assert all(v%2==0 for v in G.diagonal())
G2=matrix(GF(2),G);K=G2.right_kernel();assert K.dimension()==1
radical=vector(ZZ,[int(v) for v in K.basis()[0]])
norm=ZZ(radical*G*radical)
assert all(v%2==0 for v in G*radical) and norm%4==0
rows=[]
with ORBITS.open() as stream:
    for row in csv.DictReader(stream,delimiter='\t'):
        if row['category']!='rational':continue
        w=list(map(int,row['parent_MW17_w'].split()))
        assert int(row['minimum_norm'])==10
        rows.append({'orbit':int(row['orbit_mask']),'word':w,'stored_l1':sum(map(abs,w))})
assert len(rows)==40917
smallest=min(r['stored_l1'] for r in rows)
selected=sorted((r for r in rows if r['stored_l1']==smallest),key=lambda r:r['orbit'])
assert len(selected)<=3 and len(selected)==3
for row in selected:
    w=vector(ZZ,row['word']);assert w*G*w==10
    degrees=vector(ZZ,G.diagonal())-G*w
    assert all(v>=0 for v in degrees) and any(v%2 for v in degrees)
    row['section_intersection_degrees']=list(map(int,degrees))
    row['degree_one_sections']=[i for i,v in enumerate(degrees) if v==1]
    row['first_odd_section']=next(i for i,v in enumerate(degrees) if v%2)
report={'status':'PASS_UNIFORM_RATIONAL_BISECTION_INDEX_ONE',
    'classification':'new deduction from exact integral/mod2 Gram data and existing divisor theorem',
    'radical_dimension':1,'radical_word':list(map(int,radical)),'radical_word_norm':int(norm),
    'radical_norm_mod4':int(norm%4),'smallest_stored_l1':smallest,'selected':selected,
    'uniform_argument':'If G*w is even then w mod2 is0 or the radical word, so w^2 is0 mod4. A norm10 word therefore has an odd section intersection. The geometrically rational bisection has rational divisors of degrees2 and odd, hence index1 and a Q-point. All40917 known rational-bisection translation orbits are Q-rational, not merely geometrically rational.',
    'boundary':'Rationality of the whole conic does not imply rational splitting at any specified original elliptic parameter.',
    'inputs':read(OUT/'protocol.json')['inputs'],'protocol_sha256':sha(OUT/'protocol.json')}
save('index-and-selection.json',report)
print(report['status'],'radical norm',norm,'selected',[(r['orbit'],r['degree_one_sections']) for r in selected],flush=True)
