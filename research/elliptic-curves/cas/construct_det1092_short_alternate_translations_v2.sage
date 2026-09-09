#!/usr/bin/env sage-python
"""Fixed conic translated by the intrinsically shortest generic sections.

Only the completed130 degree1 original source curves, not a new lattice
census. Freeze before geometry; no target in geometry or map construction.
Each action25s. At most8 tied minimum-degree new image curves may be built.
"""
import argparse, hashlib, json, signal, time
from pathlib import Path
from sage.all import QQ, ZZ, PolynomialRing, EllipticCurve, matrix, vector, block_diagonal_matrix, identity_matrix

ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'det1092_short_alternate_translations_v1'
PARENT=ART/'curve302_recovered_mw17_parent_v1.json'
FRAME=ART/'det1092_genus1_picard_image_v1/frame.json'
SOURCES=ART/'det1092_low_degree_source_obstruction_v1/sources.json'
PENCIL=ART/'det1092_norm8_seed_cover_v2/generic.json'
CONIC=ART/'det1092_bifibration_conic_sources_v1/map-47755.json'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(name,data):
    path=OUT/name;path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists():assert read(path)==data,'immutable checkpoint changed'
    else:
        with path.open('x') as f:json.dump(data,f,indent=2,sort_keys=True);f.write('\n')
def ints(v):return list(map(int,v))
def freeze():
    save('protocol-v2.json',dict(classification='generic-only shortest alternate section construction',
        source_conic=47755,section_height='least positive among the completed130 common sections',existing_source_count=130,
        rule='From the completed original-section source set of alternate degree1, take every section of least positive alternate height. Translate the fixed bidegree2/2 conic once by each sign of each such section. Keep all distinct new image classes of minimum original degree greater than1. Build the whole tied minimum set only if at most8; otherwise retain UNKNOWN_MAP_CAP and do not choose orbit IDs.',
        limits=dict(seconds_per_action=25,maps=8,point_searches=0,parameter_sweeps=0,new_source_enumerations=0),
        inputs={str(p.relative_to(ROOT)):sha(p) for p in [PARENT,FRAME,SOURCES,PENCIL,CONIC,Path(__file__),
            ART/'det1092_low_degree_source_obstruction_v1/independent-replay.json',
            ART/'det1092_bifibration_conic_sources_v1/independent-replay.json']}))
def protocol():
    p=read(OUT/'protocol-v2.json')
    for name,h in p['inputs'].items():assert sha(ROOT/name)==h
    return p
def lattice():
    p=protocol();parent=read(PARENT);frame=read(FRAME)
    G=matrix(QQ,parent['generic_height_gram']);NS=block_diagonal_matrix(matrix(QQ,[[-2,1],[1,0]]),-G)
    D=vector(QQ,frame['D']);O=vector(QQ,frame['new_zero_section']);F=vector(QQ,[0,1]+[0]*17)
    def dot(v,w):return v*NS*w
    vertical=[vector(QQ,row['NS_coordinates']) for row in frame['vertical_old_sections']]
    roots=[v for v in vertical if dot(v,O)==0];assert len(roots)==5
    triv=matrix(QQ,[D,O,*roots]).transpose();H=triv.transpose()*NS*triv
    projection=triv*H.inverse()*triv.transpose()*NS
    perpendicular=(triv.transpose()*NS).right_kernel().basis_matrix().transpose()
    base=triv.augment(perpendicular);assert base.rank()==19
    conic_source=read(ART/'det1092_rational_bisection_index_v1/orbit-47755.json')
    C=vector(QQ,[2,4]+conic_source['word']);assert dot(C,F)==dot(C,D)==2 and dot(C,C)==-2
    allrows=[];heights={};selected=[]
    sources=[(i,row) for i,row in enumerate(read(SOURCES)['sources']) if row['degree']==1];assert len(sources)==130
    measured=[]
    for i,row in sources:
        v=vector(QQ,row['word']);S=vector(QQ,[1,v*G*v/2,*v]);phi=S-projection*S
        measured.append(-dot(phi,phi))
    minimum_height=min(h for h in measured if h>0)
    for i,row in sources:
        v=vector(QQ,row['word']);S=vector(QQ,[1,v*G*v/2,*v]);assert dot(S,D)==1
        phi=S-projection*S;height=-dot(phi,phi);heights[str(height)]=heights.get(str(height),0)+1
        assert S==O or height>=QQ(3)/2
        if height!=minimum_height:continue
        bits=[dot(S,r) for r in roots];assert set(bits)<={0,1}
        images=[D,S,*[r if not bit else D-r for r,bit in zip(roots,bits)],
                *[v-dot(v,phi)*D for v in perpendicular.columns()]]
        action=matrix(QQ,images).transpose()*base.inverse()
        assert all(x in ZZ for x in action.list()) and action.transpose()*NS*action==NS
        assert action*O==S and action*D==D
        for sign in [-1,1]:
            signed=action**sign;image=signed*C;degree=dot(image,F)
            assert degree>=0 and degree in ZZ and dot(image,D)==2 and dot(image,image)==-2
            record=dict(source_index=i,sign=sign,word=row['word'],height=str(height),component_bits=ints(bits),
                original_degree=int(degree),image_NS=ints(image),
                action=[ints(v) for v in signed.rows()])
            allrows.append(record)
            if image!=C and degree>1:selected.append(record)
    assert allrows and selected
    minimum=min(row['original_degree'] for row in selected)
    tied=[row for row in selected if row['original_degree']==minimum]
    # Retain every translating section; group duplicate image curves, if any.
    classes=[]
    for row in tied:
        if row['image_NS'] not in [r['image_NS'] for r in classes]:
            classes.append(dict(image_NS=row['image_NS'],original_degree=minimum,
                translators=[dict(source_index=r['source_index'],sign=r['sign']) for r in tied if r['image_NS']==row['image_NS']]))
    save('geometry-v2.json',dict(status='PASS_GENERIC_MINIMUM_IMAGE_SET' if len(classes)<=p['limits']['maps'] else 'UNKNOWN_MAP_CAP',
        height_histogram=heights,short_sections=allrows,minimum_new_original_degree=minimum,
        selected_classes=classes,protocol_sha256=sha(OUT/'protocol-v2.json'),
        boundary='Selection only among the complete old degree1 source set. The minimum is within these130 common sections, not asserted to be the global minimum of the full alternate Mordell-Weil group. No target incidence or new rank.'))
    print('heights',heights,'short sections',len(allrows),'minimum',minimum,'classes',len(classes),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['freeze','lattice']);args=parser.parse_args()
    signal.alarm(25);start=time.monotonic()
    freeze() if args.action=='freeze' else lattice()
    print('seconds',round(time.monotonic()-start,3),flush=True)
