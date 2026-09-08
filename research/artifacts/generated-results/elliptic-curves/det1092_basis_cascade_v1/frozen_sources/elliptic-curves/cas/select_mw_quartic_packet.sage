#!/usr/bin/env sage-python
"""Development control: primitive MW subspace -> parity packets -> quartic.

Fixed 245 recovered rank12 subspace; numerical height84, <=200000 rays;
first13 rays per class; modular primes1019,1031,1033; <=4096 classes;
one worker, 300 seconds. No control class or section list enters selection.
"""
import argparse
from fractions import Fraction
from hashlib import sha256
import gzip
import importlib.machinery
import importlib.util
import itertools
import json
from pathlib import Path
import signal
import sys
import numpy as np
from sage.all import EllipticCurve, GF, PolynomialRing, QQ, RDF, RealField, ZZ, matrix, pari, prod, vector

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(Path(__file__).resolve().parent))
from icarm_curve245 import POINTS, GENERAL_WEIERSTRASS_COEFFICIENTS
from mestre_root_tuples import SixRootMestreConstruction
from nagao_1994 import primitive_visible_points, quartic_point_to_short_jacobian

GRAPH=ROOT/'artifacts/generated-results/elliptic-curves/latent_lattice_graph_walk_calibration_v1.json'
CLOUD=ROOT/'artifacts/generated-results/elliptic-curves/low_height_mw_sublattices_v1_245_cloud.json.gz'
PAIR=Path(__file__).with_name('recover_mw_quartic_pairs.sage')
OUT=ROOT/'artifacts/generated-results/elliptic-curves/curve302_mw_packet_245_control_v1.json'
LOCAL=ROOT/'artifacts/local/elliptic-curves/mw-packet-control'
POLICY={'height_bound':84,'ray_cap':200000,'packet_size':13,'class_cap':4096,
        'primes':[1019,1031,1033],'seconds':300,'workers':1}


def digest(p): return sha256(p.read_bytes()).hexdigest()
def mask(v): return sum((int(a)%2)<<i for i,a in enumerate(v))
def pointkey(P): return None if P.is_zero() else (int(P[0]),int(P[1]))


def load_detector():
    loader=importlib.machinery.SourceFileLoader('mw_pairs',str(PAIR))
    spec=importlib.util.spec_from_loader(loader.name,loader)
    module=importlib.util.module_from_spec(spec);loader.exec_module(module)
    return module.detect_pairs


def packets():
    graph=json.loads(GRAPH.read_text())
    entry=next(e for e in graph['default_controls'] if e['family']=='Fermigier_rank12')
    B=matrix(ZZ,entry['selected_primitive_embedding_matrix_rows'])
    assert B.nrows()==12 and B.rank()==12
    cloud=json.loads(gzip.decompress(CLOUD.read_bytes()))
    H=matrix(RealField(280),cloud['height_gram']); G=B*H*B.transpose()
    U=matrix(ZZ,pari(G).qflllgram()); R=U.transpose()*G*U
    q=pari(R).qfminim(POLICY['height_bound'],POLICY['ray_cap'],2)
    V=U*matrix(ZZ,q[2])
    assert int(q[0])==2*V.ncols() and V.ncols()<POLICY['ray_cap']
    all_rows=[]
    for v in V.columns():
        if next(a for a in v if a)<0:v=-v
        all_rows.append(tuple(map(int,v)))
    all_rows.sort()
    # Ordering is numerical, not a certified global canonical-height ball.
    real=np.asarray(G,dtype=float); a=np.asarray(all_rows,dtype=np.int64)
    heights=np.einsum('ij,jk,ik->i',a,real,a)
    byclass={}
    for pos in sorted(range(len(a)),key=lambda i:(heights[i],all_rows[i])):
        key=mask(all_rows[pos]); bucket=byclass.setdefault(key,[])
        if len(bucket)<13:bucket.append(all_rows[pos])
    assert len(byclass)<=POLICY['class_cap']
    selected={k:v for k,v in sorted(byclass.items()) if len(v)==13}
    snapshot={'basis_rows':[list(map(int,v)) for v in B.rows()],
        'ray_count':len(all_rows),'ray_sha256':sha256(json.dumps(all_rows).encode()).hexdigest(),
        'occupied_parity_classes':len(byclass), 'eligible_packets':len(selected),
        'packets':{str(k):v for k,v in selected.items()}}
    LOCAL.mkdir(parents=True,exist_ok=True)
    (LOCAL/'packets.json').write_text(json.dumps(snapshot,sort_keys=True)+'\n')
    print('PACKETS',len(all_rows),len(selected),flush=True)
    return B,selected,snapshot


def modular_chart(E, public, B, p):
    F=GF(p); ep=EllipticCurve(F,list(map(F,E.a_invariants())))
    if not ep.discriminant():return None
    group=ep.abelian_group(); generators=[g.element() for g in group.gens()]
    orders=[int(g.order()) for g in generators]
    coordinates={}; points={}
    for cs in itertools.product(*(range(n) for n in orders)):
        Q=sum((a*g for a,g in zip(cs,generators)),ep(0))
        key=pointkey(Q);coordinates[key]=cs;points[cs]=key
    assert len(coordinates)==int(ep.cardinality())
    pubcoords=[]
    for Q in public:
        if Q.is_zero() or QQ(Q[0]).denominator()%p==0:key=None
        else:key=pointkey(ep([F(Q[0]),F(Q[1])]))
        pubcoords.append(coordinates[key])
    basiscoords=np.asarray(B,dtype=np.int64)@np.asarray(pubcoords,dtype=np.int64)
    basiscoords%=np.asarray(orders,dtype=np.int64)
    def xy(v):
        c=np.asarray(v,dtype=np.int64)@basiscoords%np.asarray(orders,dtype=np.int64)
        return points[tuple(map(int,c))]
    def zvalues(rows):
        v0=np.asarray(rows[0],dtype=np.int64); Q0=xy(v0)
        if Q0 is None:return None
        cx,cy=Q0[0],(-Q0[1])%p
        if not cy:return None
        zs=[]
        for v in rows:
            delta=np.asarray(v,dtype=np.int64)-v0
            assert np.all(delta%2==0)
            pt=xy(delta//2)
            if pt is None or pt==(cx,cy):zs.append(None)
            elif pt[0]==cx:
                assert pt[1]==(-cy)%p
                zs.append((3*cx*cx+int(ep.a4()))*pow(-2*cy,-1,p)%p)
            else:zs.append((pt[1]+cy)*pow(pt[0]-cx,-1,p)%p)
        return zs
    return zvalues


def modular_pairs(z,p):
    """Necessary six-disjoint-pair test, or unresolved on collisions.

    Distinct reductions ensure a common-gap solution cannot have its pole
    at any of the twelve used points. Hence all p+1 poles suffice, even
    when the rational Mobius matrix has bad reduction before cancellation.
    """
    if z is None or len(set(z))<13:return None
    offset=0
    while offset in z:offset+=1
    w=np.asarray([0 if a is None else pow((a-offset)%p,-1,p) for a in z],dtype=np.int64)
    inverses=np.asarray([0]+[pow(a,-1,p) for a in range(1,p)],dtype=np.int64)
    delta=(w[None,:]-np.arange(p)[:,None])%p
    xs=np.vstack([inverses[delta],w]); invalid=np.vstack([delta==0,np.zeros((1,13),dtype=bool)])
    ii,jj=np.triu_indices(13,1)
    gaps=(xs[:,ii]-xs[:,jj])%p; gaps=np.minimum(gaps,p-gaps)
    gaps[invalid[:,ii]|invalid[:,jj]]=p
    sorted_gaps=np.sort(gaps,axis=1)
    possible=np.any((sorted_gaps[:,:-5]==sorted_gaps[:,5:]) & (sorted_gaps[:,:-5]<p),axis=1)
    for row in np.flatnonzero(possible):
        for gap,count in zip(*np.unique(gaps[row],return_counts=True)):
            if gap>=p or count<6:continue
            edges=[(int(ii[a]),int(jj[a])) for a in np.flatnonzero(gaps[row]==gap)]
            def match(es,n):
                if not n:return True
                if len(es)<n:return False
                for pos,(i,j) in enumerate(es):
                    if match([(k,l) for k,l in es[pos+1:] if not ({i,j}&{k,l})],n-1):return True
                return False
            if match(edges,6):return True
    return False


def build():
    B,selected,snapshot=packets()
    E0=EllipticCurve(QQ,list(map(QQ,GENERAL_WEIERSTRASS_COEFFICIENTS)))
    E=E0.short_weierstrass_model(); iso=E0.isomorphism_to(E)
    public=[iso(E0(list(map(QQ,P)))) for P in POINTS]
    remaining=dict(selected); stages=[]; exclusions={}
    for p in POLICY['primes']:
        chart=modular_chart(E,public,B,p); rejected=[]; unresolved=[]
        for k,rows in remaining.items():
            decision=None if chart is None else modular_pairs(chart(rows),p)
            if decision is False:rejected.append(k)
            elif decision is None:unresolved.append(k)
        for k in rejected:del remaining[k];exclusions[str(k)]=p
        stages.append({'prime':p,'excluded':len(rejected),'unresolved':len(unresolved),'remaining':len(remaining)})
        (LOCAL/f'after_{p}.json').write_text(json.dumps({'stages':stages,'exclusions':exclusions,'remaining':list(remaining)},sort_keys=True)+'\n')
        print('MODULAR',stages[-1],flush=True)
    detector=load_detector(); exact=[]
    basispoints=[sum((n*P for n,P in zip(v,public)),E(0)) for v in B.rows()]
    for k,rows in remaining.items():
        def point(v):return sum((n*P for n,P in zip(v,basispoints)),E(0))
        C=-point(rows[0]); z=[]
        assert not C.is_zero() and C[1]!=0
        for v in rows:
            delta=vector(ZZ,v)-vector(ZZ,rows[0]);assert all(a%2==0 for a in delta)
            R=point(delta/2)
            if R.is_zero() or R==C:z.append(None)
            elif R[0]==C[0]:
                assert R==-C;z.append(-(3*C[0]**2+E.a4())/(2*C[1]))
            else:z.append((R[1]+C[1])/(R[0]-C[0]))
        result=detector(z)
        verified=[]
        for hit in result['hits']:
            cons=SixRootMestreConstruction(tuple(Fraction(a) for a in hit['roots']))
            RT=PolynomialRing(QQ,'T');T=RT.gen();RX=PolynomialRing(RT,'X');X=RX.gen()
            samples=[cons.primitive_quartic_coefficients(Fraction(i)) for i in range(1,8)]
            coeff=[RT.lagrange_polynomial([(QQ(i+1),QQ(samples[i][j])) for i in range(7)]) for j in range(5)]
            quartic=sum(coeff[i]*X**i for i in range(5))
            product=prod((X-r-T)*(X-r+T) for r in hit['roots']);g=X**6
            for j in range(5,-1,-1):g+=(product[6+j]-(g*g)[6+j])/2*X**j
            assert g*g-product==T*T*QQ(cons.quartic_content)*quartic
            t0=Fraction(hit['T']);J=EllipticCurve(QQ,list(map(QQ,cons.primitive_jacobian_coefficients(t0))))
            jpoints=[J(list(map(QQ,quartic_point_to_short_jacobian(cons,t0,P)))) for P in primitive_visible_points(cons,t0)]
            chosen=sorted({i for pair in hit['pairs'] for i in pair});images={i:point(rows[i]) for i in chosen}
            matches=[]
            for transport in J.isomorphisms(E):
                ids=[[i for i in chosen if transport(P)==images[i] or transport(P)==-images[i]] for P in jpoints]
                if all(len(v)==1 for v in ids) and len({v[0] for v in ids})==12:
                    matches.append({'u_r_s_t':list(map(str,transport.tuple())),'indices':[v[0] for v in ids]})
            assert matches
            verified.append({'roots':hit['roots'],'T':hit['T'],'generic_identity_checked':True,'twelve_images_match_up_to_sign':matches})
        exact.append({'class':k,'rows':[list(v) for v in rows], 'detector':result,'verified_families':verified})
        print('EXACT',k,'hits',len(result['hits']),flush=True)
        (LOCAL/'exact.json').write_text(json.dumps(exact,sort_keys=True)+'\n')
    return {'schema':'curve302.mw-packet-control.v1','status':'COMPLETE_DEVELOPMENT_CONTROL',
        'policy':POLICY,'input_sha256':{str(p.relative_to(ROOT)):digest(p) for p in [Path(__file__),GRAPH,CLOUD,PAIR,Path(__file__).with_name('icarm_curve245.py'),Path(__file__).with_name('mestre_root_tuples.py'),Path(__file__).with_name('nagao_1994.py')]},
        'enumeration':{k:v for k,v in snapshot.items() if k!='packets'},
        'packet_sha256':sha256(json.dumps(snapshot['packets'],sort_keys=True).encode()).hexdigest(),
        'stages':stages,'modular_exclusions':exclusions,'exact':exact,
        'boundary':'Development calibration using the previously recovered245 primitive subspace and a numerical height bound. No known section list or parity class enters packet selection. Numerical enumeration is not an interval-certified complete canonical-height ball. Generic Mestre identity and all twelve visible images are checked for each recovered family. No full generic basis,302 parent or302 exclusion is asserted.'}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--check',action='store_true');args=parser.parse_args()
    signal.alarm(POLICY['seconds']);pari.allocatemem(100000000,1000000000)
    result=build()
    if args.check:assert result==json.loads(OUT.read_text())
    else:OUT.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    print(result['status'])
