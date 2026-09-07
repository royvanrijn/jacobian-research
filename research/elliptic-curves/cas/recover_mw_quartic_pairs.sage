#!/usr/bin/env sage-python
"""Exact six-pair recovery on supplied 245 MW images, with pair labels masked.

This calibrates the geometric recognizer, not selection of the input MW points.
One 13-point control, 4290 quadratics, 1857628 parity pairs, 180 seconds.
"""
import argparse
import gzip
from fractions import Fraction
from hashlib import sha256
import itertools
import json
from pathlib import Path
import signal
import sys
from sage.all import EllipticCurve, PolynomialRing, QQ, ZZ, gcd, lcm, matrix, prod

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(Path(__file__).resolve().parent))
from icarm_curve245 import POINTS, GENERAL_WEIERSTRASS_COEFFICIENTS
from mestre_root_tuples import SixRootMestreConstruction
from nagao_1994 import primitive_visible_points, quartic_point_to_short_jacobian

INPUT=ROOT/'artifacts/generated-results/elliptic-curves/curve302_inverse_fermigier_245_mw_control_v1.json'
CLOUD=ROOT/'artifacts/generated-results/elliptic-curves/low_height_mw_sublattices_v1_245_cloud.json.gz'
OUT=ROOT/'artifacts/generated-results/elliptic-curves/curve302_mw_quartic_pairs_245_control_v1.json'


def digest(path): return sha256(path.read_bytes()).hexdigest()


def disjoint_matchings(edges, count, chosen=()):
    if not count:
        yield chosen
        return
    if len(edges)<count: return
    for pos,(i,j) in enumerate(edges):
        rest=[(k,l) for k,l in edges[pos+1:] if not ({i,j}&{k,l})]
        yield from disjoint_matchings(rest,count-1,chosen+((i,j),))


def detect_pairs(z):
    """All six disjoint equal-gap pairs, modulo affine target coordinates.

    Only projective abscissas enter here. A finite Mobius pole is a root of
    (a-b)(c-k)(d-k) +/- (c-d)(a-k)(b-k). Every qualifying configuration
    contains two disjoint pairs, so enumerating these quadratics is complete
    for the supplied points. The affine chart (pole infinity) is also tested.
    """
    assert len(z)==13 and len(set(z))==13
    offset=QQ(0)
    while offset in z: offset+=1
    w=[QQ(0) if v is None else 1/(v-offset) for v in z]
    poles={None}; equation_count=0
    for i,j,k,l in itertools.combinations(range(len(w)),4):
        for ids in [(i,j,k,l),(i,k,j,l),(i,l,j,k)]:
            a,b,c,d=[w[h] for h in ids]
            for sign in [-1,1]:
                ab=a-b; cd=sign*(c-d)
                A=ab-cd; B=-ab*(c+d)+cd*(a+b); C=ab*c*d-cd*a*b
                equation_count+=1
                if not A:
                    assert B or C # Disjoint distinct pairs cannot give identity.
                    if B: poles.add(-C/B)
                    continue
                disc=B*B-4*A*C
                if disc>=0 and disc.is_square():
                    ds=disc.sqrt()
                    poles.update([(-B+ds)/(2*A),(-B-ds)/(2*A)])
    hits=[]
    for pole in sorted(poles,key=lambda k:(k is not None,k or 0)):
        xs=[v if pole is None else (None if v==pole else 1/(v-pole)) for v in w]
        pairs={}
        for i,j in itertools.combinations(range(len(w)),2):
            if xs[i] is None or xs[j] is None: continue
            diff=abs(xs[i]-xs[j])
            if diff: pairs.setdefault(diff,[]).append((i,j))
        for diff,edges in sorted(pairs.items()):
            if len(edges)<6: continue
            for matching in disjoint_matchings(edges,6):
                centers=sorted((xs[i]+xs[j])/2 for i,j in matching)
                assert len(set(centers))==6
                den=lcm([a.denominator() for a in centers]+[diff.denominator()])
                nums=[ZZ((a-centers[0])*den) for a in centers]
                div=gcd(nums)
                roots=tuple(a//div for a in nums)
                roots=min(roots,tuple(roots[-1]-a for a in roots[::-1]))
                hits.append({'pole':None if pole is None else str(pole),
                    'pairs':[list(p) for p in matching],
                    'roots':list(map(int,roots)), 'T':str(diff*den/(2*div))})
    return {'finite_chart_offset':str(offset), 'quadratic_equations':equation_count,
            'candidate_poles_including_infinity':len(poles), 'hits':hits}


def build():
    control=json.loads(INPUT.read_text())
    E=EllipticCurve(QQ,list(map(QQ,GENERAL_WEIERSTRASS_COEFFICIENTS)))
    Es=E.short_weierstrass_model(); iso=E.isomorphism_to(Es)
    public=[iso(E(list(map(QQ,p)))) for p in POINTS]
    columns=list(matrix(ZZ,control['exact_embedding_matrix_rows']).columns())
    # The visible relation supplies the omitted twelfth image.
    columns=columns[:11]+[-sum(columns[:11])]+[columns[11]]
    # Remove pair adjacency and the distinguished extra point from selection.
    columns.sort(key=lambda v:sha256(','.join(map(str,v)).encode()).hexdigest())
    images=[sum((c*p for c,p in zip(v,public)),Es(0)) for v in columns]
    C=-images[0]; lifted=[]
    for v,Q in zip(columns,images):
        delta=v-columns[0]
        assert all(a%2==0 for a in delta)
        R=sum((ZZ(a/2)*p for a,p in zip(delta,public)),Es(0))
        assert 2*R-C==Q
        lifted.append(R)
    assert not C.is_zero() and C[1]!=0
    assert all(R.is_zero() or R==C or R[0]!=C[0] for R in lifted)
    z=[None if R.is_zero() or R==C else (R[1]+C[1])/(R[0]-C[0]) for R in lifted]
    result=detect_pairs(z)
    # Known coordinates/parameter are accessed only for post-selection checks.
    assert len(result['hits'])==1
    hit=result['hits'][0]
    assert hit['roots']==list(map(int,control['roots']))
    assert QQ(hit['T'])==QQ(control['T_canonical'])
    cons=SixRootMestreConstruction(tuple(Fraction(a) for a in hit['roots']))
    RT=PolynomialRing(QQ,'T'); T=RT.gen(); RX=PolynomialRing(RT,'X'); X=RX.gen()
    samples=[cons.primitive_quartic_coefficients(Fraction(i)) for i in range(1,8)]
    coeff=[RT.lagrange_polynomial([(QQ(i+1),QQ(samples[i][j])) for i in range(7)]) for j in range(5)]
    quartic=sum(coeff[i]*X**i for i in range(5))
    product=prod((X-r-T)*(X-r+T) for r in hit['roots'])
    g=X**6
    for k in range(5,-1,-1): g+=(product[6+k]-(g*g)[6+k])/2*X**k
    assert g*g-product==T*T*QQ(cons.quartic_content)*quartic
    t0=Fraction(hit['T'])
    J=EllipticCurve(QQ,list(map(QQ,cons.primitive_jacobian_coefficients(t0))))
    qpoints=list(primitive_visible_points(cons,t0))
    jpoints=[J(list(map(QQ,quartic_point_to_short_jacobian(cons,t0,p)))) for p in qpoints]
    chosen=sorted({i for pair in hit['pairs'] for i in pair})
    matching=[]
    for transport in J.isomorphisms(Es):
        mapped=[transport(p) for p in jpoints]
        matches=[[i for i in chosen if q==images[i] or q==-images[i]] for q in mapped]
        if all(len(ids)==1 for ids in matches) and len({ids[0] for ids in matches})==12:
            matching.append({'to_short_u_r_s_t':list(map(str,transport.tuple())),
                             'image_indices':[ids[0] for ids in matches]})
    assert matching
    # Diagnose the old selector separately, after the geometric recognition.
    cloud=json.loads(gzip.decompress(CLOUD.read_bytes()))
    def ray(v):
        v=tuple(map(int,v))
        return v if next(a for a in v if a)>0 else tuple(-a for a in v)
    wanted=set(map(ray,columns))
    coverage={key:len(wanted.intersection(map(ray,cloud[key]))) for key in ['ball','pool']}
    parity_seeds=sum(all((a-b)%2==0 for a,b in zip(v,columns[0])) for v in cloud['ball'])
    assert coverage=={'ball':0,'pool':0} and parity_seeds==0
    def parity(v): return sum((int(a)%2)<<i for i,a in enumerate(v))
    masks=[parity(v) for v in cloud['ball']]; lookup={v:i for i,v in enumerate(masks)}
    assert len(lookup)==len(masks)
    target=parity(columns[0]); pair_count=0; triple_count=0; witnesses=[]
    for i,x in enumerate(masks):
        if lookup.get(target^x,-1)>i: pair_count+=1
        for j in range(i+1,len(masks)):
            k=lookup.get(target^x^masks[j],-1)
            if k>j:
                triple_count+=1
                if len(witnesses)<3: witnesses.append([i,j,k])
    assert pair_count==0 and triple_count==1459
    result.update({'schema':'curve302.mw-quartic-pairs-control.v1',
        'status':'PASS_SUPPLIED_MW_IMAGES_TO_SIX_PAIR_FAMILY',
        'input_sha256':{str(p.relative_to(ROOT)):digest(p) for p in [Path(__file__),INPUT,CLOUD,
            Path(__file__).with_name('icarm_curve245.py'),Path(__file__).with_name('mestre_root_tuples.py'),
            Path(__file__).with_name('nagao_1994.py')]},
        'masked_embedding_columns':[list(map(int,v)) for v in columns],
        'projective_abscissas':[None if v is None else str(v) for v in z],
        'generic_mestre_identity_checked':True,
        'specialized_twelve_covariant_images_match_up_to_sign':matching,
        'old_selector_coverage':{'seed_vectors':len(cloud['ball']),
            'sparse_pool_vectors':len(cloud['pool']), 'known_image_rays_found':coverage,
            'seeds_in_known_image_parity_class':parity_seeds,
            'two_distinct_ball_ray_parity_decompositions':pair_count,
            'three_distinct_ball_ray_parity_decompositions':triple_count,
            'first_three_triple_witness_indices_zero_based':witnesses,
            'consequence':'No vector c+2z for any old ball seed c and any integral displayed-basis vector z can equal a known control image, even up to sign. Enlarging only the shift bound cannot fix this selector.'},
        'limits':{'input_points':13,'quadratic_equations':4290,
                  'parity_pair_comparisons':len(masks)*(len(masks)-1)//2,'seconds':180},
        'boundary':'Input is the known 245 section-image set. Pair labels and quartic coordinates are withheld from detection. This certifies the geometric recognizer only, not blind MW-point selection, generic saturation, or a parent or exclusion for302.'})
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check',action='store_true'); args=parser.parse_args()
    signal.alarm(180)
    result=build()
    if args.check: assert result==json.loads(OUT.read_text())
    else: OUT.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    print(result['status'])
    print(json.dumps({k:result[k] for k in ['quadratic_equations','candidate_poles_including_infinity','hits']}))
