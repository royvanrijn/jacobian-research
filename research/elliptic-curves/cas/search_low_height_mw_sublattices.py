#!/usr/bin/env python3
"""Bounded, truth-free combination search; numerical scores, exact embeddings.

Run selection before the separate --evaluate command opens calibration truth.
No point, rank, or generic-family discovery is implied by a selected space.
"""
from __future__ import annotations

import argparse
from collections import Counter
import gzip
from hashlib import sha256
import importlib
from itertools import combinations, product
import json
import os
from pathlib import Path
import sys

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / "elliptic-curves"), str(Path(__file__).parent)]
import numpy as np
from latent_lattice.integer import canonical_unoriented
from latent_lattice.pari import gp_matrix, gp_point, gp_vector, run_gp
from latent_lattice import exact_rational_ranks

OUT = ROOT / "artifacts/generated-results/elliptic-curves"
PREFIX = "low_height_mw_sublattices_v1"
PROTOCOL = {
    "ranks": list(range(8, 21)), "seed": 20260906,
    "height_bounds": {"245": 28, "302": 70}, "digits": 80,
    "short_ball_cap": 10000, "combination_support": 3,
    "combination_coefficients": [-1, 1], "paths_per_channel": 16,
    "channels": ["determinant", "relation", "random"],
    "finalists_per_rank_per_score": 2, "gp_timeout_seconds": 120,
    "proof_boundary": "Heuristic subspace proposals, exact integer embeddings; numerical height balls and determinants, no certified interval bounds. No exhaustive sublattice search.",
}


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def save(path, obj):
    data = (json.dumps(obj, sort_keys=True, indent=2) + "\n").encode()
    if path.suffix == ".gz":
        data = gzip.compress(data, mtime=0)
    if path.exists():
        if path.read_bytes() != data:
            raise ValueError(f"refusing to replace different checkpoint {path}")
        return
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_bytes(data)
    temp.replace(path)


def read(path):
    return json.loads(gzip.decompress(path.read_bytes()) if path.suffix == ".gz" else path.read_text())


def matrix_literal(h):
    return "[" + ";".join(",".join(map(str, row)) for row in h) + "]"


def shell(n, support=3):
    """Exactly one orientation, no division by content or accidental doubling."""
    for k in range(1, min(n, support) + 1):
        for indices in combinations(range(n), k):
            for tail in product((-1, 1), repeat=k-1):
                v = [0] * n
                for i, a in zip(indices, (1,) + tail):
                    v[i] = a
                yield v


def cloud(curve_id):
    path = OUT / f"{PREFIX}_{curve_id}_cloud.json.gz"
    if path.exists():
        return read(path)
    source = Path(__file__).parent / f"icarm_curve{curve_id}.py"
    data = importlib.import_module(f"icarm_curve{curve_id}")
    p = "[" + ",".join(gp_point(x) for x in data.POINTS) + "]"
    program = f"""
default(realprecision,80);default(parisizemax,1000000000);default(parisize,100000000);
E=ellinit({gp_vector(data.GENERAL_WEIERSTRASS_COEFFICIENTS)});P={p};n=#P;
if(!prod(i=1,n,ellisoncurve(E,P[i])),error("off curve"));
H=ellheightmatrix(E,P);U=qflllgram(H);R=U~*H*U;
Q=qfminim(R,{PROTOCOL['height_bounds'][str(curve_id)]},10000,2);V=U*Q[3];
print("VERSION|",version());print("COUNT|",Q[1],"|",matsize(V)[2]);
for(i=1,n,print("H|",strjoin(vector(n,j,Str(H[i,j])),"|")));
print("U|",vector(n,i,Vec(U[i,])));
for(j=1,matsize(V)[2],print("V|",Vec(V[,j])));
"""
    lines = run_gp(program, timeout=120)
    h = [x.split("|")[1:] for x in lines if x.startswith("H|")]
    u = json.loads(next(x[2:] for x in lines if x.startswith("U|")))
    count = next(x.split("|")[1:] for x in lines if x.startswith("COUNT|"))
    if int(count[1]) >= 10000 or int(count[0]) != 2 * int(count[1]):
        raise ArithmeticError("incomplete height ball")
    ball = sorted({canonical_unoriented(json.loads(x[2:])) for x in lines if x.startswith("V|")})
    combos = np.asarray(list(shell(len(h))), dtype=np.int64) @ np.asarray(u, dtype=np.int64).T
    pool = sorted(set(ball) | {canonical_unoriented(v) for v in combos})
    doc = {"curve": curve_id, "source": str(source.relative_to(ROOT)), "source_sha256": digest(source),
           "height_gram": h, "lll_transform_rows": u, "ball": ball, "pool": pool,
           "signed_count": int(count[0]), "sparse_shell_lines": len(combos),
           "pari_version": next(x.split("|",1)[1] for x in lines if x.startswith("VERSION|"))}
    save(path, doc)
    return doc


def parity(vectors):
    counts = Counter(tuple(int(a) % 2 for a in v) for v in vectors)
    return {"cosets": len(counts), "largest_coset": max(counts.values(), default=0),
            "collision_pairs": sum(n*(n-1)//2 for n in counts.values()),
            "multiplicity_histogram": dict(sorted(Counter(counts.values()).items()))}


def proposals(doc):
    """48 deterministic growing bases, selected from combinations in the ball."""
    h = np.asarray(doc['height_gram'], dtype=float)
    v = np.asarray(doc['ball'], dtype=np.int64)
    norms = np.einsum('ij,jk,ik->i', v, h, v)
    order = np.argsort(norms, kind='stable'); v = v[order]; norms = norms[order]
    lookup = {tuple(x): i for i, x in enumerate(v)}
    degrees = np.zeros(len(v), dtype=int)
    # Complete unit ternary relations in this bounded cloud.
    edges = set()
    for i in range(len(v)):
        for j in range(i):
            for w in (v[i]+v[j], v[i]-v[j]):
                if not np.any(w):
                    continue
                # Do not primitive-normalize: a+b=2c is not a unit relation.
                if w[np.flatnonzero(w)[0]] < 0:
                    w = -w
                k = lookup.get(tuple(w))
                if k is not None and k not in (i,j):
                    edges.add(tuple(sorted((i,j,k))))
    for edge in edges:
        for i in edge:
            degrees[i] += 1
    y = v @ np.linalg.cholesky(h)
    rng = np.random.default_rng(PROTOCOL['seed'])
    rows = []
    for channel in PROTOCOL['channels']:
        seeds = (np.argsort(-degrees, kind='stable')[:16] if channel == 'relation'
                 else np.arange(16) if channel == 'determinant'
                 else rng.choice(len(v), 16, replace=False))
        for seed in seeds:
            basis, q = [], []
            residual = y.copy()
            for rank in range(1,21):
                rn = np.einsum('ij,ij->i', residual, residual)
                valid = rn > 1e-8
                if rank == 1:
                    chosen = int(seed)
                elif channel == 'determinant':
                    chosen = int(np.argmin(np.where(valid,rn,np.inf)))
                elif channel == 'relation':
                    chosen = int(np.argmax(np.where(valid,(degrees+1)/(1+rn/norms),-1)))
                else:
                    weights = np.where(valid,(degrees+1)**0.5,0.)
                    chosen = int(rng.choice(len(v),p=weights/weights.sum()))
                direction = residual[chosen] / np.sqrt(rn[chosen])
                q.append(direction); basis.append(v[chosen].tolist())
                # Two projection passes keep rank screening stable; GP certifies rank.
                for _ in range(2):
                    residual -= np.outer(residual @ direction, direction)
                if rank >= 8:
                    rows.append({'rank':rank,'channel':channel,'seed_index':int(seed),
                                 'basis_rows':[x[:] for x in basis]})
    return rows, len(edges)


def score_batch(doc, proposals_):
    h = matrix_literal(doc['height_gram']); result=[]
    pool = np.asarray(doc['pool'], dtype=np.int64)
    ball_set = set(map(tuple,doc['ball']))
    for start in range(0,len(proposals_),32):
        batch = proposals_[start:start+32]
        matrices = ','.join(gp_matrix(p['basis_rows']) for p in batch)
        program = f"""
default(realprecision,80);default(parisizemax,1000000000);default(parisize,100000000);
H={h};V=[{matrices}];
for(i=1,#V,B=V[i];r=matsize(B)[1];n=matsize(B)[2];if(matrank(B)!=r,error("rank"));S=if(r==n,matid(n),matkerint(matkerint(B)~)~);S=mathnf(S~)~;K=matkerint(S);G=S*H*S~;T=qflllgram(G);L=T~*S;R=L*H*L~;print(i,"|",prod(j=1,r,abs(matsnf(B)[j])),"|",matdet(G),"|",vector(r,j,R[j,j]),"|",vector(r,j,Vec(L[j,])),"|",if(r==n,[],vector(n,j,Vec(K[j,]))),"|",vector(r,j,Vec(S[j,]))));
"""
        lines=run_gp(program,timeout=120)
        if len(lines)!=len(batch):
            raise ArithmeticError('missing GP candidate')
        for p,line in zip(batch,lines):
            i,index,det,diagonal,basis,kernel,primitive=line.split('|')
            k=np.asarray(json.loads(kernel),dtype=object)
            if k.size:
                km=np.asarray(k % 1009,dtype=np.int64)
                mask=np.all((pool % 1009) @ km % 1009 == 0,axis=1)
                ids=np.flatnonzero(mask)
                mask[ids]=np.all(pool[ids].astype(object) @ k == 0,axis=1)
            else:
                mask=np.ones(len(pool),dtype=bool)
            selected=pool[mask]
            in_ball=sum(tuple(v) in ball_set for v in selected)
            result.append({**p,'primitive_basis_rows':json.loads(primitive),
                'reduced_basis_rows':json.loads(basis),'generated_index_in_primitive_closure':int(index),
                'determinant':det,'lll_heights':json.loads(diagonal),
                'ball_support':in_ball,'pool_support':len(selected),'parity':parity(selected)})
        print(f"SCORED|curve={doc['curve']}|{min(start+32,len(proposals_))}/{len(proposals_)}",flush=True)
    # Deduplicate exact HNF, preserving every proposal source.
    unique={}
    for p in result:
        key=tuple(map(tuple,p['primitive_basis_rows']))
        if key not in unique:
            unique[key]=p
    return list(unique.values())


def finalists(doc, candidates):
    chosen={}
    for rank in PROTOCOL['ranks']:
        group=[(i,p) for i,p in enumerate(candidates) if p['rank']==rank]
        for score,key in [('determinant',lambda t:float(t[1]['determinant'])),
                          ('density',lambda t:-t[1]['ball_support']),
                          ('parity',lambda t:-t[1]['parity']['collision_pairs'])]:
            for i,p in sorted(group,key=key)[:2]:
                chosen.setdefault(i,[]).append(score)
    output=[]
    for i,labels in chosen.items():
        p=candidates[i]; b=gp_matrix(p['reduced_basis_rows']); h=matrix_literal(doc['height_gram'])
        program=f"""
default(realprecision,80);H={h};B={b};r=matsize(B)[1];n=matsize(B)[2];G=B*H*B~;
Q=qfminim(G,vecmin(vector(r,j,G[j,j])),10000,2);
if(Q[1]!=2*matsize(Q[3])[2] || matsize(Q[3])[2]>=10000,error("SVP cap"));
N=vector(matsize(Q[3])[2],j,Q[3][,j]~*G*Q[3][,j]);m=vecmin(N);j=1;while(N[j]!=m,j++);print("MIN|",m);
print("SV|",Vec(Q[3][,j])*B);
if(r<n,{{W=matsnf(B,1)[2]^-1;C=W[1..n-r,];if(abs(matdet(matconcat([B;C])))!=1,error("completion"));K=C*H*C~-C*H*B~*G^-1*B*H*C~;T=qflllgram(K);K=T~*K*T;print("COMP|",vector(n-r,j,Vec(K[j,])));print("CD|",matdet(K));}});
"""
        lines=run_gp(program,timeout=120)
        extra={'minimum':next(x[4:] for x in lines if x.startswith('MIN|')),
               'minimum_vector_public_coordinates':json.loads(next(x[3:] for x in lines if x.startswith('SV|')))}
        comp=next((x[5:] for x in lines if x.startswith('COMP|')),None)
        if comp:
            g=np.asarray(json.loads(comp),dtype=float); scale=float(np.mean(np.diag(g))/2)
            nearest=np.rint(g/scale).astype(int)
            extra['projected_quotient']={'gram':g.tolist(),'determinant':next(x[3:] for x in lines if x.startswith('CD|')),
                'root_scale':scale,'nearest_integral_gram':nearest.tolist(),
                'relative_rounding_rms':float(np.sqrt(np.mean((g/scale-nearest)**2))),
                'is_even_integral_root_lattice':'UNKNOWN',
                'warning':'Projected quotient lattice, not an exact orthogonal integer complement. Rounded Gram is only a diagnostic.'}
        output.append({'candidate_index':i,'selected_by':labels,**p,**extra})
    return output


def run(curve_id):
    doc=cloud(curve_id)
    pp=OUT/f'{PREFIX}_{curve_id}_proposals.json.gz'
    if pp.exists():
        proposed=read(pp)
    else:
        rows,edges=proposals(doc)
        proposed={'candidates':rows,'unit_ternary_relations':edges}
        save(pp,proposed)
    cp=OUT/f'{PREFIX}_{curve_id}_candidates.json.gz'
    if cp.exists():
        candidates=read(cp)
    else:
        candidates=score_batch(doc,proposed['candidates']);save(cp,candidates)
    final=finalists(doc,candidates)
    result={'curve':curve_id,'status':'BOUNDED_CANDIDATES_CALIBRATION_PENDING',
        'protocol':PROTOCOL,'cloud_lines':len(doc['ball']),'pool_lines':len(doc['pool']),
        'ball_parity':parity(doc['ball']),'pool_parity':parity(doc['pool']),
        'candidate_count':len(candidates),'finalists':final,
        'inputs':{str(p.relative_to(ROOT)):digest(p) for p in (pp,cp,OUT/f'{PREFIX}_{curve_id}_cloud.json.gz',Path(__file__))}}
    save(OUT/f'{PREFIX}_{curve_id}_selection.json',result)
    print(f"SELECTED|{curve_id}|candidates={len(candidates)}|finalists={len(final)}",flush=True)


def evaluate():
    # This is the first and only stage allowed to load the hidden generic subgroup.
    path=OUT/f'{PREFIX}_245_selection.json'; selection=read(path)
    truth_path=OUT/'latent_lattice_calibration_truth_v1.json'
    truth=next(x for x in read(truth_path)['fermigier_family_controls'] if x['label'].startswith('ICARM_245_'))
    t=truth['embedding_matrix_columns']; tr=truth['truth_rank']
    candidates=read(OUT/f'{PREFIX}_245_candidates.json.gz')
    ranks=exact_rational_ranks([p['primitive_basis_rows']+t for p in candidates])
    intersections=[p['rank']+tr-r for p,r in zip(candidates,ranks)]
    exact=[i for i,p in enumerate(candidates) if p['rank']==tr and intersections[i]==tr]
    selected=[p['candidate_index'] for p in selection['finalists']]
    output={'status':'PASS_PRIMITIVE_SPACE_RECOVERY' if set(exact)&set(selected) else 'FAIL_BLIND_CALIBRATION',
        'truth_rank':tr,'truth_generic_index':truth['generic_subgroup_index_in_primitive_closure'],
        'exact_truth_candidate_indices':exact,'selected_exact_truth_indices':sorted(set(exact)&set(selected)),
        'best_intersection_by_rank':{str(r):max(intersections[i] for i,p in enumerate(candidates) if p['rank']==r) for r in PROTOCOL['ranks']},
        'actual_generic_subgroup_recovered':False,
        'boundary':'An exact primitive-space recovery would still not recover the index-2048 generic subgroup. No target generic-family or root-complement certification follows.',
        'inputs':{str(p.relative_to(ROOT)):digest(p) for p in (path,truth_path)}}
    save(OUT/f'{PREFIX}_calibration.json',output);print(json.dumps(output,indent=2))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--curve',type=int,choices=(245,302))
    parser.add_argument('--evaluate',action='store_true')
    args=parser.parse_args()
    save(OUT/f'{PREFIX}_protocol.json',{'protocol':PROTOCOL,'selection_script_sha256':digest(Path(__file__))})
    if args.evaluate:
        evaluate()
    elif args.curve:
        run(args.curve)
    else:
        parser.error('choose --curve or --evaluate')


if __name__=='__main__':
    main()
