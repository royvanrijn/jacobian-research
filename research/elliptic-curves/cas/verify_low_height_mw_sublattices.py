#!/usr/bin/env python3
"""Read-only exact-coordinate and independent 100-digit replay of both passes."""
from __future__ import annotations
import argparse
from collections import Counter
from decimal import Decimal, localcontext
import importlib
from itertools import combinations, product
import json
from pathlib import Path
import search_low_height_mw_sublattices as base
from latent_lattice import (EllipticCurve, height_gram, exact_graph_walk_consensus,
                            exact_rational_ranks, exact_span_mask)
from mod2_reduction_independence import (find_mod2_reduction_certificate,
                                        combined_mod2_rank,
                                        find_two_torsion_certificate_prime)

PREFIXES=('low_height_mw_sublattices_v1','common_cover_mw_sublattices_v1')


def verify_curve(curve):
    cloudpath=base.OUT/f'{PREFIXES[0]}_{curve}_cloud.json.gz'
    cloud=base.read(cloudpath)
    module=importlib.import_module(f'icarm_curve{curve}')
    short=module.short_coefficients()
    codes=find_mod2_reduction_certificate(short,module.SHORT_POINTS,prime_bound=311)
    rank=combined_mod2_rank(codes,len(module.POINTS))
    assert rank==len(module.POINTS)
    torsion_prime=find_two_torsion_certificate_prime(short,prime_bound=100)
    gram=height_gram(EllipticCurve(module.GENERAL_WEIERSTRASS_COEFFICIENTS),module.POINTS,digits=100,timeout=120)
    with localcontext() as ctx:
        ctx.prec=110
        err=max(abs(Decimal(x)-Decimal(y)) for a,b in zip(gram,cloud['height_gram']) for x,y in zip(a,b))
        assert err<Decimal('1e-65')
    # Re-enumerate at a different precision and independently rebuild the sparse box.
    program=f"""
default(realprecision,100);default(parisizemax,1000000000);default(parisize,100000000);
H={base.matrix_literal(gram)};U={base.gp_matrix(cloud['lll_transform_rows'])};if(abs(matdet(U))!=1,error("ambient basis index"));
T=qflllgram(H);Q=qfminim(T~*H*T,{base.PROTOCOL['height_bounds'][str(curve)]},10000,2);V=T*Q[3];
if(Q[1]!=2*matsize(V)[2] || matsize(V)[2]>=10000,error("cloud truncation"));
for(i=1,matsize(V)[2],print(Vec(V[,i])));
"""
    lines=base.run_gp(program,timeout=120)
    rebuilt_ball={base.canonical_unoriented(json.loads(line)) for line in lines}
    assert rebuilt_ball==set(map(tuple,cloud['ball']))
    u=cloud['lll_transform_rows'];n=len(u);rebuilt_pool=set(rebuilt_ball)
    for k in range(1,4):
        for indices in combinations(range(n),k):
            for signs in product((-1,1),repeat=k):
                v=tuple(sum(u[j][i]*s for i,s in zip(indices,signs)) for j in range(n))
                rebuilt_pool.add(base.canonical_unoriented(v))
    assert rebuilt_pool==set(map(tuple,cloud['pool']))
    stats=[]; checked=set()
    for prefix in PREFIXES:
        selectionpath=base.OUT/f'{prefix}_{curve}_selection.json'
        selection=base.read(selectionpath)
        for path,digest in selection['inputs'].items():
            assert base.digest(base.ROOT/path)==digest, path
        for candidate in selection['finalists']:
            b=candidate['reduced_basis_rows']; a=candidate['basis_rows']; s=candidate['primitive_basis_rows'];r=len(b)
            v=candidate['minimum_vector_public_coordinates']
            program=f"""
default(realprecision,100);H={base.matrix_literal(gram)};B={base.gp_matrix(b)};A={base.gp_matrix(a)};S={base.gp_matrix(s)};v={base.gp_vector(v)};
if(matrank(B)!={r} || vecprod(matsnf(B))!=1,error("primitive rank"));
if(mathnf(B~)!=mathnf(S~),error("reduction changed lattice"));
if(matrank(matconcat([B;A]))!={r},error("generated span changed"));
if(abs(vecprod(matsnf(A)))!={candidate['generated_index_in_primitive_closure']},error("index"));
G=B*H*B~;if(abs(matdet(G)/({candidate['determinant']})-1)>1e-60,error("determinant"));
Q=qfminim(G,,,2);if(abs(Q[2]-({candidate['minimum']}))>1e-60,error("minimum"));
if(abs(v*H*v~-Q[2])>1e-60 || matrank(matconcat([B;Mat(v)]))!={r},error("minimum witness"));
print("PASS");
"""
            assert base.run_gp(program,timeout=120)==['PASS']
            key=tuple(map(tuple,s))
            if key not in checked:
                mask=exact_span_mask(cloud['pool'],s)
                ballmask=exact_span_mask(cloud['ball'],s)
                assert sum(mask)==candidate['pool_support']
                assert sum(ballmask)==candidate['ball_support']
                counts=Counter(tuple(x%2 for x in row) for row,inside in zip(cloud['pool'],mask) if inside)
                assert sum(n*(n-1)//2 for n in counts.values())==candidate['parity']['collision_pairs']
                checked.add(key)
            if prefix.startswith('common'):
                assert all(all((x-y)%2==0 for x,y in zip(row,a[0])) for row in a)
                assert candidate['generated_index_in_primitive_closure'] % 2**(r-1)==0
            stats.append((prefix,candidate['candidate_index']))
        print(f"VERIFIED|curve={curve}|pass={prefix}|finalists={len(selection['finalists'])}",flush=True)
    # Materialize every distinct recorded minimum independently in the elliptic group.
    vectors=sorted({tuple(p['minimum_vector_public_coordinates']) for prefix in PREFIXES
                    for p in base.read(base.OUT/f'{prefix}_{curve}_selection.json')['finalists']})
    points='['+','.join(base.gp_point(p) for p in module.POINTS)+']'
    program=f"""
default(realprecision,100);E=ellinit({base.gp_vector(module.GENERAL_WEIERSTRASS_COEFFICIENTS)});P={points};H={base.matrix_literal(gram)};V={base.gp_matrix(vectors)};
for(j=1,matsize(V)[1],T=[0];for(i=1,#P,if(V[j,i],T=elladd(E,T,ellmul(E,P[i],V[j,i]))));if(!ellisoncurve(E,T) || T==[0],error("point"));if(abs(ellheight(E,T)-V[j,]*H*V[j,]~)>1e-60,error("point height"));print(j,"|",T[1],"|",T[2]));
"""
    witness=base.run_gp(program,timeout=120)
    assert len(witness)==len(vectors)
    return {'curve':curve,'finite_mod2_rank':rank,'two_torsion_exclusion_prime':torsion_prime,
            'height_ball_reenumerated_at_100_digits':len(rebuilt_ball),'combination_pool_reenumerated':len(rebuilt_pool),
            'height_80_vs_100_max_error':str(err),'finalists_verified':len(stats),
            'distinct_primitive_spaces_verified':len(checked),'minimum_points':[{'vector':list(v),'point':line.split('|')[1:]} for v,line in zip(vectors,witness)],
            'finite_mod2_certificates':[code.to_record() if hasattr(code,'to_record') else str(code) for code in codes]}


def legacy_regression():
    """Replay the old successful selector; explicitly separate from the new searches."""
    path=base.OUT/'latent_lattice_fermigier_replay_v1.json.gz'
    doc=base.read(path)
    matrices=[c['proposal']['primitive_basis_rows'] for c in doc['candidates']]
    shapes=[c['primitive_hermite_signature']['log_hermite_invariant'] for c in doc['candidates']]
    ledger=exact_graph_walk_consensus(matrices,shapes,pool_size=64,shape_gap_threshold=.005,graph_weight=1.5,timeout=120)
    index=ledger.selected.source_index
    # Selection has completed before truth is accessed.
    truthpath=base.OUT/'latent_lattice_calibration_truth_v1.json'
    truth=next(t for t in base.read(truthpath)['fermigier_family_controls'] if t['label'].startswith('ICARM_245_'))
    rank=exact_rational_ranks([matrices[index]+truth['embedding_matrix_columns']])[0]
    assert rank==12
    return {'status':'PASS_EXISTING_PRIMITIVE_CLOSURE_REGRESSION','selected_index':index,
            'exact_sum_rank':rank,'generic_subgroup_index':2048,
            'scope':'Replay of the old frozen 128-proposal control ledger; does not validate either new proposal generator or recover the actual index-2048 subgroup.',
            'inputs':{str(p.relative_to(base.ROOT)):base.digest(p) for p in (path,truthpath)}}


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--check',action='store_true');args=parser.parse_args()
    report={'status':'PASS_ARITHMETIC_REPLAY_CALIBRATION_FAILED','curves':[verify_curve(c) for c in (245,302)],
            'legacy_positive_control':legacy_regression()}
    for prefix in PREFIXES:
        calibration=base.read(base.OUT/f'{prefix}_calibration.json')
        assert calibration['status']=='FAIL_BLIND_CALIBRATION'
    report['inputs']={str(p.relative_to(base.ROOT)):base.digest(p) for p in [Path(__file__),Path(base.__file__),
        *[base.OUT/f'{prefix}_{suffix}.json' for prefix in PREFIXES for suffix in ('protocol','calibration','245_selection','302_selection')]]}
    path=base.OUT/'low_height_mw_sublattices_replay_v1.json'
    if args.check:
        assert base.read(path)==report
    else:
        base.save(path,report)
    print(report['status'])


if __name__=='__main__':
    main()
