#!/usr/bin/env sage-python
"""Frozen, bounded MW16/MW18 fibre-population comparison.

prepare measures equations and quarantined chart costs, then freezes both arms.
search reads that immutable population; verify replays its arithmetic and hits.
No bounded negative is an exclusion or an equivalence test.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
from fractions import Fraction as F
import gzip
from hashlib import sha256
from importlib.machinery import SourceFileLoader
import json
from math import log2
from pathlib import Path
import platform
import statistics
import sys
import time

from sage.all import EllipticCurve, PolynomialRing, QQ, ZZ, gcd, lcm, prime_range, version
from cysignals.alarm import alarm, cancel_alarm
from cysignals.signals import AlarmInterrupt

ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / 'elliptic-curves/cas'
sys.path.insert(0, str(CAS))
sys.set_int_max_str_digits(0)
import pointed_quartic_search as detector
import mw16_model_size as sizes
from mod_l_reduction_independence import (combined_mod_l_rank,
    find_mod_l_reduction_certificate, find_no_rational_l_torsion_prime,
    mod_l_reduction_signature)

ART = ROOT / 'artifacts/generated-results/elliptic-curves'
LOCAL = ROOT / 'artifacts/local/elliptic-curves/fibre-height-population'
PROTOCOL = ROOT / 'elliptic-curves/data/fibre_height_population_v1.json'
POPULATION = ART / 'fibre_height_population_v1.json.gz'
RESULTS = ART / 'fibre_height_search_v1.json.gz'
SUMMARY = ART / 'fibre_height_comparison_v1.json'
TEMPLATE = ROOT / 'elliptic-curves/data/a1_mw16_family_template_v1.json'
COVERS = ROOT / 'artifacts/generated-results/elkies-k3-r17-extreme-anchored-mw18-covers-v1.json'
R = PolynomialRing(QQ, 'z')
z = R.gen()


def read(path):
    raw = path.read_bytes()
    return json.loads(gzip.decompress(raw) if path.suffix == '.gz' else raw)


def write(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    raw = (json.dumps(data, sort_keys=True, indent=2, allow_nan=False)+'\n').encode()
    if path.suffix == '.gz':
        raw = gzip.compress(raw, mtime=0)
    tmp = path.with_suffix(path.suffix+'.tmp')
    tmp.write_bytes(raw)
    tmp.replace(path)


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def content_hash(data):
    return sha256(json.dumps(data, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def source(name):
    path = ROOT / name
    return SourceFileLoader(path.stem.replace('.', '_'), str(path)).load_module()


def provenance(protocol):
    paths = [Path(__file__), PROTOCOL, TEMPLATE, COVERS, CAS/'mw16_model_size.py',
        CAS/'build_a1_mw16_target_free_parameter_candidates.sage',
        CAS/'mod_l_reduction_independence.py', CAS/'mod2_reduction_independence.py',
        ROOT/'elkies-k3/scripts/specialize_r17_extreme_anchored_mw18_finalists.sage',
        ROOT/'elkies-k3/scripts/search_h92_q12o5867_rootless_nagao.py',
        ROOT/'elkies-k3/scripts/construct_elkies_2026_bisections.sage',
        ROOT/'artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json']
    paths += [ROOT/p for p in protocol['source_hashes']]
    return {**detector.sources(), **{str(p.relative_to(ROOT)):digest(p) for p in paths}}


def root_part(n, k):
    """A verified kth-power divisor, using small valuations and a root test.

    This need not find the maximal divisor. Incompleteness only changes the
    coordinate, never an arithmetic certificate; no large factorization.
    """
    n = abs(ZZ(n))
    if not n:
        return ZZ(1)
    root = ZZ(1)
    for p in prime_range(1000):
        v = n.valuation(p)
        n //= ZZ(p)**v
        root *= ZZ(p)**(v//k)
    residual, exact = n.nth_root(k, truncate_mode=True)
    return root*residual if exact else root


def mw16_coordinate(A):
    d = A.degree()
    shift = -A[d-1]/(d*A[d])
    centered = A(z+shift)
    q2, q3 = centered[d-2]/A[d], centered[d-3]/A[d]
    if not q2 or not q3:
        raise ValueError('coordinate recipe requires nonzero degree d-2,d-3 terms')
    common = QQ(gcd(q2.numerator()**3, q3.numerator()**2))/lcm(q2.denominator()**3,q3.denominator()**2)
    scale = QQ(root_part(common.numerator(),6))/root_part(common.denominator(),6)
    return scale, shift


def mw18_coordinate(cover):
    q = R(cover['branch_quadratic_coefficients_low_to_high'])
    den = lcm([a.denominator() for a in q])
    integral = q*den**2
    square = root_part(gcd(list(integral)), 2)
    u_scale = QQ(square)/den
    par = cover['anchor_line_parameterization']
    slope_scale = QQ(par['line_slope_scale'])
    return u_scale/slope_scale, -QQ(par['tangent_slope_at_anchor'])/slope_scale


def rational_function(record, ring=R):
    return ring(record['numerator_coefficients_low_to_high'])/ring(record['denominator_coefficients_low_to_high'])


def logheight(q):
    q = QQ(q)
    return log2(int(max(abs(q.numerator()), q.denominator())))


def families(protocol, sections=False):
    """Load fixed family identities, never ranks or exceptional point lists."""
    builder = source('elliptic-curves/cas/build_a1_mw16_target_free_parameter_candidates.sage')
    spec = source('elkies-k3/scripts/specialize_r17_extreme_anchored_mw18_finalists.sage')
    template = read(TEMPLATE)
    assert template['status'] == 'PASS_TARGET_FREE_A1_MW16_FAMILY_PRESENTATIONS'
    presentation = next(p for p in template['presentations'] if p['presentation_id']==protocol['mw16_presentation'])
    A = R(presentation['pencil']['A_coefficients_low_to_high'])
    B = R(presentation['pencil']['B_coefficients_low_to_high'])
    first = dict(name='mw16', A=A, B=B, base=z, coordinate=mw16_coordinate(A), rank=16,
        presentation=presentation, spec=builder, family_id=presentation['fibration_id'])
    certificate = read(COVERS)
    assert certificate['status'] == 'PASS_EXACT_EXTREME_ANCHORED_MW18_COVERS'
    cover, model_path, _ = spec.cover_by_label(certificate, protocol['mw18_cover'])
    model = read(model_path)
    eq = model.get('weierstrass_model',model)
    second = dict(name='mw18', A=R(eq['A_coefficients_low_to_high']),
        B=R(eq['B_coefficients_low_to_high']), base=rational_function(cover['anchor_line_parameterization']['t_of_r']),
        coordinate=mw18_coordinate(cover), rank=18, model=model, cover=cover, spec=spec,
        family_id=protocol['mw18_cover'])
    if sections:
        old = read(builder.MODEL)
        ring = PolynomialRing(QQ,'t')
        oldA, oldB = (ring(old['weierstrass_model'][key]) for key in ('A_coefficients_low_to_high','B_coefficients_low_to_high'))
        E = EllipticCurve(ring.fraction_field(), [oldA,oldB])
        basis = [E(builder.polynomial_from_record(p['X'],ring),builder.polynomial_from_record(p['Y'],ring)) for p in old['sections']['records']]
        marking = presentation['source_marking']
        def combination(v):
            return sum((int(c)*p for c,p in zip(v,basis) if c),E(0))
        chord = source('elkies-k3/scripts/construct_elkies_2026_bisections.sage')
        h,nx,m0,quartic,a,b = builder.child_geometry(combination(marking['trace_section_basis_w']),oldA,oldB,ring,chord)
        if R(a)!=A or R(b)!=B:
            raise ArithmeticError('reconstructed MW16 pencil changed')
        first['section_context'] = dict(source_points=tuple(combination(v) for v in
            [marking['new_zero_source_section_basis_coordinates'],*marking['generic_source_section_basis_coordinates']]),
            base_maps=tuple(builder.polynomial_from_record(p,ring) for p in marking['base_maps_lambda_of_old_t']),
            old_ring=ring,h=h,nx=nx,m0=m0,quartic=quartic,child_a=a,child_b=b)
    return [first,second]


def specialize(family, parameter):
    if family['name']=='mw16':
        row, failure = family['spec'].specialize(parameter=parameter, **family['section_context'])
        if failure:
            raise ValueError(failure)
        return row['raw_generic_points']
    spec, cover = family['spec'], family['cover']
    t = family['base'](parameter)
    u = rational_function(cover['anchor_line_parameterization']['u_of_r'])(parameter)
    assert u*u == R(cover['branch_quadratic_coefficients_low_to_high'])(t)
    points = spec.direct_generic_values(family['model'], t)
    p = cover['eighteenth_section']
    points.append((R(p['x0_coefficients_low_to_high'])(t)+R(p['x1_coefficients_low_to_high'])(t)*u,
        R(p['y0_coefficients_low_to_high'])(t)+R(p['y1_coefficients_low_to_high'])(t)*u))
    return [dict(x=str(x),y=str(y)) for x,y in points]


def exact_fibre(family, parameter):
    t = family['base'](parameter)
    return EllipticCurve(QQ, [family['A'](t),family['B'](t)])


def weighted_size(E):
    # Integral global minimal model; translations fixed by PARI normalization.
    return max([ZZ(1)]+[abs(ZZ(a))**(12//w) for a,w in zip(E.ainvs(),(1,2,3,4,6))])


def remove_obvious_scale(E,hints=()):
    """Remove verified 4th/6th powers before PARI sees huge common factors."""
    A,B=E.a4(),E.a6()
    den=lcm(A.denominator(),B.denominator())
    a,b=ZZ(A*den**4),ZZ(B*den**6)
    u=ZZ(1)
    # Separating denominator/coordinate factors avoids factoring a huge
    # product of a square and an unrelated nonsquare twist coefficient.
    seeds=[abs(a),abs(b),den]+[abs(ZZ(v)) for v in hints if v]
    for base in ZZ.gcd_free_basis(seeds):
        va,vb=a.valuation(base),b.valuation(base)
        residual=base
        for p in prime_range(1000):
            exponent=residual.valuation(p)
            residual//=ZZ(p)**exponent
            u*=ZZ(p)**min(exponent*va//4,exponent*vb//6)
        for power in range(12,0,-1):
            root,exact=residual.nth_root(power,truncate_mode=True)
            if exact:
                u*=root**min(power*va//4,power*vb//6)
                break
    answer=EllipticCurve(QQ,[a/u**4,b/u**6])
    assert E.is_isomorphic(answer)
    return answer


def score_candidate(common, tables, pair, height):
    c = common.Candidate(int(pair[0]),int(pair[1]),int(height))
    counts=[]
    for block in tables:
        good,bad=c.good_primes,c.bad_primes
        c = common.score_block(c,block,{})
        counts.append([c.good_primes-good,c.bad_primes-bad])
    return {**common.candidate_record(c),'block_good_bad_counts':counts}


def candidate_from_row(common,row,stage=None):
    nagao = row['nagao']
    c = common.Candidate(*row['box_pair'],max(map(abs,row['box_pair'])))
    from dataclasses import replace
    return replace(c,block_score_units=tuple(nagao['block_score_units_1e12'][:stage]),
        good_primes=sum(v[0] for v in nagao['block_good_bad_counts'][:stage]),
        bad_primes=sum(v[1] for v in nagao['block_good_bad_counts'][:stage]))


def select_arms(rows, protocol):
    """Same deduplicated union for both arms; overlap is retained and shared.

    Nagao retains the existing staged capped Pareto policy and final ordering.
    Height buckets are narrowed to width one for this small H<=6 population.
    """
    common=source('elkies-k3/scripts/search_h92_q12o5867_rootless_nagao.py')
    count=protocol['arm_size']
    if len(rows)<count:
        raise ValueError('too few eligible exact classes for the frozen arms')
    survivors=[]
    stages={}
    for coordinate in ('identity','improved'):
        current=[r for r in rows if r['coordinate']==coordinate]
        lookup={tuple(r['box_pair']):r for r in current}
        stages[coordinate]=[]
        for stage,cap in enumerate(protocol['nagao_keep_per_bucket'],1):
            retained,stats=common.retain_bucketed_pareto(
                [candidate_from_row(common,r,stage) for r in current],
                bucket_width=protocol['nagao_bucket_width'],cap_per_bucket=cap)
            current=[lookup[(c.numerator,c.denominator)] for c in retained]
            stages[coordinate].append(stats)
        survivors+=current
    nagao=sorted(survivors,key=lambda r:(common.candidate_sort_key(candidate_from_row(common,r)),r['id']))
    if len(nagao)<count:
        raise ValueError('staged Nagao policy retained fewer than the frozen arm size')
    metrics=('j_log2_height','normalized_weierstrass_log2_size','sample_cost_seconds')
    scores={r['id']:0 for r in rows}
    for metric in metrics:
        # Equal-weight empirical percentile ranks with exact ties receiving
        # the same rank. Measured costs are frozen, never re-timed at replay.
        for r in rows:
            scores[r['id']]+=sum(s[metric]<r[metric] for s in rows)
    lower=sorted(rows,key=lambda r:(scores[r['id']],r['j_log2_height'],r['id']))
    return {'nagao':[r['id'] for r in nagao[:count]],
            'height_cost':[r['id'] for r in lower[:count]],
            'height_cost_rank_sums':scores,
            'nagao_stages':stages,
            'height_only':[r['id'] for r in sorted(rows,key=lambda r:(r['j_log2_height'],r['id']))[:count]]}


def prepare(protocol,output):
    hashes=provenance(protocol)
    checkpoint=LOCAL/'preparation.json.gz'
    payload=read(checkpoint) if checkpoint.exists() else dict(schema='elliptic-curves.fibre-height-population.v1',
        status='PREPARING',inputs=hashes,protocol=protocol,software=dict(sage=version(),python=platform.python_version()),families={})
    if payload['inputs']!=hashes or payload['protocol']!=protocol:
        raise ValueError('preparation checkpoint configuration drift')
    common=source('elkies-k3/scripts/search_h92_q12o5867_rootless_nagao.py')
    exclusions=set(protocol['excluded_j_invariants'])
    for family in families(protocol,sections=True):
        name=family['name']
        if name in payload['families'] and payload['families'][name].get('status')=='FROZEN':
            continue
        frow=payload['families'].setdefault(name,dict(family_id=family['family_id'],
            improved_parameter_scale_shift=list(map(str,family['coordinate'])),rows=[],excluded=[],failures=[],duplicates=[],sample_transcripts={}))
        model=common.FamilyModel(source=TEMPLATE,source_sha256=content_hash([str(family['A']),str(family['B'])]),
            a_coefficients=tuple(map(lambda x:F(str(x)),family['A'])),b_coefficients=tuple(map(lambda x:F(str(x)),family['B'])),
            a_degree=8,b_degree=12,coordinate=name,coefficient_source_keys=('A','B'))
        blocks=protocol['prime_blocks'][name]
        tables,rejected=common.build_residue_tables(model,blocks)
        frow['usable_prime_blocks']=[list(t) for t in tables]
        frow['rejected_primes']=list(rejected)
        seen={r['id'] for key in ('rows','excluded','failures','duplicates') for r in frow[key]}
        iso={}
        for r in frow['rows']:
            iso.setdefault(r['j_invariant'],[]).append((r,EllipticCurve(QQ,r['minimal_model'])))
        jobs=[(coordinate,c) for coordinate in ('identity','improved')
            for c in common.primitive_parameters(protocol['parameter_height'],protocol['parameter_height'])]
        # A fixed hash order interleaves coordinate populations for timing.
        jobs.sort(key=lambda j:content_hash([name,j[0],j[1].numerator,j[1].denominator]))
        for coordinate,c in jobs:
            identifier=f'{name}-{coordinate}-{c.numerator}d{c.denominator}'
            if identifier in seen:
                continue
            basic=dict(id=identifier,coordinate=coordinate,box_pair=[c.numerator,c.denominator])
            if not c.denominator:
                frow['failures'].append({**basic,'reason':'PROJECTIVE_INFINITY_NOT_IMPLEMENTED'})
                continue
            v=QQ(c.numerator)/c.denominator
            scale,shift=family['coordinate'] if coordinate=='improved' else (QQ(1),QQ(0))
            parameter=scale*v+shift
            basic['native_parameter']=str(parameter)
            try:
                E=exact_fibre(family,parameter)
                j=str(E.j_invariant())
                if j in exclusions:
                    frow['excluded'].append({**basic,'j_invariant':j,'reason':'KNOWN_CONTROL_OR_PRIOR_PROSPECTIVE_J_CLASS'})
                    continue
                for prior,curve in iso.get(j,[]):
                    if E.is_isomorphic(curve):
                        frow['duplicates'].append({**basic,'representative':prior['id'],
                            'isomorphism_to_representative':sizes.map_record(E.isomorphism_to(curve))})
                        break
                else:
                    alarm(protocol['minimal_model_seconds'])
                    try:
                        hints=[v for q in (*family['coordinate'],parameter) for v in (q.numerator(),q.denominator())]
                        minimal=remove_obvious_scale(E,hints).minimal_model()
                    finally:
                        cancel_alarm()
                    points=[sizes.read_point(E,p) for p in specialize(family,parameter)]
                    if len(points)!=family['rank'] or not all(points):
                        raise ValueError('specialized subgroup has missing/infinite sections')
                    short=minimal.short_weierstrass_model()
                    phi,points=sizes.transport(E,short,points)
                    row={**basic,'j_invariant':j,'j_log2_height':logheight(E.j_invariant()),
                        'minimal_model':list(map(str,minimal.ainvs())),
                        'minimal_maximum_coefficient_bits':sizes.size(minimal.ainvs())['maximum_bits'],
                        'normalized_weierstrass_size_exact':str(weighted_size(minimal)),
                        'normalized_weierstrass_log2_size':log2(int(weighted_size(minimal)))/12,
                        'search_model':list(map(str,short.ainvs())),
                        'subgroup':[sizes.point_record(p) for p in points],
                        'raw_to_search':sizes.map_record(phi),
                        'raw_to_minimal':sizes.map_record(E.isomorphism_to(minimal)),
                        'section_identities_and_round_trips':len(points),
                        'nagao':score_candidate(common,tables,
                            [family['base'](parameter).numerator(),family['base'](parameter).denominator()],c.height)}
                    samples=[]
                    timings=[]
                    for repeat in range(protocol['sample_repeats']):
                        elapsed=0
                        for index in protocol['sample_section_indices']:
                            started=time.perf_counter()
                            chart=detector.PointedQuarticSearch(curve=row['search_model'],subgroup=[],
                                centre={'point':row['subgroup'][index]},coordinate_policy='metric:16')
                            result=chart.search(protocol['sample_height'],protocol['sample_seconds'])
                            elapsed+=time.perf_counter()-started
                            samples.append(result.record)
                            if result.record['status']!='bounded_search_complete':
                                raise ValueError('COST_SAMPLE_INCOMPLETE')
                        timings.append(elapsed)
                    row['sample_repeat_seconds']=timings
                    row['sample_cost_seconds']=statistics.median(timings)
                    frow['sample_transcripts'][identifier]=samples
                    frow['rows'].append(row)
                    iso.setdefault(j,[]).append((row,minimal))
            except (ValueError,ZeroDivisionError,ArithmeticError,AlarmInterrupt) as exc:
                # Explicit failures affect eligibility and are never called misses.
                frow['failures'].append({**basic,'reason':str(exc),'exception':type(exc).__name__})
            write(checkpoint,payload)
            done=sum(len(frow[k]) for k in ('rows','excluded','failures','duplicates'))
            print(f'POPULATION|family={name}|done={done}/{len(jobs)}|eligible={len(frow["rows"])}',flush=True)
        frow['arms']=select_arms(frow['rows'],protocol)
        frow['status']='FROZEN'
        write(checkpoint,payload)
    payload['status']='FROZEN_BEFORE_PROSPECTIVE_SEARCH'
    payload['selection_digest']=content_hash({k:v['arms'] for k,v in payload['families'].items()})
    write(output,payload)


def independence(model,points,bound):
    model=tuple(F(v) for v in model)
    pts=tuple((F(p['x']),F(p['y'])) for p in points)
    attempts=[]
    for modulus in (2,3,5):
        signatures=find_mod_l_reduction_certificate(model,pts,modulus=modulus,prime_bound=bound)
        rank=combined_mod_l_rank(signatures,len(pts),modulus)
        try:
            prime=find_no_rational_l_torsion_prime(model,modulus=modulus,prime_bound=min(bound,200))
            order=mod_l_reduction_signature(model,(),prime,modulus).group_order
        except ValueError:
            prime=order=None
        attempts.append(dict(modulus=modulus,finite_rank=rank,signatures=[json.loads(json.dumps(asdict(s))) for s in signatures],
            no_torsion_prime=prime,no_torsion_group_order=order))
        if rank==len(pts) and prime is not None:
            return dict(status='CERTIFIED_INDEPENDENT',rank=len(pts),attempts=attempts)
    return dict(status='UNKNOWN',rank=None,attempts=attempts)


def search(population_path,output):
    population=read(population_path)
    protocol=population['protocol']
    if population['status']!='FROZEN_BEFORE_PROSPECTIVE_SEARCH' or population['inputs']!=provenance(protocol):
        raise ValueError('population not frozen or source drift')
    if any(v['arms']!=select_arms(v['rows'],protocol) for v in population['families'].values()):
        raise ValueError('frozen selection changed')
    fresh=dict(schema='elliptic-curves.fibre-height-search.v1',status='SEARCHING',
        population_sha256=digest(population_path),inputs=provenance(protocol),results={})
    payload=read(output) if output.exists() else fresh
    for key in ('schema','population_sha256','inputs'):
        if payload[key]!=fresh[key]:
            raise ValueError('search checkpoint drift: '+key)
    jobs=[]
    for name,family in population['families'].items():
        selected=set(family['arms']['nagao']+family['arms']['height_cost'])
        for row in family['rows']:
            if row['id'] in selected:
                jobs.append((name,row))
    jobs.sort(key=lambda j:content_hash(['execution-v1',j[1]['id']]))
    for name,row in jobs:
        identifier=row['id']
        if identifier in payload['results'] and payload['results'][identifier]['status']!='SEARCHING':
            continue
        started=time.perf_counter()
        result=payload['results'].setdefault(identifier,dict(status='SEARCHING',family=name,charts=[],elapsed_seconds=0))
        rank=len(row['subgroup'])
        if 'generic_independence' not in result:
            result['generic_independence']=independence(row['search_model'],row['subgroup'],protocol['independence_prime_bound'])
            write(output,payload)
        if result['generic_independence']['status']!='CERTIFIED_INDEPENDENT':
            result['status']='GENERIC_INDEPENDENCE_UNKNOWN'
            write(output,payload)
            continue
        centres=protocol['centres'][name]
        for index,rep in enumerate(centres):
            if index<len(result['charts']):
                continue
            chart=detector.PointedQuarticSearch(curve=row['search_model'],subgroup=row['subgroup'],
                centre={'coefficients':rep},coordinate_policy='metric:16')
            transcript=chart.search(protocol['search_height'],protocol['search_seconds'],
                checkpoint_dir=LOCAL/'search-charts').record
            result['charts'].append(transcript)
            write(output,payload)
            print(f'SEARCH|family={name}|candidate={identifier}|chart={index+1}/{len(centres)}|{transcript["status"]}',flush=True)
        # Exact points, then exact independent subgroup extension. Ambiguous
        # finite images are UNKNOWN, never an assertion of dependence.
        E=EllipticCurve(QQ,row['search_model'])
        known={tuple(p.xy()) for p in map(lambda p:sizes.read_point(E,p),row['subgroup'])}
        discovered=sorted({tuple(sizes.read_point(E,p).xy()) for c in result['charts'] for p in c['finite_curve_points']})
        basis=list(row['subgroup'])
        unresolved=[]
        admissions=[]
        for xy in discovered:
            p=E(xy)
            if xy in known or tuple((-p).xy()) in known:
                continue
            record=sizes.point_record(p)
            cert=independence(row['search_model'],basis+[record],protocol['independence_prime_bound'])
            if cert['status']=='CERTIFIED_INDEPENDENT':
                basis.append(record)
                known.add(xy)
                admissions.append(dict(point=record,certificate=cert))
            else:
                unresolved.append(record)
        result.update(certified_new_directions=len(basis)-rank,admissions=admissions,
            unresolved_points=unresolved,finite_points=len(discovered),basis=basis,
            elapsed_seconds=result['elapsed_seconds']+time.perf_counter()-started,
            status='COMPLETE' if all(c['status']=='bounded_search_complete' for c in result['charts']) else 'INCOMPLETE')
        write(output,payload)
    payload['status']='COMPLETE' if all(r['status']=='COMPLETE' for r in payload['results'].values()) and len(payload['results'])==len(jobs) else 'INCOMPLETE'
    write(output,payload)
    return summarize(population,payload)


def describe(values):
    return dict(minimum=min(values),median=statistics.median(values),maximum=max(values))


def summarize(population,results):
    answer=dict(schema='elliptic-curves.fibre-height-comparison.v1',status=results['status'],families={},
        claim_boundary=['Reparameterization changes the sampled fibres; it does not change an existing fibre height.',
            'No bounded null, uncompleted chart, or finite-reduction ambiguity proves absence or a rank upper bound.',
            'Overlapping policy arms share outcomes and are not independent observations.',
            'All-zero arms do not establish equivalence or rule out a useful selection advantage.',
            'MW18 uses a fixed singleton-centre pilot; MW16 retrospective sensitivity does not transfer to MW18.'])
    for name,family in population['families'].items():
        rows={r['id']:r for r in family['rows']}
        arm=family['arms']
        out=dict(eligible_classes=len(rows),excluded_count=len(family['excluded']),failure_count=len(family['failures']),
            duplicate_count=len(family['duplicates']),overlap=sorted(set(arm['nagao'])&set(arm['height_cost'])),
            coordinate_populations={},arms={})
        for coordinate in ('identity','improved'):
            subset=[r for r in rows.values() if r['coordinate']==coordinate]
            out['coordinate_populations'][coordinate]={key:describe([r[key] for r in subset]) for key in
                ('j_log2_height','minimal_maximum_coefficient_bits','normalized_weierstrass_log2_size','sample_cost_seconds')}
            out['coordinate_populations'][coordinate]['count']=len(subset)
        for policy in ('nagao','height_cost'):
            selected=[rows[i] for i in arm[policy]]
            observed=[results['results'][r['id']] for r in selected if r['id'] in results['results']]
            gains=sum(r.get('certified_new_directions',0) for r in observed)
            worker=sum(c['worker_seconds'] for r in observed for c in r['charts'])
            out['arms'][policy]=dict(candidate_ids=arm[policy],candidate_count=len(selected),
                complete_candidates=sum(r['status']=='COMPLETE' for r in observed),
                charts=sum(len(r['charts']) for r in observed),
                complete_charts=sum(c['status']=='bounded_search_complete' for r in observed for c in r['charts']),
                finite_points=sum(r.get('finite_points',0) for r in observed),certified_new_directions=gains,
                worker_seconds=worker,directions_per_worker_hour=gains*3600/worker if worker else None,
                unresolved_points=sum(len(r.get('unresolved_points',[])) for r in observed),
                metrics={k:describe([r[k] for r in selected]) for k in
                    ('j_log2_height','minimal_maximum_coefficient_bits','normalized_weierstrass_log2_size','sample_cost_seconds')})
        out['conclusion']='BOUNDED_COMPARISON_REQUIRES_INTERPRETATION' if any(a['certified_new_directions'] for a in out['arms'].values()) else 'NO_CERTIFIED_GAIN; SELECTION_ADVANTAGE_INCONCLUSIVE'
        answer['families'][name]=out
    return answer


def verify(population_path,result_path):
    population,results=read(population_path),read(result_path)
    protocol=population['protocol']
    assert population['inputs']==provenance(protocol)==results['inputs']
    assert digest(population_path)==results['population_sha256']
    assert population['selection_digest']==content_hash({k:v['arms'] for k,v in population['families'].items()})
    total=0
    common=source('elkies-k3/scripts/search_h92_q12o5867_rootless_nagao.py')
    expected_results=set()
    for family in families(protocol,sections=True):
        f=population['families'][family['name']]
        assert f['arms']==select_arms(f['rows'],protocol)
        expected_results.update(f['arms']['nagao']+f['arms']['height_cost'])
        expected_ids={f'{family["name"]}-{coordinate}-{c.numerator}d{c.denominator}'
            for coordinate in ('identity','improved') for c in common.primitive_parameters(protocol['parameter_height'],protocol['parameter_height'])}
        all_rows=[r for key in ('rows','excluded','failures','duplicates') for r in f[key]]
        assert len(all_rows)==len(expected_ids) and {r['id'] for r in all_rows}==expected_ids
        for row in f['excluded']:
            E=exact_fibre(family,QQ(row['native_parameter']))
            assert str(E.j_invariant())==row['j_invariant'] in protocol['excluded_j_invariants']
        lookup={r['id']:r for r in f['rows']}
        for row in f['duplicates']:
            E=exact_fibre(family,QQ(row['native_parameter']))
            assert E.is_isomorphic(EllipticCurve(QQ,lookup[row['representative']]['minimal_model']))
        model=common.FamilyModel(source=TEMPLATE,source_sha256=content_hash([str(family['A']),str(family['B'])]),
            a_coefficients=tuple(F(str(x)) for x in family['A']),b_coefficients=tuple(F(str(x)) for x in family['B']),
            a_degree=8,b_degree=12,coordinate=family['name'],coefficient_source_keys=('A','B'))
        tables,_=common.build_residue_tables(model,protocol['prime_blocks'][family['name']])
        byj={}
        for row in f['rows']:
            scale,shift=family['coordinate'] if row['coordinate']=='improved' else (QQ(1),QQ(0))
            assert QQ(row['native_parameter'])==scale*QQ(row['box_pair'][0])/row['box_pair'][1]+shift
            E=exact_fibre(family,QQ(row['native_parameter']))
            M,S=EllipticCurve(QQ,row['minimal_model']),EllipticCurve(QQ,row['search_model'])
            assert str(E.j_invariant())==row['j_invariant'] not in protocol['excluded_j_invariants']
            assert E.is_isomorphic(M) and E.is_isomorphic(S)
            assert M.minimal_model().ainvs()==M.ainvs()
            assert row['normalized_weierstrass_size_exact']==str(weighted_size(M))
            assert row['normalized_weierstrass_log2_size']==log2(int(weighted_size(M)))/12
            assert row['minimal_maximum_coefficient_bits']==sizes.size(M.ainvs())['maximum_bits']
            assert row['j_log2_height']==logheight(E.j_invariant())
            t=family['base'](QQ(row['native_parameter']))
            assert row['nagao']==score_candidate(common,tables,[t.numerator(),t.denominator()],max(map(abs,row['box_pair'])))
            assert len(row['sample_repeat_seconds'])==protocol['sample_repeats']
            assert row['sample_cost_seconds']==statistics.median(row['sample_repeat_seconds'])
            assert all(not M.is_isomorphic(prior) for prior in byj.get(row['j_invariant'],[]))
            byj.setdefault(row['j_invariant'],[]).append(M)
            points=[sizes.read_point(E,p) for p in specialize(family,QQ(row['native_parameter']))]
            # Use retained map tuple rather than an arbitrary sign choice.
            from sage.schemes.elliptic_curves.weierstrass_morphism import WeierstrassIsomorphism
            phi=WeierstrassIsomorphism(E,tuple(map(QQ,row['raw_to_search']['u_r_s_t'])),S)
            assert [sizes.point_record(phi(p)) for p in points]==row['subgroup']
            samples=f['sample_transcripts'][row['id']]
            indices=protocol['sample_section_indices']*protocol['sample_repeats']
            assert len(samples)==len(indices)
            for c,index in zip(samples,indices):
                assert c['input']['curve']==row['search_model'] and c['input']['subgroup']==[]
                assert c['input']['centre']['point']==row['subgroup'][index]
                assert c['height_bound']==protocol['sample_height'] and c['timeout_seconds']==protocol['sample_seconds']
                detector.PointedQuarticSearch(**c['input']).verify_record(c)
            if row['id'] not in results['results']:
                continue
            result=results['results'][row['id']]
            assert result['generic_independence']==independence(row['search_model'],row['subgroup'],protocol['independence_prime_bound'])
            for index,c in enumerate(result['charts']):
                assert c['input']['curve']==row['search_model'] and c['input']['subgroup']==row['subgroup']
                assert c['input']['centre']['coefficients']==protocol['centres'][family['name']][index]
                assert c['height_bound']==protocol['search_height']
                assert c['timeout_seconds']==protocol['search_seconds']
                assert c['input']['coordinate_policy']==detector.CoordinatePolicy.parse('metric:16').record()
                detector.PointedQuarticSearch(**c['input']).verify_record(c)
                total+=1
            basis=list(row['subgroup'])
            for admission in result.get('admissions',[]):
                assert any(admission['point'] in c['finite_curve_points'] for c in result['charts'])
                basis.append(admission['point'])
                assert admission['certificate']==independence(row['search_model'],basis,protocol['independence_prime_bound'])
            assert result.get('certified_new_directions',0)==len(basis)-family['rank']
            if result['status']=='COMPLETE':
                assert len(result['charts'])==len(protocol['centres'][family['name']])
                assert all(c['status']=='bounded_search_complete' for c in result['charts'])
            assert result['finite_points']==len({(p['x'],p['y']) for c in result['charts'] for p in c['finite_curve_points']})
    assert set(results['results'])==expected_results
    return dict(status='PASS_EXACT_POPULATION_AND_SEARCH_REPLAY',chart_count=total,
        population_sha256=digest(population_path),results_sha256=digest(result_path))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('stage',choices=('prepare','search','verify'))
    p.add_argument('--population',type=Path,default=POPULATION)
    p.add_argument('--results',type=Path,default=RESULTS)
    args=p.parse_args()
    protocol=read(PROTOCOL)
    for path,expected in protocol['source_hashes'].items():
        if digest(ROOT/path)!=expected:
            raise ValueError('protocol source drift: '+path)
    if args.stage=='prepare':
        prepare(protocol,args.population)
    elif args.stage=='search':
        summary=search(args.population,args.results)
        summary.update(population_sha256=digest(args.population),results_sha256=digest(args.results))
        write(SUMMARY,summary)
    else:
        write(ART/'fibre_height_verification_v1.json',verify(args.population,args.results))


if __name__=='__main__':
    main()
