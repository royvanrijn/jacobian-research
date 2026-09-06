#!/usr/bin/env python3
"""Aggregate the finite actual-carrier height audit and small-prime saturation."""
import argparse,json
from dataclasses import asdict
from pathlib import Path
import certify_compact_r17_candidates as cert
import audit_compact_r17_ambiguous as finite
from memory_rank_certificate import checked_rank
from research_runtime.store import checkpoint
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';LOCAL=ROOT/'artifacts/local/elliptic-curves';OUT=ART/'soluble_pair_carrier_gate_v1.json'

def expected():
    paths={Path(__file__).resolve()}
    def read(q):paths.add(q);return cert.read(q)
    h=read(ART/'soluble_pair_carrier_height_v1.json');sat=read(ART/'soluble_pair_generator_saturation_v1.json');anchor=read(ART/'soluble_pair_carrier_anchor_modl_v1.json')
    for data in (h,sat,anchor):
        for n,v in data['sources'].items():
            q=ROOT/n;paths.add(q)
            if cert.hashed(q)!=v:raise ArithmeticError('frozen carrier proof source differs')
    for folder in ['soluble-pair-carrier-height-v1','soluble-pair-generator-saturation-v1','soluble-pair-carrier-height-v1/anchor-modl']:
        for label in ('build','check'):
            s=read(LOCAL/folder/(label+'.supervisor.json'));log=LOCAL/folder/(label+'.log');paths.add(log)
            if s['outcome']!='completed' or s['returncode']!=0 or cert.hashed(log)!=s['log_sha256']:raise ArithmeticError('carrier stage failed/censored or log differs')
    images=[]
    for row in h['rows']:
        if 'curve' not in row:continue
        model=tuple(map(cert.F,row['curve']));points=[tuple(map(cert.F,P)) for P in row['independent_points']];proof=row['rank_certificate']
        actual=checked_rank(model,points,[q['prime'] for q in proof['signatures']],proof['no_rational_2_torsion_prime'])
        if json.loads(json.dumps(actual))!=proof or row['rank_lower_bound']!=len(points):raise ArithmeticError('constructed point proof differs')
        images.append({'word':row['word'],'parameter':row['compact_parameter'],'rank_lower_bound':len(points),'model_coefficient_bits':row['model_coefficient_bits'],'passes_height_gate':row['passes_declared_height_gate']})
    if len(h['rows'])!=12 or len(images)!=11:raise ArithmeticError('fixed twelve-word audit differs')
    small=[r for r in images if r['passes_height_gate']]
    if len(small)!=1 or small[0]['parameter']!='774/149' or small[0]['rank_lower_bound']!=18 or [r['finite_column_rank'] for r in anchor['audits']]!=[18,18]:raise ArithmeticError('known small anchor or finite ambiguity differs')
    model=tuple(map(cert.F,sat['short_model']))
    for a in sat['audits']:
        points=[tuple(map(cert.F,P)) for P in a['columns']];matrix=[]
        for q in a['signatures']:
            sig=finite.signature(model,points,q['prime'],a['modulus'])
            if json.loads(json.dumps(asdict(sig)))!=q:raise ArithmeticError('carrier finite saturation signature differs')
            matrix.extend(sig.rows)
        if len(finite.pivots(matrix,a['modulus']))!=len(points) or a['status']!='PROVED_SATURATED':raise ArithmeticError('finite subgroup injectivity not proved')
        if a['modulus']!=2 and not finite.ml.no_rational_l_torsion_reduction_certificate(model,a['no_rational_ell_torsion_prime'],a['modulus']):raise ArithmeticError('odd torsion exclusion failed')
    return {'schema':'elliptic-curves.soluble-pair-carrier-gate.v1','status':'PASS','sources':{str(q.relative_to(ROOT)):cert.hashed(q) for q in sorted(paths)},'rows':images,'parameter_height_censored':1,'rank19_images':10,'distinct_rank19_parameters':len({r['parameter'] for r in images if r['rank_lower_bound']==19}),'new_images_passing_height_gate':0,'generator_subgroup_saturated_at':[a['modulus'] for a in sat['audits']],'claim_boundary':'Twelve fixed actual-carrier words produce ten certified19-point images on seven distinct parameters, one already known small anchor with supplied-subgroup lower bound18 modulo2,3,5, and one declared height-censored image. No new image meets the360-bit recorded-model gate. The recorded rank2 carrier generators plus rational2-torsion are2,3,5-saturated. Larger-prime saturation, full integral basis, other small carrier images, exact original-fibre ranks and novelty remain open. No original-family point search or automatic expansion occurs.'}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();d=expected()
    if a.check:
        if cert.read(OUT)!=d:raise ArithmeticError('carrier aggregate differs')
    else:
        if OUT.exists():raise FileExistsError('preserve carrier aggregate')
        checkpoint(OUT,d)
    print('ACTUAL CARRIER GATE PASS;',d['distinct_rank19_parameters'],'DISTINCT19-POINT PARAMETERS; NO NEW SMALL MODEL;2,3,5-SATURATED',flush=True)
