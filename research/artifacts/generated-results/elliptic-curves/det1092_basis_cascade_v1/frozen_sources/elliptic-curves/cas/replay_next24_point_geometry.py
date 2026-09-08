#!/usr/bin/env python3
"""Independent generic transport, centre, raw-point and cloud-provenance replay."""
from pathlib import Path
import certify_compact_r17_candidates as cert
import compact_atlas_specialization as spec
import next24_r17_pari_batch as batch
import pari_pointed_backend as backend
from pointed_quartic_search import PointedQuarticSearch,point_record
from half_lattice_pointed_sieve import linear_combination
from research_runtime.store import digest
ROOT=batch.ROOT;ART=batch.ART
KEYS=('base_point','short_model','short_model_x_shift','pointed_chart','horizontal_matrix','ordinate_scale','coefficients')

def tuples(points):return tuple(tuple(map(cert.F,P)) for P in points)
def key(P):return P[0],abs(P[1])

def cloud_check(data,charts,path,input_path):
    cloud=cert.read(path);basis=tuples(data['final_state']['state']['reductions']['points']);points=list(basis);seen={key(P) for P in points}
    for row in charts:
        for raw in row['search']['finite_curve_points']:
            P=cert.F(raw['x']),cert.F(raw['y'])
            if key(P) not in seen:seen.add(key(P));points.append(P)
    if cloud['input_sha256']!=cert.hashed(input_path) or cloud['curve']!=data['curve'] or cloud['points']!=[list(map(str,P)) for P in points]:raise ArithmeticError('retained point-cloud provenance differs')

def geometry(data,maps,initial,p):
    model=tuple(map(cert.F,data['curve']));raw_points=set()
    if [m['centre'] for m in maps['rows']]!=maps['centres'] or data['centres']!=maps['centres']:raise ArithmeticError('map centre roster differs')
    if len(data['charts'])>len(maps['rows']) or data['maps_sha256']!=cert.hashed(p['maps_path']):raise ArithmeticError('map file or chart count differs')
    for index,row in enumerate(data['charts']):
        m=maps['rows'][index];r=row['search'];C=linear_combination(model,initial,m['centre']['representative'])
        if C is None or point_record(C)!=r['base_point'] or row['centre']!=m['centre'] or row['index']!=index:raise ArithmeticError('exact rational centre/roster differs')
        search=PointedQuarticSearch(curve=data['curve'],subgroup=[],centre={'point':r['base_point']},coordinate_policy=m['coordinate_policy']);backend.validate_map(search,m);actual=search.chart_record()
        if any(r[k]!=actual[k] for k in KEYS) or r['backend']!='pari_fixed_pointed_v1' or r['source_hashes']!=backend.sources():raise ArithmeticError('backend geometry differs')
        if r['height_bound']!=p['height'] or r['timeout_seconds']!=p['seconds_per_chart'] or r['gp_binary_sha256']!=p['gp_sha256'] or r['program']!=backend.program(m,p['height']) or not r['infinity_checked']:raise ArithmeticError('retained execution budget differs')
        if r['status']=='bounded_search_complete' and (r['returncode']!=0 or '***' in r['stderr']):raise ArithmeticError('failed call asserted complete')
        hits,points,ms=backend.witnesses(search,r['stdout'],r['status'],p['height'])
        if r['primitive_square_hits']!=[[str(v) for v in h] for h in hits] or r['finite_curve_points']!=[point_record(P) for P in points] or r['search_cpu_ms']!=ms:raise ArithmeticError('raw square/point witnesses differ')
        raw_points.update(key(P) for P in points)
    final=tuples(data['final_state']['state']['reductions']['points'])
    if final[:len(initial)]!=initial or any(key(P) not in raw_points for P in final[len(initial):]):raise ArithmeticError('new basis point was not returned by a chart')

def main():
    p=batch.protocol();atlas={r['family']:r for r in cert.read(spec.ATLAS)['families']};count=0
    for row in p['rows']:
        folder=batch.D/row['id'];path=folder/'result.json';data=cert.read(path);maps=cert.read(folder/'maps.json');original,generic=spec.specialize(atlas[row['family']],row['parameter']);u=cert.F(data['family_to_curve_scale_u']);model=tuple(map(cert.F,data['curve']));initial=tuples(data['initial_state']['state']['reductions']['points'])
        if model!=(0,0,0,original[3]/u**4,original[4]/u**6) or initial!=tuple((x/u**2,y/u**3) for x,y in generic) or tuples(data['generic_points'])!=initial:raise ArithmeticError('exact generic section transport differs')
        if sorted(c['mask'] for c in maps['centres'])!=sorted(p['generic_masks'][row['family']]) or len(data['charts'])!=len(maps['rows']):raise ArithmeticError('full generic chart roster differs')
        if data['protocol_hash']!=digest(p) or maps['protocol_hash']!=digest(p):raise ArithmeticError('generic protocol binding differs')
        for c in maps['centres']:
            if len(c['representative'])!=17 or any((c['representative'][j]-(c['mask']>>j))%2 for j in range(17)):raise ArithmeticError('generic parity differs')
        geometry(data,maps,initial,{**p,'maps_path':folder/'maps.json'});cloud_check(data,data['charts'],ART/('next24_r17_extended_'+row['id'].replace('-','_')+'_mod2_v1.json'),path);count+=len(data['charts'])
    for directory,seed_id,prefix in [('next24-rank27-adaptive-v1','11952-041','next24_rank27')]:
        folder=ROOT/'artifacts/local/elliptic-curves'/directory;ap=cert.read(folder/'protocol.json');data=cert.read(folder/'result.json');maps=cert.read(folder/'maps.json');seed=cert.read(batch.D/seed_id/'result.json');initial=tuples(seed['final_state']['state']['reductions']['points'])
        if data['initial_state']!=seed['final_state'] or any(data[k]!=seed[k] for k in ('family','parameter','curve','generic_points','family_to_curve_scale_u')):raise ArithmeticError('adaptive seed differs')
        if data['protocol_hash']!=digest(ap) or maps['protocol_hash']!=digest(ap):raise ArithmeticError('adaptive protocol binding differs')
        words=sorted(range(1,1024),key=lambda w:(w.bit_count(),w));expected={(r['mask'],words[i%len(words)]) for i,r in enumerate(ap['generic_pool'])}
        if len(maps['centres'])!=301 or {(c['generic_mask'],c['quotient_word']) for c in maps['centres']}!=expected:raise ArithmeticError('adaptive mask roster differs')
        for c in maps['centres']:
            if c['parity']!=(c['generic_mask']|(c['quotient_word']<<17)) or len(c['representative'])!=27 or any((c['representative'][j]-(c['parity']>>j))%2 for j in range(27)):raise ArithmeticError('adaptive parity differs')
        geometry(data,maps,initial,{**ap,'maps_path':folder/'maps.json'});cloud_check(data,seed['charts']+data['charts'],ART/(prefix+'_all_retained_mod2_v1.json'),folder/'all-retained-point-cloud-only.json');count+=len(data['charts'])
    print('REPLAYED NEXT24 GEOMETRY, EXACT CENTRES, RAW POINTS AND CLOUD PROVENANCE',count,'charts',flush=True)
    print('Admission/archive histories passed separate retained local replays; this check does not repeat them.',flush=True)
if __name__=='__main__':main()
