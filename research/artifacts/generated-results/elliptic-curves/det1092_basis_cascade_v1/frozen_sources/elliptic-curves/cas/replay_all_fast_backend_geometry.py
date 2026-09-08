#!/usr/bin/env python3
"""Replay all602 fixed adaptive PARI geometries and raw points, without admission."""
from pathlib import Path
import certify_compact_r17_candidates as cert
from pointed_quartic_search import PointedQuarticSearch,point_record
from pari_pointed_backend import validate_map,witnesses,program,sources
from continue_fixed_pari_search import paths
KEYS=('base_point','short_model','short_model_x_shift','pointed_chart','horizontal_matrix','ordinate_scale','coefficients')

def main():
    for case in ('rank26','small-conductor'):
        d,parent,maps,*_=paths(case);seed=cert.read(parent/'candidate-00/result.json');data=cert.read(d/'candidate-00/result.json');mapping=cert.read(maps)
        rows=[(i,r) for i,r in enumerate(seed['charts'])]+[(r['parent_chart'],r) for r in data['charts']]
        if [i for i,r in rows]!=list(range(301)):raise ArithmeticError('adaptive roster differs')
        for i,row in rows:
            r=row['search'];m=mapping['rows'][i]
            if row['centre']!=m['centre'] or r['backend']!='pari_fixed_pointed_v1' or r['source_hashes']!=sources():raise ArithmeticError('backend or map binding differs')
            search=PointedQuarticSearch(curve=data['curve'],subgroup=[],centre={'point':r['base_point']},coordinate_policy=m['coordinate_policy']);validate_map(search,m);current=search.chart_record()
            if any(r[k]!=current[k] for k in KEYS):raise ArithmeticError('exact chart geometry differs')
            if r['height_bound']!=100000 or r['program']!=program(m,100000) or not r['infinity_checked'] or r['status']!='bounded_search_complete' or r['returncode']!=0 or '***' in r['stderr']:raise ArithmeticError('fixed successful invocation differs')
            hits,points,ms=witnesses(search,r['stdout'],r['status'],100000)
            if r['primitive_square_hits']!=[[str(v) for v in h] for h in hits] or r['finite_curve_points']!=[point_record(p) for p in points] or r['search_cpu_ms']!=ms:raise ArithmeticError('raw square or mapped point differs')
        print('REPLAYED ALL FAST BACKEND GEOMETRIES',case,len(rows),flush=True)
    print('Geometry and raw points only; each admission history has its separate replayer.',flush=True)
if __name__=='__main__':main()
