#!/usr/bin/env sage-python
"""Run both independent real/dyadic replays and bind their exact outputs."""
import hashlib,json,runpy,signal
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def verify():
    paths=[];results={}
    for name,script,certificate in [
        ('dyadic','verify_det1092_rr_dyadic_theta.sage','det1092_rr_dyadic_theta_v1/replay.json'),
        ('real','verify_det1092_rr_real_images.sage','det1092_rr_real_images_v2/replay.json')]:
        p=Path(__file__).with_name(script);c=ART/certificate
        result=runpy.run_path(str(p),run_name='independent_'+name)['verify']()
        assert result==json.loads(c.read_text());results[name]=result;paths.extend([p,c])
    rows=[]
    for dyadic,real in zip(results['dyadic']['cases'],results['real']['cases']):
        assert dyadic['case_index']==real['case_index']
        assert dyadic['real_Kummer_dimension']==real['complete_real_Kummer_dimension']
        rows.append({'case_index':dyadic['case_index'],
                     'Q2_true_Kummer_dimension':dyadic['local_true_Kummer_dimension'],
                     'Q2_fake_Kummer_dimension':dyadic['local_fake_Kummer_dimension'],
                     'D0_in_2JQ2':dyadic['D0_locally_divisible_by_two'],
                     'complete_real_Kummer_dimension':real['complete_real_Kummer_dimension'],
                     'real_image_spanned_by_generic_divisors':True})
    return {'classification':'verified application and new deduction',
            'status':'PASS_INDEPENDENT_REAL_DYADIC_PANEL','cases':rows,
            'scope':'All real images complete. At2, exact dimensions and inherited D0 localization only, not full image generators. Full Selmer and other Sha classes remain uncomputed.',
            'wall_cap_seconds':25,'inputs':{str(p.relative_to(ROOT)):sha(p) for p in paths},
            'checker_sha256':sha(Path(__file__))}
if __name__=='__main__':
    signal.alarm(25);result=verify();out=ART/'det1092_rr_real_dyadic_panel_replay_v1.json'
    payload=json.dumps(result,indent=2,sort_keys=True)+'\n'
    if out.exists():assert out.read_text()==payload
    else:out.write_text(payload)
    print(result['status'],flush=True)
