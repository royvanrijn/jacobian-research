"""Post-seal union certificate; no feedback into either independent route."""
import argparse
from fractions import Fraction as F
from pathlib import Path
from transfer_common import *


def compare(folder):
    payload=read(folder/'input.json')
    class_out,v3=payload['class'],payload['v3']
    model,initial=seed(payload['packet'])
    base=v3.get('points',initial)
    base=[list(map(F,p)) for p in base]
    candidates=class_out.get('points',initial)[16:]
    admission_obj,base_proof=admission(model,base)
    records=[]
    for pt in candidates:
        pt=tuple(map(F,pt));row=admission_obj.consider(pt)
        records.append(dict(point=pt,admission=row))
    points=list(map(list,admission_obj.points))
    _,proof=admission(model,points)
    write(folder/'comparison.json',dict(status='PASS_POST_SEAL_UNION_CERTIFICATE',
        v3_rank_lower_bound=len(base),class_rank_lower_bound=class_out.get('rank_lower_bound',16),
        union_rank_lower_bound=len(points),certified_complementary_gain=len(points)-len(base),
        baseline_output_verified=payload['v3_verified'],class_output_verified=payload['class_verified'],
        points=points,proof=proof,admissions=records,
        boundary='No finite-character failure is a dependence proof. Complementarity is evaluated only against the retained, independently verified V3 subgroup. Missing baseline output is not evidence of efficiency.'))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--folder',type=Path,required=True);args=ap.parse_args();guard(args.folder.resolve());compare(args.folder.resolve())
