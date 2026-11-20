#!/usr/bin/env python3
"""Exact rational support for the compact Airy physical-cut addendum."""
from fractions import Fraction as F
from pathlib import Path
from hashlib import sha256
import argparse,json
if not __debug__:
    raise SystemExit('Unoptimized execution is required')
HERE=Path(__file__).resolve().parent

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,default=HERE/'compact_airy_physical_cut_verification.json');args=ap.parse_args()
    rho,eps=F(3,8),F(1,200000)
    guards={
        'm_epsilon':10**16>200000**3,
        'outer_y':(rho-eps)*F(5,6)>F(3,10),
        'p_lower':F(3,10)/18>3*eps,
        'p_upper_positive_correction':F(17,10)*2>1,
        'small_zeta':4*eps**2==F(1,10**10),
        'cubic_loss':(F(125,486)-F(1,5))*rho**2>4*eps**2,
        'saddle_loss':rho**3/5-F(16,3*10**16)>F(1,100),
        'radial_arc':F(36,25)<F(3,2),
        'tail_constant':5*36000/rho**2<2*10**6,
        'tail_decay':rho**3/5>F(1,40000),
        'right_connector_endpoint_loss':F(1,20000)<F(1,100),
    }
    assert all(guards.values()),[k for k,v in guards.items() if not v]
    # The real cubic on p+i*y: direct exact complex-pair multiplication.
    def add(a,b):
        out=dict(a)
        for k,v in b.items():out[k]=out.get(k,F(0))+v
        return {k:v for k,v in out.items() if v}
    def mul(a,b):
        out={}
        for k,v in a.items():
            for j,w in b.items():
                ij=tuple(i+l for i,l in zip(k,j));out[ij]=out.get(ij,F(0))+v*w
        return {k:v for k,v in out.items() if v}
    def neg(a):return {k:-v for k,v in a.items()}
    p={(1,0):F(1)};y={(0,1):F(1)}
    square_real=add(mul(p,p),neg(mul(y,y)));square_imag=mul({(1,0):F(2)},y)
    cube_real=add(mul(square_real,p),neg(mul(square_imag,y)))
    assert cube_real=={(3,0):F(1),(1,2):F(-3)}
    # A much larger epsilon violates the same concrete cubic-loss guard.
    assert not (F(125,486)-F(1,5))*rho**2>4*F(1,10)**2
    assert not 5*36000/rho**2<10**6
    result={'schema':1,'status':'PASS','scope':'compact canonical Airy physical-cut and weighted-tail rational support only',
            'proof_sha256':sha256((HERE/'COMPACT_AIRY_PHYSICAL_CUT_ADDENDUM.md').read_bytes()).hexdigest(),
            'checker_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
            'rational_guards':guards,'universal_identities':1,'negative_controls':2,
            'formal_verification':False,'external_review':False,'independent_replay':False}
    args.output.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    print(json.dumps({'status':'PASS','guards':len(guards),'universal_identities':1,'negative_controls':2}))
if __name__=='__main__':main()
