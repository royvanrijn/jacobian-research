#!/usr/bin/env python3
"""Independent companion-matrix derivation of the shared conic equations."""
import argparse
from pathlib import Path
import retrospective as r
import common_class_conic as run

OUTPUT=r.OUT/'rank_jump_common_class_conic_verification_v1.json'


def compute():
    import sympy as s
    out=r.read(run.OUTPUT)
    for name,sha in out['bindings'].items():assert r.digest((r.ROOT/name).read_bytes())==sha,name
    names='A B h0 h1 h2 z0 z1 z2 l m d w'.split();symbols=s.symbols(' '.join(names));env=dict(zip(names,symbols))
    A,B,h0,h1,h2,z0,z1,z2,l,m,d,w=symbols
    parse=lambda text:s.sympify(text.replace('^','**'),locals=env)
    C=s.Matrix([[0,0,-B],[1,0,-A],[0,1,0]]);I=s.eye(3)
    assert C**3+A*C+B*I==s.zeros(3)
    H=h0*I+h1*C+h2*C**2;Z=z0*I+z1*C+z2*C**2
    q=[s.expand(x) for x in H*Z*s.Matrix([z0,z1,z2])]
    assert q==list(map(parse,out['quadratic_forms']))
    N=s.expand(H.det());assert s.expand(N-parse(out['norm_h']))==0
    assert s.expand(Z.det()-parse(out['norm_z']))==0
    Mq=q[0]*I+q[1]*C+q[2]*C**2
    assert (Mq-H*Z**2).applyfunc(s.expand)==s.zeros(3)
    matrices=[s.Matrix([[parse(x) for x in row] for row in mat]) for mat in out['quadratic_matrices']]
    vector=s.Matrix([z0,z1,z2])
    assert all(s.expand((vector.T*mat*vector)[0]-v)==0 for mat,v in zip(matrices,q))
    assert s.expand(matrices[2].det()+N)==0
    det=s.expand(m*d*(l*matrices[2]+m*matrices[1]).det())
    assert s.expand(det-parse(out['pencil_determinant']))==0
    assert s.expand(det+d*N*m*(l**3+A*l*m*m-B*m**3))==0
    target=s.expand((l*I-d*w*w*C).det())
    assert target==l**3+A*d*d*w**4*l+B*d**3*w**6
    return {'schema':'rank-jump.common-class-conic-verification.v1','status':'PASS',
        'quadratic_forms_verified':3,'matrix_norm_identity_entries':9,'conic_determinant_verified':True,
        'pencil_determinant_verified':True,'twist_point_map_norm_identity_verified':True,
        'method':'SymPy companion matrices, independent of the Sage polynomial multiplication/reduction worker.',
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (Path(__file__),run.OUTPUT,Path(r.__file__))}}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);a=p.parse_args();result=compute()
    if a.mode=='build':r.write_new(OUTPUT,result)
    else:assert result==r.read(OUTPUT)
    print('PASS independent common-conic and twist-cover identities')
