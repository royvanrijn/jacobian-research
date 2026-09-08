#!/usr/bin/env sage-python
"""Sage matrix-copy compatibility correction; all mathematical inputs frozen."""
import argparse,importlib.machinery,importlib.util,signal
from copy import copy
from pathlib import Path
from sage.all import QQ,identity_matrix
path=Path(__file__).with_name('construct_det1092_rr_real_images.sage')
loader=importlib.machinery.SourceFileLoader('real_images_v1',str(path))
spec=importlib.util.spec_from_loader(loader.name,loader);module=importlib.util.module_from_spec(spec);loader.exec_module(module)
def congruence(H):
    B=copy(H);V=identity_matrix(QQ,H.nrows());n=H.nrows()
    for j in range(n):
        if not B[j,j]:
            options=[k for k in range(j+1,n) if B[k,k]]
            C=identity_matrix(QQ,n)
            if options:
                k=options[0];C.swap_columns(j,k)
            else:
                options=[(k,l) for k in range(j,n) for l in range(k+1,n) if B[k,l]]
                assert options
                k,l=options[0];C.swap_columns(j,k);C[l,j]+=1
            B=C.transpose()*B*C;V=V*C
        assert B[j,j]
        C=identity_matrix(QQ,n)
        for k in range(j+1,n):C[j,k]=-B[j,k]/B[j,j]
        B=C.transpose()*B*C;V=V*C
    assert B.is_diagonal() and V.det() and V.transpose()*H*V==B
    return V,B
module.retain(module.OUT/'failure.json',{'status':'FAILED_BEFORE_FIRST_REAL_IMAGE_CERTIFICATE',
    'case_index':0,'error':'AttributeError: Matrix_rational_dense has no attribute copy',
    'mathematical_result':None,'original_script_sha256':module.sha(path),
    'resolution':'Use standard-library copy(H); same equations, divisor order, limits and mathematics in v2.'})
module.congruence=congruence;module.OUT=module.ART/'det1092_rr_real_images_v2';module.__file__=__file__
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--case',type=int,required=True,choices=range(10));args=ap.parse_args()
    signal.alarm(25);module.construct(args.case)
