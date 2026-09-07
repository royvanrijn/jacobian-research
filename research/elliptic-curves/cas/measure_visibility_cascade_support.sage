#!/usr/bin/env sage-python
"""Recompute the finite atlas at earliest centre support, not choice time."""
from importlib.machinery import SourceFileLoader
from pathlib import Path
CAS=Path(__file__).resolve().parent
b=SourceFileLoader('support_matrix_base',str(CAS/'measure_visibility_cascade.sage')).load_module()
b.D=b.ROOT/'artifacts/local/elliptic-curves/visibility-cascade-support-matrix-v1'
original_read=b.read
original_checkpoint=b.checkpoint

def read(path):
    result=original_read(path)
    if path==b.CHAIN:
        for a in result['recovery_arms']:
            w=a['exact_centre']['representative']
            while len(w)>17 and w[-1]==0:w.pop()
            a['rank_before']=len(w)
    return result

def checkpoint(path,data):
    if path.name=='protocol.json':
        data['atlas']='Original 32 vetted M17 charts per direction plus all historical recovery centres, made available at earliest nonzero basis support, regardless of historical selection date. Previous witnesses persist. One cross-precision CVP translate per chart/stage.'
        data['inputs'][b.rel(Path(__file__))]=b.sha(Path(__file__))
    if path.name=='curve302_visibility_cascade_matrix_v1.json':
        path=b.ART/'curve302_visibility_cascade_support_matrix_v1.json'
        data['inputs'][b.rel(Path(__file__))]=b.sha(Path(__file__))
        data['status']='PASS_EARLIEST_SUPPORT_FINITE_ATLAS_MATRIX'
    original_checkpoint(path,data)

# The base driver checks output existence before writing; redirect just that
# output path while all immutable historical input paths keep their bindings.
class OutputDirectory:
    def __truediv__(self,name):
        return b.ROOT/'artifacts/generated-results/elliptic-curves'/('curve302_visibility_cascade_support_matrix_v1.json' if name=='curve302_visibility_cascade_matrix_v1.json' else name)

def checkpoint_redirect(path,data):
    if path.name=='curve302_visibility_cascade_support_matrix_v1.json':
        data['inputs'][b.rel(Path(__file__))]=b.sha(Path(__file__))
        data['status']='PASS_EARLIEST_SUPPORT_FINITE_ATLAS_MATRIX'
    checkpoint(path,data)

b.read=read;b.checkpoint=checkpoint_redirect;b.ART=OutputDirectory()
if __name__=='__main__':b.run()
