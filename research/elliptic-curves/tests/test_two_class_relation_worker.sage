#!/usr/bin/env sage-python
"""Bounded functional regression on a small cubic, not a population control."""
from importlib.machinery import SourceFileLoader
from importlib.util import module_from_spec, spec_from_loader
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory

CAS=Path(__file__).resolve().parents[1]/'cas'
sys.path.insert(0,str(CAS))
loader=SourceFileLoader('two_class_worker_test',str(CAS/'two_class_relation_worker.sage'))
spec=spec_from_loader(loader.name,loader)
worker=module_from_spec(spec)
loader.exec_module(worker)

with TemporaryDirectory(prefix='two-class-functional-') as temp:
    root=Path(temp)
    inp=root/'input.json'
    out=root/'outputs'
    out.mkdir()
    payload={'curve_key':'small-cubic-functional-test',
             'protocol':{'factor_base_bound':47,'candidate_budget':64,'candidate_prefixes':[4,16,64]},
             'integral_monic_cubic_ascending':['-1','-1','0','1'],
             'certified_bad_rational_primes':['23'],
             'field_discriminant':'-23','field_signature':[1,1]}
    inp.write_text(json.dumps(payload))
    worker.collect(inp,out)
    result=json.loads((out/'result.json').read_text())
    assert result['noncanonical_relation_count']==33
    assert result['final_matrix']['rank_gain_beyond_canonical']==8
    assert result['final_matrix']['deficiency']==0
    assert result['global_class_2rank_estimate'] is None
    assert result['global_class_2rank_upper_bound'] is None
    assert result['factor_base_generation']=='UNKNOWN_NOT_CERTIFIED'
print('PASS_SMALL_CUBIC_FUNCTIONAL_REGRESSION_NOT_A_GLOBAL_CLASS_RANK_CERTIFICATE')
