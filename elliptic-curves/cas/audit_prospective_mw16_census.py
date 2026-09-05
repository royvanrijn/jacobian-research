#!/usr/bin/env python3
"""Exact representative-norm and selection-binding audit; no CVP optimality claim."""
import argparse
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
import certify_compact_r17_candidates as cert
from research_runtime.store import digest

def audit(directory,output):
    if output.exists():raise FileExistsError('preserve census audit')
    protocol=cert.read(directory/'protocol.json');rows=[]
    for family in protocol['families']:
        path=directory/family/'generic-census.json';d=cert.read(path);scale=d['integer_gram_scale'];gram=[[int(F(x)*scale) for x in row] for row in d['gram']]
        if d['protocol_hash']!=digest(protocol) or d['status']!='COMPLETE_DECLARED_CENSUS' or len(d['records'])!=65536:raise ArithmeticError('census boundary changed')
        for mask,r in enumerate(d['records']):
            v=r['representative']
            if r['mask']!=mask or len(v)!=16 or any((v[i]-(mask>>i))%2 for i in range(16)):raise ArithmeticError('parity enumeration changed')
            if sum(v[i]*gram[i][j]*v[j] for i in range(16) for j in range(16))!=F(r['norm'])*scale:raise ArithmeticError('exact representative norm changed')
        expected=sorted(d['records'][1:],key=lambda r:(-F(r['norm']),r['mask']))[:43]
        if expected!=d['selected'] or dict(Counter(r['norm'] for r in d['records']))!=d['norm_histogram']:raise ArithmeticError('generic class selection changed')
        pop=cert.read(directory/family/'population.json')
        if pop['protocol_hash']!=digest(protocol) or pop['candidate_count']!=1275854 or len(pop['finalists'])!=4:raise ArithmeticError('parameter population changed')
        rows.append({'family':family,'classes_checked':65536,'selected':expected,'norm_histogram':d['norm_histogram'],
                     'generic_census_sha256':cert.hashed(path),'population_sha256':cert.hashed(directory/family/'population.json')})
        print('EXACT MW16 PARITY AUDIT',family,65536,flush=True)
    cert.write(output,{'schema':'elliptic-curves.prospective-mw16-census-audit.v1','checker_sha256':cert.hashed(Path(__file__).resolve()),
        'protocol_sha256':cert.hashed(directory/'protocol.json'),'rows':rows,
        'claim_boundary':'All327680 masks occur exactly once with exact representative norm and parity. No proof that numerical CVP found every true shortest representative; no covering-radius or specialized rank assertion.'})

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--directory',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();audit(a.directory,a.output)
