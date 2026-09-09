#!/usr/bin/env sage-python
"""Bounded post-discovery conductor audit; certified local data and factorization."""
import argparse,json,hashlib
from pathlib import Path
from sage.all import EllipticCurve,QQ,ZZ,prime_range,proof
proof.all(True)
def run(packet,output):
    source=json.loads(packet.read_text());E=EllipticCurve([QQ(a) for a in source['curve']])
    if not all(a.denominator()==1 for a in E.a_invariants()):raise ArithmeticError('integral model required')
    D=abs(ZZ(E.discriminant()));remaining=D;known=ZZ(1);rows=[]
    def put(status,**extra):
        data={'status':status,'input':str(packet),'input_sha256':hashlib.sha256(packet.read_bytes()).hexdigest(),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'curve':source['curve'],'rank_lower_bound':source['rank_lower_bound'],'discriminant':str(D),'local_data':rows,'unprocessed_cofactor':str(remaining),'conductor_lower_bound':str(known),'conductor_upper_bound':str(known*remaining),'claim_boundary':'Post-discovery conductor arithmetic only; no record or exact rank claim.',**extra}
        tmp=output.with_suffix('.tmp');tmp.write_text(json.dumps(data,indent=2)+'\n');tmp.replace(output)
    if output.exists():raise FileExistsError('preserve old audit')
    output.parent.mkdir(parents=True,exist_ok=True);put('LOCAL_AUDIT_STARTED')
    for p in prime_range(2,1001):
        e=remaining.valuation(p)
        if not e:continue
        local=E.local_data(p,proof=True);f=int(local.conductor_valuation());remaining//=p**e;known*=p**f
        rows.append({'prime':str(p),'raw_discriminant_valuation':int(e),'conductor_valuation':f});put('CERTIFIED_LOCAL_BOUNDS')
    put('CERTIFIED_LOCAL_BOUNDS_PENDING_FACTORIZATION')
    factors=remaining.factor(proof=True)
    if ZZ.prod(p**e for p,e in factors)!=remaining:raise ArithmeticError('factorization identity differs')
    for p,e in factors:
        if not p.is_prime(proof=True):raise ArithmeticError('unproved prime')
        local=E.local_data(p,proof=True);f=int(local.conductor_valuation());remaining//=p**e;known*=p**f
        rows.append({'prime':str(p),'raw_discriminant_valuation':int(e),'conductor_valuation':f});put('CERTIFIED_LOCAL_BOUNDS')
    if remaining!=1:raise ArithmeticError('incomplete conductor')
    put('EXACT_CONDUCTOR_CERTIFIED',exact_conductor=str(known));print('EXACT_CONDUCTOR',known,flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--packet',required=True,type=Path);p.add_argument('--output',required=True,type=Path);a=p.parse_args();run(a.packet.resolve(),a.output.resolve())
