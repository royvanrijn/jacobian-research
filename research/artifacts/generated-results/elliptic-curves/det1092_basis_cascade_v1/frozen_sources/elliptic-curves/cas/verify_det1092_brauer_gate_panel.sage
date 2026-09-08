#!/usr/bin/env sage-python
"""Full replay plus exact original302 bad-prime footprint completeness."""
import hashlib,json,runpy,signal
from pathlib import Path
from sage.all import QQ,ZZ,PolynomialRing
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
DIR=ART/'det1092_brauer_section_gate_v3'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def main():
    checker=ROOT/'elliptic-curves/cas/verify_det1092_brauer_section_gate.sage'
    runpy.run_path(str(checker))['verify']()
    parentpath=ART/'curve302_recovered_mw17_parent_v1.json'
    arithpath=ART/'rank_jump_curve302_strict_constructor_arithmetic_inputs_v1.json'
    protocolpath=ART/'det1092_seed_local_code_v5/protocol.json'
    parent,arith,protocol=map(read,[parentpath,arithpath,protocolpath])
    ai=[QQ(d['numerator'][0])/QQ(d['denominator'][0]) for d in parent['a_invariants']]
    assert ai[:3]==[1,1,1]
    R=PolynomialRing(QQ,'x');x=R.gen()
    f=x**3+5*x*x+(16*ai[3]+8)*x+64*ai[4]+16
    D=f.discriminant()/256;assert D in ZZ and D
    factored=ZZ(1);primes=[]
    for p,e in arith['discriminant_factors']:
        p=ZZ(p);assert p.is_prime(proof=True) and e>0
        factored*=p**e;primes.append(int(p))
    assert abs(D)==factored and primes+['infinity']==protocol['places']
    replay=read(DIR/'replay.json')
    assert replay['case_place_pairs']==189 and replay['rational_tree_intersections']==17
    result={'status':'PASS_FULL_BRAUER_SECTION_AND_LOCAL_BLINDNESS_WITH_COMPLETE_302_BAD_PRIME_SUPPORT',
      'classification':'new independently verified global reciprocity boundary',
      'literal302_discriminant':str(D),'all_discriminant_primes':primes,
      'conclusion':'Every normalized surface2-primary Brauer class is identically zero on each of the nine local elliptic fibres at the21 declared places. On302 these include infinity,2 and every discriminant prime. A nonzero evaluation at any rational302 point therefore needs at least two finite primes outside this set, both good for E302. No such class or detection is claimed.',
      'limits':'153 section pairs,17 exact rational intersection locations,189 existing case/place replays; no new control or search fibre.25seconds; no point search, Brauer/Selmer group, later point, V3 input or pilot mutation.',
      'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [checker,parentpath,arithpath,protocolpath,DIR/'replay.json',Path(__file__)]}}
    payload=json.dumps(result,indent=2,sort_keys=True)+'\n';output=DIR/'final-replay.json'
    if output.exists():assert output.read_text()==payload
    else:output.write_text(payload)
    print(result['status'],flush=True)
if __name__=='__main__':
    signal.alarm(25);main()
