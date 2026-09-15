"""Exact inputs for the written CM-orbit doubling construction of a non-CM period."""
from sage.all import *
import json,hashlib,runpy,argparse
from pathlib import Path
R=Path(__file__).resolve().parents[2];P=R/'artifacts/generated-results/elkies-k3-det852-cm-gate-v1'
cm=P/'certificate.json';data=json.loads(cm.read_text());assert data['rational_CM_orders']==[-67,-163] and data['rational_CM_point_count']==8
helper=R/'elkies-k3/scripts/certify_det1236_marked_shimura_curve.sage';h=runpy.run_path(str(helper));fixed={m:h['fixed_point_record'](ZZ(6),ZZ(71),ZZ(m)) for m in [2,3,6,71,142,213]}
records=[]
for m in [2,3,6]:
 count=(fixed[m]['fixed_points']+fixed[426//m]['fixed_points'])//2
 discs=sorted(set(c['order_discriminant'] for d in [m,426//m] for c in fixed[d]['cm_contributions']))
 assert not set(discs).intersection([-67,-163])
 assert count=={2:4,3:4,6:0}[m]
 records.append({'residual_label':m,'geometric_fixed_points':count,'possible_fixed_CM_discriminants':discs})
E=EllipticCurve(QQ,[1,1,0,-286,1780]);rank=E.pari_curve().ellrank();assert int(rank[0])==int(rank[1])==1
T=E.torsion_subgroup().gens()[0].element();assert T!=E(0) and 2*T==E(0)
assert gcd([E.change_ring(GF(p)).cardinality() for p in [5,7,11]])==2
# Modulo2E, translation by T and reflection Q-x give independent translations.
V=VectorSpace(GF(2),2);t=V([0,1]);q=V([1,0]);classes=[list(v) for v in V]
assert len(set(tuple(x) for x in [V.zero(),t,q,q+t]))==4
for x in V:assert len(set(tuple(y) for y in [x,x+t,q-x,q+t-x]))==4
out={'schema':'elkies-k3.det852-intrinsic-nonCM.v1','status':'PASS_INTRINSIC_RATIONAL_NONCM_POINT','residual_involutions':records,'rank_bounds':[int(rank[0]),int(rank[1])],'torsion_order':2,'rational_2torsion_point':[str(z) for z in T],'mod_2E_classes':[[int(z) for z in v] for v in classes],'construction':'Choose any rational CM(-67) point O as origin on C. Let R be the unique rational CM(-163) point whose class R-O lies in2Jac(C)(Q). Then S=O+2(R-O) is a particular rational non-CM point. All four possible choices of O work.','proof':'Residual w6 is translation by nonzero rational2-torsion T. Residual w2,w3 have four geometric fixed points each but no rational fixed points, so their reflection constants Q,Q+T are not in2E(Q). Together T,Q generate E(Q)/2E(Q), of order4. Both four-point CM orbits meet every coset once. In2E(Q) the only CM points are O,R. R is nonzero and non-torsion;2R is neither0 norR, proving S non-CM.','coordinates_on_426b1':None,'inputs':[{'path':str(p.relative_to(R)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in [cm,helper]],'theorem_inputs':['Previously proved full marked curve and complete rational CM locus','Characteristic-zero automorphisms of genus-one curves','Ogg Atkin-Lehner fixed-point formula','Exact PARI two-descent rank bounds'],'boundary':'Intrinsic rational non-CM point certified by CM-orbit and divisibility specification; numerical elliptic coordinates remain uncomputed. Primitive NS and actual rational divisor descent remain to be completed before frame work; no MW17 equation.'}
a=argparse.ArgumentParser();a.add_argument('--check',action='store_true');args=a.parse_args();path=P/'intrinsic-noncm.json'
if args.check:assert json.loads(path.read_text())==out
else:path.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
print('PASS: intrinsic rational non-CM point O+2(R-O); coordinates uncomputed, NS/divisor descent remains.')
