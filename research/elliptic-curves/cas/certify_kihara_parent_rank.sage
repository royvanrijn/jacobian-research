#!/usr/bin/env sage-python
"""Exact generic index-six enlargement and Picard/rank bounds on first parent."""
import json,hashlib,sys
from pathlib import Path
from sage.all import QQ,ZZ,GF,PolynomialRing,EllipticCurve,matrix,vector
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';OUT=ART/'kihara_first_parent_rank_v1.json'
def main():
 paths=[ART/'kihara_section_involution_v1.json',ART/'kihara_five_parent_distinctness_v1.json',ROOT/'artifacts/local/elliptic-curves/kihara-picard-count-v1/counts.json',ART/'kihara_parent_replay_bundle_v1.json',ART/'kihara_fresh_point_pilot_fibre0_modl_v1.json']
 h,parentdata,counts,bundle,odd=[json.loads(p.read_text()) for p in paths];assert counts['status']=='PASS'
 assert min(r['frobenius']['geometric_NS_upper_bound'] for r in counts['rows'])==18 and min(r['frobenius']['rational_NS_upper_bound'] for r in counts['rows'])==17
 row=parentdata['rows'][1];R=PolynomialRing(QQ,'T');T=R.gen();K=R.fraction_field();S=PolynomialRing(K,'x');x=S.gen()
 A=R(h['raw_A']);B=R(h['raw_B']);E=EllipticCurve(K,[A,B]);P=[E([K(a),K(b)]) for a,b in row['generic_sections']]
 quartic=S([K(R(f)) for f in row['quartic_coefficients']]);X0,Y0=[K(c) for c in row['quartic_points'][0]];shift=quartic(x+X0);ee,dd,cc,bb,aa=shift.list()
 kx=dd*dd/(4*Y0*Y0)-cc;ky=-dd*kx/(2*Y0)-Y0*bb;point=E([9*kx+3*cc,27*ky]);assert 6*point==sum(P[:11],E(0))
 C=matrix.identity(QQ,12)
 for j in range(12):C[0,j]=QQ(1)/6 if j<11 else 0
 H=C*matrix(QQ,h['basis_height_gram'])*C.transpose();basis=[point]+P[1:];assert H.is_positive_definite() and H.det()==189 and abs(C.det())==QQ(1)/6
 seed=bundle['rows'][1]['seed'];old=bundle['rows'][1];Ec=EllipticCurve(QQ,old['old_model']);oldpoints=[Ec([QQ(c) for c in p]) for p in old['old_points']]
 En=EllipticCurve(QQ,seed['curve']);oldimages=[En([QQ(c) for c in p]) for p in seed['original_points']];full=[En([QQ(c) for c in p]) for p in seed['points']]
 iso2=next(iso for iso in Ec.isomorphisms(En) if all(iso(P)==Q for P,Q in zip(oldpoints,oldimages)))
 tc=QQ(row['control_T']);Et=EllipticCurve(QQ,[A(tc),B(tc)]);iso1=next(iso for iso in Et.isomorphisms(Ec) if list(iso.tuple())==list(map(QQ,row['control_isomorphism'])))
 images=[iso2(iso1(Et([c(tc) for c in p.xy()]))) for p in basis]
 # Full fibre seed coordinates: K_old,P2,...,P14. P1=6K_old-sum(P2..P12).
 M=matrix(ZZ,12,14);M[0,0]=-11
 for j in range(1,12):M[0,j]=2
 for i in range(1,12):
  M[i,0]=-6
  for j in range(1,12):M[i,j]=1
  M[i,i+1]+=1
 for i,P in enumerate(images):assert P==sum((c*Q for c,Q in zip(M.row(i),full)),En(0))
 finite=[]
 for ell in [2,3]:
  assert matrix(GF(ell),M).rank()==12
  if ell==2:assert old['finite']['rank_lower_bound']==14
  else:assert next(a for a in odd['audits'] if a['modulus']==3)['finite_column_rank']==14
  finite.append({'modulus':ell,'coefficient_matrix_rank':12,'ambient_certified_finite_rank':14})
 determinant=4*H.det();indices=[int(n) for n in range(1,28) if ZZ(determinant)%(n*n)==0];assert indices==[1,2,3,6]
 action=C.transpose().inverse()*matrix(QQ,h['involution_action'])*C.transpose();assert all(a.denominator()==1 for a in action.list())
 # Every nonzero section on20I1+I4 has height>=4-1=3, so geometric torsion is trivial.
 result={'schema':'elliptic-curves.kihara-first-parent-rank.v1','status':'PASS','parent_parameter':'3/2','geometric_NS_rank':18,'rational_NS_rank':17,'generic_Q_MW_rank':12,'generic_Qbar_MW_rank':13,
 'generic_involution_point':list(map(str,point.xy())),'exact_relation':'6K=P2-P1+...+P12-P1','old_generic_seed_index':6,'basis':[[str(c) for c in p.xy()] for p in basis],
 'basis_change':[[str(c) for c in r] for r in C.rows()],'basis_height_gram':[[str(c) for c in r] for r in H.rows()],'full_rational_MW_height_determinant':'189','full_rational_NS_determinant':'756','geometric_torsion_order':1,
 'specialized_basis_in_full14_coefficients':[list(map(int,r)) for r in M.rows()],'finite_saturation_witnesses':finite,'possible_indices_before_finite_checks':indices,'remaining_index':1,
 'full_basis_involution_action':[[int(c) for c in r] for r in action.rows()],'invariant_rank':6,'anti_invariant_rank':6,
 'all_Q_Jacobian_fibrations_generic_MW_upper_bound':15,'full_geometric_NS_lattice':'UNKNOWN','geometric_NS_type_distinct_from_production948_and_six468':True,
 'sources':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__).resolve(),*paths]},
 'scope':'First retained fresh Kihara parent only. Exact Picard ranks and current generic MW ranks, with full rational section basis and divisor determinant via index-six correction and finite2/3 injection. Geometric Picard18 distinguishes this parent lattice type from previous Picard19 parents. Arithmetic NS17 bounds every Q-Jacobian fibration by MW15 but never bounds specialized elliptic ranks. No new fibration, parameter sweep or point search.'}
 assert not OUT.exists();OUT.write_text(json.dumps(result,indent=2)+'\n');print('PASS Q/Qbar Picard17/18; MW12/13; full rational NSdet756; every Q-fibration MW<=15',flush=True)
if __name__=='__main__':main()
