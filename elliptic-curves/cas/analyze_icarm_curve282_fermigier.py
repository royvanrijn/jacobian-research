#!/usr/bin/env python3
"""Recognize ICARM curve 282 and certify its 12+8 Fermigier decomposition."""
from __future__ import annotations
import argparse, hashlib, json, sys
from fractions import Fraction as Q
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / "elliptic-curves"), str(Path(__file__).resolve().parent)]
from ecsearch.fermigier_rank import specialize_fermigier_rank_sections
from ecsearch.fermigier import FERMIGIER_DISCRIMINANT_FACTOR_COEFFICIENTS
from ecsearch.fermigier_seed import homogenized_discriminant_factor
from ecsearch.crt_lattice import crt, gauss_reduce_congruence_lattice, hensel_lift_roots
from ecsearch.rank_certification import add_rational_points, select_independent_subset
from elliptic_candidate_record import WeierstrassChange, change_weierstrass_model, source_point_to_target, target_point_to_source, weierstrass_invariants
from fermigier_mestre import FermigierMestreFamily

DEFAULT_INPUT = ROOT / "artifacts/generated-results/elliptic-curves/icarm_7fff_zip_281_282_285_286.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elliptic-curves/icarm_curve282_fermigier_v1.json"
ADAPTER_PARAMETER = Q(11671, 42)
LITERAL_SHIFT = 2 * ADAPTER_PARAMETER
EXPECTED_COMPLEMENT = (12, 13, 14, 15, 17, 18, 19, 20)
QUARTIC_PREIMAGES = (
 (Q(3228660561891889808534,7008481624607014965),Q(-2270787101241256484580503228656764866554298239,114610567592393095607438771471379219525)),
 (Q(477895008887640373,640168315627701),Q(-97096507868487743969494513391905114933,2868708306335254658671469817807)),
 (Q(8805662608547538,18759734723359),Q(1694083356877340335861386357531980951,66514325262361460113185904509)),
 (Q(60919513767045713659042,25631700102204360843),Q(-2242897884103507081954515780427069424487286270899,510987594545044588288160336364784632727)),
 (Q(11095791492540779507816,15537056441453458863),Q(2714844090598482371553935864300530937116352891,241400122864910418377410998106433252769)),
 (Q(34234907015435189996896625611,46236855173002560209894835),Q(-49310356556068129841644955626478578917689467266882043236667,1662769714891610587140075510029218145475254701971175)),
 (Q(147374499138717464059324353985,208565268538529508370613373),Q(17098624442030910760306017070329152261949492050957810854983,4833274582283214085384941590461385219686137805159681)),
 (Q(551075304532045827133118,778502835412214716083),Q(-22736350969837499835533001100563049121383760305969,4242466653214005127580145139357718092040223)),
)
CRT_CONSTRAINTS = ((5,2,13),(11,4,6204),(13,3,1481),(23,2,492),(31,2,255))
EXPECTED_CRT = (282272502437288,408808451805325)
EXPECTED_REDUCED_BASIS = ((11671,42),(126054478,-35027260519))
def txt(q): return str(q.numerator) if q.denominator == 1 else f"{q.numerator}/{q.denominator}"
def valuation(n,p):
 n=abs(n); value=0
 while n and n%p==0: n//=p; value+=1
 return value
def derive_change(source):
 u=Q(1,882); s=Q(-1,2); r=(u*u-source[1]-Q(1,4))/3; return WeierstrassChange(u,r,s,-(1+r)/2)
def main():
 ap=argparse.ArgumentParser(); ap.add_argument("--input",type=Path,default=DEFAULT_INPUT); ap.add_argument("--output",type=Path,default=DEFAULT_OUTPUT); args=ap.parse_args()
 raw=args.input.read_bytes(); data=json.loads(raw); curve=next(x for x in data["curves"] if x["id"]==282)
 target=tuple(Q(x) for x in curve["ainvs"]); submitted=tuple((Q(x),Q(y)) for x,y in curve["points"])
 specialization=specialize_fermigier_rank_sections(ADAPTER_PARAMETER); source=specialization.canonical_model; change=derive_change(source)
 assert change.to_record()==["1/882","-24880960328501059/194481","-1/2","12440480164153289/194481"]
 assert change_weierstrass_model(source,change)==target
 si,ti=weierstrass_invariants(source),weierstrass_invariants(target)
 assert ti["c4"]/si["c4"]==882**4 and ti["c6"]/si["c6"]==882**6 and ti["discriminant"]/si["discriminant"]==882**12
 generic=tuple(source_point_to_target(p,change) for p in specialization.section_differences)
 selected,certificate=select_independent_subset(target,generic+submitted,relation_prime=5,maximum_reduction_prime=1000)
 assert selected[:12]==tuple(range(12)); complement=tuple(i-11 for i in selected[12:]); assert complement==EXPECTED_COMPLEMENT and len(selected)==20
 b2=si["b2"]; c2l=WeierstrassChange(Q(1,6),-b2/12,-source[0]/2,source[0]*b2/24-source[2]/2); literal=FermigierMestreFamily.coefficients(LITERAL_SHIFT)
 assert change_weierstrass_model(source,c2l)==literal
 base=FermigierMestreFamily.known_quartic_points(LITERAL_SHIFT)[0]; base_image=FermigierMestreFamily.quartic_point_to_jacobian(LITERAL_SHIFT,base)
 for i,qpoint in zip(complement,QUARTIC_PREIMAGES,strict=True):
  assert qpoint[1]**2==FermigierMestreFamily.quartic_value(LITERAL_SHIFT,qpoint[0]); direct=FermigierMestreFamily.quartic_point_to_jacobian(LITERAL_SHIFT,qpoint)
  group=source_point_to_target(target_point_to_source(submitted[i-1],change),c2l)
  assert direct==add_rational_points(literal,base_image,add_rational_points(literal,group,group))
 h=homogenized_discriminant_factor(ADAPTER_PARAMETER.numerator,ADAPTER_PARAMETER.denominator)
 local=[]
 for prime,exponent,residue in CRT_CONSTRAINTS:
  modulus=prime**exponent
  assert ADAPTER_PARAMETER.denominator%prime and ADAPTER_PARAMETER.numerator*pow(ADAPTER_PARAMETER.denominator,-1,modulus)%modulus==residue
  assert residue in hensel_lift_roots(FERMIGIER_DISCRIMINANT_FACTOR_COEFFICIENTS,prime,exponent)
  assert valuation(h,prime)==exponent
  local.append({"prime":prime,"exponent":exponent,"modulus":modulus,"parameter_residue":residue})
 residue,modulus=crt((r,p**e) for p,e,r in CRT_CONSTRAINTS); assert (residue,modulus)==EXPECTED_CRT
 reduced=gauss_reduce_congruence_lattice(residue,modulus); assert reduced==EXPECTED_REDUCED_BASIS
 assert reduced[0]==(ADAPTER_PARAMETER.numerator,ADAPTER_PARAMETER.denominator)
 out={"schema":"icarm-curve282-fermigier-recognition-v1","input_sha256":hashlib.sha256(raw).hexdigest(),"curve_id":282,"recognition":{"family":"fermigier-mestre-rank12","adapter_parameter_u":txt(ADAPTER_PARAMETER),"literal_shift_s":txt(LITERAL_SHIFT),"canonical_to_submitted_change_u_r_s_t":change.to_record(),"invariant_scale":{"c4":"882^4","c6":"882^6","discriminant":"882^12"},"q_isomorphic":True},"backward_crt_reconstruction":{"homogenized_discriminant_factor_valuations":local,"crt_residue":residue,"crt_modulus":modulus,"congruence":"a = residue*b (mod modulus)","gauss_reduced_basis":[list(v) for v in reduced],"shortest_vector_is_parameter":True,"comparison_to_frozen_repo_search":{"frozen_primes":[89,131,137],"frozen_profiles":"three-prime exponents 2..4","curve282_primes":[5,11,13,23,31],"present_in_that_frozen_population":False},"interpretation_boundary":"This exactly replays high-power local roots, CRT, and rational reconstruction. It is strong evidence of conductor engineering, but reconstruction from the finished rational parameter alone does not prove the submitter used this algorithm."},"rank_decomposition":{"generic_rank":12,"certified_total_rank":20,"certified_exceptional_quotient_lower_bound":8,"selected_submitted_points_one_based":list(complement),"relation_prime":5,"finite_reduction_certificate":certificate.to_json_object()},"selected_quartic_preimages":[[txt(x),txt(y)] for x,y in QUARTIC_PREIMAGES],"boundary":"The eight directions are certified independent modulo the generic subgroup; no exact specialization-rank upper bound or submitter-method claim is made."}
 args.output.parent.mkdir(parents=True,exist_ok=True); args.output.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n"); print("PASS curve282 = Fermigier(u=11671/42), certified 12+8=20")
if __name__=="__main__": main()
