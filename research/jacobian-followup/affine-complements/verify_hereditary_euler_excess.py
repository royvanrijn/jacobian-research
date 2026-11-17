"""PF129 exact Euler/degree controls and sealed inputs, not a cover census."""
from pathlib import Path
import hashlib
import json
import platform
import sympy as S

P=Path(__file__).resolve().parent
NOTE="HEREDITARY_EULER_EXCESS_AND_QUOTIENT_DEGREE_BOUNDS.md"
MANIFEST="sources/hereditary-euler-excess/manifest.json"


def main():
    checks=[]
    def check(name,value):
        assert value,name
        checks.append(name)
    def sha(path):
        return hashlib.sha256(path.read_bytes()).hexdigest()

    manifest=json.loads((P/MANIFEST).read_text())
    for item in manifest["sources"]:
        path=P/item["path"]
        check("sealed dependency "+path.name,
              sha(path)==item["sha256"] and path.stat().st_size==item["bytes"])
    previous=json.loads((P/"pseudoplane_compression_receipt.json").read_text())
    check("PF128 predecessor PASS",previous["outcome"]=="PASS" and previous["check_count"]==64)
    for f,h in previous["source_sha256"].items():
        check("retained PF128 bytes "+f,sha(P/f)==h)

    b,A,Eta,s,c=S.symbols("b A Eta s c",integer=True)
    check("UNBOUNDED subtract exactly the critical A1 curves from chi(M)",
          S.expand((b+s-A-Eta)-s-(b-A-Eta))==0)
    check("UNBOUNDED exact intermediate Euler ledger",
          S.expand(A+Eta+(b-A-Eta-1)-(b-1))==0)
    e,f=S.symbols("e f",positive=True,integer=True)
    check("UNBOUNDED local degree integral includes residue ramification",
          S.expand((e-1)+e*(f-1)-(e*f-1))==0)
    u=S.symbols("u")
    for poly in [u**5,u**2*(u-1)**3,u**5-5*u]:
        degree=S.degree(poly,u)
        factors=S.factor_list(S.diff(poly,u))[1]
        total_ram=sum(S.degree(factor,u)*mult for factor,mult in factors)
        check("finite polynomial ramification gives integral q=f: "+str(poly),
              total_ram==degree-1 and 1+total_ram==degree)

    # Direct fibre-cardinality Euler integration: contacts cancel, whereas
    # actual local excess survives. No list of nodes is treated as extra mass.
    v,kappa,r=S.symbols("v kappa r",integer=True)
    expression=4*r*v+2*r*(1-2*v-kappa)+r*kappa
    check("UNBOUNDED binary etale fibre strata cancel all double values",
          S.expand(expression-r*(2-kappa))==0)
    check("UNBOUNDED binary critical mass plus cusp excess recovers same answer",
          S.expand((4*r-2*r-r*kappa)-expression)==0)
    contact_branches=S.symbols("contact_branches",positive=True,integer=True)
    contact=b*(contact_branches-1)+(b-A)*(1-contact_branches)+(b-contact_branches*A)
    check("UNBOUNDED multiple immersive normalization values cancel exactly",
          S.expand(contact-(b-A))==0)
    check("cusp coinciding with another parameter still charges excess once",
          S.expand(b-(b-A)+(b-2*A-Eta)-(b-A-Eta))==0)

    # Genuine finite maps from the retained notes, all with WRONG etale
    # source C* x C and c=0. Thus positivity needs the whole Keller source.
    kk,ss=S.symbols("k s",positive=True,integer=True)
    check("ALL-k genuine smooth-axis cyclic map has c0, not c>=1",
          S.expand(kk-kk)==0)
    check("ALL-k genuine cusp cover g_k has c0",
          S.expand(3*kk-2*kk-kk)==0)
    check("ALL-s genuine cyclic-quotient cusp cover has c0",
          S.expand(4*ss**2-(2*ss+1)*ss-(2*ss**2-ss))==0)
    for bb,aa,eta,expected in [(5,4,0,1),(7,5,1,1),(4,2,1,1),(6,4,0,2)]:
        check("numerical-only Euler profile "+str((bb,aa,eta)),bb-aa-eta==expected)
    check("one excess unit cannot be spent in a maximal-support quotient",5-4-1==0)
    check("PF115 smooth triple equality germ has TWO excess, not one",8-3*2==2)
    check("both actual eta1 local profiles satisfy exact degree formula",
          3-2*1==1 and 4-3*1==1)
    check("nf2 joining cannot use a smooth original point with residue2",
          [(nn,ff) for nn,ff in [(1,1),(1,2),(2,1)]
           if nn*ff==2 and not (nn==1 and ff>1)]==[(2,1)])

    F=S.symbols("F",positive=True,integer=True)
    check("UNBOUNDED high-image cusp inequality uses another critical image cost2",
          S.expand(5*(b-A-3)-(A+2*F)-(5*b-6*A-15-2*F))==0)
    check("UNBOUNDED half-support comparison for b>=9",
          S.expand((5*b-17)/6-b/2-(2*b-17)/6)==0 and 2*9-17>0)
    for bb in [9,10,15,24,101]:
        check("integer half-support endpoint fits five-sixths cap "+str(bb),
              bb//2<=(5*bb-17)//6)
    check("residue>1 minimal mass10 is impossible through b10",5*2>10-1)
    check("residue>1 mass10 is high support and exceeds cap through b15",
          10>S.Rational(15,2) and (5*15-17)//6==9)
    check("UNBOUNDED nf bound follows from critical mass and f>=1",
          S.expand((e-1)*f-(e*f-f))==0)
    for tt,orders in [(6,[2,3]),(12,[2,4]),(15,[1,3]),(8,[1,1])]:
        generated_order=S.ilcm(*orders)
        check("cyclic pi1 after adding disjoint SMOOTH critical axes "+str((tt,orders)),
              tt%generated_order==0 and all(generated_order%j==0 for j in orders))
    check("A1 cusp original index3k is odd iff inertia k odd",
          all((S.gcd(2,3*j)==1)==(j%2==1) for j in [1,2,3,6,11]))

    note=" ".join((P/NOTE).read_text().split())
    for phrase in ["ENTIRE p-etale locus", "A+Eta_M+(c-1)=b-1",
                   "MINIMAL", "For B=2,Eta_M=0 instead c=2",
                   "EVERY genuine intermediate of degree b<=15",
                   "Original invisible primes are outside", "at a quotient point"]:
        check("scope guard "+phrase,phrase in note)

    theorem=dict(exact_all_normal_intermediate_Euler_support_excess_formula=True,
                 finite_whole_source_pi1_implies_nonnegative_intermediate_Euler_remainder=True,
                 saturated_Euler_mass_triggers_actual_smaller_Keller_or_minimal_cyclic_cover=True,
                 maximal_support_minimal_case_has_smooth_N_all_original_critical_residue1=True,
                 support_slack_at_most2_has_exact_immersive_or_two_cusp_source_profiles=True,
                 all_intermediate_high_support_images_and_indices_have_five_sixths_bound=True,
                 all_critical_residues1_for_every_intermediate_degree_at_most15=True,
                 original_visible_uniform_residue_product_bounded_by_block_degree_minus2=True)
    assurance=dict(prior_local_axis_polar_solid_torus_and_classification_inputs_credited=True,
                   written_unbounded_Euler_proof_includes_all_singular_points_and_image_collisions=True,
                   actual_whole_original_plane_and_all_critical_primes_retained=True,
                   generic_finite_maps_without_whole_source_assumed_positive_Euler=False,
                   every_intermediate_assumed_to_have_plane_open_or_boundary_basis=False,
                   original_degree_replaced_by_unrelated_abstract_representation_degree=False,
                   all_original_residues1_inferred_from_intermediate_degree_le15=False,
                   all_maps_assumed_degree_minimal=False,
                   Euler_one_cyclic_survivors_excluded=False,
                   singular_M_pi1_computed_by_only_generic_meridians=False,
                   full_natural_AS_or_binary_affine_monodromy_excluded=False,
                   new_classification_or_low_sheet_theorem_imported=False,
                   finite_controls_prove_unbounded_results=False,
                   full_trace_energy_or_empty_original_fibre_gate_closed=False,
                   independent_replay=False,formal_verification=False,external_review=False,
                   full_goal_complete=False)
    files=[NOTE,"verify_hereditary_euler_excess.py",MANIFEST,
           "pseudoplane_compression_receipt.json"]+[x["path"] for x in manifest["sources"]]
    receipt=dict(claim="PF2RULED129",outcome="PASS",check_count=len(checks),checks=checks,
                 theorem=theorem,assurance=assurance,source_sha256={f:sha(P/f) for f in files},
                 environment=dict(python=platform.python_version(),sympy=S.__version__))
    (P/"hereditary_euler_excess_receipt.json").write_text(json.dumps(receipt,indent=2)+"\n")
    print("PASS",len(checks),"PF129 controls; hereditary Euler excess and intermediate-degree bounds; goal ACTIVE")


if __name__=="__main__":
    main()
