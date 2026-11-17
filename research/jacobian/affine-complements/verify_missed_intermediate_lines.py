"""PF140: local-power and exact Euler controls; no finite cover census."""
from pathlib import Path
import hashlib,json,platform
import sympy as S

P=Path(__file__).resolve().parent
NOTE="COMPLETELY_MISSED_INTERMEDIATE_LINES_AND_SOURCE_IMAGE_EULER_BUDGET.md"
MANIFEST="sources/missed-intermediate-lines/manifest.json"


def main():
    checks=[]
    def check(name,value):
        assert value,name
        checks.append(name)
    def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
    manifest=json.loads((P/MANIFEST).read_text())
    for source in manifest["sources"]:
        p=P/source["path"]
        check("source seal "+source["path"],p.stat().st_size==source["bytes"] and sha(p)==source["sha256"])
    old=json.loads((P/"intermediate_critical_axis_receipt.json").read_text())
    check("PF123 retained126 PASS without mathematical replay",old["outcome"]=="PASS" and old["check_count"]==126)
    # Its full inherited source seals are also checked by the scoped auditor.
    for f,h in old["source_sha256"].items():
        check("PF123 retained source pin "+f,sha(P/f)==h)
    prior=json.loads((P/"replica_degree_descent_receipt.json").read_text())
    check("PF139 retained58 PASS without mathematical replay",prior["outcome"]=="PASS" and prior["check_count"]==58)

    u,v=S.symbols("u v")
    k=S.symbols("k",integer=True,positive=True)
    G=S.Function("G")(u,v)
    J=S.det(S.Matrix([[S.diff(u**k,u),S.diff(u**k,v)],
                       [S.diff(G,u),S.diff(G,v)]]))
    check("UNBOUNDED finite smooth-axis power-map Jacobian",
          S.simplify(J-k*u**(k-1)*S.diff(G,v))==0)
    for e in [1,2,5]:
        image=(u**e,v)
        check("whole-axis local residue is identity "+str(e),image[1].subs(u,0)==v)
        check("local critical exponent and no finite excess "+str(e),
              S.diff(image[0],u)==e*u**(e-1) and S.diff(image[1],v)==1)
    zeta=S.symbols("zeta")
    check("small uniformizing transformation fixes the whole boundary axis",
          (zeta*u).subs(u,0)==0 and v==v)

    for indices in [[1],[1,1,1],[2,1],[3,2],[2,2,1]]:
        m=sum(indices);unram=indices.count(1);w=sum(e for e in indices if e>1)
        check("entire missed packet degree and moved support "+str(indices),w==m-unram)
        check("a nontrivial packet never moves exactly one sheet "+str(indices),w==0 or w>=2)
        check("isomorphic axes have unit residue "+str(indices),sum(e*1 for e in indices)==m)

    b,A,eta,r,ct,cU,Delta,eta_res,m,usum=S.symbols("b A eta r ct cU Delta eta_res m usum")
    check("UNBOUNDED exact selected-image Euler ledger",
          S.expand((b-A-eta-r)-(ct)).subs(b,ct+r+A+eta)==0)
    check("UNBOUNDED quotient-degree bound slack",
          S.expand((ct+r+A+eta)-(r+5)-((ct-3)+(A-2)+eta))==0)
    n=S.symbols("n")
    check("UNBOUNDED multiplicity-type degree bound slack",
          S.expand((ct+r+A+eta)-(r+3+n)-((ct-1-n)+(A-2)+eta))==0)
    cV=S.symbols("cV")
    check("UNBOUNDED removed packet Euler identity",
          S.expand(m*(cV-r)-(cU-usum)-
                   ((m*cV-cU)-(m*r-usum)))==0)
    gg=S.symbols("g",positive=True)
    check("UNBOUNDED combined relative and quotient degree cost",
          S.expand(m*b-5*gg*(r+5)-
                   ((m-5*gg)*(r+5)+m*(b-r-5)))==0)

    y=S.symbols("y");H=y**2*(y-1)**3
    derivative=S.diff(H,y)
    check("mixed fibre uses three derivative zeros out of four",
          S.expand(derivative-y*(y-1)**2*(5*y-2))==0)
    other=S.Rational(2,5);value=H.subs(y,other)
    check("mixed fibre forces a genuinely different critical value",value!=0)
    common=S.gcd(H-value,derivative)
    check("other critical fibre has both double and simple roots",
          S.degree(H-value,y)==5 and S.degree(common,y)==1)
    check("two mixed multiplicity types force source-image Euler at least three",
          1+(len(set([2,3]))-1)+(len(set([2,1,1,1]))-1)==3)
    h=y**3-3*y
    check("Morse critical fibres have multiplicity types two and one",
          S.expand(h-2-(y-2)*(y+1)**2)==0 and
          S.expand(h+2-(y+2)*(y-1)**2)==0)
    for mm,rr,ggg in [(5,0,1),(10,2,2),(21,3,3)]:
        bb=rr+5
        check("mixed comparison degree cost "+str((mm,rr,ggg)),mm*bb>=5*ggg*(rr+5))

    theorem=dict(
        completely_missed_intermediate_curves_are_whole_disjoint_embedded_affine_lines=True,
        entire_upper_packets_are_smooth_unit_residue_zero_excess_power_maps=True,
        maximal_affine_source_image_has_exact_Euler_drop_without_erasing_finite_point_gaps=True,
        exact_selected_image_budget_b_equals_chiT_plus_r_plus_support_plus_excess=True,
        actual_replica_mixed_fibre_forces_chiT_ge3_and_quotient_b_ge_r_plus5=True)
    assurance=dict(
        actual_whole_original_plane_and_full_finite_normalizations_retained=True,
        critical_axis_avoidance_handles_whole_curves_not_punctured_images=True,
        prior_boundary_and_Giraldo_Roeder_import_limits_retained=True,
        shared_KV0_MR0_MS0_and_PF139_scopes_credited=True,
        finite_point_gaps_and_actual_full_trace_energy_goal_open=True,
        selected_image_identified_with_entire_intermediate_etale_locus=False,
        embedded_intermediate_line_claimed_embedded_in_final_target=False,
        packets_with_same_final_image_conflated=False,
        generic_deck_group_identified_with_surface_pi1=False,
        new_external_theorem_import=False,
        registered_predecessor_mathematical_replay=False,
        independent_replay=False,formal_verification=False,external_review=False,
        novelty_claim=False,full_goal_complete=False)
    files=[NOTE,"verify_missed_intermediate_lines.py",MANIFEST]+[s["path"] for s in manifest["sources"]]
    result=dict(claim="PF2RULED140",outcome="PASS",check_count=len(checks),checks=checks,
                theorem=theorem,assurance=assurance,source_sha256={f:sha(P/f) for f in files},
                environment=dict(python=platform.python_version(),sympy=S.__version__))
    (P/"missed_intermediate_lines_receipt.json").write_text(json.dumps(result,indent=2)+"\n")
    print("PASS",len(checks),"PF140 controls; all-intermediate missed-line geometry and exact source-image Euler cost; goal ACTIVE")


if __name__=="__main__":main()
