"""PF126 transcription, source seals and countercontrols; no finite-sample proof."""
from pathlib import Path
import hashlib,json,platform
import sympy as S
import binary_affine_longitude as B
import seifert_framing_cancellation as C

P = Path(__file__).resolve().parent


def main():
    checks = []
    def check(name, value):
        assert value, name
        checks.append(name)
    def sha(path):
        return hashlib.sha256(path.read_bytes()).hexdigest()
    olds = ["binary_affine_longitude_receipt.json", "intermediate_critical_axis_receipt.json",
            "simple_ramification_ledger_receipt.json"]
    for name in olds:
        r = json.loads((P/name).read_text())
        check("sealed predecessor PASS "+name, r["outcome"] == "PASS")
        for f, digest in r["source_sha256"].items():
            check("sealed predecessor "+name+" "+f, sha(P/f) == digest)
    manifest = json.loads((P/"sources/seifert-framing-cancellation/manifest.json").read_text())
    for item in manifest["sources"]:
        f = P/item["path"]
        check("retained source seal "+f.name, sha(f) == item["sha256"] and f.stat().st_size == item["bytes"])
    t,z,u,e,s = S.symbols("t z u epsilon s")
    multiplicity_controls = []
    for m in [1,2,3,4,7]:
        f = u*z-z**(m+1)/S.Integer(m+1)
        phi = t**m
        x,y = phi,S.expand(f.subs({u:phi,z:t}))
        g,G,w = C.gauss_frame(x,y,t)
        hessian = S.diff(f,z,2).subs({u:phi,z:t})
        check("actual simple-axis Jacobian "+str(m), S.expand(S.diff(f,z)-(u-z**m)) == 0)
        check("actual normalized Gauss Bezout frame "+str(m), S.expand(G[0]*w[1]-G[1]*w[0]-1) == 0)
        check("actual intrinsic Hessian has same zero divisor as derivative "+str(m), S.cancel(hessian/g) == -m)
        res = S.resultant(s**m-t**m,
                          S.expand(y.subs(t,s)-y+e),s)
        polar = S.expand(S.diff(res,e).subs(e,0))
        # The epsilon factor is the diagonal root s=t; its coefficient
        # is the product over all the m-1 off-diagonal roots.
        if m == 1:
            check("smooth branch contributes no self intersection", polar == 1)
        else:
            check("actual branch normal push intersection length "+str(m),
                  S.Poly(polar,t).terms()[-1][0][0] == m*m-1)
        data = C.framing_data([m],[])
        check("actual Seifert correction cancels critical-pair parity "+str(m),
              data["linking"] == m*m-1 and data["preferred_critical_shear"] == 0)
        multiplicity_controls.append(dict(multiplicity=m,hessian=str(hessian),
                                          derivative_gcd=str(g),**data))
    aa = S.symbols("a")
    for m in [2,3,4]:
        yy = t**(m+1)+aa*t**(m+2)
        res = S.resultant(s**m-t**m, yy.subs(t,s)-yy+e,s)
        polar = S.expand(S.diff(res,e).subs(e,0))
        check("higher analytic coefficient retained in local intersection "+str(m),
              S.Poly(polar,t).terms()[-1][0][0] == m*m-1)
    m = S.symbols("m",integer=True,positive=True)
    check("UNBOUNDED framing difference is consecutive-integer product",
          S.factor((m*m-1)-(m-1)) == m*(m-1))
    for multiplicities, contacts in [([1,1],[1]),([2],[3,8]),([2,2,3],[2,5,7]),([4,7],[13])]:
        d = C.framing_data(multiplicities,contacts)
        check("arbitrary cross contacts contribute an even framing change "+str((multiplicities,contacts)),
              (d["linking"]-d["defect"]) % 2 == 0 and d["preferred_critical_shear"] == 0)
    for x,y in [(t**2-1,t**3-t),(t**3-3*t,(t*t-1)**2)]:
        g,G,w = C.gauss_frame(x,y,t)
        check("actual global polynomial normal frame "+str((x,y)), S.expand(G[0]*w[1]-G[1]*w[0]-1) == 0)
    for a in [2,3,4]:
        b = 2**a
        for uu in range(0,b,2):
            p = B.normal_form(a,2,uu,0)
            check("actual preferred longitude localized translation "+str((a,uu)),
                  all(p[x] == x for x in range(1,b,2)) and B.compose(p,p) == list(range(b))
                  and B.cycle_profile(p,range(0,b,2)) == ({2:b//4} if uu else {1:b//2}))
    for label, uu, mask in [("order4",4,4),("mixed",0,4)]:
        p = B.normal_form(3,2,uu,mask)
        check("PF125 "+label+" affine control fails new actual Seifert condition",
              any(p[x] != x for x in range(1,8,2)))
    residual_controls = []
    for a,columns in [(2,[]),(2,[1,1]),(3,[1,2]),(4,[1,1,2,2]),
                      (4,[1,2,3,4]),(5,[1,2,4,8]),(5,[1,2,4,7])]:
        d = C.residual_data(a,columns)
        w,L,ii = d["component_degree"],len(columns),d["infinity_inertia"]
        check("actual elementary two-cover Riemann-Hurwitz "+str((a,columns)),
              2*d["genus"]-2 == -2*w+(w//2)*L+(w-w//ii)
              and w*d["component_count"] == 2**(a-1)
              and d["all_shears_zero"])
        if L == 2:
            check("single odd contact gives Cstar normalization "+str((a,columns)),
                  d["genus"] == 0 and d["infinity_points"] == 2)
        if not L:
            check("no odd contacts split residual A1 covers", w == 1 and d["infinity_translation"] == 0 and d["genus"] == 0)
        residual_controls.append(d)
    check("paired double-value columns need not agree", C.residual_data(3,[1,2])["infinity_translation"] == 3)
    def symplectic(x,y,a):
        return sum(((x>>i)&1)*((y>>(i+1))&1)+((x>>(i+1))&1)*((y>>i)&1)
                   for i in range(0,a,2)) % 2
    for a in [2,4,6]:
        preserving = []
        for uu in range(0,2**a,2):
            linear = lambda x: x ^ (uu if x & 1 else 0)
            if all(symplectic(linear(1<<i),linear(1<<j),a) == symplectic(1<<i,1<<j,a)
                   for i in range(a) for j in range(a)):
                preserving.append(uu)
        check("common symplectic form forces identical nonzero column direction "+str(a), preserving == [0,2])
    for vo in [1,2,5]:
        d = C.residual_data(4,[1]*(2*vo))
        check("actual symplectic residual hyperelliptic profile "+str(vo),
              d["component_degree"] == 2 and d["genus"] == vo-1
              and d["infinity_translation"] == 0 and d["infinity_points"] == 2)
    check("ordinary cusp residual section factorization",
          S.expand((z-t)**2*(z+2*t)-(z**3-3*t*t*z+2*t**3)) == 0)
    text = " ".join((P/"SEIFERT_FRAMING_CANCELLATION_AND_RESIDUAL_TWO_COVERS.md").read_text().lower().split())
    for phrase in ["original finite normalization", "all omitted primes", "other branch components",
                   "every moved letter", "n=2*delta_finite(c)+m_c", "a_j/g",
                   "actual seifert framing", "f*e=0", "order at most two", "every critical shear",
                   "order-four profile", "not assumed that the two columns", "galois cover",
                   "every residual normalization component", "no new external classification",
                   "nondegenerate alternating bilinear form", "u_i=v",
                   "remain open", "full research goal remains active"]:
        check("written scope "+phrase, phrase in text)
    files = ["SEIFERT_FRAMING_CANCELLATION_AND_RESIDUAL_TWO_COVERS.md",
             "seifert_framing_cancellation.py","verify_seifert_framing_cancellation.py",
             "sources/seifert-framing-cancellation/manifest.json"]+olds+[x["path"] for x in manifest["sources"]]
    result = dict(claim="PF2RULED126",outcome="PASS",check_count=len(checks),checks=checks,
                  theorem=dict(all_genuine_intermediates_involutory_componentwise_scope=True,
                               preferred_longitude_fixes_all_moved_letters=True,
                               all_peripheral_shears_zero_for_involutory_components=True,
                               binary_longitude_order_at_most_two_and_only_two_filling_profiles=True,
                               binary_order4_and_mixed_peripheral_classes_excluded=True,
                               binary_residual_components_elementary_abelian_Galois_two_covers=True,
                               odd_contact_label_span_sum_and_exact_genus_formula=True,
                               no_odd_contacts_split_A1_and_one_odd_contact_gives_Cstar=True,
                               common_symplectic_linear_form_forces_trivial_longitude_all_degrees=True),
                  assurance=dict(prior_PF60_PF125_degree4_longitude_result_credited=True,
                                 prior_authors_code_comparison_retained=True,
                                 original_all_prime_and_actual_Seifert_framing_hypotheses_preserved=True,
                                 PF125_algebraic_order4_countercontrol_declared_false=False,
                                 paired_odd_contact_columns_assumed_equal=False,
                                 all_binary_longitudes_assumed_identity=False,
                                 intermediate_plane_recognition=False,
                                 original_integral_coverage_transferred_to_intermediate=False,
                                 linking_entries_determined_by_orbit_profiles=False,
                                 entire_binary_affine_monodromy_family_excluded=False,
                                 actual_source_empty_fibre_excluded=False,
                                 actual_global_finite_algebra_closed=False,
                                 new_external_classification_import=False,new_authors_code_run=False,
                                 finite_samples_prove_theorem=False,formal_verification=False,
                                 independent_replay=False,external_review=False,research_goal_completed=False),
                  multiplicity_controls=multiplicity_controls,residual_controls=residual_controls,
                  source_sha256={f:sha(P/f) for f in dict.fromkeys(files)},
                  environment=dict(python=platform.python_version(),sympy=S.__version__))
    (P/"seifert_framing_cancellation_receipt.json").write_text(json.dumps(result,indent=2)+"\n")
    print("PASS",len(checks),"PF126 controls; actual Seifert-framing cancellation and residual two-cover constraints; goal ACTIVE")


if __name__ == "__main__":
    main()
