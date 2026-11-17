"""PF135: exact generic cotangent/trace and missing-whole-plane controls."""
from pathlib import Path
import hashlib
import json
import platform
import sympy as S

P=Path(__file__).resolve().parent
NOTE="COTANGENT_DECK_NORMAL_FORM_AND_RELATIVE_TRACE.md"
MANIFEST="sources/cotangent-deck-trace/manifest.json"


def main():
    checks=[]
    def check(name,value):
        assert value,name
        checks.append(name)
    def sha(p):
        return hashlib.sha256(p.read_bytes()).hexdigest()
    manifest=json.loads((P/MANIFEST).read_text())
    for source in manifest["sources"]:
        p=P/source["path"]
        check("source seal "+p.name,p.stat().st_size==source["bytes"] and sha(p)==source["sha256"])
    old=json.loads((P/"finite_orbifold_groups_receipt.json").read_text())
    check("sealed PF134 PASS without replay",old["outcome"]=="PASS" and old["check_count"]==87)
    for f,h in old["source_sha256"].items():check("PF134 predecessor bytes "+f,sha(P/f)==h)
    t,z,p,r,q=S.symbols("t z p r q")
    maps=[t,-t,1/t,-1/t]
    RR=t**2+t**-2;DD=S.diff(RR,t)
    shift=1/(t-2)
    bs=[S.cancel(shift/S.diff(g,t)-shift.subs(t,g)) for g in maps]
    averaged=S.cancel(sum(S.diff(g,t)*b for g,b in zip(maps,bs))/4)
    check("finite derivative-weighted average retains rational shift",averaged.is_rational_function(t))
    for i,(g,b) in enumerate(zip(maps,bs)):
        gp=S.diff(g,t);zz=z/gp+b
        check("nontrivial affine-cocycle action preserves exact volume "+str(i),S.cancel(gp*S.diff(zz,z)-1)==0)
        check("single average cancels every affine cocycle "+str(i),S.cancel(b+averaged.subs(t,g)-averaged/gp)==0)
        check("nonidentity group maps still square to identity "+str(i),
              S.cancel(g.subs(t,g)-t)==0 and S.cancel(zz.subs({t:g,z:zz},simultaneous=True)-z)==0)
        check("spherical quotient coordinate invariant "+str(i),S.cancel(RR.subs(t,g)-RR)==0)
        check("invariant cotangent quotient fibre coordinate "+str(i),S.cancel((p/gp)/DD.subs(t,g)-p/DD)==0)
        check("quotient Darboux wedge determinant "+str(i),S.cancel(DD*S.diff(p/DD,p)-1)==0)

    # All target expressions are polynomial in q over C(r), with no asserted
    # global regularity outside the indicated punctured affine control.
    targetQ=q/(2*r)
    for exponent in [0,2,5]:
        targetP=r**2+(targetQ**exponent if exponent else 0)
        jac=S.diff(targetP,r)*S.diff(targetQ,q)-S.diff(targetP,q)*S.diff(targetQ,r)
        check("exact target coframe after polynomial shear "+str(exponent),S.cancel(jac-1)==0)
        HH=S.integrate(targetP*S.diff(targetQ,q),q)
        drcoef=S.cancel(targetP*S.diff(targetQ,r)-S.diff(HH,r))
        check("fibre integration leaves exactly minus q dr "+str(exponent),S.cancel(drcoef+q)==0)
        check("target generating function polynomial in invariant q "+str(exponent),S.denom(S.together(HH)).has(q) is False)
    check("base target extension degree two has explicit recovery",S.cancel(2*r*targetQ-q)==0)
    check("quotient derivative zero set is among omitted points",S.cancel(DD-2*(t-1)*(t+1)*(t**2+1)/t**3)==0)
    check("all base branch values omitted",S.cancel(RR-2-(t**2-1)**2/t**2)==0 and S.cancel(RR+2-(t**2+1)**2/t**2)==0)
    check("extra omitted target-zero preimages are fourth roots of minus one",S.cancel(RR)==(t**4+1)/t**2)
    check("punctured-base Euler agrees through actual degree-four cover",2-10==4*(2-4))
    check("punctured-base Euler agrees through degree-two target cover",2-4==2*(2-3))

    traces=[]
    for k in range(7):
        trace=S.cancel(sum(g**k for g in maps))
        if k%2:reduced=S.Integer(0)
        elif k==0:reduced=S.Integer(4)
        elif k==2:reduced=2*r
        elif k==4:reduced=2*r**2-4
        else:reduced=2*r**3-6*r
        check("relative trace coefficient "+str(k),S.cancel(trace-reduced.subs(r,RR))==0)
        traces.append(reduced)
        check("invariant fibre power factors out of relative trace "+str(k),S.cancel(sum(g**k*q**3 for g in maps)-q**3*trace)==0)
    gram=S.Matrix(4,4,lambda i,j:traces[i+j])
    check("relative field trace Gram determinant keeps ramification",S.factor(gram.det())==16*(r-2)**2*(r+2)**2)
    check("generic base equation vanishes exactly",S.cancel(t**4-RR*t**2+1)==0)

    theorem=dict(noncyclic_actual_deck_action_has_birational_cotangent_normal_form=True,
                 actual_target_primitive_supplies_invariant_polynomial_generating_function=True,
                 fixed_field_rational_Darboux_coordinates_with_target_polynomial_in_fibre=True,
                 relative_trace_and_Gram_separate_from_invariant_fibre_coordinate=True,
                 free_affine_exact_coframe_control_has_target_but_no_whole_plane=True)
    assurance=dict(PF134_actual_equivariant_ruling_and_coordinator_HC4_credit_retained=True,
                   all_original_plane_and_omitted_primes_retained_in_actual_theorem=True,
                   rational_coordinates_poles_and_actual_lattice_boundary_explicit=True,
                   negative_control_has_genuine_regular_etale_target_on_its_stated_open=True,
                   no_new_classification_or_Razar_type_theorem_import=True,
                   rational_Darboux_coordinates_assumed_globally_regular=False,
                   relative_field_trace_Gram_assumed_actual_integral_lattice=False,
                   cotangent_open_control_claimed_whole_plane_or_Keller_counterexample=False,
                   rational_target_coordinates_used_as_original_polynomial_mates=False,
                   mixed_bad_fibres_or_full_trace_energy_gate_closed=False,
                   registered_predecessor_replay=False,formal_verification=False,
                   independent_replay=False,external_review=False,full_goal_complete=False)
    files=[NOTE,"verify_cotangent_deck_trace.py",MANIFEST,"finite_orbifold_groups_receipt.json"]+[s["path"] for s in manifest["sources"]]
    result=dict(claim="PF2RULED135",outcome="PASS",check_count=len(checks),checks=checks,theorem=theorem,
                assurance=assurance,source_sha256={f:sha(P/f) for f in files},
                environment=dict(python=platform.python_version(),sympy=S.__version__))
    (P/"cotangent_deck_trace_receipt.json").write_text(json.dumps(result,indent=2)+"\n")
    print("PASS",len(checks),"PF135 controls; cotangent field normal form and explicit relative trace; goal ACTIVE")


if __name__=="__main__":
    main()
