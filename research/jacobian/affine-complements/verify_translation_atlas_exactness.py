"""PF138: small exact atlas/primitive controls, not a surface classification."""
from pathlib import Path
import hashlib
import json
import platform
import sympy as S

P=Path(__file__).resolve().parent
NOTE="TRANSLATION_ATLAS_EXACTNESS_AND_UNBOUNDED_MOVING_LABEL_FIBRES.md"
MANIFEST="sources/translation-atlas-exactness/manifest.json"


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
    old=json.loads((P/"separated_variable_volume_receipt.json").read_text())
    check("sealed PF137 PASS without mathematical replay",old["outcome"]=="PASS" and old["check_count"]==56)
    for f,h in old["source_sha256"].items():check("PF137 predecessor seal "+f,sha(P/f)==h)

    x,w,v,y,z,t,s=S.symbols("x w v y z t s")
    for b,m,roots in [(x,2,[0]),(x*(x-1),3,[0,1]),(x**3-x,2,[-1,0,1])]:
        key=str((b,m));a=b**m;sigma=(1-m)*S.diff(b,x);R=b**(1-m)
        check("square-free prescribed base "+key,S.gcd(b,S.diff(b,x))==1)
        check("moving label has rational primitive "+key,S.cancel(S.diff(R,x)-sigma/a)==0)
        B=S.invert(sigma,b**(m-1),x)
        quotient,remainder=S.div(S.expand(B*sigma-1),b**(m-1),x)
        check("polynomial Bezout primitive data "+key,remainder==0)
        H0=b*B*w;H1=b*B*v+quotient;h=sigma/a
        Fmodel=a*y-z*(z-sigma)
        check("whole chart0 equation "+key,S.expand(Fmodel.subs(
              {z:a*w,y:w*(a*w-sigma)},simultaneous=True))==0)
        check("whole chart1 equation "+key,S.expand(Fmodel.subs(
              {z:sigma+a*v,y:v*(sigma+a*v)},simultaneous=True))==0)
        check("exact overlap z coordinate "+key,S.cancel(sigma+a*(w-h)-a*w)==0)
        check("exact overlap y coordinate "+key,S.cancel((w-h)*(sigma+a*(w-h))-w*(a*w-sigma))==0)
        check("regular polynomial primitives on both charts "+key,
              all(q.is_polynomial(x,w,v) for q in [H0,H1]))
        check("UNBOUNDED-form primitive cocycle control "+key,S.cancel(H1.subs(v,w-h)-H0+R)==0)
        A0=-w+S.diff(H0,x);C0=S.diff(H0,w)
        A1=-v+S.diff(H1,x);C1=S.diff(H1,v)
        check("global primitive overlap dx coefficient "+key,
              S.cancel(A1.subs(v,w-h)-C1.subs(v,w-h)*S.diff(h,x)-A0)==0)
        check("global primitive overlap fibre coefficient "+key,S.cancel(C1.subs(v,w-h)-C0)==0)
        check("chart0 primitive differentiates to unit volume "+key,S.expand(S.diff(C0,x)-S.diff(A0,w)-1)==0)
        check("chart1 primitive differentiates to unit volume "+key,S.expand(S.diff(C1,x)-S.diff(A1,v)-1)==0)
        for root in roots:
            check("two distinct reduced roots over base zero "+str((b,m,root)),sigma.subs(x,root)!=0)
            check("moving-label transition has zero residue "+str((b,m,root)),S.residue(h,x,root)==0)
        check("ruling Euler and actual omitted count "+key,1+S.degree(b,x)==1+len(roots))

    # General equation/flow identities, with no bound on the coefficient degrees.
    a,c=S.symbols("a c")
    F=a*y-z*(z-c)
    delta=lambda f:S.expand(a*S.diff(f,z)+(2*z-c)*S.diff(f,y))
    check("UNBOUNDED additive derivation preserves defining ideal",delta(F)==0)
    check("UNBOUNDED nilpotence on y",delta(delta(y))==2*a and delta(delta(delta(y)))==0)
    check("UNBOUNDED nilpotence on z",delta(z)==a and delta(delta(z))==0)
    moved=F.subs({z:z+t*a,y:y+t*(2*z-c)+t**2*a},simultaneous=True)
    check("UNBOUNDED polynomial flow preserves equation",S.expand(moved-F)==0)
    yt=y+t*(2*z-c)+t**2*a
    check("UNBOUNDED additive group law",S.expand(
          yt.subs({z:z+s*a,y:y+s*(2*z-c)+s**2*a},simultaneous=True)-
          (y+(t+s)*(2*z-c)+(t+s)**2*a))==0)

    # Compatible integration constants on a three-chart cocycle.
    R0,R1,R2=S.Integer(0),1/x,1/(x-1)**2
    check("rational primitive triple cocycle is exact, not just constant",
          S.cancel((R1-R0)+(R2-R1)-(R2-R0))==0)
    check("derivative triple cocycle remains exact",
          S.cancel(S.diff(R1-R0,x)+S.diff(R2-R1,x)-S.diff(R2-R0,x))==0)
    check("nonzero simple-pole residue remains an obstruction",S.residue(1/x,x,0)==1)
    check("moving labels evade the separated-variable obstruction",
          S.residue(1/(x*(x-1))**3,x,0)!=0 and
          S.residue(-2*S.diff(x*(x-1),x)/(x*(x-1))**3,x,0)==0)

    theorem=dict(
        affine_whole_translation_atlas_volume_exact_iff_all_transition_residues_zero=True,
        actual_MQ0_translated_source_cover_satisfies_the_precise_atlas_hypotheses=True,
        arbitrary_prescribed_base_roots_and_exponents_have_explicit_exact_free_additive_models=True,
        whole_plane_scalar_units_simply_connected_Picard_and_Euler_controls_in_all_degrees=True,
        exact_volume_and_free_source_translation_do_not_bound_degenerate_fibre_count=True)
    assurance=dict(
        all_degree_written_proof_not_bounded_scan=True,
        primary_affine_cohomology_statements_imported_with_assurance_limits=True,
        PF137_constant_label_scope_and_prior_atlas_credited=True,
        shared_MP0_MQ0_actual_target_distinctions_retained=True,
        actual_full_finite_order_trace_energy_gate_open=True,
        controls_given_etale_plane_targets=False,
        arbitrary_affine_line_surfaces_claimed_to_have_this_atlas=False,
        original_normalization_replaced_by_selected_opens=False,
        new_external_surface_classification=False,
        registered_predecessor_mathematical_replay=False,
        independent_replay=False,formal_verification=False,external_review=False,
        novelty_claim=False,full_goal_complete=False)
    files=[NOTE,"verify_translation_atlas_exactness.py",MANIFEST]
    files += [q["path"] for q in manifest["sources"]]
    result=dict(claim="PF2RULED138",outcome="PASS",check_count=len(checks),checks=checks,
                theorem=theorem,assurance=assurance,source_sha256={f:sha(P/f) for f in files},
                environment=dict(python=platform.python_version(),sympy=S.__version__))
    (P/"translation_atlas_exactness_receipt.json").write_text(json.dumps(result,indent=2)+"\n")
    print("PASS",len(checks),"PF138 exact controls; whole-atlas exactness and unbounded moving-label fibres; goal ACTIVE")


if __name__=="__main__":main()
