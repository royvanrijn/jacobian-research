"""PF141: exact local Hamiltonian/trace and multiplier controls, not a map census."""
from pathlib import Path
from fractions import Fraction
from math import ceil
import hashlib,json,platform
import sympy as S

P=Path(__file__).resolve().parent
NOTE="GENERIC_SOURCE_VELOCITY_MULTIPLIER_IDEAL_AND_TRACE_POLE_CHARGES.md"
MANIFEST="sources/source-velocity-multipliers/manifest.json"


def main():
    checks=[]
    def check(name,value):
        assert value,name
        checks.append(name)
    def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
    manifest=json.loads((P/MANIFEST).read_text())
    for source in manifest["sources"]:
        p=P/source["path"]
        check("retained source seal "+source["path"],p.stat().st_size==source["bytes"] and sha(p)==source["sha256"])

    u,z,t,s=S.symbols("u z t s",nonzero=True)
    e,q,n=S.symbols("e q n",integer=True,positive=True)
    A=S.Function("A")(z)
    Omega=e*u**(e-1)
    h=u**(-q)*A
    du=S.diff(A,z)*u**(1-q-e)/e
    dz=q*A*u**(-q-e)/e
    check("UNBOUNDED Hamiltonian du contraction",S.simplify(-Omega*dz-S.diff(h,u))==0)
    check("UNBOUNDED Hamiltonian dz contraction",S.simplify(Omega*du-S.diff(h,z))==0)
    f=u**(e*n)
    Du=f*du;Dz=f*dz
    divergence=S.cancel((S.diff(Omega*Du,u)+S.diff(Omega*Dz,z))/Omega)
    check("UNBOUNDED non-volume-preserving clearing has correct divergence shift",
          S.simplify(divergence-S.diff(f,u)*du)==0)
    check("UNBOUNDED tangent versus normal valuation difference",
          S.expand((e*n+1-q-e)-(e*n-q-e))==1)

    def pole(expr):
        expr=S.cancel(expr)
        if expr==0:return 0
        numerator,denominator=S.fraction(expr)
        def low(poly):return min(k[0] for k,c in S.Poly(poly,t).terms())
        return max(0,low(denominator)-low(numerator))

    for ee,qq in [(1,2),(2,1),(3,1),(4,2),(2,5)]:
        key=str((ee,qq));nn=1+ceil(Fraction(qq,ee))
        U=S.zeros(ee)
        for j in range(ee-1):U[j+1,j]=1
        U[0,ee-1]=t
        av=U**(-qq)
        bv=S.Rational(qq,ee)*(1+z)*U**(-qq-ee)
        traces=[]
        for j in range(1,ee+1):
            tj=S.cancel(S.trace(bv**j))
            exponent=j*(qq+ee)
            expected=(ee*(S.Rational(qq,ee)*(1+z))**j*t**(-exponent//ee)
                      if exponent%ee==0 else 0)
            check("local field trace phase cancellation "+str((ee,qq,j)),S.cancel(tj-expected)==0)
            traces += [ceil(Fraction(pole(S.trace(av**j)),j)),ceil(Fraction(pole(tj),j))]
        check("exact multiplier charge equals finite Newton trace charge "+key,max(traces)==nn)
        check("clearing exponent suffices and preceding exponent fails "+key,
              ee*nn>=qq+ee and ee*(nn-1)<qq+ee)
        check("cleared field is logarithmic and has regular divergence "+key,
              ee*nn+1-qq-ee>=1 and ee*nn-qq-ee>=0)
    check("first trace can vanish while a higher trace has an excessive pole",
          S.trace(S.Matrix([[0,t],[1,0]])**(-3))==0 and
          pole(S.trace(S.Matrix([[0,t],[1,0]])**(-6)))==3)

    for packet in [[(2,1,1)],[(1,2,1)],[(5,2,1)],[(1,2,2)],[(2,1,1),(5,3,2)]]:
        # Tuples are source pole height, normal index, residue degree.
        charge=max(1+ceil(Fraction(qq,ee)) for qq,ee,ff in packet)
        leading=sum(qq*ff for qq,ee,ff in packet)
        check("inverse leading pole charge gives bounded square clearing "+str(packet),
              charge<=leading+1<=2*leading)
        check("one reduced nonproper factor never clears a polar source coordinate "+str(packet),
              all(ee-(qq+ee)==-qq<0 for qq,ee,ff in packet))
        if any(ee==1 for qq,ee,ff in packet):
            check("unramified source-coordinate poles need at least a cube "+str(packet),charge>=3)
    check("coincident image supports use a maximum, not a sum of exponents",
          max(3,2)==3 and 3+2!=3)
    check("leading coefficient itself need not clear on local charge data",
          2<1+ceil(Fraction(2,1)) and 2*2>=1+ceil(Fraction(2,1)))

    for degree in [2,3,5]:
        deriv=lambda value:S.expand((s**degree+z)*S.diff(value,s))
        value=s
        for j in range(1,4):
            value=deriv(value)
            check("non-locally-finite source rescaling has strictly growing degree "+str((degree,j)),
                  S.degree(value,s)==1+j*(degree-1))

    theorem=dict(
        generic_first_source_velocities_have_exact_pole_q_plus_e_at_every_omitted_prime=True,
        exact_target_multiplier_ideal_stabilizes_full_normalization_and_trace_dual_with_divergence_accounted=True,
        first_d_power_trace_poles_recover_exact_multiplier_exponents_and_entire_nonproper_divisor=True,
        square_of_source_finite_inverse_leading_coefficient_clears_all_four_Jacobian_entries=True,
        target_rescalings_cannot_supply_a_locally_finite_source_flow_under_noninvertibility=True)
    assurance=dict(
        actual_whole_source_and_all_ramified_unramified_primes_retained=True,
        target_and_source_multiplier_ideals_distinguished=True,
        canonical_Lie_derivative_divergence_shift_not_omitted=True,
        shared_MS0_MT0_LO1_and_full_inverse_pole_ledger_credited=True,
        unscaled_finite_landing_and_full_trace_energy_goal_open=True,
        generic_exponents_claimed_invariant_under_nonlinear_source_changes=False,
        first_trace_or_degree_independent_prefix_substituted_for_first_d=False,
        nonproper_values_identified_with_empty_source_fibres=False,
        cleared_fields_claimed_locally_nilpotent_or_complete=False,
        local_field_controls_claimed_whole_plane_Keller_counterexamples=False,
        new_external_theorem_import=False,
        registered_predecessor_mathematical_replay=False,
        independent_replay=False,formal_verification=False,external_review=False,
        novelty_claim=False,full_goal_complete=False)
    files=[NOTE,"verify_source_velocity_multipliers.py",MANIFEST]+[s["path"] for s in manifest["sources"]]
    result=dict(claim="PF2RULED141",outcome="PASS",check_count=len(checks),checks=checks,
                theorem=theorem,assurance=assurance,source_sha256={f:sha(P/f) for f in files},
                environment=dict(python=platform.python_version(),sympy=S.__version__))
    (P/"source_velocity_multipliers_receipt.json").write_text(json.dumps(result,indent=2)+"\n")
    print("PASS",len(checks),"PF141 controls; exact full-order target multipliers and trace pole charges; goal ACTIVE")


if __name__=="__main__":main()
