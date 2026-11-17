"""PF137: exact transcription controls; the note proves the all-degree result."""
from pathlib import Path
import hashlib
import json
import platform
import sympy as S

P = Path(__file__).resolve().parent
NOTE = "SEPARATED_VARIABLE_VOLUME_AND_SINGLE_ROOT_OBSTRUCTION.md"
MANIFEST = "sources/separated-variable-volume/manifest.json"


def main():
    checks = []

    def check(name, value):
        assert value, name
        checks.append(name)

    def sha(p):
        return hashlib.sha256(p.read_bytes()).hexdigest()

    manifest = json.loads((P / MANIFEST).read_text())
    for source in manifest["sources"]:
        path = P / source["path"]
        check("retained source seal " + source["path"],
              path.stat().st_size == source["bytes"] and sha(path) == source["sha256"])

    x, w, z = S.symbols("x w z")
    models = [(x**2 * (x-1)**3, z*(z-1), [0, 1]),
              ((x+2)**3, z**3-z, [-1, 0, 1]),
              (x*(x-1), z*(z-1), [0, 1])]
    for model, (a, p, roots) in enumerate(models):
        check("square-free root polynomial " + str(model), S.gcd(p, S.diff(p,z)) == 1)
        for alpha in roots:
            pi = S.div(p, z-alpha, z)[0]
            zi = alpha+a*w
            yi = w*pi.subs(z, zi)
            check("whole polynomial chart equation " + str((model,alpha)),
                  S.expand(a*yi-p.subs(z,zi)) == 0)
            check("regular inverse near retained root " + str((model,alpha)),
                  pi.subs(z,alpha) != 0)
            check("canonical volume chart coefficient " + str((model,alpha)),
                  S.cancel(S.diff(zi,w)/a) == 1)
        alpha, beta = roots[:2]
        shift = (alpha-beta)/a
        check("exact overlap root coordinate " + str(model),
              S.cancel(beta+a*(w+shift)-(alpha+a*w)) == 0)
        p_beta = S.div(p,z-beta,z)[0]
        p_alpha = S.div(p,z-alpha,z)[0]
        check("exact overlap y coordinate " + str(model),
              S.cancel((w+shift)*p_beta.subs(z,beta+a*(w+shift))-
                       w*p_alpha.subs(z,alpha+a*w)) == 0)

    # Finite controls of the residue implication, not a census or its proof.
    residue_cases = [(x, [0]), (x**2-1, [-1,1]),
                     (x**2*(x-1)**3, [0,1]),
                     (x**2*(x**2-1)**2, [0,-1,1])]
    for a, roots in residue_cases:
        residues = [S.residue(1/a,x,b) for b in roots]
        check("at least one nonzero residue in excluded example " + str(a), any(v != 0 for v in residues))
        check("finite residue theorem control " + str(a),
              sum(residues) == (1 if S.degree(a,x)==1 else 0))
    even_a = x**2*(x**2-1)**2
    check("multi-root example may have an individual zero residue",
          S.residue(1/even_a,x,0) == 0 and S.residue(1/even_a,x,1) != 0)

    for c, b, m in [(1,0,2), (3,-2,3), (-2,1,5)]:
        a = c*(x-b)**m
        primitive = (x-b)**(1-m)/S.Integer(c*(1-m))
        check("shifted scaled pure-power rational primitive " + str((c,b,m)),
              S.cancel(S.diff(primitive,x)-1/a) == 0)
        check("pure-power finite residue vanishes " + str((c,b,m)),
              S.residue(1/a,x,b) == 0)

    N, r = S.symbols("N r", integer=True, positive=True)
    check("UNBOUNDED finite pole degree deficit", S.expand((N-1)-(N-r)) == r-1)
    t, c = S.symbols("t c", nonzero=True)
    check("infinity primitive leading derivative for arbitrary degree",
          S.simplify(S.diff(t**(1-N)/(c*(1-N)),t)-t**(-N)/c) == 0)

    # Reuse the prior monomial formula algebraically, without executing its checker.
    m, E, Ep, p, pp = S.symbols("m E Ep p pp")
    D = (1-E*pp)/p
    Hp = (D-Ep)/(m-1)
    check("UNBOUNDED retained monomial primitive identity under affine rescaling",
          S.cancel(E*pp+((m-1)*Hp+Ep)*p-1) == 0)

    a = x**2*(x-1)**3
    check("degree-one p boundary: surface is plane and volume is exact",
          S.cancel(S.diff(a*w,w)/a) == 1)
    # In the PF130 mixed-fibre control the right side genuinely depends on x.
    mixed = (z**2-x)*((z-1)**3-x)
    check("prior exact mixed two-three fibre is outside separated-variable scope",
          S.diff(mixed,x) != 0)

    theorem = dict(
        arbitrary_base_polynomial_has_whole_plane_atlas_and_constant_units=True,
        exact_volume_forces_rational_primitive_of_reciprocal_base_polynomial=True,
        rational_primitive_pole_degree_forces_single_multiple_base_root=True,
        exact_volume_classification_constant_or_single_root_exponent_at_least_two=True,
        no_etale_plane_target_for_multiple_base_roots_or_linear_base=True)
    assurance = dict(
        all_degree_proof_written_not_inferred_from_controls=True,
        prior_classical_atlas_and_monomial_primitive_credited=True,
        same_natural_tower_and_shared_MP0_MQ0_scopes_retained=True,
        x_dependent_mixed_fibres_not_excluded=True,
        surviving_single_root_target_and_full_trace_gate_open=True,
        arbitrary_affine_line_surfaces_classified=False,
        exact_primitive_claimed_to_produce_Keller_pair=False,
        all_individual_residues_claimed_nonzero=False,
        registered_predecessor_replay=False,
        independent_replay=False, formal_verification=False, external_review=False,
        new_external_theorem_import=False, novelty_claim=False, full_goal_complete=False)
    files = [NOTE, "verify_separated_variable_volume.py", MANIFEST]
    files += [source["path"] for source in manifest["sources"]]
    result = dict(claim="PF2RULED137", outcome="PASS", check_count=len(checks),
                  checks=checks, theorem=theorem, assurance=assurance,
                  source_sha256={f:sha(P/f) for f in files},
                  environment=dict(python=platform.python_version(),sympy=S.__version__))
    (P/"separated_variable_volume_receipt.json").write_text(json.dumps(result,indent=2)+"\n")
    print("PASS", len(checks), "PF137 exact controls; all-degree separated-variable exclusion; goal ACTIVE")


if __name__ == "__main__":
    main()
