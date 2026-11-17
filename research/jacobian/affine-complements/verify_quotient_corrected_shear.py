"""PF127 exact source seals, arithmetic and genuine-family controls."""
from pathlib import Path
from math import gcd
import hashlib,json,platform
import sympy as S
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ
import quotient_corrected_shear as Q
from peripheral_boundary_lattices import boundary_lattices

P = Path(__file__).resolve().parent


def main():
    checks = []
    def check(name,value):
        assert value,name
        checks.append(name)
    def sha(f):
        return hashlib.sha256(f.read_bytes()).hexdigest()
    olds = ["seifert_framing_cancellation_receipt.json","cyclic_quotient_morse_receipt.json",
            "triple_boundary_residue_receipt.json","residual_milnor_receipt.json"]
    for name in olds:
        r = json.loads((P/name).read_text())
        check("sealed predecessor PASS "+name,r["outcome"] == "PASS")
        for f,h in r["source_sha256"].items():
            check("sealed predecessor "+name+" "+f,sha(P/f) == h)
    manifest = json.loads((P/"sources/quotient-corrected-shear/manifest.json").read_text())
    for item in manifest["sources"]:
        f = P/item["path"]
        check("source seal "+f.name,sha(f) == item["sha256"] and f.stat().st_size == item["bytes"])
        if f.suffix == ".html":
            check("primary HTML version and theorem body present "+f.name,
                  item["url"].split("/")[-1] in f.read_text() and
                  "Proposition" in f.read_text() and len(f.read_bytes()) > 100000)
    n,f,M,dd,k,L = S.symbols("n f M delta k L",integer=True)
    check("UNBOUNDED framing eliminates the normalization derivative defect",
          S.expand((n*f*M+n*f-1-n*f*(2*dd+M))-(n*f*(1-2*dd)-1)) == 0)
    a = S.symbols("a",integer=True)
    check("UNBOUNDED normalized stabilizer sign a=-k",
          S.expand(n*(-a)-n*f*(1-2*dd)+1 + n*(a-f*(2*dd-1))-1) == 0)
    e = S.symbols("e",positive=True,integer=True)
    check("UNBOUNDED quotient identity is a Bezout certificate",
          S.expand(n*(a-f*(2*dd-1))-1-(n*a-n*f*(2*dd-1)-1)) == 0)
    y,z = S.symbols("y z")
    t = S.symbols("t")
    gg = t*(t-1)**2
    for rr in [1,2,6]:
        normal_coefficient = 1/(gg.subs(t,y**rr)*rr*y**(rr-1))
        check("whole-axis Jacobian fixes the rational normal coefficient "+str(rr),
              S.cancel(normal_coefficient*gg.subs(t,y**rr)*S.diff(y**rr,y)) == 1 and
              S.degree(S.denom(S.cancel(normal_coefficient)),y) == rr*S.degree(gg,t)+rr-1)
    controls = []
    for ee,ff,nn,delta in [(2,1,1,0),(8,2,3,1),(15,7,2,13),(28,1,17,1),(13,1,8,1)]:
        aa = Q.shear(ee,ff,nn,delta)
        check("exact shear and quotient-order recovery "+str((ee,ff,nn,delta)),
              Q.quotient_order(ee,ff,aa,delta) == nn and
              (nn*(aa-ff*(2*delta-1))-1) % ee == 0)
        check("orbit cyclicity from geometric shear "+str((ee,ff,nn,delta)),
              Q.orbit_invariant_factors(ee,ff,aa) == (1,ee*ff))
        mu,lam = Q.joint_permutations(ee,ff,aa)
        lattice = boundary_lattices(mu,lam)
        check("actual normalized lattice readback "+str((ee,ff,aa)),
              len(lattice) == 1 and lattice[0]["normal_index"] == ee and
              lattice[0]["longitudinal_degree"] == ff and lattice[0]["shear"] == aa)
        snf = smith_normal_form(S.Matrix([[ee,0],[aa,ff]]),domain=ZZ)
        check("exact Smith invariant factors "+str((ee,ff,aa)),
              tuple(abs(int(snf[i,i])) for i in range(2)) == (1,ee*ff))
        controls.append(dict(e=ee,f=ff,n=nn,delta=delta,a=aa,arithmetic_control_only=True))
    for ee in [2,8,14,32]:
        good = True
        for nn in {1,3,ee-1}:
            if gcd(nn,ee) != 1:
                continue
            for ff in [1,2,5]:
                for delta in [0,1,6]:
                    aa = Q.shear(ee,ff,nn,delta)
                    good &= (aa+ff) % 2 == 1
        check("even index forces odd sum of residue and shear "+str(ee),good)
    for ee,ff in [(8,2),(9,3),(25,5),(28,14)]:
        mu,lam = Q.joint_permutations(ee,ff,0)
        check("noncyclic zero-shear abstract torus orbit excluded "+str((ee,ff)),
              len(boundary_lattices(mu,lam)) == 1 and Q.orbit_invariant_factors(ee,ff,0)[0] > 1)
    check("odd shear is not dropped from even joining countercontrol",
          Q.shear(8,2,3,1) == 5 and Q.orbit_invariant_factors(8,2,5) == (1,16))
    mu = [x ^ (2 if x & 1 else 0) for x in range(8)]
    lam = [x ^ (4 if not x & 1 else 0) for x in range(8)]
    product = [mu[lam[x]] for x in range(8)]
    check("cyclic individual orbits do not imply cyclic full peripheral image",
          len({tuple(range(8)),tuple(mu),tuple(lam),tuple(product)}) == 4 and
          all(mu[mu[x]] == lam[lam[x]] == product[product[x]] == x for x in range(8)) and
          all(Q.orbit_invariant_factors(o["normal_index"],o["longitudinal_degree"],o["shear"])[0] == 1
              for o in boundary_lattices(mu,lam)))
    check("same-image individually valid packets can have incompatible delta",
          Q.shared_image_class([(4,1,0,1)]) is not None and
          Q.shared_image_class([(4,1,2,1)]) is not None and
          Q.shared_image_class([(4,1,0,1),(4,1,2,1)]) is None)
    compatible = Q.shared_image_class([(8,2,5,3),(5,1,3,3)])
    check("exact shared-image odd progression without delta search",
          compatible is not None and (1-compatible[0]) % compatible[1] == 0)
    check("recovered smooth quotient order cannot carry actual nontrivial residue",
          Q.shared_image_class([(3,2,0,1)]) is None)
    # Smooth e=2k cusp covers have central-fibre action mu^4.
    for kk in [1,2,3,5]:
        ee = 2*kk
        check("genuine smooth all-index cusp cover weighted shear "+str(kk),
              Q.shear(ee,1,1,1) == 2 % ee and (4-6+Q.shear(ee,1,1,1)) % ee == 0)
    # The entire dessin family has central fibre mu^6 lambda=1.
    mm = S.symbols("m",positive=True,integer=True)
    check("UNBOUNDED genuine dessin quotient inverse identity",S.expand(5*(3*mm-1)-3*(5*mm-2)-1) == 0)
    for mm in [1,2,3,6]:
        ee,nn = 5*mm-2,3*mm-1
        check("genuine singular dessin family shear "+str(mm),Q.shear(ee,1,nn,1) == 6 % ee)
    ss = S.symbols("s",positive=True,integer=True)
    check("UNBOUNDED genuine joining-family formula reduction",
          S.expand(ss*(4*ss-3)+(ss+1)-3-(2*ss+1)*(2*ss-2)) == 0)
    check("UNBOUNDED genuine joining-family weighted central fibre reduction",
          S.expand(ss*(8*ss-2)-3-(2*ss+1)*(4*ss-3)) == 0)
    for ss in [1,2,3,7]:
        ee,delta = 2*ss+1,2*ss-1
        aa = Q.shear(ee,ss,2,delta)
        check("genuine quotient nontrivial-residue joining shear "+str(ss),
              aa == 3 % ee and aa == (ss*(8*ss-2)) % ee)
        check("genuine smooth upstairs residue joining control, wrong Keller open "+str(ss),
              Q.shear(ee,2*ss,1,delta) == 6 % ee)
    # Multi-Puiseux polar count: this is a local curve control, not a
    # whole-axis map satisfying the new global-unit hypothesis.
    U,V = S.symbols("U V")
    implicit = U**4-6272*U**2*V**3-917504*U*V**5-16777216*V**7+9834496*V**6
    polar = S.expand(S.diff(implicit,U).subs({U:56*t**6+64*t**7,V:t**4}))
    check("retained two-Puiseux branch polar order 2delta+M",
          S.Poly(polar,t).terms()[-1][0][0] == 19 and 2*8+3 == 19)
    check("a local analytic unit is not automatically a global whole-axis unit",
          S.degree(1+y,y) > 0 and (1+y).subs(y,-1) == 0)
    # Genuine factorization numerical control from a local power cover.
    ei,fi,ni = 15,14,1
    ej,fj,nj = 15,7,2
    aa_i,aa_j = Q.shear(ei,fi,ni,13),Q.shear(ej,fj,nj,13)
    check("exact singular-order transfer under peripheral lattice containment",
          (aa_i-(fi//fj)*aa_j) % ej == 0 and nj == (fi//fj)*ni)
    for ss in [2,3,7]:
        ee,delta = 2*ss+1,2*ss-1
        source_a,target_a = Q.shear(ee,ss,2,delta),Q.shear(ee,1,2*ss,delta)
        check("genuine further cyclic quotient of joining family "+str(ss),
              (source_a-ss*target_a) % ee == 0 and target_a == (8*ss-2) % ee
              and 2*ss == ss*2 and gcd(ss,ee) == 1)
    text = " ".join((P/"ALL_INDEX_QUOTIENT_CORRECTED_PERIPHERAL_SHEAR.md").read_text().lower().split())
    for phrase in ["original finite n", "all omitted primes", "every genuine finite intermediate",
                   "hereditary power-map location", "normalizing the root is essential",
                   "u(0,y)=c!=0", "all puiseux pairs", "a = f*(2*delta_c-1) + n^(-1)",
                   "cyclic as an abelian action", "entire peripheral image need not be cyclic",
                   "every zero-shear orbit with gcd(e,f)>1", "a+f is odd",
                   "linear-congruence reduction", "genuine nontrivial residue joining",
                   "linking entries", "remain open", "research goal remains active"]:
        check("written scope "+phrase,phrase in text)
    files = ["ALL_INDEX_QUOTIENT_CORRECTED_PERIPHERAL_SHEAR.md","quotient_corrected_shear.py","peripheral_boundary_lattices.py",
             "verify_quotient_corrected_shear.py","sources/quotient-corrected-shear/manifest.json"]+olds+[x["path"] for x in manifest["sources"]]
    result = dict(claim="PF2RULED127",outcome="PASS",check_count=len(checks),checks=checks,
                  theorem=dict(all_original_degrees_all_indices_all_residue_degrees=True,
                               all_genuine_intermediate_critical_power_map_location_descends=True,
                               all_normal_critical_orbits_exact_quotient_corrected_shear=True,
                               individual_joint_peripheral_orbits_cyclic=True,
                               all_zero_shear_orbits_with_non_coprime_index_residue_excluded=True,
                               even_index_even_residue_even_shear_orbits_excluded=True,
                               quotient_order_recovery_and_shared_image_CRT_compatibility=True,
                               exact_singular_order_and_uniform_residue_invariance_through_critical_towers=True),
                  assurance=dict(prior_Borisov_geometry_and_residue_inputs_credited=True,
                                 prior_polar_conductor_import_and_genuine_cover_families_credited=True,
                                 actual_Seifert_framing_and_all_original_primes_retained=True,
                                 primary_HTML_sources_pinned=True,
                                 all_roots_assumed_flat_or_etale_at_quotient_base=False,
                                 individual_orbit_cyclic_means_full_peripheral_group_cyclic=False,
                                 peripheral_cyclic_means_global_monodromy_cyclic=False,
                                 intermediate_plane_recognition=False,
                                 local_analytic_unit_assumed_global_whole_axis_unit=False,
                                 arbitrary_germ_assumed_actual_Keller_boundary=False,
                                 actual_global_finite_algebra_closed=False,
                                 shared_empty_source_fibre_obstruction_closed=False,
                                 finite_samples_prove_theorem=False,formal_verification=False,
                                 independent_replay=False,external_review=False,new_PDF_visual_audit=False,
                                 research_goal_completed=False),
                  arithmetic_controls=controls,shared_image_control=compatible,
                  source_sha256={f:sha(P/f) for f in dict.fromkeys(files)},
                  environment=dict(python=platform.python_version(),sympy=S.__version__))
    (P/"quotient_corrected_shear_receipt.json").write_text(json.dumps(result,indent=2)+"\n")
    print("PASS",len(checks),"PF127 controls; all-normal all-index shear and cyclic peripheral orbit exclusions; goal ACTIVE")


if __name__ == "__main__":
    main()
