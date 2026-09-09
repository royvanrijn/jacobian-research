# Generic two-fibration seed construction and its exact height mechanism

## Results and scope

**New construction, independently verified.** The existing norm-eight
elliptic pencil produces two explicit rational base changes over **Q**,
using only its generic inherited section and its elliptic group law:

| Alternate multiple | Base | Degree over original parameter | Smooth ramification points | Extra direction's perpendicular height | Certified generic rank |
|---|---|---:|---:|---:|---:|
| `2B` | `P1_z` | 20 | 38 | `399/10` | at least18 |
| `3B` | `P1_z` | 58 | 114 | `3363/29` | at least18 |

These are two separate rank18 constructions, **not** a rank19 construction.
No exceptional point, later302 point, search output, or catalogue rank enters
either construction. The two multiples were fixed before computation.

**New deduction from standard intersection theory.** For any smooth
multisection `C` of degree `d` and genus `g` on this rootless parent, provided
its map is unramified above every singular elliptic fibre, the induced
section has squared height perpendicular to the full inherited subgroup

\[
\boxed{\operatorname{Schur}(C)=2d+\frac{2g-2}{d}.}
\tag{1}
\]

This recovers the previously certified values3 for the rational bisection
and4 for the genus-one bisection. It proves independence for the new
degree20 and58 rational multisections without numerical heights or searches.

**Verified literature application.** The rational rank-jump parameters of
this determinant1092 parent are non-thin, by Pasten--Salgado's theorem and
the already established second elliptic fibration. Finitely many fixed
degree-at-least-two covers therefore cannot exhaust them. This gives no
effective height bound and no prediction of amplification beyond one seed.

## Audit and permitted inputs

**Previously verified applications, reused rather than reproved.** The
parent has full geometric and rational MW17, full rational
`NS=U+(-MW17)` of rank19, determinant1092, and24 `I1` fibres. The
[norm-eight construction](DET1092_GENUS_ONE_FIRST_SEED_COVER_2026-09-08.md)
already supplies a second elliptic fibration over Q, with fibre class
`D=2O+4F+phi(w)`, `D.F=2`, and inherited sections `P14` and `-P15`.
The existing
[Picard-image obstruction](DET1092_GENUS_ONE_PICARD_SPECIFICITY_2026-09-08.md)
proves that its generic sections cannot reach the historical first302 seed.
That obstruction is not contradicted or bypassed here.

The earlier non-thin theorem application in
[the K3 notes](../../elkies-k3/PASTEN_SALGADO_NONTHIN_RANK_JUMPS_2026-08-31.md)
concerns the different published determinant948 surface. We apply the same
established theorem to the independently identified determinant1092 surface.

**Verified input boundary.** Construction reads only the generic parent and
the generic pencil certificate. Although those historical filenames contain
`302`, the fields used are rational functions, the generic trace and inherited
sections, not specialized exceptional coordinates. The exact hashes and
selection are frozen in the
[protocol](../../artifacts/generated-results/elliptic-curves/det1092_pencil_multiples_v2/protocol.json).
The old nine-address control roster is read only after the covers are built.

## Explicit equation-only construction

**Verified setup.** Write the original short model as
`Y^2=X^3+a(t)X+b(t)`. The generic centre is `Z=P14-P15`, with short
coordinates `cx=nx/h^2`, `cy=ny/h^3`. Its degree-two polynomial `h`,
polynomials `nx,ny,shift`, and the branch polynomial

\[
F_z(t)=\frac{(h^2z-\mathrm{shift})^4
 -6n_x(h^2z-\mathrm{shift})^2-8n_y(h^2z-\mathrm{shift})
 -3n_x^2-4ah^4}{h^6}
\]

are fixed by the
[generic certificate](../../artifacts/generated-results/elliptic-curves/det1092_norm8_seed_cover_v2/generic.json).
The quotient is a polynomial, degree four in `t`. Let `(t0(z),s(z))`
be the inherited `P14` point, and expand

\[
F_z(t_0+u)=s^2+q_1u+q_2u^2+q_3u^3+q_4u^4.
\]

All these expressions lie in `Q(z)`. The pointed model is

\[
J:\ y^2=x^3+q_2x^2+A_4x+A_6,
\quad A_4=q_1q_3-4s^2q_4,
\quad A_6=s^2q_3^2+q_1^2q_4-4s^2q_2q_4.
\]

The conjugate inherited point represents

\[
B=(b_x,b_y),\qquad
b_x=q_1^2/(4s^2)-q_2,\quad
b_y=-(q_1b_x+2s^2q_3)/(2s).
\]

**New exact construction.** Compute `2B=(x2,y2)` by

\[
\ell=\frac{3b_x^2+2q_2b_x+A_4}{2b_y},\qquad
x_2=\ell^2-q_2-2b_x,\quad y_2=-b_y+\ell(b_x-x_2).
\]

Compute `3B` by the ordinary chord through `B` and `2B`. For either `n`, set

\[
u_n=\frac{2sy_n+q_1x_n+2s^2q_3}{x_n^2-4s^2q_4},\qquad
T_n=t_0+u_n,\qquad
W_n=\frac{x_nu_n^2}{2s}-s-\frac{q_1u_n}{2s}.
\tag{2}
\]

Then `W_n^2=F_z(T_n)` exactly. Full rational coefficients are retained in
[degree20](../../artifacts/generated-results/elliptic-curves/det1092_pencil_multiples_v2/multiple2.json)
and
[degree58](../../artifacts/generated-results/elliptic-curves/det1092_pencil_multiples_v2/multiple3.json).
The base change is simply `P1_z -> P1_t`, `t=T_n(z)`.

At `t=T_n(z)`, put `m=h*z-shift/h` and

\[
X_n=(hW_n-c_x+m^2)/2,\qquad Y_n=m(X_n-c_x)-c_y,
\]
\[
x_n^{\rm old}=X_n-b_2/12,\qquad
y_n^{\rm old}=Y_n-(a_1x_n^{\rm old}+a_3)/2.
\tag{3}
\]

These factored formulas define the new point on the original elliptic
family. The verifier checks the universal elliptic identity, the trace's
elliptic equation, and the exact quartic identity. It also certifies the
nonexceptional image and recovery of `z` by the forward pencil map.
Thus `P1_z` is the normalization of the actual multisection, not an
artificial redundant parameter cover. Its image is an alternate-fibration
section on the smooth K3, hence a smooth rational curve.

## Independent ramification and rank certificates

**Verified exact geometry.** Write `T_n=N/D` in lowest terms. Its critical
polynomial `R=N'D-ND'` is squarefree of degree `2d-2`. If `Delta(t)` is the
original degree24 discriminant, the certificates prove

\[
\gcd\bigl(R,\ D\,D^{24}\Delta(N/D)\bigr)=1.
\tag{4}
\]

For `n=2` this is certified modulo127; for `n=3`, modulo131. Degrees
are preserved and the reductions are squarefree; explicit Bezout identities
are verified. Consequently these are proofs over Q, not probabilistic tests.
Infinity is unramified because the critical polynomial has full degree
`2d-2`. All ramification is over smooth original fibres; in particular
the covers are unramified above every original bad fibre.

**Established literature, verified application.** A multisection ramified
above a smooth fibre gives an independent section after base change; see
[Garbagnati--Salgado, Lemma2.9](https://www.cambridge.org/core/journals/forum-of-mathematics-sigma/article/rank-jumps-and-multisections-of-elliptic-fibrations-on-k3-surfaces/602807FB00055D106E3CEA5418DE08F7).
Indeed, if a nonzero multiple descended to the old base, the multisection
would be a component of its finite etale division-point scheme over the
smooth locus, contradicting ramification.

**New height derivation.** There are no reducible-fibre correction terms
after these base changes. Write

\[
[C]=dO+aF+\phi(w),\quad C^2=-2d^2+2da-\langle w,w\rangle=2g-2.
\]

The pulled-back arithmetic genus is `chi=2d`. If `Q` is the induced section,
`Q.O=C.O=a-2d`, hence `h(Q)=4d+2(C.O)=2a`. For an old section `S_v`,

\[
[S_v]=O+\tfrac12\langle v,v\rangle F+\phi(v).
\]

Computing `C.S_v` and the section pairing gives
`<Q,S_v^*>=<w,v>`. The inherited Gram is `dG`, so

\[
\operatorname{Schur}(Q)=2a-\frac{\langle w,w\rangle}{d}
=2d+\frac{C^2}{d}=2d+\frac{2g-2}{d}.
\]

This calculation uses the established intersection-height formulas and
degree scaling in
[Schuett--Shioda, sections11.8 and11.13](https://arxiv.org/pdf/0907.0298).
For the two new covers the rank18 regulators are exactly

```text
degree20: 571091189760000000000000000
degree58: 120455711812198450108549951693258752.
```

The symbolic
[height replay](../../artifacts/generated-results/elliptic-curves/det1092_pencil_multiples_v2/height-identity.json)
checks these determinants and the genus0/1/2 degree-two regressions3/4/5.
For a singular multisection or ramification above bad fibres, formula(1)
cannot be reused without addressing the changed intersection/correction terms.

## Rational points, explicit exclusions, and comparison with the old route

**Verified prospective source.** Both bases are `P1` over Q: every rational
`z` is available without a rational-point search or number-field extension.
Equations(2)--(3) explicitly exclude their denominator zeros, `D(z)=0`,
`h(T_n(z))=0`, and `Delta(T_n(z))=0`, and any poles of the inherited
section coordinates used for an affine packet. This is an equation-defined
finite chart-exclusion set; the projective maps can extend beyond it.

**Established specialization theorem applied.** Since `j(T_n(z))` is
nonconstant, the generic rank18 subgroup stays independent at all but
finitely many rational addresses. This additional exceptional dependence
set has not been computed. Thus an arbitrary chart-valid rational `z` is
not automatically assigned rank18: use the existing exact
[halving-or-cycle classifier](DET1092_SPLIT_SEED_HALVING_DICHOTOMY_2026-09-08.md)
or a fresh independent rank certificate on each exported specialization.
No specialized point has been exported or submitted to V3 in this task.

**New obstruction, exact modular certificate.** These parameter maps are not
rational-function reparametrizations through the old orbit8044 conic.
For its exact branch quadratic `q(t)`, the polynomial

\[
D(z)^2q(N(z)/D(z))
\]

is squarefree of full degree40 or116 and coprime to `D`; the respective
certifying primes are127 and131. The normalized fibre products with that
conic therefore have genera19 and57. No rational map from `P1` can lift
either new parameter map through a curve of positive genus. This also
identifies an obstruction to combining these constructions while retaining
a rational parameter base. It is not a rank19 construction or a claim
about the individual rational points on those higher-genus fibre products.
The [certificate](../../artifacts/generated-results/elliptic-curves/det1092_pencil_multiples_v2/conic-lift-obstruction.json)
retains every attempted prime within the unchanged12-prime cap.

**Verified fixed-control comparison, now completed.** Initially, on302
and the eight old V3/V4-null addresses, the degree20 map had eight exact
rational-preimage exclusions and one UNKNOWN; the degree58 map had five
exclusions and four UNKNOWN. The initial degree58 test left302 UNKNOWN.
The [13 root-free certificates](../../artifacts/generated-results/elliptic-curves/det1092_pencil_multiples_v2/controls.json)
are replayed independently using `gcd(f,z^p-z)` instead of direct root
enumeration. A separate [exact follow-up](DET1092_GENERIC_TRANSLATION_ORBIT_OBSTRUCTION_2026-09-09.md)
now closes all five remaining cases using rational-root bounds, Hensel lifts
and exact two-dimensional lattice certificates, without enlarging the prime
pool. Both maps therefore miss all nine rational control parameters. These
are cover-incidence statements, not elliptic-rank upper bounds or an
explanation of302's later cascade. The initial bounded history is retained.

## Why finitely many covers cannot be the whole mechanism

**Established literature.**
[Pasten--Salgado, Theorem1.1](https://www.math.rug.nl/algebra/uploads/Main/PastenSalgado_2024.pdf)
proves non-thin rank jumps for a nonisotrivial elliptic K3 with reduced
fibres, a second distinct elliptic fibration over the same number field,
and Zariski-dense rational points.

**Verified application to determinant1092.** The parent has24 `I1` fibres
and nonconstant `j`. Its generic norm-eight pencil is a second elliptic
fibration over Q; `D.F=2` proves that it differs from the original.
Infinitely many multiples of an old nontorsion section give infinitely
many rational section curves, proving Zariski density. Therefore

\[
\{t\in\mathbf P^1(\mathbf Q): E_t\text{ smooth},
\operatorname{rank}E_t(\mathbf Q)\ge18\}
\quad\text{is non-thin.}
\]

**New implication for this programme.** Every fixed generically finite
cover of degree at least two has thin rational image. Finite unions stay
thin. The two covers constructed here, the old conic, and any finite list
of fixed genus-one carriers therefore cannot exhaust all seed-producing
parameters. The theorem proves more routes exist; it does not locate them,
give positive density, or distinguish a rank18 fibre from a rank31 fibre.

## What this explains in the completed experiments

**New synthesis of verified results.** The geometry separates several events
which can look alike to a point-search routine:

- An original generic section is a rational curve with `d=1`, `g=0`.
  Formula(1) gives zero perpendicular height over **full MW17**. This is
  consistent with the nine blinded reconstructions: each recovers a direction
  removed from M16, but none creates a direction beyond the full parent.
- The rational conic has `d=2`, `g=0` and perpendicular height3. Its generic
  seed is therefore genuinely new. The known dependent split is a
  specialization at which independence is lost, not a contradiction to the
  generic height theorem. The earlier exact cycle certificate identifies it.
- In the alternate pencil, `B` is inherited, but taking its group multiples
  is **not** the original elliptic group operation. The curves `2B` and `3B`
  cross the original base with degrees20 and58. Their positive perpendicular
  heights prove original independence after those base changes.
- None of these generic height calculations proves the existence of further
  residual directions on a particular specialization. The rank21 and302
  rank31 cascades remain additional arithmetic phenomena. Bounded rank18
  nulls are not full-rank upper bounds.

Thus the constructive mechanism is the change of elliptic fibration, with
ramification and intersection heights certifying when it really creates a
new seed. It is not rational-point existence alone, and not a new search score.

## Replay and retained limitations

**Verified resource record.** All completed commands use Sage10.9 with
25-second caps. Construction takes about1.1 and9.2 seconds; independent
replay about1.9 and15.0 seconds. The first expanded-coordinate attempt
timed out at25 seconds and is retained under `det1092_pencil_multiples_v1`;
version2 uses the same two multiples with factored maps. A checker polynomial
type error and correction are retained separately. No long/detached job,
point search, new control address, class group, subgroup census, or production
policy change was started.

From the repository root, run each command with `timeout 25s sage -python`:

```text
research/elliptic-curves/cas/verify_det1092_pencil_multiples.sage --multiple 2
research/elliptic-curves/cas/verify_det1092_pencil_multiples.sage --multiple 3
research/elliptic-curves/cas/verify_det1092_multisection_height_identity.sage
research/elliptic-curves/cas/check_det1092_pencil_multiple_controls.sage
research/elliptic-curves/cas/verify_det1092_pencil_control_certificates.sage
research/elliptic-curves/cas/check_det1092_pencil_conic_lifts.sage
```

**Unresolved, not conjecturally closed.** Neither construction explains why
the particular302 first seed exists, selects its calibrated member
prospectively, proves any full specialized rank, or predicts amplification.
The constructive advance is a generic two-fibration operation producing
new rational rank18 covers, together with its exact independence mechanism.
