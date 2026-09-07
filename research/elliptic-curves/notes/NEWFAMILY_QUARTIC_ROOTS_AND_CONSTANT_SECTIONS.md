# New six-root quartic construction and constant-section classification

Status: exact symbolic identities, exact finite-reduction independence certificates,
and exact arithmetic calculations, with one final degree-10 rational-point closure
still marked for a dedicated standalone replay.  This note does **not** claim
literature novelty, a record curve, or an exact Mordell--Weil rank for the full
elliptic surfaces.  Rank statements below refer to the explicitly displayed
subgroups unless stated otherwise.

## 1. Construction

Let

\[
f(z)=\prod_{i=1}^6(z-r_i)
\]

be a monic sextic with rational roots.  Set

\[
P_T(x)=f(x-T)f(x+T).
\]

Let `S(x,T)` be the monic degree-six polynomial whose square agrees with
`P_T(x)` in degrees `x^12` through `x^6`, and write

\[
R(x,T)=S(x,T)^2-P_T(x).
\]

For a general monic sextic

\[
f(z)=z^6+a_1z^5+a_2z^4+a_3z^3+a_4z^2+a_5z+a_6,
\]

the condition that `R` have degree at most four in `x` is

\[
a_1^5-6a_1^3a_2+7a_1^2a_3+8a_1a_2^2-8a_1a_4
-12a_2a_3+24a_5=0.
\]

After translating the roots so that their sum is zero, this becomes the compact
symmetric condition

\[
2e_5=e_2e_3.
\]

The resulting quartic

\[
y^2=R(x,T)
\]

has twelve automatic rational sections

\[
x=r_i\pm T,\qquad y=S(r_i\pm T,T).
\]

The binary-quartic covariant map sends these sections to the Jacobian

\[
Y^2=X^3-27IX-27J.
\]

All subgroup-rank certificates discussed below use exact finite reduction,
not numerical height rank alone.

## 2. Integer root search and automatic rank

A C++ search over distinct integer roots in `[-50,50]`, primitive up to the
normalizations used by the search, tested `7,080,124` tuples and found
`20,488` quartic-condition hits, reducing to `17,001` unique root sets.

The automatic-section scanner found many generic subgroups of rank at least
9 and several with rank at least 10.  Three root sets reached the theoretical
maximum rank 11 available from the twelve automatic sections after choosing
one section as origin:

- `(-47,-43,-31,30,45,46)`;
- `(-47,-37,-18,27,29,46)`;
- `(-42,-35,-15,22,29,41)`.

For example, `(-47,-43,-31,30,45,46)` certifies rank at least 11 at
`T=11` with relation prime `5`, and again at `T=11` with relation primes
`7,11,13,17`.  A good specialization proving independence of the displayed
sections gives the corresponding generic lower bound for that displayed
subgroup.

Low-degree extra-section searches on the high-automatic-rank skeletons were
negative in the tested boxes:

- affine-linear `x=A*T+B`: about two million pairs per family in the large
  search box, no new sections for the rank-10/11 families tested;
- quadratic `x=A*T^2+B*T+C`: roughly 1.4--1.7 million triples per tested
  family after modular filtering, no exact new sections.

These are bounded experiments, not nonexistence theorems.

## 3. The section-rich seed

The asymmetric root set

\[
(-6,-5,1,2,3,5)
\]

produces

\[
\begin{aligned}
y^2={}&(T^2+50)x^4-90x^3-(2T^4+589)x^2\\
&-150(T^2-10)x+T^6-50T^4+733T^2-675.
\end{aligned}
\]

Besides the twelve automatic sections, exact square-polynomial search finds
seven affine-linear/constant sections:

\[
\begin{array}{c|c}
x(T)&y(T)\\ \hline
\frac75T+\frac{31}{5}&\frac{24}{25}T^3+\frac{434}{25}T^2+\frac{2811}{25}T+196\\
-\frac75T+\frac{31}{5}&\frac{24}{25}T^3-\frac{434}{25}T^2+\frac{2811}{25}T-196\\
T+\frac65&\frac{12}{5}T^2+\frac{36}{25}T+15\\
-T+\frac65&\frac{12}{5}T^2-\frac{36}{25}T+15\\
\frac35&T^3-\frac{634}{25}T\\
\frac15T-\frac{19}{5}&\frac{24}{25}T^3+\frac{38}{25}T^2-\frac{1011}{25}T+22\\
-\frac15T-\frac{19}{5}&\frac{24}{25}T^3-\frac{38}{25}T^2-\frac{1011}{25}T-22
\end{array}
\]

Thus there are nineteen explicit sections.

### Exact rank of the displayed 19-section subgroup

Let the covariant images of the nineteen sections be `P[0],...,P[18]`, with
`P[0]` chosen as origin, and put

\[
D_i=P[i+1]-P[0],\qquad 0\le i\le17.
\]

Exact finite reduction at good specializations certifies `D_0,...,D_8`
independent.  Exact function-field group-law verification over `Q(T)` proves
that all remaining displayed differences are combinations of these nine:

\[
\begin{aligned}
D_9={}&D_0+D_1+2D_2-2D_3-D_4-D_6-D_8,\\
D_{10}={}&2D_0+2D_1+D_2-D_3-2D_4-D_5-D_7,\\
D_{11}={}&-D_1+D_3+D_4,\\
D_{12}={}&D_0+D_1+2D_2-2D_3-2D_4,\\
D_{13}={}&-D_0-D_1+D_3+D_4+D_5+D_7,\\
D_{14}={}&-D_0-D_2+D_3+D_4+D_6+D_8,\\
D_{15}={}&D_2-D_3,\\
D_{16}={}&D_0+D_1-D_4-D_5+D_8,\\
D_{17}={}&D_0+D_1-D_4+D_6-D_7.
\end{aligned}
\]

Therefore

\[
\operatorname{rank}\langle P_0,\ldots,P_{18}\rangle=9
\]

for this explicitly generated subgroup over `Q(T)`.  This is **not** an upper
bound on the full Mordell--Weil group of the elliptic surface.

## 4. A structured root locus

Impose roots

\[
(-a,a,b,c,d,e),\qquad e=-b-c-d.
\]

The quartic condition factors as

\[
(b+c)(b+d)(c+d)
\left(a^2-b^2-bc-bd-c^2-cd-d^2\right)=0.
\]

The asymmetric component used below is

\[
\boxed{a^2=b^2+bc+bd+c^2+cd+d^2.}
\]

A search up to root bound 250 produced `6237` primitive distinct root sets on
this locus.  Under the small slope signature used in the experiment, the only
asymmetric section-rich equivalence class encountered was represented by
`(-6,-5,1,2,3,5)`; its negated root set is equivalent.

## 5. Constant-x section equations

Set

\[
C=c/b,\qquad D=d/b,\qquad U=u/b
\]

on the chart `b != 0`.  Substitution `x=u` into the normalized quartic gives
an even degree-six polynomial

\[
F(T,u)=C_6T^6+C_4T^4+C_2T^2+C_0.
\]

For this to be the square of an odd cubic `v_3 T^3+v_1 T`, the two eliminated
conditions are

\[
C_0=0,
\qquad
\Delta=C_4^2-4C_6C_2=0.
\]

The resultant in `U` factors into ten irreducibles.  The `S_4` action
permuting `b,c,d,e` groups them into exactly three geometric orbits:

- three factors of degree 2;
- four factors of degree 6;
- three factors of degree 10.

Thus only three component types need arithmetic analysis.

## 6. Degree-2 orbit: rational, but rank-dependent

A representative is

\[
CD+D^2-C+D=0.
\]

It is rational:

\[
C=\frac{D(D+1)}{1-D}.
\]

On this component the root-square condition is automatically a square, and a
convenient one-parameter form is

\[
\begin{aligned}
a&=s^2+1,\\
b&=s-1,\\
c&=-s(s+1),\\
d&=s(s-1),\\
e&=s+1,
\end{aligned}
\]

with forced constant section

\[
u=\frac{s(s^2-1)}{2(s^2+1)}.
\]

The seed `s=-1/3` gives the root set proportional to
`(-6,-5,1,2,3,5)` and `u=3/5` after clearing the common scale.

A combined automatic-plus-forced scan over 508 rational parameters stayed at
certified rank 9 in the fast common specialization.  Height relation recovery
at `s=2,T=17` produced the exact small relation

\[
2D_0-D_2+2D_3+2D_5-D_7-D_9-3D_{11}=0.
\]

The same integer relation verified exactly at nine unrelated specializations,
and was then proved identically over `Q(s,T)`.  In named automatic sections,

\[
\boxed{
2P_{-a,+}+2P_{b,-}+2P_{c,-}
=P_{a,+}+P_{d,-}+P_{e,-}+3Q.
}
\]

Hence the forced constant section `Q` can change saturation/index information
but cannot add a new rational rank direction.

The other two quadratic factors are in the same `S_4` orbit, so this closes
the degree-2 component type.

## 7. Degree-6 orbit: no admissible rational root configuration

For a representative, put

\[
S=C+D,\qquad P=CD.
\]

The component becomes quadratic in `P`:

\[
(S+10)^2P^2
-2(S+1)(17S^2+26S+18)P
+S^2(S+1)^2=0.
\]

Its discriminant is

\[
144(S+1)^4(8S^2+8S+9).
\]

After rationally parametrizing the conic and imposing that `C,D` themselves
be rational, the problem is birational to

\[
z^2=8m^4-72m^3+249m^2-360m+144.
\]

The associated Jacobian is

\[
y^2=x^3+52245x+7714278,
\]

with global minimal model

\[
y^2+xy+y=x^3-x^2+40x+155.
\]

Exact Sage calculations give conductor `126`, rank `0`, and torsion
`Z/6Z`.  The six rational points on the genus-one quartic are exactly

\[
(m,z)=(0,\pm12),(3,\pm3),(6,\pm42).
\]

The `m=0,6` points produce repeated roots; `m=3` does not lift to rational
`C,D`.  Therefore this component type has no admissible distinct rational
six-root configuration.

## 8. Degree-10 orbit

For the symmetric representative, the quotient by `C <-> D` again uses

\[
S=C+D,\qquad P=CD.
\]

The irreducible equation is cubic in `P` and can be written compactly as

\[
\begin{aligned}
0={}&(S^2+36S+36)^2P^3\\
&-9(S+1)(S^2+12S+12)(11S^2+12S+12)P^2\\
&+9S^2(S+1)^2(11S^2+8S+8)P\\
&-S^4(S+1)^3.
\end{aligned}
\]

The normalized quotient curve has genus 1.  A Riemann--Roch reconstruction
from the rational place `(S,P)=(0,1)` gives minimal elliptic curve

\[
\boxed{522i1:\quad y^2+xy+y=x^3-x^2+x+7,}
\]

with trivial torsion, rank `1`, and generator `(1,2)`.

Rational `C,D` require the quadratic lift

\[
W^2=S^2-4P.
\]

Function-field divisor arithmetic shows that this double cover has branch
degree 6 and therefore genus 4.  A search through `nG`, `-40 <= n <= 40`,
found only three quotient records satisfying this square condition; all map to
obviously degenerate repeated-root configurations.

The genus-4 curve has a stabilizer of order 8 under root permutations; the
stabilizer is dihedral `D_4`.  The central genus-2 quotient admits the even
model

\[
Y^2=20736z^6+340704z^4-428415z^2-19307236.
\]

Its Jacobian splits through the two elliptic quotients

\[
E_+=174c1:
\quad y^2+xy+y=x^3+x^2-5x-7,
\]

with trivial torsion and rank `0`, and

\[
E_-=522j1:
\quad y^2+xy+y=x^3-x^2-509x+4677,
\]

with trivial torsion and rank `1`, generated by `(11,-24)`.

The explicit `E_+` quotient has no finite rational points because
`E_+(Q)={O}`.  Tracing its two points at infinity back through the rational
parameterization gives four rational points on the genus-4 lift, all with
repeated zero roots.  This strongly closes the degree-10 orbit as having no
admissible distinct rational root configuration.

**Replay status:** the ingredients above were each checked exactly during the
interactive calculation, but the final two-points-at-infinity -> four
lift-points -> repeated-roots chain should be consolidated into one standalone
repository verifier before this last sentence is promoted to the same replay
status as the degree-2 identity and degree-6 rank-zero argument.

## 9. Current constant-section conclusion

Modulo the replay caveat in the last paragraph, the constant-`x` classification
on the structured root locus is:

| resultant orbit | arithmetic status | rank consequence |
|---|---|---|
| degree 2 | rational family | forced `Q` satisfies an exact `3Q` relation; no rank gain |
| degree 6 | finitely many rational points, all non-admissible | no distinct rational root configuration |
| degree 10 | genus-1 quotient, genus-4 lift; elliptic-factor analysis leaves only degenerate rational lifts | no admissible configuration found; final replay pending |

Thus the constant-section route does not currently provide a new
Mordell--Weil direction on this root locus.

## 10. Recommended next branch: general affine-linear sections

The natural next ansatz is

\[
x(T)=\alpha T+\beta.
\]

The section-rich seed already has six nonconstant affine-linear extras, so
this locus is known to be nonempty.  The next structural program is:

1. substitute `x=alpha*T+beta` in the general quartic on the structured root
   locus;
2. impose that the resulting degree-eight polynomial is a square of a quartic
   in `T`;
3. eliminate the ordinate coefficients recursively rather than by a giant
   Groebner basis;
4. factor the projected parameter locus;
5. quotient components by the `S_4` root symmetry;
6. on each rational component, test the new section against the automatic
   subgroup by exact function-field group law.

The objective is specifically to find a component where the affine-linear
section is not forced into the saturation of the automatic subgroup.

## 11. Ephemeral scripts and logs to preserve

Most of this branch was developed under `/tmp`.  Important scripts include:

- `newfamily_quartic_roots.cpp`
- `newfamily_rank_scanner.py`
- `newfamily_linear_generic.py`
- `newfamily_quadratic_generic.py`
- `newfamily_linear_signature.py`
- `newfamily_locus_search.cpp`
- `newfamily_rich_rank.py`
- `newfamily_exact_basis.py`
- `newfamily_relation_lattice_t17.py`
- `newfamily_function_field_probe.py`
- `newfamily_verify_relations_generic.py`
- `newfamily_constant_section_locus.sage.py`
- `newfamily_forced_section_rank_sweep.py`
- `newfamily_forced_combined_cert.py`
- `newfamily_forced_generic_relation.py`
- `newfamily_recover_forced_relation.py`
- `newfamily_prove_forced_relation_generic.py`
- `newfamily_sextic_factor3_genus1.py`
- `degree10_functionfield_genus.py`
- `degree10_lift_branch.py`
- `degree10_quotient_rr_weierstrass.py`
- `degree10_multiples_lift_search.py`
- `degree10_stabilizer.py`
- `degree10_central_genus2.py`
- `degree10_genus2_elliptic_factors.py`

The corresponding raw logs belong under
`artifacts/local/elliptic-curves/newfamily/` and should remain ignored.  Once
scripts are promoted into the repository, add compact deterministic replay
entry points under `elliptic-curves/scripts/` or `elliptic-curves/cas/` and
pin only concise results/manifests under
`artifacts/generated-results/elliptic-curves/`.

## 12. Claims deliberately not made

- No literature-novelty claim has been checked for this construction.
- No full Mordell--Weil rank equality is claimed for any surface here.
- The rank-11 automatic examples are lower bounds for explicitly generated
  subgroups.
- Bounded failed affine/quadratic searches are experiments, not
  nonexistence proofs.
- The final degree-10 rational-point closure should receive a dedicated replay
  verifier before being treated as a repository-level theorem/status entry.
