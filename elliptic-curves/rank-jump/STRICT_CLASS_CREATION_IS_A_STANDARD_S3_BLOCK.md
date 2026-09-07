# Strict class creation requires the standard S3 block

The required extra strict classes occupy a specific Galois component.
Let K be the cubic two-division field, M its S3 Galois closure, and
D=Q(sqrt(disc K)) its quadratic resolvent. For S containing 2, infinity
and every bad elliptic prime, put

\[
 c_S=\dim_{\mathbb F_2}\mathrm{Cl}_S(K)/2,
 \qquad g_S=\dim_{\mathbb F_2}\mathrm{Cl}_S(D)/2.
\]

Writing V for the standard two-dimensional F2 representation of S3,
the ordinary S-class quotient of the sextic field has the exact form

\[
\boxed{\mathrm{Cl}_S(M)/2\simeq\mathbb F_2^{\,g_S}\oplus V^{\,c_S}.}
\tag{1}
\]

The first summand is trivial as an S3 representation. It contributes
**zero** elliptic strict classes. Class creation for the large jumps must
increase the required standard-module multiplicity beyond the inherited
pool, not merely supply more quadratic genus characters.

The equation-only audit computes g_S exactly on twelve retained complete
cases, including the supplemental +6 reference. It is zero on **all seven
completed high-jump panel rows**, the reference, and three of the four
completed lows. The remaining low has g_S=1. The cubic multiplicities
c_S are still UNKNOWN independently of exceptional points. Formula (1)
and these measurements locate the missing structure; they do not construct
it or explain a condition on t that forces its creation.

## Strict elliptic classes select the standard component

Let U be the strict 2-Selmer group used in the
[cubic S-class identification](STRICT_SELMER_AND_ARTIN_BLOCKS.md), so
U≅Hom(Cl_S(K),F2). At S the local class is zero, not merely soluble.
The two-division representation is V and becomes trivial over M.

For every subgroup H⊂S3,

\[
 H^1(H,V)=H^2(H,V)=0.
\tag{2}
\]

For C3 this is averaging in characteristic two. For S3, the normal C3
has no fixed vector on V, so the same assertion follows by
inflation–restriction (or its spectral sequence). For C2, 1+sigma has
rank one and its kernel equals its image; the cyclic cohomology formulas
give (2). The trivial subgroup is immediate. These exhaust the subgroups.

Global inflation–restriction therefore identifies

\[
 H^1(\mathbb Q,V)
 \simeq\operatorname{Hom}_{S_3}(G_M,V).
\]

Imposing the strict conditions gives

\[
\boxed{U\simeq
 \operatorname{Hom}_{S_3}(\mathrm{Cl}_S(M)/2,V).}
\tag{3}
\]

Indeed, a strict cocycle restricts to an everywhere-unramified character
over M, zero on decomposition groups above S. Conversely, a character
with those properties descends globally by (2). At a place in S its local
restriction over M is zero; the local inflation kernel is H1(H,V)=0 for
the corresponding decomposition subgroup H. Hence the descended local
class is zero. Outside S, E has good odd reduction, M/Q is unramified,
and the character descends to an unramified local class, which is in the
elliptic Kummer image. Ordinary class fields impose the real-place
splitting as well. This proves both directions of (3).

Thus no extra local descent kernel has been omitted at 2 or at infinity.
The [finite module certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_standard_s3_class_module_v1.json)
also checks normalized H1 and H2 cochain complexes on all six subgroups,
including the three transposition subgroups.

## Why the quadratic-resolvent part cannot contribute

Choose a three-cycle tau and work in F2[S3]. The two central orthogonal
idempotents are

\[
 e_0=1+\tau+\tau^2,\qquad e_1=\tau+\tau^2.
\]

For A=Cl_S(M)/2, extension and norm of ideals between M and D satisfy

\[
 \operatorname{Ext}\circ N=1+\tau+\tau^2=e_0,
 \qquad N\circ\operatorname{Ext}=3=1\quad\text{on mod-two quotients}.
\]

Consequently e0 A≅Cl_S(D)/2. The involution of D/Q acts by inversion
on ideal classes, hence trivially modulo two. Thus this entire summand
is a trivial S3 module, not just a C3-invariant vector space.

The other group-algebra block is

\[
 e_1\mathbb F_2[S_3]\simeq M_2(\mathbb F_2).
\]

Every module in this block is a sum of copies of V, and End_S3(V)=F2.
The finite certificate verifies the algebra isomorphism and its
multiplication, not just dimensions. Hence e1 A≅V^r and (3) gives
r=dim U=c_S. This proves (1), including the absence of a hidden
nonsemisimple C2 extension in its trivial component.

In particular

\[
 \dim\mathrm{Cl}_S(M)/2=2c_S+g_S.
\tag{4}
\]

There is no assumption on class-group exponents or on GRH. The quotient
is Cl/2; the two-torsion subgroup is not substituted for it as a Galois
module.

## Exact genus calculation from the equations

The [protocol](RESOLVENT_GENUS_COMPONENT_PROTOCOL.json) retains the
eleven complete rows of the sixteen-fibre panel and the existing
MW16-05 reference at 3/17. Five rows stay UNKNOWN because their factor
coverage is incomplete. There is no new integer factorization or
class-group computation on a research field.

The complete elliptic-discriminant factorization determines the quadratic
resolvent's fundamental discriminant d. Factor it into prime discriminants

\[
 d=d_1\cdots d_n,\qquad
 d_i\in\{-4,8,-8\}\ \text{or}\ d_i=(-1)^{(p-1)/2}p.
\]

The finite-unramified quadratic genus characters are products of the d_i,
modulo their total product d. They exhaust the quadratic unramified
characters; see Lemmermeyer,
[§1, especially Proposition 3](https://www.numdam.org/item/JTNB_1997__9_1_51_0.pdf).
For real D one must additionally impose positivity. The computation
imposes all infinite and finite S conditions explicitly.

For a rational product delta and a rational place v, it is locally square
over every factor of D⊗Q_v exactly when

\[
 [\delta]\in\langle[d]\rangle
 \quad\text{inside }\mathbb Q_v^*/\mathbb Q_v^{*2}.
\tag{5}
\]

If D splits locally, [d]=0 and (5) requires delta already square in Q_v.
If D is a quadratic field locally, the restriction kernel is precisely
the line generated by d. Applying (5) to the prime-discriminant products
gives a small F2 matrix. Its kernel, modulo the total-product relation,
is Hom(Cl_S(D),F2). No quadratic unit or regulator calculation is needed.

| Case | Ordinary quadratic genus dimension | S-split genus dimension g_S |
|---|---:|---:|
| 074d9 low, 2824/885 | 1 | 0 |
| 103b2 high, 3726/881 | 3 | 0 |
| 103b2 low, −1049/2296 | 4 | 0 |
| 11952 high, −2448/11 | 9 | 0 |
| 11952 low, −1171/1683 | 3 | 0 |
| 11952 high, 110314/102227 | 8 | 0 |
| 11952 low, 130349/28916 | 5 | **1** |
| 11952 high, 2828/2015 | 10 | 0 |
| Historic ICARM356 | 5 | 0 |
| Historic ICARM385 | 2 | 0 |
| Historic ICARM398 | 5 | 0 |
| Supplemental MW16-05, 3/17 | 3 | 0 |

For example, the ten ordinary genus characters on the 11952 high at
2828/2015 all disappear under the strict bad-place conditions. They do
not supply its required extra strict block. The low at 130349/28916 has
one surviving generator, represented by

\[
 \delta=1172785373546842861.
\]

Its full discriminant, local tests and nontrivial coefficient mask are
in the [arithmetic artifact](../../artifacts/generated-results/elliptic-curves/rank_jump_resolvent_genus_component_v1.json).
This is an explicit quadratic-resolvent class, not a cubic elliptic
strict class. After pullback to M it lies entirely in e0 A's character
space and contributes zero to (3).

The ordinary/S distinction is substantial even in small fields. Independent
certified class groups reproduce the genus calculation for d=−56,−84,12,136.
For d=−56, the ordinary class group is C4 and inverting the ramified primes
leaves S-class 2-dimension one. Thus the tempting shortcut “inverting all
ramified primes kills the quadratic genus quotient” is false in general.
The local matrix, including the dyadic condition, is necessary.

## What kind of new cover is actually required?

For each strict class, (3) gives an equivariant quotient A→V. Its class
field is a V4 extension of M. The normal closure over Q has group

\[
 V\rtimes S_3\simeq S_4,
\]

and over D it has group V⋊C3≅A4. The extension over M is unramified and
split at S. In particular it introduces no new ramified rational prime;
at S its completions over M are trivial.

For r independent strict classes the joint quotient is V^r: this follows
from the semisimple matrix block and End_S3(V)=F2. The corresponding
normal cover has groups

\[
 \operatorname{Gal}(N/M)\simeq V^r,\quad
 \operatorname{Gal}(N/D)\simeq V^r\rtimes C_3,\quad
 \operatorname{Gal}(N/\mathbb Q)\simeq V^r\rtimes S_3.
\tag{6}
\]

The splitting over Q also follows from H2(S3,V^r)=0. Its degree over M
is 4^r. This describes the common normal **incidence carrier** of an
independent block. It is not a minimal variety carrying rational points
on the elliptic 2-covers, nor evidence that several rational directions
arise from one additional irreducible module. One copy of V accounts for
one strict direction; the reason many copies exist remains open.

Equivalently, a strict class is a global lift of the given S3 representation
to S4 whose restrictions at S are conjugate to the original S3 complement.
The local class is a coboundary there. Class creation requires globally
distinct lifts with the same prescribed local behaviour, not merely a
change in the number of available local characters.

There is a stronger exclusion of the resolvent-tower shortcut. If T/D is
a finite Galois 2-extension and M/D is cyclic of degree three, then they
are linearly disjoint and

\[
 \operatorname{Gal}(TM/D)=\operatorname{Gal}(T/D)\times C_3.
\]

The three-cycle acts trivially on the pulled-back 2-group. Such extensions,
including genus, cyclic quartic and quaternion stages in a Galois 2-tower
of D, cannot supply a standard V quotient by base change alone. The
necessary A4 structure in (6) has a **nontrivial** C3 action. A construction
must address that interaction between degrees two and three.

## Relation to the forced large-jump block

Joining existing rank labels only after the equation-only calculation
gives the following required *new standard-module* dimensions beyond the
specialized global pool. These are existence necessities, not a newly
computed basis of c_S.

| Successful fibre | Required new standard copies ≥ | Corresponding standard vector dimensions ≥ |
|---|---:|---:|
| 103b2, 3726/881 | 7 | 14 |
| 11952, −2448/11 | 3 | 6 |
| 11952, 110314/102227 | 4 | 8 |
| 11952, 2828/2015 | 6 | 12 |
| ICARM356 | 9 | 18 |
| ICARM385 | 2 | 4 |
| ICARM398 | 3 | 6 |

The [comparison artifact](../../artifacts/generated-results/elliptic-curves/rank_jump_resolvent_genus_component_comparison_v1.json)
keeps all incomplete and censored control rows. It also records the
formula (4) and rank-derived lower bounds separately from the exact g_S.
No full cubic or sextic class-group dimension is inferred from g_S alone.

The next constructive target is an unramified, S-split V quotient over M
with nontrivial C3 action, independent of the inherited block. This is the
**incidence** gate. Ordinary genus counts, even when large, cannot certify
it; the entire pulled-back quadratic 2-tower is excluded as a direct source.
Once the extra block is constructed, CT and the remaining Sha obstruction
are the separate **solubility** gate. No visibility calculation occurs here.

The missing implication is still a condition on t that constructs enough
of the nontrivial quotients in (6), followed by a rational-solubility
argument. Relabelling c_S as a Galois multiplicity does not supply that
condition. The progress here is the exact decomposition, the completed
genus calculation, and exclusion of a whole class of ramification-based
constructions from the component that the jumps require.

## Reproducibility

The independent verifier uses Hilbert-dual constraints instead of the
worker's local squareclass coordinates and Sage kernels instead of its
bit elimination. It checks 1933 Hilbert bits, all complete prime products,
cubic irreducibility, and four small certified class groups. No research
field class group is run. The finite module calculation checks the two
group-algebra blocks and all six subgroup cochain complexes.

```sh
timeout 30 python3 elliptic-curves/rank-jump/standard_s3_class_module.py check
timeout 30 python3 elliptic-curves/rank-jump/resolvent_genus_component.py check
timeout 60 sage -python elliptic-curves/rank-jump/verify_resolvent_genus_component.py check
```

All workers complete within the frozen bounds. New files are rank-jump
specific; only the earlier class-creation note receives a navigation link.
Active searches, candidate populations and mathematical-status entries
remain outside this change.
