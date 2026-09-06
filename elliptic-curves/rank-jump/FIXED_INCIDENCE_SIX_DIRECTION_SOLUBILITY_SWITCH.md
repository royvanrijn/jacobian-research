# A six-direction solubility switch with fixed strict incidence

Transfer test: [fixing the cubic inside the actual families requires a
degree6 carrier of genus28 or31](FIXED_CUBIC_TRANSFER_REQUIRES_HIGH_GENUS.md).
All eight frozen matched pairs already have nonisomorphic cubic fields.
The twist mechanism cannot be transferred by a rational or elliptic base
change preserving that entire field.

For the MW16-05 reference at t=3/17 and its minus-one quadratic twist, the
equation-defined strict Selmer group is the same. Nevertheless, a
six-dimensional block of generic rational classes on the original becomes
a nondegenerate subgroup of **Sha[2] on the twist, unconditionally**.
Every one of its63 nonzero classes loses rational solubility.

Joining the completed class bound afterward gives, under its stated GRH
assumption, original rank **22**, twist rank **at most16**, shared strict
Selmer dimension **12**, and a drop in strict rational dimension from12
to at most6. Exact twist rank is UNKNOWN.

This is an explicit simultaneous-solubility mechanism with fixed strict
incidence. It is a **quadratic-twist control using a generic block**, not
an explanation or predictor of the additional +10/+11 directions in the
original K3 specialization families.

## Equation, class transport and input boundary

The original curve has the minimal model

```text
y^2 + x*y = x^3 - 182451976602578656424609725499140*x
                  + 710003150253794219215652666162794189038512805392.
```

The [masked reference](CONDITIONALLY_EXACT_SOLUBLE_BLOCK_REFERENCE.md)
transports its equation and sixteen generic sections to y²=f(x),
f(x)=x³+Ax+B. The twist is y²=g(x), g(x)=x³+Ax−B. Exactly

\[
 g(z)=-f(-z),\quad\operatorname{disc}(g)=\operatorname{disc}(f),\quad
 \theta_+=-\theta_-.
\]

The [new calculation](../../artifacts/generated-results/elliptic-curves/rank_jump_fixed_cubic_minus_reference_v1.json)
replays both equations at every discriminant prime. Their complete sets
S of bad primes together with2 and infinity agree. Their finite part is

```text
2, 3, 5, 13, 17, 19, 29, 71,
2465779087453622652131442949,
519784438179112504122441050306814600881.
```

Thus their strict class groups are canonically identified:

\[
 U_+=U_-=U=\operatorname{Hom}(\operatorname{Cl}(\mathcal O_{K,S_K}),\mathbf F_2).
\]

This equality is unconditional, even before its dimension is known. The
workers use only the masked equation, generic sections and their previously
verified six strict classes/half ideals/Artin data. Every transported strict
representative is locally square above S, positive at all real places, and
has square norm. Unramifiedness outside S is preserved by the explicit field
isomorphism. No exceptional class or coordinate enters these computations.

The [conditional comparison](../../artifacts/generated-results/elliptic-curves/rank_jump_fixed_cubic_minus_reference_comparison_v1.json)
separately joins the existing class/rank completion theorem. Its proof uses
exceptional points and GRH for the specified quadratic ordinary class
characters; it is an outcome label here, not an independent feature.

## The precise event: four-torsion changes while strict two-classes remain

Let V=G∩U be the six-dimensional strict part of the original generic
Kummer subgroup. Its original CT form is zero because its classes are
rational. The [scalar-cup comparison](INDEPENDENT_SCALAR_CUP_AND_TWIST_BLOCKS.md)
identifies the minus-one change of four-torsion action with the scalar
quadratic character chi_(−1) times the identity. On strict classes the
local correction terms vanish, giving

\[
 \mathrm{CT}_{E^-}|_V=A+A^\mathsf T
 =\begin{pmatrix}
0&1&1&0&0&1\\
1&0&0&0&0&1\\
1&0&0&1&0&0\\
0&0&1&0&1&1\\
0&0&0&1&0&0\\
1&1&0&1&0&0
\end{pmatrix}
 \simeq H^{\oplus3}.
\]

Here H is the nondegenerate alternating plane. In the frozen strict basis
e0,...,e5, the verified symplectic pairs are represented by bit masks
(1,2), (6,8), (22,45). The six basis-change vectors are independent.

Rational Kummer classes annihilate the entire Selmer group under CT. Since
the displayed restriction is nondegenerate, V intersects the twist's
rational Kummer image trivially. The map

\[
 V\longrightarrow\Sha(E^-/\mathbb Q)[2]
\]

is consequently injective. This proves **dim Sha(E^-)[2]≥6 without GRH**
and without assuming that the whole Sha group is finite. No nonzero class
of this block is divisible by2 in Sha: any such class would pair trivially
with all2-torsion, contradicting the nondegenerate restriction.

The controlled structural chain is therefore

\[
 \boxed{\text{minus-one change in four-torsion action}
 \Longrightarrow\text{rank-six CT obstruction on unchanged strict classes}
 \Longrightarrow\text{six independent rational-to-Sha switches}.}
\]

For the usual twist parameter d, this also gives exact conditions on two
squareclasses: V is rational on twists with d in Q*², and no nonzero V-class
is rational on twists with d in −Q*², using their standard identifications.
Other squareclasses of d are not classified. This d is **not** the original
family parameter t; no condition on that t is obtained.

## Whole-curve bounds and what is held fixed

Both equations have local point-product dimension ell=11. The derivative
class −disc(f)f′(theta), respectively −disc(g)g′(theta), has norm disc(f)^4
and signs(−,+,−). Its real pairing with the nontrivial point-image class is
nonzero. Global reciprocity gives boundary dimension at most10 for each.

Write c_S=dim U. Without assuming its exact value, the twist satisfies

\[
 \dim\operatorname{Sel}_2(E^-)\le c_S+10,\qquad
 \operatorname{rank}E^-(\mathbb Q)\le c_S+4,\qquad
 \dim(U\cap\delta E^-(\mathbb Q))\le c_S-6.
\]

The completed reference theorem gives c_S=12 under its stated GRH
assumption. Hence:

| Quantity | Original | Minus-one twist |
|---|---:|---:|
| Strict Selmer dimension, conditional | 12 | 12 |
| Full Selmer dimension, conditional | 22 | between12 and22 |
| Rational strict dimension, conditional | 12 | at most6 |
| Mordell–Weil rank, conditional | 22 | at most16 |
| Rational dimension inside marked V, unconditional | 6 | 0 |
| Sha[2] dimension | 0 under GRH | at least6 unconditionally |

The **full** Selmer groups need not be identical; the strict group is the
fixed incidence object. The marked V is generic on the original family,
and the twist does not preserve that family's rational generic subgroup.
The rank difference here is therefore not a specialization quotient over
the same marked generic group.

## Explicit carriers and a minimality result

The [carrier certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_minus_reference_carriers_v1.json)
gives six pairs of symmetric4-by-4 quadric matrices. For each transported
strict representative beta_i, write

\[
 \beta_i(u_0+u_1\theta_-+u_2\theta_-^2)^2
   =q_{i,0}(u)+q_{i,1}(u)\theta_-+q_{i,2}(u)\theta_-^2.
\]

Its2-cover C_i is the complete intersection in P³

\[
 q_{i,2}(u)=0,\qquad q_{i,1}(u)+v^2=0.
\]

On v≠0, its map to E^- has
x=q_(i,0)/v² and y=sqrt(N beta_i)N(u)/v³. The norm identity verifies the
curve equation; this is the usual degree-four2-cover. Each determinant
pencil has three distinct finite roots and a simple root at infinity.
Thus all six intersections are smooth, genus1 and projective degree4.
Their Jacobian is E^- through the2-cover construction; their computed
pencil j-invariants also agree. Each has points over every completion of Q,
but none has a rational point, by the CT block above.

There is a precise restricted answer to the minimal-carrier question here.
Let a carrier have Q-morphisms to all six C_i.

* Among carriers which are torsors under abelian varieties over Q, the
  **minimal dimension is exactly6**. The product of the C_i achieves it.
* Any smooth projective geometrically integral **curve** carrier has
  **genus at least6**. Existence of a genus-exactly6 carrier is not asserted.

Here is the argument, which applies to r independent order-two torsor
classes with a common elliptic Jacobian. For a curve Y, choose a geometric
point y0. Every map Y→C_i induces a Q-homomorphism
phi_i:Jac(Y)→E^-; its torsor class is the pushforward by phi_i of the
cocycle [sigma(y0)−y0]. This follows by taking differences of the images
of y0 and sigma(y0). The same argument applies to a torsor under an abelian
variety A, replacing Jac(Y) with A. The Jacobian universal property and
finite freeness of Hom are standard; see
[Milne, III§6 and I Theorem10.15](https://www.jmilne.org/math/CourseNotes/AV.pdf).

In characteristic zero over Q, differentiation embeds
Hom_Q(A,E^-) tensor Q into Hom_Q(Lie A,Lie E^-), a Q-vector space of
dimension dim A: a homomorphism with zero differential is zero.
Consequently Hom_Q(A,E^-) has rank at most dim A. Its pushforward of one
torsor class generates an abelian group with at most that many generators,
so its order-two subgroup has F2-dimension at most dim A. Six independent
classes [C_i] therefore force dim A≥6, or genus(Y)≥6. The product carrier
under (E^-)^6 attains the abelian-torsor bound.

This minimality concerns the **obstructed fibre** and morphisms to all
six fixed torsors. It is not an absolute lower bound for every way to encode
a rationality test. On the original soluble fibre the torsor classes vanish,
so this lower-bound argument does not apply. However, a family of abelian
torsor carriers with regular maps to these covers and good specialization
at the minus-one fibre cannot have dimension below6. A smooth curve family
with such specialization cannot have genus below6. A proposed smaller
uniform carrier would have to degenerate or lose the relevant maps there.

One must also distinguish a product over Q from a fibre product over E.
The latter demands a common image point, hence a common Kummer class, and
cannot test simultaneous solubility of distinct classes. The product over
Q asks for separate rational lifts and is the relevant carrier here.

## Lessons for the rank-jump programme

1. **Solubility mechanism established in this control:** a change of
   four-torsion action produces a rank-six obstruction on a fixed strict
   incidence space. This is an actual block event, not point visibility.
2. **Carrier restriction:** six independent obstructed torsors with one
   elliptic Jacobian cannot all be reached from an auxiliary elliptic curve
   or an abelian surface over Q. A low-genus uniform explanation must account
   for its failure or degeneration on the obstructed fibre.
3. **Still missing for the fresh jumps:** an independently constructed
   additional class block modulo G, and a comparable variation of its
   obstruction with the original family parameter t. Here the tested block
   belongs to G, and twisting changes the rational generic subgroup.
4. **No selector follows:** the −1 switch on inherited blocks already occurs
   on observed-low fibres. Its presence is not a high-jump discriminator.
   No visibility or candidate-scoring feature is proposed for Agent1.

The next transfer test should ask whether an analogous fixed additional
class block can be defined over an auxiliary parameter cover of the original
family and whether its four-torsion obstruction changes on that cover.
Another computation of a generic switch alone would not answer that question.

## Verification

The [verification certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_fixed_cubic_minus_reference_verification_v1.json)
replays all upstream strict half ideals and repaired Artin entries, checks
the two complete local datasets using PARI versus Sage, verifies the real
reciprocity constraints,192 quadric coefficients by rational cubic algebra,
30 independent pencil determinants, all six smoothness tests and the
symplectic basis change. All pass within the declared limits. Carrier
minimality is the mathematical argument above; the computation verifies its
independent-torsor hypothesis, not the universal theorem by enumeration.

```sh
timeout 60 sage -python elliptic-curves/rank-jump/verify_fixed_cubic_minus_reference.py check
```

The [protocol](FIXED_CUBIC_MINUS_REFERENCE_PROTOCOL.json) permits no point
search, class-group computation or new parameter sweep. Source computations,
active search protocols and mathematical-status entries are left unchanged.
