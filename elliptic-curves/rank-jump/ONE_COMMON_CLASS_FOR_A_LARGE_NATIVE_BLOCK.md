# A large native block would require one class soluble on many twists

The published R17 global cover pool has dimension **17 or 18**, so at
most **one** class remains beyond the marked generic subgroup. Two exact
odd Frobenius counts improve the previous upper bound of 19.

This sharpens the entire retained native-cover experiment:

| Native supports | Full fibre-product genus | New generic rank, at most |
|---:|---:|---:|
| 1 | 0 | 2 |
| 2 | 1 | 5 |
| 3 | 5 | 10 |
| 4 | 17 | 19 |

There is a stronger structural consequence. **A generic gain of at least
14 on any four of these supports requires one common extra global Kummer
class to be rational on at least ten of the fifteen character twists.**
A three-support gain of ten requires one common class rational on all
seven twists. These are necessary conditions, not observations that the
required class or points exist.

For an explicitly constructed class, its twist covers have a common
conic quotient. Once that conic is parametrized, simultaneous solubility
means that a single binary quartic represents the requisite distinct
twist squareclasses. The conic alone does not solve the covers.

This is a restricted carrier mechanism with a precise missing step. It
does not yet explain the class creation of the fresh +10/+11 or historic
+12/+14 fibres. The 37 supports are the same retrospective native
selection as before; no additional support, parameter or exceptional
point was selected.

## One global incidence dimension remains unresolved

Write F=Q(t), G for the 17-dimensional generic Kummer subgroup, and L
for the original global cover pool. The
[root-curve theorem](ROOT_CURVE_TORSION_AND_REAL_CAPACITY.md) identifies
L with the rational two-torsion of the genus-ten cubic root curve's
Jacobian. The preceding native calculation directly established dim L<=19
on the published model and dim G=17.

The [frozen Frobenius protocol](NATIVE_ROOT_FROBENIUS_SLACK_PROTOCOL.json)
counts the root curve at its first twelve good primes at most 1009:

| p | #C(Fp) | p | #C(Fp) |
|---:|---:|---:|---:|
| 131 | 112 | 191 | 186 |
| 137 | 146 | 197 | 186 |
| 151 | 158 | 199 | 218 |
| 157 | 156 | **211** | **205** |
| 167 | 150 | 227 | 222 |
| 181 | 156 | **229** | **225** |

On J[2], Frobenius has reciprocal characteristic polynomial of degree
20 and fixes G pointwise. Its polynomial is divisible by (T+1)^17.
The reciprocal residual cubic is either (T+1)^3 or T^3+1. An odd point
count makes the Frobenius trace odd, forcing the second alternative:

\[
 \operatorname{charpol}(\mathrm{Frob}_p\mid J[2])
 =(T+1)^{18}(T^2+T+1).
\]

The fixed space has dimension at most the algebraic multiplicity 18.
Thus

\[
 \boxed{17\le\dim L\le18,\qquad \dim(L/G)\le1.}
\tag{1}
\]

The [count certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_native_root_frobenius_slack_v1.json)
and [independent replay](../../artifacts/generated-results/elliptic-curves/rank_jump_native_root_frobenius_slack_verification_v1.json)
check all 2190 finite and infinity fibres. The worker enumerates roots;
the replay uses finite polynomial factorization and independently selects
the good primes. It explicitly enumerates the two reciprocal polynomial
possibilities. Equality dim L=18 is **UNKNOWN**. Even counts do not prove
that Frobenius acts trivially, and no fixed-space dimension is inferred
by equating it with algebraic multiplicity.

These counts concern the fixed published root curve. They do not assert
the sharper bound for every other fibration chart of the same K3.

## The native branch kernels are exactly known

For a native quadratic q_i, let Z_i be the kernel of specialization of
L at its quadratic branch place. The earlier
[branch theorem](BRANCH_FIBRES_BOUND_NATIVE_CLASS_CREATION.md) proves
that every rational Kummer class on a twist with that branch support
lies in Z_i. It measured rank 16 on G, giving dim(G intersect Z_i)<=1.

The [new half-section certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_native_branch_half_section_v1.json)
identifies the kernel exactly. If the native section on u^2=q_i(t) is
P_i=(x0+x1*u,y0+y1*u), its branch specialization is (x0,y0). Direct
arithmetic in Q[b]/(q_i) verifies

\[
 2P_i(b)=w_i(b),
\]

where w_i is the supplied integer word in the 17 generic sections.
Every word modulo two equals its previously verified one-dimensional
fingerprint kernel and is nonzero. Consequently

\[
 G\cap Z_i=\langle w_i\rangle,
 \qquad \dim Z_i\le18-16=2.
\tag{2}
\]

All 37 identities are independently replayed with rational arithmetic
modulo a quadratic and explicit chord/tangent formulas, rather than
Sage's number-field elliptic-curve implementation. The
[verification](../../artifacts/generated-results/elliptic-curves/rank_jump_native_branch_half_section_verification_v1.json)
checks 310 group operations. The points used here are generic native
section terms from the earlier equation-only atlas, not exceptional
points on a research fibre.

The 37 nonzero words are pairwise distinct, and none of their **7770
triples** sums to zero modulo two. This is an exact constraint on these
marked generic directions; it is not a rank-prediction statistic.

## All character capacities improve together

Let d_I=c_I*product(q_i : i in I), with c_I a nonzero rational scalar.
For a common multiquadratic cover, c_I is the product of its chosen
singleton scalars; the individual bounds also hold for arbitrary c_I.
Let T_I be the rational point Kummer image of E^(d_I)(F). The branch
cubics are irreducible, so the local vanishing theorem gives

\[
 T_I\subseteq\bigcap_{i\in I}Z_i.
\]

For a singleton the right side has dimension at most two. For two or
more distinct supports, the joint generic branch map has rank 17,
because its individual kernel lines are distinct. Its full kernel in
L has dimension at most one. Therefore

\[
 \dim T_{\{i\}}\le2,\qquad
 \dim T_I\le1\quad(|I|\ge2).
\tag{3}
\]

All these curves have no rational two-torsion over F, so their point
Kummer dimensions equal their ranks. Character decomposition over the
full multiquadratic cover gives

\[
 \boxed{J_{\rm generic}\le2k+(2^k-1-k)=k+2^k-1.}
\tag{4}
\]

This includes all sections, and integral gluing does not change the
rational character ranks. It is uniform in the rational scalar choices.
The genus is 1+2^(k-1)(k-2), since all branch supports are disjoint.
Thus this full native fibre-product construction needs at least three
supports for a generic gain of 6 through 10, and at least four for a
generic gain of 11 through 19. It is not a lower genus bound for every
possible auxiliary carrier.

There is an informative dichotomy. If L=G, singleton ranks are at most
one and every multiple-support character rank is zero. For the original
native scalar choices, the known native sections then give **exactly k**
new generic directions on every k-support fibre product. Hidden native
character blocks require L/G to be nonzero.

## Why +14 on four supports forces a common class

Suppose L/G has dimension one. The part of Z_i outside G is either
empty or an affine line with precisely two points:

\[
 A_i=Z_i\setminus G=\{a_i,a_i+w_i\}.
\]

A nonzero multiple-support class must be a common point of all the A_i
in its support. Two different A_i cannot share two points, because their
directions w_i differ. Represent each nonempty A_i by an edge joining
its two points. A triangle of three edges would give w_i+w_j+w_k=0,
which the generic-word certificate excludes.

For four supports this is a simple triangle-free graph with at most four
edges. A vertex shared by r edges can support at most 2^r-1-r
multiple-support characters. If all four edges do not meet at one
vertex, their total capacity is at most five. To see the bound: with
maximum degree two it is at most four; a degree-three vertex contributes
four, and the fourth edge can add at most one without forming a triangle.

The four singleton characters contribute at most eight. Hence

\[
 \boxed{J_{\rm generic}\ge14
 \ \Longrightarrow\ A_1\cap A_2\cap A_3\cap A_4=\{h\},
 \quad h\notin G.}
\tag{5}
\]

Every multiple-support point class is then either zero or h. A singleton
point space has dimension at most 1 plus the indicator that h is rational
on that twist. It follows that

\[
 J_{\rm generic}\le4+
 \#\{\varnothing\ne I\subseteq\{1,2,3,4\}:h\in T_I\}.
\]

Thus +14 requires h to be rational on at least **ten** distinct twists.
This is valid even for the arbitrary scalar choices; the native scalar
choices additionally supply each w_i as a rational singleton class.

For three supports, attaining the bound ten in (4) forces every one of
the seven character spaces to attain its bound. A nonzero full-product
class h lies in every Z_i. It then spans every multiple-support space
and lies in both-dimensional singleton spaces. Therefore h is rational
on **all seven twists**.

The [common-class capacity certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_native_common_class_capacity_v1.json)
checks all trace triples and all 22702 triangle-free labelled graphs on
at most eight vertices with at most four edges. The graph calculation
verifies the bound five independently of unknown arithmetic values of
the a_i. It does not construct those values.

## The common solubility carrier has explicit equations

An extra class must first be constructed as
h=h0+h1*theta+h2*theta^2 in K*=F[theta]*, where
theta^3+A*theta+B=0 and N(h)=n^2. For z=z0+z1*theta+z2*theta^2 write

\[
 hz^2=Q_0(z)+Q_1(z)\theta+Q_2(z)\theta^2.
\]

The Q_i are ternary quadratic forms over F. The genus-one 2-cover
representing h on the d-twist is

\[
 \boxed{Q_2(z)=0,\qquad Q_1(z)+d w^2=0\quad\text{in }\mathbb P^3.}
\tag{6}
\]

For w nonzero its point map is

\[
 x=Q_0(z)/w^2,\qquad y=nN(z)/w^3,
 \qquad y^2=x^3+d^2Ax+d^3B.
\]

Conversely a rational point with Kummer class h supplies z with
x-d*theta=h z^2. For nontrivial h there is no rational point of (6)
with w=0: such a point would make h z^2 a nonzero scalar of square norm,
forcing both that scalar and h to be squares.

All d share the **same conic** Q_2=0. Its symmetric matrix has determinant
-N(h). The determinant of the pencil for (6) is

\[
 -dN(h)\,\mu(\lambda^3+A\lambda\mu^2-B\mu^3).
\]

For d*N(h)*(-4A^3-27B^2) nonzero the pencil is separable and the
intersection is a smooth genus-one curve. The
[universal equations](../../artifacts/generated-results/elliptic-curves/rank_jump_common_class_conic_v1.json)
are checked by exact polynomial identities and an
[independent companion-matrix derivation](../../artifacts/generated-results/elliptic-curves/rank_jump_common_class_conic_verification_v1.json).

If the conic has no F-point, none of these twists can represent h.
If it has an F-point, parametrize it by quadratic forms z_i(s0,s1).
Then -Q_1(z(s0,s1)) is a single binary quartic H_h, and (6) becomes

\[
 d w^2=H_h(s0,s1).
\]

The required event is that this one quartic represent at least ten of
the fifteen prescribed squareclasses d_I, using rational functions over
F. These are separate points on a common conic. One conic point cannot
serve two distinct d_I: their ratio would be square. Conic solubility
alone is therefore not the missing simultaneous-solubility theorem.

## Where the extra incidence dimensions would actually enter

This mechanism would not place a fourteen-dimensional jump inside the
one-dimensional quotient L/G. Let M be the full free Mordell-Weil group
over the multiquadratic function field, and let Lambda be the sum of
its invariant and quadratic-character subgroups. If its generic gain
J exceeds k, some character point class lies outside G; all character
classes still lie in L. Restriction of H1(F,E[2]) is injective under this
2-power extension, because the irreducible cubic gives no rational
two-torsion there. Thus the Kummer image of Lambda has dimension exactly
18, while rank M=17+J.

Elementary lattice algebra gives

\[
 \boxed{\dim_{\mathbb F_2}(M/\Lambda)[2]
 =\dim(\Lambda\cap2M)/2\Lambda=J-1.}
\tag{7}
\]

A +14 block would require thirteen independent halving relations between
character subgroups, with [M:Lambda] divisible by 2^13. This is a
necessary integral structure, not thirteen additional free generators
created by saturation. The new mod-two classes live over the parameter
cover, outside the restriction of L. Their possible specialization
therefore does not contradict the earlier capacity theorem for L.

The known native points alone cannot construct the missing class merely
by saturation. Their 17+k classes are already independent modulo two:
sigma_i(P_i)=w_i-P_i, while sigma_i fixes every other native P_j. Applying
sigma_i-1 to a proposed mod-two relation kills its P_i coefficient,
since w_i is nonzero and restriction on G is injective. All coefficients
then vanish. This uses the retained global native trace identities.

Equation (7) supplies no visibility score and does not claim to explain
the actual 43-chart recovery data. That connection would require an
explicit block and a separate chart-to-lattice comparison.

## The condition on t0 and the remaining implications

Once the required rational twist sections have been constructed, their
simultaneous rational specialization has the explicit carrier condition

\[
 q_i(t_0)c_i\in\mathbb Q^{\times2}\quad\text{for every chosen }i,
 \qquad\Delta_E(t_0)\ne0.
\tag{8}
\]

For four supports this is a rational point on the genus-seventeen fibre
product, away from its branch divisor. It produces rational evaluations;
independence modulo the original generic subgroup still needs proof.
No such section block or successful specialization is asserted here.

The next decisive incidence calculation is **whether L/G is zero or
one-dimensional**, followed, in the latter case, by an explicit h.
That would turn (6) into actual coefficient-defined equations. The
strongest solubility test would then ask whether the shared quartic
represents the required character values. These two gates precede any
new parameter search.

The ranked outcome is: a single common-class mechanism is now necessary
for a large four-native-cover block; geometric capacity and isolated
native lifts are insufficient; the extra class and its simultaneous
solubility remain missing. This concerns the retained published-R17
supports, not a replacement measurement of the fresh panel's class
groups or CT matrices. No candidate scores or active search files change.

```sh
timeout 30 sage -python elliptic-curves/rank-jump/verify_native_root_frobenius_slack.py check
timeout 30 python3 elliptic-curves/rank-jump/verify_native_branch_half_section.py check
timeout 30 python3 elliptic-curves/rank-jump/native_common_class_capacity.py check
timeout 30 sage -python elliptic-curves/rank-jump/verify_common_class_conic.py check
```
