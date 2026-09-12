# The fixed-quintic target stabilizer: all-degree stable marked rigidity

This note continues the rank-five Tschirnhaus transition theorem by
determining the formerly open fixed-map orbit

\[
 \operatorname{Stab}_{\mathrm{st}}^t(F_R)\cdot y_R.
\]

It proves that the standard marked orbit is a point in every degree and
after any number of identity stabilizations.  The decisive all-degree
argument combines the exposed Newton intruder `P^2*B^5*C` with Kuroda's
stable-invariant theorem.  The earlier exact recursive Newton-face
certificate through target degree twenty-eight is retained as an independent
calculation and as information about the full stable group away from the
standard zero section.

The reusable descent, boundary-faithfulness, and deck-group argument is now
stated once in the
[stable intruder descent criterion](STABLE_INTRUDER_DESCENT_CRITERION.md);
this note retains the quintic-specific boundary and topology calculations.

Work over an algebraically closed characteristic-zero field.  Use the fixed
presentation

\[
 R(T)=\prod_{i=1}^5(T-i)
\]

and the compiler target coordinates `(P,B,C)`, where `B=-2b` relative to
the `(pi,b,c)` convention.  The selected target is

\[
 y_R=
 \left(
 \frac{85}{274},
 \frac{225}{137},
 \frac{120}{137}
 \right).                                               \tag{1.1}
\]

The compiler seeds are

\[
 u_4=-\frac{61712472}{10440125},\qquad
 u_5=\frac{5636405776}{4437053125}.                    \tag{1.2}
\]

## 1. The intrinsic torus action is trivial

For a stable self-equivalence of a fixed rank-five coefficient point, the
stable-moduli theorem gives

\[
 \beta=\alpha^{-2},\qquad
 \alpha^5=1,\qquad
 \alpha^6=1.                                           \tag{1.3}
\]

Since `gcd(5,6)=1`,

\[
 \boxed{\alpha=\beta=1.}                              \tag{1.4}
\]

Thus every stable target self-equivalence induces the identity on the
intrinsic `(P,r)` normalization of the ramified boundary.  The remaining
kernel is invisible to the coefficient-torus quotient.

## 2. The fixed quintic discriminant

The normalized inverse equation is

\[
 E(S)=
 u_5P^5S^5+u_4P^4S^4+PS^3-\frac B2S^2+S-\frac C2.
                                                               \tag{2.1}
\]

Exact elimination gives

\[
 \operatorname{disc}_S(E)=P^8H(P,B,C),                \tag{2.2}
\]

up to a nonzero rational scalar, where `H` is irreducible,

\[
 \deg H=16,\qquad \#\operatorname{supp}H=59,           \tag{2.3}
\]

and its unique degree-sixteen monomial is a nonzero multiple of

\[
 P^{12}C^4.                                            \tag{2.4}
\]

Moreover `H(y_R)` is nonzero.  The selected complete fibre is therefore
away from the ramified discriminant.

Let `T` be an unstabilized target component of a self-equivalence.  By
(1.4), its restriction to the normalization of `H=0` is the identity.
The normalization map is dominant onto that prime hypersurface, so `T`
fixes `H=0` pointwise.  Consequently

\[
 \boxed{
 T=\operatorname{id}+HV
 }                                                       \tag{2.5}
\]

for a polynomial vector `V`.

The same conclusion applies to the first three components of a stabilized
target automorphism on `A^(3+s)`: their differences from `(P,B,C)` are
divisible by `H` in `k[P,B,C,t_1,\ldots,t_s]`.

## 3. The logarithmic kernel

Since a target automorphism preserves the prime ideal `(H)`, write

\[
 H\circ T=\rho H,\qquad \rho\in k^\times.              \tag{3.1}
\]

Reducing (3.1) modulo `H^2` after (2.5) gives the necessary equation

\[
 V(H)-QH-\kappa=0,\qquad \kappa=\rho-1.                \tag{3.2}
\]

In fact `kappa=0` in every degree.  Set `P=1` and choose a root `s` of

\[
 15u_5s^4+8u_4s^3+3s^2-1=0.
\]

Then set

\[
 B=20u_5s^3+12u_4s^2+6s
\]

and choose `C` so that `E(s)=0`.  The resulting quintic has a triple root
at `s`.  Its coefficient point is therefore a singular point of the
discriminant `H=0`.  Evaluating (3.2) there gives `kappa=0`.

For `deg(V)<=m`, exact characteristic-zero coefficient matrices have
nullities

\[
\begin{array}{c|rrrrrrrrrrrrr}
m&0&1&2&3&4&5&6&7&8&9&10&11&12\\ \hline
\dim\mathcal L_m&0&0&0&0&0&0&0&2&21&54&104&174&267.
\end{array}                                             \tag{3.3}
\]

The two degree-seven solutions have independent evaluations at `y_R`:

\[
 \operatorname{rank}
 \left(V_1(y_R),V_2(y_R)\right)=2.                    \tag{3.4}
\]

Thus degree seven is a genuine possible marked-motion frontier, not merely
a pair of logarithmic fields vanishing at the selected target.

There is also an all-degree structural calculation.  Homogenize `H` to
degree sixteen with a variable `Z`, put

\[
 R=\mathbb Q[P,B,C,Z],
\]

and let `M` be the kernel of the homogeneous map

\[
 R(-15)^3\oplus R(-16)
 \xrightarrow{(H_P,H_B,H_C,-H)}
 R.                                                     \tag{3.5}
\]

An exact Macaulay2 resolution is

\[
 0\longrightarrow R(-25)^6
 \longrightarrow R(-24)^{18}
 \longrightarrow R(-22)^2\oplus R(-23)^{13}
 \longrightarrow M\longrightarrow0.                   \tag{3.6}
\]

After removing the common derivative shift by fifteen, the filtered
Hilbert series is

\[
 \boxed{
 \operatorname{Hilb}_{\mathrm{fil}}(M;t)
 =
 \frac{2t^7+13t^8-18t^9+6t^{10}}{(1-t)^4}.
 }                                                       \tag{3.7}
\]

In particular, the entire logarithmic module is generated by two fields in
degree seven and thirteen new fields in degree eight; higher degree does not
introduce a new generator type.  Formula (3.7) gives the nullities in (3.3).
It is an all-degree statement about the necessary linearized equation, not
yet an all-degree solution of the nonlinear equation `H(x+HV)=H(x)`.

Let `K` be the submodule of `M` generated by the six pairwise Koszul
syzygies among `(H_P,H_B,H_C,-H)`.  A second exact characteristic-zero
calculation determines the quotient

\[
 Q=M/K.
\]

Its minimal resolution is

\[
\begin{aligned}
0\longrightarrow{}&R(-45)
\longrightarrow R(-44)^5\\
\longrightarrow{}&R(-42)^4\oplus R(-43)^3\\
\longrightarrow{}&R(-25)^6\oplus R(-30)^3\oplus R(-31)^3\\
\longrightarrow{}&R(-22)^2\oplus R(-23)^{13}
\longrightarrow Q\longrightarrow0.                 \tag{3.8}
\end{aligned}
\]

Equivalently,

\[
 \operatorname{Hilb}(Q;t)=
 \frac{
 2t^{22}+13t^{23}-18t^{24}+6t^{25}
 -3t^{30}-3t^{31}
 +4t^{42}+3t^{43}-5t^{44}+t^{45}
 }{(1-t)^4}.                                         \tag{3.9}
\]

It has Krull dimension two and degree 296.  Koszul homology vanishes away
from the singular scheme `(H,H_P,H_B,H_C)`, so all non-Koszul logarithmic
classes are confined to the cone over the singular curves of the projective
discriminant.  Thus the all-degree wall problem splits into a generic
Koszul calculation and finitely many singular-support charts.  This is a
structural reduction, not a proof that either part contains no exact
symmetry.

### 3.1. Exact singular-scheme decomposition

The projective singular scheme has exactly four minimal components.  Two are the
lines at infinity

\[
 (Z,P),\qquad (Z,C).                                  \tag{3.10}
\]

The other two are the projective closures of affine prime curves.  They are
most economically specified as contractions of prime parameter ideals.
The triple-root curve is the elimination of `s` from

\[
\begin{aligned}
0={}&15u_5P^5s^4+8u_4P^4s^3+3Ps^2-1,\\
B={}&20u_5P^5s^3+12u_4P^4s^2+6Ps,\\
C={}&2\left(u_5P^5s^5+u_4P^4s^4+Ps^3-\frac B2s^2+s\right).
\end{aligned}                                         \tag{3.11}
\]

Its affine target dimension is one and its projective degree is seventeen.
For the two-double-root curve, write the root partition as

\[
 E(S)=u_5P^5(S^2-xS+y)^2(S-u).
\]

Eliminating `(x,y,u)` from

\[
\begin{aligned}
0={}&u_5P^5(u+2x)+u_4P^4,\\
0={}&u_5P^5(2xu+x^2+2y)-P,\\
0={}&u_5P^5(2xyu+y^2)-1,\\
B={}&2u_5P^5\bigl((x^2+2y)u+2xy\bigr),\\
C={}&2u_5P^5y^2u
\end{aligned}                                         \tag{3.12}
\]

gives an affine prime target curve of projective degree nineteen.  Exact
Macaulay2 calculations over `QQ` prove that both parameter ideals are prime,
that their contractions are prime, and that both contractions contain the
affine singular ideal.

Here is the exhaustion argument, which is geometric rather than a
`minimalPrimes` computation.  On `PZ != 0`, equation (2.1) is a genuine
quintic.  At a polynomial with exactly one double root, the discriminant
hypersurface is smooth and its normal is evaluation at that root.  The
`C`-direction changes the constant coefficient by `-1/2`, so the
three-parameter coefficient slice is transverse there.  Its singular points
therefore have root partition `3+1+1` or `2+2+1`, precisely (3.11) and
(3.12).  An exact saturation verifies that the chart `P=0, Z!=0` is empty.
Finally the unique degree-sixteen term (2.4) gives, and an independent exact
radical calculation verifies,

\[
 \sqrt{(H,H_P,H_B,H_C,Z)}=(Z,PC)
 =(Z,P)\cap(Z,C).                                    \tag{3.13}
\]

Consequently

\[
\boxed{
\operatorname{Min}(H,H_P,H_B,H_C)
=
\{I_{3+1+1},I_{2+2+1},(Z,P),(Z,C)\}.
}                                                     \tag{3.14}
\]

This replaces the previously running blind characteristic-zero
`minimalPrimes` attack.  It classifies the four irreducible charts on which
the support of `Q` can lie, but does not prove that every chart occurs in
that support, determine the module's scheme-theoretic multiplicities there,
or solve the exact polynomial symmetries supported there.

## 4. Stable point-orbit theorem through degree twenty-eight

Consider a stable self-equivalence after adjoining `s` identity variables,
and first suppose its target automorphism has total degree at most
twenty-two.
If it carries a standard marked target `(y_R,0)` to another standard target
`(z,0)`, specialize the first-three-coordinate quotient vector in (2.5) at
the identity variables `t=0`.  Its degree is at most

\[
 22-\deg H=6.
\]

Equation (3.2) still holds after this specialization.  The first seven
zeroes in (3.3) force the specialized quotient vector to vanish.  Hence the
first three target coordinates do not move:

\[
\boxed{
 \deg T\le22
 \quad\Longrightarrow\quad
 z=y_R.
}                                                       \tag{4.1}
\]

It remains first to justify degree twenty-three.  In degree
twenty-three, specialize the first three target components at the identity
variables `t=0`.  Equation (3.3) gives

\[
 T_{P,B,C}(x,0)
 =
 x+H(x)(\lambda V_1(x)+\mu V_2(x)).                  \tag{4.2}
\]

Both logarithmic fields have `kappa=0`, so exact preservation of the prime
boundary gives

\[
 H(T_{P,B,C}(x,0))=H(x).                              \tag{4.3}
\]

The degree-seven leading-component maps

\[
 (\lambda,\mu)\longmapsto
 (\lambda V_1+\mu V_2)_{P,7},
 \qquad
 (\lambda,\mu)\longmapsto
 (\lambda V_1+\mu V_2)_{C,7}                         \tag{4.4}
\]

both have rank two.  Therefore every nonzero pair `(\lambda,\mu)` gives
nonzero leading `P`- and `C`-components.  Since the unique degree-sixteen
term of `H` is a multiple of `P^12C^4`, the left side of (4.3) then has the
nonzero degree-368 term

\[
 \operatorname{in}_{16}(H)^{16}
 (\lambda V_1+\mu V_2)_{P,7}^{12}
 (\lambda V_1+\mu V_2)_{C,7}^{4}.                   \tag{4.5}
\]

Every other monomial of `H` has total degree at most fifteen and cannot
cancel (4.5).  This contradicts (4.3), so `lambda=mu=0`.  The first three
coordinates of the standard marked target do not move.

For degrees twenty-four through twenty-eight, the logarithmic spaces are
too large for the single leading-component test (4.4).  The following
recursive exact test replaces it.  Let `L` be a rational linear subspace of
one of the displayed logarithmic spaces and put

\[
 d_i(L)=\max_{V\in L}\deg V_i,\qquad
 w_i(L)=
 \begin{cases}
 16+d_i(L),&L_i\ne0,\\
 1,&L_i=0.
 \end{cases}                                           \tag{4.6}
\]

These are the componentwise degrees of `x+HV` at a generic point of `L`.
Expose the Newton face of `H` with the weight vector `w(L)`.  At every
subspace reached in the calculation, that face consists of one monomial

\[
 cP^aB^bC^c.                                           \tag{4.7}
\]

If the relevant leading components of `V` are all nonzero, (4.7) gives a
nonzero unique top term of `H(x+HV)`.  It cannot occur in `H(x)`.
Consequently an exact solution of

\[
 H(x+HV)=H(x)                                          \tag{4.8}
\]

must lie in the kernel of at least one participating leading-component
map.  Intersect `L` with each such kernel and repeat.  This branches over
all possible component-degree drops, so it is exhaustive; it does not
choose sample coefficients.

The exact rational recursion gives

\[
\begin{array}{c|rrrrr}
\deg V&8&9&10&11&12\\ \hline
\deg T&24&25&26&27&28\\
\text{logarithmic nullity}&21&54&104&174&267\\
\text{visited subspaces}&10&20&33&56&81\\
\text{unresolved faces}&0&0&0&0&0.
\end{array}                                            \tag{4.9}
\]

Thus no nonzero logarithmic field in these five layers satisfies exact
boundary preservation.

> **Bounded stable target-orbit theorem.**
>
> The orbit of `y_R` under stable polynomial target self-equivalences of
> total target degree at most twenty-eight is a point.

This statement allows arbitrary identity stabilization and arbitrary mixing
with the identity variables.  It is a statement about the standard
zero-section marked targets; it does not classify the action away from that
zero section.

## 5. All-degree stable marked rigidity

The second Newton vertex

\[
 D=(2,5,1)
\]

is an **intruder**: all three coordinates are positive.  The positive weight
`(1,3,1)` exposes it uniquely in the Newton polytope of `H`.

We use Corollary 3.3 of
[Derksen--Hadas--Makar-Limanov, *Newton polytopes of invariants of additive
group actions*](https://sites.lsa.umich.edu/hderksen/wp-content/uploads/sites/614/2018/05/A.I.a.8.pdf):
every vertex of the Newton polytope of a coordinate polynomial lies on a
coordinate hyperplane.

Let an unstabilized target automorphism fix the boundary normalization.
By (2.5), each component is

\[
 T_i=x_i+HV_i.
\]

Suppose `V_i` is nonzero.  Newton polytopes multiply by Minkowski sum:

\[
 N(HV_i)=N(H)+N(V_i).
\]

Choose a generic exposing functional in the open normal cone of `D`, and
let `v` be its unique maximizing vertex on `N(V_i)`.  Then `D+v` is a
vertex of `N(HV_i)`.  It has three strictly positive coordinates because
`D` does and `v` has nonnegative coordinates.  The monomial `x_i` cannot
cancel it.  Thus `D+v` is an intruder vertex of the coordinate polynomial
`T_i`, contradicting the cited theorem.  Hence every `V_i=0` and

\[
 \boxed{T=\operatorname{id}}
\]

without a degree bound.

The generic degree-five cover has monodromy `S_5`; its point stabilizer
`S_4` is self-normalizing, so its deck group is

\[
 N_{S_5}(S_4)/S_4=1.
\]

The source component is consequently also the identity.

> **Unstabilized all-degree stabilizer theorem.**
>
> Every polynomial left--right self-equivalence of `F_R` without target
> stabilization is the identity, in every degree.

This argument alone does not settle stabilization.  In
`k[P,B,C,t_1,\ldots,t_s]`, the vertex `(D,0,\ldots,0)` lies on every new
coordinate hyperplane, so the coordinate-polynomial theorem permits it.
The stable-invariant theorem below supplies the missing descent.

There is nevertheless a complete one-variable continuation.  Suppose
`s=1`, with identity variable `t`.  If some first-coordinate quotient
`V_i(P,B,C,t)` has positive `t`-degree, expose a vertex of maximal positive
`t`-degree while keeping `D` exposed in `N(H)`.  The corresponding vertex
of `H V_i` has all four coordinates positive, again contradicting the
coordinate-polynomial theorem.  Hence

\[
 T_{P,B,C}=f(P,B,C)
\]

is independent of `t`.  Write the fourth component as `g(P,B,C,t)`.
The Jacobian is block triangular:

\[
 \det DT=\det Df\;\frac{\partial g}{\partial t}\in k^\times.
\]

Both factors are therefore constants, so `g=a t+b(P,B,C)` with
`a != 0`.  If `f(x)=f(x')`, choose `t'` so that
`a t+b(x)=a t'+b(x')`; injectivity of `T` gives `x=x'`.  Thus `f` is an
injective polynomial endomorphism of affine three-space, hence an
automorphism in characteristic zero.  The unstabilized intruder argument
now gives `f=id`.

> **One-stabilization marked-orbit theorem.**
>
> After adjoining one identity variable, the first three coordinates of
> every stable target self-equivalence are the identity in every degree.
> Consequently its orbit on standard marked targets is a point.

The arbitrary-stabilization step uses Theorem 1.1 of
[Kuroda, *Initial forms of stable invariants for additive group
actions*](https://arxiv.org/abs/1304.0313).  In the terminology of that
paper, a polynomial with an intruder cannot be a stable `G_a`-invariant.
Equivalently, if a `G_a`-action on a larger polynomial ring fixes `H`, then
it fixes the whole base ring `k[P,B,C]`.

Let

\[
 {\cal A}=k[P,B,C,t_1,\ldots,t_s]
\]

and let `T` be a stable target self-equivalence.  Boundary preservation and
the same statement for `T^{-1}` give

\[
 T(H)=\rho H,\qquad T^{-1}(H)=\rho^{-1}H.             \tag{5.5}
\]

For each standard translation `tau_j` in `t_j`, conjugate it by `T`:

\[
 \sigma_j=T^{-1}\tau_jT.
\]

Equation (5.5) gives `sigma_j(H)=H`.  Since `H` has the intruder `D`,
Kuroda's theorem forces `sigma_j` to fix `P,B,C`.  Applying `T` to this
identity shows

\[
 \tau_j(T(P))=T(P),\quad
 \tau_j(T(B))=T(B),\quad
 \tau_j(T(C))=T(C).
\]

Thus the first three components of `T` are independent of every identity
variable.  Apply the same argument to `T^{-1}`.  Both `T` and `T^{-1}`
preserve `k[P,B,C]`, so their restrictions are inverse polynomial
automorphisms of affine three-space.  The unstabilized theorem makes that
restriction the identity.

The source restriction is also the identity because the generic cover has
deck group `N_{S_5}(S_4)/S_4=1`.  Vertical automorphisms of the identity
variables may remain; they are stabilization gauge and do not move the
standard marked `(P,B,C)` target.

> **All-degree stable marked-orbit theorem.**
>
> For every `s>=0`, the first three target components of every stable
> polynomial self-equivalence of `F_R times id_(A^s)` are the identity.
> Consequently the orbit of `y_R` on standard zero-section marked targets
> is a point, without any degree bound.

This closes the fixed-quintic marked stabilizer gate.  It does not claim
that the full stable automorphism group is trivial: arbitrary vertical
automorphisms can act on the identity factors without changing the marked
Keller target.

It also closes the earlier generic-fibre formulation.  Put

\[
 K=k(h),\qquad A_H=K[P,B,C]/(H-h).
\]

If a `K`-locally nilpotent derivation of
`A_H[t_1,\ldots,t_s]` moves an element of `A_H`, clear its finitely many
denominators by a nonzero polynomial in `H`.  That polynomial is invariant,
so the cleared derivation is still locally nilpotent on
`k[P,B,C,t_1,\ldots,t_s]`, fixes `H`, and moves the base ring.  Kuroda's
theorem forbids it.  Hence

\[
 \operatorname{LND}_K(A_H[t_1,\ldots,t_s])
 =
 \operatorname{LND}_{A_H}(A_H[t_1,\ldots,t_s])
\]

for every `s`: the generic fibre is stably rigid.

### 5.1. The generic-fibre topology route is unnecessary

For completeness, an exact Newton calculation checks all forty-six
nontrivial faces of `H-h` and its coordinate restrictions with Singular.
They are torus nondegenerate.  Their normalized-volume contributions are

\[
 (8,2,0),\qquad -(38,52,2),\qquad 328,
\]

so Khovanskii's torus-stratification formula gives

\[
 \chi(H=h)=246.                                       \tag{5.6}
\]

The generic fibre is a smooth connected affine surface and hence has the
homotopy type of a real CW complex of dimension at most two.  Therefore

\[
 b_2=245+b_1\ge245.
\]

In particular the `H^2=0` sufficient criterion of
[Bandman--Makar-Limanov](https://arxiv.org/abs/math/9807146) is unavailable.
This negative topology certificate is not used in the stable marked proof;
Kuroda's theorem acts directly on the ambient polynomial ring.

### 5.2. Independent bounded Jacobian certificate

In the unstabilized target, the preceding argument already makes the target
map the identity.  There is also an independent Jacobian certificate for
the degree-twenty-three logarithmic frontier.  It permits exactly

\[
 V=\lambda V_1+\mu V_2.                               \tag{5.1}
\]

Both fields have `kappa=0`, so (3.1) has `rho=1`.  Because `T` fixes `H=0`
pointwise, its tangent action there is the identity and its normal action is
`rho`; hence a polynomial automorphism must satisfy

\[
 \det DT=1.                                            \tag{5.2}
\]

Evaluate `det DT-1` at

\[
 (1,1,0),\qquad(1,0,1),\qquad(1,1,1).                 \tag{5.3}
\]

The three exact cubic equations in `(\lambda,\mu)` have reduced
characteristic-zero Groebner basis

\[
 \boxed{(\lambda,\mu).}                               \tag{5.4}
\]

Therefore `T=id`.  The source component is also the identity by the deck
calculation above.

Combining the same argument in degrees twenty-four through twenty-eight
with that deck calculation gives the independent bounded conclusion:

> **Unstabilized bounded stabilizer theorem.**
>
> Every polynomial left--right self-equivalence of `F_R` whose target
> component has degree at most twenty-eight is the identity.

## 6. The all-degree Newton reduction

The first unchecked stable degree is twenty-nine.  It corresponds to
`deg(V)=13` in (3.2).  The degree-eight logarithmic space already has
dimension twenty-one and evaluation rank three at `y_R`, so infinitesimal
motion alone does not predict an exact symmetry.  The Newton recursion
eliminates that space and the next four layers before any source-lift
condition is needed.

For orientation, at degree eight the leading
degree-eight `P`- and `C`-component maps each have rank seven, so each kernel
has dimension fourteen; their combined map has rank twelve, so the
intersection has dimension nine.  The first Newton step puts an exact
candidate in the union

\[
 \ker(V\mapsto V_{P,8})
 \ \cup\
 \ker(V\mapsto V_{C,8}).                              \tag{6.1}
\]

The recursion proves that every branch inside (6.1) also dies.  The same
statement holds through quotient degree twelve.  Across all five exact
trees, every exposed face is one of the two vertices

\[
 P^{12}C^4,\qquad P^2B^5C,                             \tag{6.2}
\]

and it always contains a component whose leading form is active on that
branch.  A simpler conjecture that the `B`-component is always maximal is
false: quotient degree twelve already has two branches with vanishing
`P`-component and `C`-degree respectively one and two above the `B`-degree.

Infinitesimal liftability alone will not settle this.  Since `det DF_R=1`,
target vector fields pull back to polynomial source vector fields, and an
identity variable can compensate divergence at first order.  The obstruction
must detect algebraization, finite polynomial complexity, or a discrete
feature of the full decorated normalization.

There is now an exact all-degree reduction of the Newton part.  Put

\[
 A=(12,0,4),\qquad D=(2,5,1).                         \tag{6.3}
\]

Every exponent `u=(i,j,k)` in the support of `H` is coordinatewise dominated
by a point

\[
 D+t(A-D)=(2+10t,5-5t,1+3t),\qquad 0\leq t\leq1.
                                                               \tag{6.4}
\]

Indeed the required interval is

\[
 \max\left(0,\frac{i-2}{10},\frac{k-1}{3}\right)
 \leq t\leq
 \min\left(1,\frac{5-j}{5}\right),                   \tag{6.5}
\]

and the 59 exact inequalities all hold.  No support point other than `A`
and `D` lies on the segment.  It follows that for every strictly positive
weight vector `w`, the exposed face is exactly `A`, exactly `D`, or their
common edge.  Thus there are no hidden higher-degree Newton vertices.

For the component degrees

\[
 w_P=\deg(P+HV_P),\quad
 w_B=\deg(B+HV_B),\quad
 w_C=\deg(C+HV_C),                                   \tag{6.6}
\]

the wall is the single equation

\[
 \boxed{10w_P-5w_B+3w_C=0.}                          \tag{6.7}
\]

Off this wall the old monomial argument is valid and kills every nonzero
candidate.  On the wall, if `p,b,c` are the leading forms of the three
components and `h_A,h_D` are the nonzero coefficients of the two vertex
monomials, the top equation is

\[
 h_Ap^{12}c^4+h_Dp^2b^5c=0,
\quad\text{or}\quad
 b^5=-\frac{h_A}{h_D}p^{10}c^3.                     \tag{6.8}
\]

Unique factorization makes this much more rigid than an arbitrary
two-monomial cancellation.  Over the algebraic closure it forces

\[
 c=\gamma r^5,\qquad b=\beta p^2r^3,\qquad
 \beta^5=-\frac{h_A}{h_D}\gamma^3                  \tag{6.9}
\]

for a homogeneous polynomial `r` and nonzero constants `beta,gamma`.

The wall cannot simply be declared absent from the logarithmic module.
It is already reachable inside the Koszul submodule.  Consider the
`P`-zero ladder

\[
\begin{aligned}
 V_P&=0,\\
 V_B&=L H_C,\\
 V_C&=-L H_B-\eta G H,\\
 Q&=-\eta G H_C.
\end{aligned}                                        \tag{6.10}
\]

This is a logarithmic field identically.  Since

\[
 \deg H_B=13,\qquad \deg H_C=15,
\]

take

\[
 \deg L=1+3n,\qquad \deg G=18+5n.
\]

The component weights are

\[
 (w_P,w_B,w_C)=(1,32+3n,50+5n),                     \tag{6.11}
\]

so every rung lies on (6.7).  At `n=0`, condition (6.8) first forces

\[
 G=PC^2R^5,\qquad R^3=P^7CL.                         \tag{6.12}
\]

Here `L` is linear.  If `L` is proportional to `P` or `C`, the exponent
pairs are respectively `(8,1)` and `(7,2)`, neither divisible by three; any
other linear `L` introduces an irreducible factor to exponent one.  Thus
the first tie, at target degree fifty, cannot cancel in this ladder.

At the next rung take

\[
 L=P^2C^2,\qquad R=P^3C,\qquad G=P^{16}C^7.          \tag{6.13}
\]

The leading `B`- and `C`-components of `T=x+HV` are

\[
 4h_A^2P^{26}C^9,\qquad
 -\eta h_A^2P^{40}C^{15}.                            \tag{6.14}
\]

Both Newton vertices then give \(P^{172}C^{60}\), and their coefficients
cancel when

\[
 \eta^3=4^5h_Dh_A^3.                                 \tag{6.15}
\]

This is an exact leading-order Koszul wall at target degree fifty-five.
It is not an identity `H(T)=H`, and it does not prove rigidity in degrees
twenty-nine through fifty-four: other Koszul combinations and the
singular-support quotient still have to be excluded.

The first nonlinear continuation of this wall can nevertheless be closed
exactly.  Keep the normalized two-generator `P`-zero slice

\[
\begin{aligned}
 V_B&=LH_C,\\
 V_C&=-LH_B-\eta GH,
\end{aligned}
\qquad
\eta^3=4^5h_Dh_A^3,                                  \tag{6.16}
\]

and allow

\[
 L=P^2C^2+L_{\leq3},\qquad
 G=P^{16}C^7+G_{\leq22}.                              \tag{6.17}
\]

An exact sparse homogeneous recursion expands `H(T)-H` from total degree
232 downward.  Each new homogeneous piece of `G` enters linearly by a
nonzero scalar times

\[
 P^{156}C^{53}.                                      \tag{6.18}
\]

It is therefore uniquely determined whenever the current residual is
divisible by (6.18); a residual outside that monomial ideal is a genuine
equation on the lower coefficients of `L`.  Projecting to `B`-degree zero
is exhaustive for these equations, since coefficients containing `B`
cannot contribute to that projection.

The first four layers force the coefficient `l_(3,0)` of `P^3` in
`L_{\leq3}` to satisfy

\[
 l_{3,0}^4=0,
\]

hence `l_(3,0)=0` over the ground field.  The eighth layer leaves one
explicit equation in the coefficients `l_(2,0)` and `l_(2,1)`.  After
reducing by that equation, the ninth layer supplies five projected
equations; their exact rational Groebner basis is `(1)`.  Consequently:

> **Degree-55 normalized-slice obstruction.**
>
> No field in the two-generator family (6.16)--(6.17) satisfies
> `H(x+HV)=H(x)`, despite its exact leading cancellation.

This closes the most economical cancellable wall, including every lower
term of `L` and `G` in that slice.  It does not cover the third `P`-zero
Koszul coefficient, which adds an independent multiple of `H` to `V_B`,
nor the singular-support quotient (3.8).  Those are the next two exact
charts.

For comparison, a simpler but higher-degree example comes from the
`H`-multiple submodule.  Let

\[
 W=(1,kP^{39}C^{14},PC^2),\qquad
 k^5=-h_A^{17}/h_D,\qquad V=HW.                     \tag{6.19}
\]

Then `V(H)=W(H)H`, so `V` is logarithmic, and

\[
 T=x+HV=x+H^2W
\]

has component degrees `(32,85,35)`.  The `A`- and `D`-contributions have
the same monomial \(P^{388}C^{136}\), and their coefficients cancel
exactly.  This does **not** give a symmetry: lower terms remain, and no
identity `H(T)=H` or polynomial source lift is claimed.  It does prove that
the former unrestricted monomial-avoidance gate was false.

Thus full marked non-descent remains open, but the correct all-degree target
is now narrower.

> **Binomial-wall termination gate.**
>
> Show that no nonzero logarithmic field satisfying (6.7)--(6.9) can satisfy
> the full identity `H(x+HV)=H(x)` together with the polynomial-automorphism
> and source-lift conditions.

The bounded theorem through degree twenty-eight is unchanged: none of its
exact branches reaches the wall.  Increasing the bound alone is therefore
secondary.  Resolution (3.8) isolates the remaining non-Koszul work on the
singular curves of `H=0`; the generic part can now be attacked directly
through the six Koszul coefficients and the power equation (6.8).

## Exact regression

Run

```bash
.venv/bin/python scripts/verify_rank_five_stable_target_stabilizer.py
```

The checker verifies the discriminant factorization and irreducibility, its
degree and support, the all-degree positive-upper-hull certificate and the
explicit cancellable wall, the logarithmic nullity table through quotient
degree twelve, rank-two and rank-three evaluations at the selected target,
the exact recursive Newton-face certificates, and the three-point Jacobian
Groebner basis.  The optional command

```bash
.venv/bin/python scripts/verify_rank_five_stable_target_stabilizer.py \
  --module-resolution
```

requires Macaulay2 and verifies the all-degree graded resolution (3.6).
The longer command

```bash
.venv/bin/python scripts/verify_rank_five_stable_target_stabilizer.py \
  --research-koszul-homology
```

verifies the exact quotient resolution (3.8), Hilbert series (3.9),
dimension, and degree.

The degree-55 normalized-slice obstruction is reproduced by

```bash
.venv/bin/python scripts/verify_rank_five_stable_target_stabilizer.py \
  --research-koszul-hensel --research-depth=4
.venv/bin/python scripts/verify_rank_five_stable_target_stabilizer.py \
  --research-koszul-hensel --research-depth=9 \
  --research-zero-l=l_3_0 --research-continue-constraints
```

The first command gives `l_3_0^4`; the second carries the resulting reduced
branch through layers eight and nine and returns the unit constraint ideal.
