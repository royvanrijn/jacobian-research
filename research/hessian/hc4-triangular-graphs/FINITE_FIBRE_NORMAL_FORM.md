# Finite-fibre normal forms for every remaining cotangent obstruction

14 November 2026. The singleton row-orbit argument extends to every
finite reduced plane fibre. Combined with explicit plane shears, it
gives templates indexed by the full critical-point count n, independently
of graph degree and source automorphism word length.

For n>=1 define the rational polynomials

\[
 A_n(s)=\prod_{i=1}^n(s-i),\qquad
 B_n(s)\equiv A_n'(s)^{-1}\pmod {A_n},\quad\deg B_n<n,       \tag{1}
\]

and the four-variable potential

\[
 T_{n,D}=D(s)+A_n(s)P+B_n(s)UQ,\qquad\deg D<n.             \tag{2}
\]

Every hypothetical nonautomorphic unit-Hessian polynomial with n>=1
full critical points is polynomially right-equivalent to one of (2).
The degree of (2) is exactly n+1. Right equivalence here identifies the
polynomial functions by an actual source automorphism; it does **not**
transport a constant Hessian determinant off the critical scheme.

For each n>=2, acquiring unit Hessian anywhere in the full source orbit
of (2) is equivalent to a concrete SL2(C[s,U]) row problem. In particular
failure of HC4/JC2 is equivalent to solving one of those row problems.
This is an exact reduction, not a solution of that row problem
and not an unconditional embedding into the Meng--Yang hypersurface.

## 1. Finite-point row equivalence is determined by orientation (HC4MYHH0)

Let a=(a1,a2) and b=(b1,b2) generate the same proper radical ideal I
of a nonempty finite set Z in C^2. For a fixed delta in C nonzero,

\[
 b=aC\text{ for some polynomial }C,\ \det C=\delta
 \quad\Longleftrightarrow\quad
 J(b)(z)=\delta J(a)(z)\quad\hbox{for every }z\in Z.          \tag{3}
\]

The local Jacobians are nonzero because both ordered pairs generate
the reduced point ideals locally. For the reverse implication, choose
any polynomial matrix M with aM=b. At Z, differentiation kills the
terms containing a, so det M(z)=delta. Radicality gives
det M-delta in I=(b1,b2); write it mu*b1+nu*b2. Put

\[
 k_a=(-a2,a1)^T,\qquad C=M+k_a(\nu,-\mu).                  \tag{4}
\]

Then aC=b and adj(M)k_a=(-b2,b1)^T, so the rank-one determinant
identity gives det C=det M-mu*b1-nu*b2=delta. Its adjugate divided
by the field unit delta is an actual polynomial inverse. The forward
implication follows by differentiating at a=0. No inverse to the
plane polynomial map b is assumed.

The relevant general orbit invariant is known: see
[Luc Guyot, Theorem 3.5](https://arxiv.org/html/2012.03056v2),
on equivalence of generating pairs by their determinant class. The
[retained primary context](sources/finite-fibre-literature-context.json)
records that source. Equations (3),(4) prove the finite reduced plane
specialization directly, extending the already retained constructive
[HC4MYOR1 point-ideal calculation](../hc4-relative-row-orbit/PROOF.md#2-exact-classification-of-the-whole-orbit).
The general orbit theorem is not asserted new and is not an additional
unproved classification input to the argument here.

Equality of Jacobian orientations only at Z is sufficient for (3).
It does **not** assert J(aC)=det(C)J(a) away from Z: derivatives of
C contribute terms with coefficients a. This distinction is precisely
why the unresolved unit-Jacobian row problem survives.

## 2. Normalize every finite fibre and its potential (HC4MYHH1)

Suppose

\[
 f=D_0(z_1,z_2)+P F(z_1,z_2)+Q G(z_1,z_2),
                             \qquad J(F,G)=c\in C^\times,  \tag{5}
\]

and n=#V(F,G)>=1. Its full critical scheme is that finite reduced
plane fibre, by [HC4MYHG0](GLOBAL_VALUE_ORBITS.md). We construct
an actual source automorphism carrying f to (2).

### Three determinant-one shears place all points at prescribed positions

Choose a determinant-one linear base change (s,U) such that the s
coordinates a_i of the n points are distinct. This is possible by
avoiding the finitely many covectors annihilating a point difference.
Write the transformed points (a_i,b_i) and label them i=1,...,n.
Interpolate polynomials J,E of degrees below n with

\[
                    J(a_i)=i-b_i,\qquad E(i)=i-a_i.
\]

The three successive plane shears

\[
 U\mapsto U+J(s),\qquad s\mapsto s+E(U),\qquad
 U\mapsto U-s                                                \tag{6}
\]

send (a_i,b_i) to (i,0). Every shear has determinant one and a
polynomial inverse. Using the inverse of this base map in (5) keeps
the plane pair Keller with Jacobian c. Its zero ideal is now exactly

\[
                              I_n=(A_n(s),U).                \tag{7}
\]

This construction proves the needed finite-point transitivity without
importing an automorphism-group classification.

### Match orientations, then remove the entire off-fibre offset

Rescale the Q source coordinate by 1/c, giving a coefficient pair
k=(F,G/c) with Jacobian one and the same ideal I_n. Define

\[
                            r_n=(A_n,B_n U).                 \tag{8}
\]

It generates I_n because gcd(A_n,B_n)=1. Its Jacobian is A_n'*B_n,
which equals one at every (i,0). Thus HC4MYHH0 supplies C in
SL2(C[s,U]) with k=r_n C. The relative-linear source change in P,Q
given by C, or its inverse according to the direction of equivalence,
identifies the two coefficient rows.

Let D(s) be the degree-below-n interpolation of the critical values:
D(i)=D_0(i,0), after the preceding base changes. The difference
D_0-D is in I_n=(A_n,B_n U); write it a*A_n+b*B_n U. Translating
P,Q by -a,-b absorbs this entire difference. This proves (2), with
every change an actual polynomial source automorphism and with no
plane inverse assumed. For a supplied candidate the construction is
effective given polynomial ideal-membership certificates; no uniform
degree bound on these certificates or source automorphisms is claimed.

Tang's retained alternative supplies (5) for every hypothetical
nonautomorphic unit-Hessian polynomial. Therefore the result applies
to every remaining graph or mixed embedding that could supply such
a polynomial. The n=0 case is the coordinate orbit already identified
in HC4MYHG2; n=1 gives the nondegenerate quadratic orbit, after
subtracting its critical value.

### Exact properties of the templates

One explicit expression for (1) is

\[
 B_n(s)=\sum_{i=1}^n
                \frac{A_n(s)}{(s-i)(A_n'(i))^2}.             \tag{9}
\]

Hence A_n'*B_n-1 is divisible by A_n. The full critical points of
T_{n,D} are exactly

\[
 \left(s,U,P,Q\right)=
             \left(i,0,-\frac{D'(i)}{A_n'(i)},0\right),      \tag{10}
\]

with values D(i) and full Hessian determinant one. Its full native
Hessian determinant is nevertheless

\[
                     (A_n'(s)B_n(s))^2,                     \tag{11}
\]

nonconstant for every n>=2. The degree is n+1, because the monic term
A_n P cannot cancel, while deg B_n<n and deg D<n.

Its whole value profiles are precisely

\[
 E(T_{n,D}=\lambda)=t^3-nt+
                      \#\{i:D(i)=\lambda\}\,t^2.            \tag{12}
\]

Thus every template has reduced Morse points, uniform critical
Hessian value and exactly the whole-value Hodge profile required by
HC4MYHG0, even though its native Hessian is not a unit for n>=2.
These necessary invariants alone cannot settle the remaining problem.
This does not make the templates counterexamples.
Linear perturbations in future unit-Hessian coordinates pull back
through the unknown source automorphism to generally nonlinear terms;
the native reference row's perturbations are not substituted for them.

## 3. Full source unitization is exactly one relative-linear row problem (HC4MYHH2)

Fix n>=2 and any D of degree below n. The following are equivalent:

1. Some polynomial source automorphism Phi makes T_{n,D} composed
   with Phi have nonzero constant full Hessian determinant.
2. Some plane Keller pair has exactly n points over one target.
3. There is C in SL2(C[s,U]) satisfying

\[
                        J(r_n C)=1.                         \tag{13}
\]

For (1) to (2), the full critical count remains n>=2, so the resulting
gradient cannot be an automorphism. Tang and HC4MYHG0 give the plane
Keller pair with n origin preimages. For (2) to (3), translate that
target to zero, apply (6) to its reduced fibre and normalize the
constant plane Jacobian to one; (3),(8) give the required C.
For (3) to (1), the actual relative-linear source automorphism

\[
                   (s,U,w)\mapsto(s,U,C(s,U)w)              \tag{14}
\]

gives the literal potential

\[
                         D(s)+P(r_nC)_1+Q(r_nC)_2,
\]

whose full Hessian is one by the cotangent determinant identity.
Its n full critical points also prove that its gradient is not
invertible. All matrix degrees are unrestricted. The equivalence
holds for **every** D, although changing its critical-value list need
not preserve the full source orbit of the polynomial function.

Consequently JC2 fails if and only if (13) has a solution for some
n>=2. In the forward direction, a nonautomorphic plane Keller map
has generic degree d>=2 by the retained birational Keller case in
[HC4MYOR2](../hc4-relative-row-orbit/PROOF.md#6-properness-and-the-exact-missing-length);
translate a generic target and take n=d. The converse follows from
its n-point zero fibre. HCTANG2 gives the corresponding HC4 equivalence.
For existence of a counterexample one may therefore fix D=0 and the
explicit sequence of rows (8), rather than increase arbitrary graph
or source-word degrees without an invariant target.

This equivalence does not bound n, solve (13), or construct a matrix
C with unit plane Jacobian for n>=2. The rows r_n themselves have
nonconstant Jacobian. Applying a low native-degree HC4 theorem to a
template does not control the degree after (14). Likewise one cannot
transport the native Hessian zeros through a nonlinear source change.

For a supplied graph embedding the normal-form changes in Section 2
are source automorphisms and keep that actual embedded image. Starting
from an abstract template, however, this note supplies no independent
Meng--Yang graph embedding, no marked gradient-covector match and no
automatic transfer of a known ambient collision. It is a necessary
universal normal form for the graph route and an exact unrestricted
HC4/JC2 row reduction, with these distinctions preserved.

## Verification, sources and remaining work

[verify_finite_fibre.py](verify_finite_fibre.py) and its
[receipt](receipt-finite-fibre.json) check the determinant correction,
three-shear inverses, several exact templates, point ideals, full
critical points, Hessians, value-profile identities and retained inputs.
The arbitrary n and source-orbit arguments are written proofs, not
finite extrapolations. Multi-point template controls are non-Keller
reference rows, not realized plane counterexamples.

The present normal-form statements are over C. A polynomial gradient
inverse established over C for a stipulated characteristic-zero candidate
descends by uniqueness after embedding its finitely generated coefficient
field. No descent of the ordered
finite-point normalization over an arbitrary field is claimed.

The next substantive gate is a global restriction on (13), its plane
properness defect, or a graph-specific embedding obstruction that the
templates do not already encode. Finite critical/Hodge conditions that
all the templates satisfy do not replace that gate. Full HC4/JC2 and
both distinct unit-ideal and singleton-zero endpoints remain open.
No formal verification, independent replay, external review or novelty
priority is claimed.
