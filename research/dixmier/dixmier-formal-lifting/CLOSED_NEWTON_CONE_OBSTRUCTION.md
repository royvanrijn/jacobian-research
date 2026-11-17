# Closed Newton support and regular quantum parameters

**QPN2.** Parameter regularity closes the equality boundary in
[QPN1](MIXED_NEWTON_P_POWER_OBSTRUCTION.md). The resulting obstruction
allows genuinely mixed terms of arbitrarily high momentum order in a
closed support cone. It includes every regular exact-R scalar Long
proposal in this cone, including ramified quantum parameters.

## Statement

Let K be a characteristic-zero field, R=K[[t]], L=K((t)), and
eta=t^e*u with e>=1 and u in R a unit. Work in
\[
 D=L(v)[a;\eta\delta],\qquad \delta=v^{-1}\partial_v.
\]
Suppose S is a finite normal polynomial in v,a, all its scalar
coefficients lie in R, and n>=2. Write
\[
 S=g a^n+\mu v+\sum c_{j,m}v^m a^j,                         \tag{1}
\]
where the sum omits the two displayed monomials. Assume g and mu
are units in R. For every nonzero term require
\[
 \begin{array}{ll}
 m=2k:& j-k\le n,\\
 m=2k+1:&j\le k,
 \end{array}\qquad j,k\ge0.                                 \tag{2}
\]
Then S has no finite scalar T with [S,T]=kappa for nonzero
kappa in L, over L(v) or any finite algebraic coefficient cover.
The same exclusion holds for the fixed scalar image S*I_d in every
finite matrix size, by normalized matrix trace.

There is a second version in which g need only be nonzero in R:
replace the even inequality by j-k<n for every even term other
than g*a^n, retaining the weak odd inequality j<=k. The odd
coefficient mu must still be a unit. This version permits a
parameter-small top constant-coefficient momentum term.

The parameter must have positive order, and all scalar coefficients
of S must be regular. No bound on their t expansions, on finite
momentum order, or on a hypothetical partner's order is imposed.

## The coefficient and its leading parameter order

First work in characteristic p>2n+1 over a residue field k, with
all coefficients in k[[t]] and eta of positive order e. Put
\[
 J=(n-1)p-(2n+1),\qquad
 C_p(t)=[v^{-2p}a^J]S^p,\qquad
 \gamma_n=\frac{n^2-1}{24n}\in k^\times.                     \tag{3}
\]
For the closed cone (2), C_p is divisible by eta^(p+1), and
\[
 \left.\frac{C_p}{\eta^{p+1}}\right|_{t=0}
       =\gamma_n g(0)^{p-2}\mu(0)^2\ne0.                   \tag{4}
\]
Thus C_p is nonzero. For the second version,
\[
 C_p=g^{p-2}\eta^{p+1}
                  \bigl(\gamma_n\mu^2+\eta A_p(t)\bigr),
                       \qquad A_p(t)\in k[[t]],              \tag{5}
\]
which is again nonzero even when g has positive t valuation.

Here is a complete selection argument. Give a monomial v^m*a^j
weight w=j-m/2 and deficit n-w. Ore contractions preserve weight.
The target has total deficit 2n+1. Every even term has
nonnegative integer deficit; every odd term has deficit at least
n+1/2. A word with no odd terms stays in k[[t]][v^2][a;eta delta]
and has no negative v powers. Parity of the v exponent excludes
an odd number of odd terms. Consequently every contribution to
C_p has exactly two odd terms, both on their boundary j=k,
and all p-2 even terms lie on j=n+k.

Index each selected term by this nonnegative k. If Ksum is the
sum of those p indices, their raw v exponent is 2*Ksum+2.
Reaching -2p requires exactly
\[
                           p+1+Ksum                         \tag{6}
\]
derivative contractions. Its scalar contribution is therefore
divisible by eta^(p+1+Ksum). Each word and each normal-order
coefficient is a finite sum, so no convergence statement is needed.

If Ksum=0, all indices are zero. The word has p-2 copies of
g*a^n and two copies of mu*v. The universal computation in
[QPC1](UNIVERSAL_P_POWER_OBSTRUCTION.md) gives precisely
gamma_n*g^(p-2)*mu^2*eta^(p+1). Every other contribution is
divisible by eta^(p+2), with regular scalar factors. This proves
(4). In the second version the strict even inequality forces all
p-2 even terms to be g*a^n, so this common scalar factor can be
removed from every contribution, proving (5).

## Why the formal coefficients admit suitable good reductions

Arbitrary formal series over K do not come with a direct reduction
modulo p. The proof replaces a hypothetical pair by another pair
with the same finite supports and regularity properties before
reducing. It does not reduce the given infinite scalar expansions.

Trace first reduces a finite algebraic-cover mate to L(v), as in
QPC1. Scale its nonzero constant commutator to one. A finite
normal-order rational mate has a common nonzero polynomial
denominator D(v) and can be written
\[
       T=t^{-N}D(v)^{-1}\sum_{i=0}^{M}N_i(v)a^i,              \tag{7}
\]
with finite polynomials D,N_i in v whose scalar coefficients are
in K[[t]], and some N>=0. Multiplying numerator and denominator
by a scalar normalizes at least one D coefficient to a unit.

Fix these finite supports. Treat their scalar coefficients and
those of S as finitely many unknowns. Write eta=t^e*u, with u
a unit. Include a scalar inverse equation for u, mu, the selected
denominator coefficient, and g in the first version. In the
second version write g=t^r*g_unit using its actual finite order r
and include the inverse equation for g_unit instead.

Clear t powers and all powers of D(v) in [S,T]=1. The Ore
Leibniz rule uses only finitely many derivatives, so coefficient
comparison produces a finite system of polynomial equations over
K[t]. This system has a solution in K[[t]], by the hypothetical
pair and its regular scalar data.

The local ring K[t]_(t) is a G-ring. Apply the directly inspected
[Stacks Theorem 16.13.2](https://stacks.math.columbia.edu/tag/07QX),
with approximation order one, to obtain a solution in an etale
neighborhood with residue field K. Supports, the powers t^e and
t^r, and all unit equations are encoded in the system. The new
solution therefore retains every condition needed in (4) or (5).
Its common denominator remains nonzero.

An etale neighborhood has a finite presentation. The Jacobian
presentation and finite-presentation descent used here are given in
[Stacks Lemmas 10.143.2 and 10.143.3(8)--(10)](https://stacks.math.columbia.edu/tag/00U0).
Choose a finite
presentation containing this solution and its point over t=0;
spread its equations, its etale chart, that residue point and the
finitely many inverse equations over a finitely generated
characteristic-zero Z-algebra B contained in K. Clear the finite
denominators needed in these data. The commutator and inverse
equations remain identities in this model. A finitely generated
characteristic-zero ring admits finite-field reductions in
arbitrarily large characteristic after excluding finitely many
units, by QPC1's specialization argument. The spread residue
point is rational over the resulting finite field k, after finite
extension if necessary.

The reduced chart is etale at that point. Its formally etale
lifting property gives a compatible lift of the point to k[[t]],
by successive lifting modulo t^q; see
[Stacks Definition 10.150.1 and Lemma 10.150.3](https://stacks.math.columbia.edu/tag/00UP).
Thus the reduced polynomial
system has regular power-series solutions with u, mu and the
required leading g coefficient still units. The resulting finite
rational operators over k((t))(v) satisfy [S,T]=1. This construction
preserves orders by their encoded factorizations, rather than by
reducing an unbounded number of original coefficients.

The reduced constant field C=k((t)) is not perfect. The required
scalar center and central-degree statements nevertheless hold
over this field; the perfect-field formulation in QPC1 is not
being applied without this additional argument. For any field C
of odd characteristic p and eta in C nonzero, delta^p=0 on C(v):
on v^j the p successive factors j,j-2,...,j-2(p-1) contain zero
modulo p, and delta^p is a derivation. Its vanishing on v and C
therefore proves the assertion on C(v). Also ker(delta)=C(v^p).
PBW and commutation with v and a give the center
\[
                         C(v^p)[a^p].                       \tag{9}
\]
Indeed commutation with v removes every a order not divisible by
p, beginning with its greatest such order, and commutation with a
puts each remaining scalar coefficient in C(v^p).

For central degree, normalize b=(v/eta)*a; then [b,v]=1 and
the Ore ring is the rational-position localization of the usual
rank-one Weyl algebra over C. Its generic central localization
has the p^2-element PBW basis v^i*b^j, 0<=i,j<p, over
C(v^p,b^p). Over an algebraic closure Omega of that center,
choose roots alpha^p=v^p and beta^p=b^p. On Omega[z]/(z^p),
let v act by alpha+z and b by beta+partial_z. These operators
satisfy all the relations and generate the full p-by-p matrix
algebra: z*partial_z has the distinct eigenvalues 0,...,p-1,
its spectral projectors isolate the monomials, and z and
partial_z supply the matrix units. Equality of dimensions then
identifies the central localization after scalar extension with
M_p(Omega). Thus it is central simple of degree p without any
perfection hypothesis on C.

In this splitting representation a nonzero constant mate forces
S^p scalar by QPC1's characteristic-polynomial translation
argument, hence central in the original algebra. Its center (9)
has a orders divisible by p, whereas J is positive
and not divisible by p. The nonzero C_p from (4) or (5) contradicts
centrality. This proves both versions in characteristic zero.

The imported approximation theorem is Stacks 16.13.2. Its G-ring
hypothesis and etale conclusion, and the finite-presentation and
formal-lifting facts just used, were inspected; the sources and exact uses are pinned in
[the source audit](closed_newton_source_audit.json).
The finite-characteristic scalar degree, QPC1 coefficient identity,
trace descent and finite-type specialization remain dependencies
on QPC1. No arithmetic center map is prescribed.

## Long and the equality boundary

Over K(r)[[t]], the classical Long operator
\[
                     S_0=2v+4a-3r a^2
\]
has n=2, g(0)=-3r and mu(0)=2, both units. Hence every finite
regular deformation with even terms satisfying j-k<=2 and odd
terms satisfying j<=k is excluded, provided an actual finite
scalar partner is required. For example the closed cone admits
v^2 a^3, v^3 a, v^(2k) a^(k+2), and v^(2k+1) a^k.
It therefore covers unbounded generic momentum orders and the
odd equality terms missing from QPN1. Ramification eta=t^e*u
does not change the conclusion.

The cancellation at QPN1's odd equality boundary requires a
parameter pole in a simple quadratic control. Indeed, with
E=a+c*v^3/2, exact Ore multiplication gives
\[
 a^2+c v^3a+\mu v
       =E^2+(\mu-3\eta c/2)v-c^2v^6/4.                    \tag{8}
\]
Cancelling the odd coefficient uses c=2mu/(3eta), which has
a parameter pole when mu is a unit. This lies outside the
regularity hypotheses. It explains why a cancellation at a
single fixed nonzero quantum parameter does not contradict (4).
No mate is asserted for this boundary control.

The degree restriction n>=2 remains necessary: QPC1's explicit
positive pair for S=a+v has a finite rational mate. Negative-v
coefficients, higher even weights, odd terms above the boundary,
nonpositional first images and genuinely nonscalar carriers remain
outside the theorem. No bounded normalization puts every Long
representative into (2).

[Exact controls](check_closed_newton_obstruction.py) check parameter
orders in finite sparse Ore powers, both versions, closed even and
odd boundary terms, a ramified parameter, a parameter-pole
cancellation and positive finite pairs. They do not replay Artin
approximation, spreading or the universal word-selection proof.
No bounded full Long lift, DC1/DC2 solution, formal verification,
external review, independent whole-proof replay or novelty
priority is claimed.
