# Bounded formal lifting and automatic nonsurjectivity

14 November 2026. **QFL2.** Let f be a noninvertible polynomial symplectic
endomorphism of affine 2n-space in characteristic zero. If f has a formal
Weyl lift whose generator images have one common finite bound on spatial
degree, then there is a non-surjective endomorphism of A_n over C.

The hbar series need not terminate, and no evaluation at hbar=1 is needed.
Nonsurjectivity follows from the classical reduction, even if the generic
images no longer preserve its displayed reducible factorization. This is
a proved sufficient criterion. No such bound is established for Long's
map or the retained degree-five seed.

## The algebra that bounded degree supplies

For a characteristic-zero field K, put R=K[[t]], L=K((t)), and let
a(t) be a nonzero element of tR. Consider the finite PBW algebra
\[
 {\cal W}_R(a)=R\langle x_1,\ldots,x_n,p_1,\ldots,p_n\rangle/
 ([p_i,x_j]-a(t)\delta_{ij},[x_i,x_j],[p_i,p_j]).
 \tag{1}
\]
An element has finitely many spatial PBW monomials; their scalar
coefficients may be arbitrary power series in t. Its t-adic completion
is larger: spatial degrees may grow with t order.

Suppose the images F_i of a formal lift satisfy
\[
                \deg_{\rm spatial}F_i\le D
 \quad\hbox{for every image and every t coefficient}.       \tag{2}
\]
There are finitely many spatial monomials of degree <=D, so all F_i
belong to the uncompleted algebra (1). They therefore define an
endomorphism phi_R of (1). After inverting t, rescaling each p_i by
a(t)^{-1} identifies W_L(a) with A_n(L). This rescaling preserves spatial
degree.

No convergence or algebraicity of the scalar coefficient series is
required. This is the first useful improvement over demanding termination
of the formal series.

## A valuation lemma makes nonsurjectivity automatic

Assume the commutative reduction
\[
        f^*:K[x_1,\ldots,p_n]\longrightarrow K[x_1,\ldots,p_n]
\]
of phi_R is injective. For any nonzero g in W_L(a), define v(g) as the
minimum t valuation among its finitely many PBW coefficients. Then
\[
                         v(\phi_L(g))=v(g).         \tag{3}
\]

To prove this, write g=t^m g_0 with g_0 in W_R(a) and nonzero reduction
bar(g_0). Reducing phi_R(g_0) gives f^*(bar(g_0)), which is nonzero by
injectivity. Thus phi_R(g_0) has valuation zero, proving (3).

If phi_L were surjective, take a finite PBW preimage g_i of each standard
generator. Equation (3) gives v(g_i)=0, so every g_i already belongs to
W_R(a). Reduction supplies a polynomial preimage of every generator
under f^*. Hence f^* is surjective. We have proved
\[
 \boxed{\quad f^*\text{ injective and not surjective}
             \ \Longrightarrow\ \phi_L\text{ not surjective}.\quad} \tag{4}
\]

A polynomial symplectic map has invertible Jacobian, hence is dominant
and induces an injective polynomial-ring map in characteristic zero.
Thus (4) applies to every noninvertible f in the opening statement.
The argument also permits ramification a(t)=t^e u(t), u(0) nonzero.

The injectivity hypothesis matters. In the commutative example x -> t*x,
the generic inverse x/t has a pole and the special reduction is zero.
There is no such pole cancellation when f^* is injective.

For completeness, the counterexample over L yields one over C. The
finitely many PBW coefficients of its generator images generate a
finitely generated characteristic-zero field K_0 over Q. The restricted
A_n(K_0) endomorphism is not surjective, since surjectivity would persist
over L. Nonsurjectivity also persists under any field extension: the
nonzero vector-space cokernel tensors faithfully. Embed K_0 into C.
This proves the opening criterion.

The field embedding is an existence argument, not an explicit choice of
complex coefficients for a counterexample. Establishing (2) for a genuine
noninvertible classical symplectic map remains the unsolved step.

## Finite type geometry of the same condition

Fix n and D. Let M_(n,D) be the affine finite type scheme over K[h]
whose coordinates are the coefficients of 2n spatial polynomials of
degree <=D and whose equations are
\[
 { [F_i,F_j]_{\star_h}\over h}=J_{ij}.                \tag{5}
\]
These are polynomial equations in h and the finitely many coefficients.
The numerator is universally divisible by h. At h=0, equation (5)
is the full set of classical Poisson relations. At h nonzero it is
the full set of Weyl relations.

Write V_e for the spatial polynomials of degree <=e, and put
\[
                         E=D^{2n-1}.
\]
For each point m of M_(n,D), substitution through its assignment gives
\[
       T_m:V_E\longrightarrow V_{DE},\qquad g\mapsto\phi_m(g).
 \tag{6}
\]
Its matrix entries are polynomial in the parameters. At every geometric
point in characteristic zero, T_m is injective: at h=0 this follows
from the symplectic Jacobian, and at h nonzero from simplicity of the
Weyl algebra. Its rank is therefore the constant dim(V_E).

The inverse-degree theorem for polynomial automorphisms and its Weyl
analogue give degree(phi_m^{-1})<=E whenever phi_m is an automorphism.
For the Weyl statement see
[Lauritzen--Thomsen, Theorem 2.5](https://arxiv.org/html/2308.09384v2);
their Theorem 4.3 states the polynomial bound, and they reference
[Bavula's extension](https://arxiv.org/abs/math/0512215).

Consequently phi_m is an automorphism exactly when every standard
generator belongs to im(T_m). With the rank already constant, these
are closed determinantal conditions on [T_m | x_i] and [T_m | p_i].
The inverse-degree bounds are used only on field-valued geometric points.
They are not being asserted over truncated h rings with nilpotents; the
reduced-ring hypothesis in Lauritzen--Thomsen's Theorem 2.5 matters.
Thus:

> Within M_(n,D), the locus of automorphisms is Zariski closed.
> Its complement is a Zariski open locus of non-surjective assignments,
> including the classical noninvertible points at h=0.

This assertion concerns the normalized relation scheme (5), not
arbitrary polynomial maps without the injectivity condition.

In particular, around a noninvertible classical point f there is a
Zariski open neighborhood U such that every point of U with h nonzero
is a Weyl counterexample. This supplies a second proof that a
bounded-degree formal arc through f has noninvertible generic fiber.

For n=2 the inverse bound is E=D^3. No determinant of that large matrix
is computed here; the inverse bound is used to prove the all-degree
geometric statement.

## What Artin approximation does and does not provide

Suppose a degree-D formal lift exists. Its finitely many scalar
coefficient series solve the finite polynomial system (5).
[Artin approximation, Stacks Theorem 16.13.2](https://stacks.math.columbia.edu/tag/07QX)
gives a solution over an etale neighborhood of the local ring K[h]_(h)
with the same prescribed finite jet. Taking the approximation order
one preserves f. After shrinking and taking the component through that
point, this is an algebraic curve family with h nonconstant.

For algebraically closed K, the curve meets U and h nonzero, giving
actual non-surjective Weyl endomorphisms over K. Over general K the
coefficient-field argument above already gives a counterexample over C.

Approximation replaces the original scalar series by algebraic ones
with the same selected jet. It does not assert that the original series
is algebraic or has a value at h=1.

Without (2), QFL1 has infinitely many potentially nonzero spatial
coefficients. There is no single finite system (5) to which this
application of Artin approximation applies. Applying approximation
separately at each growing degree does not supply a common bound.

Equivalently, one may seek a finite-dimensional spatial carrier for
one full formal lift, rather than a terminating h expansion.

## Horizontal components and the all-degree question

At a fixed D, let H_D be the scheme-theoretic closure in M_(n,D) of
the open set h nonzero. In coordinate-ring terms its ideal is the
saturation I_D:h^infinity of the equations (5).

Any degree-D formal lift puts f in H_D. Conversely, over an algebraically
closed characteristic-zero field, a point of H_D can be reached by an
algebraic curve through f with h not identically zero. Normalize the
curve at a point above f and complete: h=t^e u(t) with u(0) nonzero.
This gives a bounded-degree ramified formal lift of f of the kind
allowed in (1). This is the usual curve-selection argument for the
closure of a dense open in an affine finite type scheme. More generally,
[Stacks Lemma 32.15.1](https://stacks.math.columbia.edu/tag/0CM2) supplies
a discrete valuation ring reaching a point in the closure of a finite
type image; the valuation proof above works directly over such a ring.

Thus the decisive construction question is
\[
 \boxed{\quad
  \text{Does a known noninvertible polynomial symplectic map belong to }
               H_D\text{ for some finite }D?
       \quad}                                      \tag{7}
\]
A positive answer already gives a Weyl counterexample by (4).
It does not need the stronger condition of a lift with integral h
powers or a terminating h series.

If DC_n is true, every noninvertible classical symplectic point is
vertical in every finite M_(n,D). QFL1 is compatible with this: its
formal solution may use increasing spatial degrees and never factor
through a single finite M_(n,D).

This is a structural reformulation, not a proposal to compute the
saturations for increasing D. The remaining unproved condition concerns
horizontal membership or its failure for all D simultaneously.

The property "some bounded-degree formal lift exists" is preserved by
source and target changes admitting exact quantum automorphisms with
bounded spatial degree, and with inverses having the same property.
Composition multiplies finite degree bounds. Positional polynomial
Weyl shears are included. An arbitrary formal gauge need not preserve
boundedness of a chosen representative.

## Two controls prevent false boundedness conclusions

First, in one Weyl pair with [p,x]=h,
\[
                  X=x,\qquad P=p+{x^2\over1-h}
 \tag{8}
\]
is an exact canonical pair over K[[h]], with spatial degree at most two.
Its h series does not terminate and its displayed formula has a pole
at h=1. Nevertheless it is a finite Weyl assignment over K((h)), and
its coefficients are regular at many other nonzero parameters. This
shows why h termination and evaluation at one designated parameter
are unnecessarily strong conditions.

Second, in symmetric symbols,
\[
 X_h={x\over1-h^2x},\qquad P_h=(1-h^2x)^2p            \tag{9}
\]
satisfies [P_h,X_h]_star=h exactly. The normalized bracket is
(1-h^2x)^2 times partial_x(X_h)=1; higher Moyal terms vanish because
the combined p degree is one. The classical reduction is the identity,
which also has the bounded lift (x,p). But the coefficient of h^(2r)
in X_h is x^(r+1), so this representative has unbounded spatial degree.

In normal PBW order the second image of (9) is
\[
 (1-h^2x)^2p-h^3(1-h^2x).
\]
These identities are checked by
[the exact controls](check_realization_controls.py).
Growth in one homotopy or gauge representative cannot prove that the
formal gauge class has no bounded representative.

## Classical reduction does not prescribe arithmetic support

The bounded polynomial quantum shear
\[
                     x\mapsto x,\qquad p\mapsto p+h^2x^2       \tag{10}
\]
also reduces to the identity at h=0. At h=1 it is a nontrivial Weyl
automorphism. In every prime characteristic ell>3,
\[
                  (p+h^2x^2)^\ell=p^\ell+h^{2\ell}x^{2\ell}.
 \tag{11}
\]
Thus its specialized Frobenius center map is a nontrivial shear, while
the identity lift of the same classical map has identity center action.

Here is an elementary all-prime proof of (11). Let A=p+h^2x^2.
The constant commutator [A,x]=h shows A^ell commutes with x.
It also commutes with A, hence with p=A-h^2x^2. Therefore A^ell is
central. The difference A^ell-p^ell has p order <ell, and commutation
with x forces it to belong to K[h,x]. A word containing j copies of
h^2x^2 contributes to its p-order-zero part only with exponents
h^(ell+j) x^(3j-ell). Centrality forces 3j-ell divisible by ell.
Since ell>3 and 0<=j<=ell, the only nonnegative surviving case is
j=ell. Its coefficient is one, proving (11). Finite PBW checks at
ell=5,7,11 are implementation controls, not the all-prime proof.

Kanel-Belov and Fedorov's
[Theorem 1.4(i)](https://arxiv.org/html/2610.33919v2)
concerns Long's prescribed Frobenius pullback on the entire center,
with equality on a Zariski-dense set of closed points of an arithmetic
model. Applying it to a representative in (7) would require proof of
that center comparison. Bounded spatial degree gives a Weyl assignment
whose arithmetic center map can differ from its classical h=0 map, as
(11) illustrates. The separate
[plane Keller graph filtration](PLANE_KELLER_GRAPH_BW_FILTRATION.md)
uses the actual graph modules with their specified center maps.

## Exact remaining gap

QFL1 supplies unrestricted all-orders formal existence. QFL2 shows
that a common spatial degree bound would already supply a Weyl
counterexample for a noninvertible classical map, with nonsurjectivity
automatic. Long's verified polynomial symplectic counterexample is
an available classical input. The existence of a bounded representative,
or even a bounded ramified representative, remains unknown.

Nothing here reduces arbitrary Weyl endomorphisms to Long's map.
An all-degree nonmembership result for this one map would close its
bounded-quantization route; it would not prove DC2.

No uniform bound, Weyl counterexample, DC2 solution, external review,
formal verification or novelty priority is claimed.
