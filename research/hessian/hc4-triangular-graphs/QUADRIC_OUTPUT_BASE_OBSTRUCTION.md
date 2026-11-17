# A quadratic observable cannot repair a nonautomorphic Keller output base

14 November 2026. The stable hypersurface obstruction extends from an
ambient coordinate pair to a linear observable paired with a nondegenerate
quadratic observable. Their generic ambient fibres need not be affine
spaces. No classification or rectification of affine-plane embeddings in
a quadric is needed.

## 1. Statement (HC4MYHQ0)

Let K:C^2_z->C^2 be a plane Keller map, r>=2, and
j:A^(r+2)_(z,w)->A^(r+3) a polynomial closed embedding, with r free
w variables. Let ell be a nonconstant affine-linear ambient function
and Q any ambient quadratic polynomial with nondegenerate Hessian.
If

\[
                         (\ell,Q)\circ j=K(z),                 \tag{1}
\]

then **K is a polynomial automorphism**. No Hessian condition on a source
potential is required. The result covers arbitrary polynomial j, all
source degrees and every r>=2. Closed embedding means global generation
of the source algebra; a substitution or generic immersion does not suffice.

The two functions in (1) need not be an ambient coordinate pair. The
nonisotropic normal case has a smooth quadric generic fibre of dimension
r+1, rather than affine space. Both normal types are included below.

## 2. The two ambient algebras and their filtrations

Center Q by an affine ambient change and translate ell by its value at
that center. Complete the quadratic form relative to ell. Up to an
invertible affine change and an invertible triangular change of the two
target functions, (1) has one of the forms

\[
\begin{array}{ll}
\text{I:}& \Pi=(T,q_{r+2}(D)),\\
\text{II:}& \Pi=(T,TY+q_{r+1}(D)),
\end{array}                                                     \tag{2}
\]

where each q_n is nondegenerate in n variables. Type I occurs when the
normal of ell is nonisotropic in the inverse Hessian metric; subtracting
a constant multiple of T^2 from the second target component gives the
displayed form. Type II is the isotropic case. The target changes have
nonzero constant Jacobian and preserve invertibility of the Keller map.

Put A=C[T,U], F=Frac(A), and let R be the ambient coordinate algebra,
viewed as the following A-algebra:

\[
\begin{array}{ll}
\text{I:}&R=A[D_1,\ldots,D_{r+2}]/(q_{r+2}(D)-U),\\
\text{II:}&R=A[Y,D_1,\ldots,D_{r+1}]/(TY+q_{r+1}(D)-U).
\end{array}                                                     \tag{3}
\]

As a C-algebra R is still a polynomial ring in r+3 variables: eliminate
U. In each q choose a hyperbolic term D1*D2 with coefficient one.
The monomials not divisible by D1*D2 form an A-basis for R. Reduction
by this one monic relation gives existence and uniqueness of that normal
form, including after arbitrary coefficient field extensions or reductions.

Give all D variables and Y degree one and A degree zero. The associated
graded algebras are respectively

\[
                A[D]/(q_{r+2}),\qquad A[Y,D]/(q_{r+1}).       \tag{4}
\]

These remain domains after extending F to its algebraic closure and after
reducing A at any height-one prime and passing to that residue field.
Indeed a homogeneous quadratic of rank at least three is irreducible
over every algebraic closure, and the ranks here are at least four and
three. One can also see (4) directly in the normal basis: reduction by
D1*D2 never raises degree, and its top reduction is exactly q_n.
Thus filtered degrees of nonzero elements add. Degree zero is precisely
the coefficient field in every field fibre.

Two consequences replace the polynomial ambient fibres in the
[earlier stable hypersurface proof](STABLE_HYPERSURFACE_OBSTRUCTION.md).

**Generic factoriality.** Over an algebraic closure Omega of F, type II
eliminates Y since T is nonzero, so its algebra is Omega^[r+1]. In type I,
write q=XY+q_r(V). Its equation is XY+q_r(V)=u with u nonzero. Localizing
at X gives Omega[X,X inverse,V], a UFD. The element X is prime before
localization: its quotient has equation q_r(V)=u, irreducible since
r>=2 and u is nonzero. Nagata's criterion therefore makes the original
Noetherian quadric ring a UFD as well. Its only units are Omega nonzero,
also immediate from the additive filtered degree. The
[primary Nagata criterion](https://stacks.math.columbia.edu/tag/0AFU)
is basic ring algebra, not an embedding classification.

**Primitive composition in R.** Suppose h in R has zero coefficient
on the normal-basis monomial 1 and the gcd of its other A-coefficients
is one. Then

\[
                          F[h]\cap R=A[h].                   \tag{5}
\]

For an irreducible pi in A, localize A at (pi). At least one nonconstant
normal-basis coefficient of h survives modulo pi. Its reduction has
positive filtered degree, so its powers have distinct degrees by (4).
It is transcendental over the residue field. If coefficients of a
polynomial B(S) in F[S] had negative minimum pi-valuation, multiply by
the compensating power of pi. Its nonzero coefficient reduction could
not vanish at that transcendental element. Thus B(h) would retain that
negative valuation and could not lie in R. Applying this at every pi
proves (5). This is the global denominator step, including the singular
quadric levels; it is not a generic-chart assertion.

## 3. The disjoint components and the smoothness argument

Let L=C(z1,z2), using the stipulated Keller map for the embedding of
F=C(T,U) into L. The extension L/F is finite separable of some degree d.
As in the retained proof, C[z] tensor_A F=L: it is an algebraic finitely
generated F-domain, hence a field. The generic source fibre is L[w].

The image of j is a smooth prime hypersurface (E) in the polynomial
ambient algebra R. After base change from F to Omega its fibre is the
disjoint union of d copies of A^r. Since R tensor_A Omega is a UFD,
E factors there into d distinct irreducible factors e_i, one for each
component. Galois acts transitively, so the factors have the same
filtered degree and the same leading normal-basis monomial. Normalize
its coefficient to one.

For i!=j, e_j is a unit modulo e_i, because the components are disjoint.
Their quotient rings are Omega[w], whose units are constants. Thus

\[
                          e_j-c_{ij}=b_{ij}e_i,
                          \qquad c_{ij}\in\Omega^*.           \tag{6}
\]

Equal degrees and degree additivity force b_ij to be a scalar; leading
normalization makes it one. Consequently every nonconstant normal-basis
coefficient is shared by all factors and fixed by Galois. Write

\[
                         e_i=h_0-a_i,\qquad h_0\in R\otimes_A F,
\]

with distinct a_i. Subtract the constant normal-basis coefficient from
h0, clear denominators and divide the gcd of its remaining coefficients.
This gives a primitive h in R with F[h]=F[h0]. Factorization of E over
Omega first gives E in F[h], and (5) upgrades it to

\[
                      E=P(T,U,h),\qquad P\in A[S],\quad\deg_S P=d. \tag{7}
\]

The substitution A[S]->R is injective, since h has positive filtered
degree over F. A factorization of P would factor E, and a substituted
nonconstant factor cannot be a unit in the polynomial C-algebra R.
Thus P is irreducible over A and, by Gauss's lemma, over F.

On every geometric generic source component h is its constant label.
Hence globally

\[
                         h\circ j=a(z)\in C[z].              \tag{8}
\]

The positive-w coefficients vanish generically and therefore identically.
This also shows that P is the minimal relation of a over F and has degree d.

The element P_S(K,a) has no source zero. If it vanished at a source
point, differentiate the literal ambient identity (7). At that point,

\[
                            dE=P_T\,d\Pi_1+P_U\,d\Pi_2.      \tag{9}
\]

Pull back to the source. Since E composed with j=0 and Pi composed
with j=K(z), invertibility of DK gives P_T=P_U=0 there. Equation (9)
then gives dE=0 at an image point, contrary to smoothness of the image
hypersurface. Therefore P_S(K,a) is a unit of C[z,w], a nonzero
constant c. The lower-degree polynomial P_S-c vanishes at (K,a).
Minimality forces P_S=c identically, so P=c*S+b(T,U) and d=1.
The plane Keller map is birational, and the retained
[birational Keller endpoint](../hc4-relative-row-orbit/PROOF.md#6-properness-and-the-exact-missing-length)
gives its polynomial inverse. This proves HC4MYHQ0.

No rectification of individual A^r images inside a quadric has been
used. The proof only needs their constant units, the generic ambient
factoriality, the filtered denominator lemma and the smooth image.

## 4. HC4 and the exact Meng--Yang observables (HC4MYHQ1)

In Tang's nonautomorphic HC4 form f=D(z)+w dot K(z), r=2 and K is
the intrinsic cotangent gradient pair. The theorem implies that this
pair cannot be the restriction of any ambient (ell,Q') as above to
a genuine A4 embedding in A5. Q' is any nondegenerate ambient quadratic;
it need not be the potential itself. Such a realization would force
K invertible, and then the full cotangent gradient is invertible.

The same obstruction holds after any certified ambient polynomial
automorphism carrying the two observables to this linear-quadratic
form. This is an embedding statement about the stipulated output base,
so it requires no preservation of the potential or nonlinear source
Hessian. It does not say every extension of two source polynomials
has that ambient form.

There is also a direct exact-encoding consequence. Suppose a genuine
Meng--Yang embedding, in some polynomial source coordinates (z,w), has
native A=F(z), B=G(z), with J(F,G) a nonzero constant. In the retained
normal ambient coordinates T=A+M/(2L), the observable

\[
                  R=LT^2+NB=LT^2+NyZ+2NxW
\]

is a nondegenerate quadratic. Its output pair is
(F+M/(2L),L(F+M/(2L))^2+NG), again Keller. The theorem forces (F,G)
invertible. The source potential LA^2+MA+NB is then a polynomial
coordinate. Its entire source orbit has exactly the retained unit-ideal
endpoint for unit-Hessian existence. Thus this direct A/B encoding
cannot conceal a nonautomorphic plane Keller base; it does not resolve
the coordinate endpoint or arbitrary mixed physical-base data.

The new ambient recognition class is not just a disguised coordinate
pair. For an exact positive control take

\[
 j(s,t,u,v)=(T,x,y,z,w)=(s,1,t+uv,u,v),\quad
 \Pi=(T,T^2+xy-zw),\quad \Pi\circ j=(s,s^2+t).              \tag{10}
\]

This is a closed embedding with literal source inverse. Generic fibres
of Pi are xy-zw=c, c nonzero. Their Hodge polynomial is
(tau-1)*tau^2+(tau-1)*tau=tau^3-tau, whereas an ambient coordinate
pair has affine-three-space fibres with polynomial tau^3. Thus Pi
is not an ambient coordinate pair. The only Hodge input here is the
already retained additive/multiplicative invariant. For the independent
nondegenerate potential Q0=Tz+yw+x^2, this same control has

\[
                 Q0\circ j=1+su+tv+uv^2,
\]

with actual Hessian determinant one and explicit invertible gradient.
It illustrates the hypotheses and their compatibility with HC4; it
constructs no nonautomorphic Keller map.

## Scope and evidence

The threshold r>=2 is used in generic factoriality and the irreducible
leading quadrics. For r=1 a nonisotropic generic fibre can be the surface
XY-Z^2=1. Here XY=(Z-i)(Z+i) gives nonunique factorization into degree-one
irreducibles, so the proof must not assert factoriality. For r=0 the
generic binary quadric has nonconstant units. No lower-dimensional
extension or counterexample to that extension is asserted.

The [checker](verify_quadric_output_base.py) and
[receipt](receipt-quadric-output-base.json) retain exact algebraic controls
and source pins. Arbitrary r, arbitrary embedding degrees and the global
primitive-composition step are written proofs, not finite exhaustion.
The [primary algebra record](sources/quadric-output-base-import.json)
keeps Nagata's criterion separate from any geometric classification.

General embeddings whose intrinsic output pair has no such ambient
extension, unaligned two-affine cases, both plane endpoints, the row
equation and unrestricted HC4/JC2 remain open. No formal verification,
independent replay, external review or novelty priority is asserted.
