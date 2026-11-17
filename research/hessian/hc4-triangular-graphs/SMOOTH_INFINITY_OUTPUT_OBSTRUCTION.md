# Arbitrary-degree output observables with a smooth boundary at infinity

14 November 2026. The linear/quadratic observable restriction can be
replaced by a geometric condition on the highest homogeneous forms.
This note isolates the filtered-algebra mechanism and supplies a
complete-intersection certificate in every degree. It does not require
rectifying any individual affine-space embedding in an ambient fibre.

## 1. A filtered output-base criterion (HC4MYHV0)

Let K:A2_z->A2 be a complex plane Keller map, r>=2, and let
j:A^(r+2)_(z,w)->A^(r+3) be a polynomial closed embedding. Let
Pi=(P1,P2) be ambient polynomials with Pi composed with j=K(z).
Write A=C[t1,t2], F=Frac(A), and R=C[ambient variables], with
t_i=P_i. Suppose R has an exhaustive nonnegative filtration with:

1. A fixed A-module basis B with assigned degrees, containing 1 as its sole degree-zero
   basis element. The filtration is the span of basis elements of bounded
   degree and is compatible with multiplication.
2. After passing to an algebraic closure Omega of F, or to
   kappa=Frac(A/(pi)) for every irreducible pi in A, the associated
   graded algebra is a domain with degree-zero part the coefficient field.
3. R tensor_A Omega is a UFD.

Then **K is a polynomial automorphism**, for every polynomial degree
of j, P1 and P2. These are sufficient hypotheses; generic factoriality
alone is not substituted for the height-one fibre condition.

Here and below a closed embedding means that the ambient pullbacks
generate the entire source polynomial algebra. Generic injectivity,
an immersion or a substitution is not enough.

### The global denominator step

Filtered degrees add in each field fibre from condition 2, and units
there have degree zero, hence are constants. If h in R has zero
coefficient on 1 and gcd one among its other A-basis coefficients, then

\[
                         F[h]\cap R=A[h].                 \tag{1}
\]

Indeed localize A at an irreducible pi. Some positive-degree basis
coefficient of h remains nonzero modulo pi. Its reduction has positive
degree, so its powers have distinct degrees and it is transcendental
over kappa. If B(S) in F[S] has negative minimum coefficient valuation,
multiply by the compensating power of pi. Its nonzero reduction cannot
vanish at that transcendental reduction of h. Consequently B(h) retains
the same negative valuation in the free module R localized at pi, and
cannot lie in R. Repeat at every irreducible pi to prove (1).

### The component pencil and the unit derivative

Set L=C(z1,z2) and d=[L:F]. The Keller condition makes this a finite
separable extension; C[z] tensor_A F=L, since this localization is a
finitely generated algebraic F-domain and therefore a field. The smooth
image of j is a prime hypersurface E=0 in the polynomial ambient ring.
Its geometric generic fibre is a disjoint union of d copies of A^r.

By condition 3, E factors into d distinct irreducibles e_i over Omega.
They are conjugate, so have equal filtered degree. Choose an ordering
on B respecting degree, and normalize the coefficient on their common
last nonzero basis element to one. Disjointness makes e_j a unit modulo
e_i, whose quotient is Omega[w]. Hence

\[
                         e_j-c_{ij}=b_{ij}e_i,
                         \qquad c_{ij}\in\Omega^*.        \tag{2}
\]

Degree additivity makes b_ij a scalar, and the normalization makes it
one. Thus all positive-degree basis coefficients agree. Galois descent
gives e_i=h0-a_i, with h0 in R tensor_A F and distinct labels a_i.
Subtract h0's coefficient on 1, clear denominators and divide the gcd
of its other coefficients. This gives a primitive h in R, with
F[h]=F[h0]. Factorization and (1) give the literal global identity

\[
                    E=B(t1,t2,h),\quad B\in A[S],\quad
                    \deg_S B=d.                          \tag{3}
\]

The substitution A[S]->R is injective by positive filtered degree.
A nonconstant factor of B substitutes to a nonunit of the polynomial
C-algebra R: a factor depending on S has positive filtered degree,
and a nonconstant factor in A remains a nonconstant ambient polynomial.
Thus B is irreducible over A and, by Gauss, over F.

On each geometric generic source component h is its constant label.
All positive-w coefficients of h composed with j therefore vanish
identically, giving h composed with j=a(z) in C[z]. The polynomial
B is the minimal relation of a over F and has degree d.

If B_S(K,a) vanished at a source point, differentiate (3) there to get
dE=B_t1 dP1+B_t2 dP2. Pullback to the source and the invertible DK
force B_t1=B_t2=0 at that point. This would make dE=0 at an image
point, contrary to smoothness. Thus B_S(K,a) is a nonzero constant c.
The lower-degree polynomial B_S-c also vanishes at (K,a), so minimality
gives B_S=c identically, B=cS+b(t), and d=1. The retained birational
Keller endpoint gives a polynomial inverse for K.

This repeats the mechanism of [HC4MYHJ0/1](STABLE_HYPERSURFACE_OBSTRUCTION.md)
and [HC4MYHQ0](QUADRIC_OUTPUT_BASE_OBSTRUCTION.md) with their exact
ambient hypotheses isolated. It uses no HC4 Hessian assumption.

## 2. Smooth complete intersections at infinity certify the criterion (HC4MYHV1)

Let P1,P2 in C[Z1,...,Z_(r+3)] have positive degrees d1,d2 and
highest homogeneous forms p1,p2. Suppose their common projective zero
scheme

\[
                    D=V(p1,p2)\subset P^{r+2}             \tag{4}
\]

is a smooth complete intersection of codimension two. If a closed
embedding j as above satisfies (P1,P2) composed with j=K(z), then
**K is invertible**. All d1,d2 and all lower-degree coefficients are
allowed. The condition is scheme-theoretic: smoothness of the reduced
support of a nonreduced intersection is insufficient.

### A common free basis in every target fibre

The forms p1,p2 are a homogeneous regular sequence. The smooth projective
complete intersection D, of dimension r>=2, is connected and integral.
For completeness, connectedness follows by taking H0 in its Koszul
resolution: the intermediate cohomology of line bundles on P^(r+2)
vanishes, leaving H0(O_D)=C. Smooth components would be disjoint, so
connectedness implies irreducibility. The homogeneous ideal is saturated
because its complete-intersection quotient has depth r+1>=3. Thus
C[Z]/(p1,p2) is geometrically a domain, including after field extension.

Give Z degree one and t_i degree zero in

\[
 R=A[Z]/(P1-t1,P2-t2)\simeq C[Z].                        \tag{5}
\]

The leading regular sequence gives

\[
                       gr(R)=A[Z]/(p1,p2).               \tag{6}
\]

Here is the filtered step explicitly. A cancellation between highest
terms a*p1+b*p2 is a first syzygy of this regular sequence, hence has
the form (a,b)=(c*p2,-c*p1). Subtract the corresponding multiple of
the exact relation (P2-t2,-(P1-t1)) from the coefficients. This lowers
their maximal contribution degree. Repeating removes every leading
cancellation and proves that the initial ideal of (P1-t1,P2-t2) is
(p1,p2). The argument works over A and every field fibre, since the
same leading regular sequence remains regular.

Choose any homogeneous C-basis B of C[Z]/(p1,p2), with 1 in degree
zero, and lift it homogeneously to C[Z]. Formula (6) shows by degree
induction that these lifts span R over A. A relation would give a
nonzero highest-degree relation in the free A-module on the chosen
graded basis, impossible. They therefore form precisely the basis
required in Section 1. In particular t1,t2 are algebraically independent.
After every field specialization (6) is a domain, so both the generic
degree argument and the global height-one denominator argument apply,
even when an affine special fibre is singular.

### Generic factoriality from the projective compactification

Over Omega=algebraic closure of C(t1,t2), homogenize the two equations
in (5), using a new coordinate h0. Their projective intersection
Y in P^(r+3) is a complete intersection of dimension r+1>=3.
At h0=0 it meets the hyperplane in exactly D, and the independent
gradients of p1,p2 make Y smooth there.

The affine geometric generic fibre is smooth as well. One can use
generic smoothness in characteristic zero directly: for each irreducible
component of the critical locus of Pi, its image has dimension at most
one, by the rank of the differential at its generic smooth point.
There are finitely many such components, so their images miss the
generic point of A2. Thus Y is smooth everywhere. The same Koszul
connectedness argument makes Y integral.

The classical Grothendieck--Lefschetz Picard theorem now gives

\[
                         Pic(Y)=Z[O_Y(1)].                \tag{7}
\]

This is the explicit new [primary import](sources/smooth-infinity-import.json):
the smooth complete-intersection statement in the introduction of
[Ravindra--Srinivas](https://arxiv.org/html/math/0511134v1), over an
algebraically closed field of characteristic zero. We do not use their
singular/general-member extension or claim to reprove (7).

The boundary D is an integral prime divisor of class O_Y(1).
Restriction Cl(Y)->Cl(Y minus D) is surjective: close every prime
divisor from the open set in Y. Its kernel is generated by D, since
a divisor principal on the open set differs from that principal divisor
only along D. Since Y is smooth, Cl(Y)=Pic(Y), and (7) therefore gives
Cl(Y minus D)=0. The smooth affine ring of Y minus D is normal with
trivial class group, hence a UFD. This supplies condition 3, finishing
the proof without a degree bound or an embedding-rectification theorem.

## 3. Arbitrary-degree examples and the denominator trap

There are two-nonlinear-observable certificates in every common degree
m>=2: in r+3 variables take

\[
          p1=\sum_i Z_i^m,\qquad p2=\sum_i a_i Z_i^m,
                      \quad a_i\ne a_j\ (i\ne j).          \tag{8}
\]

The common projective zero set is nonempty of codimension two. For
example choose a nonzero vector in the kernel of the two-by-(r+3)
matrix with rows (1) and (a_i), then choose m-th roots of its entries.
If the two gradients were dependent at a common zero, distinct a_i
would allow at most one nonzero coordinate. The first equation then
forces that coordinate zero too, impossible projectively. Thus (8)
satisfies (4) for every m, with arbitrary lower-degree corrections.
This certifies ambient observables, not the existence of an embedding
carrying a prescribed nonautomorphic Keller pair.

Generic factoriality by itself misses a real denominator obstruction.
In five ambient variables (x,y,z,u,v) take

\[
               Pi=(x,xy+z^2),\qquad A=C[t1,t2].           \tag{9}
\]

Its geometric generic fibre is affine three-space: eliminate y since
x=t1 is nonzero. The ring has an A-basis consisting of monomials in
y,u,v times 1 or z, using z^2=t2-t1*y. Thus h=z has primitive
nonconstant coefficient one. Nevertheless

\[
                    (h^2-t2)/t1=-y\in R\setminus A[h].    \tag{10}
\]

Over F, z is transcendental and the unique representing polynomial has
coefficient 1/t1, proving the nonmembership. At t1=0 the reduction
of h instead satisfies h^2=t2. The leading scheme V(x,xy+z^2) has
equations x=z^2=0 and is nonreduced; its reduced support is a smooth
projective plane. This is why neither geometric generic factoriality
nor smooth reduced support can replace the conditions above. It is
not a counterexample to Keller invertibility or to the theorem.

For a genuine positive embedding control, the retained linear/quadric
example j(s,t,u,v)=(T,x,y,z,w)=(s,1,t+uv,u,v) and
Pi=(T,xy-zw) has Pi composed with j=(s,t). Its boundary is the smooth
projective quadric V(T,xy-zw). The literal source inverse is
(T,y-zw,z,w). Thus the geometric gate includes actual closed embeddings,
as well as the arbitrary-degree observable family (8).

## 4. Consequence and remaining HC4 boundary

For a hypothetical nonautomorphic unit-Hessian HC4 potential in Tang's
form f=D(z)+w dot K(z), the intrinsic pair K cannot be the restriction
of any ambient polynomial observable pair whose highest forms satisfy
(4), in **any** genuine A4 embedding in A5. The observables need not
include the potential, be linear/quadratic, or be an ambient coordinate
pair. Every proposed extension of K must fail this smooth scheme-theoretic
complete-intersection condition at infinity. The same statement holds
after any certified ambient polynomial automorphism admitting such a
frame for the two observables; no potential-preservation requirement is
needed for this output-base argument.

Arbitrary closed embeddings do supply polynomial extensions of K, but
they need not supply extensions with the required leading forms. Adding
large multiples of an image equation E to both observables commonly
introduces the common factor top(E), violating codimension two. There
is no assertion that extensions can always be improved to satisfy (4).

Images with singular, nonreduced or excess-dimensional observable
intersection at infinity remain outside this new gate. Both plane
endpoints, actual-source unit-Hessian parameterizations in their full
orbits, and unrestricted HC4/JC2 remain open. The
[checker](verify_smooth_infinity.py) and [receipt](receipt-smooth-infinity.json)
give finite algebraic controls, distinct from the written general proof.
No formal verification, independent replay, external review or novelty
priority is asserted.
