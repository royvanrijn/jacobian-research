# Divisor-class monodromy replaces generic factoriality

14 November 2026. The generic ambient fibre need not be factorial.
A torsion-free class group with abelian Galois action suffices, because
the function field of a Keller map has no nontrivial abelian Galois
intermediate extension. We also replace the height-one graded-domain
condition by geometric integrality of those affine fibres. These two
changes close the previously excluded quadric-surface case, permit
degenerate ambient quadratics with sufficiently large restricted rank,
and cover cyclic Danielewski observable fibres in arbitrary degree.

## 1. No abelian intermediate extension (HC4MYHW0)

Let K:A^n_z->A^n_t be a complex Keller map in any positive dimension.
Put F=C(K1,...,Kn) and L=C(z1,...,zn). If

\[
                 F\subset E\subset L,\qquad E/F
                 \text{ finite abelian Galois},            \tag{1}
\]

then E=F. Degrees of K are unrestricted. This statement concerns
abelian intermediate extensions; it does not say L/F is Galois or
trivial, nor settle any Jacobian conjecture.

First, every source prime divisor maps dominantly to a target prime
divisor. At its generic smooth point the invertible DK has rank n-1
on its tangent space, so the closure of its image has dimension n-1.
Equivalently use that a Keller morphism is etale and quasi-finite.
Consequently two distinct target irreducibles have relatively prime
pullbacks: a common source divisor could not map into their
codimension-at-least-two intersection.

For an irreducible target polynomial b, b(K) is nonconstant and
square-free. If a source irreducible g occurred twice, differentiating
would give

\[
           DK^T(\nabla b)(K)=0\pmod g,
                 \qquad (\nabla b)(K)=0\pmod g.             \tag{2}
\]

The source divisor would then map into the singular locus of b=0,
which has codimension at least two in characteristic zero. This
contradicts the preceding dimension statement. The same conclusion
in dimension one is immediate. The needed forward square-free fact
is classical; it is reproduced here, with
[Jedrzejewicz's Theorem 5.1](https://arxiv.org/html/1109.2113v2)
retained as [primary context](sources/abelian-classgroup-context.json).
We do not import the converse characterization.

Suppose (1) were nontrivial. A finite nontrivial abelian group has a
cyclic quotient of prime order p, giving a degree-p cyclic subextension.
Its order-p automorphism is diagonalizable on that field as an
F-vector space, since F contains all p-th roots of unity. A nontrivial
eigenvector theta satisfies theta^p=a in F, with a not a p-th power.
Indeed an eigenvalue other than one exists, and a p-th-power a would
make theta an F-multiple of a constant root of unity.

Factor a in the rational function field of the target. After multiplying
theta by an element of F and extracting a constant p-th root, write

\[
                         a=\prod_i b_i^{e_i},\qquad
                         0<e_i<p,                         \tag{3}
\]

with distinct irreducible b_i and at least one factor. Each b_i(K)
has a source prime factor, with multiplicity one, and the different
pullbacks have no common factor. The valuation of a(K) at such a
factor is e_i. It cannot be divisible by p, whereas every valuation
of theta^p in the source rational function field is divisible by p.
This contradiction proves HC4MYHW0. No surjectivity, properness,
purity theorem or general Galois-case theorem was used.

## 2. A class-group and geometric-fibre output criterion (HC4MYHW1)

Let K be plane Keller, r>=0, and let j:A^(r+2)_(z,w)->A^(r+3)
be a polynomial closed embedding. Suppose two ambient polynomials
Pi=(P1,P2) satisfy Pi composed with j=K(z). Set A=C[t1,t2],
F=Frac(A), Omega=algebraic closure of F, and R=C[ambient variables]
with t_i=P_i. Assume:

1. R has a free A-basis with assigned nonnegative degrees, with 1 its
   only degree-zero basis element. Their bounded-degree spans form an
   exhaustive multiplicative filtration. The associated graded ring
   of R_Omega is a domain with degree-zero part Omega.
2. For every irreducible pi in A, R tensor_A Frac(A/(pi)) is
   geometrically integral. These fibres may have nonconstant units
   or a nondomain associated graded ring.
3. R_Omega is normal; its Weil divisor class group is torsion-free,
   and the image of the absolute Galois action on that group is abelian.

Then **K is a polynomial automorphism**, in arbitrary source and
observable degrees. The third condition includes factoriality as
class group zero, but also allows nonfactorial generic fibres.

### Why the geometric image components become principal

Write L=C(z1,z2) and d=[L:F]. The prime smooth hypersurface image
E=0 has geometric generic fibre a disjoint union of d copies of A^r.
Let D0 be the prime divisor corresponding to the chosen embedding
L into Omega, and alpha its divisor class in Cl(R_Omega).
Its divisor stabilizer is H=Gal(Omega/L), hence its class stabilizer
S contains H. Since the action on the class group has abelian image,
S is normal and G/S is abelian, where G=Gal(Omega/F). It has finite
index because H does. Thus Omega^S/F is finite abelian Galois and
is contained in L. HC4MYHW0 forces S=G.

All d conjugate component divisors therefore have the same class
alpha. The image equation is reduced on this fibre, so

\[
                       0=[div(E)]=d\,alpha.                \tag{4}
\]

Torsion-freeness gives alpha=0. Each component is principal. More
explicitly a rational function with divisor that single prime divisor
has nonnegative valuations at every height-one prime; normality puts
it in the ring, and its ideal is exactly the component prime ideal.
The residual factor in E is a unit and hence constant by condition 1.

The equal-degree disjoint-component pencil argument of
[HC4MYHV0](SMOOTH_INFINITY_OUTPUT_OBSTRUCTION.md#1-a-filtered-output-base-criterion-hc4myhv0)
now applies: choose normalized generators e_i, deduce
e_i=h0-a_i, and descend their common nonconstant part to F.
No assertion that every divisor of R_Omega is principal was needed.

### Why geometrically integral height-one fibres suffice

Make the descended element primitive in the fixed A-basis and
subtract its coefficient on 1, obtaining h in R. At every target
prime pi at least one nonconstant basis coefficient survives, so
the reduction of h is not in kappa=Frac(A/(pi)). A geometrically
integral kappa-algebra in characteristic zero has no element algebraic
over kappa outside kappa: otherwise extension to an algebraic closure
would split a nontrivial separable minimal polynomial inside a domain.
Thus this reduction of h is transcendental even if its filtered degree
does not behave additively. The same height-one valuation argument
gives the global identity F[h] intersection R=A[h].

Consequently the image equation is E=B(t,h) globally, with B irreducible
of label degree d. Smoothness and the invertible DK make B_S(K,a)
a nonzero constant. The minimal-relation argument in HC4MYHV0 forces
d=1, and the retained birational Keller endpoint gives the inverse.
The rest of that proof uses only generic degree additivity, preserved
by condition 1. This proves HC4MYHW1, without a source Hessian condition.

The two new hypotheses have distinct roles. Torsion in the class group
could permit d*alpha=0 without alpha=0. A nonabelian action can have
a nonnormal stabilizer and a nonzero orbit with sum zero, without
giving an abelian intermediate field. Neither case is silently included.

## 3. Quadric surfaces and degenerate quadratics are covered (HC4MYHW2)

In A^(r+3), r>=1, allow additional unused affine variables and either
of the following observable pairs:

\[
 \begin{array}{ll}
 \mathrm{I}: &(T,q_m(D)),\quad m>=3,\\
 \mathrm{II}:&(T,T Y+q_m(D)),\quad m>=2,
 \end{array}                                               \tag{5}
\]

where q_m is nondegenerate and the displayed variables fit in the
ambient dimension. Every genuine closed embedding carrying a plane
Keller pair through (5) forces that pair invertible, in all embedding
degrees. We can subtract any polynomial in T from the second target
component without changing this conclusion.

In type I give D and unused variables degree one. Use the monic
hyperbolic relation D1*D2 plus a quadratic in the other variables
minus t2; its standard monomials, not divisible by D1*D2, give the
common free basis. The generic associated graded ring has an irreducible
quadratic of rank m>=3. Every height-one affine fibre q_m=t2 remains
geometrically integral, including t2=0, by the same rank bound.

For m>=4 the generic fibre is factorial by the retained Nagata argument.
For m=3, after extending to Omega, its ring has the form

\[
       B=Omega[x,y,z,U]/(xy-(z-a)(z+a)),\qquad a\ne0,       \tag{6}
\]

where U denotes unused affine variables. It is smooth but not factorial.
Localizing at x gives the UFD Omega[x,x inverse,z,U]. Hence its class
group is generated by the two divisors D_plus=(x,z-a) and
D_minus=(x,z+a). The only relation is D_plus+D_minus=div(x).
To check there is no further relation, a rational function with divisor
n_plus*D_plus+n_minus*D_minus is a unit after localizing x, hence
c*x^k. Its divisor has n_plus=n_minus=k. Thus Cl(B)=Z.
This localization argument is just restriction and closure of prime
divisors as in HC4MYHV1; it imports no class-group classification.
The Galois action on Z has image in {1,-1}, which is abelian. In the
unsplit model the involution exchanging the roots a and -a exchanges
the two divisor classes. HC4MYHW1 applies despite nonfactoriality.

For type II give Y degree two and D and unused variables degree one.
The same monic hyperbolic relation gives a common A-basis. The generic
graded relation is T*Y+q_m(D); T is a nonzero field element, so Y
eliminates and this graded ring is a domain. The generic affine ring
is itself polynomial, hence normal and factorial. At every target
prime other than T, Y still eliminates. At T=0 the residue field is
C(t2), and the ring has q_m(D)=t2. It is geometrically integral even
for m=2: the split binary equation is xy=t2 with t2 nonzero. This
fibre has nonconstant units, and its associated graded ring can fail
to be a domain. The geometric-integrality version of the denominator
argument is exactly what includes it.

### Coordinate-free consequences

Let ell be a nonconstant affine-linear observable and Q any ambient
quadratic on A^(r+3), r>=1. If the Hessian of Q restricted to ker(dell)
has rank at least three, every closed embedding carrying a plane
Keller pair through (ell,Q) forces that pair invertible. Q may be
degenerate on the whole ambient space.

Indeed put T=ell and split that restricted quadratic into a nondegenerate
q_m(D) and radical variables Y. Completing squares removes all T-linear
and constant-linear terms in D. The remaining expression is

\[
                  q_m(D)+a(T)\cdot Y+c(T),                \tag{7}
\]

where the entries of a(T) are affine-linear in T. Over C[T] a
unimodular change of the Y variables reduces this row to its gcd g(T).
If it is zero we are in type I; if a unit, (T,Q) is an ambient
coordinate pair and HC4MYHJ1 applies; otherwise g has degree one,
and an affine change of T gives type II. These are literal polynomial
ambient changes with polynomial inverses. With no radical variables
only type I occurs.

In particular every nondegenerate ambient Q is now covered for all
r>=1. If its normal is nonisotropic, the restricted rank is r+2>=3;
if isotropic it is r+1>=2, with the required nonzero T*Y term of
type II. Thus the earlier restriction r>=2 in HC4MYHQ0 is removed.
For HC4 (r=2), degenerate quadratics with restricted rank at least
three are now additional forbidden ambient extensions of a hypothetical
nonautomorphic intrinsic cotangent pair. The explicit type-II rank-two
class is covered too. The observables need not include the potential.

## 4. Cyclic Danielewski observable fibres in every degree (HC4MYHW3)

For every integer m>=2, both following pairs force a carried plane
Keller pair invertible in a genuine stable hypersurface embedding:

\[
 \begin{array}{ll}
 \mathrm{III}:&(T,xy-z^m),\quad r>=1,\\
 \mathrm{IV}:&(T,T Y+xy-z^m),\quad r>=2.
 \end{array}                                               \tag{8}
\]

Allow unused affine variables to reach A^(r+3). These observables
are often called Danielewski equations. No classification of those
surfaces or their affine-line embeddings is imported.

Assign x,y degree m, z degree two, and unused variables degree one.
In type IV give Y degree 2m. The monic relation xy=z^m+t2 (or
xy=z^m+t2-T*Y) again gives a common free A-basis consisting of
monomials not divisible by xy. In type III the generic associated
graded ring has equation xy=z^m, a domain. In type IV it has equation
T*Y+xy=z^m, also a domain since T is generically invertible.

All type-III height-one affine fibres xy=z^m+t2 are geometrically
integral, including t2=0: the polynomial is primitive linear in y,
and x has no common factor with z^m+t2. In type IV, eliminate Y
unless the target prime is T; at T=0 the same integral equation
remains over C(t2). Thus the global denominator criterion applies.

For type III the geometric generic ring is smooth, since its gradient
could vanish only at x=y=z=0, incompatible with t2 nonzero. Its m
roots a_j of z^m=-t2 are distinct. Localizing x gives a Laurent
polynomial ring, so the identical divisor calculation from (6) gives

\[
                Cl(R_{Omega})=Z^m/Z(1,\ldots,1),           \tag{9}
\]

generated by (x,z-a_j), with the sole relation div(x)=sum_j D_j.
There is no torsion: the vector (1,...,1) is primitive. Since C
contains every m-th root of unity, adjoining one root of -t2 splits
the polynomial, with cyclic Galois group of order m (Eisenstein at
t2 also proves the degree). The action on (9) therefore has abelian
image. HC4MYHW1 proves the assertion despite the generic ring being
nonfactorial. Type IV has polynomial geometric generic ring after
eliminating Y, so it is even factorial.

For m>2 the ordinary highest forms in type III are T and -z^m.
Their projective leading scheme is nonreduced. Type IV has the same
defect for m>2. Thus (8) lies genuinely outside the smooth ordinary
infinity certificate of HC4MYHV1 in the displayed frame. Moreover,
for r>=2 type III has nonzero geometric generic class group in every
ambient frame, whereas that older certificate would force generic
factoriality. Thus no ambient automorphism can put type III into that
older certificate. The supplied weighted filtration
and abelian class action give a different certificate; no good ordinary
infinity frame is assumed.

There are genuine embeddings in every degree, for example
j(s,t,w)=(s,1,t+w^m,w) for type III, with inverse (T,y-z^m,z).
Adding an unused coordinate v gives HC4-size embeddings. The independent
nondegenerate ambient potential T*z+y*v+x^2 then pulls back to
1+s*w+t*v+w^m*v, whose actual Hessian determinant is one in all m.
This is a positive compatibility example with an invertible gradient,
not a construction of a nonautomorphic Keller pair.

## 5. Boundaries and negative controls

The denominator trap (T,T*Y+z^2) remains excluded: at T=0 the
geometric fibre z^2=t2 splits, and h=z becomes algebraic. Its
identity (h^2-t2)/T=-Y still prevents primitive saturation. Generic
factoriality has not become sufficient by itself.

For r=0 and type-I binary quadratics, the generic ambient curve xy=t2
has nonconstant units. The nonnegative additive degree argument fails,
so no extension to that case is claimed. Nor are all low restricted
ranks, nonabelian class actions, torsion class groups, or arbitrary
singular infinity frames classified. There is no assertion that arbitrary
observable extensions admit frames satisfying this criterion.

The [checker](verify_abelian_classgroup.py) and
[receipt](receipt-abelian-classgroup.json) retain exact algebraic controls,
including nonabelian and torsion orbit countercontrols. Finite group
examples do not prove the general Galois argument. Both plane endpoints
and unrestricted HC4/JC2 remain open. No independent replay, formal
verification, external review or novelty priority is asserted.
