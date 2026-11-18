# Alternating Hurwitz tuples and square-discriminant tangent pencils

This note proves an all-degree bridge from the explicit alternating tuples to
Laurent covers and sharp collinear-node configurations. Exact degree-eight
algebraic and rational coefficient models accompany it. The rational examples
have geometric A8 and arithmetic S8; the literal rational-square obstruction
remains valid. The underlying passport existence is already implied by prior
theorems, as documented below. No polynomial affine-space Keller realization or
priority for the explicit formulas is asserted.

The tuple and generation input is
[the alternating-target proof](../meridional-rank/representations/ALTERNATING_TARGET_PRIME_FAMILIES.tex),
equation generators and Theorem targets. The geometric comparison is with
[the tangent-pencil manuscript](../../unpublished/discriminant-pencils/main.tex),
Proposition normalization. Its global generic cusp-node theorem is not needed
here: the degree and normalization argument is repeated below.

## Lemma 1. The ordered product

Use permutations of {0,...,2b-1}, b>=4, and products pq=p composed with q. Put

\[
x_0=(0\,1)(2\,3),\quad x_1=(1\,4)(2\,3),\quad
x_2=(0\,5)(1\,6),\quad x_3=(1\,2)(6\,7),
\qquad x_j=(1\,2j)(6\,2j+1)\quad(j\ge4).
\]

They generate A_(2b), by the input theorem. Their ordered product p_b has the
following single nontrivial cycle, with 3 fixed:

\[
p_b=(0\ 5\ 1\
 [2b-2,\ 2b-4,\ldots,8]\
 2\ 6\
 [2b-1,\ 2b-3,\ldots,9]\
 7\ 4).
\]

Both bracketed blocks are empty for b=4. Indeed direct multiplication gives
p_4=(0 5 1 2 6 7 4)(3). On passing from b to b+1 extend p_b by fixing the
new symbols 2b,2b+1. In p_b x_b the arrow out of 1 becomes
1 -> 2b -> p_b(1), and the arrow out of 6 becomes
6 -> 2b+1 -> p_b(6). All other arrows are unchanged. Induction proves the
formula and cycle type (2b-1)(1) for every b>=4.

The tuple (x_0,...,x_(b-1),p_b^-1) generates A_(2b) and has product one.
Its finite entries have index two; the last entry has index 2b-2.

## Theorem 2. Covers and square discriminants in every even degree at least eight

For every b>=4, put n=2b. There are a number field L contained in C,
a polynomial H in L[z] of degree n with H(0)=0, a nonzero a in L,
and distinct algebraic numbers s_1,...,s_b such that

\[
 f(z)=\frac{H(z)+a}{z},\qquad
 \operatorname{Mon}_{\overline{\mathbb Q}}(f)=A_n,
\]
\[
 \Delta_H(s,a):=\operatorname{Disc}_z(H(z)-sz+a)
   =C\prod_{i=1}^{b}(s-s_i)^2,\qquad C\ne0.
\]

There are precisely two simple critical points over each finite branch value
s_i and no other finite ramification. Above infinity the two points have
indices n-1 and 1. The discriminant is a square in Qbar[s], and also in
L(sqrt(C))[s]; it need not be a square over L.

Proof. Choose distinct finite s_i in Qbar and prescribe the branch tuple from
Lemma 1, with its last entry at infinity. The product-one relation defines a
permutation representation of the punctured sphere group. Transitivity makes
the corresponding topological covering connected. Fill the punctures using
the cycle lengths of their local monodromy. Riemann existence realizes this
cover algebraically; algebraic branch points give a number-field model of the
same cover. The precise imported statements are Bilu--Strambi,
Theorems 1.1 and 1.2(a), pages 2--3 of arXiv:0809.0345v4, recorded in the
[source ledger](alternating_hurwitz_source_ledger.json).

Riemann--Hurwitz gives
2g-2=-2n+2b+(n-2)=-2, so the source has genus zero. The two points above
infinity have different indices and are individually rational over a field of
definition (enlarge that number field if needed). A genus-zero curve with a
rational point is P1. Send the pole of order n-1 to source infinity and the
simple pole to source zero. The function is P_(n-1)(z)+a/z with a nonzero.
Set H=zP_(n-1). Branch cycles give full monodromy A_n and precisely the
asserted simple critical points. Each pair contributes discriminant order two.
Lemma 3 below computes the nonzero leading coefficient and degree n, completing
the displayed factorization. No rational descent was used.

## Lemma 3. Exact discriminant factorization and leading coefficient

Let n=2b, H(0)=0, h be the leading coefficient of H, and a!=0.
Let r_1,...,r_n, counted with multiplicity, be the roots of

\[
 R(z)=zH'(z)-H(z)-a.
\]

Then, with the usual root-product convention for the discriminant,

\[
 \operatorname{Disc}_z(H(z)-sz+a)
 =(-1)^{b+1}(n-1)^{n-1}h^{n-2}
   \prod_{j=1}^{n}(s-H'(r_j)).
\tag{1}
\]

In particular its degree in s is n, with leading coefficient
(-1)^(b+1)(n-1)^(n-1)h^(n-2).

Proof. Write F=H-sz+a and note R=zF'-F. At an F-root alpha,
R(alpha)=alpha F'(alpha). Since n is even, the product of its roots is a/h.
The definitions of the two resultants give

\[
 \operatorname{Res}(F,R)
   =a\operatorname{Res}(F,F')
   =(-1)^b a h\,\operatorname{Disc}(F).
\]

On the other hand, F(r)=r(H'(r)-s) at an R-root and
product r_j=-a/((n-1)h). Thus

\[
 \operatorname{Res}(R,F)
 =-a((n-1)h)^{n-1}\prod_j(s-H'(r_j)).
\]

The two resultants agree because n is even. Division proves (1).
These equalities include multiplicities and are polynomial identities.
For arbitrary H(0), replace a by a+H(0) after subtracting H(0) from H;
the leading-coefficient formula also holds when that sum is zero, by the
polynomial identity in the constant coefficient.

## Theorem 4. The maximum number of collinear ordinary nodes

For H,a furnished by Theorem 2, the reduced projective discriminant curve
D_H of H(z)-sz+t has degree n. The line t=a contains exactly b ordinary
nodes of D_H and no other intersections. It meets each of the two local
branches at each node transversely. No line can contain more than
floor(n/2)=b ordinary nodes of D_H.

Proof. The affine discriminant has parametrization

\[
 \nu(r)=(H'(r),rH'(r)-H(r))=(S(r),T(r)).
\]

Its projective extension is given by the degree-n homogeneous polynomials

\[
 S_h=\sum_{j=1}^{n}jh_j R^{j-1}U^{n-j+1},\quad
 T_h=\sum_{j=0}^{n}(j-1)h_j R^jU^{n-j},\quad Z=U^n.
\]

They have no common zero, including at U=0 since h_n!=0.
At a generic smooth point,
dT/dS=r because (S',T')=H''(r)(1,r). In characteristic zero this recovers r
in the function field of the image. The map is birational and finite, hence
is the normalization. Pullback of a general line has degree n, proving that
the integral image curve has degree n. Its only point at infinity is
[0:1:0]. With q=U/R the coordinate S_h/T_h has linear term n q/(n-1);
this point is smooth. It is outside the projective line T=aZ.

For each branch value s_i, the polynomial H-s_i z+a has exactly two
distinct double roots r_i,u_i and all its remaining roots are simple.
These two roots are nonzero since a!=0. Equivalently they are the two
points with S=s_i and T=a. At a critical point of f,

\[
 f'(r)=0,\qquad f''(r)=H''(r)/r.
\]

Simple ramification therefore implies H''(r_i)H''(u_i)!=0.
The two branches of the normalized discriminant at (s_i,a) are smooth,
with tangent vectors (1,r_i) and (1,u_i). Their determinant is
u_i-r_i!=0. They are exactly the two branches there, so the singularity is
an ordinary node. Since both r_i,u_i are nonzero, the horizontal line is
transverse to each branch and has local intersection number two.

The b nodes already contribute 2b=n intersections. The line is not a
component of the integral degree-n curve, so Bezout exhausts its
intersection divisor. Conversely every ordinary node on any line consumes
at least two of the n intersections. This proves the bound and its
attainment. The argument makes no claim about singularities away from the
chosen line or membership in the manuscript's globally generic locus.

## Corollary 5. Rational degree-eight square-discriminant descent is impossible

For H in Q[z] of degree eight and a in Q, the polynomial
Disc_z(H-sz+a) is never a square in Q[s]. Its s^8 coefficient is
-7^7 h_8^6, a negative rational number.

If a rational-coefficient member of Theorem 2 exists, it necessarily has

\[
 \Delta_H(s,a)=-7 B(s)^2,\qquad B\in\mathbb Q[s].
\]

Its geometric group is A8, its arithmetic group is S8, and the algebraic
constant field of its Galois closure is Q(sqrt(-7)).

Proof. The first assertion follows from Lemma 3, since a nonzero rational
square has positive leading coefficient. For the second, normalize the
degree-four square root of the monic discriminant factor to be monic.
Its uniqueness and Galois invariance show that it lies in Q[s]. The leading
coefficient is -7(7^3 h_8^3)^2. The arithmetic group contains the geometric
A8; the nonsquare discriminant forces an odd permutation, hence S8.
The quotient S8/A8 has order two and is detected by sqrt(Delta), giving
exactly the indicated constant field.

More generally, in degree n=2b the rational square class of the leading
coefficient is (-1)^(b+1)(2b-1). A necessary condition for a literal
Q[s]-square is that this be a rational square. This holds exactly when
2b-1 is an odd integer square. It is only a necessary condition.
This obstruction does not disprove rational coefficients with geometric
A8 and discriminant a nonsquare constant times a square.

## Further consequences and limits

The standard (n-1)-dimensional complex representation of A_n has no common
fixed vector. Each double transposition has rank(x-I)=2. If k such elements
generate, the intersection of their fixed spaces is zero, while its
codimension is at most 2k. Thus k>=ceil((n-1)/2)=b; the tuple attains equality.
The ordinary generating rank two and the class rank b are already part of
the alternating-target argument. Since A_n is perfect, every complex
representation has determinant one. A nonidentity finite-order rank-one
reflection is diagonalizable and has determinant its unique
eigenvalue different from one. Thus no representation of A_n contains such a reflection in
its image. This statement concerns representations of the alternating
quotient itself.

For a finite etale chart set
V=P1 minus {s_1,...,s_b,infinity} and U=f^-1(V).
On U the derivative f' is a unit and

\[
 U\times\mathbb A^1\longrightarrow V\times\mathbb A^1,\qquad
 (z,u)\longmapsto(f(z),u/f'(z))
\]

is finite etale of degree n with Jacobian determinant one.
Its function-field extension is the original extension after adding the
independent target coordinate v, since u=v f'(z). Its geometric monodromy
is A_n. The finite inertia double transpositions have n-4 fixed sheets
and generate A_n; the point stabilizer A_(n-1) normally generates A_n.
These are boundary inertia groups of the completed cover; all divisor inertia
inside the finite etale target chart is trivial. The abstract action and this
proposed inertia collection pass the fixed-sheet and normal-generation necessary
group tests in
[the Keller monodromy theory](../jacobian/monodromy-theory/THEORY.md),
Sections 3--5. The source and target are affine charts with deleted
divisors. Their nonconstant units persist after adjoining variables;
the stronger [full-target chart and normalization-unit argument](../jacobian/monodromy-theory/ALTERNATING_LAURENT_FEASIBILITY.md)
show that the same extensions have surjective etale Jacobian-one charts over
the entire affine target. Their finite normalizations contain the nonconstant
unit z, which excludes every polynomial affine-space realization of these
specific extensions over the retained affine base, even without a Jacobian
condition. The abstract alternating action through another extension remains
a separate realization question.

The braid W_b in the alternating-target construction fixes the ordered
finite tuple under its stated Hurwitz convention. Embed B_b in B_(b+1)
by leaving the final branch point fixed. It fixes the completed tuple,
since Hurwitz moves preserve the product. This supplies a stabilizer in
the full braid action on the appropriate Nielsen class. The corresponding
loop is in the configuration space of unordered finite branch points;
for ordered branch points take a pure stabilizing power of W_b.
This uses the usual Hurwitz-space interpretation and makes no new
correspondence claim.

## Explicit rational models and the degree-eight braid orbit

There are now [infinitely many rational-coefficient examples](RATIONAL_ALTERNATING_PENCIL_REALIZATION.md).
A simple member is

\[
 H(z)=z^8+z^6+\frac9{28}z^4+\frac{233}{5292}z^2,\qquad
 a=-\frac1{5488},
\]
\[
 \Delta_H(s,a)=-7^7
 \left(s^4+\frac{5960}{66706983}s^2+
                    \frac{2685619}{1483273860320763}\right)^2.
\]

Its four finite branch values have double-transposition inertia, its geometric
group is A8, its arithmetic group is S8, and its constant field is
Q(sqrt(-7)). The two-parameter coefficient family and infinite pairwise
nonisomorphic full-A8 subfamily are proved in the linked note. These examples
settle rational-coefficient geometric realization positively. They are
consistent with Corollary 5, which excludes a literal square over Q.

The earlier [quadratic-field model](DEGREE_EIGHT_ALTERNATING_PENCIL.md),
with e=(-1009+94sqrt(-329))/10584, retains its separate no-rational-form
statement: its coefficient invariant b_1 b_7^2/b_5^3=e is nonrational.

The [complete degree-eight Hurwitz computation](DEGREE_EIGHT_HURWITZ_ORBIT.md)
finds 93,296 double-transposition quadruples with a fixed seven-cycle product
and exactly one full-A8 braid orbit, of size 47,334. Consequently both the
rational models and the earlier quadratic-field model lie in the knot tuple's
absolute braid orbit. This is deformation equivalence of branch-cycle data,
not isomorphism of individual covers.

## Priority and remaining questions

Bilu--Strambi, *Quantitative Riemann Existence Theorem over a Number Field*,
[arXiv:0809.0345v4](https://arxiv.org/pdf/0809.0345v4), Theorems 1.1 and
1.2(a), pages 2--3, supply complex realization and a number-field model of
the same topological cover for algebraic branch points. They are imported
theorems, not independently reproved here.

The [exact primary-source comparison](ALTERNATING_HURWITZ_PRIORITY.md) identifies
the underlying passport as a direct consequence of Edmonds--Kulkarni--Stong,
Proposition 5.3 (1984), and Pakovich, Theorem 1.1 (2009). For even degrees
at least ten, classical primitive-group results also force full alternating
monodromy. Thus broad cover existence is established precedent, and the
maximal-collinear-node existence follows geometrically from that passport.
The elementary geometric interpretation, explicit knot tuple and stabilizer,
and exact coefficient certificates supply concrete connections; no new
abstract existence or worldwide priority claim is made.

The following questions remain distinct: rational-coefficient realizations in
every even degree, polynomial Keller realization of the abstract alternating
action through a different embedded extension, and exact historical priority
for the displayed rational family and explicit stabilizers. The specific
Laurent extensions are excluded from polynomial realization over their
displayed affine bases by the normalization-unit proof.

The [tuple checker](../meridional-rank/representations/verify_alternating_tuple_product.py)
and its [receipt](../meridional-rank/representations/alternating_tuple_product_verification.json)
supply bounded controls for b=4,...,64; the insertion proof establishes the
all-b theorem. The rational-family checker proves the displayed polynomial
identities in all its parameters, while the braid checker exhausts its entire
finite degree-eight domain. These are different evidence scopes.
Independent external replay, formal verification and external human review
are not claimed.
