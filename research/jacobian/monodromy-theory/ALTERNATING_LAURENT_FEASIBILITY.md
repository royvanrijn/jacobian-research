# Alternating Laurent covers: full-target charts and a normalization-unit obstruction

The natural alternating actions in the Hurwitz construction meet the retained
affine-inertia group tests. Their exact Laurent function-field extensions,
over the chosen affine target, cannot come from any dominant polynomial
self-map of affine space. The obstruction is a nonconstant unit in the finite
normalization. An explicit etale chart nevertheless reaches the entire affine
target, including the branch-value divisors.

Work over an algebraically closed field k of characteristic zero unless a
different field is specified. The
[Hurwitz construction](../../arithmetic/ALTERNATING_HURWITZ_DISCRIMINANT_BRIDGE.md)
supplies, for n=2b>=8,

\[
 f(z)=P(z)+a/z=(H(z)+a)/z,\qquad H=zP,\quad \deg H=n,\quad a\ne0.
\]

There are b distinct finite branch values s_i. Each has two simple ramification
points and n-4 unramified points. The geometric monodromy in its natural
degree-n action is A_n. The two points above infinity have indices n-1 and 1.

## Proposition 1. All the natural actions pass the affine-inertia group tests

For every b>=4 the natural degree-2b action of A_(2b), with the prescribed
double-transposition finite inertia, satisfies:

1. Each nontrivial affine inertia group fixes a sheet.
2. Those inertia groups normally generate the full group.
3. The point stabilizer normally generates the full group; consequently there
   is no nontrivial intermediate Galois field contained in the degree-n field.

Indeed a double transposition fixes exactly n-4 sheets. The literal finite
branch cycles already generate A_n, so their normal closure is A_n.
The natural point stabilizer A_(n-1) contains a double transposition, and all
double transpositions are conjugate in A_n. To check the latter assertion,
start with a conjugator in S_n. If it is odd, multiply it by the transposition
forming either one of the two cycles of the source double transposition;
that odd permutation centralizes the source element and makes the conjugator
even. Since the double-transposition class contains a generating tuple, any
subgroup containing one such element normally generates A_n. This applies
to the point stabilizer. The equivalence with exclusion of an intermediate
Galois field is the Galois correspondence.

More generally, in any faithful transitive action A_n/H with these prescribed
double-transposition inertia classes, the fixed-sheet test holds exactly when
H contains a double transposition. For a generator d,

\[
 d\text{ fixes }gH
 \iff g^{-1}dg\in H.
\]

All d belong to the same class, so one class intersection suffices for all the
prescribed inertia groups. Once it holds, normal generation of H also holds
by the preceding class-generation argument. This is an exact group-test
criterion, not a realization theorem for all such coset actions.

The (n-1)-cycle at infinity also has a fixed sheet in the natural action.
Its fixed sheet is not required by the affine-inertia theorem, which concerns
prime divisors of the chosen affine target.

These statements apply the necessary tests of
[the monodromy theory](THEORY.md), Sections 3--5. They do not classify
polynomial Keller realizations of the abstract natural A_n action.

## Theorem 2. A Jacobian-one chart over the entire affine target

Put

\[
 R(z)=zH'(z)-H(z)-a,\qquad
 U_0=\operatorname{Spec}k[z,z^{-1},R(z)^{-1}].
\]

The morphism

\[
 \Phi:U_0\times\mathbb A^1_u\longrightarrow\mathbb A^2_{s,v},
 \qquad
 (z,u)\longmapsto
 \left(\frac{H(z)+a}{z},\frac{uz^2}{R(z)}\right)
 \tag{1}
\]

is surjective, etale and quasi-finite, has generic degree n and Jacobian
determinant one, and preserves the geometric function-field monodromy A_n.
Its geometric fibers have n points when s is not a branch value, and n-4
points when s=s_i. Its nonproperness locus is exactly

\[
 \bigcup_{i=1}^b \{s=s_i\}\times\mathbb A^1_v.
 \tag{2}
\]

Proof. Since f'(z)=R(z)/z^2, this derivative is a unit on U_0. Therefore

\[
 \det D\Phi=f'(z)\frac{z^2}{R(z)}=1.
\]

More explicitly, changing the fiber coordinate to v=uz^2/R(z) is an
isomorphism U_0 x A1_u = U_0 x A1_v. In these coordinates (1) is just
(f|_(U_0), identity). The curve map has no critical point on U_0, so it is
etale. It is quasi-finite because each fiber is a subset of a fiber of the
nonconstant degree-n map f:P1->P1.

For an ordinary finite target value, all n preimages are finite, nonzero and
unramified. For s_i the two ramification points are excluded from U_0, while
the n-4 unramified points remain. Since n>=8, these fibers are nonempty.
This proves the stated geometric fiber counts and surjectivity. The
function-field extension is k(z,v)/k(f(z),v); it is the original extension
after adjoining the independent variable v, hence has the same geometric
monodromy.

For properness, first use the finite completion

\[
 \overline X=\mathbb G_m\times\mathbb A^1_v
       \longrightarrow\mathbb A^2_{s,v},\qquad
 (z,v)\longmapsto(f(z),v).
 \tag{3}
\]

Its finiteness is proved in Proposition 4 below. The source of (1), in its
v coordinates, is the open subset obtained from (3) by deleting
R(z)=0. This boundary maps exactly onto the lines (2). Outside those lines
the entire finite fiber is retained, so (1) is finite there. Conversely,
for any critical point r above s_i and any v_0, the formal arc
z=r+t, v=v_0 has generic point in U_0 x A1 and maps to an arc through
(s_i,v_0). It cannot extend in that source at t=0 because R(r)=0.
The valuative criterion shows nonproperness at each such target point.
This proves (2).

This chart keeps the whole target A2. In particular the branch-value divisors of its finite completion are present
in the affine target, and its source retains all their unramified branches.
The deleted-target finite-etale chart in the initial bridge is recovered by
also deleting the lines (2) from the target.

## Lemma 3. Units in normalization are necessary for polynomial realization

Let X=Spec A and Y=Spec C be integral affine k-varieties, with A integrally
closed and A^*=k^*. Suppose a dominant generically finite morphism X->Y
embeds C into A. Let B be the integral closure of C in Frac(A). Then

\[
 B\subseteq A,\qquad B^*=k^*.
 \tag{4}
\]

Proof. An element of B satisfies a monic equation with coefficients in C,
therefore with coefficients in A. Since it lies in Frac(A), integral
closedness gives membership in A. A unit of B and its inverse both belong
to A, so it is a unit of A and lies in k^*. The reverse inclusion is immediate.

In particular (4) is a necessary condition for the exact function-field
extension of any dominant generically finite polynomial map A^m->A^m.
It uses neither an etale hypothesis nor a Jacobian condition. It is a
condition on the embedded extension over the full affine target ring,
rather than on its abstract monodromy action.

## Proposition 4. No polynomial realization of the fixed Laurent extensions

Let k be any field. Let H in k[z] have degree n>=2, H(0)=0, nonzero leading
coefficient h, and let a in k^*. Fix the extension

\[
 K=k(s,v_2,\ldots,v_m)\ \subset\
 L=k(z,v_2,\ldots,v_m),\qquad
 s=\frac{H(z)+a}{z},\quad m\ge1.
 \tag{5}
\]

There is no dominant generically finite polynomial map A^m_k->A^m_k whose
induced extension is isomorphic to (5) over the displayed affine target
coordinate ring k[s,v_2,...,v_m].

Proof. Its finite normalization is

\[
 B=k[z,z^{-1},v_2,\ldots,v_m].
 \tag{6}
\]

To verify this directly, the equation H(z)-sz+a=0, divided by h, is monic
in z, so z is integral over k[s]. Also

\[
 z^{-1}=\frac{s-P(z)}{a},\qquad H(z)=zP(z).
 \tag{7}
\]

Consequently k[z,z^-1] is generated as a k[s]-algebra by the single
integral element z, so it is finite over k[s]. The same holds after adjoining
the v_j. The ring (6) is integrally closed: it is a localization of a
polynomial unique-factorization domain. Its fraction field is L. Any element
of L integral over k[s,v_2,...,v_m] is also integral over B, hence belongs
to B; thus (6) is the entire integral closure. Its element z is a
nonconstant unit. Lemma 3 excludes the alleged polynomial map.

Equivalently, suppose a normal affine k-variety X with O(X)^*=k^* admitted
a nonconstant regular function s and a rational function z satisfying this
equation. Monicity first makes z regular. Equation (7) then makes z a unit,
hence z is a nonzero constant. The equation forces s to be constant,
a contradiction.

For a direct elementary polynomial-ring proof of the monicity step, write
z=p/q with p,q relatively prime polynomials. Multiplication of the monic
equation by q^n shows q divides p^n. Unique factorization makes q a
constant. Equation a=z(s-P(z)) then makes z a polynomial unit and hence
constant.

The obstruction survives adjoining any number of independent source and
target variables. It also survives arbitrary polynomial coordinate changes
and any proposed choice of source birational coordinates that makes a
polynomial map over this same affine target: the root z remains a rational
function in the source field, and the monic and unit arguments still apply.

A target birational change that makes s a nonregular rational function
changes the affine target ring. Such changes are not excluded by this
argument. Nor does the argument exclude the natural A_n action as the
monodromy of some other embedded extension. It excludes every polynomial
realization of the specific extensions (5) over their displayed affine
base, with or without the Keller condition.

## Corollary 5. No nontrivial rational one-coordinate factor of a Keller map

Let k be algebraically closed of characteristic zero, and let
F=(F_1,...,F_m):A^m->A^m have nonzero constant Jacobian. If a rational
function z on the source and a nonconstant rational function f in k(T)
satisfy F_1=f(z), then deg f=1.

Proof. F_1 is nonconstant, so z is nonconstant. Choose a Mobius coordinate
w on k(z) that sends one pole of f to infinity. Write

\[
 F_1=\frac{P(w)}{Q(w)},\qquad (P,Q)=1,\qquad\deg P>\deg Q.
\]

The relation P(w)-F_1 Q(w)=0 is monic after dividing by the nonzero leading
coefficient of P, so source normality makes w a polynomial. Bezout gives
A(T)P(T)+B(T)Q(T)=1. Substitution shows

\[
 Q(w)(A(w)F_1+B(w))=1.
\]

Thus Q(w) is a polynomial unit. Since w is nonconstant, this forces Q to
be constant. Therefore F_1=p(w) for a polynomial p. The first row of
DF is p'(w) times the gradient row of w, so p'(w) divides det DF in the
source polynomial ring. It is a unit. Characteristic zero and
nonconstancy of w then force deg p=1, proving the claim.

Consequently a nontrivial rational P1-cover, after adjoining independent
target variables, can never be the exact function-field extension of a
polynomial Keller map over the affine base containing its original cover
coordinate. The two-pole case in Proposition 4 is stronger: it excludes
every polynomial map, without using the Jacobian. The corollary still
does not exclude an abstract permutation action from occurring through
a different embedded extension.

## An explicit rational chart

The [rational degree-eight model](../../arithmetic/RATIONAL_ALTERNATING_PENCIL_REALIZATION.md)
has

\[
 H(z)=z^8+z^6+\frac9{28}z^4+\frac{233}{5292}z^2,\qquad
 a=-\frac1{5488}.
\]

Its four finite inertia generators are double transpositions and its
geometric group is A8. Formula (1) is therefore defined over Q, with

\[
 R(z)=7z^8+5z^6+\frac{27}{28}z^4+\frac{233}{5292}z^2+\frac1{5488},
\]
\[
 \Phi(z,u)=
 \left(
 z^7+z^5+\frac9{28}z^3+\frac{233}{5292}z-\frac1{5488z},
 \frac{uz^2}{R(z)}
 \right).
\]

Over an algebraic closure it has eight points in a general fiber and four
points over each of the four nonproperness lines. It reaches the full
affine target, has Jacobian one, and has geometric A8. Its source is
Spec Q[z,z^-1,R^-1,u], which contains the nonconstant unit z. The entire
embedded extension is excluded from polynomial affine-space realization
by Proposition 4, not only from extending this particular displayed chart.

## Scope and verification

The arguments above are written proofs. They require no finite search,
new symbolic certificate, or new group classification. The cover and
rational-model inputs retain their own exact proof and checker records.
The cotangent determinant and normalization identities have also been
checked directly in the displayed formulas.

All-b generation is imported only from the literal alternating tuple
theorem already retained in the bridge. The fixed-sheet and normal-generation
necessary conditions are those in the existing Keller monodromy theory.
No independent mathematical replay, formal verification, external review
or novelty is claimed. A different embedded extension now realizes every
even natural alternating action A_n, n>=6, by the polynomial
[cyclic-cone maps](UNIFORM_EVEN_ALTERNATING_CONE_KELLER_MAPS.md) in dimension
three. Their varying pole coefficient and invariant-product target avoid
the fixed affine-base hypotheses used above. The plane question remains open.
