# Cyclic quotient orders, conductor balance and the scalar-Morse residue budget

PF2RULED106, 16 November 2026. Work over C with the ACTUAL normal finite
Keller normalization N, its WHOLE marked U=A2, and all omitted primes.
The field degree is d, a boundary transverse index is e, its residue
degree is f, and a finite cyclic quotient order is n. These integers
are distinct data.

**Unconditional all-index quotient theorem.** Every finite cyclic quotient
point on an actual index-e omitted prime has, as a MARKED quotient with
the lifted boundary X=0, type

    1/n(1,-e),  gcd(n,e)=1,  2<=n<=e-1.              (1)

In particular n<=d-2. At the largest possible order n=e-1 the point
is A_(n-1), its boundary residue map is unramified and hence f=1, and
its image branch has one Puiseux pair (a,n*a-1) and local surface degree
n*a for some a>=2. No single-prime or finite-image smoothness hypothesis
is imposed.

**All-degree residue corollaries.** Every omitted prime with

    e in {1,2,3,4,6}

has f=1. Consequently ALL omitted residue degrees are one whenever
d<=10, without an all-odd-index or smooth-N hypothesis. This is not an
automorphism theorem for all maps of those field degrees.

**All-index scalar-Morse theorem.** At a finite quotient point suppose
SOME scalar linear combination of the two target coordinates, pulled
back to its C2 uniformization, has a nondegenerate quadratic term.
This is the scalar-Morse condition. Then, putting k=e-1,

    n*f divides k,
    the quotient is A_(n-1),
    ell=k*f*a,  eta=f*(k*(a-1)-1),  a>=2,            (2)
    d>=2*k*f+1>=4*f^2+1.                            (3)

Here a is the multiplicity of the chosen image branch and ell is the
ambient local degree. All other omitted primes and finite excesses
remain in the global budget. In particular any actual residue-f>1
point with d<=4*f^2 must FAIL the scalar-Morse condition. A failure
has both quadratic forms zero or multiples of one linear square;
for f>1 that square, if present, is the normal boundary coordinate
squared. The scalar-Morse condition is NOT proved automatic.

The exact all-index conductor identity in Section4 also applies outside
the scalar-Morse sector. It keeps the ordinary unramified companion
curves and does not replace them with omitted primes.

## 1. Retained geometry and equivariant coordinates

Use the original affine target, not a blowup. The primary boundary and
normalized power-map LOCATION inputs are those fully audited in
[PF105](TRIPLE_BOUNDARY_RESIDUE_ONE_AND_NORMAL_ITERATED_EXCLUSION.md).
Every boundary normalization is A1, and any nontrivial residue map is a
power map ramifying at its finite cyclic quotient point; its local
degree there equals the global f. Without a quotient point it is an
isomorphism. These published facts are attributed to Borisov, including
the stronger smooth-boundary consequence already used by HC4.

At a singular boundary point let rho:C2->N have cyclic degree n>1.
The all-index endpoint-complement input of PF2OMIT7 makes the reduced
inverse image A of the boundary ONE smooth branch. This is essential;
we do not assert it for arbitrary curves through arbitrary quotients.
The group is small and free away from zero, so both eigencharacters
of a generator are primitive n-th roots.

There are coordinates simultaneously straightening A and diagonalizing
the action. To see this directly, average a local defining function
of A with the character of its linear part; it still vanishes on A
and has the same nonzero differential. Average a complementary local
coordinate with its tangent character. The pair has invertible linear
part and is equivariant. Choose a group generator so the normal
character is zeta. Then

    A={X=0}, (X,Y)->(zeta*X,zeta^alpha*Y),
    gcd(alpha,n)=1.

For h=pi composed with rho, the two target pullbacks P,Q are invariant,
have no linear terms, and

    Jac(P,Q)=epsilon*X^(e-1), epsilon(0)!=0.         (4)

This follows in codimension one from the generic ramification index
and quasi-etaleness; a holomorphic Jacobian cannot have an extra
isolated zero. Other omitted primes are disjoint and supply no local
critical branch at this point.

## 2. Read quotient order from ONE actual Jacobian coefficient

Invariance of P,Q gives

    Jac(P,Q)(g(X,Y))*det(g)=Jac(P,Q)(X,Y).

Compare the nonzero coefficient epsilon(0)*X^(e-1). It yields
e+alpha=0 modulo n, hence the marked type and coprimality in (1).

Now expand the ACTUAL invariant P,Q as convergent Taylor series. A
contribution to the coefficient of X^(e-1)*Y^0 in their Jacobian must
come from one pure monomial X^i and one monomial X^j*Y, with

    i>=1, j>=0, i+j=e.                              (5)

At least one such contribution is nonzero because the full coefficient
is nonzero. Invariance of the pure monomial forces n|i. Thus n<=i<=e.
Coprimality rules out n=e, so n<=e-1. This uses neither a coefficient
cutoff nor a search: (5) lists ALL monomial pairs that can contribute
to that single coefficient, including arbitrarily high other terms.

The original full Euler budget has d-1>=e*f>=e, giving n<=d-2.
This bounds all quotient orders at a fixed field degree without
classifying deeper boundary jets.

## 3. All-index analytic normal form under the scalar-Morse condition

Choose a generic target scalar v whose quadratic part is nondegenerate.
It may simultaneously be chosen transverse to the chosen image branch.
Let u be a complementary target coordinate. Analytic factorization at
the two distinct tangents of v=0 gives source coordinates v=X*Y.
The smooth critical branch A is not an axis, because v was chosen
nonconstant on its image. After exchanging axes write A={X=beta(Y)},
ord beta=b>=1. The v-preserving coordinate change used in PF105 makes
A={X=Y^b}. Put m=b+1. Equation (4) becomes

    (X*d_X-Y*d_Y)u=unit*(X-Y^b)^k,  k=e-1.          (6)

Subtract a convergent target function of v to remove all diagonal
monomials X^iY^i. Use weights(b,1). The first right-hand weight is k*b.
If m divides k, the nonzero binomial coefficient at j=k/m multiplies
the diagonal monomial X^(k-j)*Y^(b*j), which the left-hand operator
annihilates. This is impossible. Therefore

    m does NOT divide k.                            (7)

Otherwise the first weighted part of u is exactly

    epsilon0*sum_(j=0)^k [(-1)^j binom(k,j)/(k-m*j)]
                               *X^(k-j)*Y^(b*j).   (8)

Its restriction to A has order k*b and a NONZERO coefficient. Indeed
the coefficient sum is

    C(k,m)=-(1/m)*k! /
                 product_(j=0)^k(j-k/m).            (9)

The finite partial-fraction identity behind (9) follows by matching
residues and the value at infinity. Condition (7) makes every denominator
nonzero. All higher analytic terms have higher weighted order and
cannot cancel this term. Thus the normalized local degree q_h of
A->normalization(h(A)) divides

    gcd(m,k*b)=gcd(m,k).                            (10)

In particular q_h divides k, not merely q_h<=k. Also (8) gives orders
k and k*b on the axes v=0, so the local surface degree is

    deg_0 h=k+k*b=k*m.                              (11)

In the actual singular situation the primary power-map/location input
gives q_h=n*f. Let a be the image-branch multiplicity. Then
m=n*f*a, k=n*f*t for an integer t>=1. By (7), a does not divide t,
so a>=2. In particular b>=2: the possible transverse-to-both-axes
case m=2 would force q_h=1 and cannot occur here.

Multiplicativity of finite local degrees and (11) give ell=k*f*a.
Subtract the generic cost e*f to obtain eta in (2). The FULL global
budget, including all other primes and excesses, gives

    d-1>=e*f+eta=ell=k*f*a>=2*k*f.

Since k>=n*f>=2f, (3) follows. No assumption on the unknown higher
coefficients or on the immersion of the finite branch image enters.

Finally n|k and the marked type 1/n(1,-e) make alpha=-1 modulo n.
The quotient is therefore A_(n-1), with invariant ring generated by
X^n,Y^n,XY and relation X^n*Y^n=(XY)^n. The divisibility and the
quadratic field-degree bound are not asserted outside this sector.

## 4. A conductor identity WITHOUT the scalar-Morse hypothesis

For ANY smooth finite germ h=(P,Q) with sole smooth critical branch
A={x=0} and Jac(h)=unit*x^(e-1), let F(P,Q)=0 be the reduced equation
of its chosen irreducible LOCAL image branch C. This is not the product
of other branches of a globally singular image curve. Factor

    F composed with h=x^e*H, x does not divide H.

The other factors of H are generically unramified preimage curves;
none is an additional local critical branch. Put q=deg_0(A->C^nu),
and let c_C be the conductor of the image plane branch. Then

    I_0(A,H)=q*(c_C-1)+1.                           (12)

Proof: choose P generically on C, with order a on its normalization.
The characteristic-zero Dedekind/polar identity is

    ord_Cnu(F_Q)=c_C+a-1.

Differentiate F composed with h=x^e H and wedge with dP. After dividing
by x^(e-1) and restricting to x=0, the equality, up to sign, is

    F_Q(P(0,y),Q(0,y))*epsilon(0,y)
                       =e*H(0,y)*(d/dy)P(0,y).

The orders are q*(c_C+a-1) on the left and
I_0(A,H)+q*a-1 on the right, proving (12).

The precise primary conductor/polar statements are Garcia Barroso--Ploski,
Proposition2.1(iii), its proof, and Corollary3.4, full PDF3--4,10
(printed121--122,128), freshly rendered and read. Only characteristic
zero and an irreducible chosen branch are used. See the
[manifest](sources/cyclic-quotient-morse/manifest.json).

For a singular actual point with q=n*f, (12) becomes

    I_0(A,H)=n*f*(c_C-1)+1 = 1 modulo n*f.           (13)

In particular c_C=0 would give a negative intersection, so the chosen
image branch is nonimmersive. Its unramified companions MUST meet A.
Companion existence was already stated by Borisov Theorem3.6; (12)--(13)
give the exact intersection and congruence. These companions are ordinary
preimage curves, not automatically omitted divisors. This identity
remains available when the scalar quadratic pencil is degenerate.

## 5. Extremal quotient orders and complete residue corollaries

Suppose n=e-1>=2. In (5), the only possible pure power is i=n,
and the other monomial is X*Y. Thus at least one target component has
an XY term. For n>2 the only invariant quadratic monomial for the
type 1/n(1,-1) is XY, so this already gives the scalar-Morse condition.
For n=2,e=3 the nonzero quadratic Jacobian forces two independent
quadratic initial forms and hence a generic nodal scalar, as in PF105.

Section3 now gives n*f|k=n, so f=1. Write m=n*a, a>=2. Dividing the
two restricted orders in (8) by q_h=n gives (a,n*a-1), which are
coprime. This is one Puiseux pair and ell=n*a. In particular an actual
extremal point requires d>=2e-1. Singular points themselves are not
excluded; their exact local type and degree are restricted.

For e=4, (1) leaves only n=3; for e=6 it leaves only n=5. Both are
extremal, so their residue maps are unramified at their quotient point.
The published power-map/location input then forces global f=1. Combined
with PF105 and the unramified case this proves the set {1,2,3,4,6}.

If d<=10 and f>1, the first remaining e is at least5, so e*f>=10>d-1.
This proves the full through-ten residue theorem, without excluding
every degree-ten or lower Keller map.

More generally, if f>1 and d<=4*f^2, (3) proves that EVERY scalar
quadratic pullback at that point is degenerate. A pencil of binary
quadratics all of whose members are degenerate is zero or is spanned
by one linear square: two independent squares have a nondegenerate
generic sum. Since the restrictions to A have orders at least n*f>=4,
this square must vanish on the tangent of A; it is the normal coordinate
squared. For n>2, such a square is not invariant, so both quadratic
parts must be zero. For n=2, a common normal square is also possible.
This is a coordinate-invariant necessary escape, not its realization.

## 6. Keep genuine countercontrols and the unsolved case

The unbounded finite quotient family in PF105 has

    n=2, k=2s, f=s, a=2, ell=4s^2,
    generic degree4s^2 and complement C* times C.

It saturates BOTH local inequalities k>=2f and ell>=4f^2 for every
s>=1. Thus local quotient/ramification information alone cannot demand
a larger ell. The actual Keller budget requires one more sheet because
its marked open is A2; none of these maps supplies that open.

The scalar-Morse assumption must also be retained. For any integer
m>=3 the genuine finite polynomial map

    P=X^2+Y^m,
    Q=X^2*Y+(m/(m+1))*Y^(m+1)

has Jac=2X^3 and one smooth critical line, but its quadratic pencil
is spanned only by X^2. It is finite of degree2(m+1), since
Y^(m+1)-(m+1)*P*Y+(m+1)*Q=0 and X^2=P-Y^m. Its critical-line
restriction is primitive. This is a smooth non-plane-complement control,
not a singular Keller counterexample; it demonstrates that finiteness
and one smooth critical curve do not force the scalar-Morse hypothesis.

We have NOT proved the general inequality q_h<=e-1 outside the scalar-
Morse sector, or an automatic scalar-Morse theorem at actual residue
points. Candidate non-Gorenstein quotients permitted by (1), such as
1/3(1,1) at e=5, and common-square quadratic pencils remain unresolved.
The exact conductor balance (12) applies there; deeper parameter or
coefficient searches are not substituted for that missing argument.

[The checker](verify_cyclic_quotient_morse.py) and
[receipt](cyclic_quotient_morse_receipt.json) retain exact monomial
coefficient, resonance, partial-fraction, quotient, conductor, local
degree and global-budget controls. Written unbounded proofs and primary
imports, not finite examples, establish the scopes above. The full
Wright goal, all f1 packets and actual Laurent-terminal finite algebra
remain OPEN. Published quartic retirement and closure-before-review
priority persist; formal/external verification remains false.
