# Degree-six conic equations exclude one cubic escape and restrict the other

14 November 2026. PF2R7V1 and PF2R7G1 continue the
[three-escape reduction](../jc2-r7-universal-contact/PROOF.md), under its
explicit degree-six Keller and entire single ramified affine-prime
(e,f)=(2,1) hypotheses. Every omitted unramified prime is retained.
The arguments below do not use Borisov's type1 source-finality assertion,
singularity at the common infinity point, or any conclusion of his
Theorem 3.8. The earlier boundary, positivity and affine smoothness inputs
remain as specified in PF2R7INF1. No formalization, external review or
novelty priority is claimed.

We exclude the third infinity escape, of indices (2,4). In the first
escape, (1,5) with infinity fibre degrees (1,0), we force a multiple zero
of p'. In particular that first escape is impossible for both retained
cusp classes and every immersed specialization with squarefree p'.
The quartic escape remains open for the general family. The subsequent
[PF2R7EX2 proof](../jc2-r7-quartic-contact/PROOF.md) excludes that escape
for both retained cusp classes, completing their scoped exclusion without
source finality. The two-row window below is retained as its dependency.

## 1. A primitive conic equation remembers the whole field degree

Let z be the rational pencil coordinate and suppose the family image
conic has a defining equation over C(z). Clear its denominators and
divide out the common polynomial factor of all coefficients to obtain
F(z;x,y) in C[z,x,y]. The conic is irreducible over C(z), so Gauss's
lemma makes F irreducible in the polynomial ring. If the source field
is C(x,y,z), the primitive polynomial F is the minimal equation of z
over C(x,y), up to a scalar in that field. Consequently

    deg_z F = [C(source):C(x,y)] = 6.                  (D)

This is an equality for the actual irreducible conic, not a bound on an
arbitrary eliminant. Both cubic escapes below have source field
C(z,x) or C(z,w), with w=y-z*x, so C(x,y,z) is indeed the full field.

Along the normalization of the ramified curve write x=p(t),y=q(t), with
the retained monic degrees (6,10), and s=q'/p'. The point t=infinity is
the only pole of p or q; their pole orders there are six and ten. The
restricted pencil coordinate z(t) has degree three in these two escapes.
At infinity let k be its local degree over z0=z(infinity); 1<=k<=3.
If z0=infinity, z has pole order k. In either case

    w=q-z*p has pole order 10,
    s-z has pole order 4 and leading coefficient 5/3.   (P)

The inequalities k+6<=9<10 and k<=3<4 justify both assertions when
z0=infinity. These estimates use exact orders, not generic coefficients.

## 2. The varying-direction parabola is impossible

PF2R7INF1's (2,4) escape has conics

    x=A(z)*w^2+B(z)*w+C(z), y=z*x+w, A!=0.

Choose polynomials a,b,c,d with no common factor and d!=0 such that
(A,B,C)=(a/d,b/d,c/d). Its primitive conic equation is

    F=d*x-a*(y-z*x)^2-b*(y-z*x)-c.

The coefficients of x^2,y,y-independent constant and x respectively
give from (D)

    deg(a)<=4, deg(b)<=6, deg(c)<=6, deg(d)<=7,
    deg(d)=7 implies deg(b)=6.                        (B)

The last implication follows because the x coefficient is d+z*b,
whose degree must be at most six. Zero b or c is allowed; a is nonzero.

Tangency along the simple ramified curve gives

    2*a(z)*w+b(z)=d(z)/(s-z).                         (T)

For completeness, differentiating incidence along E and using
x_z=x*x_w yields p'=x_w*(q'-z*p'), hence x_w=1/(s-z).
Thus (T) holds as a rational identity on the entire normalization.

### 2.1 The limiting base point cannot be finite

Suppose z0 is finite. Let alpha,beta,delta be the vanishing orders of
a,b,d at z0. If b is zero the argument is even immediate. By (B),
0<=alpha<=4, while beta,delta>=0. On E the three terms in (T) have orders

    k*alpha-10, k*beta, k*delta+4.

The first is at most two, whereas the last is at least four. Therefore
the first two must cancel their leading terms, requiring

    k*(alpha-beta)=10.

No k in {1,2,3} and alpha in {0,...,4} satisfies this with beta>=0:
k=1 needs alpha>=10, k=2 needs alpha>=5, and k=3 cannot divide10.
Thus z0 is not finite.

### 2.2 The limiting base point cannot be infinity

Write n=deg(a), m=deg(b), l=deg(d), treating b=0 separately. The pole
orders in (T) are

    P1=k*n+10, P2=k*m, P3=k*l-4.

For k=1, P1>=10 exceeds both P2<=6 and P3<=3, impossible. For k=3,
these numbers have distinct residues 1,0,2 modulo3; the first two
cannot cancel, and their maximum cannot equal P3.

For k=2, a sum without leading cancellation would require P1=P3=10,
n=0,l=7. But (B) would then give m=6, so P2=12, a contradiction.
Leading cancellation is the only remaining alternative. It requires
m=n+5, leaving precisely (n,m)=(0,5) or (1,6). It also gives
a*w~(-b/2), so the leading term of a*w^2+b*w is (b/2)*w. In the
incidence equation

    d*p=a*w^2+b*w+c

that term has pole order 2*m+10>=20, exceeding c's possible pole
order12. Comparison with d*p forces 2*l+6=2*m+10, or l=m+2.
The pair m=5,l=7 violates (B); m=6 would require l=8>7.
This exhausts z0=infinity too.

Thus the varying-direction parabola escape is impossible, proving
PF2R7V1. The argument actually needs only (D), the critical contact
identity, polynomial degrees (6,10), and degree-three z on E. It does
not require the gap equations or squarefreeness of p'.

## 3. Reciprocal conics force a finite multiple root of p'

PF2R7INF1's first escape, now of ramified projection degree three, has

    y=z*x+alpha(z)+beta(z)/(x-r(z)), r'!=0, beta!=0.

Its primitive equation has the form

    F=a(z)*(x*y-z*x^2)+b(z)*x+c(z)*y+d(z), a!=0.

Here r=-c/a. Equation (D) gives

    deg(a)<=5, deg(b),deg(c),deg(d)<=6.                (G)

The tangency relation F_x+s*F_y=0 and its combination with incidence are

    a*(q+p*s-2*z*p)+c*s+b=0,                          (S)
    a*p^2*(s-z)-c*(q-p*s)-d=0.                      (I)

At infinity the leading coefficients of q+p*s, s, p^2*(s-z), and
q-p*s are respectively 8/3,5/3,5/3,-2/3, with pole orders
10,4,16,10. Terms involving z*p are lower order. All leading
coefficient comparisons below refer to these fixed nonzero numbers.

### 3.1 A finite limiting value

Let alpha0,beta0,gamma0,delta0 be the orders of a,b,c,d at finite z0.
The a-term in (I) has order k*alpha0-16<=-1 by (G), so it must cancel
the c-term; d is regular. Thus

    k*(alpha0-gamma0)=6.

If A0,C0 are the corresponding nonzero leading coefficients after
substitution on E, cancellation in (I) requires C0=-(5/2)*A0.
In (S), those same terms then have leading coefficient
(8/3)*A0+(5/3)*C0=-(3/2)*A0, which is nonzero. The b-term must match:

    k*(gamma0-beta0)=4.

Together these require k to divide both6 and4, and
alpha0=beta0+10/k<=5. The only possibility is

    k=2, ord_(z0)(a,b,c)=(5,0,2).                    (F)

In particular a is a constant multiple of (z-z0)^5. Since z has degree
three, the fibre over z0 consists of infinity with multiplicity two and
one finite parameter u with multiplicity one.

At u, p and q are regular. The ratios of the two target partials are

    F_y=a*p+c, F_x=a*q-2*z*a*p+b.

By (F), F_y has order exactly two and F_x order zero. Tangency gives

    ord_u(p')-ord_u(q')=2.                            (M)

Here the derivatives are taken in the original polynomial parameter t;
z-z0 is a uniformizer at u. Thus p' has a root of multiplicity at least
two. If the image is immersed at u, its multiplicity is exactly two.

### 3.2 An infinite limiting value

Set n=deg(a), j=deg(b), m=deg(c), l=deg(d). In (I) the first two
pole orders are k*n+16 and k*m+10. The last is k*l<=6*k.
For k=1 or2 it cannot cancel either leading a-term. For k=3 the a/c
orders have residue1 modulo3 while d's order has residue0, so the
largest a/c term again must cancel with the other one. Hence
k*(m-n)=6 and C0=-(5/2)*A0. Equation (S) then has a nonzero residual
leading coefficient, requiring k*(j-n)=10. Under (G) this leaves exactly

    (k,deg(a),deg(c),deg(b))=(2,0,3,5) or (2,1,4,6).  (J)

There is again one finite parameter u over z=infinity, with local degree
one. Regularity of p,q in F also forces l<=j: otherwise d would dominate
every other term. At u the b term dominates F_x and c dominates F_y,
with degrees differing by two. Thus (M) holds here as well.

Consequently every reciprocal-conic escape forces k=2 at infinity,
one finite u satisfying (M), and exactly one of the coefficient patterns
(F) or (J). This proves PF2R7G1. The coefficient patterns are necessary,
not realized conic covers or sufficient existence conditions.

## 4. Consequences for the retained normal-r7 problem

The balanced infinity pair (3,3) was already excluded without finality.
PF2R7V1 now excludes (2,4) too. The only remaining index pair is (1,5):

| Infinity fibre degrees | Ramified projection degree | Further necessary condition |
|---|---:|---|
| (1,0) | 3 | p' has a finite multiple root; (F) or (J) holds |
| (2,0) | 4 | two moving infinity points; quartic contact gate remains |

The two exact nonimmersed PF2R7CUSP1 classes have squarefree p'. Their
retained literal Bezout certificates verify this without root sampling.
They therefore cannot realize the first row either: **both cusp classes
are reduced to the one quartic infinity escape without source finality.**
The same holds for every immersed pair with squarefree p'. An immersed
pair with a double root of p' is not automatically excluded; the retained
witness p=t^6+t^3 has such a root, so a squarefreeness assumption cannot
be silently imposed on the whole family.

This does not by itself remove the finality dependency from the older
whole-family exclusion. Independently closing the quartic row, and the
special multiple-root reciprocal row, remains necessary. The complete
single ramified-prime hypothesis must still be justified in any intended
application to a terminal passport. Normal r5, nonnormal slices,
(75,125), and JC2 remain open.

## 5. Reproduction and assurance

[verify.py](verify.py) checks the primitive coefficient templates,
tangency eliminations, all bounded integer order/degree patterns, the
incompatible leading coefficients, retained squarefreeness certificates,
and genuine rational cubic-contact controls outside the excluded degree
hypotheses. Mutated coefficient bounds are used as negative controls.
[receipt.json](receipt.json) binds the sources and outputs.

The finite enumerations exhaust precisely the explicitly proved degree
and local-degree bounds; they are not a graph census or a general search
for Keller maps. Irreducibility/field identification, orders on P1 and the
geometric inheritance from PF2R7INF1 are written proof steps. The checker
does not formalize them or independently reprove the cited geometry.
