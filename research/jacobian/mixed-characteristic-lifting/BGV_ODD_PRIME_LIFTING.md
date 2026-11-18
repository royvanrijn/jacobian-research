# Compatible Witt lifts of odd-prime BGV plane maps cannot algebraize polynomially

For every odd prime p, a published family of noninvertible plane maps has
compatible polynomial determinant-one lifts over Z/p^n for **every n**.
A radial normalization constructs this tower uniformly. Its restricted
p-adic limit is algebraic but is not polynomial; consequently the degrees
of its polynomial Witt reductions are unbounded.

The family is due to Borisov, Gabber and Vasiu. The normalization is an
explicit Hensel and Euler construction; no novelty or major breakthrough
is claimed. It supplies a precise failure of the inference from all finite
Witt lifts to a polynomial characteristic-zero lift. It does not resolve
JC2 or prohibit other lifting architectures.

## 1. Published input and the actual unresolved target

For a prime p, in F_p[x,y] put

    z = x + x^p y,
    P = z + (xz)^p,
    Q = y + x^(p(p-1)) z^(p(p-1)+1).                     (1)

Theorem 28.6(1), (2.e) of A. Borisov, O. Gabber and A. Vasiu,
[*On endomorphisms of affine spaces and the Jacobian problem*,
arXiv:2609.05746v1](https://arxiv.org/html/2609.05746v1), supplies:

* Jac(P,Q)=1;
* geometric degree p^2-p+1;
* surjectivity over an algebraic closure when p is odd.

The specialization is l=m=r=j=1, s=q=p-1, delta=-1,
h(t,w)=tw, g(t)=(-t)^(p-1). All hypotheses of (2.e) hold:
q=ms, s=pm-1, and m=1 <= (p-1)l. The statement is on PDF pages
131-132; the proof, including the degree calculation and surjectivity
criterion, is on PDF pages 132-133. The
[retained primary PDF](sources/bgv-2609.05746v1.pdf) is this version.
The older author-hosted BGV160.pdf has a different statement and does not
supply this input.

The Jacobian assertion can also be read directly. Write a=(xz)^(p-1).
In characteristic p,

    dP=dz,       dQ=dy+a^p dz,       dz=dx+x^p dy,

so dP wedge dQ=dx wedge dy. Degree and surjectivity are imported from the
published preprint, not certified by the checker here.

JC2 asks whether, for every characteristic-zero field k and all
P,Q in k[x,y] with Jac(P,Q) in k*, the pair has a polynomial inverse.
A lift of (1) to a polynomial determinant-one pair over Z_p with suitable
noninvertibility evidence would disprove JC2 over Q_p. Passing finitely
many Witt digits neither produces such a pair nor proves that its degree
can be bounded. A generic-degree statement after reduction alone must
not be promoted to generic-degree preservation in characteristic zero.

## 2. Complete first-obstruction criterion

Let F=(f,g) in F_p[x,y]^2 have Jac(f,g)=1. Choose integer coefficient
lifts f0,g0, and define

    H = (Jac(f0,g0)-1)/p mod p.

For h=sum h_ij x^i y^j define the polynomial coefficient extractor

    C_p(h) = sum_(a,b>=0) h_(pa+p-1,pb+p-1) x^a y^b.     (2)

Then F has a polynomial determinant-one lift modulo p^2, with arbitrary
new support allowed, **if and only if C_p(H)=0**.

Here is a complete elementary proof. A polynomial two-form h dx wedge dy
is exact precisely when every coefficient selected in (2) vanishes.
If i+1 is nonzero modulo p, the monomial is the differential of
x^(i+1)y^j dy/(i+1). Otherwise, if j+1 is nonzero, use
-x^i y^(j+1) dx/(j+1). Both partial derivatives have zero coefficient
at every excluded exponent pair, proving necessity as well.

Changing the lifts by pA,pB changes H by

    Jac(A,g)+Jac(f,B),

whose two-form is d(A dg-B df). Thus the class in top polynomial de Rham
cohomology is independent of the initial lifts. Conversely suppose
H dx wedge dy is exact. Choose eta=u dx+v dy with d eta=-H dx wedge dy.
Because df,dg form a basis of polynomial one-forms, there are polynomial
A,B with eta=A dg-B df. Explicitly,

    A=-f_y u+f_x v,       B=-g_y u+g_x v.                (3)

Then f0+pA,g0+pB have determinant one modulo p^2. There is no support
bound in this argument.

The same argument applies to the extension of any *fixed* determinant-one
pair modulo p^n to modulus p^(n+1), with H the next normalized determinant
digit and (f,g) the initial reduction. It does not identify the obstructions
arising from different choices at earlier digits.

This is ordinary polynomial de Rham obstruction theory. Stronger integral
Keller equations in positive characteristic have substantial prior
literature, including S. Maubach and A. Rauf,
[*A new formulation of the Jacobian Conjecture in characteristic p*,
arXiv:1507.02946](https://arxiv.org/abs/1507.02946). No result from that
paper is needed for the elementary proof above, and no novelty is inferred
from this comparison.

## 3. Pth-power additions preserve the first class

Adding A^p to one coordinate does not change the first obstruction class,
provided the residue pair has determinant one. For literal integer lifts,
the normalized change of d f wedge d g modulo p is

    A^(p-1) dA wedge dg,

or its analogue in the other coordinate. The one-form A^(p-1)dA is closed,
so its wedge with dg is exact. Changes in both coordinates have a cross
term divisible by p^2 and give the same conclusion. This statement is
about the class, not equality of the normalized error polynomials.

In particular arbitrary coordinatewise square additions to Mondello's
characteristic-two map preserve its nonzero first class. That direction
cannot bypass its modulo-four obstruction.

## 4. Uniform odd-prime vanishing

For (1), omit the pth-power summand (xz)^p from P when calculating the
class. This leaves the pair (z,y+a^p z), where a=(xz)^(p-1). Over the
integers,

    dz wedge d(y+a^p z)
      = dz wedge dy + p a^(p-1) z dz wedge da.

Consequently its normalized error modulo p is

    H' = x^(p-1)y + x^(p^2-1) z^(p^2-p+1).             (4)

For the second term use

    Jac(z,a) = -(p-1)x^p z (xz)^(p-2)
             = x^p z (xz)^(p-2) mod p.

This proves (4) including the coefficient and sign.

The first monomial in (4) is selected by C_p only when 1=p-1,
that is p=2. Expand the second using z=x(1+x^(p-1)y).
A summand with y exponent j has x exponent

    (p^2-1)+(p^2-p+1)+(p-1)j = 2p^2-p+(p-1)j,

which is -j modulo p. If j=-1 modulo p, that x exponent is 1 modulo p,
and it can be -1 modulo p only when p=2. Therefore C_p(H')=0 for
every odd prime p. Sections 2-3 give the unrestricted conclusion:

**Proposition.** For every odd prime p, the pair (1) has a polynomial
lift over Z/p^2 with Jacobian determinant one. No degree bound is imposed.

At p=2 the extractor is instead

    C_2(H)=1+x^3+x^4y,

the previously known nonzero Mondello obstruction. The current odd-prime
calculation changes that particular obstruction; it does not solve
characteristic-zero lifting.

## 5. Exact p=3 witnesses and their scope

The checker makes a deterministic choice in Section 2: integrate x first
when possible, integrate y otherwise, and represent all F_3 coefficients
by 0,1,2. It retains the complete correction digits in the
[receipt](bgv_odd_prime_lifting_verification.json). Those data reconstruct
actual integer polynomials, not just solution flags.

The checker verifies determinant one modulo 9 and modulo 27 using exact
integer sparse arithmetic. It computes the next obstruction for the
resulting fixed modulus-27 pair. Its nonzero value excludes extension
of that exact pair to modulus 81 at any polynomial degree. Other
modulus-27 pairs can have different third obstruction classes.
There is no assertion that unrestricted lifting stops at modulus 81.

Finite controls at p=2,3,5,7 check the normalized errors, simplification
(4), and coefficient extraction, with p=2 as a required negative control.
The proposition for every odd p is the written proof in Section 4.
Sections 6-7 give a different, compatible radial tower at every Witt
length and prove that its degrees are unbounded. The canonical primitive
choice above is retained as a negative control: a poor correction choice
can create a later obstruction despite existence of another infinite tower.

## 6. An Euler primitive gives a compatible tower at every length

Let p be odd, R=Z_p<x,y> the p-adically complete restricted power-series
ring, and F0=(f0,g0) an integral polynomial pair with J=Jac(F0) in
1+p Z_p[x,y]. Put D=x partial_x+y partial_y. Assume that, for each
nonconstant monomial x^i y^j appearing in J, the integer i+j+2 is a
p-adic unit. The constant denominator 2 is also a unit.

Define the integral polynomial

    A = sum_(i,j) 2 J_ij/(i+j+2) x^i y^j.               (5)

Then A=1 mod p, and (2+D)A=2J. Hensel's lemma in the complete ring R
applied to

    E(T)=T^2 A(Tx,Ty)-1

gives a unique lambda in 1+pR with E(lambda)=0: E(1) is divisible by p,
and E'(1) is 2 mod p. This application can also be proved directly by
successive digit corrections, because at every stage the next coefficient
is uniquely determined by division by 2. Each digit is a polynomial;
therefore lambda belongs to the restricted completion R.

Set phi=(lambda x,lambda y). This is a restricted formal change of
variables congruent to the identity mod p. It is algebraic over
Q_p(x,y), because its multiplier satisfies the displayed polynomial
E(T). For v=D(lambda)/lambda, applying D to E(lambda)=0 gives

    (1+v)(2A+DA)(phi)=2A(phi).

The determinant of phi is lambda^2(1+v). Consequently

    Jac(F0 composed phi)
       =J(phi) lambda^2(1+v)
       =lambda^2 A(phi)=1.                              (6)

Differentiation is continuous on R, so these are identities in the
restricted ring. Reduction modulo p^n gives a polynomial pair for every
n; the pairs are compatible and retain the original reduction modulo p.
The Hensel equation proves existence at all lengths, independently of
any finite computation.

The BGV pair (1), lifted by the same integer formulas, satisfies every
hypothesis uniformly for odd p. Its P monomials have total degrees 1 or
0 modulo p; its Q monomials have degrees 1 modulo p. Differentiation and
multiplication show that J has total degrees 0 or -1 modulo p. Thus
all denominators i+j+2 are 2 or 1 modulo p. Also J=1 mod p by Section 1.

**Theorem.** For every odd prime p, the BGV plane map (1) has the
compatible determinant-one polynomial Witt lifts produced by (5)-(6)
over Z/p^n for every n>=1. Their restricted limit has algebraic coordinate
functions over Q_p(x,y).

This is a general volume-normalization principle under the explicit Euler
denominator hypothesis. It does not assert the hypothesis for arbitrary
Keller reductions. In particular p=2 is excluded both in Hensel division
and in Section 4's negative control.

## 7. An algebraic-degree obstruction forbids polynomial algebraization

This argument also rules out cancellation of the algebraic radial multiplier
inside the output coordinates.

For odd p, the highest monomials of the original integer pair F0 are

    P_top=x^(p^2+p)y^p,      Q_top=x^(p^3)y^(p^2-p+1).

Their Jacobian has leading coefficient

    (p^2+p)(p^2-p+1)-p*p^3=p.

Thus J and the Euler primitive A have the same total degree

    D=p^3+2p^2+p-1,

which is odd. Put M=D+2. In independent variables s=y/x and w=x lambda,
the radial equation becomes

    w^2 A(w,sw)=x^2.                                     (8)

The left side has odd w-degree M over Q_p(s). Therefore it is not a square
in Q_p(s,w), and x^2-w^2A(w,sw) is irreducible over that field as a
polynomial in x. Gauss's lemma first in Q_p(s)[w,x], and then with x in
the coefficient field, proves that (8) is irreducible of degree M in w
over K=Q_p(s,x)=Q_p(x,y). Primitivity in the second application holds
because the w^2 coefficient is 1 and the constant coefficient is -x^2.
Consequently

    [K(lambda):K]=M>p^2.                                 (9)

In contrast, the original integer polynomial map F0 has characteristic-zero
generic degree exactly p^2. To see this, use source coordinates

    u=xz,      z=P-u^p,
    x=u/z,     y=(z-x)/x^p.

These are field identities. The second target coordinate gives

    Q u^p=(P-u^p)^(p+1)-u(P-u^p)^(p-1)+u^(p^2)(P-u^p).  (10)

The numerator on the right has u-degree exactly p^2, with leading
coefficient -pP and constant coefficient P^(p+1). It is coprime to u^p.
The rational function in (10) therefore has degree p^2 over Q_p(P).
The standard rational-function field-degree identity follows, for example,
by clearing the parameter in N(u)-Q u^p: coprimality makes that polynomial
irreducible as a polynomial linear in the independent parameter Q, and
Gauss's lemma gives irreducibility over Q_p(P,Q). The reconstruction
identities recover both source coordinates. This proves the asserted
generic degree, without importing a characteristic-zero analogue of the
characteristic-p degree theorem.

Suppose both normalized outputs Fhat=F0 composed phi belonged to K.
Put E=Q_p(lambda x,lambda y) and K0=Q_p(Fhat). The coordinates lambda x,
lambda y are algebraically independent: the determinant of phi is a
unit modulo p and hence nonzero in characteristic zero. By (10),
[E:K0]=p^2. But

    E K=K(lambda),

so the elementary compositum inequality gives

    [K(lambda):K] <= [E:K0]=p^2,

contradicting (9). Thus the two coordinates of Fhat cannot both be
rational functions, and in particular cannot both be polynomials.

If the polynomial reductions of Fhat modulo all p^n had degrees bounded
by one D0, every coefficient of degree greater than D0 would be zero
modulo every p^n and hence zero in Z_p. That would make Fhat polynomial.
Therefore the degrees of the entire compatible radial tower are unbounded.
This is about this tower, not a lower bound for every possible lifting
architecture of the original residue map.

An independent short obstruction uses the source line. Here J(0,y)=1,
so A(0,y)=1 and the unique odd-p Hensel branch has lambda(0,y)=1.
Hence Fhat(0,y)=(0,y), Fhat(0)=0 and JFhat(0)=I. If Fhat were polynomial,
Gwozdziewicz's [*Injectivity on one line*, Theorem 1.1 and Lemma 2.3]
(https://arxiv.org/html/alg-geom/9305008v1), applied over an algebraic
closure of Q_p, would make it an automorphism. Its formal inverse has
Z_p coefficients by the degree recursion for a map with identity linear
part. Its polynomial inverse must equal that formal inverse. Reduction
would therefore make (1) an automorphism, contradicting generic degree
p^2-p+1>1. This primary theorem is already imported in the
[plane-completion result](../universal-keller-plane-completion/PROOF.md).
The algebraic-degree proof above does not depend on this additional import.

For p=3, the exact checker verifies compatible radial reductions modulo
9,27,81, their full polynomial Jacobians and the fixed line. These finite
checks control the formulas; uniform all-length existence and unboundedness
are proved above. Prime controls also verify the leading coefficient in
(10) and the leading Jacobian monomial.

No polynomial characteristic-zero counterexample, general no-lifting
result, independent implementation, formal verification, external human
review, or JC2 solution is claimed. The remaining structural task would
be an architecture permitting polynomial algebraization; increasing the
computed Witt length cannot satisfy that task.

Replay from the repository root:

    python3 research/jacobian/mixed-characteristic-lifting/verify_bgv_odd_prime_lifting.py

Use --write-receipt PATH only to retain an intentional new receipt. An
ordinary replay prints its evidence without modifying the retained one.
