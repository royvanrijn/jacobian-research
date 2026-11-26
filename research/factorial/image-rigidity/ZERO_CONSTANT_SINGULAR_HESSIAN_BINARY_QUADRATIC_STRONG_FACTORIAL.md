# Zero-constant singular-Hessian binary quadratics: a residue proof

Let K be a characteristic-zero field and Lambda_2(X^rY^s)=r!s!.
For a nonzero polynomial

 f=A X^2+BXY+C Y^2+dX+eY,     B^2=4AC,

let s be its actual number of nonzero monomials.

**Theorem 1.** No four consecutive positive factorial moments of f are
all zero. Moreover, for every n>=1, at least one of
Lambda_2(f^n),...,Lambda_2(f^(n+s-1)) is nonzero.

Thus every zero-constant binary quadratic with singular Hessian satisfies
its actual-support Strong Factorial window. Affine and axis supports are
included. The four-window is stronger than the required window when s=5;
smaller supports are proved separately below.

The central observation is a native residue generating function.
Its exceptional coefficient patterns are incompatible with B^2=4AC.
No analytic integral, coordinate rescaling, generic-rank assumption,
irreducibility theorem or finite enumeration is used.

## 1. A native residue identity

Work first in a finite valued extension L of Q_2. Write v(2)=1, let O
be its valuation ring, and use bars for reduction into its residue field.
Put

 f=(alpha/2)X^2+beta XY+(gamma/2)Y^2+dX+eY,
 P_n=Lambda_2(f^n)/n!,   E(T)=sum_(n>=0) P_n T^n.

Assume alpha,beta,gamma,d,e belong to O. No Hessian condition is needed
in this section. All series are formal power series.

**Lemma 2 (native residue series).** Every P_n belongs to O, and

 Ebar(T)=
 [1+(dbar+ebar+betabar)T+dbar ebar T^2]
 /[(1+alphabar T+dbar^2 T^2)(1+gammabar T+ebar^2 T^2)].       (1)

The denominators have constant coefficient1, so the quotient is defined
formally, including when some parameters reduce to zero.

To prove integrality, let i,j,h,l,m be the occurrence counts of the five
monomials, so i+j+h+l+m=n. The coefficient of
alpha^i beta^j gamma^h d^l e^m in P_n is exactly

 (2i+j+l)!(2h+j+m)!/[2^(i+h)i!j!h!l!m!]
       =j! C_(i,j,l) C_(h,j,m),                         (2)

where

 C_(i,j,l)=(2i+j+l)!/[2^i i!j!l!]
   =(2i-1)!! binom(2i+j+l,2i) binom(j+l,j)

is an integer; (-1)!!=1. Hence all coefficients are integral, and every
term with j>=2 vanishes modulo the maximal ideal, because j! is even.

For q(X)=(alpha/2)X^2+dX, define

 Q(T)=sum_(r>=0) Lambda_1(q^r)T^r/r!,
 I(T)=sum_(r>=0) Lambda_1(Xq^r)T^r/r!.

Their coefficients are integral by the same factorial formula.
Writing q_r=Lambda_1(q^r)/r!, coefficient comparison gives

 q_0=1, q_1=alpha+d,
 q_r=alpha(2r-1)q_(r-1)+d^2 q_(r-2), r>=2.

Consequently, with D=1+alphabar T+dbar^2 T^2,

 Qbar=(1+dbar T)/D.                                    (3)

The coefficient of alpha^i d^(r-i) in Lambda_1(Xq^r)/r! is r+i+1
times its coefficient in q_r. Thus, as a universal identity of formal
series over F_2[alpha,d], before evaluating the residue parameters,

 Ibar=Qbar+T partial_T Qbar+alphabar partial_alphabar Qbar=1/D. (4)

For an explicit check, multiplication of the middle expression by D^2
gives precisely D. This uses polynomial differentiation of the universal
parameter series, not differentiation of a fixed residue-field element.

By (2), only cross powers zero and one survive. Literal expansion of
exp(Tf), coefficient by coefficient and followed by Lambda_2, therefore
gives

 Ebar=Qbar_X Qbar_Y+betabar T Ibar_X Ibar_Y.

Equations(3)--(4) prove (1). Each coefficient involved is a finite sum,
so no convergence or interchange of analytic integrals is required.

## 2. A finite recurrence criterion

Suppose a formal series p(T)=sum_(j>=0)p_j T^j over a field has p_0=1
and satisfies D(T)p(T)=N(T), where D(0)=1, deg(D)=r>=1,
the degree-r coefficient a_r is nonzero, and deg(N)<r.
For every j>=r,

 p_j+a_1 p_(j-1)+...+a_r p_(j-r)=0.                      (5)

No r consecutive positive coefficients can vanish. A block
p_n,...,p_(n+r-1)=0, n>=1, makes (5) at j=n+r-1 force p_(n-1)=0.
Repeating reaches p_0=0 by applying the relation at j=r, contrary to
p_0=1. This argument uses actual coefficients and their initial value;
it does not replace them by an arbitrary locally prescribed state.

## 3. The singular-Hessian relation removes the exceptions

Now impose beta^2=alpha gamma and assume at least one of
alpha,beta,gamma,d,e is a unit of O. Let D,N be the denominator and
numerator of (1), without canceling any common factor.

If d and e are both units, D has degree4 with nonzero leading coefficient
dbar^2 ebar^2, and deg(N)<=2. Section2 excludes four consecutive zero
residue coefficients.

If exactly one linear coefficient is a unit, say d, then ebar=0.
The factor 1+alphabar T+dbar^2 T^2 has degree2, and the other factor
is 1+gammabar T. Thus D has degree3 or2, with nonzero leading
coefficient, whereas deg(N)<=1. Section2 excludes three or two
consecutive zero residue coefficients. The case with e a unit is
symmetric.

Suppose neither linear coefficient is a unit. Then

 N=1+betabar T,   D=(1+alphabar T)(1+gammabar T).          (6)

If both diagonal parameters alpha,gamma are units, D has degree2
and deg(N)<=1, so two consecutive zeros are impossible.
If exactly one is a unit, the relation betabar^2=alphabar gammabar
forces betabar=0. Then N=1 and D has degree1, so every coefficient
is nonzero. If neither is a unit, the same relation forces beta
to be a nonunit too, contradicting the assumption that at least one
of the five parameters is a unit.

These cases are exhaustive. They establish the four-zero exclusion
over every finite extension of Q_2 under the stated normalization.

For arbitrary nonzero coefficients in L satisfying B^2=4AC, consider
the five natural parameters

 (alpha,beta,gamma,d,e)=(2A,B,2C,d,e).

Choose a nonzero member c of minimum valuation and divide the ENTIRE
phase by c. All five normalized parameters are integral, at least one
is a unit, and beta^2=alpha gamma is preserved. The j-th factorial
moment is multiplied by c^(-j), so its vanishing is unchanged.
This proves the four-window for every such phase over L. No variable
is independently rescaled.

## 4. Transfer to every characteristic-zero field

The coefficients generate a finitely generated field F over Q.
Such an F embeds into an algebraic closure of Q_2. Choose a finite
transcendence basis and map it to algebraically independent elements
of Q_2; then extend the embedding through the finite algebraic
extension into the algebraically closed Q_2bar. The independent finite
set exists by induction: a countable subfield has only countably many
algebraic elements in Q_2, whereas Q_2 is uncountable.
The finitely many coefficient images lie in a finite extension of Q_2.
Every polynomial identity and every nonzero coefficient is preserved.
Thus any hypothetical four-zero block over the original field would
give one in the local setting excluded above. This proves the first
assertion of Theorem1 without asserting that the original coefficients
are algebraic over Q.

## 5. Exact smaller supports

It remains to prove the actual-support window when s<=3.

If B!=0, the relation B^2=4AC implies AC!=0. Hence all three quadratic
monomials occur, and s<=3 forces d=e=0. Put a=C,t=B/(2C), so the phase
is exactly a(tX+Y)^2 with a,t nonzero, over the original field.
Direct expansion gives

 Lambda_2(f^j)=a^j(2j)! h_(2j)(t,1),
 h_m(t,1)=1+t+...+t^m.

Since h_(2j+2)-t^2 h_(2j)=1+t, two adjacent zero moments would force
t=-1, where h_(2j)(-1,1)=1. Thus a stronger two-window holds, and
in particular the required three-window holds.

If B=0, then AC=0. If the quadratic part is nonzero, a coordinate
swap puts the phase in the form qY^2+dX+eY, q!=0.
Define H_j=Lambda_1((qY^2+eY)^j), H_0=1. Native coefficients give

 H_j=2q j(2j-1)H_(j-1)+j(j-1)e^2 H_(j-2), j>=2.          (7)

For e!=0, adjacent positive H zeros propagate backward to H_0=0
and are impossible. For e=0, every H_j=(2j)!q^j is nonzero.
The monomial Stein identity Lambda_2(partial_X P)
=Lambda_2(P)-Lambda_1(P(0,Y)) gives

 M_j=d j M_(j-1)+H_j,  M_j=Lambda_2(f^j), j>=1.          (8)

When d,e are both nonzero, three consecutive M zeros would force
two adjacent H zeros, contradicting (7). This is the actual
three-window. When d!=0,e=0, two M zeros force one zero H value,
again impossible; this is the actual two-window. When d=0,e!=0,
the phase is univariate and (7) directly gives its actual two-window.
When d=e=0 it is a single monomial with every moment nonzero.

Finally, if the quadratic part is zero, the nonzero phase is dX+eY.
Its normalized moments are the complete sums h_j(d,e), with h_0=1
and h_j=(d+e)h_(j-1)-de h_(j-2). If de!=0, the nonzero backward
coefficient excludes adjacent zeros and gives the actual two-window.
If exactly one coefficient is nonzero, every moment is nonzero.

These cases exhaust all s<=3. For s>=4, the four-window already
proved gives the required s-window by taking its initial four indices.
This proves Theorem1 in full.

## 6. Scope and classical context

The source gives a self-contained proof for zero-constant singular-Hessian
binary quadratics and also records the stronger uniform four-window.
General nonsingular-Hessian binary quadratics, unrestricted nonzero
constants, full SFC(2), and SFC(3) remain outside this theorem.

Homogeneous powers of linear forms and independent-variable linear
adjunction have classical Strong Factorial proofs: Rocks, On the
incompatibility of diophantine equations arising from the Strong Factorial
Conjecture, author2015 preprint, Theorems3.3 and3.12/Corollary3.13,
https://people.math.sc.edu/rocks/BradyRocksSFC2015.pdf.
The exact statements and adjunction recurrence were inspected.
Those results explain classical subcases; every necessary calculation
is supplied above, so no additional sector theorem is a hidden input.

All mathematical inputs are proved in the preceding sections. Separate
dense, sparse and exceptional-orientation proofs are not premises of this
self-contained result. No formal-verification, external-review or priority
claim is made.
