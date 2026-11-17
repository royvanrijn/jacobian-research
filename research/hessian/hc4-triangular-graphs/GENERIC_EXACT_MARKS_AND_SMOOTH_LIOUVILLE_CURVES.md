# Generic exact marks can hide raw torsion; smooth exact curves have every genus

16 November2026. HC4MYMF0. Over C, ALL parameters below. Written
parameter proofs, not bounded searches. HC4/JC2 and the original full
trace-pairing/full marked-energy implication remain OPEN. External
review DEFERRED.

## 1. What is tested, and what is already known

An ACTUAL plane Keller pair makes lambda=y dx exact on each smooth
inverse image of a target line: d(lambda+P dQ/J)=0, so polynomial
Poincare gives lambda+P dQ/J=dH globally. On a line P=aQ+b its second
term is also exact. This uses the ACTUAL second component.

LP0/1 already distinguishes the RAW mark [1] from unmarked connection
acyclicity, and LG1 already describes point-supported raw torsion.
The new TWO-PARAMETER family gives stronger concrete controls:

* smooth CLOSED embedded plane curves of EVERY genus have lambda exact
  with a POLYNOMIAL primitive on the curve;
* whole-plane polynomial SUBMERSIONS have trivial generic connection,
  generic lambda=0 AND generic volume mark=0, yet NO polynomial Keller
  mate; the raw volume class has EXACT annihilator (P^2).

Neither is an HC4 counterexample. The second disproves replacement of a
global raw mark by its generic localization. LR0's ACTUAL all-irreducible
target frame removes this torsion mechanism for a supplied Keller pair;
its remaining higher-genus problem is not solved here.

## 2. One identity gives smooth exact curves of arbitrary genus

For EVERY k>=1,n>=2 put

    T=x^(k+1)*y, A(T)=1+T^n,
    P_(k,n)=x^k*A(T), G(T)=T+T^(n+1)/(n+1).

There is the GLOBAL POLYNOMIAL one-form identity

    P*y dx + n/(k*(n+1))*d(T^(n+1)) = (xy/k)*dP.       (1)

Consequently on EVERY fiber C_c={P=c}, c!=0,

    y dx=d(-n*T^(n+1)/(k*(n+1)*c)).                    (2)

The primitive is the restriction of a polynomial in x,y with fixed
scalar coefficient 1/c. It is not merely analytic or defined only on
the normalization. Also x dy=d(xy)-y dx is exact.

These are smooth CLOSED curves in the WHOLE plane. On C_c, x!=0 and
(x,y)<->(x,T) is an algebraic change there. Its equation x^k A(T)=c
cannot have both derivatives zero: a zero of A' has T=0, where the
x derivative is nonzero. The valuation -1 at any simple root of A
makes X^k-c/A(T) irreducible over C(T), by Eisenstein using 1/X
at that pole; Gauss primitivity gives irreducibility of the original
plane equation. The k=1 case is immediate.

The degree-k cover of P1_T on the smooth compact completion has n fully
ramified points over roots of A. At infinity it has d=gcd(k,n) points
of index k/d, since c/A has order n. There are no other branch points.
Riemann-Hurwitz and the explicit affine chart give

    2g=(k-1)*(n-1)-d+1,
    punctures=n+d, b1=k*(n-1)+1, chi=-k*(n-1).         (3)

Exactly the n pole points of x and the d points over T=infinity are
omitted. Every other finite T is present. In particular n=2,k=2g+1
realizes EVERY g>=0 with THREE punctures. Smoothness, embeddedness and
exact Liouville periods DO NOT imply rationality. No priority is claimed
for these elementary cyclic curves or the exact-differential method.

## 3. The k=1 specialization is a whole-plane submersion

Set k=1, retaining EVERY n>=2:

    P=x+x^(2n+1)*y^n=x*(1+T^n), T=x^2*y,
    P_x=1+(2n+1)*T^n, P_y=n*x^(2n+1)*y^(n-1).

The explicit unit-ideal identity

    (1-(2n+1)*T^n)*P_x
      +((2n+1)^2/n)*x^(2n-1)*y^(n+1)*P_y=1           (4)

proves NO critical point on A2. The zero fiber has disjoint reduced
components x=0 (an A1) and x^2*y=zeta for the n roots zeta^n=-1
(each a Gm). The ENTIRE family over c!=0 is the product

    D(P) <-> Gm_c times (A1_T minus {A(T)=0}),
    x=c/A(T), y=T*A(T)^2/c^2.                         (5)

Its inverse is (c,T)=(P,x^2*y). Thus the generic curve is rational with
n+1 punctures and b1=n. The genuinely FINITE generic Gauss-Manin
connection is TRIVIAL of rank n: the fixed punctured line gives a
trivialization, with basis dT/(T-zeta). This is stronger than finite
monodromy. Equation(2) also makes its generic Liouville class zero.

## 4. Generic volume vanishes; its raw class has order two

Let D=Jac(P,-), M=R/D(R), R=C[x,y]. This is the RAW quotient, on which
P multiplication descends. The GLOBAL polynomial identity

    D(G(T))=P^2                                      (6)

gives the rational mate Q=G(T)/P^2 on D(P). Therefore [1] becomes zero
in M tensor_C[P] C(P). Nevertheless

    Ann_C[P]([1])=(P^2), P[1]!=0.                     (7)

Here is an all-degree proof. Give x,y weights1,-2; P has weight1 and
D increases weights by2. If D(H)=P, the weight -1 part would do so.
Every weight -1 polynomial is xy*h(T), h in C[T], and

    D(xy*h(T))=x*(T*A*h'+(A+T*A')*h).                (8)

For nonzero h of degree s, the bracket has degree s+n with nonzero
leading multiplier s+n+1. Equality to A forces s=0. Its constant and
T^n coefficients would then require h=1 and (n+1)h=1, impossible.
So P is not in D(R), with NO bound on H's degree. Equation(6) proves
(7): any annihilator ideal containing P^2 is generated by its divisor.
In particular [1]!=0 and NO polynomial Keller mate exists in ANY degree.

The direct mate test extracts weight -2, whose general form is y*h(T):

    D(y*h(T))=(1+(2n+1)*T^n)*h+T*(1+T^n)*h'.         (9)

Its highest multiplier s+2n+1 is nonzero, so it cannot be a nonzero
constant. These are written parameter proofs, not finite censuses.

There is an explicit raw Liouville representative as well. The relative
one-form identification from a polynomial lift of partial_P sends lambda
to tau=[lambda(D)]=[-n*x*T^n]. Identities (1) and

    D(xy)=P+n*x*T^n

give tau=P[1]!=0 and P*tau=0. Hence the raw Liouville class has exact
annihilator (P), and the raw volume class has exact annihilator (P^2).
BOTH vanish generically and BOTH survive globally. With LP's corrected
connection nabla(tau)=-[1], consistent with d(lambda)=-dx^dy.

This agrees with RETAINED LG1: the zero fiber has n+1 components, so
its point torsion has n socle directions. We do NOT assert that the
single [1] generates the entire torsion. Existence of a connection
does NOT alone imply torsion-freeness: point Weyl modules disprove it.

## 5. A genuine open coframe for every k,n

For ALL k>=1,n>=2, on X=D(P) in the original plane,

    Q=G(T)/(k*P^2), Jac(P,Q)=1,
    dx^dy=x^(-k-1) dx^dT,
    dP^dQ=(A(T)/(k*P^2))*dP^dT=dx^dy.                (10)

This is a SURJECTIVE etale map X->Gm_P times A1_Q, of generic degree
k*(n+1). Indeed solve

    G(T)=k*c^2*q, x^k=c/A(T), y=T/x^(k+1).

G has degree n+1, critical points the simple roots zeta of A, and
DISTINCT critical values n*zeta/(n+1). At a critical value exactly
one root is double; the other n-1>=1 roots are not roots of A and remain.
Every target fiber is finite/nonempty: generic size k*(n+1), and size
k*(n-1) at these critical values.

Its smooth full finite normalization over this PUNCTURED target is

    S=C[c,c^-1,T,z]/(c*z^k-A(T)),
    Q=G(T)/(k*c^2), X=D(z), z=1/x.                   (11)

It is finite free of rank k*(n+1): first adjoin the monic degree n+1
equation for T, then the monic degree k equation for z. The dense
D(z) chart is the integral X; the monic relation has no extra component
on z=0. Smoothness holds using the c derivative when z!=0, and A'(T)
when z=0. There are n omitted Gm_c curves z=0,T=zeta. There T-zeta
has order k and Q-Q(zeta) order 2k, so each ramification index is 2k.
Other sheets remain. For n=2 the two omitted points map to DISTINCT
critical values, NOT a full-sheet empty double value.

The source and target have units and omitted divisors. Q has poles
and is NOT a polynomial on the whole source. For k>=2, P has critical
points on x=0; for k=1 it is a submersion but (7) forbids a global mate.
These are NOT actual whole-plane Keller normalizations. They do NOT
refute PF123/124's empty-double-branch/even-zero profile or reopen the
published quartic plane realization queue.

For k=1 there is a stronger extension: Q is regular on U=D(A(T)),
including x=0, since

    Q=y*(1+T^n/(n+1))/A(T)^2.

The map U->WHOLE A2_(P,Q) is SURJECTIVE etale: (5) handles P!=0,
and at P=0 the retained x=0 line maps by Q=y. Its full finite
normalization is the familiar generalized DANIELEWSKI surface

    Z={G(T)=P^2*Q} subset A3_(P,Q,T).                 (12)

This is finite free of degree n+1 and smooth: at P=0, G(T)=0 has only
simple roots; at P!=0, the Q derivative is nonzero. Besides the removed
ramification curves A(T)=0, it adds n UNRAMIFIED A1 curves P=0 at the
nonzero roots T^n=-(n+1) of G. The branch curves are Gm's
P^2*Q=n*zeta/(n+1), not polynomial A1-normalized images. Source units
and missing curves still prevent whole-plane realization. The standard
Danielewski geometry/exact-volume constructions already appear in
LA1 and UADANVOL1; they are NOT claimed new. The new control is the
specified whole-plane submersion and its exact raw mark orders.

## 6. The Lie-algebra and generic-period shortcuts are separate gaps

Regeta, [arXiv1311.0232v1](https://arxiv.org/html/1311.0232),
Propositions3.3/3.6, Theorem4.1 and Corollary4.3, recovers Keller pairs
from affine Lie algebras and relates algebraicity to JC2. The complete
pertinent proofs in Sections2--4 were read in primary HTML. This is a
ROUTE AUDIT, not a new imported classification or equivalent-criterion
programme. Polynomial brackets do not prove algebraicity. Remark4.2
also gives a polynomial sl2 subalgebra with a non-locally-finite element.
No algebraic group action or complete flow follows merely from brackets.

A literature search also prompted the exact-curve question; general
meromorphic projections permit singular images. Our smooth embedded
family is proved directly, without importing a secondary discussion.
No new external mathematical theorem or primary PDF/cache is adopted.

These controls do NOT contradict LP's global [1]=0 iff a polynomial
mate exists. They show generic [1]=0, even with generic lambda=0 and a
trivial FINITE connection, cannot replace that RAW condition. They do
NOT contradict LS1, which assumes an ACTUAL Keller pair, or LR0/1's
all-irreducible frame, where raw torsion vanishes. For actual Keller
fibers the raw Liouville mark is already zero; this does not make all
generic H1 vanish. No whole-plane submersion with global mate and
nonzero H1 is established here.

## 7. Remaining objective and assurance

Bare exactness of a smooth algebraic curve cannot prove rationality.
Generic exactness cannot globalize a coframe across a reducible special
fiber. These are all-degree negative closures of weaker arguments.

The full ACTUAL polynomial plane source ==> Tr(Lambda*Lambda) subset B,
or full marked-energy integrality, remains OPEN. Binary affine and
natural alternating/symmetric quotient compatibility remains OPEN,
including PF124's end/even-zero cases. No intermediate plane or square
root of a function follows from an even zero divisor.

New symbolic controls transcribe (1),(4),(6),(8)--(11), weight/exponent,
cyclic genus, and raw torsion identities. They do NOT replay registered
mathematical checkers or certify universal geometry by a finite census.
Independent replay, formal verification and external review remain false.
