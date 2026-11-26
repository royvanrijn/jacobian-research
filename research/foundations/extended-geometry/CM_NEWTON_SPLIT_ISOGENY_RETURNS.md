# Short isogeny returns at a split failed pair

Let \(A/\mathbb F_p\) be an abelian threefold with integral real
multiplication by the integers of a totally real cubic field \(F\).
Suppose that \(p\) is odd and split
in \(F\), and exactly one of the three coefficient pairs is nonordinary.
Write \(\mathfrak q\) for its coefficient prime, \(\pi\) for the Frobenius,
and \(X^2-aX+p\) for the relative Frobenius polynomial. Put
\[
 D=a^2-4p,\qquad E=F(\sqrt D),\qquad (D)=\mathfrak r\mathfrak s^2,
\]
with \(\mathfrak r\) squarefree as an ideal of \(\mathcal O_F\). The
[one-failed-pair ramification calculation](CM_NEWTON_FROBENIUS_RAMIFICATION_AND_DENOMINATORS.md#inertia-and-the-frobenius-weyl-group)
gives absolute simplicity and a quadratic CM endomorphism field \(E/F\).
Indeed the Frobenius powers have degree six, so Honda–Tate gives the
same commutative rational endomorphism field over every finite residue
extension. In particular every geometric endomorphism is rationally
defined over \(\mathbb F_p\). It is then defined over \(\mathbb F_p\)
by torsion-freeness of the endomorphism group. The ideal \(\mathfrak r\) contains
\(\mathfrak q\) exactly once and no other coefficient prime over \(p\).
The prime \(\mathfrak q\) is ramified in \(E\), with unique prime
\(\mathfrak Q\) above it of norm \(p\). [Tate, Theorem 1](https://numdam.org/item/SB_1968-1969__11__95_0/).

The height-two, dimension-one \(\mathfrak q\)-part of the
\(p\)-divisible group has slopes \((1/2,1/2)\). It has a unique
subgroup \(H\) of rank \(p\). Over an algebraic closure it is
\(\alpha_p\), since a rank-\(p\) subgroup can be neither étale nor
multiplicative. It lies in the kernel of relative Frobenius, which
has rank \(p\) on this part. Uniqueness also gives descent of \(H\).
Write \(\phi:A\to A/H\) for the quotient. Consider isogenies
\(\psi:A/H\to A\) respecting \(\mathcal O_F\) and of degree prime
to \(p\). Polarizations and auxiliary level structures are omitted
from the isogeny condition. [Tate–Oort, Theorem 2](https://www.numdam.org/item/ASENS_1970_4_3_1_1_0.pdf).

## The short-return bound

**Theorem.** Any such \(\psi\) of degree \(m<p/64\) has trace-zero
composite \(\gamma=\psi\phi\) over \(E/F\). Conversely a nonzero
trace-zero endomorphism of degree \(mp\) with \(p\nmid m\) factors
in this way, with return degree \(m\). For a short return, the
positive squarefree part of \(m\) is
\[
 c_v=\operatorname{sf}\!\left(\frac{N_{F/\mathbb Q}\mathfrak r}{p}\right).
\]
Here \(\operatorname{sf}\) denotes the positive squarefree part of a
positive integer. This is also the rational discriminant cofactor in
the [split Frobenius square sieve](CM_NEWTON_FROBENIUS_DISCRIMINANT_SIEVE.md).

**Proof.** The degree of \(\gamma\) is \(N_{E/\mathbb Q}(\gamma)\),
since a rational Tate module away from \(p\) has rank one over the
corresponding coefficient algebra. The kernel of \(\gamma\) has
\(p\)-primary part exactly \(H\). In particular \(\gamma\) belongs
to \(\mathfrak Q\): otherwise \(\gamma\bar\gamma\) is a unit at
\(\mathfrak q\), and \(\bar\gamma/(\gamma\bar\gamma)\) inverts its
action on the \(\mathfrak q\)-part. Here \(\bar\gamma\) is integral
as an endomorphism because \(\bar\gamma=\operatorname{Tr}_{E/F}(\gamma)-\gamma\)
and the trace belongs to \(\mathcal O_F\). Both \(\gamma\) and its CM conjugate
then reduce to zero at \(\mathfrak Q\), giving
\(\operatorname{Tr}_{E/F}(\gamma)\in\mathfrak q\). At infinity
\[
 \left|N_{F/\mathbb Q}\operatorname{Tr}_{E/F}(\gamma)\right|
 \ \le 2^3\sqrt{N_{E/\mathbb Q}(\gamma)}=8\sqrt{mp}<p.
\]
An algebraic integer in \(\mathfrak q\) with norm of magnitude less
than \(p\) is zero. This proves the trace-zero claim.

A nonzero trace-zero \(\gamma\) is \(u\sqrt D\) for some \(u\in F\).
The valuation of \(\gamma^2\) at \(\mathfrak q\) is odd. It is at
least one because \(\gamma\) is integral. When its degree has
\(p\)-valuation one, the entire \(p\)-primary kernel is therefore
in the \(\mathfrak q\)-part and is precisely \(H\). Thus \(\gamma\)
factors through the quotient by \(H\), with a prime-to-\(p\) second
isogeny. Finally \(m=N_{E/\mathbb Q}(\gamma)/p\) and the negative
norm of \(D\) have the same rational squareclass after removing
\(p\), since \(\gamma^2=u^2D\). The norm relation
\(-N_{F/\mathbb Q}(D)=N_{F/\mathbb Q}(\mathfrak r)N_{F/\mathbb Q}(\mathfrak s)^2\)
identifies that squareclass with \(c_v\). \(\square\)

The integer \(N_{F/\mathbb Q}(\mathfrak r)/p\) can have rational square
factors even though \(\mathfrak r\) is a squarefree ideal. This occurs
when several of its coefficient prime ideals lie over the same rational
prime, or when a coefficient prime has even residue degree. A loop,
meaning an \(\mathcal O_F\)-linear isomorphism \(A/H\cong A\), forces
\(c_v=1\) for \(p>64\). More strongly the norm of \(\mathfrak r\) is
then exactly \(p\), by the degree formula below.

## Exact minimization over the endomorphism order

Let \(\mathcal O=\operatorname{End}_{\mathbb F_p}(A)\subseteq\mathcal O_E\),
and denote absolute ideal norms in \(F\) by \(N\). The fractional
\(\mathcal O_F\)-ideal
\[
 \mathfrak b=\{u\in F:u\sqrt D\in\mathcal O\},\qquad
 \mathfrak I=\mathfrak b\mathfrak s
\]
is nonzero, and \(\mathfrak I\) is integral. Indeed
\(2\pi-a=\sqrt D\) belongs to \(\mathcal O\). For any \(u\in F^\times\),
the element \(u\sqrt D\) is integral over \(\mathcal O_F\) if and only
if \(u^2D\) is integral, by the monic quadratic for the converse.
Since \(\mathfrak r\) is squarefree, this is equivalent to
\((u)\mathfrak s\) being integral. The possible ideals
\(\mathfrak a=(u)\mathfrak s\) for trace-zero endomorphisms are therefore
exactly the integral ideals contained in \(\mathfrak I\) and in the
class \([\mathfrak s]\). The degree of the associated endomorphism is
\[
                 N(\mathfrak r)N(\mathfrak a)^2.
\]

For an ideal class \(C\) of \(F\) let \(\mu(C)\) be the minimum norm
of an integral ideal in that class. Writing
\(\mathfrak a=\mathfrak I\mathfrak c\), the minimum degree of a
nonzero trace-zero endomorphism is
\[
        d_0=N(\mathfrak r)N(\mathfrak I)^2\mu([\mathfrak b]^{-1})^2.
\]
The finite constant \(C_F=\max_C\mu(C)\) depends only on \(F\).
The index of \(\mathcal O_F[\pi]\) in \(\mathcal O_E\) is prime to
\(p\). At the failed coefficient place the discriminant \(D\) has
valuation one. At each other coefficient place the local quadratic
splits with distinct roots modulo \(p\), so the discriminant is a
unit there. These valuations exclude an index at all places over
\(p\). Thus \(\mathfrak I\) is prime to \(p\) because \(\mathcal O\)
is already maximal at those places and \(\mathfrak b\) for
\(\mathcal O_E\) is \(\mathfrak s^{-1}\).

**Theorem.** For \(p>C_F\), there is a return isogeny of degree
\(m_0=d_0/p\) with trace-zero composite, and this is the minimum
return degree among trace-zero composites. There exists a return of
degree less than \(p/64\) exactly when \(m_0<p/64\), in which case
\(m_0\) is the minimum degree of any prime-to-\(p\) return. On a
maximal-order representative of the isogeny class one has
\[
 m_0=\frac{N(\mathfrak r)}p\mu([\mathfrak s])^2.
\]

The formula follows because every smallest ideal representative in
\(F\) has norm at most \(C_F<p\), so \(d_0\) has \(p\)-valuation
one. The short-return theorem applies in both directions. For the
maximal order \(\mathfrak I=\mathcal O_F\). A maximal-order
representative can be reached by an isogeny over \(\mathbb F_p\) of
degree prime to \(p\). In Waterhouse's construction take an integer
\(n\) with \(n\mathcal O_E\subseteq\mathcal O\) and \(p\nmid n\),
possible since the order is already maximal at \(p\). The kernel
associated to the ideal \(n\mathcal O_E\) lies in \(A[n]\), and
the endomorphism ring of the quotient contains \(\mathcal O_E\).
Thus the claimed prime-to-\(p\) passage needs no change of residue
field. [Waterhouse, Proposition 3.9 and Theorems 3.13 and 6.1](https://www.numdam.org/item/ASENS_1969_4_2_4_521_0.pdf).

For a maximal-order representative there is also an ideal-class
description without a degree cutoff. The minimum return degree is
the smallest norm of an integral ideal \(J\) of \(\mathcal O_E\)
prime to \(p\) with \([J]=[\mathfrak Q]^{-1}\). In one direction any
return composite has principal ideal \(\mathfrak QJ\). In the
other direction a generator of \(\mathfrak QJ\) gives a self-isogeny
with \(p\)-primary kernel \(H\). This exact description allows the
generator to have nonzero relative trace when the degree is large.

## Consequences for a fixed threefold

**Corollary.** For the reductions of a fixed real-cubic threefold with full geometric
endomorphism algebra \(F\), the locus of split one-failed-pair primes
admitting a return of degree at most \((\log p)^\delta\) has relative
density zero for every \(0<\delta<1/18\). Under the GRH assumptions of
the split square sieve, the same holds with \(p^\delta\) for
\(0<\delta<1/38\). These conclusions follow because such returns are
short for large \(p\) and force \(c_v\le m\). They apply even if the
return exists only on a maximal-order representative. In particular
the special case of an isomorphism between the two endpoints of the
labelled \(p\)-isogeny has relative density zero.
[Unconditional cofactor bound](CM_NEWTON_UNCONDITIONAL_COFACTOR_SIEVE.md)
· [Conditional split cofactor bound](CM_NEWTON_FROBENIUS_DISCRIMINANT_SIEVE.md).

The known class-group expansion gives a larger degree bound. Under
GRH, the bound of Jao–Miller–Venkatesan for prime ideal
generators of norm at most \(x=(\log T)^B\), with \(B>2\) and
\(T=|\operatorname{Disc}(E)|\), bounds the ratio of nontrivial to
trivial adjacency eigenvalues by
\((\log T)^{1-B/2+o(1)}\). Their resulting sufficient walk length
to a prescribed vertex is of order
\(\log h(E)/((B/2-1)\log\log T)\). Bounding \(h(E)\) by
\(T^{1/2+o(1)}\) and taking the product of the degree bounds on those
edges gives at best \(T^{B/(B-2)+o(1)}\) from this argument. For large
\(p\), the generators in this bound are prime to \(p\), since
\(T\ll_F p^3\). That archimedean discriminant bound permits this
return estimate to exceed \(p\) by a power. The input still needed
for an ordinarity conclusion is arithmetic control on this particular
return class across the Frobenius fields of one fixed threefold.
[Jao–Miller–Venkatesan, Theorem 1.1 and Lemma 2.1](https://arxiv.org/pdf/0811.0647).
