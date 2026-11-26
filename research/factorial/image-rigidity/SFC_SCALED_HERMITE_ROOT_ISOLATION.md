# Rational isolation of scaled Hermite roots

Let \(\operatorname{He}_j\) be the monic probabilists' Hermite
polynomial, specified by the coefficients of \(\exp(xu-u^2/2)\).
The normalized nodes are the roots of \(\operatorname{He}_j(\sqrt j\,x)\).
Use \(E_j\) for the splitting field over \(\mathbb Q\) of the polynomial
of the nonzero squared nodes, and \(F_j\) for the splitting field of
the nodes themselves.

**Theorem.** If \(j\ge10\) and \(L\) is the compositum of the fields
\(F_i\) for any finite collection of distinct indices \(i\ne j\), then
every nonzero squared node at \(j\) has degree \(\lfloor j/2\rfloor\)
over \(L\). Moreover, any choice of one nonzero node at each of
finitely many distinct positive indices is linearly independent
over \(\mathbb Q\).

The degree assertion gives isolation of a polynomial in a single
squared node, including when successive Gaussian moment conditions
have different means. The linear assertion gives the corresponding
test for relations among the means in fixed shifted windows.

## Square-node groups

Write \(j=2k+\varepsilon\) with \(\varepsilon\in\{0,1\}\).
By the Hermite coefficient formula,
\[
 \operatorname{He}_{2k+\varepsilon}(x)
     =x^\varepsilon\,(-2)^k k!\,L_k^{(\varepsilon-1/2)}(x^2/2).
\]
Thus \(E_j\) is also the splitting field of the displayed Laguerre
polynomial. The scaling \(x^2=jz\) relating the unscaled squares to
the squared nodes is rational. Further, \(F_j/E_j\) is obtained by
adjoining square roots, so its Galois group is an elementary abelian
\(2\)-group.

[Schur, §1(D) and §4](https://doi.org/10.1515/crll.1931.165.52)
gives irreducibility of the two square-variable polynomials.
[His theorem II on pp. 52–53](https://gdz.sub.uni-goettingen.de/id/PPN243919689_0165)
gives the full symmetric group at \(k>12\). For \(5\le k\le12\)
the following squarefree factorization degrees give alternating
group containment. In each entry \(p:(d_1,\ldots,d_q)\) means
factorization degrees of the monic integral polynomial
\(x^{-\varepsilon}\operatorname{He}_{2k+\varepsilon}(x)\) in
the variable \(x^2\), reduced modulo \(p\).

| \(k\) | \(\varepsilon=0\) | \(\varepsilon=1\) |
|---:|:---|:---|
| 5 | \(37:(2,3)\) | \(17:(2,3)\) |
| 6 | \(13:(1,5),\ 31:(1,2,3)\) | \(31:(1,5),\ 19:(1,2,3)\) |
| 7 | \(19:(3,4)\) | \(29:(2,2,3)\) |
| 8 | \(19:(1,2,5)\) | \(29:(3,5)\) |
| 9 | \(41:(1,3,5)\) | \(41:(1,3,5)\) |
| 10 | \(29:(3,7)\) | \(61:(3,7)\) |
| 11 | \(59:(4,7)\) | \(37:(1,1,2,7)\) |
| 12 | \(97:(1,1,3,7)\) | \(43:(1,2,2,7)\) |

The squarefree Frobenius cycle rule and taking powers give single
prime cycles of lengths \(2,3,3,5,5,7,7,7\), respectively.
The transitive groups in degrees \(5,7\) are primitive by primality.
At degree \(6\) the extra \(5\)-cycle gives primitivity, and at
degrees \(8,\ldots,12\) the displayed prime cycle does so.
Indeed, a prime cycle of length exceeding \(k/2\) would have to act
trivially both on any nontrivial block system and on each of its
blocks, since both their count and their sizes are at most \(k/2\).
[Jordan's theorem as stated in Jones, Theorem 1.1](https://arxiv.org/abs/1209.5169)
now gives \(A_k\) in each case, since every selected prime cycle
fixes at least three points. The table is replayed by
[the exact small-degree verifier](verify_sfc_scaled_hermite_isolation.py).

For every \(k\ge5\), the two square fields of indices \(2k\) and
\(2k+1\) are distinct. In fact
[Schur's discriminant formula, stated in Hajir, (1.1)](https://doi.org/10.4153/CJM-2009-031-6)
for the monic Laguerre normalization gives the discriminant ratio
of the two parities, modulo squares, as
\[
                    \prod_{\substack{2\le h\le k\h\text{ even}}}
                          (2h-1)(2h+1).
\]
Rescaling in the variable \(x^2\) changes a monic discriminant by
a rational square. The product comprises the odd integers from
\(3\) to \(4\lfloor k/2\rfloor+1\). A Bertrand prime strictly
between \(k\) and \(2k\) occurs exactly once, so their square
classes differ. Consequently the groups are either \(A_k\)
and \(S_k\), or both \(S_k\) with different quadratic subfields.
In either case the fields differ and neither contains the other:
in the case of different groups this also uses that \(S_k\)
has no quotient \(A_k\) for \(k\ge5\).

## Isolation in a compositum

Fix a target \(j\) with \(k=\lfloor j/2\rfloor\ge5\).
Build the compositum of the rival squared fields first,
beginning with \(E_{2k+1-\varepsilon}\) if the same-degree
rival is present, or with \(\mathbb Q\) otherwise.
Denote this starting field by \(B\). Since
\(E_j\not\subset B\) and \(E_j\cap B\) is Galois over
\(\mathbb Q\), \(\operatorname{Gal}(E_j B/B)\) is a nontrivial
normal subgroup of \(A_k\) or \(S_k\), and contains \(A_k\).
In the rest of the tower, each relative Galois group is a
normal subgroup of the Galois group of the field being added,
by restriction to that field and its Galois intersection with
the previous compositum. The remaining squared fields thus
contribute no composition factor \(A_k\): for different degrees
at least five their only nonabelian simple factor is an
alternating group of different order, and the squared fields
of degrees at most four are solvable. Finally adjoining all
the rival signed nodes adds an elementary abelian \(2\)-group.
Therefore \(\operatorname{Gal}(L/B)\) has no composition
factor \(A_k\), ruling out \(E_j\subset L\).
The group \(\operatorname{Gal}(E_jL/L)\) is a nontrivial
normal subgroup of \(\operatorname{Gal}(E_j/\mathbb Q)\), hence
contains \(A_k\) and is transitive on the \(k\) squared nodes.
This proves the degree assertion. A rational linear relation
with nonzero coefficient at a nonzero node of some index
\(j\ge10\) would put that node, and thus its square, in \(L\).

## Low indices

For \(j\le9\), use the following possible valuations of a
nonzero node, with \(v_p(p)=1\). The symbol \(\mathbb Z\) denotes
integral valuations. A dash denotes an index already removed.
Index \(1\) has no nonzero node.

| \(j\) | \(p=7\) | \(p=5\) | \(p=3\) |
|---:|:---:|:---:|:---:|
| 2 | \(\mathbb Z\) | \(\mathbb Z\) | \(0\) |
| 3 | \(\mathbb Z\) | \(\mathbb Z\) | \(0\) |
| 4 | \(\mathbb Z\) | \(\mathbb Z\) | \(1/4\) |
| 5 | \(\mathbb Z\) | \(-1/4\) | — |
| 6 | \(\mathbb Z\) | \(1/6\) | — |
| 7 | \(-1/3\) | — | — |
| 8 | \(1/8\) | — | — |
| 9 | \(0,1/6\) | \(0,1/4\) | \(-1/2,-2/3\) |

For \(p=7,\ j\le6\) and \(p=5,\ j\le4\) the unscaled
nonzero Hermite polynomials have \(p\)-unit discriminant, and
\(\sqrt j\) is a unit in an unramified extension. The other
entries are obtained from the lower Newton polygon of the
monic squared-node polynomial and division of its slopes by
\(-2\). Both checks are in the same verifier.

At each prime, the integral-valuation terms can be combined in an
unramified field, so their sum still has integral valuation or is
zero. Apart from the unramified discriminant tests above, the only
unit cases requiring this check are \(j=9\) at \(7\) or \(5\).
The unscaled squared polynomial at index \(9\) is
\[
            Z^4-36Z^3+378Z^2-1260Z+945.
\]
Its nonzero reductions are simple: the reductions are
\(Z^3(Z-1)\) at \(7\) and \(Z^2(Z^2-Z+3)\) at \(5\).
After passing to the unramified closure of the residue field and
lifting these simple roots, their square roots, including the
unit scaling by \(9\), also lift by Hensel's lemma because the
primes are odd.

Now the nonintegral valuation classes within each column of the
table are distinct, and rational coefficients shift valuations by
integers. At \(7\), a proposed sum with a nonzero coefficient at
\(7\) or \(8\) has a unique lowest valuation after grouping the
integral terms, so these coefficients must be zero. The same
argument at \(5\) eliminates indices \(5,6\); at \(3\) it eliminates
\(4,9\), with the remaining \(2,3\) terms unramified at \(3\).
Any coefficient multiplying a zero node may be disregarded. The
only remaining possible nonzero nodes are \(\pm1/\sqrt2\) and
\(\pm1\), which are linearly independent over \(\mathbb Q\).
This proves the low-index part of the theorem.

## Application to mixed ternary common translations

Take three distinct ordinary offsets \(b_1,b_2,b_3\) with
\(|b_1|=|b_2|=d\), \(|b_3|=d+k\), where \(k\) is a nonzero integer.
Fix a positive integral direction \(a\) and a starting order
\(m_0\ge1\). For all sufficiently large integers \(T\), the
translated support \(Ta+b_j\) has the actual three-term window
at orders \(m_0,m_0+1,m_0+2\), over every characteristic-zero
field and uniformly over its coefficients. The threshold is
per prescribed starting order.

Use [the tilted Gamma normalization](SFC_TERNARY_COMMON_TRANSLATION_CUBIC_OBSTRUCTION.md#exact-tilted-gamma-normalization)
with \(t=1/T\) and coefficients \(C_j\). A hypothetical sequence
of projective common zeros can be normalized by \(\max|C_j|=1\).
Write \(P=C_1+C_2\), \(Q=C_3\). After division of each order-\(m\)
polynomial by its common order factor \(m^d\), the mean leading
constant is \(P+m^k Q\). The other terms are \(O(\sqrt t)\) in
every fixed \(L^p\) norm, since the Gamma coordinates concentrate
at one and the exponents are fixed. If \(|P+m^kQ|/\sqrt t\)
were unbounded, dividing the vanished order-\(m\) moment by that
constant to the \(m\)-th power would give a limit of one.
Thus \(P+m^kQ=O(\sqrt t)\) at all three orders. Since \(k\ne0\)
and the orders are distinct, \(P,Q=O(\sqrt t)\). At the
coefficient limit \(C_1=-C_2\ne0\).

Passing to limits of \(P/\sqrt t\) and \(Q/\sqrt t\), the
rescaled polynomials converge in moments to a scalar mean
\(A+B m^k\) plus a centered Gaussian linear form of bilinear
variance \(W/m\), with
\[
          W=C_1^2\sum_\rho (b_{\rho1}-b_{\rho2})^2/a_\rho\ne0
\]
at the coefficient limit. This follows from the Gamma cumulant
formulas in the cited normalization and the linear Taylor
part at the vector of ones. The quadratic and higher remainders
vanish in each fixed norm after division by \(\sqrt t\).
Therefore the vanishing limit at \(m\) gives
\(\operatorname{He}_m(\sqrt m\, r_m)=0\), where
\(r_m=(A+B m^k)/\sqrt{-W}\) with one fixed square root.
These three nodes satisfy a rational linear relation with
all coefficients nonzero, obtained by eliminating \(A,B\).
At least one node is nonzero, because at least one order is
even and the even Hermite polynomial has nonzero constant term.
The theorem excludes the relation. It suffices to work over
\(\mathbb C\), since the original homogeneous coefficient
system is over \(\mathbb Q\), and a common projective point
over a characteristic-zero extension gives a point over
\(\overline{\mathbb Q}\) by affine charts and the Nullstellensatz.
