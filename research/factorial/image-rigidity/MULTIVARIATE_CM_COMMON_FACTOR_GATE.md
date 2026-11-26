# An exact CM-field gate for both-private physical common factors

Consider an ordinary two-monomial profile with private variables
on both sides, with arbitrary shared rows. Suppose a monic
irreducible rational polynomial \(H\) is common to the physical
moment polynomials at adjacent orders \(n,N=n+1\). Let
\(\alpha\) be a root, \(K=\mathbb Q(\alpha)\), and \(D=[K:\mathbb Q]\).
Assume \(K\) is a CM field, meaning a totally imaginary
quadratic extension of a totally real field. Assume also that
at every rational prime the roots of \(H\) have one common
valuation in an algebraic closure of \(\mathbb Q_p\).

**Theorem.** The possible field \(K\) is \(\mathbb Q(\zeta_M)\)
for one of the following orders. In the last column,
\(v_p(p)=1\), and the listed field discriminant is signed.

| \(M\) | \(D\) | Necessary lower order \(n\) | Field discriminant and valuation restriction |
|---:|---:|---|---|
| \(4\) | \(2\) | \(n\equiv5\pmod 8\) | \(-4\), with \(v_2(\alpha)\in 1/2+\mathbb Z\) |
| \(6\) | \(2\) | \(n\equiv7,8,9\pmod {12}\) | \(-3\), with \(v_3(\alpha)\in 1/2+\mathbb Z\) |
| \(10\) | \(4\) | \(n\equiv11,17\pmod {20}\) | \(125\), with \(v_5(\alpha)\in 1/4+(1/2)\mathbb Z\) |
| \(12\) | \(4\) | \(n\equiv17\pmod {24}\) | \(144\), with local ramification and residue degrees \(e=f=2\) at \(p=2,3\) |

There is no curvature threshold. In particular if the local
tests force equal valuations and exclude every row of the table,
a common factor must generate a non-CM field. The theorem
uses no bound on the number of shared rows or on the order.
It applies only to candidate factors with the indicated
field property and valuation hypothesis.

The uniform valuation hypothesis can be checked separately
from the table. If at each prime any two distinct allowable
shared valuations have reduced denominators with sum greater
than \(D\), a factor of degree \(D\) cannot use both. The
degree spent by each slope is a positive multiple of its
reduced denominator. The
[local field-discriminant ledger](MULTIVARIATE_COMMON_FACTOR_FIELD_DISCRIMINANT.md)
also imposes the field discriminants and local degrees in
the table. In that comparison one can use the actual
primitive squarefree arrays.

## Why the argument field is the whole field

By the [both-private angle boxes](MULTIVARIATE_PRIVATE_UNIT_ASYMMETRIC_ANGLE_BOXES.md),
every conjugate of a possible common root has argument in the
open right half-plane. Their arguments are distinct, from
simplicity and the one-root-per-argument-box statement.
For an angle \(\theta>0\) the common label \(k\) must satisfy
\[
                    \frac{(2k-1)\pi}{n}<\theta
                                      <\frac{2k\pi}{n+2}.
\]
The geometric-sum case has no common root.

Let \(c\) be the nontrivial involution of \(K\) over its
totally real subfield. For every complex embedding \(\sigma\)
of \(K\), \(\sigma(c\alpha)=\overline{\sigma(\alpha)}\).
Set \(q=\alpha/c\alpha\). Its complex conjugates all have
absolute value one. Its valuations at all finite primes
are zero by the equal-valuation assumption on the roots of
\(H\). Thus \(q\) is an algebraic unit. The polynomials whose
roots are powers of its conjugates have bounded integral
coefficients and hence range over a finite set. Two powers
of \(q\) are equal, so \(q\) is a root of unity.

Distinct arguments of conjugates in the right half-plane
give distinct doubled arguments modulo \(2\pi\). Thus the
\(D\) embeddings of \(K\) give \(D\) distinct values of
\(\sigma(q)\). Consequently \(K=\mathbb Q(q)\). If \(q\) has
order \(M\), its conjugates are exactly the primitive
\(M\)-th roots, and the positive arguments of the common
factor are \(j\pi/M\) for \(0<j<M/2\) coprime to \(M\).
For a primitive index \(j\) this overlap condition becomes
\[
                 M<n j\bmod 2M<2M-2j.                       \tag{1}
\]
For odd \(M\), \(j=(M-1)/2\) leaves no integer between the
two bounds. For \(M=4m\), the index \(j=2m-1\) forces the
residue \(4m+1\). If \(m\) is even, this gives
\(n\equiv2m-1\bmod8m\), which fails (1) at \(j=1\).
If \(m\) is odd, it gives \(n\equiv6m-1\bmod8m\).
The primitive indices \(j<2m\) then pass (1) exactly when
\(j\equiv1\bmod4\). For \(m\equiv1\bmod4,\ m\ge5\), the
primitive \(j=m-2\) fails. For \(m\equiv3\bmod4,\ m\ge7\),
the primitive \(j=m-4\) fails. This leaves \(M=4,12\).

For \(M=2m\) with \(m\) odd, \(M=6\) gives its three
listed residues by \(j=1\). Otherwise \(m\ge5\), and
\(j=1\) requires \(n\equiv2m+x\bmod4m\) with
\(1\le x\le2m-3\). The primitive index \(m-2\) requires
\(x(m-2)\equiv r\bmod4m\) with \(1\le r\le3\).
The primitive index \(m-4\) requires
\(r-2x\bmod4m\in\{1,\ldots,7\}\).
Since \(7-4m\le r-2x\le1\), either \(x=1,r=3\)
or \(x=2m-3,r=1\). Substitution in the first congruence
forces \(m=5\) in both cases. At \(M=10\) the two
residues are \(11,17\bmod20\), both passing the test at
\(j=1,3\). This recovers the
[primitive angle classification](MULTIVARIATE_MATCHED_TANH_TORSION_NEWTON_GATE.md#the-cyclotomic-orders-allowed-by-the-angle-boxes)
without assuming that all moduli are equal.

## Finite prime restrictions

The identity
\[
                         \alpha=(1+q)t,\qquad t\in K^+
\]
follows by complex conjugation and \(q=\alpha/c\alpha\).
Here \(K^+\) is the totally real subfield. All embeddings of
\(t\) are positive since both \(\sigma(\alpha)\) and
\(1+\sigma(q)\) have positive real part and the same argument
in \((-\pi/2,\pi/2)\).

For \(M=4,6\), \(t\) is rational. The norm of \(1+q\)
is respectively \(2,3\), and that prime is totally ramified
with index two. This gives the two half-integer conditions.
For \(M=10\), the norm is \(5\), the unique prime over \(5\)
has index four, and its restriction to
\(K^+=\mathbb Q(\sqrt5)\) has index two. Thus the valuation
of \(1+q\) is \(1/4\) and that of \(t\) is in
\((1/2)\mathbb Z\). For \(M=12\), \(1+q\) is a unit.
The cyclotomic local degrees at \(2,3\) are both \(e=f=2\).
These follow by splitting the conductor into its prime
power and prime-to-\(p\) parts. The cyclotomic discriminant
formula gives the signed values in the table. Equivalently
the defining polynomial tests use
\(\Phi_4,\Phi_6,\Phi_{10},\Phi_{12}\) in those four cases.

Thus the table remains effective without a common-modulus
hypothesis. For example its order-ten field is subject to
a quarter-slope test, whereas the separate
[equal-modulus torsion theorem](MULTIVARIATE_MATCHED_TANH_TORSION_NEWTON_GATE.md#the-remaining-norm-and-valuation-restrictions)
excludes order ten outright using equality of all moduli.
