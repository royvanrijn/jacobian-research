# Primitive power differences after lifting two targets

Lift two primitive triples to any common multiple of their targets.
Divide the opposite summands by their exact shared gcd before taking
the difference. This shift is independent of the common multiple.
Its relation to the primitive cross products involves the gcd of
the original targets. The resulting pair count applies even with
one representation per repaired target.

## Exact cancellation at a common target

Let \(A_i+B_i=C_i\) for \(i=1,2\) be distinct ordered positive
primitive triples with \(C_i\le X^{B_0}\). Put
\[
 g_A=\gcd(A_1,A_2),\quad g_B=\gcd(B_1,B_2),\quad
 g_C=\gcd(C_1,C_2),\qquad a_i=A_i/g_A,\quad b_i=B_i/g_B.
\]
The three \(g\)'s are pairwise coprime. Write \(P=a_1b_2\),
\(Q=a_2b_1\). The cross products are coprime. A shared prime
would divide the two reduced \(a_i\), the two reduced \(b_i\),
or both summands of one original row. Also \(P\ne Q\) since
a primitive positive row is determined by its odds \(A/B\).

**PD1.** Let \(H\) be any common multiple of \(C_1,C_2\),
and put \(B_i'=HB_i/C_i\). Then
\[
 \Delta=\frac{B_2'-B_1'}{\gcd(B_1',B_2')},\qquad
 \Delta g_C=g_A(P-Q),\qquad
 g_A\mid\Delta,\qquad g_C\mid P-Q.                           \tag{PD1}
\]
One may choose \(H\) as the least common multiple of the
two least-repaired targets of any order \(M\).

First take \(H_0=\operatorname{lcm}(C_1,C_2)\) and
\(J_i=H_0/C_i\). At every prime dividing \(H_0\), take a row
whose target has the maximal valuation. Both \(J_i\) and
\(B_i\) are units there. At other primes both lifts preserve
valuations, so \(\gcd(J_1B_1,J_2B_2)=g_B\). Now
\[
 J_2 B_2-J_1 B_1
   =\frac{C_1 B_2-C_2 B_1}{g_C}
   =\frac{g_Ag_B(P-Q)}{g_C}.
\]
Changing \(H_0\) to \(H\) scales the numerator and normalizer
by the same integer. The divisibilities follow from
\(\gcd(g_A,g_C)=1\).

Equivalently, if \(U_i=J_i B_i/g_B\) and \(A_i'=J_i A_i\),
the weighted difference on the common target is
\[
 U_2 A_1'-U_1 A_2'=\Delta H_0,\qquad
 \gcd(U_2 A_1',U_1 A_2')=J_1J_2 g_A.
\]
The quotient target is \(\Delta g_C/g_A\), since
\(H_0/(J_1J_2)=g_C\). This accounts for primes in just one
of the original targets. A least-repair lift induces the
same normalized equation because the total lift factor
in row \(i\) is \(H/C_i\).
[Least target repair](POWERFUL_TRIPLES_COMMON_DIVISOR_BOUND.md#a-reverse-transfer-by-repairing-only-the-target)
· [Descent of homogeneous divisibility](ABC_MINIMAL_REPAIR_PROJECTIVE_CONTENT.md#exact-descent-of-multihomogeneous-congruences)

## Counting with a shift dictionary

For \(q\ge2\) an integer and \(t\ne0\) rational, define
\[
 d_q(t;\mathcal P)=
     \prod_{\substack{p\notin\mathcal P\\
          v_p(t)\not\equiv0\bmod q}}p,
\]
where the prescribed pool \(\mathcal P\) has prime product
at most \(X^{B_1}\). Write \(\operatorname{rad}_{\overline{\mathcal P}}\)
for the ordinary radical outside this pool.

**PD2.** Suppose \(\Delta\) belongs to a dictionary of size \(F\),
and impose
\[
 d_q(P/Q;\mathcal P)\le D,\qquad
 \operatorname{rad}_{\overline{\mathcal P}}(g_B)\le G_B,\qquad
 \operatorname{rad}_{\overline{\mathcal P}}(g_C)\le G_C .
\]
The number of ordered primitive row pairs is
\[
                   O_{B_0,B_1,\varepsilon}
                       (F D G_B G_C X^\varepsilon).          \tag{PD2}
\]
The pool may depend on \(\Delta\). The order may be prescribed or
adaptive. A joint dictionary of \(\Delta\) and the exact
\(q\)-free coefficient pair may replace \(FD\) by its cardinality
at a fixed \(q\).

Fix \(\Delta\ne0\), then \(g_A\mid\Delta\) at divisor cost,
and \(g_B,g_C\) at the cost of their two radical caps and
supported-power enumeration. This fixes \(P-Q=\Delta g_C/g_A\).
Write \(P=d_1 x^q,\ Q=d_2 y^q\) with \(q\)-free \(d_i\).
Coprimality implies that outside the prescribed pool
the support product of \(d_1d_2\) is exactly
\(d_q(P/Q;\mathcal P)\). Enumerate that radical and then
the coefficients on the resulting fixed support at cost
\(O(DX^\varepsilon)\), with the inner allowances split.

For fixed \(q,d_1,d_2,P-Q\), the coefficient-uniform weighted
power-difference bound gives \(O(X^\varepsilon)\) root pairs.
It uses a binary norm and uniform unit spacing, including the
split case. As \(P,Q\le X^{2B_0}\), for \(q\) above a fixed
cutoff in terms of \(\varepsilon,B_0\) one can instead
enumerate \(x,y\le X^{2B_0/q}\) directly. Thus the estimate
is uniform in \(q\), and taking the order adaptive costs
at most a logarithmic factor. Above \(2B_0\log_2 X+1\)
all orders give the same coefficient pair and both roots
equal to one. Finally factor the cross products in
divisor-many ways to recover \(a_1,b_2,a_2,b_1\), and use
\(g_A,g_B\) to restore the original rows. Their exact repairs
are determined by those rows.
[Fixed support enumeration (EB5)](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)
· [Zeng, Lemma 3.1](https://hexagonmath.org/pdf/2611.00621v1#page=5)

## Modulus and real window completion without the shift record

**PD3.** Prescribe \(F_0\mid g_C\), \(P\le Y_1\), \(Q\le Y_2\),
and an interval of length \(\ell>0\) for \(\log(P/Q)\).
Keep the defect cap \(D\) and the gcd cap \(G_B\), and impose
\(\operatorname{rad}_{\overline{\mathcal P}}(g_A)\le G_A\).
Then the pair count is
\[
 O_{B_0,B_1,\varepsilon}\!\left(
 D G_A G_B X^\varepsilon\left[
  1+\frac{(Y_1Y_2)^{1/q}}{F_0}\min\{1,\ell/q\}\right]\right).
                                                               \tag{PD3}
\]
The pool may depend on \(F_0\), which can be one. The interval
may depend on fixed records. Without an interval replace
\(\min\{1,\ell/q\}\) by one.

Fix the coefficient class as before and set \(R_i=(Y_i/d_i)^{1/q}\).
A nonempty branch has \(R_i\ge1\). Primitivity makes \(P,Q\)
units modulo \(F_0\), so \((x/y)^q=d_2/d_1\bmod F_0\).
There are at most \((2q)^{\omega(F_0)}\) unit root residues
in the prime power unit groups. Use this at each fixed order
below the cutoff above. In one residue distinct primitive
pairs have nonzero cross determinant divisible by \(F_0\).
Scale the coordinates by \(R_1,R_2\). The scaled vectors
have length at most \(\sqrt2\), so their angles differ by
at least \(F_0/(2R_1R_2)\). The interval for \(\log(x/y)\)
has length \(\ell/q\), and the derivative of the rescaled
angle with respect to this logarithm is at most \(1/2\)
in absolute value. The angular range is therefore
\(O(\min\{1,\ell/q\})\). This bounds the number of root pairs
in each class by the bracket in (PD3). Restore \(g_A,g_B\)
by their radical caps and factor as in (PD2). Neither
\(g_C\) nor \(\Delta\) is needed as a record. For adaptive \(q\)
one can take the largest bracket over admissible orders.
[Angular packing for a fixed unit root residue](ABC_TWO_TARGET_SUMMAND_LOCKS.md#the-mixed-target-defect-count)

For example, if each pair has such a divisor \(F_0\ge X^\omega\)
of radical at most \(n^\gamma\), list \(F_0\) at cost
\(n^\gamma X^\varepsilon\), where \(n=X^L\).
With the actual height bounds \(Y_i\le X^{2B_0}\), the bracket
is bounded when \(q\omega\ge4B_0\). It is also bounded with
\(F_0=1\) in a prescribed interval of length at most
\(X^{-\omega}\) under that same numerical inequality.

## Consequence for a putative full primitive count

Take \(X\to\infty\) at each fixed \(L>0\). If a family has
\(T\ge n^d\) original rows and at least \(n^{-\tau}T^2\)
ordered distinct pairs satisfying (PD2) with
\[
 F\le n^\alpha,\qquad D\le n^\beta,\qquad
 G_B\le n^{\zeta_B},\qquad G_C\le n^{\zeta_C},
\]
then
\[
                  2d-\tau\le\alpha+\beta+\zeta_B+\zeta_C.     \tag{PD4}
\]
This follows by taking the inner allowance smaller than any
positive gap times \(L\). The same comparison uses
\(\beta+\zeta_A+\zeta_B+\gamma\) when (PD3) has bounded
bracket and a modulus dictionary of cap \(n^\gamma\).
The pairs may have distinct targets or repaired targets.

For a pair cohort of full exponent \(2/3-o(1)\) in a balanced
endpoint, one can take \(\zeta_B,\zeta_C\) arbitrarily small
at fixed margins by the prescribed-pool overlap estimate.
Fix the first row and use the product of its three letters
as the pool for the second row. At any fixed positive
overlap threshold, the discarded pair count has a fixed loss.
Thus a joint normalized shift and coefficient record below
\(2/3\) by a fixed power suffices in a full-exponent cohort.
Closure requires extracting such a cohort, or one of the
modulus or window tests with a small power-class defect,
from the original three radical caps.
[Prescribed-record overlap (DRR4)](ABC_DIVISOR_RECORD_RECIPROCITY.md#anisotropic-intersections-for-an-arbitrary-record)
· [Private primes in primitive secants](ABC_SECANT_ENDPOINT_TRANSFER.md#small-intersections-and-output-multiplicities)

## Why a diffuse degree law cannot supply the defect cap

Consider formal input degree data in which the primes carrying
\(A_1,B_1,A_2,B_2\) are all disjoint, and every exact positive
valuation on those supports is a prime integer. Suppose each of
the four radicals has normalized mass \(1/3+o(1)\) and the
largest exact slot on any of them has mass \(o(1)\). Then
\[
          \log_n d_q(P/Q;\varnothing)=4/3-o(1)
           \qquad\text{uniformly in integers }q\ge2.
\]
Indeed \(g_A=g_B=1\). At every prime of the four supports
the valuation of \(P/Q\) is one original valuation up to
sign. A prime integer is divisible by \(q\ge2\) only
when it equals \(q\). Thus an adaptive choice of the
order can omit at most one exact slot on each support.
A pool removing only \(o(1)\) radical mass has no effect
on the conclusion. Prime-index profiles with this slot
property can have the same diffuse law as fine-grid
profiles. This is a counterexample to an inference from
degree data and support separation. It does not construct
additive triples with these data.
[Prime-index obstruction to valuation grids (FV2)](ABC_FINE_VALUATION_DIVISIBILITY.md#prime-degree-indices-also-evade-sparse-signed-grids)
