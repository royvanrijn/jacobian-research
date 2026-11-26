# Conductor-18 reconstruction for a moving branch ratio

**Theorem (quartic odd polynomial and radial even term).**

Let \(q_\theta(z)=z^2-\theta z-1\) have distinct roots, and put
\[
 f=t^{18}z^{-1}q_\theta(z)^4
     +t^{17}(b+b_1z+b_2z^2+b_3z^3+b_4z^4)+t^{16}(-z^{-1}+c),
       \qquad z=x/t^2.
\]
For every fixed algebraic \(\theta,b,b_1,b_2,b_3,b_4,c\) and every
sufficiently large prime \(p\equiv7\pmod{18}\), at least one of the
unmarked moments at \(m+kp\), \(0\le k\le11\), is nonzero, where
\(m=(13p-1)/18\). If \(b_4=0\), shifts \(0\le k\le10\) suffice.
If \(b_2=b_3=b_4=0\), shifts \(0\le k\le7\) suffice.

This normalization covers every unordered pair of distinct nonzero
roots. For roots \(a,d\), choose \(\rho^2=-ad\), set \(z=\rho w\),
and put \(\theta=(a+d)/\rho\). The quartic powers of the two
root factors become \(\rho^8 q_\theta(w)^4\).
The angular coordinate scaling \(x\mapsto\rho x,\ y\mapsto\rho^{-1}y\)
preserves the sphere functional, and multiplication of \(f\) by
\(\rho^{-7}\) preserves the vanishing of its moments. Under these
two operations the prescribed negative lift has coefficient
\(-1\) at weight \(16\), and the other permitted coefficients stay
in the stated family.

**Corollary (characteristic-zero coefficient fields).**
The nonvanishing windows of the theorem apply to every fixed phase
over a characteristic-zero field, with a prime threshold allowed to
depend on the phase. For each sufficiently large rational prime in
the stated class, at least one moment in the window at that prime
is nonzero in the original field.

Indeed, generate a finitely generated \(\mathbb Q\)-algebra by the
parameters and invert the discriminant \(\theta^2+4\). A maximal
ideal gives a number-field specialization by Zariski's lemma.
The theorem applies to that one specialization. A moment that
vanished in the original field would still vanish after
specialization, since the sphere functional takes rational values
on monomials. If the roots were not normalized, one can work
first in the algebraic extension adjoining \(\rho\), or specialize
the root data while inverting \(ad(a-d)\). The resulting
nonvanishing descends to the original field.

## The coherent line after the first row

Use the [two-pole residue formula](CM_NEWTON_FULL_WIDTH_STRIP_MOMENTS.md#a-reconstruction-matrix-for-unmarked-moments).
Work at a good place above \(p\), with bars indicating reduction.
Write \(4m=2p+\ell\), where \(\ell=(8p-2)/9\). The coherent
coefficients are extracted from
\[
 q_{\bar\theta}(z)^{4m}
     = q_{\bar\theta^{\,p}}(z^p)^2q_{\bar\theta}(z)^\ell.
\]
Since \(2\ell<m+2p\) and \(0<m<p\), there are only two
potential contributions from the latter factor to these
coefficients, at \(m\) and \(m+p\). Call them \(a_0,a_1\).
The coherent generating polynomial is therefore
\[
 \sum_j A_j z^j = (a_0+a_1 z)\,q_{\bar\theta^{\,p}}(z)^2 .
\]
For the two roots of multiplicity four the coherent source
dimension is \(1+\lfloor91/18\rfloor-2\lfloor52/18\rfloor=2\).
The [full coherent source criterion](CM_NEWTON_TWO_INPUT_STRIP_MOMENTS.md#a-full-coherent-input-space-of-arbitrary-rank)
gives a nonzero coherent input. Row \(k=0\) has normalized residue
\(A_0/13=a_0/13\). If it vanishes, \(a_1\ne0\), and the
coherent profile is \(a_1 z q_{\bar\theta^{\,p}}(z)^2\).

For deficit \(p\), let \(e_1,e_2\) count the factors of deficits
one and two. Their angle losses from the maximum \(7m\) are at
least \(3e_1+7e_2\). Since \(e_1+2e_2=p\) and \(e_1+e_2\le m\),
one has \(e_1\le2m-p\). The angle is at most
\(7m-7p/2+e_1/2\le8m-4p<2p\), so only \(C_0,C_1\) remain.
If \(\bar b=0\), row 1 has no odd term at nonpositive angle
and equals \(-2a_1/899\), which is nonzero at all sufficiently
large primes.

## Exact polynomial identity

On \(\bar b\ne0\) write
\[
 u=\bar b^{\,p},\quad l_j=\bar b_j^{\,p}/u,\quad
 v=\bar c^{\,p},\quad \eta=\bar\theta^{\,p},\quad
 r=u^2,\quad U=uC_0/a_1,\quad V=uC_1/a_1 .
\]
Divide the residues by \(a_1\). Their simultaneous vanishing
in rows 1 and 2 fixes
\[
 U=2/31,\qquad V=(-2l_1+16v)/31+3364\eta/1519 .
\]

In the quartic family, rows 3 through 5 solve for \(l_2,l_3,l_4\)
in that order. In the cubic family, rows 3 and 4 solve for
\(l_2,l_3\). A term in \(l_j\) in row \(j+1\) requires
one factor of angle \(j\) and \(j\) negative lifts to reach
angle zero. Thus it occurs only in the \(C_0\) term, and
its coefficient after substituting \(U=2/31\) is
\[
 \frac{(2/31)(j+1)(-2)^j j!}{\prod_{i=0}^j(29+16j+2i)} .
\]
This nonzero rational constant allows substitution without any
new restrictions on the parameters. After substitution the quartic
family uses rows 6 through 11 of weighted degrees 5 through 10.
The cubic family uses rows 5 through 10 of weighted degrees 4
through 9. In both cases the ring is \(\mathbb Q[r,l_1,v,\eta]\)
and the weights are \((2,1,1,1)\).

When \(b_2=b_3=b_4=0\), row 3 instead fixes \(r\) as a polynomial
in \(l=l_1,v,\eta\)
\[
 r=\frac{-6l^2+48lv}{427}
    +\frac{10092\eta l+3936\eta v}{20923}
    +\frac{1309476-865554\eta^2}{1401841}.
\]
After these substitutions rows 4 through 7 yield polynomials
\(S_4,\ldots,S_7\) in \(\mathbb Q[l,v,\eta]\) of respective
degrees \(3,4,5,6\). These definitions may include nonzero
rational scalings to make the polynomials primitive integral.

The [exact verifier](verify_cm_newton_moving_branch_ratio_certificates.py)
generates these rows by expanding the Laurent phase with odd part
\(u(1+\sum_j l_jz^j)\), even part \(-z^{-1}+v\), and top
\(z^{-1}q_\eta(z)^4\). For the coherent column it uses the
angular coefficients of \(zq_\eta(z)^2\) as weights.
For each of the two pole columns it divides its odd part by \(u\).
It then substitutes \(r=u^2\) into the resulting powers of \(u\).
Each reduced row has parity \((-1)^{k-1}\) under simultaneous
sign reversal of the weight-one variables.

In each case take all monomial multiples of the residual rows
with even parity up to the indicated weighted degree.

| Odd degree | Degree bound | Matrix dimensions | Selected rows and columns |
|---:|---:|---:|---:|
| 1 | 12 | \(340\times251\) | 245 |
| 3 | 16 | \(1912\times1484\) | 1427 |
| 4 | 18 | \(2548\times2199\) | 2078 |

Reduction modulo \(65521\) selects the rows and columns for a
rational linear solve. Multiplication by the selected rows in
the full integer matrix verifies a rational combination equal
to the constant \(1\) in each case. Thus the residual polynomials
have no common zero for any value of \(\eta\). All rational
constants in the calculation exclude only finitely many residue
characteristics. Reduction of the identity at all other
sufficiently large good primes proves the stated moment assertion.

The replay passed with SageMath 10.9 from the repository root.

```sh
sage -python research/foundations/extended-geometry/verify_cm_newton_moving_branch_ratio_certificates.py
```
