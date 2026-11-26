# Simultaneously passing all odd-corner divisor tests

Fix an even order \(m=2q\ge6\), let \(n=m+1\), and take any \(\tau>0\). There are ordinary two-monomial profiles in the simultaneous envelope corner at orders \(m,n\) with signed degree difference zero, relative gaps tending uniformly to zero, and curvature tending to \(\tau\), with the following two properties.

1. Every profile satisfies physical SFC(2) in every degree.
2. Every linked quadratic anchor in the [odd-corner divisor test](MULTIVARIATE_ODD_CORNER_LINKED_QUADRATIC_ANCHORS.md) passes its necessary divisibility, simultaneously and at full prime-power depth. This remains true when, before forming the anchored matrices, each translated polynomial has every common factor in the translation parameter removed and is made primitive over the integers.

Thus even the joint anchored divisibility cannot be contradicted uniformly in the moving row data. The normalized anchor polynomials are unchanged up to a nonzero scalar at each specialization. At least one determinantal content in the joint test is therefore still nonzero. Here primitive normalization is performed on the translated polynomials over the integer polynomial ring. Division by the numerical coefficient content separately at each integer translation is a different operation.

## Choosing the heights and moduli

Take an integer \(J>\max(n,3)\) at least the [absolute angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md). Write \(\mathcal F\) for the finite list
\[
       (s,h),\qquad s=m,n,\qquad 1\le h\le 3q-s .
\]
For \(B\) a large positive integer, the corresponding primitive anchor modulus is
\[
       V_{s,h}=(sB+h)/g_{s,h},\qquad g_{s,h}=\gcd(s,h).
\]
Choose a distinct fixed prime \(\lambda_{s,h}\) for each anchor, all greater than \(J\) and a fixed constant \(L\) to be chosen. The Chinese remainder theorem gives an infinite progression of integers \(B\) for which
\[
                 V_{s,h}\equiv\lambda_{s,h}\pmod{\lambda_{s,h}^2}
                            \qquad((s,h)\in\mathcal F).
\]
The coefficient \(s/g_{s,h}\) is a unit in each of these congruences. Every prime power \(p^E\) exactly dividing one of the moduli is then at most \(2JB/L\) for all sufficiently large \(B\). For \(p\ne\lambda_{s,h}\), use \(p^E\le V_{s,h}/\lambda_{s,h}\); for the designated prime the exponent is one and the prime stays fixed.

Start the profile with \(qR\) increasing unit-gap rows and \(R\) decreasing gap-\(q\) rows at smaller exponent \(B-1\), where
\[
                         R=\left\lfloor\frac{\tau B}{q(q+1)}\right\rfloor .
\]
These rows have signed degree zero and curvature tending to \(\tau\). The other rows below have bounded gaps and smaller exponents between fixed positive fractions of \(B\), strictly below \(B-1\). Their increasing upper exponents are at most \(B\) and their decreasing gaps satisfy \(2d_i\le m(B-a_i)\) for large \(B\). Thus they preserve the simultaneous corner and preserve \(B-1\) as the maximum smaller exponent.

## Private primes certify the actual profiles

Choose one witness row of gap one at each terminal \(2\le N<J\), with smaller exponent \(a_N\in[B/8,B/6]\) and prime \(p_N=N a_N+1\). Require that no other row chosen so far has a factor divisible by \(p_N\) in either full stripped ramp at orders \(N-1,N\). Also require the new row to avoid the earlier witness primes in the respective pairs of orders.

These choices exist successively for large \(B\). The progression \(N a+1\) over the given interval contains \(\gg_J B/\log B\) primes. This follows from the fixed-modulus prime number theorem with error \(O_N(x/(\log x)^2)\), or from the corresponding ordinary estimate when \(N=2\). The already chosen rows have a bounded number of distinct types and bounded gaps. Their full stripped factors are \(O_J(B)\). Each such factor has at most one prime divisor greater than \(B/8\) for large \(B\), so the new prime has only \(O_J(1)\) blocked values from these rows. Each previously chosen witness prime has size at least \(B/8\) and exceeds \(J\), so each factor of the new row in the tested orders excludes only one residue class for its smaller exponent. The interval contains \(O(1)\) integers in each such class. Keep the chosen smaller exponents separated by more than \(2J\) as well, which excludes only a bounded number of further bases. [Bennett, Martin, O'Bryant and Rechnitzer, Theorem 1.3 and (1.15)](https://personal.math.ubc.ca/~gerg/papers/downloads/EBPAP.pdf)

Assign the witness orientations with counts differing by at most one. If needed, add one opposite unit-gap row in the same interval, avoiding all witness primes and keeping the same spacing. This makes their aggregate signed degree zero. All the subsequent rows are required to avoid the witness primes in their respective pairs.

At the pair \(N-1,N\), the factor \(p_N=N a_N+1\) is the only nonunit of the witness row in order \(N\). It appears at every coefficient except one endpoint. At the neighboring order, a factor \((N-1)a_N+h\) cannot be divisible by \(p_N\), because
\[
                  N((N-1)a_N+h)-(N-1)p_N=Nh-(N-1)
\]
is nonzero and has bounded magnitude. The other rows avoid the prime by construction. The binomial coefficients are units since \(p_N>J\). The marked Newton polygon is a single edge with slope \(1/N\) or \(-1/N\), and the other polygon is horizontal. The Newton polygon product rule rules out any nonconstant rational common factor. The absolute angular theorem covers \(N\ge J\), and the resulting coprimality holds over every characteristic-zero field. [Filaseta and Rocks, Theorem 2.1](https://www.impan.pl/shop/publication/transaction/download/product/91656)

## Rows enforcing every anchor divisibility

For each \((s,h)\in\mathcal F\) and each \(p^E\) exactly dividing \(V_{s,h}\), let \(t\) be the other order of \(m,n\). Add two rows of gap \(t\), one per orientation, at smaller exponents in \([B/4,B/3]\). In each row require a factor in the first block at order \(t\) divisible by \(p^{E+v_p(t)}\). Keep the smaller exponents separated by more than \(2J\), and avoid all witness primes as above.

To verify that these requirements can all be met, put \(Q=p^{E+v_p(t)}\). A multiple \(y\) of \(Q\) in the range of order-\(t\) factors determines the smaller exponent \(a=\lfloor(y-1)/t\rfloor\), with \(y\) in the first block \(ta+1,\ldots,ta+t\). If \(Q\le t\), every candidate base works. Otherwise the resulting bases are distinct and the candidate count is at least \(t(B/12-O(1))/Q-1\). We have \(Q\le 2J^2 B/L\), so this lower bound can be made larger than any fixed number depending on \(J\) by choosing \(L\) large.

There are \(O_J(\log B)\) tasks, since each \(V_{s,h}\) is \(O_J(B)\). At most two prime powers above \(\sqrt B\) occur in any one modulus for large \(B\), since they are powers of distinct primes. Place those tasks first. Their spacing and all witness-prime exclusions cost a bounded number of bases, and the choice of \(L\) leaves enough bases. The rest have \(Q\le J\sqrt B\), hence have \(\gg_J\sqrt B\) candidate bases. At most \(O_J(\log B)\) bases are excluded by spacing to the already placed rows, plus a bounded number for the witnesses. This proves existence. All added pairs have zero signed degree, and their combined curvature is \(O_J(\log B/B)\). The witness curvature is \(O_J(1/B)\). This proves the assertions about the limiting real invariants.

## Divisibility after removing polynomial content

Fix the resulting profile and translate it with the parameter \(X\) of the linked anchor theorem, so the actual profile is \(X=B\). At order \(j\), write each linear factorial factor as
\[
                 j(X-t_i)+r=\gcd(j,r)\,\ell_{j,i,r}(X),
\]
where the factor \(\ell_{j,i,r}\) is a primitive integer linear polynomial. Every stripped coefficient is an integer constant times a product of such primitive lines. Remove the common primitive lines with their minimum multiplicities over all coefficients. Then divide by the integer content of the remaining polynomial in \(X,z\). By Gauss's lemma, each coefficient is still an integer constant times the product of its remaining primitive lines. Indeed the product of the lines is primitive, so the integer content comes only from the integer constants.

The two rows added for a particular \(p^E\mid V_{s,h}\) each have a first-block primitive line of value divisible by \(p^E\) at \(X=B\). The extra \(v_p(t)\) used in choosing their factors suffices because \(\gcd(t,r)\mid t\). These lines cannot be removed as common polynomial factors. Each is absent at one endpoint of its own row and cannot occur in another row: the spacing of the added heights makes their full integer ramps disjoint from all other rows at every tested order. Distinct literal factors at \(X=B\) give nonassociate primitive lines, by normalizing their common slope to one. This reasoning does not require the top repeated types to have disjoint ramps.

At least one of the two opposite rows samples its first block at every coefficient index of the unmarked order \(t\). Thus every scalar-free unmarked coefficient at \(X=B\) is divisible by \(p^E\). At the anchor \(X_0=-h/s\), reduction modulo \(p^E\) is unchanged because \(p^E\mid (sB+h)/g_{s,h}\) and the reduced denominator \(s/g_{s,h}\) is a \(p\)-unit. The corresponding unmarked coefficients at \(X_0\) all lie in \(p^E\mathbb Z_{(p)}\).

In the padded Sylvester map of the two scalar-free parents, every column from the unmarked order thus lies in \(p^E\mathbb Z_{(p)}^{s+t}\). A minor of size \(s+t-1\) uses at least \(s-1\) such columns. Therefore their reduced numerators are all divisible by \(p^{(s-1)E}\). Writing \(\widetilde D_{s,h}\) for their gcd at this anchor, we have proved the stronger simultaneous condition
\[
                       V_{s,h}^{\,s-1}\mid \widetilde D_{s,h}
                             \qquad((s,h)\in\mathcal F).
\]
At each anchor the original marked polynomial has nonzero constant and the original unmarked polynomial has no zero coefficient. Consequently the removed common polynomial factors have nonzero value there. The two specialized scalar-free parents are nonzero constant multiples of the original specializations. The linked comparison still makes a selected \(D\) nonzero. This confirms that the divisibility passes at genuine nonzero contents even though private witness primes prove full ordinary coprimality.

## Exact checks

The [finite checker](verify_multivariate_odd_anchor_divisibility_obstruction.py) constructs a 48-row corner profile and independently checks its seven private-prime pairs below the test cutoff. It verifies all sixteen prime-power divisors of the five anchor moduli, 121 primitive coefficient divisibilities, 80 anchored minor divisibilities, and five exact Chinese remainder depths. Its minor computation reuses the [linked anchor checker](verify_multivariate_odd_corner_linked_quadratic_anchors.py). The asymptotic family and certification at all remaining orders follow from the proof above.
