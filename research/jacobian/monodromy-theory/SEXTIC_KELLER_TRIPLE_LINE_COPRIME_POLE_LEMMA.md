# A Laurent pole bound for centered mates of coprime degrees

**Lemma KellerCenteredCoprimeLaurentPoleBound.** Let \(k\) be a field
of characteristic zero and \(D=-x^2\partial_x\) on \(k[x,x^{-1}]\).
Let \(P,Q\in k[x,x^{-1}][Z]\) be monic of coprime degrees
\(m>n\ge2\) in \(Z\), with \([Z^{n-1}]Q=0\). Suppose that
\[
                    D(P)Q_Z-P_ZD(Q)=\kappa/x,\qquad \kappa\in k^\times .
\]
Then \(m=n+1\) and
\[
      Q=Z^n+\sum_{i=0}^{n-2}q_iZ^i+\eta/x,\qquad
      q_i,\eta\in k,\qquad \eta^2=-n\kappa/(n+1).
\]

Work at zero or infinity with uniformizer \(t=x\) or \(t=1/x\).
Divide the bracket by \(x\), so the derivation is
\(-x\partial_x=\pm t\partial_t\), and the right side is
\(\kappa/x^2\). This derivation preserves each Laurent
exponent, and acts by a nonzero scalar on every nonzero
exponent. If \(Q\) has any pole in its coefficients, put
\[
       \omega=\max_{1\le i\le n}
                 \frac{\operatorname{poleord}_t[Z^{n-i}]Q}{i}>0 ,
\]
where zero coefficients are omitted. Then \(\omega\ge1/n\).
For each \(j=1,\ldots,m\), the coefficient at
\(Z^{m+n-j-1}\) in the bracket contains \(n\) times the
Euler derivative of \([Z^{m-j}]P\). Every other term
uses a coefficient of \(Q\) of drop \(i>0\) and a
coefficient of \(P\) of smaller drop \(j-i\), along
with one Euler derivative. Induction gives pole order at
most \(j\omega\) for \([Z^{m-j}]P\).
Here the exponent of the equation is positive for all
\(j\le m\), and the Euler derivative cannot annihilate
the leading term of a coefficient with a pole.

To read the terms at pole order \((m+n-1)\omega\),
substitute \(Z=t^{-\omega}S\) and multiply the bracket
by \(t^{(m+n-1)\omega}\), allowing a common root of
\(t\) for fractional exponents. The leading normalized
polynomials \(P_*,Q_*\) are monic and \(Q_*\) has no
\(S^{n-1}\) term. Each coefficient of drop \(i\)
contributes to them only if its pole order is exactly
\(i\omega\). Thus a leading Euler derivative is, up to
sign, \(\omega(mP_*-S P_*')\) or
\(\omega(nQ_*-S Q_*')\), respectively. The leading bracket
is consequently
\[
                  \pm\omega(mP_*Q_*'-nP_*'Q_*).
\]
If \((m+n-1)\omega\) exceeds the pole order of the right
side, this expression vanishes. Then
\(P_*^n=Q_*^m\) by differentiation and monicity.
Coprimeness makes \(Q_*\) an \(n\)-th power of a monic
linear polynomial. Centering makes it \(S^n\), contrary
to the maximum defining \(\omega\).

The right side has no pole at infinity and has order two
at zero. Hence \(Q\) has no coefficient poles at infinity.
At zero a pole is possible only if \(m=n+1\) and
\(\omega=1/n\). If \(Q\) has no pole there either, its
coefficients are constant in \(x\). The same top
coefficient induction now gives no poles at either end
for \(P\), so its coefficients are constant as well.
This contradicts the nonzero bracket.

At the remaining threshold, only the drop \(n\) can
contribute a pole to \(Q\), and it is simple. Thus
\(Q_*=S^n+\eta\) with \(\eta\ne0\). Only drop \(n\)
can reach equality in \(P\) as well, so
\(P_*=S^{n+1}+\delta S\); all other terms, including
any constant term of \(P\) of pole order at most one,
are strictly below the leading threshold.
At zero the leading bracket is
\[
                     (n\delta-(n+1)\eta)S^n-\delta\eta .
\]
It equals \(\kappa\), giving
\(\delta=(n+1)\eta/n\) and the stated equation for \(\eta\).
The absence of poles at infinity and the bounds at zero
give exactly the Laurent support of \(Q\) in the statement.

## Sextic application

**Theorem SexticTripleLineShearSexticQuinticExclusion.** A
three-variable Keller map of ordinary total degree at most
six in characteristic zero with a quartic component \(g(x)+x^2(z+xy)\)
cannot have free degrees \((6,5)\) in the shear to \(r=x+x^2v\).

Use the [triple-line support and leading equation](SEXTIC_KELLER_TRIPLE_LINE_CUBIC_MATE.md#leading-coefficients-and-equal-degrees).
Both leading coefficients are univariate in \(x\) of degree
at most six, and their ratio \(p^5/q^6\) is constant because
no full fibre factor can occur. Unique factorization gives
scalar multiples of \(\ell^6,\ell^5\) with \(\ell\) linear.
The formal residue identity for the fifth root of \(q\)
requires a rational primitive of \(dx/(x^2\ell)\) over
\(k\), precluding a nonzero root of \(\ell\).
Rescale the mates to \(p=x^6,q=x^5\) and set
\(Z=xy+[y^4]Q/(5x^4)\).
Then the hypotheses of the lemma hold with \(m=6,n=5\),
and \(\kappa\) is a nonzero scalar. Since the centered
coefficients of \(Q\) belong to
\(\mathbb C[r,x,x^{-1}]\), the invariant \(q_i\) are
polynomials in \(r\), and \(\eta\) is a nonzero scalar.
At the original \(y=0\), if \(Z\) is regular at \(x=0\)
the simple pole \(\eta/x\) persists. If \(Z\) has a pole,
its fifth power dominates all other terms.
Both contradict polynomiality of \(Q\).
