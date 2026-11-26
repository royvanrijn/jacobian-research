# Prime-power row packing at every factorial moment order

Let
\[
         f=\sum_\nu (u_\nu/\nu!)X^\nu,\qquad
         E_f(T)=\sum_{n\ge0}\Lambda(f^n)T^n/n!,
\]
where the distinct exponent columns are nonzero and the naturally
normalized coefficients \(u_\nu\) are integral over a finite extension
of \(\mathbb Q_p\), with \(v(p)=1\). A unit singleton is a column of the
form \(p^b e_i\) whose coefficient \(u_\nu\) is a unit.
Fix a precision \(p^a\), \(a\ge1\).

Choose two distinct rows of the exponent array and write their entries
as \(A_\nu,B_\nu\). Suppose \(B_\nu=0\) for every unit singleton and
\(p^{a-1}\mid A_\nu\) for every column. Replace the two rows by the one
row \(p^e A_\nu+B_\nu\), leaving the other rows unchanged. In the new
phase \(g\), use the same natural coefficients \(u_\nu\) with respect
to the **new** exponent factorials.

**Theorem.** For a computable \(e\) depending on the normalized phase,
the prime, and the precision, the packed support remains distinct and
\[
                      E_g(T)\equiv E_f(T)\pmod{p^a\mathcal O[[T]]}
\]
at every moment order simultaneously. The theorem includes all \(a=1\)
profiles with the stated condition on the low row. It applies to
arbitrary support size and needs no terminal cutoff.

The [finite-window packing theorem](SFC_PADIC_ROW_MERGER_INITIAL_FORMS.md)
separates the two rows by a digit gap sufficient for a chosen order
bound. Here the valuation of the natural factorial coefficient gives
an order-independent bound on the contribution of the low row at
the chosen precision.

## Uniform bound on the low row

Let \(w_\nu=\sum_i s_p(\nu_i)\), the total exponent digit sum.
For a count vector \(m=(m_\nu)\), its normalized coefficient
satisfies the [digit-weight estimate](FACTORIAL_MOMENT_SERIES_PRIME_POWER_RATIONALITY.md#degree-bounds-in-characteristic-p)
\[
 v\left(C(m)\prod_\nu u_\nu^{m_\nu}\right)
    \ \ge\ \sum_\nu\left((w_\nu-1)v_p(m_\nu!)
                                   +m_\nu v(u_\nu)\right).        \tag{1}
\]
One can omit zero coefficients from the support. If \(B_\nu>0\),
the hypothesis on the low row ensures that \(w_\nu>1\) or
\(v(u_\nu)>0\). The corresponding summand tends to infinity with
\(m_\nu\). Choose a positive integer \(L_\nu\) for which it is
at least \(a\) whenever \(m_\nu\ge L_\nu\). Put
\[
                  Y=\sum_{B_\nu>0}B_\nu(L_\nu-1).
\]
So every count vector with possibly nonzero contribution modulo
\(p^a\) has \(\sum_\nu B_\nu m_\nu\le Y\).

Choose \(e\ge a\) with \(p^e>\max_\nu B_\nu\) and
\(p^{e-a+1}>\max(Y,\max_\nu B_\nu)\). Every packed column then has
the original total digit weight, because the digits in the two
rows do not overlap within a column. Thus the same count bounds
apply to both profiles. Distinctness follows by reading the low
digits and the high digits of the packed row.

## The common natural normalization

Set \(x=\sum_\nu A_\nu m_\nu,\ y=\sum_\nu B_\nu m_\nu\) for a
count vector within the low-row bound. Define
\[
        Q(x,y)=\frac{(p^e x+y)!}{x!\,y!},\qquad
                       D_e=1+p+\cdots+p^{e-1}.
\]
The same factorial comparison as for finite windows gives
\[
 Q(x,y)=p^{D_e x}
     \prod_{s=1}^e\prod_{\substack{1\le h\le p^s x\\p\nmid h}}h
                              \prod_{h=1}^y(1+p^e x/h).
\]
The last product is one modulo \(p^a\). Indeed
\(v_p(h)\le e-a\) for all \(1\le h\le y\).
Since \(p^{a-1}\mid x\), each complete factor interval \(p^s x\)
is a multiple of \(p^a\). If \(\epsilon\) is the product of the
unit classes modulo \(p^a\), the first product has residue
\(\epsilon^{D_e x/p^{a-1}}\). These same statements hold for
each pair \(A_\nu,B_\nu\), by the choice of gap. After natural
normalization the ratio between the two native coefficients is
\[
                 Q(x,y)\Big/\prod_\nu Q(A_\nu,B_\nu)^{m_\nu}.
\]
All powers of \(p\) and \(\epsilon\) cancel because
\(x=\sum_\nu A_\nu m_\nu\). This is a unit congruent to one modulo
\(p^a\). Counts outside the bound contribute zero in both profiles
by (1). If the low row carries no singleton column, the same bounds can
be taken universally in the natural coefficient variables. With
nonunit singleton coefficients on the low row, the bounds instead
use the specified local valuations.

Consequently every certificate for all starts of a fixed positive
window length, obtained at this modulus from the rational moment
series of either support, transfers to the other support.
For \(a=1\) one can successively pack rows carrying no unit
singleton. This reduces to one row per variable carrying at least
one unit singleton, or to one row if there is no unit singleton.
When no unit singleton exists, (1) itself implies that the reduced
series modulo any fixed \(p^a\) is a polynomial.
Packing is a transfer of the specified modular data, so it leaves
characteristic-zero lifting of any surviving residue windows open.

## An obstruction when two unit singleton rows meet

For each prime \(p\), take the two linear terms \(uX+vY\)
with \(u,v\) units and pack their rows to the exponents \(p^e,1\),
where \(e\ge1\). In the natural coefficient of \(u v^{p^e}\)
at moment order \(p^e+1\), the original normalized coefficient is
one while the packed one is
\[
                         \binom{2p^e}{p^e}\equiv2\pmod p .
\]
All coefficients through order \(p^e\) agree modulo \(p\):
a mixed count there has low exponent \(<p^e\), and both pure terms
have the unit singleton coefficient one. Thus increasing the digit
gap postpones the first universal discrepancy but does not eliminate it.

There is also a discrepancy for any fixed residue values
\(\bar u,\bar v\ne0\), without treating them as indeterminates.
The original series is \(1/((1-\bar uT)(1-\bar vT))\).
The [explicit residue formula](FACTORIAL_MOMENT_SERIES_PRIME_POWER_RATIONALITY.md#degree-bounds-in-characteristic-p)
gives the packed series
\[
                       \frac{(1-\bar vT)^{p^e-1}}
                          {1-\bar uT-(\bar vT)^{p^e}} .
\]
At \(T=1/\bar v\) the packed denominator has value \(-\bar u/\bar v\),
whereas the original rational function has a pole.
This rules out equality of the naturally normalized series for
these unit rows at any finite digit gap.
