# Centering higher signed tanh moments about the matched second moment

Let \(\epsilon_i\in\{-1,1\}\), \(0\le T_i\le b<1\),
\(U=\sum_i\operatorname{artanh}T_i\), and
\(M_p=\sum_i\epsilon_iT_i^p\). The \(T_i\) may include all Gauss factors
of an ordinary exponent profile, including private rows.

**Lemma.** For every \(q>0\) and \(0\le c\le1\), set
\[
 \rho_q(c)=\max\left(1-c,\,
           \frac{q}{q+1}c\left(\frac{c}{q+1}\right)^{1/q}\right).
\]
Then
\[
   |M_{q+2}-c b^q M_2|\ \le\ \rho_q(c)b^{q+1}U.
\]
In particular the rational pairs \((c,\rho)\) in the table give the
same inequality with \(\rho\) in place of \(\rho_q(c)\).

| \(q\) | \(c\) | \(\rho\) |
| ---: | ---: | ---: |
| \(2\) | \(3/4\) | \(1/4\) |
| \(4\) | \(2/3\) | \(1/3\) |
| \(6\) | \(5/8\) | \(3/8\) |

For \(b=0\) everything vanishes. Otherwise put \(x_i=T_i/b\).
The maximum of \(x^{q+1}-cx\) on \([0,1]\) is \(1-c\). Its minimum
is attained at \(x=(c/(q+1))^{1/q}\), with magnitude equal to the
second term in the maximum. Consequently
\[
  |T_i^{q+2}-c b^q T_i^2|
      \le b^{q+1} T_i \rho_q(c)
      \le b^{q+1}(\operatorname{artanh}T_i)\rho_q(c).
\]
Sum these inequalities. For \(q=2,c=3/4\) the two extrema have
the same magnitude. For \(q=4,c=2/3\) and \(q=6,c=5/8\),
the required upper bounds on the minimizing argument are \(5/8\)
and \(7/10\), respectively. They follow from
\(2/15\le(5/8)^4\) and \(5/56\le(7/10)^6\).

In a grouped logarithmic model the signed second moment is
\(P=M_2/2\). For \(q=2l-2\), the exact contribution from \(M_{2l}/(2l)\)
can use \(c b^{2l-2}P/l\) in the \(P\)-direction with error at most
\(\rho_q(c)b^{2l-1}U/(2l)\). The bound may be applied at each member
of a consecutive pair using its own cap. The matched transport of
\(P\) continues to use the [excess outer loss theorem](MULTIVARIATE_TANH_EXCESS_LOSS_TRANSPORT.md).
