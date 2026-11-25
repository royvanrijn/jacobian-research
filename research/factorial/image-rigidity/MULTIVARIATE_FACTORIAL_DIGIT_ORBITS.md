# Coprimality under complete base-\(p\) exponent orbits

At any terminal \(N\), a complete orbit of common exponent shifts has
an exact reduction modulo a prime \(p>N\). The reduction depends on
the exponent gap modulo \((p-1)\operatorname{lcm}(N-1,N)\).
In particular it gives consecutive coprimality for fully shared
profiles with both orientations and with unbounded exponents, gaps,
and numbers of unequal rows.

## Exact orbit reduction

Fix \(N\ge2\), a prime \(p>N\), and a finite list of nonnegative
integer pairs \((A_g,B_g)\). Put \(d_g=B_g-A_g\). For each group \(g\)
choose \(L_g=p^{h_g}\), with \(h_g\ge1\), such that
\[
 L_g\equiv1\pmod{\operatorname{lcm}(N-1,N)},\qquad
 N\max(A_g,B_g)<L_g .
 \tag{1}
\]
These powers exist because \(p\) is coprime to both orders.
Use the exponent rows \((A_g+i,B_g+i)\) for every \(0\le i<L_g\)
in group \(g\), and let \(M_j(t)\) be the ordinary factorial moment
polynomial of the resulting profile. Define
\[
 H_j(u)=(j!)^{\lfloor u/j\rfloor}(u\bmod j)!\ \in\mathbb F_p^\times
 \quad(u\in\mathbb Z),
 \qquad
 D=\sum_g d_g\frac{L_g-1}{p-1}.
\]
The remainder \(u\bmod j\) is in \(\{0,\ldots,j-1\}\) and negative
powers are taken in \(\mathbb F_p^\times\).

**Theorem.** For \(j=N-1,N\), the polynomial
\[
 \widehat M_j(t)=\frac{M_j((-p)^{-D}t)}{M_j(0)}
\]
has \(p\)-adically integral coefficients, with coefficientwise reduction
\[
 \overline M_j(t)=\sum_{k=0}^j\binom jk
               \prod_g H_j(kd_g)\,t^k\quad\text{in }\mathbb F_p[t].
 \tag{2}
\]
Both reductions retain their degrees. If they are coprime, the
original moment polynomials are coprime over every
characteristic-zero field.

For the proof, fix one group, write \(L=p^h\), and let
\(J=jA+kd\). Its factorial arguments are \(x_i=ji+J\) for
\(0\le i<L\), with \(0\le J<L\) and \(0\le x_i<pL\).
For each \(1\le e\le h\), the residues \(x_i\bmod p^e\) are
uniformly distributed, since \(j\) is coprime to \(p\). Thus
\[
 \sum_i \left\lfloor\frac{x_i}{p^e}\right\rfloor
 =\frac{\sum_i x_i-\sum_i(x_i\bmod p^e)}{p^e}
\]
is affine in \(J\), with slope \(L/p^e\). All floors for \(e>h\)
vanish. By Legendre's formula the change in total factorial valuation
from \(k=0\) to \(k\) is consequently \(kd(L-1)/(p-1)\), with
the same slope in \(k\) at both orders.

Write \(V_p(x)=v_p(x!)\), and let \(x_e\) be the base-\(p\)
digits of \(x\). The factorial unit is
\[
 p^{-V_p(x)}x!\equiv (-1)^{V_p(x)}\prod_{e\ge0}x_e!\pmod p .
 \tag{3}
\]
Indeed, in \(x=pq+r\) the product of factors prime to \(p\)
contributes \((-1)^q r!\) by Wilson's theorem, and the remaining
unit is the unit part of \(q!\). Iterate this identity. The unit
factorials for the lower \(h\) digits of all \(x_i\) have product
independent of \(J\), by the same uniform residue count.
The top digit is \(q_i=\lfloor (ji+J)/L\rfloor\).
As \(J\) rises by one from \(J-1\), precisely one \(q_i\) increases.
For that index \(ji+J=Lq_i\), so \(1\le q_i\le j\) and
\(q_i\equiv J\pmod j\) by (1). Hence the product of the top digit
factorials is multiplied by the representative of \(J\bmod j\)
in \(\{1,\ldots,j\}\). From \(J=0\) to any \(J<L\), the ratio of top
digit products is \(H_j(J)\). The change from \(J=jA\) to \(J=jA+kd\)
is \(H_j(jA+kd)/H_j(jA)=H_j(kd)\), with the same formula also
when \(d<0\).

Multiplying over groups gives total valuation difference \(kD\)
and the sign change \((-1)^{kD}\) in (3), both removed by the
common substitution \((-p)^{-D}t\). Division by \(M_j(0)\)
removes the constant coefficient and its unit. Since \(j<p\)
the binomial coefficients are \(p\)-adic units, and this proves (2)
including preservation of degree. Coprimality of (2) then proves
that the resultant of the two normalized polynomials is a
\(p\)-adic unit. Their rational Bézout identity transfers to
every characteristic-zero field.

## Consequences for unbounded shared profiles

Put \(T=(p-1)\operatorname{lcm}(N-1,N)\). For both orders
\(H_j(u+T)=H_j(u)\) by Fermat's theorem. Thus (2) gives a
finite-field criterion involving only the multiset of \(d_g\bmod T\).
For example, one may allow any number of groups whose gap is zero
modulo \(T\), since each contributes one to every coefficient of
both reductions. For fixed \(N,p\), the multiplicity of every
gap class can also be reduced modulo \(p-1\).

Suppose there is exactly one further group of gap \(1\bmod T\).
Then \(H_j(kd)=k!\) for \(0\le k\le j\). Dividing (2) by \(j!\)
gives \(R_j(t)=\sum_{k=0}^j t^k/(j-k)!\). The identity
\(R_N(t)-t R_{N-1}(t)=1/N!\) proves coprimality.
Reversing both exponent columns gives the same conclusion with
the single marked gap congruent to \(-1\).

There is also coprimality when there are two further groups
of gaps \(1\) and \(-1\) modulo \(T\). At order \(j\) their
combined contribution is \(k!(j-k)!/j!\). The reductions in (2)
are then \(\sum_{k=0}^j t^k\), which satisfy
\(\overline M_N-t\overline M_{N-1}=1\). These families allow
all \(A_g,B_g>0\) and arbitrarily large positive and negative gaps.
All conditions are compatible with (1) by taking the powers
\(L_g\) large. To make every exponent row distinct across the
groups, choose their intervals \([A_g,A_g+L_g-1]\) successively
disjoint, increasing both \(A_g,B_g\) by the same amount before
choosing \(L_g\). Increasing both exponents in a group leaves
its gap unchanged.

More generally, take \(r\) groups of gap \(1\) and \(s\) groups of gap
\(-1\) modulo \(T\), and any number of gap-zero groups. Their reduced
pair is the normalized reduction for two monomials using disjoint
variable sets of cardinalities \(s,r\), with each exponent one.
For \(r,s\ge1\), the ordinary moments of the seed profile are
coprime at both orders by the
[disjoint-support result](MULTIVARIATE_TWO_MONOMIAL_NONDIVISIBILITY.md#all-orders-when-both-private-exponent-sums-are-at-least-two)
at every terminal. For fixed \(N,r,s\) the modular criterion therefore
holds for every sufficiently large prime \(p\), by excluding primes
dividing the nonzero rational seed resultant or its denominators.

## Independent shifts and dilations of the exponent rows

The orbit theorem is stable under large common shifts of the two
exponents in each row, even when the shifts differ between rows.
More generally fix any two-monomial profile and a prime \(p\).
Independently for each row \((a_i,b_i)\), choose integers
\(e_i,H_i,C_i\ge0\) and replace the row by
\[
 (p^{e_i}a_i+p^{H_i}C_i,\ p^{e_i}b_i+p^{H_i}C_i),\qquad
 p^{e_i}N\max(a_i,b_i)<p^{H_i}.
 \tag{4}
\]
Let \(M_j^\star\) be the new moment polynomial and put
\[
 D^\star=\sum_i (b_i-a_i)\frac{p^{e_i}-1}{p-1}.
\]
Then the coefficients of
\[
 \frac{M_j^\star((-p)^{-D^\star}t)}{M_j^\star(0)}
 \quad\text{and}\quad
 \frac{M_j(t)}{M_j(0)}
 \tag{5}
\]
have identical \(p\)-adic valuations and identical unit parts
modulo \(p\), at both orders \(j=N-1,N\).
In particular this operation preserves a coprime full-degree
reduction after any additional common substitution and separate
constant scalars already used on the original pair.

To prove (5), set \(x=(j-k)a_i+kb_i\) in a single row. The transformed
argument is \(p^{e_i}x+p^{H_i}jC_i\). Its two displayed summands
occupy disjoint digit positions by (4). Legendre's digit-sum formula
and (3) give
\[
 \begin{aligned}
  V_p(p^{e_i}x+p^{H_i}jC_i)
  &=V_p(x)+\frac{p^{e_i}-1}{p-1}x+
       \frac{p^{H_i}-1}{p-1}\,jC_i+V_p(jC_i),\\
  \prod_{\ell\ge0}
  \big(\operatorname{digit}_{p,\ell}(p^{e_i}x+p^{H_i}jC_i)\big)!
  &=
  \left(\prod_{\ell\ge0}\big(\operatorname{digit}_{p,\ell} x\big)!\right)
  \left(\prod_{\ell\ge0}\big(\operatorname{digit}_{p,\ell}(jC_i)\big)!\right).
 \end{aligned}
\]
After dividing by the \(k=0\) factorial coefficient, the valuation
and sign changes in (3) are exactly the \(i\)-th contributions
to \(kD^\star\) and \((-1)^{kD^\star}\). The substitution in (5)
removes them.

Applying (4) to the coprime orbit families permits arbitrary
rowwise \(p\)-power scale factors and independent common offsets.
Taking those shifts large preserves full shared support and lets
any prescribed subset of rows contribute arbitrarily little to
the curvature \(\sum_i (b_i-a_i)^2/\max(a_i,b_i)\), with unchanged
gap in each row when \(e_i=0\). Thus the modular criterion does
not require a curvature lower bound.

The finite-field test concerns the actual consecutive factorial
polynomials before any real normalization. Its hypothesis of complete
orbits restricts the support. The argument proves no exhaustion of
arbitrary supports by those orbits.
