# Real and simultaneous modular counts for joint root dictionaries

Use the shell boxes, raw root caps of total \(O(\lambda)\), arbitrary
fixed root records of cost \(f\), and prime-position colouring from
[the modular dictionary count](ABC_CROSS_LETTER_ROOT_DICTIONARY_ENTROPY.md#one-divisor-and-two-free-letters).
Take the height limit at fixed budget and finite degree list. Suppose
the short letter is \(T\le X^{1-\tau}\) for fixed \(\tau>0\) and the
maximum is comparable to \(X\). The tall pair \(U,V\) then satisfies
\(\log(U/V)=O(X^{-\tau})\). Let \(W\) bound the sum of the conditional
dictionary exponents over all records. The count exponent satisfies
\[
                               D\le f+W/2+o(\lambda)                 \tag{RG1}
\]
if the groups partition the remaining positions on \(U,V\). The same
bound
\[
                               D\le f+W/2+o(\lambda)                 \tag{RG2}
\]
holds for groups partitioning the remaining positions on all three
letters if the record includes a full divisor \(M\) of \(U\) with
fixed positive numerical weight and fullness at least
\(\lfloor\kappa/\lambda\rfloor\) for fixed \(\kappa>0\). In this case
use the prescribed signed ratio \(T/V\) modulo \(M\).

Groups may meet different letters. Each either has bounded size or
lies in boundedly many scaled degree windows per letter with diameter
\(\rho_\lambda\to0\) uniformly in the groups. There is no lower scaled
degree condition on free positions. The error means any fixed positive
multiple of \(\lambda\) at sufficiently small budgets with the stated
uniform margins, after the height limit.

## Logarithmic moments

Fix a record. Write \(\phi_J\) for the signed monomial of the \(U,V\)
positions in group \(J\), with positive exponents for \(U\).
For (RG2) let \(\psi_J\) be the corresponding contribution to \(T/V\)
modulo \(M\). For (RG1) use the trivial character and \(\psi_J=1\).
Frozen factors enter the real and modular targets. Use the
nonnegative kernel
\[
 K_H(s)=\left(\frac{\sin(\pi Hs)}{\pi Hs}\right)^2
       =\int e^{2\pi i ts}\,d\mu_H(t),\qquad H=X^{\tau/2},
\]
where \(\mu_H\) has probability density
\(H^{-1}(1-|t|/H)\) on \([-H,H]\). This kernel is bounded below by a
positive absolute constant on the actual log strip. Multiply it
by the modular indicator in (RG2). Fourier expansion and, when
used, character orthogonality bound the count for this record by
a constant times \(\int\mathbb E_\chi \prod_J|A_J(t,\chi)|\,d\mu_H(t)\),
where
\[
 A_J(t,\chi)=\sum_{\boldsymbol x\in\mathcal D_J}
 e^{2\pi i t\log\phi_J(\boldsymbol x)}\chi(\psi_J(\boldsymbol x)).
\]
The targets contribute phases of modulus one. In (RG1) the two
tall values determine the triple. In (RG2) all unfrozen values
are counted.

For raw group caps \(\alpha_J\) of total \(A\), set
\(h_J=\lceil A/\alpha_J\rceil\). As in the modular count, freeze
zero-cap positions at inner subpower cost, and take fixed orders
using the positive limiting caps at fixed budget. Then
\(h_J\alpha_J\le 2A\) and \(\sum_J(2h_J)^{-1}\le1\).
In the \(2h_J\)-th moment, choose the first \(h_J\) vectors at cost
\(h_JW_J\) and collapse the second \(h_J\) to one product in each
coordinate. Their summed raw height exponent is at most \(2A+o_X(1)\).
Their ordered factorizations are divisor-bounded at fixed moment
order, and distinct original positions use disjoint prime colours.

In the expansion of the moment, take \(s\) to be the log ratio of
the two tall products. When \(|s|>X^{-\tau/4}\), the kernel is
\(O(X^{-\tau/2})\). The raw count of the entire second half is
\(X^{2A+o_X(1)}\), so this tail is negligible for sufficiently small
fixed budgets. In the near part, count the collapsed \(U,V\) products
in the real strip around the first half's tall ratio. For a bounded
group use (RUW1) with correction equal to one. For a window group
compress the original collapsed products in each window, with the
signed rational correction of height \(X^{O(\rho_\lambda)+o_X(1)}\)
and boundedly many roots of height \(X^{O(\lambda)+o_X(1)}\).
The exact signed monomial recovers every original product by prime
colour and original signed degree, so the exact compressed fibre has
size one. For any fixed chart parameter the wider-height condition
holds eventually, and the chart parameter can then tend to zero.
This costs \(o(\lambda)\). The empty \(U,V\) part gives one choice.
[Unequal-height real strip count](ABC_SHORT_HEIGHT_AND_REAL_GAP.md#a-real-strip-with-unequal-root-heights)

In (RG2) the tall products are now fixed. The collapsed \(T\)
products occupy a fixed unit residue modulo \(M\). Apply (FC2) for
bounded groups, or the overlap-compatible wide cofactor count and
the same exact-product recovery by colour for window groups.
This costs a further \(o(\lambda)\). Again an empty part gives one
choice. [Modular window compression](ABC_CROSS_LETTER_ROOT_DICTIONARY_ENTROPY.md#groups-of-roots-in-boundedly-many-shrinking-degree-windows)

The group norm is thus at most
\(X^{W_J/2+o(\lambda)/(2h_J)+o_X(1)}\), uniformly in groups and
records. Hölder and the reciprocal-order bound give (RG1) and (RG2).

## Entropy consequence

Binning group surprises as in (CG3), at a balanced saturation box
with a fixed positive numerical gap the tall pair has joint group
entropy \((L_U+L_V+o(\lambda))\log X\) for every allowed partition
of its positions. This uses (RG1) with \(f=0\), including when the
short letter has no cheap full divisor. Under a prescribed record,
the conditional binning of (CG5) gives
\(H(\mathcal R)+\frac12\sum_JH(\boldsymbol x_J\mid\mathcal R)
\ge (D-o(\lambda))\log X\) with the respective free positions and
hypotheses. [Joint and conditional entropy binning](ABC_CROSS_LETTER_ROOT_DICTIONARY_ENTROPY.md#joint-entropy-at-saturation)
