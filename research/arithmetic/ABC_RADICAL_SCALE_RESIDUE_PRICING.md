# Absolute residue pricing at radical-scale depths

The layerwise diagonal moment has a quantitative
version when the modulus is only a fixed power of
the radical scale. The remaining price depends on
the actual modular interval available for products
of the squarefree layer roots.

## RD1. Moment and list prices at finite depth

Use \(n=X^L\), with the height limit at fixed \(L\).
Let \(n^\kappa\le D\le X^C\) with fixed \(\kappa,C>0\).
Count integers \(T\le X^C\), units modulo \(D\),
with radical outside a prescribed polynomial-height
prime support at most \(n^r\). Prescribe a list
\(S\subseteq(\mathbb Z/D\mathbb Z)^\times\) of at
most \(n^b\) unit classes. The support may vary with
\(X\) subject to its uniform height bound. Fix its
exact supported part at height-subpower cost. Fix
also the exact part on valuations above \(J\), for
\(J\) large depending on \(L\) and a prescribed
positive normalized loss. The latter costs at most
\(X^{C/(J+1)+\varepsilon}\), by its radical and
supported-power enumeration. Both operations translate
the prescribed unit residue list by a known unit.

On one of the remaining exact valuation boxes write
\[
              T/\text{fixed part}=\prod_e x_e^e,
       \qquad x_e\le n^{\theta_e},\qquad R=\sum_e\theta_e .
\]
At fixed \(J,L\), use a finite mesh for the
normalized base-\(n\) logarithms of the roots.
Round each nonunit slot up to a positive mesh
endpoint, with mesh size at most \(\eta'/J\).
Then \(R\le r+\eta'\), and the maximum slot price
increases by at most one mesh step. Bin unit
slots separately and omit them from the moments.
There are finitely many moment lists when they
are chosen in terms of these boxes, independent
of \(X,D,S\). For any choice of
positive integers \(h_e\) with
\(\sigma=\sum_e(2h_e)^{-1}\le1\), the count per
box is at most \(X^\varepsilon\) times
\[
 n^{\,R/2+\sum_e(h_e\theta_e-\kappa)_+/(2h_e)
                         +b\max(1/2,\sigma)}.          \tag{RD1}
\]
The support and tail prices are added separately.
In the singleton case \(b=0\). For an arbitrary
nonempty auxiliary set \(\Theta\) in the unit group,
set \(B_1=\log_n |S\Theta|\) and
\(\gamma=\log_n(|S\Theta|/|\Theta|)\). Then the list
price in (RD1) can also be replaced by
\[
                  \sigma B_1+(1-2\sigma)_+\frac{\gamma}{2}.
                                                        \tag{RD2}
\]
All constants can depend on the fixed moment orders and
degree cutoff at the height limit. The statement needs
no distribution hypothesis on the original prime
supports modulo \(D\).

Indeed, average characters with probability counting
measure on the dual unit group. For
\(F_e=\sum_{x_e}\chi(x_e^e)\), expand its \(2h_e\)-th
moment. Each equality in the moment equates the
\(e\)-th powers of two unit products. For a fixed first
product the second lies in at most \((2e)^{\omega(D)}\)
unit classes, including primes dividing \(e\), by the
prime-power unit group structure. Its height is at
most \(n^{h_e\theta_e}\), so each class has
\(O(1+n^{h_e\theta_e}/D)\) possible positive integers.
Each such product has height-subpower ordered
factorizations, uniformly as the exponent list and
\(h_e\) are fixed. The first product and its ordered
factors cost at most \(n^{h_e\theta_e}\). Consequently
\[
        \|F_e\|_{2h_e}
               \ll X^\varepsilon n^{\theta_e/2+
                                  (h_e\theta_e-\kappa)_+/(2h_e)}.
\]
Overcounting the squarefree roots by the full individual
ranges and forgetting disjointness across slots is
permitted, since Fourier inversion is applied to that
larger product domain and counts its actual residues.
For the residue indicator use Fourier norm exponent
\(1/(1-\sigma)\). Its \(L^2\) norm is at most \(n^{b/2}\),
and its \(L^\infty\) norm is at most \(n^b\). Hölder gives
(RD1) by that \(L^2\) bound when \(\sigma\le1/2\) and by
interpolation otherwise. For (RD2) use the nonnegative
convolution majorant of RX1. Its Fourier norm at
\(L^1\) is at most \(n^{\gamma/2}\), and at \(L^\infty\)
it is at most \(n^{B_1}\). The majorant is at most one
pointwise, with sum \(|S\Theta|\), so Parseval bounds
its Fourier \(L^2\) norm by \(n^{B_1/2}\).
Interpolation on either side of \(L^2\) gives (RD2).
[Relative residue majorant (RX1)](ABC_RESIDUE_EXPANSION_AND_JOINT_RAY_COUNTS.md#the-relative-expansion-price-for-a-deep-absolute-modulus)
· [Absolute residue count for complete targets (AR1)](ABC_ABSOLUTE_RESIDUE_TARGET_TAILS.md#one-full-divisor-and-the-multiplicity-threshold)

## RD3. Fine spectra and growth inside a finite unit group

If the maximum outside exact-layer radical price tends
to zero in the height-first boxes, the limiting price
from (RD1) is at most the minimum of \(r\) and
\[
 {\cal C}_{\kappa}(r,b)=
 \begin{cases}
   (r+b)/2,                   &\kappa\ge r,\\
   r/2+rb/(2\kappa),          &r/2\le\kappa\le r,\\
   r-\kappa+b,                &0<\kappa\le r/2 .
 \end{cases}                                             \tag{RD3}
\]
A precise fixed-allowance statement is that for every
\(\eta>0\) there is a positive upper bound on the
maximum individual exact-layer price, depending on
\(r,\kappa,b,\eta\), that gives count
\(O_L(n^{\min(r,{\cal C}_\kappa(r,b))+\eta}
X^\varepsilon)\). Here take \(0\le b\le\kappa\) if
the list bound is on an exactly radical-scale
modulus \(D=n^{\kappa+o(1)}\). No such upper bound on
\(b\) is needed for the estimate when \(D\) is larger.

For the calculation on the nonempty box with \(R>0\),
set \(a_0=\min(1,R/(2\kappa))\) and choose
\(h_e=\lceil R/(2a_0\theta_e)\rceil\). Then
\(\sigma\le a_0\) and \(h_e\theta_e\le R/(2a_0)+\omega\)
where \(\omega=\max\theta_e\). The extra product loss
in (RD1) is at most
\(a_0(R/(2a_0)+\omega-\kappa)_+\). This gives (RD3)
after the arbitrarily small boxing and tail losses.
For \(r=0\) the fixed-support count suffices. The
plain bound \(r\) follows from enumerating the
squarefree outside radical and its supported powers.
One can also choose the moment scale with \(r+\eta'\)
in place of \(R\) to apply the same list price
uniformly across all boxes, using the trivial bound
if needed.

For an application, assume \(D=n^{\kappa+o(1)}\),
\(\kappa\ge r/2\), and let a dictionary of integers
as above have count \(n^{a-o(1)}\), with \(a>r/2\).
Let \(S\) be its actual unit residue list. Write
\(s=2a-r\) and \(a_*=r/(2\kappa)\). For every
fixed \(m\ge1\), let \(b_m=\log_n |S^m|\), with
\(b_0=0\). By (RD2) and the moment choice above,
\[
       2a_* b_{m+1}+(1-2a_*)_+(b_{m+1}-b_m)
                                    \ge s-o(1).
\]
For \(\kappa\ge r\) this is
\(b_{m+1}-(1-r/\kappa)b_m\ge s-o(1)\).
For \(r/2\le\kappa\le r\) it already gives
\(b_1\ge\kappa s/r-o(1)\). In both ranges,
\[
 b_m\ \ge\ \frac{\kappa s}{r}
                  \left(1-(1-r/\kappa)_+^m\right)-o(1).
                                                        \tag{RD4}
\]
All losses are in the fixed-moment height-first
convention. If the lower bound on \(\log_n D\)
only tends to \(\kappa\), first use a fixed lower
bound \(\kappa'<\kappa\) and then decrease the
margin \(\kappa-\kappa'\) in each fixed-order
inequality. At \(\kappa=r/2\) the choice of
\(a_0\) stays at most one and the product loss
tends to zero. To verify the recurrence, take
\(\Theta=S^m\)
in (RD2), so \(\gamma=b_{m+1}-b_m\) and
\(B_1=b_{m+1}\). Compare that bound on the complete
integer dictionary with its assumed size and
iterate the resulting linear inequality.
In particular at the occurrence exponent \(a=r\),
for every fixed \(\epsilon>0\) some fixed product
order yields \(|S^m|\ge n^{\kappa-\epsilon}\)
up to the fixed-budget height allowance. For
\(0<\kappa\le r\), (RD3) already forces the
original residue dictionary to have exponent
\(\kappa\) if the complete dictionary reaches \(r\).
For \(\kappa>r\) that conclusion concerns bounded
products as in (RD4). A residue list whose bounded
products stay in a smaller power-size subgroup
is therefore incompatible with occurrence in this
fine-layer sector.

## RD5. Radical-scale full locks in balanced cells

There is also a version independent of the fine-layer
assumption, useful for choosing original locks.
Choose \(h_e=\lceil R/\theta_e\rceil\). Then
\(\sigma\le1/2\) and \(h_e\theta_e\le R+\max\theta_e\).
Consequently the original half-price list inequality
\((r+b)/2\) holds when
\[
                  \kappa>r+\max_e\theta_e
\]
by a fixed margin, including the boxing and tail
allowances. In particular \(\kappa>2r\) suffices for
arbitrary original exact-root distributions.
These inequalities are on numerical depth in base
\(n\) and do not require absolute depth in base \(X\).

In a regular balanced family, the RB2 residue list
input in RB3–RB7 can therefore use pairwise coprime
original full divisors of modulus at least
\(n^\kappa\) for any fixed \(\kappa>2/3\), with small
cap tolerances. The chart price is charged as usual.
For example a role with a fixed positive normalized
mass on sufficiently high exact orders and with
individually cheap primes can give a lock of this
depth at an arbitrarily small fixed radical price:
choose a packet of mass \(p\) with small overshoot,
all its orders large enough that \(e p>\kappa\).
It is enough that those orders escape to infinity in
the unscaled exact index, even if they have vanishing
scaled numerical degrees \(Le\). On one chart the
usual cell-size and conditional entropy conclusions
follow because their proofs use only the uniform
half-price count, original degree caps and CRT.
This replacement is for the plain half-list price.
Relative expansion has the quantitative loss in (RD2).
Consequences involving a positive absolute
power \(X^\nu\), such as long numerical
reconstruction, still require the stated numerical
conditions for those consequences. This distinction
also applies to the original positive-multiplicity
regime, where all three locks can instead be deep.
[Balanced dictionary and register inequalities (RB3–RB7)](ABC_BALANCED_DICTIONARY_RESIDUE_CELLS.md)
