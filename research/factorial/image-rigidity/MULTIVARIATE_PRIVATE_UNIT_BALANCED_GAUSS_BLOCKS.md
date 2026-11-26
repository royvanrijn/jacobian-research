# Balanced Gauss blocks near a mirrored private-unit profile

Consider an ordinary two-monomial profile with one private
unit row in each orientation, an arbitrary shared mirrored
core, and additional shared rows. A shared row has a smaller
exponent \(\ell\ge1\), a gap \(d>0\), and an orientation.
Its \(d\) Gauss slots are \((u,x)=(\ell/d,h/d)\) for
\(1\le h\le d\), with base \(mu+x\) in degree \(m\).

Partition the additional rows into blocks, each consisting
of two opposite lists of complete rows. In every block
assume the lists have the same total gap, number of rows,
and sum of their smaller exponents. A block may have one
common slope \(u\) or may contain several slopes.
Let \(J\ge4\) be any integer at least the
[computable absolute angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).

**Theorem.** Define
\[
 C_0=\sum_{\rm common\ slope\ blocks}
          \frac{\max_{\pm}\sum_{d\text{ in the }\pm\text{ list}}d^{-1}}
               {u^3}.
\]
For every other block let \(u_{\min}>0\) be its smallest slope,
and let \((u_0,x_0)\) be its common slot centroid within either
orientation. Define the total square cost and its maximum by
\[
 T(m)=\sum_{\rm other\ blocks}u_{\min}^{-3}
        \sum_{\rm slots\ in\ both\ lists}
                     (m(u-u_0)+(x-x_0))^2,\qquad
 C_2=\max\{T(2),T(J-1)\}.
\]
If either
\[
    C_0+4C_2\le\frac{1}{10^5J^4}
          \qquad\text{or}\qquad
    C_0+4C_2\le\frac{1}{2\cdot10^9J^2},                          \tag{1}
\]
then consecutive moment polynomials have disjoint root
arguments in every degree, and SFC(2) holds for the profile
over every characteristic-zero field. The reference can have
arbitrary curvature. Its two orientations are obtained by
keeping the mirrored core and copying one complete list of
each block to the other orientation.

The balance hypotheses mean that the two slot lists have
equal cardinality and equal sums of slopes and offsets.
Indeed a row contributes \(d\) slots with total slope
\(\ell\) and total offset \((d+1)/2\).
The comparison uses cancellation of the first coefficient
variation within each block, including blocks with
different slopes.

## A second derivative for the normalized slot chord

For an increasing slot of base \(A>0\), write
\[
     f_{m,k}(A)=\log(A)_k-(k/m)\log(A)_m .
\]
A decreasing slot contributes \(f_{m,m-k}(A)\).
These are the exact normalized log contributions from
the Gauss multiplication formula. The two private unit
factorials cancel the binomial coefficient.
For \(0<k<m\), twice differentiating gives
\[
 |f_{m,k}''(A)|=
 \frac{k(m-k)}{m}
 \left(\frac1k\sum_{i=0}^{k-1}\frac1{(A+i)^2}
       -\frac1{m-k}\sum_{j=k}^{m-1}\frac1{(A+j)^2}\right)
               \le\frac{k(m-k)}{A^3}.                          \tag{2}
\]
Indeed the difference inside the parentheses is an average
of \((A+i)^{-2}-(A+j)^{-2}\), at most \(2(j-i)/A^3\)
for each \(i<k\le j\). The average of \(j-i\) over these
pairs is \(m/2\).
Both endpoint chords vanish, and (2) also applies after
reversing the orientation.

## Common slope and offset grid cancellation

At a common slope \(u\), let \(\mu\) be the signed offset
counting measure of the two lists, one minus the other.
The cumulative count \(F(x)\) on \(0\le x\le1\) is
\[
                  F(x)=\sum_+\lfloor dx\rfloor-\sum_-\lfloor dx\rfloor .
\]
The total gaps agree, so the linear terms \(dx\) cancel.
The equal numbers of rows then give an exact first primitive
\[
 G(x)=\int_0^x F(y)\,dy
    =\frac12\sum_+\frac{\{dx\}(1-\{dx\})}{d}
      -\frac12\sum_-\frac{\{dx\}(1-\{dx\})}{d}.
\]
Indeed the integral of \(\lfloor dy\rfloor-dy\) from zero
to \(x\) is \(-x/2+\{dx\}(1-\{dx\})/(2d)\).
In particular \(G(0)=G(1)=0\), and
\(\lvert G\rvert\le(1/8)\max_\pm\sum_\pm 1/d\).
For a slot chord with base \(mu+x\), integrate twice by
parts against \(\mu\). The mass is zero and there is no
atom at zero, so the first integration gives
\(-\int_0^1 f'(mu+x)F(x)\,dx\), with the atom at one
included in the counting measure. The second gives
\(\int_0^1 f''(mu+x)G(x)\,dx\).
Using (2) and summing gives the common-slope part of the
normalized log coefficient error at most
\[
                         \frac{k(m-k)}{8m^3}\,C_0 .              \tag{3}
\]

## Blocks at different slopes

Within one of the other blocks, the equal slot counts and
coordinate sums give
\[
       \sum_+(mu+x-(mu_0+x_0))
                         =\sum_-(mu+x-(mu_0+x_0))=0.
\]
Taylor expansion of the appropriate slot chord around
\(mu_0+x_0\) therefore cancels the constant and linear
terms in the difference between the lists. All intervening
bases in the remainder are at least \(mu_{\min}\).
Thus (2) bounds the error by
\[
    \frac{k(m-k)}{2m^3 u_{\min}^3}
       \sum_{\rm slots\ in\ both\ lists}(m(u-u_0)+(x-x_0))^2.
\]
Sum this bound over the blocks to obtain
\(k(m-k)T(m)/(2m^3)\).
The function \(T(m)\) is a sum of squares of affine
functions with positive weights, so \(T(m)\le C_2\)
on \(2\le m\le J-1\).
Together with (3) we have
\[
                  |\Delta\log c_{m,k}|
                   \le \frac{k(m-k)}{8m^3}(C_0+4C_2).          \tag{4}
\]

## Direct application of the angular reserve

Both profiles have normalized coefficients in \((0,1]\)
after the private factorial cancellation, so (4) bounds
the absolute coefficient differences too. The endpoints
agree exactly. Use either reference disk radius in the
[uniform private-unit stability proof](MULTIVARIATE_PRIVATE_UNIT_MIRRORED_UNIFORM_ANGULAR_RESERVE.md#stability-and-a-matched-slot-application),
\(\rho_4=1/(2\cdot10^5J^4)\) or
\(\rho_2=1/(4\cdot10^9J^2)\).
For \(m<J\), the reference modulus on the disk boundaries
exceeds \(\rho/8\), and \((1+\rho)^m<3\).
On a boundary the evaluation error is strictly below
\[
       \frac{3(C_0+4C_2)}{8m^3}
                         \sum_{k=1}^{m-1} k(m-k)
                        < \frac{C_0+4C_2}{16}
\]
if \(C_0+4C_2>0\), and is zero otherwise.
For the sum use \(m(m^2-1)/6\).
Condition (1) therefore makes the evaluation error
smaller than \(\rho/8\) for the corresponding radius.
Rouché pairs the roots in these disks, and the two
argument errors total at most \(4\rho\), strictly below
the corresponding reference cross-degree gap.
Terminal two follows from real positivity in even degree,
and terminals \(N\ge J\) from the absolute theorem.

For a common-slope example, take gaps \(1,3\) in one
orientation and \(2,2\) in the other at integer slope \(u\).
No row in the block has its exact mirror. Here \(C_0=4/(3u^3)\)
and \(C_2=0\), so \(u^3\ge(8/3)10^9J^2\) suffices
by the quadratic reserve.
For a distinct-slope example, take gap-one rows at slopes
\(u-1,u+1\) in one orientation and both at \(u\)
in the other, where \(u\ge2\) is integral. Again there are
no exact mirror rows within the block.
Here \(C_0=0\) and \(C_2=2(J-1)^2/(u-1)^3\), so
\((u-1)^3\ge16\cdot10^9J^2(J-1)^2\) suffices.
Either block may be adjoined to any mirrored private-unit
core.
