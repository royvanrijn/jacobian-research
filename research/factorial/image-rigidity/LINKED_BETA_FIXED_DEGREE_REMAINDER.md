# An order-independent remainder for linked binary beta moments

Fix nonnegative integers \(a,b\), and put \(d=a+b+1\).
For the ordinary factorial phase
\[
        X^aY^b((s-1)X+sY)
\]
the linked moment polynomial at order \(n\), up to a
nonzero scalar, is
\[
 P_n(s)=\int_0^1[((s-t)t^a(1-t)^b)]^n\,dt .
\]
This follows by the beta change of variables, or by expanding
and integrating the monomials. Put
\(J_n=P_n/A_n\), \(A_n=\mathrm B(na+1,nb+1)>0\), so \(J_n\)
is monic. The beta change of variables is the one in the
[physical linked problem](canonical-density/moment_pair_obstructions.md#1-what-would-constitute-a-factorial-counterexample).

For \(n\ge1\), define polynomials \(C_k\) by
\(C_0=0\), \(C_1=1\), and
\[
 (k+3+nd)C_{k+2}
 =\big[(k+2+n(a+b))s+k+2+n(a+1)\big]C_{k+1}
         -(k+1+na)s C_k.                                     \tag{1}
\]
Set
\[
 E_{n;a,b}(s)=
      \sum_{h=0}^{b}(-1)^h\binom bh
                      \big[s C_{a+h}(s)-C_{a+h+1}(s)\big].
                                                               \tag{2}
\]
The common roots of \(J_n,J_{n+1}\) are exactly the common
roots of \(J_n,E_{n;a,b}\), and
\(\deg E_{n;a,b}\le a+b\), independently of \(n\).
The exact resultant identity is
\[
 \operatorname{Res}(J_n,J_{n+1})
    =\left(\frac{A_n}{(nd+2)A_{n+1}}\right)^n
       J_n(0)J_n(1)\operatorname{Res}(J_n,J_n')
                          \operatorname{Res}(J_n,E_{n;a,b}).
                                                               \tag{3}
\]
The prefactors in (3) are nonzero for every indicated profile.
This reduces coprimality to one polynomial with degree
depending on the fixed exponent pair.

## Proof of the elimination

Write \(p_s(t)=(s-t)t^a(1-t)^b\) and
\(I_k(s)=\int_0^1 t^k p_s(t)^n\,dt\). The integral of the
derivative of \(t^k t(1-t)(s-t)p_s(t)^n\) is zero. Its
quadratic coefficient equation is (1), with \(I_k\)
in place of \(C_k\). At a root of \(I_0=P_n\), this
recurrence gives \(I_k=C_k I_1\) for all \(k\).
If \(I_1\) were also zero, every moment against a
polynomial in \(t\) would vanish. Taking the polynomial
with conjugate coefficients to \(p_s(t)^n\) contradicts
\(\int_0^1|p_s(t)|^{2n}\,dt>0\). Thus \(I_1\ne0\) there.
Expanding \((1-t)^b\) in the expression
\(P_{n+1}=\int p_s(t)p_s(t)^n\,dt\) gives
\(P_{n+1}=E_{n;a,b} I_1\) at these roots. Equation (1)
shows \(\deg C_k\le k-1\) for \(k\ge1\), so (2)
has the stated degree bound.

Integrating the derivative of \(t(1-t)p_s(t)^n\) gives
\[
 (nd+2)I_1=(na+1-n(s-1))P_n-s(1-s)P_n'.
\]
For the term involving \(P_n'\), use
\(t(1-t)/(s-t)=s(1-s)/(s-t)+s-1+t\).
This identity holds also at \(s\in[0,1]\), by clearing
the apparent denominator against the integral's
\((s-t)^n\) factor. At a root of \(P_n\) it proves
\(I_1=-s(1-s)P_n'/(nd+2)\). Neither endpoint is a root,
and \(I_1\ne0\) there also proves simplicity of all \(J_n\)
roots. Multiplying the identity for \(J_{n+1}\) at these
roots yields (3), where the product of the roots is
\((-1)^nJ_n(0)\).

## An unequal positive exponent application

For \((a,b)=(1,2)\), the formula obtained from (1)--(2) is
\[
 E_{n;1,2}(s)=
 \frac{(n+1)(s-2)
       ((9n+6)s^2-(4n+2)s+4n+2)}
      {4(4n+3)(4n+5)} .
\]
Under the common-root chart \(z=s/(s-1)\), its transform
\((z-1)^3 E_{n;1,2}(z/(z-1))\) is a nonzero scalar times
\[
                (z-2)((9n+6)z^2-(4n+2)z+4n+2).
\]
The positive real possibility \(z=2\) is excluded because the
linked factorial moment in this chart has positive coefficients.
The complex roots of the quadratic have
\[
 \sin^2\frac{\arg z}{2}
   =\frac12-\frac14\sqrt{\frac{4n+2}{9n+6}}
   \le \frac13+\frac1{4(9n+6)} .
\]
Here the excess over \(1/3\) follows by rationalizing the
difference of the square roots and bounding the resulting
denominator from below by \(2/3\).
For \(n\ge2\) the
[exact linked pencil sector bound](LINKED_BETA_VISHNYAKOVA_PENCIL.md#3-order-uniform-angle-and-compactness)
requires at every native order-\(n\) root
\[
                \sin^2\frac{\arg z}{2}>
                      \frac{n+1}{3n-1}.
\]
This is larger than \(1/3+1/(4(9n+6))\). Thus (2) cannot
vanish at an order-\(n\) root in this chart. The case \(n=1\)
is immediate from the real first root and even positivity.
Interchanging \(X,Y\) gives the same all-order result for
\((a,b)=(2,1)\). Equal exponent rows and characteristic-zero
scalar extension are again allowed. General \(a,b>0\) still
require excluding the roots of (2) on the actual linked orbit.
