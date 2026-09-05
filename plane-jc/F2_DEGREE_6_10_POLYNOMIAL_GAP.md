# A sharp polynomial gap bound and the F2 normal r=9 exclusion

**Status: proved, exact algebra plus written Puiseux transport (PF2D6O1).**
For a birational polynomial parametrization of degrees `(6,10)`, the first
odd Puiseux gap at infinity is at most 21. An explicit parametrization
attains 21. Consequently the normal `r=9` row of the
[degree-six Stein reduction](F2_GEOMETRIC_DEGREE_SIX_STEIN_REDUCTION.md)
is impossible. Its remaining normal odd rows are `r=5,7`; the nonnormal
conductor families remain open. This does not exclude geometric degree six,
the coordinate-degree pair `(75,125)`, or JC2.

The missing condition in the earlier numerical ledger was polynomial
parametrizability. A valid Abhyankar delta sequence need not be realized by
polynomials of the prescribed degrees. This distinction already appears
in Sathaye–Stenerson's
[*On Plane Polynomial Curves*, Section 4](https://www.ms.uky.edu/~sohum/ma561/affgeom/papnew.pdf).
The argument below gives its own exact certificate for degrees six and ten;
it does not infer an exclusion from that paper's different example.

## 1. The polynomial gap theorem

Work over an algebraically closed field of characteristic zero. Let
`p,q in k[t]` have degrees six and ten and satisfy `k(p,q)=k(t)`.
Scale their leading coefficients to one and expand at infinity using
\[
 T=p(t)^{1/6},\qquad T/t\longrightarrow1.
\]
Write `ell` for the largest odd exponent with nonzero coefficient in
the Laurent expansion `q(T)`. Then this exponent exists and
\[
 \boxed{\ell\ge-11,\qquad L:=10-\ell\le21.}
\]

Affine translation of `t` removes the coefficient of `t^5` in `p`;
translation of `p` removes its constant. These changes preserve the
first odd gap: translation of `p` changes `T` by an odd Laurent series
in `T` with only even relative exponents. Thus write
\[
 p=t^6+a t^4+b t^3+c t^2+d t.
 \tag{1}
\]

Suppose for contradiction that the odd terms of `q(T)` vanish through
exponent `-11`. Its polynomial part is then necessarily
\[
 q=\left[p^{5/3}+B p^{4/3}+C p+D p^{2/3}+E p^{1/3}+F\right]_+,
 \tag{2}
\]
where brackets mean the polynomial part in `t`. Indeed the polynomial
part of `T^j` is monic of degree `j`; removing the positive even
terms successively leaves a polynomial of negative degree, hence zero.
The free summands `Cp+F` do not affect any odd term, so discard them.

Put `V=p^(5/3)+B*p^(4/3)+D*p^(2/3)+E*p^(1/3)`. Starting with `V`,
cancel its coefficients of `t^-2,t^-4,...,t^-10` successively by adding
constant multiples of `p^(-1/3),...,p^(-5/3)`. Denote the intervening
odd coefficients at `t^-j` by
\[
 \rho_j,\qquad j=1,3,5,7,9,11.
 \tag{3}
\]
All six must vanish. To see this without confusing the two parameters:
the added terms are even powers of `T`. Once earlier odd coefficients
vanish, the next odd leading residual has the same order in `t` and
`T`. Its coefficient in `q(T)` is `-rho_j`. This triangular
argument also identifies the first nonzero obstruction.

These are rational polynomials in `a,b,c,d,B,D,E`. They are generated
without series inversion. For `u=a*z^2+b*z^3+c*z^4+d*z^5` and
`(1+u)^alpha=sum h_n*z^n`, use
\[
 h_0=1,\qquad
 n h_n=\sum_{\substack{i\in\{2,3,4,5\}\\i\le n}}
       ((\alpha+1)i-n)u_i h_{n-i}.
 \tag{4}
\]
Only coefficients through order 21 are needed.

Let `f_j` be `rho_j` multiplied by its positive common denominator,
as recorded in the certificate. Exact rational polynomials `A_j` satisfy
\[
 \boxed{b^{12}=
 A_1 f_1+A_3 f_3+A_5 f_5+A_7 f_7+A_9 f_9+A_{11} f_{11}.}
 \tag{5}
\]
The six multipliers have respectively `142,101,52,34,21,12` terms.
Their complete coefficients are retained in
[the exact certificate](../artifacts/generated-results/f2-degree-6-10-gap-v1.json).
Equation (5) is checked by literal multiplication and addition, independent
of a Gröbner basis membership decision. In particular (3) forces `b=0`.

The remaining parameter branch has the much smaller identity
\[
 \left.\left(\rho_5+\frac{2a}{3}\rho_3+
                         \frac c3\rho_1\right)\right|_{b=0}
       =-\frac5{81}d^3.
 \tag{6}
\]
It forces `d=0` without division by `d`. Equation (1) is now even,
and every polynomial part in (2) is even as well. Hence
`k(p,q) subset k(t^2)`, contradicting birationality. This proves the
theorem, including the possibility that no odd coefficient ever occurs.

## 2. Sharpness and the surviving r=7 target

The pair
\[
 p=t^6+t^3,\qquad
 q=t^{10}+\frac53t^7+\frac79t^4+\frac7{81}t
 \tag{7}
\]
has normalized parameters
`(a,b,c,d,B,D,E)=(0,1,0,0,0,2/9,0)`. Direct substitution gives
\[
 (\rho_1,\rho_3,\rho_5,\rho_7,\rho_9,\rho_{11})
       =(0,0,0,0,0,1/19683).
 \tag{8}
\]
The first odd term of `q(T)` is therefore
`-T^-11/19683`, so `L=21`. Exact polynomial gcd computation gives
\[
 \gcd_{k(t)[s]}(p(s)-p(t),q(s)-q(t))=s-t,
 \tag{9}
\]
which verifies that this is a normalization parametrization.
Its finite derivative gcd is one, so this particular target is immersed
on the affine line. No connected six-sheet cover or Keller filling is
asserted for it.

## 3. Application to the terminal cubic

The normal odd rows in the Stein reduction have projective infinity
characteristic `(4;10,10+3r)` and normalization degrees `(6,10)`.
The first odd gap is consequently `L=3r`.

For completeness, this gap agrees with the one in Section 1.
In `tau=1/T` write
`p=tau^-6, q=tau^-10*H(tau)`, where the first odd relative term of
`H` has order `L`. Projective coordinates are
`x=p/q=tau^4/H` and `y=1/q=tau^10/H`. Normalizing `x=s^4`
gives `y=s^10*H^(3/2)` with its argument reexpanded in `s`.
The even part of this reparametrization preserves parity; the first odd
coefficient is multiplied by `3/2`, a nonzero scalar. Thus the first
odd relative order is still `L`.

The theorem gives `3r<=21`. It excludes `r=9`, with former
delta sequence `(10,6,3)`, while (7) confirms that polynomiality alone
cannot exclude `r=7`. The live normal odd ledger is now

| r | Local infinity semigroup | Delta sequence | Affine delta |
| --- | --- | --- | --- |
| 5 | `<4,10,35>` | `(10,6,15)` | 15 |
| 7 | `<4,10,41>` | `(10,6,9)` | 12 |

Their global cover and affine-plane filling conditions remain necessary.
The positive-conductor terminal families are unaffected.

There is now a finite exact coefficient interface for those two normal
target rows. In (1)–(2), their loci are respectively
`rho_1=rho_3=0, rho_5!=0` and
`rho_1=rho_3=rho_5=rho_7=rho_9=0, rho_11!=0`.
The latter forces `b!=0` by (6), so scaling the parameter permits `b=1`
over the algebraic closure. These conditions fix the odd infinity gap,
not the affine singularities or the covering monodromy. Restore `Cp+F`
and the initial translations/scalings when comparing to fixed F2 target
coordinates. The checker exposes the exact obstruction polynomials for
that next step; no new broad Laurent search is needed to recover them.

## Reproduction and assurance

```bash
.venv/bin/python scripts/verify_f2_degree_6_10_gap.py
# Optional: regenerate the exact multipliers with Singular, then replay.
.venv/bin/python scripts/verify_f2_degree_6_10_gap.py --regenerate-certificate
```

Default replay reconstructs (3) from (4), verifies (5) and (6) over the
rationals, and checks (7)–(9). It also rejects a changed multiplier.
No bounded coefficient search is used to infer the theorem. The passage
from an arbitrary parametrization to (2), and the Puiseux transport in
Section 3, are explicit written proof steps. Full formal verification,
independent end-to-end review, and novelty priority are not claimed.
