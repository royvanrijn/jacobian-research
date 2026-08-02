# Maximal shifted-power failure and the next minimum frontier in GVC(3)

## 1. Setup

Work over a characteristic-zero field. Put

\[
\rho=t^2+xy,\qquad A=\rho+x^2,
\qquad C=y\rho^2-2xt^2\rho-x^3t^2,
\]

\[
P_6=AC^2,\qquad \Delta=4\partial_x\partial_y+\partial_t^2.
\tag{1.1}
\]

For every integer \(k\ge6\), define

\[
P_k=\rho^{k-6}P_6,\qquad \Lambda_k=\Delta^k.
\tag{1.2}
\]

Then \(P_k\) is homogeneous of degree \(2k\), and the companion notes prove

\[
\Lambda_k^n(P_k^n)=\Delta^{kn}(P_k^n)=0
\qquad(n\ge1).
\tag{1.3}
\]

They also prove, for \(1\le \ell\le n\),

\[
\begin{aligned}
D_{k,n,\ell}
:=\Delta^{kn+\ell}(x^{2\ell}P_k^n)
={}&2^{(k+2)n+\ell}(kn+\ell)!(2n)!\\
&\times\frac{(2kn+2\ell+1)!!}{(4n+1)!!}
\binom{n-1}{\ell-1}\ne0.
\end{aligned}
\tag{1.4}
\]

## 2. Fischer transfer

Let

\[
q(X,Y,T)=4XY+T^2.
\]

For homogeneous symbols \(S\) and polynomials \(F\), Fischer adjointness gives

\[
S(\partial)(x^rF)\big|_0
=(\partial_X^rS)(\partial)F\big|_0.
\tag{2.1}
\]

Since \(q\) is linear in \(X\),

\[
\partial_X^{2d}q^N
=4^{2d}\frac{N!}{(N-2d)!}Y^{2d}q^{N-2d}.
\tag{2.2}
\]

This elementary identity converts the positive-phase multiplier ladder into
statements about the shifted powers \(P_k^{m+d}\).

## 3. Every shifted power tail fails

Take \(n=m+d\) and \(\ell=d\) in (1.4). Applying (2.1)--(2.2), and factoring
\(q^{km}\), gives the exact identity

\[
\boxed{
\begin{aligned}
&\Delta^{(k-1)d}\partial_y^{2d}
 \left(\Lambda_k^m(P_k^{m+d})\right)\\
&\quad=
2^{(k+2)m+(k-1)d}
\bigl(km+(k-1)d\bigr)!(2m+2d)!\\
&\qquad\times
\frac{\bigl(2km+2(k+1)d+1\bigr)!!}
     {(4m+4d+1)!!}
\binom{m+d-1}{d-1}.
\end{aligned}}
\tag{3.1}
\]

Every factor on the right is nonzero in characteristic zero. Hence:

> **Theorem 3.1 — maximal power-tail failure.** For every \(k\ge6\), every
> \(d\ge1\), and every \(m\ge1\),
> \[
> \boxed{\Lambda_k^m(P_k^{m+d})\ne0.}
> \tag{3.2}
> \]

Thus the counterexample does not merely violate GVC for one specially chosen
multiplier. It violates every shifted-power conclusion in the power
formulation of GVC. In particular, already for \(d=1\),

\[
\Lambda_k^m(P_k^{m+1})\ne0
\qquad(m\ge1).
\tag{3.3}
\]

For \(k=6\), (3.1) becomes

\[
\begin{aligned}
\Delta^{5d}\partial_y^{2d}
 \left(\Delta^{6m}P_6^{m+d}\right)
={}&2^{8m+5d}(6m+5d)!(2m+2d)!\\
&\times\frac{(12m+14d+1)!!}{(4m+4d+1)!!}
\binom{m+d-1}{d-1}.
\end{aligned}
\tag{3.4}
\]

This is a direct explicit counterexample to the equivalent
"power-of-the-radical-polynomial" formulation of the generalized vanishing
conjecture.

## 4. Exact polyharmonic depth

Taking \(\ell=1\) in (1.4) and applying (2.2) only once gives

\[
\boxed{
\partial_y^2\Delta^{kn-1}(P_k^n)
=
2^{(k+2)n-3}(kn-1)!(2n)!
\frac{(2kn+3)!!}{(4n+1)!!}\ne0.
}
\tag{4.1}
\]

Together with (1.3), this proves:

> **Corollary 4.1 — exact trace length.** The polynomial \(P_k^n\) has exact
> polyharmonic index \(kn\):
> \[
> \Delta^{kn}(P_k^n)=0,
> \qquad
> \Delta^{kn-1}(P_k^n)\ne0.
> \tag{4.2}
> \]

So the pure cancellation occurs at the final possible trace, rather than
through an accidental earlier degree deficit.

## 5. A self-multiplier rank-one SIC counterexample

Introduce dual variables \(\xi=(\xi_x,\xi_y,\xi_t)\) and put

\[
f_k=q(\xi)^kP_k(x,y,t).
\tag{5.1}
\]

Under the contraction map \(\mathcal E_3\),

\[
\mathcal E_3(f_k^m)=\Lambda_k^m(P_k^m)=0.
\tag{5.2}
\]

For every fixed \(d\ge1\), choose the coordinate-only multiplier

\[
g_d=P_k^d.
\tag{5.3}
\]

Then Theorem 3.1 gives

\[
\boxed{
\mathcal E_3(g_df_k^m)
=\Lambda_k^m(P_k^{m+d})\ne0
\qquad(m\ge1).
}
\tag{5.4}
\]

Thus rank-one SIC fails in three pairs even when the multiplier belongs to
the one-generator coordinate subalgebra \(k[P_k]\). Every positive power of
that same coordinate polynomial detects the failure at every order.

## 6. Why the current endpoint architecture starts at \(\Delta^6\)

The construction can be written on the affine quadric \(\rho=1\) in the
one-profile form

\[
\Phi_{r,s,a}
=x^{-2r}(1+ax^2)^r
 \left(1-t^2(1+ax^2)^2\right)^s,
\qquad a\ne0.
\tag{6.1}
\]

In the quadric coordinate ring,

\[
1-t^2(1+ax^2)^2=xC_a,
\qquad
C_a=y\rho^2-2axt^2\rho-a^2x^3t^2.
\tag{6.2}
\]

Modulo \(x\), \(C_a\equiv yt^4\), so the endpoint factor has exact
\(x\)-adic order one. Consequently (6.1) is polynomial exactly when

\[
s\ge2r.
\tag{6.3}
\]

Its homogeneous lift is

\[
P_{r,s,a}=x^{s-2r}(\rho+ax^2)^rC_a^s,
\tag{6.4}
\]

and has degree

\[
\deg P_{r,s,a}=6s.
\tag{6.5}
\]

The matching differential operator is therefore \(\Delta^{3s}\). Since
\(r\ge1\) and \(s\ge2r\),

\[
3s\ge6.
\tag{6.6}
\]

> **Theorem 6.1 — scoped minimum.** Within the complete phase-square,
> one-profile endpoint-contact architecture (6.1), \(\Delta^6\) is the
> minimum possible Laplacian power. The minimum is attained uniquely at the
> contact/winding pair \((r,s)=(1,2)\), up to the profile torus and scaling.

Therefore a counterexample for \(\Delta^k\), \(k\le5\), must use a genuinely
different mechanism: several phase profiles, a different homogeneous power
map, or cancellation not generated by one endpoint-contact factor.

## 7. The complete linear parity completion of Long's quartic dies

There is another natural attempt at \(k=2\). On \(\rho=1\), Long's sphere
polynomial in the present normalization is

\[
L=(1+x)(y-(2+x)t^2)=E+O,
\tag{7.1}
\]

where

\[
E=xy-2t^2-x^2t^2,
\qquad
O=y-3xt^2
\tag{7.2}
\]

are respectively antipodally even and odd. For an arbitrary linear form
\(H=ax+by+ct\), the complete linear parity repair is the homogeneous quartic

\[
F_H=\rho(xy-2t^2+Hy)-x^2t^2-3Hxt^2.
\tag{7.3}
\]

Its restriction to the sphere is \(E+HO\). Exact normalized spherical
moments give first

\[
5a-3b=0.
\tag{7.4}
\]

Put \(b=5a/3\) and \(z=c^2\). Up to nonzero rational factors, moments two
and three give

\[
28a^2-52a+27z-63=0,
\tag{7.5}
\]

\[
7164a^3-36868a^2-81341a-24453=0.
\tag{7.6}
\]

After using (7.5), moment four gives

\[
77776a^4+137224a^3-745076a^2-1119246a-198747=0.
\tag{7.7}
\]

The resultant of (7.6) and (7.7) is

\[
-6466167050191094761727778002592000\ne0.
\tag{7.8}
\]

Hence no member of (7.3) has even its first four pure moments zero. This
closes the entire linear parity-homogenization of Long's quartic, not all
homogeneous quartics.

## 8. Updated minimum frontier

The new result separates three notions cleanly.

1. The dimension problem for homogeneous GVC is finished: the first failure
   is dimension three.
2. The power-tail form is maximally false: every shift \(d\ge1\) fails for
   every \(m\).
3. The remaining internal complexity problem is now
   \[
   \boxed{\text{find a homogeneous GMC/GVC witness of degree }2k<12,}
   \]
   equivalently a counterexample for \(\Delta^k\) with \(1\le k\le5\), or
   prove a lower bound in a declared architecture.

The one-profile endpoint mechanism cannot lower \(k\), and the complete
linear parity repair does not produce \(k=2\). The next calculation should
therefore be a multi-profile homogeneous sphere search in degrees six,
eight, and ten, with the antipodal parity quotient imposed before moment
elimination.

## 9. Reproduction

Run

```bash
python3 scripts/verify_gvc3_power_tail_and_minimum.py
```

The checker verifies the exact shifted-power detector for several
\((k,m,d)\), exact polyharmonic depth, the coordinate-self-multiplier SIC
identity, and the characteristic-zero resultant (7.8). The all-order claims
are the Fischer-transfer calculations above, not extrapolations from the
bounded replay.
