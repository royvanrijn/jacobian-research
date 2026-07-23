# Higher-gcd strata of the degree-42 binary-cubic pencil

This note carries out the subresultant test suggested by the
singularity-theoretic interpretation of the degree-forty-two normal cone.
It concerns the branch

\[
 w_0=w_1=0,\qquad w_2\ne0,
\]

where the first terminal equations are two binary cubics \(p_3,q_3\).
Their resultant is a nonzero scalar times \(AB\), and their degree-one
subresultant is a nonzero scalar times

\[
 S_1(u,v)=\alpha u+\beta v.
\]

Thus

\[
 \mathcal H=(AB,\alpha,\beta)
\tag{1}
\]

is the locus where the common factor has degree at least two, subject to
the usual homogeneous leading-coefficient charts.

The exact certificate is
[`verify_degree42_higher_gcd_strata.py`](../scripts/verify_degree42_higher_gcd_strata.py).

## 1. Rational decomposition

Put

\[
\begin{aligned}
P_1={}&(e_2^4-6e_1e_2^2+12e_1^2,\;
        e_2^3-4e_1e_2+6t),\\
P_2={}&(e_1,e_2),\\
P_3={}&(t,4e_1-e_2^2),\\
P_4={}&(e_1,t).
\end{aligned}
\tag{2}
\]

Exact radical and intersection calculations over \(\mathbb Q\) give

\[
\begin{aligned}
\sqrt{(A,\alpha,\beta)}&=P_2\cap P_3\cap P_4,\\
\sqrt{(B,\alpha,\beta)}&=P_1\cap P_3\cap P_4,\\
\sqrt{\mathcal H}&=P_1\cap P_2\cap P_3\cap P_4.
\end{aligned}
\tag{3}
\]

Hence the higher-gcd locus has four rational weighted curves.  The two
shared curves \(P_3,P_4\) are the intersections of the \(A\)- and
\(B\)-phenomena; \(P_2\) belongs only to \(A\), and \(P_1\) only to \(B\).

## 2. Geometric splitting

Let \(\rho^2=-3\).  The rational curve \(P_1\) splits over
\(\mathbb Q(\rho)\) into

\[
\begin{aligned}
e_2&=r,\\
e_1&=\frac{3\pm\rho}{12}r^2,\\
t&=\pm\frac{\rho}{18}r^3.
\end{aligned}
\tag{4}
\]

These are the two Galois-conjugate higher-collision branches carried by
the rational factor \(B\).

The cubic gcd at a nonvertex point of each curve has the following type:

| curve | parametrization | generic cubic gcd |
|---|---|---|
| \(P_2\) | \(e_1=e_2=0,\ t=r^3\) | \(uv\) |
| \(P_4\) | \(e_1=t=0,\ e_2=r\) | \((u-rv)^3\) |
| \(P_3\) | \(t=0,\ e_2=2r,\ e_1=r^2\) | \((u-rv)^2\) |
| \(P_1^\pm\) | (4) | a squarefree quadratic |

Thus the subresultant locus consists of familiar double- and triple-root
collisions.  It does not introduce another generic divisor.

## 3. Quartic closure on the punctured curves

After eliminating the three pivot variables through cubic order, let
\(p_4,q_4\) be the terminal quartic Kuranishi forms.  At a nonzero point
of each weighted curve, the checker proves

\[
 \gcd\bigl(\gcd(p_3,q_3),p_4,q_4\bigr)=1.
\tag{5}
\]

The corresponding homogeneous envelopes have:

| curve | Hilbert vector | length | nilpotence |
|---|---|---:|---|
| \(P_2\setminus\{0\}\) | \((1,2,3,2)\) | 8 | \(\mathfrak m^4=0\) |
| \(P_4\setminus\{0\}\) | \((1,2,3,3,1)\) | 10 | \(\mathfrak m^5=0\) |
| \(P_3\setminus\{0\}\) | \((1,2,3,2)\) | 8 | \(\mathfrak m^4=0\) |
| \(P_1\setminus\{0\}\) | \((1,2,3,2)\) | 8 | \(\mathfrak m^4=0\) |

Each geometric punctured branch is one weighted \(\mathbb G_m\)-orbit
after the indicated scalar extension.  The witness calculation therefore
describes its generic homogeneous algebra, not an isolated accidental
point.

The transported defect is already in the fifth normal power.  The
Kuranishi nilpotence cutoff consequently makes the synchronization exact
on every punctured higher-gcd curve.

## 4. The surviving vertex

All four curves meet at

\[
 e_1=e_2=t=0.
\tag{6}
\]

With \(w_2=1\), both cubic forms vanish there, while the quartics are

\[
 p_4=\frac5{16}uv^3,\qquad
 q_4=\frac5{64}v^4.
\tag{7}
\]

They retain the common cube \(v^3\), so the quartic envelope is not
Artinian.  This is the genuine contact-five/quintuple core.  It is distinct
from the sevenfold monomial collision, which additionally requires
\(w_2=0\).

Therefore the subresultant calculation has a sharp outcome:

\[
\boxed{\text{all punctured components of }\mathcal H\text{ close at
quartic order; only their common vertex survives.}}
\tag{8}
\]

## 5. Remaining limitation and next test

Equation (8) analyzes the higher-gcd locus \(S_1=0\).  It does not yet prove
that the residual quartic scalars \(\rho_{L_A}\) and \(\rho_{L_B}\) vanish
only on \(\mathcal H\).  The next exact calculation should therefore:

1. compute cleared representatives of \(\rho_{L_A}\) modulo \(A\) and
   \(\rho_{L_B}\) modulo \(B\);
2. factor their zero divisors in the corresponding coordinate rings;
3. compare those factors with \(P_1,\ldots,P_4\);
4. isolate any additional quartic-zero curve before computing a fifth
   Kuranishi term.

If those zero divisors are supported on \(\mathcal H\), the entire
\(A/B\) exceptional locus with \(w_2\ne0\) reduces to the single vertex
(6).
