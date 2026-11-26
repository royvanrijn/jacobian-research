# Exclusion of a third-order tail at a quadruple leading root

Work over a field of characteristic zero on \(2xy+t^2=1\).
Let \(\Lambda\) take the angular constant term and integrate in
\(t\) over \([-1,1]\) with measure \(dt/2\). Write \(U=1+x\)
and \(v=Ut\).

## Proposition (Third-order tail exclusion)

A regular winding-one cubic of least sphere degree seven of the form
\[
 f=\frac Ux\left(P_0(v)+\sum_{j=3}^7x^j(a_jv+b_j)\right),
 \quad P_0(v)=(v^2-1)(v+d),\qquad a_6=a_7=0,
\]
has infinitely many nonzero positive moments.

Suppose that its positive moments vanish eventually. The
[quadruple-root boundary corollary](WINDING_ONE_CUBIC_LOGARITHMIC_BOUNDARY.md#consequences-in-degree-seven)
gives vanishing at every positive order. We prove that the moments
through order ten, together with two inverse-period equations, force
all eight free tail coefficients to vanish. That would leave the
constant cubic, whose sphere degree is at most six by the
[high-tail argument](WINDING_ONE_CUBIC_QUADRUPLE_ROOT_HIGH_TAIL.md).

### Moment and inverse-period equations

Set \(B_j(T)=a_jT+b_j\) and \(B(X,T)=\sum_{j=3}^7X^j B_j(T)\).
The [fixed interval formula](WINDING_ONE_CUBIC_QUADRUPLE_ROOT_SPARSE_COEFFICIENTS.md#the-fixed-interval-formula)
gives the following exact expression for \(M_m=\Lambda(f^m)\):
\[
 M_m=\sum_{\ell=1}^{\lfloor m/3\rfloor}\binom m\ell
       \sum_{k=3\ell}^m\binom{m-1}{k-1}
       \frac12\int_{-1}^1 P_0(T)^{m-\ell}
                        [X^k]B(X,T)^\ell\,dT .
 \tag{1}
\]
The term with no factors from the tail vanishes by angular degree.
To evaluate the interval integrals use
\[
 I_{n,e}=\frac12\int_{-1}^1 T^eP_0(T)^n\,dT
 =\sum_{\substack{0\le j\le n\\j+e\ {\rm even}}}\ \sum_{h=0}^n
 \frac{(-1)^{n-h}\binom nj\binom nh}{2h+j+e+1}d^{n-j}.
 \tag{2}
\]
We write \(I_n=I_{n,0}\) and \(J_n=I_{n,1}\).

The first two inverse deficits give the same equations as in the
[fourth-order tail calculation](WINDING_ONE_CUBIC_QUADRUPLE_ROOT_FOURTH_ORDER_TAIL.md#the-linear-rows).
With \(\beta_j=\binom{-1/3}{j-1}\), they now read
\[
 N_2=a_3-\frac89 a_4+\frac{22}{27}a_5=0,\qquad
 N_3=\sum_{j=3}^7\beta_j(b_j-da_j/3)=0.                 \tag{3}
\]
Indeed the inverse-period condition takes the constant terms at
\(x=\infty\) of \(-p/(3Us)\) and \(-q/(3Us^2)\), where
\(s=(x/U)^{1/3}\) and \(p,q\) are the two depressed coefficients.
The resulting coefficient functionals are
\([x^1](1+x^{-1})^{-2/3}p\) and
\([x^1](1+x^{-1})^{-1/3}q\). At \(j\ge3\), the respective
coefficients of \(x^j\) in \(p,q\) are \(a_j,b_j-da_j/3\).
This gives (3), after dividing the first row by its coefficient
of \(a_3\).

### Certificates over the rational function field

On the locus \(I_2 I_3 I_4 I_5\ne0\), solve (3) and the first
four nontrivial moment equations using \(u=a_4\) and \(w=a_5\)
as the remaining variables. Specifically,
\[
 \begin{aligned}
 a_3&=\frac89u-\frac{22}{27}w,\\
 b_3&=-J_2a_3/I_2,\\
 b_4&=-3b_3-J_3(3a_3+u)/I_3,\\
 b_5&=-6b_3-4b_4-J_4(6a_3+4u+w)/I_4 .
 \end{aligned}
\]
The second inverse row solves for \(b_7\) in terms of \(b_6,u,w\)
because \(\beta_7=728/6561\). In \(M_6\) the coefficient of
\(b_6\) is \(6I_5\): it occurs only in the linear term of (1),
and \(b_7\) has yet to enter. This solves \(b_6\), and thereby
\(b_7\), in terms of \(u,w\). These substitutions require only
the indicated four interval moments and nonzero rational scalars
to be invertible.

Use \(M'_m\) for the moment polynomials after substitution.
The [exact checker](verify_winding_one_cubic_third_order_tail.py)
constructs and verifies three matrix identities in
\(\mathbb Q(d)[u,w]\), expressing the targets \((u,w)\) through
the respective generator rows
\[
 (M'_7,M'_8,M'_9),\qquad (M'_7,M'_8,M'_{10}),\qquad
 (M'_7,M'_9,M'_{10}).
\]
For each identity, let its denominator polynomial be the
squarefree part of the least common multiple of the denominators
in \(d\) of the entries of its coefficient matrix. The checker
obtains degrees \(209,219,291\) and computes their monic common
gcd as
\[
                    \frac{d(d-1)(d+1)(7d^2+1)}7 .       \tag{4}
\]
The checker also verifies that all denominator factors on the
substituted phase and moment side are among the guard factors
given below. If \(d\) satisfies the guard and is outside the root
set of (4), at least one complete matrix identity specializes.
Vanishing of its three moment generators then gives \(u=w=0\).
The solved coefficients all vanish as well.

### Exceptional certificates

Up to nonzero rational factors, the four guard polynomials from
(2) are
\[
 I_2\sim7d^2+1,\quad I_3\sim d(3d^2+1),\quad
 I_4\sim143d^4+78d^2+3,\quad I_5\sim d(13d^4+10d^2+1).
\]
Together with (4), this leaves the roots of
\[
 d,\quad d-1,\quad d+1,\quad 7d^2+1,\quad 3d^2+1,\quad
                     143d^4+78d^2+3,\quad 13d^4+10d^2+1.
\]
The checker verifies irreducibility of the nonlinear factors over
\(\mathbb Q\). In the field defined by each listed factor,
it constructs the moments from (1) in all eight original tail
variables, using no guard division. It verifies exact lift
identities expressing the square of each of the eight variables
through \(N_2,N_3,M_3,\ldots,M_9\). Every putative pure phase at
these values of \(d\) therefore has zero tail, completing the
exclusion.

### Replay

Run the checker from the repository root with Python 3 and Singular 4
on the path.

```sh
python3 research/foundations/extended-geometry/verify_winding_one_cubic_third_order_tail.py
```

The replay uses Python rational numbers for (2) and exact
arithmetic in Singular. It verifies each lift by multiplying the
matrix and checking its target before using the identities for
specialization. The degree-seven boundary corollary supplies the
moment vanishing under eventual purity.

## Corollary (Third-order depressed profiles)

In the normalized quartic-leading family with least sphere degree
seven and \(\deg d\le2\), eventual purity is impossible if both
depressed coefficients agree with their constant terms through
order \(x^2\).

Use the [fixed-endpoint translation](WINDING_ONE_CUBIC_CONSTANT_DEPRESSED_COEFFICIENTS.md#reduction-to-one-polynomial-of-the-base-variable)
to make the quadratic coefficient constant \(d(0)\), preserving
all moments. The translated coefficients have degree bounds five
and seven, since \(\deg d\le2\). Their differences from the
regularity constants are divisible by \(x^3\) under the stated
hypothesis. The least sphere degree stays seven. If \(b\) had
degree five, its leading coefficient is unchanged, since \(d^2\)
has degree at most four. Otherwise \(c\) had degree seven and
its leading coefficient is unchanged, since \(bd\) and \(d^3\)
have degree at most six. The translated phase therefore meets
the third-order tail proposition.
