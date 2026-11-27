# Core control at a moving Gaussian triple junction

For aligned planar balanced global maximizers near one, the
[bounded-region graph theorem](THREE_CELL_NEAR_ONE.md#tangent-alignment-across-bounded-physical-regions)
gives a chart about the actual vertex $p$. In that chart, the
noise-scale squared deformation and the displacement of $p$
obey the relative estimate
\[
 Z_e(h)+e(|p|^2+\|h\|_2^2)\le C_R\bigl(
  D_+(P_{\rm st}h)+D_-(P_{\rm sg}h)+e\|f\|_2^2+
        e^{-c_R/e^2}d^2\bigr).
 \tag{MJ0}
\]
Here $h$ is the displacement from the translated tripod,
$f$ includes the translation, $D_\pm$ are local pieces of the
interface squared-difference forms, and $d$ is the Gaussian
area of the defect from the translated tripod. The definitions
follow below. The estimate is uniform when $|p|/e$ diverges.

Put $e=\sqrt{2(1-\rho)}$, $\sigma=\sqrt{1-\rho^2}$, and $\alpha=1-\rho$. All constants are uniform as $\rho\uparrow1$ with $R>0$ fixed. Let $T_p=p+T$ be a translated central simplex with vertex $p$, allowing a common rotation. Choose a cyclic orientation of its rays. Write $\tau_k,n_k$ for the unit tangents and normals and $d_k$ for the signed label-transfer vectors. Let $P_{\mathrm{st}}$ project onto the orthogonal complement of $(1,1,1)$, and let $P_{\mathrm{sg}}=I-P_{\mathrm{st}}$. The map $Bv=\sum_k d_kv_k$ satisfies $B^*B=3P_{\mathrm{st}}$.

Assume $|p|$ is small and that the interfaces of $A$ near $p$ are normal graphs $t=h_k(r)$ over the rays of $T_p$, with $h_k(0)=0$ and $\operatorname{Lip}(h_k)\le\delta$ on $0\le r\le R$. More explicitly, the indicator difference from $T_p$ is the sum of the signed strips from $t=0$ to $t=h_k(r)$ over this range, plus an error supported in $\{|x-p|\ge R\}$. The signed-strip representation is exact on the strips, so the error is the indicator difference restricted to the omitted region. A normal-graph chart on a disk slightly larger than $B_R(p)$ gives this representation with each signed strip an actual transfer. The strips are disjoint for small $\delta$. The labels outside the chart can be arbitrary.

Suppose the centered score prices are $\beta$ and the adjusted scores tie at $p$. Write $d=\gamma_2\{A\ne T_p\}$, $z_\sigma=\int_0^R\phi_\sigma(r)\,dr$, and
\[
 w_\sigma(r)=z_\sigma^{-1}\phi_\sigma(r),\qquad
 a_h=P_{\mathrm{st}}\int_0^R w_\sigma(r)h(r)\,dr,\qquad
 Z_e(h)=e^{-1}\int_0^R e^{-c r^2/e^2}|h(r)|^2\,dr .
\]
Here $c>0$ is a fixed sufficiently small absolute constant. The comparison has no bound on $p/e$.

## Score estimate

For small $e$,
\[
 |a_h|^2\le C_R\left(e^2|\beta|^2+e^4|p|^2
       +(\delta^4+e^2|p|^2)Z_e(h)+e^{-c_R/e^2}d^2\right).
 \tag{MJ1}
\]

The score vector of $A$ at $p$ is $\frac13\mathbf 1+\beta$, since the scores sum to one. The score vector for $T_p$ at $p$ differs from $\frac13\mathbf 1$ by at most $C e|p|$. The conditional mean in coordinates centered at $p$ is $-\alpha p$, and convolution of the unit tripod with the standard Gaussian is Lipschitz. The vector score difference coming from the strips is exactly
\[
 \sum_k d_k\int_0^R \phi_\sigma(r+\alpha p\cdot\tau_k)
                      \int_0^{h_k(r)}\phi_\sigma(t+\alpha p\cdot n_k)\,dt\,dr .
\]
Its main term is $z_\sigma\phi_\sigma(0)B a_h$, with $z_\sigma$ bounded away from zero and $\phi_\sigma(0)\asymp e^{-1}$.

The exponent of $\phi_\sigma(r+\alpha p\cdot\tau_k)/\phi_\sigma(r)$ has modulus at most $C(|p|r+e^2|p|^2)$. Since $r\le R$, the ratio differs from one by at most $C_R(|p|r+e^2|p|^2)$. The ratio $\phi_\sigma(t+\alpha p\cdot n_k)/\phi_\sigma(0)$ differs from one by at most $C(t^2/e^2+e^2|p|^2)$. Use $|t|\le |h_k(r)|\le\delta r$. Cauchy–Schwarz bounds the strip error by $C_R e^{-1}(\delta^2+e|p|)Z_e(h)^{1/2}$, since the factors $r/e$ and $(r/e)^2$ are absorbed into a wider Gaussian.

For $|x-p|\ge R$, the conditional density at $p$, divided by the standard Gaussian density at $x$, is at most $C_R e^{-c_R/e^2}$. Indeed $|x-\rho p|\ge|x-p|/2$ for small $e$, and the negative conditional exponent absorbs $\exp(|x|^2/2)$ and the power of $e^{-1}$. The score error from the omitted region is at most $C_R e^{-c_R/e^2}d$. Divide the score identity by $z_\sigma\phi_\sigma(0)$, use $\|B a_h\|=\sqrt3|a_h|$, and square to obtain (MJ1).

## Local quadratic control of the core

For $u=P_{\mathrm{st}}h$, $v=P_{\mathrm{sg}}h$, and $H_e=\{(r,s)\in[0,R]^2:|r-s|\le2e\}$, define
\[
 D_+(u)=e^{-2}\iint_{H_e}|u(r)-u(s)|^2\,dr\,ds,\quad
 D_-(v)=e^{-2}\iint_{H_e}\min(rs/e^2,1)\,rs\,|v(r)/r-v(s)/s|^2\,dr\,ds .
\]
Values on the axes do not affect the integral. These are bounded
by $C_R$ times the corresponding simplex squared-difference forms
for any extension past $R$. The full same-ray Mehler interaction
kernel with its facet weights written against $dr\,ds$ is bounded
below by $c_R e^{-2}$ on $H_e$. In the sign sector, the difference
of the same-ray and cross-ray kernels equals the same-ray kernel
times $1-\exp(-3\rho rs/(2\sigma^2))$, comparable to
$\min(rs/e^2,1)$. [Simplex sector formulas (S18) and (S31)](SIMPLEX_INTERFACE_SYMMETRY.md)

The estimates are
\[
 Z_e(u)\le C_R(D_+(u)+|a_h|^2),\qquad
 Z_e(v)\le C_R(D_-(v)+e^2\|v\|_{L^2(0,R)}^2).
 \tag{MJ2}
\]
Partition $[0,R]$ into intervals of a common length $\ell\in[e/2,e]$. Each interval square and consecutive interval rectangle lies in $H_e$. For the first estimate write $\mu_j$ for the averages of $u$. Jensen's inequality and the variance identity show that $D_+$ controls the sum of squared consecutive mean differences and interval variances, each variance divided by $\ell$. The squared difference between $\mu_j$ and $\mu_0$ is at most $C(j+1)D_+$. Sum with Gaussian weights to control $Z_e(u-\mu_0)$, which in turn bounds $|a_h-\mu_0|^2$ by Cauchy–Schwarz.

For the second estimate set $\kappa(r)=v(r)/r$. From interval $j=1$ onward, the consecutive mean differences and interval variances of $\kappa$ enter $D_-$ with weights at least $c e^2 j^2$. The sum of their inverse weights is $O(e^{-2})$. In the rectangle between intervals zero and one the coefficient $\min(rs/e^2,1)rs$ is bounded below by $c r^2$ for $r$ in interval zero. The rectangle controls $e^{-1}\int_{I_0}|v(r)-r\mu_1|^2\,dr$, where $\mu_1$ is the mean of $\kappa$ on interval one. Telescoping the subsequent means from $\mu_1$ and summing with Gaussian weights bounds $Z_e(v-r\mu_1)$ by $C_R D_-$. Also $|\mu_1|^2\le C_R(e^{-2}D_-+\|v\|_{L^2(0,R)}^2)$ by comparing the interval means to their average on a fixed compact subinterval of $(0,R)$. Finally $Z_e(r\mu_1)\le C e^2|\mu_1|^2$. This proves (MJ2).

Combine (MJ1) and (MJ2) and absorb the term $C_R(\delta^4+e^2|p|^2)Z_e(h)$ for small $\delta,e$:
\[
 Z_e(h)+|a_h|^2
 \le C_R\left(D_+(u)+D_-(v)+e^2\|v\|_2^2+
            e^2|\beta|^2+e^4|p|^2+e^{-c_R/e^2}d^2\right).
 \tag{MJ3}
\]

## Displacement in the moved-reference chart

With the cyclic normal signs, $q_k=p\cdot n_k$ satisfies $P_{\mathrm{st}}q=q$ and $|q|^2=\frac32|p|^2$. Put $f=h+q$, the total linear normal speed from the central rather than translated simplex. The standard core mean of $f$ is $a_h+q$, its standard edge form is unchanged from $u$, and $P_{\mathrm{sg}}f=v$. The unweighted chain to a fixed physical subinterval of $(0,R)$ has $O_R(e^{-1})$ links with conductances bounded below. Therefore
\[
 e|a_h+q|^2\le C_R\bigl(D_+(u)+e\|P_{\mathrm{st}}f\|_2^2\bigr).
\]
Combine with (MJ3) and absorb the $e^5|p|^2$ term to obtain
\[
 e\bigl(|p|^2+\|h\|_2^2\bigr)
 \le C_R\bigl(
      D_+(u)+D_-(v)+e\|f\|_2^2+
      e^3|\beta|^2+ e^{-c_R/e^2}d^2
      \bigr).
 \tag{MJ4}
\]

Choose the two fixed volume bumps from (P4) inside $B_{R/8}(0)$,
on distinct open simplex rays. For an actual maximizer the
volume derivative matrix on these bumps is uniformly invertible
on the centered label space by convergence to the simplex.
To sharpen the price bound, subtract the first variations of
$N_\rho(A)$ and $N_\rho(T)$ against either bump as in (P15).
The pair integrand requires a changed label at some endpoint
and a point in the support of the bump. If the changed
endpoint is outside $B_R(p)$, its other endpoint must be in
$B_{R/8}(0)$. The conditional Gaussian bounds this integral,
including the derivative weight, by $C_R e^{-c_R/e^2}$ times
the Gaussian area of the far defect. On the remaining
changed points the conditional integral of the derivative weight
is uniformly bounded. The area of the defect from $T$ inside
$B_R(p)$ is at most $C_R(|p|+\|h\|_2)$ by normal strip
integration and the Gaussian area of the swept rays between
$T$ and $T_p$. Globally the area of the defect from $T$ is
at most $d+C|p|$. Using that the simplex first variation
vanishes therefore gives

\[
 |\beta|\le C_R(|p|+\|h\|_2)+C_R e^{-c_R/e^2}d .
\]

Insert this in (MJ4). The term $e^3(|p|^2+\|h\|_2^2)$ is
absorbed for small $e$. Then substitute the resulting bound
in (MJ3), using $v=P_{\rm sg}f$ and $e^2\le e$, to obtain
(MJ0). All local norms here are on $[0,R]$.

This estimate accounts for the translation-sized part of the
core error directly from the actual adjusted score tie. An
exact comparison still requires a remainder estimate for a
single competitor whose changes extend from the moving
junction through the exterior tubes and possibly their
complement.
