# Proper volume for winding-one cubic sphere phases

Let \(K\) have characteristic zero. On the sphere \(2xy+t^2=1\),
write a regular phase in the literal Laurent chart as
\[
 f=x^{-1}\bigl(C(x)+B(x)t+D(x)t^2+E(x)t^3\bigr),
 \qquad C,B,D,E\in K[x].
\]
The normalized sphere functional \(\Lambda\) takes the angular constant
term in \(x\) and integrates in \(t\) over \([-1,1]\) with measure
\(dt/2\).

**Theorem.** Suppose the angular exponent \(-1\) occurs, \(\deg E\ge 2\),
and every nonzero geometric root of \(E\) has multiplicity at most two.
Then \(\Lambda(f^m)\ne0\) for infinitely many positive integers \(m\).
In particular the horizontal degree-seven case with squarefree \(E\) of
degree five admits no eventually pure phase.

The proof also reduces the general case with negative winding and
\(\deg E\ge2\) to local area estimates at the nonzero roots of \(E\) of
multiplicity at least three. These estimates concern the polynomial
\[
 G(u,x)=\frac{E(x)+uD(x)+u^2B(x)+u^3C(x)}{x},\qquad u=t^{-1}.
\]
For each such root \(r\), it suffices to prove, for every fixed \(M>0\),
\[
 \int_{|u|<\epsilon}|u|^{-4}
    \operatorname{Area}\{x: |x-r|<\epsilon,\ |G(u,x)|\le M|u|^3\}
       \,d^2u<\infty.                                             \tag{1}
\]
One may choose the neighborhood separately at each root and for each
\(M\). For instance an area bound \(O_M(|u|^{2+\eta})\) for some
\(\eta>0\) suffices.

## From finite volume to nonvanishing

The coefficients generate a finitely generated field over \(\mathbb Q\).
Embed it in \(\mathbb C\). The moments are rational polynomials in the
coefficients, and the geometric root multiplicities persist. On the
Laurent chart the sphere's regular volume form is a nonzero constant
multiple of \(dx/x\wedge dt\). We show that this form has finite squared
integral on the part of the chart where \(|f|\le M\), for every \(M>0\).
Here squared integral means the integral of the positive density defined
by the form and its complex conjugate.

The sphere embeds into a smooth projective compactification obtained
from \(\mathbb P^1_x\times\mathbb P^1_t\) by blowing up \((0,1)\) and
\((0,-1)\). On the affine part these are the two centers of the ideal
\((x,1-t^2)\), and its blowup chart with \(x\) as the ideal generator
is exactly the sphere with \(2y=(1-t^2)/x\). Any toric blowups used below
are supported at infinity and preserve this open subset. Resolve the
rational pencil \(f\) by further blowups away from the sphere, where
\(f\) is already regular. On the resulting smooth surface \(Y\) over
finite phase values, \(f:Y\to\mathbb A^1\) is proper.

On every relatively compact chart in \(Y\), the function \(f\) is
bounded. The squared integral on the dense Laurent chart is unchanged
by the birational modifications there. A pole of order \(k\ge1\) of a
meromorphic top form along a divisor would force divergence of the
integral of \(|z|^{-2k}\,d^2z\) transversely to that divisor. Thus finite
volume on all bounded phase regions makes the pulled-back volume form
regular on \(Y\).

The [proper-volume lemma](../../../results/primitive-balanced-image-spectrum/proper-volume.tex)
then excludes a nonzero polynomial exponential period on a compact
closed cycle in \(Y\). The original real sphere
\(\Gamma=\{y=\overline x,\ t\in\mathbb R,\ 2|x|^2+t^2=1\}\)
lies in the retained affine sphere, including its two points with
\(x=0\). Its normalized volume period is \(\Lambda\), as in the
[sphere residue formula](GLOBAL_SPHERE_PURITY_RECOGNITION.md#1-a-closed-rational-period).
The period of \(e^{-\tau f}\) is entire and its coefficient of
\(\tau^m\) is \((-1)^m\Lambda(f^m)/m!\) after normalization. Its constant
term is one. Eventual purity would make it a nonzero polynomial.

## Boundary away from the two corners at infinity

We use a local sublevel estimate. If \(H(u,z)\) is holomorphic and
\(H(0,z)\) vanishes to order \(r>0\) at \(z=0\), Weierstrass preparation
writes \(H\) as a unit times a monic degree-\(r\) polynomial in \(z\) on
a smaller bidisc. Thus \(|H(u,z)|\le\delta\) implies that \(z\) is within
\(O(\delta^{1/r})\) of some root of the monic polynomial, and that
sublevel region has area \(O(\delta^{2/r})\), uniformly in \(u\).

At \(x=0\) or \(x=\infty\) with \(t\) finite, write \(f=u^{-N}H(u,t)\)
where \(u\) is a transverse coordinate and \(N>0\). At zero \(N=1\) by
the negative-winding hypothesis. At infinity \(N\ge \deg E-1\ge1\).
At a root of order \(r\) of the nonzero polynomial \(H(0,t)\), the
sublevel estimate gives area \(O_M(|u|^{2N/r})\) on \(|f|\le M\). The
volume form has pole order one, so the integral converges there. When
\(H(0,t)\ne0\), the bounded phase region stays away from the divisor.

At \(t=\infty\) and finite nonzero \(x\), the volume form has pole order
two and \(f=u^{-3}G(u,x)\). At a root \(r\) of \(E\) of multiplicity
\(m_r\), the area in (1) is \(O_M(|u|^{6/m_r})\). This proves (1) for
\(m_r=1,2\). At a point where \(E(x)\ne0\), bounded phase values again
stay away from the divisor.

## The toric corners

Near either corner the original coordinates are \((v,u)\), with
\(v=x\) or \(v=1/x\) and \(u=1/t\). Refine the first quadrant by the
rays on which the set of maximizers in the literal support of \(f\)
can change. Take a regular subdivision of that fan, so neighboring
primitive lattice rays form unimodular pairs. Such a subdivision in two
dimensions is obtained by the continued fraction subdivisions of the
rational cones. This gives smooth monomial charts over the corner.

On a ray \((a,b)\), the orders of \(v,u\) are \(a,b\), respectively.
The pole orders of \(f\) and the volume form along this divisor are
\[
 N=\max_{(i,j)\in\operatorname{supp}(f)}(\pm ai+bj),\qquad b+1,
\]
where the plus sign is for \(v=1/x\). For the latter order, write the
volume up to sign as \(u^{-1}d\log v\wedge d\log u\), and use the
unimodular monomial coordinates. At an interior ray the leading
coefficient is a Laurent polynomial in a coordinate on the open
divisor. The multiplicity of any nonzero root is at most the span of
its exponents. To see the span, the lattice points on the maximizing
line differ by integer multiples of a primitive perpendicular to
\((a,b)\); in a unimodular chart one step changes the divisor coordinate
exponent by one. When \(a>0\), one step changes the original height
\(j\) by \(a\) in absolute value. The span is therefore at most three.
If a root of order \(r\) occurs and \(N/r>b\), the sublevel estimate gives
local convergence by
\[
 \int |s|^{-2(b+1)+2N/r}\,d^2s<\infty
\]
in a transverse coordinate \(s\). A nonvanishing leading coefficient
again keeps the bounded phase region away from the divisor.

At \(x=\infty,t=\infty\), the top monomial
\(x^{\deg E-1}t^3\) gives \(N\ge(\deg E-1)a+3b>3b\) on every
interior ray. This handles every root of every initial polynomial there.

At \(x=0,t=\infty\), the coefficient of \(x^{-1}\) is divisible by
\(1-t^2\). In sphere normal form, the negative angle \(-1\) comes
from terms \(yt^j\), since the relation removes every product \(xy\).
The highest height \(q\) at that angle is therefore either two or three.
If \(q=3\), the monomial
\(x^{-1}t^3\) strictly dominates every other term on every interior
ray. If \(q=2\), all terms of nonnegative angle and height at most two
are dominated by \(x^{-1}t^2\). Among the height-three terms only the one
with least angular exponent can maximize. A tied face thus has two
support points with height difference one, so its Laurent polynomial
has simple nonzero roots. Here \(N\ge a+2b>b\), which proves convergence.

At a toric crossing choose the two neighboring rays. Some support
monomial maximizes on both rays because the fan refines the regions of
linearity of the maximum. It is unique, since the two independent
linear weights determine both of its exponents. In the corresponding
coordinates \((s_1,s_2)\), multiplying \(f\) by
\(s_1^{N_1}s_2^{N_2}\) gives a holomorphic function with nonzero value
at the crossing. Both orders are positive, also on the two boundary
rays, since a negative angular term and both a positive angular term
and a height-three term are present. Thus bounded phase values avoid
the crossing. These estimates and (1) give the required finite volume.
