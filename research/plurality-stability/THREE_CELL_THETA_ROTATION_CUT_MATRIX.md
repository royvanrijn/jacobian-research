# Cutting the rotation on a compact Gaussian theta network

Consider a planar three-cell stationary partition whose interface is an
embedded compact theta network. Its two vertices are denoted by \(P,Q\)
and its three arcs have distinct pair types. Cutting the rotational
speed gives an exact finite matrix test on the two mass constraints.
For a fixed shape with at least four connected pieces after cuts at
regular sign changes, the matrix has a negative constrained direction
for correlations sufficiently near one, provided the shape satisfies
the stationary equations at the tested correlation. The range
depends on the geometry. Backtracking in every linear projection
does not itself force four connected cut pieces, even for balanced
smooth theta networks with transverse vertices.

The related argument in the local double-bubble problem cuts a Killing
Jacobi field at zeros into four separating pieces. The cut functions
have zero mutual terms, so volume constraints and stability give
a new zero field, to which unique continuation applies. Here the
cuts have signed mutual terms. Their sign is the content of the test.
[Hutchings, Morgan, Ritoré and Ros, Proposition 5.2](https://math.berkeley.edu/~hutching/pub/bubble/bubble.pdf)

## The exact matrix and its mass constraints

Orient every arc from \(P\) to \(Q\), with unit tangent \(\tau\)
and normal \(N=J\tau\), where \(J\) rotates the plane by \(90^\circ\).
The label vector \(d_j=e_i-e_k\) of arc \(j\) uses the labels from
which and into which \(N\) points. The three vectors form a cyclic
triple, so \(\sum_jd_j=0\) and any two are independent. Put
\(\zeta(x)=Jx\cdot N(x)=x\cdot\tau(x)\) and
\(d\nu=\phi_2\,ds\), with \(\phi_2\) the standard planar Gaussian
density. Use the normalization
\[
 Q_0(f)=\int a f^2\,d\nu-\iint K_r(x,y)f(x)f(y)\,d\nu(x)d\nu(y),
 \qquad K_r(x,y)=M_r(x,y)d_{j(x)}\cdot d_{j(y)} .
\]
The coefficient \(a\) is the negative outward normal derivative of
the pair score difference. Rotating the entire stationary partition
keeps the pair pressures fixed. Differentiating the score equations
on the moving arcs gives \(K_r\zeta=a\zeta\).
[Stationary scores and interface kernel](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#translation-equations-on-an-unconstrained-space-of-speeds)

Let \(D_1,\ldots,D_m\) partition the arcs up to arclength-null sets.
For \(f=\sum_i\alpha_i1_{D_i}\zeta\), put
\[
 B_i=\sum_j d_j\int_{D_i\cap\Gamma_j}\zeta\,d\nu,\qquad
 C_{ik}=\iint_{D_i\times D_k}K_r(x,y)\zeta(x)\zeta(y)\,d\nu(x)d\nu(y).
\]
The constraint and the value of the second variation are
\[
             \sum_i B_i\alpha_i=0,\qquad
             Q_0(f)=\sum_{i<k}C_{ik}(\alpha_i-\alpha_k)^2.       \tag{TC1}
\]
Indeed the Jacobi equation replaces the local term on \(D_i\) by
\(\alpha_i^2\sum_k C_{ik}\). Symmetry gives the difference form.
Bounded piecewise speeds can be approximated by smooth speeds
supported on regular interiors. Two fixed bumps on independent
pairs correct the flux errors, and the local and Mehler terms
converge on the compact interfaces.

If a piece contains an oriented subarc from \(u\) to \(v\), its
contribution to \(B_i\) is exactly
\[
                         d_j\,[\phi_2(u)-\phi_2(v)].              \tag{TC2}
\]
This follows from \(d\phi_2(x(s))/ds=-\zeta(x(s))\phi_2(x(s))\).
For the cut into whole arms the flux matrix is
\((\phi_2(P)-\phi_2(Q))[d_1,d_2,d_3]\). If the vertex radii differ,
only the constant coefficients satisfy the mass constraints, giving
the global zero mode. If the vertex radii are equal, every whole-arm
test is balanced and its sign is given by (TC1).

For one cut of an arm into pieces \(D,E\), denote their scalar fluxes
by \(u,v\), so \(u+v=\lambda=\phi_2(P)-\phi_2(Q)\), and let \(R\)
be the union of the other two arms. If \(\lambda\ne0\), mass balance
forces a common coefficient on those two complete arms. Subtract
the corresponding multiple of the global rotation. The remaining
constrained coefficient line is represented by
\(f=v1_D\zeta-u1_E\zeta\), with
\[
 Q_0(f)=\lambda^2 C_{DE}+v^2 C_{DR}+u^2 C_{ER}.            \tag{TC5}
\]
Here \(C_{XY}\) is the same interaction integral for unions of pieces.
For one cut on each arm at the same intermediate radius, each
arm has the same \(u,v\). After subtracting the global rotation,
the three functions \(g_j=v1_{D_j}\zeta-u1_{E_j}\zeta\) span the
entire constrained cut space modulo that rotation. The compressed
\(3\)-by-\(3\) matrix has off-diagonal entries
\(-\langle g_i,K_r g_j\rangle\). Its diagonal entries are (TC5)
for each arm, with \(R\) the union of the other arms, by the
row sums of the underlying signed Laplacian (TC1).

At a vertex \(p\), the cyclic pair gradients sum to zero, giving
\[
                           \sum_j a_j(p)\zeta_j(p)=0.           \tag{TC3}
\]
At a noncentral noncollinear vertex the rotation values thus have
both signs. If \(|P|<|Q|\), each arc has positive scalar rotation
flux by (TC2). When both vertices are noncentral and noncollinear,
there are at least two interior sign changes, counted across the
arcs. If the negative values at the vertices belong to one arc it
must turn positive and then negative again. Otherwise two arcs each
have a sign change. If \(|P|=|Q|\), every noncircular arc has both
signs by (TC2) and continuity.
[Transverse and collinear junctions](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#regular-interface-arcs-and-isolated-triple-junctions)

## Sign control when cuts outnumber the constraints

Select a finite set of interior sign-changing zeros of \(\zeta\),
each on a regular analytic arc and of finite order. Take the
\(D_i\) to be the connected components after removing these points.
The three germs at a true vertex belong to a single component.
Suppose \(m\ge4\). There are at most two mass constraints, and
the constants have zero flux by rotation invariance. Choose
\(\alpha\) of Euclidean norm one in the flux kernel perpendicular
to the constants.

The graph joining pieces across the selected cuts is connected.
At some cut between pieces \(i,k\),
\(|\alpha_i-\alpha_k|\ge 1/(\sqrt m(m-1))\). Indeed some coefficient
has magnitude at least \(m^{-1/2}\), the zero mean gives another
of opposite or zero sign, and a path of at most \(m-1\) edges
joins them.

For distinct pieces, every positive contribution of
\(K_r(x,y)\zeta(x)\zeta(y)\) has \(|x-y|\ge\delta>0\) for a
constant of the fixed cut geometry. Otherwise positive pairs at
vanishing separation would limit to a common point of two piece
closures. This point must be a selected cut. Its two open sides
have the same pair type and opposite rotation signs near the cut,
making their integrand negative.

Set \(\sigma=\sqrt{1-r^2}\), and let \(\ell\ge1\) be the largest
order of the selected zeros. At a cut with the stated jump bound,
retain on each side the arclengths from \(\sigma/2\) to \(\sigma\)
from the cut. Their products have opposite signs and
\(|\zeta(x)\zeta(y)|\ge c\,\sigma^{2\ell}\) for small \(\sigma\).
Their distances are at most \(2\sigma\). The kernel against
arclength is
\[
 M_r(x,y)\phi_2(x)\phi_2(y)=
 \frac{1}{(2\pi)^2\sigma^2}
 \exp\left[-\frac{r|x-y|^2+(1-r)(|x|^2+|y|^2)}
                   {2\sigma^2}\right].
\]
On this patch product it is at least \(c'\sigma^{-2}\).
The product arclength is \(\sigma^2/4\), so the negative term
is at most \(-c''\sigma^{2\ell}\). All positive terms together
are at most \(C\sigma^{-2}e^{-r\delta^2/(2\sigma^2)}\) by
separation, bounded speeds and finite arclength. Thus
\[
 Q_0(f)\le -c''\sigma^{2\ell}
                  + C\sigma^{-2}e^{-r\delta^2/(2\sigma^2)}.     \tag{TC4}
\]
This is strictly negative for \(r\) sufficiently close to one.
For varying geometries the proof requires control of the
separation, cut orders and leading coefficients, and neighboring
arc lengths. At a fixed intermediate correlation the sign does
not follow from the piece count alone.

## A balanced geometric theta with only two rotation cuts

Choose \(0<a<b\) so \(\int_a^b s e^{-s^2/2}ds>2/3\). For large
\(k\), put \(w(s)=(1-e^{-k(s-a)})(1-e^{-k(b-s)})\) and
\(I=\int_a^b s e^{-s^2/2}w(s)\,ds>2/3\).
Let \(\delta(s)=2\pi w(s)/(3I)\), which is less than \(\pi\).
Take a smooth \(\theta(s)\) going from \(0\) to \(2\pi\),
constant near each endpoint. In polar coordinates the three arcs
\((s,\theta(s)+j\delta(s))\), \(j=-1,0,1\), form a theta graph
between \(ae_1\) and \(be_1\). The two bounded angular strips
have Gaussian masses \(1/3\) by the definition of \(\delta\).
Every arm contains points in directions \(\pi/6,5\pi/6,3\pi/2\).
Their convex hull contains the origin in its interior, so every
linear projection takes both signs on each arm. The endpoint
projections have the same sign or are both zero. Thus no arm is
monotone under any nonzero linear projection.

The initial vertex tangents lie in an open half-plane. The
following local modifications make all three opening angles at
each vertex less than \(\pi\) and introduce two rotation zeros.
Near \(P=ae_1\), leave the middle and lower arms fixed and replace
the upper germ by the polar curve
\[
    s(t)=a+t^2-\epsilon t,\qquad \psi(t)=ct,\qquad
    c=(\epsilon/a)\tan(\beta/2),\quad
    \beta=\arctan(a\delta'(a)),\quad 0\le t\le\epsilon.
\]
For small \(\epsilon>0\) it stays at positive radius. It has one
radial minimum and its outgoing tangent has angle
\(\pi-\beta/2\). The other two tangent angles are \(0,-\beta\),
so the three cyclic opening angles are less than \(\pi\).
At \(t=\epsilon\) the curve is at radius \(a\) and positive
angle, with strictly positive radial derivative. Extend it
as a radial graph over a short interval just beyond \(a\),
and blend its angle with \(\delta(s)\) using a smooth cutoff
equal to one near the join and zero beyond that interval.
Both angles are positive there, since \(\theta=0\) near the
endpoint. Hence it meets neither of the other arms and joins
the original one smoothly. Its support tends to \(P\).

Near \(Q=be_1\), leave the middle and upper arms fixed and
parametrize a replacement lower germ outward from \(Q\) by
\(s(t)=b+\epsilon t-t^2,\ \psi(t)=-c't\) for \(0\le t\le\epsilon\),
where \(c'=(\epsilon/b)\tan(\beta'/2)\) and
\(\beta'=\arctan(-b\delta'(b))\).
This has one radial maximum and initial angle \(-\beta'/2\).
The remaining two outgoing angles are \(\pi\) and \(\pi-\beta'\),
so again every opening angle is less than \(\pi\). Rejoin on
the interior radial interval as at \(P\), retaining a negative
angle relative to the endpoint direction. The polynomial germs
have strictly monotone angular coordinates, so their radial
turns remain embedded.

The two surgeries change the bounded masses by amounts tending
to zero. Correct them separately by angular bumps on the upper
and lower radial graphs away from the endpoints. The Gaussian
mass changes there are exactly the signed integrals of those
bumps against \(s e^{-s^2/2}/(2\pi)\), so two small coefficients
restore the masses. The angle ordering and the three-point
backtracking witnesses persist. These corrections leave the
radius parametrizations and their signs unchanged.

The balanced smooth theta graph now has exactly two interior
zeros of the rotation, one on each modified arm. Removing them
leaves a connected graph via the uncut middle arm. Thus
backtracking and the transverse tangent configuration do not
give the four connected pieces required by (TC4). This example
has no global stationary score or Jacobi equation imposed and
supplies no claim of a stable theta partition.
