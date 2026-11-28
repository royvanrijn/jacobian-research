# A compact theta cannot have a thin folded middle interface

Let \(A\) be a balanced planar Gaussian noise-stability maximizer at
\(0<r<1\), with two bounded strict-winner cells and one unbounded
strict-winner cell. Suppose its entire actual interface is a compact
theta with two junctions joined by three regular pair arcs. Choose
orthonormal coordinates \((x,y)\), with origin at the Gaussian origin,
such that the two exterior arcs are respectively a lower graph and an
upper graph over the open interval of the two junction abscissae.
The middle arc can backtrack arbitrarily. If for almost every vertical
fiber with middle crossings their maximum minus minimum height is at
most \(2\sqrt{2\log 2}\sqrt{1-r^2}\), a mass-preserving planar speed
has negative second variation, contradicting maximality. The
width restriction concerns only the crossings of the middle arc, so
both bounded cells can extend far beyond the crossing interval.

## Signed intersections on a vertical fiber

Call the bounded cell touching the lower graph \(A_1\), the bounded
cell touching the upper graph \(A_2\), and the exterior \(A_3\).
Orient the arcs in the cycle \(3\to1\to2\to3\), with a unit normal
pointing out of the first label of each ordered pair. Put \(d_e=e_i-e_j\)
on a pair of orientation \(i\to j\), so \(d_e\cdot d_f=-1\) on
distinct arcs. The vertical translation speed \(s=e_y\cdot N\) is
nonnegative on the exterior graphs. Since they join distinct points
with distinct abscissae, their vertical speeds have positive Gaussian
surface integral. Their Euclidean measures \(s\,d\mathcal H^1\) push
forward to \(dx\) on their graph intervals.

The region between the two exterior graphs is the union of the two
bounded phases, by the Jordan separation property of the network.
In particular the middle arc stays inside the strip between the two
junction abscissae. On almost every interior vertical fiber the
middle crossings are finite and transverse, by the arclength area
formula and regularity of the arcs. The bottom cell is \(A_1\), the
top cell is \(A_2\), and between them only these labels occur. Thus
in order of increasing height the signs for \(s\) are \(+,-,+,\ldots,+\).
These signs are the fiber pushforward of the signed
arclength measure \(s\,d\mathcal H^1\) on the middle arc. No sign
claim is needed on the fibers of zero measure containing tangencies.

## An exact threshold for alternating Gaussian pairings

Write \(\sigma=\sqrt{1-r^2}\) and \(p_\sigma\) for the centered
one-dimensional Gaussian density of variance \(\sigma^2\). Consider
an increasing sequence \(v_1<\cdots<v_{2m+1}\) and the signed measure
\(\lambda=\sum_{\ell=1}^{2m+1}(-1)^{\ell+1}\delta_{v_\ell}\).
For a real normal observation mean \(\mu\), put \(k(y)=p_\sigma(y-\mu)\).
If \(\mu\) is outside \((v_1,v_{2m+1})\), then \(k\) is monotone
on the crossing interval and the alternating pairing is strictly
positive. This follows by retaining one endpoint value and pairing
the other terms in nonnegative consecutive differences. The case
\(m=0\) directly has positive pairing.

If \(\mu\) is inside the interval, the cumulative function of
\(\lambda\) takes values zero and one. Integrating by parts and
bounding the positive variation of \(k\) by \(k(\mu)-k(v_1)\) gives
\[
  \int k\,d\lambda\ \ge\ k(v_1)+k(v_{2m+1})-p_\sigma(0).
                                                               \tag{TF1}
\]
For \(u,v\ge0\) with \(u+v=w\), the minimum of
\(e^{-u^2/2}+e^{-v^2/2}\) occurs at an endpoint or at \(u=v=w/2\).
To see this, take \(0<u\le v=w-u\). The sign of its derivative in
\(u\) equals the sign of
\[
   H(d)=\log\frac{w+d}{w-d}-\frac{wd}{2},\qquad d=v-u .
\]
For \(0<d<w\), \(H'\) is strictly increasing, with
\(H'(d)=2w/(w^2-d^2)-w/2\), while \(H(0)=0\) and
\(H(d)\to\infty\) as \(d\uparrow w\). If \(H\) has a nonzero root,
its sign changes from negative to positive in increasing \(d\).
Since \(d=w-2u\), this root gives a maximum in \(u\). Thus the
minimum has the stated form. At \(w\le 2\sqrt{2\log 2}\), the endpoint
value is strictly greater than one and the midpoint value is at least
one. Apply this to \(u=(\mu-v_1)/\sigma\) and
\(v=(v_{2m+1}-\mu)/\sigma\). This proves nonnegativity in (TF1) under
the asserted width bound. Equality in the estimate also requires
that the cumulative function be one before \(\mu\) and zero after
\(\mu\) in the crossing interval. The zero case can occur at the
threshold only, with the three signs \(+,-,+\) at the endpoints
and midpoint. Thus the threshold is sharp for an individual Gaussian
pairing with arbitrary alternating crossing pattern.

## A strictly negative mass-preserving second variation

On one exterior arc let \(f=s\), and let \(f=0\) elsewhere.
For an observation \(p=(x,y)\) there, the middle contribution to
the scalar Mehler interaction of \(f\) with \(s-f\) is nonnegative.
Indeed \(M_r(p,q)d\nu(q)\) is the conditional Gaussian arclength
density at \(q\) with independent coordinates of means \(rx,ry\)
and common variance \(\sigma^2\). The signed middle fiber has the
pairing of (TF1) with mean \(\mu=ry\), integrated against the
positive horizontal conditional density. The interaction with the
other exterior arc is strictly positive in the scalar kernel.
Multiplying the cross interactions by the cyclic incidence factor
gives
\[
                 \langle f,K(s-f)\rangle<0 .               \tag{TF2}
\]
Since \(a s=rKs\) by translation, this is exactly
\(r^{-1}Q^{(1)}(f)\), with \(Q^{(1)}=\int a f^2d\nu-r\langle f,Kf\rangle\).
Thus the first extra-coordinate stability inequality fails.
[Translation and stability on normal speeds (G1)–(G4)](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#translation-equations-on-an-unconstrained-space-of-speeds)

The same calculation supplies a balanced planar direction. With
\(P\) the weighted orthogonal projection onto the two translation
speeds, put \(h=f-Pf\). Then \(Q^{(1)}(h)=Q^{(1)}(f)<0\). By
positivity of the Mehler quadratic form, \(Q_0(h)\le Q^{(1)}(h)\),
where \(Q_0(g)=\int a g^2d\nu-\langle g,Kg\rangle\). Let \(B\)
be the two-dimensional volume flux map and write \(s_v=v\cdot N\).
At an actual maximizer \(C:v\mapsto B s_v\) is invertible by full
active moment rank. The speed \(h-s_{C^{-1}Bh}\) is balanced.
Its two terms are \(Q_0\)-orthogonal by the translation identity.
Translation has nonpositive \(Q_0\), so the resulting balanced
direction is strictly negative. Cutting off near the finitely many
junctions and correcting the small flux error on fixed regular arcs
gives compact smooth balanced speeds with the same strict sign.
[Weighted index and translation flux](THREE_CELL_GLOBAL_INTERFACE_BOUNDS.md#a-global-matrix-bound-on-volume-flux)

There is also a quantitative lower bound on thick folds. Suppose
the network and phases it encloses are in a fixed disk of radius
\(R\), and \(r\) belongs to a fixed compact subinterval \(I\) of
\((0,1)\). Let \(D\) be the set of base abscissae where the middle
crossing spread exceeds \(2\sqrt{2\log2}\sigma\). The projection
interval of the lens has Gaussian mass at least \(2/3\), the mass of
its two bounded cells. Hence each exterior vertical flux integral
\(\int s\,d\nu\) is at least \(2\phi_1(R)/3\), where \(\phi_1\) is
the standard one-dimensional Gaussian density. The scalar
exterior/exterior pairing is uniformly positive, by the positive
minimum of \(M_r\) on the fixed disk product and this flux bound.
The signed middle pairing is nonnegative outside \(D\), and
everywhere it is bounded below by \(-p_\sigma(0)\) by (TF1) and
the monotone case. On the disk the horizontal conditional density
is bounded by a fixed multiple of \(\phi_1\), uniformly in \(I\).
Therefore the scalar pairing of one exterior arc with the whole
middle is bounded below by \(-C(R,I)\gamma_1(D)\), since the
exterior vertical flux is bounded above by \(\phi_1(0)\).
Nonnegative \(Q^{(1)}\) on that exterior arc then requires
\(\gamma_1(D)\ge c(R,I)>0\).

Thus any maximizer of this graph type needs thick middle folds.
Backtracking on the exterior arcs also remains a possible geometry.

## Affine centers for rotational cuts

For the rotation about the Gaussian origin write \(w=Jx\cdot N\),
where \(J\) is a right-angle rotation. Then \(Bw=0\) and \(Kw=aw\).
Moving the rotational center adds a translation \(s_v\). Write
\(g=w+s_v\), and let \(f_\ell=1_{D_\ell}g\) for a measurable
partition of the interfaces. Set \(I_{\ell m}=\langle f_\ell,Kf_m\rangle\).
With \(c=r^{-1}-1\), the row-sum identity
\(ag-Kg=-ca s_v\) gives
\[
 Q_0\!\left(\sum_\ell u_\ell f_\ell\right)
     =\sum_{\ell<m} I_{\ell m}(u_\ell-u_m)^2
           -c \sum_\ell u_\ell^2\langle s_v,f_\ell\rangle_a,
 \qquad
    B\!\left(\sum_\ell u_\ell f_\ell\right)=\sum_\ell u_\ell Bf_\ell .
                                                               \tag{TF3}
\]
Here \(\langle \cdot,\cdot\rangle_a\) uses \(a\,d\nu\). One exact
two-piece repair uses a measurable union \(D\) of pieces. Set
\(C_D v=B(1_Ds_v)\), \(b_D=B(1_Dw)\), and \(g_D=1_Dg\). If
\((C-2C_D)v=2b_D\), then \(g-2g_D\) is balanced and (TF3) gives
\[
               Q_0(g-2g_D)=4\langle g_D,K g_{D^c}\rangle
                                -c\|s_v\|_a^2 .           \tag{TF4}
\]
Thus a choice of shifted center meeting the two mass constraints
requires a cross interaction at least \(c\|s_v\|_a^2/4\) at a
stable network. Equation (TF3) makes no sign assumption on the
cross interactions.
[Rotation and translation identities](THREE_CELL_REFLECTION_REDUCTION.md#total-positivity-in-the-odd-reflection-sector)

## Analytic continuation for balanced null cuts

The double-bubble quadratic form is local, including the junction
terms. Splitting a rotation speed on nodal cuts
preserves zero energy for each piece, and volume balance followed by
unique continuation forces the rotation to vanish on the cut sheets.
The Mehler cross terms in (TF3) persist after splitting.
[Hutchings–Morgan–Ritoré–Ros, Proposition 5.2](https://math.berkeley.edu/~hutching/pub/bubble/bubble.pdf)

There is a nonlocal continuation statement after zero energy has
been established. On an actual compact maximizing theta, suppose
\(f\in L^2(a\,d\nu)\) satisfies \(Bf=0\) and \(Q_0(f)=0\).
Nonnegativity on balanced speeds gives \(Q_0(f,h)=0\) for every
balanced \(h\) by polarization. Hence the surface equation is
\[
                        af-Kf=\lambda\cdot d_e             \tag{TF5}
\]
for a fixed label vector \(\lambda\), by factoring the linear
functional \(Q_0(f,\cdot)\) through \(B\). The reciprocal-gradient
bound (G20) gives \(f\in L^1(d\nu)\) by Cauchy–Schwarz. The Gaussian
transform \(Kf\) is real analytic, since it is a Gaussian kernel
integrated against a compactly supported finite signed surface
measure. On every connected regular pair arc, \(a\) is positive and
analytic. Equation (TF5) then makes \(f\) analytic on that arc.

For a cut of \(w\) or any affine rotation \(g\), this forces all
coefficients on open subarcs of the same analytic edge to coincide
unless the original speed vanishes identically on the edge. Indeed
the analytic representative of \(f\) equals the stated multiple on
one subarc and therefore on the whole edge. A genuinely internal
balanced cut therefore has strictly positive \(Q_0\) at a stable
theta. This identifies zero energy, as well as negative energy, as
enough for a cut obstruction. The sign of the unrestricted
rotational cross matrix on genuinely folded compact theta networks
remains open.
[Reciprocal-gradient control of regular arcs (G20)](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#regular-interface-arcs-and-isolated-triple-junctions)
