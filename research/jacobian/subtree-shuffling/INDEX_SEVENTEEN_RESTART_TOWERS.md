# Explicit cubic restart growth at every index at least seventeen

Use the [fixed-coordinate scaffold and restart definitions](SUCCESSIVE_RESTART_UNIVERSALITY.md)
and the amalgamation and saturation arguments proved there.

## Compressed tower construction

For any fixed-coordinate scaffold with \(p\ge17\) and every \(M\ge0\), there is a
tower of \(M\) rounds adding \(17M\) coordinates with final inverse degree
\(\max\{D_*,2M+1\}\). If the scaffold has exact Jacobian index \(p\), the
final map has exact index \(p\) as well.

In the [compressed jet construction](CUBIC_INDEX_SEVENTEEN_INVERSE_SERIES.md#2-a-compressed-triangular-jet-family),
specialize its four fixed coordinates to
\((u_1,u_2,u_3,\tau)=(u,u,0,u)\). Equivalently, start its jet at \(uv\),
where \(v=(1,1,0,\ldots,0,1)\) is the vector of the closed inverse slice.
The specialized branch keeps \(u\) fixed and appends its 17 coordinates at
level \(j\) in round \(j\). The derivative of the nonlinear branch is the
restriction of \(tJH(Y(t))\) to the specialized tangent subspace. Its
seventeenth power is zero, and the same calculation applies to every
shorter jet. Fixed-coordinate amalgamation gives admissibility and the
exact-index statement when the scaffold supplies that index.

The positive jet weights still give inverse degree at most \(2M+1\).
At positive target levels zero, the inverse coefficient at level \(M\)
is \(G_M(uv)=u^{2M+1}G_M(v)\), nonzero by the closed inverse slice.
For \(M=0\), the branch is the identity on the shared \(u\). The inverse
degree of the union is therefore the stated maximum. The general
all-round equivalence also rules out a uniform ceiling for complete-basis
saturation at each \(p\ge17\). The index-four question still requires
control of weighted classes in its own Engel algebra.
