# Low-branch collisions at the high pole of the trigonal pencil

Use the normalized quartic-leading winding-one cubic of exact least
sphere degree seven from the
[high-pole Lefschetz exclusion](WINDING_ONE_CUBIC_QUADRUPLE_ROOT_GENERIC_EXCLUSION.md),
with \(U=1+x\), \(v=Ut\), and \(d_3\ne0\). Put \(u=x/U\),
\(\tau=u-1\), and
\[
 D(\tau)=\tau^3d(x),\qquad B(\tau)=\tau^5b(x),\qquad
 C(\tau)=\tau^7c(x),\qquad x=-\frac{1+\tau}{\tau}.
\]
For the finite value parameter \(\zeta\), set
\[
 H(\tau,T;\zeta)=\tau T^3+D(\tau)T^2+B(\tau)T+C(\tau)
                  -\zeta(1+\tau)\tau^7.
\]
Let \(k\) be the order at zero of \(\operatorname{disc}_T H(\tau,T;0)\),
with \(k=\infty\) if that discriminant vanishes identically.

**Theorem (Finite low-branch collision exclusion).** Suppose that
the conditions of the
[Lefschetz exclusion theorem](WINDING_ONE_CUBIC_QUADRUPLE_ROOT_LEFSCHETZ_EXCLUSION.md)
at \(u=0,\infty\) and on the interior Morse critical values hold.
If \(k<7\), then infinitely many positive pure moments are nonzero.
If \(k\ge7\), write the low quadratic in centered coordinates as
\(S^2=\tau^7 K(\tau,\zeta)\). There is a unique finite \(\zeta_*\) with
\(K(0,\zeta_*)=0\). The same nonvanishing conclusion holds if
\(K_\tau(0,\zeta_*)\ne0\), including when \(\zeta_*\) coincides with
zero or with an interior critical value.

**Proposition (Reduced coefficients in the remaining collision chart).**
For \(k\ge7\), there is a polynomial \(S_0(x)\) of degree at most two
such that translating \(v=z+S_0(x)\) makes the cubic equal to
\(z^3+\widetilde d(x)z^2+\widetilde b(x)z+\widetilde c(x)\), with
\(\deg\widetilde d=3\), \(\deg\widetilde b\le2\), and
\(\deg\widetilde c\le1\). The leading coefficients obey
\(\widetilde b_2^2=4\widetilde d_3\widetilde c_1\).
The quadratic part of this translation can be removed by a sphere
shear, leaving an affine translation.

## The quadratic at infinity

All the analytic arguments can be made over \(\mathbb C\) by
embedding the finitely generated coefficient field. The change
\(W=\tau^3v\) puts the fiber equation near \(u=1\) in the form
\[
 W^3+D W^2+\tau B W+\tau^2 C=\zeta(1+\tau)\tau^9.
\]
Its root \(W_h\) with value \(-D(0)\) is simple and analytic near
every finite value of \(\zeta\). For the remaining roots put
\(W=\tau T\). Then \(H=(\tau T-W_h)Q\) for a monic quadratic \(Q\)
with analytic coefficients. Its linear coefficient is
\((D+W_h)/\tau\), which is regular, and the constant coefficient
is regular because \(W_h\) is a unit. Centering the quadratic
writes it as \((T-A)^2-M\), where \(M\) is one quarter of its
discriminant.

The cubic discriminant differs from the quadratic discriminant
by a unit with value \(D(0)^4\) at \(\tau=0\). Indeed, if \(T_1,T_2\)
are the two bounded roots, the other contribution is
\(((W_h-\tau T_1)(W_h-\tau T_2))^2\). This is analytic and a unit by
symmetry in the roots. Differentiating the cubic discriminant with
respect to \(\zeta\) gives order at least seven in \(\tau\), with
coefficient \(4D(0)^3\) at order seven. This follows directly from
the term \(-4D^3(C-\zeta(1+\tau)\tau^7)\); every other derivative term
has an additional factor of \(\tau\).

For \(k<7\), \(M=\tau^k K\) where \(K(0,\zeta)\) is a nonzero
constant. If \(k\) is even, the bounded roots give two smooth
sections \(T=A\pm\tau^{k/2}\sqrt K\) in the normalized model. If
\(k\) is odd, they give one smooth chart with coordinate \(\eta\),
where \(\tau=\eta^2\) and \(T=A+\eta^k\sqrt K\). The two square-root
choices give the same chart by changing \(\eta\) to \(-\eta\).
In both cases the projection to the finite value line is smooth.
These charts are finite over \((\tau,\zeta)\) and their meromorphic
roots define maps to the projective \(Y\)-line, where \(Y=v+d/3\).
The map is regular across \(\tau=0\) because the high and low
depressed roots have nonzero leading pole coefficients there. The
charts are finite and birational onto the corresponding graph germs
after taking their images, hence give the normalization.

The inner part of a large regular fiber stays connected by the
[same scaled limit](WINDING_ONE_CUBIC_QUADRUPLE_ROOT_GENERIC_EXCLUSION.md#connectivity-of-the-inner-part-of-a-large-fiber)
\(W^2(W+D(0))=s^9\). Its irreducibility follows from the simple zero
at \(W=-D(0)\) in the equation for \(s\). Thus finitely many loops
avoiding the branch points in the scaled coordinate give a transitive
permutation group on the three roots, and these permutations persist
for sufficiently large values. The outer part is connected by the
totally ramified point over \(u=\infty\). Since the marked sections
at \(u=1\) are smooth and avoid finite critical points, the
[finite-node lemma](WINDING_ONE_CUBIC_QUADRUPLE_ROOT_LEFSCHETZ_EXCLUSION.md#2-the-vanishing-cycles-span-compact-homology)
applies and proves the first assertion.

## A finite value with two boundary nodes

If \(k\ge7\), the centered low quadratic has \(M=\tau^7 K\), with
\(K(0,\zeta)\) affine of nonzero slope by the discriminant calculation.
The unique zero is \(\zeta_*\). Dividing \(T-A\) by \(\tau^3\) gives
the normalized chart
\[
          L^2=\tau K(\tau,\zeta).
\]
It is smooth except for one ordinary surface double point at
\((\tau,L,\zeta)=(0,0,\zeta_*)\). Indeed \(K\) is a transverse
coordinate there and the equation is \(L^2=\tau K\). This chart
is normal and finite over \((\tau,\zeta)\), and as above maps
birationally to the corresponding graph germ. Resolve the point
using the charts \(L=\tau h,\ K=\tau h^2\) and
\(L=K\widetilde h,\ \tau=K\widetilde h^2\). The exceptional curve
is a projective line. The strict marked section \(\tau=L=0\)
meets it in the second chart at \(\widetilde h=0\), where the map
to the value line is smooth because \(K_\zeta(0,\zeta_*)\ne0\).

Put \(a=K_\zeta(0,\zeta_*)\) and \(b_*=K_\tau(0,\zeta_*)\). If
\(b_*\ne0\), at \(\zeta_*\) the exceptional curve \(\tau=0\) meets
the strict fiber in two Morse nodes, at \(h=\pm\sqrt{b_*}\). In
fact the fiber equation is \(\tau G(\tau,h)=0\) with
\(G(0,h)=(h^2-b_*)/a\), so both crossings are transverse. For
nearby values the two lifts of a fixed positive small \(\tau\)-circle
with \(h_\pm=\pm\sqrt{K(\tau,\zeta)/\tau}\) are their two vanishing
cycles. Their sum is zero in compact homology, by the boundary of
the smoothed part of the exceptional curve. The two orientations
agree in this relation because \(\zeta-\zeta_*\) is the product of
\(\tau\) and a local coordinate along that curve at either node.

The period form is \(\alpha=Y\,d\tau/(1+\tau)\). Since
\(T=A+\tau^3L=A+\tau^4h\), on the two lifts it has a common term
depending only on \(\tau,\zeta\), and the remaining term is
\(\tau^2 h\,d\tau/(1+\tau)\). Their period difference is thus
\[
 I_+(\zeta)-I_-(\zeta)=2\oint\frac{\tau^2}{1+\tau}
                 \sqrt{K(\tau,\zeta)/\tau}\,d\tau.
\]
In its expansion at \(\zeta_*\), the terms of orders zero, one and
two have no residue. At order three the only contribution of pole
order three before the \(\tau^2\) factor comes from cubing the
linear variation of \(K\). The coefficient in the period is
\[
             4\pi i\binom{1/2}{3}\,a^3 b_*^{-5/2}\ne0,
\]
with the square root determined by the chosen lift. Thus the joint
positive monodromy twist of the two nodes detects their nonzero
cycle direction: the two intersection numbers of any transported
cycle are opposite, so its period variation is a multiple of the
nonzero difference just computed.

The exceptional nodes are the only new finite singularities. If
their value equals zero, the three nodes over \(u=\infty\) give an
additional summand in the period variation. Its two independent
intersection coefficients are detected by the constant and linear
terms of their periods, as in the
[three-node calculation](WINDING_ONE_CUBIC_QUADRUPLE_ROOT_LEFSCHETZ_EXCLUSION.md#4-independence-at-the-three-nodes-over-zero).
The new difference vanishes to order three, so those two tests
first give orthogonality to the old cycles. If the exceptional value
coincides with an interior critical value, the existing hypotheses
allow only one interior node there. Its vanishing period has zero
constant term and nonzero linear term by the local Gelfand--Leray
calculation. That linear term first detects the intersection with
its cycle. In each collision case, the remaining order-three test
detects the exceptional cycle. The local twists act simultaneously
along disjoint annuli, including when their values coincide.

The resolution remains rational, and all marked sections avoid
the nodes. The argument of the
[Lefschetz proof](WINDING_ONE_CUBIC_QUADRUPLE_ROOT_LEFSCHETZ_EXCLUSION.md#2-the-vanishing-cycles-span-compact-homology)
that the finite vanishing cycles span compact homology is unchanged:
the rational compactification has zero restriction in first
cohomology to a smooth fiber, and each node contributes the
two-handle along its vanishing cycle. Two nodes over one value
have disjoint local Morse charts and can be separated by the same
smooth perturbation of values used for the three nodes over zero.
The joint period test above supplies the needed orthogonality to
both new cycles. The inner and outer parts of a large fiber are
connected as before. Therefore the inverse-branch cycle is nonzero
and cannot have an identically zero period. The
[inverse-period condition](WINDING_ONE_CUBIC_QUADRUPLE_ROOT_INVERSE_PERIOD.md)
excludes eventual purity.

## Coefficients after translating the common low part

At \(\zeta=\zeta_*\) the two low roots \(v=(A\pm\sqrt M)/\tau^2\)
share the polynomial part \(S_0(x)\) of \(A/\tau^2\). Indeed
\(\sqrt M/\tau^2\) is \(O(\tau^{3/2})\), and \(\tau^{-1}=-(1+x)\).
After the translation their two roots are \(O(x^{-1})\) and the
heavy root is \(O(x^3)\). Viète's identities for the translated
fiber polynomial give \(\widetilde b=O(x^2)\) and
\(\widetilde c-\zeta_*x/U=O(x)\). Since the value term is bounded
at \(x=\infty\), these are exactly the asserted degree bounds.

The discriminant of the polynomial in \(v\) is \(\tau^{-16}\)
times that of \(H\), by scaling its three roots and leading
coefficient. If \(k\ge7\), it has growth at most \(x^9\) at every
finite value. In translated coordinates the coefficient of \(x^{10}\)
is \(\widetilde d_3^2(\widetilde b_2^2-4\widetilde d_3\widetilde c_1)\).
This proves the asserted relation. A quadratic term of \(S_0(x)\)
can be absorbed by adding a multiple of \(xU\) to the translation,
which is the volume-preserving sphere shear \(t\mapsto t+\lambda x\)
with \(x\) fixed. Its trace invariance is the
[winding-one shear identity](MIXED_CUBIC_ZERO_CENTER_AND_CONDITIONAL_OFFCENTER.md).

The second case requires \(K_\tau(0,\zeta_*)\ne0\). Its replacement
by a higher-order zero on the special fiber remains open.
