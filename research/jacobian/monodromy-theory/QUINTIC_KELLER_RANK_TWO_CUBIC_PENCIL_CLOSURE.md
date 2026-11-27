# Invertibility for a cubic factor in the rank-two quintic pencil

**Theorem QuinticKellerRankTwoCubicPencilClosure.** Let \(F=(P,Q,R)\) be a
three-variable Keller map over a characteristic-zero field, with independent
quintic leading components \(P_5,Q_5\) and a quartic component \(R\).
If \(R_4\) has a geometrically irreducible cubic factor, then \(F\)
is invertible.

We work over \(\mathbb C\) and descend the inverse by uniqueness of the
formal inverse over the finitely generated coefficient field. By the
[rank-two quartic restriction](QUINTIC_KELLER_RANK_TWO_QUARTIC_PENCIL.md),
we may take \(P_5=z^5,\ Q_5=z^2T,\ R_4=zT\), with \(T\) irreducible
cubic. The [zero-\(\gamma\) proposition](QUINTIC_KELLER_RANK_TWO_CUBIC_PENCIL_ZERO_GAMMA.md)
covers that value of \(\gamma\), and the rank-two restriction excludes
\(L=0,\gamma\ne0\). Thus use \(x=L\ne0\), with \(T(x,y,0)\) a nonzero
multiple of \(x^3\). The quartic and cubic formulas of the rank-two
restriction put \(P\) in the form
\[
  P=c x^3+A(z)x^2+B(z)x+C(z)+y(D(z)x+E(z))+H(z)y^2
  \tag{1}
\]
with \(c\ne0\) and
\[
 \deg A,\deg D,\deg H\le1,\quad \deg E\le2,\qquad
 B(z)=5z^3+O(z^2),\quad C(z)=z^5+O(z^4).
\]
In particular \(P\) has no critical point. We show that the only
possible case of (1) has \(D=H=0\) and \(E\) a nonzero constant.
Then \(P\) is a graph coordinate with inverse completion of degree
at most five, so the [bounded coordinate-completion theorem](LOW_DEGREE_COORDINATE_KELLER_MAPS_IN_THREE_VARIABLES.md)
gives invertibility.

## A nonzero coefficient of \(y^2\)

Write \(P_0=c x^3+A x^2+B x+C\). If \(H\ne0\), at \(H(z)\ne0\)
the critical-point equations reduce by \(P_y=0\) to those of
\[
               \widetilde P(x,z)=P_0(x,z)
                                 -\frac{(D(z)x+E(z))^2}{4H(z)}.
\]
Write \(\widetilde A,\widetilde B,\widetilde C\) for its coefficients
of \(x^2,x,1\), respectively. Consider the rational function of \(z\)
given by the product of \(\widetilde P_z\) at the two roots,
with multiplicity, of \(\widetilde P_x=0\). It is rational by
symmetry and has no finite pole outside the zero set of \(H\):
the quadratic \(\widetilde P_x\) has leading coefficient \(3c\),
so its roots stay bounded at any point with \(H\ne0\).
A zero at such a point gives a critical point of \(P\).

If \(H\) is constant and \(D\) has degree one, then
\(\widetilde A\sim a_2 z^2\) with \(a_2\ne0\),
\(\deg\widetilde B\le3\), and \(\widetilde C\sim z^5\).
The two roots of \(\widetilde P_x\) have one leading term
\(-2a_2 z^2/(3c)\) and the other of order at most \(z\).
At the first root \(\widetilde P_z\) has order \(z^5\) with
nonzero leading coefficient from \(\widetilde A' x^2\).
At the second it has leading term \(5z^4\).
If \(H\) is constant and \(D\) has degree at most zero, then
\(\widetilde A=O(z)\), \(\widetilde B\sim5z^3\),
\(\widetilde C\sim z^5\); the roots are \(O(z^{3/2})\) and both
give \(\widetilde P_z\sim5z^4\).
Thus in either constant case the product has a pole at infinity
and no finite pole, so it has a finite zero.

If \(H\) is linear, then again
\(\widetilde A=O(z)\), \(\widetilde B\sim5z^3\),
\(\widetilde C\sim z^5\), and the product has a pole at infinity.
Translate the root of \(H\) to \(s=0\), so \(H=\lambda s\).
If \(D(0)\ne0\), the equation
\[
              s(P_0)_x-\frac{D(Dx+E)}{2\lambda}=0
\]
has one holomorphic root \(x_0(s)\) through \(x=-E(0)/D(0)\).
Its derivative with respect to \(x\) at \(s=0\) is nonzero.
The corresponding \(y=-(D x_0+E)/(2H)\) extends holomorphically,
so the value of \(\widetilde P_s\) there is finite, and is
nonzero at zero if \(P\) has no critical point. The second root
has a simple pole in \(s\), by the root sum of this quadratic
in \(x\). At this root the term \((Dx+E)^2/(4\lambda s^2)\)
in \(\widetilde P_s\) gives a pole of order four. All other terms
have smaller pole order. Thus the product has a pole rather
than a zero at \(s=0\).

If \(D(0)=0\) and \(E(0)\ne0\), the two roots stay bounded at
zero, because \(\widetilde A,\widetilde B\) are regular there.
Both values of \(\widetilde P_s\) have poles of order two
from \((Dx+E)^2/(4\lambda s^2)\), so the same conclusion holds.
If \(D(0)=E(0)=0\), then at \(s=0\) the \(y\) derivative
vanishes. Solve the quadratic \(P_x=0\) in \(x\) and then
\(P_s=0\) in \(y\), which has nonzero quadratic coefficient
\(\lambda\). This gives a critical point directly.

In the linear \(H\) cases without a critical point already found,
the rational product therefore has poles at infinity and at
the only possible finite pole. It must have a zero elsewhere.
All possibilities with \(H\ne0\) are impossible.

## A nonzero \(xy\) coefficient

Now \(H=0\). Suppose \(D\ne0\). At \(D(z)\ne0\) set
\(x=\xi(z)=-E(z)/D(z)\), and let \(K(z)=P_0(\xi(z),z)\).
A zero of \(K'\) gives a critical point of \(P\) by choosing
\(y\) from \(P_x=0\), since \(D'\xi+E'=-D\xi'\).

If \(D\) is constant, then \(K\) is a polynomial of degree
six when \(\deg E=2\), with top term from \(c\xi^3\), and is
monic of degree five otherwise. Thus \(K'\) has a root.
If \(D\) is linear and does not divide \(E\), then \(K\) has a
pole of order three at the root of \(D\) and has leading
term \(z^5\) at infinity. Its derivative has no zero at that
pole or at infinity and has a zero elsewhere. Finally,
if \(D\) is linear and divides \(E\), then \(K\) is monic
of degree five. If \(K'\) has a root away from the root \(z_0\)
of \(D\), that root again gives a critical point. Otherwise
\(K'(z_0)=0\). At \(z_0\), both \(D,E\) vanish, and any root
of the quadratic \((P_0)_x\) different from \(\xi(z_0)\)
permits solving \(P_z=0\) for \(y\), since its \(y\) coefficient
is \(D'(z_0)(x-\xi(z_0))\). If every root equals \(\xi(z_0)\),
then \((P_0)_x=0\) there and \(P_z=K'(z_0)=0\).
In all cases a critical point occurs, so \(D=0\).

## The remaining affine coefficient of \(y\)

With \(H=D=0\), if \(E=0\) the constant-\(H\) binary resultant
argument applies to \(P_0\) itself: at either root of \((P_0)_x\),
\((P_0)_z\sim5z^4\) at infinity, so their product is a nonconstant
polynomial with a zero. Thus \(E\ne0\).
At any simple root of \(E\), one can solve the quadratic \(P_x=0\)
and then solve \(P_z=0\) for \(y\). Therefore if \(E\) is
nonconstant it must have a repeated quadratic root.

Translate that root to \(s=0\), with \(E=\lambda s^2\).
Use the same exact-residue argument as in the
[zero-\(\gamma\) proposition](QUINTIC_KELLER_RANK_TWO_CUBIC_PENCIL_ZERO_GAMMA.md#a-repeated-root-of-the-coefficient-of-y),
now on the three holomorphic arcs \(y=0,\ x=x_i(s)\)
of a generic fibre. They arise from the three distinct roots of
\(P_0(x,0)=u\). The regular primitive on the fibre and comparison
with \(x\,ds/(\lambda s^2)\) show that all three derivatives
\(x_i'(0)\) coincide, say with value \(\mu\). Differentiating the
fibre equation gives \((P_0)_s(x_i(0),0)+\mu(P_0)_x(x_i(0),0)=0\).
Both derivatives are polynomials of degree at most two in \(x\)
because the coefficient \(c\) of \(x^3\) is constant.
The three distinct roots force this polynomial identity for all \(x\).
At a root of \((P_0)_x(x,0)\), the derivatives of \(P\) therefore
all vanish, a contradiction.

Thus \(E\) is a nonzero constant and \(P\) is the graph coordinate
claimed above. This completes the cubic-factor case.
