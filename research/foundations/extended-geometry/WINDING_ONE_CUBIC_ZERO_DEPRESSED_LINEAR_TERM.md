# The zero depressed linear term in the degree-four boundary family

Let \(f=x^{-1}(C+Bt+Dt^2+Et^3)\) be a regular winding-one phase on
\(2xy+t^2=1\) over a field of characteristic zero. Suppose its least
sphere degree is seven and \(E=e(1+\lambda x)^4\), with \(e\lambda\ne0\).
Write \(U=1+\lambda x\) and assume the
[boundary divisibilities](WINDING_ONE_CUBIC_LOGARITHMIC_BOUNDARY.md#a-quadruple-root)
\(D=U^3d\), \(B=U^2b\), and \(C=Uc\) with polynomial \(d,b,c\).

**Proposition.** If the depressed linear coefficient of
\(e v^3+d(x)v^2+b(x)v+c(x)\) is zero, then the normalized sphere
moments \(\Lambda(f^m)\) are nonzero at infinitely many positive orders.

Scale \(x\) and \(y\) inversely and multiply the phase by a nonzero
scalar to take \(U=1+x\) and \(e=1\). These operations preserve sphere
degree and the zero or nonzero status of every moment. Write
\[
 p=b-\frac{d^2}{3},\qquad q=c-\frac{bd}{3}+\frac{2d^3}{27}
\]
for the depressed coefficients. The
[sphere degree formula](MIXED_CUBIC_ZERO_CENTER_AND_CONDITIONAL_OFFCENTER.md)
gives \(\deg b\le5\) and \(\deg c\le7\). Since \(p=0\),
\(b=d^2/3\) and \(\deg d\le2\). The \(D,B,E\) terms then all have
sphere degree at most six. Therefore \(\deg c=7\) and \(\deg q=7\).
Regularity gives \(b(0)=-1\) and \(c(0)=-d(0)\), so
\(d(0)^2=-3\) and \(q(0)=-8d(0)/9\ne0\).

Embed the coefficient field in \(\mathbb C\) after restricting to a
finitely generated subfield. Put \(s=(x/U)^{1/3}\) with \(s=1\) at
\(x=\infty\), and use \(v=Ut\). If \(f(x,T(x,w))=w^3\) for the inverse
branch at \(t=\infty\), then
\[
 T(x,w)=\frac{s}{U}\bigl(w^3-P(s)\bigr)^{1/3}-\frac{d(x)}{3U},
 \qquad P(s)=R(s^3),\qquad R(u)=\frac{q(u/(1-u))}{u}.
\]
The cube root has leading term \(w\). Every coefficient of
\(w^{1-3k}\), \(k\ge1\), is a nonzero rational multiple of
\((s/U)P(s)^k\). Since \((s/U)\,dx/x=3\,ds\), the
[inverse-period condition](WINDING_ONE_CUBIC_QUADRUPLE_ROOT_INVERSE_PERIOD.md)
would give \(\operatorname{Res}_{s=1} P(s)^k\,ds=0\) for every \(k\ge1\)
if the original moments vanished eventually.

The rational map \(R\) has poles of orders seven at \(u=1\) and one at
\(u=0\). It is regular at \(u=\infty\). For a large regular value \(z\)
of \(R\), let \(S(z)\) be the sum of the seven inverse branches of
\(P(s)=z\) near \(s=1\). They assign the cube root near one to each of
the seven inverse branches of \(R(u)=z\) near \(u=1\). Their sum is a
single-valued holomorphic function near target infinity and tends to
seven there. For sufficiently large \(L\), the local inverse image of
\(|z|=L\) is one closed contour about \(s=1\), covering the target circle
seven times. Change variables around these seven lifts. With the target
circle oriented clockwise to match the positive local contour, the
residue of \(P(s)^k ds\) at one equals
\(-[z^{-1}]z^k S'(z)\). The vanishing of all these residues makes every
negative Laurent coefficient of \(S(z)\) zero, hence \(S(z)=7\).
This is the local inverse-branch trace calculation for the closed
rational moment problem, as in
[Pakovich, Proposition 3.8](https://arxiv.org/pdf/0910.2105v1).

The monodromy of \(R\) over its regular values is transitive because its
domain after removing finitely many points is connected. We may also
remove the target value zero so that all cube-root lifts can be
continued. Continue the constant relation \(S(z)=7\) along a loop based
at a large regular value that sends one of its selected \(u\)-branches
to the unique near-zero branch. The continued relation is still a sum
of one cube root for each of seven distinct \(u\)-branches.

Now take seven turns around target infinity. In the \(u\)-map, these
turns fix the near-zero branch and each of the seven near-one branches.
The cube-root lifts of the near-one branches return as well, since
those branches lie in a disk about one disjoint from zero. The near-zero
branch winds seven times about zero, up to orientation, because its pole
at target infinity is described by \(u(z)\sim q(0)/z\). Its cube root
is therefore multiplied by a nontrivial cube root of unity. Subtracting
the two continued constant relations forces that cube-root branch to
vanish, contradicting \(u(z)\ne0\). Thus the residue condition and an
eventually vanishing moment tail are both impossible in this case.
