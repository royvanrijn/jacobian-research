# Exact sphere moments beyond every single angular shear

Polynomial changes of sphere coordinates can preserve every moment
while making one-sided Laurent support impossible in every linear
angular chart after a single shear. The family below has divergent
bounded-phase volume and explicit prime-power moment poles. Its
moments come from a polynomial of one variable after a composition
of shears. The obstruction to using a single shear is the pole order
of every linear coordinate on the two ends of a generic fiber.

## The phases and their moments

Let \(Q\) be the affine sphere \(2xy+t^2=1\) over a characteristic-zero
field. Write \(\Lambda\) for its normalized functional in the Laurent
chart, taking angular constant term in \(x\) and then the interval
integral in \(t\). On points of \(Q\) define the two automorphisms
\[
 \begin{split}
 S_x(x,y,t)&=(x,\,y-xt-x^3/2,\,t+x^2),\\
 S_y(x,y,t)&=(x-yt-y^3/2,\,y,\,t+y^2).
 \end{split}
\]
They are the time-one maps of the standard sphere shear vector
fields multiplied by their invariant coordinates. Their inverses
replace the shifts \(x^2,y^2\) by \(-x^2,-y^2\), respectively.
[Complete algebraic vector fields on Danielewski surfaces, Main Theorem](https://www.numdam.org/article/AIF_2016__66_2_433_0.pdf)

Put \(S=S_x\circ S_y\circ S_x\) and \(F=t\circ S\). For every
nonconstant polynomial \(R\) let \(f=R(F)\). Then
\[
                  \Lambda(f^M)=\frac12\int_{-1}^{1}R(t)^M\,dt
                       \qquad(M\ge0).
\]
Indeed each shear preserves the sphere and its volume form,
which is \(dx/x\wedge dt=-dy/y\wedge dt\) on the chart overlap.
One can scale its shift
continuously to zero, giving an isotopy over \(\mathbb C\). Thus
the shear takes the compact real sphere cycle to a homologous
cycle, and integration of a regular top form times a polynomial
is unchanged. The moment identity follows for arbitrary
characteristic-zero coefficients by embedding their finitely
generated field into \(\mathbb C\).
[Sphere residue formula](GLOBAL_SPHERE_PURITY_RECOGNITION.md#1-a-closed-rational-period)

In particular take \(R(T)=T^h+T\), \(h\ge2\). For every odd
prime \(p\nmid h\) and each \(r\ge1\) with \(p^r\equiv1\pmod h\),
put \(M=(p^r-1)/h\). In the univariate integral the leading
height is \(p^r-1\), even, with coefficient one. Every other
denominator has lower \(p\)-valuation. Hence
\[
                         p^r\Lambda(f^M)\equiv1\pmod p .
\]
This gives infinitely many exact prime-power poles at each such
prime. The phase has least sphere degree \(18h\). In the forward
word, the first transformed \(y\) has degree three and leading
term \(-x^3/2\); the next transformed \(x\) has degree nine with
leading term \(x^9/16\); finally \(F\) has highest term
\(x^{18}/256\). This term cannot be removed using the sphere
relation, since the resulting phase grows with that degree along
\(t=0,\ y=1/(2x),\ x\to\infty\). Each \(R(F)=F^h+F\) is an
absolute nonpower because \(R\) has a simple zero at zero and
\(F=0\) is a reduced irreducible fiber, carried by \(S\) to \(t=0\).

## Pole orders at the two ends of a fiber

Over \(\mathbb C\), fix \(\lambda\ne\pm1\). The fiber \(F=\lambda\)
is parametrized by applying \(S^{-1}\) to
\[
                   (X,Y,T)=(z,(1-\lambda^2)/(2z),\lambda),
                   \qquad z\in\mathbb C^\times.
\]
Its two ends are \(z=0\) and \(z=\infty\). The pole orders of the
original coordinates at the two ends are
\[
\begin{array}{c|rrr}
   & x&y&t\\ \hline
 z=\infty&9&27&18\\
 z=0&3&9&6\,.
\end{array}
\]
Here \(S^{-1}=S_x^{-1}\circ S_y^{-1}\circ S_x^{-1}\).
At infinity, after the first inverse step the orders are \((1,3,2)\).
In the second inverse step the cubic in the old \(y\) dominates
the new \(x\), and the square in \(y\) dominates the new \(t\),
giving \((9,3,6)\). The last inverse step gives \((9,27,18)\)
by the same strict comparisons. At zero, after the first step
\(x=z\) vanishes, \(y\) has a simple pole and \(t\) is bounded.
The remaining steps give \((3,1,2)\) and \((3,9,6)\).
The orders at each end are distinct and positive. Therefore every
nonzero linear form in \(x,y,t\) has pole order at least three
at both ends, since among its terms the largest order is unique.

## Neither angular sign can be removed

**Theorem.** For every nonconstant \(R\), the phase \(R(F)\) has
both positive and negative literal Laurent angles in every linear
sphere chart, after any single polynomial shear preserving that
chart's angular coordinate. This includes shears preceded by an
arbitrary linear orthogonal change over the algebraic closure.

Use an arbitrary such chart \((L,Y_0,T_0)\) with
\(2LY_0+T_0^2=1\), so \(L\) is a nonzero isotropic linear
coordinate in the original variables. A polynomial shear replaces
\(T_0\) by \(v=T_0-P(L)\) with \(P\in L\,\mathbb C[L]\),
adjusting \(Y_0\) to preserve the sphere equation. The moments
are unchanged under the linear chart change as well as the shear.
For the linear assertion, in ordinary real sphere coordinates
the even moment of a linear form of coefficient vector \(a\)
is \((a\cdot a)^k/(2k+1)\) in degree \(2k\), and odd moments
are zero. This follows by real orthogonal invariance and integration
of one coordinate. The polynomial formulas extend to complex
coefficient vectors. Polarization gives invariance of all homogeneous
moments under the complex orthogonal group.

Suppose that the transformed expansion is one-sided, and let \(R_0(v)\)
be its coefficient at angle zero. Its moments are the interval
moments of \(R_0\). The polynomial \(R_0\) cannot be constant. If
it were the constant \(c_0\), every positive moment of \(f-c_0\)
would vanish. By the moment identity these are the interval moments
of the nonconstant polynomial \(R(t)-c_0\). A nonconstant univariate
polynomial cannot have all interval moments zero: after one
number-field specialization preserving its degree, choose a good
odd prime \(p\) prime to the degree, and take the moment order
\((p^r-1)/\deg(R)\) with \(p^r\equiv1\bmod\deg(R)\). The unique
interval denominator divisible by \(p^r\) has nonzero leading
coefficient. This is also the univariate step in the
[one-sided shear criterion](SPHERICAL_SUPERCRITICAL_FACE_SHEAR.md#a-purity-criterion-after-an-angular-shear).

Choose a complex value \(s\) outside the critical values of \(R_0\),
outside \(R_0(1),R_0(-1)\), and outside \(R(1),R(-1)\). If the
angular support is nonpositive, set \(u=1/L\). If it is nonnegative,
set \(u=L\). In the selected case the Laurent expression extends
holomorphically to \(u=0\) as a polynomial in \(u,v\), with value
\(R_0(v)\) there. Choose a root \(v_0\) of \(R_0(v)=s\). It is
simple and \(v_0\ne\pm1\). The implicit function theorem gives
a holomorphic \(v=v(u)\) on the level \(f=s\), with \(v(0)=v_0\).
For \(u\ne0\) the three original coordinates along this branch
are meromorphic in \(u\). The nonpositive case has the pole
\(L=u^{-1}\). The nonnegative case has \(L=u\) and the companion
\((1-T_0^2)/(2L)\) has a pole, since \(T_0\to v_0\). In either
case the branch leaves the affine sphere at \(u=0\).

The value \(F\) on this punctured branch is constant, equal to some
\(\lambda\) with \(R(\lambda)=s\). In particular \(\lambda\ne\pm1\).
The fiber parameter \(z=X\) on the output under \(S\) is meromorphic
in \(u\) and has a zero or pole at \(u=0\). If it had a finite
nonzero limit, the inverse parametrization would give a finite
limit in the original sphere. Thus the order in \(u\) of every
linear coordinate is the order along the appropriate fiber end
multiplied by \(|\operatorname{ord}_u z|\). In the nonpositive
case this makes a simple pole of \(L\) impossible. In the
nonnegative case \(L\) would have a pole, contrary to \(L=u\).
This excludes both one-sided possibilities.

## Comparison with finite-volume obstructions

For the radial phase \(R(t)\), take \(t\) in any bounded complex
disk and let \(|x|\) tend to infinity. The phase stays bounded
and the squared volume density in the Laurent chart is a constant
times \(|x|^{-2}d^2x\,d^2t\), with divergent integral. Since
\(S\) preserves the volume form, its inverse image gives the same
divergence for \(R(F)\). Thus the bounded-phase integrability step
of the [weighted proper-volume criterion](SPHERICAL_WEIGHTED_FACE_PROPER_VOLUME_OBSTRUCTION.md)
cannot establish these moment identities. The univariate moment
formula uses a composition of shears, although no single angular
shear in any linear chart yields one-sided support.

The [exact replay](verify_spherical_bilateral_automorphic_moments.py)
computes the generic fiber orders by symbolic Laurent substitution
and checks the degree and first two moments of \(F\) directly in
the original Laurent chart.
