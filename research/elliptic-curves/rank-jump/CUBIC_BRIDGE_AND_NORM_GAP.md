# The cubic-field bridge and its norm-solubility gap

The frozen three-class test from the [CT block report](CT_VARIATION_AND_BLOCKS.md)
has a negative explanatory result: **the first shared norm conics are soluble
automatically throughout the inherited class space**. They do not explain
the simultaneous CT block or distinguish rational solubility. There is,
however, an exact way to simplify the Galois action without losing the
rational CT information: use the degree-three cubic field, not the full
degree-six splitting field.

This note proves the field-change and cover identities needed to take that
route correctly. The [protocol](CUBIC_BRIDGE_PROTOCOL.json) fixes the three
chain masks and seven existing parameters. The
[certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_cubic_bridge_v1.json)
contains 23 explicit norm witnesses, 161 transported norm-conic checks and
140 coordinate-compatibility checks. Replay:

    python3 elliptic-curves/rank-jump/cubic_bridge.py check

No new conic solving, descent, point search or parameter enumeration was
needed. All anchor point use is retrospective.

## An odd-degree field preserves the obstruction

Let \(K=\mathbf Q(\theta)\), \(f(\theta)=0\), for the irreducible \(S_3\)
cubic \(f(T)=T^3+AT+B\) in the preceding reports. Write \(L\) for its
degree-six splitting field.

For an elliptic curve defined over a number field \(F\) and a finite
extension \(F'/F\), restriction and corestriction preserve the elliptic
Kummer local conditions. They satisfy

\[
\langle\operatorname{res}c,\operatorname{res}d\rangle_{F'}
 =[F':F]\langle c,d\rangle_F.
\]

See [Morgan–Smith, Corollary 4.9(3) and Remark 4.10](https://link.springer.com/article/10.1007/s40993-024-00545-2).
For classes of order dividing two, degree three preserves the value while
degree six makes it zero. Thus the certified rational CT block survives
unchanged on the restricted classes over \(K\). Over \(L\), all pairings
between those restricted rational classes vanish.

Vanishing after the even-degree extension does **not** prove that those
classes acquire points, nor that they pair trivially with every class over
\(L\). It means that their mutual pairings alone have lost the information we
want. A full-splitting-field governing construction would need extra descent
or corestriction data to recover it.

There is a stronger odd-degree statement for the *same rational genus-one
torsor*. If a torsor of period dividing two acquires a \(K\)-point, then its
class \(c\in H^1(\mathbf Q,E)[2]\) has
\[
0=\operatorname{cor}\operatorname{res}c=3c=c.
\]
It already has a rational point. Hence testing the full torsor over \(K\)
does not weaken rational solubility. Replacing it by its isogeny quotient
cover is a different operation and can lose information.

## A rational 2-isogeny over the cubic field

Put
\[
\alpha=\theta+u\theta^2,\qquad \gamma=1-u\theta,\qquad
D(u)=1+Au^2+Bu^3.
\]

Translate the rational curve \(E_u:y^2=F_u(x)\) by \(X=x-\alpha\) over
\(K\). It becomes
\[
y^2=X(X^2+a_uX+b_u),
\]
where
\[
\begin{aligned}
a_u&=3\theta+u(3\theta^2+2A),\\
b_u&=(3\theta^2+A)\{1+u\theta+u^2(\theta^2+A)\}\\
   &=(3\theta^2+A)\frac{D(u)}{1-u\theta}.
\end{aligned}
\]

Its rational \(K\)-point of order two is \((0,0)\). The quotient by this
point is
\[
E'_u:\quad y'^2=X'\bigl(X'^2-2a_uX'+a_u^2-4b_u\bigr),
\]
with the usual degree-two map, away from \(X=0\),
\[
X'=X+a_u+b_u/X,\qquad y'=y(1-b_u/X^2).
\]
Substitution directly verifies the curve identity.

Set
\[
\delta=-3\theta^2-4A.
\]
The key factorizations are
\[
a_u^2-4b_u=\delta\gamma^2,\qquad
\delta(3\theta^2+A)^2=\operatorname{disc}f,\qquad
N_{K/\mathbf Q}(b_u)=-\operatorname{disc}f\,D(u)^2.
\]
Consequently the quadratic field splitting the *remaining* two-torsion
points of \(E_u/K\) is the constant field
\(K(\sqrt\delta)=L\). The isogenous curve's other two-torsion points instead
depend on the squareclass of \(b_u\). This locates a varying layer beyond
the fixed original \(E[2]\).

These identities are exact in \(K(u)\). The replay also checks every
coefficient for the seven retained specializations, without constructing a
new number field or computing its class group.

## Why the first norm conics always pass

For an anchor point \(P=(p,q)\), with \(q^2=f(p)\), put
\(\beta=p-\theta\). Then
\[
s=\frac{\beta(2p+\theta)}{2q},\qquad
t=-\frac{\beta}{2q}
\]
satisfy
\[
s^2-\delta t^2=\beta.
\]

Indeed,
\[
(2p+\theta)^2-\delta=4(p^2+p\theta+\theta^2+A)
 =\frac{4q^2}{p-\theta}.
\]
Thus \(s+t\sqrt\delta\) is an explicit norm witness for \(\beta\).
Multiplying such witnesses gives a witness for every product of the twenty
anchor classes. Changing a representative by a square does not affect
existence. This proves the assertion on the whole \(2^{20}\)-element
inherited space without enumerating it.

At parameter \(u\), use \(s+(t/\gamma)\sqrt{\delta\gamma^2}\).
Its norm is still \(\beta\). When \(t\ne0\), it also supplies the point
\[
\mathcal X=-\gamma s/t,\qquad \mathcal Y=\gamma/t
\]
on the conic
\[
\mathcal X^2-\beta\mathcal Y^2=\delta\gamma^2.
\]
All 23 retained witnesses have the required nonzero denominators.

In particular the three previously fixed masks
\[
e_3=317529,\quad f_0=491700,\quad e_4=631775
\]
have soluble norm conics at both \(u=-1\) and \(u=1\), although their retained
pairings are
\[
\begin{array}{c|cc}
u&\langle e_3,f_0\rangle&\langle e_4,f_0\rangle\\ \hline
-1&1&0\\
1&1&1
\end{array}
\]
in \(\mathbf F_2\). More strongly, every one of the three classes is
certifiably obstructed on both curves: at \(u=-1\) each pairs to 1 with
anchor mask \(1\); at \(u=1\), witnesses for \(e_3,f_0,e_4\) are respectively
masks \(163,260,260\). These are rational CT certificates, not bounded
point-search misses.

The **shared quadratic norm field and existence of these first norm
witnesses therefore carry no additional solubility information on \(W\)**.
They are an already satisfied class-lifting prerequisite, not a prospective
rank feature. This does not exclude *pair-specific secondary* norm data.

## The constraints omitted by a norm-only explanation

For a fixed \(\beta\in K^\times\), the dual-isogeny covering of \(E_u/K\) has
the quartic model
\[
\mathcal D_{\beta,u}:\quad v^2=\beta z^4+a_uz^2+b_u/\beta,
\]
mapping to the translated curve by \(X=\beta z^2,\ y=\beta zv\).
Completing the square gives
\[
(2\beta z^2+a_u)^2-\beta(2v)^2=\delta\gamma^2.
\]

Thus a point on the norm conic still has to satisfy
\[
\frac{\mathcal X-a_u}{2\beta}=z^2\quad\text{in }K.
\]
That is a new square constraint. The conic's rational parametrization does
not solve it.

Even a point on this isogeny quartic does not automatically give a rational
point in the original two-cover class. In the affine chart, a sufficient
and exact coordinate compatibility condition is
\[
x=\alpha+\beta z^2\in\mathbf Q.
\]
Its two nonconstant power-basis coefficients must vanish. Then
\(F_u(x)=N(\beta)N(z)^2\) has a rational square root, and
\(x-\alpha=\beta z^2\) is precisely the original rational Kummer identity.
Conversely, a rational point in this class gives such a \(z\) and hence the
isogeny-quartic point, away from the usual projective charts. A \(K\)-valued
square root of a rational number is rational because a cubic field contains
no quadratic subfield.

For the explicit transported conic point above and a single anchor basis
class, \(\mathcal X=\gamma(2p+\theta)\). Its implied value of \(z^2\) is
\[
1-u\frac{p\theta+2\theta^2+A}{p-\theta}.
\]
Even **if** this happens to be a square in \(K\), the mapped coordinate is
\[
x=p-u(p\theta+\theta^2+A).
\]
Its \(\theta^2\) coefficient is \(-u\). It is rational only at \(u=0\).
This is a universal algebraic failure of this particular transported norm
point; no extra square tests or searches are needed. It does not rule out
other points on the same conic, quartic or two-cover.

## What this changes next

The previous chain experiment has not obtained its desired CT variation
formula. It has eliminated the first shared norm-conic layer as the missing
explanation and exposed two ways of accidentally weakening the problem:
even-degree base change erases the mutual restricted pairing, and replacing
the full cover by a norm conic omits square and rationality conditions.

The useful remaining target is **secondary descent over the odd-degree
cubic field**, or a rational auxiliary construction imposing the square and
coordinate conditions for several independent classes together. Its
arithmetic inputs must go beyond the constant norm field and the existence
of its conic points. A formula explaining the chain operator is still
UNKNOWN.

Before investing in that more expensive construction, the next cheap check
should address actual high-gain fibres: determine whether their certified
quotient independence already forces maximal independent half-point
extensions over the common two-division field. That separates a genuine
shared construction of rational points from a misleading assertion that
their halving fields collapse. It can use the pinned rank-jump panel and
finite \(S_3\)-module algebra, with no additional point searches.
