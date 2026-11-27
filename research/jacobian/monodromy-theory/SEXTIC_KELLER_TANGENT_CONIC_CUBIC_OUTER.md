# Cubic outer maps on tangent conic pencils

**Theorem SexticKellerReducedTangentConicCubicExclusion.**
Let \(P=xz+y^2,\ Q=xy\). In characteristic zero there is no
three-variable Keller map of ordinary degree six with three
independent highest forms \(\Phi_i(P,Q)\), where the \(\Phi_i\)
are cubics in two variables with no common projective zero.

**Theorem SexticKellerSquareTangentConicImmersiveCubicExclusion.**
Let \(P=xz+y^2,\ Q=x^2\). In characteristic zero there is no
three-variable Keller map of ordinary degree six with three
independent highest forms \(\Phi_i(P,Q)\) as above if
\([p:q]\mapsto[\Phi_1:\Phi_2:\Phi_3]\) is immersive at every
projective parameter over the algebraic closure.

The first theorem uses the comparison of poles over the two
components of \(Q=0\). The second treats the doubled divisor
\(Q=x^2\) by comparing the two conjugate expansions when the
first perturbation is algebraic over the fixed scalar
parameter field \(k(P)\).

## Descent at the cusp for \(Q=xy\)

Work over an algebraically closed characteristic-zero field \(k\).
The [leading-pencil regularity restriction](SEXTIC_KELLER_INDEPENDENT_LEADING_PENCIL_RESTRICTIONS.md#regularity-restrictions-in-the-tangent-conic-cases)
forces the outer map for \(Q=xy\) to be nonimmersive at
\([1:0]\). Two independent target combinations vanish there
to order at least two and thus span \(pq^2,q^3\).
The remaining form has nonzero \(p^3\) coefficient.
A shear of \(p\) by a multiple of \(q\) removes the \(p^2q\)
term, and the other two components remove the \(pq^2,q^3\)
terms from that form. For \(P=xz+y^2,\ Q=xy\) the parameter
shear is realized by a shear of \(z\) by a multiple of \(y\).
We can therefore use highest tuple
\[
                         \Phi(P,Q)=(P^3,PQ^2,Q^3).
\]

The rational kernel of \(\det J(P,Q,\cdot)\) is \(k(P,Q)\):
at generic fixed \(P,Q\), \(x\) is a free variable and \(y,z\)
can be solved rationally. Also \((P,Q)\) is onto the affine
plane, using \(x\ne0\) for every prescribed pair of values.
Thus the polynomial kernel is \(k[P,Q]\) by polynomial
contraction. In fact if a reduced fraction in \(k(P,Q)\)
with nonconstant denominator pulled back to a polynomial,
one could lift a denominator zero off the numerator and
contradict the cross-multiplied identity.

Use a polynomial reference \(G(\widetilde P,\widetilde Q)\)
with \(\widetilde P,\widetilde Q\) of highest forms \(P,Q\) and
\(G\) having cubic part \(\Phi\).
For a first error \(E_j\) of positive ordinary degree \(j<6\),
the leading determinant equation with one error column is
\(\det J(P,Q,3Q\,n\cdot E_j)=0\), where
\[
                         n=(Q^3,-3P^2Q,2P^3).
\]
Here the cross product of the two parameter derivative
columns of \(\Phi\) is \(3Q n\).
The other terms with two errors have smaller degree while
the displayed determinant degree \(9+j\) is positive.
Thus \(n\cdot E_j\in k[P,Q]\).
The ideal of \(n\) in the parameter ring contains all
homogeneous parts of parameter degree at least four.
For odd \(j\) the invariant vanishes by ordinary degree,
and for \(j=2,4\) it can be absorbed by adding a vector of
polynomials in \(\widetilde P,\widetilde Q\) to \(G\) at
that degree. This leaves a tangent error of the form
\[
          E_j=(3P^2 U,\ Q^2U+2P V,\ 3Q V).
 \tag{1}
\]
Indeed the first entry is divisible by \(P^2\) and the
third by \(Q\), by reducing \(n\cdot E_j=0\) modulo
\(P^2\) and \(Q\) respectively and using coprimeness.
The middle entry then follows from the contraction.
Here \(U\) and \(V\) have degrees \(j-4\) and \(j-2\).
If \(j=1\), both vanish.

We prove that \(V\) has the same restriction \(f(P)\) on
both \(x=0\) and \(y=0\). This will give \(V=f(P)+Q V_1\)
with \(f\in k[P]\), since \(Q=xy\) is reduced and the
restriction \(P=xz\) on \(y=0\) is onto. Homogeneity
permits \(f=0\), a scalar, or a multiple of \(P\).
The tangent part with \(f(P)\) is absorbed into \(G\)
by the vector \((0,2P f(P),3Q f(P))\) at weight \(j\).
The remaining error is absorbed by adding \(U,V_1\)
to \(\widetilde P,\widetilde Q\), since
\(\partial_q\Phi=Q(0,2P,3Q)\).
Thus the reference would eventually agree with the whole
map, after also absorbing constants, and have determinant zero.

## Comparing the first two possible poles

Use ordinary rescaling with \(\epsilon\), so the determinant
has order fifteen. Write \(r=6-j\) for the first error
loss after the normal absorption. We need only \(1\le r\le4\).
Fix the rescaled first component to \(A^3\), where \(A\)
is an indeterminate, and put \(K=k(A)\).
In the reference, solve \(G_{1,\epsilon}(\pi,\theta)=A^3\)
formally for \(\pi\) with leading \(\pi=A\). Its derivative in
\(\pi\) at \(\epsilon=0\) is \(3A^2\).
Define \(b(\theta)=G_{2,\epsilon}(\pi,\theta)/A\) and
\(c(\theta)=G_{3,\epsilon}(\pi,\theta)\). Their leading
terms are \(\theta^2,\theta^3\).
There is a unique critical point \(\xi(\epsilon)\) of \(b\)
with leading term zero, because the initial second derivative
is two. Take the formal coordinate \(s\) with initial
\(s=\theta\) and \(b=b(\xi)+s^2\).
This amounts to taking the square root of the unit in
\((b(\theta)-b(\xi))/(\theta-\xi)^2\).
Use \(c(s,\epsilon)\) to denote the centered reference third
component in this coordinate.

On either divisor \(x=0\) or \(y=0\) at generic \(P=A\),
use \(s\) obtained from \(\widetilde Q_\epsilon\), fix the
actual first component to \(A^3\), and choose a free
variable along the divisor. These are regular formal charts:
at \(x=0\), \(dP,dQ\) are independent because \(y\ne0\);
at \(y=0\) the derivative in \(z\) of \(P\) is \(x\ne0\).
Write the actual functions on the chart as
\[
 F_{2,\epsilon}/A-b(\xi)=s^2+U_*,\qquad
 F_{3,\epsilon}=c(s,\epsilon)+W_*.
\]
Both errors start at order \(r\). Their leading coefficients
satisfy
\[
                   (U_*)_r=2V,\qquad (W_*)_r=3sV,
 \tag{2}
\]
where \(V\) in this identity is evaluated on the leading
level with \(P=A,\ Q=s\). For the first error in (1), fixing
the actual first component shifts the reference \(\pi\) by
\(-\epsilon^r U\) at order \(r\). Its contribution to
\(G_{2,\epsilon}/A\) cancels \(Q^2 U/A\), and its contribution
to the third component is zero at that order.

There is one scalar Laurent series \(R(s,\epsilon)\), with
coefficients in \(K((s))\), such that on both charts
\[
 c+W_*=\sum_{\nu\ge0}\frac{U_*^\nu}{\nu!}\delta^\nu R
       \pmod{\epsilon^{15}},\qquad
 \delta=\frac{1}{2s}\partial_s.
 \tag{3}
\]
To see this, over the generic source field use the first
integrals \(F_{1,\epsilon}^{1/3}\) and
\((F_{2,\epsilon}/A-b(\xi))^{1/2}\), with leading \(P,Q\).
Subtract successive rational functions of these integrals
from \(F_{3,\epsilon}\). At each unmatched order below
fifteen the coefficient is killed by \(\det J(P,Q,\cdot)\),
by the scaled determinant equation. It therefore belongs to
\(k(P,Q)\) and can be lifted to the formal integrals.
The resulting series on the fixed-\(A\) charts gives (3) by
expanding at the new squared coordinate \(s^2+U_*\).
The two charts use the same scalar series, since the generic
normalization was made before restricting to either chart.
By (2) and \(\delta(s^3)=3s/2\) one has
\(R_i=c_i\) for every \(i\le r\), with subscript the coefficient
of \(\epsilon^i\).

Write \(\alpha_i\) for the coefficient of \(s\) in \(c_i\), and
\(u=(U_*)_r|_{s=0}\) on either chart. If \(p<r\) is the first
index with \(\alpha_p\ne0\), the coefficient of
\(s^{-1}\) at order \(p+r\) in (3) is
\[
                     [s^{-1}]R_{p+r}+\alpha_p u/2.
\]
All other terms are regular: higher Taylor terms have order
at least \(2r\) in \(\epsilon\), and the lower linear
coefficients of \(c_i\) vanish. It follows that \(u\) is
the same scalar on both charts.

Otherwise \(\alpha_i=0\) for \(i<r\). The series coefficients
\(R_i\) are regular for \(i<2r\), since the only possible
negative power from a single derivative of a regular series
is from its linear term. At \(2r\), \(R_{2r}\) has at most
a simple pole and its \(s^{-1}\) coefficient \(\rho\) satisfies
\[
                       \rho+\alpha_r u/2+3u^2/8=0.
 \tag{4}
\]
At \(3r\), the terms with pole order three come from
\(\delta R_{2r}\), the linear part of \(\delta^2 R_r\),
and \(\delta^3(s^3)\). The coefficient of \(s^{-3}\) is
\[
 [s^{-3}]R_{3r} - \rho u/2-\alpha_r u^2/8-u^3/16
               =[s^{-3}]R_{3r}+(\alpha_r u^2+u^3)/8.
 \tag{5}
\]
All remaining Taylor terms have at worst a simple pole,
apart from terms of order at most two from multiplying the
listed cubic poles by higher powers of \(s\).
Both orders \(2r,3r\) are below fifteen.

Equation (4) first makes \(u\) algebraic over \(K\) on each
chart. If two distinct values \(u,v\) were possible, comparing
(4) would give \(\alpha_r=-3(u+v)/4\). Comparing (5) would give
\(\alpha_r(u+v)+u^2+uv+v^2=0\). These imply
\((u-v)^2/4=0\), a contradiction. On \(y=0\) the field of
the chart at \(s=\epsilon=0\) is a rational function field in
one variable over \(K\), so the common value is in \(K\).
By (2), \(V\) therefore has the same restriction \(f(P)\)
on both divisors. This completes the reference descent.

## Immersive cubics and the square member \(Q=x^2\)

Write \(t=x\) in this case. The polynomial constant ring of
\(\det J(P,Q,\cdot)\) is \(T=k[P,t]\): the generic curve at
fixed \(P,t\) is rational by solving for \(z\), and \((P,t)\)
is onto, including at \(t=0\) where \(P=y^2\).
Let \(S=k[P,Q]\), with \(T=S\oplus tS\), and let
\(n=\partial_p\Phi\times\partial_q\Phi\).
The three binary quartics \(n_i\) have no common factor
because the projective parametrization is immersive and the
radial derivative is nonzero. In the source ring the
substituted minors still have gcd one, since \(P,Q\) are
coprime and no projective point kills all the minors.

If \(n(P,Q)\cdot E=0\), then \(E\) is a polynomial tangent
vector
\[
                     E=\partial_p\Phi(P,Q)\,U
                       +\partial_q\Phi(P,Q)\,V .
 \tag{6}
\]
Over the fraction field this follows by rank two.
At each prime divisor, some minor is a unit and gives
regularity of the two coefficients. The polynomial ring is
the intersection of its height-one localizations, so the
coefficients are polynomials. For an error of degree \(j\)
they are homogeneous of degree \(j-4\).
The same syzygy calculation in \(k[p,q]\) shows why a
contraction with \(n\) can realize all binary forms of
parameter degree at least five. For coefficient degree
\(h\ge1\), the kernel has dimension \(2(h-1)\) from the
two tangent columns of parameter degree two. The image has
dimension \(3(h+1)-2(h-1)=h+5\), the dimension of the binary
forms of degree \(h+4\). For \(h=0\) the kernel is zero,
so the three \(n_i\) are linearly independent over \(k\).

As before a first error \(E_j\) for
\(G(\widetilde P,\widetilde Q)\) satisfies \(n\cdot E_j\in T\)
by the leading determinant equation. If this lies in \(S\),
it can be absorbed by lower terms in \(G\), since the
ordinary degree \(8+j\) is even and its parameter degree
is at least five. The tangent error is then removed by (6).
Thus any Keller completion must reach a first non-\(S\)
normal error, of loss \(r=6-j\). Its invariant lies in \(tS\),
so \(r=1,3,5\).

Choose two highest components \(\Phi_1,\Phi_2\) whose
parameter Jacobian \(n_3\) is nonzero generically on \(q=0\).
This is possible by immersion at that projective parameter.
Frame the actual rescaled first two components by solving
\[
       \Phi_1(A,B)=F_{1,\epsilon},\qquad
       \Phi_2(A,B)=F_{2,\epsilon},\qquad
       (A,B)|_{\epsilon=0}=(P,Q).
\]
This formal construction in the generic field is regular
along \(t=0\), as \(n_3(P,0)\ne0\) generically.
Fix \(A\) as a parameter on a formal chart, take \(z\)
free, and solve the equation for \(A\) for \(y\), starting
with \(y=y_0,\ y_0^2=A\) at \(t=\epsilon=0\).
In this chart \(B=t^2+O(\epsilon)\).

Subtract the reference third component expressed in \(A,B\):
solve \(G_{\epsilon,i}(\pi,\theta)=\Phi_i(A,B)\), \(i=1,2\),
with initial \((\pi,\theta)=(A,B)\) and use
\(G_{\epsilon,3}(\pi,\theta)\). This is regular at \(B=0\)
over \(K=k(A)\).
The difference \(C\) is regular on the \(t\)-chart and
begins at order \(r\) with coefficient
\((n\cdot E_j)/n_3(P,Q)\). As in (3), the determinant gives
\[
             C=R(B^{1/2},\epsilon)\pmod{\epsilon^{15}},
             \qquad R_i(s)\in K((s)),\quad R_i=0\ (i<r),
 \tag{7}
\]
where the formal root has leading \(t\). The generic
constant field is \(k(P,t)\) and constrains the coefficients
before we pass to the chart containing \(y_0\).
In \(R_r\) there is a nonzero odd power of \(s\) of exponent
\(l\le14-r\): the numerator is an element of \(tS\) of
ordinary degree \(8+j=14-r\), and the denominator is an
even unit at \(t=0\).

We claim that (7) is impossible if \(v=B_1|_{t=0}\) is
transcendental over \(K\), or if \(v^2\in K^\times\) and
\(v\notin K\). Give \(t,\epsilon\) weights one and two.
For \(v\ne0\) the minimum-weight part of \(B\) is \(t^2+\epsilon v\).
There is a least weight \(W\le14+r\) of an offending term
\(\epsilon^i s^l\) with \(l\notin 2\mathbb Z_{\ge0}\).
Remove all lower-weight harmless terms and those at weight
\(W\), by subtracting their evaluations as powers of \(B\).
If \(W\) is even, the first remaining \(\epsilon\)-index
at weight \(W\) already gives a negative \(t\)-power.
So \(W\) must be odd. At the smallest index \(n\ge i\)
with \(2n>W\) the coefficient at weight \(W\) from the
chosen term is its scalar coefficient times
\(\binom{l/2}{n-i}v^{n-i}\). This is nonzero, \(n<15\),
and all other terms at the same index carry distinct
powers of \(v\). This excludes transcendental \(v\).

If \(v\) is quadratic, put \(D=(W-1)/2\).
The weight-\(W\) part after the forward substitution,
by regularity, is a polynomial in \(\epsilon\) with
terms \(\epsilon^i t^{W-2i}\) for \(r\le i\le D\).
Its coefficients belong to \(K(v)\), by the forward
binomial formula involving only \(v\) and scalars.
Invert the edge substitution with \(t^2=s^2-\epsilon v\)
and put \(Z=\epsilon/s^2\). After factoring out \(s^W\),
the weight-\(W\) part of the normalized series is
\[
             (1-vZ)^{1/2} A_v(Z)\pmod{Z^{15}},
             \qquad \deg A_v\le D,\quad Z^r\mid A_v.
\]
Indeed a term with index \(i\le D\) gives
\(Z^i(1-vZ)^{D-i+1/2}\).
The normalized coefficients are over \(K\), so the
conjugate expression with \(-v\) agrees modulo \(Z^{15}\).
Squaring and using the divisibility by \(Z^r\) gives
\[
               (1-vZ)A_v^2-(1+vZ)A_{-v}^2
                        =0\pmod{Z^{15+r}}.
\]
The polynomial on the left has degree at most
\(2D+1=W\le14+r\), so it vanishes identically. The
multiplicities at \(Z=1/v\) have different parity unless
both polynomials are zero. That contradicts the choice
of the offending term, by invertibility of the formal
edge substitution. This finishes the claim.

It remains to relate \(v\) to (6) and finish the tail of the
descent. Even if the first degree-five invariant is in \(tS\),
we can absorb it for purposes of that decomposition by a
vector \(t M(P,Q)\), with \(M\) quadratic in the parameters.
This uses the binary degree-six image calculation above.
The remaining degree-five tangent vector has parameters
\(U,V\) linear. On \(t=0\) the normal vector vanishes and
the first-order change in \(B\) is therefore \(V|_{t=0}\).
Fixing \(A\) causes no extra term, since \(B_0=t^2\).
If \(V|_{t=0}\) involves \(z\) it is transcendental over
\(K\); otherwise any nonzero value is a multiple of \(y_0\)
and has degree two over \(K\). Both are impossible by (7).
Thus \(V\in k[t]\).

Now use references
\(\Phi(\widetilde P,\widetilde Q)+L(\widetilde P,t)\) with
\(\widetilde Q\in k[t]\) and \(L\) of lower weight.
The degree-five correction is absorbed as just described.
At degree four the tangent parameters are constant,
so the second parameter stays in \(k[t]\). At degrees
three and two the odd or even invariants are absorbed
in \(T\) by the binary image calculation and there is
no tangent part. At degree one write
\(E_1=a x+b y+c z\) with constant vectors \(a,b,c\).
If \(n\cdot E_1\in k(P,t)\), with \(t=x\) one can
substitute \(z=(P-y^2)/t\) and vary \(y\) freely.
It follows that \(n\cdot b=n\cdot c=0\).
Linear independence of the \(n_i\) gives \(b=c=0\);
the remaining \(a t\) is absorbed directly. Constants
can be included in the reference as well, yielding
a two-parameter factorization and a zero determinant.
