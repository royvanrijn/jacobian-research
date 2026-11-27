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

## Normal forms for the nonimmersive square-member cubics

**Proposition SexticKellerSquareTangentConicCuspidalNormalForms.**
Let \(P=xz+y^2\), \(Q=x^2\) over an algebraically closed field
of characteristic zero. Suppose three independent binary cubics
\(\Phi_i(p,q)\) have no common projective zero and their
projective map is nonimmersive at a point. That point is unique.
Up to linear source changes preserving the pencil and an invertible
constant target change, \((\Phi_i(P,Q))_i\) is one of
\[
\begin{array}{ll}
 (P^3,\ P Q^2,\ Q^3) & \text{at the square member }Q=0,\\
 (Q^3,\ Q P^2,\ P^3),\
 (Q^3+Q^2P,\ Q P^2,\ P^3) & \text{away from }Q=0 .
\end{array}
\]

Choose binary parameters \(u,v\) so that the nonimmersive point
is \(v=0\). In the three-dimensional space spanned by the cubics,
the two-dimensional subspace of forms vanishing at that point
vanishes to order at least two: otherwise the affine target
ratios would have a nonzero derivative there. This subspace is
exactly \(\langle uv^2,v^3\rangle\). Modulo it, the remaining
form can be normalized to \(u^3+\lambda u^2v\).
Set \(U=u+\lambda v/3\). The cubic space then has basis
\(U^3,Uv^2,v^3\). In \(s=v/U\) the map is
\([1:s^2:s^3]\), immersive for all finite \(s\ne0\),
and at infinity it is immersive by the coordinate \(U/v\).
This proves uniqueness.

If \(v\) is proportional to \(q\), the change that replaces
\(P\) by \(P+\gamma Q\) is the linear shear \(z\mapsto z+\gamma x\).
This gives the first normal form. Otherwise normalize \(v\)
to \(p+\gamma q\); the same shear realizes this choice as the
new \(P\). Then \(U\) can be taken as \(q+\mu p\) after scaling,
and constant target combinations put the leading tuple in the
form \((Q^3+3\mu Q^2P,QP^2,P^3)\).
If \(\mu\ne0\), rescale \(P\) independently of \(Q\) by sending
\((x,y,z)\) to \((x,a y,a^2 z)\) for a suitable nonzero \(a\);
constant target scalings give the stated coefficient one.
The case \(\mu=0\) gives the other tuple.

## Degree-five constraints when the cusp is away from the square

**Theorem SexticKellerSquareTangentConicAwayCuspFirstCorrection.**
Let \(P=xz+y^2\), \(Q=x^2\), and \(H=Q+\mu P\) over an
algebraically closed characteristic-zero field \(k\). If a
three-variable Keller map has highest tuple
\[
                   (H^3,HP^2,P^3),
\]
its homogeneous degree-five vector \(E\) necessarily has the form
\[
 E=\bigl(3H^2U,\ P^2U+2H(PL+\gamma x^3),\
                  3P(PL+\gamma x^3)\bigr),
 \qquad U-\mu L=b x+a y,\qquad a\ne0,
\]
where \(U,L\) are linear forms and \(a,b,\gamma\in k\).
This applies to both away-from-the-square normal forms.

### The quadratic parameter on the square divisor

Write \(t=x\), and use ordinary rescaling with \(\epsilon\), so
the determinant is \(\kappa\epsilon^{15}\). In the generic field
put
\[
 A=F_{3,\epsilon}^{1/3},\qquad B=F_{2,\epsilon}/A^2-\mu A,
 \qquad s'=B^{1/2},
\]
with leading terms \(P,Q,t\), respectively. In coordinates
\(A,s',y\) this is an invertible formal change at a generic point:
its leading map is the rational coordinate system \(P,t,y\).
Since \(F_{2,\epsilon}=A^2((s')^2+\mu A)\) and
\(F_{3,\epsilon}=A^3\), the determinant shows that the first
component in these coordinates has coefficients in \(k(A,s')\)
at every order \(i<15\). Denote them by \(R_i(s')\). In
particular \(R_0=((s')^2+\mu A)^3\) is even and regular at zero.

Near \(x=0\), \(A,t,z\) are regular formal coordinates over
\(K_1=k(A,y_0)\), where \(y_0^2=A\). Indeed at order zero the
equation for \(y\) is \(y^2=A-tz\), and its derivative at
\(t=0,y=y_0\) is a unit. Let \(K=k(A)\) and
\(v=[\epsilon]B|_{t=0}\). On this chart
\[
              v=\frac{(E_2-\mu E_3)(0,y_0,z)}{A^2}.
 \tag{8}
\]
It is polynomial in the free variable \(z\) over \(K_1\).
If it has \(z\)-degree at most one, it is \(a y_0+b_1 z\)
for constants \(a,b_1\in k\), by ordinary homogeneity.

Use also order fifteen in the same normalization. Now take
\(A,s',z\) as coordinates after inverting \(s'\) formally.
The coefficients \(R_i\) for \(i<15\) are unchanged, because
they were independent of the free variable \(y\).
The Jacobian of the two fixed components with respect to
\(A,s'\) is, up to sign, \(6A^4s'\), and the leading Jacobian
of the ambient coordinates \(A,s',z\) is, up to sign,
\(2\sqrt{A-s'z}\). Thus the \(z\)-derivative of \(R_{15}(s',z)\)
at fixed \(A,s'\) has a simple pole whose coefficient is
\(\kappa_0/(A^4 y_0)\) with \(\kappa_0\in k^\times\).
Consequently \(R_{15}\) has a scalar Laurent part over \(K_1\)
and a nonscalar principal part consisting of
\[
                       \frac{\kappa_0 z}{A^4 y_0 s'}.
 \tag{9}
\]
The remaining nonscalar terms are regular in \(s'\).

Suppose \(v\) depends on \(z\), so it is transcendental over
\(K_1\). Give \((\epsilon,t)\) weights \((2,1)\). The
minimum part of \(B\) is \(t^2+\epsilon v\). Call a scalar
term \(\epsilon^i (s')^\ell\) offending when \(\ell\) is
negative or odd. Include the pole (9), of weight \(29\),
in the weight comparison. Harmless terms of weight at most
the minimum offending weight can be removed by subtracting
their evaluations as nonnegative integer powers of \(B\).
These are regular on the original chart.

If the minimum weight is less than \(29\), it comes from
scalar terms alone. At that weight a term with \(\ell<0\)
already contributes a pole at order \(i\); a term with
\(\ell\) positive odd contributes one at order
\(i+(\ell+1)/2<15\). In either comparison the contributions
at the same order and weight have distinct powers of \(v\)
and nonzero binomial coefficients for the selected term.
The residue cannot vanish for transcendental \(v\).
Thus the minimum weight is \(29\). At that weight, the
coefficient of \(\epsilon^{15}t^{-1}\) is a scalar polynomial
in \(v\), plus \(\kappa_0 z/(A^4 y_0)\) from (9). The scalar
coefficient at \(i=15,\ell=-1\) supplies the constant
term of the polynomial. The terms at indices \(i<15\)
have distinct powers \(v^{15-i}\), with the linear
coefficient coming only from \(i=14,\ell=1\).
Regularity would make that polynomial in \(v(z)\) affine
in \(z\) of nonzero slope. Both the polynomial and \(v(z)\)
would therefore be linear. But by (8) the slope \(b_1\)
of \(v\) is a constant in \(k\), and the coefficient
at \(i=14,\ell=1\) is in \(K\). Their product times
\(\binom{1/2}{1}\) cannot cancel \(\kappa_0/(A^4 y_0)\),
which is outside \(K\). We conclude that \(v=a y_0\).

### The first normal correction

Put \(n=(P^3,-3H^2P,2H^3)\). These are the normal
coefficients of the cusp tuple after its derivative
minors' common factor is removed. The highest linear
determinant equation puts \(n\cdot E\) in \(k[P,t]\).
Indeed the rational kernel of \(\det J(P,Q,\cdot)\)
is \(k(P,t)\), using \(t\ne0\) to solve for \(z\).
The pair \((P,t)\) is onto, also at \(t=0\) where
\(P=y^2\), so the polynomial kernel is \(k[P,t]\).
The invariant \(n\cdot E\) has degree eleven, hence
it is \(t\) times a form in \(P,Q\). The order-one
normalized coefficient of the first component is
\[
                       R_1(s')=\frac{n\cdot E}{P^3}
\]
evaluated at \(P=A,\ t=s'\). This follows from the
linearization of \(F_3=A^3,\ F_2=(B+\mu A)A^2\) and
subtraction of the first change of \(R_0\).
If \(n\cdot E\ne0\), \(R_1\) is a nonzero odd
polynomial of \(s'\)-degree at most eleven.

If also \(v\ne0\), then \(v^2=a^2 A\in K^\times\) and
\(v\notin K\). The preceding weight comparison now has
an offending term of weight at most thirteen. Use
coefficients only through order fourteen and subtract
the harmless terms up to minimum offending weight \(W\).
An even \(W\) would give a negative power already at
the smallest index there, so \(W=2d+1\) is odd.
At weight \(W\), regularity after forward substitution
gives a sum of terms \(\epsilon^i t^{W-2i}\) for
\(0\le i\le d\) with coefficients in \(K(v)\).
Invert that substitution with
\(t^2=(s')^2-\epsilon v\). On dividing by \((s')^W\)
and putting \(Z=\epsilon/(s')^2\), its evaluation is
\((1-vZ)^{1/2}L_v(Z)\) modulo \(Z^{15}\), where
\(\deg L_v\le d\). The normalized coefficients lie in
\(K\), so the conjugate expression with \(-v\) agrees.
After squaring the two polynomial expressions have
degree at most \(W\le13\), hence are equal as polynomials.
At \(Z=1/v\) they have different multiplicity parity
unless both vanish identically. This would remove the
entire chosen minimum weight, a contradiction.
Therefore \(v\ne0\) implies \(n\cdot E=0\).

### Descent when the quadratic perturbation vanishes

We prove that \(v=0\) precludes a Keller completion.
Use references \(G(\widetilde H,\widetilde P,t)\), starting
with the displayed cusp tuple and keeping
\(\widetilde H-\mu\widetilde P\in k[t]\) with leading
\(t^2\). Such a reference depends on two source
polynomials, \(\widetilde P,t\). Suppose the first error
has degree \(j\) with \(2\le j\le5\), and write \(r=6-j\).
The same linear determinant equation puts its normal
contraction in \(k[P,t]\). The ideal of the three normal
forms in the parameter ring \(k[H,P]\) contains every
form of parameter degree at least four. In degree
\(6+j\) the even part of \(k[P,t]\) lies in \(k[H,P]\)
and the odd part in \(t k[H,P]\). For the displayed
\(j\), in either parity the required form in \(H,P\)
has parameter degree at least four. Absorb it by
evaluating a correction vector at \(\widetilde H,
\widetilde P,t\). The remaining error is tangent,
and coprimeness of \(H,P\) gives
\[
               (3H^2U,\ P^2 U+2H V,\ 3P V),
 \tag{10}
\]
with \(\deg U=j-4\), \(\deg V=j-2\).

At \(P=0\) fix the actual first component to \(A_0^3\)
using the local branch with \(t_0^2=A_0\) and \(t_0\ne0\).
Take \(\widetilde P\) and \(y\) as the other coordinates,
so \(y\) is free. Both the actual and reference charts
are regular: the leading coordinate system is
\(H,P,y\), with nonzero Jacobian there, and the
reference solves for \(t\) as a function of
\(\widetilde P,\epsilon\) on the fixed level.
Center the reference second component divided by \(A_0\)
as a scalar plus \(s^2\), with \(s\) leading \(P\).
By (10) the order-\(r\) difference in this component
on the fixed first level is \(2V\), since fixing the
first component subtracts the radial term \(P^2U/A_0\).
The third reference component has leading \(s^3\).
Below order fifteen its normalized Laurent coefficients
are independent of \(y\) by the determinant.
If \(V|_{P=0}\) were nonconstant in \(y\), its restriction
would be transcendental over the level field \(k(t_0)\).
The minimum-weight pole comparison for a quadratic
critical point with odd exponent three would contradict
regularity, since \(2r<15\). This is the transcendental
case of the
[binary pole argument](SEXTIC_KELLER_BINARY_SIMPLE_CUSP_DIVISORS.md#quadratic-critical-points-with-higher-odd-terms),
which uses transcendence of the leading error on the
free coordinate for this part of the proof.

Thus \(V|_{P=0}\) belongs to \(k(t)\) in the fraction
field of that irreducible divisor. Its polynomiality
on the divisor, restricted to the line \(y=z=0\),
forces it to be a polynomial in \(t\). Write
\(V=f(t)+PL\). Absorb the contribution of \(f(t)\) in
(10) directly into the reference. The rest is absorbed
by changing \(\widetilde H,\widetilde P\) by \(U,L\),
respectively. For \(j\le4\) the change \(U-\mu L\)
in \(\widetilde H-\mu\widetilde P\) is constant or
zero. For \(j=5\) the change is linear and its
restriction to \(t=0\) is precisely \(v\) in (8):
both the normal correction of degree five and \(f(t)\)
vanish at \(t=0\) by parity and homogeneity. Thus
\(v=0\) makes this change a polynomial in \(t\) as
well. The reference remains a function of only
\(\widetilde P,t\), so the descent continues.

At degree one the linear determinant equation gives
\(n\cdot(\mathbf a t+\mathbf b y+\mathbf c z)\in k[P,t]\),
with constant target vectors. Substituting
\(z=(P-y^2)/t\) forces \(n\cdot\mathbf b=n\cdot\mathbf c=0\).
The three entries of \(n\) are linearly independent,
so \(\mathbf b=\mathbf c=0\). The term \(\mathbf a t\)
and remaining constants are absorbed directly. The
map would factor through two polynomials, contradicting
the Keller condition. Therefore \(v\ne0\).

Finally \(v\ne0\) makes the first normal correction
zero by the preceding comparison. Now (10) at \(j=5\)
applies to the actual degree-five vector without any
normal subtraction. The same cusp comparison gives
\(V=PL+\gamma t^3\), and (8) gives \(U-\mu L=b x+a y\)
with \(a\ne0\). This proves the stated form.

## Degree-five divisibility when the cusp is at the square

**Theorem SexticKellerSquareTangentConicAtCuspFirstCorrection.**
For a characteristic-zero three-variable Keller map with highest
tuple \((P^3,PQ^2,Q^3)\), where \(P=xz+y^2\) and \(Q=x^2\),
write \(E_i\) for the degree-five forms. Then
\[
                 E_2=b xP^2+x^2 M_3,\qquad
                 E_3=c xP^2+x^3 M_2
\]
for scalars \(b,c\) and homogeneous forms \(M_3,M_2\)
of degrees three and two, respectively.

Work over the algebraic closure, put \(t=x\), and rescale
ordinarily in \(\epsilon\). This time take
\[
 A=F_{1,\epsilon}^{1/3},\qquad
 B=F_{2,\epsilon}/A,\qquad s'=B^{1/4}
\]
with leading terms \(P,t^4,t\). On the generic formal
chart \(A,s',y\), the two fixed components are \(A^3,A(s')^4\).
The remaining component has scalar coefficients
\(R_i(s')\in k(A)((s'))\) through order fourteen by the
same determinant argument, with \(R_0=(s')^6\).
On the chart \(A,t,z\) near the square divisor take
\(y_0^2=A\) as before. Write \(K=k(A)\).

First put \(v=[\epsilon]B|_{t=0}\). If \(v\ne0\), use
weights \(4,1\) on \(\epsilon,t\), so the minimum
part of \(B\) is \(t^4+\epsilon v\).
Terms \(\epsilon^i(s')^\ell\) with \(\ell\) a nonnegative
multiple of four may be subtracted by regular evaluations
as powers of \(B\). The least weight \(W\) of the
remaining terms below order fifteen is at most six
because of the leading term \((s')^6\), which cannot
be so removed. If \(v\) depends on \(z\), it is
transcendental over \(K(y_0)\). At weight \(W\) a term
with \(\ell<0\) contributes a pole at order \(i\).
For \(\ell\ge0\) not divisible by four the binomial
series first contributes a pole at
\(i+\lfloor\ell/4\rfloor+1\le2\). At the selected
order the contributions from different terms have
distinct powers of \(v\), so regularity is impossible.

Otherwise \(v=E_2(0,y_0,z)/A\) is independent of \(z\).
Homogeneity gives \(v=a y_0^3\) for some nonzero
\(a\in k\). Hence \(v\) is quadratic over \(K\) with
conjugate \(-v\). If \(W\) is divisible by four,
the first term at that weight already has a negative
power, contrary to regularity. Otherwise write
\(W=4d+j\), \(1\le j\le3\). Regularity after forward
substitution gives terms \(\epsilon^i t^{W-4i}\)
only for \(0\le i\le d\). Invert using
\(t^4=(s')^4-\epsilon v\) and divide by \((s')^W\).
In \(Z=\epsilon/(s')^4\) this gives
\((1-vZ)^{j/4}L_v(Z)\) modulo \(Z^{15}\), with
\(\deg L_v\le d\). It agrees with its conjugate,
because the normalized coefficients lie in \(K\).
Take fourth powers. The resulting polynomial equality
holds exactly since the degrees are at most \(W\le6\).
At \(Z=1/v\) one side has multiplicity \(j\) modulo
four and the other multiplicity zero modulo four,
unless both vanish identically. This is incompatible
with the choice of \(W\). Hence \(v=0\), so \(t\mid E_2\).

Next write \(u=(E_2/t)(0,y_0,z)/A\). If \(u\) depends
on \(z\), give \((\epsilon,t)\) weights \((3,1)\).
Now the minimum part of \(B\) is \(t^4+\epsilon t u\):
indeed its first-order part is \(E_2/P-t^4
[\epsilon]A/P\), and all terms of order at least two
have weight at least six. The same minimum-weight
comparison uses the binomial series of
\(1+\epsilon u/t^3\). There is again an offending
term of weight at most six from \(R_0=(s')^6\).
One with negative exponent gives a pole at its
own order, and one with nonnegative exponent not
divisible by four gives one at order
\(i+\lfloor\ell/3\rfloor+1\le3\).
The coefficient at that order cannot cancel for
transcendental \(u\). Thus \(u\) is independent of
\(z\); homogeneity gives \(u=b A\).
It follows that \(E_2-btP^2\) is divisible by \(t^2\).

For the last component, the normal contraction with
\((Q^3,-3P^2 Q,2P^3)\) belongs to \(k[P,t]\), by the
linear determinant equation and polynomial contraction
for \((P,t)\). Its degree is eleven, so it is congruent
to a constant multiple of \(tP^5\) modulo \(t^3\).
In the contraction, the first summand is divisible
by \(t^6\), and the second by \(t^3\) since \(t\mid E_2\).
Thus \(2P^3 E_3\) is congruent to that multiple of
\(tP^5\) modulo \(t^3\). Coprimeness of \(P,t\)
gives \(E_3=c tP^2+t^3 M_2\). This proves both
divisibilities. The formulas descend to the field
of the original coordinates, since the scalars
and polynomial quotients are unique there.

## Exclusion of cusps away from the square member

**Theorem SexticKellerSquareTangentConicAwayCuspExclusion.**
In characteristic zero there is no three-variable Keller map
with highest tuple \((H^3,HP^2,P^3)\), where \(P=xz+y^2\)
and \(H=x^2+\mu P\). Consequently neither of the two
nonimmersive cubic normal forms away from the square member
has a Keller completion.

Work over the algebraic closure. Put \(Q=x^2\), and let \(B,C\)
be the actual components \(F_2-\mu F_3,F_3\). Their ordinary
highest parts are \(QP^2,P^3\), respectively. By the
[first correction theorem](#degree-five-constraints-when-the-cusp-is-away-from-the-square),
their degree-five parts are
\[
 P^2(a y+b x)+2QP L+(2Q-\mu P)\gamma x^3,\qquad
 3P(PL+\gamma x^3),
 \tag{11}
\]
with \(a\ne0\), \(L\) linear, and \(b,\gamma\) scalars.

Give \(x,y,z\) weights \(1,2,w\), starting with \(w=1\) and
increasing to three. For \(1\le w<3\), the highest part of
\(P\) is \(y^2\). The contributions of ordinary degrees six
and five in (11) give the weight-ten part
\(y^4(x^2+a y)\) of \(B\) and the weight-twelve part
\(y^6\) of \(C\). The lower-degree terms in \(C\) stay below
weight twelve. A term of ordinary degree at most four with
positive \(z\)-exponent can tie in \(B\) before \(w=3\)
only if it is \(z^4\) at \(w=5/2\), or \(yz^3\) at
\(w=8/3\). Indeed the exponent of \(z\) must be at least
three by the weight bound, and these are the only monomials
of that degree attaining ten before three.

At the first such tie the weighted leading pair is
\[
        y^4(x^2+a y)+d z^4,\ y^6
        \quad\text{or}\quad
        y^4(x^2+a y)+d yz^3,\ y^6,\qquad d\ne0 .
\]
The pair is independent. On a geometric generic fibre,
\(y\) is fixed nonzero, and the remaining equation is
\(x^2=\alpha+\beta z^m\), with \(\alpha,\beta\ne0\) and
\(m=4\) or \(3\). The relative volume differential is a
nonzero scalar multiple of \(dz/x\). The polynomial on the
right has simple roots, and the differential is regular
on the complete normalization, including at infinity.
It cannot be exact, contradicting the
[weighted leading-pair slice](KELLER_LEADING_PAIR_RATIONAL_SLICE.md).
Thus no tie occurs before \(w=3\).

In weights \((1,2,3)\) put \(p=xz+y^2\) and \(W=x^2+a y\).
Now the leading pair has the general form
\[
 B_* =p^2 W+(b_1 p+b_2 y^2)z^2,\qquad C_*=p^3+c z^4 .
 \tag{12}
\]
Here the extra terms of \(B_*\) come from the only
ordinary-degree-at-most-four monomials of weight ten,
\(xz^3,y^2z^2\), and the extra term of \(C_*\) comes from
\(z^4\). In particular (11) contributes no further terms
at their leading weights. The pair in (12) is independent,
as is seen on \(z=0\) where its restriction is
\(y^4(x^2+a y),y^6\). Its geometric generic fibre must again
have exact relative volume differential. Take independently
generic values \(B_*=\beta_0,C_*=\gamma_0\) and work over the
algebraic closure of their field.

### The case \(c=0\)

The value \(p\) is fixed and nonzero on each component.
The locus \(x=0\) has at most finitely many points on a
generic fibre, since it fixes \(y^2=p\) and the remaining
\(B_*\) equation is either nonconstant in \(z\) or incompatible
with generic \(\beta_0\). On \(x\ne0\) substitute \(z=(p-y^2)/x\), put \(X=x^2\),
and set \(S(y)=(b_1p+b_2y^2)(p-y^2)^2\).
Equation (12) becomes
\[
                 p^2X^2+(a p^2y-\beta_0)X+S(y)=0 .
\]
The ambient volume after fixing \(p\) is proportional to
\(dx\wedge dy/x=dX\wedge dy/(2X)\). Since the displayed
equation is \(X(B_*-\beta_0)=0\), division by the derivative
in \(X\) gives, up to a nonzero scalar, the differential
\(dy/V\), where
\[
             V=2p^2X+a p^2y-\beta_0,\qquad
             V^2=(a p^2y-\beta_0)^2-4p^2S(y).
 \tag{13}
\]
If \((b_1,b_2)\ne(0,0)\), the right side of (13)
is a squarefree polynomial
of degree four or six. For squarefreeness, a zero cannot
occur at \(S=0\) for generic \(\beta_0\). Off \(S=0\),
on the curve \(v^2=S(y)\), the zero corresponds to the
generic value \(\beta_0\) of \(a p^2y+2pv\). On any
component where this function is nonconstant, its generic
value has a simple fibre, and constant components miss
\(\beta_0\). Since \(y\) is a local parameter off \(S=0\)
and the two signs of \(v\) give distinct values there,
the polynomial in (13) has no repeated zeros.
The form \(dy/V\) is nonzero and regular on the complete
normalization. Its pullback to the finite cover adjoining
\(x\) is likewise regular and nonzero. This contradicts
exactness on the original generic fibre.

If \(b_1=b_2=0\), fixing the pair fixes \(p,W\) and then
\(y=(W-x^2)/a\). The remaining generic fibre has \(x\)
as a nonzero free parameter and \(z=(p-y^2)/x\).
The relative form is proportional to \(dx/x\): the
Jacobian of \(p,W,x\) with respect to \(x,y,z\) is
a nonzero scalar times \(x\). The simple poles of
\(dx/x\) rule out exactness.

### The case \(c\ne0\)

The locus \(z=0\) is finite on a generic fibre, by the
restriction of the pair there. On \(z\ne0\), use \(x=(p-y^2)/z\). Then \(p,z\) satisfy
\(p^3+c z^4=\gamma_0\), and the equation for \(y\) is
\[
  B_*=\frac{p^2}{z^2}(p-y^2)^2+a p^2 y+(b_1p+b_2y^2)z^2
      =\beta_0 .
 \tag{14}
\]
In these coordinates the relative differential is, up to
sign, \(dz/(3z p^2\partial_y B_*)\), and one can use
\(dp/(z(\partial_z C_*)\partial_y B_*)\) instead where
\(z\ne0\) and \(C_*=p^3+c z^4\) is differentiated in
the \(p,z\) chart. On the affine generic fibre it is regular
by generic smoothness. To check the points outside it,
use the normalization of the curve \(p^3+c z^4=\gamma_0\)
as the base of the finite extension given by (14).

At a point over \(z=0\), \(p\) is a unit. A pole of \(y\)
would give an uncancelled pole from the first term of (14).
For \(y\) finite, that term also forces \((p-y^2)/z\)
to be finite. Hence all three original coordinates are
finite there. At finite nonzero \(z\) with \(p\ne0\), the
quartic equation in \(y\) has unit leading coefficient and
\(y,x\) are finite.

At a point over \(p=0\), \(z\) is a unit and \(p\) is a
base parameter. Only branches where \(y\) has a pole
need checking. If \(b_2\ne0\), the two potentially
dominant terms are \(p^2 y^4/z^2\) and \(b_2 y^2z^2\).
Their cancellation forces \(p y\) to have a finite nonzero
limit. In the equation multiplied by \(p^2\) the resulting
nonzero root for \(p y\) is simple, since the leading
polynomial is \(Y^4/z_0^2+b_2 z_0^2Y^2\), where \(z_0\)
is the nonzero limit of \(z\). Thus \(\partial_y B_*\)
has order \(-1\) in \(p\), and \(dp/\partial_y B_*\)
is regular. If \(b_2=0\), then \(p^2y^4/z^2\) must
balance \(\beta_0\ne0\), since the other terms have
higher order. In a base change with parameter \(\sqrt p\),
the limiting \(Y=\sqrt p\,y\) is a simple nonzero root
of \(Y^4/z_0^2=\beta_0\). Hence \(\partial_y B_*\)
has order \(1/2\) in \(p\), and the same differential
is regular after the base change. This checks the points
at \(p=0\), because the other factors of the \(dp\) formula
are units there.

At the point over \(z=\infty\) on the normalized base
choose a parameter \(\tau\) with \(z=\tau^{-3}\), using
that the pole order four of the right side of
\(p^3=\gamma_0-cz^4\) is coprime to three. The
equation for \(p\) gives \(p=\tau^{-4}p_1(\tau)\), with
\(p_1=p_0+O(\tau^{12})\) and \(p_0^3=-c\ne0\).
A pole order of \(y\) larger than two in \(\tau\) would
make the quartic-in-\(y\) term alone dominant in (14).
Thus \(Y=\tau^2 y\) is finite. Equation (14) takes the form
\[
 g(Y)+O(\tau^{12})=\beta_0\tau^{10},\qquad
 g(Y)=p_0^2(p_0-Y^2)^2+a p_0^2Y+b_1 p_0+b_2Y^2 .
\]
The error notation means a polynomial in \(Y\) whose
coefficients vanish to order at least twelve. If a finite
limit \(Y_0\) is a root of \(g\) of multiplicity \(m\),
then \(Y-Y_0\) has order \(10/m\) in \(\tau\), since
\(\beta_0\ne0\). The derivative in \(Y\) of the left side
has order \(10(m-1)/m\). Consequently
\(\partial_y B_*\) has order \(-8+10(m-1)/m\).
In the formula with \(dz\), the factor \(dz/z\) is a
nonzero multiple of \(d\tau/\tau\). The resulting
differential has order \(15-10(m-1)/m>0\) in \(\tau\)
before accounting for ramification. It is therefore
regular on the normalized fibre. These are all places
outside the affine fibre. The nonzero relative
differential is regular on its complete normalization
and cannot have a rational primitive. This final
contradiction proves the theorem.
