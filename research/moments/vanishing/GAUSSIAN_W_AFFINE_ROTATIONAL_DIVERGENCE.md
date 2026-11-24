# Gaussian divergence with rotationally related dominant singularities

Let \(m\ge2\). Consider the actual Gaussian unit for the phase
\[
P=W\bigl(h(Z)+c(Z)Y\bigr)+aY^2+u(Z)Y+k(Z)
\]
in the coordinates and Wick convention of the [original-source completion](GAUSSIAN_W_AFFINE_SINGLE_ROOT_DIVERGENCE.md).
Suppose
\[
c(z)=C A(z/\zeta),\qquad C\zeta\ne0,\qquad
A(z)=B(z^m),\qquad B(0)>0,
\]
where \(B\) is a nonconstant polynomial with nonnegative real coefficients
and \(\gcd\{j\ge1:[v^j]B(v)>0\}=1\). For every \(m\)-th root of unity
\(\omega\), assume
\[
h(\omega z)=\omega h(z),\quad u(\omega z)=\omega u(z),\quad
k(\omega z)=k(z),\quad u(\omega z)^2=u(z)^2.
\]
The last condition is automatic for \(m=2\) and requires \(u=0\) for
\(m>2\). The complex constant \(a\) is arbitrary. Then the original
unit \(f=\mathbb E e^{tP}\) diverges at zero, its actual image has rank at
least two over \(\mathbb C(t)\), and the phase is neither pure nor safe.
In particular the result includes \(c(z)=1+z^2\) with arbitrary odd
\(h,u\), even \(k\), and arbitrary \(a\). Its two baseline dominant points
are \(X=\pm 1/2\). A multiple-root example is \(c=(1+z^2)^4(1+2z^2)\) and \(h=z(1+z^2)^2\), where each of \(z=\pm i\) has orders \((4,2)\) in \((c,h)\).

The proof addresses the [multiple-dominant-point obstruction](GAUSSIAN_W_AFFINE_SELECTED_DOMINANT_DIVERGENCE.md#4-the-concrete-bottleneck-beyond-these-criteria)
on the rotationally invariant locus. The quotient by the common rotations
has a unique selected singularity to which singularity analysis applies.
The argument uses the actual unit and the original fixed amplitude \(P\).

## The quotient of the actual source kernel

Write \(D_a=1-2at\), and use its square root equal to one at zero. The
original completion gives
\[
\begin{split}
z&=t h(z)+t^2 c(z)u(z)/D_a+Xc(z)/\sqrt{D_a},\\
J&=\partial_z\bigl(z-t h(z)-t^2 c(z)u(z)/D_a-Xc(z)/\sqrt{D_a}\bigr),\\
K(t,X)&=\frac{\exp(tk(z)+t^2u(z)^2/(2D_a))}{\sqrt{D_a}\,J},\qquad
f(t)=\mathbb E_S K(t,tS),
\end{split}
\]
where \(S\) is standard real Gaussian and the last identity is formal.
Write \(v=z^m\), \(U=X^m\), \(c(z)=\widehat c(v)\), and
\(t h(z)+t^2 c(z)u(z)/D_a=z b(t,v)\). The assumptions give \(b(0,v)=0\)
and an analytic numerator \(E(t,v)=\exp(tk(z)+t^2u(z)^2/(2D_a))\).
Consequently the selected \(v\) solves
\[
v=U F(t,v),\qquad
F(t,v)=\left(\frac{\widehat c(v)}{\sqrt{D_a}(1-b(t,v))}\right)^m.
\]
On the selected graph,
\[
J=(1-b(t,v))(1-U F_v(t,v)),\qquad
K(t,X)=\mathcal K(t,X^m)
=\frac{E(t,v)}{\sqrt{D_a}(1-b(t,v))(1-U F_v(t,v))}.
\]
For the Jacobian identity one can differentiate the original residual
in \(z\), then use \(X\widehat c/\sqrt{D_a}=(1-b)z\) and
\(v=UF\). Any intermediate division by \(\widehat c\) can be cleared
from that identity. The formula holds at the origin by analyticity.

At \(t=0\), rescale \(v=\zeta^m w\) and \(U=\zeta^m C^{-m} V\).
The equation is \(w=V B(w)^m\), and the kernel is
\((1-V(B^m)'(w))^{-1}\). The polynomial \(B^m\) has nonnegative
coefficients, positive constant coefficient, degree at least two,
and the same support gcd as \(B\). The equality of gcds follows because
the support contains the support of \(B\), and every new exponent is a
sum of exponents of \(B\). The positive implicit-function argument in
[Section 3 of the selected-dominant theorem](GAUSSIAN_W_AFFINE_SELECTED_DOMINANT_DIVERGENCE.md#3-an-explicit-positive-coefficient-family-with-multiple-roots)
gives a unique selected dominant point of the quotient, located at a
nonzero \(U_0\). More explicitly, for \(F_0=B^m\) the positive critical
point \(\tau\) is the unique root of \(F_0(\tau)-\tau F_0'(\tau)=0\).
The inverse with nonnegative coefficients takes values at most \(\tau\)
in absolute value on its convergence disk. Equality in the triangle
inequality at a boundary critical point forces every supported monomial
to have the phase of the constant. Support gcd one gives the unique
positive critical point, and \(F_0''(\tau)>0\) makes it simple.

## The moving singularity and the Wick series

The critical point persists analytically for small \(t\). Let its
value be \(U_*(t)\) with \(U_*(0)=U_0\). Put
\(\eta(t)=(U_*(t)/U_0)^{1/m}\) with \(\eta(0)=1\) and
\(T=t/\eta(t)\), an invertible analytic parameter change with inverse
\(t(T)\). Then
\[
f(t(T))=\mathbb E_S \widetilde K(T,TS),\qquad
\widetilde K(T,X)=\mathcal K(t(T),\eta(t(T))^m X^m).
\]
In the quotient \(U=X^m\), the dominant point is now fixed at \(U_0\).
At \(t=0\) the selected graph is compact in a slightly extended slit disk
with one simple quadratic inverse point. Away from a small critical
neighborhood the implicit derivative is a unit. The implicit-function
theorem on finitely many neighborhoods and the local quadratic inverse
construction give the corresponding joint analytic slit continuation for
all small complex \(T\). On the bounded graph, \(1-b\) stays a unit after
shrinking the parameter disk. The kernel has the uniform local form
\(P_0(T,U)+Q_0(T,U)(1-U/U_0)^{-1/2}\), with \(Q_0(0,U_0)\ne0\), since its
numerator is analytic and nonzero and \(1-UF_v\) has a simple zero in
the square-root coordinate.

The slit-contour estimate of the [selected-dominant proof, (S3)–(S5)](GAUSSIAN_W_AFFINE_SELECTED_DOMINANT_DIVERGENCE.md#1-a-selected-dominant-singularity-criterion),
applied in \(U\), gives on a closed \(T\) disk an upper bound
\(O(|U_0|^{-j}(j+1)^{-1/2})\) on the coefficient of \(U^j\). At \(T=0\)
this coefficient is
\[
Q_0(0,U_0)U_0^{-j}\binom{2j}{j}4^{-j}
+O\bigl(|U_0|^{-j}(j+1)^{-3/2}\bigr).
\]
Thus \([X^n]\widetilde K(T,X)=0\) unless \(m\mid n\). Set
\(R_X=|U_0|^{1/m}\). For \(m\mid n\) the upper bound in \(X\) is
\(O(R_X^{-n}(n+1)^{-1/2})\). For such \(n\) at \(T=0\), the coefficient
has magnitude bounded below by a positive multiple of that quantity
for large \(n\).

Take \(q\to\infty\) with \(m\mid 2q\). In the coefficient of \(T^{2q}\)
of \(f(t(T))\), the \(X^{2q}\) term has magnitude comparable to
\[
L_q=(2q-1)!! R_X^{-2q}(2q+1)^{-1/2}.
\]
Cauchy's estimate on the same parameter disk bounds the sum of lower
\(X^{2j}\) contributions by \(C\sum_{j<q}L_j R_T^{-2(q-j)}=O(L_q/q)\).
Indeed \(L_{j+1}/L_j\) is comparable to \(j+1\). For \(j\ge q/2\) this
gives a geometric bound with ratio \(O(1/q)\) on each backward step
from \(q\). For \(j<q/2\) the factorial ratio absorbs the fixed
exponential and the number of terms. The surviving subsequence grows
faster than every exponential. This proves divergence, which survives
the analytic change of time. A finite actual image, purity or safety
would make the unit analytic by the original-source theorem. If \(f'\)
were a rational multiple of \(f\), the logarithmic derivative would be
regular at zero since \(f(0)=1\), and the resulting regular first-order
equation would make \(f\) analytic. Hence \(f,f'=\Phi[P]\) have rank at
least two.

This proof uses the same selected-sheet contour estimates as the
aperiodic case. Once the rotations are broken, the quotient equation
and the single normalization may both fail. The argument also makes
no determination about a regular submodule independent of the unit.

## Generic coefficients within the rotation class

The same conclusion holds for \(c(z)=A(z^m)\) with arbitrary complex
coefficients under the following finite sufficient condition. Assume
\(A(0)\ne0\), \(\deg A=k\ge1\), \(A\) has simple roots, all roots \(w_i\)
of \(A(w)-mw A'(w)\) are simple, and the values
\[
V_i=w_i/A(w_i)^m
\]
have distinct absolute values. The \(w_i\) are nonzero and disjoint from
the zeros of \(A\). They are precisely the critical points of \(w/A(w)^m\)
with nonzero finite image. The selected quotient inverse \(w=U A(w)^m\)
has finite convergence radius. Otherwise it would be entire and bounded
on \(|U|\ge1\), since the polynomial \(A(w)^m\) has degree at least two
and the defining equation bounds any root uniformly there. Liouville
would make \(w\) constant, while \(w'(0)=A(0)^m\ne0\).
The same polynomial bound at any nonzero \(U\) rules out divergence of
the selected inverse to infinity at its convergence circle. At a
boundary singularity the inverse must land at a critical point of
\(w/A(w)^m\). Distinct moduli give a unique such boundary value,
with a simple quadratic local inverse. The slit collar now follows
by compact continuation at the regular boundary points and local
quadratic inversion at the unique singularity, as in the reference
single-point proof. The quotient argument above applies.

This finite condition is open and dense in the complex coefficients
of each fixed degree \(k\) with \(A(0)\ne0\). The simple-root conditions
hold on a Zariski-open dense set, as the example \(A=1+w^k\) shows.
On that set the derivative of each local critical value under a
variation \(\delta A\) is
\[
\delta V_i=-mw_i A(w_i)^{-m-1}\delta A(w_i).
\]
Since the \(k\) points are distinct, interpolation of degree at most
\(k-1\) makes these \(k\) variations independently arbitrary. The
critical-value map is a local submersion and a tie of absolute values
has empty interior. This gives a Euclidean-open dense family in every
rotation class, even though conjugate-by-rotation critical values in
the original \(X\) plane always tie.

## Exact control

For \(a=h=u=k=0\) and \(c(z)=1+z^2\), the quadratic relation
\(z=X(1+z^2)\) gives \(K(0,X)=(1-4X^2)^{-1/2}\), and \(K\) is
independent of \(t\). The original unit therefore satisfies
\[
[t^{2j}]f=(2j-1)!!\binom{2j}{j},\qquad [t^{2j+1}]f=0.
\]
The original Wick expansion gives the same identity directly.
In \(P^{2j}=W^{2j}(1+Z^2)^{2j}Y^{2j}\), the \(Z^{2j}\) term gives
the factor \(\binom{2j}{j}\), and the \(Z,W\) contraction gives
\((2j)!\), which cancels the exponential coefficient denominator.

## A Hamburger versus Stieltjes test on the actual Wick coefficients

For \(P_m=W(1+Z^m)Y\), \(m\ge1\), define the normalized coefficients
\[
s_q=\frac{[t^{2q}]\mathbb E e^{tP_m}}{(2q-1)!!},\qquad s_0=1.
\]
Exactly
\[
s_q=\begin{cases}\binom{2q}{2q/m},&m\mid 2q,\\0,&m\nmid2q.\end{cases}
\]
This follows by expanding \(P_m^{2q}\): only the monomial with equal
\(Z,W\) exponents survives, and the \(Z,W\) Wick contraction cancels
the exponential coefficient denominator. The sequence is Stieltjes
exactly for \(m=1,2\), and Hamburger exactly for \(m=1,2,4\).

For \(m=1\) use the point mass at 1. For \(m=2\), the random variable
\(4U\) for \(U\) with beta parameters \((1/2,1/2)\) has moments
\(4^q(1/2)_q/q!=\binom{2q}{q}\). For \(m=4\), take independent beta variables
\(U_i\) with parameter pairs
\((1/4,1/12),(1/2,1/6),(3/4,1/4)\). The rising factorial product gives
\[
\mathbb E\left[\left(\frac{256}{27}U_1U_2U_3\right)^j\right]
 =\left(\frac{256}{27}\right)^j
 \frac{(1/4)_j(1/2)_j(3/4)_j}{(1/3)_j(2/3)_j(1)_j}
 =\binom{4j}{j}.
\]
The last equality follows by grouping the factors of \((4j)!\) and
\((3j)!\) according to residues modulo 4 and 3. With an independent
fair sign, the square root of this product represents the full
sequence, whose odd moments are zero. A measure supported on
\([0,\infty)\) with first moment zero must be supported at zero, and
here \(s_2=4\). This excludes Stieltjes at \(m=4\).

For all other \(m\), \(s_2=0\) but \(s_m=\binom{2m}{2}>0\).
A positive measure on the real line with second moment zero is
supported at zero, which is a contradiction. This is a finite
obstruction for the stated sequence and does not require analytic
continuation. It illustrates why nonnegative normalized coefficients
of an actual Wick unit do not certify positive-measure moment
representability.

The distinction between real support and nonnegative real support
was prompted by the external
[Apéry Hamburger classification](https://github.com/1337isnot1337/solved-conjectures/blob/758066055ef49c2c66cb1087918ae62ad75e7ac1/results/sokal-apery-moment-threshold/HAMBURGER.md).
The elementary beta product here supplies its own proof and uses no
Apéry classification as a hypothesis. The distribution is for the
normalized even coefficients and does not alter the divergence of the
unscaled original unit.

## Moment-transform provenance and hypotheses

The [Apéry threshold proof](https://github.com/1337isnot1337/solved-conjectures/blob/758066055ef49c2c66cb1087918ae62ad75e7ac1/results/sokal-apery-moment-threshold/PROOF.md)
uses a compact-support obstruction. If a real sequence with a uniform
bound \(|s_{2n}|\le C R^{2n}\) has a positive representing measure, any
mass on \(|t|\ge R+\epsilon\) contradicts that bound for large \(n\).
The generating function germ is then \(\int(1-tz)^{-1}d\mu(t)\).
It extends at every nonreal point for real support, and also at every
negative real point for nonnegative support. Checking the actual
selected germ is essential because a singularity of another algebraic
branch does not obstruct this continuation.

For the external Apéry application, the Legendre addition formula
and a torus change of variables give the argument
\(Z=1+8x\cos^2 a+8\sqrt{x}\cos a\cos b\). The exceptional level
\(1+8x-8\sqrt{x}\) lies in \((-1,1)\) for \(0<x<1\) except at
\(x=1/4\), where it equals \(-1\). For \(x<1/4\) it is the minimum,
and the left Abel transform is zero just before the level and
positive just after. For \(1/4<x<1\), the derivatives
\(-8\sqrt x\sin a(2\sqrt x\cos a+\cos b)\) and
\(-8\sqrt x\cos a\sin b\) give saddle eigenvalues
\(-8\sqrt x(2\sqrt x-1)\) and \(8\sqrt x\). The positive local
logarithmic density yields a one-sided Abel cusp with nonzero
constant, using
\[
\int_0^\infty v^{-1/2}\log(1+1/v)\,dv
 =2\int_0^\infty\frac{dv}{\sqrt v(1+v)}=2\pi.
\]
The inverse Joukowski map sends the level in \((-1,1)\) to a
nonreal point of the unit circle. At \(x=1/4\) the local form
\(Z+1\asymp a^4+(b-\pi)^2\) gives sublevel mass of order
\(h^{3/4}\). The Abel difference has order at least \(h^{1/4}\),
which gives a square-root-scale change along \(z=-r\) as \(r\uparrow1\). The obstruction there excludes
nonnegative support and permits real support. The
[Hamburger companion](https://github.com/1337isnot1337/solved-conjectures/blob/758066055ef49c2c66cb1087918ae62ad75e7ac1/results/sokal-apery-moment-threshold/HAMBURGER.md)
uses a separate two-cut sign calculation and a positive mixture for
the exceptional real representation. None of those elliptic-period
representations is a premise of the Gaussian theorem. The actual
Gaussian application uses the selected inverse of the Wick kernel,
followed by the exact factorial weighting, with no assumed measure.

The [Ferreira–Simon catalogue proof](https://github.com/1337isnot1337/solved-conjectures/blob/758066055ef49c2c66cb1087918ae62ad75e7ac1/results/ferreira-simon-beta-prime/PROOF.md)
addresses Conjecture 1 of
[Ferreira and Simon](https://arxiv.org/abs/2108.09244v2). It factors
the Mellin transform of a sum of independent beta-prime variables
of parameters \(a,b\), for \(a>0\) and \(0<b<1/2\). In the key
conditional integral, putting \(v=(1-y)/(1+y)\) in
\[
\int_0^1 v^{2b-s-1}(1-v^2)^{-b-1/2}(u+v)^s\,dv
\]
cancels the powers of \(1+y\), since
\(1-v^2=4y/(1+y)^2\). The remaining beta parameters are
\(1/2-b\) and \(2b-s\), positive on \(-2a<s<0\), and the
hypergeometric denominator \(b+1/2-s\) equals one plus the first
upper parameter minus the second. Thus the quadratic substitution
\(x=(1-u)/(1+u)\mapsto4x/(1+x)^2=1-u^2\) is valid on the
real branch as written. The product law comparison uses an open
Mellin strip and uniqueness for the log-transformed law. This
mechanism is different from assuming all integer moments exist:
the sum has no finite first moment, because each beta-prime
density has tail of order \(x^{-b-1}\) and \(b<1/2\).
Consequently its raw moment law supplies no positive infinite
moment sequence for the Gaussian DAG. The beta factors in the
normalized Wick example above have bounded support and are
verified by their own elementary beta integrals.
