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
h(\omega z)=\omega h(z),\qquad u(\omega z)=\omega u(z).
\]
The polynomial \(k\) and the complex constant \(a\) are arbitrary.
Then the original unit \(f=\mathbb E e^{tP}\) diverges at zero, its
actual image has rank at least two over \(\mathbb C(t)\), and the phase
is neither pure nor safe.

For \(c(z)=1+z^2\) this includes arbitrary odd \(h,u\) and arbitrary
\(k\), including the even \(k\) subfamily where the full kernel is
rotation invariant. Its two baseline dominant points are \(X=\pm 1/2\).
A multiple-root example is \(c=(1+z^2)^4(1+2z^2)\) and
\(h=z(1+z^2)^2\), where each of \(z=\pm i\) has orders \((4,2)\) in
\((c,h)\). For \(m=3\), the choices \(c=1+z^3\), \(h=u=z\), and
\(k=z+z^2\) give a noninvariant numerator: neither \(k(z)\) nor
\(u(z)^2\) is invariant under multiplication of \(z\) by a primitive
cube root. The earlier invariant-numerator subclass, with
\(k(\omega z)=k(z)\) and \(u(\omega z)^2=u(z)^2\), is included.

The proof addresses the [multiple-dominant-point obstruction](GAUSSIAN_W_AFFINE_SELECTED_DOMINANT_DIVERGENCE.md#4-the-concrete-bottleneck-beyond-these-criteria)
when the selected inverse equation has rotational symmetry. Its quotient
has a unique selected singularity. Decomposing the numerator by its
powers modulo \(m\) retains that singularity bound for the original
unit and the original fixed amplitude \(P\).

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
\(t h(z)+t^2 c(z)u(z)/D_a=z b(t,v)\). Then \(b(0,v)=0\) and
\[
z=X\ell(t,v),\qquad v=U F(t,v),\qquad
\ell(t,v)=\frac{\widehat c(v)}{\sqrt{D_a}(1-b(t,v))},\quad F=\ell^m.
\]
The numerator is entire in \(z\) for small \(t\). Group its Taylor
series by exponent modulo \(m\):
\[
\exp(tk(z)+t^2u(z)^2/(2D_a))=\sum_{r=0}^{m-1} z^r E_r(t,z^m).
\]
The functions \(E_r\) are jointly analytic in \(t,v\) for small \(t\)
and bounded \(v\), by Cauchy estimates on arbitrarily large fixed
\(z\) circles. At \(t=0\), \(E_0=1\) and every other \(E_r=0\).
On the selected graph,
\[
J=(1-b(t,v))(1-U F_v(t,v)),\qquad
K(t,X)=\sum_{r=0}^{m-1} X^r\mathcal K_r(t,X^m),\qquad
\mathcal K_r(t,U)=
\frac{\ell(t,v)^r E_r(t,v)}{\sqrt{D_a}(1-b(t,v))(1-U F_v(t,v))}.
\]
For the Jacobian identity, write the residual as
\((1-b(t,z^m))(z-X\ell(t,z^m))\). On the selected graph the
derivative of its first factor is multiplied by zero. The derivative
of its second factor is \(1-Xm z^{m-1}\ell_v=1-U F_v\),
using \(z=X\ell(t,v)\). Only \(1-b\) and \(\sqrt{D_a}\) need be
units for this computation.

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
\widetilde K(T,X)=K(t(T),\eta(t(T))X)
=\sum_{r=0}^{m-1}X^r\widetilde{\mathcal K}_r(T,X^m),\qquad
\widetilde{\mathcal K}_r(T,U)=\eta(t(T))^r
 \mathcal K_r(t(T),\eta(t(T))^m U).
\]
In the quotient \(U=X^m\), the dominant point is now fixed at \(U_0\).
At \(t=0\) the selected graph is compact in a slightly extended slit disk
with one simple quadratic inverse point. Away from a small critical
neighborhood the implicit derivative is a unit. The implicit-function
theorem on finitely many neighborhoods and the local quadratic inverse
construction give the corresponding joint analytic slit continuation for
all small complex \(T\). On the bounded graph, \(1-b\) stays a unit after
shrinking the parameter disk. Every quotient term has the uniform
local form \(P_r(T,U)+Q_r(T,U)(1-U/U_0)^{-1/2}\), since its numerator
is analytic and \(1-UF_v\) has a simple zero in the square-root
coordinate. For \(r=0\) at \(T=0\) the numerator is nonzero, so
\(Q_0(0,U_0)\ne0\).

The slit-contour estimate of the [selected-dominant proof, (S3)–(S5)](GAUSSIAN_W_AFFINE_SELECTED_DOMINANT_DIVERGENCE.md#1-a-selected-dominant-singularity-criterion),
applied to each quotient term in \(U\), gives on a closed \(T\) disk
the uniform bound \(O(|U_0|^{-j}(j+1)^{-1/2})\) on its coefficient of
\(U^j\). At \(T=0\), the \(r=0\) term has coefficient
\[
Q_0(0,U_0)U_0^{-j}\binom{2j}{j}4^{-j}
+O\bigl(|U_0|^{-j}(j+1)^{-3/2}\bigr).
\]
For any \(n=mj+r\), exactly the \(r\)-th quotient term contributes to
the \(X^n\) coefficient. With \(R_X=|U_0|^{1/m}\), this gives
\([X^n]\widetilde K(T,X)=O(R_X^{-n}(n+1)^{-1/2})\) uniformly for
**all** \(n\), including those that occur only because the numerator
breaks the rotation. At \(T=0\), the coefficients with \(m\mid n\)
have magnitude bounded below by a positive constant times the same
quantity for large \(n\).

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
aperiodic case. If rotation symmetry is lost in the inverse equation,
the quotient and the single normalization may both fail. The argument also makes
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

The noninvariant \(m=3\) example can be checked at \(a=1/2\).
Put \(D=1-t\) and \(b=t+t^2(1+z^3)/D\). Expanding the exact
circular residue of the completed kernel in \(X\) gives
\[
[X^n]K(t,X)
=[z^n]\frac{(1+z^3)^n
     \exp\bigl(t(z+z^2)+t^2z^2/(2D)\bigr)}
 {D^{(n+1)/2}(1-b)^{n+1}} .
\]
Indeed this is the \(X^n\) term in the geometric expansion of the
residual denominator \(z(1-b)-X(1+z^3)/\sqrt D\).
In particular \([tX]K=1\), so the original kernel really does have a
nontrivial residue class modulo three. Independent exact contractions
of \(P^j/j!\) and of the displayed residue after \(X=tS\) gave the
same coefficients through order eight:
\[
1,\quad 3/2,\quad 23/8,\quad 107/16,\quad 2155/128,\quad
14085/256,\quad 436435/1024,\quad 6371163/2048,\quad
608546355/32768 .
\]
The finite calculation checks the source and its asymmetric amplitude.
Divergence follows from the uniform coefficient majorant above.

## Exact actual rank for \(P_m=W(1+Z^m)Y\)

Let \(m\ge2\), let \(K\) be a characteristic-zero field, and use the
normalized original Wick map \(\Phi[A]=\mathbb E[A e^{tP_m}]\) over \(K\).
Put \(M=K[t]\Phi[K[Z,W,Y]]\), \(V=K(t)M\), and
\[
\ell=\operatorname{lcm}(2,m),\quad r=\ell/2,\quad d=\ell/m,\quad
e=\ell-d,\quad x=t^\ell,\quad \vartheta=\frac{t}{\ell}\partial_t .
\]
Then \(f=\Phi[1]=g(x)\), where
\[
g(x)=\sum_{q\ge0}a_q x^q,\qquad
a_q=\frac{(\ell q)!^2}{2^{r q}(r q)!(d q)!(e q)!}.
\]
Set \(p=3m/2-1\) for even \(m\), and \(p=3m-3\) for odd \(m\).
The differential subspace generated by the unit has dimension \(p\)
over \(K(t)\). Bases of the full actual space \(V\) are
\[
\begin{array}{ll}
(1,g,\vartheta g,\ldots,\vartheta^{p-1}g),& m\text{ even},\\[2mm]
(1,g,\vartheta g,\ldots,\vartheta^{p-1}g,\kappa),& m\text{ odd},
\end{array}
\qquad
\kappa(x)=\sum_{q\ge0}\frac{a_q}{q+1/2}x^q .
\]
Thus the full rank is \(3m/2\) for even \(m\), and \(3m-1\) for odd
\(m\). With \(R=K[t,t^{-1}]\) and either displayed basis \(\mathcal B\),
the exact polynomial-time lattice is
\[
M[1/t]=R\mathcal B,\qquad M=(R\mathcal B)\cap K[[t]].
\]
Over \(\mathbb C\), the elements of \(V\) that converge as Laurent
germs are exactly \(\mathbb C(t)\), and the convergent elements of
\(M\) are exactly \(\mathbb C[t]\). The maximal globally regular
actual differential submodule is the constant connection
\(\mathbb C(t)\cdot1\). These are results for the explicit phase \(P_m\).
They do not address the general higher-rank realization gap
\(\texttt{GAP-GAUSSIAN-REGULAR-TRANSCENDENTAL}\).
The case \(m=1\) is outside this derivation because \(e=0\). It is
the same phase as the \(m=1\) case of the
[binomial full-image theorem](GAUSSIAN_BINOMIAL_MIXED_FULL_IMAGE.md#6-the-m1-exception-is-completely-classified),
which gives the global lattice with basis \((1,f_1)\).

### The exact order of the unit operator

Write \(G_s=\{1/s,2/s,\ldots,1\}\). The factorial ratio gives
\[
\frac{a_{q+1}}{a_q}=C\frac{\mathcal P(q)}{\mathcal Q(q)}
\]
for a nonzero constant \(C\), where the numerator parameter multiset
before cancellation is two copies of \(G_\ell\), and the denominator
is \(G_r+G_d+G_e\). A parameter \(\alpha\) denotes the factor \(q+\alpha\).
Both \(G_r\) and \(G_d\) lie in \(G_\ell\) and cancel first.
Since \(G_e\cap G_\ell=G_{\gcd(e,\ell)}=G_d\), precisely
\(d-\gcd(r,d)\) further parameters cancel. This gives
\(\deg\mathcal P=3r-2d+\gcd(r,d)=p\) and
\(\deg\mathcal P-\deg\mathcal Q=r>0\). Three denominator parameters
equal one and two numerator parameters equal one, so \(\mathcal Q(-1)=0\).
The roots of \(\mathcal Q\) are simple. For even \(m\), its parameter
set is \(G_e\); for odd \(m\), it is \(G_e\setminus\{1/2\}\), and
\(1/2\) is also absent from the numerator parameters. The reduced
numerator and denominator parameters lie in \((0,1]\) and no cross
pair differs by an integer.

The coefficient recurrence gives the operator
\(\mathcal Q(\vartheta-1)-Cx\mathcal P(\vartheta)\), including the
constant coefficient by \(\mathcal Q(-1)=0\). For the lower bound,
clear rational denominators in any equation over \(K(x)\) with a
rational right side. In Euler form the equation then has a polynomial
right side and a left side \(\sum_{j=0}^J x^j L_j(\vartheta)g\)
with \(L_J\ne0\). At index \(n\) large enough, divide the coefficient
equation by \(a_{n-J}\) to obtain a rational identity in \(n\). For
\(j<J\), the coefficient ratio
\[
\frac{a_{n-j}}{a_{n-J}}
 =C^{J-j}\prod_{h=0}^{J-j-1}
   \frac{\mathcal P(n-J+h)}{\mathcal Q(n-J+h)}
\]
vanishes at \(n=J-\alpha\) with at least the full multiplicity of
each numerator factor \(q+\alpha\). No denominator in the product
vanishes there. The term at \(j=J\) therefore forces
\(\mathcal P(q)\mid L_J(q)\). Every such equation has order at least
\(p\). A polynomial-coefficient equation over \(K(t)\) splits by
the exponents of \(t\) modulo \(\ell\), because \(g\) and its Euler
derivatives are series in \(x\). Each residue gives an equation over
\(K(x)\) of no greater order. This proves the exact cyclic dimension
over \(K(t)\) and independence of the constant line. The dimension
is an actual-image statement because
\(\partial_t^j\Phi[1]=\Phi[P_m^j]\) for every \(j\).

### All original markers and the odd extra period

With \(c(z)=1+z^m\), the original conditional circular identity gives
\(\Phi[1-t c'(Z)Y]=1\). Conditional on \(Y\), the selected inverse
\(z=t c(z)Y\) has Jacobian \(1-t c'(z)Y\), which is exactly this
time-dependent amplitude on the selected graph. In particular all
polynomials in \(t\) belong to \(M\).

For \(Z^aW^bY^c\), with \(a,b,c\) now nonnegative integer exponents,
the coefficient at \(t^n\) vanishes unless
\(j=(n+b-a)/m\) is an integer with \(0\le j\le n\), and \(n+c\)
is even. In that case it equals
\[
\frac{(n+b)!}{j!(n-j)!}(n+c-1)!! .
\]
This follows from the exact \(Z,W\) Wick contraction and the
ordinary \(Y\) moment. The two congruences give at most one
residue \(s\) modulo \(\ell\), \(0\le s<\ell\). Write \(n=\ell q+s\)
and \(j=dq+v\), \(v=(s+b-a)/m\). At all sufficiently large \(q\),
the coefficient is \(a_q R(q)\) for a rational function \(R\).
Its poles are simple and lie at \(-\beta-j'\) for \(\beta\in G_e\)
and a nonnegative integer \(j'\). Indeed, the shifted Gaussian
double factorial contributes a polynomial, as does the numerator
factor \((\ell q+s+b)!/(\ell q)!\). If the first denominator shift
\(v\) is positive, its factors \(dq+i\), \(1\le i\le v\), cancel
against the numerator factors \(\ell q+mi\), since
\(mv=s+b-a\le s+b\). A negative first shift contributes a
polynomial. Only positive shifts of the second denominator
factor \((e q+s-v)!\) can leave poles, giving exactly the stated grid.

Every rational multiplier having only simple poles on
\(-\beta-j'\) for \(\beta\) in the \(\mathcal Q\)-parameter set is
available in the unit cyclic space up to a finite Laurent polynomial.
In fact \(x^{-J}g\) has tail multiplier
\(C^J\prod_{h=0}^{J-1}\mathcal P(q+h)/\mathcal Q(q+h)=N_J/D_J\).
The two polynomials are coprime and \(D_J\) has only simple roots.
Choose \(J\) large enough that \(D_J\) clears the desired poles of
\(R\). By Bézout there are polynomials \(U,V_0\) with
\(R=U N_J/D_J+V_0\). Apply \(U(\vartheta)\) to \(x^{-J}g\) and
\(V_0(\vartheta)\) to \(g\). Any finite Laurent difference is in
\(K(t)\cdot1\). This already accounts for every even-\(m\) marker,
including its factor \(t^s\).

For odd \(m\), the only remaining grid is the half-integers
\(-1/2-j'\). The actual marker \(Z^{m(m-1)}\), with residue \(s=0\),
has tail multiplier
\[
\frac{\prod_{i=0}^{m-2}(2q-i)}
     {\prod_{i=1}^{m-1}(e q+i)} .
\]
Its pole at \(-1/2\) is simple with nonzero residue, and its
remaining poles are on the \(\mathcal Q\)-grid. Subtract the cyclic
part and the finite correction, then divide by that residue. This
produces exactly \(\kappa\) in the full actual space. The shift
\(x^{-J}\kappa\) has a nonzero simple pole at \(-1/2-J\) and no
other half-integer pole in its rational multiplier. All its other
poles are on the \(\mathcal Q\)-grid. Together with the cyclic
terms these shifts generate all marker tails, and finite Laurent
corrections belong to the constant line.

For independence, if a relation used \(\kappa\), clearing
denominators and choosing a nonzero residue class in \(t\) would
give an equation over \(K(x)\) with nonzero polynomial coefficient
on \(\kappa\). Write that coefficient as \(\sum_k A_k x^k\).
Divide the sufficiently large coefficient identities by \(a_{n-J}\)
where \(J\) is the largest shift in the whole equation. Each
\(A_k\ne0\) contributes a pole at \(n=k-1/2\) from the factor
\(1/(n-k+1/2)\). The factorial ratios have neither zeros nor poles
there, and no other term can have that pole. This is impossible.
The displayed bases and the generic rank assertions follow.

Moreover \((\vartheta+1/2)\kappa=g\). The quotient of the span of
\(\kappa\) and the unit cyclic space by the latter is a rank-one
trivial connection over \(K(t)\), with horizontal class
\(t^{\ell/2}\kappa\). This extension is nonsplit: a horizontal
lift would be a constant Laurent germ, while the span of \(\kappa\)
and the unit cyclic space has zero intersection with the constant
line by the independence just proved.

### The polynomial-time lattice

The generation proof used \(t\) and \(x\) monomial shifts, polynomial
Euler operators, and finite Laurent corrections. The differential
operator for \(g\) has leading Euler coefficient \(-Cx\), a unit in
\(R=K[t,t^{-1}]\), so every higher Euler derivative reduces into the
displayed basis over \(R\). The extra period \(\kappa\) for odd \(m\)
was obtained from an original marker by subtracting elements of
this same Laurent span of the constant and cyclic terms. Therefore
the basis elements belong to \(M[1/t]\), every original marker
belongs to their \(R\)-span, and \(M[1/t]=R\mathcal B\).

Origin saturation follows from the same backward shift as in the
[binomial full-image theorem](GAUSSIAN_BINOMIAL_MIXED_FULL_IMAGE.md#5-why-the-generic-frame-is-also-the-exact-global-r-frame).
Here is the argument over \(K\). Let \(X=(Z,W,Y)\) and let \(C\)
be the normalized Wick covariance. For a polynomial \(A_i\), put
\(T_i A_i=X_i A_i-\sum_j C_{ij}\partial_j A_i\). These polynomials
have expectation zero by Wick Stein. Conversely every centered
polynomial is a sum of terms \(T_i A_i\): subtract such terms to
remove the leading monomials, using their leading terms \(X_i A_i\),
and descend in degree. The remaining constant has expectation
zero. Stein with the exponential amplitude gives
\[
\Phi[T_i A_i]=t\sum_j C_{ij}\Phi[A_i\partial_j P_m].
\]
Thus \(h\mapsto (h-h(0)f)/t\) preserves original marker periods, and
it preserves the polynomial-time image because
\[
\frac{a h-a(0)h(0)f}{t}
 =a\frac{h-h(0)f}{t}+h(0)\frac{a-a(0)}t f
\qquad (a\in K[t]).
\]
If \(h\in R\mathcal B\) is regular as a formal germ at zero, some
\(t^N h\) belongs to \(M\). Iterated backward shifts divide by \(t\)
while the initial coefficient is zero, and give \(h\in M\).
Consequently \(M=(R\mathcal B)\cap K[[t]]\).

For explicit membership, write \(h=\sum_{b\in\mathcal B}r_b(t)b(t)\)
in its unique Laurent coordinates and take \(N\ge0\) such that
\(t^N r_b(t)\in K[t]\) for all \(b\). Membership in \(M\) is equivalent
to the finite equations
\[
[t^j]\sum_{b\in\mathcal B}\bigl(t^N r_b(t)\bigr)b(t)=0,
\qquad 0\le j<N .
\]
The coefficients of the basis series in these equations come from
the displayed factorial formula, Euler derivatives and \(\kappa\).
This is a finite principal-part test on the Laurent coordinates.

### Convergence and regularity in the full marked space

In every residue modulo \(\ell\), any fixed polynomial-time
combination has tail coefficients \(a_qR(q)\): shifts in that
residue change \(a_q\) only by rational ratios. The multiplier
\(R\) is a rational function. If nonzero, it has polynomial-size
asymptotics along the positive integers, whereas
\(a_{q+1}/a_q\asymp q^r\). That residue alone forces divergence.
Hence convergence of an element of \(M\) is equivalent to the
vanishing of all its rational tails, which gives a polynomial.
Clearing denominators gives the corresponding assertion in \(V\).
Every formal Laurent period in a globally regular differential
submodule converges at zero by the regular-singular Frobenius
theorem, as used in the [regular-core audit](GAUSSIAN_G_FUNCTIONS_AND_GLOBAL_REGULAR_CORE.md).
Thus the only such subspace is contained in the constant line,
and that line itself is globally regular.

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
