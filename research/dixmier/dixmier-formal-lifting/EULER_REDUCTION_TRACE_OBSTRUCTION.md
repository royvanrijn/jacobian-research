# Euler reduction of the cotangent example has a trace at every level

**QER1.** The retained rank-three cotangent endomorphism preserves the
weighted Euler moment exactly and descends to every quantum Hamiltonian
reduction at that moment. These reduced algebras have a normalized trace
at every scalar level. Consequently neither a reduction nor any finite
matrix algebra over it contains a Weyl pair with scalar nonzero
commutator. No full finite matrix corner is a Weyl algebra.

This prevents direct rank-three-to-rank-two descent through this Euler
reduction, including an alternative identification of the entire reduced
algebra with a Weyl algebra or a Morita equivalence. It does not exclude a
nonfull corner, an additional localization, or a different reduction.
The natural position chart gives a proper subalgebra of a rank-two Weyl
algebra; extending the cotangent map to that entire Weyl algebra requires
spatial denominators.

## Equivariance and exact quantum moment preservation

Use [QNI1's normalized retained Keller map](IDENTITY_LIFTS_AND_UNBOUNDED_INVERSES.md),
with \(q=1+xy\):
\[
 \begin{split}
 f_1&=(2x-3x^2y-x^3z)/2,\\
 f_2&=y+3xq^2z+3xy^2(4+3xy),\\
 f_3&=q^3z+y^2q(4+3xy).
 \end{split}                                                \tag{1}
\]
It is equivariant for weights \(w=(1,-1,-2)\); its three outputs
have the same respective weights. In the rank-three Weyl algebra
\(A_3(K)\), with \([p_i,x_j]=h\delta_{ij}\), \(h\in K^\times\),
write
\[
              E=xp_x-yp_y-2zp_z,\qquad
              B_w=A_3(K)^{\mathbb G_m}.
\]
The invariant algebra is the weight-zero part, equivalently the
centralizer of \(E\). It has central element \(E\).

For the exact cotangent endomorphism
\[
       \Phi(x_i)=f_i,\qquad
       \Phi(p_i)=\sum_j(Df)^{-1}_{ji}p_j,
\]
weighted Euler's identity gives
\[
 \sum_iw_if_i(Df)^{-1}_{ji}=w_jx_j,\qquad
                             \Phi(E)=E.                    \tag{2}
\]
There is no ordering correction: each coefficient is on the left
of \(p_j\). The map is equivariant, so induces an endomorphism
of every reduction
\[
                  H_\lambda=B_w/(E-\lambda)B_w,\qquad
                            \lambda\in K.                  \tag{3}
\]
The description equals
\((A_3/A_3(E-\lambda))^{\mathbb G_m}\): the torus weight
decomposition shows that the invariant part of the left ideal is
\(B_w(E-\lambda)\). No reductive-exactness inference beyond this
explicit grading is needed.

## Fourier transformation reduces the trace problem to positive weights

Fourier-transform the second and third pairs by
\(y\mapsto p_y,\ p_y\mapsto-y\) and
\(z\mapsto p_z,\ p_z\mapsto-z\). Then
\[
                    E\longmapsto E_++3h,\qquad
                    E_+=xp_x+yp_y+2zp_z.                  \tag{4}
\]
Thus (3) is isomorphic to the positive-weight invariant algebra
for weights \((1,1,2)\), reduced at \(E_+=\lambda-3h\).
It suffices to construct a normalized trace on
\[
 H^+_\eta=B_+/(E_+-\eta)B_+,\qquad
 B_+=A_3(K)^{(1,1,2)},\quad \eta\in K.                     \tag{5}
\]

## Polynomial continuation of finite-dimensional traces

Realize \(p_i=h\partial_{x_i}\). For each integer \(N\ge0\), the
polynomial space \(V_N\) of weighted degree \(N\), for weights
\((1,1,2)\), is finite dimensional and stable under \(B_+\).
Every element of \(B_+\) therefore has the ordinary trace
\(\operatorname{Tr}_{V_N}\).

For every fixed \(D\in B_+\), its trace is eventually a
quasipolynomial in \(N\), of period two. Here is a direct proof
covering every finite operator, without a degree cutoff.
In normal PBW order a monomial \(x^\alpha p^\beta\) has zero
trace unless \(\alpha=\beta\), since otherwise it shifts the
ordinary exponent multidegree of a basis monomial. For
\(\alpha=\beta\), its trace generating series is
\[
 \sum_{N\ge0}\operatorname{Tr}_{V_N}(x^\alpha p^\alpha)t^N
  =h^{|\alpha|}\alpha!\,
     \frac{t^{\alpha_1+\alpha_2+2\alpha_3}}
      {(1-t)^{\alpha_1+\alpha_2+2}(1-t^2)^{\alpha_3+1}}.    \tag{6}
\]
This follows by multiplying the three identities
\(\sum_{n\ge0}(n)_a t^{wn}=a!t^{wa}/(1-t^w)^{a+1}\).
The only poles in (6) are \(1,-1\). Partial fractions give
the claimed eventual quasipolynomial.

Let \(\mathcal T_0(D)(s)\) and \(\mathcal T_1(D)(s)\) be the
unique polynomials agreeing with these traces for all sufficiently
large even and odd \(N\), respectively. They depend \(K\)-linearly
on \(D\). For all \(U,V,D\in B_+\),
\[
 \mathcal T_\epsilon([U,V])=0,\qquad
 \mathcal T_\epsilon(E_+D)(s)=hs\,\mathcal T_\epsilon(D)(s).
                                                               \tag{7}
\]
Indeed both identities hold on every \(V_N\) of the specified
parity. Polynomial identity over a characteristic-zero field
then proves them universally.

Evaluation at \(s=\eta/h\) consequently gives a trace on (5).
No evaluation of an \(h\)-adic series at a nonzero parameter is used:
the functions in \(s\) here are finite polynomials.

For the identity, direct counting gives
\[
 \mathcal T_0(1)(s)=(s+2)^2/4,\qquad
 \mathcal T_1(1)(s)=(s+1)(s+3)/4.
\]
Their difference is the constant \(1/4\). Therefore
\[
 \boxed{\quad
 \tau_\eta(D)=4\bigl(\mathcal T_0(D)-\mathcal T_1(D)\bigr)(\eta/h),
 \qquad \tau_\eta(1)=1,\quad \tau_\eta([U,V])=0 .
 \quad}                                                     \tag{8}
\]
This works at every level \(\eta\) over every characteristic-zero
constant field containing \(h\ne0\), including Laurent-series fields.
Equation (4) transports it to every original \(H_\lambda\).

For example, with \(H_i=x_ip_i\), the trace satisfies
\[
 \tau_\eta(H_1)=\tau_\eta(H_2)=-h/2,\quad
 \tau_\eta(H_3)=(\eta+h)/2,\quad
 \tau_\eta(H_1^2)=0,\quad \tau_\eta(H_1H_2)=h^2/4.
\]
These are exact controls on the continuation; no positivity is asserted.

## Weyl pairs, finite matrices and full corners

A relation \([P,Q]=h\) in \(H_\lambda\) contradicts its normalized
trace. On \(\operatorname{Mat}_m(H_\lambda)\), use
\(\tau_\lambda\circ\operatorname{tr}_{\rm matrix}\).
It vanishes on commutators and takes value \(m\ne0\) on the identity.
Thus no relation \([P,Q]=hI_m\) exists there either.

There is a separate full-corner consequence, which does not assume
that the trace of an arbitrary idempotent is nonzero.
Let \(A=\operatorname{Mat}_m(H_\lambda)\), let \(e\) be a full
idempotent, and suppose \(eAe\) were a Weyl algebra of positive rank.
Every PBW element of a Weyl algebra is a commutator: integrate its
normal coefficient in the first position variable and commute
with the first momentum, dividing by \(h\).
Hence every trace on \(eAe\) is identically zero.

Fullness supplies a finite identity \(1=\sum_i a_i e b_i\).
For any \(z\in A\), cyclicity of its nonzero matrix trace gives
\[
 \tau_A(z)=\sum_i\tau_A(a_i e b_i z)
         =\sum_i\tau_A(e b_i z a_i e)=0,
\]
a contradiction. No full finite matrix corner can therefore be
a Weyl algebra; equivalently no Morita identification with a
Weyl algebra follows from this reduction.

For a nonfull corner the above fullness argument is unavailable.
In particular the proof does not claim that every corner has
nonzero trace on its identity, or exclude a mere Weyl subalgebra
in an arbitrary corner.

## The natural position quotient and its spatial boundary

On \(x\ne0\), take invariant coordinates \(u=xy\), \(v=x^2z\).
The localized reduction is \(D(\mathbb A^2_{u,v})\). Indeed the
source coordinates become
\[
 y=u/x,\quad z=v/x^2,\quad E=xp_x^{\rm chart},
\]
and the weight-zero localized algebra has an independent
rank-two Weyl algebra and the central Euler operator.

The global \(H_\lambda\) injects into this localized reduction.
One way to see faithfulness is to use differential order:
\[
 \operatorname{gr}\bigl(A_3/A_3(E-\lambda)\bigr)
    =K[x,y,z,\xi,\upsilon,\zeta]/(x\xi-y\upsilon-2z\zeta).
\]
The linear primitive polynomial in the denominator is irreducible,
so this associated graded ring is a domain. In particular the
module has no \(x\)-torsion and localization is injective.
Reduction by \(E-\lambda\) is strict for differential order because
its leading symbol is a nonzerodivisor, and taking the weight-zero
part commutes with associated graded.
Localization at \(x\) kills no nonzero leading class in this domain,
so the global-to-chart embedding is strict for differential order;
a chart-order-zero element therefore already belongs to
\(F_0H_\lambda=K[xy,x^2z]\).
Every invariant monomial acts on \(x^{\lambda/h}K[u,v]\) by a
finite differential operator polynomial in \(u,v\); the powers
of \(x\) cancel by its weight. Thus the global image lies in
the finite rank-two Weyl algebra. It is proper by (8), since
the entire Weyl algebra has no nonzero trace.

The induced position map is
\[
 \begin{split}
 A&=(2-3u-v)/2,\\
 B&=u+3v(1+u)^2+3u^2(4+3u),\\
 C&=v(1+u)^3+u^2(1+u)(4+3u),\\
                    (u,v)&\longmapsto(AB,A^2C).
 \end{split}                                                \tag{9}
\]
Its Jacobian determinant is \(A^2\), not a scalar unit.
This follows either by direct differentiation or from
\(\det Df=1\): the changes \((x,y,z)\mapsto(x,u,v)\)
have determinant \(x^3\), while \(x'=xA\).

In particular these two quotient positions have no polynomial
differential conjugates forming a full rank-two Weyl tuple.
For completeness, suppose commuting polynomial positions \(U,V\)
had arbitrary finite mates \(P_1,P_2\). The highest-order symbol
of a commutator with \(U\) or \(V\) is the contraction of the
momentum derivative of that symbol with \(dU\) or \(dV\).
Their nonzero Jacobian forces each mate to have order one;
its vector coefficients must be a column of the inverse
Jacobian. Those columns cannot both be polynomial when the
Jacobian determinant is the nonunit \(A^2\).
Order-zero terms and the Euler level cannot repair this
highest-order condition.

This failure to extend to the entire finite Weyl chart persists under
finite source and target Weyl automorphisms: these carry that Weyl
algebra onto itself, so they cannot turn an image outside it into an
element inside it. A different rational Darboux chart can change the
admitted algebra and needs its own global descent argument.

There is also a direct localization warning:
\[
                              \Phi(x)=xA .
\]
Extending the original cotangent map to the entire \(x\)-chart
requires \((xA)^{-1}=x^{-1}A^{-1}\).
The factor \(A^{-1}\) is a spatial denominator, not an Euler
scalar denominator removed by choosing a special level.

The induced endomorphism of \(H_\lambda\) itself is nonsurjective.
In the rational \(x\)-chart it is a first-order differential
substitution with position map (9) and leading momentum map its
inverse Jacobian; the Euler parameter changes only order-zero
connection terms. It preserves differential order, because the
position substitution is injective and its rational Jacobian is
invertible. If the endomorphism of \(H_\lambda\) were surjective,
preimages of \(u,v\) would have order zero. The global chart
image of every order-zero element is polynomial in \(u,v\),
and all such polynomials occur. Their preimages would therefore
supply a polynomial inverse to (9), contradicting its nonunit
Jacobian determinant. This is a nonsurjective endomorphism of
a reduced invariant algebra, not a Weyl endomorphism.

## Scope and next construction requirement

Euler reduction supplies an endomorphism of a quantum reduction
of classical dimension four but does not produce
an actual finite scalar Long lift or a finite Weyl carrier.
The normalized trace excludes replacing the entire reduction
by a Weyl algebra, including through a full matrix corner.
The natural Weyl chart exposes the precise missing global
admission through \(A^{-1}\) and the ramified quotient positions.

The next positive construction would need a different reduction,
a carefully justified nonfull corner or localization with finite
global descent, or another actual mixed finite pair.
None is constructed here. Formal completion may change trace
continuity and is outside the finite-algebra obstruction.
No DC1/DC2 resolution, bounded Long lift, formal verification
or external review is claimed.

