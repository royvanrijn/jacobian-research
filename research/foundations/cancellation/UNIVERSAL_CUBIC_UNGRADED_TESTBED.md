# A universal flat cubic ungraded testbed

This note constructs a large coefficient family of actual
geometric-degree-three Keller maps.  It contains the foundational map,
allows every denominator-free polynomial deformation of the quadratic and
cubic coefficients in the root-engineered chart, admits an arbitrary
polynomial `GL_2` Tschirnhausen gauge, and records a separate unramified
nonproperness divisor.  On this family all data in the cubic-normalization
frontend are explicit.

The result also identifies a sharp limitation.  The 24-dimensional
order-four kernel in
[`CUBIC_NORMALIZATION_FRONTEND.md`](CUBIC_NORMALIZATION_FRONTEND.md)
lives on the nonfree reduced Koszul trace module.  The family below has a
finite free normalization.  No coefficient deformation or polynomial
Tschirnhausen gauge inside this family can therefore realize one of those
24 directions.  A single testbed containing both sectors would have to
deform the underlying normalization module, not only its binary-cubic
coefficients.  No compatible distinguished Keller open is currently known
in that nonflat sector.

Work over a characteristic-zero field `k`.  The normality and geometric
boundary statements are made after passing to an algebraic closure when
necessary.

## 1. The coefficient family

Choose arbitrary polynomials

\[
 \alpha(P),\gamma(P)\in k[P]
\]

and put

\[
 A(P)=\alpha(P),\qquad
 h(P)=1+P^2\gamma(P),\qquad
 a(P)=P h(P).                                      \tag{1.1}
\]

For source coordinates `(x,y,z)`, define

\[
 t=1+xy,\qquad
 q=t^2z+y^2(1+3t),\qquad P=tq.                     \tag{1.2}
\]

The testbed map is

\[
\boxed{
\begin{aligned}
 F_1={}&P,\\
 F_2={}&y+3xq+2\alpha(P)
              +3t^2xq^3\gamma(P),\\
 F_3={}&x(5-3t)-x^3z-x^3q^3\gamma(P).
\end{aligned}}
                                                               \tag{1.3}
\]

At `alpha=gamma=0`, this is the foundational root-engineered
quadratic-gauge presentation.  It is linearly left--right equivalent to the
foundational map in the main README.

The bounded quartic coefficient cell is already nontrivial:

\[
 \alpha=\alpha _0+\alpha _1P+\alpha _2P^2
              +\alpha _3P^3+\alpha _4P^4,\qquad
 \gamma=\gamma _0+\gamma _1P.                         \tag{1.4}
\]

It has seven independent coefficient parameters, contains the foundational
point, and gives `deg A<=4`, `deg a<=4`.  Formula (1.3) is valid without
these degree bounds.

The actual normalization testbed used below is the squarefree locus
`gcd(a,a')=1`.  Formula (1.3) remains a Keller map off that locus, but the
finite Deligne--Faddeev order can then be nonnormal and is not asserted to
be the canonical normalization without a further integral-closure
calculation.

### Theorem 1.1 -- universal denominator-free cubic lift

Every map (1.3) satisfies

\[
 \boxed{\det DF=-2}
                                                               \tag{1.5}
\]

and has geometric degree three.

#### Proof

On `t!=0`, put

\[
 S={x\over t},\qquad Q=y+xq,\qquad
 D=1-S(Q-PS)={1\over t}.                               \tag{1.6}
\]

For target coordinates `(P,B,C)`, define

\[
 E(P,B,C;S)
 =a(P)S^3+\left(A(P)-{B\over2}\right)S^2+S-{C\over2}.
                                                               \tag{1.7}
\]

The short-chart form of (1.3) is

\[
\begin{aligned}
 B&=Q+2A(P)+(3a(P)-P)S,\\
 C&=2\bigl(a(P)S^3+A(P)S^2+S\bigr)-BS^2.             \tag{1.8}
\end{aligned}
\]

The quadratic coefficient `A(P)` cancels from the intercept in (1.8), so
it is completely arbitrary.  The divisibility
`a(P)-P in (P^3)` is exactly what makes the remaining cubic correction
polynomial in `(x,y,z)`.  Conversely, it is forced term by term in this
coefficient-only lift ansatz: an added `P^mS^3` term pulls back through
`S=x/t`, so it requires `m>=3`.

Equations (1.7)--(1.8) give

\[
 E(P,B,C;S)=0,\qquad
 \partial_SE=D                                      \tag{1.9}
\]

on the source incidence.  At fixed `P`,

\[
 \det{\partial(B,C)\over\partial(S,Q)}=-2D,
\qquad
 \det{\partial(P,S,Q)\over\partial(x,y,z)}=D^{-1}.
                                                               \tag{1.10}
\]

Their product proves (1.5).  The polynomial (1.7) is a generic irreducible
cubic: over `k(P,B)`, it is a nonconstant polynomial in `S` minus the
independent parameter `C/2`.  Every generic simple root reconstructs by
(1.6), proving geometric degree three.  QED

The family in
[`UNIVERSAL_CUBIC_GAUGE_MULTIPLICITY.md`](../verified/UNIVERSAL_CUBIC_GAUGE_MULTIPLICITY.md)
is the specialization

\[
 \alpha=0,\qquad \gamma=P^{n-3}-1.                    \tag{1.11}
\]

Thus (1.3) is the coefficient-universal version of that construction, not
a second mechanism.

## 2. Exact finite normalization

Put

\[
 R=k[P,B,C],\qquad
 b=A(P)-{B\over2},\qquad d=-{C\over2}.                \tag{2.1}
\]

The intrinsic binary cubic in the displayed adapted gauge is

\[
\boxed{
 f_{P,B,C}(U,V)
 =a(P)U^3+bU^2V+UV^2+dV^3.
}
                                                               \tag{2.2}
\]

Its marked-root incidence

\[
 \bar X=
 \{((P,B,C),[U:V]):f_{P,B,C}(U,V)=0\}
 \subset\mathbb A^3\times\mathbb P^1                 \tag{2.3}
\]

is finite of degree three over the target.  Its affine coordinate algebra
is the Deligne--Faddeev algebra

\[
 \mathcal B=R\oplus Ru\oplus Rv                       \tag{2.4}
\]

with multiplication

\[
\boxed{
\begin{aligned}
 u^2&=-a-bu+av,\\
 uv&=-ad={aC\over2},\\
 v^2&=-bd-du+v={bC\over2}+{C\over2}u+v .
\end{aligned}}
                                                               \tag{2.5}
\]

Indeed, on the chart `V!=0`, if `S=U/V`, then

\[
 u=aS,\qquad v=aS^2+bS+1,
                                                               \tag{2.6}
\]

and reduction modulo (1.7) gives (2.5).

Assume from now on that

\[
 \boxed{a(P)=P h(P)\text{ is squarefree}.}             \tag{2.7}
\]

This includes the foundational point `h=1`.  The finite incidence (2.3) is
then smooth.  On `V!=0`, its equation has
`partial E/partial C=-1/2`.  On the infinity chart `U!=0`, with
`T=V/U`, its equation is

\[
 a+bT+T^2+dT^3=0.                                     \tag{2.8}
\]

At `T=0`, smoothness is exactly `a'(P)!=0`.  Condition (2.7) supplies it.
Hence `mathcal B` is normal and is the canonical finite normalization of
(1.3), not merely a finite cubic order.

In particular,

\[
 \operatorname{Fitt}_3^R(\mathcal B)=R.               \tag{2.9}
\]

There is no closed-point flatness defect anywhere in this coefficient
cell.

## 3. Branch and nonproperness divisors

The branch equation is the binary-cubic discriminant

\[
\boxed{
\Delta
 =b^2-4a+2b^3C-9abC-{27\over4}a^2C^2.
}
                                                               \tag{3.1}
\]

It is irreducible in `k[P,B,C]`: over `k(P)`, it is the ordinary
irreducible discriminant of the cubic whose last two coefficients vary
independently, and it is primitive over `k[P]`.

The ramification surface has the explicit normalization

\[
\boxed{
\begin{aligned}
 B&=2A(P)+3a(P)S+S^{-1},\\
 C&=S-a(P)S^3,
\end{aligned}
\qquad (P,S)\in\mathbb A^1\times\mathbb G_m .
}
                                                               \tag{3.2}
\]

Substitution gives `E=partial_SE=Delta=0`.  In particular, the
ramification support is the smooth Laurent plane.

The reduced nonproperness equation is

\[
\boxed{
 j_F=\Delta\cdot\operatorname{rad}(h).
}
                                                               \tag{3.3}
\]

Here and below equations are understood up to a nonzero scalar.  At a
simple nonzero root `rho` of `h`, the Newton polygon of (1.7) over
`P=rho` has blocks

\[
 (0,0)\longrightarrow(2,0)\longrightarrow(3,1).       \tag{3.4}
\]

It gives two finite affine roots and one escaping unramified root with
label `(e,f)=(1,1)`.  The entire hypersurface `V(h)` is therefore the
candidate phantom boundary.  At `P=0`, the two `q=0` sheets and the
`t=0` sheet are all affine, exactly as in the fiber-invisible cubic
calculation; `P=0` is not a boundary component.  Away from
`V(Delta h)`, (1.6) reconstructs every root regularly, so (3.3) exhausts
the divisorial nonproperness locus.

Consequently

\[
 u_F={j_F\over\Delta}=\operatorname{rad}(h).           \tag{3.5}
\]

For squarefree `h`, this is simply `h`.

## 4. The two closed-point saturation modules

Let

\[
 Q_{\mathrm{cot}}=\Omega_{\mathcal B/R},\qquad
 T=\mathcal B/\operatorname{Ann}_{\mathcal B}(Q_{\mathrm{cot}}).
                                                               \tag{4.1}
\]

All ramification lies in the finite-root chart, because
`partial_SE(0)=1`.  On that chart,

\[
 Q_{\mathrm{cot}}
 =\mathcal B/(\partial_SE)\,dS.                        \tag{4.2}
\]

The total algebra is a domain and `partial_SE` is nonzero, so

\[
 \boxed{
 T=\mathcal B/(\partial_SE),\qquad
 Q_{\mathrm{cot}}=T\,dS.
}
                                                               \tag{4.3}
\]

Formula (3.2) identifies

\[
 T\simeq k[P,S,S^{-1}].                                \tag{4.4}
\]

Thus `T` is already its `S_2` hull, `Q_cot` has no closed-point torsion,
and the two saturation quotients are

\[
 C/T=0,\qquad H_Z^0(Q_{\mathrm{cot}})=0.               \tag{4.5}
\]

Equivalently, the two frontend obstruction modules vanish:

\[
\boxed{
 \operatorname{Ext}_R^2(T,R)=0,\qquad
 \operatorname{Ext}_R^3(Q_{\mathrm{cot}},R)=0.
}
                                                               \tag{4.6}
\]

This calculation holds for every member satisfying (2.7), including every
collision on the triple-root curve.

## 5. Intrinsic coefficient map and general Tschirnhausen gauge

In the adapted gauge, the coefficient morphism is

\[
 \kappa_{\alpha,\gamma}(P,B,C)
 =\left(a(P),\,A(P)-{B\over2},\,1,\,-{C\over2}\right)
 \in\operatorname{Sym}^3(k^2).                        \tag{5.1}
\]

Let

\[
 \Gamma(P,B,C)=
 \begin{pmatrix}r&s\\ \ell&u\end{pmatrix}
 \in GL_2(R).
                                                               \tag{5.2}
\]

Since the only units of `R` are scalars,
`ru-s ell in k^*`.  The completely general polynomial Tschirnhausen
representative is

\[
 f^\Gamma(U,V)
 =f(rU+sV,\ell U+uV).                                  \tag{5.3}
\]

Writing `f=aU^3+bU^2V+cUV^2+dV^3`, its four coefficients are

\[
\begin{aligned}
 a^\Gamma={}&ar^3+br^2\ell+cr\ell^2+d\ell^3,\\
 b^\Gamma={}&3ar^2s+b(r^2u+2rs\ell)
              +c(2r\ell u+s\ell^2)+3d\ell^2u,\\
 c^\Gamma={}&3ars^2+b(2rsu+s^2\ell)
              +c(ru^2+2s\ell u)+3d\ell u^2,\\
 d^\Gamma={}&as^3+bs^2u+csu^2+du^3.
                                                               \tag{5.4}
\end{aligned}
\]

The optional Deligne--Faddeev determinant twist only multiplies all four
coefficients by a scalar unit here.  In either convention,

\[
 \operatorname{Disc}(f^\Gamma)
 =\det(\Gamma)^6\operatorname{Disc}(f).                \tag{5.5}
\]

Therefore (2.3)--(4.6), the branch and nonproperness divisors, and the
verdict below are independent of the displayed polynomial gauge.  Raw
coefficient degree or affine span in (5.4) is not intrinsic.

## 6. Exact answer for the base-change map (G)

The adapted coefficient hyperplane is `c=1`.  After deleting that constant
coordinate, Proposition 2.10 of
[`CUBIC_GAUGE_STRAIGHTENING.md`](CUBIC_GAUGE_STRAIGHTENING.md)
uses the base-change map

\[
\boxed{
 G_{\alpha,\gamma}(P,B,C)
 =\left(
   a(P),\ A(P)-{B\over2},\ -{C\over2}
  \right).
}
                                                               \tag{6.1}
\]

It is triangular, with

\[
 \det DG_{\alpha,\gamma}={a'(P)\over4}.                \tag{6.2}
\]

A polynomial self-map of the affine line is an automorphism exactly when
it has degree one.  Since `a(P)=P(1+P^2 gamma(P))`,

\[
\boxed{
\begin{aligned}
 G_{\alpha,\gamma}\text{ is a polynomial automorphism}
 &\Longleftrightarrow \gamma=0\\
 &\Longleftrightarrow h=1\\
 &\Longleftrightarrow u_F\in k^*.
\end{aligned}}
                                                               \tag{6.3}
\]

There is no restriction on `alpha`: when `gamma=0`, the inverse is

\[
 P=a,\qquad
 B=2(A(P)-b),\qquad C=-2d.                             \tag{6.4}
\]

Thus this large cell proves the precise coupling anticipated in the cubic
closure protocol:

> **Flat-cell base-change theorem.**  In the universal denominator-free
> cubic lift family, the base-change map `(G)` is forced to be an
> automorphism exactly when the phantom unramified divisor is absent.

On the boundary-minimal subcell, Certificate P therefore forces `(G)` to
be an automorphism and every map is foundational up to polynomial
left--right equivalence.  Outside that subcell, `(G)` fails for exactly the
phantom-boundary reason: `V(h)` is a distinct unramified nonproperness
divisor.  Theorem 8.1 below shows that two such components can certify
genuine ungradedness.

## 7. Why the 24-dimensional quartic kernel is a different sector

The order-four kernel in the cubic-normalization frontend starts with

\[
 M_K=\operatorname{coker}\left(
 R\mathop{\longrightarrow}^{(z,-y,x)^{\mathsf T}}R^3
 \right),\qquad \mathcal B_K=R\oplus M_K.              \tag{7.1}
\]

Its underlying module has

\[
 \operatorname{Fitt}_3^R(\mathcal B_K)=(x,y,z),        \tag{7.2}
\]

whereas every member of the actual Keller family (1.3) has the unit ideal
(2.9).  The 24 order-four tensors deform multiplication on `M_K`; they do
not deform a binary cubic on the free module in (2.4).  Fitting ideals of
the underlying normalization module are invariant under coefficient
changes and polynomial `GL_2` gauge.  Hence:

\[
\boxed{
\text{No specialization of }\alpha,\gamma,\Gamma
\text{ realizes a nonzero Koszul quartic-kernel direction.}
}
                                                               \tag{7.3}
\]

This is a separation theorem, not an exclusion of the Koszul sector from
all Keller maps.  Existing exact computations show:

- every one of the 24 basis-axis families over each squarefree cubic symbol
  retains saturated cotangent presentation and a length-six support-hull
  defect;
- all squarefree coordinate planes do the same;
- all smooth coordinate three-spaces do the same.
- one sum/alternating-sum parameter plane does the same for every
  squarefree symbol, although its generic tensor has nonzero coordinates in
  all 24 fixed kernel directions.
- after translation by the deterministic generic quartic lift, the same
  affine plane has pure support and the constant length-six defect for all
  ten cubic-symbol strata, including the double-line, triple-line, and zero
  symbols.
- the six complete singular-squarefree quartic nongauge complements are
  cotangent-saturated on every geometric parameter fiber, and their
  intrinsic differents require six local generators;
- every compatible formal tail with squarefree leading symbol remains
  cotangent-saturated with a six-generated non-Cartier different, by the
  all-orders formal normal-form and strict-Rees theorem `SSADPALL`.

For squarefree symbols, arbitrary compatible formal tails are therefore no
longer open at the cotangent-saturation or different-generator level.  What
remains open is normality and algebraization of nonhomogeneous lifts, the
double-line/triple-line/zero rows behind their generic-etaleness gate, and
the existence of a distinguished `A^3` Keller open.  Until those bridges are
constructed, calling the Koszul tensor family a family of Keller maps would
be incorrect.

## 8. An explicit genuinely ungraded cubic

The phantom divisor does more than obstruct the boundary-minimal
classification.  With two distinct phantom components it kills every
connected algebraic-torus symmetry.

Specialize

\[
 \alpha(P)=0,\qquad \gamma(P)=1,\qquad
 h(P)=1+P^2,\qquad a(P)=P(1+P^2).                    \tag{8.1}
\]

Thus

\[
\boxed{
\begin{aligned}
 t&=1+xy,\qquad q=t^2z+y^2(1+3t),\qquad P=tq,\\
 F_1&=P,\\
 F_2&=y+3xq+3t^2xq^3,\\
 F_3&=x(5-3t)-x^3z-x^3q^3.
\end{aligned}}
                                                               \tag{8.2}
\]

Theorem 1.1 gives `det DF=-2` and geometric degree three.  A scalar
target change normalizes the determinant to one without changing any
claim below.

> **Theorem 8.1 -- cubic algebraic-torus exclusion.**
> Over an algebraically closed characteristic-zero field, (8.2) admits no
> nontrivial algebraic \(\mathbb G_m\) source--target equivariance.
> Consequently no polynomially left--right equivalent representative is
> graded for a nontrivial weight signature.  Hence the minimum geometric
> degree of a genuinely ungraded noninvertible Keller map is exactly three.

Here is the intrinsic proof.  The two roots
\(\rho_+,\rho_-\) of \(h=P^2+1\) give two distinct unramified canonical
boundary target components \(P=\rho_\pm\).  On the normalization of the
prime ramified discriminant, use the coordinates of (3.2):

\[
 \widetilde Z_\Delta
 =\operatorname{Spec}k[P,S,S^{-1}].
                                                               \tag{8.3}
\]

The intersections with the two unramified boundary components are the
divisors \(D_\pm=(P-\rho_\pm)\).  The relative-differential Fitting divisor
of the normalization map is

\[
\boxed{J=1-3P(1+P^2)S^2.}                              \tag{8.4}
\]

Indeed, at fixed \(P\),

\[
 \partial_S B=-S^{-2}J,\qquad \partial_S C=J,
                                                               \tag{8.5}
\]

so
\(\Omega_{\widetilde Z_\Delta/Z_\Delta}
 \simeq k[P,S,S^{-1}]/(J)\,dS\).

Let a connected algebraic torus act on the source and target making
\(F\) equivariant.  Functoriality of the canonical normalization-boundary
package lifts the target action to (8.3), preserving
\(\{D_+,D_-\}\) and \(V(J)\).  Connectedness prevents permutation of
\(D_+\) and \(D_-\).  Delete these two divisors.  The unit group modulo
scalars is freely generated by

\[
 S,\qquad P-\rho_+,\qquad P-\rho_-.
                                                               \tag{8.6}
\]

A connected group acts trivially on this discrete lattice.  Hence, for
characters \(\chi_\pm,\chi_S\),

\[
\begin{aligned}
 g^*(P-\rho_\pm)&=\chi_\pm(g)(P-\rho_\pm),\\
 g^*S&=\chi_S(g)S.
\end{aligned}                                                \tag{8.7}
\]

The two expressions for \(g^*P\) and
\(\rho_+\ne\rho_-\) force
\(\chi_+=\chi_-=1\), so \(g^*P=P\).  Preservation of \(V(J)\) gives

\[
 1-3P(1+P^2)\chi_S(g)^2S^2
 =u_g\bigl(1-3P(1+P^2)S^2\bigr),                    \tag{8.8}
\]

where \(u_g\) is a unit on the deleted-divisor chart.  Comparing the two
\(S\)-exponents first makes the \(S\)-exponent of \(u_g\) zero.  Evaluating
at \(P=\rho_+\) and \(P=\rho_-\) makes its two remaining boundary
exponents zero.  The constant term then gives \(u_g=1\), and
\(\chi_S(g)^2=1\).  Connectedness forces \(\chi_S=1\).  Thus the torus acts
trivially on the decorated ramified normalization.

The discriminant hypersurface is prime and nonnormal.  Primality is the
irreducibility from Section 3, while (8.4) is a nonunit and therefore the
finite birational normalization has nonzero relative differentials.
The pointwise-fixed hypersurface lemma from
[`NO_ALGEBRAIC_TORUS_EQUIVARIANCE.md`](NO_ALGEBRAIC_TORUS_EQUIVARIANCE.md)
now forces the target action on \(\mathbb A^3\) to be trivial.  The source
action is then a connected family of automorphisms of the finite separable
generic cubic extension over the target function field.  That deck group
is finite, so the connected source action is trivial as well.

Polynomial left--right changes merely conjugate the two actions.  Thus the
exclusion applies to every representative of the polynomial left--right
class.  Finally, noninvertible Keller maps of geometric degree one or two
do not exist by the exact geometric-degree spectrum.  The upper example
(8.2) and this lower bound prove the asserted minimum.

This theorem resolves the former `OP-UG3` conjectural direction negatively:
not every geometric-degree-three Keller map is graded.  It does not classify
all cubic classes, compute every discrete self-equivalence of (8.2), or give
a stable symmetry-exclusion statement after adjoining identity variables.

## 9. Consequences for the cubic classification program

The testbed closes all requested invariants on a large actual-Keller cell:

| datum | value |
|---|---|
| finite normalization | the finite free Deligne--Faddeev algebra (2.4)--(2.5) |
| branch divisor | `V(Delta)` from (3.1) |
| nonproperness divisor | `V(Delta h)` |
| support-hull module | `Ext_R^2(T,R)=0` |
| cotangent point module | `Ext_R^3(Q_cot,R)=0` |
| intrinsic coefficient map | the polynomial `GL_2(R)` orbit of (5.1) |
| base-change map `(G)` | automorphism iff `h=1`, equivalently iff the phantom factor is a unit |

For `OP-SUSP`, this proves foundational base-change rigidity throughout the
flat denominator-free cell once boundary minimality removes `V(h)`.
For the former `OP-UG3` direction, Theorem 8.1 shows that the extra
unramified divisor can instead certify genuine ungradedness.  Thus the
arbitrary-cubic classification cannot aim to prove that every cubic is
graded.  The surviving useful target is the boundary-minimal cubic
classification, where the phantom factor is absent.

It does not prove that every cubic Keller map enters this cell.  The
remaining excluded possibility is exactly the one the frontend already
identifies: a module-level degeneration of the finite normalization,
potentially carrying an arbitrary combination of the 24 quartic-kernel
tensors, together with a compatible `A^3` Keller open.

## 10. Exact verification

Run

```bash
.venv/bin/python scripts/verify_universal_cubic_ungraded_testbed.py
```

The checker verifies the full seven-parameter quartic coefficient cell by a
direct three-variable Jacobian expansion, the inverse and derivative
identities, the reciprocal chart, the Deligne--Faddeev multiplication
table, the discriminant and Laurent ramification parametrization, the
universal `GL_2` discriminant transformation, the triangular base-change
criterion, and the Fitting-ideal separation from the reduced Koszul cell.
For (8.2) it also verifies squarefreeness of \(P(1+P^2)\), the exact Fitting
generator (8.4), the two phantom-boundary character equations, and the
remaining order-two Fitting character.
