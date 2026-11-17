# Free discriminants and Saito matrices: first experiment

This note tests whether free divisors can organize the marked-root determinant
ledgers.  It records an exact calculation, one corrected general proposition,
and a candidate search.  It does **not** claim a new Keller family or a
classification theorem.

Work over a characteristic-zero field.  For an inverse equation

\[
 E(T;y)=0
\]

write \(\Gamma\) for its reduced branch discriminant.  The object relevant to
the existing determinant ledgers is sometimes larger:

\[
 \Delta_{\mathrm{led}}
 =\Gamma\prod_i L_i,
\]

where the \(L_i=0\) are divisorial reconstruction boundaries.  They must be
distinguished from powers in the raw polynomial discriminant.  In particular,
a degree-drop factor can record roots at infinity without being a ramified
component of the affine marked-root chart.

All calculations below are reproduced by

```bash
.venv/bin/python scripts/verify_free_discriminant_saito.py
Singular -q scripts/verify_free_discriminant_saito_nonfree.sing
```

The second command takes about one minute on the reference machine.

## 1. Summary

The first experiment gives a mixed answer.

| construction | raw polynomial discriminant, up to a unit | divisor tested | result |
|---|---|---|---|
| foundational cubic | \(\Gamma_3\) | \(\Gamma_3\) | free |
| canonical weighted quartic | \(C^3\Gamma_{\mathrm w}\) | \(C\Gamma_{\mathrm w}\) | free |
| cancellation \((m,r)=(2,1)\) | \(P^2\Gamma_{\mathrm c}\) | \(P\Gamma_{\mathrm c}\) | free |
| first quadratic-gauge quartic | \(P^2\Gamma_{\mathrm q}\) | \(\Gamma_{\mathrm q}\) and \(P\Gamma_{\mathrm q}\) | neither is free |

The branch divisors \(\Gamma_{\mathrm w}\) and
\(\Gamma_{\mathrm c}\) alone are not free.  Adjoining exactly the second
boundary component used by the determinant ledger makes them free.  This is
the strongest positive evidence for the proposed viewpoint.

The quadratic-gauge row is the obstruction.  Its fixed-\(P\) plane
discriminant is free, as every reduced affine plane curve is, but the full
three-dimensional branch divisor and the full reduced ledger divisor are not
free.  Since the quadratic-gauge quartic is nevertheless a polynomial Keller
map, full-target freeness is not necessary for the present construction
engine.

Here “minimal Saito matrix” means a basis with the rank-minimal number of
columns.  No degree-minimality assertion is made for these nonhomogeneous
matrices.

## 2. Foundational cubic

Use the normalized quadratic-gauge inverse equation

\[
 E_3(T)=PT^3-\frac B2T^2+T-\frac C2.
\]

Its discriminant is \(-\Gamma_3/4\), where

\[
 \Gamma_3
 =B^3C-B^2-18BCP+27C^2P^2+16P.
\]

A Saito matrix, with columns interpreted as logarithmic derivations in the
order \((P,B,C)\), is

\[
 S_3=
 \begin{pmatrix}
 2P&-2B&0\\
 B&3BC-16&B^2-12P\\
 -C&9C^2&3BC-4
 \end{pmatrix},
 \qquad
 \det S_3=8\Gamma_3.
\]

Consequently

\[
 \operatorname{Der}(-\log\Gamma_3)
 =k[P,B,C]\delta_1\oplus
  k[P,B,C]\delta_2\oplus
  k[P,B,C]\delta_3,
\]

where the \(\delta_i\) are the columns of \(S_3\).

Each column lifts regularly to the marked-root incidence.  The corresponding
root velocities are

\[
 \tau_1=-T,\qquad
 \tau_2=T(3C+4T),\qquad
 \tau_3=BT-2.
\]

Thus

\[
 \delta_i(E_3)+\tau_i\partial_TE_3\in(E_3).
\]

The first column is the weighted scaling direction and moves the marked root
by \(-T\).  The transverse reconstruction direction is \(\partial_C\):

\[
 \partial_CE_3=-\frac12,\qquad
 \partial_C(T)=\frac1{2\partial_TE_3}
\quad\text{on }E_3=0.
\]

This is the reciprocal marked-root derivative in the usual reconstruction
formula.  It is deliberately not logarithmic: it crosses the discriminant.

## 3. First weighted quartic

Take \(H(T)=T^3(1-T)\).  The inverse equation is

\[
 E_{\mathrm w}(T)=T^3(1-T)-BCT+AC^2.
\]

Its raw discriminant is \(-C^3\Gamma_{\mathrm w}\), with

\[
 \begin{aligned}
 \Gamma_{\mathrm w}={}&
 256A^3C^3-192A^2BC^2+27A^2C\\
 &-6AB^2C+27B^4C-4B^3.
 \end{aligned}
\]

The branch divisor \(\Gamma_{\mathrm w}=0\) alone is not free.  The reduced
ledger divisor

\[
 \Delta_{\mathrm w}=C\Gamma_{\mathrm w}
\]

is free, with

\[
 S_{\mathrm w}=
 \begin{pmatrix}
 2A&-2B^2&-14B^2\\
 B&8B^2C-8AC-B&128ABC^2-64AC+B\\
 -C&16BC^2-3C&256AC^3-48BC^2+3C
 \end{pmatrix}
\]

and

\[
 \det S_{\mathrm w}=-16\Delta_{\mathrm w}.
\]

The first logarithmic field

\[
 2A\partial_A+B\partial_B-C\partial_C
\]

fixes the marked root: its lift has \(\tau_1=0\).  The other two columns also
have polynomial root velocities; the exact expressions are checked by the
script.

The marked-root derivative and reconstruction pole are

\[
 \partial_TE_{\mathrm w}=3T^2-4T^3-BC,
\qquad
 \partial_A(T)=-\frac{C^2}{\partial_TE_{\mathrm w}}.
\]

After passing from \(A\) to the plane-incidence coordinate \(t=AC^2\), this
is exactly \(-1/\partial_TE_{\mathrm w}\); the weighted source formula
\(x=-C/\partial_TE_{\mathrm w}\) supplies the corresponding scaled
reconstruction direction.

The new structural point is that the extra component \(C=0\) is not noise:
it is what completes the nonfree branch surface to the free ledger divisor.

## 4. First cancellation quartic

For \((m,r)=(2,1)\), set the harmless scalar \(C=1\).  The inverse equation is

\[
 E_{\mathrm c}(T)
 =T-\frac{Q^2T^2}{2}+\frac{2PQT^3}{3}
  -\frac{P^2T^4}{4}-R.
\]

Its raw discriminant is \(P^2\Gamma_{\mathrm c}/432\), where

\[
 \begin{aligned}
 \Gamma_{\mathrm c}={}&
 1728P^4R^3-3456P^3QR^2+288P^2Q^4R^2\\
 &+1656P^2Q^2R-729P^2-288PQ^5R+136PQ^3\\
 &+12Q^8R-6Q^6.
 \end{aligned}
\]

Here \(P=0\) is the second reconstruction-boundary image, not the generically
ramified branch surface.  Again the branch surface alone is not free, while

\[
 \Delta_{\mathrm c}=P\Gamma_{\mathrm c}
\]

is free.  Put

\[
 \begin{aligned}
 a={}&-384P^2QR^2-8PQ^2R-3Q^3+54P,\\
 b={}&144Q^4R^2+192P^2R^3-608PQR^2+24Q^2R-45,\\
 c={}&-96P^2Q^2R-2PQ^3+81P^2,\\
 d={}&36Q^5R+48P^2QR^2-152PQ^2R-18Q^3+54P.
 \end{aligned}
\]

Then

\[
 S_{\mathrm c}=
 \begin{pmatrix}
 3P&0&0\\
 Q&a&c\\
 -2R&b&d
 \end{pmatrix},
\qquad
 \det S_{\mathrm c}=-27\Delta_{\mathrm c}.
\]

The Euler column lifts with root velocity \(\tau_1=-2T\).  The remaining two
root velocities are polynomial and are verified exactly.

The controlled divisor is

\[
 \partial_TE_{\mathrm c}
 =1-T(Q-PT)^2=D,
\]

and the transverse target direction gives

\[
 \partial_RE_{\mathrm c}=-1,\qquad
 \partial_R(T)=D^{-1}.
\]

Thus the Saito frame contains the regular tangent motions, while
\(\partial_R\) is the distinguished transverse reconstruction direction
whose root lift has the reciprocal pole cancelled by the source chart.

## 5. First quadratic-gauge quartic

Use the repository's small seed

\[
 G(T)=T(T-1)(T+1)(T-2)
 =T^4-2T^3-T^2+2T.
\]

The full inverse equation is

\[
 E_{\mathrm q}(T)
 =P^4T^4-2PT^3-(P+B)T^2+2T-C.
\]

Its raw discriminant has the form

\[
 \operatorname{Disc}_T(E_{\mathrm q})=-16P^2\Gamma_{\mathrm q}.
\]

The expanded \(\Gamma_{\mathrm q}\) is recorded in both verification scripts.
The exact Jacobian-ideal calculation gives

\[
 \operatorname{pd}(\Gamma_{\mathrm q},
 \partial_P\Gamma_{\mathrm q},
 \partial_B\Gamma_{\mathrm q},
 \partial_C\Gamma_{\mathrm q})=2,
\]

and the same result for \(P\Gamma_{\mathrm q}\).  Equivalently, each
codimension-two Jacobian ideal has a minimal resolution one step longer than
the Hilbert--Burch resolution required for a free surface divisor.  Hence

\[
 \boxed{\Gamma_{\mathrm q}\text{ and }P\Gamma_{\mathrm q}
 \text{ are not free divisors}.}
\]

There is therefore no full-target Saito matrix to find.

On the canonical section \(P=1\), put \(x=B,\ y=C\).  The reduced plane
discriminant is

\[
 \begin{aligned}
 f_{\mathrm q}={}&x^4y+5x^3y-x^3+8x^2y^2+29x^2y-4x^2\\
 &+52xy^2+29xy-23x+16y^3+23y^2-2y-9.
 \end{aligned}
\]

It is free.  An exact two-column Saito basis is

\[
 S_{\mathrm q,1}=[\,v_1-8v_2,\ 11v_0+2v_2\,],
\qquad
 \det S_{\mathrm q,1}=-1500f_{\mathrm q},
\]

where the three shorter syzygy columns \(v_0,v_1,v_2\) are written explicitly
in `scripts/verify_free_discriminant_saito.py`.  Both columns lift regularly
to the root incidence.

The reconstruction direction remains completely transparent:

\[
 \partial_CE_{\mathrm q}=-1,\qquad
 \partial_C(T)=\frac1{\partial_TE_{\mathrm q}},
\qquad
 \partial_TE_{\mathrm q}=2D.
\]

Thus the determinant ledger survives even though full-target freeness does
not.

## 6. The exact Saito--incidence statement

The calculation supports the following proposition.

### Proposition 6.1 -- logarithmic frame and transverse root pole

Let \(Y=\mathbb A^n\), let \(E(T;y)\) define a finite generically separable
marked-root incidence \(X\to Y\), and let \(\Delta=0\) be a reduced free
divisor containing its branch divisor.  Let

\[
 S=(\delta_1,\ldots,\delta_n),\qquad \det S=u\Delta
\]

be a Saito matrix.  Suppose each \(\delta_i\) lifts regularly to \(X\).
Then there are regular functions \(\tau_i\) on \(X\) such that

\[
 \widetilde\delta_i=\delta_i+\tau_i\partial_T,
\qquad
 \delta_i(E)+\tau_iE_T\in(E).
\]

If a target derivation \(\eta\) is transverse in the sense that
\(\eta(E)\) is a unit on the generic ramification divisor, its root lift is

\[
 \eta(T)=-\frac{\eta(E)}{E_T}.
\]

Hence its only generic ramification pole is the reciprocal marked-root
derivative.  If a rational source chart has volume factor
\[
 J_{\mathrm{src}}=vE_T^{-1}
\]
and the incidence projection has volume factor
\[
 J_{\mathrm{inc}}=wE_T,
\]
then their composite has constant Jacobian \(vw\) wherever it is defined.
It extends to a polynomial Keller map exactly when all reconstructed source
and target coordinates extend regularly to the chosen affine-space chart.

#### Proof

The first displayed relation is the definition of a regular lift tangent to
\(E=0\).  Implicit differentiation gives the formula for \(\eta(T)\).
Multiplication of the two Jacobian factors gives \(vw\).  Polynomial
extension is a separate regularity statement and is neither implied by
Saito's determinant identity nor by the local implicit-function calculation.
QED

This proposition is a dictionary, not a new construction theorem.  Its last
Jacobian sentence is the existing boundary-cancelled incidence mechanism.
The Saito matrix packages the tangent directions; the reconstruction pole
comes from a direction transverse to the free divisor.

## 7. Assessment of the proposed flow theorem

The proposed statement

> free discriminant + primitive logarithmic derivation + unimodular
> completion produces a Keller suspension exactly when an associated
> logarithmic flow algebraizes

is not correct without substantially stronger definitions and hypotheses.

1. A logarithmic derivation is tangent to the discriminant.  The reciprocal
   reconstruction pole is produced by a **transverse** derivation.
2. Primitivity and a Saito completion do not imply that a derivation is
   complete, locally nilpotent, semisimple, or equipped with a global slice.
   Those are the conditions that turn an infinitesimal field into a useful
   algebraic \(\mathbb G_a\)- or \(\mathbb G_m\)-flow.
3. Saito's determinant identity controls a divisor but not polynomiality
   after the reciprocal rechart.
4. Most decisively, the quadratic-gauge quartic is already a Keller
   suspension although its full branch and ledger divisors are not free.

A valid sufficient version is:

> If a logarithmic derivation lifts regularly, is complete with a global
> affine slice, its completion supplies a reciprocal volume chart, and every
> reconstruction function extends regularly on that chart, then the
> incidence suspension is Keller.

The converse becomes true only if “the flow algebraizes” is defined to
include the existence and regularity of precisely that reciprocal chart.  In
that form it is essentially a restatement of the polynomiality gate, not a
consequence of freeness.

## 8. First external candidate

The type-\(A_3\) reflection discriminant gives the cleanest control outside
the repository families.  Write

\[
 E_{A_3}(T)=T^4+pT^2+qT+r.
\]

The incidence is finite flat and its source is affine three-space, since the
equation solves for \(r\).  Its discriminant is

\[
 \Delta_{A_3}
 =16p^4r-4p^3q^2-128p^2r^2+144pq^2r-27q^4+256r^3.
\]

An exact Saito matrix is

\[
 S_{A_3}=
 \begin{pmatrix}
 2p&-16r&-6q\\
 3q&2pq&2p^2-8r\\
 4r&3q^2-8pr&pq
 \end{pmatrix},
\qquad
 \det S_{A_3}=2\Delta_{A_3}.
\]

Its root velocities are

\[
 T,\qquad -4T^3-4pT-3q,\qquad -2T^2-p.
\]

So this candidate passes four gates:

- finite flat root cover;
- free discriminant and explicit Saito matrix;
- affine-space incidence normalization;
- regular lifts of the logarithmic frame.

It does **not** yet pass the decisive fifth gate: no reciprocal
\((\partial_TE_{A_3})^{-1}\) chart has been shown to extend polynomially.
This \(A_3\) marked-root cover has the full group \(S_4\), so it is a
calibration rather than a proper-subgroup example.  The same test on a
non-\(A\) reflection orbit, beginning with \(B_3\), is the route to a proper
subgroup in its natural permutation representation.

## 9. Search order outside the current families

The next search should be gated rather than broad.

1. **Reflection discriminants.**  Continue from \(A_3\) to \(B_3\) and
   \(H_3\).  Their finite quotient maps and Saito data are built in.  Search
   only for reciprocal charts whose valuation vector cancels the reflection
   Jacobian and whose reconstruction ring is polynomial.
2. **Finite-flat free discriminants with fast normalization.**  These are
   closer to the marked-root problem than an arbitrary free divisor because
   the finite cover and normalization are already present.  Test the
   zero-dimensional families first and reject any example whose normalized
   reconstruction open has nontrivial unit rank or nontrivial class group.
3. **Quiver and prehomogeneous linear free divisors.**  Their Saito columns
   are infinitesimal algebraic group actions, so flow algebraization is easy.
   The missing datum is usually the finite root cover.  Require a finite
   primitive-element presentation before doing any Keller calculation.
4. **Three-component control.**  The star-quiver example in which
   \(h=\Delta_1\Delta_2\Delta_3\) is the product of three maximal minors is
   the smallest nontrivial multi-boundary test.  Compute its normalization
   and units before searching coefficients.  The normal-crossing divisor
   \(xyz=0\) should be retained as a rejection control: its diagonal Saito
   matrix is perfect, but it supplies no nontrivial connected root cover.

For every candidate, the stop/go ledger is:

\[
\begin{array}{c|c}
\text{gate}&\text{required certificate}\\ \hline
\text{finite cover}&\text{monic primitive equation or finite algebra}\\
\text{free divisor}&\det S=u\Delta\\
\text{root lift}&\delta_i(E)+\tau_iE_T\in(E),\ \tau_i\text{ regular}\\
\text{transverse pole}&\eta(T)=-\eta(E)/E_T\\
\text{affine chart}&\text{polynomial reconstruction ring}\\
\text{Keller ledger}&J_{\mathrm{inc}}J_{\mathrm{src}}\in k^*
\end{array}
\]

This ordering prevents an explicit Saito determinant from being mistaken for
a polynomial Keller construction.

## 10. External sources

- Buchweitz--Ebeling--von Bothmer,
  [*Low-dimensional Singularities with Free Divisors as Discriminants*](https://arxiv.org/abs/math/0612119),
  includes finite flat maps and fast normalizations.
- Buchweitz--Mond,
  [*Linear free divisors and quiver representations*](https://arxiv.org/abs/math/0509221),
  supplies quiver discriminants and explicit freeness tests.
- Granger--Mond--Schulze,
  [*Free divisors in prehomogeneous vector spaces*](https://arxiv.org/abs/0912.0626),
  gives structural restrictions and component results for linear free
  divisors.
- Granger--Mond--Schulze,
  [*Partial normalizations of Coxeter arrangements and discriminants*](https://arxiv.org/abs/1108.0718),
  is directly relevant to reconstruction opens.
- Antoniou--Feigin--Strachan,
  [*The Saito determinant for Coxeter discriminant strata*](https://arxiv.org/abs/2008.10133),
  provides explicit determinant structure across Coxeter types.

The external literature supports the candidate supply.  It does not supply
the reciprocal affine-space reconstruction chart, which remains the
construction-specific bottleneck.
