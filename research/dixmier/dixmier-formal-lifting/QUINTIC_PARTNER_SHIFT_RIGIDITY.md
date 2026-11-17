# Rational shifts of the quintic canonical partner preserve spectral rigidity

**QPS1.** The exact scalar spectrum in QQL1 remains \(5ch\mathbb Z\)
after every rational shift of the lower canonical partner:
\(K\mapsto K+g(\kappa)\).
This includes parameter-small shifts with arbitrary negative powers
of \(\kappa\), poles at other positions, and arbitrary nonzero rational
source Euler functions.

Consequently these shifts do not give an actual compatible Euler pair
or a rational-target fifth root. The result concerns scalar coefficient
shifts of this one partner; it does not treat genuinely mixed
\(\kappa,K\) corrections, a deformed first image, or nonscalar carriers.

## The shifted operator and statement

Let \(k\) have characteristic zero and \(h,c,c_v\in k^\times\). Use
\[
 \mathcal D=\operatorname{Frac}\bigl(k(v)[b;h\partial_v]\bigr),\quad
 [b,v]=h,\quad H=vb,\quad
 \kappa=c_vv^5F_s(H)^{-1},\quad F_s\in k(H)^\times .
\]
The partner \(K=-\kappa^{-1}H/5\) satisfies \([\kappa,K]=h\).
For any \(g\in k(\kappa)\), keep
\[
 U=\frac{(\kappa-3)^3}{(\kappa-2)^2},\qquad
 R=-\frac{(\kappa-3)(\kappa-2)}{\kappa^2}
\]
and put
\[
 J_g=cR\bigl(H-5\kappa g(\kappa)\bigr).                     \tag{1}
\]
This is exactly the transported Euler element obtained by replacing
\(K\) with \(K+g(\kappa)\) in QCS1's scalar cotangent/Mobius cubic model.
Its derivation on the coefficient field remains
\[
              [J_g,f(\kappa)]=\delta f,\qquad
              \delta=5ch\kappa R\partial_\kappa .
\]
The claim is
\[
 \boxed{\quad
 \{\lambda\in k:\exists\,0\ne X\in\mathcal D,\
                         [J_g,X]=\lambda X\}=5ch\mathbb Z
 \quad}                                                     \tag{2}
\]
for every rational \(g,F_s\).

One may enlarge the scalar field algebraically during the negative
proof; an eigenoperator over the original field remains an actual
operator after that extension. The realizers \(U^n\), for every
integer \(n\), belong to the original field. Thus scalar-spectrum
equalities obtained after extension descend.

## Shifts with at most a simple pole

Use the \(v\)-adic skew-series valuation from QQL1. If \(g\) has
at most a simple pole at \(\kappa=0\), let \(a_{-1}\) be its
\(\kappa^{-1}\) coefficient, zero when \(g\) is regular, and put
\(\theta=5a_{-1}\). The least \(v\)-weight of \(J_g\) is \(-10\),
with coefficient
\[
 -6c\,c_v^{-2}
       F_s(H-5h)F_s(H-10h)(H-\theta).                       \tag{3}
\]
For an eigenoperator with leading term \(v^mA(H)\), restrict divisors
to the orbit \(\theta+h\mathbb Z\). The last factor in (3) contributes
one singleton. The orbit divisor is therefore
\(1+(T^5+T^{10})f(T)\).
The primitive tenth-root argument of QQL1 gives \(10\mid m\).

Fifth-root deck projection consequently makes every scalar
eigenoperator deck-fixed. The fixed division ring is
\[
 E=\langle\kappa,H\rangle_{\rm div}
   =\operatorname{Frac}\bigl(k(\kappa)[J_g;\delta]\bigr);
\]
the equality still holds because the added term in \(J_g\) is rational
in \(\kappa\). A rational leading coefficient \(f(\kappa)\) of an
eigenoperator then satisfies \(\delta f=\lambda f\).
At infinity \(\delta\sim-5ch\kappa\partial_\kappa\), so
\(\lambda\in5ch\mathbb Z\). The opposite inclusion follows from
\([J_g,U^n]=-5ch nU^n\). This proves (2) in this case.

## Higher poles and the source divisor character

Suppose \(g\) has pole order \(\ell\ge2\), with leading term
\(a_{-\ell}\kappa^{-\ell}\), and set \(q=\ell+1\).
The least weight of \(J_g\) is now \(-5q\). Its coefficient is
\[
 30c\,a_{-\ell}c_v^{-q}
                   \prod_{j=1}^{q}F_s(H-5jh).              \tag{4}
\]
For a scalar eigenoperator with leading term \(v^mA(H)\), comparison
at weight \(m-5q\) gives the same shifted-divisor equation as in QQL1,
with translation length \(5q\).

For any orbit \(\alpha+h\mathbb Z\) of zeros and poles of \(F_s\), let
\(f_\alpha(T)\) encode its finite signed divisor. The divisor of (4)
on that orbit is
\[
            d_\alpha(T)=(T^5+T^{10}+\cdots+T^{5q})f_\alpha(T).
\]
If \(f_\alpha(\zeta_5)\ne0\) for some orbit and a primitive fifth
root \(\zeta_5\), then
\(d_\alpha(\zeta_5)=q f_\alpha(\zeta_5)\ne0\).
The lowest-coefficient equation forces \(\zeta_5^{-m}=1\), hence
\(5\mid m\). This is enough for deck projection: each nonzero
deck component of an eigenoperator has weight zero modulo five.
The same rational-leading-coefficient argument in \(E\) proves (2).

For example a source function of degree at infinity not divisible
by five must have a nonzero fifth-root divisor character on some
orbit. In particular quadratic polynomial source functions fall
in this case. Degree alone is not needed for the theorem.

## The zero-character case is an actual rational norm coordinate

It remains to treat source functions for which
\(f_\alpha(\zeta_5)=0\) on every orbit.
Since these are Laurent polynomials with integer coefficients,
this is equivalent to divisibility by
\(\Phi_5(T)=1+T+T^2+T^3+T^4\).

There is an explicit rational norm criterion, including poles:
after an algebraic scalar extension,
\[
                  F_s(H)=C\prod_{i=0}^{4}B(H+ih),
                  \qquad C\ne0,\quad B\in\bar k(H)^\times. \tag{5}
\]
To prove it, on each orbit divide \(f_\alpha\) by \(\Phi_5\)
in \(\mathbb Z[T,T^{-1}]\). The resulting finite signed divisor
\[
                         b_\alpha(T)=T^4f_\alpha(T)/\Phi_5(T)
\]
is the divisor of a rational function \(B\) on that orbit.
There are only finitely many nonzero orbits. Their rational factors
combine to one \(B\), and its five shifted product has the same
finite divisor as \(F_s\); the ratio is a nonzero scalar \(C\).
Only finitely many algebraic scalar coefficients are needed.
No positivity of any divisor multiplicity was used.

Set, with the indicated ordering,
\[
                   w=vB(H)^{-1},\qquad p=w^{-1}H.
\]
Then
\[
 w^5=Cv^5F_s(H)^{-1},\quad [H,w]=hw,\quad[p,w]=h,\quad
                   v=wB(H),\quad b=v^{-1}H .
\]
Thus \(w,p\) are actual rational Weyl coordinates generating the
entire source division ring, and \(\kappa=(c_v/C)w^5\).

In these coordinates the shift extends to an actual ambient
division-ring automorphism:
\[
           w\longmapsto w,\qquad
           p\longmapsto p-5\kappa g(\kappa)w^{-1}.
\]
The correction is rational in \(w\), so the Weyl relation is retained
and the inverse translation exists. This fixes \(\kappa\) and sends
\(H=wp\) to \(H-5\kappa g(\kappa)\). Consequently it sends
\(J_0=cR(\kappa)H\) to \(J_g\).
QQL1's exact spectrum for \(J_0\), now with constant source Euler
function, proves (2). This completes all rational source functions
and all rational partner shifts.

## Consequences and exact scope

There are no actual \(P,Q\in\mathcal D\) with
\[
                 [P,Q]=ch,\qquad J_g=QP+\lambda,\quad\lambda\in k,
\]
since they would give the forbidden eigenvalue \(-ch\).

Likewise no fifth root of \(U\,G(J_g)\), for any
\(G\in k(J_g)^\times\), belongs to the actual source division ring.
Indeed commutation would force the root to have eigenvalue \(-ch\);
the ordinary-order linearized-power injectivity proved in QQL1
applies unchanged. Thus rational coefficient shifts of the lower
partner do not reopen the homogeneous Euler compatibility route.

For deformation applications a parameter-small \(g\) preserves the
classical Long symbols. The theorem covers these shifts regardless
of their pole order, but supplies no global polynomial normalization
and no theorem about arbitrary mixed \(\kappa,K\) corrections.
Its actual division-ring spectrum does not constrain completed
series with unbounded spatial tails.

[Exact controls](check_quintic_partner_shifts.py), with
[the retained receipt](quintic_partner_shift_controls.json), verify
the shifted lowest coefficients, signed norm characters, actual
norm-coordinate identities, rational ambient momentum translation,
and canonical eigenrelations. Universal exhaustion of all shifts
uses the written pole split, divisor norm criterion, deck projection,
and QQL1 spectrum proof.

No bounded Long lift, DC1/DC2 solution, formal verification,
independent whole-proof replay or external review is claimed.
