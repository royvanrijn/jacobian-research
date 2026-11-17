# Rational mate obstructions, spectral torsion and projective frames

15 November 2026. **QST1.** The universal spectral class from QSC1
can distinguish finite operators, rational differential operators
and formal pseudodifferential operators. For the retained family
G(a)+mu*v with nonlinear G and nonzero mu, it is nontorsion:
no scalar spectral localization kills it. In fact the Adler
obstruction excludes a mate in the entire rational operator skew
field, with arbitrary operator denominators.

This conclusion passes to every finite algebraic coefficient cover
and every rational conjugate. In an exact-R chart it therefore
excludes rank-one projective carrier realizations in the same
rational pair orbit, without assuming a scalar free frame.

There is a second, independent rational residue: the constant
term of a rational torus action. It excludes rational mates for
every rational function of Euler operators, including harmonic
oscillator Hamiltonians and all their polynomial reparametrizations,
in every Weyl rank and after identity stabilization. The harmonic
example passes all Adler periods and has a PDO mate, so this
rational obstruction is strictly additional.

The ingredients are established operator and invariant theory
methods. The proofs and spectral-class calculations used here are
given directly; no novelty priority or identification with QRF2
is claimed. For related context see
[Bavula, Section5](https://arxiv.org/html/math/0403037) on semisimple
Weyl elements and
[Compoint--van der Put--Weil](https://arxiv.org/abs/1001.0153)
on differential-operator descent. Our fixed linear equation
descent below is not a claim of arbitrary differential-module descent.

## Keep the three operator algebras distinct

Let F/K be a characteristic-zero algebraic coefficient curve with
nonzero derivation delta, h a nonzero constant, and set
\[
 D=F[a;h\delta],\qquad Q=\operatorname{Frac}D,\qquad
                           \Psi=F((a^{-1})).                \tag{1}
\]
The ring D consists of finite differential operators. Every
nonzero differential operator has an inverse in Psi: remove its
leading monomial and invert its order-negative remainder by a
geometric series. Hence Q embeds faithfully in Psi.

A rational mate in Q therefore is a bounded-above PDO mate.
The converse is false. Neither algebra contains every h-adic
series with unbounded positive momentum orders. Scalar fields
such as K0((h)) are permitted throughout.

For a positive-order S in D, the universal spectral module from
[QSC1](LONG_SPECTRAL_CURVATURE_AND_SCALARIZATION.md) has
\[
 {\cal H}_S=D/[S,D],\qquad
 \alpha_S=[h],\qquad
        u\cdot[T]=[T S].                                  \tag{2}
\]
This is a K[u]-module: right multiplication by S commutes with
ad_S. Its cyclic annihilator is
\[
 I_S=\{q(u)\in K[u]:[S,T]=h\,q(S)
                                  \text{ for some }T\in D\}. \tag{3}
\]
Thus alpha_S is torsion exactly when such a nonzero q exists.
A torsion class gives the rational mate T*q(S)^-1 in Q.
The converse follows from the later
[QPD1 spectral-denominator theorem](SPECTRAL_DENOMINATORS_AND_PROJECTIVE_CENTRALIZERS.md):
for this fixed linear mate equation a rational denominator can be
cleared by a nonzero polynomial in S. This strengthens the admission
boundary of QST1; its Adler and torus proofs below are unchanged.

Localization in scalar spectral polynomials is distinct from
allowing arbitrary moving coefficient denominators in (v,u).
The nontorsion statements concern (2)--(3) and their flat scalar
spectral extensions. They do not silently identify these different
localizations or all fiberwise self-Ext groups.

## A first nonzero Adler moment separates every spectral polynomial

Let S be monic of order n>=2, R=S^(1/n), and let
Tr(P)=ell(res_a(P)*omega) be an Adler trace as in
[QAD1](ADLER_PERIODS_AND_FINITE_ADMISSION.md).
Assume some moment Tr(R^j) is nonzero. There is a smallest
such integer J: moments below -1 vanish by order.
Moreover J is not divisible by n, since nonnegative integral
powers of S are differential and negative integral powers
have order at most -n.

For each integer d define the explicit functional
\[
       \Lambda_d(T)=\operatorname{Tr}(R^{J-nd}T),\qquad
                         \Lambda_d([S,T])=0.                \tag{4}
\]
The second identity is cyclicity, since R commutes with S.
For a polynomial q of degree d and leading coefficient q_d,
\[
                         \Lambda_d(q(S))
                                =q_d\operatorname{Tr}(R^J)\ne0.
                                                               \tag{5}
\]
All lower powers have moments with index below J.
The same argument works for a nonzero rational q(S), using
its highest power d at infinity and the bounded-above formal
Laurent expansion. Only finitely many terms can affect a residue.

Consequently K(S) meets [S,Psi] only in zero. In particular
I_S=0, and S has no mate even in Psi or Q.
Equivalently one may multiply a proposed relation
[S,T]=q(S) by the invertible q(S)^-1 in Psi and obtain a
forbidden PDO mate. Formula (4) supplies explicit dual
functionals, rather than just that implication.

The global class remains nonzero after localization to K(u)
and after any field extension of that scalar spectral field.
This is flat base change of the complex defining (2).
It is not a statement that every closed spectral fiber class
is nonzero; QSC1's p^2 control shows why fiberwise reasoning
needs separate justification.

## Evaluation for the retained Long-family operator

Take F=K(v), delta=v^-1*partial_v, omega=v*dv, and
L=(G(a)+mu*v)/g_n, where G has degree n>=2 and leading
coefficient g_n nonzero, and mu is nonzero. QAD1 gives
\[
 \begin{split}
 \operatorname{Tr}(L^{j/n})&=0\quad(j<2n+1),\\
 \operatorname{Tr}(L^{2+1/n})
  &=-\frac{(2n+1)(n^2-1)}{24n^2}
                           (\mu/g_n)^2h^2\ne0,
 \end{split}                                                \tag{6}
\]
where Tr is the local double residue Res_(v=0)*res_a with
one-form v*dv.

The lower-moment vanishing follows from the same root homogeneity,
not a finite list. A term with ell potential factors, c coefficient
derivatives and lower-G order loss d>=0 can contribute only if
n*ell+c+d=j+1 and ell-2c=-2. The first possibility is
ell=2,c=2,d=0, at j=2n+1. The ell=0 possibility differentiates
constants and vanishes. This proves both the first moment
and all the zeros required for (4)--(5).

Therefore the actual global spectral class alpha_L has zero
annihilator, for every finite polynomial G and every nonzero mu.
It cannot be repaired by a polynomial spectral denominator.
More strongly there is no rational mate with any operator
denominator at all, because Q is contained in Psi.

This conclusion allows arbitrary scalar coefficient fields,
including Laurent-series scalars and h-small highest terms.
It concerns this fixed operator, not all quantum corrections
to Long's classical S. A different mixed operator can have
different periods and a different spectral class.

## The sharp torsion comparison: constant-coefficient polynomials

Let D=K(x)[a;h*partial_x] and S=G(a), with G nonconstant.
Then the exact answer to (3) is
\[
 I_G=(\Delta_G(u)),\qquad
 \Delta_G(u)=\prod_{\text{distinct critical values }c}(u-c).
                                                               \tag{7}
\]
Over K this is the monic squarefree part of
Res_a(G'(a),u-G(a)). For an affine G, Delta_G=1.

Here is an all-partner proof. The finite centralizer of G(a)
is K[a]: its highest commutator coefficient forces the leading
coefficient of any commuting operator to be constant, and
subtracting that constant multiple of a power of a lowers
the order. If [G,T]=h*q(G), differentiate in x. Then T_x
commutes with G, so every rational coefficient of T has a
constant derivative. Thus
\[
 T=x H(a)+J(a),\qquad
                      [G,T]=h G'(a)H(a).                  \tag{8}
\]
Hence q belongs to I_G exactly when G' divides q(G).

At a critical point alpha of multiplicity m in G',
G-G(alpha) has order m+1. One zero of q at the critical value
therefore supplies enough divisibility for the entire multiplicity
m, and a nonzero value supplies none. This proves (7), including
coincident and degenerate critical points.

Every nonlinear G has a nonzero finite obstruction class,
but it is torsion. A rational mate always exists:
\[
                         T=x\,G'(a)^{-1},\qquad [G,T]=h.
                                                               \tag{9}
\]
Thus adding mu*v changes the rational problem qualitatively:
at mu=0 there is rational completion, while every nonzero mu
in (6) has a nontorsion class and no rational completion.

## Fixed rational mate equations descend through finite coefficient covers

Let E/F be finite algebraic in characteristic zero and extend
delta uniquely. Write Q_E=Frac(E[a;h*delta]).
After passing to a finite Galois closure E'/F, coefficient
automorphisms commute with delta and act on Q_(E').

Their fixed division ring is exactly Q_F. To prove it, for an
invariant U use the nonzero right denominator ideal
\[
            {\cal I}_U=\{B\in E'[a]:U B\in E'[a]\}.         \tag{10}
\]
The Ore polynomial ring is right Euclidean. This ideal has a
unique monic generator B. Invariance of U makes the ideal
Galois stable, so uniqueness fixes B and also A=U B.
Both have coefficients in F, and U=A*B^-1 belongs to Q_F.

Consequently the Reynolds average
\[
        {\cal R}(U)=\frac1{|\Gamma|}
                           \sum_{\sigma\in\Gamma}\sigma(U) \tag{11}
\]
is a Q_F-bimodule retraction Q_(E')->Q_F. Any fixed linear
equation [S,U]=b with S,b in Q_F descends by this average.
This applies to rational mates with arbitrary operator denominators,
not only to coefficientwise trace of a finite differential operator.

No descent of a nonlinear simultaneous tuple, idempotent or
projective frame is asserted. The result uses the one-derivation
Euclidean operator ring; it is not silently extended to a
multivariable denominator ideal.

A nontrivial quadratic control illustrates rational descent.
For w^2=v and delta=v^-1*partial_v, put b=h/(2v^2) and
Q0=a^2+b*a-v-h^2/v^4. Direct Ore multiplication gives
\[
 \begin{split}
 Q_0&=(a+w)(a+b-w)=(a-w)(a+b+w),\\
 \tfrac12\big((a+w)^{-1}+(a-w)^{-1}\big)
                                     &=(a+b)Q_0^{-1}.
 \end{split}                                                \tag{12}
\]
The noncommutative rational average has an actual base-field
fraction, without treating ordinary coefficient trace as a
multiplicative operation.

## A torus residue detects rational obstructions missed by Adler periods

Let Q be a division algebra with a rational algebraic torus
action sigma. For one torus variable z, sigma_z(U) is a rational
function in the central z with coefficients in Q and has a
Laurent expansion at zero. Define
\[
          {\cal P}(U)=\operatorname{CT}_{z=0}\sigma_z(U).
                                                               \tag{13}
\]
For several torus variables take iterated constant terms.
The coaction identity shows P(U) is invariant; on the invariant
division subalgebra B=Q^sigma, P is the identity and a B-bimodule
map. This is a linear retraction, not an algebra homomorphism.

If S belongs to the center of B, it follows that
\[
       {\cal P}([S,U])=[S,{\cal P}(U)]=0,\qquad
                            [S,Q]\cap B=0.                 \tag{14}
\]
Every nonzero invariant right hand side, including h and any
nonzero rational function of S, is excluded.

In the Weyl division algebra use independent scalings
x_i->z_i*x_i, p_i->z_i^-1*p_i and H_i=x_i*p_i.
Invariance implies commutation with every H_i by differentiating
the torus action. Thus every rational function of the commuting
H_i is central in B and satisfies (14). This covers every
polynomial degree and every Weyl rank.

The action extends trivially on spectator Weyl pairs, so the
obstruction persists under identity stabilization. Actual
division-algebra automorphisms transport S, sigma and P by
composition, including Fourier changes and rational inner gauges.
No extension to the same momentum completion is required.
Classical Poisson coordinate changes without an actual quantum
automorphism are not being applied.

For a sharp control, on [a,x]=h set H=a^2+x^2. Over K(i), use
\[
        X=a+i x,\quad P=(a-i x)/(2i),\quad
             [P,X]=h,\qquad H=2i\,X P+i h.                 \tag{15}
\]
Hence H and every nonconstant rational function of H have
no rational mate, even after identity stabilization.
In the one-variable coefficient setting, (11) also excludes
rational mates after every finite algebraic coefficient cover.

Yet QAD1 proves every Adler class of H is exact and constructs
a bounded-above PDO mate. Therefore rational operators form a
strictly smaller admission category than PDO series. A zero
Adler hierarchy cannot certify a rational mate.

The obstruction applies when S itself requires a conjugate.
A valid Weyl tuple can of course contain Euler Hamiltonians as
composite observables; the identity map is an immediate example.
No exclusion follows merely from such a composite energy.

## Finite differential, rational and projective admission

The primitive-divisor gate of QAD1 is a finite differential
obstruction and must not be promoted to rational nonexistence.
For example, on [a,v]=h/v put
\[
 S=v a^2+\frac{h}{2v}a,\qquad A=\frac{v^2a}{3},\qquad
                     [S,A]=hS.                            \tag{16}
\]
The squarefree leading coefficient v excludes every finite
differential mate. But A*S^-1 is a rational mate.
On the cover w^2=v this is the same construction from
E=w*a, S=E^2 and the algebraic coordinate u=2w^3/3.
It also illustrates the descent in (11).

More generally a rational mate U of a finite S has an exact
finite Darboux certificate. Take the monic generator B of
its right denominator ideal, and A=U B. That ideal is stable
under left multiplication by S. Hence for some finite C,
\[
                       SB=BC,\qquad SA-AC=hB.              \tag{17}
\]
Conversely (17), with B nonzero, gives [S,A*B^-1]=h.
A denominator polynomial in S is a special case C=S.
The later QPD1 clearing theorem shows that a spectral polynomial
denominator always suffices for rational mate existence in this
coefficient-curve setting. Formula (17) remains a general finite
Darboux certificate. This does not identify arbitrary moving
coefficient localization of the universal spectral module with
scalar spectral localization.

Every rational conjugate of the obstructed family (6) has
the same no-mate conclusion. Its Adler functionals transport
by composition with the inverse conjugacy, using cyclicity
in Psi. Restricted to Q, a dual also transports through any actual
division-algebra isomorphism by composition with its inverse;
the original momentum completion need not be preserved.
In particular multiplication gauges b(v) give the
whole mixed family
\[
       b^{-1}(G(a)+\mu v)b
             =G(a+h\,\delta(b)/b)+\mu v,\qquad b\in K(v)^*.
                                                               \tag{18}
\]
They can introduce negative-v coefficients and leave the pure
separable form. No rational partner appears. More general
rational Darboux conjugates are excluded as well, even when
their displayed operator coefficients change radically.

There is a direct full-carrier consequence in the rational
Long chart. Let Q2 be the rank-two Weyl division ring and
Q1=Frac(K(r,v)[a;h*v^-1*partial_v]). Then
\[
                              \operatorname{Cent}_{Q_2}(r)=Q_1.
                                                               \tag{19}
\]
Indeed expand a rational operator in bounded-above powers
of d over Q1, where [d,r]=h. Commutation with r is
-h times formal differentiation in d, so every nonzero
d power vanishes. This argument allows rational denominators.

By QMC1 a generic rank-one projective carrier endomorphism
algebra embeds in Q2. If, in one fractional realization or
after an actual division-algebra change, its first two right
images become r and a rational conjugate inside Q1 of (6),
then its required mate would commute with r and belong to Q1.
That is impossible. No free rank-one frame is needed.

Changing the fractional realization conjugates the entire
right action rationally. Therefore this is an invariant
obstruction for that actual rational pair orbit, including
its nonpositional rational conjugates. It does not prove
that every possible Long deformation lies in that orbit.
Arbitrary formal gauges may be nonrational.

Similarly a carrier coordinate rationally equivalent to
the torus-obstructed functional class in (14) cannot have
its required rational conjugate. General nonscalar matrix
actions of higher left rank remain outside these scalar
division-ring arguments.

## Controls and remaining target

[Exact controls](check_rational_spectral_obstructions.py) verify
the first full quadratic Adler moment and lower zeros, universal
triangular spectral separation on declared degrees, critical-value
resultants and their Ore witnesses, the quadratic rational
Reynolds denominator identity, mixed rational escape (16),
finite Darboux certificates, harmonic canonical normalization
and nonhomogeneous rational torus Laurent expansions.

The universal statements use the written proofs. The torus
and Adler residues have different scopes; neither is being
identified with the old QRF2 cocycle. The all-cover rational
descent is a fixed linear equation result, not a general
algebraization theorem.

The next meaningful carrier test is to locate its actual
rational pair orbit and evaluate its global classes there.
Passing finite differential tests, passing every Adler period,
and admitting a formal PDO mate are distinct milestones.
A surviving full carrier still needs every relation, one
correct classical graph, global admission and finite-left
generation. No such Long carrier, bounded full lift, universal
Long nonlifting, scalar or matrix counterexample, independent
theorem replay, formal verification, external review, novelty
priority or DC2 solution is claimed.
