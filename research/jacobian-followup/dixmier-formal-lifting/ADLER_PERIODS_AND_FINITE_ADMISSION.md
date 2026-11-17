# Adler periods, microlocal mates and the finite differential admission gate

15 November 2026. **QAD1.** A characteristic-zero residue trace gives
an order-independent mate test for arbitrary scalar differential
operators on algebraic coefficient curves, including mixed coefficients.
For the already excluded G(a)+mu*v family, one fractional-power period
is a nonzero universal multiple of QPC1's arithmetic coefficient.

A separate leading-divisor gate excludes arbitrary mixed lower terms
whenever the leading coefficient multiplicities are coprime to the
operator order. In particular, every nonconstant squarefree polynomial
leading coefficient is excluded for finite differential mates.

The period conditions exactly characterize existence in the larger
formal pseudodifferential algebra after normalization. They do not
characterize finite differential mates. The latter have an exact
constant-linear admission criterion through the almost-commuting basis.
A harmonic-oscillator control passes every period condition and still
has no finite differential mate, so this distinction is essential.

The Adler trace and almost-commuting basis are established machinery,
not new priority claims. We give the proofs used here and the explicit
universal Long-family residue. The basis theorem is reviewed in
[Jimenez-Pastor et al., Theorem3.5](https://arxiv.org/html/2405.05421).
For the trace identity see
[Sleigh--Nijhoff--Caudrelier, Lemma1.1](https://eprints.whiterose.ac.uk/178183/7/rnab288.pdf).
The residue computation below is direct.

## The algebraic trace on a coefficient curve

Let F/K be an algebraic function field of one variable in characteristic
zero, delta a nonzero K-derivation, and h a nonzero differential constant.
Let omega be the unique rational one-form with df=delta(f)*omega.
Use finite differential operators D=F[a;h*delta] and the larger algebra
\[
 \Psi=F((a^{-1})),\qquad
 a^i f=\sum_{k\ge0}\binom{i}{k}h^k\delta^k(f)a^{i-k},
 \qquad \operatorname{res}_a P=[a^{-1}]P.                    \tag{1}
\]
The notation means series with an upper bound on a order and possibly
infinitely many negative powers. It does not allow unbounded positive
order, even if that order grows only with h-adic precision.

For monomials f*a^i and g*a^j put N=i+j+1. If N>=1,
\[
 \operatorname{res}_a[f a^i,g a^j]
 =h^N\binom{i}{N}\delta\left(
    \sum_{k=0}^{N-1}(-1)^k\delta^k(f)\delta^{N-1-k}(g)
                                      \right).              \tag{2}
\]
Indeed binom(j,N)=(-1)^N*binom(i,N), and the interior sum
telescopes on differentiation. If N<=0 the residue is zero.
Only finitely many monomial pairs in two bounded-above series can
contribute. Thus res_a[P,Q]*omega is an exact rational differential.

For any linear functional ell on Omega_(F/K)/dF, including a local
residue or a period, the map ell(res_a(P)*omega) is a trace.
In particular define the de Rham-valued fractional residues
\[
 {\cal A}_j(S)=
       [\operatorname{res}_a(R^j)\omega]\in\Omega_{F/K}/dF,
       \qquad R^n=S.                                      \tag{3}
\]
Here S is initially monic of order n>=2, and R is its unique
monic formal root of order one. Its coefficients are obtained
recursively: each next root coefficient is determined with a
nonzero scalar factor n, so no integration is needed.

## Mate necessity and transport

Suppose [S,T]=kappa with nonzero differential constant kappa,
where T may already be any element of Psi. Then
\[
 [R,T]=\frac{\kappa}{n}R^{1-n},\qquad
 [R^k,T]=\frac{k\kappa}{n}R^{k-n}\quad(k\in\mathbb Z).        \tag{4}
\]
To prove the first identity, subtract the proposed right hand side.
For the resulting C, the identity [R^n,T]=kappa gives
sum_(i=0)^(n-1) R^i*C*R^(n-1-i)=0.
If C is nonzero, its highest coefficient in this sum is multiplied
by n and cannot cancel. Thus C=0. The second identity follows
by the product rule and then inversion for negative k.

Taking the trace in (4) gives
\[
                 {\cal A}_j(S)=0\qquad(j\ne-n).             \tag{5}
\]
In particular all positive-index classes and A_-1 vanish.
This conclusion permits every finite differential mate order
and even infinitely many negative powers in a proposed mate.

Finite coefficient covers do not repair a nonzero class:
normalized field trace commutes with differentiation and retracts
Omega/dF. More directly, coefficientwise trace is a D-bimodule
retraction on both finite operators and bounded-above Psi series;
each output coefficient uses finitely many terms. A mate upstairs
therefore descends. This is the linear fixed-operator mechanism
of [QPC1](UNIVERSAL_P_POWER_OBSTRUCTION.md), not trace descent
of a nonlinear simultaneous tuple or an idempotent.

The classes are invariant under invertible pseudodifferential
conjugation: W*P*W^-1-P=[W,P*W^-1], so (2) applies. If a new
first-order generator is E=c*a+b, c nonzero, then
res_E(P)=c*res_a(P) and its dual one-form is omega/c.
Their product is exactly unchanged: differential powers contribute
no negative order, powers <=-2 cannot reach order -1, and
a^-1=(E-b)^-1*c has leading coefficient c at E^-1.
This includes algebraic position changes and momentum shifts.

A Fourier transformation does not automatically preserve this
chosen order completion. Transport across another polarization
needs a compatible microlocal/module comparison; it is not supplied
by the trace formula alone. No identification with QRF2 is asserted.

## The exact finite differential admission criterion

Put E=a+r0, where R=a+r0+negative powers, so that S in the
generator E is monic with zero coefficient of E^(n-1).
Let C_F=ker(delta), and define
\[
       B_j=(R^j)_+,\qquad
       C_j=[S,B_j],\qquad \operatorname{ord}C_j\le n-2.
                                                               \tag{6}
\]
The positive part is the finite differential part. Since
[S,R^j]=0 and (R^j)_- has order <=-1, its commutator with
the monic order-n S has order <=n-2.

Every finite T with ord[S,T]<=n-2 is a finite C_F-linear
combination of the B_j. If m=ord T, its highest commutator
coefficient is n*h*delta(t_m) at order n+m-1, so t_m is
constant. Subtract t_m*B_m and induct on the order.
Conversely every such combination has the required order bound.
The B_j are linearly independent because their orders are distinct.

Consequently the finite Weyl mate problem is exactly
\[
 [S,T]=\kappa\text{ for a finite }T
 \quad\Longleftrightarrow\quad
       \kappa\in\operatorname{span}_{C_F}\{C_j:j\ge0\}.
                                                               \tag{7}
\]
The span means finite combinations, not a completion or an F-linear
span. Unknown rational coefficient functions and partner poles have
been eliminated; the remaining unknowns are finitely many constants
at some, still unspecified, finite order. Checking a finite prefix
does not prove failure of (7). A useful all-degree application needs
a structural annihilator or a uniform statement about the whole span.

This is the classical almost-commuting-basis theorem, with its direct
proof given above. It is not a claim of a new general terminating
decision algorithm.

For a nonmonic S, first adjoin an n-th root c of its leading
coefficient. Its formal root has leading part c*a+r0.
Use E=c*a+r0 as generator; then S is monic and normalized in E,
with derivation c*delta. The finite coefficient change and its inverse
preserve differential orders and finiteness. Thus (6)--(7) apply to
arbitrary scalar differential operators after this finite normalization.
If a mate is constructed on the cover, coefficient trace descends it.

## What vanishing periods actually construct

In the normalized generator E, the following are equivalent:
S has a bounded-above pseudodifferential mate, and
\[
 \omega_E\text{ is exact},\qquad
 [\operatorname{res}_E(R^j)\omega_E]=0
                         \quad\text{for every }j\ge1.       \tag{8}
\]
Necessity is (5): res_E(R^-1)=1, and n>=2 makes -1!=-n.

For sufficiency, construct an order-zero dressing W=1+negative
powers with W^-1*R*W=E. Suppose the first surviving term of the
current conjugate is r_m*E^-m. The residue of its m-th power
is m*r_m, since all other possibilities have lower order.
Conjugation preserves its de Rham class, so (8) makes r_m a
derivative. Choose w with h*(c*delta)(w)=-r_m, and conjugate
by 1+w*E^-m. This removes that term without changing earlier ones.
The infinite product converges in the negative-order filtration.
The form omega_E is dual to c*delta, and exactness supplies
u in the coefficient field with (c*delta)(u)=1.

Then the explicit bounded-above mate is
\[
 T_0=\frac{\kappa}{nh}\,
             W u E^{1-n}W^{-1},\qquad [S,T_0]=\kappa.        \tag{9}
\]
Its order is 1-n. It can have infinitely many negative powers.
Obtaining a finite differential operator still requires (7);
the dressing construction does not make its negative tail finite.
The construction is over the coefficient field with h invertible;
the dressing itself can have negative h valuations and is not
being asserted to be a regular h-adic or bounded quantum gauge.
This separates the period obstruction from the remaining finite
lattice/admission problem.

## A divisor gate that retains finite differential admission

The finite differential criterion has a further leading-symbol
consequence on any coefficient curve, independent of all lower terms.
Let H be the order-n leading coefficient of S and define
\[
          g=\gcd\bigl(n,\{\operatorname{ord}_P(H)\}_P\bigr).
                                                               \tag{10}
\]
Use all zero and pole orders on the base coefficient curve.
If H is constant, this gcd is n. A finite differential mate
with nonzero constant commutator requires g>1.

Indeed a finite mate of positive order m with leading coefficient b
has zero highest commutator coefficient, hence
\[
 nH\delta(b)-m b\delta(H)=0,\qquad
       \delta(b^n/H^m)=0,\qquad (n/g)\mid m.                \tag{11}
\]
Differential constants have zero divisor, proving the divisibility
at every point. If g=1, then n divides m, and
b/H^(m/n) is a differential constant gamma. Subtract
gamma*S^(m/n) from the mate. This preserves its commutator and
strictly lowers its order. Repeat. An order-zero mate would
have leading condition delta(b)=0 and therefore zero commutator,
a contradiction. This finite order reduction proves the result
for every possible initial partner order.
In fact the same argument proves that, on the base coefficient curve,
every finite T with ord[S,T]<=n-2 belongs to C_F[S]. Thus the
entire low-order commutator image is zero in this divisor class,
not merely a checked finite prefix. The base centralizer statement
is not asserted for all coefficient covers.

Apply normalized trace before this argument for any proposed mate
over a finite algebraic coefficient cover. Ramification can change
the divisor gcd upstairs, but cannot repair the fixed base operator's
linear mate equation.

In particular any nonconstant squarefree polynomial H(v) has a
simple zero and g=1. Thus no order-n>=2 scalar operator with that
leading coefficient has a finite differential mate, with arbitrary
mixed lower terms, partner order or finite coefficient cover.
The polynomial degree is unrestricted. More generally the same
exclusion holds whenever the zero and pole multiplicities together
are coprime to n.

This is stronger for finite differential admission than the
leading period test below. Linear and centered quadratic squarefree
leading coefficients can have a PDO mate when the remaining periods
vanish; they still cannot have a finite differential mate.
Nonprimitive divisors g>1 are only a necessary escape condition.
The positive mixed family (22) has leading coefficient v^(3n),
so g=n and does not contradict the gate. Nor does this condition
exclude a constant leading coefficient.

On F=K(v), if the operator order n is prime, g>1 forces g=n.
Every zero and pole order is then divisible by n, so H is a scalar
times an n-th power of a rational function. This follows by factoring
the numerator and denominator; it is not asserted on arbitrary
positive-genus coefficient curves with possible divisor-class torsion.

## A leading-coefficient obstruction for arbitrary mixed lower terms

There is an explicit consequence of A_-1 that does not constrain any
lower coefficient. In the quadratic chart delta=v^-1*partial_v,
let a finite scalar operator S of order n>=2 have leading coefficient
H(v) in K(v), nonzero. Adjoin c with c^n=H(v).
A mate requires omega_E=v*dv/c to be exact, since
res_a(R^-1)=1/c.

After a finite extension of constants, the Kummer cover has its
cyclic root-of-unity action. Project any putative primitive to
the character carried by 1/c. It has the form Q(v)/c,
with Q rational. Hence the leading necessary condition is exactly
\[
  Q'-\frac{H'}{nH}Q=v,\qquad Q\in K(v).                    \tag{12}
\]
Finite further covers cannot help, by normalized trace of a primitive.
The same reduction applies to every rational H: use its actual cyclic
Kummer degree, including degree one. No polynomial squarefreeness
assumption is needed for the rational equation itself. Any auxiliary finite
constant extension can also be traced down in this linear equation.

For nonconstant squarefree polynomial H of degree D, any rational
solution Q has no finite poles. At a root alpha of H, a pole of
integer order exponent r<0 would have leading multiplier r-1/n
in (12), which is nonzero. At any other finite point, differentiation
of the pole gives an uncancellable highest pole.
Furthermore (12) forces Q(alpha)=0 at every root of H.
Thus Q=H*U with U polynomial, and
\[
 H U'+\frac{n-1}{n}H'U=v,\qquad
                         \deg U=2-D.                     \tag{13}
\]
The degree equality follows from the nonzero highest multiplier
deg(U)+(n-1)*D/n. It proves that D>=3 is impossible.

This excludes every finite scalar S with squarefree polynomial
leading coefficient of degree at least three, regardless of all
mixed lower coefficients, their poles, the partner order or any
finite algebraic coefficient cover. It excludes even a bounded-above
PDO mate. This is not limited to the separable v*H(a)+G(a) form.

The leading test has the following complete boundary classification.
For constant H it passes. For H=A*v+B of degree one it passes, with
\[
 U=\frac{n}{A(2n-1)}v
            -\frac{n^2B}{A^2(n-1)(2n-1)}.
                                                               \tag{14}
\]
For H=A*v^2+B*v+C squarefree of degree two, U must be constant.
Equation (13) is solvable exactly when B=0, in which case
U=n/(2*A*(n-1)). Passing this leading test is not sufficient
for a finite mate; lower spectral classes and (7) remain required.

More generally, on delta=v^(1-d)*partial_v, d>=1, the same argument
replaces the right side of (12)--(13) by v^(d-1) and forces
\[
                              \deg U=d-D.                 \tag{15}
\]
Thus squarefree polynomial leading coefficients with D>d
are excluded in every operator order n>=2 and with arbitrary
lower terms. At D=d, the leading test passes exactly for
H=A*v^d+B, subject to squarefreeness. In the ordinary coordinate
chart d=1, every squarefree leading polynomial of degree>=2
is excluded. These are written pole/divisibility proofs, not
bounded rational-ansatz searches.

## The characteristic-zero residue behind QPC1

Return to F=K(v), delta=v^-1*partial_v, omega=v*dv.
Take S=G(a)+mu*v, deg G=n>=2, leading coefficient g_n nonzero.
Normalize L=S/g_n, so the potential coefficient is nu=mu/g_n,
and let R=L^(1/n). There is an exact all-degree local period:
\[
 \operatorname{Res}_{v=0}
        \left(\operatorname{res}_a R^{2n+1}\,v\,dv\right)
   =-\frac{(2n+1)(n^2-1)}{24n^2}\,
                         \nu^2h^2\ne0.                    \tag{16}
\]
Thus the fixed S has no mate even in the bounded-above
pseudodifferential algebra over any finite coefficient cover.
This strengthens the finite-mate conclusion of QPC1 without
using good reduction.

The coefficient is independent of every lower coefficient of G.
Here is a selection argument and direct computation.
A term in L^beta, beta=2+1/n, containing ell potential factors,
c derivative contractions and total lower-G order loss d>=0
has a order at most n*beta-n*ell-c-d and v exponent ell-2c.
To contribute to both residues one needs
\[
       n\ell+c+d=2n+2,\qquad \ell-2c=-2.                   \tag{17}
\]
These order rules follow inductively from the root recursion,
or from the homogeneous symbol expansion. They force ell even.
For ell>=4 the first left side is too large; ell=0 would require
a derivative of a constant coefficient. Thus ell=2,c=2,d=0.
Only the leading a^n in G matters. This also proves that the
quadratic-in-potential second variation below computes the
whole residue, with no omitted higher-potential contribution.

Use Tr(P)=Res_(v=0)(res_a(P)*v*dv). It is a trace by (2).
Writing U=G(a)/g_n and introducing a scalar parameter epsilon,
cyclicity of root powers gives
\[
 [\epsilon^2]\operatorname{Tr}(U+\epsilon v)^\beta
   =\frac{\beta}{2}\operatorname{Tr}
          \left(v\,d(U^{\,\beta-1})[v]\right).              \tag{18}
\]
For example this follows by differentiating R^n=U+epsilon*v,
using cyclicity in Tr(R^k), and comparing the n product-rule terms.
No analytic functional calculus is required.

Put q=beta-1=1+1/n and f(z)=z^q as its formal expansion at
infinite a order. The second-contraction term in the linear
variation of f(U) is
\[
 \frac{h^2}{2}\delta^2(v)
       \left(\frac{f''(U)}2\,U''
                     +\frac{f'''(U)}3\,(U')^2\right).      \tag{19}
\]
Primes on U mean derivatives in the momentum symbol a.
One algebraic derivation expands the linear resolvent term
(z-U)^-1*v*(z-U)^-1 in the Ore product. Its second symbol
derivative is U''/(z-U)^3+2*(U')^2/(z-U)^4.
Formal divided differences replace these two terms by
f''(U)*U''/2 and f'''(U)*(U')^2/3 respectively.
This identity can equivalently be checked on integral powers
and extended by the unique monic root recursion.

By (17) set U=a^n in the relevant coefficient. The a^-1
coefficient in the parentheses of (19) is
\[
 \frac{q(q-1)}2 n(n-1)
       +\frac{q(q-1)(q-2)}3n^2
                         =\frac{n^2-1}{6n}.               \tag{20}
\]
Also v*delta^2(v)=-v^-2. Combining (18)--(20) gives
-beta*(n^2-1)/(24n)*h^2, proving (16) after restoring nu^2.

Writing gamma_n=(n^2-1)/(24n), the bridge to QPC1 is
\[
       \operatorname{Tr}(L^{\,2+1/n})
           =-(2+1/n)\gamma_n(\mu/g_n)^2h^2.                \tag{21}
\]
It identifies this explicit scalar-family residue with the
same rational coefficient as the arithmetic calculation.
It does not identify either class with the old restricted
order-seven cocycle QRF2.

## Positive controls and a decisive insufficiency control

There are finite mates with genuinely mixed coefficients.
For example E=v^3*a, u=-v^-1 gives [E,u]=h, so for every n>=2,
\[
                  S=E^n+u,\qquad T=-E,\qquad [S,T]=h.
                                                               \tag{22}
\]
Its leading normalization has dual form omega_E=v^-2*dv=du.
In the u,E chart, all fractional-root coefficients are polynomial
in u. Every resulting rational differential is exact.
The leading polynomial in v,a is a proper power and includes
negative-v coefficients; this positive family prevents extending
the separable obstruction to arbitrary mixed coefficients.

Conversely take S=a^2+u^2 on F=K(u), [a,u]=h.
Every coefficient of its monic square root, and of every integral
root power, is polynomial in u. Thus omega=du and every positive
Adler class vanish, proving (8) for the entire hierarchy.
It has the pseudodifferential mate (9).

It has no finite differential mate, even on finite coefficient covers.
[QSO2](QUADRATIC_OPERATOR_MATE_CLASSIFICATION.md) reduces any finite
mate to a first-order one after subtracting a polynomial in S.
For Q=b(u)*a+c(u), comparing the coefficients of [S,Q]=kappa gives
\[
       b'=0,\qquad c'=0,\qquad -2h\,u\,b=\kappa,
                                                               \tag{23}
\]
which is impossible for nonzero constant kappa. Trace descends
any algebraic-cover mate first. Thus even the complete vanishing
period hierarchy cannot replace finite differential admission.

## Controls, applicability and the remaining research target

For n=2, a full root calculation gives the explicit density
\[
 \operatorname{res}_a(a^2+\mu v)^{5/2}
 =\frac{5\mu(-3h^2+2\mu v^5)(h^2+\mu v^5)}{32v^7},
                                                               \tag{24}
\]
whose period is -5*mu^2*h^2/32. The separate graded sparse
controls for n=2,3,4,5 retain exactly the mu^2*h^2 terms and
check the respective residues -5/32,-7/27,-45/128,-11/25.
The projection is legitimate for these declared coefficients;
it is not a full arbitrary-order root computation.

[Exact controls](check_adler_admission.py) also check the
universal rational identity, monomial commutator primitives,
finite-cover trace, universal leading primitive, divisor and pole/degree identities,
the normalized almost-commuting basis and
constant-linear reconstruction in declared mixed examples,
positive mixed mates and the harmonic insufficiency control.
All-order coverage, the complete pseudodifferential criterion
and finite admission equivalence use the written proofs.

Operationally, a proposed mixed scalar S should first pass the
leading differential and spectral-period tests. If they vanish,
finite existence remains the constant-linear span problem (7).
A finite prefix or a pseudodifferential dressing is not a finite
Long carrier. A structural annihilator for the entire defect span,
or a theorem controlling that span for all bounded full carriers,
would be useful progress beyond another order search.

These scalar reductions do not normalize every deformed first
image and do not apply to arbitrary nonscalar matrix operators.
Full module constructions still require every right relation,
a single correct classical graph, global admission and finite-left
generation. No finite Long carrier, bounded full Long lift,
all-representative nonlifting, scalar or matrix counterexample,
independent theorem replay, formal verification, external review,
novelty priority or DC2 solution is claimed.
