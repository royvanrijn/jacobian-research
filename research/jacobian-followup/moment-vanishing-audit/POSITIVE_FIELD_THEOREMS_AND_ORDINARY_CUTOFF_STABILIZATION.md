# Positive field theorems and ordinary cutoff stabilization

15 November 2026. Internal written deductions, independent review pending.
The [nilpotent lifting theorem](NILPOTENT_LIFTING_AND_FORMAL_SAFE_COMPLETION.md)
supplies the topology behind the result. Section2 adds a direct coefficient
argument for varying GVC operators, with phase degree rather than the degree
of an auxiliary Gaussian product. The [small checker](verify_positive_cutoff_stabilization.py)
tests identities and algebraic controls, not the universal proofs. No novelty,
formal verification or implemented general producer is claimed.

## 1. Set-up and the field input

Work in one of these three models at bounded spatial dimension and degrees:

* Gaussian moments E(QP^m), with phase degree at most D;
* constant-coefficient GVC, Lambda^m(QP^m), with phase degree at most D
  and a finite polynomial symbol of bounded order r;
* full SIC, Ecal(g f^m), with total phase degree at most D in (xi,x).

Let R=Q[t_1,...,t_q] be a finite parameter polynomial ring. All phase and
symbol coefficients may be polynomial functions of t. Auxiliary parameters
and algebraic relations H in R are allowed; work in S=R/H. A localization
can be presented by adding an inverse parameter and its defining relation.
No flatness or reducedness assumption is imposed on S.

Let I be the full ordinary pure ideal in S, let I_N be the full ordinary
output-coefficient ideal for monomial markers of total degree delta and
m>N delta, and set

    B=S/I,    K_N=I_N/I.

Constants are markers of degree0, so I is contained in every I_N. The pure
equations and all output coefficients are kept literally, including embedded
and nilpotent structure. Let C>=1 be a common integer cutoff on the pure-safe
complex field points of this specified family, and let

    J=sqrt(H+I_C^R) in R,    j=J/(H+I^R) in B.                 (1)

The superscript R denotes the corresponding ideals before quotienting by H.
Thus V(j) is the reduced safe locus in Spec(B). A sufficient C is supplied
by the existing bounded field-uniformity theorem. If a particular family
already has a proved common field cutoff, that result can be used directly:
the deductions below then do not need the general holonomic-stratification
route. Empty families and empty safe loci are allowed with the usual zero-ring
interpretation. Take q>=1 by adding an unused parameter if necessary.

## 2. Direct coefficient derivatives for varying operators

The Gaussian derivative identity is already retained. For SIC,

    partial_c^beta Ecal(g f^m)
      =(m)_v Ecal(g M_beta f^(m-v)),  v=|beta|,                (2)

where M_beta is the product of the inserted phase basis monomials and has
degree at most Dv. These are polynomial identities including their outputs.

For GVC, first take independent symbol coefficients a_eta and phase
coefficients b_theta. For coefficient multi-indices alpha,beta put

    u=|alpha|, v=|beta|, k=max(u,v),
    D_alpha=product_eta (partial_x^eta)^alpha_eta,
    M_beta=product_theta (x^theta)^beta_theta.

Commutativity of constant-coefficient operators gives the exact identity

    partial_a^alpha partial_b^beta Lambda^m(QP^m)
      =(m)_u(m)_v D_alpha Lambda^(m-u)(Q M_beta P^(m-v)).       (3)

If u>m or v>m the derivative is zero. Otherwise align the two powers:

    RHS=(m)_u(m)_v D_alpha Lambda^(k-u)
         {Lambda^(m-k)[Q M_beta P^(k-v) P^(m-k)]}.              (4)

The multiplier inside braces has degree at most

    delta+Dv+D(k-v)=delta+Dk.                                 (5)

At a safe field point the braces therefore vanish when

    m>C delta+(CD+1)k.                                       (6)

The outer operators act on zero. Their order has no effect on (5). This
is the missing power alignment: directly replacing both powers in (3) by
the same exponent without inserting P^(k-v) would be false. Variable
coefficient operators would invalidate the commutations and are excluded.

Every coefficient derivative of total order l has k<=l. Thus, in all three
models, all derivatives of a marked output coefficient of orders less than K
vanish on V(J) if

    m>C delta+(CD+1)(K-1).                                   (7)

For polynomial parameterizations c=c(t), the repeated chain rule expresses
each t-derivative as a sum of coefficient derivatives of order at most l,
times parameter polynomials. Hence (7) also holds in the actual parameter
ring R, without assuming that the parameterization is linear or immersive.
Output coefficient extraction commutes with parameter differentiation.

Apply radical Zariski--Nagata and Ein--Lazarsfeld--Smith in the regular
complex polynomial parameter space, as in the retained nilpotent proof, and
contract ordinary powers. One obtains the direct all-model containment

    m>C delta+(CD+1)(qL-1)
      implies each marked output coefficient belongs to J^L. (8)

This argument is performed before quotienting by H+I^R. It is valid even
when that quotient is singular or nonreduced. For J=0 the relevant output
coefficients vanish identically; for J=R containment is vacuous. A sharper
nontrivial bound replaces q by the largest complex component codimension
of V(J). We keep q to avoid computing it.

Put

    B_L=C+(CD+1)(qL-1).                                     (9)

After quotienting by all pure equations, (8) gives K_N subset j^L for
N>=B_L. Conversely, sqrt(K_N)=j for N>=C, so Noetherianity gives j^s(N)
subset K_N for a finite s(N). This proves cofinality for arbitrary polynomial
coefficient families, using their own parameter count and phase degree.
In particular, varying-operator GVC no longer needs Gaussian degree r+D
for this thickness estimate.

## 3. Exactly when the ordinary cutoff ideals stabilize

**Theorem.** The following are equivalent for this family:

1. The ordinary ideals K_N, equivalently I_N in S, eventually stabilize.
2. The powers of j eventually stabilize.
3. V(j) is both open and closed in Spec(B).
4. There is one integer N bounding every pure-null safe point over every
   commutative Q-algebra compatible with this coefficient family.
5. The safe functor is represented by a closed subscheme of Spec(B).

Here compatible points mean coefficient maps from B; no finiteness or
nilpotence bound is imposed on the target algebra.

For (1)->(2), let K be the stable ideal. Choose s with j^s subset K.
Cofinality gives K subset j^(s+1), hence

    K=j^s=j^(s+1).

The powers then stabilize. For (2)->(1), if j^s=j^(s+1)=L, sufficiently
large K_N is contained in L by cofinality, while each such K_N contains
a power of j, hence contains L. Thus K_N=L eventually. L is a finitely
generated idempotent ideal; it is (e) for an idempotent e in B. This elementary
fact follows from Nakayama's determinant trick: if L=L^2, there is f=1+i,
i in L, with fL=0; then f^2=f and e=1-f generates L. Its zero locus is open
and closed. See [Stacks04PN, Lemmas10.21.3--10.21.5](https://stacks.math.columbia.edu/tag/04PN).

Conversely, if V(j) is open and closed, write it as V(e), e^2=e. Split
B=eB times (1-e)B. On the e=1 factor j is the unit ideal; on the e=0
factor j is contained in the nilradical and is nilpotent by finite generation
and Noetherianity. Hence j^s=(e) for large s. This proves (2)<->(3).

The safe functor is the union of Spec(B/K_N) by the retained lifting
theorem. If all safe points have one bound N, apply this assertion to the
universal safe point B/K_M for each M>=N. It forces K_N subset K_M,
and monotonicity gives equality. Conversely stabilization gives the common
bound. If the union is represented by Spec(B/L), each B/K_N is safe and
B/L is safe; these universal points force L=K_N for all sufficiently large
N. This proves (4)<->(5) as well. When these conditions hold, the closed
representative is precisely the open-and-closed component Spec(B/(e)),
including its original nilpotent structure, not merely its reduction.

This is a criterion for a closed representative inside the pure-null
coefficient scheme. It is not a claim about arbitrary noncontinuous maps
from a completed ring or a classification of every possible scheme embedding.

## 4. A positive field theorem becomes an ordinary ideal theorem

Suppose every pure-null complex field point in the family is safe, with
the supplied common cutoff C. Then j is the nilradical of B. It is nilpotent;
let nu>=1 satisfy j^nu=0, equivalently

    J^nu subset H+I^R.                                      (10)

The direct containment above gives

    I_N=I for every N>=C+(CD+1)(q nu-1).                      (11)

In particular, every pure-null point over every Q-algebra is safe, and there
is a common cutoff independent of the target's nilradical or length. A
particular target may have arbitrarily large nilpotence exponent: its
coefficient map still kills the fixed family ideal J^nu by the pure equations.

The finite certificate is ordinary ideal equality I_N=I, not just equality
of radicals. The full pure ideal and every I_N are theoretically computable
by the retained algorithms. Compute J, test positivity by J=sqrt(H+I^R),
find nu by exact ideal membership, then (11) gives a terminating search for
the least integer N with I_N=I. A sufficient finite pure prefix generating
I is also theoretically computable. No universal implementation, numerical
nu, least N or expanded coefficient cofactors are produced here.

For a connected pure-null coefficient scheme with at least one safe point,
the stabilization theorem has a converse: stabilization implies that every
pure field point is safe. A disconnected family is different. Restrict a
retained failure scaling line to e^2=e. Its cutoff ideals are all (e), and
the safe component e=0 is open and closed, although e=1 is still a field
counterexample. Connectedness cannot be omitted.

## 5. The universal bounded-family dichotomy

For full universal Gaussian and SIC families, I is homogeneous in phase
coefficients. For GVC, it is bihomogeneous of degree(m,m) in symbol and
phase coefficients for each pure output, hence homogeneous of total degree2m.
In all cases B is a nonnegative graded algebra with B_0=Q. It has no
nontrivial idempotents: an idempotent has degree-zero term0 or1; after
subtracting it or replacing e by1-e, its least positive homogeneous piece
cannot equal the corresponding piece of its square. Therefore Spec(B)
is connected. The origin is pure and safe.

Consequently, in each fixed universal bounded size, the following are
equivalent:

* the corresponding positive moment/vanishing assertion holds over fields;
* it holds over every Q-algebra;
* the full ordinary mixed-cutoff ideals eventually equal the pure ideal;
* the mixed-cutoff ideals stabilize at all;
* one common cutoff works for every pure-null safe point over every
  Q-algebra, with no imposed nilpotent-thickness bound;
* the full safe functor is a closed subscheme of the pure-null scheme.

In the positive case that subscheme is the entire pure-null scheme. In a
negative size, no such stabilization or uniform all-ring bound exists.
The retained counterexample scaling lines also prove this directly: a
nonzero marked tail yields arbitrarily thick nilpotent safe scales.
This dichotomy does not decide a previously open positive conjecture. It
identifies exactly which ordinary scheme conclusion follows from deciding it.

## 6. Consequences for the already positive families

**Ternary Gaussian cubics.** The retained
[ternary theorem](TERNARY_CUBIC_GMC_CLOSURE.md) supplies C=3 for every
complex phase of degree<=3 in three Gaussian coordinates. There are q=20
universal phase coefficients, D=3. If nu is the nilpotence exponent of the
radical of its full ordinary pure ideal modulo that ideal, then

    N>=200 nu-7 implies I_N=I_pure.                          (12)

Thus every all-order pure-null ternary cubic over every Q-algebra is safe,
with a common theoretically computable integer bound. This deduction uses
the ternary field theorem and direct symbolic containment; it does not need
the general Gaussian holonomic-stratification theorem. The internal analytic
review boundary of the ternary theorem remains inherited. Neither nu nor a
numerical ring bound is computed; the sharp field slope3 is not asserted
over nonreduced targets. The numerical pure-prefix gate94 and least-prefix
intervals concern field/radical testing and do not replace the ordinary ideal.
All protected prefix and Laplacian frontiers remain unchanged.

**Quadratic GVC phases.** For every fixed n and finite symbol order r,
the [quadratic-phase theorem](QUADRATIC_PHASES_FOR_ALL_OPERATOR_ORDERS.md)
supplies C=3(n+1), independently of r. Put D=2 and

    q=binomial(n+r,r)+binomial(n+2,2).

For a nilpotence exponent nu of the universal pure quotient, a sufficient
all-ring cutoff is

    3(n+1)+(6(n+1)+1)(q nu-1).                               (13)

The full ordinary mixed-cutoff ideal eventually equals the full ordinary
pure ideal, including the zero-phase and constant-symbol branches. This
uses the direct power-alignment lemma; it needs no higher-order scalar heat
reduction and no general Gaussian field stratification. The reduced two-jet
classifier remains a field classifier: the ordinary pure ideal need not equal
its reduced equations, and (13) is not an order-independent ring bound.
The new statement covers each finite r. It does not cover the universal
infinite formal-symbol coefficient ring by a finite-parameter argument.

## 7. Source and assurance boundary

The direct identities, stabilization criterion, connected-cone argument
and ideal consequences are written proofs. The new computation checks
small varying-symbol/phase derivatives, power alignment, a complete linear
Gaussian ideal example, nilpotent and idempotent controls. It does not
compute the universal ternary or quadratic-GVC ordinary ideal.

[Long2607.18186](https://arxiv.org/abs/2607.18186) was rechecked: the latest
listed arXiv version is still v1,20July2026. Its existing quartic/cubic
counterexamples are retained without a new replay or changed assurance;
no Jacobian announcement is promoted to an independently verified theorem.
This targeted source check is not an exhaustive fresh literature review.
Independent specialist review, general producer implementation, exact
individual negative-size cutoff ideals and structural normal forms remain
open. The goal remains active; quantitative cubic refinement is secondary.
