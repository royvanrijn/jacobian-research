# Polynomial submersions, bijective Gauss-Manin action and the volume mark

2026-11-16. Claims HC4MYLP0/1.

FULL THEOREM TARGET: prove HC4/JC2 or construct an actual counterexample.
CURRENT LAST UNSOLVED IMPLICATION: ACTUAL polynomial plane source forces
Tr(Lambda_L*Lambda_L) subset B, or full marked-energy integrality; OPEN.
External review/outreach is DEFERRED. The previous goal turn made PROGRESS
by classifying all source additive actions extending to the full finite
normalization (LO0/1). None was shown to exist automatically.

This note returns to LG1's raw derivative quotient. Its BIJECTIVE
transverse connection is already automatic for EVERY polynomial
SUBMERSION of the whole plane, before any Keller mate is supplied.
We write the divergence correction and an elementary proof, and give
unbounded-degree whole-plane geometric controls of arbitrarily large
generic genus. These replace an abstract-only control by an actual
polynomial submersion. They are NOT Keller pairs or HC4 counterexamples.

## 1. LP0: the canonical corrected connection for any submersion

Let R=C[x,y], Omega=dx wedge dy and let P in R have NO critical point.
Its partial derivatives generate the unit ideal. Choose a polynomial
vector field E with E(P)=1, and put

    D=P_y*partial_x-P_x*partial_y,
    b=div(E), eta=i_E(Omega), M_P=R/D(R).

The letter b here is a SOURCE polynomial, not a target coefficient.
Then

    i_D(Omega)=dP, dP wedge eta=Omega, eta(D)=-1,
    [E,D]=-b*D, d(eta)=b*Omega.                        (1)

Define an operator on the ACTUAL RAW quotient, without localization:

    nabla[a]=[E(a)+b*a].                              (2)

It is well-defined because (E+b)D=D E. It satisfies
[nabla,P]=1. Another lift of partial_P has the form E'=E+hD for a
polynomial h, since (E,D) is an invertible polynomial frame. Its
divergence is b'=b+D(h), and

    (E'+b')(a)-(E+b)(a)=D(h*a).

Thus (2) is INDEPENDENT of the chosen lift. Identifying [a] with the
relative one-form class [a*eta] identifies it with the usual relative
Gauss-Manin connection. A change of lift changes eta by h*dP, so the
identification is also independent in relative forms. No claim that
arbitrary SOURCE multiplication acts on M_P is made.

### Bijectivity uses only the two polynomial Poincare lemmas

For any r in R, choose a polynomial one-form alpha with
d alpha=r*Omega. Since (dP,eta) is a polynomial coframe, write
alpha=a*eta+c*dP. Then

    r=(E+b)(a)+D(c),

so nabla is surjective on M_P. Conversely if (E+b)(a)=D(c), the
polynomial one-form a*eta-c*dP is closed. It is dH for some H in R.
Evaluation on D gives -a=D(H), so [a]=0. Therefore

    nabla:M_P -> M_P is BIJECTIVE                     (3)

for EVERY polynomial submersion P:A2->A1, with NO second component Q.
This is a classical Gauss-Manin/Brieskorn phenomenon; we do not claim
it or the Poincare lemmas as a new classical theorem. The explicit
divergence-corrected proof clarifies the scope of LG1.

For a Keller mate Q with Jac(P,Q)=j in C*, take
E=(Q_y*partial_x-Q_x*partial_y)/j. Then b=0, eta=dQ/j and
D=-j*D_Q in LG1's notation. Thus (3) specializes EXACTLY to LG1's
bijective E on its raw quotient. It does not provide additional
evidence that this quotient must vanish.

## 2. The mark which remembers the missing Keller component

The class

    kappa_P=[div(E)] in M_P

is independent of E, since b'-b=D(h). It obeys the exact identity

    kappa_P=nabla[1].                                (4)

Its vanishing means there is a polynomial h with D(h)=-b. The corrected
lift E+hD is then divergence-free. Conversely any divergence-free lift
gives kappa_P=0. A divergence-free lift has a CLOSED eta, hence an
EXACT polynomial eta=dQ by H1_dR(A2)=0. Since dP wedge eta=Omega,
Jac(P,Q)=1. Conversely a Keller mate gives a divergence-free lift.
Using bijectivity (3) gives the complete equivalences

    kappa_P=0
      iff [1]=0 in M_P
      iff 1 belongs to D(R)
      iff P admits a polynomial Keller mate.          (5)

The sign of a mate can be reversed to change D(Q)=-1 to +1.
These are elementary exact equivalences, not a proof that such a mate
exists for a general submersion.

For an ACTUAL Keller component, the mark already vanishes. JC2 is the
still-open implication

    polynomial submersion P AND [1]=0 in M_P
         ==> M_P=0.                                 (OPEN)

Here the forward conclusion gives generic fibre H1=0, hence an affine
line fibre and degree one for the supplied mate, as in retained LG0/1.
This is a reformulation, not a claim of a new independent sufficient
condition. The useful distinction is that unmarked acyclicity (3)
forgets precisely the extra component data in (5). The uniform ideal-core
and image-multiplier criteria of LG0 retain their full quantified scope.

## 3. LP1: whole-plane geometric nonzero controls in every n>=2

For EVERY integer n>=2 let

    P_n=x+x^2*y^n, u=x*y^n,
    E_n=(1-2*u)*partial_x+(4/n)*y^(n+1)*partial_y.

Direct polynomial identities give

    E_n(P_n)=1,
    b_n=div(E_n)=(2+4/n)*y^n.                         (6)

Thus P_n is a SUBMERSION ON THE WHOLE A2. Its raw quotient M_n has
the canonical BIJECTIVE connection from LP0. Its generic connection is
NONZERO and comes from actual algebraic curves; there is no abstract
module-realization assumption.

### The generic genus and finite monodromy are explicit

For c!=0, the curve C_c={P_n=c} has x!=0 and

    y^n=(c-x)/x^2.

The simple zero at x=c makes this cyclic function-field extension have
degree n and be irreducible: its valuation1 precludes a proper power.
Over C_c the polynomial submersion ensures smoothness. The projection
to Gm_x is an n-sheet cover off x=c with a SINGLE point over c, of
ramification index n. Consequently

    chi(C_c)=n*chi(Gm minus {c})+1=1-n,
    dim H1(C_c,C)=n.                                 (7)

This Euler calculation is ordinary finite-cover additivity. It is not
the degree of a plane Keller map. The smooth projective completion has
gcd(n,2) points above x=0, because the cyclic valuation there has order
-2, and one point above infinity, whose order is -1. Those and only
those points are removed from the affine curve. Therefore

    s=gcd(n,2)+1,
    g=(n-gcd(n,2))/2.                                (8)

The genus is positive for EVERY n>=3 and unbounded. This is a whole
polynomial-plane submersion with no finite critical points, not a
punctured plane source or an unspecified geometric differential module.

After the finite etale base change c=t^n on Gm, the ENTIRE family is
explicitly trivial:

    x=t^n*X, y=t^(-1)*Y,
    X+X^2*Y^n=1.

The inverse is X=x/t^n,Y=t*y. Its monodromy therefore has finite order
dividing n. In particular the generic rank-n connection is geometric,
regular singular and has quasi-unipotent (indeed finite) monodromy at
0 and infinity. This conclusion follows from the displayed finite-cover
trivialization; no blanket semisimplicity theorem for affine-fibre
monodromy is imported.

The zero fibre is the disjoint union

    {x=0}=A1 and {1+x*y^n=0}=Gm.

It is smooth and reducible. Raw relative cohomology, its torsion-free
image and its rational connection must still be distinguished. None is
silently replaced by another in (3) or (7).

### There is NO polynomial Keller mate in ANY of these controls

Give x weight n and y weight -1. P_n is homogeneous of weight n, so
Jac(P_n,-) sends weight w to weight w+1. If Jac(P_n,Q) were a nonzero
constant, the weight -1 component Q_-1 would have that same Jacobian.
Every polynomial of weight -1 is

    Q_-1=y*H(u), H in C[u].

An exact symbolic calculation gives

    Jac(P_n,y*H(u))=(1+2*u)*H(u)+n*u*(1+u)*H'(u).      (9)

For H nonzero of degree d>=0, its highest coefficient is multiplied
by 2+n*d, and the result has degree d+1 in u. It cannot be a nonzero
constant in characteristic zero. H=0 does not help. Thus NO Q exists,
in EVERY n, regardless of Q's degree or the number of weight components.

By (5), [1] and kappa_(P_n) are NONZERO in the raw quotient despite
its bijective connection. The family is NOT a Keller counterexample:
the second component is proved nonexistent. It disproves the weaker
route from whole-plane submersion, geometricity, regular singularities,
finite monodromy, or transverse acyclicity alone to quotient vanishing.
It does NOT disprove the pointed implication in Section2.

## 4. Prior literature and scope

Dimca--Saito, *Algebraic Gauss-Manin Systems and Brieskorn Modules*,
[arXiv:math/9906129v1](https://arxiv.org/html/math/9906129v1), describes the classical Brieskorn quotient
Omega^2/(dP wedge dR), its relation to the Gauss-Manin system, and the
finite-critical-point correction. Its Introduction, Theorem0.3,
Section1.3 and complete Theorem1.7 proof were read for comparison.
In the nonsingular case the positive Koszul pieces in that proof vanish;
our elementary forms computation supplies the needed connection and
bijectivity directly. No general classification, regular-holonomicity
theorem for every raw quotient, or additional external theorem is
adopted here. The primary HTML is cached and pinned. No novelty priority
is claimed for the classical exactness or these elementary submersion
examples. The new use is a fully geometric, whole-plane filter on the
LG1-to-inversion inference, with the missing volume mark explicit.

The exact checker transcribes (1)--(2) in the new family, changes lifts,
the finite-cover trivialization, cyclic-cover Euler/genus arithmetic and
the all-degree leading-coefficient argument (9). It does not replay LG's
registered frame, abstract-module, point-projector or de Rham controls.
The universal proofs above are not inferred from finite samples.

No formal verification, independent whole-proof replay, external review
or full-goal completion is claimed. Full HC4/JC2 and the actual-plane
ordinary trace-dual multiplication/full marked-energy implication remain
OPEN. The next input must use the ACTUAL Keller mate/vanishing volume
mark, or a genuine uniform source landing; mere geometric acyclicity
of a one-component quotient is now explicitly insufficient.
