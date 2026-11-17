# Hamiltonian stabilizers and boundary conditions invisible to finite source jets

2026-11-16. HC4MYLT0/1.

FULL THEOREM TARGET: prove HC4/JC2 or construct an actual counterexample.
CURRENT LAST UNSOLVED IMPLICATION: ACTUAL polynomial plane source forces
Tr(Lambda_L*Lambda_L) subset B, or full marked-energy integrality; OPEN.
External review/outreach is DEFERRED. LS0/1 made PROGRESS by separating
automatic raw p-curvature from the genuinely finite arithmetic condition.
Here we return to the finite normalization and finite trace dual. The
uniform trace-transporter criterion is ALREADY LD1; it is not reproved.

## 1. LT0: every divergence-free preserving field, not only additive ones

Let X=Spec R=A2_C be open in a normal affine surface Z=Spec S with the
same function field L. Assume the source volume Omega=dx wedge dy has

    div_Z(Omega)=sum_i a_i E_i,  a_i>=0,            (1)

where E_i are ALL omitted prime divisors. A codimension-two omitted
set, if present in this general setting, imposes no additional regularity
test on functions or derivations of normal S. In the actual Keller
application S is the FULL finite normalization, and the retained KN0
identity is a_i=e_i-1. Under failure at least one a_i is positive.

Fix the convention

    {H,K}=H_x*K_y-H_y*K_x,
    delta_H=H_y partial_x-H_x partial_y,
    i_(delta_H) Omega=dH.                          (2)

The sign delta_H(K)=-{H,K} has no effect on preservation. Define

    A_H(S)={H in R:delta_H(S) subset S}.            (3)

We prove the EXACT global description

    A_H(S)={H in S: for EVERY i with a_i>0,
                  H-c_i has v_Ei >= a_i+1
                  for SOME scalar c_i in C}.      (4)

For a_i=0 only regularity H in S is required. Constants on intersecting
positive boundary curves must agree because H is a single regular
function; formula(4) does not disregard finite intersection points.

Here is the local proof. At a generic smooth point of E_i choose a
uniformizer u and a separating coordinate v along the curve. Write

    Omega=u^a*w*du wedge dv,  w a unit.

If delta_H is regular, then dH=i_delta Omega has coefficients divisible
by u^a. H cannot have a pole there: a pole of order m>0 would give a
nonzero du term of order -m-1 in characteristic zero. Thus H is regular
at every boundary prime, hence belongs to S by normality and source
regularity. When a>0, the residue of H on the curve must have zero
differential, hence must be a scalar c. If the first nonzero normal
term of H-c has order r>0, its du derivative has order r-1, forcing
r>=a+1. Conversely that inequality makes BOTH coefficients of dH
divisible by u^a and gives a regular delta_H. When a=0 every regular
H works. Intersect the height-one tests for delta_H(s), s in S, to
obtain global preservation, including all singular finite points.

Every divergence-free polynomial SOURCE derivation is Hamiltonian:
its contraction with Omega is a closed polynomial one-form and hence
is dH by the polynomial Poincare lemma. Thus(4) classifies ALL such
preserving source fields. No completeness or locally nilpotent
hypothesis is used.

For an ACTUAL plane Keller map, retained KU1/LO1 give

    omega_S=Lambda*Omega,   (Lambda:Lambda)=S.

If delta_H preserves Lambda, the Leibniz identity
delta_H(s)*lambda=delta_H(s*lambda)-s*delta_H(lambda) makes it preserve
S. Conversely a derivation regular on S preserves regular top forms
on its smooth locus and their reflexive extension. Since Lie_delta
Omega=d(dH)=0, it preserves Lambda. Consequently

    delta_H(Lambda) subset Lambda
        iff H belongs to A_H(S).                   (5)

This extends LO1's equality of additive stabilizers to EVERY
divergence-free source field. It does NOT assert target-field descent.
Taking H to be a polynomial in an integral source coordinate recovers
the credited KX1 replica inequality r_i(ord p+1)>=e_i.

The set A_H(S) is a commutative subalgebra of S and a Poisson
subalgebra. If a_i>0 and H,K belong to it, then locally
H=c+u^(a+1)h and K=d+u^(a+1)k. Their differential wedge has order
at least 2a+1, and division by Omega yields

    v_Ei({H,K})>=a_i+1.                            (6)

Boundary evaluation H->c_i is therefore a NONZERO Poisson trace on
A_H(S): it annihilates the C-linear span of all its brackets and
sends1 to1. In particular

    1 is NOT a finite C-linear sum of brackets in A_H(S). (7)

There is NO pair H,K in A_H(S) with {H,K} in C*. For an actual Keller
map, producing such a pair would force the full inverse, including
ramified boundary. Its existence is not automatic.

## 2. LT1: arbitrary finite source jets still pass

For every boundary divisor let I_i be its prime ideal in S and put

    J=product_i I_i^(a_i+1).                       (8)

Use J=S if there are no omitted divisors. This NONZERO ideal has
zero locus contained in the omitted boundary and lies in A_H(S), by
(4), with all positive-boundary constants zero. In particular JR=R.

Take ANY FINITE set of closed points z_j in the source and ANY finite
orders n_j. If m_j are their maximal ideals in S, then J is comaximal
with intersection_j m_j^n_j. Chinese remainders give a SURJECTION

    J -> direct_sum_j S/m_j^n_j
       = direct_sum_j R/(m_j R)^n_j.               (9)

The equality uses the given open immersion at the source points.
Thus a preserving Hamiltonian can match the prescribed finite jets
of ANY polynomial Hamiltonian at ANY finite set of source points.
To match vector-field jets of order n, match Hamiltonian jets of
order n+1. Approximants may have unbounded degrees as the requested
orders grow. A single polynomial matching an entire infinite jet is
NOT supplied, nor is unrestricted interpolation on curves whose
closures meet the boundary. The same argument DOES interpolate on any
closed subscheme of Z supported entirely in X.

In particular the preserving fields span the source tangent space at
EVERY point. Their R-linear span is EXACTLY Der_C(R). Indeed their
pointwise spanning makes the quotient of this finitely generated
R-module vanish at every closed point, so it is zero. One can choose
FINITELY MANY preserving Hamiltonians for this R-module spanning.
This is NOT an S-module or B-module spanning assertion; coefficients
from R need not preserve the finite lattice Lambda.

There is even simultaneous source integration:

    R = sum_(H in A_H(S)) delta_H(R).              (10)

To prove it, write any polynomial top form r*Omega=d eta. Express
eta=sum_j b_j*i_(delta_Hj)Omega using a finite R-spanning collection
of the fields above. Since d(i_delta Omega)=0, differentiation gives
r=sum_j delta_Hj(b_j). These primitives lie in R, not necessarily S
or A_H(S). Equation(10) does NOT contradict(7), where both arguments
must be preserving Hamiltonians. Nor does a sum over several fields
give LG0's uniform integration for ONE actual Keller lift.

## 3. The same algebra is necessarily non-Noetherian under failure

Suppose one a_i>0 and write e=a_i+1>=2. The boundary evaluation has
maximal ideal

    m={H in A_H(S):H|E_i=0},  A_H(S)/m=C.

Every element of m has valuation at least e, so every element of
m^2 has valuation at least2e. Choose h in J from(8) with valuation
EXACTLY e at E_i. Such an element exists: localizing I_i at its
generic point gives the maximal ideal of the DVR, while each other
I_j becomes a unit ideal there. Choose also s in S with nonconstant
residue on E_i. For EVERY n>=0,

    h*s^n in J subset m,

and their leading coefficients at order e are linearly independent
over C. Hence their classes in m/m^2 are linearly independent. This
cotangent space is infinite dimensional. It follows that m is not
finitely generated and A_H(S) is NOT NOETHERIAN, in particular NOT
finitely generated as a C-algebra.

For the actual Keller normalization this gives a further exact
recognition statement:

    inverse iff A_H(S) is Noetherian
            iff A_H(S) is finitely generated over C. (11)

The inverse case is A_H(S)=R. The failure case has a positive a_i by
KN0 and the preceding proof. Finiteness of S over the Noetherian B
does NOT make this Hamiltonian SUBALGEBRA Noetherian. It is generally
not a B-module and need not contain B. No automatic finite-generation
argument is asserted.

## 4. Explicit globally spanning fields on the retained proper controls

Use LO1's ALREADY PROVED whole-plane enlargements S_(a,q), ALL
a>=2,q>=2. Their source is R=C[x,f], their source volume is
Omega=dx wedge df, and their ONE boundary has volume order a.
Recall ONLY the retained data

    G=x*f^q-1, H=x*G in S,
    v_E(f)=1, v_E(G)=q+a, v_E(H)=a.

Put e=a+1 as a LOCAL VOLUME THRESHOLD, not a ramification index of
an unconstructed target. The three NEW Hamiltonians

    A=f*G,   B=H*G=x*G^2,   C=f^e                (12)

belong to S and all have boundary order at least e. Therefore their
fields preserve the proper S by(4). Direct source differentiation,
with bracket in the source coordinates (x,f), gives

    {A,B}=-G^2*((q+3)*x*f^q-1),
    {A,C}=e*f^(e+q).                              (13)

If f!=0, the second bracket is nonzero. If f=0, the first equals1.
Thus THESE THREE preserving fields span EVERY source tangent space,
although no TWO preserving Hamiltonians can have constant nonzero
bracket, by the positive boundary. This is an explicit global control
of the difference between spanning and a polynomial symplectic frame.

There is no completeness inference: delta_A(f)=-f^(q+1), so for
f(0)!=0 its trajectory obeys f(t)^(-q)=f(0)^(-q)+q*t and escapes at
a finite complex time. C supplies the old additive direction, but A
is not holomorphically complete. The old geometry, normality, affine
chart gluing, exact volume and class group are retained without replay.
NO finite Keller target is supplied for these proper controls.

There is an explicit ALL-ORDERS illustration, stronger than pointwise
jets. Put z=x*f^q and, for EVERY n>=0, define

    L_n=x*G^2*sum_(j=0)^n (j+1)*z^j.

Since v_E(z)=0 and v_E(x)=-q, its boundary order is at least
2(a+q)-q=2a+q>=e. Thus L_n is in A_H(S). The finite geometric-series
identity gives

    L_n-x=x*(-(n+2)*z^(n+1)+(n+1)*z^(n+2)).        (14)

It vanishes along the ENTIRE retained source line f=0 to order at
least q(n+1). Yet x has a boundary pole and is NOT in A_H(S). Hence
arbitrarily deep finite source-line thickenings also fail to detect
this global preservation condition. The degree of L_n is unbounded;
no single preserving polynomial matching x to every order is claimed.

## 5. Literature, meaning and remaining gap

The terminology Poisson trace agrees with Etingof--Schedler,
[*Poisson traces, D-modules, and symplectic resolutions*](https://pmc.ncbi.nlm.nih.gov/articles/PMC5818674/),
Introduction1.1 and Section2, which were read for context. No theorem
is imported. Their finite-type Poisson-variety statements cannot
simply be applied to the non-Noetherian A_H(S). Nor is S automatically
a Poisson order over the symplectic target: at ramification Omega^-1
has poles. Assuming the target Hamiltonian lifts preserve S would
already remove the ramified boundary.

The new output is an EXACT classification of ALL divergence-free
source stabilizers, their unavoidable boundary traces and infinite
cotangent space, together with arbitrary finite-jet interpolation and
global source spanning. It explains why local field matching, finitely
many source jets, pointwise transitivity and several-field source
integration do not produce a finite-lattice contracting action or a
constant-Jacobian pair of preserving Hamiltonians.

The open task is to force an actual GLOBAL finite-algebra constraint
from the polynomial plane and Keller target, such as the original
ordinary trace-pairing landing. No Hamiltonian Noetherianity, global
preserving symplectic frame, source completeness or target-field
descent has been proved automatic. Full HC4/JC2 remains OPEN. All
assurance flags and the published quartic retirement are unchanged.
No registered mathematical checker is replayed; external review DEFERRED.
