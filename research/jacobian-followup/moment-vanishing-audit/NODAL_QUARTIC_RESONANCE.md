# Separating a resonant pair by its first scalar amplitude

14 November 2026. Internal written analytic proof with exact auxiliary
algebra. No independent replay, external review or formal verification.

The subsequent [full formal-frame theorem](NODAL_FORMAL_FRAME.md)
closes every irreducible nodal-leading cubic, including all finite
degeneracies, coincident critical values and quartic resonances. Earlier
method boundaries below describe the scope of this retained argument.

## Result and exact boundary

Use the full affine Gaussian convention

    Q=(A*u^2+B*v^2+m*T^2+2*C*u*v+2*D*u*T+2*E*v*T)/2
        +l_u*u+l_v*v+l_T*T+q0,   m!=0, det Hess Q!=0.

In the quartic escape chart of the [nodal normal form](NODAL_GAUSSIAN_CRITICAL_PROFILES.md),

    P=Tuv+a*u^3+d*v^3+b*u^2+c*v^2+e*u+k,
    a*d*c*e!=0,                                       (1)

assume all affine critical points of P are Morse and their P-values
are pairwise distinct. Impose the resonance left open by the
[finite Morse theorem](NODAL_MORSE_DOMINANCE.md):

    gamma=d*e*m/(2*c)-l_T=0.

Define

    Xi=-4*D*c^4+4*b*c^3*d*m+2*c^3*l_v+c^2*d+d^3*e^2*m. (2)

**G3NODALRESONANCE20261114.** If Xi!=0, some positive pure Gaussian
moment of P is nonzero. The full Gaussian means and cross terms are
retained. For the swapped quartic chart, first exchange u,v and all
corresponding coefficients in P and Q.

Combining this result with the preceding finite Morse theorem, the only
remaining quartic escape boundary on its Morse/distinct-P-value locus
is the simultaneous condition gamma=Xi=0. The latter is not excluded
here. Degenerate finite critical points and coincident affine critical
values remain open, as do other nonhomogeneous singular-leading cubics.
No uniform numerical pure prefix is asserted.

## 1. The exact rescaled phase and its parity

Put r=t^(-1/4), u=r^2 U, v=r V, T=K/r. After removing k*t,

    S_t-k*t=r^(-2) F_r(U,V,K),   S_t=tP-Q,
    F_r=F0+rF1+r^2 F2+r^3 F3+r^4 F4+r^5 F5+r^6 F6,
    F0=KUV+eU+cV^2-mK^2/2,
    F1=dV^3-d*e*m*K/(2*c),
    F2=bU^2-EVK-q0,
    F3=-DKU-l_v V,
    F4=aU^3-BV^2/2-l_u U, F5=-CUV, F6=-AU^2/2.        (3)

The leading critical points are

    z0(K)=(2*c*e/K^2,-e/K,K), K^4=-2*c*e^2/m,
    det Hess F0(z0)=4*m*K^2!=0, F0(z0)=-m*K^2.         (4)

F1(z0)=0 precisely because gamma=0. The exact involution

    F_r(U,-V,-K)=F_(-r)(U,V,K)                         (5)

pairs opposite roots. By the analytic implicit function theorem, each
point in (4) has a unique analytic critical branch z_K(r) nearby, and
z_(-K)(r)=sigma(z_K(-r)), sigma(U,V,K)=(U,-V,-K).
Consequently the stationary values have the same even coefficients
and opposite odd coefficients within each pair. In particular their
common positive action is k*t-m*K^2/r^2; their constant stationary
action term is also the same. A common growing exponential alone
does not distinguish the two periods.

## 2. The first odd scalar coefficient

Let H=Hess F0(z0), g_i=grad F_i(z0), and z1=-H^(-1)g1. Exact reduction
using K^4=-2*c*e^2/m gives

    z1=(-d*e^2/K^3,-d*e^2/(2*c*K^2),-d*e/(2*c)).      (6)

Write the stationary action as

    S_t(z_K(r))=k*t-m*K^2/r^2+c0(K)+c1(K)*r+O(r^2).

Stationarity cancels the unknown second and third displacements from
the coefficient c1. Thus

    c1=F3(z0)+g2^T*z1+(1/2)*z1^T*(Hess F1(z0))*z1
                      +(1/6)*D^3F0[z1,z1,z1]
      =e*(Xi-c^2*d)/(2*K*c^3).                       (7)

For H(r)=Hess F_r(z_K(r)),

    d/dr log det H(r) at0
      =tr(H^(-1)*(Hess F1(z0)+D^3F0[z1]))
      =-d*e/(K*c).                                  (8)

The volume factor in this chart is r^2, and the three-dimensional
Morse scale for r^(-2)F_r is r^3, so each scalar period starts in r^5.
The first ordinary stationary-phase correction starts in r^2 and does
not contribute at order r. Expanding the stationary exponential and
the inverse square root determinant in (7)-(8) therefore gives the
normalized first coefficient

    chi(K)=c1(K)-(1/2)*(d/dr log det H(r) at0)
          =e*Xi/(2*K*c^3), chi(-K)=-chi(K).            (9)

Dropping the determinant correction would give the wrong condition.
The c^2*d summand in (2) comes from that correction.

## 3. Actual full cuts, finite expansion and pair independence

We reuse the polynomial-in-|t| cut carriers and normalized tail bounds
proved in the [finite Morse note](NODAL_MORSE_DOMINANCE.md), Section3.
They depend on m!=0, not on gamma!=0. Hence a logarithmic phase cutoff
L=H0*log|t| makes the normalized tail o(r^N) for any prescribed N by
choosing H0 large enough. This statement applies to the actual flat
rank7 cut basis, not only formal local solutions.

The four limiting points (4) have disjoint fixed spatial neighborhoods
in (U,V,K). The parameter-dependent holomorphic Morse lemma is uniform
in each. The local initial cut is confined to radius
O(|r|*sqrt(log|t|)) about its analytic critical point. In Morse
coordinates, expand the analytic volume density to degree two and
integrate its Taylor terms against the Gaussian. Terms odd in the
Morse coordinates integrate to zero, and the first remaining local
error is O(|r|^2) after normalization, uniformly on a closed small
angular sector. Analytic parameter dependence of the density gives
the order-r determinant contribution (8). Exponentially small Gaussian
tails outside this local region can be made o(r^N).

The Gelfand--Leray slicing argument retained in the preceding note
identifies this local integral with the initial part of the full cut;
the global normalized tail is o(r^N). Nearby or even equal critical
values of the opposite spatial neighborhoods do not merge their local
Milnor generators: the two disjoint Morse balls each supply one.
Cuts can be chosen in a fixed parameter sector avoiding other critical
values, with local generators retained if endpoints agree. No lower
bound on the difference of the paired critical values is used for the
local Morse chart. This is the reason a first-order amplitude expansion
remains available on the resonance.

Choose the orientation constants consistently in a pair. Their nonzero
leading constants can be absorbed into a fixed rescaling of the two
flat generators. The determinant of sigma is +1, so (5) identifies their
leading Morse data and their common constant action c0. There is a
nonzero constant C_K, shared after that normalization, such that

    J_K  =exp(B_K(t))*r^5*C_K*(1+chi(K)*r+o(r)),
    J_-K =exp(B_K(t))*r^5*C_K*(1-chi(K)*r+o(r)),
    B_K(t)=k*t-m*K^2*t^(1/2)+c0(K).                  (10)

Only a finite-order asymptotic expansion is claimed; no convergence of
the full stationary-phase series is needed. The local Morse calculation
is consistent with [Pham, PartII Section1.1, formulas(1.5)'--(1.7),
printed page26](https://numdam.org/item/AST_1985__130__11_0.pdf).
The parameter and global matching estimates are the written arguments
above and in the preceding note, not an extra hypothesis imported from
the source's general discussion of parameters.

If Xi!=0 then chi(K)!=0 for all four roots. For constants A_+,A_-,
the leading term of A_+J_K+A_-J_-K has order r^5 unless
A_++A_-=0. Under that equality its next term is

    exp(B_K(t))*r^6*C_K*chi(K)*(A_+-A_-).

It can be smaller than both orders only when A_+=A_-=0. Equivalently,
the pair sum and difference have nonzero leading orders r^5 and r^6.
This two-step elimination is valid despite equal positive actions.

## 4. Grouped dominance excludes a constant period

Group the four escape generators into the two opposite-root pairs.
Their positive actions k*t-m*K^2*t^(1/2) differ between pairs.
The three finite Morse generators have pairwise distinct linear
actions alpha*t, by hypothesis. No finite action equals an escape
group's positive action, even if alpha=k. Choose a ray on which all
nonzero positive actions and their differences have nonzero real
leading coefficients. The resulting growing groups are strictly
ordered by exponential growth.

Express a putative nonzero constant period in the flat rank7 basis.
For a largest growing finite generator, its nonzero t^(-3/2) leading
amplitude forces its coefficient to vanish. For a largest growing
escape pair, first divide by its exponential times r^5 and then,
after the leading cancellation, by its exponential times r^6.
Equation(10) forces both coefficients to vanish. Lower exponential
groups are negligible after either division. Repeat for every growing
group. All remaining negative-action periods decay exponentially;
an action with no positive power can only be the at-most-one finite
Morse point with alpha=0, whose t^(-3/2) prefactor also decays.
This contradicts a nonzero constant.

The retained formal-Wick/actual-period bridge on the differential
seed submodule now excludes all-zero positive pure moments. No seed
cyclicity assumption or general formal direct-image splitting is added.

## Evidence and the retained boundary

The [exact checker](verify_nodal_quartic_resonance.py) and
[receipt](nodal_quartic_resonance_verification.json) verify the full
chart, parity, leading critical points and determinant, first critical
displacement, coefficients (7)-(9), and direct Taylor-substitution
regressions on representatives of both K^2 classes. They do not certify
the analytic estimates or dominance argument.

A nonempty example on gamma=0 is

    a=c=d=e=1, b=k=0,
    Q=(u^2+v^2-2*T^2)/2-T,
    Xi=-1, det Hess Q=-2.

Its finite critical-value polynomial is
(27*s^2+4)*(27*s-4), squarefree of degree3, so the affine critical
points are Morse with distinct values. The exact coefficient chi is
-1/2 at K=1 and i/2 at K=i. These are checks, not a sampled substitute
for the all-coefficient proof.

Changing only l_v to1/2 yields gamma=Xi=0 while preserving the same
nondegenerate metric and finite Morse/distinct-value properties.
It is an unresolved method-boundary example, not a pure-null witness.
Thus the simultaneous equations cannot be dismissed through Gaussian
nondegeneracy or by continuity. The next gate is a higher odd scalar
coefficient or a justified analysis of the entire equal-action group.

The [source audit](nodal_quartic_resonance_source_audit.json) distinguishes
the exact algebra, the retained cut mechanism and the local Morse import.
General cubic GMC, lower-power GVC and the retained coefficient-scheme
questions remain open.
