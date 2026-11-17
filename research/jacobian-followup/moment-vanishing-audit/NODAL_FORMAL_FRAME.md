# Completing the irreducible nodal-leading class by a full formal frame

Current quantitative status,15 November2026: the [global pure-prefix proof](GLOBAL_TERNARY_CUBIC_PURE_PREFIX.md)
now gives a sufficient1855 pure-only gate and least-prefix interval[10,1855].
Earlier absence-of-global-upper or unquantified-node statements below are
historical. Their earlier class/certificate scopes remain unchanged; the
new proof remains internal and independently unreviewed.

14 November 2026. Internal written proof with exact auxiliary algebra.
No independent specialist replay, external review or formal verification.
The checker does not certify the formal-module argument.

## Theorem

**G3NODALALL20261114.** Every ternary cubic whose homogeneous leading
curve is irreducible nodal has a nonzero positive pure Gaussian moment.
All lower terms are allowed. The affine coordinate changes retain the
full nondegenerate complex Gaussian quadratic form and its means.

The [isotropic-node/conic exclusion](NODAL_CUBIC_REDUCTION.md) and
[coercive-gradient exclusion](COERCIVE_GRADIENT_MOMENTS.md) already handle
the corresponding alternatives. In the remaining nonisotropic-node chart,
use the [complete normal form and seven-point count](NODAL_GAUSSIAN_CRITICAL_PROFILES.md):

    P=T*(uv+h)+a*u^3+d*v^3+b*u^2+c*v^2+e*u+f*v+k,
    a*d!=0,
    Q=(A*u^2+B*v^2+m*T^2+2*Cuv+2*DuT+2*EvT)/2
        +l_u*u+l_v*v+l_T*T+q0,
    m!=0, det Hess Q!=0.                              (1)

The positive-dimensional critical-set case h=e=f=0 is closed by the
[complete critical-line theorem](NODAL_MORSE_DOMINANCE.md). We prove
every finite-critical member below. Finite degeneracies, coincident
critical values and all quartic resonance conditions are included.

No explicit numerical uniform pure prefix follows. General nonhomogeneous
cuspidal-leading and reducible-leading cubics, coefficient schemes and
lower-power GVC remain separate open questions.

## 1. A constant seed would supply a trivial subquotient

Let K=C(t), omega=du dv dT, S_t=tP-Q, and let M be its top twisted
de Rham module, with connection induced by partial_t+P. In (1), each
fixed nonzero t has a linearly coercive complex gradient at infinity.
The retained tame comparison and seven-point intersection computation
give dim_K M=7. This rank statement is inherited; it is not inferred
from the new interpolation examples.

Suppose all positive pure moments vanish. Formal Wick integration with
the full Gaussian Q defines a horizontal functional at t=0 with

    F(omega)=1, F(nabla_t^j omega)=0 for every j>=1.   (2)

It annihilates the twisted relations coefficient by coefficient by
Gaussian integration by parts. It takes coefficients in K to C((t));
poles at0 cause no problem. On the differential submodule U generated
by omega, assign a rational functional L by the values in (2).
If sum_j a_j(t)*nabla_t^j omega=0, applying F gives a_0(t)=0.
Thus L is well defined and horizontal, with nonzero image in the
trivial rank1 differential module K. Consequently

    U is a submodule of M and U -> K is surjective.   (3)

So the trivial module must be a composition factor of M. This argument
needs neither a cyclic seed for M nor an extension to an actual period.
Formal extension to C((1/t)) is exact, so it suffices to exclude that
composition factor at infinity. We now construct the full formal frame,
including its zero-exponential part and descent, rather than judging
only the scalar seed's leading amplitudes.

## 2. Every finite singularity is a one-variable power plus a quadratic pair

At every isolated finite critical point, uv+h=0 and at least one of u,v
is nonzero. If h!=0 this is immediate. If h=0, its center would be
critical only when e=f=0, which is the already closed critical-line case.
Use the chart u!=0, or swap u,v. Put s=uv+h and

    R(u)=G(u,-h/u),
    T'=T+[G(u,(s-h)/u)-G(u,-h/u)]/s,
    G=a*u^3+d*v^3+b*u^2+c*v^2+e*u+f*v+k.

The quotient is analytic, since the numerator is divisible by s. Exactly,

    P=s*T'+R(u),
    det d(u,s,T')/d(u,v,T)=u!=0.                    (4)

If p is critical with P(p)=alpha and R has critical multiplicity mu,
an analytic univariate coordinate y gives

    P=alpha+s*w+y^(mu+1).                          (5)

This is an exact local coordinate identity, not a generic splitting of
the deformation tP-Q. The volume Jacobian is an analytic unit, as is
exp(-Q) at p. For h!=0, the retained sextic gives total finite Milnor
number6; for the other finite profiles it is4,3,2. In particular
1<=mu<=6. Multiplicities and repeated critical values are permitted.

Write d0=mu+1. The local formal period functionals for (5) have powers

    exp(alpha*t)*t^(-nu_j),
    nu_j=1+j/d0,  j=1,...,d0-1.                    (6)

All nu_j are nonintegers. The quadratic pair contributes t^-1;
the remaining exponents are those of a one-variable power, not those
of a presumed Morse perturbation.

Here is an explicit all-order construction. For character j0=j-1 and
monomial y^a s^b w^c, the normalized functional is zero unless

    b=c and a=j0+q*d0 for an integer q>=0.

Otherwise, after factoring exp(alpha*t)*t^(-nu_j), its value is

    (-1)^(b+q)*b!*((j0+1)/d0)_q*t^(-b-q),          (7)

where (z)_q is the rising factorial. An analytic amplitude has a
well-defined formal evaluation: at any order in1/t only finitely many
Taylor coefficients occur. Apply this to the transformed polynomial
form times its analytic volume Jacobian and exp(-Q).

The recurrences

    I_(a+d0)=-((a+1)/d0)*t^-1*I_a,
    J_(b+1,b+1)=-(b+1)*t^-1*J_(b,b)

prove integration by parts for y,s,w. They also prove that differentiation
in t equals insertion of P. Thus these are formal horizontal functionals
on M. The coefficients after (6) are formal series in1/t for every
amplitude whose coefficients are unramified at infinity. No convergence
or full analytic sector decomposition is required here.

## 3. The escaping formal functionals

The retained profiles have n=1,3,4,5 escaping points and finite Milnor
rank F=7-n. Every escape is eventually Morse. They have nonzero
positive-power action: for n=1 a nonzero quadratic term in t; otherwise
the common k*t and a nonzero multiple of K_i^2*t^(2/n), with K_i!=0.
In particular no escape has zero exponential factor. Opposite quartic
roots can share their positive action; this does not change that fact.

For n=3,4,5 put r=t^(-1/n), u=r^(n-2) U, v=r V, T=K/r.
The phase is k*t+r^-2 F_r with nondegenerate limiting Hessian at each
of the n distinct K roots. Formal implicit functions give z_i(r).
Put z=z_i(r)+r*Y. The exponent becomes its stationary action plus
a nondegenerate Gaussian in Y and a cubic perturbation divisible by r.
The original volume factor is r^(n+1). Formal Wick contraction therefore
gives a horizontal functional with form

    exp(q_i(t))*t^(-1-1/n)*(a formal Laurent series in r), (8)

where q_i is the entire positive-power part of the stationary action.
Its leading value on the constant amplitude is nonzero. The stationary
constant and negative-power action terms are absorbed in a formal unit.

For n=1 use u=U/t,v=V/t,T=(h/m)*t+W. The limiting exponent after removal
of its scalar action is a nondegenerate quadratic in U,V,W; the higher
terms are divisible by1/t. Formal Gaussian evaluation gives (8) with
prefactor t^-2. The exact limiting quadratic was checked in the
[finite Morse proof](NODAL_MORSE_DOMINANCE.md), Section4.

Formal Wick integration by parts proves the twisted relations and
parameter horizontality coefficient by coefficient. A parameter-dependent
coordinate change adds only total derivatives, whose formal Gaussian
evaluations vanish. The determinant normalization is retained throughout.
These are full vector-valued formal functionals. Independence of the
scalar seed's paired first, third or later coefficients is unnecessary.

## 4. A polynomial interpolation frame proves completeness

Choose a constant complex linear form L with nonzero T coefficient,
distinct values L_p at the distinct finite critical points, and nonzero
derivative on the univariate direction in (5) at each such point.
These are finitely many proper linear conditions to avoid. Write the
escape values as L_i(r)=kappa_i/r+O(1), with distinct nonzero kappa_i.
For n=1 take r=1/t. Let

    E(X,r)=product_(escapes i) (X-L_i(r)),
    D(X)=product_(finite p) (X-L_p)^mu_p, deg D=F.

Let H_(p,j)(X), 0<=j<mu_p, be the degree<F confluent Hermite basis
modulo D, with prescribed jets delta_(p,j) at the finite nodes. Use
the following seven polynomial amplitudes (their coefficients are
meromorphic formal/algebraic functions of the parameter):

    finite rows: H_(p,j)(L)*[E(L,r)/E(L_p,r)]^8,
    escape row i:
      [D(L)/D(L_i(r))]
          *product_(escape l!=i) (L-L_l(r))/(L_i(r)-L_l(r)). (9)

There are F finite rows and n escape rows. For a formal escape Gaussian,
the fluctuation of L about its critical value is O(1); its separation
from every other escape and from each finite node is order r^-1.
For a finite local model all fluctuations have positive1/t-order.
Consequently the normalized evaluation matrix has these leading blocks:

| Rows / columns | Finite character columns | Escape columns |
|---|---|---|
| Finite rows | Invertible finite jet matrix plus positive orders | O(r^(9-F)) |
| Escape rows | O(r^F) | Nonzero diagonal constants plus O(r) |

For the first off-diagonal estimate, H has degree at most F-1 and
can grow at most as r^(-(F-1)) near an escape. At escape i, one factor
of E(L,r) is O(1), its other n-1 factors are O(r^-1), and E(L_p,r)
has order r^-n. The eighth power supplies r^8. For the second estimate,
the F factors in D(L_i) supply r^-F while D(L) is bounded near a finite
point. The remaining escape ratios stay bounded. Both orders are
strictly positive because F<=6.

In the finite diagonal block, restricting to s=w=0 makes L-L_p an
invertible local coordinate in y. Multiplication by the Jacobian and
exp(-Q) is multiplication by an analytic unit. The resulting jet map
is triangular with diagonal unit*coordinate_derivative^j. It is
invertible. The characters (7) have nonzero normalized leading values
on their corresponding jets. Escape diagonal values are nonzero Gaussian
constants. Thus, after factoring the scalar column factors (6),(8),
the entire seven-by-seven matrix is invertible over a Puiseux Laurent
field, including when positive actions coincide.

Since dim M=7, these seven horizontal formal functionals are a full
frame. In the basis (9), their period matrix is

    H(r) * diag(exp(q_i(t))*t^(-beta_i)), det H(r)!=0. (10)

Here the diagonal includes every finite character (6). Horizontality
and inversion of H identify the formal connection with the corresponding
rank1 factors after ramification. This is a direct construction of the
formal frame; no stationary-phase theorem for an unidentified global
direct-image module is invoked. Distinct critical values are unnecessary.

## 5. The zero-exponential block descends over the base field

This step is essential: an even ramification can turn a noninteger
power into an integer. For example t^-3/2 becomes r^6 when t=r^-4.
One must not infer the existence or absence of a base trivial factor
from those ramified integer exponents alone.

The escape set is invariant under continuation around infinity, being
the set of unbounded critical branches. Hence E(X,r) has coefficients
in the unramified base Laurent field C((epsilon)), epsilon=1/t.
Its constant coefficient has a nonzero simple pole in epsilon; all
other coefficients are bounded. This follows either from the root
orders or symmetric descent: the n leading roots are distinct nonzero
roots of K^n=constant, so their product is nonzero and the other leading
elementary symmetric sums vanish. The n=1 statement is immediate.
Thus E(X,r)/E(L_p,r) is unramified and tends to1 for bounded X.
All finite rows of (9) have coefficients in C[[epsilon]].

Every escaping q_i is nonzero. The q=0 block is therefore exactly the
finite characters at points with alpha=P(p)=0. It is invariant under
ramification descent. Projecting the corresponding finite rows to
this block gives a base-field basis: their pairing matrix is invertible
at leading order by the same jet argument. On these rows and columns,
(7) gives an unramified formal matrix times

    diag(t^(-1-j/(mu_p+1))), j=1,...,mu_p.          (11)

So the zero-exponential block is a direct sum of these regular singular
factors over C((epsilon)), with no integer exponent modulo Z. Its
formal monodromy eigenvalues are exp(-2*pi*i*j/(mu_p+1)); none equals1.
Taking the dual reverses these eigenvalues and still gives no1.
Coincident values and repeated spectral exponents do not create a
new eigenvalue: the full local character frame has already been constructed.

For clarity about canonicity, after ramification the factors with
different q cannot have a nonzero formal horizontal morphism: the
logarithmic derivative of a nonzero Laurent series has at most a
simple logarithmic pole and cannot equal the derivative of a nonzero
negative-power exponential difference. Hence the q=0 projector is
unique, commutes with descent, and the preceding base-field restriction
does identify the entire q=0 component.

The nonzero exponential rank1 factors cannot supply a trivial
composition factor. Neither can (11), whose noninteger residues exclude
the trivial regular factor over the base Laurent field. Direct sums
remain free of that composition factor; exact formal extension preserves
any putative one from (3). This contradiction proves the finite-critical
exclusion. Together with the complete critical-line theorem and the
previous isotropic/coercive alternatives, it proves G3NODALALL20261114.

## Evidence and what this supersedes

The [checker](verify_nodal_formal_frame.py) and
[receipt](nodal_formal_frame_verification.json) verify (4), its coordinate
Jacobian, the sextic derivative, the recurrences and parameter identities
in (7) for d0=2,...,7, general-covariance Wick identities, confluent
Hermite interpolation for all21 multiplicity partitions of F=6,4,3,2,
the degree/valuation bounds in (9), and noninteger base residues.
The universal proofs of those constructions and the formal-module
completion/descent argument are written mathematics, not a machine proof.

The [source audit](nodal_formal_frame_source_audit.json) identifies the
retained rank/critical-profile inputs, the formal Gaussian and power-phase
construction, the elementary-model vocabulary checked against Sabbah,
and the targeted latest-version recheck. No broader external theorem
is imported to bypass the construction or the base-field descent.

The first-amplitude [quartic theorem](NODAL_QUARTIC_RESONANCE.md) remains
valid in its original scope. During this continuation a next odd scalar
coefficient was also computed in the simultaneous gamma=Xi=0 chart.
Its exploratory input/output hashes are retained in the source audit.
It supplies a further open-locus scalar separation, but another vanishing
coefficient would leave another boundary. The full-frame argument now
closes all of these resonances and finite degeneracies at once. The earlier
boundary example with gamma=Xi=0 is therefore excluded by the new theorem;
it was never an all-zero-moment witness.

The subsequent [affine cuspidal theorem](CUSPIDAL_AFFINE_EXCLUSION.md)
closes every cuspidal-leading cubic for GMC. The remaining arbitrary-cubic
frontier is a nonhomogeneous reducible leading curve, including repeated-line
degenerations.
All smooth-leading, irreducible nodal-leading and homogeneous ternary
cubics are covered by the retained theorems. Numerical global pure
prefixes, full coefficient schemes and lower-power GVC remain open.
