# Eventual squarefreeness at every fixed ordinary binomial support

This consequence uses the individual complete periods and complex collars;
it does not infer simplicity from consecutive coprimality.

For fixed integers a>=0,d>=1 put

    A_m(x)=sum_(k=0)^m binom(m,k)*(am+1)_(dk)*x^(m-k).

**Theorem.** There is an unevaluated integer m_0(a,d) such that A_m has
only simple complex roots for every m>=m_0(a,d), equivalently
gcd_(Q[x])(A_m,A_m')=1, also after every characteristic-zero field
extension. No effective or support-uniform cutoff, nor all-order claim
for a>0,d>1, is supplied.

## 1. Exact derivative criterion for complete packets

For a>0 let D=a+d and

    I_m(w)=integral_0^infinity [S^a*(w+S^d)*exp(-S)]^m dS
          =(am)! A_m(m^d*w)/m^(Dm+1).

The scalar and affine substitution are nonzero, so multiplicities are
unchanged. Write Phi=log[S^a*(w+S^d)*exp(-S)] locally. At two
distinct stationary points S_L,S_R, set DeltaPhi=Phi_R-Phi_L. The exact
stationary envelope identity is

    d*w*DeltaPhi_w=-(S_R-S_L).                         (1)

On a controlled two-packet neighborhood the ACTUAL integer period is

    I_m=n_L P_L+n_R P_R,
    P_i=exp(m Phi_i)*sqrt(2*pi/m)*c_i*(1+R_i),

where the active tree coefficients are fixed signs and c_i are analytic
and nonzero. A single active packet is zero-free at large order. If both
are active, choose local logarithms and put

    E=m*DeltaPhi+log(c_R/c_L)+log(1+R_R)-log(1+R_L).

Exactly I_m=n_L P_L*(1+(n_R/n_L)*exp(E)). At a zero,

    partial_w I_m=-n_L P_L*E_w.                       (2)

Thus E_w!=0 proves simplicity. Local logarithm branches do not affect
w-derivatives. No noninteger original-cycle identity is required.

## 2. Regular parameters, including regular Stokes walls

The complete regular proofs give Gaussian bounds on a larger complex
w-disc with fixed end labels. The exact R_i are holomorphic there.
Cauchy's estimate sends the uniform R_i=O(1/m) bound to
R_(i,w)=O(1/m) on a smaller fixed disc. Thus E_w=m*DeltaPhi_w+O(1).
By (1), its leading derivative is nonzero at every regular two-saddle
parameter with w!=0. Shrinking the neighborhood bounds its modulus
away from zero, proving (2) nonzero at all large orders.

This argument includes equal-height Stokes parameters. Disjoint local
Morse discs have a fixed exit drop. When the critical heights are closer
than one third of that drop, every descending exit is below BOTH critical
heights. The complete arms therefore cannot end at either saddle; their
fixed labeled end caps and uniform Gaussian bounds persist across the
wall. Unequal-height Stokes points are zero-free by the finite-sublevel
detour. On the positive boundary the strict-subcritical coexistence point
has TWO distinct nondegenerate maxima: split the actual ray at a fixed
point in their intervening valley. Both exact partial periods have the
holomorphic Gaussian/Cauchy bounds above, so (1)-(2) apply. Other positive
points and w=0 are one-packet zero-free. The negative boundary also uses
(1)-(2) at its two maxima. Inactive spinodal inflections are dominated.

## 3. Growing scaled fields at the only coalescences

At a simple turn use its analytic coordinate eta, with eta_w!=0, and set
rho=|eta|^(1/2), k=m*rho^3. On each selected two-packet cone the complete
singular Gaussian estimates give

    |S_R-S_L|>=c*rho,
    |partial_w log(c_R/c_L)|<=C/rho,
    |partial_w log[(1+R_R)/(1+R_L)]|<=C/(k*rho^2).

The Cauchy collar has radius c*rho^2 with CENTER rho fixed; no unspecified
O-term is differentiated. Near this nonzero turning value (1) gives
|m*DeltaPhi_w|>=c*m*rho, and the relative error is at most

    C*(rho/k+1/k^2).                                  (3)

Choose fixed K large and the physical neighborhood small. For k>=K
this ratio is below one half. Complementary angular cones are zero-free
by their complete one-packet estimates and scaled Stokes detours.

At a quartic cusp, eta_w!=0, rho=|eta|^(1/3), k=m*rho^4. The selected
two-packet cone instead has leading amplitude derivative O(rho^-2) and
remainder derivative O(1/(k*rho^3)), on a w-collar c*rho^3. The saddle
separation is again bounded below by c*rho. Dividing by m*rho gives
EXACTLY (3). The same choice proves simplicity for k>=K; the remaining
outer-angle regions are zero-free. Constants may depend on the support.
Only the modulus of saddle separation is needed here.

## 4. Bounded scaled fields and exhaustion of all roots

For k<=K the complete normalized periods converge uniformly on compact
scaled-field sets to a nonzero Airy solution or the quartic model Z(h,0).
Their zeros are simple. On a slightly larger compact collar choose
disjoint discs about its finitely many model zeros. Rouche gives ONE
actual root COUNTED WITH MULTIPLICITY in each disc and no other roots
on the smaller compact set, hence simplicity. The collar handles model
zeros on an originally chosen boundary.

The regular proofs cover strict subcritical supports; their stationary
inflections are strictly inactive and have zero-free neighborhoods.
Equality 4a=(d-1)^2 is the integer quartic cusp. Above-critical supports
have just the two conjugate simple turns. The preceding arguments cover
each physical coalescence and every compact regular complement.
Finally the OLD Gamma upper bound contains every root in the fixed disc

    |x/m^d|<=(a+d)^d.

A finite cover and maximum of its cutoffs prove the theorem for a>0,d>1.
No lower annulus or finite-order enumeration is needed.

## 5. Elementary boundary supports and exact source scope

For d=1, A_m=G_(m,am+1) satisfies

    xG''+(1-gamma-m-x)G'+mG=0, gamma=am+1>0.

Its constant term is nonzero. A repeated nonzero root contradicts the
lowest-order term in this ODE, so all orders are simple. For a=0 the
exact Appell identity A_m'=m*A_(m-1) is available. For d>=6, the retained
constant-plus-power analytic theorem gives simplicity once
m>=max(10000,d^2). For d=2,3,4,5 the prime-power theorem gives all-order
simplicity; d=1 was handled by the ODE. No all-exponent finite prefix is
a premise. Its independently known all-order a=0 theorem is context only.

The premises are the complete estimates, collars and labeled cycles in
ABOVE_CRITICAL_REGULAR_COMPACT_COMPLEMENT.md,
SIMPLE_TURN_SHRINKING_WINDOW_CONSECUTIVE_EXCLUSION.md,
CUSP_SHRINKING_WINDOW_LOG_RATIO_NONCOLLISION.md, the strict-subcritical
proof, and the bounded quartic-cusp/Lee-Yang model proof. Root capture
uses GAMMA_PRODUCT_ROOT_MODULUS.tex; a=0 uses
the analytic part of ALL_UNIVARIATE_TWO_TERM_SFC.tex and the small cases
of PRIME_POWER_TWO_TERM_SFC.tex. Eventual consecutive coprimality alone
is not used to infer squarefreeness. The no-vertical-pair lemma is not a
logical premise: only |S_R-S_L|>0, not its real part, enters (1)-(3).

Reproduce the exact supporting controls with
`python3 verify_ordinary_binomial_eventual_squarefree.py --repo-root REPO
--support-root SUPPORT --output review.json`. SUPPORT contains the two
new regular/simple-turn proofs; in a combined installation it may be this
directory. The guard authenticates all sixteen required sources, the
packet derivative algebra, both Laurent error reductions, Gamma ODE and
Appell coefficient identities, literal small rational polynomials, and
false-scope controls. These finite checks support the written proof;
they are not formal verification or a finite-prefix argument.
