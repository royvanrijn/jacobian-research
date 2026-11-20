# One eventual Strong Factorial cutoff for every integer quartic cusp

For every integer s>=1 put a=s^2, d=2s+1, C=s(s+1), L=sqrt(C),
and w_c=C^d*(s+1)/s. Define the actual monic moment polynomial

    A_(m;s)(x)=sum_(k=0)^m binom(m,k)*(am+1)_(dk)*x^(m-k).   (U1)

The rising factorial is (z)_j=z(z+1)...(z+j-1), with (z)_0=1.
These are integer-coefficient polynomials, of degree m and nonzero constant term.

**Theorem UniformIntegerCuspEventualStrongFactorial.** There is ONE integer
N0>=4 such that, for every integer s>=1 and every N>=N0,

    gcd_Q[x](A_(N;s),A_(N-1;s))=1.                         (U2)

The same N0 works for every characteristic-zero field K and every nonzero
binomial f=X^(s^2)*(A+B X^(2s+1)) over K: the two factorial moments
Lambda_K(f^N) and Lambda_K(f^(N-1)) are never both zero at N>=N0.
The cutoff is independent of s, K, A and B. It is not numerically evaluated.

This is an eventual theorem. It does not assert that every terminal order
works, identify the exceptional pairs, or establish the full Strong Factorial
Conjecture for ordinary binomials.

## 1. Exact actual identities

Let Lambda_K(X^j)=j!. For B!=0, direct expansion gives

    Lambda_K(f^m)=B^m*(am)!*A_(m;s)(A/B).                  (U3)

For the original positive real cycle,

    I_(m;s)(w)=integral_0^infinity [S^a(w+S^d)]^m*exp(-mS)dS
             =(am)!*A_(m;s)(m^d*w)/m^((a+d)m+1).          (U4)

Consequently a common actual root x compares the TWO own-order parameters

    q_N=x/(N^d*w_c),
    q_n=x/(n^d*w_c)=(N/n)^d*q_N, n=N-1.                 (U5)

No common-shape Gamma replacement is made.

In normalized coordinates t=(S-C)/L and p=1+t/L, the complete integer
period is, up to the positive nonzero scalar
L*[C^(a+d)*exp(-C)]^m,

    J_(s,m)(q)=integral_(-L)^infinity Psi_s(t,q)^m dt,
    Psi_s(t,q)=p^a*exp(-Lt)*[(s+1)q/s+p^d].              (U6)

This entire-function normalization remains meaningful at q=-1. Dividing
globally by q+1 would be invalid. Local logarithms and a central normalization
are used only near q=1 or near nonzero actual saddle factors.

## 2. Uniform local disc and regular complement

The two accompanying complete proofs are the new analytic inputs.

1. UNIFORM_GROWING_QUARTIC_CUSP_NEIGHBORHOOD.md, Theorem
   UniformPhysicalCuspNeighborhood, gives r>0 and N_local, BOTH independent
   of s, excluding every actual common root with |q_N-1|<r once N>=N_local.
   It uses the compact 402 phase, complete tails and fixed bounded-H theorem.
   Its growing-field argument has uniform actual packets and Cauchy collars.
   A separate positive-real one-packet catch first puts the LOWER endpoint
   in the local disc. Both endpoint cone conditions then keep the entire
   actual radial path central; no assumption m/s small, or d*rho/k small,
   is substituted. Its normalized error/main ratio is
   O(1/N+rho/k+1/k^2).

2. UNIFORM_CUSP_CURIE_WEISS_REGULAR_COMPLEMENT.md, statement (G21), gives
   for each fixed r>0 one N_regular(r), independent of s, excluding actual
   common roots with |q_N-1|>=r once N>=N_regular(r).
   Its only bounded-s antecedent is the independent narrow 405 fixed-cusp
   eventual theorem. It does NOT require the 407 all-support theorem.
   For s tending to infinity, the exact entire phase converges on fixed
   t-rectangles to

       Psi_0(t,q)=exp(-t^2/2)*[q*exp(-t)+exp(t)].          (U7)

   The proved sharp Gamma annulus puts all actual q-roots in
   1/16<|q|<8. Complete real tails are controlled by an exact derivative
   inequality, including the receding endpoint -L. Off the limiting unit
   circle a dominant-saddle horizontal has one strict simple global maximum.
   On the circle away q=1 it has two strict simple maxima on one horizontal.
   Those strict finite-contour gaps transfer to the ACTUAL Psi_s and its
   ACTUAL critical values. The proof never approximates exp(m*Phi_s) by
   exp(m*Phi_0), so m/s is unrestricted.

   Individual actual roots in the joint s,m->infinity regime are captured
   asymptotically by the unit circle, not exactly on it at finite s. At a
   putative common-root sequence both own-order radii tend to1, forcing
   (N/n)^d->1 before the adjacent path is used. Two actual packets then give
   an exact envelope derivative and a strict log-modulus drop on that path.

The conclusions in these two inputs overlap and together cover every complex
q. Choose

    N0>=max(4,N_local,N_regular(r)).                       (U8)

This one choice is made before s or the coefficients are varied. A common
root at N>=N0 would either satisfy |q_N-1|<r and violate Input1, or satisfy
|q_N-1|>=r and violate Input2. Thus the actual polynomials have no common
complex root. A nonconstant rational gcd would have a complex root, so (U2)
follows. No limiting-zero approximation is used as an all-root substitute.

## 3. Coefficients and arbitrary characteristic-zero fields

For each s and N>=N0, the Euclidean algorithm over Q provides polynomials
U,V in Q[x] such that

    U*A_(N;s)+V*A_(N-1;s)=1.                              (U9)

For any characteristic-zero field K its prime-subfield injection Q->K
preserves this identity. If B!=0 evaluate at A/B and use (U3). The two moments
cannot vanish. If B=0,A!=0, each positive moment is A^m*(am)!, hence nonzero;
A=0,B!=0 is also a nonzero monomial. This proves the uniform coefficient and
field assertion without embedding K into C.

## 4. Genuine reverse equivalence and the unevaluated finite box

Define the GLOBAL bad set

    BadCusp={(s,N): s>=1,N>=2 integers,
                     gcd_Q[x](A_(N;s),A_(N-1;s))!=1}.       (U10)

The following four properties are equivalent.

(i) BadCusp is finite.

(ii) There is one N0 such that (U2) holds for EVERY s and N>=N0.

(iii) There is one N0 such that over EVERY characteristic-zero field and for
every s and nonzero cusp coefficient pair, the consecutive moments do not
both vanish for N>=N0.

(iv) The same uniform coefficient assertion holds over C.

The theorem proves all four. Here are the reverse implications, including
the extra input needed for the global finiteness assertion.

If (i) holds, take N0=max(4,1+max{N:(s,N) in BadCusp}), using N0=4 if the
set is empty. This proves (ii). Property (ii) implies (iii) by (U9), and
(iii) implies (iv). Conversely, a nonconstant rational gcd at some N>=N0
has a complex root z; f=X^a*(z+X^d) then contradicts (iv) by (U3). Thus
(iv) implies (ii).

To deduce (i) from (ii), a uniform order bound alone does not suffice.
Use the sharp soft criterion of
GAMMA_PRODUCT_SHARP_LOWER_ANNULUS_AND_SHAPE_SEPARATION.md:
actual consecutive coprimality holds whenever

    a>=d(N-2).                                            (U11)

At a bad cusp pair this yields

    s^2<(2s+1)(N-2).                                      (U12)

In particular N>=3 and

    s<(N-2)+sqrt((N-2)(N-1))<2N-3.                       (U13)

Together with N<N0, this puts EVERY bad pair in the finite rectangle

    3<=N<N0, 1<=s<=2N0-6.                                (U14)

This is the explicit SHAPE of an unevaluated box, not an evaluated
computational box: N0 is unknown. Therefore (ii) implies (i), completing
all directions of the equivalence under the retained sharp soft criterion.

Finiteness of BadTerminal_s at each fixed s would be weaker: for example the
abstract set {(s,s+2):s>=1} has one point at every support and no uniform
terminal cutoff. No fixed-support quantifier exchange is being assumed here.
Nor does finiteness of BadCusp imply that it is empty. An all-order assertion
would require precisely that stronger emptiness statement.

## 5. Logical antecedents and scope

The source ledger separates mathematical proof antecedents, finite exact
identity support, and classical context.

- The compact normalized tails, the bounded-H theorem, and their model/jet
  mechanisms are the retained 402 sources.
- The fixed C1--C15 complete construction is revalidated uniformly in the
  new local proof; its old constants are not simply declared uniform.
- The bounded-s fallback is the narrow fixed-cusp 405 route with its local
  and global complete-contour proofs. No broad square-root-gap theorem,
  407 theorem, finite terminal prefix, or lower-one finite certificate is
  a premise.
- The old Gamma soft criterion bounds the lower-endpoint radial factor in
  the local catch. The sharp 404 annulus supplies the global compact q-set;
  the sharp soft criterion gives (U12)--(U14).
- The exact Curie--Weiss/unitary-Hermite correspondence is classical context.
  It is NOT a premise for an exact Gamma circle law.

The limiting polynomial is
sum_k binom(m,k)q^(m-k)exp[(2k-m)^2/(2m)].
In Kabluchko's convention (1.5) it equals
(-1)^m exp(m/2) H_m(-q;4/m), the critical beta=1 model in (1.8)--(1.9).
The earlier scratch GAMMA_CURIE_WEISS_CRITICAL_BRIDGE.md already recorded
this correspondence and why relative coefficient transfer is nonuniform
in m/d. The present proof instead transfers strict contour geometry while
retaining the actual exponent. See UNIFORM_CUSP_PRIMARY_CONTEXT.md.

No cutoff, finite bad set, or numerical bound on its elements is computed.
No all-order cusp theorem, full ordinary-binomial SFC, effective all-support
cutoff, or exact finite-Gamma Lee--Yang circle law is asserted. The checker
verifies exact algebra, moment normalizations, scale identities, explicit
negative controls and source bytes. It does not replace analytic estimates
by finite samples. All independent-replay, external-review and
formal-verification assurance flags remain false.

## Reproduction

Run from any working directory:

    python3 PATH/verify_uniform_cusp_eventual.py --repo-root REPO --output RECEIPT

The checker is standalone Python using only the standard library. Repository
proof imports are resolved through --repo-root; local files use the checker
directory. Its receipt records every actual source hash and the checker
ledger hash. Optimized Python execution is refused.

