# Ordinary Strong Factorial for large lower exponents at each fixed gap

Every fixed exponent gap admits an all-terminal Strong Factorial theorem for sufficiently large lower exponents. The constants below depend on the FIXED integer gap d; no uniform C*d^2 threshold is asserted.

## 1. An iterated circular-domain root lemma

Let N>=1, h=N-1 and gamma_1,...,gamma_d>=h be positive reals. For each gamma define

G_(N,gamma)(x)=sum_(k=0)^N binom(N,k)*(gamma)_k*x^(N-k).

Its Laguerre tridiagonal matrix gives all roots beta in both disks

|beta|<=gamma+3h, |beta+2gamma|<=gamma+3h,

and Re beta<=-gamma+h. These are the same exact matrix bounds as in [the retained gap-two localization](GAP_TWO_GRACE_GLOBAL_LOCALIZATION.md), independent of any positive-measure orthogonality.

Let C_r be the unsigned binomial Grace convolution of the first r factors. Define M_r>0 and nonnegative E_r,U_r recursively by

M_1=gamma_1, E_1=3h, U_1=h,
M_(r+1)=M_r*gamma_(r+1),
F_r=E_r+2U_r,
E_(r+1)=3M_r*h+gamma_(r+1)*E_r+3h*E_r,
U_(r+1)=5M_r*h+gamma_(r+1)*F_r+3h*F_r.                    (G1)

Every root zeta of C_r satisfies |zeta|<=M_r+E_r and Re zeta<=-M_r+U_r.

Indeed, these two bounds imply a second containing disk

|zeta+2M_r|^2<=(M_r+E_r)^2+4M_r*U_r<=(M_r+E_r+2U_r)^2.

Apply the exact closed-disc Grace representation zeta_new=-alpha*beta first with the origin disk and then with the disk centered at -2M_r. The first gives the stated E recurrence. The second gives

Re zeta_new<=2M_r*(-gamma+h)+(M_r+F_r)*(gamma+3h)
           =-M_r*gamma+5M_r*h+gamma*F_r+3h*F_r,

which is exactly the U recurrence. The two applications need not use the same beta. No such simultaneous representation is assumed.

## 2. Explicit all-order strip for an arbitrary fixed gap

Now let d>=2 be fixed and a>=d, and set gamma_j=(aN+j)/d. Every gamma_j>=h. The actual reversed monic moment polynomial is

A_N(x)=sum_(k=0)^N binom(N,k)*(aN+1)_(d*k)*x^(N-k)
      =d^(dN)*C_d(x/d^d).                                (G2)

This is the exact rising-factorial multiplication formula, not an asymptotic approximation. A root x=d^d*zeta has common-N-scale coefficient w=x/N^d. Normalize

w=a^d*(-1+delta*z), delta=d/sqrt(a).

For explicit bounds define integers

e_1=3, u_1=1,
e_(r+1)=3+3e_r, u_(r+1)=5+3e_r+6u_r.                     (G3)

Since G=2a/d>=2, gamma_j<=G*N and h<=N, induction in(G1) gives

E_r<=e_r*G^(r-1)*N^r, U_r<=u_r*G^(r-1)*N^r,
(a/d)^r*N^r<=M_r<=G^r*N^r.

For example (e_r,u_r), r1through6, are (3,1),(12,20),(39,161),(120,1088),(363,6893),(1092,42452).

For a>2^d*d*u_d one has U_d/M_d<1/2. In particular the negative real-part bound can be squared safely. Since

0<=d^d*M_d/N^d-a^d<=(a+d)^d-a^d<=d^2*(2a)^(d-1),

both sides of the real strip give the literal estimate

|Re z|<=R_d/sqrt(a), R_d=2^(d-1)*(d+max(e_d,u_d)).          (G4)

Also

|Im zeta|^2<=(E_d+U_d)*(2M_d+E_d).

The preceding bounds imply the uniform-in-N height estimate

|Im z|^2<=T_d^2,
T_d^2=2^(2d-2)*(e_d+u_d)*(4+e_d)/d.                       (G5)

Thus every upper-order moment root lies in a fixed bounded-height strip whose real width tends to zero as a grows, for every fixed d. These estimates allow different constants for each d; they do not provide a d-uniform threshold.

## 3. The actual fixed-d phase has the same Hermite endpoint

Use v=(s-a)/sqrt(a), epsilon=1/sqrt(a)=delta/d and

h_delta(v)=((1+delta*v/d)^d-1)/delta.

For fixed d this is a polynomial family analytic at delta0, with h_0(v)=v. The actual normalized phase is

F_delta(v,z)=d^2*delta^(-2)*[log(1+delta*v/d)-delta*v/d]
             +log(z+h_delta(v)).                         (G6)

Its removable delta0 limit is the same F_0=-v^2/2+log(z+v). The exact stationary polynomial is z*v+(v-delta)*h_delta(v)-1. Hence the limiting endpoints remain (v,z)=(-i,2i) and its conjugate, with F_vvv=2i and F_vz=1. All local normal-form, leading-amplitude, Airy-contour envelope, divided-amplitude and Robin-band arguments in [the complete gap-two endpoint proof](GAP_TWO_LARGE_LOWER_EXPONENT_ENDPOINT_THEOREM.md) apply to this bounded holomorphic family. Their constants may depend on d, but never on N or on a after delta has been restricted to a fixed small interval.

The complementary complete contours also remain available on the bounded height interval supplied by(G5). On the right real tail,

h_delta'/h_delta<=(1/v+delta)/(1+epsilon*v), h_delta>=v.

On the left tail v=-r, 0<=1-epsilon*r<=1,

0<=h_delta'<=1,
|h_delta(-r)|>=(1-exp(-delta*r))/delta.

Here is an explicit order of choices for the whole tails. First require |Re z|<=1/100 and choose a fixed R>=3 depending only on T_d. For v>=R, h>=v and

Re F_v<=[-v+(11/10)*(1/v+delta)]/(1+epsilon*v)<-1

when delta is sufficiently small; multiply by the positive denominator to verify the sign uniformly for v>=R. On the left real tail, delta*R<=1 gives |h(-r)|>=R/2 for every r>=R, and therefore

Re F_v(-r)>=r-1/(R/2-1/100)>1.

Choose R with -R^2/2+log(T_d+R+3)<-4; for example any sufficiently large integer above4*(T_d+3) works. At the real endpoints the limiting phase is below that bound. The far vertical heights are at most1+small: the bulk saddle heights tend to -y/2 with |y|<=2, and the stable outside saddle height has modulus below1. Hence the limiting vertical phase is bounded by -(R^2-(1+small)^2)/2+log(T_d+R+3), also strictly below the peak minimum -1/2 after enlarging R if needed. Finally choose delta small so the Taylor remainders of U_delta and h_delta on this fixed compact set preserve those strict margins. The derivative inequalities just proved integrate the remaining unbounded real tails. This gives complete nonlocal estimates for both bulk packets and the outside one-packet contour.

Extra stationary roots of the degree-(d+1) equation are not ignored: these estimates control the ENTIRE original contour, so no contribution from them is left unassigned.

## 4. The coupled actual orders retain the same transverse drift

The [retained exact linear lower-exponent criterion](GAMMA_PRODUCT_ROOT_MODULUS.tex) covers a>=d*(2N-3). In its complement,

N>a/(2d)+3/2=d/(2delta^2)+3/2.

For n=N-1, kappa=n/N, the actual identity is

J_N(w)=kappa^((a+d)*n+1)*I_n(kappa^(-d)*w),
z_n=z+alpha*(z-1/delta), alpha=kappa^(-d)-1>=d/N.          (G7)

For fixed d, alpha=O_d(1/N)=O_d(delta^2), so the segment from z to z_n remains inside the chosen endpoint polydisc for sufficiently small delta. In the canonical endpoint coordinate, |zeta_z+i|<=1/4. The same exact estimate as in the gap-two proof gives

Im[zeta(z_n)-zeta(z)]>=alpha/(2delta)
                    >=d/(2delta*N)=sqrt(a)/(2N).         (G8)

The two actual zero bands are O_d(1/N), while(G8) is larger by sqrt(a). They cannot both hold for sufficiently large a. Away from the endpoint windows the actual action-drop argument and complete Gaussian contours are unchanged, with fixed-d constants.

## 5. All-terminal fixed-gap conclusion

**Theorem FixedGapLargeLowerExponent.** For every fixed integer d>=2 there exists an integer a_0(d) such that for every integer a>=a_0(d), every N>=2, and every characteristic-zero field, the actual consecutive ordinary factorial moments of X^a(A+B*X^d) cannot both vanish unless A=B=0.

Proof. The all-root bounds(G4)--(G5) localize every potential upper-order root. The fixed-d analytic phase(G6), complete complementary contours, and the same compact endpoint proof yield the zero bands. The actual shift(G7)--(G8) excludes both bands once a is large. The already retained large-lower boundary covers all other N. Nonzero rational resultants transfer the complex assertion to every characteristic-zero field.

This conclusion is qualitative: no value, effectivity, or d-uniform growth rate of a_0(d) is asserted. For each d the finitely many lower exponents below a_0(d) remain outside this theorem, each still carrying an all-terminal obligation. It is neither a proof of full ordinary binomial SFC nor a finite global support reduction.


## 6. A genuine finite-exception equivalence

**Corollary FiniteBadLowerExponents.**

For each fixed integer d>=2, define the actual reversed ordinary
factorial polynomial
\[
 A_{m;a,d}(x)=\sum_{k=0}^{m}\binom mk(am+1)_{dk}x^{m-k},
\]
and let
\[
 \operatorname{Bad}_d=
 \{a\in\mathbb Z_{\ge0}:\ \exists\,N\ge2,\
       \gcd_{\mathbb Q[x]}(A_{N;a,d},A_{N-1;a,d})\ne1\}.
\]
Then Bad_d is FINITE if and only if there exists an integer
a_0(d)>=0 such that every integer a>=a_0(d) and every N>=2
has actual consecutive coprimality.

If the tail holds, Bad_d is contained in the finite set
{0,...,a_0(d)-1}. Conversely, if Bad_d is finite, choose a_0(d)=0
when it is empty, and a_0(d)=1+max Bad_d otherwise. This proves
both directions. The coprimality formulation is equivalent to the
ordinary moment statement: all coefficients are rational, the
polynomials have nonzero constant and leading terms, and a common
complex root gives a nonzero two-term phase; monomial cases never
have zero positive factorial moments.

Consequently the quantified fixed-gap tail theorem
\[
 \forall\,d\ge2\ \exists\,a_0(d)\ \forall\,a\ge a_0(d)\
       \forall\,N\ge2\quad\text{actual consecutive coprimality}
\]
is equivalent to
\[
 \forall\,d\ge2\quad \operatorname{Bad}_d\text{ is finite}.
\]
Both are proved by the theorem above. Neither statement says
Bad_d is EMPTY or classifies its possible members. No bound on those
members or on any bad terminal orders is supplied, and no finite
global computational box follows. Full ordinary binomial Strong
Factorial remains open.

## Exact reproduction and source scope

    python3 research/factorial/image-rigidity/verify_fixed_gap_large_lower_extension.py \
      --repo-root . --output /tmp/fixed_gap_large_lower_extension_review.json

The exact guards verify the universal circular-domain closure algebra,
coefficient recurrences and literal rising-factorial multiplication,
majorant and scaling identities, together with the stated sign and
scope controls. The all-parameter induction, compact contour bounds
and qualitative analytic threshold are the written proof. The fact
that the displayed majorants depend on d does NOT disprove a possible
stronger d-uniform threshold theorem.

The mandatory ledger pins all TWELVE listed proof/checker/context
inputs, including the complete prior endpoint method. This is a
checked reuse after reverifying every fixed-d phase, contour and
coupled-order hypothesis; a gap-two theorem is not simply applied
to a different exponent gap. No additional classical theorem is
imported beyond the precisely scoped prior sources.

Updated 20 November 2026.
