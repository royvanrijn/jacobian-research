# Constant perturbations away from critical dyadic radii

This note proves a valued-field transfer of all-positive-start moment
windows away from the critical dyadic radii. The statement is local in the
specified valuation; it does not assert unrestricted quadratic SFC.

Let K be a characteristic-zero field with a real-valued nonarchimedean
valuation v extending the 2-adic valuation, with v(2)=1. Write O for its valuation ring, m for its
maximal ideal, and F=O/m for its characteristic-two residue field. A finite
extension of Q2 suffices. Let Lambda be the factorial functional in any
fixed number of variables. Let q have zero constant term and suppose

             R_r = Lambda(q^r)/r! belongs to O for every r>=0.
             R_0=1.                                      (1)

For k in K*, put delta=v(k), eta=1-delta, and

             M_m=Lambda((q+k)^m),   m>=0.

Call a coefficient R_r a unit when v(R_r)=0. A positive unit window of
length w means that for every n>=1 at least one of R_n,...,R_(n+w-1)
is a unit. A moment window refers to every positive starting order.

## 1. Exact integral coefficients, without an analytic integral

One sufficient condition for (1), valid beyond quadratics, is

       q=sum_(alpha!=0) b_alpha X^alpha/alpha!,  b_alpha in O,     (2)

with finite support. In the expansion of Lambda(q^r)/r!, the coefficient
of product b_alpha^(h_alpha), where sum h_alpha=r, is

       product_i (sum_alpha h_alpha alpha_i)!
       / [ product_alpha h_alpha!
           product_(alpha,i) (alpha_i!)^(h_alpha) ].              (3)

This is an integer. To see it, take separately labelled sets for each
variable colour, of the indicated total sizes, and divide their union
into unordered blocks of the prescribed vector sizes alpha. Internal
orders contribute the factors alpha_i!, and permutations of blocks of
the same type contribute h_alpha!. Every block is nonempty because
alpha!=0, so these permutations act freely. The displayed quotient is
the number of such coloured set partitions. Formula (3) therefore proves
(1) coefficient by coefficient. It involves no convergence assertion.

The zero-constant hypothesis in this argument matters: q=1 would give
R_r=1/r!, which is not integral in general.

For a binary quadratic in natural coordinates,

       q=(alpha/2)X^2+beta XY+(gamma/2)Y^2+dX+eY,                (4)

condition (1) follows whenever alpha,beta,gamma,d,e belong to O.

## 2. The unique minimizing index in an open band

Ordinary binomial expansion gives the exact finite identity

       M_m=sum_(r=0)^m (m)_r k^(m-r) R_r,
       (m)_r=m!/(m-r)!.                                      (5)

For j=m-r, Legendre's elementary formula v(n!)=n-s2(n) gives

       v((m)_r k^(m-r))
            =m-s2(m)+s2(j)-eta j.                            (6)

Here s2 is the binary digit sum, and the formula also covers r=0,m.
The additional factor R_r has nonnegative valuation by (1).

Fix a>=1 and L=2^a. Suppose

       1-2/L < delta < 1-1/L,
       equivalently 1 < eta L < 2.                           (7)

Write m=L u+t with 0<=t<L. For j=L b+c, 0<=c<L, the varying part of
(6) is

       s2(j)-eta j
         =[s2(b)-eta L b]+[s2(c)-eta c].                     (8)

For c>0, the second bracket is strictly positive: each nonzero low bit
i<a contributes 1-eta 2^i>0. It is zero at c=0. The first bracket is
strictly decreasing with b, since

       s2(b+1)-s2(b)-eta L <=1-eta L<0.                       (9)

Under j<=m one has b<=u. Thus (8) has its unique minimum at b=u,c=0,
namely j=L u, or equivalently

                       r=t=m mod L.                         (10)

There is no missing leading-bit, endpoint-index, or carry case: for
b<u the first bracket is strictly larger and the second is nonnegative;
for b=u,c>0 the low bracket is strictly larger. The candidate c=0 is
always allowed, including t=0 and m<L.

Consequently, if R_t is a unit, the r=t term in (5) is the unique term
of least ACTUAL valuation, so M_m!=0. In fact

       v(M_m)=delta(m-t)+t-s2(t)  when R_t is a unit.          (11)

This uses s2(m)=s2(u)+s2(t), with no assumption on the other coefficients
except their integrality. They may be zero or have arbitrary positive
valuation. Unique minimal valuation forbids cancellation in any finite
sum, independently of completeness of K.

## 3. Circular-prefix transfer of the same window

**Theorem 1 (off-critical unit-window transfer).** Under (1) and (7),
if the R_r have a positive unit window of length w, then the M_m have
a nonzero window of the SAME length w at every positive start.

Indeed consider any w consecutive positive integers m and their residues
modulo L. If w>=L, one of these residues is0, for which R_0=1. If w<L
and the residue block wraps or begins at0, it also contains0. Otherwise
it is an ordinary block of w consecutive POSITIVE indices in the finite
prefix R_1,...,R_(L-1), so the assumed unit-window property supplies a
unit R_t. In each case (10)--(11) give a nonzero M_m in the desired block.
No periodicity of the infinite coefficient sequence R_r is asserted.

Even without any unit-window hypothesis, every multiple of L has
R_(m mod L)=R_0=1, so every L consecutive positive moments contain a
nonzero one. This weaker L-bound applies also when the residue series
is polynomial; it need not be a support-size SFC bound.

## 4. The remaining elementary ranges and the boundary warning

If delta<0, every positive M_m is nonzero. For r>=1, subtracting the
valuation m delta of the r=0 term in (5) gives a lower bound

       r eta+s2(m-r)-s2(m)
          >=r eta-s2(r)>=r(eta-1)>0,

where s2(m)<=s2(m-r)+s2(r). Thus the constant term is uniquely minimal.

If delta=0, every positive EVEN M_m is a unit: every (m)_r with r>=1
contains the even integer m, while k^m is a unit. Hence a universal
two-moment window holds. One may NOT infer transfer of w=1 here: for
q=X and k=-1 all R_r=1 are units but M_1=0.

If delta>=1 and R_m is a unit, formula (6) again makes r=m uniquely
minimal. For j>0 the extra cost s2(j)+(delta-1)j is strictly positive.
Therefore a positive unit window of length w for R transfers to the
same moment window for M. Equivalently, the formal series exp(kT)
has residue1, since v(k^j/j!)=j(delta-1)+s2(j)>0 for j>0. This is
a formal coefficient statement, not evaluation of a series at a
negative-valued argument.

The open bands (7) cover 0<delta<1 except the critical numbers
theta_a=1-2^-a, a>=1. At theta_a with L=2^a, write m=2L u+epsilon L+t,
epsilon in{0,1}, 0<=t<L. The same digit-cost proof shows that the
minimizing r is uniquely t if epsilon=0, and consists of t,t+L if
epsilon=1. Thus the proof of Theorem1 does NOT extend across the tied
case merely by assuming the residue series is nonpolynomial.

At delta=1/2, multiples of4 have epsilon=0,t=0, so every four consecutive
positive moments contain a nonzero one, for every q satisfying (1).
This closes an actual-support SFC window when the phase has at least
four monomials. It is not a two- or three-support endpoint theorem.
More generally every 2L consecutive positive moments contain a nonzero
one at theta_a; this is not a general support-size conclusion.

## 5. The exact binary-quadratic residue and its unit windows

For (4), let bars denote residues. The native integral generating series
R(T)=sum R_r T^r has residue

       Rbar(T)=N(T)/D(T),
       N=1+(dbar+ebar+betabar)T+dbar ebar T^2,
       D=(1+alphabar T+dbar^2 T^2)
         (1+gammabar T+ebar^2 T^2).                         (12)

Here is a coefficient-level derivation. A monomial with i quadratic-X,
h quadratic-Y, j cross, l linear-X and b linear-Y factors has coefficient

       [(2i+j+l)!/(2^i i! j! l!)]
       [(2h+j+b)!/(2^h h! j! b!)] j!.                       (13)

Both brackets are integers by labelled matchings and coloured singleton
counts. Thus j>=2 contributes0 modulo2. For j=0 the two edges factor.
Writing u_n=Lambda(((alpha/2)X^2+dX)^n)/n!, their exact native recursion is

       u_0=1, u_1=alpha+d,
       u_n=alpha(2n-1)u_(n-1)+d^2 u_(n-2), n>=2.            (14)

It follows directly from the coefficient
(n+i)!/[2^i i!(n-i)!]. After a common factorial is removed, the recursion
is the identity
2i(2n-1)+(n-i)(n-i-1)=(n+i)(n+i-1), including its zero endpoint terms.
Modulo2 its series is U=(1+dbar T)/(1+alphabar T+dbar^2 T^2).
The shifted edge A=sum Lambda(X((alpha/2)X^2+dX)^n)T^n/n! has coefficient
(n+i+1) times the preceding native coefficient. Hence over F

       A=U+T partial_T U+alphabar partial_alpha U
         =1/(1+alphabar T+dbar^2 T^2).

This derivative identity is first an identity over universal residue
indeterminates and then specializes. The j=1 contribution is beta T
times the product of the two shifted edges, giving (12). All series
identities are formal.

Suppose Rbar is NONPOLYNOMIAL. Reduce (12) to N0/D0 with constant
terms1, put w=deg D0>=1, and note deg N0<=w. To check this last
inequality before cancellation: if D is nonconstant with neither
linear residue nonzero, deg N<=1<=deg D; with exactly one linear
residue nonzero, deg N<=1 and deg D>=2; with both nonzero, deg N=2
and deg D=4. Cancellation preserves the degree difference.

The positive coefficients of N0/D0 cannot have w consecutive zeros.
If a zero block starts at n>=2, its last index is at least w+1, where
the recurrence is homogeneous. The nonzero leading coefficient of D0
then forces the preceding coefficient to vanish. Repeating moves the
block back to indices1,...,w. Those zeros would make N0=D0 by comparison
through degree w, contradicting nonpolynomiality. A block starting at1
gives that contradiction directly. Thus R_r has a positive unit window
of length w. Properness is not assumed when deg N0=w; the preceding
argument deals explicitly with that case.

## 6. Actual support size in the binary application

Let t=N(q) count the nonzero monomials in (4) over K, and let s=t+1
be the actual support of q+k, where k!=0. In the nonpolynomial case
the denominator degree w in Section5 satisfies

                             w<=s.                         (15)

One can check this without any Hessian or rank assumption. Let t0<=t
be the number of nonzero natural coefficient RESIDUES in (4). If D has
degree1 it uses a diagonal residue; degree2 uses either a linear residue
or two diagonal residues; degree3 uses a linear residue and the opposite
diagonal residue. These give degD<=t0+1. Degree4 requires both linear
residues. If t0>=3 it again gives4<=t0+1. The only remaining degree4
case is t0=2, the two linear residues alone. There beta=alpha=gamma=0
in F, and N=(1+dbar T)(1+ebar T) cancels from
D=(1+dbar T)^2(1+ebar T)^2, leaving denominator degree at most2.
This also covers dbar=ebar. Since reduced degree cannot increase,
(15) follows. Nonzero original coefficients of positive valuation may
increase actual support and do not invalidate the inequality.

**Corollary 2 (binary off-critical actual-support sectors).** For a
binary quadratic q+k in natural normalization (4), with k!=0, if
Rbar is nonpolynomial and delta is in any open band (7), or delta>=1,
then at every positive start at least one of s=N(q+k) consecutive
factorial moments is nonzero. The stronger w-window of Section5 holds.
For delta<0 every positive moment is nonzero; for delta=0 or
0<delta<1/2 the universal two-window suffices for actual supports s>=2,
even if Rbar is polynomial. At delta=1/2 the universal four-window
suffices for s>=4. These last endpoint/low-band statements do not
assert the same w for every residue series.

For an arbitrary given binary phase, natural normalization means divide
the WHOLE phase by a coefficient of minimum valuation among
2A,B,2C,d,e of its nonconstant part. This makes the five natural
coefficients integral and at least one a unit. The relative constant
valuation delta is computed after that same division. Whole-phase
scaling preserves all moment zero sets and actual monomial support.
The constant-only case is separate and has no zero positive moments.

No statement says that every complex phase has a place/embedding meeting
the successful sector hypotheses. The corollary applies at a specified
2-adic embedding and normalization. It therefore does not close full
binary-quadratic SFC or SFC(2). In particular, polynomial residue-series
exceptions and the higher critical dyadic radii are not disposed of.
For dense actual six-support phases, the first critical radius1/2 is
already covered by the four-window; theta_a with a>=2 still needs a
separate tied-anchor argument.

## 7. Exact limits of the bridge

The proof uses only finite coefficient identities, coloured partitions,
integer factorial valuations, formal rational series and unique-minimum
valuation. It requires neither Gaussian analytic continuation nor a
holonomic-rank-to-window inference. It transfers a proven residue UNIT
window, not merely a characteristic-zero nonzero window for q: a
nonzero R_r of positive valuation is not sufficient for (11).

At a critical radius the two minimal terms may cancel. Multiplying a
residue generating series by a Frobenius-power linear factor may help
describe those ties, but nonpolynomiality alone is not a proof that
every desired cross-boundary block has a surviving term. No endpoint
closure, beyond the explicit elementary windows in Section4, is claimed.
