# Numerator and denominator records of ordinary abc odds sums

Keep positive primitive rows \(A+B=C\), \(X/2<C\le X\)
with original total cap \(\operatorname{rad}(ABC)\le X^L=n\)
and individual caps \(n^{1/3+\tau}\). Let \(N\ge n^{1/3-\tau}\)
be the cohort size. In each fixed slot use an original ratio
\(U_i/W_i\) with distinct roles. Allow fixed nonzero rational
weights \(e_i\) whose numerators and denominators have height
\(\le X^K\) for fixed \(K\). Fixed here means prescribed
uniformly over tuples, and allows dependence on \(X,L\).
Whenever required a specified shift may have the same height
bound. All small fixed margins and the tuple length \(j\ge2\)
are chosen before sufficiently small fixed \(\tau,L\);
the height limit is taken at each such \(L\).

Remove a power-small class of tuples whose shared input
radical outside a prescribed polynomial-height pool has
normalized mass above a fixed small \(\sigma>0\). This is EB5
per anchor row, followed by a union over pairs of slots.
Record the remaining shared support at price \(\sigma\).
Let \(S\) include that support, the prescribed pool, and all
prime divisors of numerators and denominators of the weights
and shifts. These added prescribed supports have polynomial
height. Write \(W_i=W_i^S d_i\), where \(W_i^S\) is the full
supported part, and let each \(d_i\) be the private complement.
The \(d_i\) are pairwise coprime and are units on all other
input rows.

## One exact shifted numerator and the larger half of the denominators

**ON1.** Write \(Y=A_*/B_*=\sum_i e_i U_i/W_i+s\) in lowest terms.
At every \(p\mid d_i\) only one summand has a pole, to exactly
the order of \(p\) in \(d_i\). Hence
\[
                    B_*=B_0\prod_i d_i
\]
with \(B_0\) supported on \(S\). Its size is polynomial in \(X\)
and it is an exact supported-power record at subpower cost.
Select the \(h=\lceil j/2\rceil\) largest \(d_i\)'s, resolving
ties by index, and expose those original rows. The set of
their slot positions is a finite record. Let \(Q\) be their
product, and \(V\) the product of the remaining \(d_i\).
The selecting rule gives \(1\le V\le Q\) for every candidate
with the same exposed-position record and exposed rows.
No numerical margin, cutoff, or dyadic height estimate is
needed.

For an exposed position \(i\), multiplying by \(B_*\) and
reducing \(p\)-adically at the full powers in \(d_i\) gives
\[
 A_*\equiv \left(
       e_i U_i\, \frac{B_0}{W_i^S}
                 \prod_{\substack{\ell\ \mathrm{exposed}\\ell\ne i}}
                       d_\ell\right)\,V\pmod{d_i}.
\]
The displayed coefficient is known and a unit modulo \(d_i\).
All other summands vanish locally modulo those powers.
The equation need not have been cleared integrally at the
primes of \(S\), which are not used in these congruences.
The exact numerator and CRT therefore recover \(V\bmod Q\),
and \(1\le V\le Q\) gives only one possible \(V\), even at
equality. For \(Q=1\), the inequality already gives \(V=1\).
Thus \(B_*\), and then the exact rational output, are recovered.

In a naturally mixed orientation, the private numerator of
each \(U_i/W_i\) is at most \(d_i\). After fixing \(B_*\),
enumerating \(d_i\mid B_*\), and fixing all exact \(S\)-powers
at subpower cost, the usual private-denominator CRT determines
each private numerator. This is also true in a pure
orientation with at most one unfavorable position, since the
sum and the other fractions fix the last. The same recovery
applies to the NR1 functions with all-role denominators
without an orientation, since those denominators recover
every original valuation.

Discard zero full sums at \(O(N^{j-1})\) cost by fixing
any \(j-1\) rows in a fixed orientation. For the ordinary
numerator count its old support outside \(S\) is small
off a power-small class. At a private
prime on the denominator role, the numerator is a unit.
At a private prime on one of the other two roles, the
original odds reduces to \(\zeta\in\{0,1,-1\}\), depending
on the role and signs. Fix the other slots. At such a
prime the numerator can vanish only if it divides the
numerator of the anchor
\(s+\sum_{\ell\ne i}e_\ell U_\ell/W_\ell+e_i\zeta\).
A nonzero anchor has polynomial height, so EB5 on the
varying row discards the power-small class where that
intersection exceeds any fixed small normalized mass.
A zero anchor fixes one of the other rows after the rest
are fixed, since its exact odds is injective on primitive
rows and its weight is nonzero. This costs \(O(N^{j-1})\)
full tuples. The argument is applied over all finitely
many orientation patterns before restricting to any chosen
pattern. No extra old-support step is needed for NR1,
whose private original primes are all poles.

Consequently for every \(\delta>0\), naturally oriented
ordinary odds have fresh numerator radical greater than
\(n^{\lfloor j/2\rfloor/3-\delta}\) outside a power-small
class. The count for a proposed fresh cap \(n^r\) is
\[
                  n^{h/3+r+O_j(\sigma)+\eta}X^\epsilon ,
\]
where \(\eta\) absorbs freely small original count allowances.
The exact numerator record costs \(r\) plus the old-support
allowance, the exposed rows cost \(h/3+\eta\), and the CRT
completion needs only the already paid shared record.
This also sharpens the even-\(j\) threshold in NR2 on NR1
assignments, where sorting exact \(d_i\)'s replaces the
strict numerical margin of the previous proof.

## Complementary fixed affine shifts and the two cheapest caps

**ON2.** Use naturally oriented odds as in ON1, or a pure
orientation with at most one unfavorable slot. Let
\(Z=A/B=\sum_i e_iU_i/W_i\) be reduced. Choose two
distinct prescribed rational shifts \(s_1,s_2\) of
polynomial height and clear their denominators with a
prescribed positive \(k\), putting
\(m_i=k s_i\in\mathbb Z\). Use the raw forms
\[
             L_i=kA+m_i B,\qquad
             L_1-L_2=(m_1-m_2)B .
\]
The reduced numerator of \(Z+s_i\) is \(L_i/\gcd(L_i,kB)\).
This gcd is supported on \(k\) and divides \(k^2\), since
\(\gcd(L_i,B)=\gcd(k,B)\). The gcd of the three auxiliary
forms divides \(k|m_1-m_2|\), using
\(\gcd(L_1,B)\mid k\). In particular all reductions here
involve only prescribed coefficient primes. Enlarge the
pool by the primes of \(k(m_1-m_2)\). Discard \(L_i=0\)
at \(O(N^{j-1})\) tuple cost for each orientation and shift.
After dividing by the common gcd, taking signs, and
reordering, the three forms are a primitive positive
auxiliary additive triple. The exact value of \(Z\)
can be recovered from this auxiliary triple, up to the
finite role/sign record, by the fractional linear ratio
\(L_1/L_2\).

For the height floor, fix all other rows of a tuple and
let \(P_i\) be their full prime support together with the
prescribed pool and coefficient primes. It contains \(S\).
If any particular original role on the varying row has
\(S\)-private part \(<X^{h_0}\), it has \(P_i\)-private part
\(<X^{h_0}\). Fixing its exact \(P_i\)-part costs subpower,
and CF1 for one integer feature excludes this class with
a uniform strict saving if \(h_0>0\) is sufficiently small.
The other rows are unrestricted. Thus, outside a power-small
class, every possible \(d_i\) exceeds \(X^{h_0}\).
The primitive auxiliary triple has positive fixed numerical
height. Its normalized \(B\)-role contains every exact
\(d_i\), since both \(m_1-m_2\) and the common gcd are
units at those primes. Thus even for growing coefficient
heights one has the required numerical floor.

Fix all complete supported powers of the auxiliary triple
on the prescribed coefficient pool at subpower cost.
The remaining radical cap of the \(B\)-role is
\(n^{j(1/3+\tau)}\) since its non-coefficient support comes
only from one original role per row. By the old-support
step above, the other two residual caps are
\(n^{r_1+O_j(\sigma)},n^{r_2+O_j(\sigma)}\) if their fresh
numerator masses are at most \(r_1,r_2\).
Assume both \(r_1,r_2\) are bounded by a fixed \(C\), which
is chosen before the small-budget threshold. MCB2 applies
with the fixed supported factors as arbitrary coefficients
and the original auxiliary numerical shell. It gives
for the output exponent
\[
 \frac12\min\{r_1+r_2,\ r_1+j/3,\ r_2+j/3\}
                    +O_j(\sigma+\tau)+\eta .
\]
The exact odds sum then recovers the original tuple at
\(n^{\sigma+o_X(1)}\) cost by the private-denominator completion of ON1. In particular, for every
fixed \(C,\delta>0\), except on a power-small class, either
both shifted fresh numerator masses are at least
\(j/3-\delta\), or an offending small one is accompanied
by a complementary mass above \(C\). This is a
two-shift output condition. MCB2 cannot be invoked with the
other auxiliary numerator's raw cap unbounded.

## Two denominator supports without any numerator record

**ON3.** Take two prescribed odds sums on the same tuple, possibly
with different polynomial-height nonzero weights and shifts.
Require that their denominator roles \(W_i,W_i'\) are
different in every slot. Enlarge \(S\) by both coefficient
pools and the shared input primes. For the first reduced
denominator \(B\), its component outside \(S\) contains
exactly the private primes and exact powers of all \(W_i\).
The second \(B'\) similarly carries all \(W_i'\). Given
merely \(\operatorname{rad}(B),\operatorname{rad}(B')\),
there are subpower choices for their full powers, since
both heights are polynomial. Allocate private primes
to slots at cost \(j^{\omega(BB')}=X^{o_X(1)}\).
Fixing exact \(S\)-parts is again subpower after paying the
shared support. This reconstructs two original letters in
each slot, hence the entire tuple. It needs no private-size
orientation, no unrestricted product set, and no assumption
on complete-letter multiplicity.

With one denominator support, the same argument fixes one
complete role per slot. One may then combine it with the
numerator method in ON1 in an allowed orientation for a
different denominator role and expose only \(\lceil j/2\rceil\) genuine fixed-letter
fibres rather than unrestricted rows. In a regular cohort
the exposed exponent for these is the sum of the respective
original degrees \(u_T=1/3-\alpha_T\). This conditional
improvement disappears at occurrence if obtaining the
first denominator dictionary already costs its full
complete-role entropy.

## Two shifted numerators in different roles for pairs

**ON4.** For two rows let \(x_i=A_i/C_i\), \(\Delta=A_1C_2-A_2C_1\),
\(P=A_1A_2\), \(Q=C_1C_2\), and \(J=A_1C_2+A_2C_1\).
Choose prescribed nonzero rational shifts \(\theta,\phi\)
of polynomial height with \(\theta\phi\ne4\). The raw shifted numerators are
\[
             F=\Delta+\theta Q,\qquad G=-\Delta+\phi P .
\]
They correspond to \(x_1-x_2+\theta\) and
\(x_1^{-1}-x_2^{-1}+\phi\), respectively. Clear the fixed
shift denominators. The exact cancellation from each raw
numerator to its reduced numerator is supported only on
fixed coefficient primes and the shared prime support of
\(C_1,C_2\) or of \(A_1,A_2\), respectively. At a private
prime of such a denominator, the unshifted determinant
is a unit, and adding a fixed rational multiple of the
denominator cannot cancel it away from coefficient primes.
Thus the two reduced numerators, the shared support, and
subpower coefficient/valuation records give \(F,G\).

As \(J^2=\Delta^2+4 P Q\), with \(\kappa=\theta\phi\),
\[
       \kappa J^2=(\kappa-4)\Delta^2
                              +4(F-G)\Delta+4 FG .
\]
Completing the square gives
\[
 [(4-\kappa)\Delta-2(F-G)]^2+\kappa(4-\kappa)J^2
                     =4[F^2+G^2+(2-\kappa)FG].
\]
For \(0<\kappa<4\) the right side is nonzero since \(F,G\)
cannot both vanish. For general \(\kappa\ne0,4\), the zero
right side either fixes \(F=G=0\) or \(F/G=\lambda\) for
one of at most two constants. The first case fixes
\(x_1-x_2=-\theta\) and \(x_1 x_2=-\theta/\phi\), so has
boundedly many pairs. In the second case divide \(F=\lambda G\)
through by \(Q\) to obtain
\[
                 (1+\lambda)(x_1-x_2)+\theta
                            -\lambda\phi x_1x_2=0.
\]
After fixing \(x_1\), this has at most one \(x_2\) except
at at most one \(x_1\) where the coefficient vanishes.
It therefore gives \(O(N)\) actual pairs. The uniform
quadratic norm estimate on the remaining nonzero right
sides gives height-subpower possibilities for \(\Delta,J\),
including in the indefinite case. All denominators in
that application are polynomial in height and can be
cleared exactly. These give \(P,Q,A_1C_2,A_2C_1\),
then the primitive rows by divisor enumeration. Thus
the two exact shifted numerators, with no denominator
record, have generic fibres of exponent at most \(\sigma\).
Their joint fresh radical union is larger than
\(n^{2/3-\delta}\) outside a power-small set, since a small
joint support pays for both exact numerators at once,
using the old-support removal in ON1. The argument uses the
exact numerator records.

## The critical shift product and affine slices

For prescribed nonzero rational shifts with \(\theta\phi=4\),
put \(\lambda=\theta/2\). The two raw numerators satisfy
\[
     F-G=2(A_1+\lambda C_1)(C_2-\lambda^{-1}A_2).
\]
Except for the \(O(N)\) pairs at one of the two zero
factors, this gives the two exact affine values
\(P_1=A_1+\lambda C_1,\ P_2=C_2-\lambda^{-1}A_2\)
by divisor enumeration after clearing fixed denominators.
Now \(F=P_1 C_2+\lambda C_1 P_2\), so either original
row fixes the other.

The affine slice at fixed \(P_2\ne0\) has at most
\(n^{1/6+\tau/2+\eta}X^\varepsilon\) remaining original rows.
Indeed clear the prescribed direction denominators to write
\(b C_2-a A_2=d\), with \(a,b\) nonzero prescribed integers.
By original primitivity the gcd of the two terms divides
\(|ab|\) and can be chosen at divisor cost. At fixed shifts the normalized auxiliary has
height at least a fixed multiple of \(X\). For growing
coefficient heights, first discard the power-small class
where the part of \(C_2\) outside the prescribed coefficient
primes has height less than \(X^{h_0}\). As in ON2, fix the
exact supported part and use the one-feature bound CF1.
The remaining private part survives the normalization,
giving the fixed numerical height floor. Fix the exact
supported parts on the primes of \(abd\) at subpower cost. The
third auxiliary role is then entirely in its prescribed
coefficient, and the two variable roles retain their
original caps. MCB2 gives the asserted exponent by
charging the zero residual cap and the cheaper variable
cap. Thus at the critical shift product the joint fresh
numerator radical is generically at least \(n^{1/2-o(1)}\)
up to arbitrary small tolerances. When \(\theta=2\),
\(P_2=B_2\), and when \(\theta=-2\), \(P_1=-B_1\).
Use the actual complete \(B\)-degree \(n^{u_B+o(1)}\)
on a regular cohort. The resulting joint floor is
\(2/3-u_B\) in normalized units, reaching \(2/3\)
at occurrence.

The critical product permits large fibres in the absence
of the original radical caps. For \(M=2^k\) with \(k\ge3\),
set \(X=3M/2\) and choose odd \(C_1\) in \((3M/4,4M/5)\)
coprime to \(3\). Put
\[
 A_1=M-C_1,\quad B_1=2C_1-M,\quad
 C_2=9M/4-C_1,\quad A_2=C_2-M,\quad B_2=M.
\]
The two rows are positive primitive rows in the same shell.
At \(\theta=\phi=2\) their exact records are
\(F=9M^2/4\) and \(G=M^2/4\). These are already reduced
numerators since \(C_1C_2\) is coprime to \(6\) and
\(A_1A_2\) is odd. This gives \(\gg M\) distinct pairs,
with only a bounded shared prime support. For any two
nonconstant distinct-row linear forms defining the
letters, eliminating \(C_1\) leaves a fixed nonzero rational
multiple of \(M\). The fixed letter \(M\) has only prime
factor \(2\). Thus bounded overlap alone does not give
subpower completion at the critical product.

## Simultaneous shift dictionaries on all three actual locks

**ON5.** Suppose \(N\ge n^d\) genuine primitive rows in a
fixed sign convention are assigned a triple of actual
locks \((D_A,D_B,D_C)\), each dividing its indicated
role, with at most \(n^c\) exact joint records and
\(d>c\) by a fixed margin. In each actual record take
\(j\ge2\) rows and two fixed odds maps whose denominator
roles are different. Use nonzero integer weights of
sum zero, fixed before the height limit. The joint
dictionary of the two exact signed sums has exponent
at least \(jd-(j-1)c\), up to arbitrary fixed allowance.
If after any two exact locks are fixed the possible
third lock costs at most \(n^\xi\), then the union of
the one-map outputs over the two denominator choices
has exponent at least \(jd-(j-1)c-\xi\).

Indeed in one fixed denominator role the odds in
each term has a fixed residue \(0,1\) or \(-1\) at
each of the two complementary locks. This uses
primitivity and the actual signed addition. Thus the
two complementary locks divide the reduced numerator
of a nonzero weighted sum. They are coprime, since
they divide different primitive roles in one row.
The pair of sums with different denominator roles
recovers all three exact locks at divisor cost from
their numerators. For one orientation its two locks
are recovered at the same cost and the remaining
catalogue contributes \(n^\xi\).

Apply the rational sumset bound and the different-pole
dichotomy within each whole record fibre. The joint
two-role image has cardinality at least \(N_{\mathbf D}^j\)
up to the arbitrary fixed \(X^\varepsilon\) loss. One
orientation per fibre suffices in the one-map statement,
using the choice guaranteed by the different-pole
statement. Each zero-output
exception contributes at most \(N_{\mathbf D}^{j-1}\)
tuples because any one original odds determines the
remaining row. In large fibres delete these tuples
from the large selected set of OS1, so both outputs
can be required nonzero. Finally sum over actual lock
records and use the same convexity argument as OS4.
All exponents concern the original shell and the
actual catalogue. The value of the third lock is
charged only in the one-map conclusion.
[Exact divisor-fibre shifts (OS4)](ABC_SIGNED_RATIONAL_ODDS_SUMS.md#exact-shifts-across-actual-divisor-fibres)

For two distinct original rows and an unweighted signed difference,
this has an exact extra identity. Write \(a_W\) for the
absolute reduced numerator in denominator role \(W\), and
\(\Delta\) for the primitive cross determinant. Then
\(a_W=|\Delta|/\gcd(\Delta,W_1W_2)\). The removed gcd is
supported only at primes shared by the two \(W\)-values.
These supports are disjoint for distinct original roles,
by primitivity. Thus for distinct \(W,T\) one has
\(\operatorname{lcm}(a_W,a_T)=|\Delta|\), including at
prime powers. The two reduced numerators together retain
the entire original determinant.

## Proof inputs

[Private factors and exact odds recovery (OS1–OS4)](ABC_SIGNED_RATIONAL_ODDS_SUMS.md)
· [Numerator floors for all-role denominators (NR1–NR3)](ABC_PROJECTIVE_SUM_NUMERATOR_RADICALS.md#statements-and-the-height-first-setting)
· [Effective prime pricing (EB5)](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)
· [Complete-letter affine features (CF1)](ABC_COMPLETE_LETTER_SMALL_FEATURES.md)
· [Arbitrary fixed-coefficient counts (MCB2)](ABC_MULTIPLE_SHORT_COFACTOR_BOUNDARY.md#bounded-group-dictionaries-with-arbitrary-fixed-coefficients)
· [Quadratic norm count, Browning–Verzobio, Lemma 2.3](https://arxiv.org/html/2608.24512v1#S2.SS1)
