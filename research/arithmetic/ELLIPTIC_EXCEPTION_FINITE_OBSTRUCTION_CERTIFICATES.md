# Finite obstruction certificates on the elliptic slice

Every quadratic class passes the combined baseline `2,3,5`. For a
larger finite set, failure can be certified by a finite cover of the
actual joint congruence image by locally insoluble cells. At a large
largest prime, such a failure has an especially restrictive form:
every previously surviving reduction must be one of the two fatal
zero-branch residues. No proof that this configuration never occurs
in the two unresolved classes is supplied.

## 1. A feasible combined baseline in all four classes

Retain `Gamma=E_0(Q)`, `A={1,2,5q_0,10q_0}`, and the finite-place
counts and coefficients of
[EFA2](ELLIPTIC_EXCEPTION_FINITE_LOCAL_ASYMPTOTICS.md). Use

\[
 E_0:y^2=x^3-10656x^2+27091730x,\qquad q_0=2709173,
\]
\[
 P=\left(\frac{373248}{49},\frac{59459616}{343}\right),\quad
 Q=\left(\frac{3442688}{529},\frac{304670016}{12167}\right),
 \qquad T=(0,0).
\]

The exact coefficients in
[J235173](ELLIPTIC_EXCEPTION_JOINT_2_3_5_17_SIEVE.md) already prove
positivity at the stronger set `2,3,5,17` in all four classes. In
particular, for `S_0={2,3,5}` and every `a in A`,

\[
 \mu_{S_0,a}>0.
                                                               \tag{1.1}
\]

The next lemma records concrete witnesses rather than a new
positivity classification.

**Lemma EFO1 (concrete combined-baseline witnesses).** In the two
classes not yet known feasible for every finite set, admissible full
field-fibre witnesses passing `2,3,5` are `Q` for class 2 and
`U=2P+3Q+T` for class `5q_0`, where

\[
 x(U)=\frac{48325966141496822395419633600628716597292868897208265}
            {124529622664493793813505034507747191688681459286016},
\]
\[
 y(U)=-\frac{131593797705730303510630822117943515388925347362658832076971985222495261005832965}
             {1389661443981606892683015531872483448667279251997885870671390295580752216064}.
                                                               \tag{1.2}
\]

**Proof.** For `Q`, its two-adic
component is translated by `T`, its three-adic reduction is `X=2`
with `v_3(y)=2`, and its five-adic `x` coordinate is integral.
ETF2--ETF3 and ETC1 give all three local points.

Write the reduced coordinates of `U` as `x=a/c^2,y=b/c^3`, `c>0`.
Exact arithmetic gives

\[
 v_2(c)=9,\quad v_5(c)=v_3(c)=0,\quad
 a=1\pmod3,\quad v_3(b)=4.                         \tag{1.3}
\]

Thus its formal two-adic depth is nine, its three-adic point is in
the `X=1` disc at depth four, and it is integral at five. The same
local theorems apply. The exact checker verifies class `5q_0`,
admissibility, and irreducibility of the cubic: its reduction is
irreducible at seven. For `Q` the corresponding irreducibility prime
is 29. The quadratic discriminants are negative and nonsquare.
These are full field fibres of the existing map.

Each witness lies in the open class-specific soluble set. Hence its
Haar measure is positive by EFA1; EQA2 and ESQ3 then give a positive
rational leading coefficient and the sharp finite-place count.
\(\square\)

The baseline certificate is finite arithmetic. It does not assert that
either witness passes all the other primes.

## 2. Complete finite-cell obstruction and minimality criteria

For a finite set of places `S`, write `S_f` for its finite primes;
the real place is automatic and cannot occur in a minimal obstruction.
Use the supplied-input setting of
[FQS2--FQS3](ELLIPTIC_EXCEPTION_FINITE_QUOTIENT_SATURATION.md).
Thus a certified open product `W subset H_S`, its finite ambient
quotient, and its exact rational image

\[
 C=H_S/W
\]

are given. This includes the saturation obligations proving that `C`
is the image of the entire rational group. For each cell `f in C`
and prime `p in S_f`, let `r_{f,p} in[0,1] cap Q` be its conditional
soluble local Haar weight, computed by EQA2. Write `C_a` for its
cells with global class `a` and put

\[
 B(f)=\{p\in S_f:r_{f,p}=0\}.                       \tag{2.1}
\]

All four class sets `C_a` are nonempty, since the quotient retains
the global four-class map.

**Theorem EFO2 (finite cover certificates and minimal obstructions).**
Under these certified-quotient hypotheses, the following are equivalent:

1. no admissible full field target of class `a` passes `S`;
2. `mu_{S,a}=0`;
3. `B(f)` is nonempty for every `f in C_a`.

The third condition is a finite congruence cover: assign to each
class-`a` cell at least one prime whose local conditional weight is
zero. On that cell every admissible point is locally insoluble at
the assigned prime.

Furthermore `S_f` is a minimal obstruction, under deletion of places,
if and only if condition 3 holds and

\[
 \text{for every }p\in S_f\text{ there exists }f_p\in C_a
                          \text{ with }B(f_p)=\{p\}.          \tag{2.2}
\]

These conditions are decidable finite rational arithmetic once the
certified quotient and local weights have been supplied.

**Proof.** FQS3 gives

\[
 \mu_{S,a}=\frac1{|C|}\sum_{f\in C_a}
                           \prod_{p\in S_f}r_{f,p}.            \tag{2.3}
\]

All terms are nonnegative. Hence the sum vanishes precisely when
every cell has a zero factor. A zero local weight really certifies
insolubility throughout that local ball: if its intersection with
the admissible soluble open set were nonempty, that intersection
would be open and have positive conditional Haar measure. This
proves the stated covering interpretation. The equivalence of
conditions 1 and 2 is the finite-local dichotomy EQA2a.

The same quotient computes the count after removing any one prime
`p`, by replacing `r_{f,p}` by 1 in (2.3). Under the zero-cover
condition, that count is positive exactly when some cell has no
other zero factor, which is (2.2). Feasibility is downward closed,
so positivity of every one-prime deletion is equivalent to positivity
of every proper subset. This proves minimality. The empty finite
set has positive coefficient and is never an obstruction.
\(\square\)

The complete input is the finite group, its exact full rational image,
its class labels, and certified rational local weights. This criterion
does not compute those inputs for arbitrary `S` or replace their
partial-saturation certificates by a sampled subgroup. It also does
not claim that an actual obstruction has been found.

By EOP1 no singleton is an obstruction. The existing J235173
positivity shows that no subset of `{2,3,5,17}` is an obstruction.
By ECO1 and the certified-class result,
there are no finite obstructions in classes 1 and `10q_0` at all.
Thus any remaining minimal obstruction has at least two active
primes, includes a prime outside `2,3,5,17`, and belongs to class 2 or
`5q_0`.

## 3. A large last prime must be a fatal zero branch

Use the cutoff `p_*`, formal image depth `e_p`, branch offset depth
`k_p`, and residue hazard table of ETP1 and ETP4 in
[the triangular prime sieve](ELLIPTIC_EXCEPTION_TRIANGULAR_PRIME_SIEVE.md).
For a good active prime `p>p_*`, define

\[
 Z_p^{\rm fatal}=
 \begin{cases}
 \{\overline{R_+},\overline{R_-}\},&
 (51730/p)=1,\ k_p<e_p,\ 3\nmid k_p,\\
 \varnothing,&\text{otherwise}.
 \end{cases}                                                   \tag{3.1}
\]

Fix a feasible finite `S^-` consisting of primes smaller than `p`
and a class `a`. Let

\[
 \beta_r=\nu_\infty\bigl(
    \{\alpha=a\}\cap\bigcap_{q\in S^-}\mathscr G_q
                     \cap\{\overline P_p=r\}\bigr).           \tag{3.2}
\]

This is a rational finite-place weight with an extra finite reduction
condition. Its finite-quotient computation is (2.3) with the new
prime's local solubility weight omitted and its reduction label fixed.

**Theorem EFO3 (exact terminal obstruction criterion).**

\[
 \boxed{\mu_{S^-\cup\{p\},a}=0
       \quad\Longleftrightarrow\quad
       \{r:\beta_r>0\}\subseteq Z_p^{\rm fatal}.}    \tag{3.3}
\]

In particular an active largest prime `p>p_*` of a minimal obstruction
in class `a` must satisfy all of

\[
 p\in\mathcal P_{a,d},\qquad e_p\ge2,\qquad
                     k_p<e_p,\quad3\nmid k_p,\qquad d=51730.
                                                               \tag{3.4}
\]

Here `P_{a,d}` is the numerator-eligible prime set of EPE2.
If `e_p=1`, every feasible smaller-prime prefix in every class
extends positively at `p`. More generally this holds whenever
`Z_p^fatal` is empty, or whenever one nonfatal reduction has
positive prefix weight.

**Proof.** The prefix and global class are measurable before the
fresh `p`-primary component is exposed, by ETP1. ETP4 therefore gives

\[
 \mu_{S^-\cup\{p\},a}
       =\sum_r\beta_r\,(1-w_p(r)).                 \tag{3.5}
\]

The formal soluble mass at the pole is positive for every `e_p`.
It is also positive in a finite zero residue if `k_p>=e_p`.
If `k_p<e_p`, the whole residue coset is bad precisely when
`3` does not divide `k_p`; otherwise it is good. All nonbranch
residues are good. Thus the exact zero set of the nonnegative
survival factors in (3.5) is (3.1), proving (3.3).

A nonzero feasible prefix has nonempty support. A fatal branch
therefore requires `e_p>k_p>=1`, proving `e_p>=2`. A supported finite
zero branch has a nonzero positive valuation of `F`, so EPE1--EPE2
force both `a` and `d` to be local squares. This is numerator
eligibility. A minimal obstruction has feasible deletion of its
largest prime, so the preceding argument applies and proves (3.4).
\(\square\)

The positive-density Kummer primes of
[EKX3](ELLIPTIC_EXCEPTION_KUMMER_BRANCH_EXCLUSION.md) have no supported
finite zero reductions at all. At every sufficiently large such
prime, the right side of (3.3) is impossible. Thus the largest prime
of a remaining large-prime minimal obstruction must also lie outside
that Kummer exclusion set.

EOP2 gives a complementary finite restriction: after choosing a
certified soluble profinite coset for `S^-,a`, every killing new prime
lies in its explicit finite exceptional set. That set changes with
the prefix and coset. These results do not prove that fatal covers
cannot occur, that there are only finitely many simultaneous
elliptic-Wieferich primes, or that all sufficiently large primes
have `e_p=1`.

## 4. Exact evidence and remaining question

The [baseline certificate](elliptic_exception_finite_obstruction_baseline.json)
and [standard-library checker](verify_elliptic_exception_finite_obstruction_baseline.py)
verify (1.2)--(1.3), the class, both curve equations, admissibility,
the exact irreducibility primes and the small local tests. The
checker also exposes an exact arithmetic routine for classifying a
supplied finite local-weight table as in EFO2; certification of that
table's rational group image and local weights remains an input.
The [receipt](elliptic_exception_finite_obstruction_baseline_checks.json)
records the actual finite arithmetic scope. No obstruction table is
asserted to exist or certified by the baseline checks.

The unresolved issue for classes 2 and `5q_0` is whether an actual
finite joint image and its prior local conditions can force the
support in (3.3) entirely onto fatal zero branches. The available
independence controls the new formal coordinate, not those earlier
reduction choices. No covering-avoidance theorem or all-finite-set
positivity theorem for these two classes is claimed.
