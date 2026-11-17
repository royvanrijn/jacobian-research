# Which primes can obstruct a fixed elliptic quadratic class?

The global quadratic class makes many of the local tests on the
existing elliptic slice automatic. A denominator obstruction requires
that class to be locally square. A numerator obstruction requires one
additional quadratic splitting condition. These restrictions refine
the all-primes sieve without resolving its limiting measure or tail.

## 1. A local restriction on the varying function

Retain the curve and rational function of
[ELV1--ELV3](ELLIPTIC_EXCEPTION_EXPLICIT_LOCAL_CRITERION.md):

\[
 E_0:y^2=xQ(x),\quad Q(x)=x^2-10656x+27091730,
 \qquad F(x)=128x-25865.
\]

Put

\[
 q_0=2709173,\quad
 A=\{1,2,5q_0,10q_0\},\quad
 x_0=\frac{25865}{128},\quad d=51730=2\cdot5\cdot7\cdot739.
                                                               \tag{1.1}
\]

The [global descent classification](ELLIPTIC_EXCEPTION_ARITHMETIC.md),
Theorem 2.1, says that every admissible rational point has
`x in a Q*2` for exactly one `a in A`. Its quadratic field is
`K_a=Q(sqrt(-3a))`. The exact branch identity is

\[
 Q(x_0)=\left(\frac{639735}{128}\right)^2,
 \qquad 639735=3\cdot5\cdot42649,\qquad
 x_0=d\left(\frac1{16}\right)^2.                    \tag{1.2}
\]

**Proposition EPE1 (local squareclass restriction).** Let
`p notin {2,3,5,42649}` and let `(x,y) in E_0(Q_p)` have
`x F(x) != 0`. If `v_p(F(x)) != 0`, then `x` is a square in
`Q_p`. If, in addition, `v_p(F(x))>0` and `p` does not divide `d`,
then `d` is a square in `Q_p`.

**Proof.** If `v_p(F(x))<0`, the constants are integral and `128`
is a unit, so `v_p(x)<0`. Rewrite the curve equation as

\[
 (y/x)^2=x\left(1-\frac{10656}{x}
                         +\frac{27091730}{x^2}\right).
\]

The parenthesis belongs to `1+p Z_p` and is a square by Hensel's
lemma, because `p` is odd. Thus `x` is a square.

If `v_p(F(x))>0`, then `x` is integral and `x=x_0 mod p`.
By (1.2), `Q(x)` is congruent to a nonzero square modulo `p`.
It is a local square, so `y^2=xQ(x)` again makes `x` a square.
When `p` does not divide `d`, the number `x_0` is a unit and
`x/x_0 in 1+p Z_p`. It follows that `x_0`, and hence `d`, is a
local square. This last step is not asserted at `7,739`.
\(\square\)

This uses the same strict Hensel input as ELV1: Milne,
*Algebraic Number Theory*, version 3.00, Proposition 7.31 and
Theorem 7.32, printed pp.113--114,
[primary text](https://www.jmilne.org/math/CourseNotes/ANT300.pdf),
[retained copy](sources/milne-ant300.pdf).
No good reduction hypothesis on `E_0` is required.

## 2. An exact smaller sieve in each global class

For `a in A`, define the prime sets

\[
 \mathcal P_a=\{p\ge11:p\equiv2\pmod3,
                 \ a\in\mathbb Q_p^{*2}\},\qquad
 \mathcal P_{a,d}=\{p\in\mathcal P_a:
                 \ d\in\mathbb Q_p^{*2}\}.          \tag{2.1}
\]

Equivalently, `P_a` requires `p` not to divide the squarefree
representative `a` and `(a/p)=1`; the second set also requires
`(d/p)=1`. No prime in (2.1) divides `d`: its odd prime factors
`5,7,739` are either excluded or are one modulo three.

**Corollary EPE2 (class-specific Hasse criterion and redundant primes).**
Fix an admissible full quadratic-times-cubic field fibre with class
`a in A`, and write its reduced elliptic coordinate as

\[
 x=M/C^2,\qquad C>0,\qquad \gcd(M,C)=1,\qquad
 N=128M-25865C^2.
\]

Its fibre is Hasse-failing if and only if it passes the three local
tests at `2,3,5` and

\[
 \boxed{v_p(C)\in3\mathbb Z\quad(p\in\mathcal P_a),
 \qquad v_p(N)\in3\mathbb Z\quad(p\in\mathcal P_{a,d}).} \tag{2.2}
\]

More precisely, the following vanishing is automatic:

\[
 \begin{aligned}
 p\equiv2\pmod3,\ p\ge11,\ p\notin\mathcal P_a
       &\ \Longrightarrow\ v_p(F(x))=v_p(C)=v_p(N)=0,\\
 p\in\mathcal P_a\setminus\mathcal P_{a,d}
       &\ \Longrightarrow\ v_p(N)=0.
 \end{aligned}                                                   \tag{2.3}
\]

In particular the ramified prime `q_0` is automatic in the two classes
containing `q_0`; it is not an omitted Legendre-symbol case. It is
also automatic for `a=2`, since `q_0=5 mod8`.

For any fixed finite set of places `S`, let `S_a` retain the places
of `S` among `2,3,5` and among `P_a`. Then the finite-place counts
and their Haar coefficients of
[EFA2](ELLIPTIC_EXCEPTION_FINITE_LOCAL_ASYMPTOTICS.md) satisfy

\[
 N_{S,a}(B)=N_{S_a,a}(B),\qquad
                 \mu_{S,a}=\mu_{S_a,a}.              \tag{2.4}
\]

**Proof.** A prime `p=2 mod3`, `p>=11`, is outside all the
exceptional primes of EPE1. If `a` is not a local square, a nonzero
valuation of `F(x)` would contradict `x in a Q_p*2`. Hence that
valuation is zero. If `a` is a square but `d` is not, a positive
valuation is impossible by EPE1. These arguments include a class
ramified at `p`, where `a` cannot be a square.

For these odd primes, `gcd(N,C)` is a unit. Thus `p|C` gives
`v_p(F)=-2v_p(C)`, while `p|N` gives `v_p(F)=v_p(N)>0`.
This proves (2.3). ELV1 and ELV3 now give exactly (2.2).
The same results say that the erased primes in (2.4) are locally
soluble for every full field target of this class. Primes one modulo
three and the real place are automatic as well. This gives equality
of the counts for every `B`. Comparing their EFA2 leading terms
gives the equality of measures. It does not identify their ambient
joint groups with a product of independent groups. \(\square\)

By ELV2, this also restricts the maximal-order conductor: outside
`{2,3,5,42649}`, every prime dividing `f` in `Disc(L)=d_K f^2`
must have `a` locally square. This uses the field conductor, rather
than the squarepart of a defining-polynomial discriminant.

For example, `11` can obstruct only the classes `1,5q_0`, and only
through the denominator. The prime `17` can obstruct only the
classes `1,2`, through either numerator or denominator. In particular
`17` is automatic in the certified class `10q_0`, consistent with
the fixed-class coefficient `1/8` in
[J2174](ELLIPTIC_EXCEPTION_JOINT_2_3_17_SIEVE.md).

## 3. Densities of eligible primes

**Proposition EPE3 (prime eligibility densities).** The natural
densities among all rational primes are

| Global class `a` | Denominator-eligible `P_a` | Numerator-eligible `P_(a,d)` |
| --- | --- | --- |
| `1` | `1/2` | `1/4` |
| `2`, `5q_0`, or `10q_0` | `1/4` | `1/8` |

These count primes where the indicated obstruction is permitted by
the squareclass restrictions. They do not assert that every such
prime occurs for a rational elliptic point.

**Proof.** The four rational squareclasses

\[
 -3,\qquad 2,\qquad 5q_0,\qquad d
\]

are independent. Indeed the factors `7,739` first force the
coefficient of `d` to vanish in any square relation; the prime
`q_0` then removes `5q_0`, the prime `2` removes `2`, and a
negative sign removes `-3`. Their multiquadratic extension therefore
has degree sixteen. Its Galois group consists of all four sign
choices. Chebotarev gives each sign pattern density `1/16`.

For `p!=2,3`, the condition `p=2 mod3` is `(-3/p)=-1`.
For `a=1`, this is the only condition in `P_a`; the sign of `d`
halves the density for `P_(a,d)`. For each other `a`, its sign is
an independent nontrivial linear combination of the signs of
`2,5q_0`. Requiring it to be positive halves the density once, and
requiring the sign of `d` to be positive halves it again. Omitting
finitely many ramified or small primes does not change these
densities. \(\square\)

The exact density import is Milne, *Algebraic Number Theory*,
version 3.00, Theorem 8.31 and Corollary 8.32, printed p.140
(PDF p.148), in the retained primary text linked above. Theorem 8.31
states the conjugacy-class density; here all classes are singletons.

## 4. Scope and exact controls

The smaller prime sets leave the all-primes obligation of
[EFA3](ELLIPTIC_EXCEPTION_FINITE_LOCAL_ASYMPTOTICS.md) intact.
In a fixed class, its large-prime tail need only inspect denominator
valuations at `P_a` and numerator valuations at `P_(a,d)`, as in
(2.2). The two small-prime tests at `2,5` and the test at `3`
remain separate.
Prime eligibility densities are not elliptic-point densities. They
give neither independence of the local conditions nor an Euler
product for the limiting Haar measure, and they do not prove
Hasse infinitude on this slice.

[EKX1--EKX2](ELLIPTIC_EXCEPTION_KUMMER_BRANCH_EXCLUSION.md) make
this distinction concrete: a fixed positive-density Chebotarev set
passes the numerator eligibility test in every class, but neither
finite branch residue belongs to the rational reduction image.
For every rational point, the reduced numerator `N` in (2.2) is
prime to each prime in that set. These numerator tests can therefore
be erased for the whole slice; the denominator condition remains.
This uses the independent nonrational branch point and Kummer theory,
in addition to the quadratic character restrictions above.

Run from the repository root:

```sh
.venv/bin/python research/arithmetic/verify_elliptic_exception_prime_eligibility.py --write-receipt
```

The [checker](verify_elliptic_exception_prime_eligibility.py) verifies
the exact branch identities and squarefree constants, the independence
and sign-count arithmetic, the examples at `11,17,q_0`, and complete
finite-field branch scans at primes `11<=p<=251`, `p=2 mod3`.
It also checks (2.3) on reduced coordinates of a bounded rational
elliptic subgroup. These finite controls do not prove the unrestricted
local or density statements; their proofs are given above, with the
stated Hensel and Chebotarev imports. The
[receipt](elliptic_exception_prime_eligibility_checks.json) records
the exact domains and source hashes. No formal verification or
external review is claimed.
