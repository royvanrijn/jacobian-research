# Lower exponent three: closure at gaps five and seven

Internal proof, 21 November 2026. Base repository revision
`19eb9912fb60142f8b620dac24896931d4ff3e46`; the additional ordinary-digit
lemma from the current investigation is an explicit input below.

## 1. Statements and scope

**Theorem A (ordinary support three/eight).** Let `K` be any field of
characteristic zero, with the ordinary factorial functional
`Lambda(X^j)=j!`. For every `N>=2` and every `A,B in K` not both zero,
the two moments of

\[
 f(X)=X^3(A+BX^5)
\]

at orders `N-1,N` are not both zero.

**Theorem B (ordinary support three/ten).** Under the same hypotheses,
the moments of `X^3(A+BX^7)` at orders `N-1,N` are not both zero for
any `N>=2`.

The proof first establishes the following intermediate confinement:
a possible gap-seven consecutive zero pair would have to occur at

\[
 N\equiv446\pmod {2401}.
 \tag{1}
\]

At such a terminal, if the two monic reciprocal moment polynomials
have a nonconstant gcd over `Q`, that gcd is an irreducible quadratic
over `Q`. Its two roots are simple. Writing it as `u^2+b u+c`, one has

\[
 b\equiv-7\pmod {49},\qquad c\equiv147\pmod {343}.
 \tag{2}
\]

Section 9 excludes that entire class by a uniform seven-case exact
certificate. Thus (1)–(2) describe an intermediate obstruction, which
is fully removed; they are not a remaining exception set.

These are ordinary **univariate binomial** statements. They do not
prove unrestricted Strong Factorial, multivariate two-term Strong
Factorial, or the corresponding conclusions for arbitrary Gamma shapes.

All positivity used below is the exact rational sum-of-squares identity
in Section 7, after descent to a rational common root. No determinant
positivity is assumed, and no ordering of an arbitrary coefficient
field is used. There is no Riemann-hypothesis input.

## 2. The retained actual-jet identities

For an odd prime `p>=5`, put `N=n+1`, `gamma=3n+1`, and

\[
 T_{n,\gamma,p}(t)=\sum_{k=0}^n\binom nk(\gamma)_{pk}t^k,
 \quad C=\gamma+pn+p=(p+3)N-2,\quad D=C+1=(p+3)N-1.
\]

The actual pair is `T_(n,gamma,p),T_(n+1,gamma+3,p)`. Its monic
lower reciprocal is

\[
 P_n(u)=u^nT_{n,\gamma,p}(u^{-1}).
\]

All roots of `P_n` are nonzero, and all its coefficients are integers.
For nonzero `A,B`, the exact connection to the phases is

\[
 \Lambda((X^3(A+BX^p))^m)
       =A^m(3m)!T_{m,3m+1,p}(B/A).
\]

A common phase ratio gives `u=A/B`. If either coefficient vanishes,
the nonzero phase is a monomial and all its positive factorial moments
are nonzero. These observations handle all coefficient pairs in the
statements of Theorems A and B.

We use the exact mixed concomitant from
[DIAGONAL_PRIME_GAP_SFC.tex](DIAGONAL_PRIME_GAP_SFC.tex),
`thm:offset-rank`, with its full-jet justification in
[LOWER_TWO_PRIME_GAP_SFC.tex](LOWER_TWO_PRIME_GAP_SFC.tex),
proof of `thm:mixed`, and
[HIGHER_GAP_ADJOINT_RANK.tex](HIGHER_GAP_ADJOINT_RANK.tex),
proof of `thm:rank`. Define, always at the same lower shape `gamma`,

\[
 B_0=1,\qquad
 B_{m+1}=(z+pm+\gamma+p-1)B_m-zB_m',
 \qquad L_j=z^j(D_z-1)^j B_{n+1}.
\]

For a polynomial in `z`, its cyclic column is its remainder coefficient
column modulo `z^p-u`, in the basis `1,z,...,z^(p-1)`. Let `Z(u)` be
the `p` by `p` matrix whose column `m`, `1<=m<p`, is the cyclic column
of `-L_(p-1-m)`, and whose last column is that of `B_n`.

At a nonzero lower root `u`, the actual nonzero rising jet
`(K_1,...,K_p)` lies in the kernel of `Z(u)`. Here
`K_j=(E)_j R`, `R=z^(pn)T_(n,gamma,p)(z^(-p))`,
`E=gamma+pn-zD_z`, and `z^p=u`. The full jet cannot vanish at a
nonzero ordinary point of its monic differential equation. No lower
root simplicity is assumed.

For offset three, the upper equation at a lower root is exactly

\[
 K_p=\kappa u\bigl((C+2)K_1+K_2\bigr),\qquad
 \kappa=\frac{p}{(p+3)C D}.
 \tag{3}
\]

The canceled denominator identity follows from
`C+2=(p+3)N`. All factors removed in deriving (3) are nonzero at
the actual positive integer shapes.

The retained pointwise Hasse criterion in
[PRIME_PLUS_ONE_EXPONENT_GAP_SFC.tex](PRIME_PLUS_ONE_EXPONENT_GAP_SFC.tex),
`thm:criterion`, together with the extra ordinary digit in
[STEPANOV_NATIVE_PRIME_OFFSET_COROLLARIES.md](STEPANOV_NATIVE_PRIME_OFFSET_COROLLARIES.md),
implies that any common root for `p>=5` has

\[
 p^3\mid D.
 \tag{4}
\]

Indeed `sigma_3(x)=3(x+1)`, and only the central zero-Hasse class
survives. The rational-root exclusion in that retained criterion can,
for these ordinary pairs, use Section 7's exact identity in place of
the real Gamma integral. Thus this ordinary specialization requires
no assumed positivity principle. Formula (4) also implies `n>=2`.

## 3. A fixed unit minor and the actual normalized scalar

**Lemma 1 (fixed normalization).** Suppose `gamma=-1 mod p` and `n>=1`.
Take rows `1,...,p-1` of `Z(u)`, omitting its constant-coefficient row,
and take its first `p-1` columns. The resulting square matrix `M(u)`
has `det M(0)` a `p`-adic unit.

Consequently there is a unique rational-function vector

\[
 w(u)=(w_1(u),\ldots,w_{p-1}(u),1)
\]

annihilated by those fixed rows. Every component belongs to
`Z_p[[u]]`, with convergence at every `v_p(u)>0`. At an actual lower
root in that disk, this is the actual jet normalized by `K_p`.

**Proof.** Write `A_0=z-zD_z`. On the original polynomial space,

\[
 B_n=\Pi_n(A_0)1,\qquad
 \Pi_n(X)=\prod_{r=0}^{n-1}(X+\gamma+p-1+pr),
\]

and the recurrence for `L_j` gives

\[
 L_j=(-1)^j(A_0)_j(A_0+C-1)B_n.
 \tag{5}
\]

At `u=0`, `A_0` is well-defined on `Z_p[z]/(z^p)`. The map
`f(X) -> f(A_0)1` is a module isomorphism

\[
 \mathbb Z_p[X]/\chi_0(X)\longrightarrow
                   \mathbb Z_p[z]/(z^p),\qquad
 \chi_0(X)=\prod_{r=0}^{p-1}(X+r).
 \tag{6}
\]

Its basis matrix is triangular with diagonal entries one, since
`A_0^j1` is monic of degree `j` for `j<p`.

For a homogeneous vector in the nullspace of `M(0)`, append `w_p=0`
and set

\[
 F(X)=w_p-(X+C-1)
       \sum_{j=0}^{p-2}(-1)^j w_{p-1-j}(X)_j.
 \tag{7}
\]

Its degree is at most `p-1`. The retained rows being zero says
`F(X)Pi_n(X)=epsilon` modulo `chi_0`, for a constant `epsilon`.
Reduce modulo `p`: `Pi_n(X)=(X-2)^n` and `1-C=2`. Evaluation at
`X=2` gives `epsilon=0`. Evaluation at every other element of `F_p`
gives `F=0` there; also `F(2)=w_p=0`. A polynomial of degree at most
`p-1` with these `p` distinct roots is zero. Formula (7) and the
triangular rising basis then force all `w_m=0`. Thus `M(0)` is
invertible modulo `p`.

The actual matrix has integer polynomial entries in `u`. Its unit
constant determinant gives the claimed integral power series and
convergence on `v_p(u)>0`. If an actual nonzero lower jet had `K_p=0`,
the same invertible minor would force all its other components to
vanish. Hence `K_p!=0`, and the normalized actual jet is the unique
solution just constructed.

**Characteristic-zero safeguard.** `A_0` does **not** descend to the
characteristic-zero quotient `z^p=u` at nonzero `u`: its Euler term
differentiates `z^p` by `p`. Formula (5) is used on the original
polynomials and then projected into cyclic columns. Argument (6) is
only at `u=0`. The later argument with variable `u` is modulo `p`,
where the operator does descend. No characteristic-zero cyclic
operator identity at nonzero `u` is asserted.

Define the actual normalized scalar

\[
 S(u)=(C+2)w_1(u)+w_2(u)=\sum_{j\ge0}S_j u^j.
 \tag{8}
\]

Every `S_j` belongs to `Z_p`. At a common root, (3) becomes

\[
                 1=\kappa u S(u).
 \tag{9}
\]

## 4. Uniform first coefficients of the scalar

**Lemma 2.** Suppose `n>=2`, `gamma=-1 mod p`, and `p^2 | D=C+1`.
Then

\[
 S_0\equiv p\pmod {p^2},\qquad S_1\equiv-1\pmod p.
 \tag{10}
\]

In particular, for every integral `x` in a complete algebraically
closed extension of `Q_p`,

\[
                 S(px)/p\equiv1-x\pmod{\mathfrak m}.
 \tag{11}
\]

**Proof of the constant coefficient.** Work at `u=0` over `Z_p`, using
(6). For the normalized vector, `F(1-C)=1`. The fixed-row equations give

\[
 F(X)\Pi_n(X)=\epsilon\pmod {\chi_0(X)}.
\]

At the eigenvalue `X=2-p`, every factor of `Pi_n` is divisible by `p`,
so `epsilon` is divisible by `p^n`. At every other eigenvalue `-r`,
`0<=r<p`, the value of `Pi_n` is a unit. Hence `F(-r)=0 mod p^n`
at all these other eigenvalues. Put

\[
 H(X)=\prod_{\substack{0\le r<p\\r\ne p-2}}(X+r).
\]

The points `1-C` and the `p-1` retained eigenvalues are distinct
modulo `p`. Their Vandermonde determinant is a unit, so interpolation
gives the coefficientwise identity

\[
 F(X)\equiv H(X)/H(1-C)\pmod {p^n}.
 \tag{12}
\]

The denominator is a unit. Let `f_(p-1),f_(p-2)` be the top two
coefficients of `F`. Expanding the rising basis in (7) gives exactly

\[
 S_0=\left(2C+1+\frac{(p-2)(p-3)}2\right)f_{p-1}-f_{p-2}.
 \tag{13}
\]

The next coefficient of the monic polynomial `H` is
`p(p-1)/2-(p-2)`. Equations (12)–(13) therefore give

\[
 S_0\equiv\frac{2D-p}{H(2-D)}\pmod {p^n}.
 \tag{14}
\]

Since `H(2-D)=-1 mod p`, `p^2|D`, and `n>=2`, this is
`S_0=p mod p^2`.

**Proof of the linear coefficient.** Now work in
`F_p[u]/(u^2)`. In characteristic `p`, `A_0` does descend to the
cyclic quotient, with characteristic polynomial

\[
 \chi_u(X)=X^p-X-u.
\]

The same cyclic-basis isomorphism is valid, and
`Pi_n(X)=(X-2)^n`. The root of `chi_u` reducing to two is
`rho=2-u`; all other roots are the other elements of `F_p` minus
`u`. They remain pairwise distinct units apart. Since `n>=2`,
`Pi_n(rho)=(-u)^n=0`. The fixed-row relation again forces its
constant right side to be zero. At the other roots `Pi_n` is a unit.
Consequently

\[
 F(X)=H_u(X)/H_u(2),\qquad
 H_u(X)=\chi_u(X)/(X-\rho).
\]

The division is monic polynomial division, not division by a value
that might vanish. Explicitly

\[
 H_u(X)=\sum_{i=0}^{p-1}X^{p-1-i}\rho^i-1.
\]

For odd `p`, direct substitution of `rho=2-u`, `u^2=0`, gives
`H_u(2)=-1`: the linear term contains the factor
`sum_(i=0)^(p-1)i=p(p-1)/2`. Thus the top coefficients of `F`
are `-1` and `-rho`. Formula (13), now modulo `p`, gives
`S(u)=-u mod u^2`. This proves the second assertion of (10).
Integral higher coefficients give (11), with convergence justified
by Lemma 1.

## 5. The lower unit-root residues and valuation balance

Every coefficient at degree `k` of `T_(n,gamma,p)` has valuation at
least `k`, by its `k` complete blocks of `p` integer factors. Hence
every lower reciprocal root has `v_p(u)>=1`, and

\[
 \widetilde P_n(x)=p^{-n}P_n(px)
\]

is a monic polynomial over `Z_p`.

Assume (4). Put `n_0=n mod p`, with `0<=n_0<p`, and
`q=(gamma+1)/p=D/p-N`. In each block of `p` factors, the unit
factors multiply to `-1 mod p`, while the divided multiple of `p`
is the corresponding factor of `(q)_k`. Therefore

\[
 p^{-k}(\gamma)_{pk}\equiv(-1)^k(q)_k\pmod p.
\]

Here `q=-N=-n_0-1 mod p`. Terms with `k>=p` vanish modulo `p`.
For `k<p`, reduction of the binomial coefficient gives
`binom(n,k)=binom(n_0,k) mod p`. It follows that

\[
 \overline{\widetilde P_n}(x)
    =x^{n-n_0}Q_p(x),\qquad
 Q_p(x)=\sum_{k=0}^{n_0}\binom{n_0}{k}
                  (n_0+1)_{\underline{k}}x^{n_0-k}.
 \tag{15}
\]

For the two primes of interest,

\[
 \begin{array}{c|c|c}
 p&n_0&Q_p(x)\\\hline
 5&1&x+2\\
 7&4&(x-3)(x+3)(x^2-x+3).
 \end{array}
 \tag{16}
\]

Both polynomials are squarefree and neither vanishes at one.
The quadratic factor in the second row has discriminant three,
a nonsquare modulo seven.

A unit root `x` of the scaled lower polynomial reduces to a root of
`Q_p`; in particular its residue is not one. Equation (11) therefore
gives `v_p(S(u))=1` when `v_p(u)=1`. If `v_p(u)>1`, put `x=u/p`;
then `x` reduces to zero, and (11) gives the same conclusion.
Thus **every actual lower root** in the present central classes has

\[
                         v_p(S(u))=1.
\]

Set `t=v_p(D)>=3`. Since `p+3` and `C` are units, (9) forces

\[
                         v_p(u)=t-2
                         \quad\text{at a common root}.
 \tag{17}
\]

No global upper bound on the valuations of all lower roots is used.

## 6. Uniqueness on the deep disk

**Lemma 3.** The function `u -> uS(u)` is injective on the entire disk
`v_p(u)>1`.

**Proof.** For distinct `u,v` in that disk, its divided difference is

\[
 S_0+(u+v)S_1+
       \sum_{j\ge2}S_j\sum_{r=0}^j u^{j-r}v^r.
\]

The first term has valuation exactly one. Every other term has
strictly larger valuation; the convergent tail cannot cancel the
first term. Thus the divided difference has valuation one and is
nonzero. This proves injectivity.

When `t>=4`, (17) places every common root in this disk. The exact
upper equation (9), whose scalar `kappa` is fixed by the parameters,
therefore permits at most one common root. All roots of every
irreducible rational common factor obey the same constraints in the
chosen algebraic closure of `Q_p`. Characteristic-zero separability
forces that irreducible factor to be linear. Section 7 excludes the
resulting rational root. Consequently **all strata `t>=4` are excluded
at both prime gaps five and seven**.

This argument uses a root-count bound, not a claim that a generic
Wronskian is nonzero at a particular root. A repeated irreducible
factor would not alter the rational-root exclusion; no unproved
squarefreeness of the full lower polynomial is required.

## 7. Exact rational-root exclusion

For a polynomial `g(X)=sum_j g_j X^j` over `Q`, Vandermonde's identity
gives the exact Gram decomposition

\[
 \Lambda(g^2)
  =\sum_{k\ge0}
       \left(\sum_{j\ge k}g_j j!\binom jk\right)^2.
 \tag{18}
\]

Indeed `(i+j)!=i!j! sum_k binom(i,k)binom(j,k)`, and substituting
this identity into `sum_(i,j)g_i g_j(i+j)!` gives (18). If `g` has
degree `m` and is nonzero, the summand at `k=m` is `(g_m m!)^2!=0`.
A sum of rational squares with this summand is strictly positive.

If a common reciprocal root `u` were rational, the nonzero rational
phase `f=X^3(u+X^p)` would have both consecutive moments zero. One
of their orders is even, say `2r`. Taking `g=f^r` contradicts (18).
This is an explicit consequence of a proved algebraic identity over
`Q`, not an assumed positivity property over the user's coefficient
field. Once rational coprimality is proved, a rational Bezout identity
remains valid over every characteristic-zero field.

For clarity, positivity of the same identity over arbitrary complex
coefficients would be false: `g=X-1+i` is nonzero but `Lambda(g^2)=0`.
The rational descent preceding (18) is indispensable.

## 8. The first valuation stratum and completion of the theorems

It remains to treat `t=3`. Equation (17) says every common root has
`v_p(u)=1`. At any root of `Q_p`, the reduction of the derivative of
the scaled lower polynomial is
`x^(n-n_0) Q_p'(x)`, which is nonzero. The polynomial divided
difference therefore shows that each residue root of `Q_p` has at
most one lower root above it, and that root is simple.

For `p=5`, `Q_5=x+2` supplies at most one possible common root.
Rational conjugacy and (18) exclude it, exactly as in Section 6.
Together with the retained reduction (4) and the deep-disk argument,
this proves Theorem A at every terminal. No finite terminal scan is
an input to this closure.

For `p=7`, write `D=p^3 d`, with `d` a unit. Reducing (9) at `u=px`
and using (11) gives

\[
                        x^2-x-3\bar d=0.
 \tag{19}
\]

The exact gcds over `F_7` of (19) with the polynomial in (16) are

\[
 \begin{array}{c|c}
 \bar d&\gcd(Q_7,x^2-x-3\bar d)\\\hline
 2&x-3\\
 4&x+3\\
 6&x^2-x+3,
 \end{array}
\]

and the gcd is one for the other nonzero residues. In the first two
rows there is at most one common root, excluded by rational descent
and (18). Hence `d=6 mod 7` is necessary. This is

\[
 10N-1\equiv6\cdot343\pmod {2401},
\]

equivalently (1). Only the two simple roots reducing to
`x^2-x+3` remain possible. Any rational common factor is therefore
an irreducible quadratic; rational linear factors have already been
excluded. Gauss's lemma makes its monic coefficients integers. Its
two roots have valuation one, so after writing it as `u^2+b u+c`,
the normalized factor `x^2+(b/7)x+c/49` is integral and reduces to
`x^2-x+3`. This proves the intermediate confinement (1)–(2), without
using the certificate in the next section.

## 9. A complete certificate for the last gap-seven branch

The full written exclusion, explicit lifts, and exact coefficient
witness are retained in
[STEPANOV_SEVENTH_GAP_HENSEL_CERTIFICATE.md](STEPANOV_SEVENTH_GAP_HENSEL_CERTIFICATE.md).
Here is its precise connection to the preceding confinement.

Put `z=7/u=1/x`. The actual normalized moments

\[
 F_m(z)=T_{m,3m+1,7}(z/7)
       =\sum_{k=0}^m\binom mk\frac{(3m+1)_{7k}}{7^k}z^k
\]

have integer coefficients. After one factor of seven is removed
from each block, the product still contains `k` consecutive integers.
Thus the coefficient of `z^k` has valuation at least `floor(k/7)`.
Modulo `7^4`, all terms of degree at least 28 vanish at integral `z`.

Each of the remaining 28 coefficients is periodic in `m` with period
`7^5`, modulo `7^4`. For the binomial coefficient this follows from
Vandermonde and `v_7(binom(7^5,j))>=4` for `1<=j<=27`. For each
normalized rising block, shifting `m` by `7^5` changes its divided
multiple of seven by `3*7^4` and its other factors by `3*7^5`.
Therefore the entire class `N=446 mod 2401` is covered by exactly
seven representatives `N_i=446+2401i`, `0<=i<=6`.

In the unramified quadratic ring with `theta^2-theta+3=0`, the
remaining residue is `x=theta`, hence `z=theta^(-1)`. The lower
polynomial has unit derivative there. The companion certificate
gives a lower zero `z_i` modulo `2401` and verifies that the upper
value divided by `343`, modulo seven, has the following coordinates
in the basis `1,theta`:

| `i` | `N_i` | Upper value divided by `343`, modulo seven |
|---:|---:|---|
| 0 | 446 | `(4,6)` |
| 1 | 2847 | `(0,4)` |
| 2 | 5248 | `(3,2)` |
| 3 | 7649 | `(6,0)` |
| 4 | 10050 | `(2,5)` |
| 5 | 12451 | `(5,3)` |
| 6 | 14852 | `(1,1)` |

None is zero. For a hypothetical exact lower root in the same
residue disk, the polynomial divided difference is a unit, so its
distance from `z_i` has valuation at least four. The upper value is
therefore the same nonzero residue modulo `7^4`. This excludes a
common root, even in an arbitrarily ramified extension. The conjugate
disk is excluded by conjugating the same identities.

The coefficient period proves that these seven cases cover every
terminal in (1). No other common roots remain after Sections 2–8.
This completes Theorem B. Both rational coprimality conclusions
extend to every characteristic-zero field by Bezout identities.

## 10. Verification, failed shortcuts, and remaining scope

The accompanying
[uniform-jet checker](verify_stepanov_lower_three_prime_closure.py)
constructs the **actual** integer `B_m`
recurrence and actual cyclic columns, selects the stated fixed minor,
and solves its truncated series by modular matrix inversion. It checks
the unit denominator and recomposes every retained-row equation. It
also verifies the symbolic top-coefficient calculations, the
finite-characteristic interpolation ingredients, the residue gcds,
the intermediate terminal class, and the exact Gram identity with its
complex countercontrol. Finite series examples guard the formulas;
the all-index assertions follow from Lemmas 1–3, not extrapolation.
The separate
[Hensel verifier](verify_stepanov_seventh_gap_hensel.py) checks the
complete seven-case certificate used in Section 9.

The [uniform-jet receipt](stepanov_lower_three_prime_closure_verification.json)
records a passing internal replay with 13 actual series cases, including
both primes' deep valuation strata, together with the symbolic and
residue checks. A separately implemented
[matrix-ring verifier](verify_stepanov_seventh_gap_matrix.py) reconstructs
the 28 coefficients by direct rising products, evaluates them in a
two-by-two integer matrix representation of the quadratic extension,
and checks both conjugate residue disks. Its
[receipt](stepanov_seventh_gap_matrix_verification.json) also passes.
These are internal exact replays. The uniform scalar proof and
uniqueness argument were additionally reconstructed by the parallel
determinant reviewer and the root reviewer; this is internal AI review,
not external refereeing or formal verification.

The investigation ruled out several shortcuts:

* Plain Newton-polygon separation of the actual lower and upper
  polynomials did not occur in the tested central progressions.
* A hypothetical uniform bound `v_p(u)<2` is unnecessary for the
  proof and is not asserted. The deep-disk injectivity replaces it.
* Residue compatibility alone leaves a further class even at gap
  five. The one-root argument closes it; a numerical residue scan
  alone would not.
* The variable-`u` cyclic operator may be used modulo `p`, but cannot
  silently be used over characteristic zero. The fixed minor is
  always formed from actual projected polynomial columns.

The possible irreducible quadratic common factor at gap seven is
removed by Section 9. Other binomial supports and unrestricted
Strong Factorial remain outside these two closures. No global status
or existing proof source is edited by this note.
