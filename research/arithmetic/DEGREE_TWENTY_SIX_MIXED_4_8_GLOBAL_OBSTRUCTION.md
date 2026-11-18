# A global obstruction to the reciprocal degree-26 split (4,8)

There is no normalized good rational degree-26 reciprocal pencil with
aggregate sign degrees `(4,8)`, in either orientation. The proof excludes
every point over `Q_29` of a stronger necessary nine-equation closure,
including degenerate parameters and arbitrary denominators.

The complete finite-field census has four points in a weak seven-equation
closure. The necessary square divisibility removes three. At the sole
remaining residue point, a literal left-null vector proves that no lift
modulo `29^2` exists. Uniform denominator control is essential: a residue
test alone would not exclude models with bad reduction or denominators.

Together with the [global (2,10) obstruction](DEGREE_TWENTY_SIX_MIXED_2_10_GLOBAL_OBSTRUCTION.md),
this leaves only the balanced `(6,6)` case among rational reciprocal models.
General nonreciprocal degree-26 models are outside these exclusions.

## 1. A seven-variable necessary system

Use the normalized notation of the
[coupled-operator theorem](MIXED_RECIPROCAL_COUPLED_ODE.md): `h=6`, `L=25`,
`r=1/rho`, and

\[
 X+2(w+2)X'=25(D_6-rD_5)P Q,\qquad
 Y+2(w-2)Y'=25(C_6-rC_5)P Q.
\tag{1}
\]

Here `C_0=D_0=1`, `C_1=w-1`, `D_1=w+1`, and each sequence satisfies
`Z_(j+1)=w Z_j-Z_(j-1)`. The polynomials `X,Y` are monic of degree twelve.
Every good model satisfies `P^2|X`, `Q^2|Y` and `G=F(r)/(2r)=0`.
The apparent division in `G` cancels identically.

Unequal sign degrees make `lambda` a square in the coefficient field.
Changing the sign of `rho` if necessary, write

\[
 P=w^2+cw+d,\qquad Q=w^4+q_3w^3+q_2w^2+q_1w+q_0,
 \qquad Y=Q^2 V,
\]

\[
 V=w^4+v_3w^3+v_2w^2+v_1w+v_0.
\]

The equation for `Y` becomes

\[
 4(w-2)VQ'+\{V+2(w-2)V'\}Q
                  =25(C_6-rC_5)P.
\tag{2}
\]

Its leading coefficient agrees automatically. Matching degrees seven,
six, five and four determines `q_3,q_2,q_1,q_0` uniquely, dividing only
by `21,17,13,9`, respectively. The corresponding triangular linear
system has constant determinant `21*17*13*9`. Thus the four `q_j` are
polynomials over `Q` in

\[
                      x=(v_3,v_2,v_1,v_0,c,d,r),
\tag{3}
\]

with no variable denominators. The four remaining coefficients of (2),
in degrees zero through three, are the `Y` compatibility equations.

Construct `F,X,Y` from `P Q` by the exact Euler formula in the coupled
note. For the weak system impose the four compatibility equations,
the two coefficients of `rem_P X`, and `G=0`. For the strong system
replace `rem_P X=0` by all four coefficients of `rem_(P^2) X=0`.
Each equation is retained as its primitive integer numerator.

The [complete integer input](degree26_mixed_4_8_system.json) fixes all
coefficients and their order. Both systems have the same seven variables;
the strong system has nine equations. The weak equations have degrees
`5,5,5,5,11,11,12` and term counts `73,90,88,60,652,712,240`.

No residue squarefreeness is assumed. Square divisibility is first
deduced for an actual good characteristic-zero model and then retained
as a polynomial identity through reduction. Discarding degenerate
residue points without this argument would be invalid.

## Lemma 1. Every strong point over Q_29 is integral

Every solution of the strong nine-equation system over `Q_29` has all
seven coordinates (3) in `Z_29`.

**Proof.** The four compatibility equations make `Q^2 V` the unique
polynomial solution of the second operator equation (1). Indeed, its
operator is (2), and the eigenvalues `2j+1`, for `0<=j<=12`, are all
nonzero. Thus `Y=Q^2 V` without any nondegeneracy assumption. Strong
`P^2|X` implies weak `P|X`, so the coefficients of `P,Q` and `r`
satisfy the raw normalized system of the
[uniform finiteness and integrality theorem](RECIPROCAL_PARAMETER_FINITENESS_AND_INTEGRALITY.md),
Theorem 3. Since `29>25`, that theorem makes these coefficients integral
and makes `r` a unit.

All Euler and operator coefficients have denominators supported on
primes at most 25, hence `Y` is integral. Dividing it by the monic
polynomial `Q^2` gives an integral quotient `V`; its remainder is zero
because `Y=Q^2 V`. Therefore every coordinate in (3) is integral. This
also applies to degenerate strong points. ∎

## Lemma 2. The complete F29 census

The weak system has exactly the following four points over `F_29`:

\[
\begin{split}
 &(10,13,16,4,20,7,15),\\
 &(18,3,2,17,6,9,12),\\
 &(20,27,2,16,15,20,9),\\
 &(22,0,8,3,4,4,16).
\end{split}
\tag{4}
\]

Imposing the four strong square-divisibility equations leaves only

\[
                       x_*=(20,27,2,16,15,20,9).
\tag{5}
\]

**Coverage proof.** For fixed `v_3,v_2,v_1,v_0`, the four compatibility
equations are affine linear in

\[
                         (c,d,r,cr,dr).
\]

This follows directly from the linear determination of `Q` in (2), and
the exact integer input verifies every monomial. The
[C++ census](verify_degree_twenty_six_mixed_4_8_finite_field.cpp)
visits all `29^4=707281` outer tuples and enumerates every affine solution
of this four-row, five-variable system by reduced row echelon form.
Every rank-deficient consistent branch is retained. It then imposes
the identities `cr=c*r`, `dr=d*r`, recovering exactly, and once each,
the original triples `(c,d,r)`.

There are 36 consistent outer systems of rank three and 706443 of rank
four; the other outer systems are inconsistent. The first four original
equations have exactly 25764 full seven-coordinate solutions. Checking
the remaining three weak equations gives (4). Direct square-remainder
evaluation gives (5). A separate complete five-outer-variable scan,
solving the linear equations in `c,d` and retaining ranks zero, one and
two, independently gave the same candidate count and four-point set.
The retained fast algorithm includes a redundant direct evaluation of
all four compatibility equations at every candidate.

The new [checker](verify_degree_twenty_six_mixed_4_8_global.py)
independently reconstructs `Q` by solving the symbolic four-by-four
linear system, reconstructs `F,X,Y`, verifies every retained integer
coefficient, and executes the complete field census. ∎

## Lemma 3. The residue point has no lift modulo 29 squared

Order the strong equations `f_0,...,f_8` as the four compatibility
coefficients in ascending `w` degree, the four coefficients of
`rem_(P^2) X` in ascending degree, and `G`. At (5), let

\[
 J=(\partial f_i/\partial x_j)(x_*)\pmod{29},\qquad
 b=-\bigl(f_i(x_*)/29\bigr)_{i=0}^8\pmod{29}.
\]

The exact data are

\[
 J=\begin{pmatrix}
20&3&9&27&8&2&25\\
24&1&3&0&5&14&27\\
4&24&11&23&19&0&14\\
5&17&17&14&11&10&15\\
28&18&14&25&0&18&16\\
13&2&27&24&18&14&6\\
2&1&13&10&7&25&19\\
1&4&19&24&12&12&1\\
17&16&20&15&9&10&25
\end{pmatrix},\quad
 b=\begin{pmatrix}1\\14\\4\\15\\11\\23\\20\\3\\22\end{pmatrix}.
\tag{6}
\]

For the row vector

\[
                  \ell=(18,8,27,21,0,0,23,17,11)
\]

one has `ell J=0` and `ell b=1` in `F_29`. Any lift `x_*+29z` would
require `Jz=b`; the two identities contradict that requirement. Thus
there is no lift modulo `29^2`. The ranks `rank J=6` and
`rank(J|b)=7` provide an additional check, but the displayed left-null
identity itself proves inconsistency.

The [literal certificate](degree_twenty_six_mixed_4_8_global_certificate.json)
and checker recompute (6) from the integer equations, rather than
trusting the displayed ranks. ∎

## Theorem 4. No rational model of reciprocal type (4,8)

The entire strong parameter closure has no point over `Q_29`.
Consequently neither orientation of the good reciprocal degree-26
sign split `(4,8)` has a rational model.

**Proof.** Lemma 1 makes every hypothetical strong point integral.
Lemma 2 forces its residue to be (5), and Lemma 3 excludes even its
first lift. This proves emptiness over `Q_29`.

Every good model of this raw sign type supplies a strong point by
square divisibility. Unequal sign degrees force the reciprocal
parameter to be a square in the coefficient field, since a nontrivial
quadratic conjugation would exchange the sign factors. Over `Q` its
square root is rational, and changing its sign selects `deg P=2`,
`deg Q=4`. The resulting rational seven-tuple would be a point over
`Q_29`, a contradiction. ∎

## Replay and limits

```sh
.venv/bin/python research/arithmetic/verify_degree_twenty_six_mixed_4_8_global.py \
  --write-receipt research/arithmetic/degree_twenty_six_mixed_4_8_global_verification.json
```

The [receipt](degree_twenty_six_mixed_4_8_global_verification.json) pins
the symbolic checker, separate C++ source, complete coefficient input,
and literal no-lift certificate. Its complete census uses exact modular
arithmetic and includes all singular linear-system branches. The
denominator passage depends on the written uniform integrality theorem,
whose proof and imports are retained separately.

This excludes reciprocal `(4,8)`, rather than all rational degree-26
pencils. Balanced `(6,6)` and general nonreciprocal models require
separate arguments. No formal verification, independent external replay,
external human review or historical novelty is asserted.
