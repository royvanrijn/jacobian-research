# A complete-chart safe union with31 radical targets of degree at most four

17 November2026. Claim `SIC2B33GLOBALSAFEGATE20261117`.
**The safe-union presentation is proved; every moment radical membership
below remains OPEN.** This is a sufficient global route to Image, not
a claimed exhaustion theorem.

Use the complete positively graded highest-cubic chart from
[the global dimension theorem](HIGHEST_CUBIC_GLOBALIZATION_AND_GRADED_GATE.md).
Its twelve variables retain every phase with nonzero cubic harmonic,
after one of seven fixed rational rotations, amplitude and one unipotent
normalization. In particular no \(D\ne0\), eigenvalue, linear-vector,
base-discriminant or coefficient-support condition is imposed here.
The lower-degree complement is already one-sided when pure, by the
retained full(2,2) theorem, and hence has the coefficient-one marked cutoff.

Let \(Y_N\) be the normalized one-sided nullcone and \(Y_R\) the
explicit Rodrigues graph. Every point of either family is pure and safe
for every original multiplier. The full marked cutoffs are respectively
\(m>b\) and \(m>3b\), where \(b\) is the original coordinate degree.
We compute the radical ideal of their **union**

\[
 J_{\rm safe}=I(Y_N)\cap I(Y_R).
\]

It has an explicit31-element generating list with ordinary coefficient
degrees at most four and positive weights between five and ten. Thus the
following is a sufficient global certificate, for any finite \(N\):

\[
 f_i\in\sqrt{(R_3,\ldots,R_N)},\qquad i=1,\ldots,31.       \tag{1}
\]

If certified, (1) would prove full field-valued(3,3) Image, including
every \(D=0\) component, with the known sharp uniform cutoff \(m>3b\).
It is not necessary for Image: additional pure components outside this
particular union would be acceptable if their marked safety were proved.

## 1. Complete equations for the nullcone image

The [previous note](HIGHEST_CUBIC_GLOBALIZATION_AND_GRADED_GATE.md),
Section4, gives its entire polynomial incidence in one parameter \(s\).
There are eight equations in the twelve chart variables and \(s\):
the monic equation \(s^3+c_1s+c_0=0\), and the seven expressions for
\(l_w,l_t,b_{-2},b_{-1},b_0,c_{-2},c_{-1}\). Each dependent variable
has a constant pivot. Hence the incidence quotient is the polynomial
ring \(\mathbb Q[l_z,b_1,b_2,c_1,s]\), so its ideal is prime.
The projection is finite because of the monic cubic. Its image is closed
and its ideal is the elimination ideal, with no saturation needed.

This affine parameter covers the entire nullcone in this chart. The
only missing isotropic flag would be the \(w\) flag. A positive-charge
cubic for that flag has zero \(z^3\) coefficient and cannot satisfy the
normalization \(\gamma_3=1\). For every other flag, the lowering
coordinate \(s\) is finite. The two prescribed cubic coefficients solve
triangularly for the seed coefficients, and all remaining positive
linear/quadratic/cubic coefficients stay arbitrary. No denominator
case is discarded, even when the monic cubic has repeated roots.

Exact QQ Singular elimination gives67 generators for \(I(Y_N)\),
retained in [machine-readable form](highest_cubic_nullcone_targets.json).
The parameter variable is then explicitly removed from the ambient ring
before dimensions and ideal intersections are computed. Its image
dimension is five.

## 2. Rodrigues and the union

The other ideal has the ten generators

\[
\begin{gathered}
 c_{-2},c_{-1},c_0,c_1,b_{-2},b_{-1},b_0,\\
 l_w+b_1^2/2,\quad l_t-b_1b_2/2,\quad l_z-b_2^2/4.
\end{gathered}
\]

Its quotient is \(\mathbb Q[b_1,b_2]\), so this is prime as well.
The graph is the explicit family
\(z^3-2uzt-2vz^2-2u^2w+2uvt+v^2z\).
The existing Rodrigues shear theorem and its uniform polynomial-tail
extension cover every point, including \(u=0\).

The exact ideal intersection initially returns77 generators. Sorting by
positive weight and removing generators already in the preceding ideal
leaves the preferred31-element list. The
[input and full output](highest_cubic_safe_union_output.log),
[preferred targets](highest_cubic_safe_union_gate.json), and
[generator-change certificate](highest_cubic_safe_union_generator_change.json)
are retained. The latter gives both containments by direct polynomial
identities:31 nonzero coefficients express the preferred list in the77,
and333 express the77 in the31. These identities are independently
multiplied over \(\mathbb Q\). The number31 need not be optimal for a
set-theoretic gate; no arithmetical-rank minimality is asserted.

The preferred weight distribution is

| Weight | Number of targets |
|---:|---:|
| 5 | 1 |
| 6 | 5 |
| 7 | 11 |
| 8 | 10 |
| 9 | 3 |
| 10 | 1 |

For example, its initial targets include

\[
\begin{aligned}
 f_1&=b_2c_{-1}+3b_1c_0-3b_0c_1-6b_{-2},\\
 f_2&=2l_zc_{-1}+9l_tc_0+6l_wc_1,\\
 f_3&=2b_2c_{-2}+2b_1c_{-1}+9b_0c_0+2b_{-1}c_1,\\
 f_4&=3b_0^2+2b_{-1}b_1+4b_{-2}b_2,\\
 f_5&=l_t^2+2l_wl_z.
\end{aligned}
\]

The last two are scalar multiples of the quadratic norm of \(B\) and
the norm of the linear vector. They vanish on the known safe union.
It is **not proved here** that they vanish on every pure phase.

## 3. Exact finite ansatz for each proposed power

For a target \(f_i\) of weight \(d_i\) and a chosen exponent \(e\),
the homogeneous identity

\[
 f_i^e=\sum_{n=3}^{\lfloor ed_i/3\rfloor}A_nR_n
\]

uses only multiplier monomials of weight \(ed_i-3n\). This is a
finite rational linear system, with no unknown witness-degree cutoff
once \(e\) is chosen. If a smaller moment prefix is imposed, truncate
the displayed range accordingly. This positive grading is stronger than
the mod12 filter of the normalized cyclic quotient.

Neither the necessary exponents nor a sufficient finite prefix are
known. Separate tests \(1\in I_N+(1-v_if_i)\) are equivalent radical
tests; they must not be combined into a single test inverting all targets.
The inverse variable has negative weight, so positivity does not itself
give a finite multiplier ansatz for an unrestricted inverse-variable test.

The certificate records **zero proved moment memberships**. Constructing
this ideal of known safe phases supplies the target side of the global
gate, not the reverse containment or an upper bound on the pure quotient.

There is already an exact lower bound on the first target's possible
power: **\(e\ge7\)**, even if arbitrarily late moments are allowed.
This guards against treating low target degree as a small certificate
bound. Restrict to \(F=z^3+(A+Bt)w^2\). Circular balance forces every
nonzero moment index to be \(5s\), with coefficient degree \(3s\).
The primitive fifth and tenth moments are

\[
 5A^3+AB^2,\qquad
 8091A^6+4495A^4B^2+465A^2B^4+5B^6.
\]

Their ideal has the verified standard basis
\(5A^3+AB^2,\ 2759A^2B^4-125B^6,\ AB^6,\ B^8\), length18
and top ordinary degree7. Hence all later moments already belong to it.
The normal forms of \(A^6,A^7\) are respectively \(5B^6/2759\ne0\)
and zero. Since \(f_1\) restricts to \(-6A\), no exponent at most six
can give a global ordinary moment membership. The
[power-obstruction certificate](highest_cubic_safe_union_power_obstruction.json)
includes explicit ideal-generator changes; the checker verifies all six
S-pairs by exact division. This is a certificate lower bound, not a new
pure-point family or proof that the seventh power works globally.

## Evidence and boundary

Run `python verify_highest_cubic_safe_union.py`. It regenerates the exact
primitive integer incidence input, replays QQ elimination/intersection
and the generator reduction, checks all31/67/77 polynomial vectors and
weights, evaluates them on the full parameterized families, and directly
multiplies both generator-change matrices. The
[receipt](highest_cubic_safe_union_verification.json) records hashes.
Elimination/intersection remain an explicit Singular dependency; there
is no msolve dependency.

The initial rational-pretty-print source failed because Singular parsed
a rational denominator as part of an exponent. It is preserved with its
error output and replaced by primitive integer equations. The first
correct elimination worked in a ring still containing the now-free
parameter \(s\), and consequently printed dimension six; that extended
ring checkpoint is preserved. The accepted source explicitly maps to
the twelve-variable ring and prints dimension five. Neither failed or
extended-ring output is used as proof of a moment membership.

The global safe-cover test, all actual pure components and their generic
marked safety remain open. The dimension theorem and this safe-union
presentation do not complete the full Image goal. No independent
specialist, external or formal verification is claimed.
