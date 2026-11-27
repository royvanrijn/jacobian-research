# Five torsion directions on a septic Pell model

The target is a squarefree polynomial \(F=P^2+1\), with
\(P\in\mathbb Q[t]\) of degree seven, and five independent rational
Jacobian three-torsion classes represented by
\[
 2PZ_i+Z_i^2-1=H_i^3,\qquad \deg Z_i\le5,\quad\deg H_i\le4.
 \tag{S1}
\]
The smooth curve \(C:y^2=F\) has genus six. The integral noncube unit,
its regulator bound, the infinity difference of order seven, and the
conditional specialization requirements are in
[LY9](LEVIN_YAN_WILJANEN_DEGREE_TEN_MIRROR_SLICE.md#a-unit-built-into-the-replacement-model).
The [quartic construction CP1](QUARTIC_PELL_TWO_CLASS_DIRECTIONS.md)
already supplies two independent specialized class directions on a
different Pell curve. The septic example below is a control, not an
improvement of that rank.

The degree-seven five-direction target remains open. The results here
reduce its coefficient systems and exclude every model of the displayed
form with degree-preserving smooth reduction at five. They do not exclude
models with bad reduction there, higher-degree models, or other Pell curves.

## Literature comparison

[Levin–Yan–Wiljanen, Theorem 9 and its proof in Section 6](sources/levin-yan-wiljanen-large-three-rank.pdf)
construct five rational three-torsion directions in genus fourteen using
a degree-thirty cubic substitution. Section 1, equations (3)–(4), describes
the square-minus-cube descent, and explicitly distinguishes generalized
Picard groups when the degree is two or four modulo six. Their Theorem 3
assumes degree divisible by six and loses one class direction in the real
case. Their Theorem 2 assumes a rational Weierstrass point. None of these
statements directly supplies (S1), a compatible integral Pell unit, or
the desired real class-rank-five specialization. Here the divided divisor
is checked directly, including both infinities. No novelty or exhaustive
literature claim is made.

## Exact degree and factorization restrictions

**Proposition SP1 (the two degree profiles).** In characteristic zero,
a nontrivial solution of (S1) has
\[
 (\deg Z,\deg H)=(2,3)\quad\hbox{or}\quad(5,4).
 \tag{S2}
\]
Every such solution on a squarefree \(P^2+1\) gives a nonzero rational
Jacobian three-torsion class. This does not assert independence of
different solutions.

If \(Z=0\), then \(H=-1\), the trivial row. Otherwise put
\(d=\deg Z\le5\). The leading term of \(2PZ\) cannot cancel with
\(Z^2-1\), so \(3\deg H=7+d\). This gives exactly (S2).
Set \(G=P+Z\) and
\[
 \alpha=\frac{y-G}{y-P}.
\]
Its norm is \(-H^3\). The denominator is an affine unit because
\((y-P)(y+P)=1\). The two numerator factors cannot have a common
finite zero: that would give \(F=G=H=0\) and, by differentiation,
\(F'=0\). Thus the finite divisor orders of \(\alpha\) are
nonnegative multiples of three. At the infinity where \(y\sim P\)
its pole order is \(3e\), with \(e=\deg H\); at the other infinity
its order is zero. Hence \(D=\operatorname{div}(\alpha)/3\) is rational.

If \([D]=0\) over the algebraic closure, then \(\alpha\) is a
constant times a cube. Its cube root is regular on the affine curve,
with a single pole of order \(e\le4\). A nonconstant affine function
has the form \(A+B y\), with polynomial \(A,B\). If \(B=0\), a
nonconstant polynomial has poles at both infinities. If \(B\ne0\),
cancellation of its pole at one infinity forces at the other a pole
of order \(7+\deg B\ge7\). Neither permits the required cube root.
This proves nonzeroness.

There are also exact factorization charts. Since
\[
 H^3+1=(H+1)(H^2-H+1),\qquad
 \gcd(H+1,H^2-H+1)=1
 \quad\text{in }\mathbb Q[t],
\]
write \(Z=A B\), up to a nonzero scalar, with
\(A\mid H+1\) and \(B\mid H^2-H+1\). The latter polynomial is
strictly positive at every real argument, so every rational divisor
has even degree. The possibilities are

| \(\deg Z,\deg H\) | \((\deg A,\deg B)\) |
| --- | --- |
| \(2,3\) | \((2,0),(0,2)\) |
| \(5,4\) | \((3,2),(1,4)\) |

In the first quadratic chart, all rows are parametrized by a quadratic
\(Z\) and a linear \(L\):
\[
 H=ZL-1,\qquad
 P=\frac{Z^2L^3-3ZL^2+3L-Z}{2}.
 \tag{S3}
\]
In the second chart, the remainder \(R=H\bmod Z\) is linear,
\(Z=c(R^2-R+1)\), and \(H=LZ+R\), with \(c\ne0\) and \(L\)
linear. Then \(P=((H^3+1)/Z-Z)/2\) is polynomial. These charts
are for one row. Equating their \(P\)'s, or equating them with
quintic-row charts, remains necessary for a simultaneous construction.

## A complete obstruction at five

**Proposition SP2 (no smooth septic square-plus-one model at five).**
For every degree-seven \(p\in\mathbb F_5[t]\) such that \(p^2+1\)
is squarefree, the smooth curve \(c:y^2=p^2+1\) satisfies
\[
 \dim_{\mathbb F_3}\operatorname{Jac}(c)(\mathbb F_5)[3]\le4.
 \tag{S4}
\]
Consequently, if \(P\in\mathbb Z_5[t]\) has leading coefficient a
five-adic unit and \(P^2+1\) has squarefree reduction, its rational
Jacobian three-rank is at most four. Thus an integral degree-seven
candidate with unit leading coefficient at five must have
\(\operatorname{disc}(P^2+1)\equiv0\pmod5\).

Here (S4) is a finite exact computational assertion with a complete
domain reduction. An affine substitution makes \(p\) monic and kills
its degree-six coefficient: seventh powers permute \(\mathbb F_5^*\)
and seven is invertible. All possibilities are therefore
\[
 p=t^7+\sum_{j=0}^5 c_jt^j,\qquad c_j\in\mathbb F_5,
 \tag{S5}
\]
a set of exactly \(5^6=15625\) polynomials. The involution
\(p(t)\mapsto-p(-t)\) gives isomorphic curves. It leaves 7875
representatives, with weights one or two restoring the full set.
The [checker](verify_septic_pell_construction.py) computes squarefreeness
and the characteristic polynomial \(L_p(X)\) of Frobenius using PARI/GP.
The weighted counts are:

| Case | Number |
| --- | ---: |
| Nonsquarefree \(p^2+1\) | 5704 |
| \(v_3(L_p(1))=0\) | 5689 |
| \(v_3(L_p(1))=1\) | 2670 |
| \(v_3(L_p(1))=2\) | 1019 |
| \(v_3(L_p(1))=3\) | 396 |
| \(v_3(L_p(1))=4\) | 110 |
| \(v_3(L_p(1))=5\) | 33 |
| \(v_3(L_p(1))=6\) | 4 |

The last 37 polynomials require more than an order test. The multiplicity
of \(X-1\) in \(L_p(X)\bmod3\) is one, two or three in respectively
20, 8 and 9 cases. For an abelian variety over \(\mathbb F_5\),
the rational three-torsion is the kernel of Frobenius minus one on its
three-torsion vector space. Its dimension is bounded both by
\(v_3(L_p(1))\) and by that root multiplicity. The recorded integers
therefore prove (S4) for every case. The checker also checks the
reciprocal functional equation, direct \(\mathbb F_5\) first traces,
and divisibility of every good Jacobian order by seven, as required
by the infinity difference. These controls do not independently derive
every Frobenius coefficient.

Finally, under the stated integral hypotheses the affine chart and
the two infinity charts give a smooth proper model at five. Rational
three-torsion injects under good reduction, since three is prime to
five. This proves the characteristic-zero consequence. A bad reduction
of this particular equation need not mean the Jacobian lacks some other
good model; the consequence concerns exactly the displayed hypotheses.

One concrete permitted chart is
\[
 P(0)\equiv2\pmod5,\qquad P'(0)\equiv0\pmod5,
 \tag{S6}
\]
which makes the reduction of \(P^2+1\) have a multiple zero at zero.
This is compatible with characteristic-zero squarefreeness. It does
not force or establish five torsion directions.

## An integral control in the bad-at-five chart

Taking \(Z=2t^2+2t\), \(L=2t+8\) in (S3) gives
\[
\begin{split}
 H&=4t^3+20t^2+16t-1,\\
 P&=16t^7+224t^6+1168t^5+2740t^4\\
  &\hspace{10mm}+2708t^3+735t^2-190t+12.
\end{split}
\tag{S7}
\]
The polynomial \(P^2+1\) is squarefree over \(\mathbb Q\), and
(S6) holds. For every \(n\equiv0\pmod9\),
\[
 P(n)\equiv3,\quad Z(n)\equiv0,\quad H(n)\equiv-1\pmod9.
\]
Thus LY9 supplies an integral unit noncubic at three, while the section
is locally cubic at both places above three. This verifies compatibility
of this chart with those local requirements for one section. It does
not assert all-prime virtuality or specialized independence.

SP1 proves that its geometric divided divisor is nonzero. At the good
prime seven the Frobenius polynomial is
\[
 X^{12}+4X^9-17X^8+68X^7+392X^6+476X^5
       -833X^4+1372X^3+117649.
\]
Its value at one is \(119112\), divisible by three exactly once.
Prime-to-seven torsion specialization therefore proves that this
control has rational Jacobian three-rank exactly one. It is rejected
as a five-direction candidate.

## Coefficient systems and remaining construction

Up to permutation, the five rows have \(k\) quadratic profiles and
\(5-k\) quintic profiles, for \(0\le k\le5\). The checker emits
each exact system over \(\mathbb Q\):

```sh
.venv/bin/python research/arithmetic/verify_septic_pell_construction.py --system 0
```

Replace zero by any integer through five. For fixed \(k\), equating
coefficients gives \(65-3k\) equations in \(63-4k\) coefficient
variables. Nonzero leading coefficients and squarefree \(P^2+1\)
are necessary open conditions. These counts are a description of the
system, not a dimension theorem or an obstruction to exceptional solutions.

The next symbolic target is a simultaneous solution in a chart such
as (S6), using the two factorization profiles for quintic rows and
allowing mixed quadratic/quintic rows. The finite-field obstruction
forbids accepting a good-at-five degree-seven model regardless of how
many rows an algebraic solver returns. A characteristic-zero solution
still needs five independent divided divisors; repeated or dependent
rows do not qualify. It then needs a common integral specialization
slice with all-prime virtuality, independence modulo the unit line,
local separation at three, and infinitely many distinct fields.
Neither the older degree-thirty field counts nor CP1's quartic count
transfers automatically to this model.

## Replay and assurance

Run from the repository root in its existing Linux environment:

```sh
.venv/bin/python research/arithmetic/verify_septic_pell_construction.py
```

The deterministic input generator is (S5); no random seed, height cutoff
or incomplete search is used for SP2. The replay uses Python 3.12.3,
SymPy 1.14.0 and PARI/GP 2.15.4. Its good-representative trace, hashing
the ascending normalized coefficients and the ascending Frobenius
coefficients as compact JSON, is
`49837a583da277362346f09e77b0388112b813f55557254ec03c489d4cfabceb`.
It also checks (S7), its squarefreeness, residues and good-seven order,
a second quadratic-chart identity, all six coefficient-system sizes,
and rejection of a perturbed section and Frobenius control. This is
internal exact computational evidence with a written finite-coverage
argument. Independent replay, formal verification, external mathematical
review and historical-priority review are not claimed.
