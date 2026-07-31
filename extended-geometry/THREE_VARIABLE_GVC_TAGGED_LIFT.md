# Tagged-channel attempts to lift the two-pair witness to GVC(3)

## 1. Outcome and status

Work over a characteristic-zero field unless a finite field is displayed
explicitly.  This note asks whether the known rank-five two-pair
Image-Mathieu counterexample can be converted into a separated
three-variable GVC counterexample

\[
 A(\partial_t,\partial_z,\partial_y)^m(P^m)=0
 \quad(m\geq1),
 \qquad
 A(\partial)^m(QP^m)\ne0
 \quad\text{for infinitely many }m.
 \tag{1.1}
\]

No GVC(3) counterexample is obtained.  The calculations do produce four
concrete conclusions.

1. The two-pair witness already has a **coordinate-only** detecting
   multiplier.  Thus the GVC restriction on \(Q\) is not the obstruction.
2. The canonical five-channel auxiliary-exponent lift has empty normalized
   full-support fiber through pure moment four over \(\mathbf F_{101}\).
3. The minimal degree-tagged rank-two lift is impossible over
   characteristic zero: for its complete binary cubic operator jet, the
   first five pure moments generate the unit ideal.
4. Allowing the complete factor-compatible quadratic profile leaves a
   one-dimensional fiber through moment five over \(\mathbf F_{101}\), but
   moment six makes that normalized fiber empty.  In 200 deterministic
   exact fibers of the general normalized cubic profile, 198 die at moment
   five and the other two die at moment six.

Items 1 and 3 are characteristic-zero statements.  Items 2 and 4 are exact
bounded computations over the displayed finite field, not
characteristic-zero theorems.  In particular, they do not prove GVC(3).

## 2. The SIC witness has a coordinate-only detector

Use the notation of
[`TWO_PAIR_IMAGE_MATHIEU_COUNTEREXAMPLE.md`](TWO_PAIR_IMAGE_MATHIEU_COUNTEREXAMPLE.md):

\[
 F=(R+Z)\left(R^2W-\frac12(2R+Z)T^2\right).
 \tag{2.1}
\]

The proved mixed identity there is

\[
 \mathcal E_2(\xi _1z_2F^m)
 =\frac{(4m+2)!\,m!}{(2m+1)!!}.
 \tag{2.2}
\]

For every polynomial \(H\),

\[
 \mathcal E_2(\xi_iH)
 =\partial_{z_i}\mathcal E_2(H).
 \tag{2.3}
\]

Because \(z_2F^m\) has coordinate degree one greater than its dual degree,
its contraction is linear:

\[
 \mathcal E_2(z_2F^m)=A_mz_1+B_mz_2.
 \tag{2.4}
\]

Equation (2.3) and (2.2) give

\[
 A_m=\frac{(4m+2)!\,m!}{(2m+1)!!}.
 \tag{2.5}
\]

For the other coefficient, use the Hopf coordinates in the counterexample
proof.  On the unit sphere,

\[
 \xi_2z_2=|U_2|^2=\frac{1-t}{2}.
 \]

After phase averaging, \(p^m\) is even in \(t\).  Its zeroth height moment
vanishes, while its product with \(t\) integrates to zero by oddness.
Consequently

\[
 \mathcal E_2(\xi_2z_2F^m)=0,
\]

so \(B_m=0\).  We obtain the all-order identity

\[
 \boxed{
 \mathcal E_2(z_2F^m)
 =\frac{(4m+2)!\,m!}{(2m+1)!!}\,z_1\ne0.
 }
 \tag{2.6}
\]

Thus the rank-five witness fails even with a multiplier depending only on
the coordinate variables.  Its failure to be a GVC point is entirely the
failure of separability.

## 3. The direct five-channel ancilla

Let \(M\) be the \(5\times5\) coefficient matrix of \(F\).  Decompose \(F\)
by rows.  In the binary monomial bases, put

\[
\begin{aligned}
 A(\xi_1,\xi_2,u)
   &=\sum_{i=0}^4 a_i\,\xi_1^i\xi_2^{4-i}u^{4-i},\\
 P(z_1,z_2,v)
   &=\sum_{k=0}^4 b_k
       \left(\sum_{j=0}^4M_{kj}z_1^jz_2^{4-j}\right)v^{4-k}.
\end{aligned}
\tag{3.1}
\]

Contracting the auxiliary pair can raise binary tensor rank from one to
five, so the rank obstruction alone does not exclude (3.1).  Powers create
the difficulty: equal total auxiliary exponents couple different channel
multisets.

The checker normalizes the two global scalings by \(a_0=b_0=1\), keeps the
other eight channel scalars free, and forms every coefficient of

\[
 \mathcal E_3\big((AP)^m\big),
 \qquad 1\leq m\leq4.
 \tag{3.2}
\]

Over \(\mathbf F_{101}\), their Gröbner basis is \((1)\).  Hence this
canonical full-support affine chart has no point through moment four in
that characteristic.  This is a bounded modular obstruction to the most
literal rank-five purification.  It does not classify other exponent
assignments, other tensor factorizations, or the characteristic-zero
fiber.

## 4. Degree tagging gives exact diagonal channel matching

The smaller three-pair counterexample

\[
 \tau(t-y)(wz+vt)
 \tag{4.1}
\]

has tensor rank two.  Instead of adding an auxiliary pair, tag its two
channels by adjacent differential and polynomial degrees.  The natural
complete family is

\[
\boxed{
\Lambda=\partial_t\partial_z+B(\partial_t,\partial_y),
\qquad
P=zL(t,y)+C(t,y),
}
\tag{4.2}
\]

where \(L\) is linear and \(B,C\) are homogeneous binary cubics.

Expand both \(m\)-th powers.  Suppose \(j\) copies of \(B\) are chosen from
\(\Lambda^m\), and \(k\) copies of \(C\) are chosen from \(P^m\).
The \(z\)-derivatives require \(k\leq j\).  The remaining binary operator
has degree \(m+2j\), while its binary polynomial has degree \(m+2k\), so
survival also requires \(j\leq k\).  Therefore

\[
 j=k.
 \tag{4.3}
\]

This is the desired power-compatible diagonalization.  It gives the exact
all-order formula

\[
\boxed{
\Lambda^m(P^m)
=\sum_{k=0}^m
 \binom{m}{k}^2(m-k)!\,
 \partial_t^{\,m-k}
 B(\partial_t,\partial_y)^k
 \left(L^{m-k}C^k\right).
}
\tag{4.4}
\]

Every summand is a scalar.  Formula (4.4) is the useful new reduction: a
three-variable nonhomogeneous GVC search becomes a factorially weighted
binary apolar problem.

## 5. Exact obstruction for the minimal Long tag

Take the profile obtained directly by raising the second channel of (4.1)
by one degree:

\[
 L=t-y,\qquad C=t^2(t-y).
 \tag{5.1}
\]

Write the complete binary cubic operator as

\[
 B
 =a_0\partial_t^3+a_1\partial_t^2\partial_y
  +a_2\partial_t\partial_y^2+a_3\partial_y^3.
 \tag{5.2}
\]

The first pure moment is

\[
 6a_0-2a_1+1.
 \tag{5.3}
\]

Substituting (5.1)--(5.2) into (4.4) gives one exact polynomial
\(\mu_m(a_0,a_1,a_2,a_3)\) for every \(m\).  Exact rational Gröbner bases
for the successive ideals have sizes

\[
\begin{array}{c|ccccc}
N&1&2&3&4&5\\ \hline
\#\operatorname{GB}(\mu_1,\ldots,\mu_N)&1&2&4&7&1.
\end{array}
\tag{5.4}
\]

At \(N=5\), the single basis element is \(1\).  Hence:

> **Theorem 5.1 — minimal tagged-lift obstruction.**
> For (4.2), (5.1), and an arbitrary binary cubic \(B\), the five
> equations
> \[
> \Lambda^m(P^m)=0,\qquad1\leq m\leq5,
> \]
> have no solution over any characteristic-zero field.

This closes the complete operator jet for the literal degree-tagged lift,
not merely a monomial or integer-coefficient subfamily.

## 6. Wider exact modular calculations

### 6.1 Complete factor-compatible quadratic profile

Keep \(L=t-y\) and take

\[
 C=(t-y)(t^2+q_1ty+q_2y^2),
 \tag{6.1}
\]

with the \(t^2\)-coefficient normalized to one.  Over
\(\mathbf F_{101}\), using the complete four-parameter binary cubic \(B\),
the first five moments have a Gröbner basis of size \(158\) and dimension
one.  Adding moment six gives the unit ideal.

This is exact over \(\mathbf F_{101}\).  An empty affine fiber at one prime
is not by itself a characteristic-zero certificate, so no theorem over
\(\mathbb Q\) is asserted.

### 6.2 General normalized cubic profiles

Take

\[
 C=t^3+c_1t^2y+c_2ty^2+c_3y^3.
 \tag{6.2}
\]

For each fixed triple \((c_1,c_2,c_3)\in\mathbf F_{101}^3\), formula
(4.4) gives equations only in the four coefficients of \(B\).  The checker
solves 200 such fibers exactly.  The sample contains six fixed boundary
profiles and 194 profiles generated with seed `20260730`.

\[
\begin{array}{c|cc}
\text{first empty moment}&5&6\\ \hline
\text{number of fibers}&198&2.
\end{array}
\tag{6.3}
\]

No sampled fiber survives through moment six.  The two delayed profiles are

\[
 (96,38,76),\qquad(47,49,36)
 \quad\text{in }\mathbf F_{101}^3.
\tag{6.4}
\]

This is a bounded fiber search, not an exhaustive calculation on the
seven-parameter total space.

## 7. The three-channel weighted-quartic continuation

The minimal positive grading containing ordinary degrees two, three, and
four is

\[
 2\nu_z+\nu_t+\nu_y=4.
 \tag{7.1}
\]

Its nine monomials split into \(z^2\), the three terms
\(z\{t^2,ty,y^2\}\), and the five binary quartics.  Thus

\[
\begin{aligned}
 \Lambda&=a\partial_z^2+
   \partial_zB_2(\partial_t,\partial_y)+D_4(\partial_t,\partial_y),\\
 P&=\alpha z^2+zC_2(t,y)+E_4(t,y).
\end{aligned}
\tag{7.2}
\]

Positivity of the grading makes every surviving pure contraction a
scalar.  Unlike (4.4), matching weighted totals does not match channel
counts: two cubic selections can balance one quadratic and one quartic
selection.  This is the desired coupled convolution.

### 7.1 Exact Dvorsky-lattice compression

The four Dvorsky--Long exponent vectors have linear span three.  On the
weighted plane, enumerate rank-three parallelograms
\(q_1+q_4=q_2+q_3\) and normalize

\[
 A=X^{q_1}+X^{q_4},\qquad
 P=x^{q_1}+d x^{q_2}+e x^{q_3}
   -\frac{q_1!}{q_4!}x^{q_4}.
\tag{7.3}
\]

There are \(28\) unoriented parallelograms.  Exact rational Gröbner bases
for their \(56\) orientations make \(54\) unit at moment three and two at
moment four.  The delayed pair is one orbit:

\[
\Lambda=\partial_z(\partial_t^2+\partial_y^2),\qquad
P=z(y^2-t^2)+d\,t^2y^2+e\,z^2.
\tag{7.4}
\]

Here

\[
 \mu_1=\mu_3=0,\quad \mu_2=32(de+2),\quad
 \mu_4=165888(3d^2e^2+8de+8).
\tag{7.5}
\]

Thus \(\mu_2=0\) forces \(de=-2\), after which the last parenthesis is
\(4\).  The rank-three lattice compresses, but its factorial weight does
not.

### 7.2 The first persistent branch is terminal

Adding the middle point gives

\[
 P=z(y^2-t^2+h\,ty)+d\,t^2y^2+e\,z^2.
\tag{7.6}
\]

The exact moment ideal through order ten has reduced basis

\[
 (de,\ h^2+4).
\tag{7.7}
\]

This is also an all-order pure-zero locus.  When \(h^2=-4\), the cubic
bracket is an isotropic square \(L^2\).  If \(e=0\), quartic selections
have too little \(z\)-degree; if \(d=0\), \(z^2\)-selections have too
little binary degree.  A fixed multiplier bridges only bounded tag depth,
so in either branch the required transverse derivative degree eventually
exceeds its fixed supply.  The equation \(de=0\) is exactly what prevents
the defect from migrating linearly.

### 7.3 Complete quartic repair on the polynomial side

In isotropic coordinates take

\[
\Lambda=\partial_z\partial_L\partial_M,\qquad
P=zL^2+z^2+\sum_{j=0}^4a_jL^{4-j}M^j.
\tag{7.8}
\]

Through moment seven the exact basis is

\[
 \bigl(a_4(a_0+1),a_1a_4,a_2,a_3\bigr).
\tag{7.9}
\]

Moment eight adds \(a_4^2\), giving radical

\[
 (a_2,a_3,a_4).
\tag{7.10}
\]

On this radical, each quartic selection supplies at most one \(M\), while
tag balance allows only about \(m/2\) quartic selections.
\(\partial_M^m\) therefore gives an all-order mixed cutoff.

### 7.4 Activating operator endpoints

The smallest simultaneous repair tested is

\[
\begin{aligned}
 \Lambda&=\partial_z\partial_L\partial_M
   +A\partial_z^2+B\partial_L^4+C\partial_M^4,\\
 P&=zL^2+z^2+FL^4+GM^4.
\end{aligned}
\tag{7.11}
\]

The exact moment ideal through order ten has basis

\[
 AB^2,\ A(3A+B),\ A+12BF+12CG,\ AG,\ BG,\ G(F+1),\ G^2.
\tag{7.12}
\]

Its radical is \((A,G,BF)\).  If \(F\ne0\), every operator factor has
positive \(M\)-degree; if \(F=0\), the simultaneous \(z\)- and \(M\)
inequalities become incompatible at linear depth.  Hence these branches
are terminal for every fixed multiplier.

### 7.5 Complete odd-quartic jet

The remaining parity chart is

\[
\begin{aligned}
 \Lambda={}&\partial_z\partial_L\partial_M+A\partial_z^2
   +U\partial_L^3\partial_M+V\partial_L\partial_M^3,\\
 P={}&zL^2+z^2+RL^3M+SLM^3.
\end{aligned}
\tag{7.13}
\]

Thus the degree-two endpoint is retained on both sides, and all odd binary
quartics occur on both sides.  The exact moment ideal through order six
has reduced Gröbner basis

\[
\begin{split}
 A^2,\quad A+3RU+3SV,\quad AS,\quad S(10U+1),\\
 RS,\quad S^2.
\end{split}
\tag{7.14}
\]

Its radical is

\[
 \sqrt{I_6}=(A,S,RU)
  =(A,S,R)\cap(A,S,U).
\tag{7.15}
\]

Indeed, \(A^2,S^2\in I_6\), and the linear relation in (7.14) then puts
\(RU\) in the radical.  Conversely every element of (7.14) lies in
\((A,S,RU)\), which is the displayed intersection of two prime ideals.
This also shows what happens to the degree-two operator channel: the
moment scheme may retain it nilpotently, but every field-valued pure-zero
point has \(A=0\).

Both components are terminal to all orders.  On \(R=0\), \(P\) is
\(M\)-free while every monomial of \(\Lambda\) has positive
\(\partial_M\)-degree.  Hence for a fixed multiplier \(Q\),
\(\Lambda^m(QP^m)=0\) once \(m>\deg_M Q\).

On \(U=0\), choose \(k\) copies of
\(\partial_L\partial_M^3\) in a term of \(\Lambda^m\), and let \(c\) be
the number of \(RL^3M\) selections in \(P^m\).  For a monomial of \(Q\)
with exponents \((q_z,q_L,q_M)\), the available \(M\)-degree must satisfy

\[
 c+q_M\geq m+2k.
\tag{7.16}
\]

Since \(c\leq m\), this forces \(2k\leq q_M\), while the number of
remaining polynomial selections satisfies
\(m-c\leq q_M-2k\).  Those selections are \(zL^2\) or \(z^2\), so their
\(z\)-degree is at most \(2(q_M-2k)\).  Supplying the operator's
\(z\)-degree \(m-k\) therefore requires

\[
 m\leq 2q_M+q_z-3k\leq2q_M+q_z.
\tag{7.17}
\]

Taking the maximum over the finitely many monomials of \(Q\) proves the
fixed-depth mixed cutoff.  In particular, (7.13) contains no GVC(3)
counterexample.

These are exact characteristic-zero calculations.  They close every
four-term Dvorsky compression, the three minimal repairs above, and the
complete odd-quartic operator/polynomial jet.

## 8. What remains

The direct lift has therefore not broken GVC(3).  Formula (4.4) identifies
the remaining targets much more sharply:

1. compute the characteristic-zero radical of the general cubic-profile
   moment ideal, rather than sampling its finite-field fibers;
2. compute the simultaneous complete even-and-odd quartic radical in
   (7.2); the separate even endpoint repairs and the complete odd chart
   are terminal, but their mixed total space has not been eliminated;
3. allow \(B\) to contain further \(z\)-derivative layers, which introduces genuine
   positive-depth equations instead of the exact diagonal rule (4.3);
4. seek a different purification of the rank-five tensor rather than its
   row decomposition.

The most economical next attack is item 2, using the two prime components
in (7.15) as branch constraints before restoring the even quartic
endpoints.  The odd-only torus chart itself is now closed.

## 9. Reproduction

Run

```bash
.venv/bin/python scripts/research_three_variable_gvc_tagged_lift.py
.venv/bin/python scripts/research_three_channel_gvc_lift.py
```

The scripts write
[`three_variable_gvc_tagged_lift.json`](../artifacts/generated-results/three_variable_gvc_tagged_lift.json)
and
[`three_channel_gvc_lift.json`](../artifacts/generated-results/three_channel_gvc_lift.json).
They require SymPy and Singular.  The all-order proofs are Sections 2, 4,
and 7; the scripts replay their bounded consequences and perform the exact
Gröbner calculations described above.
