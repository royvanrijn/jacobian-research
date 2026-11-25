# Primitive two-pair Image witnesses in every odd composite balanced degree

Let K be a field of characteristic zero. Use the two contraction pairs
\((\xi_1,z_1),(\xi_2,z_2)\), the contraction functional
\(\mathcal E_2(\xi^\alpha z^\beta)=\partial_z^\alpha z^\beta\), and

\[
 R=\xi_1z_1+\xi_2z_2,\quad Z=\xi_1z_2,\quad
 W=2\xi_2z_1,\quad T=\xi_1z_1-\xi_2z_2.
 \tag{O1}
\]

Thus \(T^2=R^2-2ZW\). Balanced degree d means bidegree \((d,d)\)
in \((\xi,z)\), or ordinary total degree \(2d\). Call a polynomial
R-primitive if it is not divisible by R.

**Theorem OddCompositeImage.** Every odd composite integer d admits
an explicit R-primitive, non-proper-power polynomial \(F_d\) of
balanced degree d, and a fixed polynomial \(Q_d\), such that

\[
 \mathcal E_2(F_d^m)=0\quad(m\ge1),\qquad
 \mathcal E_2(Q_dF_d^m)\ne0\quad(m\ge2\text{ even}).
 \tag{O2}
\]

The marker is fixed as m varies. Therefore \(F_d\) supplies failure
of eventual mixed vanishing for the two-pair Special Image space.
The construction and all-order proof are rational and work over every
characteristic-zero field. They do not decide primitive witnesses in
odd prime balanced degrees, including five and seven. In particular,
nine is the smallest odd degree supplied here, not a proved minimum
odd failing degree.

The new step is a nonlinear dilation of the signed height in the
endpoint argument. It differs from the single affine height profile
classified in [HOPF_LIFT_CLASSIFICATION.md, Section 7](HOPF_LIFT_CLASSIFICATION.md#7-real-profiles-with-mixed-height-powers).
The proof below is an internal mathematical deduction, supported by
exact finite controls. External human review, formal verification and
historical priority are not asserted.

## 1. The explicit degree-nine witness

The first case is

\[
\begin{split}
 F_9={}&T(R+Z)^3\bigl(2R^4W
             -T^2(4R^3+6R^2Z+4RZ^2+Z^3)\bigr),\\
 Q_9={}&Z^2.
\end{split}                                                     \tag{O3}
\]

For every m>=1 its pure contraction is zero. Its marked contraction is

\[
 \mathcal E_2(Z^2F_9^m)=
 \begin{cases}
 (9m+3)!\displaystyle\frac{2^m m!}
                  {\prod_{j=0}^{m}(m+2j+1)},&m\text{ even},\\
 0,&m\text{ odd}.
 \end{cases}                                                     \tag{O4}
\]

For example the normalized marked angular value at m=2 is 8/105.
The nonzero values at all even powers suffice for eventual failure;
nonvanishing at odd powers is neither needed nor claimed.

At \((\xi_1,\xi_2,z_1,z_2)=(1,1,1,-1)\), one has R=0 and
\(F_9=-8\), proving R-primitivity directly. The factor T has exact
multiplicity one: after dividing (O3) by T, its value at (1,1,1,1)
is nonzero. Since T is an absolutely irreducible bilinear form, this
also proves that \(F_9\) is not a proper power, even after a field
extension and multiplication by a nonzero scalar.

## 2. The original angular and radial functionals

Work first over Q and then extend scalars. On the affine quadric

\[
 t^2+2xy=1,
\]

define \(\Lambda\) by substituting \(y=(1-t^2)/(2x)\), taking the
constant term in the Laurent variable x, and integrating in t against
dt/2 on [-1,1]. This definition is independent of the polynomial
representative modulo the quadric relation. All operations on each
polynomial are finite and give rational numbers.

For a homogeneous polynomial H(R,Z,W,T) of balanced degree D,
the retained circular-Gaussian/Hopf identity is

\[
 \mathcal E_2(H)=(D+1)!\,
       \Lambda\bigl(H(1,x,y,t)\bigr).                         \tag{O5}
\]

Its exact definitions and proof are in
[TWO_PAIR_IMAGE_MATHIEU_COUNTEREXAMPLE.md, Section 2](TWO_PAIR_IMAGE_MATHIEU_COUNTEREXAMPLE.md#2-contraction-as-circular-gaussian-expectation).
For clarity, the circular Gaussian radius R has Gamma shape two, so
its D-th moment is (D+1)!. The angular variables are x=Z/R,
y=W/R and t=T/R. This explains the radial factor in (O5) without
replacing the two-pair contraction by a ternary Gaussian functional.

We also use the classical identity

\[
 \mathcal M_2=
  (\partial_{z_1}-\xi_1)K[\xi,z]
 +(\partial_{z_2}-\xi_2)K[\xi,z]
       =\ker\mathcal E_2.                                  \tag{O6}
\]

It has the direct algebraic proof retained in that same source:
the functional annihilates both displayed operator images; conversely
replace \(\xi_i A\) by \(\partial_{z_i}A\) modulo their sum and
induct on the total xi-degree. Thus \(H-\mathcal E_2(H)\) belongs
to the image, proving the reverse inclusion. No unproved Image
conjecture or invariant-theory nullcone assertion is an input.

## 3. The nonlinear height-dilation identity

Fix integers r>=1 and s>=0, and put

\[
 a=2r+1,\quad b=2r+1+2s,\quad d=ab,\quad A=1+x,
\]

\[
 S_s(v)=\begin{cases}1,&s=0,\\1+v^s,&s>0,\end{cases}
 \qquad H_{r,s}(u)=u(1-u^2)^rS_s(u^2).
 \tag{O7}
\]

In particular H is odd, has degree b and leading coefficient
\((-1)^r\), has roots of order at least r at both 1 and -1,
and has a simple root at zero. Define the Laurent presentation

\[
 p_{r,s}=x^{-r}A^rH_{r,s}(tA^{2r}),
 \qquad q_{r,s}=x^{2r}.
 \tag{O8}
\]

It is a polynomial on the quadric. Indeed, put

\[
 B_r=\frac{1-t^2A^{4r}}x
       =2y-t^2\sum_{j=1}^{4r}\binom{4r}{j}x^{j-1}.
\]

Then the literal polynomial representative is

\[
 p_{r,s}=tA^{3r}B_r^rS_s(t^2A^{4r}).                       \tag{O9}
\]

Let

\[
 C_m=\frac12\int_{-1}^{1}H_{r,s}(u)^m\,du.
 \tag{O10}
\]

**Lemma NonlinearHeightMoments.** For every m>=1,

\[
 \Lambda(p_{r,s}^m)=0,\qquad
 \Lambda(x^{2r}p_{r,s}^m)=C_m.
 \tag{O11}
\]

Furthermore C_m=0 for odd m and C_m>0 for even m. These numbers are
rational, so the nonzero assertion also holds after scalar extension
to any characteristic-zero field.

**Proof.** For odd m, the polynomial representative (O9) is odd in t;
its other factors and x^(2r) are even in t. Both angular means and
C_m are therefore zero.

For even m>=2 define the polynomial primitive

\[
 J_m(v)=\int_0^v H_{r,s}(u)^m\,du.
\]

Changing variable u=tA^(2r) in the symmetric polynomial integral gives

\[
 A^{rm}\frac12\int_{-1}^{1}H_{r,s}(tA^{2r})^m\,dt
          =A^{r(m-2)}J_m(A^{2r}).                           \tag{O12}
\]

This is a polynomial identity: rm-2r>=0, and it can be checked
term by term without assigning a complex integration path. Since
J'_m=H^m has a zero of order at least rm at 1,

\[
 J_m(A^{2r})-J_m(1)\quad\hbox{is divisible by}\quad x^{rm+1}.
\]

Here J_m(1)=C_m, because H^m is even. Therefore

\[
 A^{r(m-2)}J_m(A^{2r})
       \equiv C_mA^{r(m-2)}\pmod{x^{rm+1}}.                 \tag{O13}
\]

The pure angular moment extracts its coefficient of x^(rm),
which is zero because r(m-2)<rm. The marked moment extracts
the coefficient of x^(rm-2r)=x^(r(m-2)), which is exactly C_m.
For even m, H is a nonzero real polynomial, so its even power has
a strictly positive integral on [-1,1]. This proves (O11) in every
order; no finite prefix is being extrapolated. QED.

The important exponent in (O12) is r(m-2). The nonlinear dilation
A^(2r) removes two powers per winding r. Parity disposes of m=1,
while the remaining exponent is polynomial for every even m>=2.

## 4. Explicit homogeneous lifts and all odd composite degrees

Define the homogeneous polynomial of degree 4r+1 in the bilinears

\[
 \mathcal B_r=2R^{4r}W
       -T^2\sum_{j=1}^{4r}\binom{4r}{j}R^{4r-j}Z^{j-1}.
 \tag{O14}
\]

Set

\[
 F_{r,s}=T(R+Z)^{3r}\mathcal B_r^r\,
 \begin{cases}
 1,&s=0,\\
 R^{2as}+T^{2s}(R+Z)^{4rs},&s>0,
 \end{cases}
 \qquad Q_{r,s}=Z^{2r}.
 \tag{O15}
\]

At s=0 the final factor is one, not a duplicated constant two.
For s>0 each summand in that final factor has degree 2as.
Thus F has balanced degree

\[
 1+3r+r(4r+1)+2as=a^2+2as=a(a+2s)=ab=d.
\]

On setting R=1, (O15) becomes exactly (O9). By (O5) and (O11),

\[
 \mathcal E_2(F_{r,s}^m)=0,
 \qquad
 \mathcal E_2(Z^{2r}F_{r,s}^m)=(dm+2r+1)!\,C_m.
 \tag{O16}
\]

Every odd composite d can be written d=ab with odd 3<=a<=b;
take r=(a-1)/2 and s=(b-a)/2. These are admissible integers,
so (O16) supplies all the degrees claimed in the theorem.

To prove R-primitivity, reduce the explicit homogeneous expression
modulo R. Its top part is

\[
 F_{r,s}\bmod R=(-1)^rT^bZ^{2rb}.
 \tag{O17}
\]

At the point (1,1,1,-1) this is \((-1)^r2^b\ne0\), while R=0.
Hence R does not divide F in the original four-variable polynomial
ring, not merely in a chosen Laurent presentation.

Finally, T has exact multiplicity one in (O15). At (1,1,1,1),
T=0, R=2, Z=1, W=2 and R+Z=3. The quotient F/T has value

\[
 3^{3r}2^{r(4r+2)}\quad(s=0),\qquad
 3^{3r}2^{r(4r+2)+2as}\quad(s>0),
\]

which is nonzero. The coefficient matrix of the bilinear form T
is diag(1,-1), of rank two. Any nontrivial factorization of a
form of bidegree (1,1) would have rank one, so T is absolutely
irreducible. Its multiplicity one excludes a proper power of F,
even up to a nonzero scalar over any field extension.

For r=1,s=0, H(u)=u(1-u^2) and F is (O3). For even m,
polynomial beta integration gives

\[
 C_m=\int_0^1u^m(1-u^2)^m\,du
       =\frac{2^m m!}{\prod_{j=0}^m(m+2j+1)}.
\]

One can derive the last expression directly by integrating the
derivative of u^(m+1)(1-u^2)^m and iterating, as in the retained
Hopf-profile proof. This proves (O4).

**Corollary CompositeImage.** Every composite integer d>=4 admits an
R-primitive, non-proper-power two-pair Image witness with one fixed marker
whose mixed contractions fail eventual vanishing. For even d use
Theorem 5.5 of the retained quartic note; for odd d use
Theorem OddCompositeImage above. This combines the two existence families
and does not identify their formulas or assert a prime-degree obstruction.

## 5. Antecedents, verification and the remaining prime-degree question

The [source ledger](odd_primitive_two_pair_image_source_ledger.json)
records the exact inputs, inspected primary statements and assurance
limits. The all-order mechanism uses polynomial integration, endpoint
Taylor divisibility, the original radial normalization and the
image-kernel identity. It requires no new imported classification.

The primary classical image-kernel statement is van den Essen--Wright--Zhao,
[*On the Image Conjecture*, Theorem 3.3](https://arxiv.org/html/1008.3962v2#S3.Thm3),
which attributes it to Zhao's Theorem 3.1; its actual statement was
inspected, and the elementary proof needed here is given above.
Euler's beta integral is
[NIST DLMF 5.12.1](https://dlmf.nist.gov/5.12.E1), whose positive-real
parameter hypotheses hold in the displayed specialization.

Long's [Theorem 2.1 and Theorem 4.2](https://arxiv.org/html/2607.19012v1)
already give mixed torus/interval and SU(2) Mathieu counterexamples by
beta-binomial cancellation. Those exact formulas and their proofs were
inspected. They are antecedents, not a claim that unrestricted failure
or this cancellation principle is new. The theorem here concerns the
additional specified R-primitivity and odd composite balanced degrees
in the two-pair contraction space; worldwide priority for that refinement
has not been established.

The [standard-library exact checker](verify_odd_primitive_two_pair_image.py)
constructs the homogeneous polynomials in the original four-variable
monomial basis and evaluates their actual contractions. It separately
checks finite angular endpoint identities, the radial factors,
R=0 evaluation, exact T multiplicity and rejection of altered contraction
claims. The [receipt](odd_primitive_two_pair_image_verification.json)
states the finite domains actually checked. These controls support the
explicit formulas; all integers r,s and all moment orders use the
written argument (O7)--(O17), not bounded experimental coverage.

The affine-multiplier quintic obstruction and real affine-profile
classification retain their exact scopes. This nonlinear family does
not produce balanced degree five or seven. Together with the retained
even-degree family, it supplies every composite balanced degree at least
four. The separate [hyperbolic inverse-cover construction](HYPERBOLIC_ODD_HEIGHT_IMAGE_WITNESSES.md)
also supplies certain odd prime degrees, including thirteen and twenty-three.
Primitive degree five remains open. The separate [joint-Hopf construction](PRIMITIVE_IMAGE_DEGREE_SEVEN_AND_HOPF_SUBALGEBRA.md) realizes degree seven and every balanced degree at least six.
No exhaustive pure-phase classification or formal verification is claimed.
