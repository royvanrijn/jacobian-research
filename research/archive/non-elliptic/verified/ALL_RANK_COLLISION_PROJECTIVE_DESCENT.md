# All-rank collision frames and projective Tschirnhaus descent

Let `A` be a rank-`N` finite-etale algebra over a characteristic-zero field,
with `N>=3`, and let `r,u in A` be two primitive coordinates.  Put

\[
 W_r=\langle 1,r\rangle\subset A.
\]

The exact projective-presentation criterion is

\[
\boxed{
 u=\frac{ar+b}{cr+d}\text{ for some }[a:b:c:d]\in PGL_2
 \quad\Longleftrightarrow\quad
 W_r\cap uW_r\ne0.}                                  \tag{0.1}
\]

Equivalently, the four vectors

\[
 1,\quad r,\quad u,\quad ru                           \tag{0.2}
\]

have rank at most three in the `N`-dimensional vector space `A`.  Thus the
projective locus is cut out without choosing or ordering roots by the
`4`-by-`4` minors of their coefficient matrix.

On the primitive-presentation open this locus is smooth of codimension

\[
\boxed{N-3.}                                          \tag{0.3}
\]

After passing to the full collision frame, the equations are the `N-3`
residuals obtained by matching the first three roots and testing each
remaining root.  This recovers:

- no equation in rank three;
- the single cross-ratio equation in rank four; and
- `N-3` genuine projective-moduli obstructions in every higher rank.

This theorem identifies exactly where canonical `PGL_2` root transport is
available.  It does **not** say that every isomorphism of Keller incidences
must act projectively on the root coordinate.

## 1. Credit and scope

The tensor collision algebra, diagonal kernel, and off-diagonal sheet used
here are credited to Chloe van der Vlugt's *Collision Ideals and
Off-Diagonal Sheets*, with the public source attribution and limitations
recorded in the
[external audit](COLLISION_IDEALS_EXTERNAL_AUDIT.md).  The higher
configuration tower, all-rank determinantal criterion, codimension theorem,
and universal Keller-coordinate formula below are deductions made in this
repository.

The proof is stated over a characteristic-zero field to match the Keller
applications.  The linear-algebra argument itself only requires a separable
degree-`N` presentation and the displayed primitive opens.

## 2. The last collision configuration is the full frame

Let `E -> S` be finite etale of rank `N`.  Define

\[
 \operatorname{Conf}_m(E/S)
 =
 E^m\setminus\bigcup_{i<j}\Delta_{ij}.                \tag{2.1}
\]

It is finite etale of rank `N!/(N-m)!`.  Over
`\operatorname{Conf}_{N-1}(E/S)`, the `N-1` disjoint marked sections have a
rank-one open-and-closed complement.  A rank-one finite-etale cover is the
base itself, so the complementary section is canonical.  Appending it gives

\[
\boxed{
 \operatorname{Conf}_{N-1}(E/S)
 \simeq
 \operatorname{Isom}_S(\{1,\ldots,N\}_S,E).}          \tag{2.2}
\]

Hence `Conf_(N-1)` has rank `N!` and is the full `S_N` frame torsor.  This
generalizes `Off_2=Conf_2` in rank three and `Conf_3` in rank four.

The frame removes the finite permutation ambiguity.  It does not identify
two embeddings of the framed `N`-point set in the root line; that remaining
problem is measured by (0.1).

## 3. Intrinsic Schubert-incidence criterion

Fix a field `k`, a rank-`N` finite-etale `k`-algebra `A`, and primitive
elements `r,u`.  After a separable closure,

\[
 A_{\bar k}\simeq\bar k^N,\qquad
 r=(r_1,\ldots,r_N),\qquad
 u=(u_1,\ldots,u_N),                                 \tag{3.1}
\]

where both lists have pairwise distinct entries.

Because `r` is primitive, `1,r` are linearly independent.  Multiplication
by `u` is injective on `W_r`: if
`u(c r+d)=0`, then the nonzero linear polynomial `cT+d` would have to vanish
at all but at most one of the `N` distinct `r_i`, which is impossible for
`N>=3`.  Therefore `W_r` and `uW_r` are both two-planes in `A`.

They meet nontrivially exactly when there is a nonzero relation

\[
 ar+b=u(cr+d).                                       \tag{3.2}
\]

The determinant `ad-bc` cannot vanish.  If it did, the numerator and
denominator in (3.2) would be proportional.  Away from the at most one root
of `cT+d`, at least `N-1>=2` of the values `u_i` would then be equal,
contradicting primitivity of `u`.

The denominator is also a unit in `A`.  If `cr_i+d=0` for one geometric
root, (3.2) also gives `ar_i+b=0`; the two linear polynomials would share a
root and `ad-bc` would vanish.  Thus (3.2) is precisely

\[
 u=(ar+b)(cr+d)^{-1}
\]

for an element of `PGL_2(k)`.  This proves (0.1).

Now choose any basis of `A` and let `C(r,u)` be the `N`-by-`4` matrix whose
columns are (0.2).  Equation (3.2) is equivalent to

\[
\boxed{
 \operatorname{rank}C(r,u)\le3
 \quad\Longleftrightarrow\quad
 I_4(C(r,u))=0.}                                     \tag{3.3}
\]

Changing the basis of `A` multiplies `C` on the left by an invertible
matrix, so its determinantal ideal is intrinsic.  This is the Schubert
incidence condition that the two two-planes `W_r` and `uW_r` intersect.

## 4. Framed residuals and exact codimension

On the frame torsor, evaluate (0.2) at the labeled roots:

\[
 {\cal E}(r,u)=
 \begin{pmatrix}
 1&r_1&u_1&r_1u_1\\
 1&r_2&u_2&r_2u_2\\
 \vdots&\vdots&\vdots&\vdots\\
 1&r_N&u_N&r_Nu_N
 \end{pmatrix}.                                      \tag{4.1}
\]

The evaluation map from `A` to the split algebra is invertible, so (3.3) is
equivalent to `rank(E)<=3`.

The first three source and target values determine a unique projective
matrix.  With

\[
 M_{123}=
 \begin{pmatrix}
 r_1&1&-u_1r_1&-u_1\\
 r_2&1&-u_2r_2&-u_2\\
 r_3&1&-u_3r_3&-u_3
 \end{pmatrix},                                      \tag{4.2}
\]

let `(a,b,c,d)` be its signed maximal minors.  Then

\[
 \det\begin{pmatrix}a&b\\c&d\end{pmatrix}
 =
 V(r_1,r_2,r_3)V(u_1,u_2,u_3),                       \tag{4.3}
\]

which is a unit on the primitive framed open.  For every `i>=4`, put

\[
 R_i=ar_i+b-u_i(cr_i+d).                              \tag{4.4}
\]

Laplace expansion identifies `R_i`, up to the fixed row-order sign, with
the `4`-by-`4` minor on rows `1,2,3,i`.  Consequently

\[
\boxed{
 I_4({\cal E})=0
 \quad\Longleftrightarrow\quad
 R_4=\cdots=R_N=0.}                                  \tag{4.5}
\]

The coefficient of `u_i` in `R_i` is `-(cr_i+d)`, a unit on the projective
primitive locus.  The Jacobian of `(R_4,\ldots,R_N)` with respect to
`(u_4,\ldots,u_N)` is therefore diagonal and invertible.  This proves that
the projective locus is smooth of codimension `N-3`, and that the equations
in (4.5) generate its ideal on this frame chart.

Equivalently, after framing, evaluation identifies the underlying affine
space of `A` with \(\mathbb A_k^N\).  The primitive coordinates form an
open subset, and the projective locus is the image of the open subset of
`PGL_2` on which no denominator vanishes.  Three points make this map
injective and give the inverse (4.2), so the locus has dimension three.

In moduli language, `(r_1,\ldots,r_N)` and `(u_1,\ldots,u_N)` define the
same point of \(M_{0,N}\) exactly on (4.5).  Its dimension `N-3` is the
number of independent defects.  This interpretation is not needed for the
determinant proof.

## 5. Polynomial coefficient matrix

The criterion can be evaluated without adjoining roots.  Normalize a
separable relation as

\[
 E(S)=a_0+S+a_2S^2+\cdots+a_NS^N,\qquad a_N\ne0,      \tag{5.1}
\]

and write

\[
 u=q(r),\qquad
 q(S)=q_0+q_1S+\cdots+q_{N-1}S^{N-1}.                \tag{5.2}
\]

In the basis `(1,r,...,r^(N-1))`, use the matrix

\[
\boxed{
 C_N(E,q)=
 \begin{pmatrix}
 1&0&q_0&-a_0q_{N-1}\\
 0&1&q_1&a_Nq_0-q_{N-1}\\
 0&0&q_2&a_Nq_1-a_2q_{N-1}\\
 \vdots&\vdots&\vdots&\vdots\\
 0&0&q_{N-1}&a_Nq_{N-2}-a_{N-1}q_{N-1}
 \end{pmatrix}.}                                     \tag{5.3}
\]

The last column is exact because

\[
 a_NS q(S)
 -
 \sum_{j=0}^{N-1}C_N(E,q)_{j,4}S^j
 =
 q_{N-1}E(S).                                        \tag{5.4}
\]

Thus its class in `k[S]/(E)` is `a_Nru`.  Since `a_N` is a unit,

\[
\boxed{
 u\text{ is projectively related to }r
 \quad\Longleftrightarrow\quad
 I_4(C_N(E,q))=0}                                    \tag{5.5}
\]

on the primitive overlap.  The equations are polynomial in the normalized
presentation and Tschirnhaus coefficients.

For the
[universal relative Keller map](UNIVERSAL_RELATIVE_KELLER_MAP.md),

\[
 E_{{\bf u},\pi,b,c}(S)
 =
 S+bS^2+\pi S^3+\sum_{j=4}^Nu_j\pi^jS^j-\frac c2.    \tag{5.6}
\]

Substitute

\[
 a_0=-c/2,\quad a_2=b,\quad a_3=\pi,\quad
 a_j=u_j\pi^j\ (j\ge4)                               \tag{5.7}
\]

in (5.3).  The ideal `I_4(C_N)` is therefore an explicit polynomial
projective-descent ideal on the actual universal Keller
parameter--target chart.  It is invariant under relabeling because (5.3)
was constructed in the quotient algebra before passing to a frame.

The equality between its codimension `N-3` and the number of seed parameters
of the universal relative Keller map is numerically exact but does not
identify those two coordinate systems.  The seed/target split in (5.6) is
not a `PGL_2` quotient construction.

## 6. Low-rank recovery and uniform witnesses

For `N=3`, the matrix `C_3` has only three rows, so its rank is automatically
at most three.  This is the projective-interpolation input used by the
[rank-three collision-framed audit](RANK_THREE_COLLISION_DESCENT.md).

For `N=4`, (5.3) is square and

\[
\boxed{
 \det C_4
 =
 a_4(q_2^2-q_1q_3)-a_3q_2q_3+a_2q_3^2.}             \tag{6.1}
\]

Writing `e_1=-a_3/a_4` and `e_2=a_2/a_4` shows that
`\det C_4=a_4\Psi`, where `Psi` is the exact defect in the
[rank-four cross-ratio audit](RANK_FOUR_COLLISION_CROSS_RATIO.md).
Under (5.7), (6.1) becomes

\[
 u_4\pi^4(q_2^2-q_1q_3)-\pi q_2q_3+bq_3^2,           \tag{6.2}
\]

with no additional factor.

For every `N>=4`, take the split roots `(1,2,\ldots,N)` and

\[
 q(r)=r+r^2.                                         \tag{6.3}
\]

The values `i+i^2` are pairwise distinct, so `q(r)` remains primitive.
Its projective defect on the first four roots is `1`, so the corresponding
evaluation minor is nonzero and `rank C_N=4`.  Equation (6.3) is therefore a
uniform exact witness that full collision framing does not make an arbitrary
primitive change projective in any rank at least four.

Conversely, the values

\[
 u_i=\frac{2r_i+1}{r_i+N+1}                          \tag{6.4}
\]

are pairwise distinct and interpolate to a polynomial of degree less than
`N`.  Its matrix `C_N` has rank exactly three.  This gives an exact
non-affine projective witness in every rank.

## 7. What the theorem buys and where to focus

The finite-etale part of presentation descent is now separated into two
complete layers:

1. `Conf_(N-1)` removes all finite `S_N` labeling ambiguity.
2. `I_4(C_N)` measures the remaining projective embedding ambiguity.

This has three concrete consequences for Keller research.

- A proposed canonical `PGL_2` transport only needs to be constructed on the
  smooth codimension-`N-3` locus `I_4(C_N)=0`; outside it, failure is
  mathematical rather than a framing artifact.
- Any descent over the whole primitive-presentation groupoid must supply
  genuinely nonprojective Tschirnhaus transport in `N-3` independent
  directions.  Collision sheets alone cannot supply those directions.
- Rank four is the minimal test case: it has one transverse direction and
  the exact witness `q=r+r^2`.  A successful or obstructed Keller lift there
  is the smallest result capable of distinguishing projective transport
  from full presentation descent.

The next focused program is therefore:

1. extend the rank-three target-local factorization transport to the
   all-rank projective locus and record every boundary denominator;
2. use the
   [rank-four nonprojective continuation](RANK_FOUR_NONPROJECTIVE_KELLER_LIFT.md),
   which isolates the ground-field Kummer class, produces a neutral
   fixed-map endpoint problem, and constructs the exact first-order and
   all-finite-order formal lift of the single normal direction; its
   four-to-two fiber-orbit drop also excludes the straight target
   translation globally;
3. decide whether a target symmetry of degree at least nineteen realizes
   the prescribed collision-frame permutation—the prime discriminant has
   degree thirteen, every lower-degree target symmetry is in `mu_5`, and
   its endpoint orbit fails, while exact logarithmic/Jacobian equations
   exclude degrees thirteen through eighteen—then test whether such
   endpoint lifts compose to a polynomial cocycle;
4. only after that rank-four gate, treat the other `N-3` residual directions
   in (4.5).

Step 2 is now complete formally but not globally.  This theorem remains a
complete projective-descent result, not a full Keller-descent theorem.

## 8. Exact regression

Run

```bash
.venv/bin/python scripts/verify_all_rank_collision_projective_descent.py
```

The checker verifies:

- the full-frame completion and stabilizers through rank eight;
- the signed-minor interpolation residuals;
- the exact coefficient identity (5.4) through rank eight;
- polynomial evaluation and rank on exact projective witnesses;
- the uniform primitive nonprojective witness (6.3);
- the quartic specialization (6.1)--(6.2); and
- rank `N-3` of the framed residual Jacobian through rank ten.

All calculations use exact symbolic or rational arithmetic.  The bounded
rank ranges replay the uniform formulas; the all-rank claims themselves are
the written linear-algebra proofs in Sections 2--5.
