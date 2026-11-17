# Coupled triangular projective curvature is rigid

The one-curved-coordinate theorem extends to coupled curvature in any
number of coordinates. The key is an affine eigen-coordinate lemma for
a cyclic projective vector field. It turns the polynomial multiplier's
divisor into an affine coordinate, allowing descending elimination.

The resulting theorem closes an entire homogeneous architecture in every
rank and degree. It does not classify arbitrary primitive directions,
nonhomogeneous inputs, or the full Keller monodromy spectrum.

## 1. An affine eigen-coordinate lemma

**KMFI128.** Let v be a rational vector field in n variables over
algebraically closed characteristic-zero constants. Put B=Dv and W=v-x.
Assume W is cyclic for B over k(x):
\[
                   \det(W,BW,\ldots,B^{n-1}W)\ne0.
\]
If a rational function q and a constant c satisfy
\[
                     dq\,B=c\,dq,\qquad W(q)=(c-1)q,      \tag{1}
\]
then q is affine.

Write H=D^2q. Differentiating the first identity gives HB symmetric:
the other derivative term is a contraction of the symmetric second
derivatives of v, and cH is symmetric. Thus
\[
                            HB=B^tH.
\]
Differentiating the second identity and using DW=B-I and the first gives
\[
                               HW=0.
\]
It follows H B^j W=(B^t)^j HW=0 for every j. Cyclicity implies H=0.
A characteristic-zero rational function with zero Hessian is affine.

The lemma is not restricted to triangular or polynomial v. Cyclicity
is essential. For v=(x,2t+x^2), the nonlinear q=t+x^2 satisfies(1)
with c=2, but W=(0,t+x^2) is not cyclic.

## 2. The complete homogeneous architecture

Let r>=2, n=r-1, p=(u,x1,...,xn), and d>=1.
Let a_i be arbitrary constants and Q_i arbitrary polynomials in
x1,...,x_(i-1) of total degree at most d. In particular Q_1 is constant.
Consider
\[
 V_0=u^d,\qquad
 V_i=a_i u^{d-1}x_i+
                 u^d Q_i(x_1/u,\ldots,x_{i-1}/u).         \tag{2}
\]
There can be nonlinear couplings in every later coordinate.
The displayed homogeneous tuple need not be primitive.

**KMFI129.** Suppose H=g-fp=FV for a nonzero homogeneous polynomial
scalar F, the sweep S=X+gamma*(1,p), X=(f,g), satisfies detDS=gamma^r,
and rankDX=r. Then every Q_i is affine. All a_i are nonzero,
pairwise distinct and different from1. Consequently
\[
                       V=u^{d-1}Mp                       \tag{3}
\]
for an invertible constant lower triangular matrix M with distinct
eigenvalues 1,a_1,...,a_n. Every survivor belongs exactly to the
KMFI121 factorized class with scalar G=F u^{d-1}.

Thus no primitive degree>=2 direction in the entire architecture(2)
admits a full-rank normalized multiplier in any scalar degree.
Every degree11 field here is S11,A11 or regularC11 by KMFI122;
neither exceptional PSL2(11) action occurs. The result includes
all repeated parameters and resonances a_i=d.

The surviving scalar multipliers are precisely those of KMFI121
for the matrix M, with the additional polynomial divisibility condition
u^(d-1) divides G. Since the first eigen-coordinate is u, the positive
monomial exponent of u must be at least d-1. Conversely these conditions
produce normalized full-rank inputs. This does not assert existence for
every affine collection Q_i or every choice of constants.

The restriction in(2) is substantive: after dehomogenizing, the vector
v_i=a_i*x_i+Q_i is triangular with constant diagonal derivative.
General triangular polynomial maps with nonconstant diagonal,
general nontriangular primitive tuples, and nonhomogeneous fields
are not covered.

## 3. Reduction to a cyclic projective field

The same determinant and Euler argument as KMFI125 applies without
first removing the coordinate gcd. It gives
\[
 f=-Fh,\quad h=\lambda u^{d-1},\quad\lambda\ne0,\qquad
                        \prod_i a_i\ne0.
\]
Indeed detDV=d*(product a_i)*u^(r(d-1)); full rank excludes f0.
Set eta=dlogF. On u=1 write
\[
 v=(a_i x_i+Q_i),\quad B=Dv,\quad W=v-x,\quad
 A=DV|_{u=1}=\begin{pmatrix}d&0\\dv-Bx&B\end{pmatrix},
 \quad V|_{u=1}=(1,v).
\]
The characteristic identity is
\[
 \det(TI-A-V\eta)=(T-\lambda)^r,\qquad
                      \chi_A=(T-d)\prod_i(T-a_i).         \tag{4}
\]
Its right side and chi_A have constant coefficients in x.
The full-rank Krylov criterion KMFI124--125 makes
C=(V,AV,...,A^nV) invertible on this chart.

There is a useful exact determinant identity:
\[
        \det C=\det B\ \det(W,BW,\ldots,B^{n-1}W).          \tag{5}
\]
Replace the columns after V by
(A-dI)V,A(A-dI)V,...,A^(n-1)(A-dI)V.
This is a unit triangular column operation. But
(A-dI)V=(0,BW), so(5) follows immediately.
Thus B is invertible and W cyclic.

For every T the projective resolvent identity is
\[
 [(TI-A)^{-1}V]_0=(T-d)^{-1},
\]
\[
 [(TI-A)^{-1}V]_{\rm high}
       =\frac{x+T(TI-B)^{-1}W}{T-d}.                      \tag{6}
\]
Let s=degF and Fbar=F(1,x). Euler gives eta_u=s-dlogFbar*x.
Consequently (4), the determinant lemma and(6) yield
\[
 \frac{s+T\,d\log Fbar\,(TI-B)^{-1}W}{T-d}
                              =1-\frac{(T-\lambda)^r}{\chi_A}. \tag{7}
\]
All coefficients on the right are constant in x.

## 4. Descending removal of a last curved coordinate

At an induction stage let x_k be the last curved coordinate, and suppose
the later coordinates have already been changed by an affine triangular
projective coordinate change to v_j=a_j*x_j for j>k. These later a_j
are distinct from1 and all earlier a_i, by the previous stages.
The earlier prefix remains triangular with constant diagonal.

Write c=a_k, t=x_k, v_k=c*t+Q(x_1,...,x_(k-1)).
Let B_0,W_0 denote the earlier projective prefix,
chi_0(T)=det(TI-B_0), and chi_s(T)=product_(j>k)(T-a_j).
Multiplying(6) by chi_A and evaluating at T=c shows all entries
of adj(cI-A)V except the k-th high entry are zero. That entry is
\[
 L= c\,\chi_s(c)\{(c-1)\chi_0(c)t+
                 \chi_0(c)Q+dQ\,\operatorname{adj}(cI-B_0)W_0\}.
                                                               \tag{8}
\]
This is a polynomial identity, valid also when c=d or c equals
an earlier eigenvalue. Empty prefixes and suffixes have determinant1.

The vector adj(cI-A)V is nonzero because V is cyclic for A.
Otherwise its resolvent numerator would have a common T-c factor,
giving an annihilating polynomial of degree at most r-1 for V.
Thus L is not the zero polynomial.

Equation(4) at T=c gives exactly
\[
                            \eta_k L=-(c-\lambda)^r.      \tag{9}
\]
If L has no t term, its right-hand expression for eta_k is independent
of t. Polynomiality of F makes its logarithmic derivative tend to0
at t infinity, hence eta_k=0 and lambda=c. Then F is independent of t;
f and every g coordinate are also independent of t, since the k-th
linear term cancels and later coordinates are diagonal.
This contradicts full rank.

It follows the coefficient of t in L is nonzero. Hence c differs
from1 and every earlier eigenvalue; it already differs from the later
ones. If lambda=c, (9) again gives eta_k=0 and the same rank collapse,
so lambda!=c.

Set
\[
 S=\frac{Q+dQ(cI-B_0)^{-1}W_0}{c-1},\qquad q=t+S.         \tag{10}
\]
Because det(cI-B_0) is a nonzero constant, S is polynomial in the prefix.
Equation(9) and polynomiality give
\[
                  Fbar=J\,q^\nu,\qquad\nu\in\mathbb Z_{>0}, \tag{11}
\]
where J is independent of t and can depend on later coordinates.
Indeed eta_k is a nonzero constant divided by q; a polynomial's
logarithmic derivative has its positive integer multiplicity as residue.

## 5. The multiplier divisor is an affine eigen-coordinate

The left c-eigenrow of B with k-th entry1 and later entries0 is
\[
                \ell=(dQ(cI-B_0)^{-1},\,1,\,0_{\rm later}).
\]
It satisfies
\[
                    \ell B=c\ell,\qquad \ell W=(c-1)q.   \tag{12}
\]
At q=0 the pole of dlogFbar in(7) has residue nu*dq.
Since its right side has no pole in x, comparison gives, on q=0,
\[
                          dq\,B^jW=0\quad(j\ge0).        \tag{13}
\]
The cyclic span of W for B on that divisor has dimension n-1.
Here is the needed rank justification, not a generic-rank assumption
on a special divisor. Projection forgetting t is a B-equivariant
quotient onto diag(B_0,diag(a_j:j>k)), with the projected vector
(W_0,((a_j-1)x_j)_(j>k)). Since W is cyclic, this quotient vector
is cyclic. The quotient matrix and vector are independent of t,
so remain cyclic generically on the graph q=0.
By(12), the full span on q=0 lies in ker ell, which has dimension n-1;
its projection already has that dimension. The two are equal.

Thus(13) makes dq proportional to ell on q=0.
Both have k-th entry1, so dq=ell there.
Their coefficients are independent of t, hence this equality
holds identically. Equations(12) now give
\[
                         dq B=c dq,\qquad W(q)=(c-1)q.
\]
Apply KMFI128: q is affine. Use z_k=q as an affine triangular
projective coordinate change, equivalently a constant linear homogeneous
change including u. The transformed coordinate is exactly v'_k=c*z_k,
since v'=W'+z and W(q)=(c-1)q. Earlier coordinates and later diagonal
coordinates are unchanged. The induction continues.

After every stage v is affine under an affine change of coordinates.
Therefore the original v was affine as well: every original Q_i is affine.
All a_i are distinct and different from1. Equation(3) follows,
and KMFI121--122 finish the classification and degree11 obstruction.

## 6. Evidence, boundaries and remaining work

[Helper](triangular_projective.py) ·
[checker](verify_triangular_projective.py) ·
[model](triangular_projective_model.json) ·
[receipt](triangular_projective_receipt.json).

Exact controls check whole projective resolvent and Krylov identities,
coupled curvature in several coordinates, repeated and resonant spectra,
parameter consistency with isolated values retained, genuine affine
full-rank survivors, and rank-deficient nonlinear survivors.
The nonlinear q countercontrol checks why KMFI128 requires cyclicity.

The argument is unbounded in rank, direction degree, scalar degree,
number of curved coordinates and number of polynomial terms.
Finite controls supplement the written proof, not replace it.
Formal verification, independent whole-proof replay and external
review are not claimed.

The full objective remains the faithful-transitive geometric monodromy
spectrum of polynomial Keller self-maps of ordinary affine space.
KMFI129 closes an architecture, not an absolute action.
Other homogeneous primitive directions, nonconstant projective diagonal,
general nontriangular and nonhomogeneous fields, exceptional11 absolutely,
and other undecided actions remain open.
