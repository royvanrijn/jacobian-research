# Differential reconstruction and the nonlinear-flag obstruction

The remaining exceptional11 normalized input must be a nonfinite
rank4 sweep with full rankDX. This note gives an exact algebraic
parameterization of every normalized sweep, a pointwise full-rank
test, and a new obstruction to whole classes of matrix inputs.
It does not assume a triangular first-coordinate presentation.

These are necessary-and-sufficient integrability and rank statements,
not a sufficient abstract-action criterion. No PSL2(11) input or
absolute exceptional11 action obstruction is supplied.

## 1. Reconstructing every polynomial normalized sweep

Write X=(f,g), with g=(g1,...,gr), p=(p1,...,pr), and define
\[
 S(\gamma,p)=X(p)+\gamma(1,p),\qquad
                   P=Dg-p\,df .                           \tag{1}
\]
Then
\[
                  \det DS=\det(\gamma I+P).                \tag{2}
\]
Here df is a row and p a column.

**KMFI115.** The entire normalized class detDS=gamma^r has
the following polynomial coefficient parameterization.
Choose an arbitrary vector polynomial g whose linear part has
trace zero. Decompose it by ordinary homogeneous degree and set
\[
 f=f_0+\sum_{n\ge1}\frac{\operatorname{div}g_{n+1}}n,\qquad
 P_0=Dg_1,\quad
 P_n=Dg_{n+1}-\frac p n\,d(\operatorname{div}g_{n+1}).       \tag{3}
\]
Then impose just the polynomial identities
\[
                    \operatorname{tr}(P^j)=0,
                         \qquad2\le j\le r.                \tag{4}
\]
For r1 this list is empty. These are exact equations on the
polynomial coefficients, with no degree cap asserted.
They produce every normalized polynomial X, with only its
irrelevant constant target translations undetermined.

Indeed trP=divg-Ef, where E=sum p_i partial_i.
On degree n>=1, Ef_n=n f_n, giving (3).
The constant term of divg is trDg_1 and must vanish.
Conversely (3) makes trP=0. Newton identities in characteristic0
make (4) equivalent to the remaining coefficients of
det(TI-P)=T^r, hence to nilpotency and (2).
Thus this is a global polynomial reconstruction, not a formal
coordinate normal form.

There is an intrinsic test for a proposed matrix P.
For r>=2 put
\[
 b_j=\sum_i\partial_iP_{ij},\qquad
                         a=-b/(r-1).                      \tag{5}
\]
For a traceless P arising from (1), a=df and
\[
 (r-1)(\partial_jP_{ik}-\partial_kP_{ij})
                         =\delta_{ij}b_k-\delta_{ik}b_j.   \tag{6}
\]
For r>=3, (6) is also sufficient for polynomial reconstruction.
It says that the row one-form beta_i=sum_jP_ij dp_j satisfies
d beta_i=-dp_i wedge a.
Applying d again gives dp_i wedge da=0 for every i.
When r>=3, each coefficient of da can be tested by choosing a
different i, so da=0. Polynomial integration gives f with df=a,
and each beta_i+p_i a is closed, giving g_i.
Contracting the curl identity proves (5), so no independent
first component is missing. All primitives are polynomial
in characteristic0.

For r2, (6) alone is insufficient: one must additionally impose
da=0. For P=[[-q,l],[-m,q]] this is exactly
\[
                         l_{11}+2q_{12}+m_{22}=0,          \tag{7}
\]
the earlier rank2 wave equation. The distinction r2 versus r>=3
must not be omitted. For r1 the matrix P=0 loses the information
in f; (3), not (5), is its reconstruction.

For r>=2 the homogeneous traceless solution space of (6),
including da=0 for r2, has dimension
\[
 \dim V_0=r^2-1,\qquad
 \dim V_n=r\binom{n+r}{r-1}\quad(n\ge1).                   \tag{8}
\]
Before nilpotency, degree n>=1 is the image of the arbitrary
r-vector g_(n+1) under the injective linear map (3).
Its kernel is zero: P=0 gives df=0 by (5), then Dg=0.
This proves (8) without a bounded rank calculation.
In rank4 the first dimensions are15,40,80,140.
The remaining nonlinear constraints are precisely trP^2,
trP^3 and trP^4; linear trace and integrability have been removed.

## 2. The coframe and the kernel-direction derivative

**KMFI116.** Let P be nilpotent and satisfy (1), over
characteristic0. At every source point,
\[
 \operatorname{rank}DX=r
 \quad\Longleftrightarrow\quad
 \det\begin{pmatrix}
 df\\dfP\\ \vdots\\dfP^{r-1}
 \end{pmatrix}\ne0.                                      \tag{9}
\]
At any such point P is a single regular nilpotent Jordan block.
Moreover, for r>=2, if u spans kerP and v spans the left kernel,
with P^(r-1)=u v after a choice of scaling, then
\[
             v\,D_uP=df(u)\,v,\qquad df(u)\ne0.            \tag{10}
\]
The derivative treats u as a fixed tangent vector at the point;
it is not a derivative of a varying kernel section.

To prove (9), DX has rows df and P+p df. If DX is injective,
its restriction to kerP injects into the one-dimensional span
of (1,p). Hence dimkerP1 and df is nonzero on that kernel.
For one Jordan block this is exactly the cyclic-covector
condition in (9). Conversely, DXw0 implies dfw0 and Pw0,
so every row of the displayed matrix kills w; its invertibility
forces w0. This proves the equivalence pointwise, including
rank-drop loci.

For (10), differentiate trP^r=0 in p_i:
tr(P^(r-1) partial_iP)=0, hence v(partial_iP)u=0.
Also vu=trP^(r-1)=0 for r>=2.
The curl identity from (1) gives
\[
 (\partial_iP)u=(D_uP)e_i+u\,\partial_i f-e_i\,df(u).
\]
Left multiplication by v proves (10).
Nilpotency has not made its directional derivative nilpotent:
a nonzero eigenvalue is forced along the kernel.

There is also a useful family of necessary differential identities:
\[
 \sum_j\partial_j(P^m)_{ji}
       =-(r-m)(dfP^{m-1})_i,\qquad1\le m\le r.             \tag{11}
\]
For a direct proof, fix A and take the polynomial map
\[
                   T_A(p)=g(p)+(A-f(p))p .
\]
Its derivative is P+(A-f)I, of determinant (A-f)^r.
The divergence of every column of its adjugate is zero:
on expanding a cofactor, mixed partials cancel in pairs.
Expand
adj(P+gamma I)=sum_(m=0)^(r-1)(-1)^m gamma^(r-1-m)P^m,
put gamma=A-f, and compare powers of A-f.
This gives (11), with the final m=r identity zero on both sides.
It supplies exact residual tests for arbitrary, including
nontriangular, candidate matrices.

## 3. A new obstruction: a fixed nilpotent matrix space cannot work

**KMFI117.** For r>=2, a full-rank normalized polynomial sweep
cannot have every kernel-direction derivative D_uP nilpotent
on its full-rank locus. In particular it cannot have its matrix
values confined to any fixed linear space W all of whose
matrices are nilpotent. It also cannot have affine-linear P.
Consequently no normalized sweep with all X coordinates of
ordinary degree at most2 has full rankDX=r.

The first assertion follows immediately from (10): df(u) is
a nonzero eigenvalue of D_uP.

For the linear-space assertion, DP(p)[u] still belongs to W,
so would be nilpotent. The same argument works over the rational
function field when the kernel vector has rational entries.
Every polynomial identity asserting nilpotency of matrices in W
remains valid after that field extension.

If P is affine-linear and nilpotent for every p, then every
matrix in span(P0,P1,...,Pr) is nilpotent.
For a nonzero coefficient c0 it is c0*P(c1/c0,...,cr/c0);
polynomial continuation in the coefficients includes c0=0.
Thus the preceding assertion excludes full rank without
assuming simultaneous triangularizability of a nilpotent
matrix space.

Finally traceP=0 gives Ef=divg.
If all g components have degree<=2, f has degree<=1;
there can be no higher homogeneous part in Ef.
Then P is affine-linear and rankDX<r.
This proves the last assertion and also covers the apparently
more general hypothesis deg g<=2.
The obstruction does not apply to arbitrary quadratic Keller
maps: it concerns this normalized X input and full-rank hypothesis.

The rank bound is sharp in every r>=2. Set
\[
 X=(p1,\ p1^2/2,\ p3,\ldots,pr,\ 0).                       \tag{12}
\]
It has rankDX=r-1. Its normalized matrix is affine and nilpotent;
its first row is zero, and rows i>=2 have -p_i in column1
and, for i<r, 1 in column i+1.
It is a regular nilpotent block generically, but df kills its
kernel. This shows why a regular Jordan matrix alone is not
the full-rank test.

In contrast, the rank4 h=t^2 full-rank control at p*=(1,0,0,0)
has kernel vector
\[
                     u=(1,1,4/3,2/15).
\]
Here df(u)=2, while D_uP has characteristic polynomial
\[
                  (T-2)(T-1)^2(T+4).                     \tag{13}
\]
The exact check shows the required moving flag explicitly.
A constant basis making every P strictly upper triangular
would put the image in a linear nilpotent space and is therefore
impossible for any full-rank normalized input.
This does not affect the KMFI104 choice of a Jordan frame at
one branch point.

## 4. Consequences for the remaining objective

The exceptional11 rank4 gate now admits an exact coefficient
formulation without postulating an integrable Jordan frame.
Choose a polynomial vector g with traceless linear part, compute
f and P by (3), impose trP^2=trP^3=trP^4=0, and require the
coframe determinant in (9) to be nonzero.
These conditions are necessary and sufficient for a polynomial
normalized full-rank input, of arbitrary degree.
They do not prove its generic degree or monodromy.

Any successful exceptional11 input must additionally meet
KMFI105's nonfinite boundary, birational branch, omitted index5
prime and irreducible unramified-tail conditions, and prove the
exact PSL2(11) action. By KMFI117 it must have a g component of
degree>=3, matrix degree>=2, and a matrix image outside every
fixed linear nilpotent space. Its kernel-direction derivative
must have the nonzero eigenvalue in (10).
These are genuine restrictions on every remaining full-rank
normalized input, not assumptions defining a triangular subclass.

The Haantjes-tensor replacement for the previous Nijenhuis route
also fails the full-rank integrated control: an exact nonzero
component is retained in the model. Neither type of local
normal-form theorem is used in this proof.

The written arguments are unbounded. Exact checks cover the
reconstruction and curl system, dimensions of representative
homogeneous linear spaces, the coframe and Piola identities,
the derivative eigenvalue and both full-rank and rank-deficient
controls. Finite checks do not classify all nonlinear nilpotent
solutions. No formal, independent whole-proof or external
review is claimed.

[Constructor](sweep_differential.py) ·
[checker](verify_sweep_differential.py) ·
[model](sweep_differential_model.json) ·
[receipt](sweep_differential_receipt.json).

This does not decide arbitrary rank-deficient architectures,
the absolute exceptional11 action or the full faithful-transitive
Keller action spectrum. The parent goal remains open.
