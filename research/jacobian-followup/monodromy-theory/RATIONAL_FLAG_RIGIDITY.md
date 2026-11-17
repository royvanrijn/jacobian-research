# Rational constant-diagonal flags reduce to the linear architecture

The one-fiber argument extends to a nonlinear rational head: its graph
shift is linear, so an allowed linear shear decouples the fiber.
Descending a full triangular flag then removes every rational nonlinear
term. This classifies a multi-fiber perturbation architecture in every
rank and degree, including primitive cubic directions with critical
factor a product of two distinct linear forms.

These are normalized-input construction results, not a complete action
spectrum or an absolute exceptional11 obstruction.

## 1. A constant linear fiber decouples from a rational head

Let p=(q,z), r=m+1>=2. Let v(q) be rational homogeneous of degree1,
T(q) likewise, and mu,lambda nonzero constants. Put
\[
                         W=(v(q),\mu z+T(q)).
\]
Suppose G is a nonzero homogeneous polynomial of positive degree n,
the reconstructed X=(-lambda*G,G(W-lambda*p)) is polynomial,
and it is normalized with generic rankDX=r.

**KMFI166.** There is a constant row L such that
\[
                         T=\mu Lq-Lv.
\]
The simultaneous constant shear z_new=z+Lq therefore decouples
W_new=(v(q),mu*z_new). Moreover
G=J(q)(z+Lq)^nu with J a polynomial and nu a positive integer.
The result does not assert that T itself is linear when v is nonlinear.

The rational Euler/cyclicity argument of KMFI164 applies unchanged:
A=DW and W are cyclic. Also
\[
 D(GW)=A(GI+p\,dG),\qquad
             \det A=\lambda^r/(n+1)\ne0.
\]
The last equality follows from normalization and Euler.
Consequently B=Dv is invertible, and the head projection makes
v,Bv,... cyclic. Since Bq=v, q is cyclic for B.

The proof of KMFI162 over k(q) gives mu!=lambda,
det(muI-B)!=0, nu>0 and
\[
                    G=J(q)(z+S(q))^\nu,\qquad
                    S=dT(\mu I-B)^{-1}q.                 \tag{1}
\]
This remains valid for the nonconstant rational matrix B.
All quantities in(1) have the required homogeneous degree, so S has
degree1. The complete graph residues give
\[
                              dT=\mu dS-dS B.             \tag{2}
\]
Their recurrence is the same as in KMFI164, with C replaced pointwise
by B(q); powers here are matrix powers, not iterated derivations.

The nonlinear-head Hessian terms must now be retained.
Differentiating(2) gives exactly
\[
 H_T=\mu H_S-B^tH_S-\sum_i S_i H_{v_i}.
\]
All terms except B^tH_S are symmetric, hence
H_S B=B^tH_S. Homogeneity gives H_S q=0.
Pointwise matrix transport then kills H_S B^j q for every j;
cyclicity makes H_S=0. Thus S=Lq with constant L.
The zeroth graph residue is dS*v-mu*S+T=0, giving the formula asserted.

Since the shift is now linear and G polynomial, its leading z
coefficient J is polynomial too. Only this final constant shear is
used as a source/output equivalence; no nonlinear rational change
is claimed to preserve normalization.

For example v=(x,y,w+x^2/y), mu2 and T=w-x^2/y have S=w.
T is nonlinear rational, but z_new=z+w removes it exactly.
This is a coordinate identity, not an asserted normalized full-rank
input. A noncyclic head can also allow nonlinear S at the level of
the residue identities; full-rank cyclicity is essential.

## 2. Complete rigidity of a rational triangular constant-diagonal flag

Over algebraically closed characteristic-zero constants, suppose
\[
 W_i=a_i p_i+T_i(p_1,\ldots,p_{i-1}),\quad 1\le i\le r,
                                                               \tag{3}
\]
where each T_i is rational homogeneous of degree1. Thus T1=0.
Let G be a positive-degree homogeneous polynomial clearing the input
denominators so X=(-lambda*G,G(W-lambda*p)) is polynomial,
normalized and of generic rank r.

**KMFI167.** Every T_i is linear. The a_i are distinct nonzero constants
different from lambda. The input is exactly among the scalar-times-linear
KMFI121 data, with its exact generic field and degree11 consequences.

First the determinant identity in Section1 forces product a_i nonzero.
Apply KMFI166 to the last coordinate.
Its multiplicity branch requires
det(a_r I-DW_head)=product_(i<r)(a_r-a_i) nonzero,
so no repeated last eigenvalue is discarded by division.
A linear shear using only earlier variables decouples p_r.

Proceed downward. To isolate the next fiber, put already decoupled
coordinates in the head. Its graph shift is a linear combination of
earlier variables and of those decoupled coordinates.
The relation T=mu*Lq-Lv, and the distinct fiber eigenvalues, force
the coefficients of the already decoupled variables to be zero.
Thus another shear using only earlier coordinates decouples the next
fiber, preserving those already separated.
At each step repeated eigenvalues are rejected by the explicit
zero-leading-coefficient branch of the multiplicity test.
Eventually every coordinate is linear and the whole map is linear.
A constant linear conjugacy cannot turn a genuinely nonlinear
rational tuple into a linear one, so every original T_i was linear.

The remaining scalar, distinctness from lambda, exact field and
arithmetic assertions are precisely KMFI121--122.
No additional all-degree abstract group classification is inferred.
In degree11 this architecture gives only the retained symmetric,
alternating or regular cyclic actions in the stated ranks, and
cannot supply exceptionalPSL2(11).

The conclusion requires a constant-diagonal flag under allowed constant
conjugacy. It does not apply merely because a rational map can be
triangularized by a nonlinear birational substitution.

## 3. A complete multi-fiber perturbation classification

Let Q be nonconstant homogeneous of degree k over algebraically closed
characteristic-zero constants, and let U be its full space of constant
translation directions, of dimension c>=1.
Choose a basis as columns of a constant matrix B0.
Let R=(R1,...,Rc) be homogeneous of degree k+1 with every component
invariant under every direction in U. For M in GL_r put
\[
                         V=M(Qp+B0R).                    \tag{4}
\]

**KMFI168.** A normalized full-rank homogeneous scalar multiplier can
exist only if Q divides every Ri.
In that case V=Q*M'*p, and its multipliers are exactly the old KMFI121
data with scalar G=FQ divisible by Q.
Thus no primitive nonlinear direction in(4) survives.
If Q has any irreducible nonlinear factor, there are no full-rank
multipliers at all. This includes the irreducible one-kernel case
KMFI165 and allows arbitrarily many fiber directions and perturbations.

The determinant lemma, using dQ*B0=0,dR*B0=0, gives
\[
                      \det DV=(k+1)\det M\,Q^r.
\]
For a hypothetical homogeneous scalar of degree s, normalization gives
f=-lambda*FQ and lambda^r=(s+k+1)detM.
For every b in U, partial_bV=Q*M*b and the trace cofactor is
independent of b. KMFI160 forces partial_(Mb)Q=0.
Hence M maps U into itself, and invertibility makes this equality.

Choose coordinates (q,z) with U the z space. Then
\[
 M=\begin{pmatrix}C&0\\L&D\end{pmatrix},\quad C,D\text{ invertible},
\]
and, after setting G=FQ and W=V/Q,
\[
                         W=(Cq,Dz+Lq+D R/Q).
\]
Triangularize C and D by constant basis changes over the algebraically
closed constants. The whole direction is now a rational homogeneous
degree1 flag with constant diagonal: all rational perturbations depend
only on q. KMFI167 makes the entire W linear.
Since D is invertible, every Ri/Q is linear, proving the divisibility.

Write R=Q*L1*p with a constant matrix L1. Its invariance gives L1*B0=0.
Then M'=M(I+B0*L1) is invertible, with detM'=detM.
The accepted scalar is G=FQ in the factorized architecture.
Conversely any such KMFI121 data with G divisible by Q yield the
original presentation via F=G/Q.
A nonlinear irreducible factor of Q cannot divide the completely split
scalar of KMFI121, giving the final obstruction.

Requiring every Ri to be invariant under the whole U is substantive.
Invariance under one selected direction in a larger kernel does not
justify M(U) invariant or the rational-flag reduction.

## 4. A primitive cubic over two distinct critical hyperplanes

Take p=(u,x,y,z), Q=u*x, U=span(e_y,e_z),
R=(u^3,x^3), and M=diag(1,2,3,4).
Then
\[
            V=(u^2x,\ 2ux^2,\ 3uxy+3u^3,\ 4uxz+4x^3)
\]
is primitive,
\[
                        \det DV=72u^4x^4,
\]
and its Krylov determinant at p=(1,1,1,1) is1056.
Every projective base point has u=x=0, where DV is the zero matrix.
Translations along U also pass the basic directional test:
Q is invariant under both their source and target directions.

Nevertheless KMFI168 excludes every scalar degree because Q does not
divide either perturbation. This is a reducible-critical-quadric cubic
obstruction distinct from the irreducible-quadric case of KMFI165.
For the last rational fiber, T=4x^2/u and its graph shift is
S=8x^2/(3u). The gradient defect at the unit head point is(4,-8/3,0),
so the necessary full residue condition already fails.

The fully split, nonprimitive boundary does have normalized full-rank
inputs. In rank3, take beta^2+1=0,
\[
 M=\operatorname{diag}(1/2,(1+\beta)/2,(1-\beta)/2),\quad
 Q=x y,\quad R=0,\quad F=z,\quad G=xyz,\quad\lambda=1.
\]
This is KMFI121 data: its entire sweep Jacobian is gamma^3
with source order(gamma,x,y,z), and a unit-point full-rank minor is -1.
It is a normalized-input control, not a new Keller occurrence action.

## 5. Evidence and the full objective

[Helper](rational_flag.py) · [checker](verify_rational_flag.py) ·
[model](rational_flag_model.json) · [receipt](rational_flag_receipt.json).

The unbounded proof retains nonlinear-head Hessians, every graph moment,
all repeated eigenvalue rejection branches, allowed constant shears,
all kernel directions and every perturbation component.
Exact whole identities and normalized positive controls supplement it.
Formal verification, independent whole-proof replay and external review
are not claimed.

The previous goal turn was progress: KMFI164--165 gave linear-head
rational fiber rigidity and an irreducible one-kernel obstruction.
This turn gives nonlinear-head fiber decoupling, complete rational-flag
rigidity and the multi-fiber classification above.
Arbitrary rational degree1 directions without such a flag, cubic
directions outside this architecture, variable-diagonal rational flags,
nondegenerate critical quadrics, nonhomogeneous fields and the absolute
exceptional11 action remain unresolved.
The complete faithful-transitive Keller action spectrum remains active.
