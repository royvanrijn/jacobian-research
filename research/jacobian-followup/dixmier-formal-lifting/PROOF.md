# Polynomial formal lifting of canonical pairs and symplectic maps

14 November 2026. **QFL1.** The unrestricted fiber lifting equation has an
explicit polynomial right inverse. Every compatible finite quantum jet of the
retained pair extends to all orders in hbar without completing the spatial
variables or introducing X poles. A separate polynomial de Rham argument
proves the same formal existence statement for every full polynomial
symplectic endomorphism in every rank.

The distinction that remains is effectivity: a formal series with polynomial
coefficients need not define a finite polynomial differential operator at
hbar=1. QRT1--QRT3 remain valid in their declared restricted categories.
Their nonzero residue cannot extend to an obstruction in the unrestricted
fiber complex.

## An explicit right inverse for every fiber defect

Let B be a commutative Q-algebra with derivation delta, put A=B[Z], and
extend delta by delta(Z)=0. Define
\[
 \{f,g\}=f_Z\delta g-(\delta f)g_Z,\qquad
 f\star g=m\exp\!\left({\hbar\over2}
   (\partial_Z\otimes\delta-\delta\otimes\partial_Z)\right)(f\otimes g).
 \tag{1}
\]
The two derivations commute, so the product is associative. It is a
well-defined product on A[[hbar]]. In fact, on two individual polynomials
the bidifferential sum is finite because each term differentiates their
combined Z degree.

Fix any S,T in A with {S,T}=1. Write
\[
 d_1(u,v)=\{u,T\}+\{S,v\},\qquad
 I_Z\!\left(\sum_{k\ge0}F_kZ^k\right)
       =\sum_{k\ge0}{F_kZ^{k+1}\over k+1}.
\]
Then
\[
 {\cal H}(F)=\bigl(S_Z I_ZF,\ T_Z I_ZF\bigr),\qquad
                 d_1{\cal H}(F)=F
 \quad\hbox{for every }F\in A.                       \tag{2}
\]

Indeed, for any I in A the Leibniz rule and commutation of the derivations give
\[
 \begin{aligned}
 d_1(IS_Z,IT_Z)
 &=I\,\partial_Z\{S,T\}+S_Z\{I,T\}+T_Z\{S,I\}\\
 &=I\,\partial_Z\{S,T\}+I_Z\{S,T\}=I_Z.
 \end{aligned}
\]
Taking I=I_ZF proves (2). There is no degree cutoff, inverse Jacobian,
division by an X polynomial, or spatial completion in this argument.
It also applies to B=K[X,X^{-1},Q] with arbitrary finite X poles.

For the retained seed,
\[
 B=K[X,Q],\qquad
 \delta=3X^2\partial_X+(2-6XQ)\partial_Q.
\]
Thus unrestricted polynomial fiber corrections already suffice.

## Every compatible finite jet extends formally

Suppose a finite pair S_h,T_h reduces to S,T and satisfies
\[
 {[S_h,T_h]_\star\over\hbar}
       =1+\hbar^mF_m+O(\hbar^{m+1}).
\]
Replace it by
\[
 (S_h,T_h)-\hbar^m{\cal H}(F_m).                     \tag{3}
\]
The coefficient at order m changes by -d_1 H(F_m)=-F_m; lower orders
are unchanged. Iterate. At every step the new coefficients belong to B[Z].
Completeness in hbar gives
\[
 S_h=S+\sum_{m\ge1}\hbar^m s_m,\quad
 T_h=T+\sum_{m\ge1}\hbar^m t_m,\qquad
 [S_h,T_h]_\star=\hbar.                              \tag{4}
\]
This extends any already compatible finite jet, including odd branches.

If the prescribed jet is even in hbar, the normalized Moyal commutator is
even. Equation (3) can then be applied only at orders 2,4,6,..., giving an
even lift. Odd corrections are therefore unnecessary for unrestricted
formal existence in this fiber problem.

Each coefficient is polynomial. No common bound on their Z orders, total
degrees or denominators is asserted. The result concerns the S,T relation.
If R is central in (1), its commutators with these lifts stay zero; the
connection equations involving a fourth image are separate.

## Explicit escape from the QRT3 category

Use the retained polynomial lower lift s,t from
[the frozen inputs](../dixmier-residue-transport/inputs.json). Its order-five
defect is
\[
 O_5=\{s,t\}+\Pi^3(s,T)/24+\Pi^3(S,t)/24+
                        \Pi^5(S,T)/1920.
\]
QRT3 gives Lambda_plus(O_5)=1 and annihilates its allowed current
corrections. Equation (2) instead gives the unrestricted corrections
\[
 u=-S_ZI_ZO_5,\qquad v=-T_ZI_ZO_5,\qquad
                  d_1(u,v)=-O_5.                   \tag{5}
\]
In particular Lambda_plus(d_1(u,v))=-1. The same formula works above
every compatible polynomial or finite-X-Laurent lower lift; it does not
depend on a selected point of that lower-lift space.

For the frozen lower lift the exact sizes are:

| Polynomial | Nonzero terms | Z order |
|---|---:|---:|
| O_5 | 55 | 4 |
| u | 162 | 9 |
| v | 135 | 8 |

Both corrections are polynomial in X,Q,Z, with minimum X exponents 12 and
14. QRT3 requires current Z orders at most 1 and 0. Thus (5) leaves that
category explicitly. It does not contradict QRT3.

More generally, for this seed deg_Z(S)=5 and deg_Z(T)=4. If
deg_Z(s)<=3 and deg_Z(t)<=2, then deg_Z(O_5)<=4, and (5) has Z orders
at most 9 and 8. This is an all-lower-lift bound, not a finite scan.

Starting at the classical seed, the same homotopy gives a fresh
third-order correction with Z orders 11 and 10. It leaves the restricted
lower-lift category but starts an all-orders even formal solution.

The [Fraction replay](verify_fiber_homotopy.py) builds the exact
[escape witness](escape_witness.json). A separate
[SymPy replay](check_independent.py), importing no producer or legacy
algebra, verifies the general calculus identity before imposing
{S,T}=1, independently reconstructs O_5 and (5), and checks the two
residue values. The [receipts](receipt.json) and
[independent coverage](independent.json) distinguish those exact identities
from the written all-orders induction.

## Full polynomial symplectic frames are formally unobstructed

Let A=K[x_1,...,x_{2n}] in characteristic zero with a constant nondegenerate
Poisson matrix J and the corresponding Moyal product. Suppose
f=(f_1,...,f_{2n}) is polynomial and
\[
                       \{f_i,f_j\}=J_{ij}.          \tag{6}
\]
Then there are polynomials f_{i,2r} such that
\[
 F_i=f_i+\sum_{r\ge1}\hbar^{2r}f_{i,2r},\qquad
                         [F_i,F_j]_\star=\hbar J_{ij}.
 \tag{7}
\]
Consequently f lifts to a continuous endomorphism of the hbar-adically
completed Weyl algebra. Any compatible finite jet extends; an even jet
admits an even extension.

Here is an elementary proof, requiring no assumption that f is invertible.

Put D_i={f_i,-}. The D_i commute by (6). Moreover the identity
Df J (Df)^t=J makes Df an invertible matrix over A. The D_i therefore
form an A-basis of Der_K(A). Their dual one-forms theta^i are constant
linear combinations of df_j, so they are polynomial and closed.
The cochain complex of this commuting frame is exactly the polynomial
de Rham complex, expressed in the theta basis:
\[
 (dc)_{ij}=D_i c_j-D_j c_i,\qquad
 (dE)_{ijk}=D_iE_{jk}+D_jE_{ki}+D_kE_{ij}.          \tag{8}
\]
In particular H^1=H^2=0. Explicitly, a closed polynomial p-form whose
coefficients are homogeneous of ordinary degree q has primitive
i_Euler(alpha)/(q+p), for p>0. Summing its homogeneous pieces gives a
polynomial primitive. This uses the original source coordinates and the
inverse Jacobian matrix, never an inverse polynomial map f^{-1}.

Suppose the normalized commutators of a compatible jet are
J_ij+hbar^m E_ij+O(hbar^{m+1}). Their Jacobi identity at the first
remaining order gives dE=0 in (8). Vanishing of H^2 supplies polynomial
c_i with dc=-E. Adding hbar^m c_i to every image removes that defect.
Induction proves (7); parity follows from the Moyal commutator as before.

Vanishing of H^1 also shows that any two such formal lifts with the same
classical map differ by a formal Hamiltonian gauge on the target. At their
first differing order the difference is a closed one-cochain, hence
c_i=D_i g for a polynomial g. Composing with the appropriately signed
exp(hbar^m ad_star(g)/hbar) removes it. The infinite product converges
hbar-adically. This is a formal gauge; its action need not be a finite
polynomial automorphism at hbar=1.

For n=2, this proves simultaneous solvability of all six relations once
a full polynomial classical symplectic frame is supplied and all four
images may be corrected. It does not assert that the retained fiber pair
has such a global frame, that a fourth image is preserved unchanged, or
that an arbitrarily prescribed pair-only jet extends to a four-image jet.
Those admission and compatibility questions precede application of (7).

This also identifies an exact limitation of a future mixed-polarization
calculation: it cannot produce an absolute formal obstruction for a full
polynomial symplectic map when its correction spaces are unrestricted.
A nonzero obstruction must retain a substantive finiteness, filtration,
boundary or arithmetic condition.

## Comparison with graph quantization

[Baranovsky, Ginzburg, Kaledin and Pecharich, Theorem 1.1.4](https://arxiv.org/html/1403.3493v2)
give existence criteria for a line bundle on a smooth closed Lagrangian
subvariety to quantize as a complete flat module over a fixed formal
quantization. For the graph of (6) in the product with opposite symplectic
form, the graph is affine space. Its canonical bundle is trivial, and
both degree-two cohomology groups in their criterion vanish by the
polynomial Poincare lemma. Their criterion therefore permits the trivial
graph line bundle to quantize. This is consistent with (7), whose regular
left module and right action through F give a formal graph bimodule.

This comparison concerns formal deformation quantization. It is not a
theorem about a characteristic-zero D-module with a prescribed family of
positive-characteristic arithmetic supports.

## Effectivity and the remaining Dixmier questions

The hbar-adic algebra does not admit evaluation at hbar=1:
1-hbar is a unit, so its quotient by (hbar-1) is zero. A polynomial
hbar model, or another algebraic model regular at a nonzero parameter,
requires a separate argument. A formal series displaying growing degree
does not prove that every gauge representative has growing degree.

Three distinct gaps remain.

1. **Finite realization.** Prove termination, a suitable positive Rees
   bound, or algebraization in an algebra regular at a nonzero hbar.
   These must control all four images and their relations.
2. **Arithmetic realization.** If prescribed Frobenius center maps or
   arithmetic supports are required, prove that compatibility separately.
   Classical reduction modulo hbar alone does not establish it.
3. **Dixmier coverage and nonsurjectivity.** A construction needs a proof
   that its finite generator images fail to generate A_n. A proof of DC2
   needs a reduction covering arbitrary endomorphisms; one seed and its
   canonical orbit do not give that reduction.

The October 2026
[Kanel-Belov--Fedorov preprint, Theorem 4.2](https://arxiv.org/html/2610.33919v1)
excludes a finite characteristic-zero lift of Long's specified Frobenius
center map, even at matrix size m>=1. Theorem (7) nevertheless gives a
formal Moyal lift of Long's polynomial symplectic map. Taken together,
these results rule out a universal passage from formal lifting to that
prescribed arithmetic realization. They do not rule out every polynomial
hbar deformation having Long's map as its classical reduction.

For the plane graph-module route, the same preprint's Theorem 1.8
replaces arithmetic existence by bounded Berest--Wilson invariants along
faithful arithmetic families. This remains a concrete all-degree
boundedness problem. It concerns generalized DC1 and JC2 and does not
supply a proof of ordinary DC2.

[Lauritzen--Thomsen](https://arxiv.org/html/2308.09384v2)
give another useful endpoint: finite generation of A_n on either side
over an endomorphism image forces surjectivity. Their automatic
bimodule finiteness in characteristic zero is weaker. Controlling
one-sided generation, or an arithmetic transfer with the needed
uniformity, is a substantial remaining theorem.

Accordingly, further positional degree sweeps and a search for an
unrestricted formal mixed-polarization residue are poor next targets.
The next useful theorem must identify a finite-realization condition,
prove it is necessary for the intended Weyl endomorphisms, and establish
its invariance under the allowed changes of coordinates. QRT3 can then
serve as a restricted obstruction test rather than as an asserted
absolute obstruction.

No novelty priority, formal verification, external review, finite A2
endomorphism, or solution of DC1/DC2 is claimed.
