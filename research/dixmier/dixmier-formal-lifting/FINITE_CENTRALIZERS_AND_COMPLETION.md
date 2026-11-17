# Finite centralizers and finite quantum completion

15 November 2026. **QFC1.** A positive completion criterion in every rank:
once actual individual conjugates exist, module finiteness of the actual
joint centralizer forces it to be polynomial. The remaining curvature
then has an actual finite polynomial primitive. Generic module finiteness
already suffices and recovers an integral correction, preserving a flat
classical tuple.

This strengthens the polynomial-centralizer sufficient condition in
[QSC1](LONG_SPECTRAL_CURVATURE_AND_SCALARIZATION.md). It supplies no Long
pair satisfying the new finiteness hypothesis. Module finiteness is a
substantive condition, distinct from finite generation as an algebra,
finite fraction-field degree, or polynomiality after completion.

## The algebra and the positive criterion

Let K have characteristic zero, R0=K[[h]], L=K((h)), and let E be an
R0-algebra that is h-torsion-free and h-separated. Assume D=E/hE is
commutative; its Poisson bracket is the reduction of [ , ]/h.
Choose actual elements F1,...,Fm,T1,...,Tm in E such that
\[
 [F_i,F_j]=0,\qquad [F_j,T_i]=h\delta_{ij}.                 \tag{1}
\]
Write fi for the reduction of Fi. Assume the fi are algebraically
independent and their full joint Poisson centralizer in D is
\[
                   \operatorname{Cent}_{\{\,,\,\},D}(f_1,\ldots,f_m)
                         =K[f_1,\ldots,f_m].              \tag{2}
\]
Put
\[
 P=R0[F_1,\ldots,F_m],\quad B=\operatorname{Cent}_{E}(F_1,\ldots,F_m),
 \quad P_L=L[F_1,\ldots,F_m],\quad B_L=B[h^{-1}].            \tag{3}
\]
Algebraic independence and h-separatedness identify P with a polynomial
ring over R0. Central localization shows B_L is the centralizer in
E[h^-1].

**Finite centralizer theorem.** Under (1)--(2), the following are equivalent:
B is module-finite over P; B_L is module-finite over P_L; and B=P.
Under these equivalent conditions there are bi in P with
\[
 [T_i+b_i,T_j+b_j]=0,\qquad [F_j,T_i+b_i]=h\delta_{ij}.      \tag{4}
\]
If the reductions of the Ti already commute in the Poisson bracket,
the bi can all be chosen in hP. Thus the entire prescribed classical
tuple is preserved.

The scalar finite-PBW algebra and finite idempotent corners with
scalar rank-one classical reduction satisfy the separation, torsion
and commutative-reduction hypotheses. For Long's pair, (2) is QSC1.
An arbitrary partial graph module needs its own proof of these hypotheses;
a completed frame alone is not being used as an algebraic frame.

## Formal centralizer and actual translation derivations

The reduction of any b in B lies in K[f]. Subtract its polynomial lift
in P and divide by h. The result is again in B: h is a nonzerodivisor
and the Fi commute with b. Repeating gives a unique expansion in the
h-adic completion,
\[
              b=\sum_{\nu\ge0}h^\nu Q_\nu(F),\qquad
                       Q_\nu\in K[t_1,\ldots,t_m].         \tag{5}
\]
Uniqueness follows by taking the first nonzero h coefficient and using
independence of the fi. Separation makes the map from E into its
completion injective. Polynomial partial sums of two such expansions
commute, so B is commutative. The same induction works for the full
completed centralizer. Laurent localization permits only finitely many
negative h powers. In particular B_L is a domain.

Since D is commutative, every commutator in E is in hE. Division by h
is unique and remains in the actual algebra. The derivations
\[
 \nabla_i=-h^{-1}\operatorname{ad}(T_i):B\longrightarrow B,\qquad
                    \nabla_i(F_j)=\delta_{ij}              \tag{6}
\]
therefore preserve B and extend to B_L. Jacobi puts [Ti,Tj] in B, and
commutativity of B implies [nabla_i,nabla_j]=0 on B. On (5) these
derivations act coefficientwise as partial_(ti). Consequently
\[
                  \bigcap_i\ker(\nabla_i:B_L\to B_L)=L.    \tag{7}
\]
This excludes nontrivial algebraic constant extensions even though
L=K((h)) need not be algebraically closed.

## Two elementary differential algebra lemmas

Let A be any commutative Q-algebra. If an ideal I in A[t1,...,tm] is
stable under every coordinate partial derivative, then
\[
                      I=(I\cap A)A[t_1,\ldots,t_m].        \tag{8}
\]
For a polynomial q, the finite Taylor projector
\[
             q(0)=\sum_\beta\frac{(-t)^\beta}{\beta!}
                                      \partial^\beta q    \tag{9}
\]
lies in I when q does. Apply this also to all derivatives of q:
each coefficient of q belongs to I intersect A. This proves (8) in
every degree and over a ring with nilpotents as well.

If a finite module N over a noetherian A[t] has operators Di satisfying
Di(pn)=(partial_i p)n+pDi(n), all its Fitting ideals are stable under
partial_i. Indeed, over dual numbers, sigma(p)=p+epsilon*partial_i p
is a base automorphism, and 1+epsilon*Di is an invertible
sigma-semilinear module map. Functoriality and base change of Fitting
ideals give sigma(Fitt_k(N))=Fitt_k(N), hence derivative stability.
The required basic Fitting facts are
[Stacks, Lemma15.8.4](https://stacks.math.columbia.edu/tag/07ZA).

Over a characteristic-zero field, (8) says each such Fitting ideal is
zero or the unit ideal. If r is the generic rank, Fitt_(r-1)=0 and
Fitt_r=the unit ideal, so N is locally free of rank r. This follows
locally by choosing r generators: vanishing Fitt_(r-1) forces every
entry of the relation matrix to vanish. In particular a generically
zero finite module with these connections is zero. Integrability
is not needed for this Fitting argument.

## Integral module finiteness removes escaping support

Assume first that B is finite over P. The finite P-module N=B/P inherits
the connections (6). The reduction argument before (5) gives B=P+hB,
hence N/hN=0. Its zeroth Fitting ideal is, by (8),
\[
                      \operatorname{Fitt}_0(N)=JP,\qquad J\subset R0.
                                                               \tag{10}
\]
Fitting base change and N/hN=0 imply JP+hP=P. Evaluation at t=0 gives
J+(h)=R0. Since R0 is local, J=R0, so Fitt_0(N)=P and N=0. Thus B=P.

It would be invalid to apply Nakayama directly to P: h is not in the
Jacobson radical of R0[t]. For example P/(1-h*t1) is nonzero but equals
h times itself. It cannot carry a partial_1 connection, since its
annihilator would have to contain both 1-h*t1 and its derivative -h,
which generate the unit ideal. The connections remove precisely this
possible escape to t1=1/h.

## Generic module finiteness also suffices

Assume only that B_L is finite over P_L. Equations (6) and the Fitting
lemma make B_L a finite locally free P_L-module. Its relative
differentials Omega_(B_L/P_L) are finite over P_L and inherit connections
\[
       \nabla_i(b\,dc)=(\nabla_i b)\,dc+b\,d(\nabla_i c).    \tag{11}
\]
This respects d(p)=0 because partial_i p is again in P_L.
The domain B_L becomes a finite field extension after passing to
Frac(P_L), and that extension is separable in characteristic zero.
Thus Omega is generically zero and the same Fitting lemma gives
Omega=0. Hence B_L is finite etale over P_L.

We use the standard external fact that finite etale covers of affine
space over a characteristic-zero field are constant. Over an algebraic
closure this follows from Riemann existence and simple connectedness
of complex affine space, with invariance under algebraically closed
base extension; see [Milne, Remark3.3(c), Theorem3.4 and the fundamental
group exact sequence](https://www.jmilne.org/math/CourseNotes/LECc.pdf#page=27).
Descent identifies every such algebra with C tensor_L P_L for a finite
etale L-algebra C. Each element of C is killed by the derivations,
because it is separably algebraic over L. Equation (7) forces C=L.
Therefore B_L=P_L. The covering theorem is an explicit external
dependency; it is not inferred from a finite control.

Finally, B=P, not merely B_L=P_L. Express b in B as a finite polynomial
in the Fi with coefficients in L. If the minimum coefficient valuation
were negative, multiply by the opposite power of h and reduce modulo h.
A nonzero polynomial in the independent fi would become zero. Thus
all coefficients are in R0. This proves the theorem's equivalence of
the integral and generic finiteness conditions.

## Polynomial curvature has a finite integral primitive

Write
\[
                 P_{ij}(F)=h^{-1}[T_i,T_j]\in B=P.         \tag{12}
\]
Jacobi and (6) show partial_i Pjk+partial_j Pki+partial_k Pij=0.
Thus the polynomial two-form with coefficients Pij is closed.
The radial polynomial homotopy gives
\[
 b_j(t)=\sum_i t_i\int_0^1 s\,P_{ij}(st)\,ds,\qquad
                     \partial_i b_j-\partial_j b_i=P_{ij}. \tag{13}
\]
The integral is coefficientwise finite: a monomial t^alpha receives
the denominator |alpha|+2, invertible in K. Commutativity of B and (6)
now give (4). If the classical Ti commute in the Poisson bracket,
[Ti,Tj]/h is zero modulo h. Independence of the fi makes P h-saturated
in E, so all coefficients Pij and bj are h-divisible.

For rank two one may instead choose the particularly simple primitive
\[
              b_1=0,\qquad b_2=\int P_{12}(t_1,t_2)\,dt_1. \tag{14}
\]
No h-pole correction or new spatial completion is needed.

For E the actual scalar finite-PBW algebra, this is a full bounded
quantum lift, so [QFL2](BOUNDED_REALIZATION.md) applies to a noninvertible
classical map. For a finite projective corner it supplies a full corner
action; the separate hypotheses of
[QMC1](FINITE_GRAPH_CARRIERS_AND_MATRIX_REALIZATION.md) then give its
matrix conclusion. Neither case constructs the required Long input.

## A certificate from leading Poisson geometry

Here is a sufficient all-degree gate for the generic finiteness
hypothesis. In the scalar Weyl algebra over L, use the Bernstein
filtration and put di=deg Fi and gi=sigma(Fi). Suppose the gi are
algebraically independent and their entire polynomial Poisson
centralizer is module-finite:
\[
 \operatorname{Cent}_{\{\,,\,\},L[x,p]}(g_1,\ldots,g_m)
               \text{ is finite over }L[g_1,\ldots,g_m].  \tag{15}
\]
Then gr(B_L) is a graded L[g]-submodule of that centralizer. Since
L[g] is noetherian, gr(B_L) is finite over L[g]. Lift homogeneous
generators to B_L and subtract their polynomial Fi combinations in
descending filtration degree. Degree decreases at every step, so
this terminates and proves B_L finite over P_L. No degree sweep is
used. The independence assumption ensures gr(P_L)=L[g].
The same argument works for another nonnegative algebra filtration
if its associated graded algebra is commutative and the stated
Poisson and finite-module hypotheses actually hold. A projective
idealizer does not acquire these filtration hypotheses automatically.

In rank two, leading-symbol independence also gives an explicit
correction bound. If ei=deg Ti and C=[T1,T2] is nonzero, then
\[
 \deg b_2\le\deg(C/h)+d_1\le e_1+e_2-2+d_1,\qquad
 \deg(T_2+b_2)\le\max(e_2,e_1+e_2-2+d_1).                 \tag{16}
\]
There is no cancellation of distinct top weighted monomials in P(F),
because the gi are independent. Integrating in t1 raises weighted
degree by d1; the Bernstein commutator drops degree by at least two.
The bound is conditional on the actual certificate (15), not on a
finite approximation to the full Poisson centralizer.

## Positive controls with genuinely momentum dependent positions

In rank two set
\[
 F_1=x_1+h^2p_1^2,\quad F_2=x_2+h^2p_2^2,\quad
 T_1=-p_1,\quad T_2=-p_2+h^2F_1^2F_2.                    \tag{17}
\]
These are finite and have the identity full classical reduction.
Their leading symbols h^2*p_i^2 are independent over L. Their full
Poisson centralizer is L[p1,p2]: commutation with p_i^2 forces
partial_(xi)=0 in the polynomial domain. This is free of rank four
over L[p1^2,p2^2], with parity basis 1,p1,p2,p1*p2. Therefore (15)
passes in every degree. Direct PBW multiplication gives
\[
       [T_1,T_2]=-2h^3F_1F_2,\qquad
       b_2=-h^2F_1^2F_2,\qquad T_2+b_2=-p_2.             \tag{18}
\]
This checks nonzero curvature and an h-small finite correction,
not a noninvertible Long graph. The checker also tests a closed
polynomial two-form and radial primitive in three variables.

## Why the finiteness gate is substantial

The escaping-boundary algebra in QSC1 has
B=R0[x,y,(1-h*x)^-1,(1-h*y)^-1]. It is finitely generated as an
R0-algebra and has polynomial formal centralizer, with individual
conjugates, but is not module-finite over R0[x,y] or after inverting h.
For example the valuation at 1-h*x is negative on its inverse,
whereas every element integral over the polynomial ring has
nonnegative valuation. Its nonzero curvature residue prevents
finite or rational flattening. It is an Ore boundary control, not
a subalgebra example in the finite Weyl algebra.

Even actual full Weyl relations and identity classical reduction
do not force the new finiteness hypothesis in all ranks.
For any Keller map f, the exact cotangent construction has
\[
 \operatorname{Cent}_{A_n}(f_1(x),\ldots,f_n(x))=K[x],
 \quad K[x]\text{ finite over }K[f]\ \Longleftrightarrow\
                                  f\text{ is an automorphism}. \tag{19}
\]
Indeed the Jacobian is invertible, so the leading differential symbol
of a positive-order centralizer element would have all its momentum
partials zero; hence every centralizer element has order zero.
Finiteness makes f a finite etale self-map of affine space, and the
same covering theorem makes it an automorphism. The converse is
immediate. The inverse-Jacobian conjugates commute already.

For a normalized f, rescale f_t=t^-1*f(t*x) and set t=h^N as in
QNI1. The classical position images become the identity, so (2) holds,
while the generic centralizer is still L[x]. Its module finiteness
over L[f_(h^N)] is equivalent to invertibility of f. Thus a universal
centralizer-finiteness assertion even in this planar identity-reduction
cotangent subclass would imply JC2. This is a precise obstacle to
deducing it merely from the Weyl relations. The rescaled full
identity examples in [QNI1](IDENTITY_LIFTS_AND_UNBOUNDED_INVERSES.md)
supply higher-rank instances where the relations, bounded forward
images and every desired finite identity contact hold but this
centralizer finiteness fails. This is not an equivalence between
JC2 and Long-specific finiteness, nor a cotangent normalization of
Long.

The finiteness criterion is sufficient for full completion and is not
necessary: QNI1's nonfinite-centralizer tuples are already flat. Thus
failure of (15), or even proven failure of centralizer module finiteness,
is not a nonlifting theorem. Conversely, a nonzero actual curvature
class under the theorem's hypotheses proves the centralizer is not
module-finite, including on the generic fiber.

Finite fraction-field degree is weaker still. The polynomial
inclusion K[x,x*y] subset K[x,y] has identical fraction fields but
is not module-finite: on x=0 its fiber is K[y]. This is a basic
non-Keller control; its Jacobian is x.

## Evidence and the remaining construction gate

The standalone checker verifies finite Taylor coefficient extraction,
a dual-number presentation with compatible connections and Fitting
minors, the escaping-support identities, exact rank-two PBW relations
and bounds for (17), the three-variable radial homotopy, the parity
basis reduction, and the nonfinite birational control. It validates
the retained QSC1 and QNI1 source hashes, without replaying their
universal proofs. All-rank exhaustion, the connection and etale
arguments, the full leading-symbol certificate and the Keller
boundary are written proofs. No independent theorem replay, formal
verification, external review or novelty priority is asserted.

The next actual input must consist of a finite commuting Long pair,
h-regular individual conjugates preserving the prescribed classical
tuple, and either centralizer module finiteness or a complete
leading-symbol certificate. QST1 already excludes the retained
G(a)+mu*v rational pair orbit. No such new Long input is supplied
here; bounded scalar lifting, a finite full Long carrier and
all-representative nonlifting remain open.

Run the declared exact controls from the repository root:

```sh
.venv/bin/python research/dixmier/dixmier-formal-lifting/check_finite_centralizer_completion.py
```
