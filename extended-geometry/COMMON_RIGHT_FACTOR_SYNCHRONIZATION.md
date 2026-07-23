# Common-right-factor synchronization

Work over a commutative ring \(S\).  In the characteristic-zero Hessian
application \(S\) is a \(\mathbb Q\)-algebra, and composition factors are
normalized to be monic and original.  This note isolates the elementary
rigidity statement behind transported Hessian synchronization.

## 1. The top-jet theorem

> **Theorem 1.1 (common-right-degree synchronization).**  Let
> \(N=mr\), with \(m,r>1\), and suppose
> \[
>  P_i=H_i\circ R_i\in S[x]\qquad(i=1,2),
> \]
> where \(H_i\) is monic of degree \(m\), and \(R_i\) is monic and
> original of degree \(r\).  Assume that \(m\) is a unit in \(S\).  If
> \[
>  \deg(P_1-P_2)<r,
> \]
> then \(R_1=R_2\), and \(P_1-P_2\) is constant.  If in addition
> \(P_1(0)=P_2(0)\), then \(P_1=P_2\).

**Proof.**  Write
\[
 R_i=x^r+\rho_{i,r-1}x^{r-1}+\cdots+\rho_{i,1}x.
\]
For \(1\leq j<r\), every nonleading term of \(H_i(R_i)\) has degree at
most
\[
 (m-1)r=N-r<N-j.
\]
Consequently
\[
 [x^{N-j}]P_i=[x^{N-j}]R_i^m
 =m\rho_{i,r-j}+
 \Phi_j(\rho_{i,r-j+1},\ldots,\rho_{i,r-1}),                 \tag{1.1}
\]
where \(\Phi_j\) is universal.  The two composites have equal
coefficients in degrees \(N-1,\ldots,N-r+1\).  Since \(m\) is a unit in
\(S\), descending induction in (1.1) gives
\[
 \rho_{1,r-j}=\rho_{2,r-j}\qquad(1\leq j<r).
\]
The original normalization gives equality of the constant coefficients,
so \(R_1=R_2=:R\).

Now \(P_1-P_2=(H_1-H_2)\circ R\).  If \(H_1-H_2\) is nonconstant, choose
its highest nonzero coefficient \(c\), say in degree \(q\geq1\).  Since
\(R\) is monic, the term \(cR^q\) has the uncancellable leading term
\(cx^{qr}\).  Thus
\[
 \deg\bigl((H_1-H_2)\circ R\bigr)=qr\geq r,
\]
contrary to the hypothesis.  Hence \(H_1-H_2\) is constant.  Equality at
zero makes that constant zero.  \(\square\)

The proof is triangular and remains valid with nilpotents and zero
divisors in \(S\), provided \(m\) is invertible.  It is the first \(r-1\)
steps of approximate-root reconstruction, but it needs neither a field nor
a Ritt classification.

## 2. Hessian consequence

Let \(X\) be an irreducible component of the intersection of two
exact-degree Hessian-composition cuts, and let \(K=k(X)\).  Their canonical
lifts, after passing to the monic chart of the projective incidence, have
the form
\[
 P_1=F+\lambda_1x,\qquad P_2=F+\lambda_2x.                   \tag{2.1}
\]
Suppose the two generic decompositions refine, after monic-original
normalization, to
\[
 P_1=C_1\circ R_1,\qquad P_2=C_2\circ R_2,\qquad
 \deg R_1=\deg R_2=r>1.                                     \tag{2.2}
\]
Theorem 1.1 applies because \(P_1-P_2=(\lambda_1-\lambda_2)x\) has
degree below \(r\).  It first gives
\[
 R_1=R_2=:R
\]
and then gives
\[
 \lambda_1=\lambda_2.                                       \tag{2.3}
\]

This proves the proposed common-right-factor statement on every
irreducible component under the precise hypothesis that both generic
lifts have terminal normalized factors of the same degree.  It is
important that “common right factor” in the hypothesis mean a common
terminal **degree**, not an already identical polynomial: equality of the
polynomials is a conclusion of the top-jet argument.

There is also a scheme-theoretic version.  If the refinements (2.2) exist
over the coordinate ring \(S\) of the possibly nonreduced Hessian
intersection, the same proof gives
\[
 \lambda_1-\lambda_2=0\quad\hbox{in }S.
\]
Equivalently, the synchronization defect belongs to the sum of the two
Hessian residual ideals.  No radicality hypothesis is needed.

> **Corollary 2.2 (incomparable-cut synchronization).**  Let the two outer
> cut degrees \(a,c\) be incomparable under divisibility.  On every
> irreducible component of their characteristic-zero projective
> Hessian-incidence intersection, the canonical lifts synchronize
> generically.  After synchronization, their coarse right factors possess
> one normalized common terminal factor of degree
> \[
>  \gcd(N/a,N/c)=\frac{N}{\operatorname{lcm}(a,c)}.
> \]

**Proof.**  At the generic point let \(B,E\) be the coarse right factors.
If \(\lambda_1-\lambda_2\ne0\), (2.1) gives \(x\in K[B,E]\).  The
Abhyankar--Moh epimorphism degree theorem, after extending \(K\) to an
algebraic closure, forces
\(\deg B\mid\deg E\) or \(\deg E\mid\deg B\), equivalently
\(c\mid a\) or \(a\mid c\), a contradiction.  Thus the lifts synchronize.
The characteristic-zero greatest-common-right-component theorem may now be
applied noncircularly to the single common polynomial.  Monic-original
normalization, or Theorem 1.1, makes its degree-\(\gcd(\deg B,\deg E)\)
terminal factor canonical.  \(\square\)

## 3. What remains nonformal

The theorem separates the transport problem into two statements:

1. **rigidity:** two lifted terminal factors of degree \(r>1\) coincide;
2. **existence:** the reduced gcd refinement extends over the Hessian
   intersection, including its primary thickening.

Theorem 1.1 proves the first statement over every \(\mathbb Q\)-algebra.
The second is not a consequence of ordinary Engström or Ritt theory:
those theorems begin with two decompositions of one already synchronized
polynomial.  Thus one must not use Engström refinement to construct
\(R_1,R_2\) before proving synchronization.

Corollary 2.2 resolves reduced components with incomparable outer cuts
without assuming a refinement in advance.  It synchronizes first and only
then invokes ordinary gcd refinement.

For nonreduced intersections, existence of the ring-valued refinements is
the genuine all-degree frontier.  A reduced normal form does not by itself
lift through nilpotents.

### 3.1 A finite-flat criterion for the primary upgrade

There is a short commutative-algebra route when the residual scheme is
controlled over its Ritt base.  Let \(B\) be a regular integral
\(\mathbb Q\)-algebra parametrizing a normalized tame Ritt component, put
\(K=\operatorname{Frac}(B)\), and let \(A\) be the Hessian residual algebra
on the formal neighborhood of that component.

> **Theorem 3.1 (finite-flat synchronization criterion).**  Assume:
>
> 1. \(A\) is a finite maximal Cohen--Macaulay \(B\)-module;
> 2. the generic fiber \(A_K=A\otimes_BK\) is geometrically reduced; and
> 3. the relation graph synchronizes the missing linear coefficients at
>    all geometric field-valued points.
>
> Then every lift difference vanishes in \(A\).  Equivalently, all lift
> differences belong to the summed Hessian residual ideal.

**Proof.**  Work locally on \(B\).  Auslander--Buchsbaum over the regular
local ring \(B\) makes the finite maximal Cohen--Macaulay module \(A\)
free.  Hence \(A\) is \(B\)-torsion-free.  Field-valued
synchronization makes every lift difference \(\delta\) vanish at every
geometric point of the reduced algebra \(A_K\), hence \(\delta=0\) in
\(A_K\).  Torsion-freeness makes \(A\to A_K\) injective, so
\(\delta=0\) already in \(A\).  \(\square\)

The familiar geometric hypotheses imply condition 1: if \(A\) is finite
over \(B\), Cohen--Macaulay and pure of dimension \(\dim B\), and every
minimal component dominates \(\operatorname{Spec}B\), miracle flatness
makes \(A\) maximal Cohen--Macaulay over \(B\).

This criterion replaces an order-by-order Krull-intersection argument by
miracle flatness.  Its reducedness hypothesis must concern the actual
geometric generic fiber of the possibly nonreduced residual algebra, not
only its reduction.
Otherwise \(B[\epsilon]/(\epsilon^2)\) is a Cohen--Macaulay countermodel:
its reduced parametrization is etale, while \(\epsilon\) vanishes on the
reduction but not scheme-theoretically.  Likewise, Cohen--Macaulayness only
at a minimal generic point is insufficient; the required hypothesis is on
the formal neighborhood under consideration.

Once the lift differences vanish,
\[
 I_{\rm full}=I_{\rm Hess}+(c_1-\lambda).
\]
Projection is therefore the graph isomorphism of one regular function.  For
the component-adic filtrations transported by this isomorphism, the Rees
algebras are isomorphic, which is the precise Rees-strict conclusion of
Theorem 3.1.  Strictness for an unrelated filtration is stronger and does
not follow merely from synchronization.

### 3.2 The weaker support-saturation criterion

Finite flatness is more than synchronization needs.  Let \(C\) be a
Noetherian residual algebra, let \(\mathfrak a\subset C\) be an ideal, and
let \(\delta\in C\) vanish on
\(\operatorname{Spec}C\setminus V(\mathfrak a)\).  Then
\[
 \delta\in H^0_{\mathfrak a}(C)=0:_C\mathfrak a^\infty.     \tag{3.1}
\]
Hence any one of the equivalent conditions
\[
\boxed{
\begin{aligned}
H^0_{\mathfrak a}(C)&=0,\\
\operatorname{Ass}(C)\cap V(\mathfrak a)&=\varnothing,\\
\operatorname{grade}(\mathfrak a,C)&\ge1
\end{aligned}}
                                                               \tag{3.2}
\]
forces \(\delta=0\).  If \(C=R/I\), the presentation form is simply
\[
 \boxed{I:\mathfrak a^\infty=I.}                            \tag{3.3}
\]
The same statement holds after flat adic completion; it is enough to
prove (3.3) before completing.

This is the
[support-saturation principle](../verified/SUPPORT_SATURATION_PRINCIPLE.md).
Theorem 3.1 is its finite-flat shortcut: local freeness over the integral
Ritt base supplies a regular base element in every nonzero support ideal.
Conversely, (3.3) can hold even when the residual algebra is not flat.

## 4. Degrees thirty and forty-two

For outer cuts \(a,c\) of degree \(N\), the expected terminal degree is
\[
 r=\gcd(N/a,N/c)=\frac{N}{\operatorname{lcm}(a,c)}.           \tag{4.1}
\]
The nontrivial incomparable patterns in the first two frontier degrees are:

| \(N\) | cuts | coarse right degrees | \(r\) | primitive core |
|---:|---:|---:|---:|---:|
| 30 | \(\{2,3\}\) | \(15,10\) | 5 | \(2\) versus \(3\), degree 6 |
| 30 | \(\{2,5\}\) | \(15,6\) | 3 | \(2\) versus \(5\), degree 10 |
| 30 | \(\{3,5\}\) | \(10,6\) | 2 | \(3\) versus \(5\), degree 15 |
| 42 | \(\{2,3\}\) | \(21,14\) | 7 | \(2\) versus \(3\), degree 6 |
| 42 | \(\{2,7\}\) | \(21,6\) | 3 | \(2\) versus \(7\), degree 14 |
| 42 | \(\{3,7\}\) | \(14,6\) | 2 | \(3\) versus \(7\), degree 21 |

Thus the degree-thirty coincidence has a structural explanation, and the
same three decorations recur in degree forty-two with \(5\) replaced by
\(7\).  On reduced components all six rows synchronize by the
incomparable-cut argument.  On any chart where the terminal refinement
lifts scheme-theoretically, Theorem 1.1 upgrades this to exact
scheme-theoretic synchronization.

The working tree also contains three transported characteristic-zero
ideal-membership checkers:
[`verify_degree30_transported_23_synchronization.py`](../scripts/verify_degree30_transported_23_synchronization.py),
[`verify_degree30_transported_25_synchronization.py`](../scripts/verify_degree30_transported_25_synchronization.py),
and
[`verify_degree30_transported_35_synchronization.py`](../scripts/verify_degree30_transported_35_synchronization.py).
They address the primary upgrade separately; Theorem 1.1 does not replace
the refinement-lifting part of those calculations.

## 5. Counterexample search is now sharply localized

A characteristic-zero search for two normalized terminal factors of the
same degree \(r>1\) which generate distinct intermediate fields cannot
succeed: their top \(r-1\) coefficients already force the factors to be
equal.  A meaningful counterexample search must instead target one of:

- a reduced Hessian component on which the predicted terminal refinement
  does not exist;
- a nonreduced deformation on which the reduced refinement fails to lift;
- positive characteristic dividing the outer degree, where the triangular
  coefficient \(m\) in (1.1) is not invertible.

This changes the computational target from “compare two reconstructed
right factors” to “test formal liftability of the common refinement.”

The last exception is real.  Over
\[
 S=\mathbb F_2[\epsilon]/(\epsilon^2)
\]
take
\[
 H=z^2+z,\qquad R_1=x^2+\epsilon x,\qquad R_2=x^2.
\]
Then
\[
 H(R_1)-H(R_2)=\epsilon x.
\]
Thus equal-degree terminal refinements need not synchronize when the outer
degree is not invertible.  This defect is invisible on the reduced finite
field point, so a search only over fields cannot test the primary
characteristic-\(p\) obstruction.

## 6. The degree-forty-two computation to run

The first degree-forty-two primary calculation is now exact.  For the
\(\{2,7\}\) pair, transport the primitive degree-fourteen power collision
through a generic cubic.  The common cubic has two coefficients, the power
core has three coefficients and one translation parameter, and triangular
recovery leaves
\[
 5\text{ normal coordinates}\mid6\text{ base coordinates}.
\]
Let \(I\) be the nineteen pulled-back Hessian residuals, let
\(\mathfrak m\) be the ideal of the five normal coordinates, and let
\(\delta=\lambda_2-\lambda_7\).

> **Proposition 6.1 (degree-forty-two fourth normal neighborhood).**
> Over \(\mathbb Q\),
> \[
>  \delta\in I+\mathfrak m^5.
> \]
> Hence the transported \(\{2,7\}\) power component is
> scheme-theoretically synchronized through normal order four, uniformly
> over its six-dimensional base.

The exact block-order Gröbner basis has size \(88\).
[`verify_degree42_transported_27_normal_jets.py`](../scripts/verify_degree42_transported_27_normal_jets.py)
reconstructs the chart, verifies its polynomial inverse and the primitive
Ritt identity, adjoins every degree-five normal monomial, and reduces the
defect to zero over \(\mathbb Q\).  The same calculation over
\(\mathbb F_{32003}\) has the same basis size.  Exact normal orders one,
two, and three have basis sizes \(8,19,55\); they also follow formally from
the order-four membership.

The untruncated modular full-basis calculation timed out after \(900\)
seconds, and exact normal order five timed out after \(300\) seconds.  These
are performance boundaries, not failed reductions.  Thus the current
\(\{2,7\}\) characteristic-zero gap begins at order five.

A good-prime order-five calculation gives
\[
 \delta\in I+\mathfrak m^6
 \qquad\text{over }\mathbb F_{32003}.
\]
The truncated basis has size \(179\).  Hence the degree-five normal symbol
vanishes uniformly over the six-dimensional base in this fiber.  This is
strong evidence, but it is not by itself a characteristic-zero membership
certificate.  The obstruction class belongs to
\[
 \frac{I+\mathfrak m^5}{I+\mathfrak m^6},
\]
not the quotient with numerator and denominator reversed.

On the unresolved divisor, the rational specialization
\[
 (e_1,e_2,t,w_0,w_1,w_2)=(1,2,3,0,5,6)
\]
passes the stronger untruncated characteristic-zero membership test: the
full basis has size \(8\) and reduces \(\delta\) to zero.  Thus no primary
obstruction is present at this boundary point.  This specialization does
not prove generic membership on \(w_0=0\); the corresponding calculation
over \(\mathbb Q(e_1,e_2,t,w_1,w_2)\) still exceeds the current timeout.
The checker option `--w0-zero` performs exactly that remaining
function-field test, while `--base-values` supports exact rational boundary
probes.

### 6.1 Conormal and Rees interpretation

The normal-jet calculation has a sharper conceptual form.  Put
\[
 A=\mathbb Q[e_1,e_2,t,w_0,w_1,w_2],\qquad
 R=A[x_1,\ldots,x_5],\qquad \mathfrak m=(x_1,\ldots,x_5),
\]
and retain \(I\) and \(\delta\) from Proposition 6.1.  The zero section is
the transported power component, so \(I\subset\mathfrak m\) and
\(\delta\in\mathfrak m\).  Taking normal-linear terms gives the residual
conormal map
\[
 \kappa:(I+\mathfrak m^2)/\mathfrak m^2
       \longrightarrow \mathfrak m/\mathfrak m^2.           \tag{6.2}
\]
It is represented by the \(19\)-by-\(5\) normal Jacobian of the residuals
on the zero section.

> **Proposition 6.2 (dense all-order Rees synchronization).**  The ideal
> of maximal minors of \(\kappa\) is
> \[
>  \operatorname{Fitt}_0(\operatorname{coker}\kappa)=(w_0^2).
>                                                               \tag{6.3}
> \]
> Consequently, over \(D(w_0)\), the completed residual ideal equals the
> completed normal ideal:
> \[
>  I\widehat R=\mathfrak m\widehat R.
>                                                               \tag{6.4}
> \]
> In particular \(\delta=0\) in the formal completion along the entire
> \(D(w_0)\) part of the six-parameter component.  Equivalently, every
> positive Rees symbol of the synchronization defect vanishes there.

**Proof.**  Exact maximal-minor reduction gives (6.3).  After inverting
\(w_0\), (6.2) is onto, hence
\[
 \mathfrak m=I+\mathfrak m^2.
\]
Modulo \(I\), the finitely generated ideal
\(J=\mathfrak m(R/I)\) satisfies \(J=J^2\).  In the
\(J\)-adic completion, complete Nakayama gives \(J\widehat{(R/I)}=0\),
which is (6.4).  Since \(\delta\in\mathfrak m\), its completed class is
zero.  This kills all graded pieces of the \(\mathfrak m\)-adic Rees
filtration at once.  \(\square\)

This is the conceptual replacement for further jet expansion on the dense
chart: conormal surjectivity implies Rees strictness, and Rees strictness
implies synchronization at every order.  It also identifies the exact
remaining locus.  Modulo \(w_0\), all \(4\)-by-\(4\) conormal minors vanish
while the \(3\)-by-\(3\) minors generate the unit ideal.  Thus
\(\kappa\) has rank exactly \(3\) everywhere on \(w_0=0\); two normal
directions cease to be controlled at first order.  Proposition 6.1 still
kills the defect through Rees degree four on that divisor.  First-order
conormal/Nakayama alone cannot promote those four jets to all orders, so
one must inspect the quadratic normal cone.

Theorem 3.1 gives a second way to close this divisor without further normal
jets: prove that the completed residual algebra is finite
Cohen--Macaulay over the full six-dimensional Ritt base.  Its generic fiber
is already the reduced zero section over \(D(w_0)\), so finite flatness would
make the generic vanishing of \(\delta\) inject into the whole algebra,
including the ramified fiber over \(w_0=0\).

[`verify_degree42_conormal_rees_synchronization.py`](../scripts/verify_degree42_conormal_rees_synchronization.py)
checks \(I,\delta\subset\mathfrak m\), the Fitting identity (6.3), and the
rank-three degeneration on \(w_0=0\).

Set
\[
 \Delta=(t+e_1e_2)^2-4e_1^3.                               \tag{6.5}
\]
On \(w_0=0\), residuals \(5,11,17\) have a unit-triangular
normal Jacobian in \(x_3,x_4,x_5\), with diagonal
\((-1,-1,-1)\).  Formal implicit elimination leaves only \(x_1,x_2\).
The last two residuals have zero differential, and their quadratic forms
on this tangent plane are
\[
\begin{aligned}
 q={}&w_1\left[
 -\frac92e_2x_1^2+
 \left(9e_2^2+\frac32e_1\right)x_1x_2+
 \left(-\frac92e_2^3+\frac{15}{4}e_1e_2+\frac34t\right)x_2^2
 \right],\\
 r={}&w_1\left[
 -\frac34x_1^2+\frac32e_2x_1x_2+
 \left(-\frac34e_2^2+\frac34e_1\right)x_2^2
 \right].
                                                               \tag{6.6}
\end{aligned}
\]
Their binary resultant is
\[
 \operatorname{Res}(q,r)=\frac{81}{256}w_1^4\Delta.           \tag{6.7}
\]

> **Proposition 6.3 (quadratic-cone synchronization).**  On
> \(V(w_0)\cap D(w_1\Delta)\), the transported \(\{2,7\}\) power
> component is formally synchronized to every normal order.

**Proof.**  The unit pivots remove three normal variables without changing
the completed quotient.  By (6.7), \(q,r\) are a regular sequence of
binary quadrics.  The residual associated graded is therefore a quotient
of their complete intersection, whose Hilbert vector is \((1,2,1)\).
If \(J\) is the completed normal ideal, then \(J^3/J^4=0\).
Complete Nakayama applied to \(J^3\) gives \(J^3=0\).
Proposition 6.1 gives \(\delta\in I+\mathfrak m^5\), hence
\(\delta=0\) in the completed quotient.  \(\square\)

The unresolved all-order locus is thus confined to
\[
 V(w_0,w_1)\ \cup\ V(w_0,\Delta).                            \tag{6.8}
\]
On the first branch the quadratic Kuranishi map vanishes; on the second
the quadrics acquire a common tangent direction.  These are distinct
degeneration mechanisms and should be treated separately.

[`verify_degree42_divisor_rees_reduction.py`](../scripts/verify_degree42_divisor_rees_reduction.py)
certifies the unit-triangular elimination, the quadratic forms, and the
resultant identity.

Both branches in (6.8) admit one more intrinsic initial-form step.
First normalize the discriminant branch by
\[
 e_1=a^2,\qquad t=2a^3-a^2e_2.
\]
On this normalization the quadrics factor as
\[
\begin{aligned}
q&=-\frac32\bigl(x_1+(a-e_2)x_2\bigr)
 \bigl(3e_2x_1-(a^2+3ae_2+3e_2^2)x_2\bigr),\\
r&=-\frac34\bigl(x_1+(a-e_2)x_2\bigr)
 \bigl(x_1-(a+e_2)x_2\bigr).                                \tag{6.9}
\end{aligned}
\]
The determinant of the two complementary linear forms is \(a^2\).
Along the escaping tangent
\(\ell=x_1+(a-e_2)x_2=0\), implicit differentiation of the three pivot
residuals through order two gives terminal cubic coefficients
\[
\begin{aligned}
c_{19}&=-\frac58w_2a^2(e_2-2a)^2,\\
c_{18}&=(7e_2-a)c_{19}.                                     \tag{6.10}
\end{aligned}
\]

> **Proposition 6.4 (common-tangent cubic closure).**  On
> \[
> V(w_0,\Delta)\cap D(w_1w_2t),
> \]
> the transported component is formally synchronized to every order.

**Proof.**  On \(\Delta=0\), inverting \(t\) also inverts \(e_1=a^2\),
and \(e_2-2a=-t/e_1\).  Thus the complementary quadratic factors span
the tangent ideal and (6.10) is nonzero.  In coordinates \((\ell,s)\)
the initial Kuranishi ideal contains
\[
 (\ell^2,\ell s,s^3).
\]
Its maximal ideal has cube zero, so Proposition 6.1 again kills
\(\delta\) exactly.  \(\square\)

The apparent \(t=0\) exception in Proposition 6.4 is removed one order
later.  On the discriminant normalization, \(t=0\) and \(e_1\ne0\)
give \(e_2=2a\).  The common-tangent cubic vanishes, but the two terminal
quartic coefficients are
\[
 c_{18}^{(4)}=\frac{15}{16}w_2a,\qquad
 c_{19}^{(4)}=\frac{5}{64}w_2.
\]
At the cusp \(e_1=t=0\), where the quadratic ideal collapses to the
single square \(\ell^2\), they specialize to
\[
 c_{18}^{(4)}=\frac58e_2w_2,\qquad
 c_{19}^{(4)}=\frac{5}{64}w_2.
\]

> **Proposition 6.5 (discriminant quartic closure).**  The entire stratum
> \[
> V(w_0,\Delta)\cap D(w_1w_2)
> \]
> is formally synchronized to every normal order.

**Proof.**  Proposition 6.4 treats \(t\ne0\).  When \(t=0\) and
\(e_1\ne0\), the initial ideal contains
\((\ell^2,\ell s,s^4)\), whose maximal ideal has fourth power zero.
At \(e_1=0\), the specialized initial ideal contains
\((\ell^2,s^4)\), whose maximal ideal has fifth power zero.  The
coefficient \(5w_2/64\) is a unit at every point under consideration,
so the latter statement lifts over the local base, including its
nilpotent discriminant thickening.  Proposition 6.1 then makes the
synchronization defect zero.  \(\square\)

For the other branch define
\[
\begin{aligned}
A={}&4e_1^3-e_1^2e_2^2+e_2^3t-6e_1e_2t,\\
B={}&e_1^4e_2^4-2e_1^2e_2^5t+e_2^6t^2-8e_1^5e_2^2
{}+20e_1^3e_2^3t-12e_1e_2^4t^2\\
&+16e_1^6-48e_1^4e_2t+27e_1^2e_2^2t^2+9e_2^3t^3
{}+36e_1^3t^2-54e_1e_2t^3+27t^4.                          \tag{6.11}
\end{aligned}
\]
When \(w_1=0\), the quadratic Kuranishi map vanishes.  After the same
unit-pivot elimination, the two terminal binary cubics \(p_3,q_3\)
satisfy
\[
 \operatorname{Res}(p_3,q_3)
 =-\frac{15625}{262144}w_2^6AB.                             \tag{6.12}
\]

> **Proposition 6.6 (cubic-cone synchronization).**  On
> \[
> V(w_0,w_1)\cap D(w_2AB),
> \]
> the transported component is formally synchronized to every order.

**Proof.**  Equation (6.12) makes \(p_3,q_3\) a regular sequence of
binary cubics.  Their complete-intersection Hilbert vector is
\((1,2,3,2,1)\), so the fifth power of the completed normal ideal is
zero.  The membership \(\delta\in I+\mathfrak m^5\) from Proposition
6.1 is therefore exact.  \(\square\)

Consequently the remaining all-order problem is supported on
\[
 V(w_0,w_2)\ \cup\ V(w_0,w_1,AB).                           \tag{6.13}
\]
This replaces the two full Kuranishi branches by one coefficient divisor
and the two cubic-resultant divisors inside \(w_1=0\).

There is now one global algebraic target.  Let
\[
 S=\mathbb Q[e_1,e_2,t,w_0,w_1,w_2],\qquad
 C=R/I,\qquad
 \widehat C=\widehat{C}^{\,\mathfrak m},
\]
where \(\mathfrak m=(x_1,\ldots,x_5)\), and put
\[
 \mathfrak k=(w_0,\;w_1w_2,\;ABw_2)\subset S.               \tag{6.14}
\]
Then
\[
 V(\mathfrak k)
 =V(w_0,w_2)\cup V(w_0,w_1,AB).
\]
The preceding propositions say exactly that the completed defect is
supported there:
\[
 \delta\in
 H^0_{\mathfrak k\widehat C}(\widehat C)
 =0:_{\widehat C}\mathfrak k^\infty.                        \tag{6.15}
\]
Therefore the strongest efficient target is not another normal jet or
full flatness, but
\[
 \boxed{I:\mathfrak k^\infty=I.}                            \tag{6.16}
\]
By the support-saturation principle and flatness of completion, (6.16)
implies \(\delta=0\) on the entire completed residual algebra.
Equivalently, it is enough to prove that no associated prime of \(C\)
contains \(\mathfrak k\).  Operationally one should compute the saturation,
identify any associated prime of its quotient if it is nonzero, and only
then return to higher Rees symbols.

On the discriminant normalization their overlap factors further:
\[
\begin{aligned}
A|_{\Delta}
 &=a^2(2a-e_2)^2(a^2-2ae_2-e_2^2),\\
B|_{\Delta}
 &=a^4(2a-e_2)^4
 (37a^4-22a^3e_2-7a^2e_2^2+4ae_2^3+e_2^4).
                                                               \tag{6.17}
\end{aligned}
\]
Thus the next calculations are localized to \(w_2=0\) or to one of the
two displayed residual factors; no point of the discriminant branch with
\(w_1w_2\ne0\), and no generic point of the cubic branch, remains.

There is also a useful generic subresultant refinement on the cubic branch.
Write the tangent coordinates as \((u,v)\).  At the generic point of
\(A=0\), the two terminal cubics share
\[
 L_A=e_2u+(2e_1-e_2^2)v.                                   \tag{6.18}
\]
At the generic point of \(B=0\), they share
\(L_B=\alpha_Bu+\beta_Bv\), where
\[
\begin{aligned}
\alpha_B={}&8e_1^5-2e_1^4e_2^2-20e_1^3e_2t
 +4e_1^2e_2^3t+6e_1^2t^2+12e_1e_2^2t^2
 -2e_2^4t^2-9e_2t^3,\\
\beta_B={}&-4e_1^5e_2+e_1^4e_2^3+4e_1^4t
 +9e_1^3e_2^2t-2e_1^2e_2^4t-9e_1^2e_2t^2\\
&\hspace{2.9em}-5e_1e_2^3t^2+e_2^5t^2+6e_2^2t^3.
                                                               \tag{6.19}
\end{aligned}
\]
These are the linear first-subresultant factors.  The certificate works on
\(D(e_2(6e_1-e_2^2))\subset V(A)\) and
\(D(\alpha_B)\subset V(B)\); these opens contain the respective generic
points.  Since \(A\) and \(B\) are irreducible and occur to the first power
in (6.12), either generic point has exactly one common projective root.
The next calculation on each divisor is therefore again one-dimensional:
restrict the quartic Kuranishi terms to \(L_A=0\) or \(L_B=0\).  The
component \(w_2=0\), where both cubics vanish, remains the genuinely
deeper sevenfold power collision.

[`verify_degree42_kuranishi_branches.py`](../scripts/verify_degree42_kuranishi_branches.py)
certifies the implicit cubic coefficients (6.10), the binary-cubic
resultant (6.12), and the common factors (6.18)--(6.19).
[`verify_degree42_discriminant_quartics.py`](../scripts/verify_degree42_discriminant_quartics.py)
certifies both terminal quartic formulas and the discriminant closure.

For the remaining calculation, the theorem removes the right factor from
the elimination variables.  For each degree-forty-two row:

1. reconstruct the unique candidate \(R\) directly from the top \(r-1\)
   coefficients;
2. expand both coarse right factors in the free
   \(S[R]\)-module basis \(1,x,\ldots,x^{r-1}\);
3. use the nonconstant basis coefficients as refinement-normal
   coordinates;
4. test whether the Hessian residual ideal kills those coordinates, or
   directly kills the linear lift defect.

This targets the only missing statement—the scheme-valued existence of the
refinement—without carrying two redundant copies of \(R\) through a
Gröbner calculation.  Modular searches should use good characteristics not
dividing \(42/r\), and should include dual-number or higher Artin probes if
the goal is to detect a primary defect.

## References

- A. Bodin, [*Decomposition of polynomials and approximate
  roots*](https://arxiv.org/abs/0910.1676).
- R. M. Beals, J. L. Wetherell, and M. E. Zieve,
  [*Polynomials with a common
  composite*](https://dept.math.lsa.umich.edu/~zieve/papers/common_composites.pdf),
  especially Theorem 4.5.
- P. Müller and M. E. Zieve, [*On Ritt's polynomial decomposition
  theorems*](https://arxiv.org/abs/0807.3578).
- K. Ziegler, [*Tame Decompositions and
  Collisions*](https://arxiv.org/abs/1402.5945).
