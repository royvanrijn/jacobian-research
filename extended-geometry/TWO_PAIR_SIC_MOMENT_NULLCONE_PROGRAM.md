# The falsified all-degree two-pair moment--nullcone program

## 1. Status and objective

The all-degree conjecture proposed in this note is false in every degree
\(d\geq4\).  The
[bidegree-\((4,4)\) two-pair counterexample](TWO_PAIR_IMAGE_MATHIEU_COUNTEREXAMPLE.md)
gives
\[
 \mathcal E_2(F^m)=0,\qquad
 \mathcal E_2(ZF^m)=\frac{(4m+2)!\,m!}{(2m+1)!!}\ne0
 \quad(m\geq1).
\]
Thus its moment-zero point is not in the one-sided nullcone and
\(\mathrm{MN}_4\) fails.  More generally, the invariant-power family
\[
 F_d=R^{d-4}F
\]
satisfies
\[
 \mathcal E_2(F_d^m)=0,\qquad
 \mathcal E_2(ZF_d^m)=\frac{(dm+2)!\,m!}{(2m+1)!!}\ne0
\]
for every \(d\geq4\) and \(m\geq1\).  Hence \(\mathrm{MN}_d\) fails for
every \(d\geq4\), with no parity or congruence restriction.  In particular,
unrestricted
\(\operatorname{SIC}(2)\) is false.

Combining ordinary powers with invariant multiplication gives the stronger
family
\[
 G_{r,k}=R^kF^r,\qquad d=4r+k,
\]
with
\[
 \mathcal E_2(G_{r,k}^m)=0,\qquad
 \mathcal E_2(ZG_{r,k}^m)
 =\frac{(dm+2)!\,(rm)!}{(2rm+1)!!}\ne0.
\]
Choosing \(r=\lfloor d/4\rfloor\) bounds the exact \(R\)-adic order by
three.  In degrees divisible by four, \(F^{d/4}\) is \(R\)-primitive.
The endpoint profiles
\(\rho_h(z)=(1-z)(1+z)^{h-1}\) further give explicit \(R\)-primitive
witnesses \(\Phi_h\in V_{4h}\) that are not proper powers, detected by the
positive integrals
\[
 \int_0^1(1-v^2)^m(1+v^2)^{(h-1)m}\,dv.
\]
Every polynomial in the classified one-profile Hopf class has minimal
balanced degree divisible by four, so that mechanism cannot produce
primitive witnesses in the other congruence classes.
Since the invariant ring on one vector and one covector is \(k[R]\), the
invariant-multiplier strategy cannot remove the residual radial factor in
the other three congruence classes.  Their \(R\)-primitive status is open.

The first degree-five multiplicative gate is nevertheless closed.  For
the quartic seed \(F\), an arbitrary bilinear multiplier
\[
 L=aR+bZ+cW+eT
\]
satisfies the first four pure identities for \(LF\) only when \(L=aR\).
The exact elimination and coprime residual polynomials are in the
[degree-five multiplier obstruction](TWO_PAIR_DEGREE_FIVE_MULTIPLIER_OBSTRUCTION.md).
Thus a primitive point of \(V_5\), if one exists, cannot be obtained by
simply multiplying the known seed by a noninvariant bilinear.

Finite prefixes require separate care.  The primitive family
\[
 G_{d,\lambda}=R^{d-4}F+\lambda Z^d
\]
has its first \(d\) pure moments zero but an explicitly nonzero moment of
order \(d+1\) for \(\lambda\ne0\).  See the
[primitive-prefix obstruction](TWO_PAIR_PRIMITIVE_PREFIX_OBSTRUCTION.md).
Consequently no proposed degree-\(d\) primitive classification can rely
only on moments through order \(d\).
The same phase argument excludes every nonzero correction
\(\sum_{j\geq1}c_jR^{d-j}Z^j\) from all-order propagation: its least
positive phase \(s\) is detected at moment \(s+1\).
The
[opposite-monomial obstruction](TWO_PAIR_OPPOSITE_MONOMIAL_OBSTRUCTION.md)
excludes every pair \(Z^sT^{d-s},W^sT^{d-s}\) for every \(d>4\).
Odd heights and even-height phases at least three have direct all-degree
phase certificates.  The two low even-height phases reduce to exact
height-polynomial eliminations with no roots modulo \(47\) and \(29\).
Thus the one-pair two-sided monomial ansatz is closed; a continuation must
mix several absolute phases or use a nonmonomial height profile.
In the first such multi-pair sector, the four degree-five odd-height
monomials have projectively empty quadratic obstruction cone: moments
through order ten generate an ideal containing the cube of every
parameter.  The quotient has Hilbert vector \((1,4,1)\), length six, and
nondegenerate socle pairing.  After adjoining all six even-height
monomials, moments \(2,\ldots,7\) eliminate their second-order
coefficients and projected moments \(8,\ldots,11\) form a length-sixteen
quadratic complete intersection.  Thus \(RF\) is formally isolated in
the full degree-five monomial correction space, although global finite
multiweight solutions remain open.
The same exact Hopf-angular elimination proves formal isolation in the
full monomial correction spaces for degrees six and seven.  Their
projected odd-height systems are six-quadric complete intersections of
length \(64\).  Good-prime Singular certificates prove degrees eight
through eleven: the first two have eight-quadric complete intersections
of length \(256\), while degrees ten and eleven have ten-quadric complete
intersections of length \(1024\).  The verified degrees five through
eleven suggest a uniform
parity theorem with \(2\lceil d/2\rceil\) linear pivots followed by a
\(2\lfloor d/2\rfloor\)-quadric complete intersection; proving its
projected determinant without Gröbner expansion is now the main
monomial-slice continuation.
Odd height delays detection further: the primitive family
\(R^{d-4}F+\lambda Z^{d-1}T\) survives through moment \(2d-1\) and fails
explicitly at moment \(2d\).  Thus a consecutive primitive cutoff must
reach order at least \(2d\).

The representation theory, nullcone incidence geometry, Hilbert-series
tests, and exact bidegree-\((2,2)\) and \((3,3)\) calculations below remain
valid at their stated scopes.  They now serve a degree-by-degree
classification program: \(\mathrm{MN}_2\) is true, while
\(\mathrm{MN}_d\) is false for every \(d\geq3\).  In degree three the
[full-rank Rodrigues survivor](TWO_PAIR_SIC_BIDEGREE33_RODRIGUES_SURVIVOR.md)
has all pure moments zero and determinant one.  Unlike the degree-four
witness, however, it is SIC-safe: every fixed multiplier \(Q\) vanishes
against its powers once
\(m>3\deg_{Z,Y}Q\).

There is a sharper operational split.  The
[split-symbol theorem](SPLIT_SYMBOL_GVC_THEOREM.md) proves the complete
rank-one Segre cone SIC-safe in every degree.  Consequently rank-one
root partitions, generic binary-symbol orbit closure, and finite
moment--nullcone certificates below are no longer gates for SIC
counterexample discovery.  They remain relevant to the strictly stronger
scheme-theoretic question \(\mathrm{MN}_d\).  The active SIC search is the
[rank-stratified programme](RANK_STRATIFIED_MOMENT_PROGRAM.md): begin at
exact coefficient rank two, compute exact determinantal components, and
extract all-order recurrences on rank-two factorizations.

For \(d\geq 1\), let \(V_d\) be the space of two-pair forms of bidegree
\((d,d)\), and put

\[
 \mu_m(f)=\mathcal E_2(f^m),\qquad m\geq1.
 \tag{1.1}
\]

The central target is:

> **Falsified moment--nullcone conjecture \(\mathrm{MN}_d\).** The common
> zero set of all \(\mu_m\) on \(V_d\) is the pair-linear one-sided
> nullcone.

The [complete bidegree-\((2,2)\) theorem](TWO_PAIR_SIC_BIDEGREE22_FRONTIER.md)
proves \(\mathrm{MN}_2\).  The Rodrigues survivor disproves
\(\mathrm{MN}_3\), and the propagated degree-four witness disproves every
\(\mathrm{MN}_d\) for \(d\geq4\).  The structural route still organizes
the exact strata and is useful for the logically weaker SIC question, but
moment--nullcone equality itself is no longer the degree-three target.

## 2. The balanced representation

Use contraction pairs \((W,Z),(V,Y)\). The diagonal
\(\mathrm{SL}_2\)-action preserving \(WZ+VY\) identifies

\[
 V_d
 \cong \operatorname{Sym}^d(\mathbb C^2)^*
       \otimes\operatorname{Sym}^d(\mathbb C^2)
 \cong \operatorname{End}(\operatorname{Sym}^d)
 \cong \bigoplus_{r=0}^{d}\operatorname{Sym}^{2r}.
 \tag{2.1}
\]

The last equality is the Clebsch--Gordan formula: the
\(\mathrm{SL}_2\)-module \(\operatorname{Sym}^d\) is self-dual and

\[
 \operatorname{Sym}^d\otimes\operatorname{Sym}^d
 \cong
 \operatorname{Sym}^{2d}\oplus
 \operatorname{Sym}^{2d-2}\oplus\cdots\oplus
 \operatorname{Sym}^{0}.
 \tag{2.2}
\]

Thus a form has irreducible coordinates

\[
 f=(F_0,F_2,F_4,\ldots,F_{2d}),
 \tag{2.3}
\]

where \(F_{2r}\) is a binary form of degree \(2r\). Every moment
\(\mu_m\) is an \(\mathrm{SL}_2\)-invariant homogeneous polynomial of
degree \(m\). The first moment is a nonzero scalar multiple of \(F_0\), so
moment vanishing removes the scalar summand.

The standard Hilbert--Mumford criterion gives an explicit description of
the nullcone \(N_d\). For direct sums of binary forms, the criterion is
recorded explicitly by Brouwer and Popoviciu in
[*Sylvester versus Gundelfinger*](https://sigma-journal.com/2012/075/):

\[
 f\in N_d
 \quad\Longleftrightarrow\quad
 F_0=0
 \ \text{and there is a common linear form }L
 \text{ such that }L^{r+1}\mid F_{2r}
 \text{ for every }r\geq1.
 \tag{2.4}
\]

Zero components impose no condition. In binary-form language, every
nonzero \(F_{2r}\) has a root of multiplicity strictly greater than half
its degree, and all these roots agree. Such a root is unique.

Equation (2.4) is exactly the pair-linear one-sided condition. Hence
\(\mathrm{MN}_d\) would imply \(\operatorname{SIC}(2)\) for every
bidegree-\((d,d)\) form. Unequal bihomogeneous bidegrees already have
nonzero total dual-minus-coordinate weight and are one-sided. Proving
\(\mathrm{MN}_d\) for all \(d\) would therefore settle the bihomogeneous
two-pair problem.

### Geometry of the target

The common-root description also gives the nullcone dimension without
elimination. Fixing \([L]\in\mathbb P^1\), the component \(F_{2r}\) has
the form

\[
 F_{2r}=L^{r+1}Q_{r-1},
 \qquad Q_{r-1}\in\operatorname{Sym}^{r-1}(\mathbb C^2),
 \tag{2.5}
\]

so its fiber has dimension \(r\). The resulting incidence space is a
vector bundle of rank

\[
 \sum_{r=1}^{d}r=\frac{d(d+1)}2
 \tag{2.6}
\]

over \(\mathbb P^1\). Away from the zero form, the destabilizing root of
the first nonzero component is unique, so the incidence map is generically
one-to-one. It follows that \(N_d\) is irreducible and

\[
 \boxed{\dim N_d=1+\frac{d(d+1)}2},\qquad
 \boxed{\operatorname{codim}_{V_d}N_d=\frac{d(d+3)}2}.
 \tag{2.7}
\]

For \(d=2\) this gives dimension \(4\); for \(d=3\), dimension \(7\).
These agree with the independent exact eliminations in the two frontier
notes.

### Why the nullcone gives the Mathieu conclusion

This implication is worth separating from the conjectural converse.
After a pair-linear change, a nullcone point is supported in the strict
one-sided positions

\[
 M_{ij}=W^{d-i}V^iZ^{d-j}Y^j,\qquad i>j.
 \tag{2.8}
\]

Every monomial in \(f^m\) therefore has total \(V\)-exponent minus total
\(Y\)-exponent at least \(m\). A fixed multiplier \(g\) changes this
difference by a bounded amount. For all sufficiently large \(m\), every
monomial of \(gf^m\) has more \(V\)-derivatives than \(Y\)-degree, and its
contraction is zero. Thus

\[
 f\in N_d
 \quad\Longrightarrow\quad
 \mathcal E_2(gf^m)=0
 \quad\text{for every fixed }g\text{ and all }m\gg0.
 \tag{2.9}
\]

Consequently moment--nullcone equality is sufficient for the Special
Image conclusion on the balanced stratum; it is not merely an invariant
theory reformulation of the pure moments.

## 3. Invariant-ring formulation

Let

\[
 S_d=\mathbb Q[V_d],\qquad
 R_d=S_d^{\mathrm{SL}_2},\qquad
 M_d=(\mu_1,\mu_2,\ldots)\subset R_d.
 \tag{3.1}
\]

The nullcone is the zero set of the positive-degree invariants
\((R_d)_+\). Consequently the set-theoretic content of
\(\mathrm{MN}_d\) is

\[
 \sqrt{M_d}=(R_d)_+,
 \tag{3.2}
\]

or, equivalently after extending to \(S_d\),

\[
 \sqrt{S_dM_d}=I(N_d).
 \tag{3.3}
\]

This suggests two scalable proof mechanisms.

1. Find finitely many moment orders and power certificates showing that a
   homogeneous generating set of \((R_d)_+\) lies in \(\sqrt{M_d}\).
2. Prove that \(R_d\) is integral over a finite moment subalgebra and that
   the fiber over the moment origin consists only of the invariant origin.

The exact radical certificates in bidegrees \((2,2)\) and on the pure
binary-sextic slice are finite instances of the first mechanism. No
uniform finite cutoff or integrality theorem is currently proved.

There are two useful qualifications.

First, \(R_d\) is finitely generated and therefore Noetherian. The ideal
generated by the infinite moment sequence is consequently generated by a
finite subcollection:

\[
 M_d=(\mu_{m_1},\ldots,\mu_{m_s})
 \quad\text{for some finite set of orders depending on }d.
 \tag{3.4}
\]

Thus a finite cutoff exists for each fixed \(d\) without any conjecture.
What is missing is an effective choice of the orders, a useful bound, and
proof that their radical is \((R_d)_+\).

Second, Krull dimension gives a sharp lower bound on how many moment
equations a global proof can use. For \(d\geq2\), a generic binary form in
the highest \(\operatorname{Sym}^{2d}\) summand has finite
\(\mathrm{SL}_2\)-stabilizer. Hence

\[
 \dim R_d=\dim V_d-\dim\mathrm{SL}_2=(d+1)^2-3.
 \tag{3.5}
\]

If \(s\) homogeneous moments define the nullcone, their ideal has radical
\((R_d)_+\), whose height is \(\dim R_d\). Krull's height theorem therefore
forces

\[
 \boxed{s\geq(d+1)^2-3.}
 \tag{3.6}
\]

This bound explains two computations already in the repository:

- for \(d=2\), the lower bound is six, and the first six moments attain it;
- for \(d=3\), any full-space moment--nullcone proof needs at least
  thirteen moment equations, although fewer can suffice on a proper
  irreducible slice.

If exactly \(\dim R_d\) moments define the nullcone, Hilbert's criterion
makes them a homogeneous system of parameters. This is the most economical
possible invariant-theoretic certificate. The binary-form version of this
criterion is described by Brouwer, Draisma, and Popoviciu in
[*The Degrees of a System of Parameters of the Ring of Invariants of a
Binary Form*](https://doi.org/10.1007/s00031-015-9353-8).

There is an essential degree-selection test before attempting the
zero-fiber geometry. Since \(R_d\) is Cohen--Macaulay in characteristic
zero, if homogeneous invariants of degrees \(e_1,\ldots,e_{\dim R_d}\)
form a system of parameters, then

\[
 H_{R_d}(t)\prod_i(1-t^{e_i})
 \tag{3.7}
\]

is the Hilbert series of an Artinian quotient and therefore is a
polynomial with nonnegative coefficients. Thus a single negative
coefficient rules out that degree sequence independently of the chosen
invariants.

For \(d=3\), the exact weight expansion gives

\[
 [t^{63}]H_{R_3}(t)\prod_{m=1}^{13}(1-t^m)=-2186.
 \tag{3.8}
\]

Consequently \(\mu_1,\ldots,\mu_{13}\), although algebraically
independent, cannot define the nullcone. Their zero fiber necessarily has
a semistable component. The least-total-degree replacement surviving the
same Hilbert test is

\[
 \mu_1,\ldots,\mu_{12},\mu_{14}.                          \tag{3.9}
\]

These corrected moments also have exact Jacobian rank thirteen. Their
proposed Hilbert numerator is nonnegative through degree \(100\), with
last observed nonzero term in degree \(76\) and zeros through degree
\(100\); this is a necessary-test result, not a proof that they are a
system of parameters.

### The apolar-adjoint field obstruction

There is a universal restriction on the algebra generated by the moments.
In the basis \(X^iY^{d-i}\) of \(\operatorname{Sym}^d\), let
\[
 S_{i,d-i}=\frac{(-1)^i}{\binom di}
 \tag{3.10}
\]
be the invariant apolar form, let \(A=C^TD\) be the operator attached to
the coefficient matrix, and define
\[
 \tau(A)=S^{-1}A^TS.
 \tag{3.11}
\]
This is an \(\mathrm{SL}_2\)-equivariant involution.  Direct conversion
back to coefficient coordinates gives
\[
 \boxed{\tau(C)_{ij}=(-1)^{i+j}C_{d-j,\,d-i}.}
 \tag{3.12}
\]
Writing \(C(x,y)=\sum c_{ij}x^iy^j\), the contraction formula is
\[
 \mu_m(C)=\sum_{I=0}^{dm}I!(dm-I)!
 [x^Iy^I]C(x,y)^m.
 \tag{3.13}
\]
Equation (3.12) sends its \(I\)-summand to the
\((dm-I)\)-summand with sign \((-1)^{2I}=1\).  Hence
\[
 \boxed{\mu_m\circ\tau=\mu_m\quad(m\geq1).}
 \tag{3.14}
\]

On the multiplicity-free decomposition (2.2), \(\tau\) acts on
\(\operatorname{Sym}^{2r}\) by \((-1)^r\).  This action is nontrivial on
the invariant ring for every \(d\geq2\).  The first odd invariant occurs
as follows:
\[
\begin{array}{c|c|c}
d&\text{first odd degree}&\text{odd dimension there}\\ \hline
2&6&1\\
3&4&3\\
4&3&1.
\end{array}
\tag{3.15}
\]
For \(d\geq4\), the degree-three contraction of the
\(\operatorname{Sym}^4,\operatorname{Sym}^6,\operatorname{Sym}^8\)
components already gives a nonzero odd invariant.  The \(d=2,3\) rows
follow from the same weight-zero-minus-weight-two character calculation;
their odd multidegrees are respectively
\[
 (3,3),\qquad
 (0,3,1),\ (1,1,2),\ (2,1,1).
 \tag{3.16}
\]

If
\[
 {\cal A}_d=\mathbb Q[\mu_1,\mu_2,\ldots]\subset R_d,
 \tag{3.17}
\]
then (3.14)--(3.16) prove
\[
 \operatorname{Frac}{\cal A}_d
 \subseteq\operatorname{Frac}(R_d)^\tau
 \subsetneq\operatorname{Frac}R_d
 \qquad(d\geq2).
 \tag{3.18}
\]
Thus the moments never generate the full invariant field once \(d\geq2\):
they forget apolar orientation even when they have maximal transcendence
degree.  The algebra conductor vanishes:
\[
 \boxed{({\cal A}_d:R_d)=0\qquad(d\geq2).}
 \tag{3.19}
\]
Indeed, a nonzero conductor element \(a\) would imply
\(r=(ar)/a\in\operatorname{Frac}{\cal A}_d\) for every \(r\in R_d\),
contradicting (3.18).

This does not by itself obstruct integrality: a finite integral extension
may have a nontrivial fraction-field extension and zero conductor.  At
\(d=4\), however, the semistable moment-zero point \(F\) and its quadratic
separator prove nonintegrality directly.  For \(d=2,3\), (3.18) instead
shows that any successful moment parameter system has generic degree at
least two on the invariant quotient.

At \(d=2\) this degree can be sharpened for the known first-six-moment
parameter system.  Multiplying the invariant Hilbert series by
\(\prod_{m=1}^6(1-t^m)\) gives
\[
 1+t^2+t^3+t^4+2t^6+t^8+t^9+t^{10}+t^{12},
 \tag{3.20}
\]
whose coefficient sum is \(10\).  Since the first six moments are already
proved to cut out the nullcone, they form a homogeneous system of
parameters, and
\[
 [\operatorname{Frac}R_2:
   \mathbb Q(\mu_1,\ldots,\mu_6)]=10.
 \tag{3.21}
\]
Adjoining all later moments gives an intermediate field fixed by
\(\tau\).  There are fourteen weighted-degree-seven monomials in
\(\mu_1,\ldots,\mu_6\).  Exact evaluation at fifteen integral points gives
a nonzero \(15\) by \(15\) determinant modulo \(1000003\) after
\(\mu_7\) is appended.  Thus
\[
 \mu_7\notin\mathbb Q[\mu_1,\ldots,\mu_6].
\tag{3.22}
\]
Since \(R_2\) is integral over the polynomial parameter ring and that ring
is normal, its intersection with
\(\mathbb Q(\mu_1,\ldots,\mu_6)\) is the parameter ring itself.  Hence
\(\mu_7\) is not in the first-six moment field.  The
\(\tau\)-fixed intermediate field has prime degree five over that field,
so \(\mu_7\) generates it:
\[
 \boxed{
 \operatorname{Frac}{\cal A}_2
 =\operatorname{Frac}(R_2)^\tau,\qquad
 [\operatorname{Frac}R_2:\operatorname{Frac}{\cal A}_2]=2.}
 \tag{3.23}
\]
Thus in bidegree \((2,2)\) the full moment sequence recovers every generic
invariant except the choice of apolar orientation.

The first relation can also be made exact.  Put
\[
 \nu_m=\frac{\mu_m}{(2m+1)!},\qquad
 x_m=\sum_{p=0}^m\binom mp(-\nu_1)^{m-p}\nu_p.
 \tag{3.24}
\]
These are the centered factorial-normalized moments, with \(x_1=0\).
The seventh centered moment satisfies
\[
 \boxed{
 x_7^5+A_1x_7^4+A_2x_7^3+A_3x_7^2+A_4x_7+A_5=0,}
 \tag{3.25}
\]
where \(A_k\in\mathbb Q[x_2,\ldots,x_6]\) has weighted degree \(7k\)
for weights \(2,\ldots,6\).  The first coefficient is
\[
 A_1
 =-\frac{37975}{1144}x_3x_4
  -\frac{1785}{52}x_2x_5
  +\frac{154875}{1144}x_2^2x_3.
 \tag{3.26}
\]
The five coefficient ansatzes have respectively
\[
 3,\ 19,\ 49,\ 120,\ 227
 \tag{3.27}
\]
possible monomials; the exact relation has
\[
 3,\ 14,\ 37,\ 68,\ 119
 \tag{3.28}
\]
nonzero terms.  The sparse rational coefficients are stored in the
generated artifact.

The certificate uses five good primes to reconstruct the \(418\)
rational coefficients.  It then evaluates (3.25) exactly at \(418\)
fixed integral coefficient matrices.  At the same points, the ansatz
evaluation matrix has rank \(418\) modulo \(1000003\).  Since integrality
and homogeneity already prove that a relation of the form (3.25) exists,
these exact evaluations and the full-rank certificate identify its unique
coefficient vector without expanding a degree-\(35\) polynomial in nine
matrix entries.

Because \(\mu_1,\ldots,\mu_6\) are algebraically independent and \(x_7\)
has field degree five, (3.25) is irreducible over
\(\mathbb Q(x_2,\ldots,x_6)\).  After undoing the triangular centering
change, it generates the prime kernel of
\[
 \mathbb Q[t_1,\ldots,t_7]\longrightarrow R_2,\qquad
 t_m\longmapsto\mu_m.
 \tag{3.29}
\]

Field generation and polynomial generation diverge immediately after this
relation.  For each \(8\leq r\leq12\), the checker evaluates all
weighted-degree-\(r\) monomials in
\(x_2,\ldots,x_{r-1}\) at fixed integral coefficient matrices.  Modulo
\(1000003\), the ranks before and after adjoining \(x_r\) are
\[
 (6,7),\ (7,8),\ (11,12),\ (13,14),\ (20,21),
 \tag{3.30}
\]
respectively.  Each displayed rank jump is an exact characteristic-zero
nonmembership certificate.

Indeed, the triangular centering change gives
\[
 \mathbb Q[\nu_1,\ldots,\nu_{r-1}]
 =\mathbb Q[\nu_1,x_2,\ldots,x_{r-1}].
 \tag{3.31}
\]
Translation changes only \(\nu_1\) in the right-hand coordinates.  Hence,
if the translation-invariant \(x_r\) were a polynomial in the preceding
moments, it would lie in
\(\mathbb Q[x_2,\ldots,x_{r-1}]\); weighted homogeneity would put it in
the span tested in (3.30).  Consequently
\[
 \boxed{x_8,x_9,x_{10},x_{11},x_{12}
 \text{ are all necessary additional generators of }{\cal A}_2.}
 \tag{3.32}
\]
Thus \(\mu_1,\ldots,\mu_7\) generate the moment **field**, but not the
moment **algebra**.

The first post-cutoff cases can now be settled positively without a
nine-variable expansion.  In operator coordinates
\(A=C^T\operatorname{diag}(2,1,2)\), take the affine slice
\[
 A=\begin{pmatrix}
 0&1&0\\ a&b&c\\ d&e&f
 \end{pmatrix}.
 \tag{3.33}
\]
At \(A_0=\left(\begin{smallmatrix}0&1&0\\0&0&1\\1&0&0\end{smallmatrix}\right)\),
the first-row coordinates of the three infinitesimal orbit vectors form
\[
 \begin{pmatrix}
 -2&0&0\\0&0&-2\\0&-1&0
 \end{pmatrix},
 \qquad \det=4.
 \tag{3.34}
\]
Thus the \(\mathrm{SL}_2\)-saturation of (3.33) is dense, so an invariant
identity verified symbolically on this slice is universal.

For \(13\leq r\leq18\), modular interpolation selects independently
evaluated weighted-degree-\(r\) monomial sets in
\(\mathbb Q[x_2,\ldots,x_{12}]\) of sizes
\[
 23,\ 33,\ 38,\ 50,\ 57,\ 75.
\]
After rational reconstruction, the checker expands each proposed
identity exactly in \(\mathbb Q[a,b,c,d,e,f]\) and obtains zero.  Hence
\[
 \boxed{x_{13},x_{14},x_{15},x_{16},x_{17},x_{18}
 \in\mathbb Q[x_2,\ldots,x_{12}].}
 \tag{3.35}
\]
Whether the first twelve moments generate \({\cal A}_2\) in every
degree remains open; (3.35) proves the cutoff through degree \(18\).

There is also a clean normalization statement.  Set
\(S_2=R_2^\tau\).  The first six moments form a parameter system, so
\(R_2\), and hence every intermediate algebra, is finite over
\(\mathbb Q[\mu_1,\ldots,\mu_6]\).  Equation (3.23) gives
\(\operatorname{Frac}{\cal A}_2=\operatorname{Frac}S_2\), while \(S_2\)
is normal.  Conversely, an element of this common field integral over
\({\cal A}_2\) is integral over \(R_2\), hence belongs to the normal ring
\(R_2\), and therefore to \(S_2\).  Thus
\[
 \boxed{\overline{{\cal A}_2}=R_2^\tau.}
 \tag{3.36}
\]
This normalization is already nontrivial in degree two.  All three
quadratic Casimir contractions are \(\tau\)-even, whereas
\[
 \dim({\cal A}_2)_2=2
 \quad\text{and}\quad
 \dim(R_2^\tau)_2=\dim(R_2)_2=3;
 \tag{3.37}
\]
the moment side is spanned by \(\mu_1^2,\mu_2\).  Hence one even
quadratic Casimir direction is missing from the moment algebra.  Since
(3.36) is a finite birational extension, its conductor is nonzero:
\[
 0\ne({\cal A}_2:R_2^\tau),
 \qquad
 ({\cal A}_2:R_2)=0.
 \tag{3.38}
\]
The second equality is the universal orientation obstruction (3.19);
an explicit generating set for the first conductor ideal remains open.

This changes the all-degree architecture. One must first use the Molien
or weight series of \(R_d\) to select admissible moment degrees, and only
then prove that the selected zero fiber is the nullcone by the
first-component and synchronization lemmas. Consecutive initial moments
need not be the correct parameters.

## 4. Stratification by the first nonzero summand

The direct-sum nullcone condition separates into two problems: make the
first nonzero component unstable, then synchronize every higher component
with its unique destabilizing root.

For \(1\leq s\leq d\), consider the stratum

\[
 F_2=\cdots=F_{2s-2}=0,\qquad F_{2s}\ne0.
 \tag{4.1}
\]

A proof of \(\mathrm{MN}_d\) would follow from the following two lemmas.

> **First-component lemma.** On (4.1), vanishing of all moments forces
> \(F_{2s}\) to have a root \(L\) of multiplicity at least \(s+1\).

> **Synchronization lemma.** Once this root \(L\) exists, moment vanishing
> forces \(L^{r+1}\mid F_{2r}\) for every \(r>s\).

Because a root of multiplicity \(>s\) in a degree-\(2s\) form is unique,
the second lemma has a canonical flag to use. If \(F_{2s}=0\), one simply
moves to the next stratum. This avoids choosing a global normal form
before the moments have produced one.

### The quadratic anchor

The first stratum is the most useful starting point. Write

\[
 F_2=aX^2+2bXT+cT^2,\qquad
 \Delta_2=b^2-ac.
 \tag{4.2}
\]

The first-component lemma for \(s=1\) becomes the concrete target

\[
 \Delta_2\in\sqrt{M_d}.
 \tag{4.3}
\]

If \(F_2\ne0\), equation (4.3) gives \(F_2=L^2\), after which the
synchronization lemma asks successively for

\[
 L^3\mid F_4,\quad L^4\mid F_6,\quad\ldots,\quad
 L^{d+1}\mid F_{2d}.
 \tag{4.4}
\]

If \(F_2=0\), the same argument restarts with the quartic component.
Thus (4.3) is an anchor, not an assumption that every moment-zero point
has a nonzero quadratic component.

## 5. Current evidence

The evidence must be kept at its proved strength.

| locus | result | status |
|---|---|---|
| all of \(V_2\) | first six moments have the full one-sided nullcone radical | exact over \(\mathbb Q\) |
| all of \(V_3\) | \(\mu_1,\ldots,\mu_{13}\) are algebraically independent, but the degree-\(63\) Hilbert numerator coefficient is \(-2186\), so they cannot define the nullcone | exact over \(\mathbb Q\); an extra semistable zero component exists |
| all of \(V_3\) | \(\mu_1,\ldots,\mu_{12},\mu_{14}\) have exact Jacobian rank thirteen and pass the necessary Hilbert numerator test through degree \(100\) | exact over \(\mathbb Q\); corrected zero fiber still open |
| rank-one Segre cone \(\Sigma_3\subset V_3\) | moments \(1,2,3,4\) cut out exactly \(\Sigma_3\cap N_3\); every fixed coordinate multiplier has cutoff \(m>\deg Q\) | exact in characteristic zero |
| rank-one Segre points in any \(V_d\) whose operator symbol has at most two roots | all pure moments force the one-sided nullcone, with mixed cutoff \(m>\deg Q\) | exact in characteristic zero via the one-variable Laurent constant-term theorem |
| rank-one quartic \((2,1,1)\)-symbol orbit in \(V_4\) | moments \(1,\ldots,5\) cut out its three one-sided components; five nullcone generators have fourth-power certificates | exact over \(\mathbb Q\) |
| generic rank-one squarefree-quartic symbol orbit in \(V_4\) | moments \(1,\ldots,6\) cut out its four annihilator lines on a nonempty Zariski-open cross-ratio set | exact in characteristic zero; finitely many exceptional squarefree orbits remain possible |
| complete rank-one Segre cone \(\Sigma_d\subset V_d\), every \(d\) | all pure moments imply eventual vanishing for every fixed coordinate-only multiplier by complete factor polarization and Laurent Newton separation | exact in characteristic zero; proves balanced GVC but not Segre moment--nullcone equality |
| maximal-torus fixed diagonal slice in \(V_3\) | moments \(1,2,3,4\) have only the origin as a common zero, with seventh-power certificates for all four diagonal coefficients | exact over \(\mathbb Q\) |
| full non-null \(F_2\) branch in \(V_3\) | after \(F_2=2XT\), five residual-torus chart orbits cover the non-diagonal locus; \(\mu_2,\ldots,\mu_{12}\) have exact Jacobian rank eleven on every representative chart | exact over \(\mathbb Q\); all five affine zero fibers remain open |
| corrected non-null \(F_2\) chart systems in \(V_3\) | \(\mu_2\) eliminates the opposite-weight variable with a constant nonzero pivot on all five charts; on \(s_0=1\), reduced \(\mu_3\) gives two explicit principal-open pivots and the boundary where both vanish | exact over \(\mathbb Q\); the resulting smaller systems remain open |
| sparse plane in the \(s_0=1,\ A=B=0\) boundary | after the \(\mu_2\) pivot, \(\mu_3=1866240a^3\) and \(\mu_4=138240(11249-8776ab-901a^2b^2)\), with an explicit unit certificate | exact over \(\mathbb Q\); this plane is excluded |
| four-parameter family in the \(s_0=1,\ A=B=0\) boundary | moments \(3,\ldots,6\) leave a quotient of length \(372\), and moment \(7\) gives the unit ideal | exact over \(\mathbb Q\); strictly contains the sparse plane |
| three \(s_0\)-chart pivot strata | \((\mu_2,\ldots,\mu_{12})\) has exact generic differential ranks \(11,10,9\) on \(A\ne0\), \(A=0,B\ne0\), and \(A=B=0\) | exact over \(\mathbb Q\); supports finite-quotient attacks |
| pure \(\operatorname{Sym}^2,\operatorname{Sym}^4,\operatorname{Sym}^6\) in \(V_3\) | moments cut out the corresponding binary-form nullcones | exact over \(\mathbb Q\) |
| \(\operatorname{Sym}^4\oplus\operatorname{Sym}^2\) in \(V_3\), with \(F_2=2cXT\) | moments through order six imply \(c^6=0\) | exact over \(\mathbb Q\) |
| \(\operatorname{Sym}^6\oplus\operatorname{Sym}^2\) in \(V_3\), with \(F_2=2cXT\) | even moments through order fourteen imply \(c^{25}=0\) | only over \(\mathbb F_{32003}\) |
| full mixed \(V_3\) | no moment--nullcone equality yet | open |
| explicit full-rank point in \(V_4\) | all pure moments vanish but the fixed \(Z\)-mixed moment is \((4m+2)!m!/(2m+1)!!\) | exact all-order counterexample; \(\mathrm{MN}_4\) is false |

The normalized \(c\)-certificates say that the tested branches cannot have
\(\Delta_2\ne0\). They motivate (4.3), but the finite-field
\(c^{25}\) membership is not a characteristic-zero certificate and neither
slice proves the global quadratic-anchor statement.

## 6. What does and does not generalize

### Arbitrary direct sums of binary forms

The geometric half is already general. For

\[
 V=\operatorname{Sym}^{n_1}\oplus\cdots\oplus
   \operatorname{Sym}^{n_s},
 \tag{6.1}
\]

Hilbert--Mumford says that \((H_1,\ldots,H_s)\) is in the
\(\mathrm{SL}_2\)-nullcone exactly when all nonzero \(H_i\) have a common
root of multiplicity \(>n_i/2\). Thus the first-component and
synchronization architecture applies to any direct sum of binary forms.

What is special here is the invariant sequence \(\mu_m=\mathcal E_2(f^m)\).
An arbitrary direct sum has no canonical contraction moments with the
required properties. The generalization is therefore a reusable nullcone
geometry, not a universal moment--nullcone theorem.

### \(n\) contraction pairs

There is a natural higher-rank formulation. Let \(U=\mathbb C^n\) and

\[
 V_{n,d}
 =\operatorname{Sym}^d(U^*)\otimes\operatorname{Sym}^d(U)
 =\operatorname{End}(\operatorname{Sym}^dU).
 \tag{6.2}
\]

Pieri's rule gives the multiplicity-free \(\mathrm{SL}_n\)-decomposition

\[
 \boxed{
 V_{n,d}\cong
 \bigoplus_{j=0}^{d}
 V_{j(\omega_1+\omega_{n-1})}},
 \tag{6.3}
\]

where \(V_\lambda\) denotes the irreducible highest-weight module of
highest weight \(\lambda\). One derivation writes
\((\operatorname{Sym}^dU)^*\), up to a determinant twist, as the Schur
module of rectangular shape \((d^{\,n-1})\). Adding a horizontal
\(d\)-strip gives precisely

\[
 (d+j,d,\ldots,d,d-j),\qquad 0\leq j\leq d,
 \tag{6.4}
\]

which becomes \(j(\omega_1+\omega_{n-1})\) after removing the determinant
twist. When \(n=2\), formula (6.3) reduces to (2.1).
As a dimension check, Weyl's formula gives

\[
 \dim V_{j(\omega_1+\omega_{n-1})}
 =\frac{2j+n-1}{n-1}
   \binom{j+n-2}{n-2}^{\!2},
 \qquad
 \sum_{j=0}^{d}\dim V_{j(\omega_1+\omega_{n-1})}
 =\binom{n+d-1}{d}^{\!2}
 =\dim V_{n,d}.
\]

Define

\[
 \mu_{n,m}(f)=\mathcal E_n(f^m).
 \tag{6.5}
\]

Because \(f\) is balanced, this is a scalar
\(\mathrm{SL}_n\)-invariant. The formal higher-rank question is

\[
 \mathrm{MN}_{n,d}:\qquad
 V(\mu_{n,1},\mu_{n,2},\ldots)
 \stackrel{?}{=}\mathcal N(V_{n,d}).
 \tag{6.6}
\]

The nullcone is now described by an arbitrary destabilizing one-parameter
subgroup, or equivalently by a weighted flag in \(U\); there is generally
no reduction to one common point of \(\mathbb P^1\).

The easy direction still holds in every rank:

\[
 \mathcal N(V_{n,d})
 \subseteq V(\mu_{n,1},\mu_{n,2},\ldots),
 \tag{6.7}
\]

because every positive-degree invariant vanishes on the nullcone.
Moreover, a balanced nullcone point satisfies eventual mixed contraction
vanishing. Indeed, choose a one-parameter subgroup for which every weight
of \(f\) is positive. Weights in \(f^m\) then grow at least linearly in
\(m\). A fixed \(g\) contributes bounded weight, while every monomial that
survives \(\mathcal E_n\) has nonnegative residual coordinate exponents
whose total degree is fixed by the bidegree of \(g\), hence belongs to a
finite set of weights. Equivariance gives a contradiction for large
\(m\). Therefore

\[
 f\in\mathcal N(V_{n,d})
 \Longrightarrow
 \mathcal E_n(gf^m)=0\quad(m\gg0).
 \tag{6.8}
\]

### Exact obstruction in three or more pairs

The converse (6.6) is false as soon as \(n\geq3\), already for \(d=2\).
The repository's
[four-term three-pair counterexample](THREE_PAIR_IMAGE_MATHIEU_COUNTEREXAMPLE.md)

\[
 f=\tau(t-y)(wz+vt),\qquad g=y
 \tag{6.9}
\]

lies in \(V_{3,2}\) and satisfies

\[
 \mathcal E_3(f^m)=0,\qquad
 [t]\mathcal E_3(gf^m)=(-1)^{m-1}(m+1)!\,m!
 \ne0
 \quad(m\geq1).
 \tag{6.10}
\]

If \(f\) were in the \(\mathrm{SL}_3\)-nullcone, (6.8) would force the
second expression eventually to vanish. Hence \(f\) is a moment-zero
semistable point and

\[
 \boxed{\mathrm{MN}_{3,2}\text{ is false}.}
 \tag{6.11}
\]

Ignoring additional contraction pairs preserves both identities in
(6.10), so the same argument proves

\[
 \boxed{\mathrm{MN}_{n,2}\text{ is false for every }n\geq3.}
 \tag{6.12}
\]

This was the first decisive boundary of the generalization.  The later
two-pair witness shows that rank two also fails, at \(d=4\).  Rank two
still has the special common-root description of instability, but the
contraction moments do not force that instability in every degree.

### Nonhomogeneous forms and positive characteristic

A single unequal bidegree is one-sided, but the proposed balanced theorem
already fails in bidegree \((4,4)\).  A sum of several bidegrees has the
additional complication that mixed products can cancel the central
dual-minus-coordinate grading.

All statements in this program are over characteristic zero, with
nullcone tests made after scalar extension to an algebraic closure. In
small positive characteristic the displayed Clebsch--Gordan decomposition
need not be semisimple, so no characteristic-free version is asserted.

## 7. The next attack

This section concerns the residual degree-three moment--nullcone
classification, not the opening of the SIC counterexample search.  For
SIC in bidegree \((4,4)\), the next attack is exact rank two as described
in the rank-stratified programme.  The calculations below should discover
small invariant or covariant identities, rather than eliminate the full
coefficient space.

1. **Construct uniform Clebsch--Gordan coordinates.** Express the
   projections \(f\mapsto F_{2r}\) and the contractions \(\mu_m\) by
   transvectants, with \(d\) retained as a parameter.
2. **Select admissible parameter degrees.** Compute the Molien/weight
   Hilbert series before attacking a zero fiber. In \(d=3\), degrees
   \(1,\ldots,13\) are now excluded and the first corrected target is
   \((1,\ldots,12,14)\). Seek a uniform rule for choosing
   \((d+1)^2-3\) moment orders whose proposed Hilbert numerator is
   nonnegative.
3. **Test the corrected minimal zero fiber.** Determine whether
   \(V(\mu_1,\ldots,\mu_{12},\mu_{14})\) equals the nullcone, preferably
   by invariant/covariant saturation. If it fails, isolate the extra
   semistable component and move only to another Hilbert-compatible
   degree set; do not merely extend a consecutive cutoff.
4. **Prove the quadratic anchor.** Search for a finite identity
   \(\Delta_2^N\in(\mu_{m_1},\ldots,\mu_{m_k})\) in invariant coordinates,
   using the explicit global projection
   \[
   r_0=(3c_{10}+2c_{21}+3c_{32})/10,\quad
   r_1=(-9c_{00}-c_{11}+c_{22}+9c_{33})/20,\quad
   r_2=-(3c_{01}+2c_{12}+3c_{23})/10.
   \]
   The torus-fixed slice is now closed exactly; the next test must allow
   generic higher-component weights.
5. **Prove one synchronization step.** On the incidence chart
   \(F_2=L^2\), show that the moments force the forbidden coefficients of
   \(F_4\) to vanish, equivalently \(L^3\mid F_4\). Formulate the
   certificate covariantly so that \(L\) can be eliminated without a large
   global Gröbner basis.
6. **Generalize the step.** With
   \(L^{j+1}\mid F_{2j}\) for \(j<r\), use the lowest remaining
   \(\mathrm{SL}_2\)-weights in the moments to force
   \(L^{r+1}\mid F_{2r}\).
7. **Handle a zero anchor.** Repeat the same construction on the strata
   \(F_2=0\), then \(F_2=F_4=0\), so the proof does not discard components
   whose lower summands vanish.

The first decisive milestone is therefore not a rational reconstruction
of the isolated \(c^{25}\) calculation. The consecutive dimension-sized
set is now excluded and the corrected set (3.9) is the first viable
target. The next milestone is a global \(d=3\) quadratic-anchor
certificate for that corrected ideal, stated in
\(\mathrm{SL}_2\)-invariant form, followed by the first common-root
synchronization certificate. Those two identities would expose the
pattern needed for arbitrary \(d\).

### Decision tree after the present tests

The best next calculation is the localized non-null-quadratic problem.
On \(\Delta_2\ne0\), use \(\mathrm{SL}_2\) to put

\[
 F_2=2cXT,\qquad c\ne0,
\]

then use overall homogeneity to set \(c=1\). The scalar component has
already been removed by \(\mu_1\), leaving twelve variables
\((F_4,F_6)\) and the twelve corrected equations
\(\mu_2,\ldots,\mu_{12},\mu_{14}\). The residual diagonal torus gives
weights

\[
 (6,4,2,0,-2,-4,-6;\;4,2,0,-2,-4).
\]

There are already \(246354\) residual-weight-zero monomials of degree at
most thirteen, and the corrected system also contains order fourteen, so
a raw expanded Gröbner calculation remains large. Recorded 180-second
runs with Singular and `msolve` do not close even \(s_0=1\), over
\(\mathbb F_{101}\) or \(\mathbb F_{43}\). This is a bounded failure, not
evidence for a component.

The first triangular layer is now exact. On the five representative
charts, \(\mu_2\) has constant nonzero derivatives

\[
 -72,\quad432,\quad-1080,\quad336,\quad-1344
\]

in the respective opposite variables \(s_6,s_5,s_4,t_4,t_3\).
Consequently it eliminates that variable globally on each chart. On
\(s_0=1\), the reduced third moment is

\[
 \mu_3=-103680\,A\,s_5-17280\,B\,t_4+C,
\]

with \(C\) independent of \(s_5,t_4\) and explicit polynomials \(A,B\)
recorded in the bidegree-\((3,3)\) frontier. The chart therefore splits
into \(A\ne0\), \(A=0,B\ne0\), and \(A=B=0\). The first two branches
permit a second variable elimination; the third gains two boundary
equations. This branchwise triangularization, rather than a raw basis, is
the next attack.

The first exact boundary piece is also closed. On the two-parameter plane

\[
t_0=a,\quad t_3=3a,\quad t_4=b,\quad
s_6=(14ab+70)/3
\]

with \(s_0=1\) and the other higher coordinates zero, \(A=B=\mu_2=0\)
identically, while \(\mu_3=1866240a^3\) and the normalized fourth moment
has nonzero constant \(11249\) modulo \(a^3\). Hence
\((\mu_3,\mu_4)=(1)\) on that plane. This validates the triangular
strategy but does not describe the remaining common boundary.

Exact rational points on all three pivot strata also show that the
Jacobian of \(\mu_2,\ldots,\mu_{12}\) is invertible before restriction.
Since \(dA,dB\) are independent, the restricted moment maps have maximal
generic ranks \(11,10,9\). Consequently none of the three strata carries
a representation-theoretically forced positive-dimensional moment fiber.
The next algebraic target is the finite special fiber at the moment
origin, not another generic-rank calculation.

The outcomes have clear implications:

1. If the dehomogenized ideal is the unit ideal over \(\mathbb Q\), then
   \(\Delta_2\) lies in the radical of the corrected moment ideal: the
   global quadratic anchor is proved for the viable minimal set.
2. If it has a component, test whether the component survives all higher
   moments. A surviving exact point or recurrence is a bidegree-\((3,3)\)
   semistable moment-zero point; a component killed by later moments only
   shows that the corrected set is still not sufficient. Its geometry
   should determine the next Hilbert-compatible replacement.
3. Once the non-null branch is excluded, move to the incidence chart
   \(F_2=L^2\) and prove \(L^3\mid F_4\). This is the first synchronization
   lemma.
4. Do not interpolate these identities to all \(d\): the \(d=4\)
   counterexample forbids that conclusion. A full raw Gröbner basis in
   sixteen coefficients remains the least informative route for \(d=3\).

For computation, the preferred order is now: eliminate with \(\mu_2\);
split and eliminate with \(\mu_3\); decompose the \(A=B=0\) boundary into
covariant or coordinate strata; and represent the two principal opens
without expanding high powers of the inverse pivot. Run modular sparse
elimination only on those smaller systems, reconstruct any unit
certificates over \(\mathbb Q\), then repeat the triangular step on the
other four chart orbits. A finite-field unit ideal is evidence only until
the rational certificate is reconstructed.

## 8. Claim boundary

The conjecture \(\mathrm{MN}_d\) for all \(d\), and hence unrestricted
\(\operatorname{SIC}(2)\), is false.  This program does not decide
\(\mathrm{MN}_3\).  Its surviving contributions are:

- an exact all-\(d\) description of the one-sided nullcone as a common
  high-multiplicity-root condition;
- a stratified proof architecture with two explicit lemma families;
- the exact Hilbert-series obstruction to consecutive degrees
  \(1,\ldots,13\), plus the algebraically independent corrected candidate
  \(1,\ldots,12,14\);
- a precise first invariant target, \(\Delta_2\in\sqrt{M_d}\);
- a rule for using the existing exact and finite-field calculations
  without promoting experiments to theorems.

It also proves the structural statements (2.7), (3.4)--(3.8), the
corrected Jacobian independence in (3.9), and the higher-rank obstruction
(6.11)--(6.12). These are consequences of exact calculation and standard
invariant theory together with the repository's proved counterexamples.

The cited literature supplies the Image/Mathieu framework, the
Hilbert--Mumford and binary-form nullcone criteria, and the invariant-ring
tools. It does not assert the now-falsified moment--nullcone conjecture
\(\mathrm{MN}_d\); that conjecture and the proposed
anchor/synchronization proof were specific to this repository program.

## 9. Sources

- A. van den Essen, D. Wright, and W. Zhao,
  [*On the Image Conjecture*](https://arxiv.org/abs/1008.3962),
  J. Algebra 340 (2011), 211--224. This supplies the Image/Mathieu
  framework and the contraction setting.
- H. Derksen,
  [*Constructive Invariant Theory*](https://www.cse.iitb.ac.in/~sohoni/CS782/DerksenGIT.pdf),
  Sections 2.2 and 4. This is used for the invariant-theoretic nullcone,
  Hilbert--Mumford criterion, and finite invariant bounds.
- A. E. Brouwer and M. Popoviciu,
  [*Sylvester versus Gundelfinger*](https://sigma-journal.com/2012/075/),
  SIGMA 8 (2012), 075. This records the common high-multiplicity-root
  criterion for direct sums of binary forms and Hilbert's nullcone
  criterion for systems of parameters.
- A. E. Brouwer, J. Draisma, and M. Popoviciu,
  [*The Degrees of a System of Parameters of the Ring of Invariants of a
  Binary Form*](https://doi.org/10.1007/s00031-015-9353-8),
  Transform. Groups 20 (2015), 953--967.
- C. D. Long,
  [*Counterexamples to the \(xz\)-Conjecture and the Mathieu Conjecture for
  \(SU(2)\)*](https://arxiv.org/abs/2607.19012), arXiv:2607.19012v1
  (2026). The repository's three-pair witness is a bihomogeneous lift of
  Long's \(SU(2)\) seed.
- M. Müger and L. Tuset,
  [*The Mathieu Conjecture for \(SU(2)\) Reduced to an Abelian
  Conjecture*](https://arxiv.org/abs/2210.06582),
  Indag. Math. 35 (2024), 114--118.
