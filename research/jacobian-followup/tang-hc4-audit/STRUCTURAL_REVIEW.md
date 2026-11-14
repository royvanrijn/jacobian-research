# Full structural proof review of Tang v1

Reviewed 13 November 2026. This note completes the written audit left open by
[HCTANG1](README.md). It reproduces Qiyue Tang's external theorem, not a new
attribution of that theorem. The source is the byte-pinned
[31-page manuscript](source-v1.pdf), Hexagon 2611.00156v1. The independent
local calculations in HCTANG1 remain part of the evidence. The review below
is a mathematical proof audit; no finite checker certifies its geometric
classification or all its arguments. No formalization or external review
is claimed.

**Conclusion (HCTANG2).** For a complex polynomial f in four variables
with nonzero constant Hessian determinant, either its gradient is an
automorphism, or after an invertible constant linear source change

\[
 f=R(x,y)+zA(x,y)+wB(x,y),\qquad J(A,B)\in\mathbb C^\times.
\]

Consequently HC4 and JC2 are equivalent. Neither conjecture is resolved.
The quantitative and repository-specific consequences are proved separately
in [APPLICATIONS.md](APPLICATIONS.md).

## 1. Geometry behind the singular-Hessian theorem

Tang Theorem 2.1 concerns a rational V in three variables with nilpotent
N=DV, symmetrized by a generically nonsingular **polynomial Hessian** H.
Its conclusion is that the image of V lies in an affine plane. Nonsingularity
means that det H is not identically zero, not that it is a unit. Polynomiality
is essential, as the rational counterexamples in HCTANG1 demonstrate.

At rank one, write V=A(t), t=p/q reduced, using the characteristic-zero
rank criterion and Luroth. Restriction to a suitable affine line shows
that the image curve is unirational, which supplies this application of
Luroth. If its image is nonplanar, the components of A' are independent.
Nilpotence gives A'(t).grad t=0. Clearing denominators, primitive
divisibility by p-Tq and comparison of source degrees give the stronger
identity A'(T).grad(p-Tq)=0 over C(T). The coefficient matrix of this last
gradient has rank at most two and entries affine in T.

At rank two its cross-product kernel has degree at most two. Independence
of A' makes its three coefficient vectors independent, so a constant
conjugation puts that kernel at (1,T,T^2). The resulting identities are
p_x=0, p_y=q_x, p_z=q_y, q_z=0. Integrating and comparing coefficients
forces p,q affine and makes the same coefficient matrix rank one, a
contradiction. At rank one, p-Tq depends on one linear form over C(T).
If its source degree is at least two, taking two T derivatives of its
top power forces that linear form to have constant coefficient ratios.
This would give a constant relation on A'. The remaining possibility is
that p,q are affine with independent linear parts.

Use p=x,q=y. Set F(t,y,z)=f(ty,y,z) and k=A'_3/A'_2. Self-adjointness
gives F_y+k F_z depending only on t. The coefficient of y^j z^ell in F
has t-degree at most j. Its highest z coefficient is therefore constant.
If the z degree is at least two, the next coefficient forces k to be
affine in t, giving a constant relation on A'. At z degree at most one,
the same degree constraint makes f affine. Both alternatives contradict
the hypotheses. This is an all-degree argument, not a bounded search.

At rank two, the exact first-jet calculation retained in HCTANG1 proves
that im N is involutive. Assuming a nonplanar image, choose local image
coordinates (u,v,w(u,v)), allowing algebraic w, and put

\[
 M=((1,0,w_u)^T\mid(0,1,w_v)^T),\quad J_0=D(u,v),\quad N=MJ_0.
\]

The reversed product is nonzero, square-zero and rank one. After exchanging
the first two coordinate axes if necessary it is

\[
 J_0M=b\begin{pmatrix}-\beta&1\\-\beta^2&\beta\end{pmatrix}.
\]

Involutivity of the two column fields of M gives
w_uu+2 beta w_uv+beta^2 w_vv=0. At least one second derivative is nonzero,
since otherwise the image would be planar. Thus beta is algebraic over
the image field. On each geometric generic fibre, w_u,w_v,beta are
constant. Its tangent direction ker N=im N^2 is then the constant vector
(1,beta,w_u+beta*w_v). Two independent annihilating linear forms are
constant on that curve; its closure is a line.

The closure B of these projective lines in the Grassmannian is an
irreducible surface. The conditions that V be defined and constant on a
line hold on an open subset of B: in a line chart they say that the
numerator and denominator coefficient vectors are proportional, with a
nonzero denominator. The omitted at-most-one-dimensional subfamily
sweeps at most a surface. At a general source point every remaining
line must have tangent ker N, so exactly one passes through it.

Resolve B and pull back the universal family. This is the precise
congruence convention of De Poi, Section 1: a smooth P1 bundle over a
resolution of a Grassmannian surface, with evaluation to P3. Evaluation
is proper birational. Its ramification divisor has degree two on a
general incidence fibre, from the canonical degrees -2 and -4.
A proper birational map to the normal P3 is an isomorphism over the
generic point of every divisor; hence each ramification component maps
to dimension at most one. A component mapping to a point must dominate
B, because a fibre line embeds and can meet that point only once. All
general lines then pass through that point. Otherwise the focal locus
has dimension one.

The primary [De Poi Theorem 0.1](https://arxiv.org/html/math/0407339v1)
explicitly includes reducible and nonreduced focal curves. It gives
secants of a twisted cubic, pencils in planes through a fixed line with
centres on that line, or joins of a fixed line and a curve. The point
case just proved is added separately. Smoothness of the Grassmannian
surface itself is not an extra hypothesis. These are all the required
families, including their placements relative to the affine chart.

The identity ker N=im N^2 and self-adjointness imply that ker N is
H-orthogonal to im N. In particular each fibre direction r is H-isotropic.
The following table accounts for every placement of the focal endpoints.

| Line family | Nilpotence and the polynomial-Hessian obstruction |
| --- | --- |
| Two finite generic endpoints | Write X=(1-lambda)A(s)+lambda B(t). In the frame A',B',r=B-A, trace zero removes the A' component of W_s and the B' component of W_t. The second characteristic coefficient says that the other two coefficients alpha,gamma satisfy alpha*gamma=0. Independence leaves exactly one nonzero. If alpha=0, im N=span(A',r). For its normal n, Hr is proportional to n, so polynomial regularity at A integrates to grad f(X)=grad f(A)+n R(lambda), while f(X)=f(A)+lambda grad f(A).r. Differentiating in t forces R=0 because n.B' is nonzero. The gradient would depend only on s, contrary to det H nonzero. The other sign uses B. |
| Star at a finite point c | The identity (X-c)^T H(X)(X-c)=0 kills every homogeneous term of degree at least two in the polynomial expansion at c. |
| Star at infinity | Use the common direction e3. Then V depends only on x,y and f=zP(x,y)+Q(x,y). Trace zero and H-orthogonality make the nonconstant planar image parameter algebraic over P. Generic polynomial levels of P have constant tangent directions and their closures are lines. Lines in different levels are disjoint, so all are parallel. This gives a constant relation on the first two components of V. The constant planar pair already gives a relation. |
| Pencils with a finite fixed line | Isotropy and expansion at the finite centre make f affine on each pencil plane. Expanding at any common finite point of the fixed line kills every higher homogeneous part on a dense union of planes. |
| Pencils with the fixed line at infinity | Use X=(lambda,q+a(t)lambda,t). If a is constant the lines are parallel. Otherwise trace zero forces W_3,q=0; the remaining derivative W_3,t is an eigenvalue and hence zero. |
| Joins, fixed line infinite and curve finite | The directions at a general finite endpoint sweep a fixed two-dimensional vector plane. The same polynomial expansion makes f affine on each corresponding affine plane. Thus f=alpha(z)x+beta(z)y+gamma(z), whose Hessian is singular. |
| Joins, fixed line finite and curve infinite | Use X=(lambda,t lambda,s+lambda q(t)), with q possibly algebraic, and W=(a,t a+b,c). Trace gives b_t=-a and c_s=q a_s+q' b_s. Rank two forces b_s nonzero. The reversed-product determinant gives c_t=q a_t. Compatibility then gives q''=0. In an algebraic function field over C(t), zero derivative means a constant, so q is affine. The constant covector (tq'-q,-q',1) annihilates W_s,W_t. |
| Both endpoint loci infinite | Their joins stay in the plane at infinity and cannot dominate affine three-space. |

In the two-finite-endpoint row, passing to an algebraic parameter cover
does not alter the argument: it is a function-field identity, and
polynomial f is regular at both finite endpoints regardless of whether
V is defined there. The uses of polynomiality in this table cannot be
replaced by generic rational regularity. This completes Theorem 2.1.

For a four-variable polynomial h with Hessian rank three, take a transverse
slice with invertible three-by-three Hessian A, and v=A^-1 grad h_w on
that slice. HCTANG1 proves the characteristic transport identity and,
for polynomial h, nilpotence of Dv and symmetry of A Dv. Theorem 2.1
therefore puts v in an affine plane b.v=d. If R defines the gradient-image
hypersurface and C=grad R(grad h), its kernel direction on the slice is
(-v,1). Hence (b,d).C=0 on a dense image. The derivative of R in that
constant direction vanishes on R=0 and has smaller degree than R; it
vanishes identically. This proves Tang Theorem 2.3 and its singular-rank
relation corollary. No rational auxiliary potential was silently treated
as polynomial in this assembly.

## 2. Common rules for changing weights

Let E be the leading weight and W=4E-2 sum(w_i). The determinant-weight
identity gives W>=0. If W>0 the leading Hessian is singular; if W=0 its
determinant is the same nonzero constant as the full determinant.
Finite support supplies a strict gap until the first new wall.

At W=0 with positive weights, write A=Hess P(0). Both A and A^-1 pair
only complementary weight blocks. Consequently A^-1 grad f has ith
component of weight at most w_i and linear part x_i within its own
weight block. A nonlinear monomial of that weight uses only lower-weight
variables. Successive substitution supplies an actual polynomial inverse.
This also justifies stopping at W=0 when another term ties there.

The primary [de Bondt Corollary 7.3 and Theorem 7.2](https://arxiv.org/html/1609.03904v2)
give the singular forms used here. At rank three, one apex gives a
rank-one complementary Hessian over the active-variable field; two give
joint affineness in the complementary pair; three give a constant kernel
line. At rank two, the alternatives are a constant two-dimensional
kernel or one apex with a complementary affine triple. At rank one
there is a three-dimensional constant kernel. Homogeneous tops of rank
at most three depend on exactly that many variables after a constant
change. The extra algebraic-dependence conditions in de Bondt's formulas
are respected; the common binary coefficient pair is not assumed free.

Apex spaces and constant kernels split across equal-weight blocks by
equivariance. Old homogeneous minors survive as leading parts at a wall.
If a passive vector has a unique coordinate of greatest old weight and
the corresponding principal cofactor is nonzero, its adjugate contraction
cannot vanish. An apex must have zero such contraction. This remains
usable when the new weights coincide.

For an apex plane that is a graph from base to fibres, the transpose
source shear aligns it with the target base plane. It preserves new
weights. In every continuing use, the replacement variable has strictly
smaller old weight, so it also preserves the old leading form. This
direction of the shear is essential.

If nonzero positive-weight binary coefficients A,B have zero Jacobian,
weighted Euler gives p A dB=q B dA. Unique factorization yields
A=a S^r,B=b S^s with S weighted homogeneous. If S has zero constant
directional derivative it is a polynomial in a constant linear form.
With unequal positive base weights, such a homogeneous S is a power of
one base coordinate. Otherwise it is genuinely binary.

## 3. Lemmas 3.2–3.8: coverage and termination

**3.2, genuinely binary common coefficient.** The curvature polynomial
U=S_y^2 S_xx-2 S_x S_y S_xy+S_x^2 S_yy is nonzero. Indeed U=0 would
make generic irreducible polynomial level components lines, and comparison
of lines in different levels makes them parallel, contrary to genuine
binary dependence. For H=z A+w B the two fibre cofactors and every
mixed base/fibre minor are nonzero, by the displayed cofactor formulas.
Raise both fibres with slope one. E has slope one and W is constant.
At a wall, the cofactors exclude an apex in the fibre plane. A single
apex would have an annihilating old-homogeneous base direction and the
other fibre axis; their nonzero surviving minor contradicts complementary
rank at most one. Three apices also intersect that fibre plane. Thus
the apex plane is a graph over the base. The graph shear restores H
exactly: a new tied fibre-degree-j term had old weight E+(1-j)tau,
so no fibre-free term is present and the fibre-linear part is unchanged.
Each genuine continuing wall is a base/fibre weight equality. There
are at most four. Beyond them no fibre-degree-at-least-two term can
remain, since it would eventually overtake H. The full f is cotangent.

**3.3, three distinct phase entries.** For F(x,L), L=ell(x).Y, the
cofactor deleting a passive variable is
-F_LL F_L^2(ell'_i ell_j-ell'_j ell_i)^2. All three are nonzero for
phase entries proportional to 1,x^r,x^s with distinct integral exponents.
Raising all passive weights gives W slope 4m-6>0. At every wall the
apex is the active axis. The old passive-degree-m coefficients recover
exactly the same phase. A genuine continuation strictly increases m.
Ordinary polynomial degree bounds these increases.

At the final maximal passive degree D, its whole homogeneous passive
part N_D still has the three cofactors. Weights (0,1,1,1) and the old
contraction test force its only apex to be active. Thus
N_D=a(x)(ell(x).Y)^D. Its three rational phase coefficients are independent
over C, or a constant kernel would make a rank-three image a hyperplane
with three apices. Their Wronskian is nonzero. At a regular x0 choose
passive coordinates from ell,ell',ell''/2 and write x=x0+xi. Every
monomial xi^k u^i v^j t^h of N_D has k>=j+2h, with equality part
a(x0)(u+xi v+xi^2 t)^D.

The deliberately negative active weight path is
(-eta,1,1+eta,1+2eta). A lower passive degree at least two can reach the
top only by eta<=D-2; the affine part has been removed, so passive
degree one cannot do so before D-1. Such an earlier wall has W>=2,
surviving cofactors, and the same phase. Any lower phase power would
require the negative exponent (j-D)/eta in its polynomial pure-u
coefficient. Thus no such wall exists. At xi=0 all terms of passive
degree at least two depend only on u, making the full Hessian rank at
most three. This contradicts its unit determinant. Positive weights
are not used on this negative-weight path; the inverse lemma is not
invoked there.

**3.4, two phase entries plus an affine partner.** For
H=a x^ell(w+x^A y)^m+b x^p z, increase the phase weights at slope one
and the absent variable at slope m. W has slope 2m-4>=0. The exact
kernel and adjugate formulas retained by HCTANG1 give three nonzero
passive cofactors. They force the active apex at every wall, even at
new equalities. Old weights recover the same two phase coefficients
and exactly the old affine partner. A third phase entry with a new
exponent invokes 3.3. Equal exponents are removed by a constant shear
at a crossing from below. Other continuations strictly increase m.
The ordinary degree bounds those increases; the absent weight minus
each phase weight strictly increases, allowing at most two shears.
After the last wall the monomial inequality is
deg_y+deg_w+m deg_z<=m. Hence f=F(x,y,w)+z B(x), and
det Hess f=-(B')^2 det Hess_(y,w) F. Units force B' constant and p=1.
HC2 over C(x) gives a fibre inverse. Its unique formal inverse is over
C[x], because the binary linear determinant is a unit there, and
truncates over C(x); it therefore is polynomial over C[x]. The z and
x outputs finish the inverse. No specialization-only inverse is used.

**3.5, one phase entry and two affine powers.** Use slopes (0,1,m,m)
for x,y,z,w. Here W=2(p+q-1-beta), so the positive-weight endpoint is
beta=p+q-1. The two fibre cofactors survive every earlier wall. A line
apex other than x has an annihilator containing x and the other fibre;
their surviving minor excludes it. An apex plane is a graph over x,y;
only y can mix with fibres, because x has strictly smaller weight.
Its common coefficient is genuinely binary (3.2), or the nonzero old
x powers fix it to x. In the latter case the old-weight comparison
fixes the degree-m term and permits only higher powers in y. A line
apex gives a phase with one, two or three entries. Equal exponents are
removed at crossings; the other cases use 3.4, 3.3 or restart with a
strictly larger m. Both fibre weights minus beta strictly increase on
continuing segments, giving at most two shears. Finite degree and these
crossings prove termination before or at W=0.

**3.6, a minimum-weight unique apex.** Normalize its weight to one and
one phase entry per passive weight block. An absent phase variable
has a nonconstant affine coefficient, or there would be a constant
kernel and three apices. Two absent variables have distinct weights,
or their monomial coefficients would be proportional with the same
consequence. A small paired raise selects the highest phase power and
the absent affine partners. All passive weights become greater than
one, so exactly 3.3, 3.4 or 3.5 applies. Normalizing the greatest-weight
phase coefficient gives a polynomial coefficient and a nonnegative
active exponent. Laurent phase ratios were not silently made polynomial.

**3.7, three affine partners.** Raise all passive weights at slope one;
W decreases at slope two. Their distinct weights and surviving mixed
minors exclude every constant kernel. At a rank-two wall the unique
apex must be x, and its affine representation fixes every old partner,
so that no genuine wall exists. At rank three, one apex uses 3.6; a
two-apex plane must contain x, or its complementary plane would contain
an old nonzero mixed minor. Its common polynomial uses 3.2, or is x.
In the latter case a new nonlinear passive power is present and uses
3.5. The remaining endpoint is the positive-weight W=0 inverse.

**3.8, a pure power in a lower-weight variable.** Initially the weights
are (1,beta,M beta-p,M beta-q), with 0<q<p<M and p/M<=beta<1.
The possible zero fibre weight at the initial endpoint becomes positive
after any positive step. The w weight is always greater than beta and
becomes greater than one. The only crossings are
p/(M-1), (p+1)/M and 1; they are increasing for p<M-1 and simultaneous
for p=M-1. The fibre cofactors force a line or plane apex transverse
to the fibres. For a line apex, aligning its nonzero x component uses
only lower-old-weight variables. The y^M term fixes the highest phase
power, and the other phase exponents are p-(M-1)beta and q-(M-1)beta.
The second is negative, and the first is less than one, so only a
constant entry at the first crossing survives polynomiality. After its
shear the old form returns and has at least two apices, a contradiction.
For a plane apex use its graph shear. The common polynomial uses 3.2,
or is a power of x (at beta=1, of x+lambda y, rectified by an allowed
shear). Old and new weight comparisons fix R=a y^M and the two old
partners. Every continuing wall is one of the three crossings. The
path therefore exits or reaches W=0 with positive weights.

## 4. Proposition 4.1: all homogeneous top ranks

After removing affine terms, a degree-D>=3 top has singular Hessian.
The homogeneous low-rank classification gives essential rank one,
two or three. Rank zero is impossible. This partitions all potentials;
there is no generic-top or squarefree hypothesis.

At rank three, use h(x,y,z) and raise only w to its first wall. It occurs
by weight D-1, because a non-affine w monomial reaches weight D by
then. W remains at least 2D-4>0. The ternary cofactor survives. An apex
on w is excluded by complementary rank one; a single apex in the old
block uses 3.6. A two-apex plane involving w would make the old ternary
top affine in two coordinates with proportional homogeneous coefficients,
contradicting its rank. Thus the plane is in the ternary block. Its
cotangent coefficients use 3.2 or are powers of one base coordinate;
the remaining base second derivative is nonzero by rank three, so a
paired raise uses 3.5. Three apices would give a constant kernel,
excluded by the ternary cofactor and the occurrence of w.

At rank two, keep h(s,t) and raise z,u equally. Either W=0 is reached
at common weight D-1, or an earlier wall beta occurs. A constant
two-dimensional kernel would be the whole heavy block, contrary to
new terms. Rank two with one apex lies in the old block and reduces
to a s^D+b t s^(D-1)+c z s^(D-beta), with u unused. At rank three,
a unique apex is in the old block (3.6). Two apices cannot be the old
block, since the equal-degree affine coefficients would be proportional,
or the heavy block, since h has rank two. The mixed plane gives
G(s,z)+t A(s,z)+u B(s,z); common powers use 3.2 or the paired raise
3.5. Three apices leave a ternary form H(s,t,z) and unused u.

For either remaining unused-u form, raise only u. Its first wall gamma
satisfies beta<gamma<=D-1 and W>=2(D-beta-1)>0. With ternary H, its
cofactor excludes an apex on u and a constant kernel. A unique apex
in s,t uses 3.6. A unique apex on z gives H=Phi(z,L(s,t))+C(z)M(s,t);
C has positive weight and vanishes at z=0, making the old binary top
rank at most one, impossible. Two apices are s,t or z plus a base line;
the old ternary rank supplies the remaining nonlinear base derivative,
so 3.2 or 3.5 applies.

For the other unused-u form, the old mixed minors in (s,t),(s,z)
survive, and there is no constant kernel in those three variables.
A single apex must be s: for a base apex d_s e_s+d_t e_t, its
annihilator minor with z is -(d_t c q s^(q-1))^2, forcing d_t=0.
At rank three use 3.6; at rank two the three distinct affine powers
use 3.7. A two-apex plane cannot be z,u by the (s,t) minor; every
mixed plane must contain s by the same minor or the nonzero t
coefficient in the second directional derivative. The planes s,t and
s,z force the other coefficient to depend on s alone and give rank
at most two, a contradiction. The remaining plane s,u uses 3.2 or 3.5.

At rank one, retain a0 s^D and initially raise the three heavy weights
together. Any rank-three face containing that term is terminal as
follows. A minimum-weight unique apex uses 3.6. A heavier unique apex
would have a complementary phase normalized at s; its highest power
is D with coefficient a0. Coefficients of s^(D-1) times a heavier
phase variable have negative active weight and vanish. The image is
then a three-dimensional product of a plane with a curve, giving two
apices, a contradiction. A two-apex plane contains s because otherwise
its complementary space contains s, contradicting f_ss nonzero. Its
common polynomial uses 3.2, or is a power of s (3.5), or of the other
base variable t. In the last case raise s at speed one and both fibres
at speed D. Only a0 s^D and the two affine t powers have slope D.
After normalization this is exactly 3.8 with beta strictly between p/D
and one.

At a rank-one wall the constant kernel would be the whole heavy
block and contradict new terms. A rank-two wall, whether it has a
two-dimensional kernel or one apex, leaves a binary H(s,t) and two
equally weighted unused variables. A rank-three wall with three
apices leaves a ternary H and one unused variable of maximal weight.
Raising that variable gives rank three with no constant kernel, hence
one of the already terminal rank-three cases. From binary H, raise
the unused pair equally: rank two with two-dimensional kernel is
impossible, rank three either is already terminal or gives the ternary
case, and rank two with one apex gives two distinct affine partners
and one unused variable. Raise the last variable: rank three is terminal;
rank two gives the third distinct partner and 3.7. The W=0 endpoints
have positive weights throughout. This is a finite case graph, not an
indefinite search for a favourable weighting.

For degree two the gradient is already an invertible affine map; degrees
below two cannot have a nonzero Hessian determinant. All exits are actual polynomial inverses or constant-linear cotangent
forms with unit coefficient Jacobian. Restoring affine terms and undoing
the finite constant source changes preserves that assertion. In cotangent
form, replacing the fibre coordinates by grad R+DG^T fibres converts
the gradient to G times an identity map, up to an output exchange.
This proves the equivalence with JC2. All geometric classification
imports, polynomiality uses, endpoint cases and termination measures
have now been accounted for. Independent symbolic checks support the
displayed algebra; the argument in this note supplies the coverage.
