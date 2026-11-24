# Explicit Mathieu counterexamples over characteristic-zero fields

The main results concern standard rational unitary and sphere forms,
product instability over every formally real field, and every
nondegenerate quadratic complement of even dimension at least4 over
any characteristic-zero field. All phases and single fixed markers
are displayed explicitly. The interval--Laurent cancellation is
Long's known construction; the field and quadratic-complement
arguments below preserve its exact all-order identities.

## Exact assertions

Let K be ANY field of characteristic zero. Define all groups and varieties
below by scalar extension of the specified Q-models; this qualification
is part of the statement.

1. For n>=2, let E=Q[i]=Q[e]/(e^2+1), with conjugation e|->-e.
   Let U_n be the determinant-one unitary Q-group of the standard
   Hermitian form on E^n. Its real points give the compact SU(n),
   and its geometric group is SL_n. If M_11=a_11+e*b_11, put in
   O(U_n)[t,t^-1]

       u=a_11^2+b_11^2,  v=1-(1-u)^(n-1),
       f=(1-t^-1)((1-v)+v*t),  h=t^-1.                         (F1)

2. Let S=Spec Q[x,y,z]/(x^2+y^2+z^2-1), homogeneous under
   SO(x^2+y^2+z^2), and put in O(S)[t,t^-1]

       v=(1+z)/2,
       f=((1+z)*t-2*z+(z-1)*t^-1)/2,  h=t^-1.                  (F2)

3. Let H_Q=(-1,-1)_Q be the Hamilton quaternion Q-algebra. On the
   unit group G=GL_1(H_Q)=Spec Q[a,b,c,d,N^-1], put

       N=a^2+b^2+c^2+d^2, alpha=a^2+b^2, beta=c^2+d^2,
       F=(N-1)*(alpha*N+beta)/N^2,  Hmark=N^-1.                 (F3)

For their unique normalized invariant K-linear traces, after EVERY
extension Q->K, all these ACTUAL phases and SINGLE fixed markers satisfy

    Lambda_K(f^m)=0,
    Lambda_K(h*f^m)=(-1)^(m-1)/(m+1) !=0, EVERY m>=1,           (F4)

with F,Hmark replacing f,h in case3. Thus each stated product/group
has a non-Mathieu trace kernel over EVERY characteristic-zero K.

If K is FORMALLY REAL, the trace kernels on U_n,K alone and on S_K
alone are Mathieu: each has ONLY the zero pure element. The Laurent
torus K[t,t^-1] with constant term has the Mathieu property over EVERY
characteristic-zero K. Consequently, over EVERY formally real K, two
individual Mathieu trace kernels can have a non-Mathieu tensor-product
trace kernel. In particular K=Q already gives these product failures.

For formally real K, GL_1(H_Q tensor K) has no nonidentity K-rational
unipotent elements, while its trace kernel still fails Mathieu by(F3)--(F4).

No compactness of groups of K-points is asserted for a general K.
No embedding of the whole coefficient field is assumed. Only finitely
generated coefficient fields are embedded in C for the scalar-torus
comparison; no ordering of K is chosen.
The assertion is about these standard rational forms; arbitrary twists
or arbitrary real-reductive homogeneous spaces are not classified.

## Real calculations used in the rational argument

Sections1--6 fix the real field and establish the complete real
moment/positivity lemmas. Sections7--12 prove the standard rational
form and arbitrary-field extensions. Sections13--20 give the
independent arbitrary even-quadratic construction. Section21 compares
the routes without circular dependencies.

## Exact conclusion

For every n>=2, let K_n be the compact real algebraic group SU(n), and put

    A_n=R[K_n][t,t^-1]=O(K_n x G_m,R).

Let Lambda_n be the unique normalized invariant real Reynolds functional.
Write a_11,b_11 for the real and imaginary entry coordinates of a matrix in
K_n, and define the REAL regular functions

    u=a_11^2+b_11^2,
    v=1-(1-u)^(n-1),
    f=(1-t^-1)((1-v)+v t),
    h=t^-1.                                                   (R1)

Then, for EVERY integer m>=1,

    Lambda_n(f^m)=0,
    Lambda_n(h f^m)=(-1)^(m-1)/(m+1) != 0.                     (R2)

Thus ker Lambda_n is not a Mathieu--Zhao subspace, even though the
derived real group is compact. Both the compact group factor alone
and the split torus alone have the real Mathieu property.

A smaller homogeneous example is the compact real sphere

    X={x^2+y^2+z^2=1}=SO(3)_compact/SO(2)_compact.

On X x G_m, put

    f=((1+z)t-2z+(z-1)t^-1)/2,  h=t^-1.                       (R3)

The same identities (R2) hold. Hence adjoining a split torus to a
compact homogeneous variety can destroy the real Mathieu property.

There is also an explicit group example on the real algebraic group of
nonzero Hamilton quaternions H^times, an inner real form of GL_2.
Its derived group is compact SU2 and it has no nonidentity real
unipotent elements; nevertheless its real normalized trace is not
Mathieu. Section5 gives regular formulas and a scheme-theoretic
finite etale double cover proving their exact moments.

These are applications of Long's existing mixed interval--Laurent
construction to REAL coordinate algebras, with a regular polynomial
change of the compact coordinate. The interval--Laurent pair and its
cancellation are credited to Long; no novelty claim for them is made.

## 1. Regular uniform-variable criterion

Let B be any commutative unital real algebra, and let lambda:B->R be
linear with lambda(1)=1. Suppose v in B satisfies, for EVERY polynomial
p in R[s],

    lambda(p(v))=integral_0^1 p(s) ds.                         (U1)

On A=B[t,t^-1] define Lambda by applying lambda to the constant
Laurent coefficient. Set f,h as in (R1), using this v. Then (R2) holds.

Proof. The hypothesis lets us integrate in a scalar s in place of v.
For every m>=1, an exact Laurent-polynomial identity is

    integral_0^1 [(1-t^-1)((1-s)+s t)]^m ds
       =(t-1)^(m-1)(t-t^-m)/(m+1).                            (U2)

To check it, integrate (1+s(t-1))^m, multiply by (t-1)^m t^-m,
and cancel its factor t-1. Equivalently one may use the binomial
theorem, so no specialization t=1 or analytic branch is involved.

The first term on the right of (U2) has only positive t-exponents,
and the second has only negative t-exponents. Its constant term is
therefore zero. Its coefficient of t^1 is exactly
(-1)^(m-1)/(m+1), because only t(t-1)^(m-1) contributes. Taking the
constant term after multiplication by the SINGLE fixed h=t^-1
proves (R2). This proves all orders, not a finite purity prefix.

For a homogeneous B=O(Y) with normalized reductive trace lambda, this
Lambda is precisely the normalized invariant functional on
Y x G_m: it is invariant under the acting group and under G_m, has
value1 on1, and uniqueness identifies it.

More generally, Long's parameter family works over R: for any integer
d>=1 and nonzero real alpha,beta, set

    f_d=alpha(1-beta t^-d)((1-v)+beta^-1 v t^d), h_d=t^-d.

Then every pure moment is zero and

    Lambda(h_d f_d^m)=(-1)^(m-1) alpha^m beta^-1/(m+1)

for all m>=1. This is the realization of Long Proposition3.1, rather
than a separate new interval cancellation.

## 2. A polynomial beta-distribution uniformizer

Suppose a real regular u has the beta(a,b) distribution under the
normalized compact invariant measure, where a,b are POSITIVE INTEGERS.
Write

    C=(a+b-1)!/((a-1)!(b-1)!),
    V(s)=C integral_0^s r^(a-1)(1-r)^(b-1) dr.

Then V is a polynomial in Q[s], V(0)=0, V(1)=1, and V'(s)>0 for
0<s<1. For every polynomial p, substitution y=V(s) gives

    integral_0^1 p(V(s)) C s^(a-1)(1-s)^(b-1) ds
       =integral_0^1 p(y) dy.

Therefore v=V(u) is a REGULAR uniform variable satisfying (U1), and
the criterion produces a pure unsafe phase on the product with G_m.
Integer parameters matter here: they make V a polynomial and retain
regularity. No assertion for a general beta distribution is needed.

## 3. Compact SU(n): exact distribution and algebraic trace

View K_n=SU(n) as a real algebraic group by using real and imaginary
matrix-entry coordinates. Its real points are compact, and its
complexification is SL_n,C. Haar integration of a REAL regular
function is real. It is the normalized invariant Reynolds functional:
invariance under K_n(R), together with Zariski density and rational
local finiteness, is scheme-theoretic K_n-invariance; uniqueness then
applies. The same identification holds for the homogeneous sphere
in the next section.

There is an elementary uniqueness argument on these real examples.
A regular function belongs to a finite-dimensional rational module of
translates. Compact Haar averaging in that finite-dimensional module
gives its constant orbit-average. Any normalized invariant linear
functional commutes with this finite-dimensional vector integral, so
must take that same constant value. On the product with G_m, invariance
under t|->r t for real r!=0 kills every nonzero Laurent weight; on
weight zero the compact argument applies. This directly proves
existence and uniqueness of the normalized trace used here.


For completeness, density here is not a rational-point transitivity
assumption. A smooth real algebraic group has a real analytic
neighborhood of its identity of real dimension dim K_n; a regular
function vanishing on all real points vanishes on that neighborhood
and therefore in the irreducible algebraic group. Invariance of a
linear functional can be tested on one finite-dimensional rational
submodule, where its coefficient equations are regular functions.

SU(n) acts transitively on the unit sphere in C^n for n>=2. Indeed,
extend a unit vector to a unitary basis and adjust the phase of a
different column so that the determinant is1. Consequently the first
column of a Haar SU(n) matrix is the uniform unit vector: its law is
invariant, and a transitive compact action has a unique invariant
probability measure.

Here is a direct computation of u=|M_11|^2. Let Z_1,...,Z_n be
independent standard complex Gaussian variables. The direction
Z/||Z|| is uniform on that sphere by unitary invariance. The independent
variables E_i=|Z_i|^2 have density exp(-r) on r>0. Put
S=E_2+...+E_n. Its density is s^(n-2)exp(-s)/(n-2)!.
The change of variables

    u=E_1/(E_1+S), R=E_1+S

has Jacobian R. Integrating R gives the exact density

    (n-1)(1-u)^(n-2) du,  0<u<1.                              (U3)

Thus u is beta(1,n-1). Its polynomial cumulative distribution is

    V(u)=1-(1-u)^(n-1).

Section2 proves that v=V(u) is uniform, and Section1 proves (R1)--(R2).
All coefficients in the displayed formula are real; no i-valued scalar
or complex-valued regular phase is being smuggled into the real algebra.

The classical Haar column construction is also explicitly stated in
Meckes, Section1.2, printedpp19--20, including the SU(n) last-column
adjustment. The derivation of (U3) above supplies the needed density
rather than assuming an asymptotic distribution or a numerical model.

## 4. Compact S^2 and individual-factor Mathieu properties

Normalized surface area on S^2 is the normalized SO(3)-invariant
measure. In cylindrical coordinates its height z has density dz/2
on [-1,1]. Thus v=(1+z)/2 is uniform on [0,1], and the uniform-variable
criterion gives exactly (R3) and (R2).

For the compact factor alone, the real invariant functional is
positive definite on squares: a nonzero real regular function is
nonzero on an open subset of its real locus, so its squared Haar
integral is strictly positive. Real points are Zariski dense for
SU(n) and S^2. A pure real function has zero second moment, so it
must be zero. Therefore every pure real function on either compact
factor is safe.

For the split torus alone, the real constant-term functional has
the Mathieu property by the classical Laurent theorem of
Duistermaat--van der Kallen (their introduction printedp222 and the
rank-one conclusion in Remark3 printedp223; arbitrary-rank proof
in Section3). Apply the complex Laurent theorem to a real phase and
then restrict its conclusion to real markers. Thus (R1) and (R3)
give explicit failure of closure under products, even when BOTH
real factors separately have the Mathieu property.

This does not say every compact factor becomes bad after adjoining
G_m. A compact REAL TORUS times a split torus is still a real torus;
after complexification the Laurent theorem applies. Nor does this
note give a classification for all real reductive homogeneous spaces.

## 5. A familiar real reductive group: invertible quaternions

Let H be the Hamilton quaternion algebra over R, with coordinates
q=a+bi+cj+dk and norm

    N=a^2+b^2+c^2+d^2.

Its unit group is the connected reductive real algebraic group

    G=GL_1(H)=Spec R[a,b,c,d,N^-1].

For later clarity put alpha=a^2+b^2 and beta=c^2+d^2, so N=alpha+beta.
Define the REAL REGULAR phase and SINGLE fixed marker

    F=(N-1)(alpha*N+beta)/N^2,
    Hmark=N^-1.                                              (Q1)

Negative powers of N are regular on G by definition. Then, for every
m>=1, its normalized invariant real Reynolds functional satisfies

    Lambda_G(F^m)=0,
    Lambda_G(Hmark F^m)=(-1)^(m-1)/(m+1) !=0.                  (Q2)

### The cover is finite etale as a scheme

Let K=SL_1(H)={q:N(q)=1}, the compact real SU2, and consider

    p:K x G_m -> G,   (q,s) |-> s q.

Scalars are central and quaternion norm is multiplicative, so p is
a group homomorphism. At the level of coordinate algebras there is
the explicit isomorphism

    O(K)[s,s^-1]
      = O(G)[s]/(s^2-N).                                     (Q3)

On the right, s is already invertible because s^2=N is invertible.
The inverse coordinate formulas are q_i=a_i/s. The norm-one relation
becomes sum(q_i^2)=N/s^2=1, proving (Q3) without a pointwise argument.

Thus p is finite locally free of rank2 and faithfully flat. It is
etale because the derivative 2s is a unit. Its kernel is the diagonal
mu_2={(1,1),(-1,-1)}, and it is a central Galois isogeny. In particular,
this is a scheme cover over R; no assertion about quotient sets of
real points substitutes for (Q3).

### Exact trace and moment computation

The normalized trace upstairs is Haar_K times CT_s. Under p,

    p^*N=s^2,
    p^*(alpha/N)=u=q_a^2+q_b^2.

The variable u is uniform on [0,1] by the SU2 calculation in Section3.
Direct substitution into (Q1) gives

    p^*F=(1-s^-2)((1-u)+u s^2),
    p^*Hmark=s^-2.                                           (Q4)

The first is the d=2 parameter specialization of Section1. Therefore
its pure moments vanish and its marked moments are the numbers in
(Q2).

Finally the trace upstairs restricts to the trace downstairs:
p is equivariant for the central isogeny K x G_m -> G, so the
restriction is normalized and G-invariant by faithfully flat
surjectivity; uniqueness identifies it with Lambda_G. As a result the upstairs computation proves (Q2). The two functions
in (Q4) are invariant under (q,s)->(-q,-s), as their explicit
downstairs formulas already demonstrate. No phase averaging is used.

### Fixed-phase trace, purity and safety in both directions

The double cover also gives an exact statement for EVERY GIVEN
downstairs phase, not just F. Write C=O(G) and D=O(K)[s,s^-1],
so D=C direct_sum C s with s^2=N. The deck automorphism tau fixes C
and sends s to -s. For b=b0+b1 s define

    T_p(b)=(b+tau(b))/2=b0.

In the C-basis 1,s, multiplication by b has matrix

    [[b0,N b1],[b1,b0]].

Hence T_p=(1/2)Tr_(D/C), the normalized algebra trace. The upstairs
functional is tau-invariant: tau is translation by the central point
(-1,-1) of K x G_m. If g in C is ANY FIXED phase and b in D is ANY
SINGLE FIXED marker, averaging b and tau(b), and using trace restriction,
gives for ALL m>=0

    Lambda_up(b (p^*g)^m)
       =Lambda_down(T_p(b) g^m).                             (Q5)

The marker T_p(b) is independent of m, with exactly the same complete
set of surviving orders. Conversely a downstairs marker a uses p^*a,
and T_p(p^*a)=a.

Taking b=1 shows separately that g pure implies p^*g pure and that
p^*g pure implies g pure. If g is safe, every fixed upstairs b is
controlled by the fixed downstairs T_p(b), proving p^*g safe.
Conversely, if p^*g is safe, pull each fixed downstairs a up to p^*a
and use (Q5), proving g safe. The cutoff may depend on each marker.
A given unsafe marker pushes or pulls with identical surviving orders,
so unsafety is equivalent as well. These are fixed-phase implications
under this ACTUAL cover: an arbitrary upstairs phase is not claimed
to descend. No converse from purity of a power to purity of its root
is asserted.

### Group structure and the real-form distinction

The explicit complex matrix representation

    a+bi+cj+dk |-> [[a+bi,c+di],[-c+di,a-bi]]

identifies H tensor_R C with Mat_2(C), with determinant N.
Consequently G_C=GL_2,C and K_C=SL_2,C. Reductivity and connectedness
of G follow from this complexification. The derived subgroup of G is
K, as can be checked after this faithfully flat extension using
[GL_2,C,GL_2,C]=SL_2,C. Its real points are compact.

There are no nonidentity real unipotent elements in G. For the
displayed matrix, unipotence forces its trace2a=2 and determinant N=1.
Thus a=1 and b^2+c^2+d^2=0 over R, so b=c=d=0. In particular this
example does not obtain its real failure from a real unipotent root
subgroup; the split central direction and the compact uniform
coordinate already suffice.

This does not classify inner forms or all anisotropic-derived real
groups. It gives one completely explicit real counterexample and a
family arising from the general regular-uniform-variable criterion.

## 6. Scope and exact sources

The real Mathieu property is not determined by complexification alone.
For example compact real SU2 has no nonzero pure REAL regular function
by positivity, while the split real SL2 has real-coefficient bad phases.
Long's Theorem4.2 supplies F=(1+c)(ad+b), H=-c in the split coordinate
algebra R[a,b,c,d]/(ad-bc-1); uniqueness and scalar extension of the
normalized invariant trace identify its moments with the stated complex
Haar computation. The same expression is complex-valued in the compact
matrix-entry realization; this does not contradict real positivity.

The present uniform-variable criterion gives explicit real product
failures and the quaternion-group example. It is not a classification
of all real reductive homogeneous spaces, does not assert that all
compact factors become bad after adjoining G_m, and does not transfer
the global Mathieu property merely from a complexification. The
interval--Laurent construction is Long's known method, with its
all-order cancellation reproduced in Section1.

Actual primary/author sources inspected:

1. C.D.Long, Counterexamples to the xz-Conjecture and the Mathieu
   Conjecture for SU(2), arXiv:2607.19012v1, Theorem2.1 and its full
   proof; Proposition3.1; Lemma4.1 for the uniform SU2 coordinate.
   https://arxiv.org/html/2607.19012v1
   This supplies the exact credited mixed interval--Laurent construction.
   Theorem4.2 is used only for the split-versus-compact comparison.
   No claim about the Jacobian conjecture from its introduction is used.

2. E.Meckes, The Random Matrix Theory of the Classical Compact Groups,
   author PDF, Section1.2, printedpp19--20, explicit column construction
   of Haar measure for U(n),SU(n),SO(n).
   https://case.edu/artsci/math/esmeckes/Haar_book.pdf
   Only that narrow statement is used as a literature cross-check;
   Section3 gives its own transitivity and Gaussian proof of (U3).

3. J.J.Duistermaat and W.van der Kallen, Constant terms in powers of
   a Laurent polynomial, Indag.Math.(N.S.)9(1998),221--231,
   printedpp222--223 and Section3.
   https://webspace.science.uu.nl/~kalle101/powers.pdf
   Used only for the separate torus-positive comparison.

No finite purity test is used as an inference to (R2). The optional
exact controls are supporting checks only. No claim that this
real-form application is absent from all literature has been established.

## 7. Rational traces and arbitrary scalar extension

The Q-groups in the statement are connected reductive: after an
algebraic closure the unitary group becomes SL_n, the quaternion unit
group becomes GL_2, and the sphere is a homogeneous space of SO_3.
Their coordinate algebras have a rational locally finite action and
their invariant functions are Q.

Here is the exact classical trace construction. In characteristic zero,
every finite-dimensional rational module of a connected reductive group
is semisimple, over the given coefficient field. For a nontrivial simple
module V, the Lie-image gV is a G-stable subspace. If it were zero, the
connected image group would have zero differential and therefore be
trivial in characteristic zero. Thus gV=V. On trivial modules g acts
as zero. Rational local finiteness and A^G=Q consequently give

    A_Q=Q*1 direct_sum g_Q A_Q.                                (F5)

Projection defines the unique normalized invariant Q-linear Lambda_Q.
The sum g_Q A_Q is a VECTOR subspace, not an ideal.

For ANY extension field K/Q, Lie algebras and derivations base-change:

    (g_Q A_Q) tensor_Q K = g_K A_K.

Tensoring(F5) therefore gives the same normalized invariant projection
over K:

    Lambda_K=id_K tensor Lambda_Q.                            (F6)

This proves extension compatibility for arbitrary K, including fields
too large to embed in C. It is not an assertion that the GLOBAL Mathieu
property is preserved by arbitrary coefficient extension.

The narrow classical source is Milne, Algebraic Groups(2017 author PDF),
Theorem22.42/Corollary22.43 printedpp477--478 and Section12l,
Remark12.53/Reynolds paragraph printedpp248--249. These give
same-field complete reducibility and rational locally finite decomposition.
The elementary Lie-image and scalar-extension steps have been written
above. No all-automorphism theorem is required for this argument.

## 8. The explicit moment identities descend to Q, then extend to K

After Q->R, the normalized traces on U_n,R and S_R are the normalized
compact Haar traces. The real proof gives:

- u on compact SU(n) has the EXACT beta(1,n-1) law.
- v=1-(1-u)^(n-1) is a polynomial regular uniform variable.
- v=(1+z)/2 is a polynomial regular uniform variable on the sphere.
- For every m>=1, the exact Laurent identity is

    integral_0^1 [(1-t^-1)(1-s+s*t)]^m ds
       =(t-1)^(m-1)*(t-t^-m)/(m+1).                            (F7)

The constant term of(F7) is0 and its t^1 coefficient is
(-1)^(m-1)/(m+1). This is Long's credited mixed interval--Laurent
method, realized by REAL regular variables.

Every phase, marker and trace in these calculations is defined over Q.
For each fixed m, the corresponding trace value in A_Q belongs to Q.
By(F6), its image in R is exactly the computed rational number in(F4).
Faithfulness of Q->R proves the equality IN Q. Then(F6) transports it
to every characteristic-zero K, where the same nonzero rational number
remains nonzero. The argument works for every integer m; it is not
a bounded moment experiment.

For the quaternion example, the real proof uses the explicit central
finite etale isogeny

    SL_1(H_Q) x G_m -> GL_1(H_Q), (q,s)|->s*q,
    O(cover)=O(base)[s]/(s^2-N).

The ring identity and derivative2s-unit argument are over Q.
The phase/marker pull back to the d=2 Long family

    (1-s^-2)(1-u+u*s^2), s^-2,

with u uniform on compact SU2 after Q->R. The direct equivariant
isogeny trace restriction proves the real formulas. The same rational
trace argument just given then proves(F4) for(F3) over every K.
No compact measure over K, descended K-cover point, or K-point
transitivity is used.

Equivalently the word 'uniform' over an arbitrary K can be replaced
by its algebraic moment statement

    Lambda_K(v^j)=1/(j+1), ALL j>=0.

The formal functional J_K:K[s]->K with J_K(s^j)=1/(j+1) replaces the
real interval integral. Long's binomial cancellation and(F7) are then
identities over K, with all denominators invertible in characteristic zero.

## 9. A rational Gram lemma for every formally real field

A field K is formally real here if -1 is not a finite sum of squares
in K. We use this definition directly, not the theorem that such a
field admits an ordering.

**Gram lemma.** Let B_Q be a Q-algebra with Q-linear lambda_Q.
Suppose its scalar extension lambda_R is strictly positive on squares
of nonzero elements of B_R=B_Q tensor R. Then, for EVERY formally
real K and every f in B_K,

    lambda_K(f^2)=0  ==>  f=0.                                 (F8)

Proof. Write f as a finite sum

    f=sum_(i=1)^r c_i a_i,

where the a_i in B_Q are Q-linearly independent and the c_i lie in K.
Such a finite independent set is obtained by choosing a Q-basis of the
span of the finitely many coefficients from B_Q used to express f.
Field extension preserves this independence.

Form the finite rational symmetric Gram matrix

    G_ij=lambda_Q(a_i a_j).

For every nonzero real vector x, independence after Q->R makes
sum x_i a_i a nonzero element of B_R. Thus x^T G x>0.
Consequently G is a positive definite rational matrix.

An exact rational LDL^T decomposition gives

    c^T G c=sum_i d_i ell_i(c)^2,

where each d_i is a positive rational number and the ell_i form an
invertible rational triangular linear transformation. This decomposition
needs no irrational square root: eliminate the first row/column, whose
positive rational pivot leaves a positive definite rational Schur
complement, and iterate.

Write d_i=A_i/B_i with POSITIVE INTEGERS A_i,B_i. Then d_i is literally
a sum of A_i*B_i copies of the rational square (1/B_i)^2. Therefore
lambda_K(f^2)=c^T Gc is a FINITE SUM OF SQUARES in K.

In a formally real field, a sum of squares can vanish only if every
term vanishes. Otherwise divide by one nonzero term's square to
express -1 as a sum of squares. Hence every ell_i(c)=0; invertibility
gives c=0, proving f=0. This proves(F8) without ordering K or embedding
it into an Archimedean field.

For B_Q=O(U_n) and B_Q=O(S), the positivity hypothesis holds by the
compact real Haar interpretation and Zariski density of the real locus.
A nonzero real regular function is nonzero on a real open subset;
its squared Haar integral is strictly positive. The rational trace
construction identifies this Haar functional with lambda_Q tensor R,
so the Gram lemma applies.

If f is pure over formally real K, its second moment is zero.
Equation(F8) forces f=0. Thus both standard-form compact-factor
trace kernels have the Mathieu property over every formally real K.
The term 'compact-factor' refers to their REAL realizations, not to
any claimed topological compactness of K-points.

## 10. The torus factor over every characteristic-zero field

Let f,h lie in K[t_1^+-1,...,t_d^+-1], and suppose every positive
constant term of f^m vanishes. The finitely many coefficients of
this GIVEN pair f,h generate a finitely generated field K0/Q.

The field K0 embeds in C: choose an algebraically independent complex
image of a finite transcendence basis and then extend the embedding
through the finite algebraic extension, using algebraic closedness of C.
This is a statement about K0, not the whole coefficient field K.

All pure constant terms lie in K0 and are zero there. Their complex
images therefore satisfy the Duistermaat--van der Kallen theorem.
For this fixed marker h it follows that CT(h*f^m)=0 for sufficiently
large m in C. Injectivity descends those equalities to K0 and hence
to K. The cutoff may depend on f and h, exactly as safety permits.

Thus the torus trace kernel is Mathieu over every char0 K.
Combined with Section9 and the explicit all-field bad product pair,
this proves the claimed product instability for EVERY formally real K.

The exact imported Laurent theorem is Duistermaat--van der Kallen,
Constant terms in powers of a Laurent polynomial, printedpp222--223
and Section3. It is applied only to ordinary scalar Laurent coefficients
after embedding the finitely generated coefficient field, never to
function-valued Laurent coefficients.

## 11. Quaternion unipotents over formally real K

For q=a+bi+cj+dk, the usual matrix after adjoining a square root of -1 is

    [[a+ib,c+id],[-c+id,a-ib]],

with trace2a and determinant N. A K-rational unipotent element must have
a=1 and N=1. Therefore b^2+c^2+d^2=0. Formal reality forces
b=c=d=0, so q=1. This proves absence of nonidentity K-rational
unipotents without choosing an ordering or using a compactness claim.

Its derived group is SL_1(H_Q tensor K), whose real realization is
compact SU2; the precise all-field description is this rational inner
form. The all-field failure(F3)--(F4) coexists with the absence of
K-unipotents for formally real K.

## 12. Evidence, exact sources and boundaries

Sections1--6 above give the complete real moment and positivity proofs.
Sections7--11 establish the rational/base-change and formally real extensions.

Actual narrow sources already inspected:

- Milne, Algebraic Groups(2017 author PDF),22.42--22.43 printed477--478,
  12l/12.53/Reynolds printed248--249:
  https://www.jmilne.org/math/Books/iAG2017.pdf
  retained primary SHA256
  503f1a7d511632af05541c52fb2bf6fc36f5a75b0323f2f9110f03dcce11fd67.
  This is the necessary rational trace/complete-reducibility input.

- Long, arXiv:2607.19012v1,Theorem2.1 complete beta-binomial proof and
  Proposition3.1:
  https://arxiv.org/html/2607.19012v1
  retained HTML SHA256
  abeec76e95cba0e6244f4af31240634e28d657116c8afcf743b58db144dd7ebd.
  The interval--Laurent construction is known and credited.

- Duistermaat--van der Kallen, Indag.Math.(N.S.)9(1998),221--231,
  printed222--223 and Section3:
  https://webspace.science.uu.nl/~kalle101/powers.pdf
  retained PDF SHA256
  bb002d12b539a1ff3fdd099f343d3accd921ccf60f8a149f29b07cd0f7e87e1b.

The Gram proof uses only finite rational linear algebra and the
DEFINITION of formal reality. It does not import Artin--Schreier,
an ordering theorem, Lagrange's four-square theorem or an embedding
of a general K into R/C.

These results do not establish a Mathieu classification over arbitrary
fields or for arbitrary twists. The global Mathieu property is not
being propagated by base extension; only GIVEN rational trace identities
are, with the separate Gram argument proving the individual-factor
property for formally real fields.


## 13. Every nondegenerate even-dimensional quadratic complement

Let K be ANY field of characteristic zero. Let q be a nondegenerate
quadratic form on a K-vector space V of EVEN dimension N=2r>=4.
Choose an orthogonal nondegenerate decomposition

    q=A+B,  rank A=2, rank B=2r-2.

Such a decomposition exists by diagonalizing q in characteristic !=2.
Put X={x in V:q(x)!=0}, the AFFINE principal open with coordinate ring
C=K[V,q^-1]. Let Lambda_X be the unique normalized invariant trace for
the geometrically transitive action of SO(q) x G_m given by
(g,t)x=t*g*x. Define ACTUAL K-regular functions

    u=A/q,
    v=1-(1-u)^(r-1)=1-(B/q)^(r-1),
    F=(1-q^-1)((1-v)+v*q),  H=q^-1.                            (E1)

Then for EVERY integer m>=1,

    Lambda_X(F^m)=0,
    Lambda_X(H*F^m)=(-1)^(m-1)/(m+1) !=0.                      (E2)

In particular the normalized trace kernel is not Mathieu on X, for
EVERY nondegenerate q of every EVEN dimension>=4 over EVERY char0 K.
There is ONE fixed marker H, and it survives at EVERY positive order.
Every assertion persists after EVERY coefficient-field extension,
with the displayed base-changed q,A,B,F,H. No K-point of {q=1},
embedding of K in R/C, positive definiteness of q, or splitness is
assumed.

For r=2 the formula becomes

    F=(q-1)*(A*q+B)/q^2,

recovering the quaternion norm formula when q is the four-square norm.
For general r an equivalent formula is

    F=(q-1)-(q-1)^2*B^(r-1)/q^r.                               (E3)

Thus all denominators are powers of the explicitly inverted q.

This is an application of Long's credited interval--Laurent cancellation,
transported through a regular beta-CDF variable and an actual finite
etale double cover. No novelty claim for that cancellation is made.

## 14. Geometrically homogeneous spaces and normalized traces

Write W={w:q(w)=1}. It is smooth affine and geometrically integral.
The group SO(q) is connected reductive for N>=4 and becomes transitive
on W over an algebraic closure. W(K) may be empty.

The product Y=W x G_m is geometrically homogeneous under the connected
reductive group SO(q) x G_m. So is X: geometrically, scale one nonzero
q-value to another by a square root, then use transitivity on a fixed
nonzero quadric. This is a geometric/scheme statement, not an assertion
that SO(q)(K) x K^* is transitive on X(K).

For completeness, a rational locally finite action of a connected
reductive group in characteristic zero on such an algebra A has

    A=K*1 direct_sum g*A,

because A^G=K and finite-dimensional rational modules are semisimple.
Each nontrivial simple module is its Lie-image; a zero differential
would make the connected characteristic-zero image trivial.
Projection defines the unique normalized invariant trace. The
Lie-image is a VECTOR subspace, not an ideal.

This construction commutes with every field extension:

    Lambda_L=id_L tensor_K Lambda_K.

It therefore applies to W even without K-points. The narrow complete
reducibility input is Milne22.42--22.43, printed477--478, with local
finiteness/Reynolds in12l/12.53, printed248--249.

On Y the trace is Lambda_W times CT_s, by torus weight decomposition
and normalization. The source and target traces used below are thus
defined algebraically over K, including for a q having no real
norm-one point.

## 15. The exact beta moments for the rank-two summand

On W put a=A(w). We claim, for every k>=0,

    Lambda_W(a^k)=k!*(r-1)!/(k+r-1)!.                          (E4)

Here is a proof with no embedding assumption on K.

First take the STANDARD Q-model

    W_std={y_1^2+...+y_(2r)^2=1},
    a_std=y_1^2+y_2^2.

After Q->R, the normalized trace on this sphere is normalized surface
area: it is SO(2r)-invariant and normalized, and uniqueness identifies
it. Pair the 2r independent real Gaussian coordinates into r complex
Gaussian pairs. The sum of the squares in each pair has the independent
exponential law. Normalizing the Gaussian vector gives uniform sphere
direction, and the ratio

    a_std=E_1/(E_1+...+E_r)

has the EXACT beta(1,r-1) density

    (r-1)*(1-a)^(r-2) da, 0<a<1.

This follows directly from the ratio/total Jacobian with E_1 exponential
and E_2+...+E_r gamma(r-1), exactly as in the SU(n) calculation.
Integrating a^k gives(E4), a rational number.

The normalized trace on W_std is defined over Q and commutes with
Q->R. Faithful extension therefore proves(E4) over Q itself and
then over any characteristic-zero extension.

Now extend the ORIGINAL K to an algebraic closure. Orthogonal
diagonalization and square roots of the nonzero diagonal coefficients
identify the pair(q,A) with the standard pair(sum of2r squares,
first2 squares). The change of variables conjugates SO(q) to SO(2r),
so equivariance and normalized-trace uniqueness transport(E4).
The standard model over this algebraic closure already has the rational
moments just proved. Faithfulness of K->Kbar descends those identities
to K. The argument does not embed K or Kbar in C.

## 16. The polynomial uniform variable

The beta(1,r-1) cumulative distribution is the POLYNOMIAL

    V(a)=1-(1-a)^(r-1).

Over the standard real model, substitution using
V'(a)=(r-1)*(1-a)^(r-2) proves

    Lambda_W(p(V(a)))=integral_0^1 p(t)dt

for every real polynomial p. Equivalently,

    Lambda_W(V(a)^j)=1/(j+1), EVERY j>=0.                       (E5)

These are rational polynomial consequences of(E4), so the same
rational-trace descent and equivariant transport prove(E5) for the
original arbitrary K and q. Over K the interval integral is only
shorthand for the formal K-linear map J_K(t^j)=1/(j+1); no ordered
or analytic interpretation of K is used.

The integer r-1>=1 is essential to this construction: V is a
polynomial, so it gives a regular function. On X the corresponding
regular variable is exactly v in(E1).

## 17. The finite etale double cover and actual phase descent

Define

    p:Y=W x G_m -> X,  (w,s)|->s*w.

Its coordinate ring has the explicit presentation

    O(Y)=C[s]/(s^2-q),                                        (E6)

with inverse w_i=x_i/s. Because q is a unit, so is s, and
q(w)=q(x)/s^2=1. This proves(E6) as a ring isomorphism.
Consequently p is finite free of rank2 and faithfully flat;
the derivative2s is a unit, hence p is etale. Its deck action is

    tau(w,s)=(-w,-s).

For EVEN N, -I belongs to SO(q), so tau is action by the central
point(-I,-1) of SO(q) x G_m. In particular it preserves the source
normalized trace.

The pullback is equivariant for the stated SO(q) x G_m actions,
so the normalized source trace restricts to the unique normalized
target trace:

    Lambda_Y(p^*c)=Lambda_X(c), c in C.                         (E7)

This direct equivariance argument suffices; an additional general
finite-etale trace theorem is not needed to prove(E2).

The ACTUAL explicit phase and marker satisfy

    p^*u=A(w)=a,
    p^*v=V(a),
    p^*q=s^2,
    p^*F=(1-s^-2)((1-V(a))+V(a)*s^2),
    p^*H=s^-2.                                                (E8)

They are already tau-invariant: A(w) is quadratic, so A(-w)=A(w),
and s^2 is invariant under s|->-s. Thus the displayed functions
really descend; no averaging of a pure phase is used.

By(E5), Long's all-order Laurent identity applies to the source,
with the substitution t=s^2. Laurent constant term is unchanged by
this substitution. Explicitly,

    J_K([(1-t^-1)((1-v)+v*t)]^m)
       =(t-1)^(m-1)*(t-t^-m)/(m+1).

Its constant term is0 and its t^1 coefficient is
(-1)^(m-1)/(m+1). This proves the source pure and marked moments in
(E2). Equation(E7) proves the identical target moments. Every positive
order has been covered algebraically.

## 18. Given-phase marker and safety transport under this cover

The extension(E6) has basis1,s. If b=b0+b1*s, then its normalized
algebra trace is

    T_p(b)=Tr_(O(Y)/C)(b)/2=b0=(b+tau(b))/2.

For EVERY GIVEN downstairs phase g in C and SINGLE fixed upstairs
marker b, deck invariance and(E7) give for ALL m>=0

    Lambda_Y(b*(p^*g)^m)=Lambda_X(T_p(b)*g^m).                  (E9)

The pushdown marker is independent of m and has exactly the same
complete set of surviving orders. Conversely a downstairs marker c
pulls to p^*c, and T_p(p^*c)=c.

Taking b=1 proves both purity directions for this GIVEN g.
Safety downstairs implies safety upstairs because every fixed b
pushes to fixed T_p(b). Safety upstairs implies safety downstairs
by pulling each fixed c up. The cutoff may depend on the marker.
Given fixed unsafe witnesses push and pull with identical surviving
orders; infinitely many survivors force a nonzero pushed marker.

These are implications for phases that ACTUALLY descend through p.
They do not make an arbitrary upstairs phase invariant and do not
assert that purity of a power implies purity of its root.

## 19. The odd-dimensional boundary of this construction

For odd dimension N=2r+1>=3, the rank-two ratio has beta(1,r-1/2)
moments in the standard real model. Its cumulative distribution has
half-integral exponent and is not the polynomial used above.

A single sphere height z has polynomial density proportional to
(1-z^2)^(r-1). Its polynomial CDF V has V(-z)=1-V(z), so this is a
uniform variable. However it is NOT deck invariant on
W x G_m, where tau(z,s)=(-z,-s). The even-power Long phase

    (1-s^-2)((1-V(z))+V(z)*s^2)

is transformed to

    (1-s^-2)(V(z)+(1-V(z))*s^2),

which is a different regular function in general. Thus this exact
even-dimensional descent argument does not give an odd-dimensional
phase. Neither taking a norm nor averaging such phases preserves
purity automatically.

This is a limitation of the PRESENT regular-uniform-variable
construction, not a claim that the Mathieu property for all these
odd-dimensional forms is classified here or unknown in every field.
In particular the existing complex classification has a different
scope. No global open-status or novelty claim is inferred from this
failure of a candidate formula.

There is also an exact determinant-isotropy distinction. Over an
algebraic closure choose a norm-one base point. Its stabilizer in
SO(q) x G_m consists of t=+/-1 and the orthogonal transformation
on the perpendicular complement with determinant t. On tangent
space V the determinant is

    det(t*g)=t^N.

It is trivial for even N, and nontrivial on the second stabilizer
component for odd N. The even case has invariant volume
q^(-N/2) dx_1 wedge...wedge dx_N. In odd dimension the base has no
invariant volume for this action. In ODD dimension, geometrically the double cover
above is its determinant cover; over K it may be a twisted form of
a pointed determinant cover if W(K) is empty, and no pointed
identification over K is asserted.

## 20. Exact source boundary and evidence

The necessary classical algebraic input is normalized reductive trace
existence/base change from finite-dimensional complete reducibility:

- Milne, Algebraic Groups(2017 author PDF),22.42--22.43,
  printed477--478;12l/12.53/Reynolds, printed248--249.
  https://www.jmilne.org/math/Books/iAG2017.pdf
  retained primary SHA256
  503f1a7d511632af05541c52fb2bf6fc36f5a75b0323f2f9110f03dcce11fd67.

The exact known cancellation is credited to:

- Long, arXiv:2607.19012v1,Theorem2.1 and its full beta-binomial proof,
  Proposition3.1.
  https://arxiv.org/html/2607.19012v1
  retained HTML SHA256
  abeec76e95cba0e6244f4af31240634e28d657116c8afcf743b58db144dd7ebd.

The uniform real-sphere column law also appears in Meckes Section1.2,
printed19--20. The exact Gaussian ratio proof needed here is reproduced
in Section15; no asymptotic or numerical spherical law is used.

Sections7--12 above contain the rational trace and formally real Gram
arguments. The even-quadratic proof is independent of the direct SU(n)
moment calculation: its standard real-sphere Gaussian ratio argument
was given explicitly in Section15.


## 21. Independent comparison bridges

### The standard even-quadric pair pulls back to SU(n)

Over Q, regard E^n as a real-coordinate Q-vector space of dimension2n,
with standard q=sum_i(a_i^2+b_i^2) and A=a_1^2+b_1^2.
There is a regular morphism

    Phi:U_n x G_m -> X_(2n),   (M,s)|->s*(first column of M).

Unitarity gives Phi^*q=s^2 and Phi^*A=s^2*u. The real-coordinate
representation rho:U_n->SO(q) makes this map equivariant for the
SO(q) x G_m action through rho x id.

The group map rho need NOT be surjective. The trace comparison uses
the following precise replacement: U_n x G_m is itself geometrically
transitive on X_(2n). Over an algebraic closure E splits into two factors,
and a first column is the pair(v,ell) in the standard SL_n module and
its dual with ell(v)=1. SL_n acts transitively on such pairs: choose
the image vector and a basis of ker ell, then adjust the determinant
in the complementary basis, which is possible for n>=2. Scaling
handles arbitrary nonzero q-values over the algebraic closure.

Consequently the normalized SO(q) x G_m trace on X_(2n), restricted
along rho x id, is its unique normalized U_n x G_m trace. Equivariance
and normalized uniqueness therefore give

    Lambda_(U_n x G_m)(Phi^*c)=Lambda_X(c).

For the actual quadratic-complement pair in Section13,

    Phi^*F=(1-s^-2)((1-v)+v*s^2), Phi^*H=s^-2.

Thus the even-quadric theorem gives a SECOND proof of the SU(n)
degree-two Laurent family over every characteristic-zero field.
It does not use the direct SU(n) Haar calculation in Section3.

Finally the Laurent substitution t|->s^2 preserves constant term
coefficient by coefficient. Hence the moments of this degree-two
family equal those of f(v,t),h=t^-1 in(F1). This proves the original
degree-one family by a given-phase pullback identity. It is not the
false inference that a pure power has a pure root.

Only the FORWARD even-quadric-to-SU(n) comparison is used. No reverse
existence implication from a SU(n) witness to an arbitrary quadratic
complement is asserted.

### Quaternion norm specialization has the same normalized trace

For q=N on the Hamilton quaternion vector space of dimension4,
X={N!=0} is exactly GL_1(H) as an affine variety, and(E3) is exactly(F3).

Its trace defined by SO(N) x G_m and its trace defined by group left
translation agree. Indeed the common group SL_1(H) x G_m maps to
SO(N) x G_m by norm-preserving left quaternion multiplication and
scaling. It also maps onto GL_1(H) through the central isogeny of
Section5. The former trace is invariant under this common group,
hence under GL_1(H) by faithfully flat surjectivity; uniqueness
identifies it with the group trace.

Thus the quadratic-complement theorem supplies an independent route
to the quaternion example, while the direct Section5 proof is retained.
Neither route assumes an all-automorphism or general finite-etale
trace theorem.

## 22. Exact verification and assurance boundary

The single source-guarded replay is

    python3 verify_field_mathieu_product_quadratic.py --output /tmp/field-mathieu-counterexamples.json

The proof establishes its unrestricted orders and fields. The checker
authenticates the precise statement, this proof, the narrow primary
sources and its own source via the ledger. Its finite rational models
support polynomial CDF identities, marked Laurent moments, the
quaternion substitution and rational Gram factorizations; they do not
prove complete reducibility, all-field scalar descent, arbitrary
quadratic-form equivalence or an unrestricted theorem by extrapolation.

Milne's exact retained PDF is REUSED from the existing canonical
source artifact and the narrow source leaf
SRC-RHT-MILNE-COMPLETE-REDUCIBILITY. The package carries a copy only
for portable replay; canonical installation reuses those exact bytes
and that source role. It does not create a second canonical book.

Long's interval method is credited as known. His Theorem4.2 remains
a distinct mathematical import only for the separate split-real-SL2
comparison; it is not a premise of the uniform-variable, product or
even-quadratic moment cores. Meckes is a narrow cross-check of the
column law, which is proved directly. Duistermaat--van der Kallen is
used only for the scalar-torus Mathieu comparison and its finite
coefficient-field extension.

No classification of all real/forms/arbitrary-field homogeneous spaces
is claimed. No odd-dimensional starting phase is produced. Publication,
external review, independent replay, formal verification and worldwide
priority are separate from the written proof and bounded replay.
