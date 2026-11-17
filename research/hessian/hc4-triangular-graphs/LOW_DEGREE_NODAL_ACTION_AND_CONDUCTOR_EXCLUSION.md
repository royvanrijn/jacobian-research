# Low-degree nodal actions and a degree-two conductor exclusion

Let

    A_0={f in C[t]:f(1)=f(-1)}
       =C[t^2-1,t*(t^2-1)].

If p,q generate A_0 and at least one of p,q has degree at most THREE,
then

    integral_(-1)^1 p(t)*q'(t) dt != 0.                       (1)

There is NO degree bound on the other generator. The proof uses the
exact finite-order index, including a cubic norm identity. It does not
assume that all embeddings of the nodal curve are ambient-equivalent.

Combined with the actual Keller conductor theorem, (1) excludes the
entire global conductor-degree-two case whenever either output has
directional source degree at most three, even after a polynomial target
automorphism. The case in which both generators resist every such
reduction remains open.

## 1. The one-dimensional index of the nodal algebra

Put R_0=C[t]. The linear functional f->f(1)-f(-1) has kernel A_0,
so R_0/A_0 has complex dimension one. Let p be a nonconstant element
of A_0, of degree m. Then D=C[p] is a polynomial ring and R_0 is finite
free over D with basis 1,t,...,t^(m-1), after dividing the monic source
equation by the nonzero leading coefficient of p.

If A_0=C[p,q], then q is primitive over Frac(D), since Frac(A_0)=C(t).
It is integral over D. Its minimal polynomial is monic of degree m,
with coefficients in the normal D. Therefore 1,q,...,q^(m-1) is a
D-basis of A_0.

Let Delta(U) be the determinant of this basis inclusion into R_0,
where U represents p. Smith normal form over C[U] gives

    deg_U Delta=dim_C(R_0/A_0)=1.                            (2)

More precisely Delta is a nonzero constant times U-p(1), because the
one-dimensional quotient is annihilated by p-p(1).
This is an exact index statement, not a bounded coefficient search.

Subtracting p(1) from p changes neither its algebra nor its action
integral: q(1)=q(-1). We henceforth take p(1)=p(-1)=0.

## 2. Every quadratic generator: arbitrary other degree

A degree-two p in A_0 has the form p=a*(t^2-1), a!=0.
Write the other generator uniquely over D as

    q=H(p)+t*B(p), H,B in C[U].

The inclusion determinant is B(U). Equation(2) and B(0)=0 give
B(U)=beta_0*U with beta_0!=0. Thus

    q=H(p)+beta*t*(t^2-1), beta=beta_0*a!=0.

The H(p) term contributes zero action, because p*dH(p) is a polynomial
differential in p and p has the same endpoint values. Direct integration
gives the all-degree identity

    integral p*dq=8*a*beta/15 != 0.                         (3)

Every nonconstant element of A_0 has degree at least two. Thus degree
zero or one cannot supply an additional nonconstant-generator case.

## 3. Every cubic generator: the exact norm forces low remainder degree

After subtracting its node value, a degree-three p has the form

    p=(t^2-1)*(a*t+b), a!=0.

Rescale p by a nonzero constant, and first treat a=1. Its source
equation over D is

    t^3+b*t^2-t-b-U=0.

Write q uniquely as

    q=H(p)+B(p)*t+C(p)*t^2, H,B,C in C[U].                  (4)

Let t_1,t_2,t_3 be the three roots over a splitting field of C(U).
For i!=j, the quotient of the differences of q and t is
B+C*(t_i+t_j). Since the sum of the three roots is -b, the product
over the three pairs gives

    Delta(U)=Norm_(C(t)/C(p))(B(p)-b*C(p)-t*C(p))
            =B^3-2*b*B^2*C+(b^2-1)*B*C^2-U*C^3.          (5)

The additive H(p) does not affect this determinant. The same formula
follows from the companion multiplication matrix, independently of
the splitting-field transcription.

For EVERY nonzero polynomial r(t), its norm to C(p) has U-degree
deg_t r. Indeed that norm is the resultant of the monic source
polynomial p(T)-U with r(T); reversing the resultant expresses it
as a product of deg_t r nonconstant linear factors in U, counted
with multiplicity.

By(2), the polynomial

    r(t)=(B-b*C)(p(t))-t*C(p(t))

therefore has degree ONE. If C=0 its degree is a multiple of three,
contradicting (2). If C!=0, the two displayed terms have degrees in
different residue classes modulo three. Their leading terms cannot
cancel. Hence C is a nonzero constant and B-b*C is constant, so B
is constant too. The endpoint identity q(1)=q(-1) in(4) now gives B=0.

Consequently EVERY cubic-generator pair has the form

    q=H(p)+C*t^2, C!=0.                                    (6)

Restore the original leading coefficient a. The exact action is

    integral p*dq=-8*a*C/15 != 0.                          (7)

The coefficient b and the arbitrary polynomial H do not affect it.
This proves(1) with no bound on deg H or on deg q.
If it is q that has small degree, interchange the generators;
integration by parts changes the action's sign.

## 4. Action covariance under polynomial target automorphisms

For a polynomial automorphism Phi=(F,G) of the ambient plane, let
lambda=Jac(Phi) in C*. The polynomial one-form

    F*dG-lambda*U*dV

is closed, hence polynomial-exact on A2. The endpoint values of
(p,q) coincide, so its integral along the nodal normalization path
is zero. Thus

    action(Phi(p,q))=lambda*action(p,q).                     (8)

Accordingly action nonvanishing holds for every generator pair that
can be reduced to the preceding low-degree class by an actual
polynomial ambient automorphism. No such reduction is asserted for
all generator pairs.

## 5. Actual Keller exclusion in the degree-two conductor class

Use the setup and exact lowest-degree classification in Section6 of
[the coordinate-order conductor theorem](GENERIC_COORDINATE_ORDER_CONDUCTOR_AND_CLOSED_COLLISIONS.md).
Suppose an ACTUAL noninvertible plane Keller map has global degree two
in x for its coordinate-order conductor. After its polynomial source
shear and constant scale, write the source as C[u,y]. Then

    C[P,Q,y]=A_0 tensor_C C[y]
            =C[y,u^2-1,u*(u^2-1)],

where A_0 is identified with the fixed nodal algebra in u.
The normalized source remains the whole affine plane. Its Jacobian
J=Jac_(u,y)(P,Q) is a nonzero CONSTANT.

On EVERY fibre y=c, the polynomials p_c=P(u,c),q_c=Q(u,c) generate
the exact nodal algebra A_0 in C[u]. In particular

    P(1,y)=P(-1,y), Q(1,y)=Q(-1,y).

Define the polynomial in y

    I(y)=integral_(-1)^1 P(u,y)*Q_u(u,y) du.

Differentiation and polynomial integration by parts give

    I'(y)=[P*Q_y]_(u=-1)^(u=1)
           +integral_(−1)^1 (P_y*Q_u-P_u*Q_y) du
         =-2*J.                                            (9)

The endpoint term vanishes by the displayed polynomial identities.
Thus I is affine with nonzero slope and has a complex zero.

If either output has u-degree at most THREE, its restriction on each
fibre is nonconstant of degree two or three. A constant restriction
could not generate the nonnormal A_0 together with just one other
element, and no degree-one polynomial belongs to A_0. Sections2--3
therefore make I(c) nonzero for EVERY c, contradicting (9).

This proves the actual exclusion

    global conductor x-degree two
    AND min(deg_x P,deg_x Q)<=3
       ==> no noninvertible Keller map.                   (10)

The source shear used in the conductor theorem preserves these
directional degrees. By(8), the same exclusion applies if ANY
polynomial target automorphism achieves that degree bound. Target
automorphisms preserve B and hence A=B[y] and its intrinsic conductor.

No assertion about unrestricted zero-action nodal embeddings is used.
No trace-integrality or JC2 solution follows for the remaining class.

## 6. A zero-action immersion that fails the generator hypothesis

The polynomial functions

    p=t^2-1, q=t*(t^2-1)*(t^2-3/7)

have zero action between -1 and 1 and give an immersion everywhere.
They also identify the additional pair t=+sqrt(3/7),-sqrt(3/7).
Their quadratic inclusion determinant is U*(U+4/7), of degree TWO,
so C[p,q] is a proper subalgebra of A_0 and R_0/C[p,q] has dimension two.
This exact control shows why the full generator equality in(1) matters.
A zero action for arbitrary nodal functions or for an immersion with
a chosen node would not disprove the theorem.

## 7. Verification boundary

The proof uses finite free monogenic modules, the determinant index
over C[U], a written norm-degree argument, and polynomial integration.
The actual Keller application retains the proved coordinate-order
classification. It imports no classification of all nodal embeddings.

Exact controls check the full cubic norm determinant, quadratic and
cubic actions with higher-degree H, action covariance under two
nonlinear target automorphisms, the all-parameter variation identity,
and rejection of the false claim that every pair of nodal functions
with zero action generates A_0. The controls support these written
unbounded arguments; they do not replace them.
No independent mathematical review, formal verification or external
review is claimed. The general nodal-action question and full
trace-integrality/HC4/JC2 goal remain open.
