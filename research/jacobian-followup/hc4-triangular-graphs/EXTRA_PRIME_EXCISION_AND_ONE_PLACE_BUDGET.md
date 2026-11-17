# Extra source divisors, etale excision and one-place critical budgets

15 November 2026. Claims **HC4MYJA0/1/2**.
The [IZ descended-trace theorem](ONE_PLACE_PENCILS_AND_DESCENDED_TRACE_CLOSURE.md)
is retained. This note constrains the genuinely non-descended cases in every
helper degree. It does not close all rank-one embeddings or HC4/JC2.

## Setup and exact inputs

Let K:A2_z->A2_T be complex plane Keller. Write R=C[z], B=C[T], with
B embedded in R through K. Suppose f is a nonconstant source polynomial
such that

    B -> R/(f) is surjective.                              (1)

For the actual rank-one HC4 embedding

    j(z,w)=(w1,w2,X,Y1,Y2),
    X^2+w1Y1+w2Y2=D(z)+w dot K(z),

IX1 proves (1) with f=X(z,0) and D=f^2, from the closed zero-section
and its exact first free jet. No source-field or differential substitute
for that closed immersion is used.

The following new external inputs are explicit:

- [Stacks 024T, Proposition 41.6.1](https://stacks.math.columbia.edu/tag/024T):
  a section of a separated unramified morphism is open and closed.
- The one-place pencil Milnor identities in
  [Assi-Garcia-Sanchez, arXiv:1407.0490v1](https://arxiv.org/html/1407.0490v1),
  Introduction item (3) and the final genus paragraph of Section 5:
  the total Jacobian/Milnor length of a one-place plane polynomial is
  twice its smooth pencil genus, and for a fibre with normalization
  genus zero the local formula gives
  sum(mu_P)=2 sum(delta_P)-sum(branches_P-1).

These are three individual statements (the section theorem and two
Milnor identities), not conclusions of the symbolic checker.
The already explicit Vistoli one-place theorem and intrinsic pencil
genus formula of IZ0 are reused. All prior source bytes remain pinned.
The section theorem can also be seen from the diagonal: its open and
closed diagonal pulls back to the graph of the section.

## JA0: the complete image equation and disjoint source factor

There are target and source polynomials a,g, with a nonconstant, such that

    ker(B -> R/(f))=(a),
    R/(f) = B/(a),
    a(K)=f g,                 (f,g)=R.                    (2)

Thus the whole source divisor f=0 is an open and closed piece of the
scheme-theoretic inverse image of a=0. Every remaining source prime in
g=0 is disjoint from the *entire* f=0 divisor. This holds with all
multiplicities in f, not just its radical.

**The target ideal is principal.** Put I=ker(B->R/(f)). The closed image
has dimension one. The ring R/(f), and hence B/I, is a one-dimensional
Cohen-Macaulay hypersurface ring: locally f is a nonzero divisor in the
regular source plane. In particular I has no embedded closed-point
associated primes. In a primary decomposition of I all associated primes
are height-one primes (a_i) in the target UFD. Localizing a primary
component at (a_i) gives (a_i^m_i) in that DVR; contracting gives that
same principal power in B, since a primary ideal is recovered by this
localization and contraction. Their intersection is their principal
least common multiple. Thus I=(a).

The source quotient surjection gives a(K)=f g for some polynomial g.
Base change of the etale separated map K over A=Spec B/(a) gives

    W=Spec R/(a(K)) -> A.

The isomorphism Spec R/(f)->A from (1) is a section of W->A.
By 024T its image is open and closed in W. At every source point on
f=0, the two quotient ideals (a(K)) and (f) agree locally.
Cancelling the nonzero f shows g is a unit there. Hence f and g have
no common complex zero, and the Nullstellensatz gives (f,g)=R.

This proves more than gcd(f,g)=1: their zero sets do not intersect.
It also explains exactly why a unit modulo f need not be a scalar
polynomial in R.

**Etale excision.** On V=Spec R[1/g], K is etale and its entire
scheme-theoretic inverse image of a=0 is f=0, isomorphic to a=0.
With U=Spec B[1/a] the cartesian square

    Spec R[1/(f g)]  ->  Spec R[1/g]
             |                    |
             v                    v
        Spec B[1/a]      ->     Spec B

is an elementary distinguished Nisnevich square: the bottom map is open,
the right map is etale, and it is an isomorphism on the reduced
complements (indeed on the full divisor schemes).
No motivic or excision theorem is imported to infer an inverse.

The factor g is constant precisely when f is a nonzero constant multiple
of a(K). In the rank-one HC4 setting this is a descended trace and IZ2
already gives the inverse. A hypothetical remaining nonautomorphic
embedding must therefore have a genuinely nonconstant g and actual
additional prime divisors outside f=0.

### Local excision is not global polynomial invertibility

The non-Keller polynomial pair

    L(s,t)=(s,t+s t^2),        f=t,        a=T2,
    a(L)=t(1+s t)

has a closed restriction t=0 onto the target line, with
g=1+s t and 1=g-s f. Its Jacobian is J=1+2s t.
On the source open where g J!=0 it is etale, its complete inverse image
of the target line is t=0, and its generic field degree is still two.
This open source is not A2. It supplies the same elementary excision
pattern without a global polynomial inverse.

This is not a Keller counterexample and does not weaken the hypotheses
of IZ1. It forbids deriving a global inverse merely from the excision
square or from g being invertible on the radical divisor.

## JA1: an additional component must carry a nonconstant source unit

Let C be an irreducible component of the reduced f=0 curve, with one
place at infinity. By (1), K identifies C with its closed target image
A. Write p for an irreducible equation of C and b for one of A.
Then all additional irreducible components of b(K)=0 are disjoint
from the entire f=0 curve.

To verify the componentwise assertion, apply the section argument to
the isomorphism C->A: C is open and closed in K^(-1)(A), so other
components are disjoint from C. If another such component met a second
component of f=0, its image point on A would already have a preimage
on C. Injectivity of K on the entire f=0 curve would identify those
two source points, contradicting disjointness from C.

Assume K is nonautomorphic. IZ1 makes at least one additional component
necessary. Then **at least one additional component D has two or more
points at infinity in its normalization, and p restricted to D is a
nonconstant unit**. In particular a nonautomorphic realization cannot
have only scalar units on all the additional components.

Here is the all-degree proof.

A unit on an irreducible affine curve with one place at infinity is
constant: on its smooth projective normalization its divisor is
supported at one point and has degree zero, so it has no pole and
is constant.

If two one-place plane curves V(p),V(q) are disjoint, restrict each
equation to the other curve. Both are nonzero constant units, so
q-c is divisible by p and p-lambda is divisible by q.
Comparison of total degrees gives equality, and

    q = c1(p-lambda),        lambda!=0.                   (3)

Suppose every component of b(K)=0 had one place at infinity.
Each extra component is disjoint from C, so (3) gives

    b(K)=P(p),

where P has at least two distinct roots; square-freeness is HW0.

But the target polynomial b has finite critical locus. Vistoli's
one-place theorem makes every b-level integral and reduced. A critical
curve would force a squared factor in one of those levels, impossible.
The chain rule and the invertible DK make the critical locus of b(K)
the quasi-finite inverse image of that finite set, hence finite.
A nonlinear P(p), in contrast, has a whole critical curve at any
root of P'. This contradiction proves that some extra D has at
least two punctures.

On that D the function p has no zero, by disjointness. If it were a
constant lambda, D would be contained in p=lambda. Vistoli makes that
whole level integral and one-place. It would equal D, a contradiction.
Thus p is an actual nonconstant unit on D.

Equivalently, a one-place radical component together with scalar units
on every additional source component is a new inverse certificate.
It allows additional components in its hypotheses, unlike the earlier
complete-preimage criterion. No target descent of the trace is assumed.

## JA2: polynomial radical components have a bounded one-place part

Assume in addition that the affine normalization of C is A1. The
same is then true for A. Under hypothetical nonautomorphy C is singular,
by IX1; put

    delta=delta(C)>0,
    B_C=sum_(P in Sing C)(number_of_branches_P-1),
    mu_C=2 delta-B_C>0.

Let s be the number of irreducible one-place components of the entire
b(K)=0 curve, including C itself. Then

    s mu_C <= 2 delta,       so s<=2.                     (4)

If s=2, every singular point of C is an ordinary node. If C has any
non-nodal singularity, s=1. Thus there is at most one additional
one-place component, and it is possible in the necessary conditions
only for an all-nodal selected curve. JA1 still requires at least one
further component with multiple punctures and a nonconstant unit.

### Why each one-place component is the same abstract curve

Let D be one of the one-place components. Its normalization maps etale
to the normalization A1 of A. Indeed normalization commutes with the
etale base change K^(-1)(A)->A: the base change of the finite normalization
is finite and birational, is normal because it is etale over the normal
curve, and hence is the normalization. Its connected pieces are the
normalizations of the source components.

The nonconstant projective extension has a preimage of the target
puncture. Since D has one puncture, it is the only such point.
The exact boundary argument retained in IY gives

    2g(D_normal)-1 = -e,

where e is the normalized degree over A1. Therefore e=1 and
g(D_normal)=0. The normalized map is an isomorphism of A1.

By (3), D is a level p=lambda. IZ0 makes its generic pencil genus the
same as that of p, namely delta(C). Since its normalization is A1,
the intrinsic formula applied to this level gives

    delta(D)=delta(C).

The map D->A is finite birational. To see finiteness without assuming
a component of an etale curve map is automatically proper, align their
normalizations C[t]. The target ring is contained in O(D) inside C[t],
and C[t] is finite over the target ring. Noetherianity makes O(D) a
finite submodule, hence finite. The equality of the two normalization
deficits then gives O(D)=O(A). Thus D->A is an isomorphism.

### The complete critical budget

Every such D is therefore a level p=lambda whose finite singularities
have the same total Milnor number mu_C as C. The newly explicit
one-place Milnor identity gives the global critical length of p:

    dim_C C[z]/(p_z1,p_z2) = 2 G(p) = 2 delta(C).

Different levels use disjoint parts of this finite critical scheme.
Summing the s identical contributions proves (4).

For completeness delta>=B_C locally and globally: normalization allows
independent constants on r branches, whereas a local function has one
common residue, so its quotient has at least r-1 dimensions.
Hence mu_C>=delta, proving s<=2. Equality s=2 forces delta=B_C.
At every singular point the completed local ring then consists of all
branch series with one common constant. Its embedding dimension is r;
a plane curve has embedding dimension at most two. The singular case
is r=2 with completed ring C[[x,y]]/(xy), an ordinary node.
Thus a cusp, tacnode or higher multiple singularity excludes any
additional one-place component.

### The necessary nodal exception is not a constructed Keller failure

The pencil p=y^2-x^2(x+1) has generic genus one. Its two critical
values 0 and -4/27 both have one node of Milnor number one. Thus
two isomorphic one-place nodal levels really can exhaust the budget.
The identity

    p(-x-2/3,i y)=-p(x,y)-4/27

relates those two levels. This is not a Keller map identifying both
levels with one target curve, and no nonautomorphic example is asserted.

### A forced finite-image boundary on an extra component

Choose the multiple-puncture D supplied by JA1. For its normalization
over the target normalization A1, let e be the degree, g_D its
completion genus, r_D>=2 its puncture count. The same exact boundary
formula gives

    E_D=e+2g_D-2+r_D >= e >=1.                            (5)

Here E_D sums indices at source punctures whose image is a *finite*
point of the target normalization. Thus actual finite target values
are approached from infinity on this extra source component. They
belong to the nonproperness set of K. In this scope the additional
divisor is not merely an abstract source prime: it carries a
nonconstant unit and a compulsory finite-image boundary.

## Assurance and remaining work

See [source context](sources/extra-prime-excision-context.json) for the
new section and Milnor inputs and retained normalization, one-place
and Keller arguments. In JA2 the standard etale-preserves-normality
fact is the same retained fact used in HY0. No new theorem about
Nisnevich cohomology, global inverse reconstruction from an excision
square, or finite affine curve maps is imported.

Run

    .venv/bin/python research/jacobian-followup/hc4-triangular-graphs/verify_extra_prime_excision.py

for [exact controls](receipt-extra-prime-excision.json); see
[scoped validation](validation-extra-prime-excision.json).
The disjointness, unit obstruction, nonfinite open-source example,
cusp/nodal critical budgets and genuine nonlinear embedding controls
are separated from the written all-degree proof and its external inputs.

What remains is to use the full w-dependent closed embedding to exclude
these actual extra primes and their boundary behaviour, or force trace
descent or image saturation. The zero-section and local excision alone
do not do that. Genuinely non-descended traces without the one-place
component hypothesis, higher helper ranks, unaligned intrinsic ambient
projections, actual endpoints and unrestricted HC4/JC2 remain open.
