# Sharp quadratic pure-moment threshold in every covariance rank

14 November 2026. This completes the local written audit of Theorem 5.1 in
Zongmin Liu's [Sharp finite Gaussian moments for quadratic polynomials, v3](https://zenodo.org/records/21688251).
It is a credited reproduction, not a new priority claim or an external review.
The downloaded source identity is retained in the [quadratic audit](QUADRATIC_GAUSSIAN.md).
The proof below includes degree drop, repeated poles and the polynomial
transfer part. The exact checker realizes the small-rank algebra independently;
it does not enumerate the projective intersection in arbitrary rank.

For every N>=1, the least number of initial pure Gaussian moments forcing
all positive pure moments to vanish for every complex quadratic in covariance
rank N is exactly 2N+1. The [previous finite theorem](QUADRATIC_GAUSSIAN.md)
proves sufficiency. Here is a complete existence argument proving necessity.

## A rational symmetric realization with no genericity assumptions

Let D(0)=1, deg D<=N and deg B<N. There is a complex symmetric N-by-N
matrix A and a vector l such that

    det(I-2tA)=D,
    l^T*(I-2tA)^(-1)*l=B/D.                         (1)

To verify this, split the reduced rational function B/D into its polynomial
part q and strictly proper part b0/D0, with D0(0)=1 and b0,D0 coprime.
For the strictly proper part let s=deg D0, write
D0=1+d1*t+...+ds*t^s and b0/D0=sum h_j*t^j. On

    W=C[x]/(x^s+d1*x^(s-1)+...+ds)

define phi(x^j)=h_j for 0<=j<s. The denominator recurrence makes this
identity valid for every j>=0. The Hankel form <f,g>=phi(fg) is
nondegenerate: a nonzero vector in its radical would give an annihilating
recurrence of degree less than s for the whole sequence h_j, contradicting
the reduced denominator. Multiplication T by x is self-adjoint for this
form. With v=1,

    <v,(I-tT)^(-1)v>=b0/D0,       det(I-tT)=D0.

Any nondegenerate complex symmetric Gram matrix H has H=S^T*S with S
invertible. Consequently STS^(-1) is symmetric in the ordinary bilinear
form, and Sv gives the required vector for this block. This argument also
handles repeated roots of D0.

If q has degree p>=0, use C[x]/(x^(p+1)), multiplication by x and
phi(x^j)=[t^j]q. Its Gram matrix is anti-triangular with nonzero
anti-diagonal [t^p]q, so it gives a symmetric nilpotent realization of q
in p+1 dimensions, with determinant one. This block is essential when D
has degree drop; it must not be discarded as a supposedly proper transfer.

Realize the unused determinant factor D/D0 by a diagonal block, with zero
vector on it. If q=0 the total dimension is deg D. Otherwise it is
p+1+deg D=deg B+1. Either is at most N. Pad with zero blocks and set
A=T/2. This proves (1), including zero proper or polynomial parts.

## Nonzero logarithmic-gap data exists in every rank

For d=(d1,...,dN), let

    D_d=1+d1*t+...+dN*t^N,
    g_k(d)=[t^k](D_d*log D_d).

Give d_j weight j. Each g_k is weighted homogeneous of degree k; it is
nonzero since, on d2=...=dN=0,

    g_k=(-1)^k*d1^k/(k*(k-1))  (k>=2).

Consider the N-1 equations g_(N+2)=...=g_(2N)=0. Substitute d_j=x_j^j.
They become N-1 ordinary positive-degree homogeneous equations in N
variables. Krull's height theorem gives an affine zero set of dimension
at least one, so over C it contains a nonzero point, not merely the origin.
For N=1 there are no equations and any nonzero x1 works. At such a point
some d_j is nonzero, hence D_d!=1.

This argument needs neither squarefree D nor full degree N. It also makes
no claim that this intersection is reduced, finite or has a particular
number of geometric points. The source's Catalan intersection product is
not needed for the threshold theorem and is not independently certified here.

## Turning the gaps into a sharp quadratic

Choose the nonzero gap point above, write D=D_d and put

    C(t)=D(t)*(log D(t)-d1*t),
    t^2*B(t)=truncation of C to degrees 2,...,N+1.

Then deg B<N. Because d1*t*D has degree at most N+1, the gap equations
imply

    d1*t/2-log(D)/2+t^2*B/(2D)=O(t^(2N+1)).          (2)

Use (1) and let P=X^T*A*X+l^T*X-tr(A). Since d1=-2tr(A), the left
side of (2) is the formal logarithm of E(exp(tP)). Thus the first 2N
cumulants, and hence moments, are zero.

The next coefficient cannot be zero: if the first 2N+1 moments vanished,
the finite quadratic theorem would force A nilpotent and D=1, contradicting
the chosen nonzero gap point. Therefore the (2N+1)-st moment is nonzero.
This proves that the sufficient threshold is optimal in every N.
All data can be chosen over the algebraic numbers by doing the zero-set
argument over their algebraic closure and using algebraic matrix congruence.
No coefficient-height bound or explicit all-rank list is supplied.

## Exact small-rank realization checks

Rank one uses P=X^2-1+sqrt(-2)*X. Its first moments are 0,0,-4,-48,
confirming that two are insufficient and three suffice.

For rank two choose delta with 6delta^2-6delta+1=0 and set

    D=1+t+delta*t^2,       B=delta-1/2-t/6,
    h0=delta-1/2, h1=1/3-delta, h2=delta/2-1/6,
    H=((h0,h1),(h1,h2)),   T=((0,-delta),(1,-1)).

Then HT=T^T H, det H=1/18-delta/4!=0, det(I-tT)=D and
v^T H(I-tT)^(-1)v=B/D for v=(1,0). An explicit ordinary symmetric
realization has

    alpha=h1/h0,   beta^2=det(H)/h0^2,   gamma^2=h0,
    A=((alpha,beta),(beta,-1-alpha))/2,   l=(gamma,0).

Consequently

    P=(alpha*X1^2+2beta*X1*X2-(1+alpha)*X2^2)/2
         +gamma*X1+1/2.

Direct Gaussian Wick expansion in the exact field Q(sqrt(3)), with
beta and gamma eliminated by their squared values, gives the first four
moments zero and the fifth 2*(5delta-1)!=0. Either root delta and either
choice of each square root works. The independent
[checker](verify_quadratic_sharpness.py) also checks the determinant,
self-adjointness and transfer identities, so these are actual Gaussian
polynomial witnesses, not merely a formal cumulant sequence.
