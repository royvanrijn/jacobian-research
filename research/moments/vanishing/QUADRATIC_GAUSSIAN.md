# Quadratic Gaussian vanishing: finite gate and explicit consequences

Reviewed 14 November 2026. Attribution: Zongmin Liu,
[Sharp finite Gaussian moments for quadratic polynomials, v3](https://zenodo.org/records/21688251).
The downloaded source archive has SHA-256
`196fed360a70d9150709bf981af4304644021294783a525b2eafef246d23be2f`.
Its Theorems 2.1 and 2.2 supply finite determination, classification and the
optimal mixed bound. The separate [all-rank sharpness audit](QUADRATIC_SHARP_THRESHOLD.md)
now supplies the complete written realization/intersection argument for
Theorem 5.1, with independent exact Gaussian witnesses in ranks one and two.
Its general existence proof is not an enumerated all-rank certificate.
No priority claim is made for the following deductions or explicit examples.

## Audited finite gate

For standard Gaussian X in N>=1 real coordinates, write a centered complex
quadratic as P=X^T A X+l^T X-tr(A), with A symmetric. Vanishing of the first
2N+1 moments implies

    A^N=0,     l^T A^j l=0 (0<=j<N),

and these identities imply every positive pure moment vanishes and
E(QP^m)=0 for m>floor(N/2) deg(Q). A polynomial whose first moment vanishes
is already centered. Arbitrary real Gaussian laws reduce to their covariance
rank by restriction to their affine support. Rank zero is deterministic and
must be treated separately: its first moment forces P=0 on that support.

For clarity, here is the short algebraic proof used in this audit. Put
D=det(I-2tA), B=l^T adj(I-2tA)l. The formal logarithmic generating function is

    L=-t tr(A)-(1/2)log D+(t^2/2)B/D.

Thus L'=R/D^2 with deg R<=2N. The first 2N+1 zero moments give
L'=O(t^(2N+1)), hence R=0. Exponentiating 2L=0 gives exp(S)=D for a
rational S with S(0)=0. The derivative of a rational function has zero
residue at every finite pole, whereas D'/D has a nonzero positive integral
residue at a root of D. Therefore D=1. Cayley--Hamilton and the displayed
formula for L give precisely the finite matrix identities above. Conversely
those identities give L=0. This is a formal power-series argument; analytic
convergence of a complex quadratic moment generating function is unnecessary.

## A smaller per-input bound

Suppose the matrix conditions hold. Let nu be the nilpotency index of A
(take nu=1 when A=0) and s=dim span{l,Al,...} (take s=0 when l=0).
For d=deg Q the following sufficient cutoff refines the uniform one:

    E(QP^m)=0 for m>B_A,l(d),
    B_A,l(d)=s*d+floor(d/2)*max(0,nu-1-2*s).

Indeed Gaussian tilting gives polynomial mean
mu=t(I-2tA)^(-1)l of degree at most s and covariance
Sigma=(I-2tA)^(-1) of degree at most nu-1. The cyclic Krylov space is totally
isotropic since (A^i l)^T A^j l=l^T A^(i+j)l=0. Hence
s<=floor(N/2) and A^s l=0. A Wick partition of a degree-d monomial with
p pairs and d-2p singletons has t-degree at most
(d-2p)s+p(nu-1). Maximizing over 0<=p<=floor(d/2) proves the formula;
it is nondecreasing in d and therefore covers lower-degree monomials too.
Since nu<=N it is at most floor(N/2)d. Constant Q has bound zero.
This is a sufficient per-input bound, with no universal sharpness claim.

## An explicit chain attaining the uniform bound

Let h=floor(N/2)>=1 and take circular pairs (Z_i,W_i) with
E(Z_i^a W_i^b)=delta_ab a!, independently for 1<=i<=h. Ignore an extra
real coordinate when N is odd. Set

    P_h=Z_1 W_2+Z_2 W_3+...+Z_(h-1) W_h+Z_h.

Every summand has weight one for wt(Z_i)=h-i+1 and
wt(W_i)=-(h-i+1). Thus all positive pure moments vanish. In
E(W_1^d P_h^m), balancing each circular pair forces exactly d copies of
each of the h summands, so

    E(W_1^d P_h^m) = (h*d)! if m=h*d, and 0 otherwise (m>=1).

The multinomial coefficient (hd)!/(d!)^h cancels the h Wick factors d!.
This proves attainment for every d>=1, including linear Q, directly in
polynomial coordinates. The N=1 case has only P=0: a scalar nilpotent A
is zero and l^2=0 forces l=0. The checker verifies the chain for h=1..4,
d=0..3 and m<=3h+1; the balancing argument establishes all h,d,m.

## Complete three-variable quadratic normal form

For N=3 every pure-moment-zero quadratic, after a complex orthogonal
change of polynomial coordinates and a circular normalization, has

    P=Z*(a*Z+b*T+c),   E(ZW)=E(T^2)=1,
    E(Z^2)=E(W^2)=E(ZT)=E(WT)=0.

Conversely every such polynomial has all positive pure moments zero and
E(QP^m)=0 for m>deg Q. These are identities of the Gaussian Wick
functional; a complex orthogonal change need not be a real change of the
probability space.

Proof: the isotropic Krylov space has dimension at most one, so Al=0 and
l^T l=0. If A=0, choose its isotropic direction from l (or choose any if
l=0). If A has rank one, write A=alpha*v*v^T; nilpotence makes v isotropic.
The conditions v^T l=0 and l^T l=0 force l proportional to v, because
v-perp modulo v is a one-dimensional nondegenerate space. If A has rank
two, its Jordan type is (3), and ker A is an isotropic line v: it lies in
im A=(ker A)-perp. Again l lies in that line. Choose a Witt frame based
on v. A v=0 makes the quadratic independent of W; trace A=0 removes its
T^2 term. The linear part is proportional to Z, giving the displayed form.
Conversely every monomial has positive Z-minus-W weight. The weight of
any monomial in Q is at least -deg Q, proving the mixed bound. The zero
polynomial is included. Taking P=Z and Q=W^d proves the coefficient one
cannot be decreased even within this three-variable layer.

This closes the entire quadratic layer, including arbitrarily dense
quadratics in the original real coordinates. It leaves open three-variable
cubic GMC failures: the known degree-four witness and sparse Pareto theorem
do not exclude all dense cubics. It also leaves the ordinary-Laplacian
lifting and lower-power homogeneous GVC questions untouched.
