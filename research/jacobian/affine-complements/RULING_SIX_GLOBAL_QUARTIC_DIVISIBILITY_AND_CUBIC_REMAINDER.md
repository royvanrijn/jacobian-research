# Global quartic divisibility and a regular cubic remainder

PF2RULED169. 17 November 2026. Actual complex plane polynomial Keller maps,
within the complete marked quadratic source model.

Every remaining ruling-six coordinate has the stronger regular square
reduction

    P=W3^2+R3, deg_w W3=3, deg_w R3<=3,

where both polynomials belong to the actual source-image algebra.
The coefficient of the former quartic remainder is divisible by h*M.
This is a necessary condition on the entire remaining A class; W3 is
not asserted to be a target coordinate or integral over the target.

## 1. Exact setup and retained dependencies

Use PF165 notation

    w^2=N+M*z, N=u*Sigma^2, M=u*A,
    u=t^2, w=t*(Sigma+A*x), z=2Sigma*x+A*x^2,
    Omega=du wedge dw/M,
    P=M*h^2*z^3+2hK*w*z^2+b2*z^2+w*(a1*z+a0)+b1*z+b0.

Every displayed coefficient is in C[u], h!=0, and K=A2/2. Put

    T=(h*w+K)*z, r=b2-N*h^2-K^2.

[PF165](RULING_SIX_SIMPLE_ESCAPE_GRAPHS_AND_QUARTIC_FACE_REDUCTION.md)
gives the exact identity

    P=T^2+r*z^2+w*(a1*z+a0)+b1*z+b0,                         (1)

and a globally reduced actual companion of degree s=3q, q>=3 odd.
The leading coefficient relation is

    p6=(h/M)^2, q_s=beta*(h/M)^q, beta!=0.                    (2)

All original full-normalization primes remain retained.

We use two further explicitly scoped results.
[PF167](RULING_SIX_APPROXIMATE_CUBIC_ROOT_REGULARITY_AT_M_DIVISORS.md) proves

    ord_c r>=ord_c h+ord_c M at every root c of M.             (3)

[PF168](RULING_SIX_OUTER_SQUARES_AND_EVEN_SPECIALIZATIONS.md) proves that
the specialization P*(c,w), c off M, cannot have
degree three or five, and that every outer degree-six weighted initial
is a square. The simple-escape argument below uses the same cached
Chau component theorem as PF165: a closed affine-line nonproper-value
component is impossible for an actual complex plane Keller map.
These are written-proof dependencies, not conclusions from sampling.

## 2. The global divisibility theorem

We prove

    h*M divides r in C[u].                                  (4)

If r=0 this is immediate. Otherwise it suffices, in view of (3), to
prove ord_c r>=ord_c h at every root c of h off M.

Fix such c and write v=u-c, j=ord_c h>=1. Suppose, for a
contradiction, delta=ord_c r<j. All coefficients of P*, Q* are
regular at v=0. If K is nonzero put k=ord_c K; allow k=infinity
when K=0. With P*=sum p_i*w^i, its top coefficients are

    p6=h^2/M^2,
    p5=2hK/M^2,
    p4=(r+K^2-2N*h^2)/M^2,
    p3=a1/M-4N*hK/M^2.                                     (5)

The other coefficients are regular. On the whole original source
line over c, P and every reduced companion have degree at least two.

Choose the outermost positive Newton slope ell at which the degree-six
term meets a lower term. It satisfies ell>j/2. To check this even
when the leading p4 coefficient cancels, if delta<2k then p4 has
order delta<j, forcing ell>=j-delta/2>j/2. If delta>=2k, then
k<j/2 and p5 forces ell>=j-k>j/2. The case K=0 belongs to
the first alternative.

Let I(Y) be the entire initial after w=v^-ell*Y. Its degree is six.
By PF168 it is a square. Its lowest nonzero Y exponent is even.
The regular constant/linear coefficients cannot appear, since
H=6ell-2j>2ell>ell. The lowest exponent is therefore two or four.

If it were two, regularity of p2 gives H<=2ell, so ell<=j/2,
a contradiction. Consequently it is four, and

    I(Y)=C*Y^4*(aY+b)^2, a*b!=0.

Comparing p6,p5,p4 in (5) yields

    k<infinity, ell=j-k,
    ord p4=2k, initial(p4)=initial(K^2/M^2).

The nonzero leading coefficient of r cannot cancel in this last
equality if delta<=2k: for delta<2k it would dominate, and for
delta=2k it would alter that coefficient. Thus necessarily

    2k<delta<j, ell=j-k>k.                                  (6)

This argument includes the delta=2k cancellation case.

## 3. A lower escape controls the cubic coefficient

If k=0, regularity already gives ord p3>=k. Suppose k>0.
Then (6) makes h,K,r vanish at c, so deg P*(c,w)<=3.
It is at least two, and PF168 excludes degree three. Hence

    p2(c)!=0.

Let a=ord p3. If a<k, then a>=1 because the specialization has
degree two. Substitute w=v^-a*Y and multiply P*-p by v^(2a).
The full initial is

    Y^2*(gamma+eta*Y), gamma*eta!=0.                         (7)

Indeed the degree-four term has positive margin 2k-2a;
degrees five and six have positive margins j+k-3a and
2j-4a. Positivity follows from a<k and j>2k.
Degrees zero and one have positive margins as well.

The nonzero root of (7) is simple and constant, and the exact coframe
has order a>=1. All-order Hensel recursion is over C[p], with only
constant unit divisions. The finite companion residue is therefore
a polynomial in p. The original lifts have a source-coordinate pole,
so the full-normalization image is the entire closed polynomial graph.
Chau's component theorem excludes it, exactly as in PF165 Section2.
It follows that

    ord p3>=k, and hence ord a1>=k.                          (8)

The latter follows from (5), since ord(hK)=j+k>k.

## 4. The outer double cluster splits into a forbidden graph

Use the EXACT zero

    wc=-K/h

of T=(h*w+K)(w^2-N)/M. Its pole order is ell=j-k.
Set

    beta0=j-delta/2, w=wc+v^-beta0*Y.

By (6), 0<beta0<ell, so (w-wc)/wc tends to zero. Pass to
v=s^2 if necessary; this is a fixed, p-independent Puiseux base
change. The degree-six polynomial (1), including all its lower
terms, after multiplication by v^(4ell-delta), has reduction

    C1*Y^2+C0, C1*C0!=0.                                   (9)

For clarity, the orders before scaling are:

    T^2:       2k-2ell-2beta0=delta-4ell,
    r*z^2:     delta-4ell,
    a1*w*z:    >=k-3ell, gap >=j-delta>0,
    b1*z:      >=-2ell,  gap >=2ell-delta>0,
    a0*w,b0:   strictly higher.

The two roots of (9) are nonzero, simple and constant. Differentiating
the entire polynomial at either lift gives the exact order

    ord_v(w_p)=3ell-k-delta/2>0.                            (10)

After the Puiseux base change, omega=-w_p*du/M is holomorphic.
The lifted Laurent coefficients of w and Q belong to C[p]((s)).
Holomorphic dQ excludes every negative Laurent term of Q; its
residue is polynomial in p. Since wc has a pole and M(c)!=0,
both unit t lifts have x polar. These are original omitted primes,
with full closed polynomial-graph image, again forbidden by Chau.

This contradiction excludes delta<j and proves (4).
Neither algebraic residues of degree two nor omitted high-ramification
primes have been suppressed: (7) and (9) construct additional simple
constant-root branches whenever the proposed divisibility fails.

## 5. A regular cubic square with remainder of ruling at most three

Write

    B0=r/(2hM) in C[u], W3=T+B0*w.

Both are actual regular source-image polynomials. Using w^2=N+M*z
in (1) gives the exact identity

    P=W3^2+R3,
    R3=(a1-2B0*K)*w*z
       +(b1-2B0*h*N-B0^2*M)*z+a0*w+b0-B0^2*N.              (11)

Thus the previous quartic remainder is reduced to ruling at most
three, without changing the actual target field or normalization.
On t=0, W3=T=2*Sigma(0)*K(0)*x, a nonconstant affine function.
Consequently W3 is outside F=C(P,Q), by the same source-intersection
and one-line argument as PF165. R3 is also outside F: otherwise
W3^2=P-R3 is in B, and the squarefree pullback argument of PF165
would put W3 in B, a contradiction.

To lower the remainder to ruling two by adding a polynomial in u
to W3 would require the NEW divisibility

    h divides a1-2B0*K.                                     (12)

It is not proved here; the subsequent
[quartic trace argument](RULING_SIX_QUARTIC_TRACE_AND_CONSTANT_SQUARE_ROOT_CORRECTION.md)
and [nested two-cluster theorem](RULING_SIX_REGULAR_CUBIC_APPROXIMATE_ROOT_AND_QUADRATIC_REMAINDER.md)
prove it globally and reduce the regular remainder to ruling two.
A rational polynomial-part square root at
w=infinity does not by itself establish (12), regularity, target
membership or target integrality. In particular (11) closes neither
the remaining A faces nor derivative/trace integrality nor JC2.

## 6. Verification scope

The checker tests the entire lower and translated outer substitutions,
their simple roots and exact coframe orders, including odd delta
requiring Puiseux parameters. It also checks the outer square
coefficient obstruction, the algebraic reduction (11), and numerical
controls at all relevant inequality boundaries. The written order
comparisons and all-order constant-unit recursion establish the
unrestricted statement. No coefficient-degree cutoff is part of it.
The checker passes 116 exact controls.

Replay:

    .venv/bin/python research/jacobian/affine-complements/verify_ruling_six_global_quartic_divisibility.py

All earlier source-model and external component inputs retain their
recorded verification scope. Independent whole-proof replay, formal
verification and external human review remain absent.
