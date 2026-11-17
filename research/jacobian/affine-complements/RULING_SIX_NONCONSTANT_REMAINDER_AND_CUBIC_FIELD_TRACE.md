# The quadratic remainder cannot depend only on u

PF2RULED172. 17 November 2026. Actual complex polynomial plane Keller maps.

In PF171's regular normal form P=H3^2+R2, the remainder R2
cannot be a polynomial in u alone. Its w degree is therefore
one or two. A field trace through the cubic H3 gives the
obstruction uniformly for every odd companion degree.

This is an exclusion of the entire ruling-zero remainder class,
with arbitrary coefficient degrees. The remaining linear and
quadratic remainders, full ruling-six closure and JC2 remain open.

## 1. Actual setup and exact trace domain

Use the original complete quadratic source-image algebra and
actual fields of PF165--171:

    w^2=N+M*z, u=t^2, E=C(u,w), F=C(P,Q), L=C(x,t),
    dP wedge dQ=du wedge dw/M,
    P=H^2+R2,
    H=H3 in O(V), deg_w H=3,
    R2=d(u)*z+e(u)*w+f(u), d,e,f in C[u].                (1)

The polynomial coefficients in u and the actual companion degree
are unrestricted. The full normalizations and every original
prime remain those of the actual map.

[PF171](RULING_SIX_REGULAR_CUBIC_APPROXIMATE_ROOT_AND_QUADRATIC_REMAINDER.md)
establishes the regularity of H and R2. In C(u)[w] write

    H=a3(u)w^3+a2(u)w^2+a1(u)w+a0(u),
    a3=h/M!=0.                                         (2)

The actual globally reduced companion supplied by PF165 has

    deg_w Q=3k, k>=3 ODD,
    leading_w Q=beta*a3^k, beta in C*.                   (3)

Target changes used to obtain (3) preserve B, F and all original
finite normalizations.

We use the SEPARATE field

    K0=C(u,H), [E:K0]=3.                                (4)

This is the cubic function-field extension defined by (2).
It is separable in characteristic zero. The trace below is
Tr_(E/K0), not a claim about trace landing in the actual target
ring B, the full codifferent or an auxiliary normal-tail algebra.

## 2. Polynomial cubic traces retain the odd leading degree

Let s_m=Tr_(E/K0)(w^m). Dividing (2) by a3 gives the exact
monic relation

    w^3+b*w^2+c*w+d0-H/a3=0,
    b=a2/a3, c=a1/a3, d0=a0/a3.                          (5)

Newton sums, or the actual rank-three multiplication algebra, give

    s0=3, s1=-b, s2=b^2-2c,
    s_m=-b*s_(m-1)-c*s_(m-2)
                           +(H/a3-d0)*s_(m-3), m>=3.    (6)

All s_m lie in C(u)[H]. Induction in (6) proves

    deg_H s_m<=floor(m/3),
    leading_H s_(3j)=3*a3^(-j), j>=0.                   (7)

For m=3j, only the H*s_(m-3)/a3 term can contribute degree j;
the other two terms have degree at most j-1. This proves the
second assertion without a finite degree cutoff.

Apply this to the entire actual polynomial companion (3).
Every term of w degree less than 3k contributes H degree at
most k-1, whereas the top term contributes 3beta*H^k.
Therefore

    T=Tr_(E/K0)(Q) in C(u)[H],
    deg_H T=k ODD, leading_H T=3beta!=0.                 (8)

In particular this trace is not zero. Nothing in (7),(8)
asserts that T belongs to B or is regular at finite u poles.

## 3. If R2=f(u), trace of the actual volume vanishes

Suppose for contradiction that the whole remainder has w degree
zero, namely d=e=0 in (1) and R2=f(u).

The root trace in (6) is s1=-a2/a3, a rational function of u.
Taking the trace of the actual two-form gives

    Tr_(E/K0)(du wedge dw/M)
                   =du wedge d(s1)/M=0.                 (9)

This identity can be checked by summing the three conjugate
roots of the cubic. Derivations extend uniquely through a finite
separable characteristic-zero field extension, so differentiation
commutes with this FIELD trace.

Now P=H^2+f(u) belongs to K0. Thus tracing the actual Jacobian
identity in (1) gives

    dP wedge dT=0.

In the independent coordinates u,H this is the exact equation

    f'(u)*T_H-2H*T_u=0.                                 (10)

All u derivatives in (10) keep H fixed.

Assume first that f is nonconstant. Write
T=sum_(j=0)^k t_j(u)H^j, t_j in C(u), with
t_k=3beta!=0. The coefficient of H^m in (10) is

    (m+1)*f'*t_(m+1)-2*t'_(m-1)=0,
    t_(-1)=0.                                          (11)

At m=0 this gives f'*t1=0, so t1=0.
At m=2 it gives 3f'*t3=2t1'=0, so t3=0.
Inductively every odd coefficient vanishes. This contradicts
the nonzero odd leading term in (8). Rational coefficient poles
in u do not affect the argument because C(u) is a field and
f' is a nonzero element.

If f is constant, then dP=2H*dH. H is a nonconstant polynomial
on the original whole source. More precisely its restriction
to the central source line is nonconstant affine by PF171,
so it has a zero there. At that source point dP=0, impossible
for a coordinate of a Keller map. Equivalently the source
Jacobian is divisible by the nonunit H.

Both alternatives are excluded. Hence

    R2 is not in C[u],
    deg_w R2 in {1,2}.                                  (12)

The same conclusion holds if one initially allows f in C(u):
the nonconstant case uses only f'!=0, and the constant case
is unchanged. For the actual regular normal form (1), f is
already polynomial.

## 4. Consequences and exact remaining scope

The theorem excludes every constant-plus-cubic-square shape
P=H^2+f(u) satisfying the actual polynomial companion conditions,
for all k>=3 odd and all coefficient degrees. It does not
exclude rational-in-w companions of arbitrary form.

In the canonical regular remainder R2=d*z+e*w+f,
either d!=0, giving w degree two, or d=0 and e!=0,
giving w degree one. The field escapes H not in F and
R2 not in F from PF171 remain intact. The trace argument does
not promote either one to a target coordinate or an integral
element of the original full finite normalization.

The checker replays the actual cubic multiplication traces,
their complete degree-leading recurrence, the volume-trace
identity and the odd-coefficient recursion. The checker passes 89 exact controls, which
support the written inductions (6),(7),(11); bounded moment or
coefficient sampling does not prove the theorem.

Replay with

    .venv/bin/python research/jacobian/affine-complements/verify_ruling_six_cubic_trace_remainder.py

No new external classification is imported. All original primes
and source/image distinctions remain retained. Entire ruling-six
closure, derivative integrality, automatic trace perfection,
normal-tail termination and HC4/JC2 remain open. Independent
whole-proof review, formal verification and external human
review remain unclaimed.
