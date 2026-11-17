# Weighted period recovery and the global nilpotence gap

Canonical claims HC4MYKH0/1/2; mathematical authority is
[MATH_STATUS.json](../../MATH_STATUS.json). These are all-degree structural
results. They recover the full finite plinth scheme from weighted periods,
identify exactly what a Morse degeneration forces on that reconstruction,
and show why local Morse and finite-jet data do not prove reducedness.
No native action, source alignment or HC4 closure is constructed.

## 1. Hypotheses and retained inputs

Use the actual A4 source and the hypotheses of
[KF0/1](FLAT_PLINTH_FAMILIES_AND_ONE_VALUE_RECOGNITION.md):

    R=C[X1,X2,X3,X4], 0!=delta in LND(R),
    A=ker(delta)=B[h], B=C[b1,b2],
    J=pl(delta), sqrt(J)=I_Z A,
    Z={b_i:1<=i<=n} finite, reduced, nonempty.

All three kernel generators are algebraically independent. The entire
plinth ideal need not descend from B. Put D=C[h], K=C(h), U=A2 minus Z.
KF proves the following facts, which are inputs here, not new assumptions:

* One Cech class xi in H1(U,O) tensor D has exact annihilator J.
* T=A/J is finite free over D, with component ranks L_i and L=sum L_i.
* On the entire source, delta=c N for c in C*, where
  N(f)=Jac(b1,b2,h,f) in the actual source coordinates.
* On some dense value open S, the whole X_S is the actual Ga torsor
  over U times S. Its generic H3 has the constant point-link basis.

For the native cycle-exhausted functions, [KA1](GLOBAL_INTEGRAL_VANISHING_CYCLES_AND_VOLUME_PERIODS.md)
and KF give n=N, a bijection between support points and Morse points,
and primitive point-link projections of the finite vanishing cycles.
The native unweighted periods are exactly nonzero constant multiples of
h-c_i, including critical-value collisions. Use
[KE0](FINITE_PLINTH_DESCENT_AND_VANISHING_CYCLE_RECOGNITION.md)'s finite
negative-Laurent pairing at each support point. No closed-value kernel
or action-plinth base change is assumed.

## 2. Weighted periods recover the full plinth: HC4MYKH0

Let Omega=dX1 dX2 dX3 dX4 and let omega=Omega/dh, with convention
dh wedge omega=Omega. On a local-slice chart choose delta(t)=1. Since
delta=c N, the exact Jacobian calculation gives

    Omega=c db1 db2 dh dt,
    omega=c dt wedge db1 wedge db2 on a value fibre.

For any f in B[h], a relative primitive is c f t db1 wedge db2:
terms from df wedge db1 wedge db2 vanish on that fibre. On two charts
the primitives differ by c f xi db1 wedge db2. Cech--de Rham Stokes
integration over a positively oriented point link therefore gives

    Pi_i(f;h):=integral_(gamma_i(h)) f omega / (c(2pi i)^2)
              = Res_(b_i)(f xi_i db1 db2).                 (1)

Here gamma_i projects to the positive integral link of b_i. Reversing
that orientation changes both sides' common sign. The residue is the
coefficient of u^-1 v^-1 in parameters
u=b1-b_i1,v=b2-b_i2. One may verify the Stokes normalization on the
two-coordinate annular cover: the difference of the primitives is
integrated over |u|=epsilon,|v|=epsilon. The contractible affine-line
torsor identifies that boundary calculation with its three-cycle.
The same formula follows by pulling back a smooth section over the
point link and applying the two-chart Stokes formula there. This uses
the whole generic torsor, not only a dense subvariety of one fibre.

After decomposing support, xi_i has finitely many negative Laurent
coefficients, each in D. Formula (1) makes every Pi_i(f;h) a polynomial
in h, even if the cycle calculation initially used S. It is an exact
identity there, hence throughout its polynomial continuation.

The entire scheme can now be read from the weighted periods:

    J={a in A: Pi_i(a g;h)=0 for every g in B and every i}. (2)

Indeed a nonzero Laurent coefficient of a xi_i is detected by
multiplication with the appropriate nonnegative monomial u^r v^s.
Conversely all such residues vanish when a xi=0. This is the same
elementary coefficient duality as KE, now identified with actual
weighted source-volume periods, not just an abstract inverse system.

Over K put lambda_i(f)=Pi_i(f;h). The Hankel multiplication pairing

    H_i(f,g)=lambda_i(fg)

has radical J_i,K and induces a perfect pairing on T_i,K. In particular

    generic rank(H_i)=L_i.                              (3)

The rank means the supremum of ranks of finite moment matrices; a
monomial basis of the finite quotient attains it. If verified slice
powers u^r,v^s annihilate xi_i, use the rectangle
{u^a v^b:0<=a<r,0<=b<s}. Its finite moment matrix has exactly this rank.
There is no bound on r,s or the degree of the unknown quotient generators.
This is a finite symbolic certificate once the actual powers are supplied,
not an unproved uniform cutoff or a coefficient-degree search.

Consequently the following are equivalent for a candidate action:

    every generic H_i has rank one;
    L_i=1 for all i;
    J=I_Z A;
    the full cotangent source frame is recovered by KD/KF. (4)

When lambda_i(1)!=0 in K, rank one is equivalently the exact covariance
identity

    lambda_i(fg) lambda_i(1)=lambda_i(f) lambda_i(g)      (5)

for all f,g, or for all pairs in a supplied spanning rectangle. The
support is the fixed point b_i, so the resulting character is evaluation
there. For native cycle-exhausted functions, (4) supplies pure alignment
with R_(N,c); unit existence still has the separate E_N conjunct.

The criterion uses generic symbolic rank over C(h). A rank drop of a
specialized period matrix, even at a regular value, is not asserted to
be the length of the actual ideal fibre. That distinction remains the
specialization issue proved in KF.

## 3. The leading Morse moment is automatically reduced: HC4MYKH1

Now use a native Morse point p_i, with h(p_i)=c_i and b(p_i)=b_i.
Choose the local vanishing sphere and transport it to the point-link
cycle above, choosing its orientation accordingly. For any holomorphic
weight F near p_i, the elementary Morse scaling formula in four variables is

    integral_(gamma_i(h)) F omega
      = C_i (h-c_i) [F(p_i)+O(h-c_i)], C_i!=0.         (6)

To see this, make analytic Morse coordinates h-c_i=sum z_j^2 and write
Omega=phi(z) dz1...dz4 with phi(0)!=0. Set z=sqrt(h-c_i) theta on
the real unit three-sphere. The relative form scales by h-c_i. Odd
homogeneous Taylor terms of F phi integrate to zero under theta -> -theta;
the constant term is a nonzero sphere volume times F(p_i)phi(0).
The convergent even Taylor series supplies the displayed holomorphic
remainder. This proves (6) directly; no isochore classification is imported.

Apply (6) to every quotient polynomial f(b). Formula (1), and detection
of every Laurent coefficient by a monomial, give

    xi_i=(h-c_i) eta_i,
    eta_i at h=c_i = a_i/(uv), a_i in C*,              (7)
    ann(eta_i)=ann(xi_i)=J_i.

The last equality is over D: the inverse-system module is value-torsion
free, so multiplying its element by h-c_i changes no annihilator.
The first equality also holds over D, coefficient by coefficient.
For the native exact linear unweighted period, Res eta_i=a_i at every
value; its other coefficients vanish at c_i. Coincident critical values
cause no difficulty because the support components remain distinct.

Thus the normalized leading Hankel matrix at the critical value ALWAYS
has rank one. It does not test L_i=1. More precisely choose a finite
free D-module W_i of negative Laurent monomials containing all eta_i
and all their B-multiples that survive. Such a module exists because
fixed powers kill the class. The injective residue embedding is

    T_i -> W_i, f mod J_i -> f eta_i.                  (8)

After h=c_i, (7) identifies its kernel exactly as Nil((T_i)_(c_i)),
of dimension L_i-1. Its image has dimension one. This is compatible
with T_i itself being finite free of rank L_i: it is the embedding
that loses rank on specialization.

Equivalently, for Q_i=coker(8), the exact tensor sequence gives

    Tor_1^D(Q_i,D/(h-c_i)) = Nil((T_i)_(c_i)).          (9)

In particular a sufficient and, under these hypotheses, equivalent
minimum-length test is that this normalized residue embedding have
saturated image at h=c_i. It is NOT the already-proved flatness of
T_i. A new global action argument eliminating the cokernel's value
torsion would force L_i=1. Equations (7)--(9) expose the precise
specialization defect instead of mistaking its reduced leading class
for a reduced actual ideal fibre.

## 4. Arbitrary thickness with genuine local Morse periods: HC4MYKH2

The following family proves that local Morse-volume data, even to every
order, do not by themselves eliminate that defect. Fix ANY integer m>=1.
In C[u,v,P,Q] put

    q=uP+vQ, sigma(P)=P+P^(m+1),
    phi(P)=sigma'(P)=1+(m+1)P^m,
    Omega_m=phi(P) du dv dP dQ.

Near the origin, Omega_m is a nonvanishing holomorphic volume and q
has its unique nondegenerate Morse point. For all polynomial weights
f(u,v), the exact local period on the shrinking quadratic sphere is

    integral f Omega_m/dq / integral Omega_0/dq
        = f(0,0)+s^m [u^m v^0] f, q=s.               (10)

For a direct computation put P=bar(u),Q=bar(v) on the q=s sphere for
positive small s. The induced relative-volume density is a nonzero
constant times the rotationally invariant sphere measure. Separate
phase rotations of u,v kill all monomials except the constant and
u^m P^m. On the unit complex two-sphere |u|^2 is uniform on [0,1],
so its m-th moment is 1/(m+1). The coefficient m+1 in phi cancels
that denominator. Scaling supplies s^m. Analytic continuation proves
(10) as an exact holomorphic identity. In particular the unweighted
period is EXACTLY linear in s, at every m.

The normalized inverse system in (10) is

    eta_m=1/(uv)+s^m/(u^(m+1)v),
    xi_m=s eta_m,
    ann_D[u,v](eta_m)=(u^(m+1),v).                    (11)

The quotient is free of rank m+1 over D=C[s], has fixed support u=v=0,
and specializes to that SAME length at s=0. Independence of the
classes 1,u,...,u^m follows by extracting their distinct negative
Laurent coefficients over the domain D. Nevertheless eta_m at0 is
1/(uv), with reduced annihilator (u,v).

The full (m+1)-by-(m+1) moment matrix on 1,u,...,u^m is

    H_m=s^m J_(m+1)+E_00,
    det H_m=(-1)^(m(m+1)/2) s^(m(m+1)).               (12)

Here J_(m+1) is the reversal matrix. The determinant follows from
J^2=1, det J=(-1)^(m(m+1)/2), and J_00=0 for m>=1: the rank-one
update E_00 changes no determinant. It has generic rank m+1 and
special rank one. For any fixed finite Taylor order, m may be taken
larger, making these weighted data indistinguishable to that order
from the reduced control. There is no uniform finite-jet reducedness
test based solely on this local Morse-volume information.

Even requiring constant volume in some analytic Morse coordinates
does not remove the example. Set a=(m+1)/m and make the local analytic
change

    U=u exp(-aP^m), Z=P exp(aP^m), V=v, W=Q.

It preserves q=UZ+VW and its Jacobian is exactly phi(P). Thus Omega_m
is the constant volume in U,Z,V,W, with the harmless constant sign
adjusted by coordinate order. The weights u,v become analytic
functions in these new coordinates. This is an explicit local
inverse-function construction, not a polynomial source automorphism
or an import of an isochore Morse theorem.

### The precise missing global additive action

On the smooth factorial principal open

    X_m=Spec C[u,v,P,Q,phi(P)^-1]

the primitive volume-normalized field fixing u,v,q is, up to orientation,

    E=(v partial_P-u partial_Q)/phi(P).                (13)

It has polynomial embedded kernel A=C[u,v,q]. To justify that equality,
the standard field has invariant field C(u,v,q) and polynomial kernel
C[u,v,q]. An invariant rational function regular on X_m cannot acquire
new poles: each possible source prime pole would divide phi(P), but
none of its linear factors P-r is invariant under v partial_P-u partial_Q.
Factorial closure, or direct differentiation of a reduced denominator,
excludes such invariant poles.

One can check even its exact slice-value ideal, despite E not being
locally nilpotent. Define J_E=E(O(X_m)) intersect A; call this a
slice-value ideal, not the plinth of a Ga action. On v!=0 a local time is
t_v=sigma(P)/v. On u!=0 use

    t_u=[sigma(P)-sigma(q/u)]/v.

The latter is in C[u,u^-1,v,P,Q]: substitute q=uP+vQ and expand the
binomial difference to cancel v. Both have E(t)=1. Their difference is

    t_v-t_u=q/(uv)+q^(m+1)/(u^(m+1)v)=xi_m.            (14)

Their localized kernels in the ORIGINAL polynomial ring are A_u,A_v.
The two-chart intersection argument therefore gives J_E equal to the
annihilator in (11), first for polynomial primitives. It also gives
the same result on X_m: if r had a pole of order k>0 at a root of phi,
E(r) would have pole order k+2 there, with nonzero leading coefficient
on the dense set v!=0. Thus E(r) in A forces r to have no such pole.
Consequently

    J_E=(u^(m+1),v) A,

an exact fixed-cylinder, finite-free, arbitrarily thick slice-value
scheme on this actual principal open source. These are not native
Ga quotients. In particular no local-slice polynomial-ring theorem
may be inferred solely from the identities E(t)=1.

The obstruction is explicit. Write E^n(P)=v^n f_n(P), where

    f_1=1/phi, f_(n+1)=f_n'/phi.

Every root r of phi is simple. Inductively f_n has a pole of order
2n-1 at r, with nonzero leading coefficient, so E^n(P) never vanishes.
E is NOT locally nilpotent on O(X_m). At every finite jet at the
origin its induced linear operator IS nilpotent: it preserves the
maximal ideal, its associated graded is the linear nilpotent root
field, and the higher terms strictly raise degree. The finite filtered
quotients consequently all exponentiate to compatible polynomial
Ga actions on those fat points. Such jets do not give a global LND.

This example has the wrong global source and no algebraic Ga action;
it is not a counterexample to KF or to the open native minimum-length
claim. The actual A4 hypotheses include constant global units and
global local nilpotence, both absent from this attempted substitution.
It isolates the extra global information a successful proof must use.

## 5. Source audit, assurance and next direction

The inspected [Masuda preprint](https://arxiv.org/html/2512.06687v1)
has relevant quotient criteria in Proposition2.5 and Theorem4.3, with
surjectivity/equidimensionality, and in the latter a free action and an
A3-fibration, among their hypotheses. Only its introduction,Section2
through Proposition2.5 and the Theorem4.3 statement were inspected here.
No theorem from it is imported. The native generic quotient instead
misses the entire fixed point lines on a dense value open. Finite
freeness of A/J over C[h] is not faithful flatness of the source over
its invariant ring, and a native Morse fibration is not an A3-fibration.
The [Francoise bibliographic page](https://www.numdam.org/item/AST_1978__59-60__119_0/)
was revisited only as context; no full-paper local normal-form theorem
was inspected or used. The explicit change above supplies what is needed.

The [checker](verify_weighted_period_recovery.py),
[receipt](receipt-weighted-period-recovery.json) and
[source context](sources/weighted-period-recovery-context.json) bind
retained inputs and verify residue pairings,rank defects,local times,
the analytic volume Jacobian and the arbitrary-m recurrence and matrix
identities. They support the written all-degree arguments rather than
enumerating quotient degrees. Independent replay, formal verification
and external review are not claimed.

The previous turn's [KG exclusion](THREE_POINT_TRANSLATION_QUOTIENT_BARRIER_AND_MORSE_COUNTERCONTROL.md)
still requires a genuinely new embedded quotient field for S_rho. For
a candidate with KF fixed support, use either the actual one-value ideal
length test or the equivalent generic rank-one weighted moment test.
A structural proof can now target saturation of the NORMALIZED residue
embedding at each Morse value; this is the exact defect left by (9).
It must use a global algebraic action, not only local Morse-volume
geometry, formal additive jets or a vector field with local time functions.
No such saturation theorem, native action or alignment is proved here.
E2/E3, both universal endpoint inversions, full graph embedding recovery
and unrestricted HC4/JC2 remain open. The original goal is not complete.
