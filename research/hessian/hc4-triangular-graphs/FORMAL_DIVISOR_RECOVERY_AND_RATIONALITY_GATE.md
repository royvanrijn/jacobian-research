# Divisor inverse jets are automatic; rational recovery is the global gate

15 November 2026. Claim **HC4MYJB0**. This is an all-order consequence
of the closed-divisor restriction in JA0, with a quantitative global
recovery criterion. The thickening statement is a reproduction, in
this setting, of Stacks Lemma 10.150.6; it is not claimed as a new
general formal-etale theorem. HC4 and the genuinely non-descended
embedding problem remain open.

## 1. Every divisor thickening already has its inverse

Let K:A2_z -> A2_T be complex plane Keller, R=C[z1,z2], and
B=C[T1,T2], embedded through K. Suppose f is nonconstant and
B -> R/(f) is surjective. Retained JA0 gives

    ker(B -> R/(f))=(a),    a(K)=f g,    (f,g)=R.          (1)

Then for every m>=1 the natural map is an isomorphism

    B/(a^m) -> R/(f^m).                                 (2)

This keeps all multiplicities in f. Consequently the inverse limits
are isomorphic as adic rings, with the original polynomial map as
their structural map. There are compatible polynomial coordinate
recoveries q_i,m(T) with

    q_i,m(K)-z_i in (f^m),                              (3)

unique modulo (a^m). They also satisfy K(q_1,m,q_2,m)=T modulo a^m.
No extra hypothesis on the free-variable degree of the HC4 embedding
enters (2).

Here is an elementary proof using (1), so no new unproved formal
lifting theorem is required. The element g is a unit modulo every
f^m. Indeed choose u,v with u g+v f=1; a finite geometric sum gives

    u g * sum_{j=0}^{m-1}(v f)^j = 1-(v f)^m.             (4)

For surjectivity, assume it for m and take r in R. Choose b in B with
r-b(K)=f^m c. Since g is a unit modulo f, choose d in B whose image
modulo f is c/g^m. Then

    r-b(K)-a(K)^m d(K) in (f^(m+1)).

For the kernel, if b(K) is in (f^m), the case m=1 gives b=a b1.
Cancellation of f and invertibility of g modulo f^(m-1) show
b1(K) in (f^(m-1)); induction gives b in (a^m).
The converse containment follows from a(K)=f g.
The isomorphisms are the natural ones and commute with reduction;
taking inverse limits proves the completion assertion. Applying
injectivity in (2) to K_i(q_m)-T_i proves the last assertion.

This is the principal-divisor proof of the standard thickening fact
in [Stacks 00UP, Lemma 10.150.6](https://stacks.math.columbia.edu/tag/00UP).
The source statement is explicitly attributed and independently
reproduced above. It is not an additional external assumption.

The extra divisor g=0 is absent from this entire formal neighborhood:
g is invertible in the f-adic completion by (4). Existence of inverse
jets, at any finite order or compatibly at all orders, therefore adds
no constraint to JA0. This does not say the coefficients of a formal
inverse contain no global information. Rationality of that inverse
is precisely a global condition which (2) does not supply.

## 2. Rational recovery, and an explicit degree threshold

Put d=max(deg K1,deg K2), using total source degree, and e=deg f>=1.
Fix a source coordinate z_i. If target polynomials P!=0 and Q obey

    P(K) z_i-Q(K) in (f^m),    deg P,deg Q<=N,            (5)

and

    m e > 1+N d,                                        (6)

then P(K) z_i=Q(K) identically. Indeed a nonzero polynomial divisible
by f^m has degree at least m e, whereas the left side of (5) has
degree at most 1+N d. Dominance of K makes P(K) nonzero.

If (5)--(6) hold for both source coordinates, their rational recovery
makes C(z)=C(K). The retained Keller denominator-removal argument
of HW0 then gives a polynomial inverse. Equivalently, compatible
formal coordinate inverses which are represented by fixed rational
functions have a polynomial inverse in this Keller setting: cross
multiply, use (3), and use the elementary separation
intersection_m (f^m)=0. Rational functions here mean fixed pairs
P,Q with P!=0 satisfying the cross-multiplied identities at every
order; no division in a potentially reducible completion is assumed.

A uniform bound N on numerator and denominator degrees for arbitrary
large orders m therefore forces the inverse. No compatibility among
the chosen pairs P,Q is needed for this implication: take an order
large enough for (6). Conversely an actual polynomial inverse gives
uniformly bounded polynomial representatives.

Under hypothetical nonautomorphy, at least one z_i is not in C(K).
For that coordinate every rational approximation (5) necessarily has

    N >= ceil((m e-1)/d).                               (7)

Thus (2) automatically supplies polynomial inverse jets in all
orders, but at least one coordinate requires target degrees growing
at least linearly. This is an all-degree obstruction to bounded
rational reconstruction, not a proposal to compute longer jets.
The estimate is sufficient, not an optimal bound in every example.
Neither the actual closed embedding nor its fixed polynomial recovery
formulas are proved to bound these target degrees.

## 3. One explicit model realizes all the formal and boundary warnings

Retain the JA0 comparison map on the etale source open

    L(s,t)=(s,t+s t^2),   J=1+2s t,   U=D(J).
    f=t,   a=T2,   a(L)=t(1+s t)=f g,   1=g-s f.        (8)

Its entire polynomial extension on A2 is not Keller. On U the whole
selected divisor C:t=0 maps isomorphically onto the target line,
and its other inverse-image component is D:1+s t=0. They are
disjoint; D has coordinate s in Gm, with t=-1/s and J=-1.
The normalized map D -> A1 is s -> s, degree one and etale,
with two punctures, genus zero, and finite-image boundary index one.

At s=0 the target tends to (0,0), but f=-1/s has a pole.
The zero of this nonconstant unit is instead at s=infinity, which
maps to target infinity. Thus the forced unit's zero and the
finite-image boundary puncture need not coincide under the section
and etale hypotheses. The example does not satisfy the whole-A2
Keller hypothesis or give an HC4 counterexample.

Its formal inverse along T2=0 is

    s=T1,
    t=sum_{n>=1} (-1)^(n-1) Catalan_(n-1) T1^(n-1) T2^n
     =(-1+sqrt(1+4 T1 T2))/(2 T1).                       (9)

The apparent denominator in the displayed expression is removed by
the series. Truncating at n<m gives a polynomial inverse modulo
t^m, including both compositions. The recurrence for Catalan
numbers verifies t+T1 t^2=T2 coefficient by coefficient for every
order. The discriminant 1+4 T1 T2 is irreducible, hence has odd
valuation on its prime divisor and is not a square in C(T1,T2).
Thus this inverse branch is algebraic of degree two and nonrational.

For m>=2 the unique representative of (9) having T2 degree <m has
total degree 2m-3. Adding a multiple of T2^m cannot cancel its
nonzero monomial T1^(m-2)T2^(m-1); every polynomial representative
has at least that total degree. This proves its unbounded target
degree without an asymptotic experiment. Replacing f by t^r and a
by T2^r gives the same conclusion for all fat-divisor orders.

## 4. Consequence for the next HC4 step

The actual rank-one embedding gives (1) by IX1 and JA0. Formal
neighborhood recovery, conormal lifting and arbitrarily many normal
jets are already automatic; they cannot be added as an unproved
global-generation premise. A useful next implication would force
fixed rational recovery, target trace descent, or full source-prime
saturation from the actual globally closed, w-dependent embedding.

The attempted bounded-image argument does not do this: at a
finite-target escape, a pole in f makes f^2+w dot K diverge for
bounded w. The model above exhibits exactly that pole placement
under the local section hypotheses. It does not settle what the
full embedding forces.

The fibre-affine class has already been classified in HT1 and is
not a new all-degree search space. Nonconstant rank-two/rank-three
helper traces and intrinsic ambient alignment remain alternative
structural gates. No quadratic-symmetry rectification or transfer
of normality between different implicit images has been proved.

The checker verifies exact finite controls of these formulas and
pins prior proof/source bytes. The all-order statements are proved
in Sections 1--3; their validity is not inferred from the checked
orders. No prior receipt is regenerated.
