# Integer degree growth, target framing and the positive-thinness boundary

2026-11-16. HC4MYMG0. All complex plane polynomial Keller maps, all degrees.

FULL THEOREM TARGET: prove HC4/JC2 or construct an actual counterexample.
CURRENT LAST UNSOLVED IMPLICATION: ACTUAL polynomial plane source forces
Tr(Lambda*Lambda) subset B, or full marked-energy integrality. OPEN.
External review/outreach remains DEFERRED.

This gives a published-theorem corollary about iteration, an elementary
all-degree target framing, and an exact comparison with EVERY omitted
prime of the full finite normalization. No novelty priority is claimed
for the dynamical-degree corollary or the framing mechanism.

## 1. Published input and reading boundary

We use Charles Favre and Mattias Jonsson, *Dynamical compactifications
of C2*, Annals of Mathematics173(2011),211-249,
DOI10.4007/annals.2011.173.1.6:
[publisher PDF](https://annals.math.princeton.edu/wp-content/uploads/annals-v173-n1-p06-p.pdf).

The inspected published package consists of:

* the general inequality lambda_2<=lambda_1^2 and Proposition2.3's
  eigenvaluation in V1;
* Proposition2.5's integer/irrational eigenvaluation dichotomy when
  lambda_2<lambda_1^2, and Proposition2.6;
* the thinness definition in Section1.3: for v_E=ord_E/b_E,
  A(v_E)=(1+ord_E(dx wedge dy))/b_E; irrational valuations have
  irrational thinness; V1 has A<=0 and skewness alpha>=0;
* the Jacobian formula(5.2),
  A(nu)+nu(JF)=d(F,nu) A(F_bullet nu);
* Propositions5.1-5.3: in the maximum case lambda_2=lambda_1^2,
  unbounded normalized degree growth gives a polynomially conjugate
  skew product, bounded growth with lambda_1>1 gives a PROPER map,
  and bounded iterate degrees give an automorphism.

The pertinent text and proofs on printed pp216-221,224-235 were read;
printed pp220-222 were also inspected as rendered pages. The long
geometric machinery and the preceding 2007 *Eigenvaluations* paper
are NOT independently reproduced here. This is an explicit published
proved import, not certification by the elementary checker below.
The source PDF is byte-pinned in the context. All source conventions
are those of the published 2011 version.

Write lambda=lambda_1(F)=lim deg(F^n)^(1/n) and d=lambda_2(F) for the
generic sheet degree. These are DIFFERENT quantities.

## 2. Every plane Keller map has integer dynamical degree

Let F=(P,Q) be polynomial Keller. If lambda=1, then d<=lambda^2=1;
the retained birational Keller theorem gives an inverse, so this
case already has integer lambda.

Suppose lambda>1. We first rule out d=lambda^2. If normalized degree
growth deg(F^n)/lambda^n is unbounded, Proposition5.1 makes F
polynomially conjugate to a skew product (p(x),q(x,y)).
Its Jacobian p'(x) q_y(x,y) is a nonzero constant. Both factors are
therefore units: p is affine-linear and q=a y+b(x), a!=0.
This is an automorphism with lambda=1, a contradiction.
If normalized growth is bounded, Proposition5.3 makes F proper.
A proper Keller map is finite etale; the complex affine plane has
no nontrivial connected finite etale cover. Thus F is an automorphism,
with d=1, again contradicting d=lambda^2>1.

Consequently

    lambda>1 ==> d<lambda^2.                              (1)

Proposition2.3 supplies an eigenvaluation nu_* in V1 with
d(F,nu_*)=lambda. Its thinness is finite, in [-2,0]. Since JF is a
nonzero constant, the Jacobian formula gives

    A(nu_*)=lambda A(nu_*), hence A(nu_*)=0.                (2)

An irrational valuation has irrational thinness, so nu_* cannot be
irrational. Proposition2.5 now gives

    lambda is a POSITIVE INTEGER for EVERY plane Keller F. (3)

When lambda>1 the eigenvaluation is divisorial or infinitely singular.
If divisorial, it has alpha>0 and A=0: alpha=0 in the quasimonomial
part of V1 is a rational pencil valuation with A<0. Proposition2.8
explicitly identifies the alpha>0,A=0 eigenvaluation situation as a
potential Jacobian-counterexample situation. It does NOT exclude it.

For a nonautomorphic F, d>1 and (1) is strict. This is a necessary
all-degree dynamical restriction, not a new degree exclusion or a
proof that lambda=1. No claimed Lyapunov or entropy inequality beyond
the cited statements is used.

## 3. Arbitrarily fast degree growth with EXACTLY the same finite cover

This construction works for ANY dominant polynomial pair F=(P,Q);
Keller is only needed for the later application.

Set p=deg P, q=deg Q>=1, and write Q_q for Q's leading homogeneous
form. Choose b in C with Q_q(1,b)!=0. Such b exists because Q_q is a
nonzero binary form. For ANY integer M>=2 with Mq>p define the
determinant-one target automorphism

    H_M(s,t)=(s+t^M, b(s+t^M)+t).
    H_M^-1(S,T)=(S-(T-bS)^M, T-bS).

Put G=H_M composed with F and D=Mq. Its degree-D homogeneous part is

    G_D=v h,  v=(1,b), h=Q_q^M, h(v)=a!=0.

For n>=1 the highest homogeneous part of G^n is

    c_n v h^(D^(n-1)),
    c_1=1, c_(n+1)=a c_n^D.

It is nonzero for EVERY n. This proves, by induction without expanding
iterates,

    deg(G^n)=D^n, lambda_1(G)=D.                           (4)

The choice of M makes D arbitrarily large. Since H_M is a target
automorphism, however,

    C[G1,G2]=C[P,Q]=B, C(G1,G2)=C(P,Q)=K,
    d(G)=d(F), S_G=S_F, Lambda_G=Lambda_F                 (5)

as LITERAL embedded objects in the fixed source field L.
In particular Tr(Lambda*Lambda) is literally unchanged.
Also Jac(G)=Jac(F). This is postcomposition, not dynamical conjugacy;
there is no claim that dynamical degree is invariant under it.

For Keller F choose M even larger so D>d. Proposition2.6 excludes a
divisorial eigenvaluation for G; (2)-(3) exclude the irrational case.
Thus its eigenvaluation is INFINITELY SINGULAR, exactly as for the
familiar polynomial automorphism case. This applies conditionally to
every hypothetical counterexample as well. Neither high integer
growth nor this eigenvaluation type supplies an inverse certificate.

## 4. Every omitted finite-normalization prime has POSITIVE thinness

Now use the ACTUAL source X=A2 and its full finite normalization
Z=Spec(S) over the target A2. Keep ALL omitted primes E in Z-X,
including omitted unramified primes. Let e_E>=1 be the ramification
index over its target image curve. Put v=ord_E as a valuation of
L=C(x,y).

The center of v is outside the source A2. Therefore a generic affine
source function has a pole of some integer order b_E>0, and the
normalized valuation is nu_E=v/b_E in V0. Were both source
coordinates regular at its center, it would have a source A2 center;
uniqueness of its center on Z would place it in the source open,
a contradiction.

At the generic point of E the finite map is a tame DVR map.
The retained KN0 differential calculation gives

    ord_E(F^*(dP wedge dQ))=e_E-1.

Since F^*(dP wedge dQ)=Jac(F) dx wedge dy and Jac(F) is a unit,

    ord_E(dx wedge dy)=e_E-1,
    A(nu_E)=e_E/b_E>0.                                    (6)

In this formula dP wedge dQ denotes the target coordinate volume
form before pullback. The target residue extension is separable;
no wild or codimension-two formula is assumed.

Equation(6) holds for EVERY omitted finite-normalization prime.
For ramified primes e_E>1 these are precisely where the codifferent
has order -(e_E-1); unramified omitted primes have e_E=1 and are
still of positive thinness. No prime is deleted to make the statement.

Since V1 requires A<=0, NONE of these divisorial valuations lies in
the tight tree V1. In particular none is the zero-thinness eigenvaluation
from (2). The tree deliberately used to control global degree growth
does not itself contain the omitted divisors carrying the finite trace
defect. This is not a claim that dynamics could never constrain them
indirectly; an additional bridge would have to be proved.

The separation is compatible with (5): target reframing changes degree
growth while preserving the entire finite normalization and its defect.

## 5. What this establishes and what remains

MG0 supplies the integer-degree corollary, the uniform target framing
with exact embedded-ring preservation, and the positive-thinness
location of ALL omitted primes. Published Favre-Jonsson statements
are explicitly credited. KN0's canonical ramification formula is
retained, not presented as a new local different theorem.

Neither a tight dynamical compactification nor a recurrence for iterate
degrees proves finite properness of the original map. Do not promote
the partial holomorphic-at-infinity conclusions for the strict regime
to an everywhere holomorphic finite compactification. Do not assert
that postcomposition leaves lambda_1 fixed, that infinitely singular
eigenvaluations force an inverse, or that V1 covers the omitted primes.

The new checker verifies symbolic parameter identities for H_M and
its inverse, the local tame-volume transcription and finite leading-term
controls. The all-degree conclusions come from the written proofs and
the declared published import. No registered mathematical checker is
rerun. Independent replay, formal verification and external review
remain false. Full HC4/JC2, the actual full trace-dual pairing, and
full marked-energy integrality remain OPEN.
