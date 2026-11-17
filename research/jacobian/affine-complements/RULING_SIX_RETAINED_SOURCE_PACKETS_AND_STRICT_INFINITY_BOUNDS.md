# Retained source packets give strict ruling-six infinity bounds

PF2RULED166. 17 November 2026. Actual complex polynomial plane Keller maps.

Every intermediate infinity prime shares its target image with a
different, unramified prime reached by the whole source. The complete
degree formula therefore improves the former infinity upper bound
e<=b to e*f<=b-1, keeping the residue degree over the image.

For the remaining ruling-six class this forces D>=22 and b>=5,
where D=deg b3 and b=[E:F]. The original generic map degree is
2b>=10. It also forces the product of the infinity u-index and
the residue degree over its image to be at most two. If that
product is ever two, then D>=42 and b>=15.

These are necessary conditions on the entire remaining class.
The quartic A faces, derivative integrality and JC2 remain open.

## 1. Actual complete normalization and notation

Keep the complete quadratic source image and fields of PF165:

    R=C[x,t], B=C[P,Q], F=C(P,Q), E=C(u,w), L=C(x,t),
    u=t^2, w^2=N+M*z, dP wedge dQ=du wedge dw/M,
    m=deg M>=2, D=deg b3, b=[E:F].                       (1)

Let S_E be the ORIGINAL full finite normalization of B in E.
No source point, ramified prime, unramified prime or P-vertical
prime is removed from it.

[PF165](RULING_SIX_SIMPLE_ESCAPE_GRAPHS_AND_QUARTIC_FACE_REDUCTION.md)
proves that only A pole divisors remain, that

    b3=M*h^2, h in C[u],

and that every reduced companion has degree 3k for odd k>=3.
Let the n_A>=1 actual A pole divisors have slopes L_1,...,L_(n_A),
and put S=sum_i L_i. From PF164/PF165,

    b<=sum_i(2L_i-5)=2S-5n_A,
    ord_(c_i) h=2L_i, D=m+2deg h>=m+4S.                 (2)

On every smooth complete geometric generic P fibre, the coframe
is dQ. At a point over u=infinity put u=s^-alpha. The full monic
parameter estimate gives

    ord_s(dQ)+1>=ceil(alpha*(D+3m-6)/6).                 (3)

These points are holomorphic for dQ, so P and Q are finite
there, whereas u has a pole. Each belongs to a genuine height-one
prime of S_E over a target curve with nonconstant P.

## 2. A different source prime lies above every such image

Let Dinf be one of these original infinity primes and Gamma its
irreducible target image. Write

    e=e(Dinf/Gamma),
    f=[kappa(Dinf):kappa(Gamma)].                        (4)

There is a DIFFERENT prime Dsrc of S_E over Gamma, of ramification
index one, which is reached by the whole source.

To prove this without assuming finiteness of the Keller map, let
H(P,Q) define Gamma. Its source pullback H(P(x,t),Q(x,t)) is a
nonconstant polynomial: P,Q are algebraically independent. Over C
it has a nonempty curve of zeros; choose an irreducible component C.

The actual Keller map is etale and quasi-finite. Hence C cannot
map to a point and instead dominates Gamma. Its generic target
valuation has index one, because an etale pullback of a reduced
target divisor is reduced.

Every element of S_E lies in R: it is integral over B subset R,
and R is integrally closed in L. Thus C has a centre Dsrc on
S_E, necessarily a height-one prime over Gamma. The ramification
indices in the valuation tower multiply. Since the index at C
over Gamma is one,

    e(C/Dsrc)*e(Dsrc/Gamma)=1,

so e(Dsrc/Gamma)=1. Its residue degree is a positive integer.
The function u=t^2 is regular at C, hence has nonnegative
valuation at Dsrc. It has negative valuation at Dinf. Therefore

    Dsrc!=Dinf.                                         (5)

This proof uses the whole affine source, not a hypothetical
source sheet selected from a comparison normalization. It does
not assert that Dsrc is the only other prime above Gamma.

## 3. The complete divisor degree formula is strict at infinity

Localize B at the generic point of Gamma. It is a DVR with
uniformizer H. The localized S_E is finite and torsion-free
of rank b over that DVR, hence free of rank b.

Reducing modulo H and decomposing over ALL primes above Gamma
gives the complete degree formula

    b=sum_(D over Gamma) e(D/Gamma)*f(D/Gamma).           (6)

Equivalently each local DVR factor with H a unit times tau^e
contributes e times its residue degree to the quotient dimension.
Both Dinf and the distinct Dsrc of (5) occur in this sum.
Consequently

    e*f<=b-1.                                          (7)

More strongly, all infinity primes over the SAME Gamma can be
counted together:

    sum_(Dinf over Gamma) e(Dinf/Gamma)*f(Dinf/Gamma)<=b-1. (7a)

The retained Dsrc is different from every one of them because u
is regular there. Every other prime remains in (6); keeping additional primes
can only strengthen (7). In particular the estimate does not
depend on an assumption that the boundary has only one prime.

To connect e with the geometric generic coframe, note that P is
nonconstant on Gamma. On the generic P line its equation is a
separable polynomial in Q over C(P). After algebraic closure its
roots q_j(p) are distinct; near each such point Q-q_j(p) is a
target uniformizer. At a geometric point of Dinf it pulls back
with order e. Separable constant-field base change splits the
residue points and does not change this ramification index.
Therefore

    ord_s(dQ)=e-1.                                      (8)

The number f in (4) is the residue degree over the IMAGE CURVE
Gamma. It is not the degree of kappa(Dinf) over C(P).
These two degrees must not be interchanged.

Combining (3),(7),(8) gives the strict full infinity estimate

    f*ceil(alpha*(D+3m-6)/6)+1<=b.                       (9)

Since alpha>=1 and f>=1, this implies the simpler bound

    ceil((D+3m)/6)<=b.                                  (10)

## 4. Infinity index times image residue degree is at most two

Set rho=alpha*f. The elementary inequality
f*ceil(alpha*c)>=ceil(rho*c), for c>0, applies to (9).
If rho>=3, then

    b>=ceil(3*(D+3m-6)/6)+1
      >=D/2+3*(m-2)/2+1,

contradicting PF164's b<=floor(D/2)-5.
Thus every original intermediate infinity prime satisfies

    alpha*f<=2.                                        (11)

Applying the same argument to (7a) also gives

    sum_(Dinf over Gamma) alpha_D*f_D<=2                  (11a)

for each fixed infinity image Gamma. In particular no three
index/residue-one infinity primes share the same target image.

The only individual possibilities are

    (alpha,f)=(1,1),(1,2),(2,1).

In particular an index-two infinity prime is birational on its
target image, while an index-one prime has image residue degree
at most two. Nonlinear residue extensions over C(P) remain
possible in every statement here. No graph-image conclusion
follows from (11).

## 5. Nonlinear image packets of total u degree six

Put g_D=[kappa(Gamma_D):C(P)]. We have g_D>=2.

Here is a direct justification using the retained PF159 polynomial
leading-coefficient lemma, rather than a general classification of
nonproper curves. At an infinity branch put u=s^-alpha. The first
coefficient of w outside C occurs after a finite Laurent truncation
with C coefficients; write it lambda. Expanding the whole P*
polynomial there gives H(lambda)=p with H in C[Y] nonconstant.
The same full Taylor argument says that the finite Q residue is
G(lambda), G in C[Y], or zero: a negative minimum order would
give a pole because lambda is transcendental over C.

Thus the target image field is represented by C(H(lambda),G(lambda)).
If g_D were one, then G belongs to
C(H(lambda)) intersect C[lambda]=C[H(lambda)].
Indeed a finite pole of a reduced rational function of H pulls back
to a root of H-c and cannot become a polynomial in lambda.
So the target image would be Q=q0(P) for a polynomial q0.
It is the whole closed graph, as in PF165's full-normalization
argument. Every original source prime over this infinity prime
has t a pole and is omitted. The retained Chau component theorem
excludes its graph image. This proves g_D>=2.

This argument does not identify lambda with a generator of
kappa(Dinf); the embedding of the IMAGE field is sufficient.

The degree of u on the full intermediate generic P curve is six.
The complete valuation degree formula over C(P) therefore gives

    sum_D alpha_D*f_D*g_D=6.                             (11b)

Let r_D=alpha_D*f_D in {1,2}. The possible MULTISETS of (r_D,g_D)
are consequently exactly

    [(1,6)],
    [(2,3)],
    [(1,4),(1,2)],
    [(2,2),(1,2)],
    [(1,3),(1,3)],
    [(1,2),(1,2),(1,2)].                                 (11c)

Thus there are at most three closed intermediate infinity primes.
Distinct entries may have the same target image only when their
combined r is at most two, by (11a). In particular the two entries
in [(2,2),(1,2)] have different images; the last case has at least
two images. This is a field-degree packet classification, not
a claim that the original quadratic cover splits over every
non-algebraically-closed residue field.

## 6. Quantitative constraints on the entire quartic-face ledger

Insert D>=m+4S from (2) into (10), and combine with b<=2S-5n_A:

    (2S+2m)/3<=b<=2S-5n_A,
    4S>=15n_A+2m,
    S>=ceil((15n_A+2m)/4).                                (12)

Because n_A>=1 and m>=2, equation (12) gives S>=5. Hence

    D>=m+4S>=22,
    b>=ceil((D+3m)/6)>=5,
    [L:F]=2b>=10.                                       (13)

These are GLOBAL bounds on the remaining quadratic ruling-six
class. They do not identify ruling degree with generic map degree.

For example n_A=1,m=2 forces L_1>=5. Its least numerical
ledger is S=5,D>=22,b>=5. For n_A=2,m=2, S>=9,D>=38,b>=8.
These are compatible integer lower limits, not realizations
of a Keller map.

Suppose some infinity IMAGE has total r_D equal to two.
This includes one prime with rho=2 or two different primes
with rho=1 over the same image. The joint bound (7a) and
the ceiling-product estimate give

    b>=ceil((D+3m-6)/3)+1.

Combining again with (2) gives

    2S>=15n_A+4m-3.                                       (14)

Since n_A>=1,m>=2, this requires S>=10, and therefore

    D>=42, b>=15, [L:F]>=30.                             (15)

For the b bound use D>=42,m>=2 in the preceding strict
rho-two estimate, giving b>=ceil(42/3)+1=15.

Thus if D<42 OR b<15, every infinity prime has

    alpha=f=1.                                         (16)

Their target images are also PAIRWISE DISTINCT, since two primes
over the same image would have combined r equal to two.

This is unramified projection plus degree one OVER ITS IMAGE.
It does not make the image curve's P-degree one, its finite
image a graph, or the map invertible.

## 7. Scope and verification

The strict source-packet lemma uses etaleness, quasi-finiteness,
normality of the full source, a source polynomial u in E, and
the complete finite intermediate normalization. The local degree
consequences additionally use the one-dimensional plane generic
P fibre and the exact ruling-six ledger. No dimension-free
derivative-integrality implication is inferred.

The checker verifies the local coframe/index equality, exact
complete-packet degree accounting, the ceiling products, all infinity index/residue pairs, the COMPLETE nonlinear image
packet list of total degree six, and the symbolic ledger inequalities with complete residue
class controls for the stated integer minima. The checker passes 141 exact controls. The written proof,
not bounded sampling, establishes the unrestricted assertions.

Replay with

    .venv/bin/python research/jacobian/affine-complements/verify_ruling_six_strict_infinity.py

All original primes remain in scope. No new external classification
is imported; independent whole-proof review, formal verification
and external human review remain unclaimed. Full ruling-six closure,
automatic trace perfection, normal-tail termination and JC2 remain open.
