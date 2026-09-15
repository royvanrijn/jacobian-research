# Good branch reduction separates branch values from trace norms

On direct11952 Q80, suppose the reduced branch divisor at131 is squarefree,
has degree two or four, and avoids the reduced parent discriminant. The
space of branch-value patterns compatible with the seventeen inherited
reciprocity characters has dimension **zero in degree two and at most one
in degree four**. This allows zero-section values at branches and imposes
no bound on section heights or coefficient denominators.

Consequently, for an actual quadratic cover in this scope, write r for
its rank increment and n for the dimension of the inherited trace parities
realized over that cover. Then

```
genus0: r=n,
genus1: r=n+epsilon,  epsilon in {0,1}.                 (1)
```

The parent has arithmetic rank17. Thus a rank-two increment on a genus-one
cover with this reduction must realize at least one nonzero inherited
trace parity; on a genus-zero cover it must realize at least two independent
trace parities. These are function-field trace norms, not numerical heights
or specialized Selmer proxies. Equations (1) do not bound r by one.

The [positive correlated-gain objective](CORRELATED_QUADRATIC_GAINS_2026-09-12.md)
remains open. No new MW17 section or soluble norm torsor is constructed.
The general norm criterion below is standard descent made explicit for
this objective; its new application is the complete good-branch code bound.

## 1. The branch character space

Let K=Q(t), L=K(w), w^2=d, and let M=E(K), G=E(L), M^- the anti-invariant
subgroup. Assume d has geometric branch degree b=2 or4, with branches only
at smooth parent fibres. On the degree-two pullback chi=4 and all fibres
are irreducible, so every nonzero section has height at least8 and G has
no torsion. The genus of the covering base is (b-2)/2.

Normalize the binary branch form to be primitive over Z_131. Our additional
reduction hypothesis is that this form has squarefree reduction disjoint
from the reduced degree24 discriminant. This is a condition on the branch
divisor in the fixed Q80 model. The literal constant squareclass of d is
arbitrary; in particular, we do not assume the covering curve itself has
good reduction when its scalar has odd131-valuation.

The branch algebra is unramified over Q_131. At a closed residual branch
point v of degree f, let V_v=E_v(F_(131^f))[2]. It is an F2-space of
dimension zero, one or two. Its nonzero elements are the cubic roots over
that field; O is its zero element. Define

```
C_v : V_v -> F2^17,
C_v(O)=0,
C_v(e)_j = squareclass of Norm(X_(T_j)-e).              (2)
```

Use the derivative 3e^2+A at a coincident2-torsion value and the square
class at a pole, as in the retained
[regularized reciprocity proof](Q80_DEGREE_TWO_RECIPROCITY_2026-09-14.md#1-coupled-norm-conditions-including-zero-ordinates-and-poles).
The product of the three root characters is a square, so C_v is linear.
Every P in M^- specializes into the2-torsion at each branch. Its actual
specializations inject into the residual spaces by good reduction of
prime-to-131 torsion. Reciprocity implies

```
sp(P) in B(d_bar) := ker( direct_sum_v V_v --sum C_v--> F2^17 ).   (3)
```

This uses actual branch values, including O. It needs no assumption on
Gauss valuations of coordinate polynomials. The local spectral argument
and norm specialization are inherited from the cited proof.

## 2. Complete finite bounds, including zero branch values

The [new result](../artifacts/generated-results/elkies-k3-q80-good-branch-code-v1/result.json)
reuses the complete rational/quadratic roster and completed cubic/quartic
certificates. The new joins below operate on these records; no large
census is repeated.

| Residual branch degrees | Bound on dim B | Reason |
|---|---:|---|
|1+1|0|All110 rational nonzero root codes are distinct|
|2|0|Every smooth quadratic-site map is injective|
|1+1+1+1|1|Complete short-word check, allowing O at any site|
|1+1+2|1|No quadratic code plane occurs in a rational-pair span|
|2+2|1|Distinct smooth quadratic sites have distinct code planes|
|1+3|1|Both retained cubic gates, with the rank-one case included|
|4|1|The retained quartic census excludes a zero code plane|

Here is the precise linear algebra and its certificate boundary.

For rational sites, the110 codes are nonzero and globally distinct. For
each of the16 full-splitting sites its three codes form one plane. Hence
there are no words supported at one or two distinct sites. A direct
complete enumeration finds four words on three sites and forty on four
sites. No two different words have support union of size at most four.
Two independent codewords would be two such words, proving the bound.
Root-free sites may be included with the zero value; the support-union
test accounts for them.

For the1+1+2 case, enumerate all8,646 pairs of distinct rational sites,
including root-free ones. Their maps are injective, and their spans contain
11,786 distinct planes. None equals any of the1,456 smooth full-splitting
quadratic code planes. A quadratic site with only one root cannot give a
kernel of dimension two. The quadratic roster has no zero norm code at a
smooth site; its only zero is at the omitted nodal quadratic. Distinct
quadratic code planes are also unequal, proving2+2.

For1+3, if the cubic-site space has dimension at most one, the bound is
automatic. If it has dimension two and its code map is injective, a kernel
of dimension two would require equality with a full rational code plane.
The retained k0 cubic-signature census excludes this. If its code map has
rank one, a kernel of dimension two would require its nonzero image to
be a rational root code. The retained k1 census lists all seven smooth
cubic sites with a zero norm code; each has exactly one zero and its other
code is absent from all110 rational codes. No cubic code plane is zero.
The unique nodal cubic is outside the present smooth-residue hypothesis.

For4, a kernel of dimension two requires all three root norm codes to be
zero on a full-splitting quartic site. This is precisely what the completed
73,620,690-orbit producer and independent replay exclude. Their ten
gap-free checkpoint ranges are checked against the retained hashes and
against each other. Their arithmetic is inherited, not rerun here.

The four rational three-site patterns, written (base,root), are

```
[(2,54), (76,3), (103,26)],
[(27,88), (50,79), (103,26)],
[(65,11), (70,54), (84,100)],
[(98,64), (112,30), (115,68)].                           (4)
```

They are necessary residue data, not points or global sections.

## 3. Why the kernel of branch specialization is an anti-trace

The arithmetic input is the Q80 unramified Kummer theorem proved in
[the branch-stratum note, section3](COMMON_QUARTIC_BRANCH_STRATA_2026-09-14.md#3-the-unramified-arithmetic-input-on-q80):

```
every norm-square spectral class with all valuations even is delta(M). (5)
```

This includes bad fibres and infinity. It uses the retained Brauer and
Picard results on the determinant948 surface and their transfer to the
literal Q80 fibration. We inherit that written sheaf/descent argument.

Let N=(1-sigma)G, a subgroup of M^-. Every anti-trace specializes to O at
every branch, so N is contained in ker(sp). Conversely, if P specializes
to O at all branches, its twist Kummer representative d*(x(P)-theta) has
even valuations everywhere. By (5) it equals delta(T) for some T in M.
Restricting to L and using Kummer exactness gives an actual R in G with

```
2R=P+T.
```

Conjugating and using absence of2-torsion gives Tr(R)=T and
(1-sigma)R=P. This proves

```
ker(sp)=N,
H^1(<sigma>,G)=M^-/N = im(sp) subset B(d_bar).           (6)
```

No half is introduced formally. Its existence is supplied by (5) and
the Kummer map. The identification (6) is over Q; the code bound was
obtained by specializing the branch values over Q_131.

## 4. Rank accounting and an explicit sufficient norm criterion

Set W=Tr(G)/2M, a subspace of M/2M, and n=dim W. The map

```
N/2M^- -> W,     (1-sigma)R |-> Tr(R) mod2M            (7)
```

is an isomorphism. Changing R by an invariant point changes its trace
by twice that point. Its kernel consists exactly of the doubles: if
Tr(R)=2T, then R-T is anti-invariant and (1-sigma)R=2(R-T).
Conversely a doubled anti-invariant has even trace. Surjectivity is by
definition. Since M^- is free of rank r,

```
r = dim(N/2M^-) + dim(M^-/N) = n + dim im(sp).          (8)
```

Combining (8) with section2 proves (1). For a rank-two increment in the
good quartic scope, at least one primitive anti-invariant combination
is an anti-trace with nonzero inherited trace parity. If its trace parity
were zero, that combination would be twice another anti-invariant.
The exact norm subgroup, not the branch code alone, is the missing part.

For explicit construction take T=(a,b_T) in M and put

```
H_T(m)=m^4-6*a*m^2-8*b_T*m-3*a^2-4*A.
```

The trace-norm torsor for T has the affine equation

```
H_T(m)=d*s^2,    m,s in Q(t).                          (9)
```

Given a solution with s nonzero, the point

```
x(Q)=(m^2-a+s*w)/2,
y(Q)=m*(x(Q)-a)-b_T                                   (10)
```

satisfies the parent equation and Tr(Q)=T. This is the regular chord
identity without any degree restriction on m. If T has nonzero parity,
s=0 would give a rational half of T and is impossible; the two points
at infinity of the quartic require sqrt(d), so add no missing K-point
for a genuine quadratic extension. The torsor is the fibre of the norm
map Res_(L/K)E -> E over T, with Jacobian E^d.

Thus a concrete sufficient condition on the parent and d is:
**two independent classes T1,T2 in M/2M have soluble torsors (9) for the
same d.** Then n>=2 and (8) gives two independent new directions, without
having to establish independence from numerical heights. The elementary
inequality r>=n and this sufficient condition do not require the special
Q80 unramified theorem; that input is needed for the branch description
and exact small remainder in (1).

This condition is not yet verified on a low-genus cover of an MW17 parent.
Solubility means actual rational functions m,s, not local points or Selmer
classes. The base still needs a rational parametrization or a pointed
genus-one model with a certified nontorsion point to provide infinitely
many rational specializations. Equation (9) identifies the remaining
global lifting problem rather than asserting it has been solved.

## 5. Consequences for the first rational-pole layer

Use the [pole formula](CORRELATED_GAIN_POLE_CRITERION_2026-09-14.md):
height8 means no O-valued branch, and height10 means exactly one O-valued
branch. A height8 and a height10 section have different nonzero branch
patterns, hence two independent patterns over F2. Section2 excludes such
a pair in the good quartic scope, without coefficient bounds.

Two height10 sections can survive only with identical branch values,
including the same rational branch pole. If independent, their cross
pairing must be0 or+/-2. Indeed P+Q and P-Q have O at all four branches,
so each has height16+4m for an integer m>=0, and their heights sum to40.
For cross pairing0 the two anti-traces both have height20; their half-
descended bisection images have arithmetic genus2 and normalization genus1.
For cross pairing+/-2 the heights are16 and24, giving arithmetic genera1
and3. These use the retained image formula p_a(B)=h(U)/4-3, where
U=R-sigma(R). They isolate the remaining singular-image problem; they do
not exclude it. With four rational residue branches the nonzero values
must be one of (4).

In genus zero, good branch reduction forces every anti-invariant section
to take O at both branches, hence heights12+4m. This extends a necessary
height restriction beyond the retained rational-bisection atlas, without
claiming completeness of that atlas at higher heights.

## Evidence, replay and boundaries

The [producer](scripts/audit_q80_good_branch_code.py) finishes the new
finite joins and retained checkpoint checks in about0.022 seconds. The
[independent checker](scripts/verify_q80_good_branch_code.py) enumerates
triples directly and lists all planes inside rational-pair spans, rather
than using the producer's pair-bucket joins. Its
[receipt](../artifacts/generated-results/elkies-k3-q80-good-branch-code-v1/independent-replay.json)
records0.055 seconds, the44 words,11,786 rational-pair planes and1,456
quadratic planes. A synthetic rank-two configuration is detected, so the
negative result is not a hardcoded dimension cap. Both use at most20 CPU
seconds and1GiB; the imported rational audit tightens its producer to15
CPU seconds. Higher-degree arithmetic and global unramified descent are
inherited dependencies. No formal or external verification is claimed.

```
python3 research/elkies-k3/scripts/verify_q80_good_branch_code.py
```

The first [rational-only audit](../artifacts/generated-results/elkies-k3-q80-rational-branch-code-v1/audit.json)
is retained unchanged. Discovery scripts refuse to overwrite result files.
Repeated branch reductions, branches reducing to nodal fibres, other primes,
other MW17 parents, and actual solubility of the norm torsors remain outside
this theorem. The old genus-zero norm-eight replay gap is not used or
upgraded. The positive rank-at-least19 objective remains open.
