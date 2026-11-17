# Centralizer rigidity and a normalized fourth quantum image

14 November 2026. **QCR1.** In every characteristic-zero polynomial
symplectic frame, fixing all images except one conjugate leaves only
a univariate polynomial ambiguity in the paired image. The same
statement holds h-adically for formal quantum frames, with a completed
univariate algebra rather than a bounded polynomial algebra.

For Long's retained triple (R,S,T), every classical completion is
D_Long+P(R). The displayed D_Long already has the canonical
normalization defined below, and every such completion has total
degree at least 23 and momentum order exactly five. This closes
searches for another lower-degree Hamiltonian correction while
keeping that triple fixed.

When the first quantum image is the exact positional polynomial
R=x(2-3xq), the remaining formal ambiguity has an explicit projection
that never increases spatial degree. If a completing formal D exists,
a bounded one exists for that fixed quantum triple if and only if
its normalized D is bounded. This does not prove existence of D
for an arbitrary quantized triple, or boundedness for Long.

## A centralizer theorem for every full polynomial symplectic frame

Let A=K[z_1,...,z_(2n)] in characteristic zero, and suppose
\[
 \{d,r\}=1,\qquad
 (r,d,s_1,t_1,\ldots,s_{n-1},t_{n-1})
\]
is a full polynomial symplectic frame. Define
\[
 C=\{g\in A:\{g,r\}=\{g,s_i\}=\{g,t_i\}=0
                                    \text{ for all }i\}.
\]
Then
\[
                                  C=K[r].              \tag{1}
\]

Indeed the symplectic Jacobian identity makes the differentials
of the full frame a polynomial basis of Omega_A^1. Expanding dg
in this basis, the bracket constraints force
\[
                                   dg=b\,dr
\]
for some b in A. In characteristic zero, dg wedge dr=0 implies
that g and r are algebraically dependent.

We use the classical polynomial dependence theorem:
algebraically dependent nonconstant polynomials lie in one
polynomial algebra K[u]. This is
[Arzhantsev--Petravchuk, Lemma 2](https://publications.hse.ru/pubs/share/direct/219130966.pdf),
following their generative-polynomial lemma. Write r=P(u), g=Q(u).
The polynomial Hamiltonian mate d gives
\[
                       1=\{d,r\}=P'(u)\{d,u\}.
\]
Units of A are scalars, so P'(u) is a nonzero scalar. Since u is
nonconstant and the characteristic is zero, P has degree one.
Thus g belongs to K[r]. Constants cause no exception.

Consequently all polynomial choices d' completing the same
remaining 2n-1 images are exactly
\[
                                  d'=d+P(r),\qquad P\in K[u].
 \tag{2}
\]
No inverse polynomial map is assumed. The theorem applies to
noninvertible symplectic endomorphisms just as to automorphisms.
It assumes an admitted classical frame; it does not establish
admission from an arbitrary partial tuple.

## Formal quantum centralizers have exactly the completed ambiguity

Let the fixed formal images
\[
                  \mathsf R,\mathsf S_i,\mathsf T_i\in A[[h]]
\]
have the above classical reductions and satisfy their exact mutual
Weyl relations. In particular \(\mathsf R\) commutes with every
\(\mathsf S_i\) and \(\mathsf T_i\). Its formal joint centralizer is
\[
 \left\{\sum_{\nu\ge0}h^\nu P_\nu^\star(\mathsf R):
                                      P_\nu\in K[u]\right\}.     \tag{3}
\]
The notation means star-polynomial evaluation. There is no common
bound on the degrees of the P_nu in (3).

To prove this, let C commute with all the fixed images. The
coefficient of h in each commutator shows that its reduction C_0
lies in the classical centralizer (1), so C_0=P_0(r). Subtract
\(P_0^\star(\mathsf R)\), which commutes with all the fixed images, and
divide the remaining series by h. Repeat. The resulting sum
converges h-adically, and subtraction is unique: at the first
differing order, a nonzero polynomial in the nonconstant r cannot
vanish. Conversely every series in (3) commutes with the tuple.

If two formal conjugate images \(\mathsf D\) and \(\mathsf D'\) exist for
this fixed partial tuple, their difference is exactly an element
of (3). For a prescribed common classical image, its h^0 term
must vanish. This classifies the ambiguity, not existence of a
conjugate for an arbitrarily prescribed quantum tuple.

## Apply the theorem to the corrected Long triple

For
\[
 R=x(2-3xq)
\]
and the retained S,T of Long's full map, (1)--(2) give
\[
 \operatorname{Cent}_{\{,\}}(R,S,T)=K[R],\qquad
                       D'=D_{\rm Long}+P(R).           \tag{4}
\]
The [retained exact transcription](../dixmier-audit/verify_long_pbw.py)
has
\[
 \deg D_{\rm Long}=23,\qquad
 \operatorname{ord}_{p,z}D_{\rm Long}=5,\qquad \deg R=3.
\]
Here momentum order is the maximum combined exponent of p and z
in the original polarization. Exact nonzero witnesses in D_Long are
\[
                 [x^{12}q^6z^5]D_{\rm Long}=144/5,
             \qquad [x^6z^5]D_{\rm Long}=32/405.         \tag{5}
\]
The full support bounds are checked, not inferred from these two
monomials alone.

Since P(R) has momentum order zero, every completion in (4) has
momentum order exactly five. For nonconstant P of degree k,
deg P(R)=3k. No multiple of three equals 23. Thus
\[
                     \deg(D_{\rm Long}+P(R))
                          =\max(23,3k)\ge23.            \tag{6}
\]
Constants leave the degree unchanged. There is no lower-degree
or lower-momentum-order fourth image for the retained triple,
irrespective of the degree allowed in the search.

This is a fixed-triple result. Altering S,T or the source
polarization changes the problem. It does not establish an
intrinsic minimum over Long's full canonical orbit.

## An explicit projection removes the remaining R freedom

For a normal PBW polynomial F in (x,q,p,z), put
\[
 {\cal N}(F)=F-\sum_{k\ge0}
       \frac{[x^{2k}q^k p^0z^0]F}{(-3)^k}\,R^k.        \tag{7}
\]
The sum is finite for every polynomial F.

Indeed
\[
 R^k=\sum_{j=0}^k
       {k\choose j}2^{k-j}(-3)^j x^{k+j}q^j
\]
has exactly one monomial of the form x^(2l)q^l, namely l=k,
with coefficient (-3)^k. Thus (7) is a linear idempotent with
kernel K[R]. Its image consists of polynomials with all these
pivot coefficients zero.

Every subtracted R^k has degree 3k, equal to the degree of the
pivot that caused its subtraction. Therefore
\[
                            \deg{\cal N}(F)\le\deg F.  \tag{8}
\]
It also leaves every positive-momentum-order term unchanged.
The exact retained D_Long has no pivot terms, so
\[
                             {\cal N}(D_{\rm Long})=D_{\rm Long}.
\]

Extend \({\cal N}\) coefficientwise to A[[h]]. Since R depends only
on positions, its star powers are its ordinary powers. If the
fixed quantum tuple has \(\mathsf R=R\) exactly, equation (3) is
K[R][[h]], and
\[
        {\cal N}(\mathsf D)={\cal N}(\mathsf D')
\]
for any two completing conjugates. Subtracting the removed
central series preserves every required Weyl relation.

Moreover the existence of a spatially bounded completion of
that fixed quantum tuple is equivalent to boundedness of this
single normalized representative:
\[
 \begin{split}
 &\text{some completing }\mathsf D'\text{ has a common spatial bound}\\
 &\hspace{30pt}\Longleftrightarrow\
           {\cal N}(\mathsf D)\text{ has a common spatial bound}.
 \end{split}                                           \tag{9}
\]
One direction uses (8); the other uses the fact that the normalized
series is itself a completion. In particular no larger search
through the P_nu(R) ambiguity can improve its coefficient degrees.

Within the finite PBW algebra with scalar coefficients K[[h]],
the centralizer of the fixed tuple is exactly K[[h]][R].
If an element has spatial bound D, its expression in R has degree
at most floor(D/3), since each coefficient has that bound and
powers of R have distinct total degrees.

## The exact R slice is formally reachable

For Long, keeping R exact is possible at the formal level. Begin
with any full formal lift supplied by [QFL1](PROOF.md). In the
polynomial coordinates (x,q,B,D_0) from [QCC1](CUBIC_CORE_NONCOMPLETION.md),
the classical operation {-,R} is partial_(D_0), hence is surjective
on the polynomial ring by integration in D_0.

At the first order where the lifted first image differs from R,
choose a polynomial Hamiltonian whose bracket with R cancels that
coefficient. Its h-adically small inner gauge preserves all Weyl
relations and removes that term. Repeat. Every coefficient remains
a polynomial, so this produces a full formal lift with first image
R exactly. Applying (7) then normalizes its completing fourth image
without changing the classical Long map, since D_Long is already
normalized.

This formal gauge can introduce growing spatial degrees. It is not
a proved bounded algebraic normalization of every bounded lift.
Consequently a failure for all R-fixed triples, were it established,
would still need a bounded-normalization comparison to exclude all
possible lifts. No such failure or comparison is proved here.

## Why the undeformed R hypothesis matters

For a general deformed first image, formal centralizer closure
does not become a finite polynomial algebra merely because an
individual centralizer element is spatially bounded.

For example, in normal order with [p,x]=h, set
\[
                  X=x+h^2x^2,\qquad
                  P=(1+2h^2x)^{-1}p.
\]
Then [P,X]=h exactly in the h-adic algebra, with classical
identity reduction. The element x has spatial degree one and
commutes with X, but
\[
                              x\notin K[[h]][X].
\]
Over K((h)), any nonconstant polynomial in X has even positive
x degree. Nevertheless x belongs to the completed algebra:
\[
       x=\sum_{k\ge0}(-1)^k C_k h^{2k}X^{k+1},
       \qquad C_k=\frac{1}{k+1}{2k\choose k}.           \tag{10}
\]
Thus (3) genuinely needs h-adic completion. In fact this fixed X
has no spatially bounded polynomial Weyl mate at all. Over K((h)),
a finite normal PBW operator Q of p order m>0 has leading
commutator coefficient h*m*q_m(x)*(1+2h^2*x). A scalar nonzero
commutator forces m=1 and q_1(x)*(1+2h^2*x)=1, impossible for a
polynomial coefficient. The same p-order argument works with any
number of spectator pairs as coefficients.

So even a spatially bounded first image with identity classical
reduction can force its completing mate to have unbounded spatial
degrees. The identity formal gauge class still has the bounded
representative (x,p). This is failure of the specified first-image
choice, not intrinsic nonliftability of its formal class.

## The remaining constructive question

The search for another polynomial H with the same classical
R,S,T is closed by (4). For a fixed quantized triple with R
undeformed, the only remaining D ambiguity is removed by (7).
What remains is to construct bounded S,T and a completing
normalized D, or to prove an obstruction covering all such
choices. The normalization makes D-growth meaningful within
one fixed triple; it does not make it invariant under arbitrary
formal gauges or changes of the other images.

The [cubic-core obstruction](CUBIC_CORE_NONCOMPLETION.md) explains
why replacing corrected D by D_0 is not such a simplification.
The [scalar-anomaly criterion](SCALAR_ANOMALY_AND_CUBIC_LIFTS.md)
and [bounded realization theorem](BOUNDED_REALIZATION.md) remain
positive construction gates. Neither is proved for Long by
this centralizer classification.

[Exact controls](check_centralizer_normal_form.py) verify the
retained support and degree witnesses, projection identities,
polynomial R ambiguities of several degrees, the exact quantum
inverse-coordinate control and its finite Catalan truncation.
The general centralizer theorem uses the stated published
polynomial-dependence result. Formal centralizer exhaustion,
normalization optimality and all-degree conclusions use written
proofs. No independent theorem replay, formal verification,
external review, novelty priority, bounded Long lift or DC2
solution is claimed.
