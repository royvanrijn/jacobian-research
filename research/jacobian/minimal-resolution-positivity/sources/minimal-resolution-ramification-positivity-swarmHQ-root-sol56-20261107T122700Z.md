# Positivity of ramification on the actual minimal resolution

Claim ID: MINIMAL-RESOLUTION-RAMIFICATION-POSITIVITY-1. Producer: swarmHQ ROOT. Independent reviewer: gpt-5.6-sol / ultra configuration; actual hosted identity UNVERIFIED. Evidence tier **MANUAL / PROVISIONAL / CONFIRMED_AT_NAMED_FOUNDATIONS**. Novelty unassessed. JC2 remains unresolved.

The full ramification divisor and its canonical adjoint are nef and big on the actual minimal smooth model of any hypothetical degree>1 Keller map. Their zero intersections are confined to contracted smooth rational (-2)-curves. This strengthens the retained arithmetic-genus obstruction and supplies inequalities against every curve on that specific model.

## Reviewed statements

Let f:Y -> P^2 be a generically finite morphism of degree N>1 from a smooth irreducible projective surface. Let U be an open subset isomorphic to the whole A^2, with f etale on U and f(U) missing a target line L. Assume f contracts no smooth rational (-1)-curve. Write H=f^*L, K=K_Y, R=R_f for the full effective determinant divisor, including multiplicities on contracted components, and Q=2K+3H=K+R in divisor class. The reviewed conclusions are:

1. Q and R are nef and big.
2. For an irreducible curve C, Q.C=0 if and only if C is an f-contracted smooth rational (-2)-curve. The same equivalence holds with R.C=0.
3. p_a(R)>=2 and h^0(Y,Q)>=2.

These statements are intended for the SAME minimal resolution of the finite projective normalization of a hypothetical Keller map. They do not assume that normalization smooth, Q-Gorenstein or rational-singular, or that f itself is finite. Arbitrary extra boundary blowups are excluded. R is not assumed reduced, connected or SNC.

The retained MINIMAL-RESOLUTION-RAMIFICATION-GENUS-1 interface gives p_a(R)>=1 under exactly these hypotheses, at its MANUAL/PROVISIONAL, independently confirmed named-foundation scope. Its [original report](minimal-resolution-ramification-genus-swarmHQ-root-sol56-20261107T0755Z.md), full SHA2564d38a1b2b1c72e3dc9c32f38a3ad036aba9786af594438afd8e258a9546213df. No fresh proof audit or promotion of that interface is claimed.

Named foundations additionally used below: smooth-surface cone/contraction theorem; rational-surface classification and K^2+rho=10; adjunction/Riemann--Roch and Serre duality; Hodge index; Bertini and Stein factorization; divisor-class localization; proper quasi-finite implies finite; and absence of a nontrivial connected finite etale cover of complex A^2. These are named standard inputs, not newly source-proof audited. No semiampleness or vanishing theorem is required.

## 1. Nefness of Q

Y is rational, H is nef and H^2=N>0, while R is effective and R~K+3H. Suppose Q is not nef. On the K-nonnegative part of the curve cone, Q=2K+3H is nonnegative. A K-negative extremal ray therefore has Q negative, exactly as in the retained genus proof.

A birational extremal contraction has a smooth rational (-1)-curve C. If H.C=0, f contracts C, contrary to the hypothesis. Otherwise H.C is a positive integer and Q.C=-2+3H.C>=1, again a contradiction. A contraction to a point gives Y=P^2 and H=m[line]; Q negativity forces m=1, hence N=1.

For a contraction to a curve, a general fibre F has K.F=-2 and H.F>0 by Hodge index. The inequality Q.F=-4+3H.F<0 forces H.F=1. The one-ray and connected-fibre argument in the retained proof then makes every fibre irreducible and reduced and gives a P^1-bundle Y=F_n. Explicitly, a connected reducible fibre would have a component of negative square although every component is proportional to the square-zero fibre class; a multiple fibre mC is excluded by H.F=1. Arithmetic genus zero gives smooth P^1 fibres. Write H=S+bF. Then

R ~ S+(3b-n-2)F,

whose arithmetic genus is zero by adjunction on F_n, independently of effectivity, reducibility or multiplicities of a representative. This contradicts the supplied p_a(R)>=1. Thus Q is nef. The reasoning does not merely deduce nefness from the sign of Q.R: it excludes every possible Q-negative extremal contraction.

Since R=(Q+3H)/2 numerically, R is nef. Its square is (Q^2+6Q.H+9N)/4>0, so R is big.

## 2. Q is not numerically zero

If Q is numerically zero, K is numerically -(3/2)H. Thus K^2=9N/4 is an integer. Rationality gives K^2=10-rho(Y)<=9, so 4 divides N and N>1 forces N=4, K^2=9 and rho(Y)=1. A smooth rational projective surface of Picard rank one is P^2. Hence H has class2[line].

The boundary-prime group is isomorphic to Pic(Y): divisor-class localization and O(U)^*=O(Y)^*=C^*, Cl(U)=0 give this directly. There is exactly one boundary curve C, and its effective class must generate Pic(P^2), so C is a line. Both R and f^*L are supported on that line: f is etale on U and f(U) misses L. Consequently R=3C and f^*L=2C as effective divisors.

H is ample, so f has no positive-dimensional fibre and is finite. Its restriction over P^2 minus L has source P^2 minus C=A^2, and is etale because the full determinant divisor is supported on C. This would be a connected finite etale degree4 cover of complex A^2, impossible. This argument does not require excluding isolated boundary points separately. Therefore Q is not numerically zero.

## 3. Excluding Q^2=0 without a semiampleness import

Suppose Q^2=0. Since Q is nonzero nef and H is nef of positive square, Q.H>0 by Hodge index. From K=(Q-3H)/2 we obtain K.Q=-(3/2)H.Q<0. For large n, H^2(Y,O(nQ))=H^0(Y,O(K-nQ))=0 because (K-nQ).H<0. Riemann--Roch gives

h^0(Y,nQ) >= chi(O(nQ)) = 1-(n/2)K.Q,

which tends to infinity. Choose n with h^0(nQ)>=2 and remove the fixed part of |nQ|, leaving a nonzero movable divisor class M. Nefness and nQ.M+nQ.Z=n^2Q^2=0, where Z is the fixed part, give Q.M=0. Two general members of |M| have no common component, so M^2>=0. Hodge index for the isotropic class Q implies M^2=0 and M is a positive numerical multiple of Q.

Because M^2=0, two general members cannot meet, hence |M| has no base points. Its image has dimension one. After Stein factorization its general connected fibre F is smooth, has F^2=0 and is numerically a positive multiple of Q. Thus K.F<0, so adjunction forces F=P^1 and K.F=-2. But Q.F=0 then says

0 = 2K.F+3H.F = -4+3H.F,

contrary to the integral intersection number H.F. This excludes Q^2=0. Therefore Q is nef and big.

## 4. Null curves, arithmetic genus and section count

Let C be irreducible and Q.C=0. Hodge index gives C^2<0. The identity 2K.C+3H.C=0 gives K.C<=0. Adjunction gives K.C=2p_a(C)-2-C^2>=-1. If H.C>0, integrality in that identity instead gives K.C<=-2, a contradiction. Hence H.C=K.C=0. Adjunction then gives p_a(C)=0 and C^2=-2, so C is a smooth rational f-contracted curve. The converse follows immediately from H.C=K.C=0. As R=(Q+3H)/2 and both summands are nef, the same null curves characterize R.

R is nef and big and Q is nef and nonzero, so Q.R>0 by Hodge index. Cartier adjunction yields Q.R=2p_a(R)-2. Consequently p_a(R)>=2. Finally H^0(K)=H^1(K)=0 on rational Y, so 0 -> O(K) -> O(K+R) -> omega_R -> 0 identifies H^0(Q) with H^0(omega_R). Duality on the Cartier divisor gives h^0(omega_R)=p_a(R)+h^0(O_R)-1>=2, without assuming R reduced or connected.

## 5. Numerical identity, minimality and remaining gap

Let gamma be the genus of a smooth general member of |H|. Adjunction and K^2+rho(Y)=10 give the exact identity

p_a(R)=1+K^2+(9/2)(H^2+K.H)=2+9gamma-rho(Y).

The bound implies rho(Y)<=9gamma. More substantively, Q^2=16+24gamma-4rho(Y)-3N>0, and the nef inequalities hold for EVERY curve on the actual minimal model. No low-degree census or degree-bound re-audit is proposed.

For an extra boundary blowup beta, H'=beta^*H, R'=beta^*R+E and Q'=beta^*Q+2E. Thus R'.E=-1 and Q'.E=-2: nefness fails exactly as expected, and p_a(R')=p_a(R)-1. A count on an arbitrary common resolution therefore cannot be used as the opposite inequality on this minimal model.

No contrary intersection, empty adjoint space, properness, surjectivity, finite affine Keller realization or JC2 proof is obtained. The two-section conclusion does not repair the existing trace, finite-pole or saturation gaps. The remaining question is whether the whole affine-plane source forces a curve or another global condition incompatible with these nef/big and null-curve constraints on the same actual minimal model.

## Independent review and integration

The [independent review](actual-ramification-nef-big-review-sol56-20261107T1205Z.md), full SHA25683d4c42e463a1a09913e60f0154059746e6323434c73c0cc26cf079d38c59763, confirms the statement and every listed endpoint at the named foundations. ROOT read the complete159-LF sealed review, checked its mathematical reconstruction, and verified the input, command and reading records. The frozen producer required no mathematical repair. The exposition above retains its proof, with lifecycle wording updated after review. The [evidence bundle](../box/minimal-resolution-ramification-positivity-swarmHQ-20261107T1220Z/README.md) preserves exact producer bytes and compact review custody.

The review covers all Q-negative extremal contractions, the numerically-zero Q endpoint, the elementary movable-pencil argument for Q^2=0, the null-curve adjunction calculation, the two-section consequence, and sensitivity to extra boundary blowups. It explicitly checks the actual-normalization application: for finite projective normalization h:X->P^2 and its minimal resolution q:Y->X, q is an isomorphism over the smooth source open, and the curves contracted by g=hq are q-exceptional. No smoothness, rational-singularity or Q-Gorenstein assumption on X is added.

All18 charged packet files including the manifest and all tools/prompt were unchanged. The42 verified windows cover1836 non-state LF and98 state LF, including all175 historical-file LF. A declared171-LF subtotal underestimated the producer seal by four lines; the reviewer identified and charged all175 within the unchanged240-LF cap. No input or sealed report was repaired. These execution checks support custody and scope; they are not a formal proof certificate.

The statement depends on MINIMAL-RESOLUTION-RAMIFICATION-GENUS-1 at its existing named-foundation tier, and thereby retains that interface's smooth-finite dependency. It is a necessary condition, not an opposing intersection/genus theorem or a proof of JC2.
## OPEN(S) RAISED

None. The existing whole-source contradiction gap remains open; no new OPEN identifier is introduced.

## COLLISIONS

status: EMPTY

- NONE — the report contains no explicitly raised `OPEN[...]` entries.

Corpus basis: `27c4ab29317fa30501ddcecb9ba25b0ae6bca09c` (Git blobs only).

<!-- BODY-END -->

## Seal

- Body definition: every byte through the unique standalone `<!-- BODY-END -->` line,
  including its terminating newline; this seal is outside the body.
- Body bytes: `10635`.
- Body SHA-256:
  `bf55ceb608090e056759353e6714f4493e44b4fc823546fb56e0cab1d686d027`.
- Frozen basis: `27c4ab29317fa30501ddcecb9ba25b0ae6bca09c`.
