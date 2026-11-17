# Complete individual phase classification in the known cotangent algebra

16 November 2026. Internal written deductions from the retained normalized prime-separating theorem PSMT1 and the explicit joint inverse. Independent replay, external review and formal verification remain pending. No new map, smallest Gaussian endpoint or novelty is claimed.

This closes the previously open classification of EVERY individual zero-constant H(S,T) in the explicit pair S=Y1,T=-Y2. The [actual joint image](ACTUAL_JOINT_IMAGE_AND_FACTORIAL_RECOMPOSITION.md) exposed a small algebraic period. Its quadratic source marker turns out to carry a one-radial, prime-separating functional. This extra arithmetic input, rather than generic rank or two scalar source generators alone, decides the full original source-marker question.

## 1. A general prime-separating marker certificate

Let a Gaussian Wick functional have an integer charge grading that annihilates every nonzero charge, with every source variable of charge of absolute value at most w. Suppose two polynomial phases S,T have charges-1,+1, every positive joint moment vanishes, and a fixed marker Q0 satisfies

    L[Q0 S^i T^j]=delta_ij*a_i for all i,j>=0.               (1)

Assume every a_i is nonzero rational and the sequence is normalized prime-separating in the exact sense of retained [PSMT1](../../foundations/extended-geometry/FACTORIALLY_WEIGHTED_MULTITORUS_THEOREM.md), Section5. For H(s,t) with H(0,0)=0, let chi=i-j be the angular charge of a monomial s^i t^j with nonzero coefficient. Then

    H(S,T) is safe iff all its nonzero charges have
    ONE STRICT SIGN.                                      (2)

The zero polynomial is safe separately. A monomial of charge0 violates the condition. If both positive and negative charges occur, the condition fails even when special complex coefficients cause finite moment cancellations.

For necessity, use the injective algebra map

    iota(s)=z, iota(t)=r*z^(-1)
    into K[r,z,z^(-1)].

The distinct pairs (j,i-j) distinguish all monomials. Equation(1) says that the marked functional is exactly Phi_a composed with iota. Safety of H(S,T) forces eventual vanishing for this one fixed marker. PSMT1 therefore excludes0 from the convex hull of its angular support. In one angular dimension this means all charges have one strict sign. For sufficiency, if eta=min|chi|>0, a degree-delta source marker has charges in[-w*delta,w*delta], whereas H(S,T)^m has one sign and absolute charge at least m*eta. Thus

    L[Q H(S,T)^m]=0 for m*eta>w*deg Q.                      (3)

This controls ALL source markers without assuming they belong to the pullback algebra. Every H(S,T) in the zero-constant algebra is pure-null by the joint hypothesis. This is a reusable marker-certificate theorem, not a statement that every Gaussian phase admits such a marker.

## 2. The known pair has a universal quadratic certificate

For the retained six-variable cotangent pair, the joint inverse gives

    U=(1-2uv)^(-1/2), X*Y=U-1,
    L[(1+xy) exp(uS+vT)]=U.

Consequently

    L[(1+xy) S^i T^j]=delta_ij*(2i)!/2^i.                  (4)

Here Q0=1+xy has degree2 and a_n=(2n)!/2^n, including a_0=1. Since H(S,T) is pure-null, its positive powers have the SAME marked moments against xy and against1+xy. Thus a single fixed QUADRATIC marker xy detects every unsafe member of this whole algebra.

The required arithmetic is direct. For any integers0<=n<M and any odd prime p,

    a_d/a_(np)=2^(np-d)*(2d)!/(2np)! is p-integral
       whenever np<=d<=Mp.

If n<j<=M, the factorial ratio from2np to2jp contains p*(2n+1), so

    a_(jp)/a_(np) belongs to p*Z_(p).                       (5)

These are precisely the normalized PS1--PS2 conditions. They hold for all odd primes, so the exposed-face/Frobenius proof gives more: for every H whose angular support contains0 in its convex hull, there is k>=1 such that

    L[xy H(S,T)^(k*p)]!=0                                  (6)

for every sufficiently large prime p. The integers k and the exceptional-prime threshold are not computed here. The proof is the retained face theorem with this explicit normalization; this is not a numerical moment-prefix inference.

The reusable source is credited: Long's factorially weighted exposed-face argument and the repository's PSMT1 arithmetic extension. The [primary torus manuscript](sources/long_torus_719ca23022821868320fc0a5b9824c051105ff7e.tex) was read for its complete face,good-reduction,and factorial-isolation proofs. The rational-radial extension is the already retained local theorem. Neither its old assurance nor its rank-one formalization is promoted to the new application.

The source variables have charges (1,-1,-2,-1,1,2), so w=2. Combining(2)--(5) gives the complete all-degree FIELD classification:

* H=0 is safe.
* H is safe exactly when all occurring i-j are positive, or all are negative.
* In that case m*min|i-j|>2deg Q is sufficient.
* Every other H is pure-null but unsafe, with the SAME fixed marker xy and an infinite prime-index subsequence.

The common coefficient2 works for every safe H in every recomposition degree. It is sharp for this whole family because S,T already have exact linear-marker values L[p3*S^2]=-15 and L[z*T^2]=-8. This restricted-family uniformity does not contradict the retained absence of a degree-free cutoff across all Gaussian-safe phases.

All statements extend from C to arbitrary characteristic-zero fields by embedding the finitely generated coefficient subfield into C and reflecting nonvanishing. Fixed complex linear covariance conjugation gives ordinary independent Gaussian coordinates, preserving marker degrees. The theorem classifies this explicit cotangent-generated algebra; it is not a general support-only classification of specialized Gaussian polynomials.

## 3. A complete reduced safe scheme in every coefficient degree

Fix recomposition degree d>=1 and parameterize

    H(s,t)=sum_(1<=i+j<=d) c_ij*s^i*t^j.

Partition the coefficients into P for i>j,N for i<j,and Z for i=j. Let I_pos=(N,Z),I_neg=(P,Z). The reduced safe locus is the union of the two coordinate subspaces, and its radical ideal is explicitly

    J=I_pos intersection I_neg=(Z, c_p*c_n : p in P,n in N). (7)

This equality is an ordinary monomial-ideal calculation. J is radical as an intersection of coordinate prime ideals. Every specialized coefficient choice is covered; there is no generic-coefficient exception. This is the reduced SAFE ideal,not an assertion that a full ordinary mixed-cutoff ideal equals J.

Over ANY coefficient K-algebra A, the complete safety criterion is

    H(S,T) is safe iff the image of J is nilpotent.          (8)

Necessity follows by reducing to every prime fraction field and applying the field classification. The finitely many generators of J then form a finitely generated nil ideal,which is nilpotent without assuming A Noetherian or its full nilradical nilpotent.

There is a direct quantitative proof of sufficiency. Suppose (JA)^L=0. Expand a term of H^m using n_+ positive-charge,n_- negative-charge,and n_0 neutral factors. Its coefficient belongs to

    J^(n_0+min(n_+,n_-)),

by pairing opposite-sign coefficient factors. A surviving term must have n_0+min(n_+,n_-)<=L-1. Every nonzero charge has magnitude between1 andd. If n_+>=n_-,balance with a degree-delta source marker requires n_+<=d*n_-+2delta. The opposite case is symmetric. Hence every nonzero mixed contribution has

    m<=2delta+(d+1)(L-1).

Therefore

    m>2deg Q+(d+1)(L-1) implies L[Q H(S,T)^m]=0,           (9)
    C_L=2+(d+1)(L-1)

is a common sufficient INTEGER linear coefficient,using purity for constant markers. For L=1 this recovers coefficient2. These thickness bounds are sufficient,not claimed optimal. A disconnected reduced algebra may have positive coefficients on one component and negative ones on another;J=0 handles this without forcing one global sign. It is wrong to require every coefficient of H to be nilpotent,which was correct only in the neutral subfamily h(ST).

If I_C denotes the full ordinary cutoff ideal of the universal coefficient family,(9) gives I_(C_L) contained in J^L as an ordinary ideal;the family is already identically pure. For integer N>=2,radical(I_N)=J by the field classification. Noetherianity gives the opposite cofinal containment of some power of J in I_N. Thus the actual cutoff-completion functor is the formal completion along this explicit union of subspaces. Individual mixed levels are not computed or identified with powers. The previously exact lines(e^(N+1)) and(e^(2N+1)) remain unchanged and illustrate why levels matter.

## 4. Exact normalization and a positive-dimensional pure base-change obstruction

Let H_P be the exact Gaussian source for P=H(S,T),and let its polynomial-Hom normalization ideal be

    N_H={a in A:a*L[Q exp(tP)] is polynomial in t for every Q}.

The [exact formal-Hom theorem](GAUSSIAN_FORMAL_COMPLETION_AND_POLYNOMIAL_HOM.md) and(8) give

    N_H=Gamma_J(A)=union_(k>=0) Ann_A((JA)^k).               (10)

For a fixed scalar a,its marked tails vanish exactly when P is safe over A/Ann_A(a). Criterion(8) says precisely that a power of JA annihilates a. Safety requires1 in N_H;N_H!=0 alone is insufficient.

In the UNIVERSAL polynomial coefficient ring A=K[c_ij],J is nonzero for every d>=1 because both c_10 and c_01 occur and their product lies in J. The domain therefore has Gamma_J(A)=0. Yet every fiber in either positive-dimensional safe coordinate subspace has a one-dimensional polynomial Hom space. This is an entirely pure-null family with base-change failure on two nontrivial safe components,not only at the zero phase. It does not invalidate valid fiberwise Hom methods or the already closed field locus.

The all-ring functor in(8) is a formal completion,not the literal equation JA=0 at every safe point. For example nilpotent neutral coefficients give safe points with nonzero JA. The normalization ideal is not a prescribed ordinary cutoff ideal.

## 5. Implemented classifier and remaining programme

The [checker/classifier](verify_complete_cotangent_phase_classification.py) accepts H in Q(i)[s,t] with H(0)=0. It returns the EXACT field safety classification from collected charge support,an explicit sufficient rational cutoff for positive cases,and the fixed quadratic bad marker for negative cases. It does not guess safety from sampled moments. A guarded coefficient extractor computes any requested small marked moment from(4). The [receipt](complete_cotangent_phase_verification.json) records exact small comparisons with the source Wick extractor,normalized prime-ratio controls,coefficient-ideal intersections and nilpotent charge bounds. The [source audit](complete_cotangent_phase_source_audit.json) preserves the primary git BLOB identity,retained theorem scopes and all protected records.

The written proof,using the retained PSMT1 theorem,establishes completeness and all-degree/all-ring statements. Finite controls are transcription checks;they do not prove prime isolation or certify its independent review. No general Gaussian classifier,large prime campaign,expanded degree49 inverse gate,or generic ordinary-ideal/Hom solver is implemented.

The explicit H(S,T) classification is CLOSED subject to independent review. Next study other cotangent algebras for a prime-separating marked period,or prove their actual module/marker reductions separately. Full moment-state necessity without inverse,other Gaussian normal forms/modules,OP-GVC-FORMAL-LINEAR,general producers and independent specialist review remain open. The global field94,least-prefix intervals,smaller GMC endpoints,exact scaling ideals and concurrent radial bridge are unchanged.
