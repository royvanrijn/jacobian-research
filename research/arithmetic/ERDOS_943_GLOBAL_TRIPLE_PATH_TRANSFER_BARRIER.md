# Global square-full triple counts transferred to selected paths: a numerical barrier

## 1. Graph, dyadic block and notation

Fix a positive real T>=2 and an integer n. Consider a simple bipartite graph with squarefree positive kernel vertices b,d in [T,2T), and select one fixed positive representation

    a^2*b^3+c^2*d^3=n

on each edge, with a,c in [T,2T) and gcd(a^2*b^3,c^2*d^3)=1. These are coefficient-one actual representations. Do not replace them by arbitrary congruence witnesses. A nonempty graph automatically has 2*T^5<=n<64*T^5.

An ordered shared-right path is (b1,d,b2), with b1!=b2 and the two selected edge labels (a1,c1),(a2,c2). Put

    alpha=gcd(a1,a2), a_i=alpha*r_i,
    g=gcd(b1,b2),     b_i=g*u_i,
    L=(c1^2-c2^2)/(alpha^2*g^3).

Actual edge coprimality makes L a nonzero integer. The reduced kernels u1,u2 are coprime, squarefree and distinct, and

    r1^2*u1^3-r2^2*u2^3+L*d^3=0.                      (1)

These statements are the actual reduced-cubic lift in the existing Pell-recovered-kernel-gcd source. Restrict to a dyadic block alpha in[A,2A), g in[G,2G), with real A,G>=1. If nonempty, A,G<=2T and

    A^2*G^3<4*T^2.                                    (2)

Indeed alpha^2*g^3 divides the nonzero difference c1^2-c2^2, whose absolute value is less than 4*T^2.

Write B0=T^5/(A^2*G^3), and B=32*B0. Then B0>T^3/4>=2 in a nonempty block. Let P_sf(A,G) count only those ordered paths in this block for which the ADDITIONAL property

    |L|*d^3 is square-full                             (3)

holds. Property (3) is not automatic. If |L|=e*v^2 with e squarefree, it can fail at a prime having exponent one in L and not dividing d.

For H>=1 let Gamma(H) count ordered positive triples (s,t,w) of square-full integers with s+t=w<=H. It includes nonprimitive triples.

## 2. Exact transfer with subpower fibers

**Theorem 1.** For every epsilon>0,

    P_sf(A,G) <<_epsilon n^epsilon*Gamma(32*B0).        (4)

More generally, any bound Gamma(H)<<_(theta,eta) H^(theta+eta), for a fixed theta>=0 and every eta>0, implies

    P_sf(A,G) <<_(theta,epsilon)
                    n^epsilon*B0^theta.               (5)

The constants are uniform in n,T,A,G and in the selected graph. This statement is conditional on the stated bound for Gamma, not a new estimate for that global counting function.

Put m_i=r_i^2*u_i^3 and w=|L|*d^3. Both m_i are positive square-full. Equation(1) says m1+w=m2 if L>0, and m2+w=m1 if L<0. Thus the path maps to an ordered square-full triple by

    L>0: (s,t,w_total)=(m1, L*d^3,m2),
    L<0: (s,t,w_total)=(m2,-L*d^3,m1),

together with one sign bit. Here w_total is the largest of the three terms. Since r_i<2T/A and u_i<2T/G, both m_i are less than B=32*T^5/(A^2*G^3); the difference is smaller still. This proves the asserted height bound. No coprimality assumption on the square-full triple is required by Gamma.

It remains to bound the fiber uniformly at the FIXED original n. Fix the triple and the sign bit. The two known m_i have unique square-times-squarefree-cube decompositions, recovering r_i and u_i. For the known positive difference t, enumerate squarefree d with d^3 dividing t; there are at most tau(t) choices, and each fixes the nonzero signed integer L. Further caps, coprimality and agreement with actual selected edges can only discard these candidates.

For each such choice use the existing local norm recovery. The exact relation is

    r2^2*u2^3*c1^2-r1^2*u1^3*c2^2=n*L.

The integers

    x=r2*u2^2*c1, y=r1*u1*c2,
    D=u1*u2,      M=u2*n*L

therefore satisfy x^2-D*y^2=M. Here D>1 is squarefree: u1,u2 are distinct coprime squarefree integers. The signed M is nonzero. All parameters and variable heights are bounded by fixed powers of n: r_i,u_i,c_i,d<2T, |L|<4*T^2, D<4*T^2, |M|<8*n*T^3, |x|<16*T^4 and |y|<8*T^3. Since n>=2*T^5 these are uniform polynomial bounds.

The elementary norm lemma in ERDOS_943_PELL_RECOVERED_KERNEL_GCD_WEIGHT.md, Lemma1, bounds these norm solutions by n^epsilon after resizing epsilon. This applies to either sign of M, nonprimitive x,y, and all target valuations. The map from x,y back to c1,c2 is injective for the fixed r_i,u_i.

The positive quantity

    X=(n-d^3*c1^2)/(r1^2*u1^3)

then has at most one expression X=alpha^2*g^3 with g squarefree and alpha positive integral; prime-valuation parity fixes g and then alpha. Nonpositive, nonintegral or impossible values are discarded. It fixes b_i=g*u_i and all selected edge data. Hence each remaining tuple gives at most one ordered actual path.

Finally t<=B<=32*T^5<=16*n, since A,G>=1. Thus tau(t)<<_epsilon n^epsilon. Combining the sign factor, divisor choices and uniform norm count proves(4). It is important that the original n is fixed during the norm count. No independent sum over g, alpha or original targets has been introduced. Applying a bound for Gamma and using B0<<n proves(5) after allocating smaller epsilons.

This transfer does NOT require the primitive diagonal-cubic theorem as a separate completion step: the global triple count has already counted both reduced square-full values. It uses the local norm/recovery route only. In particular no primitive-triple restriction is silently imposed, and L's signed or vanishing cases are not merged: zero was excluded by actual distinct kernels.

## 3. The two existing all-path bounds

**Lemma 2.** The number P_all(A,G) of ALL ordered actual paths in the same block, without(3), satisfies

    P_all(A,G) <<_epsilon n^epsilon
          * min(T^2/G, T^4/(A^4*G^3)).                  (6)

For T>16 the general actual codegree theorem applies with all four caps2T, since 2*((2T)*(2T))^2=32*T^4<n. For a fixed g in[G,2G), ordered kernels both divisible by g number at most(2T/g)^2. Summing over the O(G) positive integers in that range gives O(T^2/G) kernel pairs; restricting to exact gcd g, squarefreeness and distinctness only decreases this. The actual codegree theorem supplies n^epsilon right neighbors for each pair. Each pair and common right vertex specify at most one selected path. This proves the first term of(6).

For the second term, count the coefficient data BEFORE completing them. There are O(T^2/A^2) choices of ordered r1,r2. The nonzero integer L obeys

    0<|L|<4*T^2/(A^2*G^3).

The number of such signed integers is O(T^2/(A^2*G^3)); if the range is smaller than one there are no paths. The existing fixed-data proposition in the Pell source supplies n^epsilon actual selected completions per(r1,r2,L). It combines the coefficient-uniform primitive cubic theorem with local norm recovery, and crucially does not sum over g. This proves the second term of(6).

When 2<=T<=16, the roots and kernels belong to an absolute bounded range, hence n<64*16^5. A nonempty block has A,G<=2T and positive right sides in(6), bounded away from zero by an absolute constant on that compact range. Increasing the constants therefore covers this finite initial range as well. Empty graphs or empty dyadic blocks are immediate. This establishes(6) uniformly for every stated T.

## 4. Numerical domination for EVERY global triple exponent at least one half

**Proposition 3.** If theta>=1/2, the majorant B0^theta in(5) cannot improve the existing minimum in(6), even on the additional square-full-difference sector. The assertion concerns numerical majorants, not an impossibility of a sharper arithmetic count.

Define the two positive majorants

    U=T^2/G, V=T^4/(A^4*G^3).

There is the EXACT identity

    sqrt(B0)=U^(3/4)*V^(1/4).                           (7)

It holds without taking logarithms or assuming power-scale A,G. By(2), B0>=1 for every nonempty block. Therefore theta>=1/2 gives

    B0^theta >= sqrt(B0)
              =U^(3/4)*V^(1/4) >= min(U,V).             (8)

Thus adding(5) to the already-valid menu(6) does not lower its pointwise exponent envelope. Constants and independently resized n^epsilon factors do not change this comparison of principal powers. This is not a lower bound for P_sf, P_all, or Gamma.

Equivalently, when A=T^a,G=T^b, the three exponent expressions obey

    (5-2a-3b)/2
       =(3/4)*(2-b)+(1/4)*(4-4a-3b).

Both endpoints and a term below one are handled by the scale identity(7) and the nonempty condition(2); no unjustified assumption that an arbitrary counting weight is at least one is made.

## 5. Exact primary scope and the remaining useful target

The currently listed official arXiv record for Heath-Brown2601.07817 has only v1. Its actual Theorem1.1 counts ALL positive square-full triples, including nonprimitive ones, with

    Gamma(H)<<_eta H^(theta_HB+eta),
    theta_HB=3/5-3/1555=186/311>1/2.

Its separate Theorem1.3 is coefficient-uniform for primitive signed diagonal cubic zeros. The latter is an indispensable input to the existing fixed-data bound used in(6); it is not the same statement as the global theorem used conditionally in(5). Source: https://arxiv.org/html/2601.07817v1 .

Consequently Theorem1.1 really does transfer to the restricted path count by(5), but Proposition3 shows why this coarse transfer provides no new power saving over(6). The sector(3) itself is additional and does not exhaust the paths or cycles.

The older 7/19 numerical exponent also needs its exact source distinguished. The repository's prescribed-sum7/19 argument in ERDOS_943_CUBE_MODULUS_BOXES.md uses Chan1407.0054's coefficient-unit box-completion theorem, whereas Chan2012's7/19 is a fixed-difference estimate. Blomer--Schoebel's stronger fixed-shift result is not a uniform prescribed-sum theorem; their uniform positive-sum Corollary2 gives2/5. The current dual graphs already give3/10. None of these distinctions is changed by this note.

A potentially useful next target is a joint ORIGINAL-n estimate on the norm-realizable reduced triples in the proof of Theorem1, rather than an unrestricted Gamma(B) bound. Such an estimate would retain the fixed norm M=u2*n*L and the selected-edge caps while summing the reduced cubic data. Alternatively a uniform averaged/weighted line estimate in Heath-Brown's Sections4--6 could be useful only after those same fixed-target conditions and its frequency cost are tracked. The numerical barrier here excludes neither approach, nor a finer rectangular triple theorem or a different map. It excludes only treating the unrestricted global triple count as an improvement to these existing path bounds.

No full-C6 estimate, global improvement below3/10, fixed-target2/7 for every integer, or novelty/priority claim is made.


## 6. Every valid unrestricted global power exponent meets this barrier

**Corollary 4.** For H>=1,

    Gamma(H)>=floor(sqrt(H)/5).

Consequently, if Gamma(H)<<_eta H^(theta+eta) for every eta>0, then theta>=1/2. Thus no admissible exponent-only upper bound for this unrestricted global Gamma can improve the majorant minimum in(6) through the coarse transfer(5).

For every positive integer j<=sqrt(H)/5, the ordered triple

    (9*j^2,16*j^2,25*j^2)

consists of positive squares, hence square-full integers, and satisfies the required equation. The triples are distinct. If theta<1/2, choose eta=(1/2-theta)/2>0; the displayed lower bound contradicts the upper bound as H tends to infinity. The quantifier over every eta is essential here.

The numerical threshold is exact for this menu: for T>1 set A=sqrt(T),G=1. Then B0=T^4 and U=V=T^2. Thus B0^theta<min(U,V) holds exactly when theta<1/2. This is a numerical example only, not an existence assertion for actual paths or triples in the fiber. In particular the square-family lower bound for unrestricted Gamma does not assert that those triples are norm-realizable at a prescribed original n. Refinements that retain selected-edge information, coprimality or the original-n norm condition are not ruled out.
