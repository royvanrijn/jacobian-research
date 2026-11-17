# The contact507 branch: geometric-degree bounds and zero seven-jet compatibility

14 November 2026. This continues the actual leading-stratum analysis in
[JC2CFN1](PROOF.md). The whole F2 problem remains open.
Claims JC2CFMV1 and JC2CFJ1 below distinguish an unconditional support
bound from a conditional target-attachment test.

## 1. Mixed volume gives a geometric-degree bound without generic coefficients

Let d be the generic inverse-sheet degree of the F2 Keller map (P,g),
where g=-(5/9)Q. The source support envelopes are

    P0=conv((0,0),(3,0),(15,60),(0,75)),
    Q0=conv((0,0),(5,0),(25,100),(0,125)).

After subtracting generic target constants both polynomials include the
origin in their support. A generic fibre misses the two source axes:
their images have dimension at most one. The classical Bernstein upper
bound for isolated torus roots therefore applies to the entire generic
fibre, even though its coefficients lie in the special Keller locus.

Use the mixed-volume convention

    MV(A,B)=area(A+B)-area(A)-area(B).

Then MV(P0,Q0)=2175, giving d<=2175 directly from the retained envelope.

The imported upper-bound theorem is Bernshtein, The number of roots of
a system of equations, 1975, DOI10.1007/BF01075595.
The [pinned Philippon--Sombra author paper](sources/bkk-manifest.json)
states the upper bound, including its distinction from generic equality,
in its abstract and introduction. Standard toric-surface interpretation
of mixed volume is an explicit geometric input below; no independent
full proof of the external BKK theorem is claimed.

## 2. Use the smaller contact polynomial and count its three sheets

Put W=P^5-g^3 and let L be its highest eta=(5,-1) band.
JC2CFN1 gives

    L in {39,40,45,50,55,60,65,70},
    deg_(5,2)(W)=13L.

Thus its ordinary total degree T_W satisfies

    T_W<=T=min(375,floor(13L/2)).                         (2.1)

The map (P,W) has generic degree3d: for a generic target (a,b),
g^3=a^5-b has three distinct nonzero roots, and each corresponding
(P,g) target has d distinct preimages. The cubic extension is not
silently treated as birational.

Let B_(L,T) be the rational polygon cut out by

    i,j>=0, 5i-j<=L, i+j<=T.

It contains the actual support of W-b. The facet formula gives

    MV(P0,B_(L,T))=12L+15T.                             (2.2)

The two nonzero facet contributions of P0 have primitive lattice lengths
12 and15, with outward degree directions (5,-1) and (1,1).
Hence already

    3d<=12L+15T, d<=4L+5T.                             (2.3)

The checker independently obtains these same numbers from Minkowski sums
of the exact INTEGER hull of all allowed W monomials, also imposing
5i+2j<=13L. It does not substitute rational vertices for integer supports.

## 3. A mandatory boundary intersection improves every row

The carrier face has

    c(v)=v(v-1)^2 R(v), v=xy^5,

with a simple root at0, a double root at1, and the roots of R distinct
from0,1 with total multiplicity two. On the toric eta boundary use the
smooth chart

    q=y, v=xy^5, q=0, v!=0.

It is a genuine lattice chart: x=v/q^5,y=q.
The cleared equations of a generic fibre of (P,W) are

    F=q^15(P-a), H=q^L(W-b).                            (3.1)

We count their obligatory local intersections, which are included in the
mixed-volume intersection number but are absent from the affine source.

At a root of R of multiplicity mu, the reductions modulo q have orders
3mu and at least3mu. Thus both belong to (q,(v-r)^(3mu)); their
intersection length is at least3mu. Summing gives at least6.

At v=1 use the local parameter s=v^(1/5)-1. Taking this root is an
unramified formal coordinate change in characteristic zero. The retained
Puiseux/terminal chart is

    x=X^5, y_old=y_new+1/X,
    s=X*y_new, q=y_old=(1+s)/X.

Every translated P monomial X^a*y_new^b has 5a-17b<=3.
Writing ell=a-b, its contribution to F in (3.1) is

    q^(15-ell)*(1+s)^ell*s^b.

For local weights wt(q,s)=(5,12), every term has weight at least72:
75-(5a-17b)>=72. The generic target constant contributes q^15,
also above that bound.

Likewise the terminal height of W=P^5-g^3 is at most15. Every term
in H has local weight at least5L-15, including the generic constant.

For completeness, if two plane germs have weighted orders at least
A,B for weights(r,t), their local intersection is at least A*B/(r*t).
Substitute q=u^r,s=v^t. This finite flat local map has degree rt, so
it multiplies the intersection length by rt. The ordinary intersection
inequality is at least the product of the two multiplicities. That
inequality follows by summing branch valuations: a function of
multiplicity n has order at least n times the branch multiplicity on
every normalized branch. This proves the stated weighted bound.

Here it yields

    I_(v=1)(F,H)>=72*(5L-15)/60=6L-18.                  (3.2)

Thus the eta boundary alone consumes at least6L-12 intersections.

The toric intersection bookkeeping is as follows. Use the normal toric
surface of P0. P0 is integral, so its divisor D_P is Cartier.
All four facet support values of P-a are attained: the constant at the
origin, the forced total-degree75 vertex, and the retained eta-leading
edge. Hence its zero divisor contains no toric boundary component.
The support values (0,L,T,0) define an effective Q-Cartier toric class
for W-b, with intersection D_P.D_W=12L+15T. It may have extra boundary
components when the ordinary degree is below T; none is shared with
the generic P divisor. The pair has no common curve in the torus
because its generic fibre is finite. Their proper intersection is
therefore the sum of nonnegative local contributions. The eta points
above are smooth and contribute their ordinary plane intersection lengths.

Subtracting only the explicitly proved contributions gives

    3d+(6L-12)<=12L+15T,
    d<=2L+5T+4.                                        (3.3)

No Newton nondegeneracy or equality in the BKK theorem is assumed.
Uncomputed boundary contributions can only improve this upper bound.

## 4. Exact eight-row geometric-degree atlas

| highest band L | cusp contact | total-degree cap T | plain mixed-volume bound | boundary-refined bound |
| ---: | ---: | ---: | ---: | ---: |
| 39 | 507 | 253 | 1421 | 1347 |
| 40 | 520 | 260 | 1460 | 1384 |
| 45 | 585 | 292 | 1640 | 1554 |
| 50 | 650 | 325 | 1825 | 1729 |
| 55 | 715 | 357 | 2005 | 1899 |
| 60 | 780 | 375 | 2115 | 1999 |
| 65 | 845 | 375 | 2135 | 2009 |
| 70 | 910 | 375 | 2155 | 2019 |

Consequently every actual pair in this retained F2 support/carrier stratum
has geometric degree d<=2019. The lowest contact507 branch has d<=1347.
The retained squarefree/double lower floors6/12 remain unchanged.

This excludes larger sheet-degree packets within each contact stratum;
it does not exclude any contact stratum entirely or establish existence
at an endpoint. These are GEOMETRIC degrees, not coordinate degrees.
The conditional coordinate-degree frontier125 is not improved.

As a local bound control, the germs s^6+q^15 and s^15+q^36 have local
intersection216: their resultant in s is q^216*(1+q^3)^3.
They meet the weighted bound at L39. This is a local control, not a
full F2 source or a prescribed terminal-normal-form realization.

## 5. Exact zero seven-jet locus and its Euler connection

The source profile L39 has, on the generic carrier,

    pi=P^3/g^2 ~ q^5/c(v), h=P^5/g^3,
    h-1=W/g^3=O(q^36).

Thus its first seven integer carrier centres are zero.
This does NOT force every nonproper-value component to share those
centres: the global attachment is an additional hypothesis, as emphasized
by the retained [puncture audit](../../foundations/plane-jc/F2_AFFINE_PURITY_PUNCTURE_ATTACHMENT.md).

We can nevertheless solve the precise compatibility problem completely.
Restore all affine target transports of a normalized k1 parametrization,
after normalizing leading scalars:

    p=t^3+a*t+mu,
    q=t^5+b*t^4+eta*t^3+c*t^2+(d+eta*a)*t+nu.            (5.1)

Here pi=p^3/q^2=t^(-1)*(1+...), h=p^5/q^3=1+....
Vanishing of its first seven carrier coefficients is equivalent to

    deg_t(p^5-q^3)<=7.                                 (5.2)

It is therefore cut out exactly by the seven coefficients at t degrees
14,13,...,8. Their complete ideal has Groebner basis

    b,
    3*eta-5*a,
    3*c-5*mu,
    9*d+10*a^2,
    9*nu-10*a*mu,
    9*mu^2-a^3,
    a^2*mu,
    a^5.                                               (5.3)

A compact certificate comes from the first five constant triangular
pivots. Substitute

    b=0, eta=5a/3, c=5mu/3,
    d=-10a^2/9, nu=10a*mu/9.

The last two high coefficients become exactly

    -(5/27)*(a^3-9mu^2), -(5/9)*a^2*mu.                (5.4)

Finally a^5=9*(mu*(a^2*mu)-a^2*(mu^2-a^3/9)).
The seven standard monomials are1,a,a^2,a^3,a^4,mu,a*mu.
All target transports are retained; no unjustified zero-translation
normalization is used.

Eliminating mu,eta,nu from (5.3) gives exactly

    (b,d+10a^2/9,c^2-25a^3/81,a^2*c),

the adapted-quintic Euler-eligibility ideal of JC2CFQ1.
Indeed mu=3c/5, eta=5a/3 and nu=2ac/3 give an explicit scheme
isomorphism. Both schemes have length7 and just one geometric point,
the monomial E8 endpoint with all seven parameters zero.

Thus Euler eligibility of this adapted quintic is equivalent to the
existence of a zero seven-jet transport, scheme-theoretically.
Any NON-E8 k1 branch that is required to follow the L39 zero-centre
profile is excluded. A mere nonproperness component is not automatically
such a branch. The E8 global endpoint and other attachments remain open.

## 6. Exact next gate

The contact507 branch now has a fixed degree51 Wronskian, geometric
degree at most1347, and a completely solved zero-seven-centre k1
compatibility test. The missing step is still to couple its actual
source coefficients, or its actual global boundary attachment, to
one of these obstructions. No missing historical circuit was rebuilt.

The all-parameter curve/toric arguments above are written proofs.
The checker verifies exact polygon arithmetic, the seven coefficient
ideal with literal triangular reductions, and local intersection controls.
It does not formally verify the external BKK theorem or supply a global
Keller realization or exclusion.
