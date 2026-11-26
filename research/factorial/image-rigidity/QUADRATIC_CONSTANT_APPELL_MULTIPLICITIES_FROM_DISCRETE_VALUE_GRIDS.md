# Critical Appell multiplicities from a discrete coefficient-value grid

This note proves a multiplicity bound and an all-constant window transfer from a discrete grid for the nonconstant coefficient field. It uses only finite polynomial reductions and factorial valuations. The [native coefficient companion](QUADRATIC_CONSTANT_NONCRITICAL_DYADIC_WINDOW_TRANSFER.md) supplies the complete divided-factorial integrality and binary residue identity used in Section4. No completeness, Hessian-rank or fixed-correlation hypothesis is imposed.

Let E be any characteristic-zero nonarchimedean valued field with v(2)=1 and

                  v(E^*) contained in (1/e)Z               (1)

for some positive integer e. Let R_0=1 and R_r in O_E for every r>=1. Put

 P_m(Z)=sum_(r=0)^m (m)_r R_r Z^(m-r),
 (m)_r=m!/(m-r)! .                                        (2)

The polynomial P_m is monic of degree m. The constant k may lie in ANY valued extension of E, even if the extension valuation initially has higher rank. If k is transcendental over E, every P_m(k) is nonzero and there is nothing to prove. Otherwise let F be its minimal polynomial. At a zero of F at least two distinct summands have equal least valuation. Hence (j-i)v(k)=v(a_i)-v(a_j) belongs to (1/e)Z for some i<j. An ordered value group is torsion-free, so v(k) belongs to the divisible hull of this rank-one grid and is a rational real number. Thus the real-valued case partition below is exhaustive for every possible zero. The case k=0 is handled separately. For native moments of a fixed nonconstant divided-power-integral polynomial q over O_E, set R_r=Lambda(q^r)/r!: these are integral and P_m(k)=Lambda((q+k)^m). The argument below applies to any integral sequence R.

## 1. Multiplicative widths at a specified valuation

For a nonzero polynomial A(Z)=sum a_i Z^i over E and a real number delta realized by an element k!=0 in a valued extension, let I_delta(A) be the nonempty finite set of indices attaining

                         min_i(v(a_i)+i delta).

Define

           width_delta(A)=max I_delta(A)-min I_delta(A).   (3)

For nonzero A,B one has

           width_delta(AB)=width_delta(A)+width_delta(B). (4)

To prove this, work in the given valued extension and replace Z by kU. Divide A(kU) by any one of its least-valuation coefficients. The resulting polynomial has integral coefficients and a nonzero reduction; its smallest and largest nonzero reduced exponents are min I_delta(A),max I_delta(A). Do the same for B. The two reductions multiply in a polynomial ring over a residue FIELD, which is an integral domain. Therefore their lowest nonzero exponents add and their highest nonzero exponents add. This also proves that the minimum valuation of a product is the sum of the two minima, including any cancellation in intermediate coefficients. Equation(4) follows. No completeness or convergence is involved.

If A(k)=0, then width_delta(A)>0: a unique least-valuation summand cannot sum to zero. If delta=(L-1)/L with L=2^a and a>=1, equality of the endpoint costs at i<j in I_delta(A), together with(1), gives

                 delta*(j-i) in (1/e)Z.

Since L-1 and L are coprime, L divides e*(j-i). Consequently, whenever A(k)=0,

          L/gcd(L,e) divides width_delta(A),
          width_delta(A)>=L/gcd(L,e).                     (5)

## 2. Exact critical multiplicity bound

**Proposition.** Let L=2^a,a>=1, theta=1-1/L, and t=m modL. Suppose R_t is a UNIT. Every root k of P_m with v(k)=theta has multiplicity at most

                         g=gcd(L,e).                     (6)

If the a-th bit of m is zero, no such root exists. If it is one and such a root exists, width_theta(P_m)=L EXACTLY.

The coefficient of Z^(m-r), evaluated at valuation theta, has valuation

 v((m)_r)+(m-r)theta+v(R_r).

Legendre's exact factorial formula, after subtracting m*theta, gives the base cost

             r/L+s_2(m-r)-s_2(m),                         (7)

with the extra v(R_r)>=0 left outside it. Write j=m-r=2L u+w,0<=w<2L. Minimizing(7) is equivalent to minimizing

                   s_2(u)-2u+s_2(w)-w/L.

The high part strictly decreases with u, since its increment is at most1-2<0. The low part is nonnegative and is zero exactly at w=0,L; all binary-place weights below place a are positive and the place-a weight is zero. Thus use the largest allowed u. If m mod2L=t<L, the unique minimizing r is t. If m mod2L=L+t, the only minimizing r are t,t+L. Every other r has a STRICTLY larger base cost. This is exact for every m, including m<L.

Since R_t is a unit, the first case cannot give a root. In the second case a root forces R_(t+L) to be a unit too: otherwise the r=t term would be the unique least term. Hence the only least-cost exponents are m-t-L and m-t, separated by L. This proves width_theta(P_m)=L.

Let F in E[Z] be the monic minimal polynomial of k. In characteristic zero F is separable. If k has multiplicity h in P_m, polynomial factorization over E therefore gives F^h dividing P_m. By(4) and nonnegativity of every polynomial width,

                       L>=h*width_theta(F).

By(5), width_theta(F)>=L/g. Therefore h<=g, proving(6). The proof does not count conjugates by their valuations, extend a valuation uniquely, or assume E is complete. QED.

## 3. A positive residue window transfers to all constants

Assume w>=1 and every block of w consecutive POSITIVE indices of R_r contains a unit. Write b=v_2(e). Then for every k in any valued extension of E, the sequence P_m(k),m>=1, has no block of

                              w+2^b                        (8)

consecutive zeros. At the particular critical valuation theta=1-2^(-a), the sharper bound is w+2^min(a,b).

First take a critical k and put g=2^min(a,b). Suppose P_n(k),...,P_(n+w+g-1)(k) all vanish, where n>=1. Among the LAST w indices m=n+g,...,n+g+w-1 some R_(m modL) is a unit. Indeed, if the residues wrap through0 or contain a full period, R_0=1 suffices. Otherwise they form a w-block of positive R-indices. This argument also covers w>=L.

The exact Appell identity is

                 P_m^(j)(Z)=(m)_j P_(m-j)(Z).             (9)

The zero block contains P_m(k),P_(m-1)(k),...,P_(m-g)(k). Thus k has multiplicity at least g+1 in P_m, contradicting(6).

The other valuations give bounds no larger than(8). If v(k)<0, k^m is the unique least-valuation term in(2). If v(k)=0, every even m gives a unit after division by k^m, because (m)_r is even for every r>=1. Thus a two-window suffices.

If 0<delta=v(k)<1 is not critical, choose L=2^a with1<L(1-delta)<2. For j=Lu+c,0<=c<L, the base cost apart from a common constant is

 s_2(u)-L(1-delta)u+s_2(c)-(1-delta)c.

The high part strictly decreases; the low part is positive unless c=0, because every low binary-place weight is positive. Hence the UNIQUE minimizing r is m modL. A unit at that R-index forces P_m(k)!=0. The same circular-prefix reasoning gives a w-window.

If v(k)>=1, the exact expansion

 P_m(k)/m!=sum_(r=0)^m R_r k^(m-r)/(m-r)!

has residue Rbar_m: every term with m-r>=1 has valuation at least (m-r)-v((m-r)!)=s_2(m-r)>0. A w-window suffices again. For k=0 use P_m(0)=m!R_m. Since w+2^b>=max(w,2), all cases prove(8).

## 4. Native binary consequence and exact scope

Take a naturally integral binary quadratic

 q=alpha X^2/2+beta XY+gamma Y^2/2+dX+eY

with alpha,beta,gamma,d,e in O_E. Its native residue generating series is

 Rbar(T)=N(T)/D(T),
 N=1+(dbar+ebar+betabar)T+dbar ebar T^2,
 D=(1+alphabar T+dbar^2 T^2)(1+gammabar T+ebar^2 T^2).       (10)

The complete native proof of(10) is in Section5 of the [coefficient companion](QUADRATIC_CONSTANT_NONCRITICAL_DYADIC_WINDOW_TRANSFER.md); its integrality proof is Section1. These identities are independent of fixed-weight or rank-one hypotheses. Suppose this series is NONPOLYNOMIAL. It has a positive unit-window w<=4. Here is the elementary rational-series bound: if original degD>=3, Rbar is proper and its reduced recurrence of order at most4 has a nonzero final coefficient, so a zero block of that order would propagate backward and forward to the zero series. If degD is1 or2, a polynomial part can have degree at most2-degD; after its contribution is removed from a block at positive indices, a suffix as long as the reduced proper recurrence remains. Thus a positive two-window suffices in these cases. DegD=0 would make the series polynomial and is excluded.

Consequently, for ALL constants k in valued extensions:

* if e is odd, q+k has no five consecutive zero factorial moments;
* if e is twice an odd integer, q+k has no six consecutive zero factorial moments.

The six-window gives actual-support SFC for a phase with all six binary-quadratic monomials. These are uniform bounds; sparse phases are not assigned a smaller actual-support window merely by restriction. Polynomial residue series remain outside this corollary. So do coefficient fields for which no indicated small grid denominator is available. Enlarging the coefficient field can increase e and weaken the bound, and no specialization is claimed that gives arbitrary coefficients an odd or twice-odd value grid.

The odd case recovers the parent's initial complete-local-field argument. The exact gcd bound extends it: every critical root with a unit anchor has multiplicity at most the power of2 already present in the COEFFICIENT value grid. The constant itself need not belong to that grid. This distinction is essential.

