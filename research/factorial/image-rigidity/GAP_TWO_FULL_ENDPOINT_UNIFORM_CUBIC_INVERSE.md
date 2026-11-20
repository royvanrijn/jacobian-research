# A fixed regular cubic inverse across the full actual endpoint window

The interval input is an outward-rounded boundary certificate on closed cells. The conclusion is a true inverse/Jacobian bound; no complete-period C_R or numerical Strong Factorial threshold is inferred.

Let |delta|<=10^(-8) and

|Re z|<=1/100000, 18999/10000<=Im z<=21001/10000.            (F1)

Set u=v+z/2 and use the physical disk |u|<47/50. Let F_delta, Phi and zeta be the exact phase and matched data from QUANTITATIVE_GAP_TWO_ENDPOINT_DATA.md. This rectangle lies inside its |z-2i|<=101/1000 domain. Normalize the cubic by v_t=i and zeta_z=-i at the limiting upper merger.

**Theorem FullActualEndpointInverse.** The equation

F_delta(v,z)-Phi(z,delta)=t^3/3-zeta(z,delta)*t              (F2)

has a selected inverse v(t,z,delta) holomorphic and regular on |t|<91/200. It satisfies |v+z/2|<47/50, matches both physical stationary points with t=+/-sqrt(zeta), and obeys

|v_t|<19 on |t|<=2/5,
|G_t|<350000 on |t|<=3/8,                                 (F3)

where h=v_t=A+B*t+(t^2-zeta)*G. The statement concerns the actual thin coefficient rectangle(F1), not the whole complex disk of radius1/10.

## 1. The complete interval input

The primary supporting computation is verify_hermite_endpoint_outer_boundary_ball.py, SHA2569041c5ba279121d17986fe4167a809a14c43ee83f228f563c7814cdba0fcac35. It evaluates

W_0(u,iy)=F_0(-iy/2+u,iy)-Phi(iy,0)

on u=(47/50)*exp(i*theta), for EVERY theta in[0,2pi] and y in[18999/10000,21001/10000]. It proves

|W_0|>9/100.                                              (F4)

This is not a mesh test. The512 exact rational root rectangles cover the full closed parameter set; every split uses exact rational midpoints and two closed children covering its parent. Each surviving leaf is accepted only when the outward Arb/Acb modulus ball is strictly above9/100. All1300 leaf boxes pass at100-bit precision, with2088 attempted boxes and maximum depth7. The leaf stream hash is ac7b501b1144b46b929a6711185955b9df2582756ab87b241e1d5630aadf6a31. Independent replays reproduce those counts and the stream. Throughout, Im(iy/2+u)>99/10000, which fixes the logarithm's analytic upper-half-plane branch. The negative control rejects a false1/10 lower bound.

## 2. Transfer to the actual complex window

At fixed u write D=z/2. The exact derivative in z of the Hermite matched phase is

d/dz [F_0(-z/2+u,z)-Phi(z,0)]
 =(u+D)/2+1/[2*(D+u)].

On(F1), |D|<1051/1000 and |D+u|>=18999/20000-47/50=199/20000. Its modulus is below52. Moving from z=iy to the allowed real part therefore changes W_0 by less than52/100000.

At the same physical boundary |v|<2, and |z+v|>=199/20000. Direct Taylor differentiation bounds the delta derivative of the first log term by

2^3/[6*(1-|delta|)^2],

and the second by 1/[199/20000-|delta|]. Their sum is below102. At the selected physical stationary points, |v_j|<3/2 and |z+v_j+delta*v_j^2/4|>3998/10000, so the matched mean has delta derivative below2, exactly as in the matched-data proof. Thus the actual normalized phase changes by less than110*|delta|. Equations(F4) and these two bounds give, on |u|=47/50,

|W_delta|=|F_delta-Phi|
 >9/100-52/100000-110/100000000
 >447/5000=0.0894.                                       (F5)

The data lemma gives |zeta|<1/8. For every |t|<=91/200,

|t^3/3-zeta*t|<(91/200)^3/3+(1/8)*(91/200)
               <447/5000.                               (F6)

The phase W_delta has exactly three physical zeros. At the base delta0,z2i it is S(u)=i*u^3*integral_0^1 s^2/(1-i*u*s) ds, whose integral factor has positive real part on |u|<1 and hence no zeros. It has a triple zero only at u0. The boundary(F5) never vanishes on the connected parameter domain, so the argument principle preserves that count. Rouche using(F5)--(F6) then gives exactly THREE physical v roots of(F2), counting multiplicity, for each t in the stated closed cubic disk.

The phase has exactly TWO physical critical points. Its stationary numerator P_delta is the exact cubic from the data lemma. On |u|=47/50, P_0=u^2-(1+z^2/4), with |1+z^2/4|<11/100. The delta terms are bounded by4*|delta|+|delta|^2, much smaller than (47/50)^2-11/100. Rouche gives two critical roots. Both are the selected pair; the third cubic root lies outside the physical disk. The logarithmic denominators remain nonzero throughout this disk.

## 3. A regular-sheet criterion with parameters

Here is the covering argument in full. Suppose first zeta!=0 and put b^2=zeta. Regard(F2) as a cubic cover t(u) of the fixed physical u disk. At either critical cubic value the equation in u has exactly three roots, by(F5); the target critical value has modulus at most(2/3)*(1/8)^(3/2)<3/100. One physical root is the corresponding selected stationary point of multiplicity two, and the remaining one is simple.

At the selected physical critical, the two quadratic factors cancel locally: phi(t)-phi(b)=(t-b)^2*(t+2b)/3 and W(u)-W(u_c)=(u-u_c)^2*a(u), with a nonzero analytic factor. Taking the local square roots gives two analytic crossing branches. In the normalization this is a NODE, not a branch point. At the unmatched simple physical root the cubic t cover has an ordinary square-root branch point. Thus there are precisely two simple branch points of the normalized cubic cover, one for each critical value. They are distinct while zeta!=0.

The boundary monodromy is the identity. At the exact merger all three cube-root branches of W=i*u^3 times a nonzero unit are single-valued around the boundary. For every parameter in(F1) and the delta disk, boundary t roots stay distinct: the phase lower bound(F5) is larger than both critical cubic values. Continuous parameter transport preserves the boundary permutation.

The boundary permutation is the product of the two interior transpositions. Consequently those transpositions are equal, and their common fixed letter gives exactly one globally holomorphic sheet T(u,z,delta). This does not yet assume injectivity or stationary matching.

For completeness this sheet is jointly holomorphic, not merely pointwise. Fix a physical boundary base point; its three distinct roots can be labelled holomorphically throughout the simply connected parameter rectangle times the delta disk. Off the curve zeta0, the unique globally unramified sheet fixes one boundary label, constant on the connected punctured parameter domain. At ordinary interior points the implicit-function theorem applies. At a selected node use the joint local factorization above, with its nonzero analytic unit and locally analytic b; the sheet chosen from the boundary selects one analytic node factor. Thus its branches agree jointly in u and parameters.

The data bound |zeta_z+i|<1/4 makes zeta a valid transverse parameter. Its unique zero curve lies within |z-2i|<10^(-5). Indeed on that circle the data lemma gives |zeta(z,delta)+i*(z-2i)| <= (49/100)/99999+(1/4)*10^(-5)<10^(-5); Rouche gives exactly one zero, and the derivative bound makes it simple. This punctured parameter domain is connected: any punctured z rectangle can be joined to the fixed safe point z=2.05i, which remains away from that small zero curve while delta is moved to0.

All roots of the monic cubic in t are uniformly bounded on the compact physical/parameter closure, since both logarithms have a uniform positive distance from their poles. Therefore the selected holomorphic sheet is bounded as zeta tends to0. Its Cauchy integral in the transverse zeta variable extends it removably across that curve, jointly in u and the other parameters.

At a merger W has one triple zero and no others in the physical disk. Any extended cube-root branch has one simple zero. The selected T therefore has zero count one there. Boundary(F5)--(F6) imply |T|>91/200 on the physical boundary, so the argument principle preserves this zero count over the full parameter domain.

For every |tau|<91/200, Rouche applied to T(u)-tau gives exactly one physical root, counting multiplicity. It is simple. This proves the full holomorphic regular inverse on the cubic disk and its bound |u|<47/50.

## 4. Matching and actual Jacobian estimates

Both canonical stationary values t=+/-b lie strictly inside the inverse disk, since |b|<sqrt(1/8)<91/200. The derivative of(F2) along the inverse is

F_v*v_t=t^2-zeta.

The inverse is regular, so v_t is nonzero. At t=+/-b its physical image is therefore a physical critical point. The two critical values are distinct when zeta!=0, hence each image is its corresponding selected physical critical. This proves matching without assuming that the unmatched values +/-2b lie in the inverse disk; they need not. The normalization v_t=i follows by continuity of the matched ratio q/b at the Hermite merger.

The closed cubic disk |t|<=9/20 lies strictly inside radius91/200. Cauchy on circles of radius1/20 around |t|<=2/5, using |u|<47/50, gives

|v_t|<= (47/50)/(1/20)=94/5<19.

Interpolate h=v_t at t=+/-b, h=A+B*t+(t^2-zeta)*G. The exact Cauchy divided-amplitude formula on |tau|=2/5 gives, for |t|<=3/8,

|G_t|<= (2/5)*19 /
 [((2/5)-(3/8))^2*((2/5)^2-1/8)]
 <350000.

These are actual inverse/Jacobian bounds over the full thin endpoint window. To certify a numerical a_0, one must still combine them with quantitatively compatible Airy contours, joint-norm constants, and complete nonlocal connecting-arc bounds. Those further estimates are not proved by the boundary certificate or by this inverse lemma alone.

## Exact reproduction

    .venv/bin/python research/factorial/image-rigidity/verify_gap_two_full_endpoint_uniform_inverse.py \
      --repo-root . --output /tmp/gap_two_full_endpoint_inverse_review.json

The verifier checks all eleven mandatory source records, every exact
closed boundary leaf and a fresh replay of the immutable boundary source.
The [matched-data theorem](QUANTITATIVE_GAP_TWO_ENDPOINT_DATA.md) and
[critical basepoint lemma](HERMITE_CRITICAL_CUBIC_INVERSE_DISK.md) retain
their separate scopes. The source applicability audit records which
analytic steps are proved directly and which finite enclosures support them.

Updated 20 November 2026.
