# Height nine: complete gcd-one source-degree and resonant chart restrictions

14 November 2026. PF2RULED17 narrows the three gcd-one height-nine pairs
on the actual marked Wright model. It proves the complete necessary rows
below; it does not exclude a whole remaining pair or construct a finite
map. The gcd-three pair(6,9) is outside this theorem.

Let k be the contact multiplicity of the primitive common cubic tau.
Thus div(tau)=(3-k)*p+k*q for k=1,2, or3q for k=3. The away point p
is arbitrary, including the marked point r=infinity. All rows retain
all lower-pole lifts.

| Pole pair | Ordinary source degrees(P,Q) | Remaining k |
| --- | --- | --- |
| (5,9) | (48,128) | 2,3 |
| (7,9) | (48,128) | 3 |
| (7,9) | (80,128) | 2,3 |
| (7,9) | (96,128) | 1,2,3 |
| (8,9) | (48,128) | 3 |
| (8,9) | (80,128) | 2,3 |
| (8,9) | (96,128) | 2,3 |
| (8,9) | (112,128) | 1,2,3 |
| (7,9) | (72,108) | 2,3 |
| (8,9) | (72,108) | 2,3 |

The degree128 rows have chart-leading forms proportional to T^j,T^8,
where T=v*(u+v)^3 and deg(P)=16j. Their cusp multiplier has slope3/4.
The last two rows have the retained resonance R^6,R^9, with
R=v*(u+v)^2, and cusp multiplier slope2/3. The reversed ordinary
orientation(108,72) is excluded, including for the pole(8,9) pair.

## Complete spaces and the actual chart

Use the [actual marked inclusion](WRIGHT_MARKED_MODEL.md)

 A=C[Y,XY,X^2Y,X^3Y+X] inside C[X,Y],
 X=u^-1, Y=u^2*(u*v-1).

For W_(m,d)=V_m intersect source-degree<=d, reconstruct the full
[PF2COORD1 Cox basis](SOURCE_COORDINATES.md), with

 alpha=1+X^2Y, beta=-alpha^2-Y, t=1+X^2*alpha, s=X,
 alpha^(m-j)*beta^j*s^i*t^(m+j-i),
 0<=j<=m, 0<=i<=m+j.

The exact source echelon procedure is the one justified in
[PF2RULED15](CONTACT_SOURCE_FILTRATION.md). Its invertible rational
triangular basis changes retain complete source polynomials, fixed
boundary symbols and every lower-pole term. The symbol in the actual
g0=beta*t^2,r=s/t frame is(-1)^(m-j)*r^(2(m-j)+i).

The [checker](verify_height_nine_frontier.py) constructs all sections for
m=5,7,8,9, reconstructs every source row, checks regularity in u,v and
then performs exact linear elimination of the contact and homogeneous
chart constraints. A triangular basis change on the filtered input rows
certifies that the zero-constraint rows span the entire kernel, not a
selection of lifts. All transformations and dimensions are retained in
the [certificate](height_nine_frontier_certificates.json).

The all-height period and cusp results force contact symbols for every
gcd-one pair here: P and Q have symbols tau^m,tau^9, tau(q)=0. The
source bound is therefore

 deg(P)<=16(m-1), deg(Q)<=128.

Imposing only contact order at least the pole gives:

| Complete contact space | Dimension | Maximum chart total degree |
| --- | --- | --- |
| W_(5,64), symbol order>=5 | 46 | 16 |
| W_(7,96), symbol order>=7 | 85 | 24 |
| W_(8,112), symbol order>=8 | 109 | 28 |
| W_(9,128), symbol order>=9 | 136 | 32 |

The full unconditioned W_(9,72) has chart degree<=21. All these are
complete finite-space facts. No all-height extrapolation is made.

## The degree128 branch

The pinned2022 degree theorem leaves maximum ordinary degree>=125 or
the exceptional pair72,108. In the first branch Q must be the larger
coordinate, since the three P bounds are64,96,112. The complete source
echelon basis has no degree125,126 or127 row in W_(9,128). Hence Q has
ordinary degree exactly128.

On every row of the contact space for Q, the entire chart projection
of total degree>=32 is

 Q_[>=32]=c*v^8*(u+v)^24=c*T^8.                         (1)

The nonzero c is exactly the ordinary X^96Y^32 coefficient and the
v^32 restriction coefficient. Thus deg_(u,v)(Q)=32.

Let d be the nonzero chart degree of P. Its leading Jacobian must
vanish. In w=u+v, the operator for Q_32=v^8*w^24 is

 8*w*partial_w-24*v*partial_v.

On w^i*v^(d-i) its eigenvalue is32i-24d. Therefore d=4j and

 P_d=a*v^j*w^(3j)=a*T^j, a!=0.

The bounds above give1<=j<=4,6,7 respectively. Converting a chart
monomial u^a*v^b to the source has highest Laurent total degree4b-a.
Consequently P_d forces ordinary source degree16j. Its coefficient of
v^(4j) is nonzero, and no lower chart component can cancel it.

Furthermore the source degrees16j-1,16j-2,16j-3 are absent. Indeed,
a monomial u^a*v^(4j-a) in T^j has source degree16j-5a; for a>=1
this drops by at least5. Lower chart components drop the bound by at
least4. The same gap holds for Q, which leads as T^8.

If j=1,2 or4, subtract a scalar multiple of P^(8/j) from Q. This
cancels source degree128 and leaves source degree<=124, while P has
degree16,32 or64. The target shear preserves the Keller identity,
field degree four and target algebra. It may raise the pole height,
which is allowed here: the published degree alternative applies to
every such target presentation. A maximum below125 with a component
below72 satisfies neither alternative, a contradiction.

Thus j is3 for m5;3,5,6 for m7; and3,5,6,7 for m8. This proves the
ordinary-degree rows in the first part of the table.

## The72,108 exception and its genuine resonance

For the normal orientation P has ordinary degree72 and Q has degree108.
The pole-five pair cannot have that orientation because its P bound
is64. For the other two pairs the full P spaces have chart degree<=20.
The contact part of W_(9,108), dimension131, has the universal high
projection

 Q_[>=27]=a*u*v^7*(u+v)^20
          +v^7*(u+v)^18*(c*v^2+b*u*v+d*u^2),            (2)

where c is the nonzero X^81Y^27 coefficient. The coefficients a,b,d
are unrestricted in this linear projection. Dropping the quadratic
term would incorrectly discard the resonance below.

If a!=0, Q has degree28 with linear-factor multiplicities1,7,20.
If a=0, Q has degree27. Its highest form is primitive except when the
quadratic in(2) is c*v^2. To check this exhaustively, factor that
quadratic over C. Its remaining multiplicity lists are

 (7,18,1,1), (7,18,2), (8,18,1), (7,19,1),
 (8,19), (7,20), or (9,18).

Only the last list has gcd>1; it is precisely the c*v^2 case. For
homogeneous forms F_d,G_e with J(F_d,G_e)=0, the elementary relation
F_d^e=constant*G_e^d follows in the variable u/v. Factor valuations
then imply e divides d when the multiplicities of G_e are coprime.
Thus a primitive form of degree27 or28 cannot commute with a nonzero
P_d of degree<=20.

The sole surviving leading form is therefore

 Q_27=c*v^9*(u+v)^18=c*R^9.

Here the chart degree of P must be divisible by3. Ordinary source
degree72 requires chart degree at least18, since source degree is at
most four times chart degree. With the upper bound20, this forces

 deg_(u,v)(P)=18, P_18=a*v^6*(u+v)^12=a*R^6.

The two displayed forms do commute. They are a retained leading control,
not a finite realization or a solution of the full Jacobian equation.

The reversed orientation is possible numerically only for pole(8,9).
Its pole-eight/source108 coordinate has contact order>=8 and the
[PF2RULED16 primitive top](POLE_HEIGHT_EIGHT_EXCLUSION.md)
c*v^7*(u+v)^20. The other coordinate lies in the complete W_(9,72),
whose chart degree is<=21. That primitive degree27 form gives the
same nonresonance contradiction. This excludes(108,72) in full.

## Contact refinements on the surviving rows

For each necessary source/chart-leading family, impose the complete
linear constraints: source degree bound, contact order>=m, chart degree
bound and the exact leading form T^j or R^6. Exact elimination gives
the following minimum possible leading-symbol orders at q:

| Pole m | Source degree | Leading chart form | Minimum symbol order |
| --- | --- | --- | --- |
| 5 | 48 | T^3 | 9 |
| 7 | 48 | T^3 | 17 |
| 7 | 80 | T^5 | 9 |
| 7 | 96 | T^6 | 7 |
| 8 | 48 | T^3 | 21 |
| 8 | 80 | T^5 | 13 |
| 8 | 96 | T^6 | 9 |
| 8 | 112 | T^7 | 8 |
| 7 | 72 | R^6 | 11 |
| 8 | 72 | R^6 | 15 |

An actual symbol tau^m has order m*k. Comparing with these exact
lower bounds gives every contact list in the first table. The minima
are attained in the full linear spaces, but this certifies neither the
complete nonlinear common-symbol condition nor the full map equations.
No retained row is asserted realizable.

## The two cusp-multiplier slopes and the remaining problem

The [actual finite-quartic cusp identity PF2OMIT9](../jc2-omitted-divisors/WRIGHT_IMAGE_RESTRICTION.md) has h=v-v* and
h*F1=lambda*F0', with lambda(v*)!=0. The leading chart forms prove:

* Degree128 branch: boundary restriction degrees4j,32; lambda is linear
  with slope24/32=3/4.
* Normal72,108 branch: boundary restriction degrees18,27; lambda is
  linear with slope18/27=2/3.

The nonzero value at the cusp and all lower coefficients still need the
full Jacobian and global finite-map constraints. These slopes do not
supply such a map.

The [receipt](height_nine_frontier_receipt.json) seals the full exact
finite-space calculation. Section-space geometry, the retained published
2022 degree alternative, homogeneous-factor argument and source-gap
shear step are the written inputs identified above. No conditional2026
upgrade, independent full reproof, formal verification, external review
or novelty priority is claimed.

The first pairs remain(5,9),(6,9),(7,9),(8,9). The three gcd-one pairs
are narrowed, not excluded. The gcd3 pair(6,9) needs an algebraic-cover
period and its own full source/jet calculation; the rational gcd1/gcd2
criterion cannot be copied to it. Higher heights, all-height pole-four
coordinates, actual arbitrary-normalization compatible completion, both
full quartic packets and the original full goal remain open.
