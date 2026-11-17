# Two constant moment relations reduce degree thirty to cubic balance

PF2RULED184. 17 November 2026. Actual complex polynomial plane Keller maps,
in the complete marked quadratic obstruction class.

Every remaining degree-D=30 survivor has intermediate degree EIGHT
and its compulsory critical quartic has

    M=constant*u^2, h=constant*(u-c)^14, c!=0,
    L=7, ord_c e=3, ord_c d=10,
    b_field=[C(u,w):C(P,Q)]=8,
    ord_v q=10, ord_v J>=7, ord_v Delta=14.              (1)

The original generic map degree in this conditional D=30 class is
sixteen. This does NOT exclude that remaining cubic-balance class,
all of D=30, or full ruling six.

The new ingredient is the FULL constant companion resonance in two
coefficient residues. Triple/simple, two-double, and even-quartic
escape limits have exact moment growth for EVERY positive odd index.
The highest term of either a scalar or a two-component constant
relation cannot cancel. No finite moment prefix is promoted to
an arbitrary-companion theorem.

## 1. Two relations with the same constant coefficients

Retain the original source, full normalizations and regular square
form of PF171--PF183:

    P=H3^2+d*z+e*w+f,
    w^2=N+M*z, dP wedge dQ=du wedge dw/M.                (2)

PF176 constructs constants c_n, with k>=3 odd and c_k!=0, such that

    Q-sum_(n=0)^k c_n*P^(n/2)=O(w^-3).                  (3)

The square root has the leading term of H3. Define COEFFICIENT
residues, using the same branch throughout,

    m_n=[w^-1]P^(n/2),
    l_n=[w^-2]P^(n/2)=[w^-1](w*P^(n/2)).                (4)

The polynomial Q, constants, and even half-powers have both
coefficients zero. Thus (3) gives TWO identities with the SAME
constant coefficients:

    sum_(1<=n<=k,n odd) c_n*m_n=0,
    sum_(1<=n<=k,n odd) c_n*l_n=0.                       (5)

Only POSITIVE odd n are used below. No negative-character identity
or comparison with an unaccounted finite contribution is used.

Fix the actual A support and put v=u-c. For any p-independent local
Laurent center a(v), define

    lhat_n=l_n-a(v)*m_n
           =[w^-1]((w-a(v))*P^(n/2)).                   (6)

It too satisfies the constant relation (5). Multiplying either
identity by a scalar Laurent series in v is allowed. This is
an operation on local identities; a(v) is not put in the target.

## 2. The complete degree-thirty local data

PF182--183 give exactly one A support, with

    L=7, b_field in {6,7,8}, all FOUR quartic sheets poles,
    ord_c e in {2,3,4,5}, t=7+ord_c e in {9,10,11,12}.    (7)

The two other generic roots are finite. With Y=v^7*w the entire
equation and separated depressed quartic are

    F0(v,Y)-p*v^14=F2*G,
    G=X^4+2J*X^2+q*X+J^2+Delta,
    ord_v J=ell>=1 or infinity,
    ord_v q=t, t+1<=delta=ord_v Delta<=14.               (8)

The depressed translation is p-independent below order fourteen.
All coefficients below that order are complex constants as functions
of v; the first p coefficient of Delta at order fourteen is nonzero
and affine. On all four geometric sheets,

    8<mu=ord_v G_X<21/2,
    sum_(four sheets) mu=32+b_field in {38,39,40}.        (9)

These statements keep ramified points with their sheet multiplicity.
The companion and all coefficient degrees are unrestricted.

For the moment calculations set p=0 in the factorization (8), since
the moments concern P itself. Every coefficient below order fourteen
is unchanged. In particular all leading data used in Sections 4--6
are unaffected. The generic parameter in the Newton and curvature
arguments remains p; it is not set to zero there.

## 3. A contour decomposition with the finite cluster retained

The coefficients of the actual polynomial P(v,w) are holomorphic
near v=0 because the support is off M. Its specialization has
degree TWO in w. On a fixed sufficiently large finite circle,
P(v,w) converges uniformly to a nonconstant quadratic, with no
zero on that circle. Exactly two roots remain inside; the other
four escape near w=a(v), where a(v) has pole order seven.

For v in a small sector, choose an escape circle after the indicated
Newton scaling, enclosing those four roots. Both clusters have EVEN
root counts, two and four. The branch of sqrt(P) selected at w
infinity is consequently single-valued outside the two cluster
discs. Cauchy's theorem splits its large-circle coefficient
integral into the two cluster integrals.

For every fixed positive odd n, the finite circle contributes

    O(1) to m_n, O(v^-7) to lhat_n.                      (10)

Indeed P is uniformly bounded there, and w-a(v)=O(v^-7).
The estimates do not omit the finite roots or assume their
contribution zero.

On the escape circle put

    w=a(v)+v^(a-7)*Z.

After clearing a finite fractional power of v, the rescaled factors
are convergent power series. The other quadratic factor is a
nonvanishing analytic unit U(v,Z), with constant nonzero reduction
U(0,Z)=tau^2, tau in C*. It has a holomorphic square root on the
whole fixed Z disc. This follows either from the analytic
implicit factor recursion or directly from the separation of
the two root clusters.

Thus a uniform leading expansion of the escape integral gives
the actual leading coefficient of a rational moment whenever
its exponent is strictly below the finite bounds (10).
Each residue below is a large Z-circle coefficient residue.
All signs are fixed by the original sqrt(P); a different choice
multiplies the displayed nonzero constants by a nonzero sign.

## 4. Triple and simple roots force strictly growing scalar moments

Suppose the first escape scale is X=v^3*Z, and the leading monic
quartic is

    R0(Z)=(Z-a0)^3*(Z+3a0), a0 in C*.                   (11)

Then P has leading v^-2*tau^2*R0(Z), and
dw=v^-4*dZ. Consequently

    m_n=tau^n*v^(-n-4)
                  *[Z^-1]R0(Z)^(n/2)+o(v^(-n-4)).       (12)

The finite contribution in (10) is strictly later.
For EVERY positive odd n the coefficient is exactly

    [Z^-1]R0(Z)^(n/2)
       =(4a0)^(2n+1)*binom(n/2,2n+1)!=0.                (13)

To verify it, put Z=a0+4a0*T. Then
R0=(4a0)^4*T^3*(T+1), and coefficient residue transforms
as a differential. The required coefficient is that of T^-1
in T^(2n)*(1+T^-1)^(n/2), giving (13).
A positive half-integer binomial coefficient at an integer
index is never zero.

Thus m_1,m_3,m_5,... have exact pole orders 5,7,9,... .
The top nonzero character c_k*m_k in (5) has a strictly larger
pole than every smaller character. This is impossible.
Every leading triple/simple face (11) is therefore excluded.

## 5. Two double roots are excluded by a nonvanishing moment vector

Suppose the first escape scale is X=v^3*Z, and its full normalized
quartic has the form

    (Z^2+j(v))^2+v*(A*Z+B)+O(v^2),
    j(0)=j0!=0, (A,B)!=(0,0).                           (14)

The O(v^2) remainder is polynomial in Z of degree at most one;
the complete square base j(v) is retained, not truncated.
Again P has the prefactor v^-2 and dw=v^-4*dZ.

Write n=2r-1, r>=1. Expand on a fixed escape circle using the
complete quadratic H=Z^2+j(v). In

    U(v,Z)^(n/2)*H^n
       *sum_j binom(n/2,j)
           [v*(A*Z+B+O(v))/H^2]^j,                     (15)

every term j<r is holomorphic inside that circle, even at the
zeros of H: the H exponent n-2j is nonnegative. The unit
U has a holomorphic square root on that entire disc.
These lower terms therefore have integral ZERO, rather than
an uncontrolled small contribution. The first possible term
is j=r; all later ones are O(v^(r+1)) on the contour.

Choose z0 with z0^2=-j0, and put

    U0=B+A*z0, V0=B-A*z0.

The leading two coefficient residues are

    M_r=[Z^-1](AZ+B)^r/(Z^2+j0)
         =(U0^r-V0^r)/(2z0),

    N_r=[Z^-1]Z*(AZ+B)^r/(Z^2+j0)
         =(U0^r+V0^r)/2.                                (16)

They cannot BOTH vanish, for any r>=1, because U0,V0 are
not both zero. This includes cancellation at either square
center and every root-of-unity ratio U0/V0.

Equations (10),(15),(16) give the vector asymptotic

    (m_n, v^4*lhat_n)
      =tau^n*binom(n/2,r)*v^(-r-3)*(M_r,N_r)
                                       +o(v^(-r-3)).    (17)

The finite contributions have orders O(1) and O(v^-3),
respectively. Since r+3>=4, they cannot affect (17).
The binomial factor is nonzero. Each successive odd index
has a strictly larger vector pole order, and the leading
vector is never zero. Both identities (5) cannot hold with
c_k!=0. Thus every case (14) is excluded.

## 6. A dominant constant quartic is excluded by the centered second moment

Suppose the escape scale is X=v^(13/4)*Z and the leading
quartic is

    Z^4+d0, d0 in C*.                                   (18)

Then P has prefactor v^-1 and dw=v^(-15/4)*dZ.
The centered second coefficient has leading term

    lhat_n=tau^n*v^(-(n+15)/2)
                     *[Z^-1]Z*(Z^4+d0)^(n/2)
                              +o(v^(-(n+15)/2)).        (19)

For every positive odd n,

    [Z^-1]Z*(Z^4+d0)^(n/2)
       =binom(n/2,(n+1)/2)*d0^((n+1)/2)!=0.             (20)

The finite contribution is only O(v^-7), while the first
displayed pole order is eight. Thus lhat_1,lhat_3,... have
strictly growing exact pole orders 8,9,... . Their constant
relation is impossible. This excludes (18).

All three moment arguments used positive indices only and
explicitly retained the finite cluster. No conclusion about
a negative half-power or a third centered moment is inferred.

## 7. Exhaustion of the remaining critical Newton regions

We now show that (11),(14),(18), elementary low-pole budgets,
and the actual affine response exhaust all cases except (1).
Use the complete arbitrary-coefficient Newton mechanisms of PF180.

### Dominant Delta and the square/Delta boundary

If delta<2ell and delta<4t/3, all four mu=3delta/4.
Equation (9) gives

    b_field=3delta-32 in {6,7,8}.

The only integral possibility is delta=13, b_field=7.
Its first quartic is (18), already excluded.

If delta=2ell<4t/3, at least two mu=3delta/4. The strict
interval (9) and even delta force ell=6,delta=12.
Its range is t=10 or11. Four simple initial roots give
b_field=4, too small. The sole repeated initial is a double
zero with two outer simple roots, each mu=9.
On the inner double cluster the whole second derivative has
order six, the third at least three, and higher derivatives
have nonnegative order. The first response index is
R=14-mu>7/2>3. The exact Taylor argument of PF181 therefore
forces inner mu<7+6/2=10. The total pole degree is then
STRICTLY below 2*(9-8)+2*(10-8)=6. This is impossible.

### Cubic/linear balance

If ell>2t/3 and delta>4t/3, all four mu=t. Hence

    b_field=4t-32 in {6,7,8}

forces t=10 and b_field=8. Its integer ranges are exactly

    ell>=7 or infinity, delta=14.                       (21)

These are the cases left in the theorem.

### Square centers

The region is ell<2t/3, delta>2ell.
Put m0=min(delta,t+ell/2). There are at least two
sheets with mu0=(m0+ell)/2. If delta!=t+ell/2,
all four have that value. If equality allows one-sign
cancellation, the guaranteed pair remains and the other
two obey mu<21/2. Thus in the equality case

    b_field<2mu0-11.                                    (22)

In particular mu0<=17/2 cannot reach b_field>=6.

Solving these bounded INTEGER order conditions in
t=9,10,11,12 and t+1<=delta<=14 gives exactly:

| Type | Possible (t,ell,delta) not excluded by degree/strict bounds |
|---|---|
| m0=delta<t+ell/2 | (12,5,14), (11,6,13), (12,6,13), (12,6,14) |
| m0=t+ell/2<delta | (10,6,14) |
| delta=t+ell/2, with possible one-sign cancellation | (10,6,13), (11,6,14), (12,4,14) |

For completeness, in the first row the common degree is
2(delta+ell)-32, hence six or eight; this gives
delta+ell=19 or20. In the second it is
2t+3ell-32, which must be six, seven or eight.
In the equality row ell is even and (22) bounds the pair.
These explicit equations, with the proved finite order ranges,
establish the table without a coefficient cutoff.

Every listed delta=14 case with m0=14 is excluded by a
quadratic first response. At a p-independent square center,
the whole F0 has constant order at least fourteen, first
derivative order at least t, exact second derivative order
ell, third derivative order at least ell/2 and higher
derivatives nonnegative. Put R=7-ell/2. Since
t+ell/2>=14 and ell<7,

    ell+2R=14,
    t+R>=21-ell>14,
    ell/2+3R=21-ell>14, 4R>14.

The whole initial is a nonzero quadratic in the first
parameter coefficient plus a constant minus p.
This contradicts PF159 at the compulsory actual pole,
exactly as in PF182. The nonzero affine-p coefficient at
order fourteen prevents a one-sign cancellation from
removing this response.

The other listed cells are

    (10,6,13), (10,6,14), (11,6,13), (12,6,13).

At X=v^3*Z their first non-square deformation is precisely
(14): A!=0 when t=10, B!=0 when delta=13.
At least one is nonzero in every cell. They are excluded
by the moment vector (17).

### The remaining boundary at t=9

The only additional boundary, not already in the preceding
regions, has t=9, ell>=6, delta>=12. The t=12 equality
ell=8 or delta=16 is already dominant Delta or outside
the parameter cap. For t=10,11 the corresponding equality
orders are nonintegral.

At X=v^3*Z, the leading quartic has the form

    H(Z)=Z^4+2j0*Z^2+q0*Z+j0^2+d0, q0!=0,             (23)

where j0 or d0 can be zero. Its cubic coefficient vanishes.
Four simple roots all have mu=9 and degree four.
A double root with two simple roots has two outer mu=9;
the two inner roots have curvature order six and thus
mu<10 by the same full Taylor argument above. Its total
pole degree is strictly below six.
Two double roots or a quadruple root are impossible with
q0!=0: a depressed quartic with all roots multiple has
zero linear coefficient.

The remaining multiplicity is triple plus simple. Depression
makes (23) equal (Z-a0)^3*(Z+3a0), a0!=0.
This is excluded by (12),(13). Thus the entire boundary
is excluded, including delta=13 or14 with ell=6 and
delta=12 with ell>6. No cell from the enlarged cap is lost.

All Newton regions and their exceptional repeated initials
have now been retained. Only the cubic/linear cases (21) remain.

## 8. Exact conclusion and remaining obstruction

The cubic/linear range gives b_field=8, t=10 and (21).
Since t=7+ord_c e and ord_c d=7+ord_c e, this gives
ord_c e=3 and ord_c d=10. PF182 already fixed M,h,L,
proving (1).

The formal degree-eight control retained in PF183 belongs to
this remaining class. Its leading quartic Z^4+q0*Z has
residues vanishing on entire congruence subsequences. The
strict scalar/vector growth arguments above do not close it.
Neither actual companion existence nor full constant resonance
for that control is asserted.

The critical D=30 cubic-balance class, larger leading degrees,
target trace perfection and HC4/JC2 remain open.

## 9. Exact evidence and replay

The checker verifies the two coefficient identities on exact
formal controls, every all-index binomial and two-vector formula,
whole sextic moments with a finite simple-root cluster, the
complete integer Newton table, and the parameter-order/curvature
inequalities. Its finite examples validate formulas; the written
contour, binomial and case arguments prove the unrestricted result.

Replay:

    .venv/bin/python research/jacobian/affine-complements/verify_ruling_six_degree30_two_moment_reduction.py

Essential proof sources are PF159, PF164, PF176, PF178, PF180,
PF181, PF182 and PF183 in this directory. Internal derivation
and replay do not constitute formal verification or external
human review. No full ruling-six or JC2 closure is claimed.
