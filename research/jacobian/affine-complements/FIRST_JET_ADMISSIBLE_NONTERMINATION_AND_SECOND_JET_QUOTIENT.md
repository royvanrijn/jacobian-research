# First-jet-admissible HC4 nontermination and the second-jet quotient

PF2RULED154. 17 November 2026.

The nilpotence certificate in [PF153](HC4_TRACE_QUOTIENT_AND_UNIFORM_CONE_NILPOTENCE.md)
is valid, but its first-jet quotient A_plane is too large for ANY uniform
nilpotent-ideal termination certificate. There are exact characteristic-zero
points of A_plane with cone-deviation ideal J=0 whose unique formal graph is
nonpolynomial. The same points admit the square-zero unit-propagation control
of the separate Rodrigues certificate, now INSIDE
A_plane/J. These are controls for the proposed transfer, not a counterexample
on the full surviving polynomial-graph locus.

An explicit next quotient A_2 is defined in Section 5. Its second-jet degree
equations remove this entire control family by the nonzero scalar -27/4.
Thus the reduced polynomial-body condition on the whole A_plane is disproved,
not merely unverified. A successful certificate must pass to further genuine
candidate equations before attempting that condition.

## 1. A polynomial target with admissible first jet

Let L,N,k be units in a commutative Q-algebra; M,tau are arbitrary. Define

    b(y)=9Ny^3+k,
    h_tau(y)=-27N^2 y^5/(2L)+tau y^2,

    F_tau=Lp^2+(M+4Ny^2)p+Nyq
          +x(b(y)p+h_tau(y))
          +x^2(b(y)^2-k^2)/(4L).                 (1)

This target potential has total degree eight over a field. Directly, since
q appears as Nyq, its Hessian determinant is N^2 times
(F_xp^2-F_xx F_pp). Here

    F_pp=2L, F_xp=b, F_xx=(b^2-k^2)/(2L),
    det Hess(F_tau)=N^2 k^2.                     (2)

This constant is independent of tau.

In fact grad(F_tau) is a polynomial automorphism over the unit-parameter
base, including nonreduced base changes. For its output (U,V,W,Z), recover

    y=Z/N,
    x=(b(y)(W-M-4Ny^2)+2L h_tau(y)-2LU)/k^2,
    p=(W-M-4Ny^2-b(y)x)/(2L),
    q=(V-8Nyp-x(b'(y)p+h_tau'(y))
                  -x^2 b(y)b'(y)/(2L))/N.        (3)

Consequently this control does not construct an HC4 or JC2 counterexample.

Retain the actual Meng--Yang parent

    u=1+xy,
    A0=u^3p+3xu^2q,
    B0=y^2u(4+3xy)p+(y+3xy^2(4+3xy))q,
    Psi0=LA0^2+MA0+NB0,
    C=N(2-3xy)-Mx^2-2Lx^2 A0,
    D_tau=(Psi0-F_tau)/x.

The numerator in D_tau vanishes at x=0, so D_tau is polynomial. The exact
graph equation Psi(x,y,p,q,R_tau)=F_tau is

    Lx^5 R_tau^2+C R_tau+D_tau=0.                (4)

Its implicit derivative at x=0 is 2N, a unit, so it has a unique formal
solution in A[y,p,q][[x]]. By the retained all-normal theorem HC4MYGJ2,
this is also the unique constant-determinant graph for its plane trace and
constant c=N^2k^2.

The first two normal coefficients R_tau=T_tau+x S_tau+... are

    T_tau= -27N y^5/(4L)+y^3p-6y^2q
       -(3L/N)(yp^2+pq)-(3M/(2N))(yp+q)
       +(k/(2N))p+(tau/(2N))y^2,                 (5)

    S_tau= -(12L/N)y^2p^2-(39L/(2N))ypq
       -(9L/(2N))q^2-(15M/(4N))y^2p
       -(21M/(4N))yq-(27/2)y^3q
       +(3k/(4N))yp+(3tau/(4N))y^3
       +(9k/(4L))y^3.                            (6)

They have plane degrees five and four. More precisely,

    (T_tau)_5=-27N y^5/(4L), (T_tau)_4=y^3p,
    J maps to zero.                             (7)

The exact plane identity F_plane(T_tau)-8LN^3 S_tau=N^2k^2 proves that
(5) defines a homomorphism from PF153's ACTUAL A_plane to the unit-parameter
ring. It factors through A_plane/J. This is not a point of a radical or
initial substitute, and it is not the old first-jet-degree-six control.

## 2. All-order nontermination at a reduced point

Over a characteristic-zero field, put

    Delta_tau=C^2-4Lx^5 D_tau.

Exact calculation gives

    degree_q Delta_tau=1,
    partial_q Delta_tau=-12LNx^3(xy+2).           (8)

Thus Delta_tau has valuation one at its linear q-prime over K(x,y,p).
It is not a square in K(x,y,p,q). Equation (4) has no rational solution,
although its unique x-adic solution exists. This proves nonpolynomiality
without extending a finite coefficient prefix.

For an explicit all-order witness, set y=p=0 after extracting the graph:

    Delta_tau|_(y=p=0)=(2N-Mx^2)^2-24LNx^3q.

Choose sqrt(Delta_tau)|_(x=0)=2N. From

    R_tau=(-C+sqrt(Delta_tau))/(2Lx^5),

for every integer j>=2 one obtains

    [x^(3j-5)y^0p^0q^j] R_tau
       =(N/L) binomial(1/2,j)(-6L/N)^j !=0.       (9)

The first q^j term of the square root has x order 3j; C has q degree at
most one. This proves the displayed identity for all j. The nonzero
coefficient is a unit times a nonzero rational number. In particular it
persists at L=1,M=13,N=2,k=1,tau=0, an ordinary rational point of
A_plane/J with determinant four.

Every nilpotent ideal of A_plane maps to zero at this field point.
Therefore there do not exist ANY nilpotent ideal K of the WHOLE A_plane
and finite integers D,c_0,N_0 with

    K^N_0=0,
    a_(D+m) in K^ceil(m/c_0) A_plane[y,p,q]
                for all m>=1.                    (10)

for its universal normal coefficients. Such containment would specialize
to a polynomial graph, contradicting (9). This conclusion allows an
arbitrary normal shift D; changing the proposed cone ideal cannot rescue
a whole-first-jet-quotient certificate. It does not apply to a further
quotient that excludes this field point.

## 3. The square-zero propagation control now factors through A_plane

Differentiating (4) in tau gives

    partial_tau R_tau=y^2/sqrt(Delta_tau),
    partial_tau Delta_tau=4Lx^5 y^2.              (11)

Use the ACTUAL Image Rodrigues length17 algebra H and eta=b_H^6, with
eta!=0 and eta^2=0, as retained in PF152. The symbol b_H is its algebra
generator, distinct from the polynomial b(y) in (1). Base change tau=eta
is valid by (5)--(7), so the full first-jet equations hold in H and the
cone ideal still maps to zero. The lifted graph difference is EXACTLY

    V=R_eta-R_0=eta*y^2/sqrt(Delta_0).             (12)

At L=1,M=13,N=2,k=1,

    [y^2p^0q^0]V=eta/(4-13x^2).                 (13)

Its projected transition is multiplication by the unit 13/4. Every late
coefficient ideal of V in H is (eta), which is nonzero and square zero;
there is no increasing-power bound in (eta), or in H's maximal ideal.
For M=0 the first q^j term still gives

    [x^(3j)y^2p^0q^j]V
       =eta/(2N) binomial(-1/2,j)(-6L/N)^j !=0.   (14)

There is also a minimal ACTUAL FINITE TRACE QUOTIENT here. At
L=1,M=13,N=2,k=1,c=4, fix the 55 trace coefficients other than [y^2]T
to their values in (5) at tau=0. The one remaining coefficient is tau/4.
The linear fixing equations eliminate all other trace variables, and the
plane equations vanish identically in tau by (5)--(6). Therefore the
whole resulting quotient is exactly Q[tau], not just an image of it.
Adjoining the relation tau^2=0 gives the genuine quotient

    B=Q[tau]/(tau^2), dim_Q B=2,
    K=(tau), K!=0, K^2=0.                        (14a)

In basis (1,tau), actual multiplication by tau is [0 0;1 0]. Its unit
vector and this matrix retain the nonzero nilpotent direction. The
coefficient variation V over B satisfies (12)--(13) with tau in place of
eta, so its late-tail ideals are K and its projected transition is the
invertible matrix (13/4)Id_B. Thus even an actual finite HC4 first-jet
quotient with its actual multiplication does not force increasing powers.
Its reduced graph is nonpolynomial by (9), so it fails the separate
reduced-body requirement as well. These finite quotients are not claimed
to satisfy the later degree-five candidate equations.

The original Rodrigues multiplication is still used to certify eta and
the specialization tau=eta, without identifying A_plane with H or using
the Rodrigues initial table. Equations (11)--(14a) are exact all-order
identities and quotient identifications, not sampled normal prefixes.

## 4. This is a first-jet control, not a surviving degree-five graph

Let a_2 be the second normal coefficient of (4). Write C=sum C_i x^i,
D_tau=sum D_i x^i. The quadratic term Lx^5R_tau^2 has not entered yet, so

    a_2=-(C_1 S_tau+C_2 T_tau+D_2)/(2N).

For this ENTIRE parameter family,

    [y^5p]a_2=-27/4.                             (15)

The plane degree is six. A graph of total degree <=5 would instead require
degree(a_2)<=3. Thus the necessary second-jet equations already exclude
this control family, including all its square-zero deformations. No new
surviving coordinate family is excluded by this observation. The role of
(15) is to identify an actual further gate that removes the obstruction
to (10); it is not a new counterexample on that smaller quotient.

## 5. An explicit finite next quotient

There is no need to guess which equations should follow I_plane. In its
coefficient algebra, set S=(F_plane(T)-c)/(8LN^3). Expand the ACTUAL parent
potential on the truncated graph R=T+xS as

    Psi(T+xS)=f_0+x f_1+x^2 f_2+x^3 f_3+O(x^4).

For source variables v=(y,p,q), define 4 by 4 matrices

    H0_xx=2 f_2, H0_xi=partial_i f_1,
    H0_ij=partial_i partial_j f_0,
    H1_xx=6 f_3, H1_xi=2 partial_i f_2,
    H1_ij=partial_i partial_j f_1.                 (16)

The forced second coefficient is the finite polynomial recipe

    B_2(T)=trace(adj(H0) H1)/(24LN^3),
    I_2=I_plane+(all source-degree>3 coefficients of B_2(T)),
    A_2=P/I_2.                                  (17)

To prove (17), the coefficient of x in det Hess(Psi(T+xS)) is
trace(adj(H0)H1). Adding x^2U to R changes that coefficient by
-24LN^3 U, precisely HC4MYGJ2 with k=2. Terms of normal order >=3 are
invisible here. Thus (17) gives exactly the next forced jet degree bound.
All derivatives in (16) are source derivatives of polynomials, not
parameter derivatives. The recipe is finite and its denominators are
units; using untruncated S is harmless modulo I_plane.

PF153's J^4096=0 descends to A_2. Equation (15) shows that the control map
to the parameter ring cannot factor through A_2: its second-jet ideal
maps onto the unit ideal. This supplies a specified next trace quotient,
without claiming that A_2 is the full surviving locus, is nonempty, or
has a polynomial reduced formal graph. Other retained normal-form and
higher-jet equations still apply.

The larger candidate quotient still lacks a uniform increasing-power
certificate and a positive termination theorem. The
[complete quintic residual proof](FULL_COMPLEX_QUINTIC_GRAPH_EXCLUSION_AND_CUBIC_NILPOTENCE.md)
closes the bounded degree-five graph locus after its additional necessary
constraints; it does not identify this partial-jet ring with that locus.
The global [PF151 coverage obligations](NORMAL_FORM_COVERAGE_AND_COMMON_LEADING_DIVISORS.md)
remain unproved. The
marked polynomial-graph collision route stays closed separately. Full
finite normalization S, its codifferent and every ramified boundary prime
are retained; no Artin parameter algebra is identified with S. The 3D
integrality controls and PF142 mixed-minimal retirement remain in force.
Full HC4/JC2 and an actual counterexample remain OPEN. Independent
whole-proof replay, formal verification and external review are unclaimed.
