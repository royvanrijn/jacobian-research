# Simple escape graphs close ruling-six faces B and C

PF2RULED165. 17 November 2026. Actual complex polynomial plane Keller maps.

Every ruling-six pole face of type B or C is excluded in the complete
quadratic source-image class. Each has another escaping branch with a
simple p-independent initial root and holomorphic coframe. Its entire
companion residue is polynomial in P, so its image is a forbidden
closed affine-line component of the nonproper-value set.

Consequently only quartic faces A remain. Combining their weighted
initials with the actual coupled leading coefficients forces every
reduced companion degree to be an odd multiple of three. Therefore
b3 vanishes at every M-root, and the central degree-three alternative
is closed. The remaining central degree-two case retains PF162's
quadratic residue extension and original ramified prime.

In fact every survivor has the global polynomial normal form
b3=M*h(u)^2. Its reduced companion s=3k, k odd, has top odd
coefficient proportional to M^((k-1)/2)*h^k.

[PF188](RULING_SIX_SQUARE_DEFORMATION_AND_MARKED_QUADRATIC_EXCLUSION.md) excludes the remaining A class within the retained marked quadratic
model. Trace integrality and JC2 remain open.

## 1. Actual setup and imported component theorem

Keep the complete quadratic source model and actual fields of PF162:

    w^2=N+M*z, N=u*Sigma^2, M=u*A,
    u=t^2, w=t*(Sigma+A*x), z=2Sigma*x+A*x^2,
    B=C[P,Q] subset O(V), F=C(P,Q), E=C(u,w), L=C(x,t),
    dP wedge dQ=du wedge dw/M.                              (1)

P has ruling degree six, with expression

    P=b3*z^3+w*(a2*z^2+a1*z+a0)+b2*z^2+b1*z+b0.             (2)

All coefficient degrees and all companion degrees are unrestricted.
Retain the ORIGINAL full finite normalizations S_E and S of B in
E and L. No target change below replaces either ring or deletes
a prime.

The imported global input is the exact COMPONENT Remark of Chau,
Two remarks on non-zero constant Jacobian polynomial maps of C2,
Ann. Polon. Math.82(2003), printed43 / full PDF page5: no irreducible
nonproper-value component of a nonsingular plane polynomial map is
isomorphic to A1. This is the retained published input and application
of [PF161, Sections 1,4,5](POLYNOMIAL_RESIDUE_RAMIFIED_LINE_AND_QUADRATIC_RULING_FIVE_CLOSURE.md).
It is not a prohibition on all singular rational curves.

For c off M, write v=u-c, P*=sum_(i=0)^6 p_i(v)w^i, with
p_i in C[[v]]. The actual companion Q* is also polynomial in w
over C[[v]], since only powers of M are inverted. On P=p,

    omega=dQ*=-du/(M*P*_w).                                (3)

Every target coordinate after Q->Q-f(P), f in C[P], has degree
at least two on the whole source line over c, by PF52's retained
line-injectivity consequence. Here w is affine-linear on each such
line. The map remains noninvertible because its field has the
nontrivial quadratic deck factor L/E.

## 2. A simple constant escape root with holomorphic coframe forces a graph

Suppose c is off M, h,H are positive integers, and substitution
w=v^-h*Y in the ENTIRE equation gives

    v^H*(P*(c+v,v^-h*Y)-p) in C[p][[v]][Y],
    its reduction I(Y), I(y0)=0, I'(y0)!=0, y0 in C*.       (4)

Assume H>=h. The derivative in (4) is a unit with nonzero CONSTANT
term I'(y0), independent of p. All-order implicit recursion gives

    Yhat in C[p][[v]], Yhat(0)=y0,
    what=v^-h*Yhat in C[p]((v)).                            (5)

Every new coefficient divides only by that same nonzero constant.
Thus no p-dependent denominator or deleted finite p-value occurs.

The exact differentiation identity in (4) gives

    ord_v(P*_w)=h-H, ord_v(what_p)=H-h,
    ord_v omega=H-h>=0.                                   (6)

Its leading coefficients are nonzero constants independent of p.
The actual Q* evaluates to a Laurent series in C[p]((v)).
Since its v derivative is the holomorphic (3), it has no negative
Laurent term. Consequently

    Qhat in C[p][[v]], q0(p)=Qhat(0) in C[p].               (7)

This produces an ORIGINAL divisorial valuation of E regular on B,
with residue C(p) and centre D over

    Gamma={Q=q0(P)}.

Indeed P has transcendental residue p and Q-q0(P) has positive
finite order, by (6). The full normalization S_E is contained in
the valuation ring. Its height-one centre has

    S_E/D=C[P], image(D)=Gamma, residue degree1.             (8)

The quotient is integral over C[P] and contained in C(P), so
integral closedness proves equality. In particular (8) is the
whole closed graph, not a punctured curve.

As c is off M, c!=0 and A(c)!=0. The two unit square-root choices
t^2=c+v lift (5) to the original source field, with coefficients
in C[p]((v)). At either choice x=(what/t-Sigma)/A has a pole of
order h. Thus the resulting original primes of S are outside
the whole plane source, have the same graph residue, and map
onto Gamma. This is an irreducible nonproper-value component,
forbidden by the imported Chau Remark.

One can also fix any p0 in C. The analytic implicit function theorem
applies to (4), with the same constant unit derivative, giving an
actual escaping source arc with P=p0 and Q tending to q0(p0).
It confirms directly that every finite point of the graph is
retained. Other primes above Gamma and every other target curve
remain in the original full normalization.

Therefore (4) with H>=h is IMPOSSIBLE for an actual target (2).

## 3. Every B face has a forbidden simple extra escape

PF162 covers every possible pole by faces A,B,C, and PF164 sharpens
the pole budgets. For B the weighted initial has the form

    I=C*Y^2*(Y-c0)^3*K(Y),
    H=2L, L>=5, deg K<=1, K(0)*K(c0)!=0.                  (9)

If deg K=1, its root d is nonzero and different from c0, hence
simple in I. Use w=v^-L*Y and H=2L in (4), at y0=d.
Then (6) gives a coframe zero of exact order L. Section2 excludes
this entire case, including the possible missing middle coefficient.

If deg K=0, the face has degree five. Put delta=ord_c p6>4L,
h=delta-3L>L. The nonzero face coefficients have orders

    ord p2=0, ord p3=L, ord p4=2L, ord p5=3L.

Now substitute w=v^-h*Y and scale by v^H with H=5h-3L.
The whole reduction is

    Y^5*(eta5+eta6*Y), eta5*eta6!=0.                       (10)

The nonzero root -eta5/eta6 is simple. All lower coefficients lie
strictly above this face: degrees2,3,4 have margins respectively
3(h-L),2(h-L),h-L; degrees0,1 have positive margins as well.
The coframe zero in (6) is H-h=4h-3L>0.
Section2 again gives the forbidden graph.

Thus NO actual ruling-six pole divisor has type B, with either
constant or linear K, arbitrary later coefficients retained.

## 4. Every C face has a forbidden simple LOWER escape

For C, PF162 gives

    p2(0)=gamma!=0, ord p3=B0,
    (ord p3,ord p4,ord p5,ord p6)
                  =(B0,B0+L,B0+2L,B0+3L), 5<=B0<L.       (11)

Use the LOWER slope h=B0, not the triple-root slope L. Substitute
w=v^-B0*Y and scale by v^(2B0). The ENTIRE reduction is

    Y^2*(gamma+eta3*Y), gamma*eta3!=0.                     (12)

The higher terms of degrees4,5,6 have positive margins
L-B0,2(L-B0),3(L-B0); degrees0,1 also have positive margins.
The nonzero root -gamma/eta3 is simple. Equation (6) gives an
exact coframe zero of order B0. The all-order graph argument
excludes the entire C class, including B0 divisible by three.

The root in (12) is the fourth escaping sheet left over from the
triple cluster. It is not an original boundary prime discarded
because a different weighted face was chosen. Including that
sheet is what produces the contradiction.

## 5. Independent weighted subduction and its proper-power boundary

There is a second exclusion of B which does not use the component
theorem. It also describes the exact surviving A alternative.

Give v weight -1 and w weight L. For a nonzero polynomial in w
over C[[v]], its weighted degree is max_i(iL-ord_v coefficient_i);
it is a finite INTEGER. Suppose

    in(P*)=v^-H I(v^L*w), in(Q*)=v^-K J(v^L*w).

For every coordinate Q-f(P), the whole-line criterion gives K>=2L.
For the B faces, H=2L; for C, H=3L-B0; for A, H=2L. In each
case H+K+1-L>0, whereas the actual Jacobian 1/M has weight zero.
Thus the highest Jacobian term vanishes:

    H*I*J'-K*I'*J=0.                                      (13)

At every root of I, equation (13) forces
H*ord(J)=K*ord(I); J has no other roots. Let g be the gcd of
the multiplicities of all roots of I. If gcd(H,g)=1, then H divides
K and J=C*I^(K/H). Subtract C*P^(K/H) from Q. This is a legitimate
polynomial target automorphism. It strictly lowers the integer
weight K while the whole-line bound K>=2L persists. Finite descent
is impossible, so that initial type is excluded.

Both B forms have g=gcd(2,3[,1])=1 and are excluded this way.
For C, g=3; the same argument would already exclude 3 not dividing
B0, but Section4 closes the remaining congruence too.

For A,

    H=2L, I=C*[Y*(Y-c0)^2]^2.

After any possible integer-power cancellations, (13) instead allows

    K=kL, J=C'*[Y*(Y-c0)^2]^k,
    k>=3 ODD.                                             (14)

An odd half-integral K/H cannot be removed by subtracting a
polynomial power of P. No false Euclidean step is made here.
Since deg I=6, each integer-power cancellation has w degree
6(K/H)<=deg_w Q*, so it never increases the companion ruling
degree. Equation (14) already forces companion ruling at least9.

## 6. Coupled leading coefficients force odd companions and central degree two

Take a globally reduced companion: s=deg_w Q*>6 and 6 does not
divide s. Its existence by elementary target reductions is retained
in PF162. Put d=gcd(6,s). The actual leading Jacobian relation is

    p6=alpha*H0^(6/d), q_s=beta*H0^(s/d), H0 in C(u)*.      (15)

Every actual P has a coframe pole. Sections3,4 force at least one
A face, where delta=ord_c p6=4L, H=2L. Relation (15) gives

    ord_c q_s=(s/6)*delta=2sL/3.

Hence the top w term of Q has weight sL/3, so K>=sL/3.
On the other hand (13), with the A root multiplicities, gives

    J=C'*[Y*(Y-c0)^2]^k, K=kL, deg J=3k<=s.

Therefore K<=sL/3. Equality holds throughout:

    s=3k, K=kL, and k is ODD because 6 does not divide s.   (16)

This proves that every globally reduced companion of a ruling-six
coordinate has odd degree divisible by three, and s>=9.

Write s=2n+1. Its unique odd leading coefficient has q_s=a_n/M^n.
Using p6=b3/M^3, equation (15) gives

    b3^(s/d)=C*M^(3/d)*a_n^(6/d).                         (17)

At EVERY M-root the right side vanishes, so b3 vanishes there as
well. In particular b3(0)=0. Thus the central restriction has degree
TWO, and the entire central degree-three alternative is excluded.

All remaining ruling-six coordinates therefore satisfy the complete
central conclusions of PF162: one original ramified prime with
quadratic residue C(lambda)/C(P), P=lambda^2; the nonconstant image
class T=lambda*H(P); its conductor and all original ramification.
Those singular rational images are not forbidden graphs and are
not excluded by this note.

## 7. Noncentral M-roots force the polynomial square normal form

Let c!=0 be an M-root, m_c=ord_c M>=2 and
delta=ord_c b3>=1, the last inequality supplied by (17).
Suppose delta<m_c. The reached values w0=+/-sqrt(N(c)) are nonzero,
and b2(c)+w0*a2(c) is a unit by the whole-line condition in PF162.

Put y=1/z. The exact equation is

    b3+y*(b2+w*a2)+y^2*(b1+w*a1)
                              +y^3*(b0+w*a0-p)=0.          (19)

Its y derivative at (v,w-w0,y)=0 is the unit b2(c)+w0*a2(c).
All-order implicit lifting over C[p] gives

    y=b3*U(v,w,p), U a unit,
    w^2=N+v^(m_c-delta)*V(v,w,p), V a unit.

The second equation has unit w derivative 2w0 at v=0, so it too
has a unique lift with coefficients in C[p][[v]]. The partial
p derivative of y at fixed w has order 3delta, since (19)
gives y_p=y^3/(unit). Differentiating the second equation, its
w derivative remains a unit and gives the EXACT orders

    ord_v(w_p)=m_c+delta,
    ord_v omega=delta>=1.                                 (20)

All leading coefficients here are nonzero constants independent
of p. Although Q* can have finite M-denominators at c, evaluation
still lies in C[p]((v)). As before (20) forces its residue to be
a polynomial q0(p) and its whole Laurent series to be regular.

All original t choices satisfy t^2=c+v. If w0/t=Sigma(c), then

    w/t-Sigma=(w^2-N)/(t*(w+t*Sigma))

has exact order m_c-delta; division by A of order m_c gives
an x pole of order delta. For the opposite sign the numerator
is a unit and x has a pole of order m_c. Thus every original
lift is omitted. The full-normalization closed graph argument
of Section2 applies and contradicts Chau's component theorem.
Consequently delta>=m_c at EVERY noncentral M-root.

At the central root PF162 already gives delta=m_0+2a>=m_0.
Together with (17), using s=3k, d=3, k odd, we have

    b3^k=C*M*a_n^2.

At each finite root this implies ord b3 congruent to ord M modulo
two. We have also proved ord b3>=ord M at every M-root.
Therefore b3/M is a POLYNOMIAL with all zero orders even.
Over C absorb its nonzero constant into a square, obtaining

    b3=M*h(u)^2, h in C[u] nonzero.                        (21)

Since p6=(h/M)^2, relation (15) and s=3k show
q_s=beta*(h/M)^k. For n=(3k-1)/2 this gives the exact
polynomial leading coefficient

    a_n=beta*M^((k-1)/2)*h^k.                              (22)

No uniform bound on h or any lower coefficient is imposed.
In particular D=deg M+2deg h, and at every A pole divisor off
M the zero order of h is 2L.

## 8. A regular cubic square completion, with an explicit field obstruction

In fact h divides a2 in C[u]. At an arbitrary root c of h put
j=ord_c h>=1, m_c=ord_c M>=0. Suppose a=ord_c a2<j and set

    ell0=2j-a>j, H=5ell0-a+2m_c.

From (21), the degree-six and degree-five coefficients of P* have
orders 2j-2m_c and a-2m_c. EVERY lower coefficient has order at
least -2m_c: after replacing b3=M*h^2 the whole expression P*
has denominator at most M^2, and N and all canonical coefficients
are regular at c.

Use w=v^-ell0*Y and scale P*-p by v^H. Its entire reduction is

    Y^5*(eta5+eta6*Y), eta5*eta6!=0.                       (23)

For degree i<=4, its margin is at least
(5-i)ell0-a>0. The nonzero root is simple. All-order Hensel lifting
over C[p] applies even though M can vanish here, since the whole
scaled equation has regular coefficients. The exact orders are

    ord_v w_p=H-ell0,
    ord_v omega=H-ell0-m_c=4ell0+m_c-a>0.                  (24)

The actual companion again evaluates in C[p]((v)), and (24)
makes its residue polynomial in p. If c!=0, take the two unit
t lifts; w is polar and so is x. If c=0, pass to t with u=t^2;
w stays polar and x=(w/t-Sigma)/A is polar as well. The original
full-normalization graph argument is therefore unchanged.
This contradiction proves h divides a2.

Put A2=a2/h in C[u]. The regular source-image polynomial

    T=h*w*z+(A2/2)*z

has ruling degree three. The relation w^2=N+M*z gives the exact
square completion

    P=T^2+R,
    R=(b2-N*h^2-A2^2/4)*z^2
                         +w*(a1*z+a0)+b1*z+b0,            (25)

where R is regular of ruling degree at most four.

T is provably outside F; membership in the full finite normalization
is not established. Indeed PF162's central
equalities give ord_0 a2=ord_0 h, so A2(0)!=0. On the original
whole source line t=0,

    T=Sigma(0)*A2(0)*x,                                   (26)

a nonconstant affine-linear function. If T belonged to F,
the source-intersection identity F intersect C[x,t]=B of PF52
would put it in B. That element would make the actual Keller
map injective on this source line, forcing an automorphism and
contradicting the nontrivial quadratic factor. Thus T is
provably outside F. The identity (25) consequently does not
supply a lower-ruling target coordinate or a false closure of
the remaining problem.

The central leading identity of PF162 also gives
A2(0)^2=4b2(0). Hence R's quadratic z coefficient vanishes at0,
and its restriction to t=0 has degree at most one.

R is outside F as well. Otherwise T^2=P-R belongs to
F intersect C[x,t]=B. For a Keller map, the pullback of an
irreducible reduced target divisor is reduced, because the map
is etale. It is nonempty since its defining target polynomial
has a nonconstant source pullback. Distinct target divisors have
no common source component, since a source curve cannot map to
their finite intersection under a quasi-finite map. Thus if a
target polynomial pulls back to a square, every exponent in its
irreducible target factorization is even. It is itself a square
up to a nonzero constant. Applied to T^2, this would put T in B,
contradicting (26). Both lower-ruling functions T and R therefore
remain outside the actual target field.

## 9. Updated complete ledger and limits

Every actual ruling-six coordinate must now have at least one A
pole divisor. If these divisors have slopes L_i and delta_i=4L_i,
PF164's complete ledger specializes to

    1<=b=[E:F]<=sum_i(2L_i-5)<=floor(D/2)-5,
    D=deg b3, D>=16,
    ceil(alpha*(D+3m-6)/6)<=b at every infinity point,
    alpha in {1,2}.                                       (18)

All six original points over t=infinity remain unramified for t.
The no-A subfamily is now CLOSED, stronger than its former necessary
D>=48 estimate. Since D>=4*sum_i L_i, the infinity lower bound also
shows that not all A slopes can equal three: at least one L_i>=4.

The [strict retained-source infinity bounds](RULING_SIX_RETAINED_SOURCE_PACKETS_AND_STRICT_INFINITY_BOUNDS.md)
strengthen the current remaining-class bounds to D>=22, b>=5
and original generic map degree 2b>=10. They also retain the
residue degree over each infinity image and the complete
nonlinear infinity packet alternatives.

The remaining A class has central degree two, arbitrary coefficient
degrees, possible nonlinear infinity residues, arbitrary higher
ramification and the full quadratic central image class. Neither
unramified source infinity nor the local residue conductor supplies
a missing trace-integrality or normal-tail termination step.

The checker controls the full reductions (10),(12), all simple-root
and exact coframe orders, weighted commutation identities, multiplicity
lattice/divisibility, finite target-subduction mechanism, coupled
degree/valuation identities, reciprocal M-root p-response orders,
the global square/odd-companion normal form, and the regular
square completion with its central source-line field obstruction.
Uniform Hensel recursion, closed full-normalization images and
integer weight descent establish the theorem. The checker passes
191 exact controls; finite controls do not replace those arguments.

Replay with

    .venv/bin/python research/jacobian/affine-complements/verify_ruling_six_simple_escape.py

No new external classification is imported. Independent whole-proof
review, formal verification and external human review remain unclaimed.
[PF188](RULING_SIX_SQUARE_DEFORMATION_AND_MARKED_QUADRATIC_EXCLUSION.md) excludes the retained marked quadratic ruling-six class.
HC4/JC2 and derivative integrality remain open.
