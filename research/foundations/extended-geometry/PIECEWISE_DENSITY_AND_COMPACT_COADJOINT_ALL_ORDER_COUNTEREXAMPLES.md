# Piecewise-polynomial height laws and all-order counterexamples on compact coadjoint products

The normalizer endpoint method, credited to Long, gives regular counterexamples on a product with a split torus when an actual height coordinate has a suitable polynomial measure law. The finite piecewise-polynomial extension and shear below give ONE fixed unit marker nonzero at EVERY positive order. Classical Duistermaat--Heckman density theory realizes the criterion on compact real coadjoint orbits. Prescribed compact Chevalley rational models give witnesses over every characteristic-zero field, with genuine two-factor instability over every formally real field.

A phase F is pure when every positive pure moment vanishes; it is safe when every fixed regular marker has eventually vanishing marked moments. A functional kernel is Mathieu when every pure phase is safe. An absolute unit-nonpower cannot be written u*g^e for any unit u and integer e>=2 over an algebraic closure. All positivity statements below concern real values or positive rational values; after scalar extension to a general field they assert nonvanishing, not an ordering on that field.

## 1. The distinguished-endpoint criterion

Let A be a geometrically integral finitely generated Q-algebra, with X=Spec(A) smooth and normalized Q-linear functional lambda. Suppose z in A satisfies, for EVERY real polynomial g,

    lambda_R(g(z))=sum_i integral_(a_i)^(a_(i+1)) rho_i(u)g(u)du,

where a_0<...<a_r=1 are finitely many real knots, each rho_i is a real polynomial of degree at most an INTEGER b>=0, and the resulting law is nonnegative with mass one. Require the last density polynomial rho_right not identically zero. The real version instead takes A over R, a normalized R-linear functional, and f with real coefficients; its conclusion is a real theorem. A rational height can be normalized to this endpoint only when that normalization is actually defined over Q; no such operation is inferred for arbitrary algebraic endpoints.

Suppose f in Q[u] has f(0)=1, a SIMPLE zero at1, and a zero of order at least2 at every other NONZERO knot. Knot0 needs no root. Such a polynomial exists when all other knots are algebraic and1 is rational: take (1-u) times the SQUARE of the product of the distinct minimal polynomials of all nonzero knots other than1, and normalize its value at0. No such minimal polynomial can contain1, so its zero at1 remains simple. Extra Galois-conjugate roots are harmless.

Let q be the vanishing order of rho_right at1. Then 0<=q<=b and

A_q=(-1)^q rho_right^(q)(1)>0.

Indeed a nonzero polynomial nonnegative on the interval immediately to the left of its support endpoint has a positive first inward Taylor coefficient.

Put L=b+1, T=1+t and

p=T t^-2 f(Tz)^2,
F=p^L (1+t^L p),
H=t^(-(q+1)).

**Theorem 1 (piecewise-polynomial distinguished endpoint).** For every m>=1, with k=Lm,

lambda_bar(F^m)=0,
lambda_bar(H F^m)
 = A_q f'(1)^(2Lm) (2Lm)!/(2Lm+q+1)! >0.                 (1)

Here lambda_bar=CT_t lambda. Thus H is ONE fixed unit marker nonzero at EVERY positive order. If z has a geometric zero on X, then F is absolutely nonpower even up to units. The same rational F,H and identities remain valid over every characteristic-zero extension of Q. The real-base-field version uses arbitrary real knots and a real polynomial f.

The proof tracks every positive and negative exponent. Appendix A supplies a broader algebraic-knot criterion without a rational distinguished endpoint; its conclusion about a finite marker family is explicitly weaker.

## 2. Negative endpoint part

Write rho_i(u)=sum_{j=0}^b c_ij u^j. For each integer k>=L let

J_jk(v)=integral_0^v w^j f(w)^(2k)dw,
C_jk=sum_i c_ij(J_jk(a_{i+1})-J_jk(a_i)),
E_k(T)=sum_{j=0}^b C_jk T^(k-j-1).

E_k is a polynomial of degree at most k-1. Direct interval integration gives

lambda_R(p^k)
 =t^-2k sum_j T^(k-j-1) sum_i c_ij
        (J_jk(Ta_{i+1})-J_jk(Ta_i)).                       (2)

Replacing each primitive value by its value at the knot produces the endpoint part

N_k(t)=t^-2k E_k(1+t).

Its largest Laurent exponent is at most -k-1. Thus it contributes no constant term, or positive coefficient, even if some of its coefficients vanish.

At any nonzero knot a other than1, the double root of f makes J_jk' vanish to order at least4k. Consequently

J_jk(Ta)-J_jk(a)=O(t^(4k+1)).

After multiplication by t^-2k and the polynomial factors T^(k-j-1), every such error starts at exponent at least2k+1. At knot0 the error is identically zero.

## 3. The unique selected-endpoint contribution

The error contributed by the right endpoint1, after summing the monomial density coefficients, is EXACTLY

R_k(t)=t^-2k T^(k-1) integral_1^T rho_right(w/T) f(w)^(2k)dw. (3)

This follows from sum_j c_j T^-j w^j=rho_right(w/T). Since1 is the outer support endpoint, it appears with positive orientation from its last interval and no adjacent nonzero-density interval lies to its right.

Set w=1+s. Taylor expansion at1 gives

rho_right((1+s)/(1+t))
 =rho_right^(q)(1)/q! * (s-t)^q/(1+t)^q
   + terms of higher total order in(s,t),
f(1+s)^(2k)=f'(1)^(2k)s^(2k)+O(s^(2k+1)).

Therefore R_k starts at exponent q+1. Its coefficient there is

rho_right^(q)(1)/q! * f'(1)^(2k)
     integral_0^1 v^(2k)(v-1)^q dv
 =(-1)^q rho_right^(q)(1) f'(1)^(2k)
       (2k)!/(2k+q+1)!.

The beta integral cancels q!. All coefficients of t^1 through t^q are zero. Because q+1<=L<=k, every other knot error starts strictly later than this selected coefficient.

Combining Sections2-3, lambda_R(p^k) has its whole negative endpoint part bounded above by exponent -k-1, its first possible positive coefficient at q+1, and that coefficient is the positive expression displayed above. It has no constant term.

These are formal local expansions at t=0 of an actual Laurent polynomial. No convergence or unrestricted shift-domain assertion is being substituted for coefficient calculation.

## 4. The nonpower modification preserves every marked order

Expand

F^m=sum_{h=0}^m binom(m,h) t^(Lh) p^(Lm+h).

Set k=Lm+h. Its endpoint part has largest exponent at most

-k-1+Lh=-Lm+(L-1)h-1 <= -m-1.

Its selected positive part starts at q+1+Lh. Its other-knot errors start at2k+1+Lh. Thus every summand has zero constant coefficient, proving purity.

For h>=1, the selected positive part starts strictly above q+1 and the other errors also start above q+1. The endpoint part remains strictly negative. Hence the coefficient of t^(q+1) in lambda_R(F^m) comes ONLY from h=0, with binomial coefficient1. Section3 with k=Lm proves(1).

This includes L=1 and q=0: k=m+h and the largest endpoint exponent is -m-1, while every h>=1 selected error starts at least2. No special exception is hidden in the inequalities.

### The greatest unsafe integer monomial exponent

Put s=q+1. The negative endpoint part of lambda_R(F^m) has largest exponent at most -m-1, and its positive part begins exactly at s with nonzero coefficient. For every fixed integer ell>-s, the marker t^ell therefore has eventually zero moments. If -s<ell<0, both support bounds give zero at every order. If ell>=0, the negative part stays negative for m>=max(1,ell), while the positive part remains positive. The marker at ell=-s survives every order by the proved coefficient identity. Hence -s is exactly the GREATEST unsafe integer Laurent-monomial exponent for this displayed F. Equivalently s is the least positive magnitude e for which t^-e detects, and any UNSAFE monomial marker must have ell<=-s. No assertion that every smaller exponent is unsafe is made. For the coadjoint phases below this greatest exponent is -D.

## 5. A simple zero of the modifying factor

At a geometric point x0 with z=0,

p(x0,t)=(1+t)/t^2,
W(x0,t):=1+t^L p(x0,t)=1+t^(L-2)+t^(L-1).

For L=1, W=2+t^-1 has the simple nonzero root -1/2. For L=2, W=2+t has the simple root -2. Both roots differ from-1, so p is a unit there.

For L>=3 set n=L-2>=1 and g(t)=1+t^n+t^(n+1). It has no root0 or-1. If g and g' had a common root alpha, then alpha=-n/(n+1), since g'=t^(n-1)(n+(n+1)t). This rational number is not an integer and cannot be a root of the monic integer polynomial g. Thus g is squarefree.

Choose any root alpha. At(x0,alpha), p is nonzero and the derivative of W in the t direction is nonzero. Smoothness gives a reduced prime hypersurface through this point, with valuation1 on W and valuation0 on p. The valuation of F=p^L W is therefore1. It cannot equal a unit times a proper power. Replay this argument after embedding Qbar into an algebraic closure of any characteristic-zero K; it proves absolute nonpower after all base changes.

When f and z and lambda are rational, the left side of(1) is rational. In particular its positive nonzero value for m=1 shows A_q is positive rational: f'(1) is a nonzero rational and all factorials are rational. Formula(1) therefore transfers faithfully to every characteristic-zero field. No embedding of that field intoR is used.

### Radial primitivity of the displayed phase

The same geometric height-zero specialization proves that the radial Laurent support has differences of greatest common divisor one. Precisely, F cannot be written u t^a G(t^r) with r>1, a an integer, u a unit of O(X_Kbar), and G in O(X_Kbar)[v,v^-1]. Indeed specialize z=0. Then

    F(t)=t^(-2L)(1+t)^L+t^(-L-2)(1+t)^(L+1).

For L=1 this is t^-3(1+3t+2t^2), and for L=2 it is t^-4(2+5t+4t^2+t^3). For L>=3 its lowest coefficient, at exponent -2L, is 1, and the adjacent coefficient is L except that for L=3 the overlap adds 1, giving 4. Thus two adjacent global Laurent coefficients are nonzero, since their specializations are nonzero. Multiplying by a Laurent unit changes support only by translation: over a domain every Laurent unit is a coefficient unit times a power of t. A pullback along t->t^r has all support differences divisible by r. This proves the stated radial primitivity; it is not a claim about arbitrary morphisms or arbitrary phase transformations.

## 6. Endpoint order on compact adjoint orbits without Darboux theory

For a compact adjoint orbit, use an invariant positive form < , >_c, and identify it with its coadjoint orbit. Let x be the basepoint, xi its varying compact conjugate, and let

w(xi)=1-(<x,x>_c-<x,xi>_c)/c, c>0.

All xi have the same norm as x. Hence

1-w(xi)=||x-xi||^2/(2c).                                 (4)

Thus w has its UNIQUE maximum1 at xi=x. Suppose the orbit is smooth real2D-dimensional and its invariant probability measure has a smooth positive local volume density, as normalized Liouville measure does.

In local coordinates y in R^(2D) centered atx, the embedding is xi(y)=x+Ay+O(|y|^2), with A injective. Consequently

1-w(xi(y))=Q(y)+O(|y|^3),

where Q(y)=||Ay||^2/(2c) is positive definite. The measure density is a(y)=a(0)+O(|y|), with a(0)>0. Compactness and uniqueness of the maximum put the entire superlevel set{w>1-epsilon} in this chart for all small epsilon. Positive definiteness bounds its coordinates by O(sqrt(epsilon)).

Rescale y=sqrt(epsilon)u. The rescaled domain is uniformly bounded, its indicator tends almost everywhere to the ellipsoid{Q(u)<1}, and a(sqrt(epsilon)u) tends to a(0). Dominated convergence therefore gives

Prob(w>1-epsilon) ~ C epsilon^D, C>0.                     (5)

If the right-endpoint density is a nonzero polynomial of degree at most D-1, its tail instead has leading order epsilon^(q+1), where q is its vanishing order. Comparing with(5) forces q=D-1.

This proof uses only the distance identity, smooth embedding, positive smooth volume and the DH polynomial degree bound. It imports no equivariant Darboux normal form or stationary-phase formula, corrected or otherwise.

## 7. All compact real coadjoint models under the stated affine hypotheses

Let G be a connected reductive real algebraic group, let K=G(R)^0 be compact, and let X=G/H be a positive-dimensional affine geometrically homogeneous coadjoint orbit through a real point, with reductive algebraic stabilizer H. Assume its compact K-orbit M is Zariski dense. One may equivalently use an adjoint orbit and a nondegenerate invariant real bilinear identification with the coadjoint representation. Write dim_R M=2D>0.

The KKS form is explicit: on an adjoint orbit, with invariant bilinear identification, put

    omega_xi([a,xi],[b,xi])=<xi,[a,b]>.

It is well-defined because adding a centralizing a changes the right side by <[xi,a],b>=0. Nondegeneracy follows in the same way from invariance and the stabilizer description. Differentiating this invariant form on fundamental vector fields gives closure by the Jacobi identity. The linear coordinate mu_y(xi)=<xi,y> satisfies dmu_y=iota_(Y_y)omega, with a consistent choice of fundamental-action sign. Thus M is Hamiltonian, not merely a homogeneous manifold carrying an unspecified measure.

Normalized Liouville measure is K-invariant and therefore is the unique invariant probability measure on this compact transitive orbit. Its integrals on regular functions equal the normalized algebraic Reynolds functional. Indeed K contains a real analytic neighborhood of the identity in the smooth real group and is Zariski dense; Haar invariance is consequently algebraic invariance. Characteristic-zero reductive complete reducibility gives uniqueness. A nonzero real regular function has nonzero pullback to G and is nonzero somewhere on the dense compact orbit, so its squared integral is strictly positive.

Choose a maximal torus of K. Integral lattice vectors in its Lie algebra avoid the finitely many root hyperplanes after a generic choice; for example choose integer coordinates outside the zero set of the finite product of their nonzero rational root linear forms. Hence there is a REGULAR periodic circle generator y. After conjugating the base point into this compact Cartan, the fixed points of this circle on M are precisely M intersect t, the finite Weyl orbit of that point: centralizing regular y forces the point into t, and compact Cartan conjugacy identifies the intersection with the Weyl orbit. The moment map is nonconstant, because a trivial circle action would fix all of positive-dimensional M.

Normalize the real regular linear moment coordinate to

    z=(mu_y-min_M mu_y)/(max_M mu_y-min_M mu_y).

It has image [0,1] and has a real zero. Replace the circle by its effective finite quotient. The moment map is proper because M is compact. The actual Duistermaat--Heckman density theorem, Section3 Lemma3.1, Proposition3.2 and Corollary3.3, gives a nonnegative absolutely continuous density polynomial of degree at most D-1 on each interval of regular values. Critical points are fixed points, so there are finitely many critical values and finitely many pieces. Finite isotropy is allowed in the primary theorem; global freeness is not assumed.

At each maximum the linearized circle generator A on the tangent space is invertible. To see this, average a Riemannian metric; its exponential map identifies a fixed-point neighborhood with the fixed vectors of the linear action. The fixed set here is isolated, so ker A=0. Differentiating dmu=iota_Y omega at a fixed point gives

    Hess(mu)(v,w)=omega(A v,w).

The Hessian is nondegenerate because A and omega are nondegenerate; at a maximum it is therefore definite. In smooth local coordinates, 1-z is a positive definite quadratic form plus O(|u|^3). Liouville measure has a smooth positive local density. The volume rescaling proof of Section6 applies at each maximum, and their finitely many positive leading constants add. Compactness puts the full top superlevel set in these charts. Thus its probability tail is asymptotic to C epsilon^D, C>0. The endpoint polynomial consequently has vanishing order q=D-1.

Take b=D-1,L=D and

    f(u)=(1-u) product_(a_i not in{0,1})(1-u/a_i)^2,
    p=(1+t)t^-2 f((1+t)z)^2,
    F=p^D(1+t^D p), H=t^-D.

**Theorem 2 (real compact coadjoint products).** On every real affine coadjoint model just specified, F is pure and

    barlambda(H F^m)=C_* f'(1)^(2Dm)(2Dm)!/(2Dm+D)!>0

for EVERY m>=1, with C_*=(-1)^(D-1)rho_right^(D-1)(1)>0. It is an absolute unit-nonpower and radially primitive as defined in Section5. The compact factor has only zero as a pure real phase by strict squared positivity, and the split torus factor is Mathieu by scalar Laurent constant-term theory. Their product therefore fails the Mathieu property. Point orbits, arbitrary nonsymplectic compact homogeneous varieties, and unspecified affine models are not included. This is a REAL theorem; the real normalization of its moment coordinate is not by itself a rational-model construction.

## 8. The actual compact Chevalley rational models

We now give actual Q-defined examples, so the constructed phase works over every characteristic-zero field.

Take a simple complex Lie algebra with a Chevalley basis H_i,E_alpha,E_-alpha normalized so that the structure constants are integers and B(E_alpha,E_-alpha)>0. Let i^2=-1. Inside its Q(i)-form, define

k_Q = span_Q { i H_i, E_alpha-E_-alpha,
                        i(E_alpha+E_-alpha) : alpha positive }.

This is a Q-Lie algebra. One way to check closure is to take the fixed algebra of conjugation onQ(i) composed with the Chevalley involution H_i->-H_i, E_alpha->-E_-alpha. Equivalently substitute the displayed basis in the integer Chevalley brackets. It complexifies to the original simple Lie algebra. The Killing form is negative definite on its real span: the Cartan block is minus the positive root-space Cartan form, each real root-plane basis vector has norm -2B(E_alpha,E_-alpha)<0, and different root planes are orthogonal.

Let G_Q=Aut(k_Q)^0, the connected identity component of the algebraic bracket-preserving subgroup of GL(k_Q). It is the adjoint connected semisimple Q-group of the chosen type. OverR its automorphisms preserve the positive definite form -B, so G_Q(R) is compact. This provides an explicit linear-algebraic Q-model; it is not an invocation of a hypothetical uniform Q-form for arbitrary real homogeneous spaces.

A primary antecedent for this construction is Kerr--Robles, Hodge theory and real orbits in flag varieties, exact arXiv1407.4507v1, whole Section6.3 printed38-39, especially formulas6.6-6.12, Definition6.9 and Remark6.10. Their construction with grading element 0 makes every root even and gives precisely this compact Q-form. Their printed39 also explicitly notes that i times the fundamental coweights are Q-defined. No Hodge-theoretic orbit or later boundary theorem is imported.

Choose a nonzero rational coweight h and set x=i h in k_Q. The real group exp(Rx) is a closed circle after scaling x by a positive integer: all root eigenvalues are i times rational numbers, hence the adjoint representation has a common period. Its algebraic closure S is a Q-torus, because the corresponding complex one-parameter torus is stable under the Q(i)/Q descent; its Galois action is inversion. Let H=C_G(S). At the full group level H=C_G(x), not just at the level of Lie algebras. Over C, Ad(g)x=x is equivalent to commuting with every exp(tx). This periodic one-parameter circle is Zariski dense in S_C, so g centralizes S_C; conversely centralizing S_C fixes x by differentiation. Characteristic-zero algebraic groups are smooth, so equality of complex point subgroups gives equality of subgroup schemes and descends to Q. Equivalently, in the faithful adjoint representation the ad_x eigenspaces are the integer cocharacter weight spaces after common scaling, so fixing x forces commutation with S. This excludes a disconnected-centralizer shortcut.

A torus centralizer in a connected reductive group is smooth connected reductive. The homogeneous quotient X_Q=G_Q/H is therefore a smooth geometrically integral affine variety (Matsushima in characteristic0). Its regular orbit map is

gH -> Ad_g x.

Thus the coordinate expressions below are actual regular functions on the affine quotient. The Q-point x supplies the base point of the compact orbit; Section9 supplies a geometric height-zero Weyl conjugate for nonpower.

Normalized Reynolds trace on X_Q is unique and commutes with every field extension. OverR the compact orbit through x is Zariski dense: G_Q(R)^0 contains a real analytic neighborhood of the identity in the smooth real group and hence is Zariski dense, and the quotient orbit map is dominant. Its Haar trace agrees with the normalized Liouville measure as in Section7.

## 9. Positive pairing conventions and the finite rational knot set

Use the compact invariant inner product

    <u,v>_c=-B(u,v)/(2h_dual),

and transport it to the rational coweight space by

    (h,h')_cw:=<ih,ih'>_c=B(h,h')/(2h_dual).

This is the POSITIVE root-Cartan pairing. In the simply laced E6/E7 examples its matrix on fundamental coweights is the inverse Cartan matrix. Minus the Killing form on the split h itself is not being called positive.

Put x=ih, r=<x,x>_c=(h,h)_cw and sigma=min_(v in W h)(h,v)_cw. Some Weyl element moves the nonzero h in a simple root system, and Cauchy--Schwarz gives sigma<r. Use the Q-regular height

    z(Ad_g x)=(<x,Ad_g x>_c-sigma)/(r-sigma).

It is a nonzero real scalar multiple and an additive translate of the circle moment map generated by x. It has unique maximum 1 at x, by equality in Cauchy--Schwarz. It has minimum 0 at a Weyl conjugate.

Differentiation in a direction [a,xi] shows that a critical point satisfies [x,xi]=0, by nondegeneracy and invariance of the form. If xi commutes with x, both lie in a maximal torus of the compact centralizer of x: choose a maximal torus of that centralizer containing the closure of exp(Rxi); its center contains exp(Rx). Conjugating that torus into the fixed maximal torus by the centralizer leaves x fixed. Two elements of a fixed compact Cartan conjugate in K are Weyl-conjugate. Therefore every critical height is

((h,w h)_cw-sigma)/(r-sigma), w in W.

Conversely each Weyl conjugate gives a critical point. These are exactly the critical values, forming a finite rational set A contained in [0,1] and containing 0 and 1. They are rational for rational h. This argument accommodates positive-dimensional fixed components; it does not use the original uncorrected isolated-stationary-point formula.

The symplectic dimension is twice the number D of positive roots with alpha(h)!=0: each such root pair gives one real tangent2-plane. Consequently the density bound is b<=D-1 and the construction may take L=D. The height satisfies

    1-z(Ad_g x)=||x-Ad_g x||_c^2/[2(r-sigma)].

Its maximum is unique, so Section6 proves the endpoint order q=D-1 using a smooth embedding and positive volume alone. No Darboux or stationary-phase theorem is needed. A Weyl conjugate with z=0 exists over Qbar as well, supplying the geometric zero for absolute nonpower after all field extensions.

## 10. The field-uniform Chevalley-orbit theorem and factor property

Define the finite rational Weyl-product polynomial

    f(u)=(1-u) product_(a in A minus{0,1})(1-u/a)^2.

It has f(0)=1, a simple selected endpoint, and double roots at every other nonzero knot. Put

    p=(1+t)t^-2 f((1+t)z)^2,
    F=p^D(1+t^D p), H=t^-D.

All coordinates and coefficients are Q-defined on X_Q x G_m. Normalized Reynolds traces commute with field extension. The real compact Haar trace equals normalized Liouville measure. Hence the rational trace identities of Theorem1 apply over Q and after Q->K.

**Theorem 3 (standard compact Chevalley rational coweight orbits).** For every simple root type, every nonzero rational coweight h in the explicit construction of Section8, and every characteristic-zero field K, the displayed regular F,H on the base change of X_Q x G_m satisfy

    barlambda_K(F^m)=0,
    barlambda_K(t^-D F^m)
      =C_* f'(1)^(2Dm)(2Dm)!/(2Dm+D)! !=0, m>=1,

where C_*=(-1)^(D-1)rho_right^(D-1)(1) is positive rational. The equality at m=1 proves its rationality from the rational trace and nonzero rational f'(1); its closed numerical value is unnecessary. F is absolutely unit-nonpower and radially primitive. The real positivity is only a way to establish that these rational scalar values are nonzero; no ordering of a general K is asserted.

For completeness the individual compact factor is Mathieu over every formally real K. For a given finite Q-independent span a_i in O(X_Q), the Gram matrix G_ij=lambda_Q(a_i a_j) is rational positive definite over R: every nonzero real combination has strictly positive squared Haar integral on the dense compact orbit. Rational LDL^T elimination gives positive rational pivots d_i and an invertible rational linear change ell_i(c), with

    lambda_K((sum c_i a_i)^2)=sum_i d_i ell_i(c)^2.

Each d=A/B>0 is the sum of AB copies of (1/B)^2. In a formally real field a zero sum of squares forces every term zero, since otherwise division by a nonzero square expresses -1 as a sum of squares. Thus the Gram form is anisotropic and only zero is pure. This uses no choice of ordering and no embedding of K into R.

The scalar Laurent factor is Mathieu over every characteristic-zero field: embed the finitely generated coefficient field of a given phase into C, apply DvK's scalar Newton-polytope theorem, and strictly separate the finite support from zero. For any fixed marker the support of high powers has strictly positive separating weight beyond the marker's finite range, so its constant term vanishes. Only that finite coefficient field is embedded, and no function-valued DvK theorem is used. Consequently both factors are Mathieu over every formally real K while their product fails by the displayed pair.

This theorem concerns the specified compact Chevalley Q-forms and their rational coweight orbits. It does not claim that all real coadjoint models descend to these rational models, that arbitrary twisted forms have the same coordinate formulas, or that new phases over an extension descend.

## 11. EIII: E6 with a D5 times circle stabilizer

Fix the simply laced Cartan diagram with chain0-1-2-3-4 and extra node5 attached to2. Let h be fundamental coweight0, characterized by alpha_j(h)=delta_{j0}. Deleting node0 givesD5. The model of Section8 therefore has compact orbit customarily denoted EIII, with stabilizer of type (Spin10 times U1)/finite center.

The finite root and weight computations give:

number of roots72; Weyl orbit of h has27 elements;
D=16; (h,h)_cw=4/3;
{(h,w h)_cw:w in W}={4/3,1/3,-2/3}.

Here r=4/3,sigma=-2/3, so the actual normalized height is z(Ad_g x)=(<x,Ad_g x>_c+2/3)/2. Its critical values are 0,1/2,1, its unique maximum is 1 at x, and D=16. Define

    f(u)=(1-u)(1-2u)^2,
    p=(1+t)t^-2 f((1+t)z)^2,
    F=p^16(1+t^16 p), H=t^-16.

Here f'(1)=-1. Theorem3 gives, for every characteristic-zero base change and every m>=1,

    barlambda(F^m)=0,
    barlambda(t^-16 F^m)=C_* (32m)!/(32m+16)! !=0.

The phase is absolutely unit-nonpower and radially primitive. Over every formally real field the separate factors are Mathieu. No arbitrary EIII form is substituted for this prescribed Q-model.

### EVII: E7 with an E6 times circle stabilizer

Use chain0-1-2-3-4-5 and node6 attached to2. Let h be fundamental coweight5. Deleting node5 leavesE6. The compact orbit is customarily EVII, with stabilizer of type (E6 times U1)/finite center.

The exact finite data are:

number of roots126; Weyl orbit of h has56 elements;
D=27; (h,h)_cw=3/2;
{(h,w h)_cw:w in W}={3/2,1/2,-1/2,-3/2}.

Here r=3/2,sigma=-3/2, so the actual normalized height is z(Ad_g x)=(<x,Ad_g x>_c+3/2)/3. Its critical values are 0,1/3,2/3,1 and D=27. Define

    f(u)=(1-u)(1-3u)^2(1-(3/2)u)^2,
    p=(1+t)t^-2 f((1+t)z)^2,
    F=p^27(1+t^27 p), H=t^-27.

Again f'(1)=-1. At every positive order on every characteristic-zero base change,

    barlambda(F^m)=0,
    barlambda(t^-27 F^m)=C_* (54m)!/(54m+27)! !=0.

The same nonpower, radial primitivity, and formally real two-factor conclusions hold. Here and in EIII, C_* is the positive rational endpoint coefficient for the particular normalized orbit. The two constants need not agree.

### A nonspherical regular G2 example

Take type G2 and h with both simple-root evaluations 1. All six positive roots have nonzero evaluation: alpha,beta,alpha+beta,2alpha+beta,3alpha+beta,3alpha+2beta. Hence H is a maximal torus and D=6. The finite Weyl-product f from Section10 specifies the rational phase without needing a closed density formula, and the one marker is t^-6 at every order. The complex affine quotient G2/T has dimension 12, whereas a Borel has dimension 8, so it is not spherical. No complex spherical classification is a premise of this construction.

## 12. Exact primary inputs and attribution

The endpoint calculations, valuations, radial-support arguments, local-volume asymptotics, KKS verification, rational-coordinate construction and scalar transfers are proved above. The following classical inputs have separate narrow roles.

### The actual Duistermaat--Heckman density theorem

J. J. Duistermaat and G. J. Heckman, *On the Variation in the Cohomology of the Symplectic Form of the Reduced Phase Space*, Inventiones Mathematicae69(1982),259--268, DOI[10.1007/BF01399506](https://doi.org/10.1007/BF01399506). Exact [original scan](https://gdz.sub.uni-goettingen.de/download/pdf/PPN356556735_0069/LOG_0035.pdf), SHA256 `e314c2033972ab39d244b8bc98ddfb37d8ade96c8c36bfa95e05d82fc1033be5`.

The import is Section3 Lemma3.1, Proposition3.2 and Corollary3.3, printed264--266, with Theorem1.1 and its complete proof in Sections1--2. A connected effective Hamiltonian torus action with proper moment map has an absolutely continuous Liouville pushforward with polynomial density of degree at most n-l on each regular-value component, where dimM=2n and dimT=l. The proof treats finite isotropy using V-manifolds; global freeness is unnecessary. Compactness supplies properness here. The exact original printed259--268, including the note added in proof, were inspected.

The printed268 correction says the simple global stationary-phase formula(4.6), Theorem4.1, is only correct as printed for isolated fixed points. That uncorrected formula is NOT used, even in applications with isolated fixed points. The endpoint-order proof uses elementary local volume rescaling, not a Darboux theorem, stationary phase or an unread later addendum.

### A specified rational compact form

Matt Kerr and Colleen Robles, *Hodge theory and real orbits in flag varieties*, exact [arXiv:1407.4507v1](https://arxiv.org/pdf/1407.4507v1), SHA256 `7d60954c11e8c22ee2353f45ad89dd26a427347e97409668d26c1621aac18a8a`.

Only Section6.3, printed38--39, is imported: the normalized Chevalley basis, formulas(6.6)--(6.12), Definition6.9 and Remark6.10 give the integral/rational real-form construction and rational i times fundamental coweights. Grading element zero makes all roots even and gives the compact rational form used here. Closure and negative Killing form are also explained in Section8. The surrounding paper's Hodge classification and later real-orbit theorems are not imported.

### Reductive projection, torus centralizers and affineness

The same-field normalized trace, scalar compatibility and generic rational Gram argument have complete local proofs in [STANDARD_COMPACT_CLASSICAL_PRODUCT_COUNTEREXAMPLES.md](STANDARD_COMPACT_CLASSICAL_PRODUCT_COUNTEREXAMPLES.md), Sections3 and13; these GENERIC propositions apply here, not that document's seven model-specific witness conclusions. Their algebraic input is J. S. Milne, *Algebraic Groups*(2017), exact physical PDF SHA256 `503f1a7d511632af05541c52fb2bf6fc36f5a75b0323f2f9110f03dcce11fd67`: complete reducibility22.42--22.43, printed477--478, and Reynolds Section12l/12.53, printed248--249. No stronger all-automorphism or finite-cover trace theorem is imported.

A DISTINCT scope of the same Milne PDF is used for H=C_G(S): connectedness17.37--17.38, printed363--364, and Corollaries17.57--17.59 with their complete proofs, printed371, give smooth, connected, reductive torus centralizers. This geometry is not attributed to complete reducibility merely because both occur in the same book.

Ivan V. Arzhantsev, *Invariant ideals and Matsushima's criterion*, exact [arXiv:math/0506430v1](https://arxiv.org/pdf/math/0506430v1), SHA256 `786bff5625883f4456f6decd82e3fc72c60ace70ced673d58aeb20566889b685`, Proposition2.2 and complete proof, printed2, gives affineness of G/H for reductive H over an algebraically closed characteristic-zero field. Affine field descent gives the Q-model. Only this sufficient direction is needed for rational examples; the real coadjoint theorem already assumes an affine model.

### Scalar Laurent factor and known method

J. J. Duistermaat and W. van der Kallen, *Constant terms in powers of a Laurent polynomial*, Indagationes Mathematicae9(1998),221--231, DOI[10.1016/S0019-3577(98)80020-7](https://doi.org/10.1016/S0019-3577(98)80020-7), exact PDF SHA256 `bb002d12b539a1ff3fdd099f343d3accd921ccf60f8a149f29b07cd0f7e87e1b`. The scalar Newton-polytope theorem and the locally proved finite-coefficient-field/support argument make the torus factor Mathieu. This is needed for genuine two-factor instability, not for the product bad phase itself. No coordinate-algebra-valued Laurent theorem is invoked.

Christopher D. Long's [arXiv:2607.19012v1](https://arxiv.org/html/2607.19012v1), exact HTML SHA256 `abeec76e95cba0e6244f4af31240634e28d657116c8afcf743b58db144dd7ebd`, Theorem2.1 and complete proof and Proposition3.1, gives the known interval--torus construction. His *A Classification of Mathieu's Second Conjecture for Reductive Spherical Pairs*, exact revision `e9edac262b4b64eb68feb0e6e80496cefe4a1dee`, Section5 Proposition5.1 and complete proof, gives the normalizer endpoint method; exact PDF SHA256 `b2bf25fbc2d6211333d4ff3e1f6de0b6562990e869c6078a2e219ef746f1026b`, matching TeX SHA256 `934ad65de5ed4358fb33b29af15728ae86a1f7c5836e72050b00fc948c78df64`. These are method credits, rederived here, not real-descent theorems. Known complex spherical existence does not supply the real coordinate, split radial torus, prescribed all-field Q-model, or nonpower argument. No priority conclusion follows from bounded literature comparison.

## Appendix A. Algebraic knots without a rational distinguished endpoint

This independent broader criterion has a weaker marker conclusion. It does NOT assert that one marker survives every order.

Let A be a geometrically integral finitely generated Q-algebra, X=Spec(A) smooth, and let lambda:A->Q be a normalized linear functional. Write lambda_R for its scalar extension. Assume z in A has the following probability law under lambda_R: there are finitely many real knots a_0<...<a_r and polynomials rho_i in R[u] of degree at most an integer b>=0, such that

lambda_R(g(z)) = sum_{i=0}^{r-1} integral_{a_i}^{a_{i+1}} rho_i(u)g(u)du

for every real polynomial g. The resulting density is nonnegative almost everywhere and has total mass1. In particular it is a non-atomic positive measure; density coefficients need not be rational.

Suppose a polynomial f in Q[u] satisfies f(0)=1 and f(a_i)=0 at every NONZERO knot. Put L=b+1, T=1+t and, in A[t,t^-1],

p = (T/t) f(Tz)^2,
F = p^L (1+t^-L p).

Extend the functional by lambda_bar=CT_t lambda. Then:

* lambda_bar(F^m)=0 for every m>=1.
* For each m>=1, at least one of the L values lambda_bar(t^ell F^m), 1<=ell<=L, is nonzero.
* Consequently at least one fixed unit marker t^ell, 1<=ell<=L, detects F^m for infinitely many m. Thus ker(lambda_bar) is not a Mathieu--Zhao subspace.
* If X has a geometric point x0 at which z=0, then F is not a constant or a unit times a proper power in (A tensor Qbar)[t,t^-1].

The same F and finite candidate markers work after EVERY characteristic-zero field extension Q->K, provided lambda_K is the scalar extension. For the absolute-nonpower conclusion, embed Qbar into an algebraic closure of K and replay the same smooth-point and simple-divisor argument there. A single ell chosen to detect infinitely many rational values works over every such K. This asserts a product counterexample; the separate Mathieu property of X requires the independent positivity argument stated in the main text.

When all nonzero knots are algebraic numbers, the required f always exists: take the product of their distinct monic minimal polynomials over Q, and divide by its nonzero constant coefficient. Rational knots are the special case f(u)=product_{a_i!=0}(1-u/a_i). Adding Galois conjugate roots does no harm.

### Endpoint identity and its precise range

Expand rho_i(u)=sum_{j=0}^b c_{ij}u^j and define polynomial primitives

J_{j,k}(v)=integral_0^v w^j f(w)^(2k)dw.

Let

C_{j,k}=sum_i c_{ij}(J_{j,k}(a_{i+1})-J_{j,k}(a_i)),
E_k(T)=sum_{j=0}^b C_{j,k} T^(k-j-1).

For k>=L this is a polynomial of degree at most k-1 and is divisible by T^(k-L).

Integrating f(Tz)^(2k) on each interval, first for real T near1 and then as a rational-function identity, gives

lambda_R(p^k)
 = t^-k sum_{j=0}^b T^(k-j-1) sum_i c_{ij}
       (J_{j,k}(Ta_{i+1})-J_{j,k}(Ta_i)).

At a nonzero knot a, f has a zero, so J'_{j,k}(w)=w^j f(w)^(2k) vanishes to order at least2k at w=a. Hence

J_{j,k}(Ta)-J_{j,k}(a) is divisible by t^(2k+1).

At a=0 the difference is identically zero. It follows that, in the formal power series ring R[[t]],

lambda_R(p^k) = t^-k E_k(1+t) + O(t^(k+1)).                 (1)

Here lambda_R(p^k) is itself a Laurent polynomial, so (1) determines its entire nonpositive Laurent part. For a shifted monomial t^d p^k the remainder begins at exponent k+1+d. In particular the same endpoint expression determines its nonpositive part whenever d>=-k.

This is a local formal expansion at t=0. No integral substitution at T=0, convergence assertion, or unrestricted extension to d<-k is used.

Since deg E_k<=k-1, the expression t^-k E_k(1+t) has only strictly negative exponents. Thus lambda_bar(p^k)=0 for k>=L.

Moreover,

E_k(1) = integral rho(u) f(u)^(2k)du > 0.                  (2)

Positivity follows because a nonzero polynomial has only finitely many real roots and the positive probability law is absolutely continuous. In particular E_k is nonzero.

### A finite set of markers detects every order

For k>=L the marked moments lambda_bar(t^ell p^k), 1<=ell<=L, are the coefficients of t^(k-ell) in E_k(1+t).

If all L were zero, E_k(1+t), which has degree at most k-1, would have degree at most k-L-1. Translation does not change polynomial degree. But E_k(T) is nonzero by(2) and divisible by T^(k-L), forcing degree at least k-L, a contradiction. This also covers k=L, where the supposed degree bound is -1 and means the zero polynomial.

Thus for each k>=L at least one of these L moments is nonzero. This is an exact polynomial argument, not inference from a finite moment prefix.

For completeness, p^a itself is pure for every a>=L, since every exponent am>=L. This conclusion does not require a purity-root implication.

### The modified phase preserves the marker vector

Expand

F^m = sum_{h=0}^m binom(m,h) t^-Lh p^(Lm+h).

Set k=Lm+h. Then k>=L and k-Lh = Lm-(L-1)h >=m. The endpoint formula applies after multiplication by t^-Lh.

For h>=0 the endpoint part has largest Laurent exponent at most -Lh-1; its remainder begins at k+1-Lh>=m+1. Therefore every summand has zero constant term, proving purity of F.

After multiplication by a marker t^ell with 1<=ell<=L, every summand h>=1 still has endpoint largest exponent at most ell-Lh-1<=-1, and remainder starts at k+1-Lh+ell>0. Thus these summands have zero marked constant term. Only h=0 remains:

lambda_bar(t^ell F^m)=lambda_bar(t^ell p^(Lm)), 1<=ell<=L. (3)

The preceding marker-vector argument proves that the vector in(3) is nonzero for EVERY m. A finite union of finite sets is finite; hence at least one fixed ell has nonzero value for infinitely many m. This gives one unsafe fixed UNIT marker, without claiming that one marker is nonzero at every order.

The case L=1 is included: k=m+h, k-Lh=m, and all inequalities remain valid.

### Absolute nonpower

Work over Qbar and choose x0 in X(Qbar) with z(x0)=0. Since f(0)=1,

p(x0,t)=(1+t)/t,
1+t^-L p(x0,t)=g_L(t)/t^(L+1),
g_L(t)=t^(L+1)+t+1.

The polynomial g_L has no root0 or-1. It is squarefree: if g_L(alpha)=g'_L(alpha)=0, then alpha^L=-1/(L+1) and substituting in g_L gives alpha=-(L+1)/L. The latter rational number has absolute value greater than1, contradicting |alpha|^L=1/(L+1).

Choose any root alpha of g_L. At (x0,alpha), p is nonzero and the partial derivative in t of 1+t^-L p is g'_L(alpha)/alpha^(L+1), which is nonzero. Since X times Gm is smooth, the zero locus of 1+t^-L p is a smooth reduced hypersurface near this point. At least one irreducible prime divisor through the point therefore has valuation1 on 1+t^-L p, and valuation0 on p. Consequently its valuation on F is1.

If F=uQ^d for a unit u and d>1, the same valuation would be divisible by d, a contradiction. This proves nonpower even up to units over the algebraic closure. It also rules out constant-times-proper-power representations.

## Appendix B. Finite Cartan support for the numerical EIII/EVII specializations

The universal standard-Q theorem uses W symbolically and needs no numerical enumeration. The explicit16/27 exponents and small knot polynomials additionally use `verify_piecewise_coadjoint_exceptional_roots.py` and its deterministic output `piecewise_coadjoint_exceptional_root_data.json`. They check ONLY finite Cartan data, not all-order moments or geometry. The output retains the complete Cartan/inverse matrices, full reflection-closed root and weight sets, positive roots, centralizer roots and pairings used in Section11. No weak-circuit or strong-circuit moment prefix is installed as a theorem checker.

From the repository root, reproduce ONLY this finite Cartan support with:

```sh
python3 research/foundations/extended-geometry/verify_piecewise_coadjoint_exceptional_roots.py --output /tmp/pwc-exceptional-root-data.json
```

The output must be byte-identical to the retained data, SHA256 `ae03576c1d869d8fd033efa9bf15d0c9f5ccf5de0bbbb41e306cb907a42e57e7`. This command does not replay the all-order mathematical theorem.


For an auditable derivation, set diagonal entries2, edge entries-1 and remaining entries0 for the E6/E7 graphs in Section11. Exact rational LDL pivots are respectively

    E6: 2,3/2,4/3,5/4,6/5,1/2;
    E7: 2,3/2,4/3,5/4,6/5,7/6,2/7.

Thus both matrices A are positive definite. In simple-root coordinates r, reflection i changes only r_i to r_i-(Ar)_i. Starting from all simple roots and repeatedly applying these reflections generates the complete root system. In simple-root-evaluation coordinates v for weights, reflection i sends v_j to v_j-v_i A_ji; start at the indicated fundamental-coweight unit vector. Each reflection preserves respectively r^T A r or v^T A^-1 v. Integral vectors of fixed norm in a positive definite lattice form a finite set, so exhaustive closure queues terminate, with no truncation or height cutoff.

The exact closures give72/126 roots and27/56 weights. Counting positive roots with nonzero selected coordinate gives16/27. Pairing a weight v with the initial coweight uses the selected row of A^-1 times v, giving {4/3,1/3,-2/3} and {3/2,1/2,-1/2,-3/2}. These are exhaustive finite algebra calculations. Withdrawing this support would affect the numerical EIII/EVII small-knot specializations; it would not affect the generic endpoint criterion or the symbolic finite-Weyl-product theorem.
