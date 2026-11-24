# Universal compact rational coadjoint product counterexamples

The actual orbit, rational height and density statements are proved independently of the endpoint-amplification theorem. The final product construction applies [the algebraic-spline amplification criterion](ALGEBRAIC_SPLINE_ENDPOINT_AMPLIFICATION_COUNTEREXAMPLES.md) only after all its hypotheses have been established. A worked example gives a closed every-order marked formula without a periodic Lie generator.

Existence of complex counterexamples on these nontrivial homogeneous products is already covered by [the published complex classification](https://hexagonmath.org/pdf/2611.00556v1), Theorem 1.1; the compact-group cases also have Long's antecedents. The statements here concern actual arbitrary compact Q-forms, one Q-defined unit marker surviving every order after all characteristic-zero scalar extensions, and the stated divisor/support properties. No new complex-existence or classification claim is made.

## 1. Statement

Let G be ANY connected semisimple algebraic group over Q with compact G(R), let g=Lie(G), and let 0!=x in g(Q). Use the invariant rational positive real form

    (u,v)=-Killing_g(u,v).

There is a maximal Q-torus T with x in t(Q), t=Lie(T), and one may choose y in t(Q) which is regular in g, separates the finite set (Ad(G_C)x) intersect t_C, and satisfies (x,y)>0. Put

    X=Ad(G)x, phi(z)=(z,y), 2D=dim X.

Then X is an actual closed, smooth, geometrically integral affine Q-orbit, canonically G/C_G(x), with connected reductive full centralizer. Its unique normalized rational invariant trace, after real extension, is normalized Liouville measure on the compact connected orbit M=Ad(G(R)^0)x. The phi-pushforward of that probability has compact support, is nonnegative and piecewise polynomial of degree at most D-1, and all its knots are real algebraic. Extended by zero, its density is GLOBAL C^(D-2) when D>=2; for D=1 no continuity is asserted or required. The finite critical points are nondegenerate. Finally phi=0 has an ACTUAL real orbit point and a geometric Qbar-point.

This statement does NOT require y or x to generate a circle, or to be a rational multiple of an integral cocharacter in the real torus. The polynomiality proof uses approximating REAL circle directions solely as an analytic device. D>=1, and zero components of x on some semisimple factors are allowed.

**Theorem 1 (compact rational coadjoint products).** For every G and x above there are a rational regular height phi, a polynomial f in Q[u] with f(0)=1, and an integer L>=max(D,4) for which p=(1+t)t^-2 f((1+t)phi)^2, F=p^L(1+t^L p) and the SAME fixed unit H_marker=t^-D satisfy Lambda(F^m)=0 and Lambda(H_marker F^m)!=0 at EVERY positive order over EVERY characteristic-zero basechange. F is absolutely not a unit times a proper power and has primitive Laurent support. Over every formally real field, the compact orbit factor and split torus factor are Mathieu while their product is not. The same conclusions pull back to every actual G/H with a nonzero rational Lie vector fixed by the FULL subgroup H, as proved in Section 9.

## 2. Semisimplicity, closed orbit, and rational torus

The real Lie algebra g_R is compact semisimple, so minus its Killing form is positive definite and every ad(x) is skew-adjoint. Thus ad(x) is diagonalizable over C. The adjoint representation of the semisimple Lie algebra is faithful; functorial Jordan decomposition in characteristic zero therefore makes x semisimple as an algebraic Lie element, not merely in one accidental representation.

Here is a direct closed-orbit argument. Let m_x be the squarefree minimal polynomial of ad(x), and let chi_x be its characteristic polynomial, both over Q. In g define the closed G-stable set

    Z={z:m_x(ad(z))=0 and charpoly(ad(z))=chi_x}.

For every geometric z in Z, ad(z) is semisimple and its zero-eigenvalue multiplicity is the same fixed integer. Hence dim ker ad(z) is constant. In characteristic zero stabilizers are smooth and Lie(C_G(z))=ker ad(z), so ALL G-orbits in Z have the same dimension. If Ad(G)x had a nonempty boundary in its Zariski closure, that boundary would contain an orbit of strictly smaller dimension: the orbit is locally closed and open dense in its own closure, and the boundary is invariant of smaller dimension. This contradicts the constant dimension in Z. Thus the geometric orbit is closed; closedness descends to Q. It is affine, smooth and geometrically integral, and its orbit map identifies it with the scheme quotient G/C_G(x).

The rational torus is supplied by the exact arbitrary-field result in Milne, Proposition 17.90 and its proof, together with 17.82--17.83. Choose a maximal Q-torus T of C_G(x)^0. The proof of 17.90 gives x in Lie(T) and says that T is maximal in G. This is not an assertion that an arbitrary compact rational Lie element is a rational coweight. The full centralizer is connected reductive as well. To see this without discarding components, work over C and take S to be the Zariski closure of exp(Cx) inside T_C. It is a connected subtorus; no periodicity is assumed. An element fixes x iff it commutes with every exp(tx), iff it centralizes S. Thus FULL C_G(x)_C=C_G(S). Torus centralizers in connected reductive groups are smooth connected reductive (Milne 17.38/17.40 and 17.59), and these geometric properties descend to Q.

For clarity, the orbit dimension is positive. If it were zero, connected G would fix x and x would lie in the zero Lie center of the semisimple group. The Kostant--Kirillov form below is symplectic, so dim X=2D for an integer D>=1.

## 3. Native rational trace, compact orbit, and Liouville probability

The real group K=G(R)^0 is connected compact and Zariski dense in G: a real analytic identity neighborhood in a smooth real group is Zariski dense, and G is connected. Its orbit M=Kx is compact, connected and Zariski dense in X by the dominant orbit map. No claim that all real points of an arbitrary rational orbit lie in M is needed.

The actual affine homogeneous Q-algebra O(X) has the unique normalized Reynolds projection to constants. This is Q-valued on O(X), is invariant under the Lie algebra, and commutes with all field extensions, by the existing characteristic-zero complete-reducibility/normalized-trace theorem. Its real extension is the Haar-orbit integral, by normalization and invariance.

Using the rational form ( , ) to identify adjoint and coadjoint orbits, the KKS form on M is

    omega_z([a,z],[b,z])=(z,[a,b]).

It is well-defined, nondegenerate, closed and K-invariant. For example closedness follows from the Jacobi identity in its invariant exterior derivative; nondegeneracy follows because a pairing with all b vanishes precisely when [z,a]=0. Its normalized positive Liouville measure omega^D/D! divided by its total mass is K-invariant and hence is exactly the Haar probability on Kx. The height phi_y(z)=(z,y) is a Hamiltonian component, with vector field [y,z] and the convention i_X omega=dphi. Thus the rational normalized trace is the actual Liouville pushforward functional required by the density theorem.

## 4. Rational generic y, finite algebraic critical scheme, and Morse Hessians

Over C the semisimple orbit meets t_C in a finite Weyl orbit. One can derive this directly from maximal-torus conjugacy. If z=Ad(g)x lies in t_C, then T_C and gT_Cg^-1 are maximal tori of C_G(z), which is connected reductive as above; conjugate them by an element fixing z to make g normalize T. This gives z=w x. Consequently P=X_C intersect t_C is a finite set of distinct algebraic points.

Choose y in t(Q) outside all root hyperplanes and all hyperplanes (z_i-z_j,y)=0 for distinct points of P, and, if desired, outside (z_i,y)=0 as well. These are finitely many proper real linear subspaces when restricted to t_R: its real span complexifies to t_C, and the pairing on t_C is nondegenerate. Their complement meets the open half-space (x,y)>0. Since t(Q) is dense in t_R, such a rational y exists. Rational points of the torus itself, and rational integral cocharacters, are irrelevant to this density assertion about its Lie vector space.

Differentiation of phi along [a,z] shows that z is critical iff [z,y]=0. Since y is regular, its Lie centralizer is t_C, so the geometric critical set is exactly P. The critical scheme is defined over Q by these equations. At any of its points z, the Hamiltonian linearization on T_z X=[g,z] is ad(y). This tangent space is invariant under ad(y), and it has zero intersection with t_C: semisimplicity of ad(z) gives im ad(z) intersect ker ad(z)=0, while t_C is in that kernel. The only kernel of ad(y) is t_C. Hence ad(y) is invertible on T_z X. Since the Hessian is omega_z(ad(y)u,v), it is nondegenerate. The critical scheme is therefore finite and reduced in characteristic zero, hence finite etale over Q. In particular every actual real critical value is REAL ALGEBRAIC, and the chosen y separates them.

For real circle approximations below, the real critical set P_real=M intersect t_R is the same fixed finite set for EVERY regular y' in t_R. The remaining complex critical points, if any, play no role in real measure integration. This avoids needing a rational Weyl representative or a claim that every geometric critical point is a rational point.

## 5. Ordinary Morse oscillatory decay gives GLOBAL C^(D-2)

We use the ordinary stationary-phase BOUND, not an exact Duistermaat--Heckman localization formula. If h is a smooth real function with finitely many nondegenerate critical points on a compact n-dimensional manifold with smooth density nu, then

    integral exp(i tau h) dnu=O(|tau|^(-n/2)).

Here is enough of the proof to isolate the input. A partition of unity reduces to neighborhoods of critical points and regions where dh!=0. On the latter, repeated integration by parts along a vector field V with Vh=1 gives arbitrary inverse powers of tau. Near each critical point the Morse lemma reduces the phase to a nondegenerate real quadratic form Q plus a constant, and the density to a compactly supported smooth amplitude a(u)du. Fourier inversion of a and the Gaussian Fourier identity bound that integral by C|tau|^(-n/2) times the L1 norm of the Schwartz Fourier transform of a. The identity can be regularized by exp(-epsilon|u|^2), and dominated convergence justifies the limit. No periodic Hamiltonian flow is a hypothesis.

A retained primary proof is Bates--Weinstein, Lectures on the Geometry of Quantization, Appendix B, printed/PDF pages 121--124: Lemmas B.1--B.3 and the complete quadratic-Fourier proof and final stationary-phase principle. The needed Morse lemma is also proved there in its parametrized form, Section 4.3 printed 51--52; its entire deformation proof has been read. Only the decay bound is used. There is no import of the more specialized exact global stationary-phase formula in DH1982 Section 4.

Apply the bound to h=phi and the normalized Liouville probability on the 2D-dimensional M. If mu=phi_*nu, its Fourier transform obeys |hatmu(tau)|<=C(1+|tau|)^(-D). For D>=2 and every 0<=j<=D-2, tau^j hatmu(tau) is integrable. Fourier inversion and differentiation under the integral therefore give a GLOBAL density rho in C^(D-2)(R). It is nonnegative and is zero off the compact support, because it represents a positive compactly supported measure.

For D=1 we do not claim this continuity. In EVERY dimension the pushforward has no atoms: a level set minus its finitely many critical points is a regular hypersurface, covered by countably many submersion charts and hence has volume zero, and the finite critical set has volume zero as well. This fact will handle the finite knots in the polynomiality argument.

## 6. Circle approximation proves polynomiality without a circle hypothesis

Let Lattice=ker(exp:t_R to T(R)^0), a full lattice in the compact real torus Lie algebra. The rational span Q Lattice is dense in t_R. Choose y_n in Q Lattice tending to y and staying in the open regular/separating chamber chosen in Section 4, with (x,y_n)>0. Each flow exp(R y_n) is a closed circle; it need not be a Q-defined subgroup of G. This statement is purely about the compact real torus and is not confused with the rational vector y in t(Q).

The nonconstant Hamiltonian phi_n(z)=(z,y_n) has exactly the same finite critical POINT SET P_real as phi. Its critical values a_j(n)=(z_j,y_n) tend to the distinct a_j=(z_j,y). Divide the real circle by its finite action kernel if necessary. The action is nontrivial because phi_n is nonconstant, and the resulting effective circle has the same Hamiltonian component up to a harmless choice of scalar coordinate.

Apply ONLY Duistermaat--Heckman 1982, Section 3 Proposition 3.2 and Corollary 3.3, using the proved linear-variation theorem from Sections 1--2. For compact connected symplectic M of dimension 2D and effective circle action, the Liouville pushforward has a density which is polynomial of degree at most D-1 on each interval of regular values. Normalizing its total mass or changing the circle coordinate does not change that degree. Denote these pushforwards by mu_n. The actual retained source permits finite stabilizers in the regular reduction; no globally free circle action is required.

Since M is compact, phi_n tends uniformly to phi. Hence mu_n tends weakly to mu=phi_*nu. Fix a relatively compact open interval J avoiding the finite critical values a_j. For all sufficiently large n the interval J avoids every a_j(n), so mu_n on J is a polynomial density P_n of degree <=D-1. A finite-dimensional polynomial space is closed under this distributional limit. Explicitly choose a nonnegative smooth chi of compact support in J, positive on an interval, and test with chi(s)s^i for 0<=i<=D-1. The coefficient vectors of P_n are recovered from these tests by the invertible Gram matrix

    M_ij=integral chi(s)s^(i+j) ds.

Weak convergence makes the test vectors converge, so the polynomial coefficients converge. Testing any other smooth compactly supported function in J then gives mu|J=P(s)ds for the limiting polynomial. Overlapping intervals show that a single polynomial works on each whole component of R minus {a_j}.

There is no residual mass at the finitely many a_j by Section 5's no-atoms argument. Thus mu has a nonnegative, compactly supported piecewise-polynomial probability density rho with degree <=D-1, whose possible knots are precisely within this finite REAL ALGEBRAIC list. If D>=2 the polynomial pieces agree with the GLOBAL C^(D-2) density from Section 5. For D1 they give a piecewise-constant density with no additional continuity requirement. This finishes the hard density gate with no torus-localization, Harish-Chandra orbital-integral or corrected/exact DH Section 4 formula needed.

Notice that no uniform stationary-phase bound in n, coefficient formula for the P_n, or rationality of the approximating y_n in the Q-Lie structure is assumed. Only weak convergence and closure of a fixed finite-dimensional polynomial space are used. This is why irrational periods of the original rational y cause no gap.

## 7. Mean zero and an actual height-zero fiber

The vector integral integral_M z dnu is invariant under K. Zariski density makes its invariant subspace the Lie center of g_R, which is zero because G is semisimple. Thus lambda_R(phi)=0. On the other hand phi(x)=(x,y)>0. The Haar probability has full support on M, so phi must also take a negative value; a continuous everywhere nonnegative function with zero integral would be identically zero. Since M is connected, phi takes zero at an actual point of M.

The fiber X intersect {phi=0} is a Q-defined finite-type scheme and is nonempty after real extension. It is therefore geometrically nonempty and has a Qbar-point. After any characteristic-zero extension one may embed Qbar into its algebraic closure, giving the same needed geometric nonempty fiber. This proves the ACTUAL point gate for simple-divisor and Laurent-support arguments, not just that zero belongs to a proposed interval of abstract moments. The optional extra genericity in Section 4 even makes zero a regular value, but that strengthening is unnecessary.

## 8. Application of the separately stated algebraic-knot criterion

The preceding arguments supply EVERY hypothesis of [the algebraic-spline amplification criterion](ALGEBRAIC_SPLINE_ENDPOINT_AMPLIFICATION_COUNTEREXAMPLES.md): actual rational affine algebra and Q-valued normalized functional, scalar compatibility, nonnegative probability density of degree at most D-1 with algebraic knots and global C^(D-2) regularity, smooth geometric integrality, and actual phi=0 point.

Consequently, applying that separate endpoint-amplification theorem, there exist f in Q[u] with f(0)=1 and an integer L>=max(D,4) such that, with T=1+t,

    p=T t^-2 f(T phi)^2,
    F=p^L(1+t^L p), H_marker=t^-D

are actual regular Laurent functions on X times split Gm. Every unmarked moment of F is zero and the SAME fixed unit marker survives EVERY positive order with a strictly positive RATIONAL value. Its actual smooth-point divisor proof gives absolute nonpower even up to units, and its phi=0 specialization gives adjacent Laurent support. No rational endpoint, rational circle generator, split Chevalley Q-form or coweight hypothesis remains.

These conclusions hold on EVERY characteristic-zero basechange of this GIVEN rational model because all functions and the normalized trace are defined over Q. The choice of f and L may depend on the actual algebraic spline; no universal sharp L or closed numerical endpoint coefficient is asserted. The endpoint-amplification theorem is a distinct mathematical premise; no density or actual-model hypothesis is inferred from that conditional theorem alone.

The compact factor is Mathieu over EVERY formally real field. Indeed on the real compact orbit lambda(a^2)>0 for every nonzero regular real a, by dense orbit pullback and positive Haar mass of a nonzero open set. Every finite rational Gram matrix is positive definite. Rational LDL and the positive-rational-as-sum-of-squares argument then make that Gram form anisotropic over every formally real extension, without embedding that whole field in R. Only zero is pure. The split scalar Laurent factor is Mathieu by the exact existing scalar DvK theorem and coefficient-field argument. Thus the displayed all-field pair gives genuine two-factor instability over every formally real field.

## 9. Homogeneous sources with a full-stabilizer-invariant vector

Let H be ANY closed Q-subgroup of G with a nonzero vector x in g(Q) fixed by the FULL H, not just by its identity component. Then H is contained in C_G(x), so there is the actual map

    G/H to G/C_G(x)=X.

The identity component of H is reductive. Its real unipotent radical would be a positive-dimensional real unipotent subgroup with unbounded real points inside compact G(R), impossible in characteristic zero. Thus Matsushima and affine descent give the affine source G/H; alternatively one may take its affineness as the usual homogeneous-source hypothesis of HQP. All quotients are smooth geometrically integral and the map is smooth surjective, commuting with all characteristic-zero scalar extensions.

Apply [the homogeneous-quotient transport theorem](HOMOGENEOUS_QUOTIENT_TRACE_AND_COUNTEREXAMPLE_PULLBACK.md) to the GIVEN downstairs pair of Section 8. Native trace restriction preserves every pure and marked value; the scalar Laurent unit remains the same unit. Injective pullback preserves literal Laurent support. For absolute nonpower, use the downstairs ACTUAL order-one divisor supplied by the smooth phi=0 proof, whose smooth-surjective pullback is a nonempty reduced Cartier divisor and hence has an order-one prime component. This is the required stronger geometric certificate; an arbitrary assertion of unit-nonpower without such a certificate is not substituted.

The resulting all-field product counterexamples therefore hold for every such rational homogeneous source. Its compact-factor Mathieu property over every formally real field follows directly from compact Haar positivity and rational Gram, independently of any general global-Mathieu descent theorem. In particular the full G itself is included by H=1. A disconnected stabilizer which negates all candidate vectors may have g(Q)^H=0 and is not covered by this gate; no phase averaging is proposed to bypass that restriction.

### 9.1. A separate rational-line and full-normalizer extension

**Corollary 2.** The conclusion also holds for ANY closed rational subgroup H preserving a nonzero rational LINE Qx under the FULL adjoint action, not merely fixing x. The argument uses the amplification theorem's EVEN polynomial and the actual line-stabilizer quotient. In particular it applies to the full normalizer of C_G(x) whenever its Lie center is one dimensional. No symmetry of the density or rational representative of a sign coset is required.

Let G/Q be connected semisimple with compact G(R), let0!=x in Lie(G)(Q), and C=C_G(x). Use the actual smooth affine coadjoint orbit X=G/C and rational generic height phi from Sections 1--7 above. In particular C is the FULL connected reductive centralizer, the native Q-trace realizes the required algebraic spline law, and phi=0 has an actual geometric point.

Define N_line to be the closed stabilizer in G of the rational line Qx under the adjoint representation. Every a in N_line satisfies

    Ad(a)x=chi(a)x.

This scalar is a Q-group homomorphism chi:N_line->Gm. The nonzero rational Killing norm (x,x) gives chi(a)^2=1, so chi factors through mu2. Its kernel is EXACTLY C, scheme-wise. Thus N_line/C is a closed subgroup of mu2, either trivial or the constant order-two group in characteristic zero. In particular N_line^0=C and N_line/C is finite etale. This conclusion requires no rational representative of a nontrivial geometric coset.

The free right action of N_line/C on the affine X has an actual affine quotient Y=G/N_line; pi:X->Y is finite etale and surjective. Therefore Y is smooth and geometrically integral. On the actual adjoint orbit the right action sends

    Ad(g)x -> Ad(ga)x=chi(a) Ad(g)x,

so phi transforms by the SAME sign. The coordinate S=phi^2 is an ACTUAL rational regular function on Y.

Use the even f from the even-polynomial amplification theorem, and write f(u)=v(u^2) with v in Q[s]. The seed is an actual function downstairs:

    p=T*t^-2*v(T^2 S)^2.

Thus F and the scalar Laurent unit H_marker=t^-D descend coefficientwise to Y times splitGm. Native trace restriction, by G-invariance and normalization, identifies EVERY downstairs pure and marked moment with the upstairs one. This proves the all-order counterexample on Y after EVERY characteristic-zero extension of this given Q-model.

The strong geometric certificates can also be proved directly downstairs. The actual phi0 point upstairs maps to S0 downstairs. At that point p=T/t^2 and the modifier is Q_L(t)=1+t^(L-2)(1+t), squarefree with no roots0,-1. Its nonzero t derivative on the smooth actual quotient gives a prime divisor of order1 with p a unit. The same specialization yields coefficients1,L at adjacent exponents -2L,-2L+1. Hence the downstairs phase is absolutely not a unit times a proper power and is radially primitive.

Now let H be ANY closed Q-subgroup of G preserving Qx as a FULL group scheme. Then H is contained in N_line. Since G(R) is compact, the identity component of H is reductive: a nontrivial real unipotent radical would have unbounded real points in compact G(R). The actual affine homogeneous source G/H exists, and the map

    G/H -> G/N_line

is smooth and surjective. [The homogeneous-quotient theorem](HOMOGENEOUS_QUOTIENT_TRACE_AND_COUNTEREXAMPLE_PULLBACK.md) pulls back the GIVEN downstairs pair, preserving every pure and marked value, the same scalar Laurent unit and literal Laurent support. The downstairs ACTUAL simple divisor pulls back smoothly to an order-one prime component, so the absolute unit-nonpower certificate is preserved as well.

For every such actual quotient or source, compact Haar positivity and the rational Gram/SOS argument make its separate compact factor Mathieu over every formally real field. The split scalar torus factor is Mathieu by the separate scalar Laurent theorem. Their product fails by the displayed pair over every formally real field. The pair itself works after every characteristic-zero extension of the given Q-model.

A useful special case is the FULL normalizer N_G(C) when the center of Lie(C) is one dimensional. That center is Qx, so every normalizer preserves Qx and N_G(C) is contained in N_line. Conversely every element taking x to +x or -x normalizes C_G(x), so N_line=N_G(C). Thus the preceding result covers those full normalizer quotients, including disconnected sign stabilizers, with no additional model-identification or sign-representative assumption.

This corollary requires actual preservation of the rational line. It does not cover arbitrary finite normalizers acting without such a line on a higher-dimensional centralizer center, arbitrary phases averaged to become invariant, or nonpower transport without an actual order-one divisor.

## 10. A fully explicit nonperiodic worked example

Let F=Q(sqrt(2)), let U=SL1((-1,-1)_F), and let G=Res_(F/Q)U. Its real group is SU2 times SU2. The coadjoint orbit of x=(1+sqrt(2))*i, using the rational invariant pairing -Killing/8, is the actual rational affine variety

X: a0^2+2a1^2+b0^2+2b1^2+c0^2+2c1^2=3,
   a0*a1+b0*b1+c0*c1=1.

Indeed a=a0+a1sqrt(2),b=b0+b1sqrt(2),c=c0+c1sqrt(2) satisfy a^2+b^2+c^2=(1+sqrt(2))^2. After algebraic closure this is a product of two smooth affine quadrics, so X is smooth geometrically integral. Its compact real orbit is the product of two spheres and is Zariski dense. The normalized Reynolds trace lambda is the product rotation probability. All further statements extend to every characteristic-zero field through this Q-model.

Take y=(1+2sqrt(2))*i and the actual height z=(point,y)=Tr_(F/Q)((1+2sqrt(2))*a)=2a0+8a1. Under the real embeddings it has the law

z=c_plus U1+c_minus U2, c_plus=5+3sqrt(2),c_minus=5-3sqrt(2),

where U1,U2 are independent uniform[-1,1]. Both c are positive. Convolution gives

rho(s)=(5-3sqrt(2))/14 for |s|<=6sqrt(2),
rho(s)=(10-|s|)/28 for 6sqrt(2)<=|s|<=10,
rho(s)=0 otherwise.

This density is globally continuous and degree-one piecewise polynomial. Its derivative jumps are +1/28 at +/-10 and -1/28 at +/-6sqrt(2). Every odd moment is zero, and for j>=0,

lambda(z^(2j))=[100^(j+1)-72^(j+1)]/[14(2j+1)(2j+2)].

This follows either by direct convolution or two integrations by parts using those four jumps. In particular lambda(1)=1 and lambda(z^2)=86/3. The values are rational despite the irrational inner knots.

Set T=1+t and

f(u)=(1-u^2/100)(1-u^2/72),
p=T*t^(-2)*f(Tz)^2,
Psi=p^4(1+t^4*p), H=t^(-2).

The exact spline-endpoint identity applies with D=2. For its two outer knots, a f'(a)=7/9; their total B=c*a^2 is50/7. For its two inner knots, a f'(a)=-14/25; their total B is-36/7. Thus their squared derivative factors are beta_outer=49/81 and beta_inner=196/625, respectively. Since49/81>196/625 and50>36, the signed endpoint sum is positive for every integer k>=0; the coefficient identity itself is used only at k>=D=2.

The all-order baseline/remainder argument, with L=4, therefore proves for EVERY integer m>=1

Lambda(Psi^m)=0,
Lambda(H Psi^m)=[50(49/81)^(4m)-36(196/625)^(4m)]/[7(8m+1)(8m+2)]>0.

These are exact rational numbers; the displayed same unit H survives every order after every characteristic-zero scalar extension. No extra amplification is needed in this example (M=0 already works).

There is an ACTUAL rational height-zero point (a0,a1,b0,b1,c0,c1)=(0,0,1,1,0,0). At it p=T/t^2 and1+t^4p=1+t^2+t^3. This modifier is squarefree: its derivative is t(2+3t), and its only nonzero candidate repeated root -2/3 gives value31/27. It has no root0 or-1. The nonzero t-derivative at any root on this actual smooth fiber gives a geometric prime divisor of order1, with p a unit, so Psi is absolutely not a unit times a proper power. Also Psi specializes to t^-8 T^4+t^-6 T^5; its coefficients at -8,-7 are1,4. Thus its actual Laurent support is primitive.

Over every formally real field the compact X factor is Mathieu by dense compact Haar square-positivity and rational Gram/LDL, while the split scalar torus is Mathieu by the scalar Laurent theorem. Their product fails via this explicit pair. These claims concern the stated rational model and its scalar extensions.

Finally neither x nor y generates a real circle: their component speed ratios are (1+sqrt(2))/(1-sqrt(2)) and (1+2sqrt(2))/(1-2sqrt(2)), both irrational. The convolution proof requires no periodicity. This example illustrates the broader algebraic-knot construction with a closed all-order marked formula; it makes no historical-priority or assurance claim.

## 11. Exact source scopes and assurance boundary

The source roles are deliberately distinct.

* Milne, Algebraic Groups (2017), retained canonical SHA256 `503f1a7d511632af05541c52fb2bf6fc36f5a75b0323f2f9110f03dcce11fd67`: 17.76 and proof give the smooth Lie centralizer; 17.82--17.83 and proofs give maximal tori over the base field; 17.90 and its COMPLETE proof give x in Lie(T) and maximality in G. Torus centralizer connectedness uses 17.37--17.40 with proof, and reductivity uses 17.59 with proof. Underlying Borel/unipotent results cited in those proofs remain classical inputs. This is not a whole-book audit. Normalized rational Reynolds and scalar compatibility retain only the exact existing complete-reducibility scopes from SCCP/RHT.
* Duistermaat--Heckman (1982), retained [original scan](sources/duistermaat-heckman-1982-original.pdf), SHA256 `e314c2033972ab39d244b8bc98ddfb37d8ade96c8c36bfa95e05d82fc1033be5`: actual Section 3 printed 264--266, especially Proposition 3.2 and Corollary 3.3 and their proofs, together with the complete Theorem 1.1 variation proof printed 262--264 previously inspected. Only actual effective REAL circle approximants use this result. No original Section 4 global stationary-phase assertion, nonisolated stationary formula or claim that rational y generates a circle is imported.
* Bates--Weinstein, Lectures on the Geometry of Quantization, retained [author PDF](sources/bates-weinstein-geometry-quantization.pdf), SHA256 `8596ed6391af5aaa34a06e1a9eddd7c21416f05fede0dff0d1d5ce8624ba45eb`: whole Appendix B printed 121--124, including B.1--B.3, full quadratic-Fourier argument and final stationary-phase principle, plus entire parametrized Morse proof in Section 4.3 printed 51--52. The only imported conclusion is the ordinary nondegenerate oscillatory DECAY BOUND; its shorter proof is written above. No Hamiltonian integrality, exactness or localization result is attributed to this bound.
* The closed-orbit argument, generic rational height choice, Hessian calculation, circle-approximation distributional limit, Fourier regularity and actual zero-fiber argument are written here. The invariant homogeneous trace and the rational Gram/DvK/HQP facts have their already canonical exact scopes. The separate algebraic-spline amplification theorem is the explicit source for constructing f,L and the final every-order phase.

These arguments concern ANY compact rational semisimple group and its actual rational orbits, then every characteristic-zero basechange of those models. They do not assert a classification of all algebraic varieties or all rational forms of a fixed complex homogeneous space; the fixed-vector source branch requires full-H invariance of a nonzero rational Lie vector, while the separate Corollary 2 requires preservation of its rational line. No novelty follows from the primary searches or from failure to locate an equivalent earlier formulation. External review, formal verification, independent replay and worldwide priority remain false.
