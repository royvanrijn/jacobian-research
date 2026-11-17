# Higher-rank return: descent, sharp escapes and marked deformations

15 November 2026. Written general proofs with selected exact controls. No independent whole-proof replay, formal verification, external review or publication-priority certificate. OP-KDESC remains open. This extends the [one-direction joint theorem](../universal-keller-joint-feedback/PROOF.md), identifies a genuine higher-rank escape from its hypotheses and gives all-degree deformations preserving the canonical filter fibres.

Unless a universal ring base is specified, work over a characteristic-zero field k. Retain N compiler blocks P_i, any finite higher profiles g_ij(U,sigma), j>=4, and first outputs

    t_i=1+x_i*y_i, q_i=t_i^2*z_i+y_i^2*(1+3*t_i),
    sigma_i=t_i*q_i.

The profiled outputs are

    Phi_i=(sigma_i,
      -y_i/2-3*x_i*q_i/2
       -sum_j (j/2)*g_ij*t_i^2*x_i^(j-2)*q_i^j,
      x_i*(5-3*t_i)-x_i^3*z_i
       -sum_j (j-2)*g_ij*(x_i*q_i)^j).

UAPARALLEL1 gives whole-source polynomial Keller maps when U is retained as identity coordinates. The source is never restricted to an inverse chart.

## 1. UAMULTIKER1: a commuting return frame in every rank

Let G:A^m->A^m be polynomial with determinant1 and let M be a constant m-by-r matrix of rank r, 1<=r<=m. Its columns are mu_j. Put

    Delta_j=(DG)^(-1)*mu_j,
    W_j=(D_P Phi)^(-1)*(D_U Phi)*Delta_j,
    D_j=Delta_j-W_j.

All are polynomial derivations. The Delta_j commute, and so do the D_j. Indeed Delta_j(G)=mu_j and D_j(Phi,G)=(0,mu_j). Their commutators annihilate the respective Keller coordinate tuples; the invertible Jacobian makes each commutator zero.

Choose an invertible constant output change L taking M to the last r coordinate columns, and write L*G=(R,t), with m-r quotient coordinates R and r time coordinates t. We do not assume det L=1; if r=m such a normalization need not preserve the specified return basis. Then

    Delta_j(t_i)=D_j(t_i)=delta_ij,    Delta_j R=0.

Let C be the FULL common polynomial kernel of the Delta_j on k[U], and let Kc be their common rational kernel on k(U). Then Kc is the relative algebraic closure of k(R) in k(U), and Frac C=Kc. For the rational assertion, the r derivations are independent over k(U) by their time matrix. A further transcendental common constant over k(R) would leave a field extension of transcendence degree r-1 but r independent derivations annihilating its base, impossible in characteristic zero. Conversely algebraic constants are annihilated. Clearing a monic minimal polynomial over k[R] and using normality of k[U] gives Frac C=Kc, as in UAGEOMFEED1. Finite generation and polynomiality of C are not assumed.

Call block i active if at least one Delta_j g_iq is nonzero. Process just those blocks, leaving the others and U unchanged; denote the resulting 3N coordinates by Z. UAPARALLEL1 makes H_I:(P,U)->(Z,U) an actual polynomial Keller map. Every original signal remains polynomially recoverable, and each D_j fixes Z. Then

    intersection_j ker(D_j on k(P,U))=Kc(Z),
    intersection_j ker(D_j on k[P,U])=C[Z].             (1.1)

### Proof of the joint kernels

Apply the exact generic-field Bertini step from UAJOINTKER1, now starting with a geometrically integral auxiliary function field of relative dimension r over Kc. For each active block the inverse equation is

    S+b*S^2+sigma*S^3+sum_q g_iq*sigma^q*S^q-c/2=0.

Some preceding horizontal derivation acts nontrivially on its coefficient sum. Extending that derivation to fix S proves the required independence of S^2 and the inverse polynomial without b,c. Stacks Lemma37.32.1 therefore gives a geometrically integral field after adjoining the independent target pair b,c. That lemma is not restricted to surface fields or one-dimensional resulting fibres. Each step adds two source and two target variables, preserving relative dimension r.

After all active blocks, the full source field is geometrically integral over Kc(Z) and the r horizontal derivations remain independent by D_j t_i=delta_ij. The transcendence-degree argument above identifies their full common rational kernel. Polynomial saturation through the already constructed etale affine map H_I then puts every polynomial common constant in k[Z,U], where the D_j act coefficientwise as Delta_j. This gives C[Z]. This proof imposes no degree bound on profiles or invariants and does not use the one-direction reciprocal-sign calculation.

For the retained cubic G=K in dimension three, C=k[R] for EVERY return rank r=1,2,3 and every constant return subspace. In output coordinates (a,b,c), the inverse equation S+bS^2+aS^3-c/2 is linear in every chosen time direction. A nonzero direction mu has nonzero coefficient mu_1*S^3+mu_2*S^2-mu_3/2. Solve for one time variable: k(U)=k(R,S,remaining r-1 times), a rational extension of relative dimension r. Its common rational constants are k(R), and smooth flat regular descent gives C=k[R].

For any auxiliary G, C=k when r=m. When r=m-1, C=k[R] as well: the single quotient polynomial has gradient generating the unit ideal because its row extends to the invertible Jacobian DG. The classical one-generator closedness theorem (Ayad, as recorded by Jedrzejewicz-Zielinski, Section1) says k[R] is algebraically closed in k[U]. The common-kernel description therefore gives exactly k[R]. No multigenerator generalized Jacobian conjecture is imported.

## 2. UABASICTIME1: all polynomial time maps that descend still factor

For an arbitrary polynomial vector selector beta(P,U) with r entries define

    T=(Phi(P,U),G(U)-M*beta(P,U)).

The exact small-matrix identity is

    det DT=det(I_r-J),    J_ij=D_j(beta_i).              (2.1)

For r>1 this is a determinant condition, not r separate first-integral conditions. UAMULTIKER1 by itself does not classify every solution of (2.1).

Say the selector descends to the processed coordinates and times when each entry is in the specified subring C[Z,t]. This is a condition to PROVE, not a consequence of (2.1). Write

    F(t;Z,C)=t-beta(Z,C,t).

On this subring the D_j are ordinary partial derivatives in t_j, so

    T is Keller with determinant c
       iff det(partial F/partial t)=c in k^times.       (2.2)

Thus these selectors are exactly polynomial families of r-dimensional Keller time maps. One must not assume these maps invertible for r>=2: taking constant parameters already includes JC(r). Invertibility is unnecessary for the following factorization.

Assume the FULL auxiliary common constant ring is polynomial, C=k[eta_1,...,eta_(m-r)]. Write R=Q(eta) and H_aux=(eta,t). The chain rule for L*G=(Q(eta),t) gives a product of two polynomial Jacobians equal to the nonzero constant det L. Hence Q and H_aux are polynomial Keller self-maps. If m=r, Q is the empty identity stage. Every descended T in (2.2) has actual whole-affine stages:

1. Process active compiler blocks, retaining U.
2. Apply H_aux, retaining Z.
3. Apply the Keller time map F, retaining Z and eta.
4. Process inactive compiler blocks using their coefficient expressions in eta and retained signals.
5. Apply Q on eta and the fixed final output change L^-1.

No intermediate variety replaces the specified affine space. If e_F is the generic degree of the parameter-retaining time stage, then

    deg T=deg G * product_i d_i * e_F.                  (2.3)

For N>=1 and deg G>1 the compiler and auxiliary stages are both nontrivial. A noninvertible time map adds a factor. It does not remove these factors, even if it were a counterexample to JC(2) or another JC(r). With N>=2 the compiler stages alone force decomposability. If C becomes polynomial after adjoining identity variables or extending constants, the same statement gives a corresponding stable obstruction, provided the selector belongs to C[Z,t] before that change. No general descent of a non-descended selector is claimed.

### Cyclic nilpotent and nonlinear time blocks are included

For coefficients p,q in C[Z], the matrix

    A=((p*q,p^2),(-q^2,-p*q))

has A^2=0. The time map F=(I-A)*t has polynomial inverse (I+A)*t, despite its two nonzero off-diagonal entries. Its apparent cycle is a descended automorphism and preserves the above factors. The same applies in every rank to nilpotent coefficient matrices with entries in the common constant ring.

A recent explicit nonlinear family from Wang-Yan, equation(1.1), gives another all-degree example. With w=y+x^2 and any univariate P, put

    F_s(x,y,z)=(x+s*P(w), y+s*z-2*s*x*P(w), z+s*P(w)^2).

For invariant parameter s, the output identity Y+X^2-s*Z=w recovers x=X-s*P(w), z=Z-s*P(w)^2, y=w-x^2. This is a polynomial time automorphism and still preserves the displayed factors when used in the descended class. Its ordinary time Jacobian must not be confused with an abstract Schur matrix on a non-descended source leaf. We do not import that paper's full nilpotent-map classification.

## 3. UAMARKEDDEF1: all-degree atomic profiles with identical marked restrictions

Let I(T) be any polynomial, and replace the cubic's higher coefficients by

    a_j=I(sigma)*h_j(sigma,parameters), j>=4,            (3.1)

with arbitrary finitely many polynomial h_j. The first output remains sigma and every other output changes by a multiple of I(sigma). The profile determinant identity UAPROFILE1/UAPARALLEL1 gives Jacobian one.

On the entire inverse image of the target divisor I(a)=0, the new map and cubic are IDENTICAL as morphisms of schemes. Equivalently their pullbacks over that target divisor coincide, including arbitrary base change. This follows by reducing modulo I(sigma); it is stronger than testing known fibre points. For every target lying in that divisor, its complete scheme fibre is unchanged. Retaining the first output prevents new preimages outside the unchanged source slice.

In particular choose

    I(T)=T*(3*T-1)*(4*T+1).

Over every Z[1/6]-algebra, the COMPLETE empty fibre at mu=(1/3,-1,2/3), singleton at0 and split triple at rho=(-1/4,0,0) are exactly those of UAFILTER1. In fact the entire target planes with first coordinates0,1/3,-1/4 have unchanged restrictions, not just those three targets. The Jacobian and ideal identities are polynomial identities over Z[1/6], so nonreduced bases are included.

For every d>=4 take only h_d=1; for d=3 take the cubic. Over characteristic-zero fields the resulting three-coordinate map G_d has generic degree d and is absolutely stably self-atomic by UAALLATOM1. Higher promoted parameters are also allowed whenever the effective top profile is nonzero. Thus fixed marked fibres, even the entire restricted map over these three planes, do not determine generic degree or stable left-right class: the G_d are pairwise stably inequivalent by degree. Generic degree can drop on parameter specializations where the top profile vanishes; it is not asserted constant over every parameter base.

### Arbitrarily high finite-order agreement does not determine degree

For every integer ell>=1 replace (3.1) by a_j=I(sigma)^ell*h_j. The same proof now identifies the TWO maps over the entire target thickening V(I(a)^ell): the defining source ideal is I(sigma)^ell for both, and every output correction is in it. This is equality of the restricted morphisms, universally under base change, including nilpotent bases. Taking only h_d=1 gives an absolute stable atom of every generic degree d>=4, with the cubic as d=3, agreeing with the cubic on this whole thickening.

Consequently no prescribed FINITE order of agreement along the three marked target planes determines generic covering degree or stable left-right class in this family. The chosen maps may depend on ell. We do not assert that two distinct polynomial maps agree on the full formal completion, or to every order simultaneously.

### Canonical trace receivers can be deformed without losing their complete fibre

In UACANFIB1 keep the first stage K_A, its trace C_A and beta=1-lambda_N(C_A). In the auxiliary block use (3.1), allowing parameters h_j to be arbitrary globally defined invariant polynomial functions of the compiler source and sigma. For example h_j=f_j(C_A,sigma) is intrinsic and glues over arbitrary permitted schemes. Form

    (P,U) -> (K_A(P), G_(P)(U)-beta(P)*mu).              (3.2)

Its Jacobian is one, by triangularity and the frozen-parameter profile identity. On the complete first-stage fibre, beta is an idempotent taking the two selector values0 and1. The auxiliary first equation at marked output0 gives sigma=beta/3. Therefore I(sigma)=0 in that fibre ring. Every profile correction vanishes, so the FULL marked scheme fibre is exactly the original Spec(A), including nonreduced base changes. This uses the whole fibre ideal, not merely vanishing at listed points.

The dimension remains3N+3, with the inherited canonical volume lift in dimension3N+4. When the effective auxiliary degree is d, the generic degree over a characteristic-zero field is d*3^N. These are flexible existence families, not improved minima. There are still actual factors

    (P,U)->(P,G_(P)(U)-beta(P)*mu), then K_A times id.

After splitting A, processing each compiler block separately gives a displayed atomic chain with degrees {d,3,...,3}; no assertion about all maximal chains over arbitrary arithmetic bases is made. Every construction (3.2) remains non-atomic.

## 4. UARETURNESC1: sharp higher-rank escapes and their limits

### Rank two can replace the cubic by an atom of another degree

For the marked-preserving degree-four map G_4 above, the first output equals K_1. Consequently

    G_4=K-e_2*(K_2-G_4,2)-e_3*(K_3-G_4,3)

is an actual rank-two constant-return selector. It is NOT a polynomial in the old output coordinates. An exact witness uses

    v0=(0,0,-1),    v1=(-1/2,3,26).

Both map under K to(-1,0,0), but G_4(v0)=(-1,0,0) and G_4(v1)=(-1,24,24). No polynomial function of the old outputs can distinguish them. Also deg G_4=4 cannot be a multiple of the old cubic degree3. The new auxiliary map is an absolute stable atom, and it retains all three filter fibres from Section3. The universal factor claim of Section2 would be false without its descent hypothesis.

This is a change of the auxiliary map, not an atomic canonical receiver or a new stable class beyond the retained profile family. Adding independent compiler blocks can still leave other factors.

### These two-dimensional leaves are not affine planes

When the retained signal sigma=a is nonzero, the entire source leaf is

    Spec k(a)[x,y,(1+x*y)^-1].                          (4.1)

Indeed t*q=a makes t=1+xy invertible, and

    z=(a-t*y^2*(1+3*t))/t^3

eliminates z. Conversely this expression satisfies the signal equation. The nonconstant unit t proves that this leaf is not affine two-space, even after identity stabilization. Since partial sigma/partial z=t^3, a determinant-one compiler restricts to

    dBout wedge dCout = t^-3 dx wedge dy.               (4.2)

Thus the rank-two example is not a plane JC(2) counterexample. Its source is the hyperbola complement with the indicated volume form; rational time coordinates are not a polynomial plane model. This distinction is essential when applying published plane or nilpotent-Jacobian results.

### Full-rank cancellation can remove the old auxiliary factor entirely

Take one compiler block with profile a_4=U_1, auxiliary G=K on three variables and return matrix M=I_3. Choose beta(U)=K(U)-U. Then

    T=(Phi_4(P;U_1),U).

This is the known absolutely stable atomic degree-four profile with identity parameters, by UAALLATOM1. The original auxiliary degree-three factor has disappeared. In the descended class (2.3) would force a degree multiple of12, so this selector does not descend. This is a positive sharpness example, not a new canonical construction or a new stable atomic class; the return does not perform the canonical fibre filter.

### The native marked-preserving upgrade plus trace filtering uses three directions

Relative to the FIXED auxiliary cubic K(U), combine a nonzero native higher-profile deformation with beta(P)*mu, where beta(P) is not identically zero and mu_1 is nonzero, as in the trace filter. Assume the higher coefficients depend polynomially on sigma and compiler parameters P. In the auxiliary inverse coordinates sigma,S the profile correction has components

    (0,-sum_j (j/2)*a_j*sigma^j*S^(j-2),
       -sum_j (j-2)*a_j*sigma^j*S^j).

If a constant linear form annihilated the combined return vector, its third coefficient would vanish by the highest S power. Its second coefficient would then vanish by the next highest power, and its first would vanish because beta*mu_1 is nonzero. Hence the combined return spans all three constant output directions. Rank three is sharp for this SPECIFIED native-profile-plus-trace-shift architecture. It is not a general lower bound excluding other rank-two feedback.

In particular (3.2) with d=4 is a fibre-preserving non-descended rank-three deformation relative to the original cubic. Its auxiliary stage has degree4; the full receiver has degree4*3^N and cannot be a composition containing the entire old receiver of degree3^(N+1). It nevertheless has the different explicit factors stated in Section3. Removing an old factor is not a primitive-monodromy or self-atomicity certificate for the full map.

## 5. Current remaining gate

Higher-rank return genuinely permits non-descended selectors; the explicit rank-two and full-rank examples demonstrate this. Common first integrals still have the uniform ring C[Z], but the higher-rank Keller condition is the determinant (2.1), not scalar affine-slice classification. Polynomial time maps built in C[Z,t], including cyclic nilpotent and nonlinear published normal forms, preserve whole-affine factors.

The next canonical atomic candidate therefore needs a non-descended higher-rank interaction retaining the exact full marked scheme, or another stated departure such as a different compiler or auxiliary constant geometry. Section3 supplies a general ideal-based method to preserve marked fibres while varying auxiliary degree. It does not prove that arbitrary perturbations in a vanishing ideal are Keller: the profile Jacobian identity and retained first output are both required. The new canonical deformation families remain decomposable. Absolute canonical atomicity, improved minima and external review/novelty remain open.
