# Relative Euler defect and ruled-chart compression beyond Euler one

PF2RULED130. 16 November 2026. Written all-degree deductions for actual
complex plane Keller maps and their genuine finite intermediates. The
full Wright/JC2 and marked trace/energy implications remain OPEN.
No novelty, formal verification or independent review is asserted.

## 1. Exact defects between ANY two genuine intermediates

Fix an actual Keller map and genuine intermediate normalizations

    M2 --h--> M1 --p1--> A2,  deg(h)=t, deg(p1)=b1,
    p2=p1*h, b2=t*b1.

The top M2 need not be the original normalization N. Let V_i be the
ENTIRE p_i-etale locus and write A_i for total moved critical support,
Eta_i for total normalized finite excess, c_i=chi(V_i). PF129 proves

    c_i=b_i-A_i-Eta_i>=1.

All critical curves are disjoint A1 axes by PF123/127. Call an upper
critical curve RESIDUAL when it does not map into the lower critical
divisor. Such a curve D lies ENTIRELY over V1. Otherwise it meets
h^-1(Crit(p1)), whose curve components are also critical for p2;
the lower quotient axes are locally Q-Cartier, so every point of this
inverse image lies on one of those curve components. This would
intersect two distinct upper critical curves. Put

    Delta=A2-t*A1=sum_(residual D) e_D*f_D,
    Eta_res=sum_(points of residual D) eta_(p2).

Include every branch image. Delta counts movement inside fixed lower
blocks, including on visible target components; it is not merely the
sum over invisible target components. Degree conservation above each
lower critical curve proves the first equality for Delta.

The normalized image of a residual D in V1 is A1: it is a finite image
of A1, so Luroth and the unique point at infinity give this assertion.
Its finite map under
p1 to the final image normalization is unramified, since p1 is etale
along the whole image. A finite unramified polynomial A1->A1 is an
isomorphism. Thus the relative ramification index/residue/local excess
of h at D are exactly e_D,f_D,eta_(p2).

Let Z=h^-1(V1). It is V2 together with the residual critical A1s, so
chi(Z)=c2+#residual_curves. The finite map Z->V1 is flat over smooth
V1 by normal surface Cohen--Macaulayness. Apply the finite Euler
integration from PF129 to this relative map, retaining source singular
points and image self-identifications. The exact result is

    t*c1-c2=Delta+Eta_res,
    Eta2=t*Eta1+Eta_res.                             (1)

The second equality also follows by substituting PF129's two absolute
Euler identities. In fact excess transforms pointwise: above any lower
critical point q, with g_p the LOCAL TOPOLOGICAL degree of h at p,

    eta_(p2,p)=g_p*eta_(p1,q),  sum_(p over q) g_p=t.

To see this, use PF127's uniformized smooth-axis cover. If k is relative
normal index and r is relative local axis degree, g_p=k*r. At smooth
axis points r=1; at quotient points r is the quotient-order ratio.
The final local degrees multiply, e2=k*e1 and q2=r*q1, so subtraction
gives the identity. g_p is NOT asserted to be the scheme length of h's
fibre at a singular lower point.

Consequently c_i/b_i decreases along field extensions, strictly exactly
when residual ramification is introduced:

    c1/b1-c2/b2=(Delta+Eta_res)/b2>=0.

Equality holds iff h^-1(V1)=V2 and V2->V1 is FINITE ETALE. Indeed a
residual critical curve costs e*f>=2; absence of all of them and purity
over smooth V1 gives the asserted equality of entire etale opens.
Along a whole tower these identities telescope with weights d/b_i.
This is conservation of the exact budget, not an extra inequality
obtained by iterating it or a proof that every defect vanishes.

## 2. A plane chart in the universal cover is enough

Return to the original N and a genuine intermediate M of degree b>1,
with relative degree t=d/b and entire etale V. PF128 gives finite
pi1(V) of order m|t. Let W->V be its universal finite etale cover.
Assume W contains an OPEN A2. No surjectivity of this chart onto V is
needed. Composing the chart, W->V and V->A2 gives an ACTUAL polynomial
plane Keller map of generic degree b*m>1. If m<t it is strictly smaller
than the original map.

For an original counterexample of MINIMAL GENERIC DEGREE, m=t. As in
PF128, the finite Z=h^-1(V) lifts algebraically to W, the lifted map
is finite of degree t/m=1, and normality makes it an isomorphism. Thus

    Z=W=U_N, h|Z universal finite etale,
    Delta=Eta_res=0, b0=t*c-1, c=chi(V).             (2)

Here U_N is the entire original etale completion and b0 counts its
omitted unramified primes. This extends the mechanism of PF128 to
ARBITRARY c when the plane-chart premise is independently verified.
It does not assert that every universal cover has a plane chart.

## 3. A uniform affine-ruling criterion for that chart

Suppose V has an affine A1-fibration rho:V->A1. For a fibre
sum_j mu_j*C_j write g=gcd_j(mu_j), its WHOLE-FIBRE multiplicity.
Assume the following precise component condition:

    EVERY fibre has SOME component of multiplicity equal to g.   (3)

This includes every reduced-fibre A1-fibration and every irreducible-
fibre affine pseudo-plane, but allows reducible multiple fibres too.
It is stronger than merely requiring gcd1 or finite pi1.

Use the standard affine-line-fibration input already adopted in
HC4MYML0: Dubouloz2002.12838v1 Section1, Proposition3 and the paragraph
after Definition4. Reduced fibre components are disjoint A1s, and a
fibration with every fibre one reduced A1 is a Zariski A1-bundle.
We reuse the explicit source/reproof limits in HC4MYML0, including its
warning that some general results in the source have proof sketches.

For the fundamental group use Miyanishi's primary Oberwolfach2005
report, printed1111/PDF39, Lemma2.3: the generic-fibre pi1 maps onto
the kernel of the orbifold base group. Since the generic fibre is A1,
pi1(V) is the free product of the cyclic whole-fibre multiplicity groups.
This is the GENERAL ruling version of the input, not only the
irreducible-fibre pseudo-plane specialization used in PF128. Finiteness
therefore allows at most one g>1, and pi1(V)=C_m where m is that g;
if none, m=1.

When m>1, shift its base value to0 and normalize the Kummer cover
u^m=rho. It is connected because C(rho) is algebraically closed in
C(V)=C(rho)(v). At every divisor, the order of rho is divisible by m,
so the normalized cover is unramified in codimension one; purity over
smooth V makes it finite etale. Its degree m makes it the universal
cover W. This does not assume flatness of an unnormalized root ring.

The induced W->A1_u has the same ordinary fibres after base change.
At u=0, each component of multiplicity mu_j=m*a_j lifts to m disjoint
A1s of multiplicity a_j, since W->V is finite etale and C_j=A1.
Condition (3) supplies a_j=1. Other fibres have g=1 and already have a
multiplicity-one component by (3). The case m=1 uses W=V directly.

Choose one multiplicity-one component in each degenerate W fibre,
and delete all other components of those finitely many fibres. The
result is an open affine morphism to A1 with EVERY fibre one reduced
A1; the generic fibre remains A1 over C(u). By the adopted structure
input it is an A1-bundle. Pic(A1)=H1(A1,O)=0, so this open is A2.
Deleting an effective Cartier divisor is an affine immersion, as in
HC4MYML0, so the affine-morphism hypothesis has been checked.

Thus (3) proves the chart premise in Section2 WITHOUT chi(V)=1.
In particular a MINIMAL counterexample has NO proper relative factor
t>1 through an entire V with a REDUCED-fibre affine A1-ruling: here
m=1 forces t=1. For a whole-fibre multiple case, minimality instead
forces t=m and the universal cyclic relative cover. No ruling or
condition (3) is assumed to exist automatically.

## 4. Full generic deck groups after this compression

In the affine-ruling case, (2) supplies a subgroup C_t of the FULL
generic rational deck group D. If the block quotient is primitive,
PF128's normalizer argument gives D=C_t and M is its fixed field,
for ANY c. No Euler-one hypothesis is needed for that argument.

Without primitivity write |D|=t*h. The free deck Euler formula in
HC4MYLA0 gives h|c. If D is noncyclic, then h>=2 and HC4MYLA1 forces

    c/h>=2, hence c>=4.                             (4)

Therefore for c=1,2,3 the FULL D is cyclic even for a NONPRIMITIVE
intermediate satisfying the chart criterion. Its order may exceed t
when c>1; no false identification of the relative and full deck groups
is made. D need not preserve the original plane chart or its boundary.

## 5. Exact volume does not supply the missing component condition

The shared HC4MYML0 already proves that an exact-volume smooth affine
surface contains no open A1-bundle over P1. It applies to EVERY entire
intermediate V, since its target coframe gives dP wedge dQ=d(P*dQ).
Thus a complete-type ruling must have some fibre whose ALL component
multiplicities are at least2. We credit this as its theorem, not a new
period obstruction proved here.

For affine-base rulings, exactness alone does NOT imply (3). Here is
one explicit countercontrol, not an actual finite Keller intermediate:

    S={x^2*z=P(x,y)}, P=(y^2-x)*((y-1)^3-x).

It is a smooth affine surface. Away from x=0 use F_z=x^2, where
F=x^2*z-P. At x=0 the two roots y=0,1 have F_x=-1,1 respectively.
The x-ruling has general fibre A1 and special fibre 2*A1 disjoint
union 3*A1. Its gcd is1 and it has NO component of multiplicity1.
The ruling fundamental-group input gives pi1(S)=1.

Its nonvanishing residue volume omega=dx wedge dy/x^2 is EXACT.
Here is an algebraic regularity certificate for a global primitive.
Put

    p0=y^2*(y-1)^3, p1=-(y^2+(y-1)^3),
    R=y*(y-1)*(4y^2-7y+2)/2,
    H=(5y-2)*(4y^2-7y+2)/2, T=(8y-7)/2,
    B=(z-1)*(2R'-2-H+x*T)-p1*T.

The following POLYNOMIAL identities hold:

    R*p0'=H*p0,
    p1*(H-R'+1)-R*p1'=T*p0,
    x*F_x*(R'-1)-R*P_y=x^2*B+(H-x*T)*F.

On D(x) set theta=-dy/x+d(R/x), so dtheta=omega. On D(F_x), the
last identity and dF=0 rewrite the SAME rational form as

    theta=(B/F_x)*dy+(R/F_x)*dz,
    omega=dy wedge dz/F_x.

These two opens cover S, proving that theta and omega are globally
regular and omega is nowhere zero. This is a single exact polynomial
certificate, not a bounded experiment or a new surface classification.
It realizes the (2,3) multiplicities in a genuine exact-volume affine
example with AFFINE base; no complete-type ruling is claimed.
No whole A2 chart, dominant etale A2
chart, closed target coordinate coframe or finite target map is claimed.

## 6. Remaining global obligation

The exact tower identities and the chart criterion reduce complete
classes of imprimitive candidates without a degree/depth census.
They do not make (3) automatic. The complete-base and affine-base
bad-fibre branches, arbitrary primitive natural A/S actions, remaining
binary affine actions, and the actual cyclic finite-target realization
remain OPEN. No failure of a selected local passport is presented as
a full closure. The shared trace/energy and empty-original-fibre gates
retain priority; external review and outreach remain deferred.
