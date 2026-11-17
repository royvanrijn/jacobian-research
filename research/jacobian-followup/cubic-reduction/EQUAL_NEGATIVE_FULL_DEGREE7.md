# Full equal-negative degree-seven closure and the sharp diagonal-weight frontier (KCUB40)

Subsequent [KCUB41](RAW13_AFFINE_ORBIT_MINIMUM.md) closes arbitrary affine
source/output sparsification of the classified112 orbit,with sharp13 and
equality rigidity. Source-frame-open descriptions below delimit this
KCUB40 stage;they are superseded for affine transformations by KCUB41.
Nonlinear coordinates and unweighted constructions remain open.

Over every characteristic-zero field, a three-variable polynomial Keller
map of weights (1,-1,-1) and ordinary degree at most seven has no
algebraic-closure collision. After tangent normalization, its positive
output is the positive source coordinate. This closes the full 77* sector
left outside [KCUB23](../equal-negative-completion/PROOF.md), by a changed
structural proof. No old timeout is increased or relabelled, and the old
large primary certificate is neither replayed nor used.

Together with [KCUB39](ZERO_WEIGHT_CLASSIFICATION.md) and the complete
[KCUB36 census](WEIGHT_RESONANCE_CENSUS7.md), every collision through degree
seven with a nontrivial displayed diagonal G_m action has weights
(1,-1,-2), up to permutation and simultaneous sign reversal. By
[KCUB15](WEIGHTED_COMPLETE_DEGREE7_CLASSIFICATION.md) it is a diagonal
transform of the established representative. Consequently, for this
three-variable diagonal-weight form, minimum ordinary collision degree is
**seven**, and minimum raw nonlinear occurrences through degree seven is
**thirteen**. Across dimensions and all degrees, its minimum collision
dimension is **three**, by KCUB36's rank-one lemma in dimensions at most two
and the verified existing upper. These matching statements now cover every
nontrivial diagonal weight, including repeated and zero weights.

This remains a restricted form. Raw counts refer to its displayed diagonal
coordinates. General linear diagonalization preserves degree and dimension
but may change monomial counts. Thus arbitrary source-frame sparsification
of the known map, unweighted constructions and higher-degree sparsity are
still open. No global positive bound changes, no new witness is constructed,
and the unrestricted strict census remains 62 excluded / 12 remaining.
The full simplicity goal remains open.

## 1. Complete reduction to a cubic with a pure highest power

Work over an algebraic closure; exclusions then descend to every base field.
Normalize the tangent blocks and the constant at the origin. The complete
ansatz, using Y=xy and Z=xz, is

    F=(x*f(Y,Z), x^-1*A(Y,Z), x^-1*B(Y,Z)),
    f(0)=1, A(0)=B(0)=0,
    deg f<=3, deg A,deg B<=4,
    J=f*{A,B}+A*{f,B}+B*{A,f}=1.                    (1)

Here {P,Q}=P_Y*Q_Z-P_Z*Q_Y. Initially A1=Y,B1=Z; later constant output
changes leave their linear block invertible without retaining that identity
normalization. The full normalized ansatz has 33 free coefficients and is
independently reconstructed, including its full three-variable determinant.

Nonconstant f automatically gives a collision: take a zero of f with x=1;
its first output is zero, whereas x=0 maps linearly and bijectively onto
that target plane. Moreover f=0 must be smooth, since an affine singular
zero makes the entire gradient of x*f(xy,xz) vanish at x=1.

KCUB23 already excludes nonconstant f of degree at most two. It also
excludes the deg f=3 branch with both negative invariant degrees at most
three: their highest equal-degree forms are proportional, giving a lower
negative combination and a strict flag below764, as covered there. Constant
f=1 has no collision by the imported plane-fibre theorem in KCUB23; the
plane degrees are at most four. We import those stated conclusions and
their written/plane-theorem boundaries, not an unproved full111 assertion.

It remains deg f=3 and at least one negative invariant degree four. The
identity Jac(fA,fB)=f*J and the homogeneous binary common-power argument in
KCUB23 show that A4 and B4 are proportional. A constant target change makes
A4 nonzero and B4=0. Put e=deg B, so 1<=e<=3. The highest forms f3*A4 and
f3*B_e have degrees 7 and e+3, with gcd(7,e+3)=1. Their zero Jacobian
therefore makes them powers of one linear form. In particular f3 is its
cube, A4 its fourth power and B_e its e-th power.

For clarity, the binary-form fact follows by setting one variable to one:
zero Jacobian gives e*U'*V-d*U*V'=0 for degrees d,e, hence U^e/V^d is
constant. Unique factorization gives common powers; coprime degrees force
the common homogeneous form to be linear. This is a written universal
argument, not a numerical leading-coefficient search.

An invertible linear change in the negative source pair thus gives

    f=Y^3+h20*Y^2+h11*Y*Z+h02*Z^2+h10*Y+h01*Z+1,
    A in span{Y,Z,Y^2,YZ,Z^2,Y^3,Y^2Z,YZ^2,Z^3,Y^4},
    B in span{Y,Z,Y^2,YZ,Z^2,Y^3}.                    (2)

Their coefficients are free. These complete spaces have 5+10+6=21
parameters. No tangent coefficient, leading parameter-zero boundary or
raw support hypothesis is removed in the models that follow. The three
cases h02!=0; h02=0,h11!=0; h02=h11=0 are exhaustive.

## 2. The quadratic-in-Z branch: a formal pole obstruction

Suppose h02!=0. Complete the square by an affine change of Z; the zero
curve has equation U^2=Q3(Y), with Q3 cubic of nonzero leading coefficient.
This is a calculation in the quotient curve, not a claimed polynomial
source-coordinate change of the three-variable map. The algebra
K[Y,U]/(U^2-Q3) has unique representatives a(Y)+U*b(Y).

There is an injective formal Laurent embedding at infinity with
Y=t^-2 and U=t^-3*S(t), after scaling U by a nonzero field constant, where
S(0)!=0 and S is the formal square root of a power series with constant
term one. Its coefficients exist recursively in characteristic zero.
The possible pole orders of a nonconstant representative are 2*deg(a) or
3+2*deg(b). Their parities differ, so they cannot cancel. Every
nonconstant regular function therefore has a positive pole order at this
one formal place, at least two; a constant has order zero. This reasoning
uses only the explicit quadratic algebra and formal series.

On f=0 put Omega=-dY/f_Z. Then dH=delta_f(H)*Omega for every restricted
polynomial H, where delta_f={f,-}. In the completed-square coordinates,
f_Z is a nonzero constant times U, so Omega has order zero at the formal
place: both dY and U have order -3. Restricting (1) gives

    A*dB-B*dA=Omega.                                 (3)

If A,B have the same pole order, subtract their constant leading ratio to
make the orders different. The new function cannot be zero, for that
would make (3) zero. For distinct pole orders m,n>=0, the leading
coefficient of A*dB-B*dA is a nonzero constant times m-n, at order
-m-n-1. At least one order is positive, since two constants give zero.
The resulting pole contradicts order zero on the right of (3).
This excludes the entire quadratic-in-Z branch without a genus theorem
or a bounded curve classification.

## 3. A polynomial Euler-Wronskian lemma on the punctured line

We use the following all-Laurent-degree algebraic fact. If a,b in
K[t,t^-1] satisfy

    a*D(b)-b*D(a)=kappa!=0,   D=t*d/dt,

then for an integer k>=1 they lie in span{t^k,t^-k}, with an invertible
coefficient matrix. Indeed they are coprime in the Laurent UFD: a common
nonunit would divide the constant Wronskian with its square. Differentiating
by D gives a*D^2(b)=b*D^2(a). Coprimality gives a common Laurent polynomial
H with D^2(a)=H*a and D^2(b)=H*b. Neither function can be constant: if one
were constant, the other would have constant Euler derivative, which is
zero in degree zero and cannot be a nonzero constant.

Thus H!=0. Comparing highest exponents in D^2(a)=H*a gives max(H)<=0;
comparing lowest exponents gives min(H)>=0. Hence H is constant. The
coefficient equation (j^2-H)*a_j=0 shows H=k^2 for an actual nonzero integer
exponent k, and only exponents k,-k occur in either function. Nonzero
Wronskian forces the coefficient matrix invertible. This proves the lemma
without a degree bound or external ramification theorem.

The same support conclusion appears in the separate
[UAAFFINER1 proof](../universal-keller-feedback-normal-form/PROOF.md) via
Riemann-Hurwitz; that theorem is not imported here. The elementary
Euler/UFD argument above is the premise used for this result. No novelty
or priority claim is made.

## 4. The punctured-line branch: only two lift types, both impossible

Suppose h02=0,h11!=0. Linear scaling and a Z-shear preserving the highest
line normalize (over the algebraic closure)

    f=1+alpha*Y+beta*Z+YZ+Y^3.

Set

    T=Y+beta,
    U=Z+Y^2-beta*Y+alpha+beta^2,
    D=1-alpha*beta-beta^3.

This triangular quotient change has Jacobian one, and f=T*U+D. If D=0,
the zero curve is singular at T=U=0, already impossible. If D!=0 its
coordinate algebra is K[T,T^-1], with U=-D/T. The kernel of substitution
is exactly (f); the quotient is explicitly isomorphic to that Laurent
ring because T^-1=-U/D. On this curve delta_f=-T*d/dT. Equation (1) and
Section3 imply that the restrictions of A,B span opposite powers T^k,T^-k.

Here (2) sharply bounds k, including every coefficient-zero case. At T=0,
Y has no pole and Z has a simple pole. Therefore A has pole order at most
three and B at most two; collectively some negative power occurs, so k<=3.
At T=infinity, Y has leading term T and Z leading term -T^2. In A's
complete space, Z^3 is the unique monomial with exponent six; all other
monomials have exponent at most five. For k=3 its coefficient must vanish,
which leaves pole order at most two at zero, a contradiction. Thus k<=2.
For k=2, Z^2 is the unique monomial of B's space with exponent four at
infinity. Its coefficient must vanish; B now has pole order at most one
at zero and cannot contain T^-2. Its restriction is a nonzero multiple of
T^2. These bounds concern actual Laurent support, not a generic branch.

For k=1 polynomial lifts of the two powers are T and U. Differences from
these lifts are divisible by f. Degree at most four permits a linear
quotient for A; its pure Y^4 highest term rules out the Z coefficient.
Degree at most three permits only a constant quotient for B. Absorbing
constants and Y=T-beta into their free coefficients, the complete forms are

    A=a*T+b*U+f*(u*T+v),
    B=c*T+d*U+w*f.                                   (4)

No condition on their constant or tangent coefficients is needed for the
contradiction. Expand J in T,U and write Cij for its coefficient of T^iU^j.
Independent exact expansion gives

    C00=D*(D*d*u+a*d-b*c),
    C11=5*D*d*u+3*(a*d-b*c),  C22=4*d*u,
    C00-(D/3)*C11+(D^2/6)*C22=0.                    (5)

A constant nonzero J would have C11=C22=0 and C00!=0, contradicting (5).

For k=2 the complete lifts are

    A=a*T^2+b*U^2+f*(u*T+v),
    B=c*T^2+w*f.                                    (6)

U^2 has degree four with only a pure Y^4 leading term, so it lies in A's
allowed space; the same division/degree argument makes (6) complete.
Here C00=0 identically. Equivalently, at T=U=0, f=D and df=0, while dB=0;
formula (1) is zero there. Again J cannot be a nonzero constant.
The two contradictions cover every k and every free parameter, including
beta=0, a=0 and vanishing leading lift coefficients. No ideal solver,
saturation query or long-certificate replay enters the proof.

## 5. The graph and parallel-line branches

It remains h02=h11=0, so f=P3(Y)+beta*Z, with P3 cubic and f(0)=1.
If beta!=0 the curve is a polynomial graph. Restricting (1) gives a
nonzero constant polynomial Wronskian in Y. Two polynomials with such a
Wronskian have degrees at most one: subtract their leading ratio if the
degrees coincide; for different degrees m,n its noncancelling degree is
m+n-1, forcing degrees zero and one. Thus both A,B restrict to affine
functions in Y. Divisibility by f and the same top-degree bounds as before
then give

    A=a*Y+b+f*(u*Y+v),   B=c*Y+d+w*f.

The original zero constants force b=-v,d=-w. In quotient coordinates Y,W=f,
whose Jacobian is beta, these are

    A=(a+u*W)*Y+v*(W-1),   B=c*Y+w*(W-1).

Let Cij now be J's coefficient of Y^iW^j. Direct exact expansion gives

    C00=-beta*(a*w-c*v),
    C01=beta*(2*a*w-2*c*v-u*w), C02=2*beta*u*w,
    C00+C01/2+C02/4=0.                              (7)

It is incompatible with a constant nonzero J. This argument uses the
original constant constraints explicitly; it does not mistake a quotient
translation for a valid source normalization.

If beta=0, smoothness forces P3 to have three distinct roots. On each
vertical root line the restricted Wronskian is a nonzero constant in Z,
so A,B restrict to affine functions of Z. In (2), A's Z^2 coefficient has
Y-degree at most one, its Z^3 coefficient is constant, and B's Z^2
coefficient is constant. Vanishing at all three roots makes every one
zero. Hence A=a(Y)*Z+b(Y), B=c(Y)*Z+d(Y).

The coefficient of Z in J is f*(a'*c-a*c'), so a,c are proportional over
K (and not both zero). A constant output change gives c=0, a!=0. Now

    J=-a(Y)*(f(Y)*d(Y))'.                            (8)

A constant nonzero J forces both factors units. But deg f=3 makes the
nonzero derivative (f*d)' have degree at least two, a contradiction.
If d=0 the determinant is zero. This closes the final branch without a
plane theorem beyond the already imported constant-f case.

## 6. Consequences, evidence and remaining work

Sections1-5 exclude every nonconstant cubic f in the remaining full111
case. Combined with the explicit lower-degree and strict imports of
KCUB23, they exclude all nonconstant f through ordinary degree seven.
The constant-f plane-fibre conclusion then gives the stated full no-collision
result. All old inputs, scripts, logs and checkpoints are pinned in
[prior-evidence.json](equal-negative-full7/prior-evidence.json) and unchanged,
including the two forty-second timeout records and the unreplayed large
chart-I primary result.

The [independent Fraction checker](verify_equal_negative_full7.py) reconstructs
all 33 native atoms and the full source determinant, the complete 21-parameter
post-leading spaces and all three short coefficient obstructions. It also
checks the exact C* triangular normalization, every allowed Laurent basis
monomial and the extreme-pole bounds, the Euler-Wronskian differentiation
identity, eight generic opposite-power controls, the completed-square
identity, symbolic formal square-root/differential controls through order
eight, all fifteen degree-four quadratic-algebra monomials and distinct-pole
leading controls, and the generic vertical derivative-jet identity.

The existing raw13 upper map is independently reconstructed coefficientwise
from its retained input, with full determinant one, weights (1,-1,-2),
row counts (2,5,6), degree seven and distinct integer points (0,1,-8),
(-1,1,4) mapping to (0,1,0). It is not a new witness. Run
`python verify_equal_negative_full7.py` from this directory. The producer
`equal_negative_full7.py` uses SymPy for saved sparse expressions and no
solver. See [artifact index](equal-negative-full7/README.md) and
[validation](equal-negative-full7/validation.json).

Finite algebraic controls do not replace the written universal coverage,
formal-place, coprimality/Euler-support, lift-divisibility and root/unit
arguments. Those proofs are not formalized or externally reviewed.
KCUB23 retains its plane/geometric trust boundary in the imported branches;
KCUB36/38/39 retain theirs in the complete weighted census and dimension
corollary. No unrelated theorem is silently strengthened.

There are now zero remaining degree<=7 raw<=12 types in displayed diagonal
weight coordinates: the only collision class is112, with exactly13 terms.
An actual smaller raw map must leave those counted diagonal coordinates
(or exceed degree seven); a degree<=7 construction with smaller invariant
flag must have no nontrivial linear torus symmetry. Source-frame changes,
unweighted maps, nonlinear graphs and changed cubic/carrier constructions
remain routes to the original global objective. Neither the raw13 orbit
under arbitrary source/output frames nor any unrestricted simplicity
minimum is settled here.
