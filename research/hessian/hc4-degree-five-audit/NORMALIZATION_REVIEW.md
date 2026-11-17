# Normalization coverage and the degree five closure

This is the repository's review of the assembly of the two version-pinned
preprints, not a new attribution of their theorem. The original proof and
programs are retained under `source-2611.00161` and `source-2611.00134`.
The programs were executed, not merely inspected or inferred from their
printed check totals. The extra complete-block ideal below was generated
independently.

The logical conclusion is injectivity and polynomial invertibility of
gradients of potentials of degree at most five. It is not an exclusion of
all constant-Hessian potentials: many such potentials are automorphic.

## Imported results and terminal arguments

We use the established Hessian theorems in dimensions two and three, the
homogeneous Gordan–Noether theorem in at most four variables, the plane
Jacobian theorem through degree four (a small part of Moh's range), and
the degree-divisibility property of plane polynomial automorphisms. These
are imports, not theorems newly reproved by the ancillary programs.
[Ni's paper](https://arxiv.org/html/2608.14217v1) states the lower-dimensional
Hessian imports and the scalar descent argument. The repository's `HC4CQ1`
already covers potentials of degree at most four.

The terminal arguments used by the case split have the following precise
coverage. Their validity does not depend on a numerical coefficient bound.

**A nonzero constant second derivative.** For
f=P(y)+sR(y)+a s^2, a nonzero constant, let t=s+R/(2a). The output f_s
determines t. The remaining outputs are the gradient in y of
P-R^2/(4a)+tR. Schur complementation gives constant nonzero Hessian
determinant for that three-variable potential for every fixed t. HC3
therefore gives injectivity on every fibre and hence injectivity of the
full gradient. A polynomial map with nonzero constant Jacobian which is
injective over C is an automorphism.

**A quadratic block with constant nonzero determinant.** Write
f(X,z)=z^T A(X)z/2+b(X)^T z+c(X), with X,z each two-dimensional and det A
a nonzero constant. The equation f_z=Y_z solves for z polynomially. The
partial Legendre potential has two-variable Hessian determinant
det Hess(f)/det A. HC2 gives injectivity for every fixed Y_z, and the
full gradient is again an automorphism. This uses the general Schur identity,
not just the concrete instance printed by one source checker.

**A missing constant direction.** For f=P(u1,u2,w)+sR(u1,u2), if P_ww is
nonzero, constancy of det Hess f makes
R2^2 R11-2 R1 R2 R12+R1^2 R22 vanish. The normal direction to every general
irreducible fibre of R is constant along that fibre. Such a component is a
line; lines from distinct polynomial fibres cannot meet, so all are parallel.
Thus R=g(l) in one constant linear form. The determinant contains g'^2 as
a factor, forcing g' to be a nonzero constant. The remaining problem is
HC2 on the fibres of that form. If P_ww=0, write P=P0+wP1. The determinant
is the square of J(R,P1); both degrees are at most four. The plane
Jacobian theorem applies and the two remaining variables follow by a
linear system with unit determinant.

For f=P(y)+sR(y)+a s^2 with deg R at most two and a=0, the coefficient of
s^2 forces the quadratic Hessian of R to have a nonzero kernel vector
orthogonal to the linear part of R. This supplies the preceding missing
direction. If that quadratic Hessian were invertible, its highest-degree
bordered determinant would already be nonzero. At ranks two, one and zero,
the kernel and orthogonality conditions give the required direction.

**Triangular leaves.** One repeatedly obtains an output equal to a nonzero
constant times a base variable. A second output determines another
variable, and the remaining equations are affine in the remaining pair
with unit determinant. These are actual inverse constructions; a
nonvanishing determinant of a leading form alone would not suffice.

## Degree and stratum coverage

Remove affine terms, which do not change the Hessian. A genuinely quintic
top has singular Hessian, since the top determinant layer must vanish.
Gordan–Noether makes it a cone. Its essential number of variables is
therefore exactly one, two or three. The rank-zero top is absent because
the degree is five. This partitions the whole quintic problem.

### Three essential top variables

The ordinary and weighted determinant layers are checked with fully
symbolic lower coefficients in the ternary verifier. The permitted
weighted top has the form phi+s h+s^2 m. The classification is made
before imposing an arbitrary specialization of phi. Nondegeneracy is
treated as an open condition: every coefficient of det Hess(phi) is
individually used in Rabinowitsch tests, so there is no restriction to
one chosen coefficient chart.

The transformations act on the three equal-weight base variables and
preserve the weighted shape. At m=0 the cubic h is zero or has singular
ternary Hessian, hence is binary. Its distinct-root, double-root and
pure-cube possibilities are all retained. At m nonzero, align m with x1.
The passive singular-Hessian condition and homogeneity put h in the form

\[
 a x_2^3+b x_1x_2^2+x_1^2(c x_2+c' x_3)+e x_1^3.
\]

Split by a nonzero, then a=0 with b nonzero, then a=b=0 with (c,c')
nonzero, and finally a=b=c=c'=0. In the second branch a shift of x2 by
a multiple of x1 removes c; in the third a passive linear change and
shift align the nonzero pair and remove e. All other coefficients of phi
remain generic. No division by a possibly zero leading coefficient is
used without retaining its zero branch.

The surviving forms end in the inverse mechanisms above, a two-fibre
plane-Jacobian endpoint, or an inconsistent completion.
In the plane endpoint the actual degrees are four and three. A plane
Keller pair in this range is an automorphism, while neither degree divides
the other; this is impossible. The potentially exceptional weighted form
is excluded with all 41 admissible lower coefficients retained. The 29
ternary checks cover the identities and ideal computations used in these
steps. The classification and the imported plane degree obstruction are
written arguments.

### Two essential top variables

The first transverse quartic Hessian has rank at most one over the base
function field. Choose its kernel vector primitive over the polynomial
base ring. Gauss's lemma makes the coefficients of powers of that primitive
linear form polynomial. Homogeneous degree four leaves either a constant
transverse direction or a degree-one primitive vector occurring quadratically.
The latter is an invertible linear pair of base coordinates; the weighted
determinant identity excludes it. This step prevents discarding a moving
quartic direction without proof.

With a constant direction chosen, the general source form is
P+sR+s^2 L+e s^3. The cubic part of R is a binary cubic and L is affine.
The relevant coordinate group preserves the binary top: GL2 on its base,
affine changes of the third variable by base linear forms, translations,
and scalings. It does not freely mix the third variable into the base.

| Flag | Exhaustive splits and admissible reductions |
| --- | --- |
| e nonzero | Scale e and remove L by an affine change of s. The cubic of R is zero, a cube, a double-root cubic, or a cubic with three distinct roots. |
| Three distinct cubic roots | Two base translations remove the two square coefficients. The other quadratic and lower coefficients remain free until forced by the equations. |
| Double-root cubic | If the x1*x3 quadratic coefficient is nonzero, normalize it and remove the other two base quadratic coefficients by allowed third-variable shears. On its zero branch, base translations remove those coefficients; the remaining square coefficient is zero or nonzero. |
| Pure cube | The transverse quadratic restriction has rank zero or one. The flag-preserving group distinguishes a square in the third variable, a square in the second variable, and cross terms with the first. The source checks retain all unnormalized linear coefficients. |
| Cubic part zero | A quadratic R with a critical point is impossible: the s row of the full Hessian vanishes there at s=0. For a submersion, rank-two, rank-one and rank-zero quadratic parts have linear parts outside their images. The flag determines whether the surviving linear direction lies inside or outside the binary base. All such cases occur in the source list. |
| e zero, L nonconstant | Its linear part lies in the binary base. Align it with x1 and remove its constant term by translation. Split in succession by the coefficients of x3^2 in R, x3^3 in P, x1*x3^2 in P, and x3^2 in P. Each zero complement is retained. |
| e zero, L constant | A nonzero value is a scalar pivot. At zero, the highest s coefficient depends on R alone, so its conclusions are valid with every transverse quartic coefficient of P present. The preliminaries include those extra coefficients where later layers need them. |

The finite lists in this table are classifications of coefficient flags,
not degree or height samples. A normalized nonzero scalar is made one by
scaling in its stabilizer over C; the complementary zero stratum is always
explicit. When a lower nonzero P coefficient is normalized, the binary top
remains arbitrary rather than being fixed to a numerical quintic.

The mixed surviving family has an explicit inverse over Q(sqrt(-3)); its
conjugate is also checked. Two additional leaves have triangular inverses,
and the linear-R leaf reduces fibrewise to HC2. The rest reach a scalar
pivot, a constant-determinant quadratic block, a missing direction, or an
inconsistent ideal. The main binary verifier has 148 checks and its
preliminaries have 19. The distinction between the restricted quartic
ansatz and the full quartic ansatz in the L=0 branch was specifically
checked: implications depending only on R use the highest s layer, and
the full P is restored before using later layers.

### One essential top variable

Translate the active variable to remove its pure fourth-power term.
Write the quartic part by descending active powers. The successive free
passive blocks have the following complete classification:

| First nonzero block | Shapes to retain |
| --- | --- |
| Passive quartic | A genuinely binary quartic: four distinct roots, one double root, two double roots, or one triple root; a fourth power; or zero. The generic even-quartic chart retains its cross-ratio parameter. |
| Passive cubic at the next active power | Two genuine binary cubic types, one cube, or zero. |
| Passive quadratic at the next active power | Rank two, rank one, or zero. |
| Passive linear part | Zero or one nonzero linear form. |

The vanishing Hessian conditions that justify these reductions are among
the checked determinant layers. Gordan–Noether excludes a genuinely ternary
passive quartic or cubic at these stages. Rank-two quadratic and cubic
forms are not treated as pure powers.

At each subsequent normalization in two passive variables, the earlier
nonlinear forms involving that pair have vanished. Thus GL2 on that pair
preserves the previous flag. Perfect-square and perfect-cube covariants
are verified before these changes; both the zero form and the nonzero
power are retained. Remaining affine partners are left free or split by
zero/nonzero after the power has been aligned. Translating the active
variable to remove a newly introduced pure fourth-power coefficient
does not change the other homogeneous quartic coefficients.

If the passive linear form is nonzero while the next cubic is genuinely
binary, the cubic's two-dimensional coefficient span either contains that
linear form or is complementary to it. These give the last two branches
of the unary tree. Lower-rank cubics fall into the containing-span branch;
they are not discarded.

The unary peeling procedure only solves equations with constant nonzero
linear pivots and imposes zeros proved by exact radical membership. Its
nonzero hypotheses are included by a product inverse. Failure to derive
more zeros would leave a failed check, not a theorem. All 80 unary checks
passed. Terminal leaves are the three inverse mechanisms above or
inconsistent systems, apart from the quadratic-passive block addressed next.

## Independent completion of the quadratic-passive block

The remaining block form is

\[
 f=x^5+P_{\le4}(x,y)+z b(x,y)+w d(x,y)
           +\tfrac12(z,w)(A_0+xA_1)(z,w)^T.
\]

Here P contains every base monomial of degrees two through four except
x^4; b and d contain every base monomial of degrees one through three;
and A0,A1 are arbitrary constant symmetric matrices. Affine terms have
no effect. This is a complete 35-coefficient class: 11+18+6 coefficients.
It is stable under GL2(z,w).

The independently generated polynomial system sets all nonconstant
coefficients of det Hess f to zero, while inverting det A1 and the constant
term of det Hess f. There are 167 equations in 35 coefficient variables
and two inverse variables. The exact Groebner basis is `[1]`. Thus no
constant nonzero full determinant is possible with rank A1=2. This result
does not require selecting a diagonal representative inside an already
normalized leaf.

If rank A1=1, a GL2 change of the full block puts it in diag(0,1).
If the complementary A0 diagonal entry is nonzero there is a scalar pivot.
If it is zero and the cross entry is nonzero, det(A0+xA1) is a nonzero
constant and the partial Legendre argument applies. If both vanish,
one passive variable is linear with a coefficient depending only on x,y,
giving the missing-direction argument. At rank A1=0 the same alternatives
follow from the constant matrix A0. These inverse arguments allow all b,d,P,
so they remain valid after a coordinate change alters earlier leaf coefficients.

This closes the normalization concern that motivated the independent
calculation. It does not turn the bulk author programs into independent
implementations. The optional attempt to generate a smaller certificate
through the unsaturated degenerate variety yielded no certificate and is
not used.

## Assembly and descent

Every quintic top has entered exactly one of the three essential-variable
branches. Every retained branch ends in a proved inverse mechanism or a
characteristic-zero inconsistent ideal. Potentials of smaller degree are
covered by HC4CQ1. Hence every potential of degree at most five has an
automorphic gradient, and no nonzero antipodal gradient collision exists.

For any characteristic-zero coefficient field, descend from a finitely
generated subfield embedded in C. The formal inverse at the origin is
unique and its coefficients lie in that subfield. Since it truncates over
C, it is already a polynomial inverse over the original field. This
does not invoke an unproved characteristic-zero lifting conjecture.

The direct-quintic collision question `OP-HC4-D5` is therefore closed by
reproduction. The abstract Fitting-module computations and all-degree
conditional theorems retained elsewhere are not declared solved or false
by this conclusion. The Meng–Yang graph degree and the degree of its
resulting potential must still be distinguished. Neither unrestricted
HC4 nor JC2 follows from this bounded-degree result.
