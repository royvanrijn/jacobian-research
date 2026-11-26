# Compatible Witt lifts of odd-prime BGV plane maps cannot algebraize polynomially

For every odd prime p, a published family of noninvertible plane maps has
compatible polynomial determinant-one lifts over Z/p^n for **every n**.
A radial normalization constructs this tower uniformly. Its restricted
p-adic limit is algebraic but is not polynomial; consequently the degrees
of its polynomial Witt reductions are unbounded.

For the prime-three seed, the minimum polynomial degree of a
determinant-one lift is \(81\) both modulo \(9\) and modulo \(27\).
Modulo \(81\), the minimum within the Frobenius shifts of the radial
first lift is \(109\). The degree-\(109\) witness gives the same upper
bound over all first lifts. For any prime, polynomial
lifts of a fixed nonautomorphic residue map modulo every \(p^n\), with
one common degree cap, would already give a characteristic-zero
counterexample.

The family is due to Borisov, Gabber and Vasiu. The normalization is an
explicit Hensel and Euler construction; no novelty or major breakthrough
is claimed. It supplies a precise failure of the inference from all finite
Witt lifts to a polynomial characteristic-zero lift. It does not resolve
JC2 or prohibit other lifting architectures.

## 1. Published input and the actual unresolved target

For a prime p, in F_p[x,y] put

    z = x + x^p y,
    P = z + (xz)^p,
    Q = y + x^(p(p-1)) z^(p(p-1)+1).                     (1)

Theorem 28.6(1), (2.e) of A. Borisov, O. Gabber and A. Vasiu,
[*On endomorphisms of affine spaces and the Jacobian problem*,
arXiv:2609.05746v1](https://arxiv.org/html/2609.05746v1), supplies:

* Jac(P,Q)=1;
* geometric degree p^2-p+1;
* surjectivity over an algebraic closure when p is odd.

The specialization is l=m=r=j=1, s=q=p-1, delta=-1,
h(t,w)=tw, g(t)=(-t)^(p-1). All hypotheses of (2.e) hold:
q=ms, s=pm-1, and m=1 <= (p-1)l. The statement is on PDF pages
131-132; the proof, including the degree calculation and surjectivity
criterion, is on PDF pages 132-133. The
[retained primary PDF](sources/bgv-2609.05746v1.pdf) is this version.
The older author-hosted BGV160.pdf has a different statement and does not
supply this input.

The Jacobian assertion can also be read directly. Write a=(xz)^(p-1).
In characteristic p,

    dP=dz,       dQ=dy+a^p dz,       dz=dx+x^p dy,

so dP wedge dQ=dx wedge dy. Degree and surjectivity are imported from the
published preprint, not certified by the checker here.

JC2 asks whether, for every characteristic-zero field k and all
P,Q in k[x,y] with Jac(P,Q) in k*, the pair has a polynomial inverse.

### Integral inversion and a uniform degree cap

Let \(F\) be a polynomial map over \(\mathbb Z_p\) whose Jacobian
has unit determinant at the origin. Subtract \(F(0)\) from its
target. The inverse linear part is integral, and the homogeneous
degree recursion for the formal inverse only uses that integral
linear inverse, addition and multiplication. If a polynomial inverse
exists over \(\mathbb Q_p\) or any field extension, it agrees term by
term with this unique integral formal inverse. Reducing the polynomial
inverse and its two composition identities modulo \(p\) would invert
the residue map. The BGV seed has residue degree \(p^2-p+1>1\), so
every polynomial determinant-one lift over \(\mathbb Z_p\) is
noninvertible even over \(\overline{\mathbb Q}_p\).
Its finitely generated coefficient field over \(\mathbb Q\) embeds
in \(\mathbb C\), and uniqueness of the formal inverse over that
field also precludes a polynomial inverse after the embedding.

Fix an integer cap \(B\) for the total degree of each coordinate.
All pairs of degree at most \(B\) over \(\mathbb Z_p\) form the
compact coefficient space \(\mathbb Z_p^{2\binom{B+2}{2}}\).
The residue condition and determinant-one congruence modulo \(p^n\)
cut out a closed subset of this space. Every admissible pair
over \(\mathbb Z/p^n\) has an integral coefficient representative
in that subset. If there is a solution at each length \(n\), these
closed subsets are nonempty and nested, so their intersection
contains an integral polynomial lift. Its determinant is identically
one by \(p\)-adic separatedness. Equivalently, the finite tree of
bounded-degree lifts has an infinite compatible branch. The
individual witnesses at each length may be chosen independently.

For such a plane lift, the four-variable potential
\(\Phi(X,Y)=Y^{\mathsf T}F(X)\) has gradient
\((DF(X)^{\mathsf T}Y,F(X))\). The variable
\(DF(X)^{\mathsf T}Y\) is a polynomial change of the \(Y\) coordinates,
by the adjugate identity. The gradient is therefore polynomially
equivalent to \(F\times\operatorname{id}_{\mathbb A^2}\) and is
noninvertible. Its Hessian determinant is
\((-1)^2(\det DF)^2=1\), by exchanging the two blocks of
columns in the Hessian. This gives the conditional implication to the
four-variable constant-Hessian-determinant problem.

In particular, since the radial tower gives finite lifts at every
length for odd \(p\), excluding every polynomial \(\mathbb Z_p\) lift
of this seed is equivalent to proving that the minimum lift degrees
over all choices of corrections tend to infinity. The radial
tower's own degree growth does not decide that minimum.

### Degree growth for lifts affine along a lift of the seed line

For every degree cap \(B\), some Witt level admits no degree-\(\le B\)
determinant-one lift of (1) whose restriction to a line lifting
\(x=0\) is affine-linear. The line may vary with the Witt level.
Here a line in this chart is \(x=a+by\), with \(a,b\in
p(\mathbb Z/p^n)\), and the affine-linear restriction is a
polynomial identity
\[
 F(a+bt,t)=c+wt.
\]
Reducing modulo \(p\) gives \(c=0\) and \(w=(0,1)\) there.
Precompose by \((x,y)\mapsto (x+a+by,y)\) and postcompose by
\[
 (P,Q)\mapsto
 \bigl(w_2(P-c_1)-w_1(Q-c_2),\, (Q-c_2)/w_2\bigr).
\]
These affine changes have determinant one, reduce to the identity,
and preserve the degree cap. The new pair fixes \(x=0\) pointwise.
Thus if the affine-line class admitted bounded-degree lifts at
every level, the closed coefficient conditions for fixing \(x=0\)
would give a polynomial limit fixing that line. Over
\(\overline{\mathbb Q}_p\), Gwoździewicz's pointwise fixed-line
lemma would make the limit an automorphism, contrary to integral
inversion above. [*Injectivity on one line*, Lemma 2.3](https://arxiv.org/html/alg-geom/9305008v1)

## 2. Complete first-obstruction criterion

Let F=(f,g) in F_p[x,y]^2 have Jac(f,g)=1. Choose integer coefficient
lifts f0,g0, and define

    H = (Jac(f0,g0)-1)/p mod p.

For h=sum h_ij x^i y^j define the polynomial coefficient extractor

    C_p(h) = sum_(a,b>=0) h_(pa+p-1,pb+p-1) x^a y^b.     (2)

Then F has a polynomial determinant-one lift modulo p^2, with arbitrary
new support allowed, **if and only if C_p(H)=0**.

Here is a complete elementary proof. A polynomial two-form h dx wedge dy
is exact precisely when every coefficient selected in (2) vanishes.
If i+1 is nonzero modulo p, the monomial is the differential of
x^(i+1)y^j dy/(i+1). Otherwise, if j+1 is nonzero, use
-x^i y^(j+1) dx/(j+1). Both partial derivatives have zero coefficient
at every excluded exponent pair, proving necessity as well.

Changing the lifts by pA,pB changes H by

    Jac(A,g)+Jac(f,B),

whose two-form is d(A dg-B df). Thus the class in top polynomial de Rham
cohomology is independent of the initial lifts. Conversely suppose
H dx wedge dy is exact. Choose eta=u dx+v dy with d eta=-H dx wedge dy.
Because df,dg form a basis of polynomial one-forms, there are polynomial
A,B with eta=A dg-B df. Explicitly,

    A=-f_y u+f_x v,       B=-g_y u+g_x v.                (3)

Then f0+pA,g0+pB have determinant one modulo p^2. There is no support
bound in this argument.

The same argument applies to the extension of any *fixed* determinant-one
pair modulo p^n to modulus p^(n+1), with H the next normalized determinant
digit and (f,g) the initial reduction. It does not identify the obstructions
arising from different choices at earlier digits.

This is ordinary polynomial de Rham obstruction theory. Stronger integral
Keller equations in positive characteristic have substantial prior
literature, including S. Maubach and A. Rauf,
[*A new formulation of the Jacobian Conjecture in characteristic p*,
arXiv:1507.02946](https://arxiv.org/abs/1507.02946). No result from that
paper is needed for the elementary proof above, and no novelty is inferred
from this comparison.

## 3. Pth-power additions preserve the first class

Adding A^p to one coordinate does not change the first obstruction class,
provided the residue pair has determinant one. For literal integer lifts,
the normalized change of d f wedge d g modulo p is

    A^(p-1) dA wedge dg,

or its analogue in the other coordinate. The one-form A^(p-1)dA is closed,
so its wedge with dg is exact. Changes in both coordinates have a cross
term divisible by p^2 and give the same conclusion. This statement is
about the class, not equality of the normalized error polynomials.

In particular arbitrary coordinatewise square additions to Mondello's
characteristic-two map preserve its nonzero first class. That direction
cannot bypass its modulo-four obstruction.

## 4. Uniform odd-prime vanishing

For (1), omit the pth-power summand (xz)^p from P when calculating the
class. This leaves the pair (z,y+a^p z), where a=(xz)^(p-1). Over the
integers,

    dz wedge d(y+a^p z)
      = dz wedge dy + p a^(p-1) z dz wedge da.

Consequently its normalized error modulo p is

    H' = x^(p-1)y + x^(p^2-1) z^(p^2-p+1).             (4)

For the second term use

    Jac(z,a) = -(p-1)x^p z (xz)^(p-2)
             = x^p z (xz)^(p-2) mod p.

This proves (4) including the coefficient and sign.

The first monomial in (4) is selected by C_p only when 1=p-1,
that is p=2. Expand the second using z=x(1+x^(p-1)y).
A summand with y exponent j has x exponent

    (p^2-1)+(p^2-p+1)+(p-1)j = 2p^2-p+(p-1)j,

which is -j modulo p. If j=-1 modulo p, that x exponent is 1 modulo p,
and it can be -1 modulo p only when p=2. Therefore C_p(H')=0 for
every odd prime p. Sections 2-3 give the unrestricted conclusion:

**Proposition.** For every odd prime p, the pair (1) has a polynomial
lift over Z/p^2 with Jacobian determinant one. No degree bound is imposed.

At p=2 the extractor is instead

    C_2(H)=1+x^3+x^4y,

the previously known nonzero Mondello obstruction. The current odd-prime
calculation changes that particular obstruction; it does not solve
characteristic-zero lifting.

## 5. Exact p=3 witnesses and their scope

The checker makes a deterministic choice in Section 2: integrate x first
when possible, integrate y otherwise, and represent all F_3 coefficients
by 0,1,2. It retains the complete correction digits in the
[receipt](bgv_odd_prime_lifting_verification.json). Those data reconstruct
actual integer polynomials, not just solution flags.

The checker verifies determinant one modulo 9 and modulo 27 using exact
integer sparse arithmetic. It computes the next obstruction for the
resulting fixed modulus-27 pair. Its nonzero value excludes extension
of that exact pair to modulus 81 at any polynomial degree. Other
modulus-27 pairs can have different third obstruction classes.
There is no assertion that unrestricted lifting stops at modulus 81.

Finite controls at p=2,3,5,7 check the normalized errors, simplification
(4), and coefficient extraction, with p=2 as a required negative control.
The proposition for every odd p is the written proof in Section 4.
Sections 6-7 give a different, compatible radial tower at every Witt
length and prove that its degrees are unbounded. The canonical primitive
choice above is retained as a negative control: a poor correction choice
can create a later obstruction despite existence of another infinite tower.

### Minimum degrees modulo \(9\) and \(27\) at \(p=3\)

Use the literal integer polynomials \(P_0=z+(xz)^3\) and
\(Q_0=y+x^6z^7\), with \(z=x+x^3y\). Let
\(H=(\operatorname{Jac}(P_0,Q_0)-1)/3\bmod 3\), and denote
the residue polynomials by \(\bar P,\bar Q\). Any lift of
degree at most \(80\) is represented by
\((P_0+3A,Q_0+3B)\), where
\(A,B\in\mathbb F_3[x,y]\) have degree at most \(80\). The
determinant condition is exactly
\[
 H+\operatorname{Jac}(A,\bar Q)+\operatorname{Jac}(\bar P,B)=0.
\]

Define \(\ell\) over \(\mathbb F_3\) by the following weights on
coefficients of \(x^i y^j\). Its other weights are zero.

| \(i\) | \(j\) | Weight |
| ---: | ---: | ---: |
| 38 | 9 | 1 |
| 40 | 10 | 2 |
| 42 | 11 | 1 |
| 67 | 16 | 1 |
| 69 | 17 | 1 |
| 73 | 19 | 1 |
| 75 | 20 | 1 |
| 85 | 25 | 2 |

Here \(\ell(H)=1\). For a fixed monomial \(x^u y^v\) in
\(\bar Q\), a selected output \((i,j)\) can only receive a
contribution from the correction exponent \((s,t)=(i-u+1,j-v+1)\).
Its multiplier in \(\operatorname{Jac}(A,\bar Q)\) is \(sv-tu\)
times the coefficient of \(x^u y^v\) in \(\bar Q\).
The rule for \(\operatorname{Jac}(\bar P,B)\) uses \(\bar P\) and
the opposite sign. Summing these contributions gives the two
coefficient functionals

| Functional | Nonzero coefficient weights \((s,t,\mathrm{weight})\) |
| --- | --- |
| \(\ell(\operatorname{Jac}(A,\bar Q))\) | \((61,20,1),(65,22,2),(67,23,1),(68,16,2),(70,17,1),(71,25,1),(73,26,2),(74,19,2),(76,20,1),(86,25,1)\) |
| \(\ell(\operatorname{Jac}(\bar P,B))\) | \((65,16,1),(67,17,1),(71,19,1),(73,20,1),(83,25,2),(85,26,1)\) |

Every listed \(s+t\) is at least \(81\). Thus both correction
terms have zero \(\ell\)-value at the proposed cap, contradicting
\(\ell(H)=1\). This single functional covers every allowed support.

For the upper bound, set \(L_{ij}=-H_{ij}/(i+j+2)\) for each
monomial of \(H\), with denominators taken modulo \(3\), and take
\(L\) with coefficients \(0,1,2\). The denominators are units by
the Euler calculation of Section 6. The pair
\[
 (P_0,Q_0)\bigl((1+3L)x,(1+3L)y\bigr)\pmod 9
\]
has determinant one. Modulo \(9\), the source scaling has
Jacobian \(1+3(2+x\partial_x+y\partial_y)L\), which cancels
the determinant error \(3H\) of the original pair.
The polynomial \(L\) has degree \(47\) with top term
\(2x^{38}y^9\). Modulo \(9\), a degree-\(j\) monomial is
multiplied by \(1+3jL\). In the first coordinate the summand
\((xz)^3\) is unchanged, and the term \(x^3y\) gives degree \(51\).
The term \(x^{27}y^7\) of the second coordinate gives degree
\(81\), with coefficient \(6\) at \(x^{65}y^{16}\).
Thus the minimum degree, measured by the larger of the two
coordinate degrees over \(\mathbb Z/9\), is exactly \(81\).

The minimum remains \(81\) modulo \(27\). Take the coefficient
representatives \(0,\ldots,8\) of the displayed modulus-\(9\) pair.
The two second digits \(A_2,B_2\) retained in the certificate checker
have \(24\) and \(17\) terms. Adding \(9(A_2,B_2)\) to those
representatives gives a pair of degrees \((62,81)\), whose
entire polynomial Jacobian determinant is \(1\bmod 27\).
Reduction modulo \(9\) gives the same lower bound \(81\).

### The radial modulo-\(9\) branch at modulus \(81\)

The displayed radial modulus-\(9\) pair has no determinant-one
extension modulo \(81\) with both degrees at most \(81\). Keep its
coefficient representatives \(P,Q\) in \(0,\ldots,8\). Every extension
with this cap has the form \((P+9U,Q+9V)\) with coefficients of
\(U,V\) in \(\mathbb Z/9\) and total degrees at most \(81\). Its
determinant condition is equivalent to the linear equation over
\(\mathbb Z/9\)
\[
 \frac{\operatorname{Jac}(P,Q)-1}{9}
       +\operatorname{Jac}(U,Q)+\operatorname{Jac}(P,V)=0 .
\]
The term \(81\operatorname{Jac}(U,V)\) disappears before division.

The checker retains a dual coefficient functional over
\(\mathbb Z/9\) with twenty nonzero weights. It vanishes on both
Jacobian correction expressions for every monomial of degree at
most \(81\) in the corresponding correction variable. Its value on
the divided error is \(3\bmod9\), proving inconsistency. The calculation
allows every choice of both subsequent Witt digits over this fixed
first lift. Other modulus-\(9\) lifts can have different linear
extension equations, and the unrestricted modulus-\(81\) minimum
remains open.

The [certificate checker](verify_bgv_first_lift_degree.py) checks
the monomial-annihilation conditions and all the retained polynomial
determinants with Python 3.12 and the standard library. Run it
from the repository root with

    python3 research/jacobian/mixed-characteristic-lifting/verify_bgv_first_lift_degree.py

The two sharp bounds concern moduli \(9\) and \(27\).
Growth of the minimum degrees across the whole tower remains open.

### Sharp degree 109 for the Frobenius shifts of the radial first lift

**Proposition BGVFrobeniusFirstLiftObstruction.** At \(p=3\), the
minimum of the larger coordinate degree of a determinant-one lift
modulo \(81\) is \(109\) among lifts that reduce modulo \(9\) to a pair
of the form
\[
 (P+3A,Q+3D),
 \qquad A,D\in\mathbb F_3[x^3,y^3],
\]
where \(P,Q\) are the radial first lift above. The minimum is attained
with \(A=D=0\).

Use integral coefficient representatives of \(A,D\) with the same
Frobenius support. Their partial derivatives are multiples of \(3\).
Every capped extension of the pair has the form
\((P+3A+9U,Q+3D+9V)\). All products involving two correction
derivatives vanish modulo \(81\) because each derivative of a
coordinate correction is a multiple of \(9\). Thus the full divided
determinant equation over \(\mathbb Z/9\) is
\[
 \frac{\operatorname{Jac}(P,Q)-1}{9}
 +\operatorname{Jac}(U,Q)+\operatorname{Jac}(P,V)
 +\frac{\operatorname{Jac}(A,Q)}3
 +\frac{\operatorname{Jac}(P,D)}3=0 .
\]
All quantities are integral. At cap \(108\), the fifteen-weight
functional **FROBENIUS_CAP108_DUAL** in the certificate checker
annihilates the correction columns and has value \(3\bmod9\) on the
first term. The checker tests all \(11990\) ordinary monomial columns
from \(U,V\) and all \(1406\) divided columns from \(A,D\). The latter
are indexed by \(x^{3i}y^{3j}\) with \(i+j\le36\) in each coordinate.
This excludes every lift of degree at most \(108\) in the family.

At degree \(109\), the checker retains two correction polynomials
with \(54\) and \(39\) terms to add with multiplier \(9\) to \(P,Q\).
It checks the literal degree pair \((109,105)\), exact reduction to
the radial first lift modulo \(9\), and the entire polynomial
Jacobian determinant modulo \(81\). The upper and lower bounds give
the stated minimum. The earlier twenty-weight certificate also
annihilates the \(812\) divided columns at cap \(81\), where the
Frobenius family has \(3^{812}\) first lifts.

A first lift with different gradients after dividing its difference
from \((P,Q)\) by \(3\) and reducing modulo \(3\) lies outside this
family. In particular the unrestricted minimum modulo \(81\) lies
in \([81,109]\). The degree-\(81\) test for the remaining branches is
the all-branch criterion below.

### All degree-\(81\) branches modulo \(81\)

**Proposition BGVCap81AllBranchCriterion.** For the \(p=3\) pair (1), the existence of a lift modulo \(81\) with both coordinate degrees at most \(81\) is equivalent to the following finite condition. Fix the correction and error spaces
\[
V=(\mathbb Z[x,y]_{\le81})^2,\qquad W=\mathbb Z[x,y]_{\le160}.
\]
In this section \(\operatorname{Jac}(P,Q)\) denotes the two-coordinate determinant. Let \(L_0(U,V')=\operatorname{Jac}(U,Q_0)+\operatorname{Jac}(P_0,V')\) for the integer pair \(P_0,Q_0\) from (1). Put \(K=\ker\overline{L_0}\), \(C=\operatorname{coker}\overline{L_0}\) over \(\mathbb F_3\), and \(\pi:W/3W\to C\). Enumerate the affine space of first corrections
\[
\overline{L_0}(A,D)=-\big((\operatorname{Jac}(P_0,Q_0)-1)/3\big)\bmod 3.
\]
For each pair use its integer coefficient lifts and put \(P_u=P_0+3A\), \(Q_u=Q_0+3D\), \(e_u=(\operatorname{Jac}(P_u,Q_u)-1)/9\). Set \(L_u(U,V')=\operatorname{Jac}(U,Q_u)+\operatorname{Jac}(P_u,V')\). Project \(\bar e_u\) to \(C\). For a branch where \(\pi(\bar e_u)=0\), solve \(\bar L_0(b)=-\bar e_u\) and use any integer lift \(\widetilde b\). The criterion is the existence of one such branch with
\[
-\pi\!\left(\frac{e_u+L_u(\widetilde b)}{3}\bmod3\right)
\ \in\ \operatorname{im}\beta_u,\qquad
\beta_u:K\longrightarrow C,\qquad
\beta_u(h)=\pi\!\left(\frac{L_u(\widetilde h)}{3}\bmod3\right).
\]
The integer lift \(\widetilde h\) in this definition is arbitrary.

Indeed, a pair lifting \(P_u,Q_u\) modulo \(81\) with the cap is obtained by adding \(9(U,V')\) with \((U,V')\in V/9V\). The term quadratic in this correction has a factor \(81\). The determinant equation is therefore exactly \(e_u+L_u(U,V')=0\) modulo \(9\), with \(L_u\bmod3=\overline{L_0}\). Once the first equation is solvable, every solution modulo \(3\) is \(b+h\) for some \(h\in K\). Write \((U,V')=\widetilde b+\widetilde h+3c\bmod9\). Dividing the remaining equation by \(3\) and reducing modulo \(3\) gives the criterion. Changing \(\widetilde h\) by \(3t\) changes the divided expression by \(\overline{L_0}(t)\), which projects to zero. This also proves additivity and \(\mathbb F_3\)-linearity, including the coefficient carries in the lifted vectors. Changing the chosen \(b\) changes the class by an element of \(\operatorname{im}\beta_u\).

The map can be assembled from one matrix over \(\mathbb F_3\) and an update linear in \((A,D)\). In fact,
\[
\beta_u(h)=\pi\!\left(\frac{L_0(\widetilde h)}{3}\bmod3\right)
 +\pi\!\big(\operatorname{Jac}(h_1,D)+\operatorname{Jac}(A,h_2)\big).
\]
The degree assumptions fix \(V,W\) for every branch. The radial branch is one member of the affine space and its obstruction above leaves the existential criterion unresolved. This degree-bounded criterion does not decide the existence of a horizontal polynomial lift over \(\mathbb Z_3\).

## 6. An Euler primitive gives a compatible tower at every length

Let p be odd, R=Z_p<x,y> the p-adically complete restricted power-series
ring, and F0=(f0,g0) an integral polynomial pair with J=Jac(F0) in
1+p Z_p[x,y]. Put D=x partial_x+y partial_y. Assume that, for each
nonconstant monomial x^i y^j appearing in J, the integer i+j+2 is a
p-adic unit. The constant denominator 2 is also a unit.

Define the integral polynomial

    A = sum_(i,j) 2 J_ij/(i+j+2) x^i y^j.               (5)

Then A=1 mod p, and (2+D)A=2J. Hensel's lemma in the complete ring R
applied to

    E(T)=T^2 A(Tx,Ty)-1

gives a unique lambda in 1+pR with E(lambda)=0: E(1) is divisible by p,
and E'(1) is 2 mod p. This application can also be proved directly by
successive digit corrections, because at every stage the next coefficient
is uniquely determined by division by 2. Each digit is a polynomial;
therefore lambda belongs to the restricted completion R.

Set phi=(lambda x,lambda y). This is a restricted formal change of
variables congruent to the identity mod p. It is algebraic over
Q_p(x,y), because its multiplier satisfies the displayed polynomial
E(T). For v=D(lambda)/lambda, applying D to E(lambda)=0 gives

    (1+v)(2A+DA)(phi)=2A(phi).

The determinant of phi is lambda^2(1+v). Consequently

    Jac(F0 composed phi)
       =J(phi) lambda^2(1+v)
       =lambda^2 A(phi)=1.                              (6)

Differentiation is continuous on R, so these are identities in the
restricted ring. Reduction modulo p^n gives a polynomial pair for every
n; the pairs are compatible and retain the original reduction modulo p.
The Hensel equation proves existence at all lengths, independently of
any finite computation.

The BGV pair (1), lifted by the same integer formulas, satisfies every
hypothesis uniformly for odd p. Its P monomials have total degrees 1 or
0 modulo p; its Q monomials have degrees 1 modulo p. Differentiation and
multiplication show that J has total degrees 0 or -1 modulo p. Thus
all denominators i+j+2 are 2 or 1 modulo p. Also J=1 mod p by Section 1.

**Theorem.** For every odd prime p, the BGV plane map (1) has the
compatible determinant-one polynomial Witt lifts produced by (5)-(6)
over Z/p^n for every n>=1. Their restricted limit has algebraic coordinate
functions over Q_p(x,y).

This is a general volume-normalization principle under the explicit Euler
denominator hypothesis. It does not assert the hypothesis for arbitrary
Keller reductions. In particular p=2 is excluded both in Hensel division
and in Section 4's negative control.

## 7. An algebraic-degree obstruction forbids polynomial algebraization

This argument also rules out cancellation of the algebraic radial multiplier
inside the output coordinates.

For odd p, the highest monomials of the original integer pair F0 are

    P_top=x^(p^2+p)y^p,      Q_top=x^(p^3)y^(p^2-p+1).

Their Jacobian has leading coefficient

    (p^2+p)(p^2-p+1)-p*p^3=p.

Thus J and the Euler primitive A have the same total degree

    D=p^3+2p^2+p-1,

which is odd. Put M=D+2. In independent variables s=y/x and w=x lambda,
the radial equation becomes

    w^2 A(w,sw)=x^2.                                     (8)

The left side has odd w-degree M over Q_p(s). Therefore it is not a square
in Q_p(s,w), and x^2-w^2A(w,sw) is irreducible over that field as a
polynomial in x. Gauss's lemma first in Q_p(s)[w,x], and then with x in
the coefficient field, proves that (8) is irreducible of degree M in w
over K=Q_p(s,x)=Q_p(x,y). Primitivity in the second application holds
because the w^2 coefficient is 1 and the constant coefficient is -x^2.
Consequently

    [K(lambda):K]=M>p^2.                                 (9)

In contrast, the original integer polynomial map F0 has characteristic-zero
generic degree exactly p^2. To see this, use source coordinates

    u=xz,      z=P-u^p,
    x=u/z,     y=(z-x)/x^p.

These are field identities. The second target coordinate gives

    Q u^p=(P-u^p)^(p+1)-u(P-u^p)^(p-1)+u^(p^2)(P-u^p).  (10)

The numerator on the right has u-degree exactly p^2, with leading
coefficient -pP and constant coefficient P^(p+1). It is coprime to u^p.
The rational function in (10) therefore has degree p^2 over Q_p(P).
The standard rational-function field-degree identity follows, for example,
by clearing the parameter in N(u)-Q u^p: coprimality makes that polynomial
irreducible as a polynomial linear in the independent parameter Q, and
Gauss's lemma gives irreducibility over Q_p(P,Q). The reconstruction
identities recover both source coordinates. This proves the asserted
generic degree, without importing a characteristic-zero analogue of the
characteristic-p degree theorem.

Suppose both normalized outputs Fhat=F0 composed phi belonged to K.
Put E=Q_p(lambda x,lambda y) and K0=Q_p(Fhat). The coordinates lambda x,
lambda y are algebraically independent: the determinant of phi is a
unit modulo p and hence nonzero in characteristic zero. By (10),
[E:K0]=p^2. But

    E K=K(lambda),

so the elementary compositum inequality gives

    [K(lambda):K] <= [E:K0]=p^2,

contradicting (9). Thus the two coordinates of Fhat cannot both be
rational functions, and in particular cannot both be polynomials.

If the polynomial reductions of Fhat modulo all p^n had degrees bounded
by one D0, every coefficient of degree greater than D0 would be zero
modulo every p^n and hence zero in Z_p. That would make Fhat polynomial.
Therefore the degrees of the entire compatible radial tower are unbounded.
This is about this tower, not a lower bound for every possible lifting
architecture of the original residue map.

An independent short obstruction uses the source line. Here J(0,y)=1,
so A(0,y)=1 and the unique odd-p Hensel branch has lambda(0,y)=1.
Hence Fhat(0,y)=(0,y), Fhat(0)=0 and JFhat(0)=I. If Fhat were polynomial,
Gwozdziewicz's [*Injectivity on one line*, Theorem 1.1 and Lemma 2.3]
(https://arxiv.org/html/alg-geom/9305008v1), applied over an algebraic
closure of Q_p, would make it an automorphism. The integral inversion
argument of Section 1 excludes this for any polynomial lift of (1).
This primary theorem is already imported in the
[plane-completion result](../universal-keller-plane-completion/PROOF.md).
The algebraic-degree proof above does not depend on this additional import.

For p=3, the exact checker verifies compatible radial reductions modulo
9,27,81, their full polynomial Jacobians and the fixed line. These finite
checks control the formulas; uniform all-length existence and unboundedness
are proved above. Prime controls also verify the leading coefficient in
(10) and the leading Jacobian monomial.

No polynomial characteristic-zero counterexample, general no-lifting
result, independent implementation, formal verification, external human
review, or JC2 solution is claimed. The remaining structural task would
be an architecture permitting polynomial algebraization; increasing the
computed Witt length cannot satisfy that task.

Replay from the repository root:

    python3 research/jacobian/mixed-characteristic-lifting/verify_bgv_odd_prime_lifting.py

Use --write-receipt PATH only to retain an intentional new receipt. An
ordinary replay prints its evidence without modifying the retained one.
