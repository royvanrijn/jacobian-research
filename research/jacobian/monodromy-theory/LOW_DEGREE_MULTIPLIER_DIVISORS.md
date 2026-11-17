# Low-degree multiplier divisors in the remaining quadratic direction class

An irreducible quadratic multiplier is not a surviving smallest case.
The obstruction is stronger: in the entire remaining matched pure-square
rank4 quadratic direction class, no non-u linear factor or quadratic irreducible
factor may occur in the multiplier, whatever its total degree.

The proofs use the first two feedback residues, polynomial automorphy,
the exhaustive two-image reduction, and two elementary highest-degree
arguments. This is still an input-construction theorem, not an absolute
exceptional11 action obstruction.

## 1. The first two residues preserve quadratic factors

Let r>=2, p=(u,x), d>=1 and
\[
                     V=(u^d,u^d v(x/u)),
\]
where v is a polynomial automorphism of degree at most d.
Suppose a nonzero homogeneous F gives normalized H=g-fp=FV.
For this lemma full rank is not required.
Put B=Dv, W=v-x and Fbar=F(1,x).
The constant determinant of v gives f=-Fh, h=lambda*u^(d-1).
This follows from the determinant argument even without full rank:
detDV and detD(FV) are nonzero, so f is nonzero.

**KMFI136.** Every irreducible factor q of Fbar of degree1 or2 satisfies
\[
                              q(v)=c q                  \tag{1}
\]
for a nonzero constant c.
Here q(v) is polynomial composition; for degree2 it requires an argument
beyond the linear-coordinate irreducibility used in KMFI135.
The assertion is unbounded in the degree of F and in the multiplicity
of q. It does not assert preservation for factors of higher degree.

Write eta=dlogF. The projective resolvent identity of KMFI129--133 gives
\[
 \frac{s+T\,d\log Fbar\,(TI-B)^{-1}W}{T-d}
       =1-\frac{(T-\lambda)^r}{(T-d)\det(TI-B)},\quad s=\deg F.
                                                               \tag{2}
\]
The right side has no pole at a source divisor q=0 when T is independent:
its denominator is monic as a polynomial in T.
The separate residue of q, with its nonzero positive integer
multiplicity, therefore yields
\[
                     dq\,W=0,\quad dq\,BW=0\pmod q.       \tag{3}
\]
In particular W(q)=Lq with polynomial L.
Differentiating along W gives W^2(q) divisible by q, while
\[
 W^2(q)=D^2q(W,W)+dq(B-I)W.
\]
Thus D^2q(W,W) is divisible by q as well.
For a quadratic q, the exact finite Taylor formula is
\[
                     q(v)=q+W(q)+\tfrac12D^2q(W,W).
\]
It is divisible by q. For linear q the first two terms already suffice.
A polynomial automorphism preserves irreducibility; hence q(v)=cq
with c a nonzero constant. This proves(1).

Discrete preservation alone is not sufficient for a normalized input.
For example v=(y,x+y^2,z+x^2-y^2) preserves q=x^2+z,
but W(q)=-(x-y)^2 is not divisible by q.

## 2. A quadratic factor contradicts cyclicity in the two-image class

Return to r4,d2,V0=u^2 and a normalized full-rank input.
The nonaffine projective direction has a two-dimensional quadratic
image by KMFI132--133 and is polynomial automorphic.
Its homogeneous quadratic part Q is dominant onto its image plane:
in the two normal forms it is an invertible linear image of
(0,x^2,y^2) or (0,x^2,xy).
The former is onto that plane over algebraically closed constants;
the latter is onto a dense open subset.

Full rank makes W cyclic for B:
\[
                         \det(W,BW,B^2W)\ne0.             \tag{4}
\]
This is the exact Krylov determinant identity
detC=detB*det(W,BW,B^2W), without a triangular assumption.

**KMFI137.** No irreducible homogeneous quadratic factor may divide F
in this remaining nonaffine direction class, in any total scalar degree.
In particular an irreducible quadratic F cannot supply a candidate.

Suppose such a factor exists. Its dehomogenization q is irreducible
of degree2: a degree drop would introduce a factor u into its
homogenization. By KMFI136, q(v)=cq. Also W(q)=Lq with degL<=1.
The Taylor formula gives
\[
                   \tfrac12D^2q(W,W)=(c-1-L)q.           \tag{5}
\]
Classify the rank of the homogeneous quadratic part q2 in three variables.

If rankq2=3, comparison of degree4 in q(v)=cq gives q2(Q)=0.
Density of Q on its image plane makes q2 vanish on that plane.
A nondegenerate three-dimensional quadratic form cannot vanish on
a plane: its polar bilinear form would make that plane lie in its
one-dimensional orthogonal complement. This is impossible.

If rankq2=1, write q2=ell^2 over the algebraically closed constants.
Equation(5) says
\[
                         (\ell W)^2=(c-1-L)q.
\]
Primality of q makes ell W divisible by q.
Both have degree at most2, so ell W=a q with a constant a.
Then a^2 q=c-1-L. Comparing degrees gives a=0.
Thus ell W=0 for a nonzero constant covector ell.
Differentiating yields ell B=ell, so ell B^j W=0 for all j,
contradicting(4).

If rankq2=2, an affine change writes q=XY+aZ+b.
All components of the conjugated automorphism v remain nonconstant
polynomials of degree at most2. Preservation gives
\[
                   v_X v_Y+a v_Z+b=c(XY+aZ+b).
\]
The product has degree at most2, so v_X and v_Y must both be affine.
Equation(5) now says q divides W_X*W_Y.
The two factors are affine, and irreducible q has degree2.
Hence one W coordinate is identically zero.
Again it gives a nonzero constant row annihilating the entire Krylov
span, contrary to(4).

These ranks exhaust q2. This proves the theorem.
The argument does not rule out irreducible factors of degree>=3.

## 3. A single linear multiplier factor is also impossible

**KMFI138.** In the same remaining nonaffine matched pure-square rank4
quadratic class, F has no linear factor not proportional to u, even
when other factors are nonlinear. After removing any power of u,
every remaining irreducible factor has degree at least3.
At least one must occur, by KMFI135; in particular degF>=3.
This strengthens the completely split multiplier obstruction to
individual non-u linear factors mixed with higher-degree factors.

Suppose a non-u linear factor dehomogenizes to q=ell*x+q0.
KMFI136 gives q(v)=mu*q with mu!=0.
In either normal form v=M(x+H_A/B)+b, comparison of quadratic terms
gives ell*M*e2=ell*M*e3=0. Comparison of linear terms gives ell*M=mu*ell,
so ell*e2=ell*e3=0. Thus ell is proportional to the first coordinate row
and
\[
               M=\begin{pmatrix}a&0&0\\b&c&d\\e&f&g\end{pmatrix},
                       \quad a(cg-df)\ne0.              \tag{6}
\]
Affine constants remain arbitrary.

If d=0, the full projective map is triangular with constant diagonal.
Its genuinely nonlinear quadratic part is excluded to full rank by
KMFI129, for every multiplier. Assume d!=0.
Use coordinates X=x,Y=y,Z=z-(g/d)y and put
\[
                           r=f-gc/d .
\]
The determinant makes r nonzero.
For H_A the transformed projective vector field W has the form
\[
 W_X=(a-1)X+k_1,
\]
\[
 W_Y=bX+(c+g-1)Y+dZ+cX^2+dY^2+k_2,
\]
\[
 W_Z=(e-gb/d)X+rY-Z+rX^2+k_3'.
\]
Its trace gate for Fbar has cofactor kappa-2dY, with kappa constant.
Let m be the degree of Fbar in Y and J(X,Z) its nonzero leading
coefficient. Degree m+1 in Y gives exactly
\[
                         r\partial_Z J=-d(m+2)J.
\]
A constant derivative cannot have a nonzero constant eigenvalue on
polynomials. This contradicts d!=0.

For H_B the corresponding expressions replace dY^2 in W_Y by dXY;
W_Z is unchanged. The trace cofactor is kappa-dX.
Take the highest total degree t in the two variables Y,Z in Fbar.
The X^2 terms are translations in that plane and lower this degree,
so do not enter the equation for its top part.
In that top part compare its highest X-degree, with nonzero
coefficient J(Y,Z). The only term raising X-degree is dXY*partial_Y,
and the equation is
\[
                              Y\partial_Y J=-J.
\]
An Euler derivative has no polynomial eigenvalue -1.
This closes every matrix in(6), with all constants b and every scalar degree.

Factors proportional to u dehomogenize to units and are outside this
divisor argument. They are not excluded by the theorem.
KMFI135 rules out a pure power of u as the entire multiplier,
so at least one non-u factor of degree>=3 is required.
The necessary factorization is therefore
\[
                  F=\kappa u^e\prod_i D_i^{b_i},
       \quad e\ge0,\quad b_i>0,\quad \deg D_i\ge3,
\]
with at least one D_i, all irreducible and not proportional to u.
This is a necessary form, not an existence assertion.

## 4. Boundaries, checks and the remaining problem

[Helper](low_degree_multiplier.py) ·
[checker](verify_low_degree_multiplier.py) ·
[model](low_degree_multiplier_model.json) ·
[receipt](low_degree_multiplier_receipt.json).

The u-factor boundary is essential. Dehomogenization forgets its divisor;
it cannot be included in the nonconstant-factor residue argument.
The written theorem is unbounded in the total degree of F, not a
bounded scalar search.

The full-rank assumption is also essential. An exact two-image nonlinear
rank3 normalized control with F=x is retained below the rank4 boundary.
Automorphism preservation of a quadratic is not sufficient on its own:
another control preserves the quadratic but fails the Darboux condition.

For the remaining matched pure-square primitive quadratic direction,
one may now seek only higher-degree irreducible multiplier divisors,
possibly together with a power of u, in the two arbitrary-M,b families.
These are necessary conditions. No surviving field, requested generic
group, ordinary-affine PSL2(11) realization or absolute action obstruction
is supplied. Directions without the matched pure-square component,
higher primitive directions, nonhomogeneous inputs and the complete
faithful-transitive action spectrum remain open.

Formal verification, independent whole-proof replay and external review
are not claimed. The goal remains active.
