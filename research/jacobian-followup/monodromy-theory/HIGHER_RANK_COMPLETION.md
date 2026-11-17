# A constructive affine completion in every tautological rank

The rank-two completion theorem KMFI59 need not remain a
dimension-specific tool. Larger weights reduce its higher-jet problem
to a finite formal curve calculation and unipotent Gaussian elimination.
The resulting source change is a polynomial automorphism with an
explicit polynomial inverse. Thus a formal chart is not being used
as a substitute for ordinary affine space.

This supplies the affine completion for a hypothetical higher-rank
exceptional PSL2(11) field. It does not construct that field or decide
the full Keller action spectrum.

## 1. A finite weighted-jet rectification lemma

**KMFI103.** Let G:A^(r+1)->A^(r+1), r>=1, be polynomial with
G(0)=0 and DG(0)=I. Write its variables as (t,u1,...,ur), with weights
\[
                   e0=1,\qquad ei=r+i\quad(1<=i<=r).
                                                               \tag{1}
\]
There is an explicit triangular polynomial automorphism
\[
                   \phi(t,u)=(t,h(t)+A(t)u)                 \tag{2}
\]
with det Dphi=1, h(0)=h'(0)=0, A lower unitriangular and A(0)=I,
such that every G_i composed with phi has weighted order at least ei.
Its inverse is (t,A(t)^-1(u-h(t))). Both h and A are computed
modulo t^(2r), and their polynomial representatives are used in (2).

Here is the finite algorithm and proof.
Put D=2r. The formal solution H(t) of G_high(t,H(t))=0 starts
in t^2. Compute h=H mod t^D by D iterations
\[
 H_0=0,\qquad
 H_{j+1}=H_j-G_{\rm high}(t,H_j)\pmod {t^D}.                 \tag{3}
\]
Indeed G_high=u+O((t,u)^2), so this iteration improves the t-adic
error by at least one order. It follows that
G_high(t,h)=0 mod t^D.

Set M(t)=D_u G_high(t,h) mod t^D. It is the identity at0.
Right-column Gaussian elimination, starting with the bottom row,
finds a lower unitriangular A over k[t]/(t^D) such that M A
is upper triangular modulo t^D. The pivots are units since they
equal1 at0. Use the polynomial representatives of the entries of A;
the resulting exact polynomial matrix is still lower unitriangular,
so has determinant1 and polynomial inverse.

Expand G_high(t,h+A u) by degree in the high variables u.
Its constant part is O(t^D). For output i, a linear term in u_j
with j<i has coefficient O(t^D), while j>=i already has weight
ej>=ei. Every term of degree at least two in the high variables
has weight at least2(r+1)>2r, greater than every required output
weight. These three cases exhaust all monomials.
The first output only requires its constant term to vanish.
This proves the lemma in every r.

This is weaker than matching an arbitrary full inverse jet.
No assumption that G has constant Jacobian is needed.
Truncating an arbitrary formal automorphism need not give a
polynomial automorphism; here the special form (2), its inverse
and the weight inequalities supply the missing global argument.
The retained [formal-orbit warning](../../archive/non-elliptic/extended-geometry/FORMAL_ORBIT_TRIVIALITY.md)
does not supply this lemma and is not used to bypass algebraization.

## 2. All-rank completion theorem

Let p=(p1,...,pr), Delta=(1,p), and let
\[
 S(\gamma,p)=X(p)+\gamma\Delta,\qquad
 X\in k[p]^{r+1},\quad\det DS=\gamma^r,\quad
                         \operatorname{rank}DX=r.          \tag{4}
\]
Work over algebraically closed characteristic-zero constants.
Let L=k(gamma,p), K=k(S), G its geometric generic closure group,
H its sheet stabilizer, and D=closure(X(A^r)).

**KMFI104.** The following are equivalent, without a finiteness assumption:

1. H normally generates G.
2. For a reduced equation d of D, the polynomial
   d(S)=gamma^(r+1) h has nonconstant h.
3. Some smooth generic branch value X(p*) has a regular preimage
   (gamma0,p0) with gamma0!=0.
4. After adjoining one common independent variable C, the exact
   field pair K(C) subset L(C) has a polynomial Keller model on
   ordinary A^(r+2).

The proof gives a finite coefficient algorithm for the map.
For a finite sweep of degree N and branch degree delta,
these conditions are exactly N>(r+1)delta.
Over an original nonclosed coefficient field, the chosen points
and coefficients may require a finite algebraic extension.
No arithmetic descent or minimum-dimension statement is asserted.

The first implications generalize the written KMFI59 proof.
At the generic point of gamma0, rank DX=r makes S quasi-finite.
The Jacobian different exponent r gives ramification index r+1.
Thus d(S) has exactly gamma-adic order r+1.
Every other prime factor is etale and has multiplicity1.

If h is constant, gamma^(r+1) lies in K but gamma does not:
the source valuation has order1 on gamma and order divisible
by r+1 on every element of K. This gives a nontrivial cyclic
Galois intermediate field, contradicting1.
If h is nonconstant, a prime of h gives an actual etale tail
dominating D. Intersect its image with the generic rank-r image
of X to choose the points in3. These are algebraic choices,
not approximate numerical roots.

For a finite S, all normalization points are present and
N=(r+1)delta+u, where u is the degree of the unramified tail.
This proves the finite numerical criterion.
The implication4=>1 is the absolute Keller normal-generation theorem.

Here is the full construction for3=>4.

### Polynomial first stage from nilpotency

Put N=[-p | I_r] and P=N DX. Expanding the determinant of (4) gives
\[
                         \det(\gamma I_r+P)=\gamma^r,
\]
so P is polynomial and nilpotent.
Translate p* and X(p*) to0 and make a matching constant linear
source/target direction change so that P(0) is strictly upper
triangular. These operations preserve (4).
At0, DX_(i+1) has no linear terms in p_j with j<=i.

Choose
\[
 a_i=r+i-1,\qquad e0=1,\quad e_i=a_i+1=r+i.
\]
Substitute p_i=gamma^(a_i) u_i and define
\[
 E_i(\gamma,u)=
       S_{i+1}(\gamma,\gamma^{a_1}u_1,\ldots,\gamma^{a_r}u_r)
                        /\gamma^{e_i}\quad(0<=i<=r).        \tag{5}
\]
All E_i are polynomials.
The constants vanish after branch centering. A surviving linear
term in X_(i+1) has source index j>i and weight a_j>=e_i.
Every nonlinear term has weight at least2r>=max e_i.
The gamma p_i terms have exactly the required weight.
These cases are the complete polynomiality check.

Let theta=(gamma,u), and let theta0 be the regular preimage in3
after the scaling. Then E(theta0)=0, gamma0!=0, and
\[
                   \det DE(\theta0)=\gamma0^{-1}.           \tag{6}
\]
Indeed sum(e_i)=r+1+sum(a_i), so the sweep/scaling/rescaling
determinants give gamma^(r+sum(a_i)-sum(e_i))=gamma^-1
at a zero, where derivatives of target rescaling contribute nothing.

### An actual polynomial source change and the whole map

Put J=DE(theta0), L0=J^-1, and
G(v)=E(theta0+L0 v). Apply KMFI103 to G with the weights e_i.
Set
\[
                 \Theta(v)=\theta0+L0\phi(v).
\]
This is a polynomial automorphism with polynomial inverse and
constant Jacobian det L0=gamma0.
On ordinary coordinates (x,y0,...,yr), substitute
v_i=x^(e_i)y_i. Write gamma for the first component of Theta(v),
and define
\[
 C=x\gamma,\qquad
 F=(C,E_0(\Theta(v))/x^{e0},\ldots,E_r(\Theta(v))/x^{er}).
                                                               \tag{7}
\]
Every quotient is polynomial by KMFI103, including on x0.

For the whole Jacobian, first pass from (x,y) to (x,theta).
Its determinant is x^(sum e_i)*det L0.
Next use (C,theta), whose determinant contribution is gamma.
At fixed C the remaining outputs equal
S_i(gamma,gamma^a u)/C^(e_i), so their determinant is
gamma^(r+sum a_i)/C^(sum e_i).
Consequently
\[
 \det DF=x^{\sum e_i}\det L0\,
       \gamma\,\gamma^{r+\sum a_i}/C^{\sum e_i}
                    =\det L0=\gamma0\ne0.                 \tag{8}
\]
This is an exact identity on a dense open and hence globally,
because all outputs in (7) are polynomials.

Finally the target equations S_i=C^(e_i)F_(i+2) recover K(C).
Conversely x=C/gamma, u_i=p_i/gamma^(a_i), the explicit inverse
of Theta, and y_i=v_i/x^(e_i) recover L(C).
Thus the exact field pair, degree and geometric action are preserved.

The construction3=>4 only needs a centered X with nilpotent P(0)
and a regular second preimage. It also works for rank-deficient X
when such points are explicitly supplied; the equivalence with1
and2 above is asserted under the full-rank hypothesis (4).

## 3. The remaining exceptional11 gate has one precise rank

**KMFI105.** Within full-rank normalized polynomial tautological
sweeps (4) of arbitrary rank, either exceptional degree11 PSL2(11)
field action can occur only for r4.
Any such rank-four field would receive an actual ordinary-A6
Keller completion by KMFI104. It must be nonfinite and satisfy
\[
 [k(p):k(X)]=1,\qquad
                     d(S)=\gamma^5 h,                    \tag{9}
\]
where h is irreducible and its divisor maps birationally to D.
There is exactly one additional index5 normalization prime over D
absent from the polynomial sweep source. D must be covered by
horizontal target straight lines.

Here are the action and geometry arguments.
The nontrivial A5 stabilizer normally generates the simple group,
so KMFI104 applies to any proposed field (4).
In its completion, D times G_m remains an affine target divisor
under the target rescaling S_i=C^(e_i)F_(i+2).
Full affine inertia must fix a sheet, so its orders can only be
2,3 or5 by KMFI82. The visible source cycle has length r+1,
forcing r in{1,2,4}.

For r1, write X=(f(w),g(w)); the determinant condition gives
g'=w f'. Eliminating gamma gives the generic linear pencil
Y2=wY1+g-wf, with full natural symmetric group, by the
birational critical-value calculation of KMFI38.
For r2, KMFI102 excludes both exceptional actions.
Thus r4 is the only remaining rank in this full-rank architecture.
This is not a lower bound on the dimension of arbitrary Keller maps.

The rank-four visible prime has index5. Its full inertia must
therefore be C5. Over algebraically closed constants decomposition
centralizes tame inertia, and C_G(C5)=C5.
KMFI82 consequently gives exactly two ramified primes of indices
5 and residue degrees1, plus one unramified prime of residue degree1.
The single visible ramification divisor gamma0 therefore has
residue degree **one**, not two: X is birational onto D.
The other ramified prime must be missing, since every source
ramification prime is gamma0. Hence S cannot be finite.
The unique unramified prime must be present by KMFI104, and is
defined by the squarefree irreducible h in (9).

Compactify the line family just as in KMFI101, now over P4
directions. Each missing normalization divisor comes from the
boundary of the universal projective line family. Every limiting
line with an affine part has direction (0,0,p1,...,p4), so is
horizontal. Its boundary line family sweeps D.
This retains all exceptional graph-boundary components.
The rank-two horizontal-ruling rigidity KMFI100 does not extend
automatically: the generic Y1 fiber now has dimension3, and a
line ruling need not make it a line. No such higher-rank rigidity
is claimed.

## 4. Concrete controls and the remaining field problem

The [checker](verify_higher_completion.py) retains:

- nontrivial finite formal-curve and lower-unitriangular matrix
  calculations through rank4, with exact inverse and Jacobian;
- strictly upper triangular companion nilpotent surfaces through
  rank4, whose group is cyclic and which correctly fail the tail gate;
- a fully expanded six-variable Keller map from (7), with whole
  Jacobian -3, defined by rational coefficients and exact divisions.
  Its geometric action is the already-known S3, not PSL2(11).

For the last control, use r4 and
\[
 X=(2p1+p1^2,\ p1^2+2p1^3/3,\ p3,\ p4,\ 0).
\]
Its sweep has determinant gamma^4 and the regular zero
(gamma,p)=(-3,-3,0,0,0).
With weights a=(4,5,6,7), the stage is
\[
 E=(1+2\gamma^3u1+\gamma^7u1^2,\
 u1+\gamma^3u1^2+2\gamma^7u1^3/3,\
 u2+u3,\ u3+u4,\ u4).
\]
The marked theta0 is (-3,-1/27,0,0,0), and its inverse Jacobian is
\[
 L0=\begin{pmatrix}
 9&324&0&0&0\\
 -11/27&-15&0&0&0\\
 0&0&1&-1&1\\
 0&0&0&1&-1\\
 0&0&0&0&1
 \end{pmatrix}.
\]
Equations (2)--(3) and (7) specify all coefficients finitely.
The retained model also serializes every expanded coordinate.
X has rank3; this control verifies the marked-point construction,
and is not falsely treated as a full-rank index5 candidate.
The sweep target equations Y3=gamma*p2+p3,Y4=gamma*p3+p4,Y5=gamma*p4 recover p2,p3,p4 rationally once gamma is known. Its root parameter p1 satisfies p1^3+3p1^2-3Y1*p1+3Y2=0, the generic cubic pencil with nonsquare discriminant. Thus the exact field group is S3. The last three outputs are identity spectators. The first three
have degrees15,125,135 and respectively9,323,388 monomials.
An exact polynomial-ring determinant of the entire expanded map
is -3, including the divisor x0.

The companion family is
\[
 X=(p1,\ldots,pr,-[T^{r+1}]
                       \log(1+p1T+\cdots+prT^r)).
\]
It has full rank and det DS=gamma^r, but its generic extension
is cyclic of degree r+1. Taking logarithmic coefficients of
(1+gamma T)(1+sum p_i T^i) expresses gamma^(r+1) in the target.
There is no unramified tail, so it does not give a Keller map. More explicitly, if L_Y is that logarithmic target coefficient, then d=Y_last+L_Y is linear in Y_last and d(S)=(-1)^r*gamma^(r+1)/(r+1). The first r target coordinates recursively recover p from gamma, proving the cyclic degree and the constant-tail assertion in every r.
Strictly upper triangular nilpotency alone is insufficient.

The concrete remaining PSL2(11) task is now to construct the
rank-four polynomial field with the group and missing-boundary
configuration in (9), or to rule it out. Its affine completion
would then be supplied by the actual finite algorithm above.
This does not decide rank-deficient higher-dimensional inputs,
other architectures, arbitrary Keller actions, or the full spectrum.

Formal verification, independent whole-proof replay and external
review are not claimed. Algorithm termination, polynomial inverses,
polynomiality at x0 and the whole Jacobian have separate written
arguments and exact controls.

[Constructor](higher_completion.py) · [model](higher_completion_model.json) ·
[receipt](higher_completion_receipt.json) ·
[retained inputs](sources/higher-completion-manifest.json).
