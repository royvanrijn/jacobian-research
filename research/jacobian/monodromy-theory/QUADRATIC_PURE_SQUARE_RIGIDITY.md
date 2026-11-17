# Complete rigidity of the matched pure-square quadratic architecture

The higher-degree multiplier frontier closes without a degree search.
Every multiplier divisor is preserved by a straight-line homotopy.
In three projective variables, a nonlinear divisor is ruled by a
single selected family of lines. A polynomial automorphism must act
affinely on each selected line. The leading quadratic base locus
then contradicts full-rank cyclicity in both remaining normal forms.

This closes the entire homogeneous rank4 quadratic architecture with
matched V0=u^2, including arbitrary linear mixing, constants and scalar
degree. It is not an absolute Keller action obstruction: quadratic
directions without that matched component, higher primitive directions
and nonhomogeneous architectures remain open.

## 1. All multiplier divisors are invariant under the line homotopy

Let r>=2, p=(u,x), d>=1 and
\[
                       V=(u^d,u^d v(x/u)),
\]
where v is any polynomial map of degree at most d.
Let nonzero homogeneous F give a normalized input H=FV.
Put W=v-x, B=Dv, J=DW=B-I and Fbar=F(1,x).

**KMFI139.** For every irreducible factor q of Fbar, of any degree,
\[
                         q(x+tW(x))\in(q)\,k[x,t].       \tag{1}
\]
Thus the affine hypersurface q=0 contains the entire straight line
through a generic point x with direction W(x).
If v is a polynomial automorphism, also q(v)=c*q for a nonzero constant c.
No full-rank assumption or automorphy is needed for the line-containment lemma.

The homogeneous determinant identity
detD(FV)=(s+d)/d*F^r*detDV=(-f)^r makes f=-Fh for a polynomial h.
This includes detDV0, when f0; otherwise it follows by UFD divisibility.
The projective resolvent identity has numerator (T-h(1,x))^r and
a denominator monic in T, hence no source-divisor pole.
Its residues, exactly as in KMFI136 but without requiring h constant
on u1, give
\[
                           dq\,B^jW=0\pmod q\quad(j\ge0).
\]
Only j=0,...,r-2 are needed; Cayley--Hamilton gives all higher powers.
Since J=B-I, this is equivalent to dq*J^jW divisible by q for all j.
The following formal argument turns these finitely generated
tangencies into exact polynomial line containment.

Work in R[[t]], R=k[x], and put Phi_t(x)=x+tW(x).
It is an automorphism over the formal power-series ring because
it is the identity modulo t. Its formal inverse G_t has polynomial
coefficients in x; no global polynomial inverse of Phi_t is asserted.
Let
\[
               Y_t=(I+tJ)^{-1}W=\sum_{j\ge0}(-t)^jJ^jW .
\]
Each coefficient is tangent to q=0, so dq*Y_t=q*a_t for a_t in R[[t]].
Differentiating Phi_t(G_t(y))=y gives
\[
                       \partial_t G_t=-Y_t(G_t).
\]
Therefore U_t=q(G_t) satisfies
\[
                    \partial_t U_t=-a_t(G_t)U_t,\quad U_0=q.
\]
In characteristic zero, the unique solution is q times a formal unit
(the exponential of an integrated formal coefficient).
Thus G_t preserves the principal ideal(q), and so does its formal
inverse Phi_t. Since q(Phi_t) is polynomial in t, coefficientwise
divisibility in R[[t]] is precisely divisibility in R[t], proving(1).
If v is automorphic, at t1 irreducibility makes the quotient
a nonzero constant.

The full residue sequence matters. For q=z-xy^2 and
v=(2x,-y,-2z), the first two residues vanish modulo q but the third
does not. The homotopy does not preserve q. This does not contradict
the quadratic-factor theorem; a cubic Taylor remainder is present.

## 2. A degree bound from the leading base locus

Use three affine variables x, so the homogenized input has rank4.
Suppose v is a polynomial automorphism of degree D>=2, and its highest
homogeneous part Q_D has base locus contained in finitely many
linear hyperplanes:
\[
                  \{Q_D=0\}\subseteq\bigcup_{i=1}^m\{\ell_i=0\},
                  \qquad\ell_i\ne0.                     \tag{2}
\]

**KMFI140.** If the corresponding homogeneous input is normalized
and full rank, every irreducible factor q of Fbar has degree at most D.
This is unbounded in total multiplier degree and in factor multiplicities.
It applies to any direction homogenization degree d>=D.

Suppose degq>D. W has degree D, so it is nonzero at a generic point of
the irreducible surface S_aff=V(q): otherwise q would divide every
component of W, forcing W0 globally and contradicting D>=2.
By KMFI139, the assignment
\[
                         x\longmapsto L_x=x+kW(x)
\]
gives a rational family of lines covering this surface.

Here is the elementary line-family fact needed in the proof.
For an irreducible nonplane surface S in P3, an irreducible algebraic
family of lines covering S is a curve, and through a general point
there is exactly one line from that selected family.

Indeed a covering family of dimension>=2 would give infinitely many
lines through a general smooth point p. All lie in the tangent plane
T_pS, and their union is dense in that plane; this would make S a plane.
Thus a covering family has dimension1.
If through a general point there were two lines from the same family,
then infinitely many members would meet a general member L.
The closed incidence condition and irreducibility make every pair
of family members intersect. A pairwise-intersecting family in P3 is
either coplanar or concurrent: take two members meeting at p;
if a third misses p, every line meeting all three lies in their plane.
The coplanar case makes S a plane. In the concurrent case a general
point off the vertex lies on just the unique joining line, contrary
to the assumed multiple incidence. This proves the fact.
It concerns the selected irreducible family, so it does not assume
that a surface has only one possible ruling component.

Apply it to the projective closure of S_aff, using the closure of
the rationally assigned lines L_x. On a general selected affine line L,
a generic y in L is assigned the same line L_y=L.
Thus v(y)=y+W(y) lies in L, and v(L) is contained in L.
Because v is an ambient polynomial automorphism, its restriction is
an isomorphism onto a closed one-dimensional subset of the affine line.
That subset is the whole line. Hence v|L is an affine automorphism.

The coefficient of t^D in v(x+tW(x)) must therefore vanish:
\[
                              Q_D(W(x))=0
\]
at a generic point of q=0.
Condition(2) and irreducibility force one ell_i(W) to vanish there.
Its degree is at most D<degq, so ell_i(W) vanishes identically.
Differentiation gives ell_i*B=ell_i, hence ell_i*B^jW=0 for every j.

But full rank makes W cyclic for B by the exact determinant identity
detC=detB*det(W,BW,B^2W), valid without triangularity.
The nonzero constant row ell_i cannot kill that span.
This contradiction proves the degree bound.

The hypotheses are essential. Without cyclicity, the quadratic
automorphism v=(x+1,y,z+y^2) preserves the irreducible cubic z-xy^2
and acts affinely on its ruling lines. Its W_y is identically zero.

## 3. The complete matched pure-square quadratic rank4 theorem

Let V be homogeneous quadratic in(u,x,y,z), V0=u^2.
No primitive or triangular hypothesis is imposed.

**KMFI141.** If a nonzero homogeneous polynomial F gives normalized
H=g-fp=FV with generic full rankDX4, then the projective
dehomogenization v is affine. Consequently
\[
                              V=uMp,                     \tag{3}
\]
with M constant, and every multiplier and generic field is exactly
in the factorized KMFI121--122 classification.
Conversely the KMFI121 data with first matrix row(1,0,0,0)
and scalar G divisible by u give all such surviving representations
via F=G/u.

In particular no primitive quadratic direction with this matched
pure-square first component admits a normalized full-rank multiplier
of any degree. Every degree11 field in the whole architecture is
S11,A11 or regularC11; neither exceptional PSL2(11) action occurs here.

The pure-square determinant gate and KMFI133 make v polynomial
automorphic. If its quadratic image has dimension<=1, KMFI132 already
excludes a nonaffine full-rank input.
Otherwise, with arbitrary invertible M and constants b retained,
\[
 v=M(x+H_A)+b\quad\hbox{or}\quad v=M(x+H_B)+b,
\]
\[
                       H_A=(0,x^2,y^2),\qquad H_B=(0,x^2,xy).
\]
Their leading quadratic part Q=M H has base locus contained in x=0,
because H's second coordinate is x^2 and M is invertible.

KMFI140 with D2 bounds every irreducible factor of Fbar by degree2.
KMFI137--138 exclude each such quadratic or non-u linear factor.
Thus Fbar is constant, and homogeneity makes F a pure power of u.
KMFI135 then forces v affine, a contradiction to the remaining
nonaffine case. This closes every M,b and every scalar degree.

The conclusion has an explicit geometric boundary: the component
V0=u^2 must be matched to that same source coordinate u under allowed
constant homogeneous conjugacy. An unrelated square linear combination
of output components does not automatically supply this hypothesis.
The theorem does not classify all homogeneous quadratic directions.

## 4. Evidence and the open action problem

[Helper](divisor_homotopy.py) ·
[checker](verify_divisor_homotopy.py) ·
[model](divisor_homotopy_model.json) ·
[receipt](divisor_homotopy_receipt.json).

Exact controls verify the formal inverse and velocity equation,
all-residue versus truncated-residue boundary, cubic line containment,
whole restriction-to-line identities, both leading base loci,
and an affine full-rank survivor and normalized nonlinear rank3 control.
The formal inverse is used only in a completed ring to prove a polynomial
divisibility identity; no formal chart is being substituted for affine space.

The full faithful-transitive Keller monodromy spectrum remains open.
The matched pure-square primitive quadratic route is now closed.
Other quadratic primitive directions, higher primitive degree,
nonhomogeneous inputs, actual exceptional11 or an absolute obstruction,
and other undecided actions remain outstanding.

Formal verification, independent whole-proof replay and external review
are not claimed. The goal remains active.
