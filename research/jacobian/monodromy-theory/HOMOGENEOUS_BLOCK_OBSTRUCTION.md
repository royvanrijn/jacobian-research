# Actual signed Keller fields with forbidden pair-block fields

Work over algebraically closed characteristic-zero constants, keeping
coefficient fields separate. Write l>=2, m=l-1, n=2l.

**KMFI62.** The homogeneous normalized sweeps
\[
 \sigma=w^{2l-1}-z,\quad b_0=(l-1)/l,\quad d_0=(l-1)/(2l-1),
\]
\[
 X=(\sigma,b_0w^{2l},d_0w^{4l-2}-\sigma^2/2),\qquad
 S=X+\gamma(1,w,z)                                          \tag{1}
\]
all admit explicit polynomial Keller completions on ordinary affine4.
Their degree is4l. The action is binary SL2(F3) on8 for l2, and
V_(2l):A_(2l) on4l signed letters for l>=3.
The coefficient field can be chosen of degree at most4l-6 over Q.
The maps are absolutely and stably atomic with their unique pair block.
The formulas below specify the marked point and every stage operation;
there is no remaining affine-completion equation.

**KMFI63.** The unique proper pair-block field of each completed
map is independently induced from the fixed rational cover
\[
 \tau=\phi_l(u)=F_l(u)^l/u,\qquad
 F_l(u)=\frac12+\frac{u}{l}
       -\frac{(l-1)^2}{2(2l-1)l^2}u^2.                       \tag{2}
\]
Its degree is2l, group A_(2l), and its three branch passports are
\[
 l^2,\qquad 3\,1^{2l-3},\qquad (2l-1)\,1.                    \tag{3}
\]
Their fiber-gcd lcm is l<2l. By KMFI9, this exact block field has
**no polynomial Keller model in any affine dimension**, even after
common independent stabilization and arbitrary birational recharts.

Thus Keller realizability of field pairs is not preserved by taking
a proper intermediate field, even with stable recharting allowed.
Both sides of every block pair are rational function fields by(8),
and its sheet stabilizer A_(2l-1) normally generates A_(2l),
including the C3 stabilizer in A4. Thus normal generation plus
rational source and target is **not sufficient for arbitrary field
pairs**. The normalized polynomial presentation in KMFI59 is an
essential additional hypothesis.
For l>=3 the quotient action A_(2l) is nevertheless realized by
different Keller fields in KMFI15. For l2 this obstructs the newly
available binary map's exact A4 block field; it does not exclude
natural A4 for every other field pair.

**KMFI64.** In the l6 member, the block normal closure has group A12.
Reducing it to either natural PSL2(11) embedded transitively on12,
or exceptional PSL2(11) embedded as11+1, gives fixed target curves of
genera respectively74038 and74226. Both target fields, with their
actual independent spectators, are not stably rational.
Consequently neither reduction produces a Keller map of ordinary
affine space, including the eleven-sheet component in the second
case. This excludes these specific targets, not either abstract
PSL2(11) action.

[Family and marked-stage helper](homogeneous_blocks.py) ·
[exact checker](verify_homogeneous_blocks.py) ·
[model](homogeneous_blocks_model.json) ·
[receipt](homogeneous_blocks_receipt.json).

## 1. Exact homogeneous family and marked affine lift

The normalization of(1) is lambda=w,rho=w^(2l-2), so det DS=gamma^2.
KMFI44 proves it finite of degree4l; KMFI45 gives its exact groups.
Its branch-surface degree is2 and its unramified degree is4l-6>0.
KMFI59 therefore guarantees the completion.

Here are explicit algebraic choices for that theorem, avoiding a
generic-point instruction. At p*=(1,0) the branch target is
\[
                  (Y1,A,Y3)=(1,b_0,d_0-1/2).
\]
With B=Y3+Y1^2/2, the block eliminant is
\[
 Q(q)=(A-b_0q^l)^2+2q^l(A-b_0q^l)+2d_0q^{2l}-2qB.
\]
At A=b0,B=d0 it has the exact triple factor(q-1)^3.
Define the polynomial of degree2l-3
\[
 T_l(q)=
 \frac{(l-1)q^{2l}-2(2l-1)q^l+2l^2q-(l-1)(2l-1)}
      {(q-1)^3}.                                           \tag{4}
\]
It has
\[
 T_l(0)=(l-1)(2l-1),\quad
 T_l(1)=l^2(l-1)(2l-1)/3\ne0.
\]
Choose any root r of T_l(r^2)=0 and put
\[
 \gamma0=b_0(1-r^{2l})/r,\quad
 w0=r,\quad z0=r^{2l-1}-1+\gamma0.                           \tag{5}
\]
These are a regular preimage of the displayed branch target.

To check nonvanishing and separability uniformly, use
\[
 qQ'-Q=-(A-b_0q^l)^2.
\]
At the chosen branch target, a repeated root other than1 would
satisfy q^l=1. But the numerator of(4) at such a root is
2l^2(q-1), so it is nonzero. Thus T_l is separable, r!=0 and
r^(2l)!=1, proving gamma0!=0.
The coefficient degree is at most deg T_l(r^2)=4l-6.

Use old w=1+w,old z=w+z and target rows
\[
 Y1'=Y1-1,\quad Y2'=Y2-Y1+1-b_0,\quad
 Y3'=Y3-Y2+Y1-d_0-1/2+b_0.
\]
The centered normalization is lambda=1+w,
rho=(1+w)^(2l-2)-1. Form the polynomial stage
E_i=S_i'(gamma,gamma*u,gamma^2*v)/gamma^i, i1,2,3, and set
\[
 \theta0=(\gamma0,(r-1)/\gamma0,
                   (r^{2l-1}-r+\gamma0)/\gamma0^2).
\]
Compute L0=DE(theta0)^-1 and
f20=-(D_(L0[:,1])^2 E3)(theta0)/2, where the column label1 means
the first column. Apply exactly KMFI59's single triangular shear,
then v1=xy,v2=x^2zsrc,v3=x^3tsrc and
F=(xgamma,E1/x,E2/x^2,E3/x^3).
These are explicit finite derivative and inverse-matrix operations
over Q(r), not an unresolved existence test.
The polynomial coordinates have full Jacobian gamma0 and the exact
field pair of(1) with one common independent C=xgamma.
The checker replays the entire marked stage for l2 and3 and the
universal algebra and marked equations in further controls.
The all-l conclusion is the written proof and KMFI59.

For l>=3 the arithmetic normal-closure constants over the chosen
coefficient field K0=Q(r) are exactly
\[
 K0\big(\sqrt{(-1)^{l-1}(2l-1)},\sqrt{-(2l-1)}\big).          \tag{6}
\]
The arithmetic group is the preimage of this Galois image under
block sign and flip parity in C2 wr S_(2l).
Indeed the leading coefficient of Q is
-(l-1)^2/[(2l-1)l^2], its root product is
-(2l-1)(lA/(l-1))^2, and its discriminant square class is
(-1)^(l-1)(2l-1), by the displayed qQ'-Q resultant identity.
The geometric group is the full common kernel V_n:A_n, so these
two characters account for every constant.
For l2 use the binary normalizer result KMFI60: Q(r)=Q(sqrt(-3))
already splits its constants.
Atomicity and stable atomicity follow from KMFI30 and61, since
the exact actions have their unique degree-two pair block.

## 2. Two-sided reconstruction of the block curve

Let L=k(gamma,w,z),K=k(S),q=w^2. The central generic deck involution
of(1) sends gamma,w to their negatives. The pair-block intermediate
field is exactly
\[
 E=k(Y1,q,A),\quad A=Y2,\quad B=Y3+Y1^2/2.
\]
The remaining relation is
\[
 B=\frac{A^2}{2q}+\frac{Aq^{l-1}}l
      -\frac{(l-1)^2}{2(2l-1)l^2}q^{2l-1}.                 \tag{7}
\]
Conversely adjoining w with w^2=q recovers gamma=(A-b0q^l)/w
and z=w^(2l-1)-Y1+gamma. Thus [L:E]=2 and [E:K]=2l.

Set the **common target spectator**
\[
                    e=A^2/B,\quad u=q^l/A.
\]
Equation(7) gives qB/A^2=F_l(u), hence q=eF_l(u). Put
tau=B^l/A^(2l-1). The exact source and target reconstructions are
\[
 q=eF_l(u),\ A=e^lF_l(u)^l/u,\ B=A^2/e,
\]
\[
 A=\tau e^l,\qquad B=\tau^2e^{2l-1},\qquad
                    \tau=F_l(u)^l/u.
\]
Therefore
\[
 K=k(Y1,e,\tau)\ \subset\ E=k(Y1,e,u)                         \tag{8}
\]
is precisely the fixed curve cover(2) with independent Y1,e.
The actual affine-four-space map adds the same independent C to both.
This is a field identity, not a conclusion drawn from a particular
pencil or a branch profile alone.

For completeness F_l has discriminant1/(2l-1), nonzero constant
term and distinct nonzero roots, while
\[
 l uF_l'-F_l=-\frac12(1-(l-1)u/l)^2.
\]
The third critical point is u0=l/(l-1), where
F_l(u0)=l^2/[(l-1)(2l-1)]!=0.
Thus(3) lists every branch point, by the derivative and the poles.
The indices sum4l-2=2(2l)-2. Their fiber gcds are l,1,1,
whose lcm is l, strictly below the degree2l.
KMFI9 now forbids every Keller presentation of(8) after common
independent stabilization. The obstruction is intrinsic to this
field pair, not to the particular affine coordinates.

The quadratic layer of the actual Keller field is especially clear:
\[
                         w^2=eF_l(u).
\]
The spectator e is coupled to that layer. Discarding it produces
the forbidden block curve. This does not prove that the full
signed field is curve-induced or non-curve-induced for all l;
the action-level KMFI31 test remains separately scoped.

## 3. PSL2(11) subgroup targets for l=6

Let N/k(tau) be the geometric normal closure of phi_6, with group
A12, and let J be either of the two order660 subgroups:
the natural twelve-letter PSL2(11), or its exceptional eleven-letter
action with one fixed letter. The index is
\[
                         [A12:J]=362880.
\]
The three inertia types in A12 are 6^2,3*1^9,11*1.
For a class g, the number of fixed cosets is
|C_(A12)(g)|*|g^(A12) intersection J|/660.
All powers are required for the order6 inertia.

For the natural embedding the fixed-power arrays are
\[
\begin{array}{c|c|r}
g&\operatorname{Fix}(1,g,\ldots)&\#\text{cycles}\\ \hline
6^2&(362880,6,162,1920,162,6)&60856\\
3\,1^9&(362880,0,0)&120960\\
11\,1&(362880,1,1,1,1,1,1,1,1,1,1)&32990.
\end{array}
\]
For the eleven-plus-one embedding the first array is
(362880,0,0,0,0,0), giving60480 cycles; the other two are the same.
Every fixed-point-free power is absent from that subgroup, and
its order3 type is3^3*1^3, not one3-cycle.

The 11*1 class splits into two alternating classes. Each order11
cyclic subgroup has five powers in each class, so each meets J
in60 elements. The checker verifies this using the signs of all
eleven-cycle multipliers and enumerates all660 matrices for both
actions. Neither A12 nor its362880 cosets needs to be enumerated.

Riemann--Hurwitz gives
\[
 g(N^J)=1+\frac{362880-\sum\#\text{cycles}}2
       =74038\quad\hbox{or}\quad74226.                       \tag{9}
\]
By(8), the actual fixed target fields have exactly these curve
fields with the independent Y1,e,C adjoined.
They are not merely generic high-genus pencils on rational surfaces.
A positive-genus curve field cannot become rational after independent
stabilization: a dominant rational map from projective space to that
curve would restrict to a nonconstant P1 map, contradicting
Riemann--Hurwitz. Thus both targets are intrinsically not stably rational.

For the transitive J the base-changed block cover has degree12.
For the eleven-plus-one J it splits into a degree1 component and
a degree11 component; the latter has the exceptional PSL2(11) action.
Its nonrational target still excludes an affine-space Keller model.
The same target obstruction applies if one chooses another cover
from that J-Galois closure. It does not decide either action from
a different normal closure or a different field construction.

## 4. Scope and assurance

The all-l family, explicit field reconstructions and branch passports
have written proofs. The checker performs exact polynomial algebra,
marked stage controls and the full660-element subgroup censuses.
No formal verification or independent whole-proof review is claimed.

The result rules out the exact natural A4 block field and two
additional PSL2(11) subgroup targets while producing an infinite
family of actual signed Keller fields. Absolute natural A4,
absolute PSL2(11) and the full action spectrum remain open.
