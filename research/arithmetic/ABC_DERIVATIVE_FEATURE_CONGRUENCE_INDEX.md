# Integral comparison of derivative congruences and feature differences

The quotient lattice below gives an integral comparison when two
complete letters are affine functions of the same integer features.
It accounts for coincident prime constraints by their actual joint
index. The computation is uniform in the feature dimension and in the
number of rows.

## DI1. The exact index and its affine coordinates

Take \(N\) positive primitive rows \(A_i+B_i=C_i\) and an integer
feature model with both \(A,B\) in its rational evaluation space
\(F\subseteq\mathbb Q^N\). Include the constant vector \(1\) in
\(F\). Let \(U\) be an integer basis matrix for \(F\cap\mathbb Z^N\)
with \(r=\dim F\) columns. This lattice is saturated. The constant
vector is primitive and can be the first basis column. Denote the
row of \(U\) at \(i\) by \(U_i\). The original integer-valued
letters have integral coefficient vectors \(\ell_A,\ell_B\),
\(\ell_C=\ell_A+\ell_B\) in this basis.

Choose a set of primes for arithmetic derivatives. As in AD1, an
integer prime-coordinate vector generates \(D(T)=T\sum_p v_p(T)z_p/p\).
For now suppose every selected prime occurs in at most one row-letter
and each row-letter is represented. Put
\[
 g_{i,T}=\gcd_{p\ {\rm selected},\,p\mid T_i}
                                 \bigl((T_i/p)v_p(T_i)\bigr).
\]
The actual derivative evaluation image in the row-letter slots is
exactly the product of the lattices \(g_{i,T}\mathbb Z\), by
independent Bezout identities on the disjoint supports. Requiring
\(D(A),D(B)\in F\) and additivity places the derivative values in
the image of the colored evaluation matrix \(\mathcal U\). Its rows
are \((U_i,0),(0,U_i),(U_i,U_i)\) in colors \(A,B,C\).
Consequently the parameter lattice is
\[
 \Gamma=\{\theta\in\mathbb Z^{2r}:
                     g_{i,T}\mid\mathcal U_{i,T}\theta
                     \text{ for every slot}\}.
\]
All implications here are integral, including existence of a
preimage under the private derivative evaluation.

For any prime \(p\), set \(m_j=v_p(g_j)\) and omit slots with
\(m_j=0\) in the following expression. Let \(W\) be their colored
row matrix and \(D_p=\operatorname{diag}(p^{m_j})\). The \(p\)-adic
index exponent is exactly
\[
 \log_p[\mathbb Z_p^{2r}:\Gamma_p]
    =\sum_j m_j-v_p\bigl(\Delta[\,W\mid D_p\,]\bigr),
\]
where \(\Delta\) is the gcd of maximal minors.
Indeed the image of the constraint map is
\((W\mathbb Z_p^{2r}+D_p\mathbb Z_p^{\# j})/
D_p\mathbb Z_p^{\# j}\); take its index by Smith normal form.
This index formula also applies to any externally prescribed
depths. When considering shared derivative prime coordinates it
computes the index of these constraints if they are given, and
requires an independent description of the actual evaluation image.

At each depth \(h\), take the spans \(V_T\) of the active
evaluation rows of color \(T\). The colored rank, over either
\(\mathbb Q_p\) or \(\mathbb F_p\), is
\[
          \dim V_A+\dim V_B+\dim V_C
                       -\dim(V_A\cap V_B\cap V_C).
\]
This follows by taking the kernel of
\((u,v,w)\mapsto(u+w,v+w)\).
In particular, for a prime which is present at a single slot
of the congruence system, the constant basis column makes the
index contribution exactly \(p^{m_j}\), with no dependence on
the feature differences. The derivative-minor and feature-minor
problems already separate on this genuine independent evaluation
case.

If each row-letter has one selected private prime, the derivative
coordinates themselves give a further exact formula. Write
\(M=\operatorname{diag}(g_j)^{-1}\mathcal U\). The restricted
derivative lattice is \(M\Gamma\), of Euclidean covolume
\[
 [\mathbb Z^{2r}:\Gamma]\left(
       \sum_{|J|=2r}\frac{\det(\mathcal U_J)^2}
                          {\prod_{j\in J}g_j^2}\right)^{1/2}.
\]
This is the determinant of the parameter sublattice times its
linear volume factor, with Cauchy–Binet for the latter. Thus
colored feature minors occur in the real factor, and their
\(p\)-adic overlaps occur in the preceding integral index.
The formula does not multiply each depth by an independent
dimension loss.

## DI2. Removing the complete-letter directions at a shared prime

Consider any one-prime constraint system as above with working
depths \(m_{i,T}\le v_p(T_i)\). Suppose \(p\) is active on two
distinct colors. The coefficient vectors \(\ell_A,\ell_B\) are then
independent modulo \(p\). In fact, if they are dependent, all
nonzero coefficient vectors among \(\ell_A,\ell_B,\ell_C\) are
proportional. Vanishing of two actual letters on a row contradicts
primitivity. Thus whenever a letter vanishes its whole coefficient
vector must be zero modulo \(p\), and only one color can occur.
Independence allows integral \(p\)-adic coordinates with
\[
                        U_i=(A_i,B_i,w_i).
\]
The transformation is unimodular since the two primitive independent
coefficient columns extend to a basis over \(\mathbb Z_p\).

For an active \(A\)-row divide the constraint by the unit \(B_i\);
for a \(B\)- or \(C\)-row use the unit \(A_i\). The first two
entries of the resulting \(U_i\), modulo its working depth, are
\((0,1),(1,0),(1,-1)\) in colors \(A,B,C\). Denote the normalized
tail by \(\bar w_i\), and choose, within each active color, one
anchor \(i_T\) of maximal depth \(M_T\). Subtract the anchor
constraint from all the constraints of its own color. This
preserves the congruence system since the anchor depth is at least
the depth being subtracted. The residual constraint at \(i,T\)
is one of
\[
 (\bar w_i-\bar w_{i_T})\cdot\alpha',\quad
 (\bar w_i-\bar w_{i_T})\cdot\beta',\quad
 (\bar w_i-\bar w_{i_T})\cdot(\alpha'+\beta')
                 \equiv 0\pmod{p^{m_{i,T}}},
\]
where the tail parameters \(\alpha',\beta'\) have total dimension
\(2(r-2)\). Each original anchor likewise has its tail but its
head, on four parameter coordinates, is one of
\[
                  (0,1,0,0),\quad (0,0,1,0),\quad(1,-1,1,-1).
\]
These rows have a unit maximal minor for any selected subset.
A triangular integral change of head coordinates absorbs the anchor
tails into those head coordinates. Therefore the exact index is
\(p^{\sum_T M_T}\) times the index of the residual tail constraints.
The same calculation applies when just one color is active if the
coefficient vectors are independent modulo \(p\).

In particular on a rank-one affine feature model for **both**
letters (\(r=2\)) in the independent-head case above, the index is exactly
\(p^{\sum_T M_T}\). Within each color only the deepest constraint
contributes. Across distinct primes the index factors, but the
head or tail directions over different primes can be the same
rational direction. The expected dimension losses cannot be
summed from the index alone.

For an elementary check, write the normalized rank-one model as
\(A_i=a_0+a_1 z_i,\ B_i=b_0+b_1 z_i\), with
\(1,(z_i)\) an integral evaluation basis. Then
\[
 A_i B_j-A_j B_i=(a_0b_1-a_1b_0)(z_j-z_i).
\]
The product of the three same-role gcds of two original rows
divides the determinant by additivity and primitivity. The same
product outside the support of \(a_0b_1-a_1b_0\) therefore divides
the true feature difference. A prime switching roles between
the rows gives a unit projective determinant, so the switch
cannot create this divisibility.


## DI3. Shared derivative coordinates and the actual feature lattice

Allow selected primes to occur in arbitrary slots. Let \(E\) be the
integer derivative evaluation matrix with slots \((i,T)\) as rows
and selected primes \(q\) as columns. Its entry is
\((T_i/q)v_q(T_i)\) if \(q\mid T_i\) and zero otherwise. Keep
the saturated colored feature matrix \(\mathcal U\) of DI1. Define
\[
  L_E=E\mathbb Z^{\#q},\quad
  K=\{\theta\in\mathbb Z^{2r}:\mathcal U\theta\in(L_E)_{\mathbb Q}\},
  \qquad
  \Gamma=\{\theta\in\mathbb Z^{2r}:\mathcal U\theta\in L_E\}.
\]
Here \(K\) is saturated and \(\Gamma\) has finite index in \(K\).
If \(B\) consists of the columns \(\mathcal U\theta_j\) for a
basis of \(K\), then
\[
 [K:\Gamma]=\frac{\Delta(E)}{\Delta([\,E\mid B\,])},
\]
where \(\Delta\) denotes the positive gcd of minors at the rank
of the matrix in question, with \(\Delta=1\) at rank zero.
Indeed \(\mathcal U\) identifies the quotient \(K/\Gamma\) with
\((L_E+\mathcal U K)/L_E\). Both lattices on the right have the
same rational span and saturation, so Smith normal form gives
the ratio of their saturation indices. Equivalently, in a Smith
presentation of \(E\), the rows beyond its rank impose exact
linear equations on \(\mathcal U\theta\) defining \(K\), and the
other rows impose divisibility by the nonzero Smith factors.
The rational equations and integral divisibility are separate
steps even when all derivative equality constraints have a
nonzero kernel.

There is a one-digit comparison to the actual joint valuation
depths. For a prime \(p\) write
\[
 W_p=\bigoplus_{i,T}p^{v_p(T_i)}\mathbb Z_p .
\]
Every non-\(p\) column of \(E\) belongs to \(W_p\). The \(p\)
column \(E_p\), taken as zero when unselected, satisfies
\(pE_p\in W_p\). On the available parameter space \(K_p\) put
\[
 \Gamma_p^0=K_p\cap\mathcal U^{-1}(W_p),\qquad
 \Gamma_p^+=K_p\cap\mathcal U^{-1}(W_p+\mathbb Z_p E_p).
\]
The actual \(\Gamma_p\) is contained in \(\Gamma_p^+\), and the
quotient \(\Gamma_p^+/\Gamma_p^0\) embeds in
\((W_p+\mathbb Z_p E_p)/W_p\), of order at most \(p\).
Consequently
\[
 \log_p[K_p:\Gamma_p]\ \ge\
 \log_p[K_p:\Gamma_p^0]-1 .
\]
The loss is zero if \(E_p\in W_p\), or if its nonzero class is
outside the image of the available parameter space modulo \(W_p\).
If the non-\(p\) columns span \(W_p\), then
\(\Gamma_p=\Gamma_p^+\) and the comparison is exact. A sufficient
condition for this span is a selected private \(q\ne p\) in every
slot whose valuation there is a \(p\)-adic unit.

When \(K_p=\mathbb Z_p^{2r}\) and at least two colors are active,
DI2 computes the coordinatewise index exactly, with the full
depths \(m_{i,T}=v_p(T_i)\). In a common rank-one affine model
(\(r=2\)), the sole possible lost digit can be tested directly.
Normalize each active slot as in DI2 by its chosen unit \(w_i\).
Its \(E_p\) image is
\[
 p^{m_{i,T}-1}m_{i,T}(T_i/p^{m_{i,T}})/w_i
                      \pmod {p^{m_{i,T}}}.
\]
For its class to be hit, within each color these residues must
be compatible with the maximal-depth anchor. The independent
head forms in DI2 show that this compatibility is also sufficient.
At all lower positive depths it requires \(p\mid m_{i,T}\).
At equal maximal depth the leading digits must coincide modulo
\(p\). Thus for \(p>\log_2 X\), variation among the positive
depths of any one color prevents loss of that digit.

For comparison, Laniewski–Müller compute the Wronskian image
for one primitive row by a pair of integer linear forms. Their
pointwise index does not impose a common feature evaluation
space on derivatives of multiple rows. The shared-coordinate
restriction above retains the simultaneous exact linear equations
before computing the integral index.
[Laniewski–Müller, Lemma 4.1 and its proof](https://arxiv.org/html/2610.26490v1#S4)

### Tests on genuine primitive rows

The two rows \((6,23,29),(7,34,41)\) in the \(X=41\) shell show
the sharp one-digit loss. Select all their primes and take
\(F=\mathbb Q^2\) on the two row positions. Every slot has a
private odd prime to order one. At \(p=2\) the depth-one slots
are the first \(A\) and second \(B\), whose derivative column
entries are \(3\) and \(17\). Thus the coordinatewise system
requires two independent evenness constraints and has index
\(4\). The actual image also permits simultaneous odd residues
in the two slots and has index \(2\).

The effect of the rational restriction is visible on
\((1,8,9),(1,9,10),(1,11,12)\), all in the \(X=12\) shell.
Select \(2,3,5,11\). Their three additive derivative equations
have primitive kernel generated by \((2,4,7,40)\), in the indicated
prime order. Its three \(A\)-derivatives are zero, and its
\(B\)-derivatives are \((24,24,40)\). With the full evaluation
space \(F=\mathbb Q^3\) this gives rank one for \(K\). With the
actual rank-one feature \(B=(8,9,11)\), use
\(F=\operatorname{span}_{\mathbb Q}(1,B)\). The same derivative
vector is not affine in \(B\), so \(K=0\). The additive derivative
matrix, including all its AD2 data, is identical for these two
models. These finite tests use genuine addition and primitivity
and do not impose the limiting small-budget caps.


## DI4. All common chart moduli act through two integral directions

Let positive primitive rows lie on one chart with fixed full divisors
\(D_A,D_B,D_C\) and fixed signed register projections modulo
\(m_T\mid D_T\). The pairwise coprimality of the \(D_T\) follows from
primitivity. Suppose \(A,B\) are affine integer functions of common
integer features in a lattice \(V\cong\mathbb Z^r\). Subtract an actual
anchor row, and write \(\ell:V\to\mathbb Z^2\) for the integral map
from feature differences to \((\Delta A,\Delta B)\). Write
\(\Lambda=\ell(V)\), \(D=\prod_T D_T\) and \(m=\prod_T m_T\).

The simultaneous position congruences define a lattice \(L_0\subseteq
\mathbb Z^2\). Its quotient has Smith factors \(m,D\), with a factor
of one allowed. At a prime of \(D_T\), take \(T\) and the signed
register operand as position coordinates. The change has determinant
\(\pm1\) and the two local factors are \(p^{v_p(m_T)}\) and
\(p^{v_p(D_T)}\). At any other prime both local factors are one. The
global statement follows by CRT.

The simultaneous congruences in the feature lattice define
\(\Gamma=\ell^{-1}(L_0)\), with the exact index
\[
 [V:\Gamma]=[\Lambda:\Lambda\cap L_0]
          =\frac{mD}{[\mathbb Z^2:\Lambda+L_0]} .
\]
The quotient \(V/\Gamma\) is isomorphic to the subgroup
\((\Lambda+L_0)/L_0\) of the position quotient. It has at most two
nontrivial Smith factors, all congruences vanish on the same
saturated sublattice \(\ker\ell\) of rank at least \(r-2\), and these
statements hold at primes where \(\ell\) is degenerate modulo \(p\).
When the affine differences span the feature space, their maximal
minors in a basis of \(\Gamma\) are their original minors divided by
\([V:\Gamma]\). This extracts the joint determinant factor at once.
Elimination of more than two directions needs information beyond the
common position congruences. JD8 takes \(\Lambda=\mathbb Z^2\) in
position coordinates.

There is a corresponding restriction on the actual derivative
columns in AD2 for the primes already fixed by the chart. If a
selected prime \(p\) divides a full \(D_T\), its valuation on every
row is the same \(e_p=v_p(D_T)\) in role \(T\) and is zero in the
other two roles. In the integer derivative evaluation matrix its
column is supported on the \(T\) slots and there it is
\((e_p/p)(T_i)_i\). All fixed chart primes span at most three rational
directions in derivative evaluation and at most two after taking
the signed additivity residual. For the latter, use
\((C_i)_i=(A_i)_i+(B_i)_i\). If at least two rows are distinct,
\((A_i)_i,(B_i)_i\) are independent by primitivity. With chart primes
selected in every role, an additive derivative formed only from
those primes therefore has equal logarithmic rates on the three
roles and its evaluation image has rational rank one. These column
equalities hold for genuine primitive rows at all radical caps.
Large local divisibility from common support can coexist with this
bounded rational rank. DI3 is still needed for the derivative
columns on varying support.

### The last derivative digit on a full register cell

There is an integral comparison for the remaining derivative
digit. Fix a selected prime \(p\) with \(p^e\parallel D_C\), and
work on the basis rows of AD2 on a full joint register cell.
Write \(r_0\) for their number and let the \(i\)-th row of \(W\)
be \(V_{A_i}-V_{B_i}\). Its rational rank is \(s\). The relation
module \(\mathcal K=\ker(W^{\mathsf T}:\mathbb Z_p^{r_0}\to
\mathbb Z_p^S)\) is saturated. Let \(M\) be the additive
derivative matrix of these rows, with each row divided by the
unit \(A_i\). Since \(B_i=C_i-A_i\), for any \(\lambda\in\mathcal K\)
every non-\(p\) column of \(\lambda^{\mathsf T}M\) is divisible by
\(p^e\). For the \(p\) column there is the exact identity
\[
 p^{1-e}(\lambda^{\mathsf T}M)_p
     =-e\sum_i\frac{\lambda_i(C_i/p^e)}{A_i}.
\]
The fixed registers give \(A_i=A_0\bmod p^e\) and
\(C_i=C_0+D z_i\) with an actual anchor. Consequently, for
\(\sum_i\lambda_i=0\) the expression on the right is
\[
            -\frac{e}{A_0}\frac{D}{p^e}\sum_i\lambda_i z_i
                                                     \pmod{p^e}.
\]
This comparison uses unit normalization and a saturated relation
module, with no denominator cost at \(p\). If \(p\nmid e\),
the coefficient is a unit, so the divisibility on the affine
relations is exactly the divisibility of the corresponding
linear form in the actual quotient features, up to depth \(e\).

In particular, if all relations in \(\mathcal K\) annihilate
\((C_i/p^e)_i\bmod p\), or if \(p\mid e\), the contribution
\(e(r_0-s)-1\) from this \(p^e\) in AD3 improves to
\(e(r_0-s)\). Extend a basis of \(\mathcal K\) to a
\(\mathbb Z_p\)-basis. The corresponding \(r_0-s\) transformed
rows now have every column divisible by \(p^e\), and the row
change is unimodular. For \(p\nmid e\) the annihilation
condition uses all of \(\mathcal K\); a check only on affine
relations suffices if \(1\) belongs to the \(\mathbb Q_p\)-column
span of \(W\). The result tests one digit of the common-depth
bound. A uniform rank saving would still need correlations on
the varying supports or a numerical bound on the quotient
features across different charts.

[Full position index (JD8)](ABC_JOINT_DISPLACEMENT_CATALOGUES.md#a-rank-two-affine-index-bound)

## DI5. Prime-by-prime restoration forces an exact feature relation

Retain the full joint register cell of DI4 and its basis rows.
Let \(s'\) be the rank over \(\mathbb Q\) of the matrix
\([1\mid W]\), where \(W\) is the valuation-difference matrix,
and write \(Z=\max_i|z_i|\). Let \(\mathcal P\) be a set of
selected primes dividing \(D_C\) for which \(p\nmid e_p\),
with \(e_p=v_p(D_C)\). Suppose for each \(p\in\mathcal P\)
and each integral affine relation \(\lambda\) among the rows of
\([1\mid W]\) that
\[
                   p^{e_p}\mid(\lambda^{\mathsf T}M)_p
                   \quad\text{in }\mathbb Z_p .
\]
The normalized derivative matrix \(M\) is the one in DI4.
If
\[
       \prod_{p\in\mathcal P}p>
       (s'+1)(\sqrt{s'}\max(1,\log_2 X))^{s'} Z ,
\]
then the actual quotient vector \(z=(z_i)_i\) belongs to the
rational column span of \([1\mid W]\).

To prove this, choose \(s'\) independent rows and \(s'\)
columns giving a nonzero square minor. Adjoin each remaining row
in turn and use the signed \(s'\)-minors on those columns to
produce an integral affine relation, supported on at most
\(s'+1\) rows. The coefficient of the newly adjoined row is
nonzero and these relations span the whole relation space.
Their coefficients have absolute value at most
\((\sqrt{s'}\max(1,\log_2 X))^{s'}\) by Hadamard. DI4 makes
each \(\lambda^{\mathsf T}z\) a multiple of every prime in
\(\mathcal P\); the displayed numerical bound forces it to vanish.
In fact the hypothesis only needs the specified relations.
The factor in the bound depends on the valuation rank, with no
dependence on the number of rows. For example, if \(Z\le n^K\)
and the prime product is at least \(n^{K+\sigma}\) for fixed
\(K,\sigma>0\), this forces exact dependence whenever
\(s'\log(s'\log X)=o(L\log X)\) at fixed \(L>0\).

The restoration input refers to the selected relation values.
The depth-defect bound on the gcd of maximal minors alone gives
no such input. All the primes here recover a single complete
letter's quotient direction on the same cell. A saving in a
growing feature dimension needs more than duplicating this
direction at many primes.

## Obstruction and the remaining extraction

Here is a genuine rank-one obstruction to deducing the numerical
size of a feature minor from derivative-minor gcd data alone.
Use the rows
\((t+h,\,2(t+h)+1,\,3(t+h)+1)\) at the shifts \(h=0,1,H\).
Let \(P\) be a product of nine fixed distinct primes greater than
100, take \(H\equiv2\bmod P\), and \(4H\le t\le5H\).
At each selected prime, the nine affine forms have distinct roots
modulo that prime, by their shifts \(0,1,2\) modulo \(p\).
Assign one selected prime privately to each form, and choose \(t\)
by CRT modulo \(P^2\) so the assigned valuation is exactly one.
For large \(H\) such \(t\) exist. The three primitive rows lie in
one shell, for example \(X=20H\). The selected derivative matrix
has disjoint row support; each row gcd is one from the coprime
coefficients on its first two letters. Thus its maximal-minor
gcd is one and all AD2 depth deficits are zero. The saturated
feature lattice has basis \(1,(h_i)\), since a difference is one,
and its largest difference is \(H\). Replacing \(H\) by the shift
2, and retaining the same large parameter range for \(t\), gives
the same derivative-minor gcd and deficits and bounded feature
minors. This counterexample imposes primitivity and the numerical
shell. It makes no claim of simultaneous small radical caps.

### The same obstruction on a deep register cell

The comparison can have a nonzero rank defect and still leave the
feature minors undetermined. Fix primes \(q>3\), \(s,t>20q\),
and \(p>\max(s,t)\) with \(p\equiv1\bmod 3\). Take arbitrarily large
\(D=p^m\equiv1\bmod 3(qst)^2\) with \(p\nmid m\). Choose
\(H\equiv2\bmod (qst)^2p\) with \(H\asymp D^\gamma\) for any fixed
\(\gamma>0\). Compare the shift lists \(z=(0,1,H)\) and \(z=(0,1,2)\)
in the primitive rows
\[
 A_i=(D-1)/3+D(T+qz_i),\quad B_i=2A_i+1,\quad C_i=3A_i+1.
\]
Use the same \(T\in[4qH,5qH]\) for both lists. Require
\(T\equiv q\bmod q^2\); all three \(A_i\) then have \(q\)-valuation
one. At \(s,t\), the nine form roots modulo the prime are distinct,
since the shifts \(0,q,2q\) and the translates by \(0,1/2,1/3\)
are distinct. Choose \(T\) modulo \(s^2,t^2\) so that \(s\) occurs
only on \(B_3\) and \(t\) only on \(C_3\), both to order one.
Avoid the at most three values of \(T\bmod p\) that would raise a
target valuation above \(m\). These constraints work for both
lists and have fixed combined modulus. Thus CRT gives the chosen
\(T\) for large \(D\). Both lists lie in \(X/2<C_i\le X=20DqH\)
and in one full register cell with full divisors \(q,1,D\).

Select \(p,q,s,t\) in the derivative matrix. The first two rows
are supported at \(p,q\); their \(2\)-minor there has absolute value
\(mD/p\), from \(C_i=3A_i+1\) and \(A_2-A_1=Dq\). The third row
has the private coordinates \(s,t\). Thus the three rows are
independent and the nonzero maximal minors are the \(2\)-minor
times \(B_3/s\) or \(C_3/t\). Those two factors are coprime by
primitivity, so the maximal-minor gcd is \(mp^{m-1}\) for both
lists. At every prime and depth, if the third row contributes to
AD2, its difference vector has a private \(s\) or \(t\) coordinate
and increments both the row count and the rank by one. Hence
all AD2 deficits agree between the two lists. In particular the
deficit sum at \(p\) is \(m\), with the one-digit loss attained.

The true saturated rank-one feature has largest difference \(H\)
versus \(2\), since both feature lists have a difference of one.
For any fixed \(L,K>0\), choosing \(0<\gamma<KL/2\) gives depth
\(Dq\ge X/n^K\) with \(n=X^L\). The example has the numerical
height, primitive rows and exact common divisors. It imposes no
small-budget radical caps.

The exact reduction in DI2 concerns a congruence system on an
already common affine model. DI3 computes the actual shared-coordinate restriction and its local
index on that model. An unrestricted application still needs to
extract the model on enough genuine rows retaining all three caps,
or a feature-level divisibility bound without that model. In either
case it needs a bound on the joint rational directions of the
residual systems at primes selectable at low radical cost. The
local index comparison holds only after restriction to \(K_p\). The index alone leaves both
obligations open. A model for only one complete letter does not
give the two complete-letter directions used in DI2.

[Arithmetic derivative lattice (AD1–AD2)](ABC_ARITHMETIC_DERIVATIVE_FAMILY_LATTICES.md)
· [Small affine minors (CF1)](ABC_COMPLETE_LETTER_SMALL_FEATURES.md)
· [Pastén, Theorems 2.5 and 3.3](https://arxiv.org/pdf/2106.16165v3)
