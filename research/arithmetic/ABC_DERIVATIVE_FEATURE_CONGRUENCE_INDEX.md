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
letters (\(r=2\)) the shared-prime index is exactly
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

The exact reduction in DI2 concerns a congruence system on an
already common affine model. An unrestricted application needs
to extract such a model on enough genuine rows retaining all
three caps, or construct an integral map for shared derivative
coordinates directly. It also needs a bound on the joint
rational directions of the residual systems at the primes that
one can select at low radical cost. The index alone leaves both
obligations open. A model for only one complete letter does not
give the two complete-letter directions used in DI2.

[Arithmetic derivative lattice (AD1–AD2)](ABC_ARITHMETIC_DERIVATIVE_FAMILY_LATTICES.md)
· [Small affine minors (CF1)](ABC_COMPLETE_LETTER_SMALL_FEATURES.md)
· [Pastén, Theorems 2.5 and 3.3](https://arxiv.org/pdf/2106.16165v3)
