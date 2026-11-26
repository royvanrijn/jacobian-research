# Uniform prime-first exclusion of the zero cross residue

For \(abck\ne0\), put \(f=ax^2+bxy+cy^2+k\) and
\(M_j=\Lambda_2(f^j)\), with \(\Lambda_2(x^iy^j)=i!j!\).
For every prime \(p\ge7\), the moments \(M_p,\ldots,M_{p+3}\)
cannot all vanish with the normalized residue \(\overline{b/2}=0\) at \(p\).
Together with the
[negative residue lift](BINARY_EVEN_QUADRATIC_PRIME_FRONT_NEGATIVE_RANK_LIFT.md)
and the [diagonal lift](BINARY_EVEN_QUADRATIC_PRIME_FRONT_DIAGONAL_LIFT.md),
this proves the four-window whenever its first index is a
prime at least seven.

The first correction eliminates the singular residue branch
and reduces the other branch to one fixed quadratic.
The second correction cuts its roots down to a fixed finite
set of residue characteristics. A modular factorial
recurrence excludes the one remaining residue after the
small-degree algebraic checks.

We use the normalization \(a=1\) of the
[prime-front reduction](BINARY_EVEN_QUADRATIC_PRIME_FRONT_REDUCTION.md).
The prime \(7\) is covered by the
[zero-residue gcd at that prime](BINARY_EVEN_QUADRATIC_PRIME_FRONT_ZERO_RANK_GATE.md#exact-finite-field-checks-and-scope).
Hence take \(p\ge11\). Write
\[
 h=(p-1)/2,\quad C=(-1)^h,\quad t=b/2,\quad d=c-t^2,\quad
 \kappa=k^p/p^2,\quad \alpha=p/k^h,\quad \beta=C\alpha,
 \quad x=d/k,\quad z=t/\beta.
\]
By the [zero-residue conditions](BINARY_EVEN_QUADRATIC_PRIME_FRONT_ZERO_RANK_GATE.md),
\(v(\beta)=1/p\), \(x,z\) are units, and \(v(\kappa+2)\ge7/p\).
Here \(v(p)=1\). Set \(\theta=\overline{(\kappa+2)/\beta^7}\).
We use \(X\) for the residue of \(x\), and in the following
algebraic formulas for its polynomial variable as well.
Define \(g_i(X)=\sum_{j=0}^i\binom ij(2j)!X^j\),
\(G=g_{h+1}(X)\), \(H=g_h(X)\), and \(\mathcal D=X\,d/dX\).
Working modulo \(p\),
\[
 \mathcal DG=(G-H)/2,\quad
 \mathcal DH=((G-1)/X-H)/2,\quad
 G_2=1+6XG-3XH,\quad G_3=1+20X G_2-15XG,\quad
 G_4=1+42XG_3-35XG_2 .
\]
Here \(G_j=g_{h+j}\). The derivative formulas and the
three consecutive replacements follow from the coefficient
identities in the
[negative residue lift](BINARY_EVEN_QUADRATIC_PRIME_FRONT_NEGATIVE_RANK_LIFT.md#prime-at-the-second-index).
All computations at a candidate have \(X\ne0\).

## The native moment expansion

For \(0\le j\le2\), and a positive integer \(\nu\), let
\[
 \mathcal A_{j,\nu}=
 2\sum_{u=0}^j\binom ju(-2)^u R!z^{R+\nu}
       \sum_{i=0}^u\binom ui X^i(R+\nu+1)_{2i},
                     \qquad R=2(j-u)
\]
inside the summation, with rising factorial notation.
Put \(G_1=G\),
\(\eta_j=-(-4)^{j+1}j!/(2j+1)!!\),
and \(\nu_j=-(-2)^{j+1}/(j+1)\).
The three moment layers after division by \(\beta^{2p+2j+1}\)
are \(E_j+\beta F_j+\beta^2 L_j\), with an extra
\(24\theta\beta^2\) in order \(p+2\), where
\[
\begin{split}
 E_j&=\mathcal A_{j,1}+\eta_jG_{j+1},\\
 F_j&=\mathcal A_{j,2}+\eta_j z(2\mathcal D+1)G_{j+1}
                                      +\nu_jg_{j+1},\\
 L_j&=\mathcal A_{j,3}
     +\tfrac12\eta_jz^2(2\mathcal D+1)(2\mathcal D+2)G_{j+1}
     +\nu_jz(2\mathcal D+1)g_{j+1}
     +\frac{\eta_j}{2j+3}G_{j+2}.
\end{split}                                                \tag{1}
\]
The errors have valuation strictly greater than \(2/p\).
All coefficients can be lifted arbitrarily from their
residue formulas at these precisions.

For the expansion write \(f=(x_1+t y)^2+d y^2+k\).
If \(L\) is the total perturbation count and \(i\) the
number of \(d\)'s, the literal summand in \(M_{p+j}\) is
\[
 \binom{p+j}L\binom Li k^L x^i R!\sum_{m=0}^R(m+1)_{2i}t^m,
                        \qquad R=2(p+j-L).
\]
In the exterior bands pair \(L=u\) with \(L=p+u\), \(0\le u\le j\).
After the coefficient congruences, \((2p+R_0)!/p^2\equiv2 R_0!\)
for \(R_0=2(j-u)\), while the high band has \(R_0!\).
The factor \(k^p=\kappa p^2\) cancels the shared low powers
in the two bands. Their differences at the next three
indices \(m=R_0+1,R_0+2,R_0+3\) give \(\mathcal A_{j,\nu}\).
The lowest uncancelled effect of \(\kappa+2\) is from
\(u=i=m=0\), of degree \(2p+7\) in \(\beta\), with
coefficient \((2j)!\theta\). It enters only for \(j=2\).

In the middle band the groups with \(L\le h+j\) have
one factor \(p\) in the factorial and one in the order
binomial. Only \(L=j+1\) can contribute through the
retained precision, at its \(m=0,1\) terms with residue
coefficient \(\nu_j\). The group \(L=h+j+1\), with
\(R=p-1\), first appears at degree \(2p+2j+1\); the
binomial divided by \(p\) times \((p-1)!\) reduces to
\(j!/[(h+j+1)!h!]\). Using \((h!)^2\equiv-C\) and
\(\beta=C\alpha\) gives \(\eta_j\), and \(m=0,1,2\)
give the three indicated derivative expressions. The
next group \(L=h+j+2\) gives \(\eta_j/(2j+3)\) by the
ratio of adjacent order binomials and the last two
factorial terms. All later middle groups have larger
\(L\) and strictly higher valuation. These bounds
include \(i>h\), for which each rising term has a factor
\(p\). The original smallest possible degree is \(2p\),
so a full factor \(p\) in a coefficient error gives a
strict gain at every displayed precision for \(p\ge11\).

## Boundary data and the first correction

Write \(B_j(a)=\Lambda_1((k+ax^2)^j)\).
For the lower weighted boundary use the polynomial
\(((c-t)B_{p+2}(c)+(1-t)B_{p+2}(1))/\beta^{2p+5}\).
Put \(w'=z^2-2X,\ P=-8(1-2w'+6w'^2)\) and \(A=128/15\).
Its first two layers are
\[
                 E_b=A-zP,\qquad F_b=8/3-zA+w'P .           \tag{2}
\]
The exterior bands of \(B_{p+2}(1)\) cancel at this
precision. The middle indices \(h,p-1\) give
\(A\beta^{2p+5}\) and \((8/3)\beta^{2p+6}\).
At \(c\) only indices \(0,1,2\) contribute to the
lowest layer \(P\beta^{2p+4}\).
At order \(p+3\) the same calculation at precision \(2p+7\)
gives the weighted boundary condition
\[
        720\theta-512/35-z\,16g_3(W)=0,\qquad
                    W=X-2G^2,\quad z=-2G.                  \tag{3}
\]
The first term is the exterior pair at index \(3\) and
\(p+3\), the second is the middle index \(h\), and at \(c\)
the indices \(0,\ldots,3\) give \(16g_3(W)\beta^{2p+6}\).
These use the [weighted boundary equations](BINARY_EVEN_QUADRATIC_SINGLE_FACE_STRONG_FACTORIAL.md#4-the-exact-reduction-for-general-b-and-its-limitation).

Since \(E_0=2z+4G\), eliminate \(z\) by its leading
expression \(z_0=-2G\) and first correction
\(z_1=-F_0(X,z_0)/2\). For \(i=1,2,b\) let
\(e_i=E_i(X,z_0)\) and \(r_i=F_i(X,z_0)+z_1\partial_z E_i(X,z_0)\).
The leading equations \(e_1=e_b=0\) give
\[
 H=\frac{12G^3-3G-6XG+2}{6X},\quad
 G=\frac{8}{15D_0(W)},\quad X=W+2G^2,\quad D_0=1+4W+24W^2 .
\]
In particular \(D_0\ne0\). Then \(e_2=0\) is automatic.
All derivatives below are taken before these substitutions.
Put
\[
 Q_1=120W^2+12W+1,\quad Q_2=360W^2+24W+1,\quad
 I=129600W^5+43200W^4+14400W^3+1800W^2+225W+128.
\]
After substitution \(X=I/(225D_0^2)\), so \(I\ne0\).
The derivative \(J_2=\mathcal D e_2\) is
\(120D_0 Q_1^2/I\). If \(Q_1\ne0\), it is a unit,
forcing the correction in \(x\) to have valuation at
least \(1/p\) by the leading linear part. The two first
corrections must then be compatible. Their determinant is
\[
             J_1 r_2-J_2 r_1=80D_0Q_1Q_2/I,
                               \qquad J_1=\mathcal D e_1 .
\]
Thus \(Q_2=0\) on this branch.

If \(Q_1=0\), the first derivatives of \(e_1,e_2,e_b\)
vanish and the following reductions hold:
\[
          \mathcal D^2 e_b=448(10W+1)/25,\qquad
                 \mathcal D^2 e_2=0,\qquad r_2=8/3 .
\]
The first value is a unit since \(Q_1(-1/10)=1\).
Take a representative of the residue of \(x\) in a
finite unramified extension. The equation for \(e_b+\beta r_b\)
forces the distance from this representative to have
valuation at least \(1/(2p)\). At a smaller valuation,
its nonzero quadratic term is uniquely least.
At the resulting distance all variation of \(e_2\) has
valuation strictly above \(1/p\), using its first and
second derivative vanishings. But \(\beta r_2\) has
valuation \(1/p\). This excludes the singular branch.

## The second correction and fixed exceptional primes

On \(Q_2=0\), \(Q_1\ne0\): their resultant has prime
support in \(2,3,5,7\). The derivative \(J_2\) therefore
remains nonzero. Use
\(z_2=-(L_0(X,z_0)+z_1\partial_z F_0(X,z_0))/2\).
The second coefficients after eliminating \(z\) are
\[
 \begin{split}
 s_j={}&L_j(X,z_0)+z_1\partial_z F_j(X,z_0)
       +z_2\partial_z E_j(X,z_0)
       +\tfrac12z_1^2\partial_z^2 E_j(X,z_0)
       +24\theta\,\boldsymbol 1_{j=2}.
 \end{split}
\]
Here \(\theta\) is determined by (3). If \(x_0\) is an
unramified lift of \(X\), write the change in \(x/x_0\)
through the next two valuation digits as
\(1+\beta U+\beta^2(V+U^2/2)\).
For the first digit use \(U=-r_2/J_2\). The nonzero
linear coefficient in the order-\(p+2\) equation then
forces the remaining correction to valuation at least
\(2/p\), giving the second digit \(V\). This is only a
quadratic Taylor change of coordinates. In this
notation the second-order equations for \(j=1,2\) read
\[
           J_j V+A_j=0,\qquad
 A_j=(\mathcal D^2 e_j)U^2/2+(\mathcal D r_j)U+s_j .
\]
As before derivatives are taken before imposing any
residue equations. The first and second correction
identities are verified by the
[exact symbolic checker](verify_binary_even_quadratic_prime_first_zero_rank.py).

Reducing \(J_2A_1-J_1A_2\) modulo \(Q_2\), at
primes other than \(83,33047233\) in the current range,
gives the necessary condition
\[
        \frac{16(81133029162 W+3369377549)}{10409878395}=0 .
\]
The resultants of \(Q_2\) with \(D_0,I,Q_1\)
have prime supports contained in
\(\{2,3,5,7,83,33047233\}\).
The resultant of \(Q_2\) with the displayed linear
numerator after removing the factor \(16\) is
\[
           4108714878796509472092
              =2^2 3^2\,7\cdot83\cdot33047233\cdot5944184339 .
\]
It remains to check the last three primes. The table
uses the original rational expressions for \(J_j,A_j\)
at the \(Q_2\) roots, with no division by the reduced
modulus resultant.

| \(p\) | \(W\) | Obstruction or remaining condition |
|---|---:|---|
| 83 | 78 | \(D_0=0\) |
| 83 | 16 | \(J_2A_1-J_1A_2=22\) |
| 33047233 | 26667925 | \(X=0\) |
| 33047233 | 21801350 | \(J_2A_1-J_1A_2=31229009\) |
| 5944184339 | 5504500288 | \(J_2A_1-J_1A_2=2163193724\) |
| 5944184339 | 43405095 | \(G=4250342174,\ X=1951318268\) required |

Both roots in each exceptional field are displayed,
so this also covers residues in algebraic extensions.
At the last row the actual value \(g_{h+1}(X)\) is
\(1059228468\), a contradiction. Its summands \(a_i\) start
at \(a_0=1\) and satisfy
\[
                   a_{i+1}=X(1-4i^2)a_i\pmod p .
\]
In particular \(a_{h+1}=0\). The
[modular recurrence checker](verify_binary_even_quadratic_first_position_big_modular_control.py)
sums \(a_0,\ldots,a_h\) by polynomial blocks. For a block
of multipliers, keep its product and its sum of prefix
products omitting the full product. Two successive
blocks with data \((P,S),(Q,T)\) combine as
\((PQ,S+PT)\). Applying this recursively with a symbolic
starting index gives polynomials that are evaluated at all
block starting indices by a product and remainder tree.
This performs the one modular calculation exactly.
