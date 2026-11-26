# Half-price counts on long affine lines

An affine line in the plane of primitive sums has a uniform
half-price count, with no restriction on its index length.
Normalizing each coordinate by the gcd of its intercept and
slope fixes the gcd in every two-coordinate difference.
The arbitrary-coefficient fixed-shift bound then applies,
and the remaining small projective quotients contribute at
most one point. This also allows the divisor-lattice argument
to use translated boxes and an explicit successive minimum.

## Uniform fixed lines

**LA1.** Fix \(h_0,H,A>0\). Put \(n=X^L\) and let
\(0\le r_T\le A\) on the three roles \(T=A,B,C\).
Count positive primitive rows \(A+B=C\) with
\(X^{h_0}\le C\le X^H\) on a prescribed real affine line
in this plane. One may place the original caps
\(\operatorname{rad}(T)\le n^{r_T}\) on all three roles
and impose any original eligibility restrictions.
With an integral line parametrization
\[
 T(j)=b_T+a_Tj,\qquad g_T=\gcd(a_T,b_T),
\]
one may instead cap \(\operatorname{rad}(T/g_T)\le n^{r_T}\)
in the count. Let \(r_*\) be the least \(r_T\) of a
nonconstant role. For every \(\eta>0\), for all sufficiently
small fixed \(L\) and every \(\varepsilon>0\), the count is
\[
                     O(n^{r_*/2+\eta}X^\varepsilon).
                                                               \tag{LA1}
\]
The constants are uniform in the line and the caps within
the stated ranges. At most one positive primitive row lies
on a line through the origin.

If fewer than two rows occur the claimed bound is automatic.
Keep any specified integral parametrization and count its
integer parameters. For a bare real line take the first
point as an intercept and the primitive integral direction
between it and a second point as the slope vector. All
integer points on that line then have an integer parameter.
In either setting translate the parameter to put an eligible
point at zero. The slopes have magnitude \(O(X^H)\) by
comparing two eligible points at distinct integer parameters.
The divisors \(g_T\) are pairwise coprime, since they divide
the corresponding coordinates of a primitive row. In a
nonconstant role write \(T'=T/g_T=b'_T+a'_Tj\), with
\(\gcd(b'_T,a'_T)=1\). In particular
\(\gcd(T',a'_T)=1\). If a role is constant, fix that
complete value and apply (MCB2) to the original additive
law, using \(g_T\) as fixed coefficients of the two
variable parts. The largest complete value has the
required height floor. The coefficient and root-record
losses from anatomic factorization fit within an arbitrarily
small prescribed allowance. This gives the desired half
cost of the cheaper varying role.
[Arbitrary-coefficient two-cheapest count (MCB2)](ABC_MULTIPLE_SHORT_COFACTOR_BOUNDARY.md#bounded-group-dictionaries-with-arbitrary-fixed-coefficients)
· [Anatomic factorization, Lemma 2.2](https://arxiv.org/html/2410.12234v2#S2)

Assume all roles are nonconstant and the line avoids
the origin. For any pair \(U,V\) their line projection
avoids the origin as well, since that coordinate map is
invertible on \(A+B=C\). Put \(h=\gcd(a'_U,a'_V)>0\).
The exact fixed-difference equation is
\[
          (a'_V/h)U'-(a'_U/h)V'
                  =(a'_V b'_U-a'_U b'_V)/h\ne0 .
\]
The two weighted terms have coprime absolute values.
In fact
\[
              \gcd(a'_V U',a'_U V')=h.
\]
For a prime dividing both slopes, the two primitive
cofactors \(U',V'\) are units. A prime dividing just
one slope cannot divide the weighted term using the
opposite slope and that slope's own cofactor. Primes
dividing neither slope could occur on both sides only
if they divide both primitive cofactors. This proves the
identity, including the full valuation at each prime.

Take the constant \(C_0\) from (MCB2) at the required
allowance and with complete polynomial height cap at least
\(3H\). If some normalized complete coordinate \(T'\) is
at least \(n^{C_0}\), use a pair containing that role and
a cheapest role. If they are the same, choose either of
the other two roles. The two weighted terms and the
nonzero shift, after permuting signs into a positive
three-term sum, have maximum at least \(n^{C_0}\).
Their two free radical caps are the original normalized
caps of the selected roles. Apply (MCB2) with the shift
as a fixed letter; once the weighted terms are fixed the
original two roles and then the third are determined.
There are only three possible large-role choices. This
counts these branches at the cost in (LA1).

There is at most one row for which all three \(T'\)
are less than \(n^{C_0}\). If there were two, put the
first at \(j=0\). The three normalized slopes would
be bounded in magnitude by \(n^{C_0}\), since they are
differences of the two normalized coordinate vectors
divided by a nonzero integer parameter difference.
The normalized intercept vector has the same bound.
These vectors are independent and their cross product
is a nonzero integer multiple of \((g_A,g_B,-g_C)\).
This vector is primitive by the pairwise coprimality
of the \(g_T\). Hence all \(g_T\le2n^{2C_0}\).
The first row has \(C<2n^{3C_0}<X^{h_0}\) at small
enough fixed budget, a contradiction. On an origin
line primitivity gives at most one point by direct
rational reconstruction of the common slope. This
proves all the cases.

In particular, under balanced individual caps
\(n^{1/3+o(1)}\), a prescribed moving affine line has
count exponent at most \(1/6+\eta\). The same holds
after conditioning on any additional divisors that
are incorporated in the coordinate gcds. A family
covered by a dictionary of lines of exponent
\(f<1/6\) by a fixed margin therefore admits a strict
saving, independent of complete multiplicity. The
general bounded-degree curve estimate (BC2) gives
an absolute gap for degree one as well. (LA1) makes
the line exponent explicit and depends only on the
individual caps and the actual changing coordinates.
[Moving polynomial curves (BC2)](ABC_BOUNDED_DEGREE_CURVES_AND_SECANT_GCD.md#a-uniform-saving-on-moving-curves)

## Arbitrarily eccentric and translated boxes

**LA2.** Let a rank-two lattice \(\Lambda\) of complete
coordinate pairs contain the eligible rows, in a fixed
real box
\[
                  \|\mathcal S(U,V)-z\|_\infty\le1 ,
\]
where \(U,V\) are two original roles and \(\mathcal S\)
is invertible. Let \(\ell_1,\ell_2\) be the lengths of
a Gauss-reduced basis of \(\mathcal S\Lambda\), ordered
with \(\ell_1\le\ell_2\), and let \(v_1\) be its first
direction pulled back to the complete coordinates.
For \(\rho\) the least radical cap exponent of the
roles varying in direction \(v_1\), the eligible count
is bounded by
\[
        O\!\left(n^{\rho/2+\eta}(1+\ell_2^{-1})
                                    X^\varepsilon\right).
                                                               \tag{LA2}
\]
The original radical caps may be replaced by upper
bounds on the quotient radical costs on each line,
provided the same \(\rho\) works on those lines.

Project perpendicular to the first basis vector.
A reduced companion has perpendicular length
\(\asymp\ell_2\), so at most \(O(1+\ell_2^{-1})\)
parallel affine lines can meet the real box. Pull
them back and apply (LA1) independently of the
index width. If each pulled-back line contains at
most one eligible integer parameter, the same bound
holds directly. When the affine index interval on
each line has length \(\ll n^{L^{-\xi}}\) for the
appropriate fixed \(\xi>0\), use (DL1) in place of
(LA1) to remove the factor \(n^{\rho/2}\). In
particular this works if
\(\ell_1\ge n^{-L^{-\xi'}}\) with \(0<\xi'<\xi\).
In the homogeneous setting the origin-line
argument of (DL2) may give the still stronger
bound one when \(\ell_2\) exceeds a fixed absolute
constant.
[Short-index lattice lines (DL1–DL2)](ABC_DIVISOR_LATTICE_LINE_OCCUPANCY.md#dl1-a-short-lattice-line-has-subpower-occupancy)

For example, in any such fixed shape, an actual
record with at least \(n^r\) rows and \(r>\rho/2+\eta\)
forces a second independent lattice direction of
length \(O(n^{\rho/2+\eta-r}X^\varepsilon)\) in the
scaled coordinates. This is a restriction on lattice
shape in addition to the determinant bound of (DL4).
For joint records of exponent \(f\), writing
\(v_2=\max(0,-\log_n\ell_2)\), the count exponent is
at most \(f+\rho/2+v_2+\eta\) when these bounds are
uniform on the record dictionary. Thus a fixed
saving can also follow from a sparse transverse
direction when the determinant itself falls far
below the area threshold. This conclusion keeps the
original additive law and all three radical caps.

## Explicit dual frequencies for divisor records

**LA3.** In a three-divisor record, set \(D_T\mid T\)
on the three roles. Any nonzero functional with
coefficient vector
\[
   \mathbf u=
        \left(\frac{i}{D_A}+\frac{h}{D_C},
              \frac{j}{D_B}+\frac{h}{D_C}\right),
          \qquad i,j,h\in\mathbb Z,
\]
takes integer values on the complete lattice. Its range on
the translated box of (LA2) has length at most
\(2\|\mathcal S^{-T}\mathbf u\|_1\). Thus the record has
count at most
\[
            O\!\left(n^{\rho/2+\eta}
                    (1+\|\mathcal S^{-T}\mathbf u\|_1)
                              X^\varepsilon\right),
\]
with \(\rho\) the least cap exponent of the changing
coordinates on the level lines of this functional.
In fact the displayed vectors make up the entire dual
lattice. To see this, start with the rectangular dual
lattice of covolume \(1/(D_A D_B)\). The class of
\((1,1)/D_C\) has order \(D_C\) modulo it by pairwise
coprimality. Thus adjoining it gives the requisite
covolume \(1/(D_A D_B D_C)\), and every vector pairs
integrally with the actual complete-value lattice.
The count follows by grouping the original rows by
the integer value of the functional and applying (LA1).

For the central coordinate scale \(A,B\ll X\), the
smallest \(\ell^1\) range bound in a nonrectangular
dual coset with \(1\le h<D_C\) is, up to an absolute
constant,
\[
       X\left(\frac{\|hD_A/D_C\|}{D_A}
                    +\frac{\|hD_B/D_C\|}{D_B}\right).
                                                               \tag{LA3}
\]
The double bars on a scalar mean distance to the
nearest integer. The rectangular coset contains the
two coordinate functionals. Consequently, a record
dictionary of size at most \(n^f\), on which some
functional has range \(O(n^s)\) and cap \(\rho\) as
above, gives exponent \(f+s+\rho/2+\eta\).
This criterion can be checked on the actual divisors
and can apply independently of their product height.

Here is a restriction on deriving such dual records from
raw prices alone. Consider the height-one balanced limit
with full raw cost on any selected prime power, total
radical mass \(d=1/3\) on each role, and all the numerical
height represented by finite scaled degrees in an interval
\([a,b]\) with \(b\le 2a\) and \(a>0\).
For a fixed divisor selection let \(t_T\) denote its raw
radical mass and \(w_T\) its numerical weight on \(T\).
Then \(w_T\le b t_T\) and \(1-w_T\ge a(d-t_T)\).
When a nonrectangular functional in (LA3) has range
\(n^{O(1)}\), neither of the two distances can vanish.
Both are at least \(1/D_C\) by coprimality, so
\(w_A+w_C\ge1\) and \(w_B+w_C\ge1\) in the limit.
It follows that
\[
             t_A+t_B\ \ge\ {2(1-w_C)\over b}
                        \ \ge\ {2a\over b}(d-t_C),
\]
and the total price is at least \(d\). A rectangular
functional of that range would fix a full divisor of
numerical weight one on some role, again costing \(d\).
This covers the formal \([2,4]\) example. Discounts in
the actual joint record price, numerical height outside
these bands, and noncentral boxes can change the test.
[Raw diffuse record models](ABC_BALANCED_FIXED_SHIFT_REVERSE_TESTS.md#surviving-conditional-profiles-across-the-full-range)

## Inverse-square content in translated boxes

**LA4.** Let \(\mathcal F\) be a family of \(M\)
eligible rows in the translated box of (LA2), and
put \(\mathcal B=|\det\mathcal S|^{-1}\).
In any two complete coordinates write \(P_i\)
for its integer vectors. Let \(r_{\rm med}\) be
the median of the three individual radical caps
in normalized exponents. If \(M\ge n^r\) with
\(r>r_{\rm med}/2\) by a fixed margin, then the rows
have affine rank two and
\[
 G_{\rm aff}:=\gcd_{i,j,k}
  |\det(P_j-P_i,P_k-P_i)|
       \ \ll\ \mathcal B\,M^{-2}n^\eta X^\varepsilon .
                                                               \tag{LA4}
\]
The quantifiers are height-first, with the small
budget chosen after the stated margin and \(\eta>0\).
The gcd can be computed in either pair of complete
coordinates.

Indeed any line has at least two varying roles, so
its cheapest varying cap is at most \(r_{\rm med}\).
By (LA1) the family cannot be collinear. Fix \(P_0\)
and generate a rank-two lattice by all differences
\(P_i-P_0\). Its determinant is \(G_{\rm aff}\).
The determinant of a basis is an integer combination
of pair determinants of the anchored differences;
conversely it divides every such determinant, and
bilinearity makes the gcd of all triangle minors
equal to that anchored gcd. All the points lie in
the corresponding affine lattice coset.

Write \(\Delta=G_{\rm aff}/\mathcal B\). Its reduced
companion length \(\ell_2\) in the transformed
coordinates is bounded above by an absolute
constant. Otherwise the perpendicular projection
shows that the box meets at most one line in the
first lattice direction, contradicting (LA1).
If \(\Delta\ge n^{-L^{-\xi_1}}\), use the same
slicing as in (DL2): there are \(O(1+\Delta^{-1/2})\)
lines and their parameter widths are
\(O(\Delta^{-1})\). Unlike on a single line of large
transverse spacing, the widths are now in the
short-index range of (DL1), after decreasing
\(\xi_1\). Counting all lines gives
\(M\ll n^{\eta'}X^{\varepsilon'}\Delta^{-1/2}\);
one may absorb the extra one since
\(\Delta\le8\) by taking any nonzero triangle
minor in the translated box. If
\(\Delta<n^{-L^{-\xi_1}}\), the same count
inequality is weaker than the trivial count
from the radical caps of two roles, for sufficiently
small budgets. Rearrange the bound and resize the
allowances to obtain (LA4).
[Homogeneous predecessor (DL4)](ABC_DIVISOR_LATTICE_LINE_OCCUPANCY.md#dl4-the-gcd-of-the-actual-projective-minors)

In particular, for actual rows in any one coset of
a rank-two congruence lattice of index \(D\), the
translated-box count is
\[
    O\!\left(n^\eta X^\varepsilon
             \max\{n^{r_{\rm med}/2},\sqrt{\mathcal B/D}\}\right).
\]
Use the half-cap floor for a small family. Otherwise
the lattice of differences is a sublattice of that
congruence lattice, so \(G_{\rm aff}\ge D\) and
(LA4) applies. This supports inhomogeneous residue
records as well as the projective classes in (LC1).
[Projective class occupancy (LC1)](ABC_PROJECTIVE_CONTENT_MOMENTS.md#projective-classes-from-one-congruence)
