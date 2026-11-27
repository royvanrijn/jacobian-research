# Two weighted tangent profiles in the same quadratic fields

A rational translation of a weighted tangent elliptic cover
changes the tangent profile and preserves the field squareclass.
For two explicit profiles the resulting six prime forms
are jointly admissible. This gives infinitely many fields
with two genuine coherent triangles on different lines.
The triangles share their maximal-order cubic class.
Their order-three ideal differences fail to remain
three-torsion in the common order with enlarged conductor.

## Translating the elliptic involution

Use the fixed-profile coordinates of
[WM1](ABC_WEIGHTED_CUBIC_TANGENT_ELLIPTIC_MARKS.md#the-double-cover-and-its-three-marks).
For \(k\ne0,-1\), put \(w=k(k+1)\).
The model \(E_w:Y^2=X^3+16w^2\) has
deck involution \(T-P\) with
\(T=(4w,-4w(2k+1))\).
For \(k=b/a\) in lowest terms with \(a>0\),
choose the integer normal \((a,b,-a-b)\)
and weights \((s_1,s_2,s_3)=(a,|b|,|a+b|)\).
Then the homogeneous prime markings
for \(t=u/h\) are
\[
 (p_1,p_2,p_3)=(u,\ u-(k+1)h,\ u-kh).
\]
Here \(h\) may be rational. The normalized
coefficients of the reconstruction are
\((x,y,z)=(|a+b|,a,|b|)\).

**Proposition XT1.** Suppose \(A\in E_w(\mathbb Q)\)
and \(T_1=T+2A\) is finite with nonzero
abscissa. Set
\[
 c=\frac{X(T_1)}{4w},\qquad
 k_1=\frac{-Y(T_1)/(4w)-1}{2},\qquad
 w_1=k_1(k_1+1).
\]
Then translation by \(A\) followed by
\((X,Y)\mapsto(c^2X,c^3Y)\) identifies the
two elliptic covers with their involutions.
It induces a rational projective
transformation of the parameter lines.
Under the induced map, the two adjacent
edge classes \([\alpha_{12}]\) agree in
the common function field modulo cubes,
with consistent choices of orientation.

The curve equation gives \(w_1=w c^3\).
The group isomorphism \(\iota_c\) from
\(E_w\) to \(E_{w_1}\) sends \(T_1\)
to its deck center and
\(D=(0,4w)\) to \(D_1=(0,4w_1)\).
Thus \(P\mapsto\iota_c(P+A)\)
commutes with the deck maps.
The resulting isomorphism of the
invariant rational function fields
is projective linear on the
parameter lines.

For the norm class, work in the
affine parameter chart \(h=1\) of
the first weighted profile, with
\(U^2=f_k(t)\) as in WM1.
Write
\[
 A_1=t^3,\quad A_2=k^2(t-k-1)^3,\quad
 A_3=(k+1)^2(t-k)^3.
\]
The adjacent norm generators in this
chart, up to a common orientation, are
\[
 \alpha_{12}=\frac{A_1+A_2-A_3+wU}{2|k|},\qquad
 \alpha_{23}=\frac{A_2+A_3-A_1+wU}{2|w|}.
\]
Indeed dividing the Heron expression
by the square of the coefficient
product gives
\(-E=f_k(t)/a^2\) in this chart.
Choose the signed square root
\(a^2 wU/(a|b(a+b)|)\);
the reconstruction formulas give
the displayed generators.
The norm of \(\alpha_{12}\)
is \(t^3(t-k-1)^3\).
At \(t=0,U=w\) the generator
vanishes. At \(t=k+1,U=-(k+1)\)
it also vanishes. The other sheets
at those points are nonzero.
The two zeros have order three.
The two poles of \(t\) each give
a pole of order three, since the
leading trace coefficient
\(1+k^2-(k+1)^2=-2k\)
is nonzero, whereas \(U\) has
pole order two. The two zero
points are \(O,T-D\), and the
two poles of \(t\) sum to \(T\).
Thus \(\operatorname{div}(\alpha_{12})/3\)
has class \(-D\) on \(E_w\).

The same description on the second
profile gives class \(-D_1\).
The isomorphism \(\iota_c\) matches
these classes and a translation is
trivial on \(\operatorname{Pic}^0\).
Consequently the quotient of the
pulled-back generators is a cube
times a rational constant, by descent
of principal divisors. Its norm to
the parameter line is a cube.
This makes the square of the
constant a rational cube, hence
the constant itself is a cube.
Specialization at nonzero norm
generators preserves the equality.
One may multiply either of the
actual homogeneous prime markings
by a rational scaling factor:
the generators scale by its cube.

## Six jointly admissible prime markings

**Theorem XT2.** There are
\(\gg X^{1/2}/(\log X)^6\) distinct
imaginary fields of absolute
discriminant in \([X,2X)\) with
genuine coherent triangles on both
fixed profiles
\[
\begin{array}{c|rrr|rrr}
 &a_1&a_2&a_3&x&y&z\\ \hline
 1&1&2&-3&3&1&2\\
 2&55566&-4913&-50653&50653&55566&4913
\end{array}
\]
for every sufficiently large \(X\).
The six prime vertices are distinct
and at a common scale \(X^{1/4}\)
up to fixed constants. Both triangles
give the same nonzero maximal-order
virtual cube class, and its associated
ideal class has exact order three.

Take \(k=2\), \(T=(24,-120)\)
on \(Y^2=X^3+576\) and translate
by \(A=T\). Then
\[
 3T=\left(-\frac{2516}{441},-\frac{182960}{9261}\right),
 \qquad k_1=-4913/55566.
\]
Evaluation on the three rational
marks of WM1 gives the parameter
map
\[
 t\longmapsto
 \frac{1186923t-1978205}{-3333960t+13724802}.
\]
Put \(c_0=3\cdot7\cdot17\cdot37=13209\)
and \(g=c_0^2=174477681\).
For the first profile the
markings are \(u,u-3h,u-2h\).
For the second they are
\[
 \frac{1186923u-1978205h}{g},\quad
 \frac{4226103u-14489496h}{g},\quad
 \frac{892143u-764694h}{g}.
\]
Indeed their first coordinate
is \(u_1\), and their step
\[
 h_1=-\frac{55566}{g}(60u-247h)
\]
makes all three target
mark relations hold.
Then \(u_1/h_1\) is the
displayed projective transform.
In homogeneous form, writing
\(F_k(u,h)=h^4 f_k(u/h)\),
one has
\[
 F_{k_1}(u_1,h_1)=F_2(u,h)/g .
 \tag{1}
\]
The constants in (1) follow
by substitution, or by the
square factor in XT1 and
evaluation at a single
ordinary parameter.

Choose a small fixed rectangle
about \((u,h)=(Z,-Z)\).
All six forms are positive
of scale \(Z\) with
pairwise distinct values
there (the six values at
the center are pairwise
distinct). The six linear
directions are pairwise
nonproportional.
The quartic \(F_2\) is negative
there and has magnitude
comparable to \(Z^4\).

Use a congruence cell with
the following projective
residues. The three target
valuations in the last
column are for their
numerators before dividing
by \(g\). The source
valuations are all zero
in each row.
\[
\begin{array}{c|c|r|r|c|c}
 q&\text{modulus}&u&h&
     v_q(\text{target numerators})&v_q(F_2)\\ \hline
 2&2&1&0&(0,0,0)&0\\
 3&81&1&15&(2,2,2)&3\\
 5&5&1&1&(0,0,0)&0\\
 7&2401&18&1&(2,2,2)&3\\
 17&4913&86&1&(2,2,2)&2\\
 37&50653&392&1&(2,2,2)&2
\end{array}
\]
The target numerator row
contents are \(17^2 37^2\),
\(3^2 7^2 37^2\) and
\(3^2 7^2 17^2\).
Fix any unit-scaled
representative of each row.
The Chinese remainder theorem
gives an affine lattice cell.
After the division by \(g\)
all six forms are integral
and prime to the primes
in the table throughout
the cell. For another
prime the six forms
exclude at most six
lines modulo that prime.
Each local density is
positive, since the
remaining primes are
at least \(11\).
The six-form theorem
of [Green–Tao](https://arxiv.org/pdf/math/0606088)
and [Green–Tao–Ziegler](https://arxiv.org/pdf/1009.3998)
therefore gives a positive
constant times \(Z^2/(\log Z)^6\)
prime pairs asymptotically.

We can also impose that
\(d=F_2(u,h)/g\) be
squarefree with the
same order of count.
The table already fixes
its required behavior
at every listed prime.
Outside the table the
four simple projective
roots of \(F_2\) and
the six prime-form
lines exclude at most
ten directions modulo
\(q\ge11\). Thus for
each fixed finite list
of primes outside the
table we may avoid
their squares with
positive local factors.
The local defects are
\(O(q^{-2})\), as in
[CT4](ABC_CUBIC_TANGENT_NORM_RANK_AND_CONDUCTORS.md#squarefree-dyadic-fields-and-rank-two-conductors).
The same squarefree
sieve applies. In the
range from a fixed \(Y_0\)
to \((\log Z)^8\)
the six-form Selberg
upper bound has sieve
weights \(6/q+O(q^{-2})\),
and the \(O(q^2)\) bad
residue pairs modulo
\(q^2\) give
\(O(Z^2/(Y_0(\log Z)^6))\)
plus a power-small error.
The pairs \((u,h)\)
are primitive when the
three source values
are distinct primes.
Above \((\log Z)^8\)
the primitive root-lattice
and Helfgott quartic
bounds of CT4 give
\(o(Z^2/(\log Z)^6)\).
Thus the total squarefree
count has the claimed
lower order.
By homogeneity, shrinking
the rectangle and fixing
a suitable constant in
\(Z\asymp X^{1/4}\)
puts \(|d|\) in \([X,2X)\).
The chosen cell gives
\(d\equiv1\bmod4\).

For either profile, the
three full trace congruences
follow from \(\sum a_i p_i=0\)
and \(\sum a_i=0\), as
in WM4. For a large
prime \(p_i\), the
remaining weighted
Heron inputs are units,
and their ratio modulo
\(p_i\) is \(-a_l/a_j\)
with \(i,j,l\) distinct.
Their difference is a
unit because \(a_i\ne0\)
is a fixed integer.
Thus the Heron parameter
is a unit at every
incident prime, and
both triangles are
genuine by exact
reconstruction.
The sign and field
agree by (1).
The class of XT1
is nonzero in the
maximal order: a
principal edge-quotient
ideal of norm \(p_i p_j\)
could not have a rational
generator, by its
one-sided supports.
A nonrational integral
generator has norm
at least \(|d|/4\),
larger than \(p_i p_j\).
The class therefore
has exact order three.

Finally the source
quartic is irreducible
modulo \(5\), where
after scaling it is
\(t^4+2t^2+3\).
It has trivial
rational branch
automorphism group
by the asymmetric
case of
[FM4](ABC_CUBIC_TANGENT_COLLISION_MOMENT.md#asymmetric-weighted-profiles).
The collision incidence
bound \(O_\varepsilon(Z^{39/23+\varepsilon})\)
therefore removes a
negligible fraction
of the simultaneous
prime pairs. Hence
the count holds for
distinct fields.

## The two order lifts cannot be combined without change

**Proposition XT3.** For the paired
triangles of XT2, the conductors
of their respective quadratic
orders are
\[
 f_1=c_0|h|,\qquad
 f_2=|60u-247h|/g .
\]
Outside a fixed finite set of
primes these conductors have
disjoint support. On all but
\(O(Z(\log Z)^C)\) of the
paired triangles at scale \(Z\),
each conductor has a prime
outside that fixed set.
The two displayed adjacent
edge differences give
nonzero order-three classes
in their own orders. After
lifting the same supporting
prime ideals to the order
of conductor \(\operatorname{lcm}(f_1,f_2)\),
neither displayed difference
has order dividing three.

For the first profile
\(-E=h^2 F_2(u,h)\).
For the second
\(-E_1=(h_1/a_1)^2 F_{k_1}(u_1,h_1)\),
where \(a_1=55566\) is its
positive first weight.
Together with (1) and
\(d\equiv1\bmod4\) this
gives the exact conductors.
As \((u,h)\) are primitive,
a prime dividing \(h\)
and \(60u-247h\) divides
\(60\). Outside fixed
primes the two conductors
are supported on disjoint
linear directions.
An integer of size \(O(Z)\)
supported on a fixed set
of primes has only
\(O((\log Z)^C)\) possible
values. Each of the two
linear forms has
\(O(Z)\) pairs per value.
Thus both conductors
have a prime outside
any prescribed fixed set
for all but the stated
number of pairs.
All six supporting primes
are coprime to the
common conductor at
large \(Z\). Indeed
the fixed determinant
of each prime form
with the other
profile's conductor
form is nonzero.
A prime above the
fixed determinant
divisors cannot
divide both values
at a primitive pair.

In the \(h=1\) chart
the cube root of
\(\Theta=\alpha_{12}\overline{\alpha}_{23}\)
has the form
\[
 \beta_k=\operatorname{sgn}(k+1)
    \frac{B_k(t)-U}{2},\qquad
 B_k(t)=t^2-2(k+1)t+k(k+1).
 \tag{2}
\]
This follows by cubing
and using the two
norm generators of XT1.
The homogeneous cube root
has the additional factor
\(h^2\) in that profile.
At a good prime dividing
its own \(h\), reduction
gives a nonzero rational
multiple of
\((1-\sqrt{-3})/2\).
The resulting class modulo
rational units has exact
order three. The
corresponding prime-ideal
difference is therefore
nonprincipal in the
order, since the units
in the maximal imaginary
order are \(\{\pm1\}\).
Its cube is principal
there by the genuine
norm generators.

At a prime dividing
the other profile's
conductor, the parameter
is instead specialized
to a fixed finite value.
It is \(247/60\) in
profile \(1\), and
\(-395641/1111320\)
in profile \(2\).
At these values, write
\(U=\pm v\sqrt{-3}\).
The respective values
of \(B_k,v\) are
\[
 \left(-\frac{6311}{3600},\frac{4403}{1200}\right),
 \qquad
 \left(\frac{278135623}{400075200},
       \frac{248858189}{1200225600}\right).
\]
In each case \(v(B_k^2-v^2)\)
is a nonzero rational
number. The coefficient
of \(\sqrt{-3}\) in
\(\beta_k^3\) is, up
to sign, \(3v(B_k^2-v^2)/8\).
It stays nonzero at
every cross-conductor
prime outside a
fixed finite set.
The homogeneous factor
\(h^6\) is a unit there
by the disjoint
conductor supports.
Thus the corresponding
\(\Theta\) is not scalar
modulo that prime.

If the lifted ideal
difference in the
common order had
order dividing three,
its cube would have
a generator in that
order. After extending
to the maximal order
this generator would be
\(\Theta\) times one of
\(\pm1\), since the
maximal ideal cube
is generated by
\(\Theta\). This is
impossible at the
cross-conductor prime.
The result concerns
the displayed supporting
prime ideals. The
common order can
contain other cubic
residue classes.

## Exact paired example

The projective congruence
conditions admit the primitive
pair
\[
 (u,h)=(8053580009,-8202927966).
\]
The prime triples are
\[
 (8053580009,32662363907,24459435941),
 \qquad (147789977,876281423,77131211).
\]
Their common fundamental
discriminant is
\(-529788943876447727940522414336123\).
The order conductors
are \(108352475502894\)
and \(14382\).
The first has the
exclusive good prime
\(23\), and the second
has \(47\).
In the exact basis
\(1,\sqrt d\) the
irrational coefficient
of \(\Theta_1\) is
\(8\bmod47\)
and that of
\(\Theta_2\) is
\(22\bmod23\).
Thus the same prime
supports exhibit the
failure in XT3 as
well as the shared
maximal-order class.
[The Sage checker](verify_weighted_tangent_shared_field_translations.py)
replays the parameter and
cube-root identities,
the full residue table,
primality, the common
discriminant, the
cube-class equality
and the two cross
residues.
