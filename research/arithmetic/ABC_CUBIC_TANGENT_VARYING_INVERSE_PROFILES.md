# Varying inverse weighted profiles with order-three ideal classes

The virtual product construction extends to varying coefficient
classes. An elementary separation argument keeps the unordered
inverse coefficient pairs distinct. A fixed elliptic quotient
then controls the number of times the fields recur, uniformly
even as those profiles vary.

## A rational source of inverse weighted profiles

**Proposition VI1.** For any rational \(t\ge2\), put
\[
 A=t^3/2,\qquad B=(A^2+1)/(2A),\qquad b=A^2-1.
\]
There are weighted profiles with rational deck halves and
coefficient classes \(b,b^2\) modulo cubes. On the doubled
genus-five model with this \(b\), the specialization
\[
 s=\frac{70}{35t^4+2}
 \tag{1}
\]
gives two positive directions. If \(t\) is even at \(2\)
and integral at \(5,7\), the field is imaginary and the
two cubic norm classes span a two-dimensional space whose
virtual-unit line maps injectively into the ideal class
group.

For a rational branch parameter \(r\) in
[GT1](ABC_RECIPROCAL_GENUS_TWO_TANGENT_OBSTRUCTION.md#which-weighted-profiles-have-a-rational-branch-point),
write \(Z=3r/(r+2)\). The coefficient of its weighted elliptic
model satisfies the exact identity
\[
 4k(k+1)=(Z^2-1)\left(-\frac{r+1}{r^2}\right)^3 .
 \tag{2}
\]
Now take \(Z=A,B\) for the two profiles, using
\(r=2Z/(3-Z)\).
Both inputs are admissible: \(A\ge4,\ B>2\), and \(B=3\)
would require \(A=3\pm2\sqrt2\).
The first coefficient in (2) has class \(b\), and the
second has class \(b^2\) because
\[
 B^2-1=\frac{(A^2-1)^2}{4A^2}=\frac{b^2}{t^6}.
\]
Neither coefficient is a rational cube, since each weighted
profile admits a rational deck half and no symmetric
profile admits one by GT5.
In particular the normalized coefficients fit DB1.

With (1) the new parameter \(w=b^2s^3\) is
\[
 w=\frac12
       \frac{(1-4/t^6)^2}{(1+2/(35t^4))^3},
 \qquad \frac13<w<\frac12\quad(t\ge2).
\]
Thus \((w^2-1)(w^2+4w-1)<0\). The two positive
directions follow from the real abscissa argument of
[DB3](ABC_CUBIC_TANGENT_DOUBLED_BRANCH_FAMILY.md#simultaneous-positive-independent-directions).

For \(t\) even at \(2\), \(b=(t^6-4)/4\) is a dyadic unit
and so is \(s\). This makes \(w\) a dyadic unit. For \(t\)
integral at \(5,7\) the coefficient \(b\) is integral and
\(s\) is divisible by each prime, so \(w\) has positive
valuation at both. At either place the quadratic
radicand is a square, since it is congruent to one modulo
that odd prime. At \(7\), choose \(Y\equiv1\bmod7\); the
value of the product function \(\gamma\) of VP1 is then
\(-2\bmod7\), a noncube. This local residue argument does
not require the prime to divide the coefficient, as it only
uses \(w\equiv0\). By
[VP1](ABC_CUBIC_TANGENT_VIRTUAL_PRODUCT_FAMILY.md#the-product-and-its-dyadic-condition),
\(\gamma\) is a noncube virtual unit. The field splits at
\(5\) and is therefore distinct from \(\mathbb Q(\sqrt{-3})\).
Its units are cubes, giving the nonzero class-group
direction. Finally \(b\) is noncubic and a dyadic unit,
so an odd prime has noncubic coefficient valuation.
Both individual norms fail there by DB2. This proves
independence and identifies the virtual-unit line.

## A count with both fields and profiles distinct

**Theorem VI2.** At least
\[
 X^{1/36}\exp\!\left(-O\!\left(\frac{\log X}{\log\log X}\right)\right)
\]
distinct imaginary quadratic fields of absolute discriminant
at most \(X\) occur in VI1 with pairwise distinct unordered
inverse coefficient-class pairs. The estimate counts both
objects in the same family, with a uniform field-multiplicity
bound that does not depend on \(b\).

Choose \(t=a/c\) in lowest terms, \(H\le a\le2H\) with
\(a\) even, and \(1\le c\le H^{1/3}/4\) with \(c\)
coprime to \(70\). There are \(\gg H^{4/3}\) such pairs.
For example counting \(a\) of the stipulated parity
coprime to \(c\) gives its interval length times
\(\phi(c)/(2c)\), with error \(O(\tau(c))\).
The mean of \(\phi(c)/c\) on denominators coprime to
\(70\) is bounded below, by the divisor expansion.

For a fixed positive cube-free integer coefficient
representative \(k\), one would have
\[
 \frac{a^6-4c^6}{4}=k v^3,\qquad v\in\mathbb Z_{>0}.
\]
Here \(b=(a^6-4c^6)/(4c^6)\) has exactly the same cube
class as the integer on the left.
The ratios \(v/a^2\) for that \(k\) lie within
\(4(1/4)^6 H^{-4}\) of \((4k)^{-1/3}\), by taking cube
roots and using \(a\ge H\) and \(c\le H^{1/3}/4\).
Two distinct such rationals would be separated by at least
\(1/(16H^4)\), since their denominators are at most
\(4H^2\). All the ratios must therefore be equal.
The displayed equation then fixes \(c/a\), and
primitivity fixes \(a,c\).
Thus each positive coefficient class occurs at most once
throughout the skew box.

The rational function \(w(t)\) has degree at most twelve,
so its logarithmic height is at most \(12\log H+O(1)\).
Clearing the fourth power of its denominator in the common
quartic in \(w\) gives \(|d|=O(H^{48})\).
To bound field multiplicity, use the point on the **fixed**
\(H_3:U^2=X^3+X^2-X\) from DB6. Its abscissa
\(X=(w+1)/(w-1)\) has logarithmic height \(O(\log H)\).
Different parameter values can yield the same \(X\)
at most twelve times, since \(w(t)\) is nonconstant.
The proof of [DM2](MORDELL_CUBIC_DOUBLED_SIMULTANEOUS_TWISTS.md#counting-distinct-simultaneous-twists)
applies to this fixed elliptic curve. The point is
anti-invariant, the cubic norm gives the same canonical
height separation \(\gg\log(2|d|)\), and the fixed
two-descent algebra is \(\mathbb Q\times\mathbb Q(\sqrt5)\).
This gives at most
\(\exp(O(\log H/\log\log H))\) parameters per field with
constants independent of the varying coefficients.
Select one parameter per field. Passing from ordered
coefficient classes to unordered inverse pairs loses at
most a factor two. Taking \(H\) a sufficiently small
constant times \(X^{1/48}\) proves the bound.

The [exact checker](verify_cubic_tangent_varying_inverse_profiles.py)
replays the coefficient and inverse identities and a
specialization, including the dyadic unit and noncube residue.
This count concerns distinct quadratic fields and coefficient
classes, with an actual order-three ideal class in each field.
The individual weighted norms still fail the local
valuation condition for genuine triangles.
