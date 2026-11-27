# A normalized cusp relation and finite indicial occupancy in degree forty-two

**Theorem PolynomialWeylFortyTwoNormalizedCuspRelation.** Take
the centered actual degree-forty-two Fourier residual and its
[regular finite pole](POLYNOMIAL_WEYL_FORTY_TWO_REGULAR_FINITE_POLE.md).
Use \(t=X^{1/3}\), \(E=tY\), \(C=XY\), and
\(H(E)=(E^3-u)^4(E^3-v)\) as in
[the root constraints](POLYNOMIAL_WEYL_FORTY_TWO_SECONDARY_ROOT_CONSTRAINTS.md).
There exist a shifted mate \(Q^\ast=Q-D(P)/2\) with
\(\deg D\le3\) and a polynomial \(F\) of degree seven and
leading coefficient \(\kappa=\beta^2/\alpha^7\) such that
\[
             S=(Q^\ast)^2-F(P),\qquad [P,S]=-2Q^\ast .
\]
It has \(t\)-degree \(17\) and lowest Euler grade \(-51\).
Its top \(t\)-coefficient is \(H^5 B\), where
\[
           B'-\frac23\frac{H'}H B=\frac{\beta}{\alpha},
              \qquad \deg_E B=10 .
                                                               \tag{1}
\]
The \(X,Y\) support of \(S\) is subrectangular with corner
\((34,85)\), of Bernstein degree \(119\).

Let \(X^{-18}f(C),X^{-63}g(C),X^{-51}h(C)\) be the lowest
Euler components of \(P,Q^\ast,S\), respectively. They
commute and \(\deg f=30,\deg g=105,\deg h=85\). For the
twelve excess roots of \(f\) after its factor \((C)_{18}\),
the multiplicity counts in every \(\rho+\mathbb Z\) orbit
are invariant under shift by three modulo eighteen. Thus
they occupy exactly two cycles of six, counted with
multiplicity. Each cycle lies in a single
\(\rho+3\mathbb Z\) progression and has one member in
each residue modulo eighteen there. In particular every
noninteger indicial class that occurs has six or twelve
roots counted with multiplicity.

## Removing the higher \(t\)-components

Start from \(Q^2-\kappa P^7\). The \(t\)-top terms of \(P,Q\)
are \(t^6\alpha H^2,t^{21}\beta H^7\), so the degree
\(42\) component cancels. The commutator with \(P\) is
\(-2Q\). At any stage of reducing this element, let \(K\)
be its largest \(t\)-degree. If \(K>17\), the bracket of
its \(t^K A(E)\) component with \(t^6\alpha H^2\)
must vanish at degree \(K+4>21\). The first contraction
formula gives
\[
                  A'=\frac K3\frac{H'}H A .
\]
Because \(H\) has simple roots, \(K=3j\) with \(j\) an
integer, and \(A\) is proportional to \(H^j\). Since
\(17<K<42\), one has \(6\le j\le13\). Write
\(j=2a+7b\) with \(a,b\ge0\), necessarily \(b\le1\).
A scalar multiple of \(P^a Q^b\) cancels the top
component. When \(b=1\), its commutator with \(P\) is
\(-P^a\), whose \(t\)-degree is at most \(18\) and
whose lowest Euler grade is at least \(-54\).
For \(b=0\) the commutator is zero.
Thus the degree \(21\) and grade \(-63\) coefficients
of \([P,S]\) are those of \(-2Q\) throughout the
reduction. The procedure terminates at \(t\)-degree
at most \(17\), which must be exactly \(17\) because
the commutator still has degree \(21\).

The removed terms have the ordered form \(D(P)Q+G(P)\)
with \(\deg D\le3\). Since \([Q,D(P)]=D'(P)\), completing
the square with \(Q^\ast=Q-D(P)/2\) gives the claimed
form of \(S\), with
\[
                   F(P)=\kappa P^7+G(P)-D'(P)/2+D(P)^2/4 .
\]
The shift in \(Q\) preserves both its \(t\)-top and its
lowest Euler component.

Let \(A=S_{17}(E)\). The degree \(21\) commutator equation
gives
\[
                  A'-\frac{17}3\frac{H'}H A
                       =\frac{\beta}{\alpha}H^5 .
\]
At a root of \(H\) of multiplicity \(\mu=1\) or \(4\),
the leading coefficient on the left for an initial order
\(j\) has factor \(j-17\mu/3\), which is nonzero.
Thus the order of \(A\) there is \(5\mu+1\).
Consequently \(A=H^5 B\) with \(B\) polynomial, vanishing
at all six distinct roots of \(H\). Its equation is (1).
Since \(\deg H=15\) and \(\deg B\ge6\), the highest degree
on the left of (1) would be \(\deg B-1>0\) with a
nonzero coefficient unless \(\deg B=10\).
Hence \(\deg_E A=85\).

## What the top \(t\)-coefficient forces at zero

The term \(t^{17}E^{85}\) in the normal form comes from
the ordinary monomial \(X^{34}Y^{85}\), since corrections
in the normal change of variables only lower \(t\)-degree.
It has Euler grade \(-51\). Let \(a\le-51\) be the true
lowest grade of \(S\), with component \(X^a h(C)\).
The grade \(a-18\) is below the lowest grade \(-63\) of
\([P,S]=-2Q^\ast\). Thus the lowest components commute:
\[
                 f(C+a)h(C)=f(C)h(C-18).
\]
Writing \(e=\deg h\) and comparing coefficients one
below the leading degree gives \(30a=-18e\).
The \(t\)-degree of the highest \(Y\)-exponent monomial
of \(X^a h(C)\) is \(3a+2e=-a/3\). The \(t\)-degree
cap \(17\) forces \(a\ge-51\), giving \(a=-51\),
\(e=85\). The bounds \(i-j\ge-51\) and \(3i-j\le17\)
on any \(X^iY^j\) in \(S\) imply \(i\le34,j\le85\).
The corner term exhibits equality.

Let \(P_0,Q_0,S_0\) denote the three lowest components.
The equality \([Q^\ast,S]=-F'(P)\) also shows that
\([Q_0,S_0]=0\), since the lowest possible grade of
the left is \(-114\), while the right has lowest grade
\(-108\). Alternatively homogeneous centralizers of
\(P_0\) commute in the rational Euler algebra. In
particular,
\[
                 S_0^3=c\,Q_0P_0^5,\qquad c\ne0.            \tag{2}
\]
Both powers at grade \(-153\) commute with \(P_0\).
The ratio of their coefficient polynomials is rational
of period \(18\), so is constant by the same displayed
commutation equation.

For root occupancy, take the divisor of
\[
                 f(C-51)h(C)=f(C)h(C-18).
\]
Within each \(\rho+\mathbb Z\), project the multiplicities
to indices modulo eighteen. The \(h\)-terms cancel, so
the \(f\)-record is invariant under addition of \(51\),
equivalently of \(3\).
The base falling factorial \((C)_{18}\) has invariant
record as well, and the remaining record has total size
twelve. Each nonzero cycle of addition by \(3\) has
length six. This proves the occupancy statement.

## Two simultaneous finite root inequalities

Remove the falling factorials \((C)_{18},(C)_{63},(C)_{51}\)
from \(f,g,h\), respectively, and discard their leading
scalars. In a fixed \(\rho+3\mathbb Z\)-orbit write
\(U(z),V(z),W(z)\) for the corresponding excess-root
multiplicity records in step-three coordinates. The
commutation equations cancel the falling factorials and
give
\[
    V(z)=\frac{1-z^{21}}{1-z^6}U(z),\qquad
    W(z)=\frac{1-z^{17}}{1-z^6}U(z).                       \tag{3}
\]
For each residue \(i\bmod6\), let \(r_{i,t}\) be the
\(t\)-th member of the sorted multiset of excess-root
indices there. The counts are equal across residues by
the occupancy statement. For any \(b=17\) or \(21\),
nonnegativity of the coefficients of
\((1-z^b)U(z)/(1-z^6)\) says that the number of roots
in residue \(i\) through any index \(n\equiv i\bmod6\)
is at least the number in residue \(i-b\) through
\(n-b\). Equivalently
\(r_{i,t}\le r_{i-b,t}+b\) for each rank \(t\).
Thus the two polynomial compositions give
\[
              r_{i,t}\le r_{i+1,t}+17,\qquad
              r_{i,t}\le r_{i+3,t}+21.                   \tag{4}
\]
Conversely these inequalities give nonnegative excess
records in (3). A packet of six indices at a given
rank \(t\) has width at most \(55\) in step-three
coordinates. To compare any two residues, move by a
sequence of offsets \(+1,+3\); the bounds for offsets
\(1,2,3,4,5\) can be taken to be \(17,34,21,38,55\),
respectively. This constrains the six-root packets
without fixing their locations in \(\mathbb C\).
