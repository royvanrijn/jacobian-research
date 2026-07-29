# Three-puncture nonlinear completion frontier

This note keeps the double-incidence core as a bounded moonshot.  It does
not claim an affine-space Keller completion.  It proves two dimension-free
rank-drop gates, identifies the exact collapse of the first nonlinear
five-variable repair, and records a finite six-variable screen.

Work first with \(a=b=1\).  Write

\[
c=L-1,\qquad
D_0=1-ru,\qquad
D_1=1-(r+c)v
\]

and

\[
R=\int_0^r(1-\xi u)(1-(\xi+c)v)\,d\xi .
\]

Then

\[
R_r=D_0D_1.
\]

The selected curve is still

\[
c=-1,\qquad u=r^{-1},\qquad v=(r-1)^{-1},
\]

so its coordinate ring is
\(k[r,r^{-1},(r-1)^{-1}]\).

## 1. Replacing \(u\) is not enough if \(R\) is retained

Suppose a completion keeps \(c,v,R\) as target coordinates and appends
arbitrarily many variables.  Expanding its determinant along the \(R\)-row
shows that the determinant belongs to

\[
(R_r,R_u).
\]

This ideal is proper.  For example,

\[
(c,v,r,u)=(0,3/2,1,1)
\]

annihilates both generators.  Therefore:

\[
\boxed{\text{No Keller completion in any padded dimension can keep
\(c,v,R\) unchanged.}}
\]

Thus the next search must modify the primitive coordinate as well as replace
\(u\) or \(v\).  This is stronger than the earlier plane-Keller gate.

## 2. The first nonlinear \(\mathbb A^5\) transfer

The smallest modification in the \(u\)-orientation is

\[
A=R+D_1z.
\]

Set

\[
\begin{aligned}
N={}&c^2uv^2-cruv^2-2cuv-2r^2uv^2+ruv+u+6v^2z.
\end{aligned}
\]

The polynomial map

\[
\Psi(c,v,r,u,z)=(c,v,A,r,N)
\]

has the exact determinant

\[
\boxed{\operatorname {Jac}(\Psi)=(1-cv)^3.}
\]

This is a useful transfer identity, but not a Keller map.  It collapses the
two incidence factors to the single residual factor

\[
q=1-cv.
\]

On the selected three-puncture curve,

\[
q=\frac r{r-1}.
\]

Its valuation vector is the difference of the \(r=0\) and \(r=1\)
punctures, so only a rank-one character remains visible.

The failure is not an artifact of the displayed \(N\).  If \(A\) and \(r\)
are retained and any number of variables is appended, every determinant
lies in

\[
(A_u,A_z)=(R_u,D_1).
\]

Both generators vanish on \(V(r,1-cv)\).  Hence no alternative polynomial
last coordinate, and no identity-padded sixth coordinate, can make this
orientation Keller while \(r\) remains exposed.

## 3. The opposite orientation adds a puncture

The symmetric-looking modification

\[
A'=R+D_0z
\]

has the rank-drop locus

\[
D_0=0,\qquad 1+3cu=0.
\]

Indeed the polynomial

\[
M=-3cru^2v+3cuv-2r^2u^2v+ruv+6u^2z+v
\]

gives

\[
\operatorname {Jac}(c,u,A',r,M)=1+3cu.
\]

On the selected curve this residual factor is

\[
1+3cu=\frac{r-3}{r}.
\]

It vanishes at \(r=3\).  Thus this orientation does not merely compress the
two desired characters: it introduces a fourth puncture.  It is excluded
from the bounded three-puncture search.

## 4. Finite coupled \(\mathbb A^6\) screen

The first two-modification primitive was tested in the two forms

\[
\widetilde R=R+D_1z+D_0w,
\qquad
\widetilde R=R+D_1z+D_0w+zw.
\]

For the second output the screen used

\[
r+z,\quad r+w,\quad r+z+w,\quad r+uz,\quad r+uw,
\]

and for the third output

\[
u,\ z,\ w,\ u+z,\ u+w,\ z+w,\ D_0+z,\ D_1+w.
\]

For each of these \(2\cdot5\cdot8=80\) exact coordinate skeletons, the
linear Jacobian-slice equation for a fourth block output has no solution of
total degree at most three over \(k(c,v)\).  Since \(c,v\) are retained, this
is exactly the six-variable determinant equation.  The calculation is over
the rational function field, so failure there also excludes polynomial
coefficients in \(c,v\).

This is a finite search result only.  It does not exclude arbitrary
degree-four outputs, different modifications, or replacing both \(c\) and
\(v\).

## 5. Affine transverse outputs are impossible

The affine-coupling system proposed below has a uniform first stratum which
can be closed exactly.  Put

\[
\begin{aligned}
\widetilde R&=R+D_1z+D_0w+zwH(r,u),\\
B&=r+zP(r,u)+wQ(r,u)+zwS(r,u),
\end{aligned}                                        \tag{1}
\]

where \(P,Q\) are arbitrary affine-linear functions of \(r,u\) over
\(K=k(c,v)\).  No degree restriction is needed on \(H,S\).  Suppose the
remaining two block outputs \(C,D\) are affine-linear in \(r,u,z,w\), again
with coefficients in \(K\).

On the slice \(z=w=0\), the first two gradient rows are

\[
\begin{aligned}
d\widetilde R&=(D_0D_1,\ R_u,\ D_1,\ D_0),\\
dB&=(1,\ 0,\ P,\ Q).                                  \tag{2}
\end{aligned}
\]

In particular \(H,S\) disappear.  Let

\[
p_{12},p_{13},p_{14},p_{23},p_{24},p_{34}
\]

be the six \(2\)-by-\(2\) minors of the constant gradient matrix
\((dC,dD)\).  They satisfy the Plücker relation

\[
p_{12}p_{34}-p_{13}p_{24}+p_{14}p_{23}=0.             \tag{3}
\]

Conversely, over the field \(K\), equation (3) is the complete
decomposability condition for an affine transverse two-plane.

Contracting \(d\widetilde R\wedge dB\) with these six Plücker coordinates
and setting the result equal to one gives eleven coefficient equations in
\(r,u\).  Together with (3), their exact Gröbner basis over
\(K=\mathbb Q(c,v)\) is

\[
\boxed{(1).}                                          \tag{4}
\]

Scaling one output converts any nonzero constant determinant to one.
Therefore:

\[
\boxed{\text{No map of the form (1) has two affine transverse outputs
and nonzero constant Jacobian.}}
\]

This conclusion allows arbitrary \(K\)-coefficients in \(P,Q,C,D\) and
arbitrary polynomial \(H,S\).  It is not a finite specialization argument.
It does not exclude a nonlinear \(C\) or \(D\).

There is also an exact first quadratic screen.  Keep \(P,Q\) arbitrary
affine and let \(D\) be a completely general polynomial of total block
degree at most two.  For

\[
C\in
\{u,z,w,u+z,u+w,z+w,D_0+z,D_1+w\},                   \tag{5}
\]

the coefficient ideal of the determinant-one equation restricted to
\(z=w=0\) is the unit ideal over \(\mathbb Q(c,v)\) in all eight cases.
This allows every coefficient of \(D\), rather than selecting a finite list
of quadratic outputs.  Notice that \(D_0+z\) is already a genuinely
quadratic third output.

The eight-row result is still a screen, not a classification: it does not
cover an arbitrary affine \(C\) outside (5), an arbitrary quadratic \(C\),
two simultaneously general quadratic outputs, or a fourth output of degree
at least three.

## 6. Surviving bounded moonshot

The next search should be confined to \(\mathbb A^6\) and impose all of the
following from the outset:

1. deform \(R\); retaining it is impossible;
2. replace \(u\) or \(v\);
3. do not retain \(r\) as an exposed output;
4. couple both modification variables before solving the determinant
   equation;
5. make at least one of the two transverse outputs genuinely nonlinear;
6. reject residual determinants whose restriction to the selected curve
   spans rank at most one, and reject any new zero such as \(r=3\);
7. after a Keller solution is found, test whether its target algebra
   recovers the source variables or is stably left-right equivalent to an
   existing \(\mathbb A^1/\mathbb G_m\) branch.

A concrete next coefficient system is

\[
\begin{aligned}
\widetilde R&=R+D_1z+D_0w+zwH(r,u),\\
B&=r+zP(r,u)+wQ(r,u)+zwS(r,u),
\end{aligned}
\]

with \(H,P,Q,S\) affine-linear at first, while the remaining two block
outputs have total degree at most four and at least one has degree at least
two.  The determinant is then imposed coefficient by coefficient over
\(k(c,v)\).  This is genuinely nonlinear in four block variables, uses two
interacting modification coordinates, and avoids every gate proved above.

One explicit noninvertible solution surviving the rank-two character and
stable-equivalence tests would give the desired mechanism beyond the
marked \(\mathbb A^1\) and \(\mathbb G_m\) branches.

## 7. Verification

Run

```bash
.venv/bin/python scripts/verify_three_puncture_nonlinear_frontier.py
```

The checker verifies both dimension-free common-zero gates, the two exact
five-variable transfer determinants, their restrictions to the selected
curve, all 80 degree-three six-variable slice systems, and the twelve
coefficient/Plücker equations giving the affine-transverse unit ideal (4).
It also checks the eight unit ideals in the general quadratic-fourth-output
screen (5).
