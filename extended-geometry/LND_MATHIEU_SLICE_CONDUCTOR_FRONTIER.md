# The LND conjecture for Mathieu--Zhao images:
# slice, finite fibers, and the conductor frontier

## 1. Status and scope

Let \(k\) be a field of characteristic zero, let
\(B=k[x_1,\ldots,x_n]\), let \(D\) be a locally nilpotent \(k\)-derivation
of \(B\), and let \(I\subseteq B\) be an ideal.  The strongest formulation
considered here is

\[
 D(I)\text{ is a Mathieu--Zhao subspace of }B.                 \tag{1.1}
\]

Thus, if \(f^m\in D(I)\) for every \(m\geq1\), then for every fixed
\(g\in B\),

\[
 gf^m\in D(I)\qquad(m\gg0).                                   \tag{1.2}
\]

There are two statements in the literature which must not be conflated.

1. The **LND conjecture** usually denotes the case \(I=B\), namely that
   \(\operatorname{Im}D\) is Mathieu--Zhao.
2. The ideal-image or **LNED** strengthening asks for \(D(I)\) for every
   ideal \(I\).  Principal-ideal results address only one part of this
   stronger statement.

As of July 2026, no counterexample to either characteristic-zero statement
was located.  The following positive results are the relevant boundary.

- Van den Essen, Wright, and Zhao proved the two-variable locally finite
  image theorem in
  [*Images of locally finite derivations of polynomial algebras in two variables*](https://doi.org/10.1016/j.jpaa.2010.12.002).
- Sun and Liu proved the rank-two and homogeneous rank-three cases in
  three variables in
  [*Images of locally nilpotent derivations of polynomial algebras in three variables*](https://doi.org/10.1016/j.jalgebra.2020.10.025).
- Sun and Wang proved the principal-ideal theorem for rank-one LNDs in
  every dimension and for a large class in three variables, including
  rank two and homogeneous rank three, in
  [*On the LND conjecture*](https://doi.org/10.1017/S000497272300059X).
- Liu, Sun, and Zeng proved the principal-ideal case in dimension two and
  further ideal classes in
  [*Images of locally nilpotent derivations acting on ideals of polynomial algebras*](https://doi.org/10.5486/PMD.2024.9435).
- Sun and Wang proved a relative two-variable theorem over a
  one-dimensional UFD, with conditional Dedekind-domain extensions, in
  [*Images of locally nilpotent derivations of bivariate polynomial algebras over a domain*](https://doi.org/10.21136/CMJ.2024.0008-24).
- Gupta, Ahuja, and Kour proved that the image of every rank-at-most-two
  linear derivation of a four-variable polynomial algebra is
  Mathieu--Zhao in
  [*Image of linear \(K\)-derivations and linear \(K\mathcal E\)-derivations of \(K[x_1,x_2,x_3,x_4]\)*](https://doi.org/10.1080/00927872.2023.2224885).
- Gupta proved a transfer result for primary ideals under linear
  derivations, together with monomial-ideal cases, in
  [*Image of ideals under linear K-derivations and the LNED conjecture*](https://doi.org/10.1016/j.jpaa.2025.108041).

The present note does not claim a counterexample or a proof of (1.1) in
full generality.  It records an exact slice reduction, closes every
carrier-free finite residual scheme, every monic carrier with finite
residual, and the normalized repeated nonmonic degree-drop chart.  It also
proves \(D(qB)\) Mathieu--Zhao for every nonzero principal carrier \(q\)
in the four-variable two-slice model
\(D=u\partial_x+v\partial_y\), as well as every ideal containing a power
of its plinth ideal \((u,v)\).  These are internal results and have not
been externally reviewed.

## 2. Slice normalization gives an exact membership test

Suppose first that \(B=A[s]\), \(D=\partial_s\), and
\(A=\ker D\).  Define the zero-constant primitive

\[
 {\cal J}\!\left(\sum_{j\geq0}a_js^j\right)
   =\sum_{j\geq0}\frac{a_j}{j+1}s^{j+1}.                     \tag{2.1}
\]

Then

\[
 \boxed{\quad
 h\in D(I)
 \ \Longleftrightarrow\
 {\cal J}(h)+c\in I\text{ for some }c\in A .
 \quad}                                                       \tag{2.2}
\]

Indeed, if \(h=D(F)\) with \(F\in I\), then
\(F-{\cal J}(h)\in\ker D=A\).  Conversely, differentiating
\({\cal J}(h)+c\in I\) gives \(h\in D(I)\).

For a local slice \(r\), put \(a=D(r)\in\ker D\) and
\(s=r/a\).  The slice theorem gives

\[
 B_a=(\ker D)_a[s],\qquad D=\partial_s
 \quad\text{on }B_a.                                         \tag{2.3}
\]

Equation (2.2) therefore solves membership away from the plinth divisor
\((a=0)\).  Returning from \(B_a\) to \(B\) is a saturation problem:
denominators and their valuations along the components of \(a=0\) must be
controlled uniformly in the power \(m\).

## 3. A safe nonprincipal class

The next theorem is elementary but useful because it excludes the most
obvious reduced finite-fiber counterexamples.

### Theorem 3.1 (radical finite fibers)

Let \(k\) be algebraically closed of characteristic zero, let
\(B=k[x,s]\), let \(D=\partial_s\), and let \(I\subset B\) be a radical
zero-dimensional ideal.  Then \(D(I)\) is a Mathieu--Zhao subspace.
More strongly, if \(f^m\in D(I)\) for every \(m\geq1\), then

\[
 gf^m\in D(I)
 \quad\text{for every }g\in B\text{ and every }m\geq1.        \tag{3.1}
\]

#### Proof

Write the finite zero set of \(I\) as

\[
 S=\{(a,b):a\in T,\ b\in S_a\},
\]

where \(T\subset k\) is finite and each \(S_a\subset k\) is finite.
By (2.2), \(h\in D(I)\) exactly when there is a polynomial \(c(x)\) such
that

\[
 {\cal J}(h)(a,b)+c(a)=0
 \qquad(a\in T,\ b\in S_a).                                  \tag{3.2}
\]

Polynomial interpolation chooses the values \(c(a)\) independently for
distinct \(a\).  Hence (3.2) is equivalent to

\[
 \int_{b_0}^{b}h(a,t)\,dt=0
 \quad
 (a\in T,\ b_0,b\in S_a).                                    \tag{3.3}
\]

Assume \(f^m\in D(I)\) for every positive \(m\).  Fix a fiber with two
distinct points \(b_0,b\).  Then

\[
 \int_{b_0}^{b}f(a,t)^m\,dt=0\qquad(m\geq1).                  \tag{3.4}
\]

The one-variable polynomial moment lemma says that (3.4) forces
\(f(a,t)=0\) as a polynomial in \(t\).  This is Lemma 2.2 in the 2024
Sun--Wang paper cited above.  Thus \(f\) vanishes identically on every
vertical fiber containing at least two points.  For any \(g\) and any
\(m\geq1\), all the integrands \(g(a,t)f(a,t)^m\) on those fibers are
zero.  The conditions (3.3) hold for \(gf^m\), proving (3.1). \(\square\)

The theorem is a statement about a nonprincipal class: a radical finite
set with a repeated \(x\)-coordinate generally has a nonprincipal ideal.
It is not merely the principal-ideal theorem in different coordinates.

### Theorem 3.2 (arbitrary finite residual scheme)

Let \(k\) be algebraically closed of characteristic zero, let
\(B=k[x,s]\), let \(D=\partial_s\), and let \(J\subset B\) be any
zero-dimensional ideal.  Then \(D(J)\) is a Mathieu--Zhao subspace.

#### Proof

First consider a primary ideal \(Q\), supported at
\(P=(a,b)\), and choose \(N\) with
\(\mathfrak m_P^N\subseteq Q\).  There are two cases.

If some \(F\in Q\) has \(\partial_sF(P)\ne0\), the formal implicit-function
argument modulo \(\mathfrak m_P^N\) gives a polynomial
\(\phi(x)\) such that

\[
 s-\phi(x)\in Q.
\]

Indeed, modulo \(\mathfrak m_P^N\), \(F\) is a unit times
\(s-\phi(x)\); that unit remains invertible modulo the primary ideal
\(Q\).  Hence \(1=D(s-\phi(x))\in D(Q)\), so \(D(Q)=B\).

Otherwise every element of \(Q\) has zero \(s\)-derivative at \(P\).
If \(f^m\in D(Q)\) for every \(m\), choose \(F_m\in Q\) with
\(D(F_m)=f^m\).  Evaluation at \(P\) gives

\[
 f(P)^m=D(F_m)(P)=0,
\]

and therefore \(f\in\mathfrak m_P\).  For fixed \(g\), put

\[
 H_m(x,s)=\int_b^s g(x,t)f(x,t)^m\,dt.
\]

Its local order at \(P\) tends to infinity with \(m\), so
\(H_m\in\mathfrak m_P^N\subseteq Q\) for \(m\gg0\).  Thus
\(gf^m=D(H_m)\in D(Q)\).  This proves the primary case.

Now write

\[
 J=\bigcap_{a\in T}J_a,
\]

where \(J_a\) is the intersection of the primary components supported on
the vertical fiber \(x=a\).  Primitive constants for distinct \(a\)'s can
be patched by Hermite interpolation, because each \(J_a\) contains a
power of \(x-a\).  Consequently membership in \(D(J)\) is equivalent to
simultaneous membership in all \(D(J_a)\).

If \(J_a\) has one support point, the primary argument applies.  If it has
at least two support points, reduction modulo \(\sqrt{J_a}\) and Theorem
3.1 show that \(f(a,s)=0\) identically.  Hence \(x-a\mid f\), and the
zero-constant primitive of \(gf^m\) lies in \(J_a\) for \(m\gg0\), since
its \((x-a)\)-adic order grows linearly.  Hermite interpolation patches
the finitely many singleton-fiber primitive constants.  Therefore
\(gf^m\in D(J)\) for all sufficiently large \(m\). \(\square\)

## 4. Principal carriers and finite residual schemes

Every nonzero ideal of the UFD \(k[x,s]\) can be written

\[
 I=qJ,                                                        \tag{4.1}
\]

where \(q\) is a greatest common divisor of its generators and either
\(J=B\) or \(J\) has height two.  In the latter case \(B/J\) has finite
length.  Thus the two endpoints are:

- \(J=B\): the principal carrier \(qB\), covered by the known
  principal-ideal theorem;
- \(q=1\) and \(J\) radical: covered by Theorem 3.1.

The next theorem closes the entire monic positive-\(s\)-degree carrier
interface, including arbitrary nonreduced residual schemes.

### Theorem 4.1 (monic carrier with finite residual)

Let \(k\) be algebraically closed of characteristic zero, let
\(B=k[x,s]\), let \(D=\partial_s\), and let

\[
 I=qJ,
\]

where \(q\in k[x,s]\) is monic of positive \(s\)-degree and \(J=B\) or
\(B/J\) has finite length.  Then \(D(I)\) is a Mathieu--Zhao subspace.

#### Proof

First suppose that \(q\), considered over \(\overline{k(x)}\), has two
distinct roots \(\alpha,\beta\).  If \(f^m\in D(qJ)\), then
\(f^m\in D(qB)\), and primitives evaluated at the two roots give

\[
 \int_\alpha^\beta f(x,t)^m\,dt=0\qquad(m\geq1).              \tag{4.2}
\]

The one-variable moment lemma over \(\overline{k(x)}\) forces \(f=0\).
The Mathieu conclusion is then immediate.

It remains to treat the case where \(q\) has one distinct root.  Since
\(q\) is monic and the characteristic is zero,

\[
 q=(s-\alpha(x))^d
\]

for some \(\alpha(x)\in k[x]\).  Translation by \(\alpha(x)\) preserves
\(\partial_s\), so write \(t=s-\alpha(x)\) and assume \(q=t^d\).  Let
\({\cal J}_0\) denote the \(t\)-primitive vanishing at \(t=0\).  Since every
element of \(t^dJ\) vanishes at \(t=0\),

\[
 h\in D(t^dJ)
 \Longleftrightarrow
 {\cal J}_0(h)=t^d A\text{ for some }A\in J.                  \tag{4.3}
\]

Assume \(f^m\in D(t^dJ)\) for every \(m\).  From the case \(m=1\),

\[
 f=D(t^dA_1)
   =t^{d-1}(dA_1+t\partial_tA_1),                             \tag{4.4}
\]

so \(t^{d-1}\mid f\).  Consequently
\({\cal J}_0(gf^m)/t^d\) is a polynomial for every \(g\) and \(m\geq1\).

Take a primary component \(Q\) of \(J\), supported at
\((x,t)=(a,b)\), and choose \(N\) with
\((x-a,t-b)^N\subseteq Q\).

If \(b\ne0\), evaluation of
\({\cal J}_0(f^m)/t^d\in J\) at \((a,b)\) gives

\[
 b^{-d}\int_0^b f(a,t)^m\,dt=0\qquad(m\geq1).
\]

Hence \(f(a,t)=0\), so \(x-a\mid f\).  The quotient
\({\cal J}_0(gf^m)/t^d\) is divisible by \((x-a)^m\) locally at \(b\ne0\),
and therefore belongs to \(Q\) once \(m\geq N\).

If \(b=0\) and \(d\geq2\), (4.4) shows that the same quotient is divisible
by

\[
 t^{m(d-1)+1-d},
\]

whose order tends to infinity.  If \(b=0\) and \(d=1\), evaluation at
\((a,0)\) gives \(f(a,0)^m=0\), hence
\(f\in(x-a,t)\); integration and division by \(t\) preserve the resulting
linearly growing local order.  In either case
\({\cal J}_0(gf^m)/t^d\in Q\) for all sufficiently large \(m\).

There are only finitely many primary components.  Thus
\({\cal J}_0(gf^m)/t^d\in J\) for \(m\gg0\), and (4.3) gives
\(gf^m\in D(I)\). \(\square\)

### Theorem 4.2 (normalized repeated degree-drop carrier)

Let

\[
 p=xs-1,\qquad I=x^c p^dJ,\qquad c\geq0,\ d\geq1,
\]

where \(J=B\) or \(B/J\) has finite length and is supported at
\((0,0)\).  Then \(D(I)\) is a Mathieu--Zhao subspace.

#### Proof

For \(d=1\), write \(h=\sum_{j=0}^n h_j(x)s^j\).  Solving

\[
 h=D(pA)=xA+p\partial_sA
\]

from the highest \(s\)-coefficient down shows that

\[
 h\in D(pB)
 \Longleftrightarrow
 T(h):=\sum_{j=0}^n\frac{h_j(x)}{(j+1)x^{j+1}}\in k[x].       \tag{4.5}
\]

The corresponding quotient \(A=\sum A_js^j\) satisfies

\[
 A_j=\frac{h_j}{(j+1)x}+\frac{A_{j+1}}x,\qquad A_{n+1}=0.    \tag{4.6}
\]

For nonzero \(f=\sum f_j(x)s^j\), define

\[
 w(f)=\min_j\bigl(v_x(f_j)-j\bigr)
\]

and let \(F(u)\ne0\) be its initial polynomial after substituting
\(s=u/x\).  Then

\[
 T(f^m)
 =x^{mw(f)-1}
   \left(\int_0^1F(u)^m\,du+O(x)\right).                    \tag{4.7}
\]

For general \(c,d\), the unique primitive vanishing at the rational root
\(s=1/x\) gives the quotient

\[
 A_m=x^{-c}p^{-d}\int_{1/x}^s f(x,t)^m\,dt.                 \tag{4.8}
\]

If \(f^m\in D(x^cp^dJ)\) for all \(m\), then \(A_m\) is a polynomial for all
\(m\).  At \(s=0\), substitution \(u=xt\) gives

\[
 A_m(x,0)
 =(-1)^{d+1}x^{mw(f)-1-c}
   \left(\int_0^1F(u)^m\,du+O(x)\right).                    \tag{4.9}
\]

Were \(w(f)\leq0\), (4.9) would force
\(\int_0^1F(u)^m\,du=0\) for every positive \(m\).  The one-variable
moment lemma would give \(F=0\), a contradiction.  Thus

\[
 w(f)\geq1.                                                  \tag{4.10}
\]

Moreover, the order of (4.8) at \(p=0\) and the case \(m=1\) give

\[
 p^{d-1}\mid f.                                              \tag{4.11}
\]

For fixed \(g\), the \(p\)-order of \(gf^m\) is therefore at least
\(m(d-1)\), while

\[
 w(gf^m)=w(g)+mw(f)\longrightarrow\infty.
\]

Thus, for \(m\gg0\), the primitive of \(gf^m\) has \(p\)-order at least
\(d\), so its quotient by \(x^cp^d\) is polynomial.  Substitution \(u=xs\)
and expansion back in \(s\) show that every coefficient of this quotient
is divisible by an \(x\)-power tending to infinity; the fixed loss \(c\)
does not affect growth.  For \(c=0,d=1\), this is also immediate from
(4.6).  Since \(x^N B\subseteq J\) for some \(N\), the quotient lies in
\(J\) for \(m\gg0\).  Hence \(gf^m\in D(x^cp^dJ)\). \(\square\)

### 4.3 Where a slice counterexample can still hide

The first unresolved interfaces are now:

\[
 \deg_s q>0\text{ with }q\text{ nonmonic},                   \tag{4.12}
\]

excluding the normalized linear chart of Theorem 4.2.  For a remaining
nonmonic carrier, its leading \(s\)-coefficient vanishes on a divisor of
the \(x\)-line.  Repeated roots then escape to infinity or collide with
the content divisor, and the monic proof works only after localization.
Returning across that divisor is a valuation/conductor saturation problem.

This leads to a sharper necessary profile for a two-variable slice
counterexample:

1. it cannot be a monomial ideal, since \(\partial_s(I)\) is then a
   monomial ideal;
2. it cannot have \(q=1\) with any zero-dimensional residual scheme,
   by Theorem 3.2;
3. it cannot have a monic carrier of positive \(s\)-degree with any finite
   residual scheme, by Theorem 4.1;
4. it cannot have carrier \(q=x^c(xs-1)^d\) with the residual scheme on
   the degree-drop point, by Theorem 4.2;
5. its obstruction must survive every power through a leading-coefficient
   valuation face or a nonreduced jet functional;
6. a bounded zero prefix is not enough: both pure membership and mixed
   nonmembership require all-order certificates.

## 5. How the five techniques fit together

### Slice normalization

Use (2.2) on a global slice and (2.3) on a local slice.  Record the plinth
element \(a=D(r)\), not just the rational coordinate \(s=r/a\).

### Support cones

Give \(a\), \(q\), the residual generators, \(f\), and \(g\) a common
multifiltration.  If every exposed face of \(f\) lies strictly on one side
of the primitive obstruction, multiplication by a fixed \(g\) cannot cross
the linearly growing gap.  Such a chart is safe and should be discarded
before elimination.

### Valuation faces

For every height-one prime \(p\mid aq\), take the \(p\)-initial face of

\[
 {\cal J}(f^m)+c_m\in qJ.                                    \tag{5.1}
\]

A counterexample must retain at least two tied faces for infinitely many
\(m\); a unique leading face gives either a contradiction or eventual
mixed membership.

### Conductor saturation

Normalize the finite residual algebra and compute its conductor.  Conditions
visible on separate normalization branches are reduced interval moments and
are safe by Theorem 3.1.  The only remaining information is conductor gluing
and nilpotent jet data.  Saturation by \(a q\) must be performed before
interpreting a component as a candidate.

### Factorial and constant-term separation

After a valuation face is fixed, encode its primitive obstruction as a
finite sum

\[
 L_m=\sum_{\nu} c_\nu\,\lambda_\nu^m R_\nu(m),                \tag{5.2}
\]

where \(R_\nu(m)\) is typically a polynomial, factorial ratio, or constant
term.  Distinct exponential bases and distinct extremal factorial slopes
separate.  A viable pure cancellation needs a tied group in (5.2), while a
counterexample needs the corresponding mixed sequence to remain nonzero
for infinitely many \(m\).

## 6. Concrete search order

The first search should stay with \(D=\partial_s\) and use

\[
 I=qJ
\]

with \(q\) nonmonic of small bidegree and \(J\) a primary ideal of length
two through six supported over the vanishing of the leading
\(s\)-coefficient of \(q\).  Carrier-free, monic-carrier, and normalized
linear \(q=x^c(xs-1)^d\) charts should be rejected before any moment
elimination.  The first surviving carrier must therefore have a more
complicated leading-coefficient divisor or conductor gluing not reducible
to that normalized one-branch model.  For each remaining chart:

1. compute the exact primitive membership functionals from (2.2);
2. split them on the normalization of \(B/J\);
3. remove support cones with a strict separator;
4. saturate by \(q\) and the conductor;
5. solve the first moments only as candidate generation;
6. promote a candidate only after deriving an all-order recurrence or
   constant-term identity.

If this slice search is exhausted, the next target is a rank-three,
nonhomogeneous LND on \(k[x,y,z]\) whose local slice has a reducible plinth
divisor.  Rank-one, rank-two in three variables, and homogeneous rank-three
principal-ideal charts are already inside the proved literature boundary.

The first plinth audit uses

\[
 D=x\partial_y+y\partial_z,\qquad
 \ker D\supset k[x,2xz-y^2],\qquad s_{\rm loc}=y/x.
\]

The exact homogeneous compiler tests five ideals crossing \(x=0\), 45
sparse seeds, and pure powers through six.  Thirteen seeds survive their
pure prefixes, all in the visible \(x\)-support cone, and none has a mixed
tail obstruction for the multipliers \(1,x,y,z\).  This is a bounded
negative search, not an all-order theorem.  The next plinth pass should use
the concrete reducible-plinth LND

\[
 D=x(x-1)\partial_y+y\partial_z,\qquad
 2x(x-1)z-y^2\in\ker D,
\]

with ideals coupling the fibers \(x=0\) and \(x=1\).  The exact
\(\mathrm{wt}(x,y,z)=(0,1,2)\) module search tests five such ideals and
205 branch-aware seeds through the sixth pure power.  Fourteen seeds
survive; every one is \(x(x-1)y^r\) for \(1\leq r\leq3\), and none has a
mixed-tail obstruction for \(1,x,x-1,y,z\).  Thus both plinth valuations
again grow together in this grid.

The nonlinear triangular follow-up

\[
 D=x(x-1)\partial_y+y^2\partial_z,\qquad
 3x(x-1)z-y^3\in\ker D,
\]

has the exact weight grading \((0,1,3)\).  Five coupled ideals and 210
branch-aware seeds through weight four leave seventeen six-power
survivors, again all \(x(x-1)y^r\), now with \(1\leq r\leq4\), and no
mixed-tail obstruction.  Thus the homogeneous triangular family still
forces both branch valuations to grow.

The first branch-asymmetric perturbation is

\[
 D=x(x-1)\partial_y+(y^2+x)\partial_z,\qquad
 3x(x-1)z-y^3-3xy\in\ker D.
\]

Write this as \(D=D_2+x\partial_z\), where
\(D_2=x(x-1)\partial_y+y^2\partial_z\).  The two summands commute.
Although the \((0,1,3)\)-grading is broken, the leading component of every
kernel element is a polynomial in \(3x(x-1)z-y^3\).  Subtracting powers of
the full invariant therefore normalizes any primitive of a target of
weight \(H\) to weight at most \(H+1\).  This gives an exact two-stage
membership test:

1. lift one normalized primitive by a finite \(\mathbf Q[x]\)-module
   standard-basis calculation;
2. reduce it modulo \(I\) and test whether its class lies in the finite
   image of \(k[x,3x(x-1)z-y^3-3xy]\).

The compiler applies this test to eight ideals, including three
non-weight-homogeneous jet charts, and 350 mixed-weight seeds.  There are
48 chartwise six-power survivors but only seven distinct forms; every one
is divisible by \(x(x-1)\), and none has a mixed-tail obstruction.

The next geometry makes the plinth components meet rather than placing them
at two disjoint points of the invariant \(x\)-line.  On \(k[u,v,y,z]\) take

\[
 D=uv\partial_y+(y^2+u)\partial_z,\qquad
 3uvz-y^3-3uy\in\ker D.
\]

The same leading-weight normalization works over \(\mathbf Q[u,v]\).
For each target, a Singular module standard basis constructs one exact
bounded primitive; its class in a zero-dimensional quotient is then tested
against the finite image of
\(\mathbf Q[u,v,3uvz-y^3-3uy]\).  Five crossing ideals and 956 seeds,
including quadratic base coefficients, invariant combinations, and images
of ideal generators, leave 177 chartwise six-power survivors representing
36 distinct forms.  Every survivor is divisible by the full plinth
equation \(uv\), and no multiplier in \(1,u,v,y,z\) gives a mixed-tail
obstruction through exponents four to six.  More sharply, every tested
seed outside the \(uv\)-support cone fails by its square: the longest
nonplinth pure prefix has length one.

Thus merely intersecting the two components does not defeat the
support-cone mechanism in a finite residual quotient.  The next
counterexample chart should test a genuinely nonprincipal plinth ideal,
where two slice denominators must be glued rather than multiplied.  The
minimal concrete target is

\[
 D=u\partial_x+v\partial_y
 \quad\hbox{on }k[u,v,x,y],\qquad
 \ker D=k[u,v,uy-vx],
\]

whose plinth ideal contains \((u,v)\).  Ideals meeting \(u=v=0\) retain a
two-chart Čech obstruction that the single equation \(uv\) cannot see.

The exact two-slice compiler again works weight by weight over
\(\mathbf Q[u,v]\), now normalizing primitives modulo
\(\mathbf Q[u,v,uy-vx]\).  Five zero-dimensional plinth ideals and 1,055
seeds leave 494 chartwise six-power survivors, representing 100 distinct
forms.  All 100 lie in the base-support cone \((u,v)A\), none gives a
mixed-tail obstruction for \(1,u,v,x,y,uy-vx\), and every tested seed
outside \((u,v)A\) already fails its first power.  Thus nonprincipality
alone is still insufficient when the residual quotient is finite.

The positive-dimensional continuation admits an all-order support
saturation, not merely another bounded search.

### Theorem 6.1 (nonprincipal plinth-power saturation)

Let

\[
 A=k[u,v,x,y],\qquad D=u\partial_x+v\partial_y,\qquad
 \mathfrak p=(u,v),
\]

over a characteristic-zero field.  If an ideal \(I\subseteq A\) contains
\(\mathfrak p^q\) for some \(q\ge1\), then \(D(I)\) is a Mathieu--Zhao
subspace.  Equivalently, the conclusion holds whenever
\(\mathfrak p\subseteq\sqrt I\).

**Proof.**
Grade \(A\) by total \((u,v)\)-degree:

\[
 A=\bigoplus_{r\ge0}A_r,\qquad D(A_r)\subseteq A_{r+1}.
\]

For every \(N\ge1\),

\[
 \boxed{\quad
 \operatorname {Im}D\cap\mathfrak p^N
 =D(\mathfrak p^{N-1}).
 \quad}                                                       \tag{6.1}
\]

The right-to-left inclusion is immediate.  For the converse, write a
primitive \(P=\sum_rP_r\), with \(P_r\in A_r\), of
\(h\in\operatorname {Im}D\cap\mathfrak p^N\).  Its components of
\((u,v)\)-degree below \(N\) vanish.  Hence \(D(P_r)=0\) for
\(r<N-1\).  Subtracting these kernel components leaves a primitive in
\(\mathfrak p^{N-1}\), proving (6.1).

Now suppose \(f^m\in D(I)\) for all sufficiently large \(m\), and fix
\(g\in A\).  The derivation \(D\) is linear of rank two, so the theorem of
Gupta--Ahuja--Kour cited in Section 1 says that \(\operatorname {Im}D\)
is Mathieu--Zhao.  Therefore \(gf^m\in\operatorname {Im}D\) for all large
\(m\).  Since
\(\operatorname {Im}D\subseteq\mathfrak p\), primality of
\(\mathfrak p\) gives \(f\in\mathfrak p\).  Thus
\(gf^m\in\mathfrak p^m\), and for \(m\ge q+1\), (6.1) gives

\[
 gf^m\in D(\mathfrak p^{m-1})
 \subseteq D(\mathfrak p^q)
 \subseteq D(I).
\]

This proves the theorem. \(\square\)

For the first free-line ideal \(I=(u,v,x)\), one can see the saturation
through a single explicit conductor functional.  Put \(w=uy-vx\) and

\[
 \Lambda(h)=[v^1]\,h(0,v,0,y)\in k[y].
\]

Since \(\ker D=k[u,v,w]\) maps to the constants in \(A/I=k[y]\), primitive
normalization gives

\[
 D(I)=\operatorname {Im}D\cap\ker\Lambda.                    \tag{6.2}
\]

Indeed, if \(D(P)=h\), then
\(\Lambda(h)=\partial_yP(0,0,0,y)\).  Theorem 6.1 proves this ideal safe,
as well as \((u,v,\ell(x,y))\) for every nonzero linear form \(\ell\), and
all nilpotent thickenings containing a power of \((u,v)\).

The positive-dimensional compiler also tests the four nilpotent variants

\[
 (u^2,uv,v^2,x),\quad (u^2,v^2,x),\quad
 (u^2,uv,v^2,x-uy),\quad (u^2,v^2,x-uy).
\]

Together with \(I\), these five charts and 1,055 seeds leave 289 chartwise
six-power survivors representing 75 forms.  Every survivor lies in
\((u,v)A\), every seed outside that cone fails its first power, and there
is no mixed-tail obstruction.  Membership is exact without truncating the
free \(y\)-direction because the images of \(u,v,w\) are nilpotent in each
quotient.  The exponent window is bounded, but Theorem 6.1 independently
proves all five charts safe to every order.

The next live counterexample locus cannot be set-theoretically supported
on the plinth locus.  The minimal target is the principal ideal \(I=(x)\).
Here

\[
 A/I=k[u,v,y],\qquad
 \operatorname {im}(\ker D\to A/I)=k[u,v,uy],
\]

which has positive module rank and a degree-drop conductor along \(u=0\).
The grading saturation no longer places a primitive back in \(I\).

Nevertheless its kernel correction has a sharp semigroup description.
For a primitive \(P\) of \(h\),

\[
 h\in D((x))
 \Longleftrightarrow
 P\bmod x\in k[u,v,uy].
\]

Equivalently, every residue monomial \(u^av^by^c\) must satisfy \(a\ge c\).
The exact valuation-face compiler tests 1,055 seeds through six pure
powers.  It leaves 48 survivors, all divisible by \(u\); none gives a
mixed-tail obstruction, and every failed seed has pure-prefix length at
most one.  Thus the first principal-conductor grid is again negative, but
only boundedly.

The strict faces admit a general all-order statement.  For a nonzero
polynomial \(F\), write \(\nu_u(F)\) and \(\nu_v(F)\) for its \(u\)- and
\(v\)-adic orders, and put

\[
 \omega(F)=\min\{a-c:u^av^bx^dy^c\in\operatorname {supp}F\}.
\]

### Theorem 6.2 (principal-conductor strict cones)

For \(D=u\partial_x+v\partial_y\) and \(I=(x)\), fix \(f\in A\).
Then \(gf^m\in D(I)\) for every \(g\in A\) and all sufficiently large
\(m\) if any one of the following holds:

1. \(\nu_u(f)>\deg_y f\);
2. \(\nu_v(f)>\deg_x f\) and \(\omega(f)>0\);
3. \(f=ub\) for some \(b\in\ker D=k[u,v,uy-vx]\).

**Proof.**
Over \(A_u\), let

\[
 J_x(h)=u^{-1}\int_0^x h(t,y)\,dt,\qquad
 S=J_xv\partial_y.
\]

Since \(S\) lowers \(y\)-degree,

\[
 P=J_x\sum_{j\ge0}(-S)^j(gf^m)
\]

is a finite sum and \(D(P)=gf^m\).  A term with \(j\) applications of
\(S\) loses \(j+1\) powers of \(u\).  The first inequality makes every
term polynomial for large \(m\), and every term is divisible by \(x\).
This proves case 1.

Similarly, over \(A_v\), put

\[
 J_y(h)=v^{-1}\int_0^y h(x,t)\,dt,\qquad
 T=J_yu\partial_x.
\]

Now \(T\) lowers \(x\)-degree and preserves \(\omega\), while the outer
\(J_y\) lowers \(\omega\) by one.  The second pair of inequalities makes

\[
 P=J_y\sum_{j\ge0}(-T)^j(gf^m)
\]

polynomial and gives \(\omega(P)\ge0\) for large \(m\).  Hence
\(P\bmod x\in k[u,v,uy]\); subtracting the corresponding invariant makes
the primitive divisible by \(x\).  This proves case 2.

Finally, if \(f=ub\) with \(D(b)=0\), the local-slice primitive

\[
 \sum_{j\ge0}
 \frac{(-1)^j}{(j+1)!}\,
 f^mD^j(g)\left(\frac{x}{u}\right)^{j+1}
\]

is finite, differentiates to \(gf^m\), and is polynomial and divisible by
\(x\) once \(m\) exceeds the nilpotence order of \(g\).  This proves case
3. \(\square\)

All three cases sit inside a simpler local-slice polynomial cone.

### Theorem 6.3 (principal-conductor local-slice polynomial cone)

Put \(B=\ker D=k[u,v,w]\), \(w=uy-vx\).  If

\[
 f=u^rF,\qquad r\ge1,\qquad F\in B[x],
\]

then \(gf^m\in D((x))\) for every \(g\in A\) and all sufficiently large
\(m\).

**Proof.**
In \(A_u=B_u[s]\), \(s=x/u\) and \(D=\partial_s\).  Write

\[
 g=\sum_{j=0}^d b_js^j,\qquad b_j\in B_u,
\]

and choose \(e\) such that \(u^eb_j\in B\) for every \(j\).  If
\(F=\sum_\ell c_\ell x^\ell\), then

\[
 F^m=\sum_k c_{m,k}u^ks^k,\qquad c_{m,k}\in B.
\]

Termwise zero-constant integration gives the primitive

\[
 P_m=\sum_{j,k}
 \frac{b_jc_{m,k}}{j+k+1}\,
 u^{rm+k}s^{j+k+1}.
\]

After substituting \(s=x/u\), every summand contains

\[
 b_j\,u^{rm-j-1}x^{j+k+1}.
\]

For \(m\) large relative to \(d,e\), this is polynomial and divisible by
\(x\).  Hence \(P_m\in(x)\) and \(D(P_m)=gf^m\). \(\square\)

The exact survivor census verifies that all 48 forms lie in \(uB[x]\).
For the five \(y\)-dependent forms, division by \(u\) gives respectively

\[
 w,\quad w+2vx,\quad w-uv+2vx,\quad u(w+2vx),\quad
 2vw+2v^2x.
\]

Thus every survivor is now proved safe to all orders.  The strict-cone
Theorem 6.2 remains useful because its hypotheses can be checked directly
in the original coordinates, while Theorem 6.3 gives the larger invariant
normal form.

The finite census suggested a square gate, but that gate is false.  Over a
characteristic-zero extension containing \(\alpha=\sqrt{-15}/3\), put

\[
 p(t)=2t-1+\alpha(6t^2-6t+1).
\]

Orthogonality of the first two shifted Legendre polynomials gives

\[
 \int_0^1p(t)\,dt=\int_0^1p(t)^2\,dt=0,
 \qquad
 \int_0^1p(t)^3\,dt=\frac{32\sqrt{-15}}{315}\ne0.             \tag{6.3}
\]

Homogenize this face as

\[
 F(z,w)=w(2z-w)+\alpha(6z^2-6wz+w^2)
\]

and set \(z=uy\), \(w=uy-vx\).  Then

\[
 f=F(uy,uy-vx)
\]

satisfies \(f,f^2\in D((x))\) but \(f^3\notin D((x))\).  Thus two powers
do not force \(f\in uB[x]\).  This is not an LNED counterexample: the
third power already fails.  The exact checker verifies all three
statements symbolically.

The full eventual-power hypothesis nevertheless closes the gap.

### Theorem 6.4 (moving-linear principal ideal)

Let

\[
 A=k[u,v,x,y],\qquad D=u\partial_x+v\partial_y
\]

over a characteristic-zero field.  Then \(D((x))\) is a Mathieu--Zhao
subspace.  In fact, if \(f^m\in D((x))\) for all sufficiently large \(m\),
then, for every fixed \(g\in A\),

\[
 gf^m\in D((x))\qquad(m\gg0).
\]

After a simultaneous linear change of the pairs \((u,v)\) and \((x,y)\),
the same conclusion holds with \((x)\) replaced by
\((\ell(x,y))\) for every nonzero linear form \(\ell\).

**Proof.**
Put \(B=k[u,v,w]\), where \(w=uy-vx\), and embed

\[
 A\hookrightarrow B_u[x],\qquad y=\frac{w+vx}{u}.
\]

On this localization \(D=u\partial_x\), with \(u,v,w\) held fixed.  For
\(h\in A\), define

\[
 {\cal T}(h)=\frac1x\int_0^x h(u,v,t,(w+vt)/u)\,dt.
\]

The unique primitive of \(h\) vanishing at \(x=0\) is
\(x{\cal T}(h)/u\).  Any element of \((x)\) also vanishes there, so
invariant constants cannot alter this primitive.  Consequently

\[
 \boxed{\quad
 h\in D((x))\Longleftrightarrow {\cal T}(h)\in uA.
 \quad}                                                       \tag{6.4}
\]

Write the finite Laurent expansion of a nonzero \(f\) as

\[
 f=u^rF(v,w,x)+\text{terms of higher \(u\)-degree},
 \qquad F\ne0.
\]

If \(r\ge1\), then \(f\in uB[x]\), and Theorem 6.3 gives the asserted
mixed-power conclusion.  Suppose instead that \(r=-q\le0\).
Writing \(z=w+vx\), the lowest-face description of
\(A=k[u,v,x,z/u]\) gives

\[
 z^q\mid F.                                                   \tag{6.5}
\]

We use the same description once more.  If \(H\in uA\), then its
coefficient of \(u^R\), for every \(R\le0\), is divisible by
\(z^{1-R}\).  Indeed, a monomial of \(uA\) having \(u\)-degree \(R\)
comes from \(u^{a+1}y^c\) with \(R=a+1-c\), and therefore contains
\(z^c\) with \(c=a+1-R\ge1-R\).

Assume now that \(f^m\in D((x))\) for every sufficiently large \(m\).
The coefficient of \(u^{-qm}\) in \({\cal T}(f^m)\) is
\({\cal T}(F^m)\).  By (6.4) and the preceding face divisibility,

\[
 z^{qm+1}\mid{\cal T}(F^m)                                   \tag{6.6}
\]

for every such \(m\), with the assertion automatic if that coefficient
vanishes.

Pass to \(K=k(v,w)\) and change from \(x\) to \(z=w+vx\).  Define

\[
 \Phi(t)=F\!\left(v,w,\frac{t-w}{v}\right)\in K[t].
\]

Equation (6.5) says \(t^q\mid\Phi(t)\), while direct substitution gives

\[
 {\cal T}(F^m)
 =\frac1{z-w}\int_w^z\Phi(t)^m\,dt.                           \tag{6.7}
\]

At \(z=0\), the denominator \(z-w\) is a unit and
\(\int_0^z\Phi(t)^m\,dt\) has order at least \(qm+1\).
Therefore (6.6) is equivalent to

\[
 \int_0^w\Phi(t)^m\,dt=0.                                    \tag{6.8}
\]

This holds for every sufficiently large \(m\).  Choose one such \(M\).
Then (6.8) for the exponents \(Mn\), \(n\ge1\), says that every positive
power of the polynomial \(\Phi^M\) has zero integral from \(0\) to \(w\).
The one-variable polynomial moment lemma used in Theorems 3.1 and 4.1,
applied after extending \(K\) algebraically if necessary, forces
\(\Phi^M=0\).  This contradicts \(F\ne0\).

Hence \(r\ge1\), and Theorem 6.3 completes the proof for \((x)\).  Given
\(\ell=ax+by\ne0\), extend \((a,b)\) to a matrix
\(M\in\operatorname {GL}_2(k)\), set

\[
 \binom{x'}{y'}=M\binom{x}{y},\qquad
 \binom{u'}{v'}=M\binom{u}{v},
\]

and observe that \(D=u'\partial_{x'}+v'\partial_{y'}\) and
\((\ell)=(x')\).  The first case therefore applies. \(\square\)

The square-gate failure explains why the sparse rational grid saw no
candidate: cancellation of the first two moments occurs only on an
algebraic isotropic face, and it still breaks at the next moment.  The
all-moment argument shows that no finite prefix can replace the eventual
hypothesis.

The moving-linear theorem bootstraps through every multiplicity.

### Corollary 6.5 (moving-linear thickenings)

For every nonzero linear form \(\ell(x,y)\) and every \(d\ge1\),

\[
 D((\ell^d))
\]

is Mathieu--Zhao.

**Proof.**
It is enough to take \(\ell=x\).  The case \(d=1\) is Theorem 6.4.
For \(d\ge2\),

\[
 D(x^dA)\subseteq x^{d-1}A.
\]

Thus eventual membership \(f^m\in D((x^d))\) forces \(x\mid f\).
Theorem 6.4 gives, for every fixed \(g\) and all large \(m\), a primitive

\[
 D(P_m)=gf^m,\qquad P_m\in xA.
\]

Also \(gf^m\in x^mA\).  If \(P_m=x^ra\), \(1\le r<d\), and
\(D(P_m)\in x^rA\), then

\[
 D(P_m)=ru\,x^{r-1}a+x^rD(a).
\]

Reduction modulo \(x^r\) gives \(x\mid ru a\).  Since
\(A/(x)=k[u,v,y]\) is a domain and the characteristic is zero, \(x\mid a\).
Starting at \(r=1\) and taking \(m\ge d\), iteration gives
\(P_m\in x^dA\).  Hence \(gf^m\in D((x^d))\). \(\square\)

Multiplicity is therefore not the next obstruction.  Nor is the first
reducible carrier: two distinct generic orbit intersections force the
Mathieu radical to vanish.

### Theorem 6.6 (generic two-root carrier)

Put \(B=\ker D=k[u,v,w]\), \(w=uy-vx\), and
\(L=\operatorname {Frac}(B)\).  Let \(0\ne q\in A\).  If, as a polynomial
in the local-slice coordinate \(x\) over \(\overline L\), the carrier
\(q\) has at least two distinct roots, then

\[
 \{f\in A:f^m\in D(qA)\text{ for all }m\gg0\}=\{0\}.
\]

In particular, \(D((xy))\) and \(D((x(x-1)))\) are Mathieu--Zhao.

**Proof.**
In \(L[x]\), \(D=u\partial_x\).  Suppose \(f^m\in D(qA)\) for every
sufficiently large \(m\), and choose primitives

\[
 P_m=qa_m,\qquad D(P_m)=f^m.
\]

Let \(\alpha\ne\beta\) be two roots of \(q\) in \(\overline L\).  Both
primitive values vanish, so

\[
 \int_\alpha^\beta f(t)^m\,dt
 =u\bigl(P_m(\beta)-P_m(\alpha)\bigr)=0
\]

for all sufficiently large \(m\).  Choosing one sufficiently large \(M\)
and restricting to exponents \(Mn\) shows that every positive power of
\(f^M\) has zero integral from \(\alpha\) to \(\beta\).  The
one-variable polynomial moment lemma forces \(f^M=0\), hence \(f=0\).
\(\square\)

For \(q=xy=x(w+vx)/u\), the two roots are \(0\) and \(-w/v\).  Equivalently,
the exact averaged-primitive criterion is especially transparent:

\[
 h\in D((xy))
 \Longleftrightarrow
 \frac{{\cal T}(h)}{z}\in A,\qquad z=w+vx=uy.                 \tag{6.9}
\]

Indeed, the zero-constant primitive is \(x{\cal T}(h)/u\), and division by
\(xy=xz/u\) leaves \({\cal T}(h)/z\).  Thus pure-power membership makes
the primitive vanish at both \(z=w\) and \(z=0\), which is exactly the
two-root moment obstruction above.

The first invariant-content controls also close.

### Theorem 6.7 (linear invariant-content carriers)

For \(a\in\{u,v\}\) and \(\ell\in\{x,y\}\),

\[
 D((a\ell))
\]

is Mathieu--Zhao.

**Proof.**
By the simultaneous exchanges \(u\leftrightarrow v\) and
\(x\leftrightarrow y\), it is enough to treat \((ux)\) and \((uy)\).
Since \(u\in\ker D\),

\[
 D(uxA)=uD(xA),\qquad D(uyA)=uD(yA).                          \tag{6.10}
\]

First suppose that \(f^m\in D(uxA)\) for every sufficiently large \(m\).
Then \(f^m\in D(xA)\), and the proof of Theorem 6.4 gives

\[
 f=uF,\qquad F\in B[x].
\]

For fixed \(g\), apply the termwise local-slice proof of Theorem 6.3 with
the fixed multiplier \(g/u\).  This merely adds one fixed power of \(u\)
to the coefficient denominators; the factor \(u^m\) in \(f^m\) still
clears them for large \(m\).  Hence

\[
 \frac{gf^m}{u}\in D(xA)\qquad(m\gg0).
\]

Multiplication by \(u\) and (6.10) give \(gf^m\in D(uxA)\).

Now suppose that \(f^m\in D(uyA)\) eventually.  Since
\(D(uyA)\subseteq uA\), primality of \(u\) gives \(u\mid f\).  Also
\(f^m\in D(yA)\), so the \(y\)-version of Theorem 6.4 gives

\[
 f=vF,\qquad F\in B[y].
\]

The reduction map

\[
 B[y]/(u)=k[v,w,y]\longrightarrow A/(u)=k[v,x,y],
 \qquad w\longmapsto-vx,
\]

is injective.  Therefore \(B[y]\cap uA=uB[y]\), and \(u\mid f=vF\)
forces \(F=uH\) with \(H\in B[y]\).  Thus \(f=uvH\).

Use the local slice \(t=y/v\), so \(A_v=B_v[t]\) and
\(D=\partial_t\).  The fixed multiplier \(g/u\) has coefficients in
\(B_{uv}\).  In the termwise primitive of

\[
 \frac{gf^m}{u}=g\,u^{m-1}v^mH^m,
\]

the powers \(u^{m-1}v^m\) clear every fixed coefficient denominator and
every bounded loss caused by substituting \(t=y/v\).  For large \(m\) the
primitive is polynomial and divisible by \(y\).  Hence
\((gf^m)/u\in D(yA)\), and (6.10) gives \(gf^m\in D(uyA)\).
\(\square\)

The first rational-root ladder also closes.  This is the point at which a
single valuation coefficient is insufficient: in the coordinates below,

\[
 y=\frac{v(q_1-w)}{u^2}+\frac wu,
\]

so, for example, the coefficient of \((q_1-w)^0\) in \(u^2y^3\) is
\(w^3/u\).  Thus an arbitrary negative \(u\)-coefficient need not be
divisible by \(q_1-w\).  What survives is divisibility of the **lowest**
\(u\)-face, and that is exactly enough for the all-order moment argument.

### Theorem 6.8 (rational-root ladder)

For every \(n\ge1\), put

\[
 q_n=u^nx+w=(u^n-v)x+uy.
\]

Then \(D(q_nA)\) is Mathieu--Zhao.  More precisely, eventual pure-power
membership

\[
 f^m\in D(q_nA)\qquad(m\gg0)
\]

forces

\[
 f\in uB[q_n],\qquad B=k[u,v,w],
\]

and this local form gives \(gf^m\in D(q_nA)\) for every fixed \(g\in A\)
and all sufficiently large \(m\).

There is also an exact untruncated membership criterion.  If \(P\in A\)
is any primitive of \(h\in\operatorname {im}D\), substitute

\[
 x=ut,\qquad y=-(u^n-v)t
\]

in the residue of \(P\) modulo \(q_n\).  Then

\[
 \boxed{\quad
 h\in D(q_nA)
 \Longleftrightarrow
 P( u,v,ut,-(u^n-v)t)\in k[u,v,u^{n+1}t].
 \quad}                                                       \tag{6.11}
\]

Equivalently, every normalized residue monomial \(u^av^bt^c\) must
satisfy \(a\ge(n+1)c\).

**Proof.**
The quotient parametrization in (6.11) kills
\((u^n-v)x+uy\).  Its kernel is exactly \((q_n)\): the image has
transcendence degree three, while the primitive linear polynomial \(q_n\)
is prime.  Moreover

\[
 w=uy-vx\longmapsto-u^{n+1}t.
\]

Thus the image of \(B\) in the quotient is
\(k[u,v,u^{n+1}t]\), whose monomial support is exactly
\(a\ge(n+1)c\).  A primitive \(P\) can be corrected into \(q_nA\) exactly
when its residue modulo \(q_n\) lies in this invariant image, proving
(6.11).

For the all-order assertion, embed

\[
 A\hookrightarrow B_u[q_n],\qquad
 x=\frac{q_n-w}{u^n},\qquad
 y=\frac{v(q_n-w)+u^nw}{u^{n+1}}.                            \tag{6.12}
\]

Write \(\delta=q_n-w\), and let \(\nu\) be the lowest \(u\)-exponent in
this Laurent presentation.  We need the following lowest-face fact:
if \(0\ne H\in A\), \(\nu(H)=R<0\), and

\[
 H=u^RH_R(v,w,q_n)+\text{higher \(u\)-terms},
\]

then

\[
 \delta^{\lceil-R/(n+1)\rceil}\mid H_R.                      \tag{6.13}
\]

Indeed, the associated graded algebra of (6.12) has

\[
 \operatorname {in}(x)=\delta u^{-n},\qquad
 \operatorname {in}(y)=v\delta u^{-(n+1)},
\]

Writing \(X=\operatorname {in}(x)\) and
\(Y=\operatorname {in}(y)\), it has the presentation

\[
 \operatorname {gr}_\nu(A)
 \cong
 \frac{k[u,v,w,q_n,X,Y]}
 {(u^nX-q_n+w,\ uY-vX)},                                    \tag{6.13a}
\]

with degrees \(1,0,0,0,-n,-(n+1)\), respectively.  To see this directly,
adjoin \(w,q_n\) to the polynomial presentation of \(A\).  The relation
\(q_n-u^nx-w\) is homogeneous, while the lowest part of
\(w-uy+vx\) is \(-uY+vX\); after eliminating \(w\), the remaining
filtered ideal is principal, so there are no further initial relations.
A homogeneous monomial of degree

\[
 R=a-nd-(n+1)c<0
\]

coming from \(u^ax^dy^c\) contains \(\delta^{d+c}\), and

\[
 -R\le(n+1)(d+c).
\]

The leading relation preserves this divisibility, so every homogeneous
sum has the factor asserted in (6.13).  Notice that (6.13) concerns only
the lowest face; the displayed \(u^2y^3\) term explains why the analogous
claim for every coefficient is false.

Now write

\[
 f=u^rF(v,w,q_n)+\text{higher \(u\)-terms},\qquad F\ne0.
\]

Since \(D(q_n)=u^{n+1}\), the unique primitive of \(f^m\) vanishing at
\(q_n=0\) is

\[
 Q_m=\frac1{u^{n+1}}\int_0^{q_n}f(s)^m\,ds.                  \tag{6.14}
\]

If \(f^m\in D(q_nA)\), its primitive in \(q_nA\) also vanishes at
\(q_n=0\), hence equals \(Q_m\).  Suppose \(r\le0\).  The lowest
\(u\)-degree of (6.14) is \(rm-(n+1)<0\), with coefficient

\[
 G_m(q_n)=\int_0^{q_n}F(s)^m\,ds.
\]

By (6.13), \(\delta\mid G_m\).  Evaluation at \(q_n=w\) therefore gives

\[
 \int_0^wF(s)^m\,ds=0                                       \tag{6.15}
\]

for every sufficiently large \(m\).  Choose one sufficiently large
\(M\) and apply the one-variable polynomial moment lemma over
\(k(v,w)\) to the polynomial \(F^M\).  Equation (6.15) for the exponents
\(Mn\), \(n\ge1\), forces \(F^M=0\), a contradiction.  Hence \(r\ge1\),
and the Laurent expansion gives \(f=u^rH\) with \(H\in B[q_n]\).

Finally fix \(g\in A\).  In \(B_u[q_n]\), the coefficients of \(g\) have
only a bounded \(u\)-denominator.  The zero-constant primitive

\[
 \frac1{u^{n+1}}\int_0^{q_n}g(s)u^{rm}H(s)^m\,ds
\]

belongs to \(q_nB[q_n]\subseteq q_nA\) once \(m\) is large enough: the
factor \(u^{rm}\) clears both the fixed denominator of \(g\) and
\(u^{n+1}\).  Its derivative is \(gf^m\), completing the proof.
\(\square\)

For \(n=1\), the exact \(a\ge2c\) checker retains seventeen of the 1,055
sparse seeds on the genuinely eventual window \(m=4,5,6\).  All
seventeen lie in \(uB[q_1]\), and none gives a bounded mixed-tail
obstruction.  The bounded census is only a regression; Theorem 6.8 is
the independent all-order conclusion.

The first tied two-prime family closes by applying the lowest-face
argument independently at both plinth primes.

### Theorem 6.9 (two-prime monomial rational-root grid)

For every \(r,s\ge1\), put

\[
 q_{r,s}=u^rv^sx+w.
\]

Then \(D(q_{r,s}A)\) is Mathieu--Zhao.  Eventual pure-power membership
forces the stronger normal form

\[
 f\in uvB[q_{r,s}].                                         \tag{6.16}
\]

If \(P\) is any primitive of \(h\in\operatorname {im}D\), there is also
the exact criterion

\[
 h\in D(q_{r,s}A)
 \Longleftrightarrow
 P(u,v,ut,-(u^rv^s-v)t)
 \in k[u,v,u^{r+1}v^st].                                   \tag{6.17}
\]

Equivalently, every residue monomial \(u^av^bt^c\) must satisfy the two
simultaneous inequalities

\[
 a\ge(r+1)c,\qquad b\ge sc.                                 \tag{6.18}
\]

**Proof.**
The substitution in (6.17) kills

\[
 q_{r,s}=(u^rv^s-v)x+uy
\]

and sends \(w\) to \(-u^{r+1}v^st\).  The carrier is a primitive linear
prime, and the parametrized image has transcendence degree three, so the
kernel is exactly \((q_{r,s})\).  The invariant image is therefore
\(k[u,v,u^{r+1}v^st]\), proving (6.17)--(6.18) by the same
primitive-correction argument as in Theorem 6.8.

Put \(q=q_{r,s}\) and \(\delta=q-w\).  In the two-prime localization,

\[
 A\hookrightarrow B_{uv}[q],\qquad
 x=\frac{\delta}{u^rv^s},\qquad
 y=\frac{\delta}{u^{r+1}v^{s-1}}+\frac wu,\qquad
 D=u^{r+1}v^s\partial_q.                                    \tag{6.19}
\]

First use the lowest \(u\)-exponent \(\nu_u\), treating \(v\) as a
coefficient.  If \(H\in A\) has lowest degree \(R<0\), then

\[
 \delta^{\lceil-R/(r+1)\rceil}
 \mid\operatorname {in}_{\nu_u}(H).                         \tag{6.20}
\]

Indeed, with \(X=\operatorname {in}(x)\) and
\(Y=\operatorname {in}(y)\), the associated graded presentation is

\[
 \frac{k[u,v,w,q,X,Y]}
 {(u^rv^sX-q+w,\ uY-vX)}
\]

with \(u\)-degrees \(1,0,0,0,-r,-(r+1)\).  A negative homogeneous
monomial of degree \(R=a-rd-(r+1)c\) contains
\(\delta^{d+c}\), and
\(-R\le(r+1)(d+c)\).

Now use the lowest \(v\)-exponent \(\nu_v\), treating \(u\) as a
coefficient.  Formula (6.19) gives \(v\)-degrees

\[
 \deg_v(x)=-s,\qquad \deg_v(y)=1-s.
\]

The same graded monomial estimate gives, for a lowest degree \(S<0\),

\[
 \delta^{\lceil-S/s\rceil}
 \mid\operatorname {in}_{\nu_v}(H).                         \tag{6.21}
\]

For \(s>1\), both negative initial generators contain \(\delta\), and
\(-S\le s(d+c)\).  For \(s=1\), \(y\) has degree zero and only unmatched
negative powers of \(x\) contribute to \(S<0\); each such power contains
\(\delta\), giving the same estimate.  As in (6.13), these are
lowest-face statements, not assertions about arbitrary Laurent
coefficients.

Suppose \(f^m\in D(qA)\) for every sufficiently large \(m\).  The unique
primitive vanishing at \(q=0\) is

\[
 Q_m=\frac1{u^{r+1}v^s}\int_0^q f(z)^m\,dz.                 \tag{6.22}
\]

It equals the assumed primitive in \(qA\).  Write the lowest \(u\)-face
of \(f\) as \(u^\rho F(v,w,q)\).  If \(\rho\le0\), the lowest face of
(6.22) has degree \(\rho m-(r+1)<0\).  By (6.20), its coefficient
\(\int_0^qF(z)^m\,dz\) is divisible by \(q-w\), so

\[
 \int_0^wF(z)^m\,dz=0\qquad(m\gg0).
\]

The one-variable polynomial moment lemma over \(k(v,w)\) gives a
contradiction.  Hence \(\rho\ge1\).

Apply the same argument to the lowest \(v\)-face
\(v^\sigma G(u,w,q)\).  Equation (6.21) and the denominator \(v^s\) in
(6.22) again produce all sufficiently large moments from \(0\) to \(w\);
hence \(\sigma\ge1\).  The finite Laurent expansion in
\(B_{uv}[q]\) now has every \(u\)-exponent and every \(v\)-exponent
positive, proving (6.16).

Finally, for fixed \(g\in A\), its expression in \(B_{uv}[q]\) has
bounded \(u\)- and \(v\)-denominators.  Substituting
\(f=uvH\), \(H\in B[q]\), in the zero-constant primitive of \(gf^m\)
supplies \((uv)^m/(u^{r+1}v^s)\), which clears those fixed losses for
large \(m\).  The resulting primitive lies in \(qB[q]\subseteq qA\).
\(\square\)

For \(q_{1,1}=uvx+w\), the exact cone \(a\ge2c,\ b\ge c\) retains eight
of the 1,055 sparse seeds on powers \(4,5,6\).  All eight lie in
\(uvB[q_{1,1}]\), and none gives a bounded mixed-tail obstruction.

The monomial coefficient, rather than the particular intercept \(w\), is
the real input to the proof.

### Theorem 6.10 (monomial slope with arbitrary invariant intercept)

Let

\[
 r\ge1,\qquad s\ge0,\qquad b\in B=k[u,v,w],
\]

and put

\[
 q=u^rv^sx+b.
\]

Then \(D(qA)\) is Mathieu--Zhao.  If \(b\ne0\), eventual pure-power
membership forces

\[
 f\in
 \begin{cases}
  uB[q],&s=0,\\
  uvB[q],&s\ge1.
 \end{cases}                                                \tag{6.23}
\]

**Proof.**
Set \(a=u^rv^s\) and \(\delta=q-b\).  In \(B_{uv}[q]\),

\[
 x=\frac{\delta}{a},\qquad
 y=\frac{\delta}{u^{r+1}v^{s-1}}+\frac wu,\qquad
 D=u^{r+1}v^s\partial_q.                                    \tag{6.24}
\]

For \(b\ne0\), the lowest-\(u\) proof of Theorem 6.9 applies verbatim
with \(q-w\) replaced by \(q-b\).  A nonpositive lowest \(u\)-order of
\(f\) would force

\[
 \int_0^bF(z)^m\,dz=0\qquad(m\gg0),
\]

contradicting the one-variable polynomial moment lemma over
\(\operatorname {Frac}(B)\).  Hence \(f\in uB_{v^{\pm1}}[q]\).
When \(s\ge1\), the independent lowest-\(v\) argument gives positive
\(v\)-order as well, so \(f\in uvB[q]\).  If \(s=0\), no \(v\)-denominator
occurs and the first conclusion already gives \(f\in uB[q]\).
The content in (6.23) clears the fixed denominators of \(g\), together
with \(u^{r+1}v^s\), in the zero-constant primitive of \(gf^m\).

It remains to take \(b=0\), so \(q=ax\) and

\[
 D(qA)=aD(xA).
\]

Eventual membership implies \(f^m\in aA\), hence \(u\mid f\), and also
\(v\mid f\) when \(s\ge1\).  It also implies \(f^m\in D(xA)\), so
Theorem 6.4 gives \(f\in uB[x]\).  Reduction modulo \(v\), using the
injection

\[
 B[x]/(v)=k[u,w,x]\hookrightarrow A/(v)=k[u,x,y],
 \qquad w\longmapsto uy,
\]

shows that \(v\mid f\) forces the \(B[x]\)-cofactor to be divisible by
\(v\).  Thus the same normal forms as in (6.23) hold.  The local-slice
proof for \(D(xA)\) accepts the fixed multiplier \(g/a\), because the
growing \(u\)- and, when necessary, \(v\)-content clears its denominator.
Multiplication by \(a\) completes the proof. \(\square\)

The same valuation argument closes a dense class of genuinely
nonmonomial slopes.  It also isolates the exact conductor alignment where
the argument can fail.

### Theorem 6.11 (coprime invariant slope)

Let \(0\ne a,b\in B=k[u,v,w]\) satisfy \(\gcd(a,b)=1\), and put

\[
 q=ax+b.
\]

Then \(D(qA)\) is Mathieu--Zhao.

**Proof.**
Work in \(B_{ua}[q]\), where

\[
 x=\frac{q-b}{a},\qquad
 y=\frac{w+v(q-b)/a}{u},\qquad
 D=ua\,\partial_q.                                          \tag{6.26}
\]

Let \(\pi\ne u\) be an irreducible factor of \(a\), of multiplicity
\(e\).  In the lowest \(\pi\)-filtration, \(x\) has degree \(-e\).
If \(\pi\ne v\), the negative face of \(y\) has the same degree; if
\(\pi=v\), it has degree \(1-e\), with degree zero when \(e=1\).
Every negative lowest face is therefore divisible by a positive power of
\(\delta=q-b\), exactly as in (6.20)--(6.21).  Since
\(\gcd(a,b)=1\), the image of \(b\) modulo \(\pi\) is nonzero.
If the lowest \(\pi\)-order of \(f\) were nonpositive, the
zero-constant primitive

\[
 \frac1{ua}\int_0^q f(z)^m\,dz
\]

would consequently give

\[
 \int_0^{\bar b}F(z)^m\,dz=0\qquad(m\gg0)
\]

over \(\operatorname {Frac}(B/(\pi))\).  The polynomial moment lemma is
a contradiction.  Thus \(f\) has positive order at every prime factor
of \(a\) other than \(u\).

If \(u\mid a\), the \(u\)-lowest face is again a positive power of
\(q-b\); coprimality makes \(\bar b\ne0\), and the same argument forces
positive \(u\)-order.  Suppose \(u\nmid a\).  Then \(x\) has \(u\)-degree
zero, while the degree-\(-1\) face of \(y\) is

\[
 w+\frac{v(q-\bar b)}{\bar a}.
\]

Its root in the \(q\)-line is

\[
 \beta=\bar b-\frac{\bar a w}{v}.
\]

If \(\beta\ne0\), divisibility of the lowest primitive face gives all
sufficiently large moments from \(0\) to \(\beta\), again impossible.
Hence \(f\) has positive \(u\)-order in this case.

It remains to suppose \(\beta=0\), or equivalently

\[
 v\bar b=\bar a\,w.                                        \tag{6.25}
\]

Reduction modulo \(u\) sends \(w\) to \(-vx\), so (6.25) says exactly
that \(q=0\) in \(A/(u)\), that is,

\[
 q\in uA.                                                    \tag{6.27}
\]

Since \(u\) is invariant,

\[
 D(qA)\subseteq uA.
\]

Eventual membership \(f^m\in D(qA)\) therefore gives \(u\mid f\), because
\(u\) is prime in \(A\).  Thus positive \(u\)-content holds on the
aligned face as well, without a moment argument.

It follows that, in the finite Laurent expansion (6.26), \(f\) contains
one copy of \(u\) and one copy of every irreducible factor of \(a\).
These growing powers clear \(ua\) and all fixed coefficient denominators
in the primitive of \(gf^m\).  The primitive then belongs to
\(qB[q]\subseteq qA\), proving the result. \(\square\)

The common invariant factor is removed by an abstract scaling argument;
no further conductor calculation is needed.

### Lemma 6.12 (Mathieu scaling)

Let \(R\) be a commutative domain, let \(M\subseteq R\) be a
Mathieu--Zhao subspace, and let \(0\ne c\in R\) satisfy

\[
 cM\subseteq M.
\]

Then \(cM\) is Mathieu--Zhao.

**Proof.**
Suppose \(f^m\in cM\) for all sufficiently large \(m\).  Then
\(f^m\in M\) eventually.  Choose one sufficiently large \(N\), and put

\[
 h=\frac{gf^N}{c}\in R
\]

for a fixed \(g\in R\); this is an element of \(R\) because
\(f^N/c\in M\subseteq R\).  Since \(M\) is Mathieu--Zhao,

\[
 h f^{m-N}\in M\qquad(m\gg0).
\]

But

\[
 h f^{m-N}=\frac{gf^m}{c}.
\]

Therefore \(gf^m\in cM\) for all sufficiently large \(m\). \(\square\)

### Corollary 6.13 (complete invariant-linear carrier theorem)

Let \(0\ne a\in B=k[u,v,w]\) and \(b\in B\), and put

\[
 q=ax+b.
\]

Then \(D(qA)\) is Mathieu--Zhao, with no coprimality assumption on
\(a,b\).

**Proof.**
Write

\[
 c=\gcd(a,b),\qquad a=ca_0,\qquad b=cb_0,
\]

where \(c=a\) and \(a_0=1,b_0=0\) when \(b=0\).  Put
\(q_0=a_0x+b_0\).  If \(b\ne0\), Theorem 6.11 proves that

\[
 M=D(q_0A)
\]

is Mathieu--Zhao; if \(b=0\), this is Theorem 6.4 for \(q_0=x\).
Because \(c\in B=\ker D\),

\[
 D(qA)=D(cq_0A)=cD(q_0A)=cM,
\]

and \(M\) is a \(B\)-module, so \(cM\subseteq M\).  Lemma 6.12 completes
the proof. \(\square\)

Degree one on the generic orbit has an intrinsic description that removes
the remaining denominator descent.

### Theorem 6.14 (irreducible invariant-affine carrier theorem)

Let \(h\in A\) be irreducible and satisfy

\[
 D^2h=0,\qquad Dh\ne0.
\]

Then \(D(hA)\) is Mathieu--Zhao.

**Proof.**
Put \(s=Dh\in B\).  The plinth identity

\[
 \operatorname {im}D\cap B=(u,v)B
\]

gives \(s=b_1u+b_2v\) for some \(b_1,b_2\in B\).  Consequently

\[
 h=b_0+b_1x+b_2y,\qquad b_0\in B,                            \tag{6.28}
\]

because the difference has derivative zero.  Conversely, every
expression (6.28) lies in \(\ker D^2\).  With \(w=uy-vx\), direct
elimination gives

\[
 x=\frac{u(h-b_0)-b_2w}{s},\qquad
 y=\frac{v(h-b_0)+b_1w}{s}.                                 \tag{6.29}
\]

Thus \(A\hookrightarrow B_s[h]\) and \(D=s\partial_h\).

Let \(\pi\) be an irreducible factor of \(s\), and work in the DVR
\(B_{(\pi)}\).  At least one of \(u,v\) is a unit there.  If \(u\) is a
unit, (6.29) becomes

\[
 x=\frac{u(h-\widetilde\beta_\pi)}s,\qquad
 \widetilde\beta_\pi=b_0+\frac{b_2w}{u},
 \qquad y=\frac vu x+\frac wu.
\]

If \(u\) is not a unit, then \(\pi=u\), \(v\) is a unit, and the analogous
formula is

\[
 y=\frac{v(h-\widetilde\beta_\pi)}s,\qquad
 \widetilde\beta_\pi=b_0-\frac{b_1w}{v},
 \qquad x=\frac uv y-\frac wv.
\]

Write \(\beta_\pi\) for the residue of
\(\widetilde\beta_\pi\) in \(\operatorname {Frac}(B/(\pi))\).
These formulas give the associated graded algebra directly: every
negative lowest \(\pi\)-face of an element of \(A\) is divisible by a
positive power of

\[
 h-\beta_\pi.                                                \tag{6.30}
\]

Suppose \(f^m\in D(hA)\) eventually.  The unique primitive vanishing at
\(h=0\) is

\[
 P_m=\frac1s\int_0^h f(z)^m\,dz\in hA.                      \tag{6.31}
\]

If the lowest \(\pi\)-order of \(f\) were nonpositive and
\(\beta_\pi\ne0\), (6.30)--(6.31) would give

\[
 \int_0^{\beta_\pi}F(z)^m\,dz=0\qquad(m\gg0),
\]

contradicting the one-variable polynomial moment lemma.

If \(\beta_\pi=0\), the displayed unit-coordinate formula, together with
\[
 v(ub_0+b_2w)-u(vb_0-b_1w)=sw,
\]
shows that both constant terms in (6.29) vanish modulo \(\pi\).
The identities

\[
 uh=sx+(ub_0+b_2w),\qquad
 vh=sy+(vb_0-b_1w)
\]

then show \(h\in\pi A\): use the first identity when
\(\pi\ne u\), and the second when \(\pi=u\).  This is impossible because
\(h\) is irreducible and noninvariant whereas \(\pi\in B\).  Thus every
\(\beta_\pi\) is nonzero, and the moment argument forces positive
\(\pi\)-order of \(f\).

Hence the order of \(f^m\) at every prime denominator in (6.29) grows
linearly.  For fixed \(g\), these powers clear \(s\) and all fixed
coefficient losses in

\[
 \frac1s\int_0^h g(z)f(z)^m\,dz.
\]

The resulting primitive has nonnegative order at every factor of \(s\),
so it lies in \(B[h]\).  It vanishes at \(h=0\), hence belongs to
\(hB[h]\subseteq hA\).  Therefore \(gf^m\in D(hA)\) for all sufficiently
large \(m\).
\(\square\)

The degree-one result bootstraps through every multiplicity.

### Corollary 6.15 (invariant-affine thickenings)

If \(0\ne h\in\ker D^2\setminus\ker D\) and \(d\ge1\), then

\[
 D(h^dA)
\]

is Mathieu--Zhao.

**Proof.**
Factor \(h\) in the UFD \(A\).  The kernel of a locally nilpotent
derivation is factorially closed, so every invariant irreducible factor
belongs to \(B\).  Since \(h\) has degree one in the generic orbit
coordinate, there is exactly one noninvariant irreducible factor, with
multiplicity one.  Hence

\[
 h=cH,\qquad c\in B,\qquad H\ \text{noninvariant irreducible}.
\]

Also \(D^2H=0\), so Theorem 6.14 proves \(D(HA)\) Mathieu--Zhao.
It is enough first to treat \(H^d\).  Since

\[
 D(H^dA)\subseteq H^{d-1}A,
\]

eventual pure membership forces \(H\mid f\).  The degree-one theorem
gives, for fixed \(g\) and large \(m\), a primitive

\[
 D(P_m)=gf^m,\qquad P_m\in HA.
\]

Write \(P_m=H^ra\), with \(1\le r<d\).  Since \(gf^m\) is divisible by
\(H^m\),

\[
 D(P_m)=rH^{r-1}D(H)a+H^rD(a)\in H^rA
\]

for \(m\ge d\).  Reduction modulo \(H^r\) gives
\(H\mid D(H)a\).  The polynomials \(H,D(H)\) are coprime: otherwise
irreducibility would give \(H\mid D(H)\), which for a locally nilpotent
derivation forces \(H\in\ker D\).  Thus \(H\mid a\).  Iteration yields
\(P_m\in H^dA\).

Finally,

\[
 D(h^dA)=D(c^dH^dA)=c^dD(H^dA),
\]

and \(D(H^dA)\) is a \(B\)-module.  Thus Lemma 6.12 removes the
invariant factor \(c^d\). \(\square\)

These results close every principal carrier for the present rank-two
linear locally nilpotent derivation.

### Theorem 6.16 (complete principal-ideal theorem for the model LND)

For every \(0\ne q\in A=k[u,v,x,y]\),

\[
 D(qA)
\]

is Mathieu--Zhao.

**Proof.**
If \(q\in B\), then \(D(qA)=qD(A)\).  Theorem 6.1 applied to \(I=A\)
proves \(D(A)\) Mathieu--Zhao; it is a \(B\)-module, so Lemma 6.12
applies.

Assume \(q\notin B\), and view \(q\) in
\(L[x]\), \(L=\operatorname {Frac}(B)\).  If it has at least two distinct
roots over \(\overline L\), Theorem 6.6 gives zero eventual-power radical.
Suppose it has exactly one root.

Factor \(q\) in the UFD \(A\).  Every invariant irreducible factor lies
in \(B\).  A noninvariant irreducible factor remains irreducible after
localizing at \(u\), because it is not associated to \(u\).  Since

\[
 A_u=B_u[x],
\]

Gauss's lemma and the one-root hypothesis force every such factor to
have degree one over \(L\).  Any two are associates in \(A_u\); the only
units introduced by the localization are scalar multiples of powers of
\(u\), and neither factor is divisible by \(u\), so they are already
associates in \(A\).  Therefore

\[
 q=cH^d
\]

for some \(0\ne c\in B\), \(d\ge1\), and one noninvariant irreducible
\(H\) of generic degree one.  Hence \(D^2H=0\), Corollary 6.15 proves
\(D(H^dA)\) Mathieu--Zhao, and its \(B\)-module structure lets
Lemma 6.12 remove \(c\). \(\square\)

The principal-ideal frontier for \(D=u\partial_x+v\partial_y\) is
therefore closed.  For the ideal-image/LNED strengthening, the remaining
model problem consists of genuinely nonprincipal ideals not covered by
plinth-power saturation or the finite-residual theorems.  For the usual
LND Conjecture \(I=A\), this model was already covered by the
rank-at-most-two linear theorem cited in Section 1; a full proof still
requires control of locally nilpotent derivations whose slice/conductor
geometry is not equivalent to this rank-two linear normal form.

## 7. Exact replay

Run

```bash
.venv/bin/python scripts/verify_lnd_radical_slice_fibers.py
```

The checker constructs the radical complete intersection of six points in
three vertical fibers, verifies reducedness, compares the generic
primitive-remainder conditions with the three interval obstructions, and
replays pure and mixed powers through exponent twelve.  It also checks two
nonreduced length-two residual charts for the carrier \(q=s^2\), one on and
one off the carrier.  The bounded powers are regression checks only.  The
all-order conclusions are Theorems 3.1, 3.2, and 4.1.  The normalized
nonmonic support-weight calculation of Theorem 4.2 has a separate search
replay.

Run

```bash
.venv/bin/python scripts/search_lnd_nonmonic_degree_drop.py
```

for the exact \(p=xs-1\) membership compiler.  It tests \(x^cp^dJ\) for
\((c,d)=(0,1),(0,2),(1,1)\), three primary residual schemes at the
degree-drop point, 256 monomial/binomial seeds, and six multipliers.  The
primitive-carrier support-weight assertion is a theorem regression.  Every
finite prefix remains only a replay, even though Theorem 4.2 supplies the
all-order exclusion.

Run

```bash
.venv/bin/python scripts/search_lnd_plinth_ideal_images.py
```

for the exact homogeneous plinth-divisor search described in Section 6.

Run

```bash
.venv/bin/python scripts/search_lnd_reducible_plinth.py
```

for the exact \(\mathbf Q[x]\)-module searches on the two plinth branches,
covering \(y\partial_z\), \(y^2\partial_z\), and the branch-asymmetric
\((y^2+x)\partial_z\) profile.  Singular is required.  Module and
quotient/kernel membership are exact, but each six-power search window is
bounded.

Run

```bash
.venv/bin/python scripts/search_lnd_crossing_plinth.py
```

for the exact \(\mathbf Q[u,v]\)-module and finite quotient/kernel search
at the crossing \(uv=0\).  Singular is required.  The individual
membership decisions are exact; the six-power pure window and
three-exponent mixed window are bounded counterexample searches only.

Run

```bash
.venv/bin/python scripts/search_lnd_nonprincipal_plinth.py
```

for the exact homogeneous two-slice search at the nonprincipal plinth
ideal \((u,v)\).  Singular is required.  The finite quotient/kernel
membership tests are exact; the pure and mixed exponent windows remain
bounded.

Run

```bash
.venv/bin/python scripts/search_lnd_positive_dimensional_plinth.py
```

for the exact free-line and positive-dimensional nilpotent-jet search.
The checker also verifies (6.2) degree by degree through total degree
eight.  That bounded identity is a regression for the free-line corollary
of Theorem 6.1.  The power windows are bounded, while the theorem proves
all five displayed ideals safe independently of those windows.

Run

```bash
.venv/bin/python scripts/search_lnd_principal_conductor.py
```

for the exact \(a\ge c\) valuation-face search at \(I=(x)\).  Singular
constructs the primitives.  Every membership decision is exact; the pure
and mixed exponent windows are bounded.  The same checker verifies the
algebraic face (6.3), including first- and second-power membership and
third-power failure, and the exact identity
\({\cal T}(D(xya))=(uy)a\) in the local slice.  It also checks
\(D(uxa)=uD(xa)\) and \(D(uya)=uD(ya)\).  For
\(q_1=ux+(uy-vx)\), it uses the exact normalized support cone
\(a\ge2c\), checks the rational-root normalization through \(q_4\), and
finds seventeen eventual-window survivors, all in \(uB[q_1]\), with no
bounded mixed-tail obstruction.  For \(q_{1,1}=uvx+(uy-vx)\), it checks
the simultaneous cone \(a\ge2c,\ b\ge c\), retains eight
eventual-window survivors, and verifies that all eight lie in
\(uvB[q_{1,1}]\).  The normalization identities are also checked for
\(1\le r,s\le3\).  For one sample
\(h=b_0+b_1x+b_2y\), it additionally verifies \(D^2h=0\) and both
inverse identities (6.29).  Theorems 6.4, 6.6, 6.7, 6.8, 6.9, 6.10,
6.11, 6.14, and 6.16, Lemma 6.12, and Corollaries 6.5, 6.13, and 6.15
are the separate all-order arguments; they are not inferred from the
bounded calculation or the coordinate-identity regression.
