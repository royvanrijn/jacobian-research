# Quartic endpoint semigroups: the finite-feasibility gate

## Status

The cusp/connector pole-semigroup data do **not** close either surviving
geometric-degree-four packet.  Even after conductor pairing is restored,
the proposed inputs do not define an unconditional finite lattice problem.

What is proved here is a no-finiteness theorem and an exact bounded
compiler.  For externally supplied bounds on connector count, displayed
pole order, and connector contact, the compiler enumerates every endpoint
and sheet pairing and exactly counts the resulting decorated lattice
records.  Geometric degree four does not currently supply those bounds.  No
quartic Keller normalization is constructed or excluded.

The computation is
[`cas/experiment_quartic_endpoint_semigroups.py`](cas/experiment_quartic_endpoint_semigroups.py).
It continues the negative result in
[`JC2_GLOBAL_COX_PACKET_ATTACK.md`](JC2_GLOBAL_COX_PACKET_ATTACK.md),
Sections 6.4--6.5.

## 1. Local semigroups and conductors

At the clean cusp, the normalization parameter has coordinate pole orders
\(2,3\).  The numerical semigroup and its conductor are

\[
S_{\mathrm{cusp}}=\langle 2,3\rangle
 =\{0,2,3,4,\ldots\},
\qquad c(S_{\mathrm{cusp}})=2.
\tag{1.1}
\]

The local odd-square multiplier has contact two.  This agrees with the
completed identity

\[
\ell=4r-9T^2
\tag{1.2}
\]

at \(r=0\).

For \(n\ge1\) and \(\lambda\ne0\), consider the smooth connector

\[
d_{n,\lambda}=xy-y^{n+1}-\lambda.
\tag{1.3}
\]

Its normalization is

\[
k[x,y]/(d_{n,\lambda})\simeq k[t,t^{-1}],
\qquad
x=t^n+\lambda t^{-1},\quad y=t.
\tag{1.4}
\]

At \(0,\infty\), the displayed pole vectors of the two generators are

\[
p(x)=(1,n),\qquad p(y)=(0,1).
\tag{1.5}
\]

Pole vectors are not additive under multiplication when a zero of one
factor cancels a pole of another.  The additive object used by the polar
filtration is therefore first the monoid of *displayed polar bounds*

\[
S_n^{\mathrm{bound}}
 =\mathbb N(1,n)+\mathbb N(0,1)
 =\{(a,b)\in\mathbb N^2:b\ge na\}.
\tag{1.6}
\]

It has no conductor translate in \(\mathbb N^2\): after fixing any proposed
conductor vector, increasing only the first coordinate eventually violates
\(b\ge na\).

Initial residues change the result.  The exact cancellation

\[
\frac{x-y^n}{\lambda}=t^{-1}
\tag{1.7}
\]

adds the missing pole vector \((1,0)\).  Together with \(p(y)=(0,1)\), it
shows that every vector in \(\mathbb N^2\) is realized as a pole vector by
a Laurent polynomial in the coordinate algebra.  The residue-completed
polar-bound monoid is therefore

\[
S_n^{\mathrm{res\mbox{-}bound}}=\mathbb N^2,
\qquad c(S_n^{\mathrm{res\mbox{-}bound}})=(0,0).
\tag{1.8}
\]

The actual additive multivaluation semigroup uses signed orders.  Since the
coordinate algebra is \(k[t,t^{-1}]\), it is

\[
\Gamma_{\mathrm{conn}}
=\{(\operatorname{ord}_0f,\operatorname{ord}_\infty f):0\ne f\in
  k[t,t^{-1}]\}
=\{(u,v)\in\mathbb Z^2:u+v\le0\}.
\tag{1.9}
\]

It is generated as a monoid by
\((1,-1),(-1,1),(0,-1)\).  Formula (1.9) follows by taking the least and
greatest Laurent exponents; conversely \(t^u+t^{-v}\) realizes every pair
with \(u+v<0\), and a monomial realizes equality.

This is the semigroup form of the automorphism reduction already seen in
the endpoint Laurent compiler.  It separates four notions which must not
be conflated:

1. \(S_n^{\mathrm{bound}}\) is the monoid of additive polar bounds from the
   displayed restrictions of \(x,y\);
2. \(S_n^{\mathrm{res\mbox{-}bound}}\) includes simultaneous residue
   cancellation;
3. \(\Gamma_{\mathrm{conn}}\) is the exact signed valuation semigroup; and
4. the connector normalization is already normal, so its algebra
   conductor is the unit ideal, while a target self-identification has
   conductor exponent one at each of its two normalized endpoints.

The conductor equivalence relation pairs the two connector endpoints.  It
does not follow from either numerical monoid.

## 2. The odd-square valuation vector

For \(r\) connectors, write their completed contact orders as
\(m_i\ge1\).  The endpoint orders carried by the odd-square multiplier are

\[
\delta_{\mathrm{odd}}
  =(2,m_1,m_1,\ldots,m_r,m_r).
\tag{2.1}
\]

The first entry is the cusp contact.  Each connector contributes the same
contact order at its two local copies, and the conductor relation pairs
those copies.  The completed chart permits every finite \(m_i\); no
degree-four bound on the \(m_i\) has been proved.

On the global source bridge

\[
C_h=k[x,y,s,a]/(h-a s^2),
\qquad \deg s=1,\quad\deg a=-2,
\tag{2.2}
\]

every graded piece is rank-one free over \(R=k[x,y]\):

\[
(C_h)_d=
\begin{cases}
Rs^d,&d\ge0,\\
Ra^{-d/2},&d<0\text{ even},\\
Ra^{(1-d)/2}s,&d<0\text{ odd}.
\end{cases}
\tag{2.3}
\]

Moreover

\[
(as)^2=ha,
\qquad
\operatorname{coker}\bigl((C_h)_{-1}^{\otimes2}\to(C_h)_{-2}\bigr)
\simeq R/(h).
\tag{2.4}
\]

Consequently factoriality, rank-one freeness, and the odd-square cokernel
are compatibility checks, not restrictive inequalities on the endpoint
lattice.

## 3. Uniform carriers with arbitrarily many connectors

Choose arbitrary \(r\ge1\), \(n_i\ge1\), and nonzero
\(\lambda_i\), with pairs \((n_i,\lambda_i)\) distinct, and put

\[
h_{\boldsymbol n,\boldsymbol\lambda}
 =(y^2-x^3)
   \prod_{i=1}^r
   (xy-y^{n_i+1}-\lambda_i).
\tag{3.1}
\]

Every connector factor is irreducible and smooth: as a polynomial in \(x\)
over \(k[y]\), it is primitive and linear, and its \(x\)-derivative is
\(y\), while \(y=0\) is incompatible with (1.3).  Distinct factors are
coprime.  The cusp does not divide a connector: on
\((x,y)=(t^2,t^3)\), the connector restricts to

\[
t^5-t^{3n_i+3}-\lambda_i\ne0.
\tag{3.2}
\]

Hence (3.1) is reduced, \(R=k[x,y]\) remains factorial, and (2.2)--(2.4)
hold.  Each factor supplies the connector semigroups (1.6)--(1.8).
The endpoint set admits any chosen perfect-matching conductor decoration,
and the completed local model admits arbitrary positive contacts \(m_i\).

This family is an endpoint carrier only.  The conductor pairing belongs to
the finite map and is added as decoration; (3.1) is not asserted to be a
rank-four finite normalization or a Keller map.  Its purpose is to prove
that the listed semigroup/module constraints do not bound \(r\),
\(\boldsymbol n\), or \(\boldsymbol m\).

### Theorem 3.1 -- conductor-decorated no-finiteness

The following data do not define a finite quartic endpoint enumeration:

1. the cusp semigroup \(\langle2,3\rangle\) and contact two;
2. the displayed connector polar-bound monoids and the exact signed
   two-ended valuation semigroups, including residue cancellation;
3. a conductor equivalence relation pairing connector endpoints;
4. the odd-square endpoint vector (2.1);
5. factoriality of \(k[x,y]\); and
6. rank-one freeness of every graded piece of \(C_h\).

Indeed, (3.1) supplies compatible endpoint carriers for arbitrary
connector count and arbitrary displayed pole parameters, while the
completed connector chart independently permits arbitrary positive
contact.  Therefore a finite search requires external bounds on all three.

This theorem strengthens the earlier four-filter insufficiency statement
by restoring the conductor pairing and by computing the
residue-completed connector bounds and signed semigroup.  It is not an
impossibility theorem for quartic Keller maps.

## 4. The finite lattice problem after bounds are supplied

Fix integers

\[
r\le R,\qquad 1\le n_i\le B,\qquad 1\le m_i\le C.
\tag{4.1}
\]

For a labeled \(r\)-connector packet, the finite feasibility variables are:

- one of the 24 nondegenerate ordered cusp braid pairs;
- one of the three \(2+2\) sheet matchings for each connector;
- one perfect matching of the \(2r\) normalized connector endpoints;
- the integers \(n_i\in[1,B]\) and \(m_i\in[1,C]\);
- the displayed monoids
  \(\mathbb N(1,n_i)+\mathbb N(0,1)\), their forced
  residue-completion vectors \((1,0)\), and the odd-square vector (2.1).

The exact finite constraints are

\[
\begin{aligned}
S_{\mathrm{cusp}}&=\langle2,3\rangle,
&c_{\mathrm{cusp}}&=2,\\
S_i^{\mathrm{bound}}&=\{(a,b):b\ge n_i a\},
&S_i^{\mathrm{res\mbox{-}bound}}&=\mathbb N^2,\\
\Gamma_i&=\{(u,v)\in\mathbb Z^2:u+v\le0\},&&\\
\delta_{\mathrm{odd}}&=(2,m_1,m_1,\ldots,m_r,m_r),
&p_{i,0}&\sim p_{i,\infty}.
\end{aligned}
\tag{4.2}
\]

All 24 cusp braid pairs combined with any one of the three connector
matchings generate \(S_4\).  Thus every bounded record passes the
sheet-transitivity filter.  Counting labeled endpoint pairings, the number
of formally feasible records at fixed \(r\) is

\[
24\cdot3^r\cdot(2r-1)!!\cdot B^r C^r.
\tag{4.3}
\]

The compiler evaluates (4.2), lists the holes of every displayed connector
monoid in a requested finite box, verifies its residue completion, replays
the 24-by-3 monodromy audit, checks the graded-piece formulas through a
requested cutoff, and reports (4.3).  Passing means only that the stated
finite lattice/module constraints are feasible.

## 5. Outcome for the quartic frontier

The experiment gives neither requested terminal outcome:

- it is not a clean impossibility theorem for the quartic packets; and
- it is not a concrete rank-four packet model satisfying the global
  finite-normalization and distinguished-\(\mathbb A^2\) conditions.

It does give a clean obstruction to the proposed method as currently
stated.  Before surface geometry can be applied to a finite list, one still
needs a degree-four theorem bounding:

1. the number of connector identifications;
2. the automorphism-minimized marked multivaluation height; and
3. the completed connector contact orders.

The ramified-plus-unramified two-boundary row additionally needs the marked
valuation data of the unramified component \(D\).  Since even the
one-boundary row has the infinite family above, the requested
all-packet finite enumeration cannot yet be formed from the available
semigroup, conductor, factoriality, and rank-one data.

## 6. Reproduction

Run

```bash
.venv/bin/python plane-jc/cas/experiment_quartic_endpoint_semigroups.py \
  --max-connectors 4 \
  --max-pole 8 \
  --max-contact 8 \
  --cutoff 12 \
  --output artifacts/generated-results/quartic_endpoint_semigroups.json
```

The bounds are experiment parameters, not theorem constants.
