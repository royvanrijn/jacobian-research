# A finite atlas for every odd-height fibre of one rational coordinate

Let \(K\) be a characteristic-zero field and fix
\[
 u\in K(x),\qquad u(0)=0,\quad u'(0)=1,\qquad N=\deg u\ge1.
 \tag{H1}
\]
At odd height \(a\ge1\), the fibre consists of normalized polynomial
pairs \(C(0)=-D(0)=1\) for
\[
 p=t^a(C+t^2D)/x,\qquad t^2+2xy=1,
 \tag{H2}
\]
whose positive pure moments all vanish and whose selected endpoint
coordinate is \(u\). The coefficient degrees are unrestricted.

The exact reciprocal-unit data stabilize at a sharp height bounded only
by \(N\). This does not say that the polynomials \(C,D\) stop changing:
their dependence on the height remains explicit. Moreover, each
individual reciprocal unit is present on one interval of odd heights.
The whole collection of fibres is therefore determined by one finite
atlas, computable over a presented number field without importing a
general Orbit decision theorem.

## 1. Reciprocal units and the exact height interval

Let \(\mathcal C_u\) be the marked irreducible component of
\(u(x)+u(y)=0\) through \((0,0)\). A normalized reciprocal unit is
\[
 \eta\in K(x),\qquad \eta(0)=1,\qquad
                 \eta(x)\eta(y)=1\quad\hbox{on }\mathcal C_u.
 \tag{H3}
\]
This definition does not depend on \(a\). Put
\[
 b=\frac{u}{xu'},\qquad X_\eta=\frac b2(1+\eta).
 \tag{H4}
\]
Then \(X_\eta(0)=1\). The retained fixed-coordinate theorem supplies
exactly the following interface: the reciprocal unit gives a pair at
height \(a\) if and only if both rational functions
\[
 C_a=\frac{xX_\eta^a}{u},\qquad
 D_a=-\frac{xX_\eta^{a+2}}{u}
 \tag{H5}
\]
are polynomials. When they are, they are normalized, pure to all orders,
and have coordinate \(u\); every pair arises from one such unit.

Write \(Z\) for the finite unmarked zeros of \(u\), and \(P\) for its
finite poles. For \(\beta\in Z\), write \(s_\beta=\operatorname{ord}_\beta u>0\);
for \(\beta\in P\), write \(m_\beta=-\operatorname{ord}_\beta u>0\).
A unit is rejected at **every** height if either of these conditions fails:

* \(X_\eta\) has no finite pole outside \(P\);
* \(r_\beta:=\operatorname{ord}_\beta X_\eta>0\) for every \(\beta\in Z\).

For a unit passing these height-independent conditions, set
\[
 L_\eta=\max\left(1,\max_{\beta\in Z}
                         \left\lceil\frac{s_\beta}{r_\beta}\right\rceil\right),
 \tag{H6}
\]
where the inner maximum over an empty set is omitted. At every finite
pole \(\beta\) of \(X_\eta\), put
\(t_\beta=-\operatorname{ord}_\beta X_\eta>0\), and set
\[
 U_\eta=\min_{\beta\text{ a finite pole of }X_\eta}
           \left(\left\lfloor\frac{m_\beta}{t_\beta}\right\rfloor-2\right),
 \qquad U_\eta=+\infty\quad\hbox{if }X_\eta\in K[x].
 \tag{H7}
\]

**Theorem ReciprocalUnitHeightInterval.** The complete set of heights at
which this unit contributes is
\[
             \{a\in\mathbb Z_{\ge1}:a\text{ odd},\ L_\eta\le a\le U_\eta\}.
 \tag{H8}
\]
It may be empty. In particular, a unit cannot disappear and then reappear.
If it contributes at infinitely many heights, \(X_\eta\) is polynomial;
conversely a polynomial \(X_\eta\) passing the two independent conditions
contributes at every sufficiently large odd height.

**Proof.** At a finite point \(\beta\ne0\), writing
\(s=\operatorname{ord}_\beta u\) and \(r=\operatorname{ord}_\beta X_\eta\),
the exact orders in (H5) are
\[
                   ar-s,\qquad (a+2)r-s.              \tag{H9}
\]
At a zero of \(u\), they are both nonnegative exactly when \(r>0\)
and \(ar\ge s\), giving (H6). At a pole of \(u\), nonnegative \(r\)
imposes no restriction. If \(r=-t<0\), the stronger inequality is
\((a+2)t\le m\), giving (H7). At a point with \(u\) regular and nonzero,
they are both nonnegative exactly when \(r\ge0\). At zero, the factors
\(x/u\) and \(X_\eta\) are units with value one; hence no additional
condition appears. A rational function with no finite poles is a
polynomial. Nothing is required at infinity. This proves both directions
and includes every local branch. \(\square\)

The assertion concerns the fixed unit, or equivalently the fixed
\(X_\eta\). Different units can enter and leave at different heights;
no general monotonicity of the complete low-height fibres is asserted.

## 2. Sharp stabilization in terms of the coordinate degree

Put
\[
 S=\max_{\beta\in Z}s_\beta,\qquad
 M=\max_{\beta\in P}m_\beta,
 \tag{H10}
\]
with maximum zero for an empty set. Let \(A(u)\) be the smallest positive
odd integer at least \(\max(1,S,M-1)\), and define
\[
 A_N=2\lceil N/2\rceil-1.
 \tag{H11}
\]
Thus \(A_1=A_2=1\); for larger \(N\), \(A_N=N-1\) when \(N\) is even
and \(A_N=N\) when \(N\) is odd. Define the stable unit set
\[
 \mathcal E_\infty(u)=
 \{\eta\text{ satisfying (H3)}:
       X_\eta\in K[x],\ X_\eta(\beta)=0\text{ for every }\beta\in Z\}.
 \tag{H12}
\]

**Theorem SharpOddHeightStabilization.** One has \(A(u)\le A_N\).
For every odd \(a\ge A(u)\), the complete fibre is parametrized by
exactly \(\mathcal E_\infty(u)\), using (H5). Its emptiness, cardinality,
finite or infinite character, and exact reciprocal-unit parameters are
independent of \(a\) in this range. The same statements hold for the
geometric fibres. The degree-only cutoff \(A_N\) is sharp for every
\(N\ge3\). The coordinate-specific \(A(u)\) is a sufficient cutoff;
it need not be the least stabilization height for that particular coordinate.

**Proof.** Since the zero fibre has total multiplicity \(N\) and the
marked zero has multiplicity one, \(S\le N-1\). Every pole multiplicity
is at most \(N\), so \(M\le N\). This proves \(A(u)\le A_N\), including
\(N=1,2\).

If a contributing \(X_\eta\) had a finite pole of order \(t\ge1\),
(H9) would imply \(a+2\le(a+2)t\le M\). But \(a\ge M-1\) gives
\(a+2>M\), a contradiction. Thus it is polynomial. At each unmarked
zero, the remaining requirement is \(ar_\beta\ge s_\beta\). Because
\(a\ge S\) and \(r_\beta\) is an integer, this is precisely
\(r_\beta\ge1\). These conditions are exactly (H12), proving the
assertion over \(K\) and over its algebraic closure. The sharp example
is given in Section 5. \(\square\)

There is also stabilization before the final jet filtering. The finite
pole cap for a reciprocal unit is exactly one, while the zero part of
\(J_{a,u}\) has disappeared. Its remaining part is the fixed divisor
\[
 J_\infty=\prod_{\beta:\,u(\beta)\ne0,\infty}
                   (x-\beta)^{\operatorname{ord}_\beta u'},
 \tag{H13}
\]
where only positive orders are used and the product is monic.
Equivalently, \(J_\infty\in K[x]\) is defined by those geometric orders.
Consequently every finite candidate list, isolated exponent, and
infinite congruence progression obtained from the fixed-coordinate
classification is the same once \(a\ge A(u)\).

When the canonical polynomial reciprocal generator \(R\) exists, write
\(u=U(R)\). In the global cyclotomic theorem's notation, the stable
pole bounds simplify to
\[
 q_\alpha=\begin{cases}1,&e_\alpha=1,\\0,&e_\alpha>1.\end{cases}
 \tag{H14}
\]
Indeed a quotient pole of order \(h\) pulls back to order \(eh\), and
the allowed order is now one. The twist intervals and the lower exponent
bound therefore stabilize exactly, as do the active set, simultaneous
cyclotomic conditions, and their modulus. The active set consists of values \(R(\beta)\) at finite critical
points \(\beta\) of \(u\) with \(u(\beta)\ne0,\infty\); zero values
of \(U\) are inactive.
This is equality of the full parameter data, not just of their number.

## 3. Reconstruction and exact degree at all later heights

**Corollary StableHeightTransport.** For a fixed \(\eta\in\mathcal E_\infty(u)\)
and odd \(a\ge A(u)\),
\[
              C_{a+2}=X_\eta^2 C_a,\qquad
              D_{a+2}=X_\eta^2 D_a.                  \tag{H15}
\]
This gives a bijection between consecutive stable fibres when the
selected \(X_\eta=\sqrt{-D_a/C_a}\), normalized at zero, is included in
the formula. The inverse recovers the same unit and uses (H5) at the
previous height; its polynomiality follows from stabilization.

Put \(d_\eta=\deg X_\eta\ge0\) and \(s_\infty=\operatorname{ord}_\infty u\).
Then throughout the stable range
\[
 \begin{aligned}
 \deg C_a&=1+a d_\eta+s_\infty,\\
 \deg D_a&=1+(a+2)d_\eta+s_\infty,\\
 \deg_S p_a&=(a+2)(d_\eta+1)+s_\infty.
 \end{aligned}                                                   \tag{H16}
\]

**Proof.** The transport identities follow directly from (H5). Orders
at infinity give the first two degree equalities; the normalized
polynomials are nonzero, so these are their actual degrees. The retained
least sphere-degree formula is
\(\max(\deg C_a+a-1,\deg D_a+a+1)\). The second entry exceeds the first
by \(2(d_\eta+1)>0\), giving (H16). No eventual leading-coefficient
exception is needed: \(X_\eta\) is already a fixed nonzero polynomial.
\(\square\)

## 4. One terminating finite atlas over a number field

**Theorem EffectiveAllOddHeightAtlas.** Over a presented number field,
one terminating computation returns a finite atlas for the complete
fibres at every odd positive height. At most \(\lceil N/2\rceil\)
height entries are required. The last entry applies to every later odd
height, with reconstruction (H5), transport (H15), and degrees (H16).
The atlas consists of explicit isolated units and explicit congruence
families, not a numerical degree truncation. The same computation can
return the geometric atlas over one finite extension, followed by exact
descent tests for the original field.

Here is an explicit construction, specifying the dependencies.

1. Compute the finite zeros and poles with multiplicities, \(A(u)\), the
   marked component, and the bounded canonical-polynomial search of
   Sections 2--3 of the retained effective note. This decides whether
   \(R\) exists and, if so, finds it. These data are independent of the
   height.
2. Work over the common splitting field of the finite support from the
   retained coefficient-field theorem. The degree bound \((2N)!\) is
   uniform in the height. That theorem uses the finite support polynomial
   for the poles and the fibre \(u=-u(\infty)\); it does not require
   splitting every jet polynomial.
3. If \(R\) does not exist, enumerate **once** the finite integer box
   of the effective note's Section 4, using the largest pole bound, at
   height one. Explicitly, if there are \(s\ge1\) support points, one
   may use
   \[
    B=1+\lfloor N/3\rfloor,\qquad H=s!B(2N)^{s-1},
             \qquad -B\le n_i\le H.
    \tag{H17}
   \]
   If \(s=0\), test only \(\eta=1\). All later pole bounds are smaller,
   so this box contains every unit that can occur at any height. Test
   reciprocity exactly on the marked component, then apply (H6)--(H8).
   This returns finitely many units with their complete height intervals.
4. If \(R\) exists, compute \(U\) and enumerate **once** the pole-twist
   box of the effective note's Section 5, again using height one. For
   each of the finitely many odd heights
   \(1,3,\ldots,A(u)\), impose its smaller pole intervals and its lower
   exponent bound. For every surviving twist, solve the full equation
   modulo \(J_{a,u}\), including all nilpotent jets, by
   **ElementaryRationalPointOrbit** from the new elementary exponent
   decision note. Multiplication by \(R\) on the finite quotient,
   restricted to rational scalars, is a rational matrix; the start is
   the bounded-twist representative and the target is \(-1\).
   Every denominator is invertible modulo \(J_{a,u}\), since its roots
   avoid poles of \(u\). The algorithm also handles a nonunit multiplier.
   It returns exactly an empty set, a singleton, or a progression.
   The case \(J=1\) has a vacuous congruence and is handled directly.
5. Attach the exact height to each early entry. Attach all odd
   \(a\ge A(u)\) to the final entry. Use (H5) to reconstruct the
   polynomials at the requested height. Since
   \((A(u)+1)/2\le\lceil N/2\rceil\), this is the claimed finite bound.

The canonical \(R\) descends to the original field. A fixed orbit
\(\eta_0R^j\) descends exactly when \(\eta_0\) does: multiplying by a
known rational function over \(K\) cannot repair failure of descent.
For isolated units, use the same direct coefficient test. The
fixed-coordinate bijection then identifies these with the pairs over
\(K\). This proves the field assertion and termination.

The algorithm uses only the bounded constructions of the older effective
note and the separately proved elementary exponent algorithm. It does
not import the older effective theorem as a whole, its
`RationalOrbitReduction`, or a general Orbit or recurrence-zero theorem.
The accompanying checker does not implement this entire atlas: it
checks the local interval, stabilization, reconstruction, and explicit
families. Universal termination is a written theorem based on the
specified finite constructions.

## 5. Sharp and infinite examples

### 5.1 The degree-only cutoff cannot be decreased

For every integer \(N\ge2\), put \(c=N/(N-1)\) and
\[
 u_N=\frac{x(1-x)^{N-1}}{(1-cx)^N},\qquad
 b_N=(1-x)(1-cx).                                    \tag{H18}
\]
These numerator and denominator are coprime, so \(\deg u_N=N\), and
\(u_N(0)=0,u_N'(0)=1\). Direct logarithmic differentiation gives
\(u_N/(xu_N')=b_N\). The constant reciprocal unit \(\eta=1\) yields
\(X=b_N\) and
\[
 C_a=(1-x)^{a-N+1}(1-cx)^{a+N},\qquad
 D_a=-(1-x)^{a-N+3}(1-cx)^{a+N+2}.                    \tag{H19}
\]
It contributes exactly when \(a\ge N-1\). For \(N\ge3\), the preceding
odd height \(A_N-2\) is positive and excludes this unit, while \(A_N\)
includes it. Thus exact unit-data stabilization cannot start at any
smaller universal odd height. Its stable least degree is \(3(a+2)\).

The polynomial-\(b\) Davenport--Zannier classification already describes
these coordinates more precisely. Only the direct identities above are
needed for sharpness. The new theorem applies also when \(b\) is rational
with finite poles and when there is no polynomial reciprocal generator.

For \(N=1\), \(u=x\) has the stable singleton \(X=1,C=1,D=-1\).
For \(N=2\), (H18) already stabilizes at the first allowed height.
These cases do not require a fictitious height zero or negative height.

### 5.2 An infinite stable fibre

For \(u=x/(1+x)\), put
\[
 R=1+2x,\qquad X_k=(1+x)(1+R^k)/2,\qquad k\ge0.
 \tag{H20}
\]
The degree-one fixed-coordinate theorem gives the complete unit set
\(\eta=R^k\) at every odd height. Thus
\[
 C_{a,k}=(1+x)X_k^a,\qquad D_{a,k}=-(1+x)X_k^{a+2},
 \qquad \deg_S p_{a,k}=(a+2)(k+2).                   \tag{H21}
\]
This is an infinite set of units even though only one height entry is
needed. The atlas does not replace an infinite fibre by a bounded scan.

## 6. Evidence and limits

The main new statements are deductions from the exact local
polynomiality conditions. No new analytic estimate, distribution
hypothesis, or external classification is needed. The closest retained
all-height result is the Davenport--Zannier theorem, whose additional
hypothesis is that \(u/(xu')\) is polynomial. The finite atlas here has
no such hypothesis. The separate uniform coefficient-field theorem
already bounded all coefficients by one extension; that bound alone
did not give stabilization, the exact height intervals, or the sharp
number of height entries.

All claims concern the fixed-coordinate, normalized two-height class
(H2). They do not classify arbitrary sphere phases or assert that the
polynomials themselves are unchanged across heights. No historical
novelty, practical complexity bound, formal verification, or external
review is claimed. Exact local controls with a rational nonpolynomial
\(X\) test polynomiality only unless reciprocity is separately verified;
they are not automatically pure examples.

The source ledger records the exact retained premises and the checker
scope. The receipt binds the proof, checker, ledger and source bytes.

Proof support: [exact checker](verify_all_odd_height_fibres.py),
[source ledger](all_odd_height_fibres_source_ledger.json), and
[verification receipt](all_odd_height_fibres_verification.json).

Run from the repository root:

    python3 research/foundations/extended-geometry/verify_all_odd_height_fibres.py --repo-root . --output research/foundations/extended-geometry/all_odd_height_fibres_verification.json
