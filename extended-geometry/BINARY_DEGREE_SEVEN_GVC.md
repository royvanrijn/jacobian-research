# Binary GVC through polynomial degree seven

## 1. Theorems

Work over a characteristic-zero field.  Let
\[
 \Lambda=\Lambda_r+\Lambda_{r+1}+\cdots
\]
be a constant-coefficient operator in two variables whose lowest nonzero
positive homogeneous part has order \(r\), and let \(\deg P=7\).

> **Theorem 1.1 (high-order septic rows).**
> If \(4\leq r\leq6\) and
> \[
>  \Lambda^m(P^m)=0\qquad(m\geq1),
> \tag{1.1}
> \]
> then, for every fixed polynomial \(Q\),
> \[
>  \Lambda^m(QP^m)=0\qquad(m\gg0).
> \tag{1.2}
> \]

Together with the already proved quadratic- and cubic-leading septic rows,
the lowest-order-one theorem, and the theorem for \(\deg P\leq r\), this
gives the complete degree statement.

> **Corollary 1.2 (binary GVC through degree seven).**
> Every constant-coefficient operator in two variables satisfies the GVC
> conclusion for every polynomial of degree at most seven.

This degree theorem did not by itself prove unrestricted
\(\operatorname {GVC}(2)\).  The later
[Hall-envelope theorem](BINARY_GVC_ENVELOPE_CLOSURE.md) supplies the
all-degree proof, so the septic computation is retained as an independent
exact regression rather than a current counterexample frontier.

The proof is a generated finite Newton-state calculation followed by an
all-order terminal-face argument.  The calculation does more than report
that no example was found: it identifies the complete reduced geometry of
every non-origin face and proves that its surviving components cannot
recur.

## 2. Uniform Hall charts

Fix \(r\in\{4,5,6\}\).  Factor \(\Lambda_r\) and \(P_7\) after scalar
extension.  At a marked operator root of multiplicity \(e\), Hall failure
is equivalent to at least \(8-e\) polynomial factors annihilating that
direction.  A deficient set cannot contain two nonparallel derivative
directions.  Hence every leading pure zero lies on
\[
 M_e=X^eY^{r-e},
 \qquad
 P_7=y^{8-e}C_{e-1}(x,y),
 \qquad 1\leq e\leq r.
\tag{2.1}
\]

Let \(t\) be the largest \(x\)-degree appearing in the nonzero form
\(C_{e-1}\).  Then
\[
 0\leq t\leq e-1,
 \qquad
 [x^ty^{7-t}]P_7\ne0.
\tag{2.2}
\]
After scaling, the marked polynomial monomial is \(x^ty^{7-t}\).  Every
other degree-seven monomial has smaller \(x\)-degree and is strict for all
weights used below.  Thus the pairs \((e,t)\), not a classification of the
remaining roots of \(C_{e-1}\), give all projective Hall charts.  There are
\[
 \sum_{e=1}^r e=10,15,21
\tag{2.3}
\]
charts at orders four, five, and six.

Formal division by \(M_e\) gives the exact two-wing support rule.  At order
\(r+h\), \(h\geq1\), a normalized operator monomial can remain only when
\[
 X^aY^{r+h-a},
 \qquad
 a<e\quad\text{or}\quad a>e+h.
\tag{2.4}
\]
At order \(r\), the marked monomial has \(a=e\), while its homogeneous
cofactor has \(a>e\).  The divided differential unit is locally invertible
on polynomials and preserves (1.1)--(1.2).

The checker exhausts Hall matching for every integer partition of
\(r=4,5,6\), all distributions of seven annihilator factors, and every
chart (2.2).

## 3. Newton states and the no-cycle invariant

Use weights
\[
 w_s(x)=w_s(X)=s,
 \qquad
 w_s(y)=w_s(Y)=1,
 \qquad s>1.
\tag{3.1}
\]
The initial operator minimum and polynomial maximum are
\[
 U=X^eY^{r-e},
 \qquad
 V=x^ty^{7-t}.
\tag{3.2}
\]
Their common threshold is
\[
 s_*=1+\frac{7-r}{e-t}.
\tag{3.3}
\]

More generally, a state consists of a marked minimum
\(U=X^aY^b\), a marked maximum \(V=x^iy^j\), and a current slope interval.
Only an operator term with \(X\)-exponent \(a'<a\) can cross below \(U\)
as the slope increases.  Only a polynomial term with \(x\)-degree \(i'>i\)
can cross above \(V\).  Their crossing slopes are
\[
 s=\frac{b'-b}{a-a'},
 \qquad
 s=\frac{j-j'}{i'-i}.
\tag{3.4}
\]
The two-wing rule bounds every possible \(a'\) and \(b'\) in a bounded
interval, while \(P\) has finite support.  Hence every state has a complete
finite crossing list.

At one crossing let \(A_s\) and \(B_s\) be the complete equality faces.
The maximal output-weight part of (1.1) gives
\[
 A_s(\partial)^m(B_s^m)=0\qquad(m\geq1).
\tag{3.5}
\]
The exact radicals in Section 4 use the coefficients of (3.5) for
\(1\leq m\leq10\).

The only non-origin radicals have the form
\[
 \sqrt I=(\text{all other face coordinates},uv),
\tag{3.6}
\]
where \(u\) is one newly tied operator coordinate and \(v\) is one newly
tied polynomial coordinate.  Thus both sides cannot migrate at the same
wall.  On the \(u\)-axis, the new operator exponent satisfies \(a'<a\); on
the \(v\)-axis, the new polynomial exponent satisfies \(i'>i\).  In either
case the positive integer
\[
 g=a-i
\tag{3.7}
\]
strictly decreases.  This is the promised no-cycle invariant.  A surviving
face is a one-sided Newton pivot, not a cancellation family, and a chain of
such pivots cannot recur indefinitely.

## 4. Exact face classification

The complete census is
\[
\begin{array}{c|r|r|r|r|r}
r&\text{Hall charts}&\text{initial faces}&\text{axis faces}
 &\text{child states}&\text{child faces}\\ \hline
4&10&97&8&16&58\\
5&15&112&7&14&40\\
6&21&78&0&0&0.
\end{array}
\tag{4.1}
\]
Every unlisted initial face and every one of the 98 child faces has radical
equal to the face-coordinate origin.  The fifteen axis faces are listed
below.  In each row, \(u\) and \(v\) denote the coefficients of the two
displayed monomials, and the radical is exactly (3.6).

\[
\begin{array}{c|c|c|c|c|c}
r&e&t&s&u\text{-monomial}&v\text{-monomial}\\ \hline
4&2&0&2&XY^4&xy^5\\
4&3&0&\frac32&XY^4&x^2y^4\\
4&3&1&2&X^2Y^3&x^2y^4\\
4&4&0&\frac43&XY^4&x^3y^3\\
4&4&0&\frac32&X^2Y^3&x^2y^4\\
4&4&0&\frac53&XY^5&x^3y^2\\
4&4&1&\frac32&X^2Y^3&x^3y^3\\
4&4&2&2&X^3Y^2&x^3y^3\\ \hline
5&3&0&\frac32&XY^5&x^2y^4\\
5&4&0&\frac43&XY^5&x^3y^3\\
5&4&1&\frac32&X^2Y^4&x^3y^3\\
5&5&0&\frac54&XY^5&x^4y^2\\
5&5&0&\frac43&X^2Y^4&x^3y^3\\
5&5&1&\frac43&X^2Y^4&x^4y^2\\
5&5&2&\frac32&X^3Y^3&x^4y^2.
\end{array}
\tag{4.2}
\]

For example, the first row means that the complete slope-two face radical
contains every face coordinate except
\(u=[XY^4]\Lambda\) and \(v=[xy^5]P\), and contains \(uv\).  Its two
components pivot to \((XY^4,y^7)\) and \((X^2Y^2,xy^5)\).  The checker
generates both child intervals from the radical, verifies that the gap
(3.7) decreases, derives every later crossing from (2.4), and proves that
every later radical is the origin.  The other fourteen rows are treated
identically.

There is therefore no surviving nonlinear component, no mixed semigroup
return, and no adelic collision on the high-order septic rows.  The only
objects which “survive” an initial radical are one-sided changes of the
marked Newton endpoint.

## 5. All-order termination

If an initial face radical is the origin, the marked pair advances to its
common threshold (3.3).  If it is one of (4.2), each of its two axes enters
the generated child interval, all of whose crossings have origin radical.
Every branch therefore reaches a slope \(s_*\ne1\) and a weight \(W\) with
\[
 w_{s_*}(\operatorname {Supp}\Lambda)\geq W,
 \qquad
 w_{s_*}(\operatorname {Supp}P)\leq W.
\tag{5.1}
\]

The unequal-weight terminal-face theorem applies to the equality faces in
(5.1), giving a coordinate-derivative deficit linear in \(m\).  Every
strict operator or polynomial selection consumes positive integral weight
defect.  A fixed multiplier \(Q\) permits only \(O_Q(1)\) strict selections,
which cannot repair the linear deficit.  Hence (1.2) holds on every branch.
This proves Theorem 1.1.

For Corollary 1.2, degrees at most six were already closed.  In degree
seven, lowest positive order one is safe; orders two and three are the
[quadratic-leading](BINARY_QUADRATIC_SEPTIC_GVC.md) and
[cubic-leading](BINARY_CUBIC_SEPTIC_GVC.md) septic theorems; orders four
through six are Theorem 1.1; and order at least seven is covered by
\(\deg P\leq r\).  These cases are exhaustive.

## 6. Reproduction

Run

```bash
.venv/bin/python scripts/verify_binary_high_order_septic_gvc.py
```

The command requires SymPy and Singular.  It exhausts Hall matching for all
root partitions of orders four through six, generates every chart and
crossing, verifies all 385 exact characteristic-zero face radicals, checks
the fifteen squarefree axis ideals and their 30 primaries, proves strict
decrease of the gap (3.7), and audits every final common threshold.  The
all-order passage from those thresholds to (1.2) is the written
terminal-face argument, not a bounded extrapolation.
