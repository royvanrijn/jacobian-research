# Fixed-endpoint radical counts for physical shifts

Fixing one input endpoint makes one letter of the normalized
difference triple a prescribed complete value. The fixed-letter
count then bounds the physical shift by its radical outside a
prescribed prime pool. The pool may be chosen from the first
endpoint or earlier differences. A large radical on the pool
causes no loss, since the powers on that pool are recorded inside
the fixed-letter count.

## A physical difference on one role

Write \(n=X^L\), with \(L>0\) fixed before \(X\to\infty\).
Let \(\mathcal E\) consist of positive primitive triples
\(a+b=c\) with \(X/2<c\le X\), \(\operatorname{rad}(abc)\le c^L\),
and each individual letter radical at most \(n^{1/3+\kappa}\).
Write \(W(x)\) for the complete value in role \(W\) of a row \(x\).
For an integer \(P\ge1\), put
\[
 R_P(m)=\prod_{\substack{p\mid m\\p\nmid P}}p\qquad(m\ne0).
\]
For \(\mathcal B\subseteq\mathcal E\), let \(F_{\mathcal B,W}\)
be the largest number of its rows at one complete \(W\)-value.

Fix \(B,K,\zeta>0\). There exist \(\kappa_0,L_0,\sigma>0\)
such that for \(0<L<L_0\), \(0\le\kappa\le\kappa_0\),
\(0\le b\le B\), each \(x\in\mathcal E\) and \(1\le P\le X^K\),
\[
\begin{split}
 &\#\{y\in\mathcal B:W(y)\ne W(x),\
                   R_P(|W(y)-W(x)|)\le n^b\}\\
 &\quad\ll_{B,K,\zeta,L}
    n^{1/3-\sigma}
     +F_{\mathcal B,W}\,
       n^{\,\min\{b,1/3+\kappa\}/2+\zeta}.
\end{split}                                                     \tag{FE1}
\]
The constants are uniform in the role and the centre. There is no
cap on the raw radical of \(W(y)-W(x)\).

For the proof use the uniform positive-scale consequence of (SH1)
as in (SC5). There are absolute \(a_0,a_1,\gamma>0\), with
\(0<a_0<a_1\), such that outside a uniform count of exponent
strictly below \(1/3\), each letter of \(y\) has normalized
radical mass at least \(\gamma\) on its exact prime powers with
degrees \(Lv_p(W(y))\in[a_0,a_1]\). Choose \(0<\xi<\gamma/2\).
With \(x\) fixed, the rows \(y\) whose \(W\)-radical mass on the
support of \(W(x)\) is at least \(\xi\) also have a uniform
exponent gap, by (EB5). These two classes are absorbed into
the first term of (FE1) by decreasing \(\sigma\), \(\kappa_0\)
and \(L_0\) as needed.
[Positive scaled degrees (SH1)](ABC_SHORT_HEIGHT_AND_REAL_GAP.md#positive-scales-in-every-endpoint-letter)
· [Fixed-support discount (EB5)](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)
· [Physical coordinate differences (SC5)](ABC_PRIMITIVE_SECANT_RADICAL_AMPLIFICATION.md#the-radical-of-a-physical-coordinate-difference)

Put \(w=W(x)\), \(w'=W(y)\), \(h=|w-w'|>0\) and
\(g=(w,w')\). The positive primitive triple with values
\(w/g,w'/g,h/g\), ordered with the larger of \(w/g,w'/g\)
last, has height at least \(X^{a_0(\gamma-\xi)}\). In fact
the indicated band of \(w'\) outside the support of \(w\)
survives in \(w'/g\) in full. First choose \(g\mid w\)
at height-subpower cost. The complete value \(w/g\) is now
fixed. Its raw radical and the raw radical of \(w'/g\)
are at most \(n^{1/3+\kappa}\). The outside-\(P\) radical of
\(h/g\) is at most \(n^b\). Apply (FL1) with the pool on
\(h/g\), the empty pool on \(w'/g\), and the appropriate
role for \(w/g\). The possible role and summand ordering
give only a bounded factor. For each output, \(g\) and the
ordered branch recover \(w'\), which has at most
\(F_{\mathcal B,W}\) input preimages.
[Fixed complete letters and prescribed pools (FL1)](ABC_FIXED_LETTER_RADICAL_CAPS.md#a-pointwise-count-with-prescribed-prime-pools)

Here and below, the shell change in (FL1) is made as follows.
If the derived heights are in \([X^{\omega},X]\), with
\(\omega>0\) fixed, split into dyadic shells of upper endpoint
\(Y\). Use the fixed budget unit \(L\) in the application of
(FL1) at height \(Y\). A normalized original cap \(e\) in base
\(n\) becomes the cap \(\ell=e\log X/\log Y\) in base \(Y^L\).
The factor \(\log X/\log Y\) is bounded in terms of \(\omega\).
The same comparison keeps the pool height polynomial in \(Y\).
The term \(\ell L\log Y\) equals \(eL\log X\) exactly. The
normalized allowance in (FL1), the divisor loss and the shell
sum fit inside \(\zeta\) by taking their inner allowances smaller
than fixed multiples of \(\zeta L\). Thus (FL1) gives the
claimed uniform estimate with the original budget fixed before
the height limit. Its two free-letter conditions concern outside
radicals, which is why a large raw radical of \(h/g\) is allowed.

## A uniform bound inside any fixed complete fibre

Let \(\mathcal E_T\) have one fixed complete value in any role
\(T\). For distinct \(x,y\in\mathcal E_T\), write \(h(x,y)\)
for the absolute difference in either other letter. The
magnitudes are equal because the remaining two signed changes
are opposite. With the parameters of (FE1), uniformly in the
fibre, its centre \(x\), and \(1\le P\le X^K\),
\[
 \#\{y\in\mathcal E_T\setminus\{x\}:
                        R_P(h(x,y))\le n^b\}
       \ll_{B,K,\zeta,L} n^{\,\min\{b,1/3+\kappa\}/2+\zeta}.
                                                                    \tag{FE2}
\]
This estimate has no exceptional set of input rows.

In the two moving roles write the positive values as \(v,v+h\)
and \(w,w+h\). Their normalization factors are
\(g_v=(v,h)\) and \(g_w=(w,h)\). They are coprime by
primitivity of the original rows and \(g_vg_w\mid h\).
Both larger values exceed \(h\), and at least one exceeds
\(X/4\). For the latter, the target and one summand of each
original row both exceed \(X/4\), so any two distinct
roles include one of these. Therefore
\[
 \frac{v+h}{g_v}\frac{w+h}{g_w}
       \ge \frac{(v+h)(w+h)}{h}\ge X/4.
\]
Choose a moving role in which the larger normalized value is
at least \(\sqrt X/2\). Once the role and \(g\) have been
chosen, the corresponding endpoint value at \(x\), divided
by \(g\), is fixed. There are divisor-many choices of \(g\)
because \(g\) divides that original value at \(x\). Use (FL1)
with this fixed value, the other endpoint quotient with its
original radical cap, and \(h/g\) with its outside-\(P\) cap.
The height comparison above applies, for example with
\(\omega=1/3\) for large \(X\). The recovered shift and its sign
on the selected role, together with \(x\) and the fixed \(T\),
determine \(y\). The constant accounts for the two signs.
[Tall primitive difference branches](ABC_POPULAR_SHIFT_RADICAL_AMPLIFICATION.md#a-tall-primitive-difference-triple)

## Dense pair relations and target fibres

If \(\mathcal A,\mathcal B\subseteq\mathcal E\) both have
exponent \(1/3-o(1)\) in a height-first endpoint sequence,
suppose that pairs of distinct \(W\)-values with
\(R_{P(x)}(|W(x)-W(y)|)\le n^b\) have exponent \(2/3-o(1)\)
for bounded \(b\). The pool \(P(x)\le X^K\) may vary with \(x\).
Then the complete \(W\)-dictionary on the \(\mathcal B\) side
has a subfamily at exponent \(1/3-o(1)\) with dictionary
exponent \(1/6+o(1)\), and \(b\ge1/3-o(1)\).
To prove this dyadically bin the \(\mathcal B\) fibre sizes
over \(W\), at subpower cost, and keep a full-exponent pair
bin. If \(f\) is the retained fibre exponent, (FL1) on the
original rows gives \(f\le1/6+o(1)\). Sum (FE1) over \(x\).
The bound forces \(f+\min(b,1/3)/2\ge1/3-o(1)\).
The projection of the retained pair relation to \(\mathcal B\)
has exponent \(1/3-o(1)\). Keep the entire bin on the
\(\mathcal B\) side to preserve the lower fibre sizes.
Its size is at most \(n^{1/3+o(1)}\), so the dictionary
has exponent at most \(1/3-f+o(1)=1/6+o(1)\). Its lower
bound follows from (FL1). The same argument gives the
conclusion on the \(\mathcal A\) side if the test also holds
with pools assigned to the second endpoint. In particular
a fixed raw radical cap gives both conclusions by taking
\(P=1\). Reversal alone does not reverse an adaptive pool.

More quantitatively, for each fixed \(F<1/6\), if
\(F_{\mathcal B,W}\le n^F\) and both families are sufficiently
close to full exponent, then for every fixed \(B,K>0\) there is
a uniform power saving in the number of pairs with
\(R_{P(x)}(|W(x)-W(y)|)\le n^B\). Choose the balance tolerance
and exponent allowances in (FE1) small compared with
\(1/6-F\). Equal-value pairs have the separate bound
\(|\mathcal A|F_{\mathcal B,W}\). In the central least-target-repair
normalization \(Z=Y^{1/M}=n^{1/3+o(1)}\), a multiplicity cap
\(Z^{u+o(1)}\) with fixed \(u<1/2\) for each occurring target
radical implies this condition for complete target values.
Thus near-full input families in this range have, on all but a
power-small fraction of pairs, outside radical exceeding \(n^B\)
on their physical target differences. The discount may use any
pool of polynomial height chosen using the first input. This
statement concerns primitive quotients of the least repairs.
[Repaired target multiplicity (FL5)](ABC_FIXED_LETTER_RADICAL_CAPS.md#size-of-a-minimally-repaired-target-radical-fibre)
· [Least-repair equivalence](POWERFUL_TRIPLES_ABC_ENDPOINT_TRANSFER.md)

The conditioning in (FE2) can be iterated inside any fibre
\(|\mathcal E_T|\ge n^r\) with \(r>0\). Fix an integer \(k\ge1\)
and \(0<\eta<2r\). For a given centre \(x\) and pool \(P_0\) of
polynomial height, all but a power-small fraction of ordered
\(k\)-tuples \((y_1,\ldots,y_k)\) in the fibre satisfy
\[
 \log_n
 R_{P_0\,abc(x)\prod_{j<i}h(x,y_j)}(h(x,y_i))>2r-\eta
 \qquad(1\le i\le k).                                             \tag{FE3}
\]
The bound is uniform in the centre and pool with a fixed height
exponent. Condition on the earlier endpoints and use (FE2) with
\(b=2r-\eta\) and \(\zeta<\eta/4\). The pool at stage \(i\) is
fixed independently of \(y_i\) and has polynomial height since
\(k\) is fixed. Terms equal to the centre or previous endpoints
are negligible. The \(k\) counted prime sets are disjoint, so
the union of the shift supports outside \(P_0\,abc(x)\) has
normalized radical mass greater than \(k(2r-\eta)\).
In particular a fibre of exponent \(1/6-o(1)\) cannot have full
pair exponent on a branch whose shift has outside radical below
\(n^{1/3-\eta}\), even with a pool chosen from the first row.

For a near-full multiplicity cohort in the central least-repair
normalization, a target-radical fibre of size \(Z^{u-o(1)}\)
has a complete-target subfibre of size \(n^{u/3-o(1)}\):
the fixed radical has height-subpower many complete targets
in the original shell. For \(u>0\), (FE3) applies to such
subfibres with \(r\) arbitrarily close to \(u/3\) from below.
The argument can also be summed over all complete-target
subfibres of exponent at least \(r\) with the same uniform
constants. Thus the fixed-target and variable-target tests
yield different constraints in the intermediate multiplicity
regime. Neither gives an upper bound for the occurrence of
those fibres by itself.

## Transfer at maximum multiplicity

There is a separate consequence of an upper bound on the
*raw* shift radical. Take a height-first balanced endpoint sequence
and an ordered common-\(T\) pair relation of exponent \(1/2-o(1)\)
with \(\operatorname{rad}(h)\le n^B\) for some fixed \(B\).
After retaining a full-exponent pair relation, both input sides
have projected exponent \(1/3-o(1)\) and all three complete-letter
dictionaries have exponent \(1/6+o(1)\) in ambient subfamilies
containing those projections. This transfer does not assume
minimum nonatomic input profiles.

For either moving role \(W\), fix \(x\) and let
\(g_W=(W(x),h)\). If \(\log_n\operatorname{rad}(g_W)\ge\xi>0\),
the second row \(y\) has at least \(\xi\) normalized radical
on the fixed pool of \(W(x)\) in that role. Apply (FL1) on
the **original** fixed \(T\), then sum over \(x\). These pairs
have exponent at most \(1/2-\xi/2+o(1)\). Thus on a retained
full-exponent pair relation the two \(g_W\) have radical
\(n^{o(1)}\). Discard rows failing the uniform compact-band
conclusion of (SH1). Their number is power smaller than
\(n^{1/3}\), and each has at most \(n^{1/6+o(1)}\) partners
by (FL1). The positive band away from the primes of \(g_W\)
then makes both normalization branches tall.

For one moving role the normalized shift triple has raw
radical at most \(n^{2/3+B+o(1)}\), and its two endpoint
quotients have individual caps \(n^{1/3+o(1)}\). The two-letter
bound (GD3), with singleton root groups, counts at most
\(n^{1/3+o(1)}\) such output triples. Apply it in dyadic
height shells with cap conversion as above. The fixed raw
bound, including the third output, keeps the relative
anatomic budget \(O_B(L)\). Record \(g_W\) at cost \(n^{o(1)}\)
by radical enumeration and fixed-support counting. For
each output and its multiplier the two original \(W\)-values
and \(h\) are fixed, up to the bounded orientations.
Completions are common \(T\)-values among the two endpoint
families. For each fixed \(T\) there is at most one completion
in each family because \(T\) and \(W\) determine the third
role. Thus the multiplicity is bounded by the smaller of
the two original \(W\)-fibre caps.
[Two-letter dictionary bound (GD3)](ABC_GROUPED_ROOT_DICTIONARY_ENTROPY.md#group-dictionary-count)

Bin all six original fibre sizes at height-subpower cost.
For either moving role the preceding count forces both
retained fibre exponents to be at least \(1/6-o(1)\).
For \(T\) this follows directly from the pair exponent and
the bound \(n^{1/3+o(1)}\) on each input family. All six
upper fibre bounds follow from (FL1). Each side projects
onto exponent \(1/3-o(1)\) since each row has at most
\(n^{1/6+o(1)}\) opposite partners. Its three dictionary
caps follow from the lower fibre-size thresholds in the
ambient family used to form the bins, even if the bin
intersection decreases individual fibres. This gives the
claimed transfer. It leaves open whether a full-exponent
relation of bounded-radical shifts must occur.

## Relation to existing counting methods

The numerical height in (FE1) comes from a positive
scaled-degree band in an original letter. In (FE2) it comes
from the two moving branches together. Once a branch is
tall, fixing one of its complete values uses a different
estimate from the two-free-letter secant bound (SC5). It
permits pools that vary with the fixed endpoint and places
no restriction on the shift radical inside a pool. The
unrestricted exceptional-set estimates of Bernert and of
Bernert, Browning, Lichtman and Teräväinen give global
bounds in terms of the original radical budget. They do
not give this fixed-endpoint estimate for physical
differences or the abundance of compressed shifts that
would be needed for a strict unrestricted improvement.
[The exceptional set in the abc conjecture, Theorem 1](https://arxiv.org/html/2506.13364#Thmtheorem1)
· [Bounds on the exceptional set in the abc conjecture, Theorem 1.3](https://arxiv.org/html/2410.12234v2#Thmtheorem1.3)
