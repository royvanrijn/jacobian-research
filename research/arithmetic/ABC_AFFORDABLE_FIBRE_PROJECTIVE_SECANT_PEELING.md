# Real and projective peeling for secants in an affordable fibre

A large fixed-record fibre admits simultaneous thinning on the
real projective line and modulo prime powers. After the thinning,
full prime powers in several different secants can be combined
as congruences for a single query. This uses a lower bound on
the prime valuations in those packets. Their valuations need
not lie on a common divisibility grid.

## Fixed-record setup

Fix \(L>0\) before the height limit, and write \(n=X^L\).
Consider a family \(\mathcal A\) of actual primitive triples
with distinct targets at height at most \(X\), in one
affordable record
\[
 U_i=F u_i,\qquad V_i=G v_i,\qquad
 T_i=\sigma F u_i+\tau G v_i=D w_i.
\]
The positive integers \(F,G,D\) are fixed full factors,
\(D=1\) is allowed, and the signs are fixed. Refine to dyadic
cofactor intervals \(u_i\in[U,2U]\), \(v_i\in[V,2V]\).
Suppose the resulting family has \(N\ge n^r\) members for a
fixed \(r>0\). For an anchor \(i\) and a distinct query \(j\)
write
\[
 d_i(j)=|u_i v_j-u_j v_i|,\qquad
 H_{ij}=\frac{d_i(j)}
                {D\gcd(u_i,u_j)\gcd(v_i,v_j)\gcd(w_i,w_j)}.
\]
All secants are nonzero because the pairs \((u_i,v_i)\)
are primitive and distinct. The quotient is integral by
(DC1) and fullness of \(D\).
[Quotient secants and common target factors](ABC_AFFORDABLE_FIBRE_DETERMINANT_COMPLEXITY.md#three-orientations-of-the-same-quotient-determinant)

Let \(\theta_j\) be the angle of \((u_j/U,v_j/V)\).
These vectors lie in the fixed positive compact box
\([1,2]^2\), so for all anchors and queries
\[
             d_i(j)\asymp UV\,|\theta_i-\theta_j|,            \tag{PS1}
\]
with absolute constants, uniformly also at very small
separations. This follows from the determinant formula
using the two normalized vector lengths and the sine of
their angle difference.

## Simultaneous peeling

Fix \(0<\gamma<r/10\), a numerical tolerance \(\eta>0\),
and a polynomial support bound \(C\ge3\). There is a
\(\delta>0\), depending only on \(\gamma,\eta\), such that
one may retain a family \(\mathcal A_*\) of size
\(N_*\ge N n^{-\gamma}\) with the following properties.
Its angular span is \(A>0\), and all its points are in one
projective residue class modulo an integer \(M\) divisible by
\(D\). Put
\[
                          h=\log_X(UV A/M).
\]
For every fixed anchor from \(\mathcal A\), the number of
queries in \(\mathcal A_*\) with
\[
                       |\theta_i-\theta_j|<A X^{-\eta}
\]
is at most \(N_* n^{-\delta}\). For every fixed anchor and
every fixed prime set \(S\) with \(\prod_{p\in S}p\le X^C\),
at most \(N_* n^{-\delta}\) distinct queries satisfy
\[
          \sum_{p\in S}\bigl(v_p(d_i(j))-v_p(M)\bigr)_+
                                  \log_X p\ge\eta.           \tag{PS2}
\]
Here \(a_+=\max(a,0)\). Both bounds are simultaneous
statements for the chosen family. The height threshold for
the subpower enumerations is taken at fixed \(L,C,\gamma,\eta\).

Start with \(M=D\), since \(D\) divides all secants. For
primitive vectors, a fixed projective class modulo \(M\)
is equivalent to pairwise vanishing of the determinant
modulo \(M\). Locally at a prime power choose a unit
coordinate of the reference vector. The vanishing
determinant makes the other vector a scalar multiple,
and the scalar is a unit by primitivity. The Chinese
remainder theorem combines the local classes.

At a current stage, let \(A\) be the actual angular span.
If an anchor's angular neighborhood at the indicated scale
contains more than the prescribed fraction, retain those
queries. Their new span is at most \(2A X^{-\eta}\).
If (PS2) fails for some anchor and pool, retain the large
violating class and fix the exact \(S\)-supported part of the
nonzero secant \(d_i(j)\). There are \(X^{o_X(1)}\) possibilities
by the uniform fixed-support estimate, since \(d_i(j)\) has
polynomial height. Call the fixed supported value \(d_S\).
The primitive anchor prescribes the projective class
modulo \(d_S\) by the same unit-coordinate argument.
The new common modulus can be taken to be
\(\operatorname{lcm}(M,d_S)\). Its numerical logarithm
increases by at least \(\eta\) by (PS2). The use of the
lcm is essential when some primes of \(S\) already divide \(M\).
[Uniform fixed-support enumeration](ABC_EXCEPTIONAL_SET_AND_POWERFUL_TRIPLES.md#a-radical-count-contradiction-to-the-claimed-smooth-number-lower-bound)

Each step raises \(\log_X(M/A)\) by at least \(\eta/2\)
for large \(X\). As long as two points remain, (PS1) and the
common modulus give \(M/A\ll UV\le X^2\).
Initially the potential is nonnegative, since the angular
span in the normalized box is less than one. Thus there
are at most \(J=\lceil 12/\eta\rceil+3\) steps.
For example choose \(\delta=\gamma/(10J)\).
Each step retains at least a fraction \(n^{-2\delta}\)
for large \(X\), absorbing the within-pool subpower
selection at fixed \(L\). Hence the size stays at least
\(N n^{-\gamma}\) through the possible steps, in particular
greater than one. The procedure stops with both simultaneous
properties. No collection of possible pools or anchors is
enumerated in this existential thinning.

Projective packing in the terminal angular interval gives
\[
                N_*\ll 1+UV A/M,
         \qquad h\ge (r-\gamma)L-o_X(1).                     \tag{PS3}
\]
Indeed the spacing of distinct angles in a common class
modulo \(M\) is \(\gg M/(UV)\) by (PS1).
For a fixed \(L\), one may choose \(\eta\) small
compared with the lower bound in (PS3).
In the iterated balanced endpoint setting there is a
stronger choice: the uniform packing gap (with projective
modulus \(M/D\) after dividing secants by \(D\)) gives
\(h\ge\nu>0\) for all sufficiently late budgets, with
an arbitrary fixed smaller margin absorbed in \(\nu\).
Indeed \(|\Delta_{ij}|/D\ll UV A/D\), and the retained
class still has exponent at least \(r-\gamma\). This allows
\(\eta\) to be an absolute small constant with respect to
\(L\), after the positive exponent and affordable margin
are fixed.
[Uniform projective modulus gap](ABC_AFFORDABLE_FIBRE_SECANT_VALUATION_WEIGHTS.md#setup-and-the-projective-packing-gap)

## Counting overlapping families of full secant divisors

Fix a real \(0<\vartheta<1\) and integer \(m\ge2\) with
\(m\vartheta>1\). Use a support bound \(C\) large enough
for the union of the supports of \(M\), the \(m\)
anchor cofactor triples, and their pairwise determinants.
For instance \(C>3+3m+2m^2\) suffices for large \(X\).
Assume
\[
                         m(\vartheta h-\eta)>h              \tag{PS4}
\]
with fixed slack at the chosen budget.
The lower bounds for \(h\) above show how to achieve
this by the choice of \(\eta\).

Let \(k\ge2\) be an integer and put \(\kappa=kL\).
For each pair \((i,j)\) take the full part of \(H_{ij}\)
on primes outside \(M\) with \(v_p(H_{ij})\ge k\).
One may also impose the cutoff \(\log_n p\le\tau\) for
some fixed \(\tau>0\).
Call the pair eligible when the numerical logarithm
of this part is at least \(\vartheta h\).
If
\[
 \begin{cases}
       h/\kappa+\tau\le r-2\gamma, &\text{with the cutoff},\\
       2h/\kappa\le r-2\gamma,     &\text{without the cutoff},
 \end{cases}                                                     \tag{PS5}
\]
then the number of eligible ordered pairs in
\(\mathcal A_*\) is
\[
                              O(N_*^2 n^{-\delta/(2m)}).
                                                               \tag{PS6}
\]
An analogous statement applies to any chosen full divisor
of \(H_{ij}\) outside \(M\) of numerical weight at least
\(\vartheta h\) and radical at most \(n^b\), if
\(m b\le r-2\gamma\). In that version no lower bound
on each individual selected valuation is needed.

Indeed, otherwise count ordered \(m\)-tuples of distinct
anchors incident to each query in the eligible-pair graph.
The ordered star count and the power mean give \(m\)
anchors with \(\gg_m N_* n^{-\delta/2}\) common queries.
The average degree tends to infinity by the lower bound
on \(N_*\), so the falling-factorial counts of distinct
anchors give this same estimate. Use the pool \(S\)
consisting of all primes of \(M\), of all the chosen anchor
\(U_iV_iT_i\), and of every nonzero determinant of two
anchors. All but at most \(mN_* n^{-\delta}\) of the
common queries satisfy the terminal bound (PS2) for
these anchors and this pool.

Discard primes in \(S\) from each selected full part.
For primes of \(S\) outside \(M\), the entire valuation
is an increment in (PS2), so this discards numerical
weight less than \(\eta\) per anchor. No primes of \(M\)
were used in the selected parts in the first place.
Outside \(S\) the valuations of the raw and quotient
secants coincide. Also a prime dividing the two raw
secants of different anchors must divide their determinant:
eliminating the query coordinates shows that the prime
divides that determinant times both coordinates, and the
query is primitive. Thus the remaining prime powers from
all \(m\) anchors have pairwise disjoint supports. Their
total numerical weight exceeds \(h\) by (PS4).

With the full-part condition, select prime powers from
these parts up to the first time their accumulated
numerical weight is at least \(h\). Before the last power
the normalized radical cost is at most \(h/\kappa\),
since each selected order is at least \(k\).
The last prime costs at most \(\tau\) if the cutoff is
imposed. Without the cutoff its cost is at most
\(h/\kappa+o_X(1)\): a prime power outside \(M\) in a
single secant divides \(d_i(j)/M\ll X^h\) by (PS1).
In the analogous arbitrary-divisor version, simply take
the product of the \(m\) stripped divisors, at cost at
most \(m b\). In either case the product of chosen moduli
\(Q\) is coprime to \(M\), and \(\log_X Q\ge h\).

Record the radical of \(Q\), its assignment to the \(m\)
anchors, and the exact powers. The latter two contributions
are \(X^{o_X(1)}\), uniformly at fixed \(L,m\) and the
polynomial height of the moduli. With (PS5) the exponent
for all records is at most \(r-2\gamma+o_X(1)\) in base
\(n\). For each record the primitive query is in one
projective class modulo \(M Q\), because each selected
modulus comes from the secant at its prescribed anchor.
Packing in the common angular interval gives
\(O(1+UV A/(M Q))=O(1)\) points. This contradicts
the lower exponent \(r-\gamma-\delta/2\) for the common
queries. This proves (PS6).

For \(m=2\), any \(\vartheta>1/2\) works once the peeling
tolerance is small enough. In particular a thick relation
consisting of \(k\)-full quotient secants with more than this
fraction of effective weight on new primes is excluded at
the cost (PS5). Fullness alone need not supply that much
weight outside the fixed modulus. The weight condition in
the eligible relation is the additional precise requirement.

## A hybrid with full factors on the original letters

Suppose a query has three selected full packets on its
actual residual \(u_j,v_j,w_j\), with product \(P\) coprime
to \(FGDM\), numerical weight at least \(a\) and normalized
radical cost at most \(b\). The packets are pairwise coprime
by primitivity. A packet on \(u_j\) or \(v_j\) specifies the
projective coordinate zero modulo that packet. A packet on
\(w_j\) specifies \(\sigma F u_j+\tau G v_j=0\) modulo the
packet, with \(F,G\) units. In each case the other primitive
coordinate is a unit in the indicated chart, so the packet
gives one projective class.

A prime of \(P\) that also divides an anchor secant belongs
to the support of the corresponding anchor letter. Indeed
the prime is outside \(FGD\) and identifies the two primitive
projective points. Any one of \(u_j,v_j,w_j\) that vanishes
there also vanishes at the anchor, using the signed linear
form for \(w_j\). The selected original packets are therefore
coprime to all secant packets after stripping the anchor
pool \(S\).

If \(a\ge h\), the number of queries admitting such original
packets is at most \(n^{b+o_X(1)}\), with no secant condition.
Enumerate their joint radical, the assignments and the
powers. For each record \(MP\) already gives packing
multiplicity \(O(1)\) in the terminal angular interval.

For \(0\le a<h\), count only eligible pairs whose query
admits these original packets. In place of (PS4) assume
\[
                        m(\vartheta h-\eta)+a>h .
\]
If secants are needed, the root-cost hypotheses replacing
(PS5) are
\[
 \begin{cases}
      b+(h-a)/\kappa+\tau\le r-2\gamma, &\text{with the cutoff},\\
      b+((h-a)+h)/\kappa\le r-2\gamma, &\text{without the cutoff}.
 \end{cases}
\]
Their conclusion is the same \(O(N_*^2 n^{-\delta/(2m)})\)
bound as in (PS6). Indeed the same star count and the
same pool stripping leave \(\gg_m N_* n^{-\delta/2}\)
queries if the bound fails. For each such query choose its
original packets first. Then add the available secant prime
powers in a fixed prime order until \(P Q\) has numerical
weight at least \(h\), or use \(Q=1\) if the original
packets already suffice. Each next prime is outside the
modulus already selected, by the disjointness just proved
and the removal of the pairwise anchor differences.
The available weights ensure that the selection reaches
the threshold. Before the final prime the secant radical
cost is at most \((h-a)/\kappa\). The last prime has cost
at most \(\tau\), or \(h/\kappa+o_X(1)\) without the cutoff,
as before. This accounts for the exact chosen packets
without counting a free placement or origin for the query.
Enumerating the combined record and packing modulo \(MPQ\)
contradicts the same common-query bound.


## A consequence without any secant power hypothesis

For each anchor \(i\), apply (PS2) with the fixed pool
\(S_i\) of primes dividing \(M U_iV_iT_i\).
For all but an \(O(n^{-\delta})\) fraction of pairs,
(PS1) and the terminal real separation give
\(\log_X(d_i(j)/M)\ge h-\eta-o_X(1)\).
The part of \(d_i(j)/M\) supported on \(S_i\) has
numerical weight less than \(\eta\) by (PS2).
Off \(S_i\) the remaining value equals the corresponding
part of \(H_{ij}\), and its primes also avoid the complete
query triple by the same projective residue argument as
in the hybrid count above.

Take \(\kappa=kL\) uniformly large enough for the
no-cutoff case of (PS5). For instance one may choose
\(k=\lceil\kappa_0/L\rceil\) with
\(\kappa_0>2h_0/(r-2\gamma)\), given any uniform upper
bound \(h\le h_0\).
For any fixed \(\vartheta>0\) choose \(m\) and then
the peeling tolerance as in (PS4). Equations (PS3)–(PS6)
then show that all but a power-small fraction of pairs
have numerical weight at orders below \(k\), outside
\(S_i\), at least
\[
                             (1-\vartheta)h-2\eta-o_X(1).
\]
The high-order cutoff \(\kappa_0\) can be chosen independently
of the desired small relative loss \(\vartheta\).
In late balanced endpoint boxes, the uniform gap for \(h\)
allows the choices of \(m,\eta,\delta\) before the
fixed-budget height limits and uniformly across sufficiently
small budgets.

The near-height modulus gap gives a positive contribution
from bounded scaled orders in the whole quotient secant.
The simultaneous peeling and thick-pair count above put
nearly all of the remaining numerical scale into bounded
scaled orders, after the common projective modulus and
the fixed support are removed. Those removed primes can
still carry high numerical weight. The bounds do not
control the radicals of the moderate-order packets and
do not exclude the large-radical determinant configurations.
