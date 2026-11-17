# Rational line Fourier transforms: an all-degree branch obstruction

14 November 2026. **QFT1.** Logarithmic amplitudes and finite singularities
do not automatically escape [QPH1](QUINTIC_PHASE_AND_BRANCH_OBSTRUCTION.md).
The obstruction extends to every rational rank-one pre-Fourier connection
whose h=0 coefficient is a polynomial of degree at least three.

The theorem concerns the simple minimal extension of that line connection,
its Fourier transform, and finite algebraic covers of its position field.
It permits arbitrary finite pole orders and numbers, poles depending on
h and on the source position, and rational quantum corrections of any
finite degree. It is not a theorem excluding every Long carrier:
a different model can have generic pre-Fourier connection rank greater
than one, or involve additional boundary constituents requiring a separate
classical-support analysis.

## Set-up and the classical Long coefficient

Let k have characteristic zero, let h be the central deformation parameter,
and initially work with constants in k((h)). Source positions and r may
be included in the constant field for the a-variable analysis. Let
chi_h(a) be rational in a, regular in h valuation, with reduction
\[
             \chi_0(a)\in k[a],\qquad \deg_a\chi_0=\ell>=3.
                                                               \tag{1}
\]
The coefficient field may itself be a function field such as K(r,s);
finite algebraic constant extensions cause no difficulty.

Consider the rank-one rational differential module represented by
h*partial_a+chi_h, and take its simple minimal extension M to the affine
a-line. Our sign convention under Fourier transformation is
\[
  a\longmapsto h\partial_q,\qquad h\partial_a\longmapsto-q.
\]
Thus the corresponding exponential-module presentation has
q_left=-h*partial_a+chi_h. For the Long phase in QPH1,
\[
             \chi_0=(s-4a+3ra^2)^2/8,\qquad \ell=4
                            \quad\hbox{over }K(r,s).          \tag{2}
\]
The geometric target chart then uses q=v^2/2.

**Theorem.** The Fourier transform of M is a nonzero connection on a
dense position open set. Its differential module over k((h))(q) remains
irreducible after algebraic closure of constants and every finite
algebraic function-field extension. Consequently no algebraic projector
can select a proper generic branch after q=v^2/2 or after any additional
finite coefficient cover.

The assertions are about a simple line connection and its minimal
extension. The algebraic full right action is still required separately
for a carrier. Additional compatible right actions cannot produce a
projector absent already for this connection.

## Remove apparent line poles before counting rank

After algebraic closure of constants, a finite simple pole has
\[
                          \chi_h/h=\alpha/(a-c)+\cdots .
\]
If alpha is an integer, replace the rational frame by (a-c)^alpha
times that frame, removing this apparent pole. The resulting gauge is
rational. Repeat at every apparent regular pole. Irregular poles remain.

A rational frame change f changes chi_h by
-h*(partial_a f)/f. It does not change the reduction in (1):
normalize f by a scalar of its h valuation. Its leading rational
coefficient is nonzero, so
\[
                     v_h((\partial_a f)/f)>=0.               \tag{3}
\]
This also holds after a finite ramified extension of the scalar
parameter field. Poles moving to infinity do not invalidate it.

Let the resulting rational coefficient be
\[
       \widetilde\chi_h=P(a)+B(a)/A(a),\qquad
       \deg B<d:=\deg A,\quad \gcd(A,B)=1.                    \tag{4}
\]
Every remaining regular pole has nonintegral residue and every irregular
pole is genuine. The polynomial part P is unchanged by rational gauges.

For this cleaned line module a cyclic presentation of its minimal
extension is
\[
 D_a/D_a L,\qquad
       L=A(a)h\partial_a+A(a)P(a)+B(a).                       \tag{5}
\]
Here D_a is the one-variable polynomial Weyl algebra over the generic
constant field.

We spell out why cleaning is necessary. At a pole of order j>=2,
the coefficient B(c) is nonzero. On a finite derivative-of-delta vector
the B(c) term has the highest delta-derivative degree, while
(a-c)^j*partial_a lowers it. Thus L has no kernel on a delta module,
so there is no point-supported quotient. Its formal adjoint has the
same nonzero constant term there and excludes point-supported submodules.
At a simple pole the corresponding leading coefficient is proportional
to alpha-n-1 on the n-th delta derivative; the adjoint replaces alpha
by 1-alpha. Nonintegrality excludes both possibilities.

Away from these poles the rational connection has rank one and is
irreducible. Any proper submodule or quotient of (5) would therefore
be supported at its finite poles, just excluded. Thus (5) is the simple
minimal extension. This local argument uses the standard equivalence
between a point-supported holonomic module and a finite direct sum
of delta modules; it does not assume a rank formula from a finite sample.

## A nonconstant polynomial part gives no finite singularities

The normal-ordered Fourier operator associated to A*h*partial_a+C is
\[
 \widehat L=-q A(h\partial_q)-h A'(h\partial_q)
                                  +C(h\partial_q).           \tag{6}
\]
If P in (4) has degree m>=1, then C=AP+B has degree d+m.
The transformed operator has order d+m and nonzero constant leading
coefficient. Consequently the entire Fourier module is a nonsingular
polynomial connection on the affine q-line. Fourier transformation is
an equivalence, so it is irreducible.

It remains irreducible after finite algebraic covers. A nonsingular
connection has a fundamental matrix of entire functions over C.
Every algebraic element of its Picard--Vessiot field is meromorphic
and unramified on the affine line, hence rational. Thus its differential
Galois group is connected. Finite algebraic base change cannot reduce
that connected group or split its irreducible representation.
The same argument and finite-coefficient descent were used in QPH1.

For the standard Picard--Vessiot finite-part theorem and the entire-function
argument, see [Singer, Corollary3.10 and Section6 after Proposition9](https://arxiv.org/html/0712.4124).
We use the standard Galois correspondence; no independent verification
of all the results in that source is asserted.

## A constant polynomial part gives one finite singularity

It remains to handle P=c, including c=0. Replace q by z=q-c.
Then (6) becomes
\[
 \widehat L=-z A(h\partial_z)-h A'(h\partial_z)
                                      +B(h\partial_z).       \tag{7}
\]
Its order is d, with leading coefficient -A_d*h^d*z.
There is exactly one possible finite singular point, z=0.
Since all normalized coefficients have at most a simple pole there,
it is regular singular.

If d=0 the original cleaned coefficient is constant. It cannot satisfy
(1), because its rational degree and its reduction degree are zero.
In fact d>=ell: the rational function c+B/A has degree d, while
reduction of a rational function cannot increase its degree.
To verify this last statement even when c or pole locations have
parameter poles, choose a coprime numerator and denominator and
normalize their h valuations. Since the quotient has valuation zero,
their reductions are both nonzero. Cancellation on reduction only
decreases the rational degree. Equation (3) preserved chi_0, so
\[
                       N:=\operatorname{rank}\widehat M
                                   =d>=\ell>=3.              \tag{8}
\]

Let beta=B_(d-1)/(h*A_d), with beta=0 if that coefficient vanishes.
The indicial polynomial of (7), up to a nonzero scalar, is
\[
 \lambda(\lambda-1)\cdots(\lambda-d+2)
                                  (\lambda+1-\beta).        \tag{9}
\]
Only the top two differential coefficients contribute: the coefficient
of partial_z^(d-1), after normalization, has residue d-beta.
Therefore finite monodromy has at least N-1 eigenvalues equal to one.
Resonant Jordan blocks do not alter this eigenvalue statement.

## One finite monodromy of this type excludes every finite-cover split

We use the following general criterion.

Let E be an irreducible connection of rank N>=3 on the punctured affine
line, with no other finite singularities, and suppose its monodromy A
at the puncture has at least N-1 eigenvalues one. Then its differential
Galois identity component acts irreducibly. E therefore stays irreducible
over every finite algebraic function-field extension.

Here is a proof. The finite algebraic part of a Picard--Vessiot field
is unramified away from the singularities. A connected finite cover of
the punctured line can ramify only at the puncture and infinity.
If its degree is e and it has r0,rinf points over these two values,
Riemann--Hurwitz gives
\[
                  2g-2=-r_0-r_{\infty}.
\]
Hence g=0, r0=rinf=1, and the cover is a cyclic power map.
The Galois component group is therefore cyclic.

If the representation split over its component cover, the usual
finite-index restriction argument makes the restriction semisimple.
One can prove this directly: the orbit of an irreducible submodule
under the finite deck group spans the original irreducible module,
and that span is a sum of simple submodules. The endomorphisms
of the restricted module consequently have dimension greater than one.

Descent along a cyclic cover decomposes those endomorphisms into
characters. A nontrivial character therefore gives an isomorphism
\[
                       E\simeq E\otimes\mathcal K_\zeta,
                 \qquad \zeta\ne1,                           \tag{10}
\]
where K_zeta is a Kummer line with finite puncture monodromy zeta.
Equivalently use the projection formula
End(pi^*E)=direct_sum_j Hom(E,E tensor K_j).
Local monodromy would have to satisfy A similar to zeta*A.

That is impossible. Its spectrum consists of N-1 copies of 1 and
one nonzero value mu. Comparing traces of the first three powers gives
\[
 \begin{split}
 \mu&=1-N,\\
 \operatorname{tr}(A^2)&=N(N-1)\ne0
                       \quad\Longrightarrow\quad\zeta=-1,\\
 \operatorname{tr}(A^3)&=-N(N-1)(N-2)\ne0 .
 \end{split}                                                  \tag{11}
\]
The last line contradicts similarity to -A. This argument works for
every possible finite character, without bounding its order.

Thus restriction to the identity component is irreducible. Any finite
algebraic field extension leaves that component unchanged, proving the
criterion. Ordinary points introduce no hidden ramification: an algebraic
Picard--Vessiot element is locally meromorphic and single-valued there.
The argument uses no regularity assumption at infinity. It applies to
(7) by (8)--(9).

The proof is first made over C, then descends by the finitely many
coefficients and algebraic relations of any hypothetical finite-cover
submodule. This covers arbitrary characteristic-zero constant fields.

## The rank-two exception is real

For c nonzero consider the pre-Fourier line coefficient
\[
             \chi=-c/a^2+3h/(2a).
\]
Its Fourier equation, after multiplication by -1, is
\[
                  h^2 q\, y''+\tfrac12h^2 y'+cy=0.
\]
On q=v^2/2 it becomes
\[
                         h^2 y_{vv}+2cy=0.                  \tag{12}
\]
After adjoining kappa with kappa^2=-2c, the first-order system
h*partial_v*(y,h*y_v)^t=J*(y,h*y_v)^t, J=[[0,1],[-2c,0]],
has the constant horizontal projectors (I+J/kappa)/2 and
(I-J/kappa)/2. It really splits algebraically.

Its local monodromy eigenvalues are 1,-1, exactly the exceptional
case left by (11). This is a rank-two control, not a Long deformation:
its classical coefficient is rational of degree two rather than
a polynomial of degree at least three.

## Consequences for amplitudes and rational phase changes

For every finite collection of nonzero rational functions f_j(a,s),
a logarithmic amplitude changes the Long coefficient by
\[
          \chi_h=\chi_0-h\sum_j\alpha_j\,\partial_a\log f_j.
                                                               \tag{13}
\]
Each logarithmic derivative is proper at a=infinity, even when
the f_j have mixed a,s divisors or parameter-dependent poles.
Thus the generic polynomial part remains the nonconstant quartic
chi_0. The first case of the theorem already excludes every finite
algebraic branch splitting. There is no degree, pole-number or
residue search to perform in this class.

More generally any finite rational chi_h with the same reduction (2)
is covered, including rational phase corrections that change its
polynomial part at generic h. This distinction matters:
chi_0/(1+h*a^5) is proper generically although its reduction is quartic.
It is handled by the one-singularity and monodromy argument, not by
assuming its behavior at infinity remains polynomial.

Tensoring the Fourier module by any rank-one connection afterwards
does not help: tensoring by a line is an equivalence and its
endomorphism connection is unchanged. The same statement holds
after every finite algebraic base extension.

The rational connection must still have a compatible s derivative,
and ultimately the full fourth-image action, to be a carrier.
Logarithmic amplitudes with constant residues have the expected
mixed-derivative compatibility. The theorem excludes splitting
already for the single a/Fourier-q connection, so adding the other
relations cannot repair it.

Minimal extension is part of the statement. Nonminimal extensions
may add puncture-supported constituents before Fourier transformation.
Their classical support, parameter lattice, global admission and
finite-left generation need a separate audit. No theorem here turns
such a constituent into a Long carrier or rules out every conceivable
different algebraization.

## Why the generic pre-Fourier rank remains an actual gap

Classical or formal position rank one does not force generic position
rank one. A small exact control is
\[
 \mathcal E_R=D_R/D_R(hA^2+A-a),\qquad [A,a]=h.                \tag{14}
\]
It is h-flat. If hQ=U*(hA^2+A-a), reduction modulo h in the
polynomial domain forces U divisible by h; cancel h to prove
h saturation of the left ideal. Its special module is the
rank-one graph A=a, and flatness gives a formally free rank-one
position completion. Its generic connection rank is two, because
the defining operator has order two when h is invertible.

This control carries no constructed full Long right action.
It demonstrates why the present theorem cannot silently be promoted
to all-representative nonlifting. A candidate may have additional
generic connection branches invisible in the formal classical frame.
One must prove or construct the relevant rank and lattice properties.

The next useful Fourier construction must leave the rational line
ansatz, or establish an actual carrier using separately controlled
boundary constituents. A higher-rank pre-Fourier connection is a
remaining structural possibility. The full bounded scalar target
and [QMC1's finite-left graph carrier](FINITE_GRAPH_CARRIERS_AND_MATRIX_REALIZATION.md)
remain open.

[Exact controls](check_rational_fourier.py) verify the Fourier normal
ordering, two distinct rational behaviors at infinity, a mixed logarithmic
amplitude and its flatness identities, apparent-pole removal, the uniform
three-trace identities, the indicial formula in declared controls,
the rank-two split and the rank-jump control. Clean minimal extension,
absolute irreducibility, all finite covers and the application to all
rational deformations use written proofs and the cited standard
Picard--Vessiot correspondence. No independent theorem replay, formal
verification, external review, novelty priority, bounded Long lift,
finite-left Long carrier, all-representative nonlifting or DC2 solution
is claimed.

For the wider local Fourier framework see
[Bloch--Esnault](https://arxiv.org/html/math/0312343).
The normal-order and monodromy arguments above do not assume that
a local Fourier model is globally left finite.
