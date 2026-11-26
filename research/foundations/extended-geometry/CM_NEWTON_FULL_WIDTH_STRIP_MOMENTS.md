# Coherent input moments on the full even-conductor strip

The wider even-conductor strip allows all three top Laurent weights
\(n,n-1,n-2\). That exact bound is enough to recover the full coherent
source by unmarked moments when the source index is at most half the
conductor. At larger indices there is one further interval pole for
each angular input. Its contribution is coupled to odd-height terms
and can cancel a nonzero coherent input. A longer unmarked window
recovers both sets of inputs for all but finitely many values of any
one odd lower coefficient when the others are fixed. In the conductor-18
example, finite identities recover the coherent component with both
radial coefficients and one positive-angle odd coefficient arbitrary.

## The full-width theorem

Let \(n\ge6\) be even, \(3\le a\le n/2\), and let \(P(z)\) have
degree \(a+1\) and \(R\ge2\) distinct nonzero roots of positive
multiplicities \(m_1,\ldots,m_R\). Use the literal Laurent chart of
\(2xy+t^2=1\). Prescribe the top angular polynomial \(z^{-1}P(z)\)
on \(j+2i=n\), and permit every smaller-weight regular term obeying
the full support strip
\[
        j+2i\le n,\qquad (n-1)j+2ni\ge n(n-3).
\]
For a unit representative \(1\le s<n\), impose the full coherent
source condition
\[
 d_s=1+\left\lfloor\frac{as}{n}\right\rfloor
             -\sum_{\nu=1}^R\left\lfloor\frac{sm_\nu}{n}\right\rfloor=R.
\]

**Theorem.** If \(s\le n/2\), for each fixed phase \(f\) of this
shape over a characteristic-zero field and each sufficiently large
rational prime \(p\equiv s^{-1}\pmod n\), at least one of
\[
          \mathcal A(f^{m+kp}),\qquad
          m=(sp-1)/n,\qquad 0\le k<R
\]
is nonzero. The same conclusion for \(s>n/2\) holds if there are no
literal terms of weight \(n-1\) at angles \(0,1,\ldots,R-2\).
In rank two the additional condition excludes just the radial
coefficient of \(t^{n-1}\).

Every such phase has least sphere degree \(n+1\) and is an absolute
nonpower. The certifying prime set in each class has natural density
\(1/\varphi(n)\). If \(s\) generates the units modulo \(n\), then
the primitive Jacobian factor of \(W^n=z^{n-1}P(z)\) has
\(p\)-rank zero in every good reduction over these primes.

### The exact support bound

Put \(d=n-j-2i\). The old strip becomes
\(0\le d\le 2(n+i)/(n-1)\). Since \(j\ge0\), also
\(i\le(n-d)/2\), giving \(d\le3\). If \(d=3\), both inequalities
must be equalities and \(i=(n-3)/2\), impossible for even \(n\).
Thus \(d\) is 0, 1 or 2. At angle \(-v\) with \(v>0\), regularity
requires a height width of at least \(2v\) from the factor
\((1-t^2)^v\). The two support bounds permit width at most
\(2(n-v)/(n-1)\), so \(v=1\) is the only possible negative angle.
Its coefficient is exactly \(P(0)x^{-1}(t^{n+2}-t^n)\).
Consequently, writing \(z=x/t^2\), the literal phase has the form
\[
   f=t^n\bigl(Q(z)+t^{-1}B(z)+t^{-2}D(z)\bigr),\qquad
   Q(z)=z^{-1}P(z),\quad D(z)=-P(0)z^{-1}+D_+(z),
\]
where \(B,D_+\) are ordinary polynomials. Each positive-angle
top monomial has nonnegative height by \(n\ge2a\). The regular
negative-angle lift \(-2P(0)yt^n\) has degree \(n+1\). The
nonnegative angles have degree at most \(n\) in the normal sphere
basis, and passage to that basis cannot raise degree. A proper
power in the rational function field would give a proper power
as its leading coefficient in the expansion at \(t=\infty\)
over \(\overline{k}(z)\). The order of \(Q(z)\) at \(z=0\) is
\(-1\), excluding such a power.

### Interval poles and coherent inputs

At a good prime use \(f^{m+kp}=f^m(f^k)^p\) modulo \(p\).
All angular-constant pairs have angles \((jp,-j)\) with
\(0\le j\le k\), since \(f^m\) has minimum angle at least
\(-m>-p\) and \(f^k\) has minimum angle at least \(-k\).
The low height at angle \(jp\) in \(f^m\) is
\((s-2j)p-1\) minus the accumulated weight deficit. The deficit
is at most \(2m<2p\). The interval denominator is a multiple of
\(p\) only when the deficit is 0 or \(p\). Define
\[
 \begin{split}
 A_j&=[z^{m+jp}]P(z)^m,\\
 C_j&=[z^{jp}\xi^p](Q(z)+\xi B(z)+\xi^2D(z))^m.
 \end{split}
\]
These are respectively the even and odd low-height coefficients
that can contribute. A coefficient at an impossible height is zero.
Write \(F_{k,j,u}=[x^{-j}t^u]f^k\) in the literal Laurent ring.
For \(p>s+kn\) the maximum interval denominator is less than
\(p^2\). This gives the exact reduction formula
\[
 p\mathcal A(f^{m+kp})\equiv
 \sum_{j=0}^k A_j\sum_{u\,\mathrm{even}}
                   \frac{F_{k,j,u}^{\,p}}{s-2j+u}
 +\sum_{j=0}^k C_j\sum_{u\,\mathrm{odd}}
                   \frac{F_{k,j,u}^{\,p}}{s-2j-1+u}.
\]
Only terms that occur contribute to the displayed sums. In the
second sum an index with \(C_j=0\) is omitted. The denominators
of contributing terms lie between 1 and \(s+kn\).

For \(s\le n/2\), all the \(C_j\) vanish because \(2m<p\).
For the larger \(s\) covered by the theorem, an odd-height term of
\(f^k\) must involve a term from \(B\), at angle at least \(R-1\).
The other \(k-1\) angles are each at least \(-1\). Thus its total
angle is at least \(R-k>0\) for \(1\le k<R\), and the second sum
again vanishes.

The remaining system is triangular in the \(A_j\). On the diagonal
the only contribution to angle \(-k\) comes from powers of the exact
negative lift. With \(L=s+kn-2k\), the coefficient of \(A_k\) is
\[
 P(0)^{kp}\sum_{v=0}^k\frac{(-1)^{k-v}\binom{k}{v}}{L+2v}
       =\frac{(-1)^k2^k k! P(0)^{kp}}{\prod_{v=0}^k(L+2v)}\ne0
              \pmod p.
\]
The full coherent source condition ensures \(s\ge 2R-1\), since
\(R-1\le\lfloor as/n\rfloor\le \lfloor s/2\rfloor\) and \(s\)
is odd. On the normalized connected cover \(W^n=z^{n-1}P(z)\),
the basis \([W^s/z^{s+j}]\), \(0\le j<R\), has Frobenius entries
\(A_j\) on the target line \([W/z]\). The selected source spans the
entire de Rham character space by \(d_s=R\), so the Cartier
factorization gives surjectivity onto the target coherent line.
The inputs therefore cannot all vanish. The same eigensheaf argument
gives \(d_{-s}=0\). If \(s\) generates the units modulo \(n\),
every primitive character lies on the same orbit and its coherent
Frobenius is nilpotent. These facts and the fixed number-field
specialization for general characteristic-zero coefficient fields are
proved in the [coherent input criterion](CM_NEWTON_TWO_INPUT_STRIP_MOMENTS.md#a-full-coherent-input-space-of-arbitrary-rank).

This proves the all-parameter conclusions in degrees 11, 15 and 19
on the full wider strips. Their full source indices are respectively
\(3,3,5\), all below half the respective conductors \(10,14,18\).
The fixed-rank unbounded-conductor families in the coherent input
criterion also satisfy this bound.

## A four-point family of uniform positive density

For \(h\ge3\) put \(N=3^{h-1}\), \(n=6N\), \(a=2N-1\), and
choose both root multiplicities equal to \(N\). The two distinct
nonzero roots and the leading scalar are arbitrary. Use all
\(1\le s<3N\) with \(s\equiv5,11\pmod{18}\). If
\([t^{n-1}]\operatorname{CT}_x f=0\), also use
\(3N<s<4N\) in those two residue classes. To verify the coherent
dimension, write \(s=6q+5\). Then
\[
 \left\lfloor\frac{sm_1}{n}\right\rfloor
 =\left\lfloor\frac{sm_2}{n}\right\rfloor=q,
 \qquad
 \left\lfloor\frac{as}{n}\right\rfloor
   =\left\lfloor\frac{s}{3}-\frac{s}{n}\right\rfloor=2q+1
\]
for \(s<4N\), giving \(d_s=2\). These \(s\) generate the units
modulo \(n\) because their reductions modulo nine are 5 or 2,
and \(v_3(s^{2\cdot3^j}-1)=1+j\). Thus the primitive factor has
\(p\)-rank zero on all these good certifying classes.

There are \(N/3\) such indices below \(3N\) and \(N/9\) between
\(3N\) and \(4N\). Since \(\varphi(n)=2N\), their disjoint prime
progressions give density \(1/6\) with all compatible lower terms
arbitrary, and density \(2/9\) if the single odd radial coefficient
vanishes. Both densities are independent of \(h\).

## Exact checks and the upper-half obstruction

At \(n=10,a=5,p=17,s=3\), take \(P=(z^2-1)^3\) and add the
strict lower terms \(4t^8+9x^4t\) to the regular top lift. Both
added terms lie in the wider strip and outside the [narrower strip](CM_NEWTON_TWO_INPUT_STRIP_MOMENTS.md#the-strip-and-its-two-moments).
The normalized residues at orders 5 and 22 are still 0 and 16.

The secondary correction at \(s>n/2\) is genuine. At
\(n=18,a=7,s=13,p=43\), take \(P=(z^2-1)^4\), \(B=6\),
and \(D=-z^{-1}\). The source has dimension two and \(m=31\).
The coherent inputs are \(A_0=0\), \(A_1=5\) modulo 43, while
\(C_0=13\). The two normalized unmarked residues are both zero.
In the second residue the contributions are
\(5(1/31-1/29)=24\) and \(13\cdot6/29=19\) modulo 43.
This example rules out recovering the two inputs from these two
residues for arbitrary subleading terms at upper-half indices.
[Exact replay of both cases](verify_cm_newton_full_width_strip_moments.py)

## A reconstruction matrix for unmarked moments

Put \(s=2q+1\). The low heights for \(A_j,C_j\) are respectively
\((s-2j)p-1\) and \((s-1-2j)p-1\), so it suffices to allow
\(A_0,\ldots,A_q\) and \(C_0,\ldots,C_{q-1}\), including any
coefficients that vanish for further support reasons. In the
\((z,t)\) chart write \(G_{k,j,u}=[z^{-j}t^u]f^k\). Define the
column indexed by \((j,\epsilon)\), with \(\epsilon=0\) for \(A_j\)
and \(\epsilon=1\) for \(C_j\), by
\[
 M_{k,j,\epsilon}
   =\sum_{u\equiv\epsilon\ (2)}\frac{G_{k,j,u}}{s-\epsilon+u}.
\]
The corresponding literal \(x\)-chart height in \(f^k\) is
\(u+2j\). At any sufficiently large good prime for the chosen finite
rows, the vector of normalized unmarked residues is \(M^{(p)}(A,C)^T\),
where \(M^{(p)}\) is the reduction followed by entrywise Frobenius.

For a prescribed window the kernel has zero coherent projection
exactly when
\[
 \operatorname{rank}(M_A\mid M_C)
   = \operatorname{rank}M_C+\#\{\text{coherent columns}\}.
\]
This follows by taking the images of the coherent columns modulo
the column space of \(M_C\). Known coherent relations can be
incorporated by restricting the domain of \(M_A\) first.

For \(B\ne0\), put \(\ell=\operatorname{ord}_z B\) and
\(b_\ell=[z^\ell]B\). The extra pole block is triangular on
rows \(\ell+1,\ldots,\ell+q\). The first odd product of
angle \(-j\) occurs in row \(J+1\), where \(J=j+\ell\).
It uses one angle-\(\ell\) odd factor and \(J\) negative lifts.
The pivot in column \(C_j\) is
\[
 \frac{(J+1)b_\ell(-2P(0))^J J!}
      {\prod_{v=0}^J\bigl(s+(n-2)(J+1)+2v\bigr)} .
\]
Rows \(0,\ldots,\min(\ell,q)\) contain only coherent entries
and recover them triangularly up to that index. In particular,
when \(\ell\ge q\), those rows and the shifted pole block recover
all inputs. In the stratum \(\ell=0\) this gives the original
elimination using rows \(E=\{1,\ldots,q\}\). Write each row as
\(U_k A+V_k C\). After these rows vanish, the remaining equations
are \((U_k-V_kV_E^{-1}U_E)A=0\).

Fix the top data, so that \(P(0)\) is a nonzero scalar.
After clearing denominators, the residual minors are polynomials
on the coefficient locus \(b_0=\cdots=b_{\ell-1}=0\). For any
finite choice of rows, full coherent rank at every geometric
point with \(b_\ell\ne0\) is equivalent to the maximal minors
generating the unit ideal after inverting \(b_\ell\). Equivalently
their ordinary ideal contains a power of \(b_\ell\) before
localization. This follows from the Nullstellensatz and clearing
the remaining denominators in a unit-ideal identity.

At a fixed \(t\) where the angular Laurent polynomial has both a
simple pole at zero and a pole at infinity, eventual vanishing of
all angular moments against a polynomial differential forces that
differential to vanish, by Pakovich's [Theorem 6.6](https://arxiv.org/pdf/0910.2105)
(arXiv:0910.2105v1, §6). Here the entries also integrate the angular
constant terms over \(t\), against \(t^{s-1-\epsilon}\).
The angular theorem gives no implication from vanishing after that
integration.

### Effective bounds for the reconstruction window

**Proposition (computable reconstruction cutoff).**
For fixed \(n,s\) and fixed support, a finite initial window
has the same kernel as all reconstruction rows, simultaneously
at every characteristic-zero specialization. A uniform cutoff
can be computed on any coefficient locus given by polynomial
equations over \(\mathbb Q\). At an individual phase over a
presented finitely generated field one can compute the cutoff
without calculating it on a whole coefficient locus.

The column \((j,\epsilon)\) has the regular marker
\[
 h_{j,\epsilon}=x^j t^{s-1-\epsilon-2j},
\]
since its \(t\)-exponent is nonnegative. In the \((z,t)\) chart
this is \(z^j t^{s-1-\epsilon}\), so \(M_{k,j,\epsilon}
=\mathcal A(h_{j,\epsilon}f^k)\). Use a series variable \(T\)
and write \(g_i=h_i/(1-Tf)\) for each marker. On the localized
smooth affine quadric over \(K(T)\) use the invariant residue
form \(\omega\). The sphere trace annihilates exact top forms by
Stokes, using the [closed sphere realization](GLOBAL_SPHERE_PURITY_RECOGNITION.md#1-a-closed-rational-period)
after an embedding of the finitely generated coefficient field
into \(\mathbb C\). The rational coefficient formula descends
the identity to the original field.
In the relative de Rham complex the coefficientwise extension
to Laurent series at \(T=0\) still kills exact forms and
commutes with \(\partial_T\).

Apply the [indicial argument of the finite-tests manuscript](../../../unpublished/finite-tests-for-mixed-moments/main.tex)
(§4, Lemma 4.1) in the direct sum of one top relative de Rham
space per column. These spaces are finite-dimensional. The
tuples \(([\partial_T^j(g_i\omega)])_i\) therefore have a
linear dependence with the same \(p_j(T)\in K[T]\) for every
component. It is computable by enumerating bounds on derivative
order, polynomial degrees and powers of \(1-Tf\) in the exact-form
certificates, and solving the simultaneous linear identity
in the relative coordinate ring. The search requires the
\(p_j\) to be nonzero as a tuple. Clearing denominators in any
de Rham dependence gives such a certificate at finite bounds.

Then \(\sum_j p_j(T)\partial_T^j\) annihilates all the
column series and every constant linear combination of them.
With \(v_j=\operatorname{ord}_0 p_j\) for \(p_j\ne0\) and
\(\delta=\min_j(v_j-j)\), the polynomial
\[
 I(N)=\sum_{v_j-j=\delta}[T^{v_j}]p_j\,(N)_{\underline j}
\]
is nonzero by the distinct degrees of the falling factorials.
The first nonzero coefficient of any such linear combination,
at index \(N\), gives \(I(N)=0\) by inspecting degree \(N+\delta\).
To bound the nonnegative integer roots over a function field
\(\operatorname{Frac}(\mathbb Q[\mathbf c]/\mathfrak p)\), clear
the parameter denominators and reduce the coefficients modulo a
Gröbner basis of \(\mathfrak p\). Any nonzero rational coordinate
polynomial then gives a rational root bound. Exact sphere moments
through that bound determine the kernel by linear algebra.
This decides whether it contains a vector with nonzero coherent
projection, also when restricted to a prescribed coherent subspace.

To compute one cutoff on a parameter locus, start with its
finitely many irreducible components over \(\mathbb Q\), using
the same effective minimal-prime operation as in
[the global sphere prefix argument](GLOBAL_SPHERE_PURITY_RECOGNITION.md#4-a-computable-degree-only-prefix).
For a component with prime \(\mathfrak p\), compute the generic
all-row rank \(r\) by the preceding cutoff over its fraction field.
If \(r=0\) then all entries vanish on the component. Otherwise
choose one nonzero \(r\)-minor from the prefix. All \((r+1)\)-minors
involving any rows vanish in the coordinate ring, since the
ring embeds in the fraction field. Wherever the chosen minor
stays nonzero, the prefix and the full matrix both have rank
\(r\), hence the same kernel. Recurse on the components where
the chosen minor vanishes. Each has smaller dimension and the
branching is finite. The maximum of the finitely many cutoffs
is a uniform cutoff. This argument computes equality of kernels
after characteristic-zero specialization. The abstract equality
of row submodules over the parameter ring also follows from
Noetherianity, without computing their generators.

Once this cutoff is computed, existence of a geometric blind
direction on a prescribed algebraic locus is decidable as well.
For each coherent coordinate \(A_j\), add the equation \(A_j=1\)
to the prefix equations in the phase and input variables.
Homogeneity in the inputs makes some such normalization possible
for every nonzero coherent projection. Nonvanishing conditions
on phase coefficients can be encoded by an inverse variable.
Testing the resulting finitely many ideals for the unit ideal
decides existence over \(\overline{\mathbb Q}\).

### Varying one odd lower coefficient

**Theorem (reconstruction on odd coefficient lines).**
In the full-width strip with \(P(0)\ne0\), vary the coefficient
\(\lambda\) of any allowed odd term \(x^\ell t^{n-1-2\ell}\)
and fix every other coefficient. Put \(L=\ell+1\). Partition the
\(s\) columns \((j,\epsilon)\) by \(j+L\epsilon\bmod 2L\).
For each nonempty block, let \((j_0,\epsilon_0)\) be its first
column in increasing \(j\) and \(d\) its number of columns.
Use the \(d\) shifts
\[
 k=j_0+L\epsilon_0+2Li,\qquad 0\le i<d .
\]
The determinant of the resulting \(s\) rows is a nonzero polynomial
in \(\lambda\) of degree \(\sum_{\text{blocks}}(d\epsilon_0+d(d-1)/2)\).
All shifts satisfy \(k\le 2q+L\). For each fixed phase over a number
field where the determinant is nonzero, these unmarked residues
recover all \(A_j,C_j\) at every sufficiently large prime
\(p\equiv s^{-1}\pmod n\), at every place above \(p\). Under the
full coherent source condition at least one corresponding moment
is nonzero. The nonvanishing assertion also holds over any
characteristic-zero coefficient field.

For the radial coefficient \(b\), the two blocks use the shifts
\[
 \mathcal K_s=
 \begin{cases}
  \{0,1,\ldots,2q\},&q\text{ even},\\
  \{0,1,\ldots,2q-1,2q+1\},&q\text{ odd}.
 \end{cases}
\]
The degree is \(q(q+1)\) for even \(q\) and \(q(q+1)+1\)
for odd \(q\).

In the \((z,t)\) chart, \(r\) copies of \(\lambda z^\ell
t^{n-1}\) in \(f^k\) force angle at least \(Lr-k\).
For angle \(-j\) this gives \(r\le(k-j)/L\). Set
\(\lambda=w^L\), and scale row \(k\) by \(w^{-k}\) and column
\((j,\epsilon)\) by \(w^j\). The constant terms at \(w=\infty\)
come from equality at \(r=(k-j)/L\), with all \(e=k-r\)
other factors taken from \(P(0) z^{-1}t^{n-2}(t^2-1)\).
The height parity requires \(r\equiv\epsilon\pmod2\), giving
the stated blocks. The leading entry when \(r,e\ge0\) is
\[
 \frac{(-2P(0))^e k!/r!}
      {\prod_{\nu=0}^e\bigl((n-2)k+r+s-\epsilon+2\nu\bigr)}.
\]
Expand \((t^2-1)^e\) against the monomial weight and take
successive finite differences of reciprocals to obtain the
denominator product.

Within a block the columns have \(j=j_0+Li\), \(0\le i<d\),
with \(\epsilon\) alternating. Write \(r_0=(k-j_0)/L\).
Relative to the first column, column \(i\) has \(r=r_0-i\)
and \(e=k-r_0+i\). For \(k\ge j_0\) of the required residue
the first column is nonzero. Divide each row by its first entry.
At each increment of \(i\), the grid of denominator factors gains
one factor at either the lower or the upper endpoint, because
\(\epsilon\) alternates. Thus the normalized column \(i\) is
\((-2P(0))^i(r_0)_{\underline i}\) divided by the product of the
\(i\) added affine factors in \(k\). This expression includes
zero entries with \(r_0<i\) through the falling factorial.
All grid factors are positive at nonnegative \(k\) where evaluated,
since even the lower endpoint is
\((n-2)k+(k-j)/L+s-\epsilon>0\).

Use the product of the \(d-1\) added factors as a common
denominator. The normalized numerators are polynomials of degree
at most \(d-1\) in \(r_0\). They are independent by triangular
evaluation at \(r_0=0,1,\ldots,d-1\), where the added factors
are nonzero and the \(i\)th falling factorial first becomes nonzero
at \(i\). They form a basis and therefore also give a nonsingular
evaluation at the chosen distinct \(r_0=\epsilon_0+2i\).
It follows that the scaled determinant has a nonzero constant
term at infinity. Undoing the scaling gives the degree in \(w\)
as \(\sum k-\sum j\). The contribution of a block is
\(L(d\epsilon_0+d(d-1)/2)\); dividing by \(L\) proves the
degree in \(\lambda\).

A fixed nonzero algebraic determinant is a unit at all but
finitely many places, and its Frobenius power remains a unit.
For a general characteristic-zero field, specialize the finitely
generated coefficient ring to a number field while preserving
the determinant, all nonzero roots, and their distinctness and
multiplicities. This is the fixed-specialization argument in
the [coherent input criterion](CM_NEWTON_TWO_INPUT_STRIP_MOMENTS.md#a-full-coherent-input-space-of-arbitrary-rank),
with one additional element inverted. The resulting nonzero
moments specialize from the original phase.

### All parameters in a conductor-18 family

**Theorem (conductor 18).** Take the phase on \(2xy+t^2=1\)
\[
 f_{b,c,d}=x^7t^4-4x^5t^8+6x^3t^{12}-4xt^{16}
             -2yt^{18}+bt^{17}+ct^{16}+dx^2t^{13}.
\]
For every fixed algebraic \(b,c,d\) and every sufficiently large
prime \(p\equiv7\pmod{18}\), put \(m=(13p-1)/18\). At least one
of the moments of orders \(m+kp\) is nonzero with \(0\le k\le7\).
If \(d=0\), shifts \(0\le k\le5\) suffice. If \(c=d=0\), the
four shifts \(0,1,3,5\) suffice.

In the \((z,t)\) chart this is \(P=(z^2-1)^4\), \(B=b+dz^2\)
and \(D=-z^{-1}+c\). At any fixed good residue place, put
\(w=z^2\) and write \(4m=2p+\ell\), where
\(\ell=(8p-2)/9<p\). For \(p\equiv25\pmod{36}\), \(m\) is even
and \(A_0=[w^{m/2}](w-1)^{4m}\) is nonzero because
\(0\le m/2\le\ell<p\). This follows by factoring
\[
 (w-1)^{4m}=(w^p-1)^2(w-1)^\ell
 \quad\text{in characteristic }p.
\]
The \(k=0\) row is \(A_0/13\), for any lower terms in the strip.
For \(p\equiv7\pmod{36}\), \(m\) is odd, and the angular bound
\(7m<6p\) and parity leave only \(A_1,A_3,A_5\). Their indices
in \(w\) are \((m+p)/2+ip\) for \(i=0,1,2\), with
\(0\le(m+p)/2\le\ell\). The same factorization gives
\((A_1,A_3,A_5)=(A,-2A,A)\) with \(A\ne0\).

For the \(C_j\), compare the maximum angle with \(7m\).
A deficit-one factor costs at least five in angle, and a
deficit-two factor costs at least seven. At total deficit \(p\)
the maximum angle is thus at most \(7m-7p/2<2p\), leaving
only \(C_0,C_1\). If \(c=d=0\), the number of \(B\) factors at
deficit \(p\) is odd. All other factors have odd angular
exponents, so since \(m\) is odd the total angle is even and
\(C_1=0\).

Write \(u=b^p\), \(v=c^p\) and \(r=u^2\) in the residue field.
If \(u=0\) the \(k=1\) row is \(-2A/899\). For \(u\ne0\), set
\[
 U=uC_0/A,\qquad V=uC_1/A,\qquad h=d^p/u.
\]
Dividing the residues by \(A\), rows 1 and 2 are
\[
 -\frac2{899}+\frac{U}{29},\qquad
 -\frac{4v}{2115}+\frac{2vU}{45}-\frac{4V}{2115}.
\]
Their simultaneous vanishing gives \(U=2/31\) and \(V=16v/31\).

When \(h=0\), after these substitutions the next three rows,
up to nonzero rational factors, are
\[
 \begin{aligned}
  &26724-28609 r,\\
  &v(33636-44737r),\\
  &(4535708-7209279r)v^2-1383599r^2+351436r+652208.
 \end{aligned}
\]
The first fixes \(r=26724/28609\), and
\(33636\cdot28609-44737\cdot26724\ne0\) forces \(v=0\) in the
second for all sufficiently large primes. With
\(L=28609r-26724\) and \(H=1383599r^2-351436r-652208\),
the exact identity
\[
 28609^2H-(28609\cdot1383599r+26921067152)L
       =185622733382800
\]
excludes a common zero with the third. When \(c=d=0\), row 1
already fixes \(U\) and \(C_1=0\), so only rows 3 and 5
are needed for the same identity. The two minors in this scalar
case are nonzero rational multiples of \(uL\) and \(uH\).

For general \(h\), row 3 now forces \(h=427r/6-4454/67\).
After eliminating \(h\), row 4 is proportional to \(v(7r-3)\).
If \(v\ne0\), substitute \(r=3/7\) into rows 5 and 6. They
give the two equations
\[
 -624512896v^2-649463671=0,\qquad
 199006080v^2+222720503=0.
\]
The determinant of their coefficient and constant columns is
nonzero. If \(v=0\), rows 5 and 7 instead give
\[
 \begin{aligned}
 Q_5(r)&=5555305r^2-9164116r-85804=0,\\
 Q_7(r)&=1417868646401r^3-1910019598132r^2
                    -1292362283932r-11534389824=0 .
 \end{aligned}
\]
These have the same degrees modulo 3 and reduce to
\(r^2-r-1\) and \(-r^3-r^2-r\). The latter has remainder
\(1-r\) modulo the former, and the former is nonzero at \(1\).
Their resultant is therefore a nonzero integer.

The row identities are obtained by iterating multiplication by
\(t^{18}Q(z)+t^{17}(u+uhz^2)+t^{16}(-z^{-1}+v)\) from \(1\),
and using the reconstruction matrix above. For the coherent
column use angles \(-1,-3,-5\) with weights \(1,-2,1\), and
for the remaining columns use angles \(0,-1\) and divide their
entries by \(u\). Parity leaves polynomials in \(r,v,h\).
The [exact checker](verify_cm_newton_full_width_strip_moments.py)
generates these polynomials by finite multiplication over the
rationals and verifies the stated eliminations. The replay used Python
3.12.3 from the repository root.
```sh
python3 research/foundations/extended-geometry/verify_cm_newton_full_width_strip_moments.py
```
All divisions
and nonzero constants in the argument exclude only finitely many
rational primes. Nonzero reduction of a normalized residue
implies nonvanishing of the original moment.

### A determinant zero locus resolved by another row

**Theorem (full reconstruction on an exceptional radial line).**
Keep \(n=18\), \(s=13\) and \(P=(z^2-1)^4\). On the line
\[
 B=b+z^3,\qquad D=-z^{-1},
\]
the unmarked rows with \(0\le k\le13\) recover every coherent
and extra pole input, for every fixed algebraic \(b\) at every
sufficiently large prime \(p\equiv7\pmod{18}\).
The original window \(0\le k\le12\) is singular at exactly
42 distinct algebraic values of \(b\), all nonzero.

In the order \((A_0,\ldots,A_6,C_0,\ldots,C_5)\), form the
two square minors
\[
 \Delta_0=\det M_{0,\ldots,12},\qquad
 \Delta_1=\det M_{0,\ldots,11,13}.
\]
The [certificate checker](verify_cm_newton_full_width_reconstruction_certificates.py)
constructs the rows by multiplying the Laurent phase and checks
over \(\mathbb Q[b]\) that
\[
 \deg\Delta_0=42,\quad \Delta_0(0)\ne0,\quad
 \gcd(\Delta_0,\Delta_0')=1,\quad
 \deg\Delta_1=39,\quad
 G_0\Delta_0+G_1\Delta_1=1.
\]
It computes \(G_0,G_1\) by extended Euclidean division and
verifies the polynomial identity by multiplication. Thus every
root of the first determinant has a nonzero replacement
determinant. The identity descends at all sufficiently large
primes, and entrywise Frobenius preserves the full rank.
The coherent source has \(d_{13}=1+\lfloor91/18\rfloor
-2\lfloor52/18\rfloor=2\), so the window contains a nonzero
moment. This certificate uses all 13 height-allowed columns,
including columns that vanish by support or by Lucas in this family.

### A unit ideal after triangular elimination

**Theorem (odd quartic and even radial terms).** With the same
\(n,s,P\), let
\[
 B=b+b_1 z+b_2 z^2+b_3 z^3+b_4 z^4,\qquad D=-z^{-1}+c.
\]
For every fixed algebraic choice of these coefficients and every
sufficiently large \(p\equiv7\pmod{18}\), one of the moments
of order \(m+kp\), \(0\le k\le9\), is nonzero. If \(b_4=0\),
shifts \(0\le k\le8\) suffice.

The \(p\equiv25\pmod{36}\) case again uses row 0. For
\(p\equiv7\pmod{36}\), the coherent vector is \((A_1,A_3,A_5)
=(A,-2A,A)\), where \(A\ne0\). At deficit \(p\), let \(e_1,e_2\)
count factors from \(B,D\). Their angle losses from \(7m\)
are at least \(3e_1+7e_2\). Since \(e_1+2e_2=p\)
and \(e_1+e_2\le m\), one has \(e_1\le2m-p\). The angle is
at most \(7m-7p/2+e_1/2\le8m-4p<2p\), leaving only \(C_0,C_1\).
If \(b=0\), row 1 detects \(A\) directly.

For \(u=b^p\ne0\), put \(r=u^2\), \(v=c^p\),
\(\ell_j=b_j^p/u\), \(U=uC_0/A\), and \(V=uC_1/A\).
Rows 1 and 2 force
\[
 U=2/31,\qquad V=(16v-2\ell_1)/31.
\]
In row \(k\), for \(k=3,4,5\), the first occurrence of
\(\ell_{k-1}\) comes from one such odd term and \(k-1\)
negative lifts at angle zero. Its coefficient is therefore
a nonzero rational constant after substituting \(U=2/31\).
Use these rows to solve for \(\ell_2,\ell_3,\ell_4\).
After those substitutions, rows \(k=6,7,8,9\) give four polynomials in
\(\mathbb Q[r,\ell_1,v]\), up to nonzero rational scalings.
Their total degrees are \(5,6,7,8\). The
[Bézout multipliers](cm_newton_conductor_18_quartic_bezout.json)
give a linear combination equal to \(1\) in this ring.
They were obtained by multiplying the four rows up to total
degree 15 and solving the rational coefficient system. Reduction
modulo \(65521\) selected the pivot rows and columns. The checker
verifies the identity by direct multiplication over \(\mathbb Q\).
For the cubic subfamily, omit the row-5 substitution and use
rows \(k=5,6,7,8\), of total degrees \(4,5,6,7\). The checker
computes and multiplies the four rational multipliers for that
case as well. This precludes simultaneous vanishing
even at a singular specialization of the original determinant.
Reversing the linear substitutions and clearing the powers of
\(u\) gives a localized unit-ideal certificate for the coherent
line singled out by Lucas. All its rational constants and
identities transfer to characteristic \(p\) away from a finite
set of primes.

Both certificates replay with SageMath 10.9 from the repository root.

```sh
sage -python research/foundations/extended-geometry/verify_cm_newton_full_width_reconstruction_certificates.py
```

The general upper-half question is whether extra unmarked rows
force the coherent projection of the kernel to vanish at every
allowed lower polynomial. The determinant above covers the
complement of a finite set on each odd coefficient line. The two
coprime minors settle full reconstruction on the line \(B=b+z^3\),
and the unit ideal settles the actual coherent direction for the
larger odd-quartic family. For arbitrary exceptional parameters the
computable cutoff decides the all-row kernel at a fixed phase. A
kernel direction with nonzero coherent component concerns the residue
system. Its relevance to moment vanishing would also require an
argument that Frobenius inputs can occupy it.
