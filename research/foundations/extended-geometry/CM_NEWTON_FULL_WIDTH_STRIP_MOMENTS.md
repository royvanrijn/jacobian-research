# Coherent input moments on the full even-conductor strip

The wider even-conductor strip allows all three top Laurent weights
\(n,n-1,n-2\). That exact bound is enough to recover the full coherent
source by unmarked moments when the source index is at most half the
conductor. At larger indices there is one further interval pole for
each angular input. Its contribution is coupled to odd-height terms
and can cancel a nonzero coherent input. A longer unmarked window
recovers both sets of inputs for all but finitely many values of any
one odd lower coefficient when the others are fixed. For the
conductor-18 top polynomial, the prime class \(25\bmod36\) already
gives a nonzero shift-zero residue throughout the full strip.
The longer windows address \(7\bmod36\), where that residue vanishes.
A weight argument also excludes Hamiltonian brackets from the marker
space and excludes generic-fibre exactness for the coherent profiles.

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
In particular the even class contains a fixed-prime witness.

**Corollary (the conductor-18 witness at 61).**
At \(p=61\) the order is \(m=44\) and \(\ell=54\).
For every phase with this top polynomial and any compatible lower
terms integral at a place above \(61\), the normalized residue is
\[
 61\,\mathcal A(f^{44})
      \equiv \frac{1}{13}\binom{54}{22}\equiv5\pmod{61}.
\]
This single moment certifies nonvanishing throughout that local
integrality class. The same prime class gives a nonzero shift-zero
moment at every sufficiently large prime for each fixed algebraic
phase on the full strip.
For \(p\equiv7\pmod{36}\), \(m\) is odd, and the angular bound
\(7m<6p\) and parity leave only \(A_1,A_3,A_5\). Their indices
in \(w\) are \((m+p)/2+ip\) for \(i=0,1,2\), with
\(0\le(m+p)/2\le\ell\). The same factorization gives
\((A_1,A_3,A_5)=(A,-2A,A)\) with \(A\ne0\).
In this harder prime class the shift-zero residue is zero regardless
of the lower terms. The later reconstruction windows address this
case.

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

### Three extra pole inputs for odd degree five and even degree one

Keep \(n=18\), \(s=13\), and \(P=(z^2-1)^4\). Consider
\[
 B=b+\sum_{j=1}^{5} b_jz^j,\qquad D=-z^{-1}+c+ez.
\]

**Theorem (three-pole coordinate slices and a fixed prime).**
If \(bb_1b_2ce=0\), for every fixed algebraic choice of the
coefficients and every sufficiently large \(p\equiv7\pmod{18}\),
one of the moments of order \(m+kp\), \(0\le k\le16\), is
nonzero. The same large-prime conclusion holds when \(b\ne0\)
and \(b^2,b_1/b,b_2/b,c,e\) are integral at one place above
\(331\). For every phase with all coefficients integral at a
place above \(331\), one of the seventeen moments
\(\mathcal A(f^{239+331k})\), \(0\le k\le16\), is nonzero.
The final tested order is \(5535\). This includes every choice
of algebraic-integer coefficients.

At \(p\equiv25\pmod{36}\) the shift-zero moment is nonzero
as before. At \(p\equiv7\pmod{36}\) use the coherent line
\((A_1,A_3,A_5)=(A,-2A,A)\) with \(A\ne0\).
A deficit-\(p\) term uses \(u_1,u_2\) factors of \(B,D\)
with \(u_1+2u_2=p\) and \(u_1+u_2\le m\).
The angle is at most \(7m-2u_1-6u_2=7m-2p-2u_2\).
Since \(u_2\ge p-m\), it is at most
\(9m-4p=(5p-1)/2<3p\). Thus only \(C_0,C_1,C_2\)
can contribute. When \(b=0\), shift 1 detects \(A\).
The third input occurs in practice. At \(p=151\), taking the
coefficients of \(B\) in increasing angle as
\((2,3,5,7,11,13)\) and \((c,e)=(17,19)\) gives \(C_2=97\)
by literal Laurent expansion.

For \(u=b^p\ne0\), let \(r=u^2\), \(\ell_j=b_j^p/u\),
\(v=c^p\), \(w=e^p\), and \(U_j=uC_j/A\).
If shifts 1 through 3 vanish, their triangular equations give
\[
 \begin{aligned}
 U_0&=2/31,\\
 U_1&=(16v-2\ell_1)/31,\\
 U_2&=(2\ell_1^2-16\ell_1v-2\ell_2+34w)/31
                 +(427/93)r-8908/2077.
 \end{aligned}
\]
Then shifts 4 through 6 solve successively for
\(\ell_3,\ell_4,\ell_5\). In shift \(k\) the first occurrence
of \(\ell_{k-1}\) uses this one odd factor and \(k-1\)
negative lifts. It contributes at angle zero to the \(U_0\)
column, so the pivot is a nonzero rational constant.
After the substitutions each remaining normalized row \(S_k\)
lies in \(\mathbb Q[r,\ell_1,\ell_2,v,w]\). It suffices to
find a power of \(r\) in the ideal generated by these rows.
Reversing the substitutions and clearing powers of \(u\)
then gives the desired localized identity in the original
coefficients.

For each of \(w=0,\ell_2=0,v=0,\ell_1=0\), the
[exact checker](verify_cm_newton_full_width_reconstruction_certificates.py)
specializes the ten rows \(S_7,\ldots,S_{16}\).
Give \(r,\ell_2,w\) weight two and \(\ell_1,v\) weight one.
The products with monomials up to total weight 18, 18, 20 and 20,
respectively, contain \(1\) in their rational span.
Each computation selects an invertible submatrix modulo \(65521\),
solves for rational multipliers, and checks every coefficient of the
linear combination over \(\mathbb Q\). The weights and the symmetry
restriction on \((\ell_1,v)\mapsto(-\ell_1,-v)\) bound which monomials
the search uses. The checked product gives the proof of the identity.
Clearing denominators excludes finitely many primes and proves
the first assertion. Each slice identity also lifts to an identity
\(1=\sum G_k S_k+xH\) over \(\mathbb Q\), where \(x\) is the
coordinate set to zero. Thus each of \(\ell_1,\ell_2,v,w\) is a unit
in the characteristic-zero residual quotient.

The weights also give a denominator bound for undoing the
odd-coefficient substitutions. Before elimination give \(u\)
weight one and each \(\ell_j\) weight \(j\). Each Laurent factor
of angle \(i\) has parameter weight at most \(i+1\). The resulting
normalized row has weight at most \(k-1\), because \(U_j\)
has weight at most \(j\) and its column has been divided by \(u\).
Each solved \(\ell_j\) therefore has weight at most \(j\).
Undoing the substitutions preserves the certificate's weight
bound by adding multiples of the linear pivot equations.
Replacing \(\ell_j\) by \(b_j/u\) and \(r\) by \(u^2\) then
clears the remaining parameter denominators after multiplication
by \(u^D\), where \(D=18\) or \(20\) is the corresponding bound.

Over \(\mathbb F_{331}\) the full unspecialized rows have a unit
combination with total weight at most 22. The calculation reduces
the rational rows before any primitive normalization and verifies
that the constant vector is a row of the reduced coefficient matrix.
Here \(331>13+18\cdot16\). In row \(k\le3\) the normalized
pole pivot is
\[
 \frac{(-2)^{k-1}k!}{\prod_{a=0}^{k-1}(11+18k-2a)} .
\]
The coefficient pivot for \(4\le k\le6\) has the same
formula times \(U_0=2/31\). Both are units at \(331\). This
proves the fixed-prime assertion over every residue field above \(331\).
Moreover, if the five normalized parameters are integral at
such a place, they cannot annul all \(S_k\) in characteristic
zero, since their reductions would annul a unit ideal. Some
evaluated \(S_k\) is therefore a nonzero algebraic number.
At every sufficiently large prime its Frobenius power is
still nonzero. This proves the large-prime assertion for the
local integrality class, even when one of the eliminated
coefficients has a pole above \(331\).
The rational unit-ideal question on the unrestricted
\(b_1b_2ce\ne0\) locus remains open.

### A quintic for the three-pole response

Use \(Q_+=(z^2+1)^4/z\), with the same lower support.
The change \(f_-(z,t)=i f_+(iz,t)\) gives the preceding
minus-top convention. Angular constant terms are unchanged by
the rotation, and the \(m\)-th moment acquires a factor \(i^m\).
The lower coefficients change by powers of \(i\) and retain their
valuations after adjoining \(i\). Thus the slice and local results
apply to either convention.

**Proposition (quintic obstruction for three corrections).**
On \(bb_1b_2ce\ne0\), failure to detect the coherent component
with rows \(0,\ldots,4\) is equivalent to \(F=0\), where
\[
\begin{aligned}
F={}&155708cb^5+1407b_1b^4+8880cb^3\\
 &+(3417b_1e+13362b_1+1608b_2c+201b_3)b^2\\
 &-(1608b_1^2c+402b_1b_2)b+201b_1^3 .
\end{aligned}
\]
Every line varying only \(b\) in this remaining coefficient locus
has at most five candidates for failure of the seventeen-shift
coherent test.

In the five-marker space use
\(A=(a_0+a_1z)(z^2+1)^2\) and \(C=c_0+c_1z+c_2z^2\).
Row zero is \(a_0/13\). After setting \(a_0=0\), the determinant
on columns \((a_1,c_0,c_1,c_2)\) in rows \(1,\ldots,4\) is
\[
       -\frac{4096}{433798734017054481075}\,F .
\]
The reconstruction checker multiplies the Laurent polynomials
with all coefficient variables free and verifies the identity over
\(\mathbb Q\). Rows \(1,2,3\) have an invertible correction block
for \(b\ne0\). Hence the zero locus of the determinant is exactly
where a kernel vector with \(a_1=1\) satisfies these four rows.
The nonzero leading coefficient when \(c\ne0\) gives the bound.

Put \(R_0=\mathbb Q[b_1,b_2,b_3,c,e,(b_1b_2ce)^{-1}]\) and
\(\mathscr R_0=R_0[b]/(F/(155708c))\).
This is free of rank five over \(R_0\) in the basis
\(1,b,\ldots,b^4\). Its constant coefficient is a unit, so \(b\)
is a unit too. Over \(\mathscr R_0\), rows 5 and 6 solve
successively for \(b_4,b_5\). Their pivots are units by the same
negative-angle calculation as for rows 4 through 6 above.
Thus over each fixed \((b_1,b_2,b_3,c,e)\) there are at most five
candidate triples \((b,b_4,b_5)\), counted with multiplicity.
One can test the remaining equations on this chart by taking the
five columns of multiplication by each residual in \(\mathscr R_0\).
Their image over \(R_0\) is the ideal generated by the residuals,
viewed as a submodule of \(\mathscr R_0\). At a geometric
specialization they span the fibre precisely when no candidate
satisfies every residual. On this chart \(b_3\) can alternatively be
solved from \(F\), since its coefficient \(201b^2\) is a unit.
This recovers the five-parameter presentation above.

The rank-five representation keeps the power of \(b\) below five
at the cost of expanding the Laurent coefficients in the other
variables. In a direct exact comparison using the normalized
minus-top rows, first substitute \(r=b^2,\ell_1=b_1/b,\ell_2=b_2/b\).
Clear the smallest sufficient power of \(b\) and reduce by the
monic quintic in that convention. The checker obtains the
following sizes.

| Row | Normalized terms | Power of \(b\) cleared | Degree before reduction | Terms after reduction | Largest pole order in \(c\) |
|---|---:|---:|---:|---:|---:|
| 7 | 71 | 6 | 12 | 353 | 8 |
| 16 | 1330 | 13 | 27 | 12816 | 22 |

These counts describe this presentation; they do not decide the
localized ideal.

**Corollary (generic detection for every radial value).**
There is a nonempty Zariski-open set in the base
\((b_1,b_2,b_3,c,e)\), with \(b_1b_2ce\ne0\), over which
rows \(0,\ldots,7\) detect the coherent component for every
\(b,b_4,b_5\).

Indeed, after imposing the earlier rows, let \(N_7\) be the
determinant of multiplication by the seventh residual on
\(\mathscr R_0\). The open set is \(N_7\ne0\), with \(b=0\)
already detected by row 1. To check that \(N_7\) is nonzero,
take \(b_1=b_2=c=e=1,\ b_3=0\) and reduce at 1009.
The normalized minus-top arguments of \(S_7\) are then
\((r,\ell,h,v,w)=(-b^2,i/b,-1/b,i,-1)\), with \(i=469\).
The determinants of multiplication by \(b^6S_7\) and by \(b^6\)
are \(801\) and \(336\) in the reduced quintic algebra.
Both are units. Thus the norm is nonzero in characteristic
zero. The full unit-ideal question on its zero locus is open.

### Eliminating the three angular coefficients first

There is a second finite free presentation with one fewer base
parameter. Use the normalized minus-top rows \(S_k\), after all
three angular pivots, and clear the smallest sufficient power
of \(b\) in \(S_k(b^2,b_1/b,b_2/b,c,e)\). For shift 7 the
cleared polynomial is \(G_7=b^6S_7\) in this substitution.
It has degree 12 in \(b\) with leading and constant terms
\[
 -\frac{23055471104}{584486340375375}\,b^{12},
 \qquad -\frac{2048}{64942926708375}\,b_1^6 .
\]
Thus \(G_7\) gives a rank-twelve algebra over the base with
coordinates \(b_1,b_2,c,e\) and \(b_1\) inverted. The constant
coefficient makes \(b\) a unit. Reducing the cleared row 8 and
row 16 in that algebra gives respectively 105 and 2753
coefficient monomials. This is smaller at the last row than
the tested rank-five normal form and introduces no inverse
powers of \(c\). The residual ideal still requires a global
calculation.

### Large-prime detection with an arbitrary radial coefficient

**Corollary (local integrality away from the radial coefficient).**
For either top convention, let an algebraic phase have
\(b_1,\ldots,b_5,c,e\) integral at one place above \(331\).
For every sufficiently large \(p\equiv7\bmod18\) at least one
moment of order \(m+kp\), \(0\le k\le16\), is nonzero.
Here \(b\) can be any algebraic number. In particular this
large-prime statement holds for every algebraic-integer phase
in the family.

Write \(\nu\) for the chosen valuation and use the plus-top
convention. If \(\nu(b)<0\), use the eight-row determinant
\(\Delta_8\) of the [full marker matrix](CM_NEWTON_SPHERE_ESCAPING_VARIATIONS.md#a-generic-eight-row-window).
Its degree in \(b\) is 12. All its coefficients evaluated at the
other parameters are integral, since the moment denominators for
these rows are below \(331\). The common leading coefficient is
\(5\bmod331\) by the same leading-entry formula. The degree-12
term has uniquely smallest valuation, so \(\Delta_8\ne0\).

Suppose now that \(\nu(b)\ge0\). The slice result handles
\(b=0\). Choose a coefficient \(\beta\) of \(B\) of minimum
valuation. If \(\nu(b)=\nu(\beta)\), then \(b^2,b_1/b,b_2/b,c,e\)
are integral, and the normalized-integrality result applies.
Otherwise let \(q\ge1\) be the first index with
\(b_q/\beta\) a unit. In the five-marker matrix, divide each
correction column by \(\beta\). Every term in such a column
has an odd number \(d\) of odd factors. Its coefficient includes
\(\beta^{d-1}\) and \(d\) factors from \(B/\beta\), so it is
integral. Modulo the place the contributing odd factors have
angle at least \(q\); the other angles are at least \(-1\).
For column \(C_j\) at row \(k\), a contribution therefore needs
\(d(q+1)\le k-j\). At \(k=q+j+1\), its first possible term
uses one angle-\(q\) factor and \(k-1\) negative lifts.
Its coefficient is the single diagonal entry
\[
  \frac{b_q}{\beta}\,
  \frac{(-2)^{k-1}k!}{\prod_{a=0}^{k-1}(11+18k-2a)} .
\]
This is a unit. Row 1 has no correction term after reduction
and its \(a_1\)-entry is \(-2/(31\cdot29)\).
Together with row zero and the correction rows \(q+1,q+2,q+3\)
this gives an invertible normalized minor modulo the place.
All these rows have index at most eight.

In each minor argument the evaluated characteristic-zero minor is
nonzero. It remains nonzero after reduction at all sufficiently
large primes, and entrywise Frobenius preserves its rank. At the
difficult primes, the coherent input lies in the indicated
\(a_1\)-line and is nonzero. The full rank on this line and the
three correction columns excludes simultaneous vanishing, and
shift zero handles the other prime class. The fixed-prime
statement above applies to phases integral in every coefficient
at the chosen place.

### A single normalized pole at 331

**Proposition (single normalized poles in the three-pole system).**
At every algebraic common zero of the rational residuals,
at least two of \(r,\ell,h,v,w\) have negative valuation at
each place above \(331\). Thus the seventeen-shift theorem
applies if at most one normalized parameter has a pole at
some place above \(331\). The sign rotation preserves the
number of poles.

All the coefficients of the residuals \(S_7,\ldots,S_{16}\)
are integral at \(331\). For a pole in one coordinate divide
each equation by its largest characteristic-zero power of that
coordinate before reducing. The resulting leading coefficients
must vanish. In \(r,\ell,h\), the leading coefficients of \(S_7\)
are respectively \(17,-121,-121\bmod331\). For the leading
coefficients \(L_{v,k}\) in \(v\) there is the reduced identity
\(-102L_{v,7}-41L_{v,8}=r\).

For the leading coefficients in \(w\) one has
\[
 56 L_{w,7}+96v L_{w,8}-120 L_{w,9}=r,\qquad
 L_{w,8}=81\ell,\qquad L_{w,7}|_{r=\ell=0}=117(h-137)
 \quad\bmod331.
\]
Thus \(r,\ell,h\) reduce to \(0,0,137\), respectively.
Set \(s=1/w\), and write \(T_k=s^{\deg_w S_k}S_k(r,\ell,h,v,1/s)\).
The checker verifies at \(s=0\), with \(h_0=1466/1139\) lifting
\(137\), that the Jacobian of \(T_7,T_8,T_9\) in \(r,\ell,h\)
at \((0,0,h_0)\) has constant unit determinant. It also verifies
the exact degree-two solution jet with those initial coordinates.
On this jet, \(T_{11}\) has zero constant and linear coefficients,
and its \(s^2\)-coefficient is
\[
       -\frac{56321516765184}{1239568379329485639766985107},
\]
a unit at \(331\).

To apply the jet at a proposed algebraic zero, evaluate its integral
coefficients at the given \(v\). The true coordinates and the
evaluated jet have the same reduction. If their difference has
minimum coordinate valuation \(a>0\), the unit Jacobian makes the
linear part of the difference in the three pivot equations have
minimum valuation \(a\). Its nonlinear terms have valuation
greater than \(a\). The pivot equations at the evaluated jet vanish
to order \(s^3\), so \(a\ge3\nu(s)\). The same error bound holds
for \(T_{11}\) by integrality of its coefficients. Its quadratic
unit term therefore gives a contradiction.

For a pole in \(v\), the second leading identity is
\(L_{v,7}|_{r=0}=109(\ell^2-h-34)\bmod331\).
Use \(s=1/v\) and normalize each row by its degree in \(v\).
The starting coordinates for \(r,h\) are
\(0,\ell^2+370/67\), giving the required reduction. The pivot
rows 7 and 8 have constant unit Jacobian in \(r,h\) there.
The checker solves their equations recursively to order \(s^6\),
over the ring with the two integral parameters \(\ell,w\).
The same Taylor comparison bounds the error valuation at a
proposed zero below by \(7\nu(s)\).

On this jet the normalized rows 9, 10 and 11 start at order
\(s^4\). Write \(p(\ell,w)\) for the coefficient of \(s^4\) in
row 9. Its leading coefficient is a unit. Successively at
orders \(4,5,6\) divide the current leading coefficient of
each of rows 10 and 11 by \(p\), subtracting the quotient
times the appropriate power of \(s\) times row 9. The
remainders at orders 4 and 5 vanish exactly. At order 6 the
two remainders reduce to \(92w+121\) and \(244w+243\).
The quotients are integral. These linear polynomials have
distinct roots, so the two rows cannot both vanish.
This excludes the last single-pole case.

### The projective boundary at 331

Use the characteristic-zero parameters
\((r,\ell,h,v,w)=(b^2,b_1/b,b_2/b,c,e)\) in \(S_7,\ldots,S_{16}\).
All coefficients of these ten rows are \(331\)-integral.
Homogenization uses the degrees over \(\mathbb Q\) before
reducing coefficients. For ordinary homogenization the
geometric boundary over \(\mathbb F_{331}\) is the union
\(\ell=0\) and \(\ell=8v\) in \(\mathbb P^4\).
For weights \((2,1,2,1,2)\) the geometric boundary is the union
of the following three weighted projective curves:
\[
 \begin{array}{c|l}
 L_1 & r=\ell=h=0\\
 L_2 & r=h=0,\quad \ell^2-8\ell v+17w=0\\
 L_3 & r=0,\quad \ell=8v,\quad h=17w .
 \end{array}
\]
The [boundary checker](verify_cm_newton_full_width_projective_boundary.py)
substitutes each component in the corresponding rational leading
forms. Thus these boundary components occur already over \(\mathbb Q\).
For the reverse inclusion in characteristic 331, it verifies
membership of \((\ell(\ell-8v))^4\) in the ordinary leading ideal.
For the weighted ideal it checks membership of \(r^9\) and the
seventh powers of
\[
 h(h-17w),\qquad h(\ell-8v),\qquad
 (\ell^2+h)(\ell-8v)-\ell(h-17w).
\]
Splitting according to \(h=0\) and then \(\ell=0\) gives the
three listed components over the algebraic closure. In particular
neither unsaturated homogeneous system has empty boundary.

**Proposition (local restrictions on a residual zero).**
Suppose an algebraic point annuls \(S_7,\ldots,S_{16}\), and
let \(\nu\) be a valuation above \(331\), normalized by \(\nu(331)=1\).
Put
\[
 \delta=\min\{\nu(r)/2,\nu(\ell),\nu(h)/2,\nu(v),\nu(w)/2\}.
\]
Then \(\delta\le-1/22\) and \(\nu(r)>2\delta\). Its weighted
reduction cannot lie on \(L_1\). A point on \(L_2\) with
\(\ell\ne0\) must have
\[
             v/\ell\in\{2/17,\,5/56,\,1/8\}
             \quad\text{in characteristic }331 .
\]
Thus any reduction outside \(L_3\) is among the two points on
\(L_2\) with \(v/\ell=2/17\) or \(5/56\).

To see the valuation bound, lift the coefficient multipliers of
the weight-22 unit identity over \(\mathbb F_{331}\) to
\(\mathbb Z_{(331)}\). The weighted leading forms of the ten
rows stay nonzero modulo 331, so the lifted combination still
has weight at most 22 and equals \(1+331H\) with \(H\)
\(331\)-integral. At a zero of the residuals,
\(\nu(H)=-1\). Each monomial in \(H\) has valuation at least
\(22\delta\) if \(\delta<0\), and valuation at least zero
otherwise. This gives the bound on \(\delta\).
In the valued field choose \(\lambda\) of valuation
\(-\delta\), extending the field if a square root is needed.
Scale the coordinates by \(\lambda\) to their respective
weights. The scaled coordinates are integral and at least one
is a unit. Reducing the homogenized equations with \(T=\lambda\)
puts the point on the weighted boundary, so \(\bar r=0\).

For the \(L_1\) exclusion with \(w\ne0\), normalize \(w=1\)
by a further unit scaling and put \(\tau=T^2\). In this chart
the transverse coordinates \(r,\ell,h\) and \(\tau\) are in
the maximal ideal. The checker solves the pivot equations
from shifts \((7,8,9)\) and \((7,8,10)\) through order three in
\(\tau\), over \(\mathbb Q(v)\). The two pivot determinants
reduce to coprime polynomials. At least one triple therefore
has unit determinant at the residue value of \(v\). For this
triple the actual transverse coordinates agree with the computed
expansion modulo \(\tau^4\). Indeed the difference of the two
pivot evaluations factors through the coordinate differences,
with matrix reducing to the invertible pivot Jacobian.
The constant, linear and quadratic residual terms vanish exactly.
The reduced cubic numerators have gcd 1, with denominators
units at every point where the pivot determinant is a unit.
Some residual consequently has a unit cubic term and cannot
vanish.

At the remaining point of \(L_1\), normalize \(v=1\) and put
\(d=\ell^2-8\ell+17w\), \(\tau=T^2\).
Shifts 7 and 8 solve for \(r,h\) with unit Jacobian.
For weights \((1,1,2)\) on \((\ell,d,\tau)\), the checker
uses exact jets through weight ten and row combinations with
\(331\)-integral multipliers. The radical of their reduced
initial ideal is \((\ell d,\tau)\). If
\(\rho=\min\{\nu(\ell),\nu(d),\nu(\tau)/2\}\), scaling by
an element of valuation \(\rho\) leaves just two possible
directions, the \(d\)-axis and the \(\ell\)-axis.

On the \(d\)-axis put \(\ell=d L\), \(\tau=d^2 Z\), with
\(L,Z\) in the maximal ideal. The equation \(S_9/d^6=0\) on the
implicit branch has unit derivative in \(L\), and its solution
has \(L(d,0)=0\). The weight-seven combination with reduced
leading form \((8\ell+d)\tau^3\), divided by \(d^7\) and
evaluated on this branch, has unit \(Z^3\) coefficient.
The first three coefficients in \(Z\) vanish exactly. To see
this, work first with \(d\) invertible and \(\tau=0\), where
\(w=d/17\). Over a field extension one may normalize \(w=1\).
The shifts \(7,8,9\) implicit solution in that chart has all
residuals of order at least three in \(\tau\). The two solutions
agree after a change of chart by uniqueness of the formal
implicit solution. The order statement therefore holds
identically in \(d\), including for the integral power series
on the \(d\)-axis. Because \(\tau\ne0\), the unit cubic
coefficient excludes this direction.

On the \(\ell\)-axis put \(d=\ell D_0,\ \tau=\ell^2 Z\)
with \(D_0,Z\) in the maximal ideal. In the
\((d,\tau)\)-jets of weights \((1,2)\), shift 9 begins
in weight five and its leading coefficients are divisible
by \(\ell\). Subtracting the exact multiples of this row
that cancel weight five in the other rows leaves weight-six
coefficients also divisible by \(\ell\). The checker verifies
integrality at the chart point for these divisions. After the
change of variables, divide shift 9 by \(\ell^6\) and the
other combinations by \(\ell^7\). Their lowest-weight forms
in \((D_0,Z)\) reduce to an ideal with radical \((D_0,Z)\).
A second scaling by \(\min\{\nu(D_0),\nu(Z)/2\}\) contradicts
that radical. This excludes the last point on \(L_1\).

For \(L_2\) with \(\ell\ne0\), normalize \(\ell=1\) and set
\(d=17w+1-8v\). The checker solves shifts 7 and 8 for \(r,d\)
in terms of \(h,\tau\). If \(J\) is their transverse determinant,
write \(h=J^2 x,\ \tau=J^2 y\). Outside the three listed
residues, \(J\) is a unit. The checker obtains the exact
expansions through weight six in \((x,y)\) of weights \((1,2)\).
The leading terms of shifts 9 and 10 after one linear cancellation
are \(x^5 Q(y/x^2)\) and \(x^6 H_0(y/x^2)\), of weights
five and six. Here \(Q\) has degree two, \(H_0\) has degree
three, and \(Q(0)\), both leading coefficients and the resultant
are units outside the three residues. All lower-weight terms
vanish over \(\mathbb Q(v)\). The truncation gives the actual
coordinates to valuation at least \(7\min\{\nu(x),\nu(y)/2\}\)
by the same unit-Jacobian difference argument.
If \(\nu(y)\le\nu(x)\), the \(y^3\) term has strictly smallest
valuation in the weight-six row. If \(\nu(x)<\nu(y)<2\nu(x)\),
the \(xy^2\) term has strictly smallest valuation in the
weight-five row. If \(\nu(y)>2\nu(x)\), the \(x^5\) term does.
The case \(x=0\) also gives a nonzero \(y^3\) term. The remaining
ratio \(\nu(y)=2\nu(x)\) would force a common root of \(Q,H_0\)
over the residue field extension, contradicting the resultant.

The [refined off-component charts](CM_NEWTON_331_OFF_COMPONENT_CHARTS.md)
exclude both directions with \(v-\alpha\) of strictly least
valuation. In the other chart each of the two points off
\(L_3\) is confined to one further residue of a transverse
ratio.

These calculations replay with SageMath 10.9.

```sh
sage -python research/foundations/extended-geometry/verify_cm_newton_full_width_projective_boundary.py
```

### Varying the conductor-18 branch ratio

Let \(P_\lambda=(z-1)^4(z-\lambda)^4\), with \(\lambda\ne0,1\).
Keep \(\deg B\le5\) and \(D=-\lambda^4 z^{-1}+c+ez\).
Put \(Q_s=(z-1)^2(z-\lambda)^2=\sum_{i=0}^4q_i z^i\).
In good characteristic \(p\equiv7\bmod18\), factor
\[
 P_\lambda^m=Q_s(z)^p(z-1)^\ell(z-\lambda)^\ell,\qquad
                 \ell=4m-2p=(8p-2)/9.
\]
The last two factors have degree \(2\ell<2p\).
Writing \(X,Y\) for their coefficients at \(m,m+p\) gives
\(A_j=q_j^p X+q_{j-1}^p Y\), with out-of-range \(q_i=0\).
Since \(q_0=\lambda^2\), the zero shift forces \(X=0\) when its
residue vanishes. The full coherent source condition gives
\(Y\ne0\). The coherent component is then proportional to
\((q_0^p,\ldots,q_4^p)\) at indices \(1,\ldots,5\), with
\[
      (q_0,\ldots,q_4)=
      (\lambda^2,-2\lambda(\lambda+1),\lambda^2+4\lambda+1,
       -2(\lambda+1),1).
\]
The angle bound above still leaves only \(C_0,C_1,C_2\).

**Proposition (moving branch ratios above 331).**
For any fixed algebraic coefficients with \(b=B(0)=0\),
shifts \(0,1\) detect the coherent input for every sufficiently
large \(p\equiv7\bmod18\) and every \(\lambda\ne0,1\).
For \(b\ne0\), suppose that at one place above 331,
\(\lambda\) reduces to \(-1\) and
\(b^2,b_1/b,b_2/b,c,e\) are integral. Then shifts \(0\)
through \(16\) detect the coherent input for every sufficiently
large \(p\equiv7\bmod18\). If all phase coefficients are integral
at such a place, the window \(239+331k\), \(0\le k\le16\),
also gives a nonzero moment at the fixed prime 331.

If shift zero vanishes and \(b=0\), the coherent pivot in
shift 1 is \(-2(\lambda^4)^p A_1/899\). The two possible
secondary terms in that shift vanish because the negative
lift is even and the odd radial coefficient is zero.
For \(b\ne0\), form the characteristic-zero rows with coherent
profile \(zQ_s\). Set \(U_j=bC_j/Y\), \(r=b^2\) and
\(\ell_j=b_j/b\). The pivots for \(U_0,U_1,U_2\) in rows 1
through 3 are the fixed-top pivots above multiplied by
\((\lambda^4)^{k-1}\). In particular \(U_0=2\lambda^6/31\).
The pivots for \(\ell_3,\ell_4,\ell_5\) in rows 4 through 6
are those first three factors with the respective shift index
\(k\), multiplied by \(U_0\). After these substitutions the
remaining rows are polynomials in \(r,\ell_1,\ell_2,c,e\),
with rational coefficients in \(\lambda\) obtained by inverting
only the pivot scalars. The divisibility and weight argument
for the fixed top applies with \(\lambda\) a scalar parameter.
At the chosen place the pivots are units. Reduction gives
exactly the fixed-top rows over \(\mathbb F_{331}\), whose
unit identity forces some evaluated residual to be nonzero.
At sufficiently large good primes that residual remains
nonzero, also after entrywise Frobenius. For the fixed-prime
statement the phase itself reduces to the fixed-top family,
and the normalized residues agree with the same seventeen tests.

The [exact moving-branch certificate](CM_NEWTON_MOVING_BRANCH_RATIO_CERTIFICATES.md)
covers every branch ratio on a smaller lower support.

For a degree-19 phase in these families, the existing two-pair
Hopf homogenization \(F\) satisfies
\(\mathcal E_2(F^j)=(19j+1)!\mathcal A(f^j)\).
Each unmarked nonvanishing statement therefore transfers to that
balanced Image family.
[Hopf normalization, (13)](PRIMITIVE_IMAGE_DEGREE_SEVEN_AND_HOPF_SUBALGEBRA.md).

An exact rational certificate for [linear odd and radial even lower
terms](CM_NEWTON_MOVING_BRANCH_RATIO_CERTIFICATES.md) covers every
branch ratio in that smaller support.

### Hamiltonian symmetries and the coherent profiles

Over a characteristic-zero field \(K\), let \(\mathscr R\) be the
coordinate ring of the affine quadric. Its regular volume form is
\(\omega=dz/z\wedge dt\) in this chart. Set
\[
             \{f,V\}=z(f_zV_t-f_tV_z).
\]
The markers have the form \(H=t^{s-1}A(z)+t^{s-2}C(z)\),
with \(\deg A\le(s-1)/2\) and \(\deg C\le(s-3)/2\).
The leading \(t\)-degree in this chart is the weighted degree
for \(w(x)=2,w(t)=1,w(y)=0\) in the normal sphere basis. A
negative-angle monomial \(y^j t^d\) has leading \(t\)-degree
\(d\) after substituting \(y=(1-t^2)/(2zt^2)\), and distinct
angular coefficients cannot cancel.
For regular markers \(h_i\), the reconstruction matrix is the
Jacobian of the unmarked moment map in the marker directions,
after scaling row \(k\) by \(1/(k+1)\).
\[
 \mathcal A(h_i f^k)=
 \frac1{k+1}
 \left.\frac{\partial}{\partial\varepsilon_i}
 \mathcal A\left((f+\sum_j\varepsilon_j h_j)^{k+1}\right)
 \right|_{\varepsilon=0}.
\]
Expanding the first-order term proves the identity. Higher
moment derivatives concern nonlinear deformations. The kernel
criterion for residue reconstruction uses these first derivatives.

For regular \(V\) and \(k\ge0\),
\[
        f^k\{f,V\}\omega=\frac{d(f^{k+1}dV)}{k+1},
\]
so these brackets have all sphere moments zero.

**Proposition (Hamiltonian exclusion in the marker space).**
If \(V\in\mathscr R\) and the \(t\)-degree of \(\{f,V\}\) in the
\((z,t)\) chart is less than \(n-1\), then \(V\in K[f]\) and
\(\{f,V\}=0\). In particular no nonzero marker is such a bracket.

Let \(G(z)t^\mu\) be the leading term of \(V\) in \(t\).
For each negative angle \(-j\) in the \((x,t)\) chart, its coefficient
is divisible by \((1-t^2)^j\). After \(x=zt^2\) the leading degree
of each nonzero angular coefficient is therefore nonnegative.
The coefficient of \(t^{n+\mu-1}\) in the bracket is
\(z(\mu Q'G-nQG')\), which must vanish. Thus \(G^n/Q^\mu\)
is constant. The simple pole of \(Q\) at zero forces \(n\mid\mu\),
and the differential identity gives \(G=\lambda Q^{\mu/n}\)
with \(\lambda\in K\). Subtract \(\lambda f^{\mu/n}\) and repeat.
Each step decreases the nonnegative top degree and leaves the
bracket unchanged.

The rational centralizer is \(K(f)\). A rational function
commuting with \(f\) is algebraic over \(K(f)\) by the differential
criterion for algebraic dependence. Multiplying by a suitable
nonzero polynomial in \(f\) makes it integral over \(K[f]\).
The quadric is smooth and its coordinate ring is integrally
closed, so this multiple is regular. The preceding argument
puts it in \(K[f]\).

For the coherent component the exclusion persists on the generic
fibre. Write the nonzero roots of \(P\) as \(\alpha_i\), with the
multiplicities \(m_i\) used above.

**Proposition (coherent profiles on the generic fibre).**
Suppose \(A\ne0\) obeys
\[
       \deg A\le\left\lfloor\frac{as}{n}\right\rfloor,\qquad
       \operatorname{ord}_{\alpha_i} A
             \ge\left\lfloor\frac{sm_i}{n}\right\rfloor
             \quad(1\le i\le R).
\]
For any nuisance polynomial \(C\) in the marker space, the
relative differential \(H\omega/df\) on the smooth affine generic
fibre of \(f\) has nonzero algebraic de Rham class.

The generic fibre is smooth in characteristic zero. On each
component of the critical locus, \(df=0\) forces \(f\) to be
constant, so its critical values form a finite set.
Exactness on the affine fibre would give a regular primitive,
since the affine de Rham complex is computed by global sections.
The fibre's coordinate ring is \(\mathscr R\otimes_{K[f]}K(f)\).
Clearing its denominator gives \(\{f,V\}=J(f)H\) with
\(V\in\mathscr R\) and \(J\) monic of degree \(d\). The case
\(d=0\) is already excluded. For \(d\ge1\), strip off centralizing
powers of \(f\) whenever the leading bracket term would be above
degree \(nd+s-1\). The remaining leading term of \(V\) has
degree \(\mu=n(d-1)+s\), since the bracket has degree exactly
\(nd+s-1\). Its coefficient \(G\) is Laurent in \(z\). With
\(T=G/Q^{d-1}\), the leading equation becomes
\[
               sz(Q'/Q)T-nzT'=A.
\]
At \(\alpha_i\), for \(T\) of order \(u\), the leading coefficient
on the left is proportional to \(sm_i-nu\).
It cannot vanish since \(s\) is a unit modulo \(n\) and
\(0<m_i<n\). Thus \(\operatorname{ord}_{\alpha_i} T
=\operatorname{ord}_{\alpha_i} A+1\).
At zero the multiplier is \(-s-nu\), and the equation excludes
a pole there as well. These are the only possible finite poles,
so \(T\) is a polynomial. At infinity the multiplier is
\(sa-n\deg T\ne0\) since \(0<a<n\), so \(\deg T=\deg A\).
But
\[
   \deg T\ \ge\ \sum_i\left\lfloor\frac{sm_i}{n}\right\rfloor+R
   \ >\ \left\lfloor\frac{as}{n}\right\rfloor .
\]
The strict inequality follows from \(\sum_i m_i=a+1\) and
the sum of the \(R\) fractional parts of \(sm_i/n\) being
less than \(R\).

The profile conditions express holomorphy of
\(Y^{-s}A(z)\,dz/z\) on the normalized cyclic cover \(Y^n=Q(z)\).
At a root of multiplicity \(m_i\) the order is
\((n(\operatorname{ord}_{\alpha_i} A+1)-sm_i)/\gcd(n,m_i)-1\).
At infinity the order is \((as-n\deg A)/\gcd(n,a)-1\),
and at zero it is \(s-1+n\operatorname{ord}_0 A\).
In good characteristic \(p\equiv s^{-1}\pmod n\), Cartier
applied to \(Y^{-1}dz/z=Y^{-sp}Q^m\,dz/z\) gives
\[
        Y^{-s}\sum_j A_j^{1/p}z^j\,\frac{dz}{z},
\]
which has this profile by preservation of holomorphy.
For conductor 18 the characteristic-zero coherent line is
\(A(z)\propto z(z^2-1)^2\). The degree bound is five and
each of \(\pm1\) requires order at least two. The generic-fibre
result applies to this line for every compatible lower strip.

The remaining period question is whether the sphere resolvent
\(\mathcal A(H/(\lambda-f))\) detects these coherent classes.
Its Laurent coefficients at infinity are the marker moments,
so their simultaneous vanishing gives the zero germ.
On the open parameter set where the fibre complements form
a locally trivial family, monodromy differences of the sphere cycle
lie in the kernel of the map from the homology of the fibre
complement to that of the quadric. By the relative homology
sequence and the Thom isomorphism this kernel consists of tubes
around fibre cycles. Integrating \(H\omega/(\lambda-f)\) on a tube
gives the corresponding period of \(H\omega/df\), up to sign and
a factor of \(2\pi i\). It would suffice to show that the periods
on the cycles obtained from the sphere distinguish a nonzero
coherent component in the presence of \(C\).

Gavrilov's result for a real vanishing cycle assumes a semiweighted
homogeneous plane Hamiltonian with Morse critical points and gives
the rank from the span of that cycle under monodromy. The monodromy
statement of Novikov and Yakovenko in §1.5 concerns generic Morse
plane Hamiltonians regular at infinity. The quadric and its fixed
sphere cycle require a separate argument.
[Gavrilov, Theorem 2 and Corollary 2](https://www.numdam.org/item/10.5802/aif.1684.pdf)
and [Novikov and Yakovenko, §1.5](https://pi.math.cornell.edu/~novikov/Papers/redundant.pdf).

The nuisance marker \(t^{11}\) illustrates the period issue
when \(n=18,s=13,B=0\). Reflection in \(t\) fixes \(f\) and
kills all sphere moments against this marker. Its generic-fibre
class is nonzero. The leading equation for a putative primitive
has \(s\) replaced by 12 and right side 1. Its multipliers at
the two roots, at zero and at infinity are nonzero for integral
orders. It forces a polynomial \(T\) of degree zero with a zero
at each of \(\pm1\).


The [sphere period calculation](CM_NEWTON_SPHERE_PERIOD_CYCLES.md)
identifies a relative path for the sphere resolvent and computes
two limiting sphere variations on the separable conductor-18 model.

The same cyclic-cover factorization yields [one-step rank bounds
and a rank-two block measured by marked moments](CM_NEWTON_BLOCK_RANK_MOMENTS.md).

## Open reconstruction question

The [seventeen-row locus calculation](CM_NEWTON_SEVENTEEN_ROW_LOCUS.md)
projects the common-zero equations through smaller finite algebras.
The remaining normalized locus has dimension at most one. In
particular the seventeen shifts suffice whenever the two even
coefficients are integral at one place above \(331\), with all
odd coefficients arbitrary.

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
