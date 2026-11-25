# Coherent input moments on the full even-conductor strip

The wider even-conductor strip allows all three top Laurent weights
\(n,n-1,n-2\). That exact bound is enough to recover the full coherent
source by unmarked moments when the source index is at most half the
conductor. At larger indices there is one further interval pole for
each angular input. Its contribution is coupled to odd-height terms
and can cancel a nonzero coherent input.

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
residues for arbitrary subleading terms at upper-half indices. It
does not decide whether additional moments suffice there.
[Exact replay of both cases](verify_cm_newton_full_width_strip_moments.py)
