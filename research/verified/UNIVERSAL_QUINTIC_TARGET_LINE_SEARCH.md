# Searching for one target line with all five quintic groups

The fixed determinant-minus-two map has inverse polynomial

\[
E_{\Pi,B,C}(S)=
\Pi^5S^5-5\Pi S^3-2BS^2+4S-2C.
\]

The next search asks for one rational affine line

\[
(\Pi,B,C)=u+tv
\]

with rational parameters realizing \(C_5,D_5,F_{20},A_5,S_5\).  This note
records a bounded first incidence search.  It does **not** claim that no such
line exists.

## Completed bounded search

The anchor catalogue exhausts primitive projective targets through height
\(30\), modulo \((B,C)\mapsto(-B,-C)\), using the exact discriminant and
Frobenius screens from the height search and PARI/GP only to classify the
survivors.  After expanding the sign involution, it contains

\[
8\ A_5,\qquad 2\ C_5,\qquad 4\ D_5,\qquad 2\ F_{20}
\]

affine points.  The only new exceptional target between heights \(22\) and
\(30\) is

\[
(\Pi,B,C)=(-1,27,28),
\]

with group \(D_5\).

For every affine line through two catalogue points of different exceptional
groups, rational parameters \(t=p/q\) are screened by exact arithmetic.
The completed prototype covered all \(81\) distinct such lines through parameter
height \(120\).  Lines involving the new height-\(30\) dihedral point were
additionally covered through parameter height \(300\).  No line acquired a
third exceptional group.  Since an \(S_5\) value cannot repair two missing
exceptional groups, this excludes a five-group witness only inside this
finite search region.

There is a clean three-group partial line:

\[
(\Pi,B,C)=
\left(1,0,-\frac7{10}\right)
+t\left(-\frac35,-\frac{21}{10},\frac{27}{10}\right).
\]

It realizes \(C_5,D_5,S_5\) at \(t=0,1,-1\), respectively.  The first two
are the existing exact rows.  At \(t=-1\), the primitive normalized
polynomial is

\[
15625T^5-78125T^3-105000T^2+256000T+1114112;
\]

the patterns \(3:(5)\) and \(19:(2,1,1,1)\) certify \(S_5\) without a
Galois-group oracle.  A deeper exact screen of this line through rational
parameter height \(600\) found no \(A_5\) or \(F_{20}\) value.

The reproducible search driver is

```bash
.venv/bin/python scripts/search_universal_quintic_target_lines.py \
  --parameter-bound 120 --screen-primes 12 --write
```

The twelve primes are only sound exclusion filters.  PARI is a discovery
classifier here, not a witness verifier.  Any positive line must receive a
separate exact certificate without `polgalois`.

## Cost function for a positive witness

To make comparisons reproducible, encode every rational as a reduced signed
numerator and positive denominator.  For each candidate line, normalize two
of the five witness parameters to \(0\) and \(1\), try all ten choices, and
minimize lexicographically:

1. maximum bit length among the coordinates of \(u,v,t_1,\ldots,t_5\);
2. total bit length of those reduced rationals;
3. largest prime used by the modular witnesses;
4. total degree of the displayed resolvents;
5. UTF-8 byte length of canonical minified JSON containing the exact
   certificate, including resolvent factors or automorphisms but excluding
   prose.

There is a useful certificate improvement independent of the line search.
For the height-\(21\) \(D_5\) row, Dummit's sextic factors as a linear
polynomial times a quintic irreducible modulo \(3\).  Together with the
square discriminant and the \(11:(2,2,1)\) pattern, this is a three-check
dihedral certificate of resolvent degree \(6\).  Using Dummit's sextic for
both \(D_5\) and \(F_{20}\) lowers the current total resolvent-degree target
from \(16\) to \(12\).

## Next search

The bottleneck is anchor diversity, not parameter height.  The constant
\((\Pi,B)=(1,0)\) line already contains the small \(A_5\) and \(C_5\)
targets, but an exact rational scan through height \(600\) found neither
\(D_5\) nor \(F_{20}\).  Likewise, 2,644,128 primitive trace vectors in the
known De Moivre \(F_{20}\) field produced only its already-known target
presentation.

The next useful computation should therefore generate new presentations of
the cyclic and Frobenius fields on the trace-descent variety, then hash
lines between those presentations and the dihedral/alternating catalogues.
Merely extending \(t\) on the same \(81\) lines has lower expected value.
