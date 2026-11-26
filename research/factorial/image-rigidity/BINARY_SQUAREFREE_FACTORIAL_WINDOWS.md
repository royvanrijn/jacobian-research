# Finite factorial windows for squarefree binary lowest parts

**Theorem.** Let \(f(X,Y)\) over a characteristic-zero field have
squarefree lowest homogeneous part \(F\) of degree \(d\ge2\).
Suppose every higher-degree term has degree at least \(2d-1\).
Then the factorial moments of \(f\) are nonzero infinitely often.
For each fixed phase there is a finite length such that every
window of that length at positive indices contains a nonzero
factorial moment. The length can depend on the phase.

This includes every squarefree homogeneous binary form of
degree at least two, and arbitrary higher terms when the
lowest part is squarefree quadratic.

## The binary coefficient matrix is invertible at good primes

Work first over a finite residue field \(k\) of characteristic
\(p\) where the reduction of \(F\) is squarefree of degree \(d\)
as a homogeneous form. We use the same symbols for the reduced
form and put \(P(x)=F(x,1)\).
The matrix on the \(d-1\) positive binary indices of total
degree \(d\) from the
[Hasse–Witt factorial witness](FACTORIAL_HASSE_WITT_NATIVE_ORDERS.md#the-finite-frobenius-test)
takes the form
\[
                     A_{t,r}=[x^{pt-r}]P(x)^{p-1},
                       \qquad 1\le r,t\le d-1 .
\]
We prove it invertible, including when the binary form has
a root at zero or infinity.

Pass to a finite splitting extension of \(k\) and let \(Z\)
be the set of its \(d\) distinct projective zeros.
Consider the Laurent quotient
\[
              L=k[x,x^{-1}]/(k[x]+x^{-d}k[x^{-1}])
\]
over that extension. To map values specified at \(Z\) into
\(L\), interpolate them by \(g_0(x)\) on the finite \(x\)
chart and by \(g_\infty(y)\), \(y=1/x\), on the finite \(y\)
chart. Then
\[
                 \frac{g_0(x)-g_\infty(1/x)}{P(x)}
\]
is Laurent, because the numerator vanishes at each common
point with \(x\ne0,\infty\) and those roots are simple. A
different choice of \(g_0\) adds a multiple of \(P\).
A different choice of \(g_\infty\) adds a multiple of
\(y^dP(1/y)\), so the change of the quotient is in
\(x^{-d}k[x^{-1}]\). This proves well-definedness. If the
class is zero, the two interpolants can be adjusted by
these multiples to agree on the overlap. Their common
Laurent polynomial is constant because it belongs to both
\(k[x]\) and \(k[x^{-1}]\). Conversely a constant tuple
gives the zero class. Thus the quotient of the \(d\)
values by constant tuples injects into \(L\). The two
spaces have dimension \(d-1\), so the map is an isomorphism.

The \(p\)-th power on each value corresponds under this
map to \(u\mapsto P(x)^{p-1}u^p\) on Laurent classes.
In the basis \(x^{-1},\ldots,x^{-(d-1)}\), its matrix
is \(A^{\mathsf T}\) acting after coefficient Frobenius.
The coordinatewise \(p\)-th power is invertible over
the finite splitting field and preserves the constants.
It follows that \(A\) is invertible there, and hence over
the original residue field.

## The consequence for moments

For number-field coefficients, choose any good prime
\(p>\deg f\) avoiding denominators and the binary
discriminant, and use its finite residue field.
The gap from the lowest degree to the next degree
is greater than \(d-2\). In the native-order witness
at \(p^s-1\), this makes the first matrix square with
only bottom-degree contributions. Invertibility gives
infinitely many unit residues of the normalized moments.
The fixed-phase polynomial recurrence then gives a
finite-width window.

For arbitrary characteristic-zero coefficients, eventual
vanishing would give a solution to the polynomial
equations for a zero tail of factorial moments at the
fixed support. Adjoin an inverse of the discriminant of
the bottom binary form. Noetherianity and the weak
Nullstellensatz give a solution over \(\overline{\mathbb Q}\),
with squarefree bottom form and the same degree-gap
condition. It would violate the number-field conclusion.
Thus the moments of the original phase are not eventually
zero, and the same recurrence argument applies.
