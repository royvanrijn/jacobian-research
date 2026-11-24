import Mathlib.Algebra.Polynomial.Laurent

namespace LaurentPolynomial

abbrev support {R : Type*} [Semiring R] (p : LaurentPolynomial R) := p.coeff.support

end LaurentPolynomial

