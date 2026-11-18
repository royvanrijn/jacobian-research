# Literature check and independent adversarial audit

17 November 2026. This records a bounded literature search and a second
reading of the derivation; it is not external peer review or a proof of
novelty.

## Literature actually located

1. **Dietmann and Elsholtz**, *Sums of two numbers having only prime factors
   congruent to one modulo four*, Banach Center Publications 118 (2019),
   177--179, DOI 10.4064/bc118-11.
   https://www.impan.pl/shop/en/publication/transaction/download/product/113105

   This proves the Gaussian result for all sufficiently large admissible
   integers, using Greaves. It also provides many representations, so
   merely asking for multiple representations is not a new objective.
   Its final remarks discuss effectiveness and the different problem of
   summands supported on primes 3 modulo 4. It does not state the
   Eisenstein two-summand theorem.

2. **Olivier Ramaré and G. K. Viswanadham**, *The trigonometric polynomial on
   sums of two squares, an additive problem and generalization*,
   Transactions of the AMS 379(5) (2026), 3425--3465,
   DOI 10.1090/tran/9598; arXiv:2509.23260.
   Author PDF: https://ramare-olivier.github.io/Maths/B2Vinogradov-Arxiv.pdf
   Author bibliography: https://ramare-olivier.github.io/Maths/contenu.html

   The introduction, page 3, explicitly includes primitive Löschian
   integers obtained by sieving integers 1 modulo 3 by primes 2 modulo 3.
   This is the set A3 used here (with the condition 3 does not divide the
   integer). The additive application stated in its abstract and main
   theorem is ternary, not the binary all-integers claim. This is relevant
   prior art for the set and analytic framework, not an identified proof
   of the present result.

3. **Robin Chapman**, *Representations of integers by the form
   x²+xy+y²+z²+zt+t²*, International Journal of Number Theory 4 (2008),
   709--714, DOI 10.1142/S1793042108001638.
   Author bibliography:
   https://empslocal.ex.ac.uk/people/staff/rjchapma/publications.html

   This concerns the unrestricted representation count of the quaternary
   form. It supplies useful normalization context; it does not by itself
   exclude inert primes from the component norms.

4. **Oliver Knill**, *Goldbach for Gaussian, Hurwitz, Octavian and
   Eisenstein primes*.
   https://people.math.harvard.edu/~knill/primes/papers/goldbach.pdf

   Its Eisenstein Goldbach problem concerns sums of prime elements of the
   Eisenstein ring, a different and substantially stronger condition than
   sums of rational integers in A3. It is not a duplicate of the theorem
   investigated here.

5. **Howe, Huffman, Kohlmeyer, Kruer**, *Euler's 1747 conjecture on sums of
   integers with restricted prime factors*, Hexagon 2611.00370v1,
   12 November 2026.
   https://hexagonmath.org/pdf/2611.00370v1

   This is the immediate inspiration for torus-normalizer averaging and
   the explicit-bound-plus-certificate strategy. The Eisenstein spectral
   note specifies its own algebra, order, units, level, and constants.

Searches covered exact prime-factor formulations modulo 3 and 6, primitive
Eisenstein norms, primitive Löschian integers, sums of two primitive binary
norms, and the quaternary form with coprimality conditions. No source
proving the exact binary A3 statement was found. **Absence from these
search results is not a novelty certificate.** A specialist literature
check is still needed before a priority claim.

## Independent checks of the unbalanced argument

I looked for failures of normalization, missing one-dimensional
representations, and invalid reductions, rather than assuming the
Gaussian argument transfers automatically. No fatal gap was identified
in the current derivation.

### Algebra and representation count

- The trace-pairing blocks in `spectral_core.md` have determinants -3 and
  3; the order's reduced discriminant is 3. The algebra is definite and
  is the quaternion algebra ramified at 3 and infinity.
- The two-coordinate nearest-point norm bound is 2/3, yielding a genuine
  left/right Euclidean algorithm and principal one-sided ideals.
- The twelve units are exactly the six units in one component or the
  other. This fact is responsible for Gamma being wholly inside the
  torus normalizer; the Gaussian factor cannot simply be reused.
- Counting all integral right ideals of norm N gives sigma(N), not the
  primitive sphere count. Twelve generators per ideal yield 12 sigma(N).
- For N=2 modulo 3, both component norms are 1 modulo 3 and positive;
  zero summands and ramified prime 3 cannot appear.

### The Hecke operator

- Replacing a representative g of Gamma\\G by gamma*g changes alpha*g
  into alpha*gamma*g; right multiplication by a unit permutes Q_N.
  Therefore the displayed operator is well-defined.
- Quaternion conjugation preserves Q_N and becomes inversion
  projectively, proving self-adjointness for counting measure.
- At the ramified prime 3, a norm-3 Atkin--Lehner element is not an
  accidental stabilizer: rational scalar rescaling changes norm
  valuations by even integers, so its odd valuation cannot become an
  order-unit valuation. The chosen projective-unit level matters here.
- In the all-ideal normalization, T_(ell²)=T_ell²-ell I and its degree is
  ell²+ell+1. This agrees with the stated Satake formula and includes the
  scalar ideal. Using primitive spheres instead would change the formula.
- Norm characters restrict through the determinant square class of
  PGL2(F_p), so the constant function and the displayed quadratic
  determinant character account for the noncuspidal one-dimensional
  part. The transfer to weight-two cusp forms requires (N,3p)=1, which
  is explicitly assumed in that part of the proof.
- The standard Jacquet--Langlands/Deligne input is substantial. The
  specified level and normalization make its application plausible and
  internally consistent; an independent specialist should still verify
  this automorphic identification before publication.

### Torus averaging

- For p=5 modulo 6, the quotient sizes are |Gamma|=6,
  |K|=2(p+1), |Gamma\\K|=(p+1)/3. Hence density is 2/[p(p-1)].
- Gamma is contained in K, so an entire twelve-element left-unit class
  is bad or good; b_p=12 T_N(f_p)(Gamma) is exact.
- The norm from F_(p²) onto F_p is surjective. Its two determinant square
  classes remain equally numerous after quotienting by Gamma, whose
  determinants are squares. This gives the required orthogonality to
  the determinant character.
- Right-K invariance turns point evaluation into the normalized quadratic
  form. This is the step that removes the p-dependent loss.

### Explicit arithmetic estimates

- For the divisor product, the only primes below 32 to check are
  2,5,11,17,23,29; for M(N), they are 7,13,19,31. The finite exponent
  checks are justified globally because the simple next-term ratio is
  <=1 from the reported endpoint onward. The additional corrected-tau
  ratio is <1. It need not itself be monotone for this proof to work.
- The lattice label L=x(x-1)/2+y+1 is injective on the unit-sector
  representatives x>=1 and 0<=y<x. Its difference from the norm is
  (y-(x+1)/2)^2+(x²-5)/4; x=1,2 are checked separately. Thus the bound
  of 6t for nonzero lattice points is valid without an extra constant.
- The bound for the tail of odd inert primes can be obtained using the
  wheel residues 11,17,23,29 modulo 30. Their counting discrepancy has
  total oscillation 8/5. Integrating the interval-count bound yields
  2/(15P)+8/(5P²), for P>=5. One cannot derive that exact leading
  constant merely by dropping the negative endpoint term from a global
  prime-count bound; the interval-discrepancy argument is needed.
- With P=2X and X=N^(1/4), small-prime errors are
  (2/3)X^(-1/5)+(17/6)X^(-6/5). Large-prime errors are
  (3/4)X^(-1/5)+(9/2)X^(-6/5). Their sums are exactly 17/12 and 22/3
  in the stated inequality. Including the 1/4 and 3/20 main losses gives
  the constant 2/5.
- At X>=100, X^(-1/5)<2/5 and X^(-6/5)<1/250, yielding exactly
  249/250<1. The very small numerical margin is nevertheless strict.

### Finite certificate

- The residue index k represents 6k+1, so it excludes 2 and 3 at the
  outset. Any integer 1 modulo 6 with a forbidden 5-modulo-6 prime has
  an even number of such prime factors with multiplicity. At least one
  lies below its square root. Thus sieving only the indicated small
  forbidden primes is exact.
- The first multiple of p=5 modulo 6 that is 1 modulo 6 is 5p; subsequent
  such multiples differ by 6p. The slice start and step are correct.
- The witness encoding sums to 6k+2, includes 1, keeps a<=b, and never
  indexes beyond the sieve. The balanced inequality is correctly
  translated to t>=ceil(k/3)=(k+2)//3.
- There are exactly 16,666,667 targets below 100 million, ending at
  99,999,998. The analytic region beginning at 100 million leaves no
  uncovered admissible target.

## Outcome

No concrete defect was found in the unbalanced proof or finite certificate
by this review. The two material qualifications remain the lack of
external specialist verification of the automorphic step and the lack of
a conclusive novelty search. The balanced statement remains computational
evidence plus a proposed route; it has not been promoted to a theorem.
