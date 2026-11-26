# Exact calibration of whole shared Newton residue blocks

In the [grouped shape test](MULTIVARIATE_GROUPED_SHAPE_CONGRUENCES.md) and the [quadratic denominator-prime test](MULTIVARIATE_DENOMINATOR_PRIME_SHAPE_CONGRUENCES.md), the required local unit factors are automatically determined to any finite precision when the candidate takes whole coprime blocks of the shared residual gcd. This applies when the gcd is squarefree or when the candidate saturates a shared residual width. It uses the two original factorial parents together and allows repeated residual factors in both.

## A common factor determines a lifted block

Let \(P_n,P_N\) be original adjacent factorial parents. Fix one proposed valuation \(a/e\). In the Eisenstein field \(K=\mathbb Q_p(\pi)\), \(\pi^e=p\), scale each parent by \(z=\pi^aX\) and divide its coefficients by their minimum power of \(\pi\). Call the resulting integral polynomials \(F_n,F_N\). Each has nonzero reduction over the residue field \(k\).

If there is a rational common factor, the monic group polynomial \(H\) of its roots at this valuation, after the same scaling, has integral coefficients and unit constant. It divides \(F_n,F_N\) over the valuation ring \(\mathcal O\). Write their monic residual gcd as \(X^\ell G\), where \(G(0)\ne0\). Since \(H(0)\) is a unit, its reduction divides \(G\).

Fix any monic \(G_0\mid G\) with
\[
                            \gcd(G_0,G/G_0)=1 .                 \tag{1}
\]
If \(G_0\) divides the reduction of \(H\), then \(H\) has a unique monic factor with reduction \(G_0\). More importantly, that local factor can be computed without knowing \(H\). Choose any residual Bézout combination of the two scaled reductions equal to \(X^\ell G\), and lift the two multipliers to \(\mathcal O[X]\). The mixture \(F\) of the two scaled parents reduces exactly to \(X^\ell G\). Then a unique monic factor \(H_0\) of \(F\) reducing to \(G_0\) exists over \(\mathcal O\), and every possible original common \(H\) whose reduction contains \(G_0\) must contain this \(H_0\).

For the factor lift, start with monic \(H_0\) modulo \(\pi\) equal to \(G_0\), and quotient \(Q_0=X^\ell G/G_0\). These are coprime because \(G_0(0)\ne0\) and (1). At stage \(L\) with \(F-H_0 Q\equiv0\bmod\pi^L\), solve
\[
 \overline{(F-H_0Q)/\pi^L}
      = \delta H\,Q_0+G_0\,\delta Q,\qquad
                        \deg\delta H<\deg G_0 .
\]
Reducing modulo \(G_0\) and inverting \(Q_0\) determines \(\delta H\) uniquely. The quotient then determines \(\delta Q\); its degree is at most \(\deg F-\deg G_0\) by the displayed identity. Lifting digits gives the claimed factor and uniqueness. The quotient of \(F\) may have leading coefficient of positive valuation. The argument uses the monic factor \(H_0\) and does not require the leading coefficient of \(F\) to be a unit.

To see why a common \(H\) contains this factor, apply the same lift inside the monic \(H\). Its residual quotient by \(G_0\) is coprime to \(G_0\) since \(\bar H\mid G\) and (1). That factor then divides the mixture \(F\) by divisibility of \(H\) in both parents, so uniqueness identifies it with the computed \(H_0\).

## Consequences for the denominator digit test

If the proposed local root count equals the degree of \(G\), test the unique lift of the entire \(G\). A failed division against either parent at any finite precision excludes that root count. A successful comparison to the chosen precision gives all unit-factor digits needed by the coefficientwise shape formula. If \(G\) is squarefree, enumerate its distinct irreducible factors over the residue field. Each has a unique lifted block and each can be tested against both parents. Every possible local group of a global common factor is composed of some of the surviving blocks. In the repeated case the same calculation applies to a full primary block \(g^r\) of \(G\), taken with its entire multiplicity.

In particular, for a common quadratic with two distinct integer valuations \(s<t\), suppose each selected unit root is simple in the appropriate shared residual gcd. The joint mixture recovers each scaled unit to the whole valuation-spread precision \(w=t-s\), even if one of the parent reductions alone is multiple at that root. The denominator-prime congruence then supplies at least the full local denominator power for the scaled rational shape. When those residues are unique at every denominator prime, this removes the denominator loss from the joint native interval test. The same implication holds in higher common degree for full calibrated groups. Ambiguous choices of which blocks occur, and partial multiplicities of repeated common blocks, remain separate constraints.

There is a direct lower-precision version from either one of the parents. If the proposed group fills its entire Newton edge with scaled indices \(l,u\), retain all coefficients between \(l\) and \(u\), including those strictly above the edge. Let \(L\) be the minimum valuation of a coefficient outside this window after scaling. Modulo \(\pi^L\) the parent is \(X^l\) times the window polynomial, whose constant and leading coefficients are units. A degree-\(u-l\) common group factor with unit constant must equal that window after monic normalization. Convexity of the lower polygon gives \(L\) at least the minimum adjacent slope separation, measured with \(w=e v_p\). This shortcut uses only the original marked window when its precision suffices.

## Digit ambiguity for a partial repeated block

A monic candidate \(H_0\bmod\pi\) need not be a whole coprime
block. At any branch already lifted through precision \(\pi^L\),
take the current monic candidate \(H\) and its two quotients
\(Q_j\) modulo \(\pi^L\). Its next monic correction \(\delta H\)
of degree below \(m=\deg H_0\) must satisfy
\[
         \delta H\,\bar Q_j
           \equiv \overline{(F_j-HQ_j)/\pi^L}\pmod{H_0},
                           \qquad j=n,N .                     \tag{2}
\]
This follows by reducing the two quotient-correction equations
modulo \(H_0\). If the pair is solvable, its solution space before
additional conditions is a translate of the kernel of multiplication
by the two \(\bar Q_j\)'s on \(k[X]/(H_0)\). Put
\(J_0=\gcd(H_0,G/H_0)\). This kernel consists exactly of multiples of
\(H_0/J_0\) of degree less than \(m\), and its dimension is
\(\deg J_0\). Indeed divisibility of \(\delta H\bar Q_j\) by \(H_0\)
requires \(\delta H\) to be a multiple of
\(H_0/\gcd(H_0,\bar Q_j)\); intersect the two conditions, using
\(\gcd(\bar Q_n,\bar Q_N)=X^\ell G/H_0\) and \(H_0(0)\ne0\).

There is a further rationality restriction. Before the scaling
\(z=\pi^aX\), the group polynomial for the roots of valuation \(a/e\)
is monic over \(\mathbb Q_p\), since that set of roots is invariant
under its Galois group. After scaling, its coefficient of \(X^i\)
lies in \(\pi^{r_i}\mathbb Z_p\), where \(r_i\) is the least
nonnegative residue of \(a(i-m)\bmod e\). Choose the current
representative with digits of this form. Only coefficient positions
\(i\) with \(L\equiv r_i\bmod e\) may then have a nonzero
\(\pi^L\) correction. Intersect (2) with this coordinate subspace.
If nonempty, its dimension is at most the dimension of the
intersection of that subspace with the preceding kernel. This is
a finite linear calculation at every level. It reduces repeated
block ambiguity without assuming simple residual roots.

For an exact example at \(p=2\), take eight \((0,3)\) rows, five
\((3,0)\) rows, and \((3,4),(5,1)\) in orders \(5,6\). At
valuation \(9/2\) the two raw scaled reductions are
\(1+X^4\) and \(X+X^5\). For the partial candidate \(H_0=1+X^2\)
both residual cofactors are multiples of \(H_0\). At the first
digit all four unconstrained corrections modulo \(\pi^2\) divide
both scaled arrays there. Rationality permits only the zero and
linear corrections at this digit, since \(a=9,e=2,m=2\).

## An ordinary paired calculation

Take four rows \((2,5)\) and four rows \((5,2)\), in orders \(8,9\), at \(p=3\) and valuation zero. The raw scaled residues are
\[
 \bar F_8=1+X^2+2X^3+2X^5+X^6+X^8,\qquad
 \bar F_9=X+X^2+X^7+X^8 .
\]
The shared residue gcd is the squarefree \(G=(X+1)(X^2+1)\). Its overlap with the residual cofactor of the first parent is \(X+1\), while its overlap with the other residual cofactor is \(X^2+1\). Thus neither parent individually separates the entire block. One may take multipliers \(1-X^2+X^3-X^4\) and \(1+X^2+X^3+X^4\) for their residual combination giving \(G\).

The unique linear mixed lift is \(X+7\bmod9\), with remainder \(3\) in each parent, so that unit branch is excluded. The quadratic mixed lift is \(X^2+3X+1\bmod9\) and divides both arrays at that precision. At the next precision it lifts to \(X^2+12X+19\bmod27\); the two remainders are \(9(X+1)\) and \(-9(X+1)\), which exclude the branch. Since zero is the only shared root valuation at \(3\), there is no rational common factor at this pair. The conclusion follows from these local block lifts and is not intended to exclude independent certificates at other primes.

The [exact checker](verify_multivariate_shared_residual_block_lifting.py) recomputes the original factorial ramps and Newton slopes, checks the crossed residual cofactors and Bézout combination, and lifts each proper block and the full block by testing all digit corrections at the two stated precisions. It also checks the fractional-slope partial-block example.
