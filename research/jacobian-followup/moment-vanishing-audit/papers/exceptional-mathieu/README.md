# Haar moments and the Mathieu property of reductive homogeneous spaces

The standalone paper is available as a [PDF](../../../../../output/pdf/mathieu-reductive-homogeneous-spaces.pdf) and an
[editable single-file LaTeX export](../../../../../output/pdf/mathieu-reductive-homogeneous-spaces.tex). The canonical source is
[main.tex](main.tex). It is 16 pages, including proofs, appendices and references.
It does not require the earlier repository notes to follow the argument.

The paper proves the exceptional fixed-multiplier counterexamples, the full
reductive classification including disconnected acting groups, the separate
moment-nullcone criterion, polynomial recomposition failures, transcendence of
ordinary marked resolvents, and the coefficient-algebra statement on R[P].
Classical inputs are stated explicitly and cited. The characteristic-p
extraction rule, group reduction and invariant-form checks are proved in the paper.

The simplification discovered while writing is the radial identity in Section5.
The common construction with f4(X,Y)=X^2Y^2 gives

    P4=x^2 y^2 (1+x z w^2)^4 (1+y z^2 w)^4,   Q4=z^2 w^2,
    beta_m=(2m/(3m+1)) (2m)! (4m)!/(6m)! >0.

For every m>=1 the normalized Haar moments are h(P4^m)=0 and
h(Q4 P4^m)=beta_m. On C^2 use M(v)=[[v0,-conjugate(v1)],
[v1,conjugate(v0)]] and compose these polynomials with M. For every
compactly supported finite positive SU2-invariant measure mu,

    integral P4(M(v))^m dmu=0,
    integral Q4(M(v)) P4(M(v))^m dmu
        =beta_m integral norm(v)^(12m-4) dmu.

Thus the same construction supplies the line-and-doublet and normalizer
witnesses needed in the classification. The proof uses a beta integral
and a partial-fraction identity. No separate special radial theorem or
old normalizer circuit is assumed. The exceptional proof also avoids
boundary geometry and rational Haar transport.

To build from this directory, run `latexmk -pdf main.tex` with a standard
TeX Live installation. All text and references are in that one source file.
The optional command `python verify.py` requires SymPy and checks the
printed polynomial identities and two radial normalizations. It is not
an all-degree theorem verifier. The [QA record](qa.json) records the build,
source/PDF hashes and rendered-page inspection.

The paper is a research draft. Independent replay, formal verification and
external review are not claimed; external review remains deferred.
