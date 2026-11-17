# Bertini supplies the joint-selector extension

Reviewed 15 November 2026. [Proof](PROOF.md), [source pins](sources.json). This records a new application with written arguments, not a publication-priority certificate.

The precise imported result is [Stacks Lemma37.32.1](https://stacks.math.columbia.edu/tag/0G4C), a generic-field Bertini lemma. With algebraically independent u,v, the proof checks its four-element independence condition and obtains geometric integrality of K(b) over F(b,v+b*u). In the compiler u=S^2 and v is the inverse polynomial without its independent quadratic tilt and constant target. An active profile makes v nonconstant along the auxiliary curve. The lemma can therefore be iterated over every active block with no degree bound.

[Poonen and Slavov, The exceptional locus in the Bertini irreducibility theorem for a morphism](https://arxiv.org/html/2001.08672v3), Theorem1.5, gives modern morphism-Bertini context and quantitative exceptional-locus bounds. Those bounds and its finite-field random-slicing proof are not needed here; the application uses the precise Stacks field lemma. This is published mathematics, not a new external AI-assisted discovery.

[Stacks Lemma10.47.12](https://stacks.math.columbia.edu/tag/0G33) supplies the relative-algebraic-closedness interpretation. The [smooth-morphism section](https://stacks.math.columbia.edu/tag/01V4), flatness and composition, supports polynomial saturation for the already explicit intermediate affine map. No arbitrary imprimitive intermediate field is promoted to a polynomial model.

UARANKONE1 supplied the earlier Schur criterion; UAALLSELECT1 and UAONECHANNEL1 treated selectors independent of auxiliary sources. The Bertini argument now determines the entire joint kernel as a polynomial extension of the auxiliary constant ring. It also permits every nonzero constant Jacobian, via an affine auxiliary slice term. The earlier source-only common-kernel description remains valid.

The factorization criterion concerns the FULL auxiliary constant ring. It does not assume all such rings are polynomial, finitely generated or stably polynomial. Stable polynomiality after any constant extension is a sufficient obstruction to absolute stable atomicity. Local nilpotence with a slice gives a direct automorphism argument; the slice alone gives no local nilpotence claim. Products of two retained cubics explicitly demonstrate why polynomial kernels extend beyond output recovery.

Higher-rank return and auxiliary invariant geometry outside this criterion remain open. The all-active natural marked-fibre obstruction requires the auxiliary source to be fixed by S_N. No new complete canonical fibre, dimension optimum, absolute atom, formal proof or external review is supplied.
