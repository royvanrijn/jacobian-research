# Quantitative Hasse-Principle Failures in the Fibers of a Fixed Keller Map

**Unpublished · draft.**

Classifies every Hasse-failing rational parameter on one line of a fixed degree-five Jacobian-one map and proves the exact target-height asymptotic `N(B) ~ C B²/log B`.

[PDF](paper.pdf) · [LaTeX](main.tex) · [Metadata](metadata.json)

The complete criterion is in [the local calculation](sections/03-local-global.tex), and the height normalization, coprimality argument, coupled 3-adic condition and Euler product are in [the rational count](sections/04-rational-counting.tex). The constant is

`C = 11/(288√3 π²) ∏_{p≡2 mod3}(1+p⁻³)/(1−p⁻³)`.

[Verification scope](VERIFICATION.md) distinguishes the written rational proof, internal finite checks and the earlier integer-family Lean statements. The theorem counts targets on this line; ambient targets and distinct fibre algebras remain separate problems.

Run the current finite audit from the repository root:

```sh
python3 unpublished/fixed-map-hasse-failures/verify_rational_line.py --maximum-height 10000 --output unpublished/fixed-map-hasse-failures/rational_line_verification.json
```

It gives `N(1000)=617`, `N(3000)=4377`, `N(10000)=41415`. These are finite audits, not the asymptotic proof.

Next: Obtain specialist review of the complete rational count and complete the priority review in [LITERATURE_AUDIT.md](LITERATURE_AUDIT.md) before preparing a submission bundle.
