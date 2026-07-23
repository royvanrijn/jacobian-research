# The Gaussian Moments Conjecture in Two Variables

Standalone source for the paper by Roy van Rijn.

Build from this directory with:

```sh
latexmk -pdf main.tex
```

The proof has four dependencies:

1. the circular-coordinate identity
   \(\mathbb E(F)=\mathcal L(\operatorname{CT}_T F)\);
2. an elementary supporting-line argument for the lowest radial orders;
3. the Duistermaat--van der Kallen constant-term theorem;
4. Frobenius and factorial isolation after reduction modulo a prime.

The computational support classifications in the surrounding repository
motivated the result but are not needed by the proof.
