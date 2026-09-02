# Heegner / CM-vector enumeration for Elkies X(6,79)

This starts from one representative of the unique-looking signature-(2,1), determinant -948 transcendental genus.

For each primitive negative vector v it computes:

- norm `v^2`;
- divisibility `gcd(<v,T>)`;
- exact positive binary complement `v^\perp`;
- the determinant identity
  `det(v^\perp) = |det(T)| |v^2| / div(v)^2`;
- the primitive positive binary quadratic form `[a,b,c]`;
- its negative quadratic-order discriminant `b^2-4ac`;
- fundamental field discriminant and conductor when recognized.

The grouping key `(norm, divisibility, reduced primitive binary form, content, order discriminant)` is an O(T)-orbit surrogate, not an orbit proof. That is deliberate: we want to identify a tiny shortlist before implementing exact indefinite-lattice automorphism/orbit certification.

Run first:

    sage elkies-k3/scripts/enumerate_heegner_orbits.sage \
      --T 0 --bound 40 --norm-max 2000 --top 100

Then summarize:

    python3 elkies-k3/scripts/summarize_heegner_orbits.py

If stable and quick, deepen only the vector bound:

    sage elkies-k3/scripts/enumerate_heegner_orbits.sage \
      --T 0 --bound 100 --norm-max 10000 --top 300 \
      --out artifacts/local/elkies-k3/heegner-orbits-deep.txt
