\\ Coarse Q_2 decompositions of the selected quartic and quintic fibers.
\\
\\ For each irreducible Q-factor, idealprimedec over the maximal order gives
\\ the primes above 2.  We record [local degree=e*f, e, f].  The output is a
\\ multiset because a Q-irreducible factor may split over Q_2.

localdata(P) =
{
  my(F = factor(P), out = List());
  for (i = 1, matsize(F)[1],
    my(q = F[i, 1] / pollead(F[i, 1]), m = F[i, 2], d = poldegree(q));
    if (d == 1,
      for (j = 1, m, listput(out, [1, 1, 1])),
      my(nf = nfinit(polredbest(q)), D = idealprimedec(nf, 2));
      for (k = 1, #D, listput(out, [D[k].e * D[k].f, D[k].e, D[k].f]))
    )
  );
  vecsort(Vec(out))
}

checkrow(label, P, expected) =
{
  my(actual = localdata(P));
  if (actual != expected,
    error(Str(label, ": expected ", expected, ", got ", actual))
  );
  print(label, " ", actual)
}

\\ Four canonical quartic arithmetic witnesses already used by the repository.
checkrow("quartic-witness-card", x^4 - 3*x^2 - 1, [[4, 2, 2]]);
checkrow("quartic-common-u1", x^4 + x^3 - 2*x^2 + x + 1, [[2, 1, 2], [2, 2, 1]]);
checkrow("quartic-small", 2*x^4 - x^3 - x^2 + x + 1, [[1, 1, 1], [3, 3, 1]]);
checkrow("quartic-optimized", 9*x^4 - 19*x^3 + 10*x^2 - 8*x - 4, [[1, 1, 1], [3, 3, 1]]);

\\ All ten rows of the generated fixed-quintic arithmetic-zoo ledger.
checkrow("quintic-split", x^5 - 5*x^3 + 4*x, [[1, 1, 1], [1, 1, 1], [1, 1, 1], [1, 1, 1], [1, 1, 1]]);
checkrow("quintic-S5-real-1", x^5 - 5*x^3 + 2*x^2 + 4*x + 2, [[2, 1, 2], [3, 3, 1]]);
checkrow("quintic-S5-real-3", x^5 - 5*x^3 - 2*x^2 - 4*x - 2, [[2, 1, 2], [3, 3, 1]]);
checkrow("quintic-S5-real-5", x^5 - 5*x^3 + 4*x + 1, [[5, 1, 5]]);
checkrow("quintic-A5", x^5 - 5*x^3 + 4*x + 4/5, [[2, 1, 2], [3, 3, 1]]);
checkrow("quintic-C5", x^5 - 5*x^3 + 4*x + 7/5, [[5, 1, 5]]);
checkrow("quintic-D5", x^5 - 5*x^3 + 42/25*x^2 + 32/125*x - 128/3125, [[1, 1, 1], [2, 2, 1], [2, 2, 1]]);
checkrow("quintic-F20", x^5 - 5*x^3 - 3/2*x^2 + 1/2*x - 1/40, [[5, 5, 1]]);
checkrow("quintic-product", x^5 - 5*x^3 + 3*x^2 + 4*x + 9, [[2, 1, 2], [3, 3, 1]]);
checkrow("quintic-Hasse", x^5 - 5*x^3 + 288*x^2 + 500*x + 376, [[1, 1, 1], [1, 1, 1], [3, 1, 3]]);

\\ The separate connected rank-five multiplicity witness.
checkrow("quintic-witness-card", x^5 + x^3 + 1, [[5, 1, 5]]);

print("PASS: exact PARI ideal decomposition reproduces every Q_2 row");
