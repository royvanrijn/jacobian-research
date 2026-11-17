\\ Exact local certificates for the three quartic S_4 extensions of Q_2.
\\
\\ Roe--Turturean count three Aut(S_4)-orbits of admissible epimorphisms.
\\ The three local fields below are the corresponding unlabelled arithmetic
\\ models in the Jones--Roberts/LMFDB classification.  In addition to their
\\ group and ramification data, this script computes the normalized relative
\\ degree-two Stiefel--Whitney obstruction of their trace forms.  That bit
\\ names the unique anisotropic Roe--Turturean lift.

allocatemem(512000000);

pcorder(pc) =
{
  my(result = 1);
  for (i = 1, #pc[2], result *= pc[2][i]);
  result
}

traceform(P) =
{
  my(nf = nfinit(P), n = poldegree(P));
  matrix(n, n, i, j, nfelttrace(nf, Mod(x^(i+j-2), P)))
}

hasse2(G) =
{
  my(diagonal = qfgaussred(G, 1)[2], result = 1);
  for (i = 1, #diagonal,
    for (j = i + 1, #diagonal,
      result *= hilbert(diagonal[i], diagonal[j], 2)
    )
  );
  result
}

\\ The zero section S_3 fixes one point and acts through the unique tame
\\ cubic.  Hence its degree-four trace form is <1> plus the trace form of
\\ Q_2[root]/(root^3-2); the <1> summand changes neither determinant nor
\\ Hasse invariant.
base_trace = traceform(x^3 - 2);
base_disc = matdet(base_trace);
base_hasse2 = hasse2(base_trace);

isq2(a) =
{
  my(v = valuation(a, 2), unit);
  if (v % 2, return(0));
  unit = lift(Mod(a / 2^v, 8));
  unit == 1
}

relative_sw2(P) =
{
  my(G = traceform(P), d0 = base_disc);
  if (!isq2(matdet(G) / d0),
    error("quartic and tame-base discriminants have different square classes")
  );

  \\ Degree two of w(Tr_P) / w(Tr_base):
  \\ w_2(P) + w_2(base) + w_1(base)^2.
  my(relative_hasse_sign =
    hasse2(G) * base_hasse2 * hilbert(d0, d0, 2));
  if (relative_hasse_sign == 1, 0, 1)
}

checkresolventkummer(label, p, q, r) =
{
  \\ For T^4+pT^2+qT+r, the cubic resolvent roots y_i satisfy
  \\ (root_1+root_2)^2 = y_i-p.  Over the tame S_3 splitting field these
  \\ three square classes span V_4 and have their sole relation in degree 3.
  my(R = x^3 - p*x^2 - 4*r*x + (4*p*r - q^2));
  my(S = subst(nfsplitting(R), x, y), nf = nfinit(S));
  my(primes = idealprimedec(nf, 2));
  if (#primes != 1 || primes[1].e != 3 || primes[1].f != 2,
    error(Str(label, ": resolvent closure is not the tame (e,f)=(3,2) field"))
  );

  my(roots = nfroots(nf, R));
  if (#roots != 3, error(Str(label, ": resolvent does not split in its closure")));
  my(classes = vector(3, i, roots[i] - p), pr = primes[1]);
  for (i = 1, 3,
    if (nfislocalpower(nf, pr, classes[i], 2),
      error(Str(label, ": a resolvent Kummer generator is a local square"))
    )
  );
  for (i = 1, 3,
    for (j = i + 1, 3,
      if (nfislocalpower(nf, pr, classes[i] * classes[j], 2),
        error(Str(label, ": two resolvent Kummer classes coincide"))
      )
    )
  );
  if (!nfislocalpower(nf, pr, classes[1] * classes[2] * classes[3], 2),
    error(Str(label, ": missing product-one Kummer relation"))
  );
  R
}

checkmodel(label, P, p, q, r, expected_disc, expected_ram_orders, expected_relative_sw2) =
{
  if (!polisirreducible(P), error(Str(label, ": polynomial is reducible")));

  my(pg = polgalois(P));
  if (pg[1] != 24 || pg[4] != "S4",
    error(Str(label, ": expected global S4, got ", pg))
  );
  if (poldisc(P) != expected_disc,
    error(Str(label, ": wrong polynomial discriminant"))
  );

  my(v = valuation(expected_disc, 2));
  my(unit_mod_8 = lift(Mod(expected_disc / 2^v, 8)));
  if (v % 2 != 0 || unit_mod_8 != 5,
    error(Str(label, ": discriminant root is not unramified quadratic"))
  );

  my(nf4 = nfinit(P), quartic_primes = idealprimedec(nf4, 2));
  if (#quartic_primes != 1 ||
      quartic_primes[1].e != 4 ||
      quartic_primes[1].f != 1,
    error(Str(label, ": quartic completion is not totally ramified"))
  );

  my(S = nfsplitting(P));
  if (poldegree(S) != 24,
    error(Str(label, ": splitting field degree is not 24"))
  );
  my(nfS = nfinit(S), gal = galoisinit(S));
  if (#gal.group != 24, error(Str(label, ": Galois model has wrong order")));

  my(split_primes = idealprimedec(nfS, 2));
  if (#split_primes != 1 ||
      split_primes[1].e != 12 ||
      split_primes[1].f != 2,
    error(Str(label, ": splitting completion does not have e=12, f=2"))
  );

  my(ram = idealramgroups(nfS, gal, split_primes[1]));
  my(actual_ram_orders = vector(#ram, i, pcorder(ram[i])));
  if (actual_ram_orders != expected_ram_orders,
    error(Str(label, ": expected ramification orders ",
              expected_ram_orders, ", got ", actual_ram_orders))
  );

  my(resolvent = checkresolventkummer(label, p, q, r));
  my(sw2 = relative_sw2(P));
  if (sw2 != expected_relative_sw2,
    error(Str(label, ": expected relative SW2 bit ",
              expected_relative_sw2, ", got ", sw2))
  );

  print(label, " disc=", expected_disc,
        " quartic=(e,f)=(4,1) closure=(24;12,2) ram=",
        actual_ram_orders, " resolvent=", resolvent,
        " relative-sw2=", sw2)
}

checkmodel("s4-local-2.1.4.4a1.1", x^4 - 2*x + 2, 0, -2, 2, 1616, [24, 12, 4], 0);
checkmodel("s4-local-2.1.4.8a1.1", x^4 - 4*x + 2, 0, -4, 2, -4864, [24, 12, 4, 4, 4, 4, 4], 0);
checkmodel("s4-local-2.1.4.8a1.2", x^4 + 4*x^2 - 4*x + 2, 4, -4, 2, 9472, [24, 12, 4, 4, 4, 4, 4], 1);

print("PASS: the three dyadic quartic S_4 models and Kummer classes are exact");
print("PASS: the unique relative-SW2-one model is 2.1.4.8a1.2");
