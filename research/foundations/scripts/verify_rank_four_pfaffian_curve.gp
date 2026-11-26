\\ Run from the repository root with PARI/GP 2.15.4 or compatible:
\\ gp -fq research/foundations/scripts/verify_rank_four_pfaffian_curve.gp

check_equal(actual, expected, label) = {
  if(actual != expected, error(label, ": expected ", expected, ", found ", actual));
};

branch_A = x^3 + x + 1;
branch_B = x^3 - 2*x^2 - 2*x + 1;
branch_invariant = poldisc(branch_A) * poldisc(branch_B) * polresultant(branch_A, branch_B);
cubic_unit = Mod(z, z^2 + z + 1);

branch_value(argument) = {
  subst(branch_A, x, argument) * subst(branch_B, x, argument)^2;
};

field_element(index, p, r, generator) = {
  my(argument = 0 * generator);
  for(j = 0, r - 1,
    argument += (index % p) * generator^j;
    index = index \ p;
  );
  argument;
};

character_counts(p, r, cube_root) = {
  my(generator = ffgen(p^r, 'u), counts = [0, 0, 0], exponent = (p^r - 1)/3);
  for(index = 0, p^r - 1,
    my(value = branch_value(field_element(index, p, r, generator)));
    if(value != 0,
      my(character = value^exponent);
      if(character == 1,
        counts[1]++,
        if(character == cube_root,
          counts[2]++,
          check_equal(character, (cube_root^2 % p) * (1 + 0 * generator), "cubic character");
          counts[3]++;
        );
      );
    );
  );
  counts;
};

newton_coefficients(traces) = {
  my(coefficients = vector(#traces + 1));
  coefficients[1] = 1;
  for(k = 1, #traces,
    coefficients[k + 1] = -sum(j = 1, k, traces[j] * coefficients[k - j + 1])/k;
  );
  coefficients;
};

squarefree_factor_degrees(polynomial, p) = {
  my(factors = factormod(polynomial, p), size = matsize(factors)[1]);
  check_equal(Vec(factors[,2]), vector(size, i, 1), "squarefree factorization");
  vecsort(vector(size, i, poldegree(factors[i,1])));
};

check_relative_factorization(polynomial, p, specialization, expected) = {
  check_equal((specialization^2 + specialization + 1) % p, 0, "coefficient field specialization");
  check_equal(
    squarefree_factor_degrees(subst(lift(polynomial), z, Mod(specialization, p)), p),
    expected, "relative factorization degrees"
  );
};

verify_split_reduction(p, cube_root, expected_counts, expected_relative, certificates) = {
  check_equal((cube_root^2 + cube_root + 1) % p, 0, "primitive cube root");
  check_equal(gcd(branch_invariant, 3*p), 1, "distinct tame branches");
  my(counts = vector(4, r, character_counts(p, r, cube_root)));
  check_equal(counts, expected_counts, "finite-field character counts");
  my(traces = vector(4, k, -1 - sum(j = 0, 2, cubic_unit^j * counts[k][j + 1])));
  my(coefficients = newton_coefficients(traces));
  my(relative = sum(k = 0, 4, coefficients[k + 1] * x^(4 - k)));
  check_equal(relative, subst(expected_relative, z, cubic_unit), "relative Frobenius polynomial");
  my(conjugate = subst(lift(relative), z, cubic_unit^2));
  if(relative == conjugate, error("the relative polynomials must be distinct"));
  check_equal(conjugate, simplify(x^4 * subst(relative, x, p/x) / polcoef(relative, 0)),
    "polarization reciprocity");
  my(absolute = simplify(lift(relative * conjugate)), middle = polcoef(absolute, 4));
  check_equal(gcd(middle, p), 1, "ordinary middle coefficient");
  my(determinant_unit = simplify(lift(polcoef(relative, 0))/p^2));
  check_equal(denominator(determinant_unit), 1, "integral scaled relative determinant");
  for(root = 1, p - 1,
    if((root^2 + root + 1) % p == 0,
      check_equal(valuation(simplify(subst(determinant_unit, z, root)), p), 0,
        "unit scaled relative determinant");
    );
  );
  check_relative_factorization(relative, certificates[1][1], certificates[1][2], [4]);
  check_relative_factorization(relative, certificates[2][1], certificates[2][2], [1, 3]);
  print("p=", p, " verified point counts, ordinary reduction, and S4 witnesses");
  absolute;
};

count_points(p, r) = {
  if(p^r % 3 == 2, return(p^r + 1));
  my(generator = ffgen(p^r, 'u), count = 3, exponent = (p^r - 1)/3);
  for(index = 0, p^r - 1,
    my(value = branch_value(field_element(index, p, r, generator)));
    if(value == 0,
      count++,
      if(value^exponent == 1, count += 3);
    );
  );
  count;
};

verify_all() = {
  check_equal([poldisc(branch_A), poldisc(branch_B), polresultant(branch_A, branch_B)],
    [-31, 125, 31], "branch discriminants and resultant");

  absolute_seven = verify_split_reduction(
    7, 2, [[3,1,2], [17,14,15], [129,93,117], [809,771,818]],
    x^4 + (-z+2)*x^3 + (-3*z+3)*x^2 + (-14*z+7)*x - 49*z,
    [[79,23], [73,8]]
  );

  absolute_thirteen = verify_split_reduction(
    13, 3, [[6,3,2], [62,46,55], [758,705,732], [9518,9631,9406]],
    x^4 + (z+5)*x^3 + 16*x^2 + (-13*z+52)*x + 169,
    [[73,8], [79,23]]
  );

  check_equal(squarefree_factor_degrees(absolute_seven, 2), [4, 4], "splitting at 2 for p=7");
  check_equal(squarefree_factor_degrees(absolute_thirteen, 2), [2, 6], "splitting at 2 for p=13");
  print("verified distinct Frobenius fields by splitting at 2");

  check_equal(gcd(branch_invariant, 3*11), 1, "distinct tame branches at 11");
  inert_counts = vector(4, k, count_points(11, k));
  check_equal(inert_counts, [12,132,1332,14172], "point counts at inert prime 11");
  inert_coefficients = newton_coefficients(vector(4, k, 11^k + 1 - inert_counts[k]));
  check_equal(inert_coefficients[5], -105, "middle coefficient at inert prime 11");
  check_equal(gcd(inert_coefficients[5], 11), 1, "ordinary middle coefficient at 11");
  print("p=11 verified ordinary reduction with nontrivial coefficient action");
};

iferr(verify_all(), failure, print(failure); quit(1));
