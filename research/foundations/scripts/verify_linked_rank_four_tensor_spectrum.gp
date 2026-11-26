\\ Run from the repository root.
read("research/foundations/scripts/verify_rank_four_pfaffian_curve.gp");

verify_tensor_pattern() = {
  my(counts, coefficients, polynomial);
  branch_B = x^3 - x^2 - x - 1;
  check_equal(
    gcd(poldisc(branch_A) * poldisc(branch_B) * polresultant(branch_A, branch_B), 21),
    1, "distinct tame branches for the tensor spectrum"
  );
  counts = vector(4, r, character_counts(7, r, 2));
  check_equal(counts, [[2,3,1],[16,13,17],[135,102,102],[832,833,733]], "tensor pattern counts");
  coefficients = newton_coefficients(
    vector(4, r, -1 - sum(i=0, 2, counts[r][i+1] * cubic_unit^i))
  );
  polynomial = sum(j=0, 4, coefficients[j+1] * x^(4-j));
  check_equal(
    polynomial, (x^2 + 2*(cubic_unit+1)*x + 7*cubic_unit) * (x^2 - 7*cubic_unit),
    "tensor pattern polynomial"
  );
  print("p=7 verified exact relative quartic with tensor-flip spectral pattern");
};

iferr(verify_tensor_pattern(), failure, print(failure); quit(1));
