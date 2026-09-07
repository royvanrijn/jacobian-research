\\ Independent PARI/GP verifier for proof-carrying arithmetic compilation.
\\ Run from the repository root. Set ARITHMETIC_CERTIFICATE to override input.

jqlines(filter) =
{
  my(command = Str("jq -r '", filter, "' '", certificatePath, "'"));
  my(lines = externstr(command));
  if(#lines == 0, error(Str("empty jq result for ", filter)));
  return(lines);
}

jqlinesOptional(filter) =
{
  my(command = Str("jq -r '", filter, "' '", certificatePath, "'"));
  return(externstr(command));
}

jqone(filter) =
{
  my(lines = jqlines(filter));
  if(#lines != 1, error(Str("expected one jq result for ", filter)));
  return(lines[1]);
}

splitFields(line) = strsplit(line, "\t");

parseCSV(text) =
{
  if(text == "", return([]));
  my(parts = strsplit(text, ","));
  return(vector(#parts, index, eval(parts[index])));
}

assertEqual(actual, expected, label) =
{
  if(actual != expected,
    error(Str(label, ": expected ", expected, ", got ", actual)));
}

assertTrue(condition, label) =
{
  if(!condition, error(label));
}

canonicalExpandedMap(inputMap) =
{
  my(payload = "[", counts = vector(3), degrees = vector(3));
  for(coordinate = 1, 3,
    if(coordinate > 1, payload = Str(payload, ","));
    payload = Str(payload, "[");
    my(firstTerm = 1, count = 0, totalDegree = 0);
    for(exponentX = 0, poldegree(inputMap[coordinate], x),
      for(exponentY = 0, poldegree(inputMap[coordinate], y),
        for(exponentZ = 0, poldegree(inputMap[coordinate], z),
          coefficient = polcoeff(polcoeff(polcoeff(
            inputMap[coordinate], exponentX, x), exponentY, y), exponentZ, z);
          if(coefficient != 0,
            if(!firstTerm, payload = Str(payload, ","));
            payload = Str(payload, "[[", exponentX, ",", exponentY, ",",
              exponentZ, "],\"", coefficient, "\"]");
            firstTerm = 0;
            count++;
            totalDegree = max(totalDegree,
              exponentX + exponentY + exponentZ)))));
    payload = Str(payload, "]");
    counts[coordinate] = count;
    degrees[coordinate] = totalDegree);
  payload = Str(payload, "]");
  return([payload, counts, degrees]);
}

verifyCertificate() =
{
iferr(externstr("command -v jq"), E,
  error("jq is required only to parse the portable JSON certificate"));

certificatePath = getenv("ARITHMETIC_CERTIFICATE");
if(certificatePath == 0,
  certificatePath = "artifacts/generated-results/arithmetic_keller_quintic.json");
if(#strsplit(certificatePath, "'") != 1,
  error("certificate path must not contain a single quote"));

assertEqual(jqone(".schema"),
  "proof-carrying-arithmetic-compilation/v1", "schema");

globalCoefficients = parseCSV(jqone(
  ".global_polynomial.monic_coefficients_ascending | join(\",\")"));
degree = eval(jqone(".degree"));
assertEqual(#globalCoefficients, degree + 1, "global coefficient count");
assertEqual(globalCoefficients[degree + 1], 1, "monicity");
T = 'T;
P = Polrev(globalCoefficients, T);

baseDenominator = eval(jqone(".coefficient_crt.base_denominator"));
crtModulus = eval(jqone(".coefficient_crt.crt_modulus"));
multiplier = eval(jqone(".coefficient_crt.multiplier"));
commonDenominator = eval(jqone(".coefficient_crt.common_denominator"));
residues = parseCSV(jqone(
  ".coefficient_crt.coefficient_residues_ascending | join(\",\")"));
assertEqual(commonDenominator,
  baseDenominator * (1 + multiplier * crtModulus), "CRT denominator");
assertEqual(#residues, degree, "CRT residue count");
for(index = 1, degree,
  assertTrue(denominator(commonDenominator * globalCoefficients[index]) == 1,
    "CRT coefficient was not integral after scaling");
  assertEqual(lift(Mod(commonDenominator * globalCoefficients[index], crtModulus)),
    lift(Mod(residues[index], crtModulus)), "global CRT residue"));

intervalLines = jqlines(
  ".coefficient_crt.real_coefficient_intervals[] | @tsv");
assertEqual(#intervalLines, degree, "coefficient interval count");
for(index = 1, degree,
  fields = splitFields(intervalLines[index]);
  assertTrue(eval(fields[1]) < globalCoefficients[index]
    && globalCoefficients[index] < eval(fields[2]),
    "global coefficient outside real box"));

modulusLines = jqlinesOptional(
  ".coefficient_crt.local_moduli[] | [.prime,.modulus] | @tsv");
localPrimes = vector(#modulusLines);
localModuli = vector(#modulusLines);
for(index = 1, #modulusLines,
  fields = splitFields(modulusLines[index]);
  localPrimes[index] = eval(fields[1]);
  localModuli[index] = eval(fields[2]));

localLines = jqlinesOptional(Str(
  ".local_models_and_precision_claims[] | ",
  "[.prime,.discriminant,.discriminant_valuation,.universal_precision,",
  ".claimed_precision,.coefficient_modulus,",
  "(.model_coefficients_ascending|join(\",\"))] | @tsv"));
localModels = vector(#localLines);
modelPrimes = vector(#localLines);
for(index = 1, #localLines,
  fields = splitFields(localLines[index]);
  localPrime = eval(fields[1]);
  modelPrimes[index] = localPrime;
  assertTrue(isprime(localPrime), "local prime is not prime");
  modelCoefficients = parseCSV(fields[7]);
  model = Polrev(modelCoefficients, T);
  localModels[index] = model;
  claimedDiscriminant = eval(fields[2]);
  assertEqual(poldisc(model), claimedDiscriminant, "local discriminant");
  discriminantValuation = valuation(claimedDiscriminant, localPrime);
  assertEqual(discriminantValuation, eval(fields[3]),
    "discriminant valuation");
  universalPrecision = 2 * discriminantValuation + 1;
  claimedPrecision = eval(fields[5]);
  assertEqual(universalPrecision, eval(fields[4]), "universal precision");
  assertTrue(claimedPrecision >= universalPrecision,
    "claimed precision below universal radius");
  modulusIndex = select(j -> localPrimes[j] == localPrime, [1..#localPrimes]);
  assertEqual(#modulusIndex, 1, "local modulus lookup");
  modulusIndex = modulusIndex[1];
  expectedModulus = localPrime^(
    claimedPrecision + valuation(baseDenominator, localPrime));
  assertEqual(localModuli[modulusIndex], expectedModulus, "local modulus");
  assertEqual(eval(fields[6]), expectedModulus, "recorded coefficient modulus");
  for(coefficientIndex = 1, degree,
    difference = globalCoefficients[coefficientIndex]
      - modelCoefficients[coefficientIndex];
    assertTrue(difference == 0
      || valuation(difference, localPrime) >= claimedPrecision,
      "global polynomial outside local stability ball");
    scaledCenter = baseDenominator * modelCoefficients[coefficientIndex];
    assertTrue(denominator(scaledCenter) == 1,
      "base denominator did not clear a local center");
    assertEqual(lift(Mod(scaledCenter, expectedModulus)),
      lift(Mod(residues[coefficientIndex], expectedModulus)),
      "local CRT residue")));
assertEqual(prod(index = 1, #localModuli, localModuli[index]),
  crtModulus, "CRT modulus product");

primitiveDenominator = eval(jqone(".global_polynomial.common_denominator"));
primitiveCoefficients = parseCSV(jqone(
  ".global_polynomial.primitive_integer_coefficients_ascending | join(\",\")"));
assertEqual(primitiveCoefficients,
  vector(degree + 1, index,
    primitiveDenominator * globalCoefficients[index]),
  "primitive polynomial");
assertEqual(gcd(primitiveCoefficients), 1, "primitive content");
assertEqual(poldegree(gcd(P, deriv(P))), 0, "squarefreeness");

realRootCount = eval(jqone(".real_roots.count"));
assertEqual(polsturm(P), realRootCount, "total real-root count");
rootLines = jqlinesOptional(Str(
  ".real_roots.isolating_intervals[] | ",
  "[.left,.right,.multiplicity] | @tsv"));
isolatedCount = 0;
previousRight = 0;
for(index = 1, #rootLines,
  fields = splitFields(rootLines[index]);
  left = eval(fields[1]); right = eval(fields[2]);
  multiplicity = eval(fields[3]);
  assertTrue(left < right, "invalid real isolating interval");
  if(index > 1, assertTrue(previousRight < left,
    "overlapping real isolating intervals"));
  assertTrue(subst(P, T, left) != 0 && subst(P, T, right) != 0,
    "isolating endpoint is a root");
  assertEqual(polsturm(P, left, right), multiplicity,
    "real isolating interval");
  isolatedCount += multiplicity;
  previousRight = right);
assertEqual(isolatedCount, realRootCount, "isolated real-root total");

irreducibilityPrime = eval(jqone(".irreducibility_witness.prime"));
assertTrue(isprime(irreducibilityPrime), "irreducibility witness is not prime");
assertTrue(denominator(content(P)) % irreducibilityPrime != 0,
  "irreducibility prime divides a denominator");
assertTrue(polisirreducible(Mod(1, irreducibilityPrime) * P),
  "reduction is reducible at irreducibility witness prime");

translationValue = eval(jqone(".selected_translation"));
S = 'S;
translated = subst(P, T, S + translationValue);
seed = translated - subst(P, T, translationValue);
linearCoefficient = polcoeff(seed, 1, S);
assertTrue(linearCoefficient != 0, "translation kills the linear coefficient");
assertTrue(polcoeff(seed, 3, S) != 0,
  "translation kills the cubic coefficient");
target = parseCSV(jqone(".target | join(\",\")"));
assertEqual(target,
  [1, 0, -2 * subst(P, T, translationValue)
    / subst(deriv(P), T, translationValue)], "target");
seedCoefficients = parseCSV(jqone(
  ".inverse_polynomial_identity.seed_coefficients_ascending | join(\",\")"));
inverseCoefficients = parseCSV(jqone(
  ".inverse_polynomial_identity.inverse_coefficients_ascending | join(\",\")"));
translatedCoefficients = parseCSV(jqone(Str(
  ".inverse_polynomial_identity.translated_input_coefficients_ascending | ",
  "join(\",\")")));
assertEqual(Polrev(seedCoefficients, S), seed, "seed coefficients");
assertEqual(eval(jqone(".inverse_polynomial_identity.linear_coefficient")),
  linearCoefficient, "linear coefficient");
inversePolynomial = seed - linearCoefficient * target[3] / 2;
assertEqual(inversePolynomial, translated, "inverse-polynomial identity");
assertEqual(Polrev(inverseCoefficients, S), inversePolynomial,
  "inverse coefficients");
assertEqual(Polrev(translatedCoefficients, S), translated,
  "translated input coefficients");

stablePresent = eval(jqone(Str(
  "if .keller_map.stable_multiplicity == null then 0 else 1 end")));
stableParameter = if(stablePresent,
  eval(jqone(".keller_map.stable_multiplicity.family_parameter")), -1);
if(stablePresent,
  assertTrue(denominator(stableParameter) == 1 && stableParameter >= 0,
    "stable parameter is not a nonnegative integer");
  stableGaugeExponent = eval(jqone(
    ".keller_map.stable_multiplicity.gauge_exponent"));
  stableInvariant = jqone(
    ".keller_map.stable_multiplicity.separation_invariant");
  stableValue = eval(jqone(
    ".keller_map.stable_multiplicity.separation_value"));
  if(degree == 3,
    assertEqual(stableGaugeExponent, stableParameter + 4,
      "cubic gauge exponent");
    assertEqual(stableInvariant, "geometric_boundary_target_components",
      "cubic separation invariant");
    assertEqual(stableValue, stableParameter + 4,
      "cubic boundary-component count");
    assertEqual(eval(jqone(Str(
      "(.keller_map.stable_multiplicity.fitting_support // []) | length"))),
      0, "cubic fitting support");
    assertEqual(jqone(Str(
      ".keller_map.stable_multiplicity.boundary_prime_count | ",
      "if . == null then \"none\" else tostring end")),
      "none", "cubic boundary-prime count");
    assertEqual(jqone(Str(
      ".keller_map.stable_multiplicity.boundary_ramification_index | ",
      "if . == null then \"none\" else tostring end")),
      "none", "cubic boundary ramification"),
    assertEqual(stableGaugeExponent, stableParameter,
      "power-shift gauge parameter");
    assertEqual(stableInvariant, "normalized_fitting_newton_area",
      "power-shift separation invariant");
    assertEqual(stableValue,
      2 * degree - 3 + (degree - 2) * stableParameter,
      "Fitting Newton area");
    stableGcd = gcd(degree - 3, stableParameter + 2);
    assertEqual(eval(jqone(
      ".keller_map.stable_multiplicity.boundary_prime_count")),
      stableGcd, "boundary-prime count");
    assertEqual(eval(jqone(
      ".keller_map.stable_multiplicity.boundary_ramification_index")),
      (degree - 3) / stableGcd, "boundary ramification index");
    for(exponent = 4, degree,
      assertTrue(polcoeff(seed, exponent, S) != 0,
        "stable power shift has a zero higher coefficient"));
    supportLines = jqlines(
      ".keller_map.stable_multiplicity.fitting_support[] | @tsv");
    assertEqual(#supportLines, degree - 1, "Fitting support size");
    supportFields = splitFields(supportLines[1]);
    assertEqual([eval(supportFields[1]), eval(supportFields[2])],
      [0, 0], "Fitting support origin");
    supportFields = splitFields(supportLines[2]);
    assertEqual([eval(supportFields[1]), eval(supportFields[2])],
      [1, 2], "Fitting support cubic point");
    for(supportIndex = 3, #supportLines,
      supportFields = splitFields(supportLines[supportIndex]);
      supportDegree = supportIndex + 1;
      assertEqual([eval(supportFields[1]), eval(supportFields[2])],
        [supportDegree + stableParameter, supportDegree - 1],
        "Fitting support higher point"))));

actionLines = jqlinesOptional(Str(
  ".local_action_certificates[] | ",
  "[.kind,.prime,((.cycle_degrees//[])|join(\",\")),",
  "([.factors[]|(.coefficients_ascending|join(\",\"))]|join(\";\")),",
  "([.factors[].witness]|join(\",\")),",
  "([.factors[]|(.ramification_index//0)]|join(\",\")),",
  "([.factors[]|(.residue_degree//0)]|join(\",\"))] | @tsv"));
for(actionIndex = 1, #actionLines,
  fields = splitFields(actionLines[actionIndex]);
  kind = fields[1]; actionPrime = eval(fields[2]);
  factorStrings = strsplit(fields[4], ";");
  factors = vector(#factorStrings,
    index, Polrev(parseCSV(factorStrings[index]), T));
  witnesses = strsplit(fields[5], ",");
  if(kind == "local_factorization",
    modelIndex = select(j -> modelPrimes[j] == actionPrime, [1..#modelPrimes]);
    assertEqual(#modelIndex, 1, "local action model lookup");
    assertEqual(prod(index = 1, #factors, factors[index]),
      localModels[modelIndex[1]], "local action factorization");
    ramificationIndices = parseCSV(fields[6]);
    residueDegrees = parseCSV(fields[7]);
    for(leftIndex = 1, #factors,
      for(rightIndex = leftIndex + 1, #factors,
        assertEqual(valuation(polresultant(
          factors[leftIndex], factors[rightIndex]), actionPrime), 0,
          "local factors do not have unit resultant"));
      assertEqual(ramificationIndices[leftIndex]
        * residueDegrees[leftIndex], poldegree(factors[leftIndex]),
        "local component degree");
      if(witnesses[leftIndex] == "eisenstein",
        assertEqual(pollead(factors[leftIndex]), 1, "Eisenstein monicity");
        for(exponent = 0, poldegree(factors[leftIndex]) - 1,
          assertTrue(
            valuation(polcoeff(factors[leftIndex], exponent), actionPrime) >= 1,
            "Eisenstein divisibility"));
        assertEqual(
          valuation(polcoeff(factors[leftIndex], 0), actionPrime), 1,
          "Eisenstein constant valuation"),
        if(witnesses[leftIndex] == "unramified_irreducible_reduction",
          assertTrue(
            polisirreducible(Mod(1, actionPrime) * factors[leftIndex]),
            "unramified factor reduction is reducible");
          assertEqual(
            valuation(poldisc(factors[leftIndex]), actionPrime), 0,
            "unramified factor discriminant"),
          error("unknown local factor witness")))),
    if(kind == "frobenius_factorization",
      assertEqual(prod(index = 1, #factors,
        Mod(1, actionPrime) * factors[index]), Mod(1, actionPrime) * P,
        "Frobenius factorization");
      cycleDegrees = parseCSV(fields[3]);
      assertEqual(cycleDegrees,
        vector(#factors, index, poldegree(factors[index])),
        "Frobenius cycle degrees");
      for(leftIndex = 1, #factors,
        assertTrue(
          polisirreducible(Mod(1, actionPrime) * factors[leftIndex]),
          "Frobenius factor is reducible");
        for(rightIndex = leftIndex + 1, #factors,
          assertEqual(poldegree(gcd(
            Mod(1, actionPrime) * factors[leftIndex],
            Mod(1, actionPrime) * factors[rightIndex])), 0,
            "Frobenius factors are not coprime"))),
      error("unknown local action kind"))));

x = 'x; y = 'y; z = 'z;
t = 1 + x * y;
g1 = polcoeff(seed, 1, S);
g2 = polcoeff(seed, 2, S);
g3 = polcoeff(seed, 3, S);
qpoly = t^2 * z + (g1 / g3) * y^2 * (1 + 3 * t);
mapFirst = t * qpoly;
mapSecond = y + 3 * (g3 / g1) * x * qpoly
  + 2 * (g2 / g1) * t * qpoly;
mapThird = x * (5 - 3 * t) - (g3 / g1) * x^3 * z;
powerShift = if(stableParameter >= 0, stableParameter, 0);
if(stableParameter >= 0 && degree == 3,
  cubicExponent = stableParameter + 4;
  mapSecond += 3 * (g3 / g1) * (
    t^(cubicExponent - 1) * x * qpoly^cubicExponent
      - t^2 * x * qpoly^3);
    mapThird -= (g3 / g1) * (
    t^(cubicExponent - 3) * x^3 * qpoly^cubicExponent
      - x^3 * qpoly^3),
  for(exponent = 4, degree,
    mapSecond += exponent * polcoeff(seed, exponent, S) / g1
      * t^(powerShift + 2) * x^(exponent - 2)
      * qpoly^(exponent + powerShift);
    mapThird -= (exponent - 2) * polcoeff(seed, exponent, S) / g1
      * t^powerShift * x^exponent * qpoly^(exponent + powerShift)));
mapping = [mapFirst, -mapSecond / 2, mapThird];
mapVariables = [x, y, z];
jacobian = matrix(3, 3, row, column,
  deriv(mapping[row], mapVariables[column]));
assertEqual(matdet(jacobian), 1, "Jacobian determinant");

mapSerialization = canonicalExpandedMap(mapping);
assertEqual(mapSerialization[3],
  parseCSV(jqone(".keller_map.coordinate_degrees | join(\",\")")),
  "coordinate degrees");
assertEqual(mapSerialization[2],
  parseCSV(jqone(".keller_map.expanded_map_hash.term_counts | join(\",\")")),
  "expanded-map term counts");
assertEqual(jqone(".keller_map.expanded_map_hash.serialization"),
  "expanded-sparse-rational-json/v1", "map serialization");
hashOutput = externstr(Str("printf '%s' '", mapSerialization[1],
  "' | shasum -a 256"))[1];
computedHash = Strchr(Vecsmall(hashOutput)[1..64]);
assertEqual(computedHash, jqone(".keller_map.expanded_map_hash.digest"),
  "expanded-map SHA-256");
assertEqual(eval(jqone(".keller_map.geometric_degree")), degree,
  "geometric degree");

print("PASS: ", certificatePath);
print("PASS: PARI local stability, CRT, real isolation, and irreducibility");
print("PASS: PARI local actions, inverse identity, Jacobian, and map SHA-256");
}

verifyCertificate();
