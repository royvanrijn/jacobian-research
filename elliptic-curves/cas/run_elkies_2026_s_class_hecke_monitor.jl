#!/usr/bin/env julia

# Monitor Hecke's exact principal-relation collector modulo the declared
# Selmer-prime ideal classes.  This is a targeted S-class experiment: reaching
# quotient dimension zero proves only that the materialized factor base is
# killed modulo 2.  A separate factor-base generation proof and unit/local
# computations are still required for a complete 2-Selmer result.

using Hecke
using Random
using TOML

const PROTOCOL = "ELKIESR17HECKESCLASS"
const DEFAULT_POLYNOMIAL_ASCENDING = (
  "4105678984643853583390832544029019669185034999158," *
  "-774250153578278482962797863407542,-1,1"
)
const DEFAULT_S_PRIMES = (
  "2,3,5,7,13,23,29,89,43207,226549," *
  "22823593909227592035983291,44013936637595415741483513793"
)

function parse_args(args)
  values = Dict{String, String}(
    "--polynomial-ascending" => DEFAULT_POLYNOMIAL_ASCENDING,
    "--s-primes" => DEFAULT_S_PRIMES,
    "--factor-base-bound" => "240",
    "--timeout-seconds" => "300",
    "--method" => "2",
    "--valuation-radius" => "0",
    "--precision" => "100",
    "--enumeration-limit" => "10",
    "--random-seed" => "20260902",
    "--include-s-in-factor-base" => "false",
    "--checkpoint" => (
      "artifacts/local/elliptic-curves/elkies-r17-hecke-s-class-v1/" *
      "control-r21-bound240.toml"
    ),
  )
  index = 1
  while index <= length(args)
    key = args[index]
    haskey(values, key) || error("unknown argument: $key")
    index == length(args) && error("missing value after $key")
    values[key] = args[index + 1]
    index += 2
  end
  return values
end

function emit(message)
  println(message)
  flush(stdout)
end

function row_mask(row)
  mask = BigInt(0)
  for (index, value) in row
    if isodd(value)
      mask = xor(mask, BigInt(1) << (index - 1))
    end
  end
  return mask
end

function insert_mask!(pivots::Dict{Int, BigInt}, row::BigInt)
  while row != 0
    pivot = ndigits(row, base = 2) - 1
    if !haskey(pivots, pivot)
      pivots[pivot] = row
      return true
    end
    row = xor(row, pivots[pivot])
  end
  return false
end

function relation_rows(context)
  return vcat(collect(context.M.bas_gens), collect(context.M.rel_gens))
end

function quotient_state(context, s_columns)
  pivots = Dict{Int, BigInt}()
  for column in s_columns
    insert_mask!(pivots, BigInt(1) << (column - 1))
  end
  masks = BigInt[]
  for row in relation_rows(context)
    mask = row_mask(row)
    push!(masks, mask)
    insert_mask!(pivots, mask)
  end
  dimension = length(context.FB.ideals) - length(pivots)
  return dimension, masks, length(pivots)
end

function s_factor_base_columns(order, context, rational_primes)
  columns = Int[]
  missing = Dict{String, Vector{String}}()
  for rational_prime in rational_primes
    decomposed = prime_decomposition(order, rational_prime)
    local_missing = String[]
    for (prime_ideal, _exponent) in decomposed
      index = findfirst(==(prime_ideal), context.FB.ideals)
      if index === nothing
        push!(local_missing, string(prime_ideal))
      else
        push!(columns, index)
      end
    end
    if !isempty(local_missing)
      missing[string(rational_prime)] = local_missing
    end
  end
  unique!(sort!(columns))
  return columns, missing
end

function checkpoint!(path, status, context, s_columns, missing, started, settings)
  dimension, masks, augmented_rank = quotient_state(context, s_columns)
  rows = relation_rows(context)
  data = Dict(
    "schema" => "elliptic-curves.elkies-2026-hecke-s-class-monitor.v1",
    "status" => status,
    "claim_boundary" => [
      "Rows are exact principal-ideal relations produced by Hecke.",
      "The displayed quotient is only for the materialized factor base modulo 2 and the declared S columns.",
      "No factor-base generation proof, unit group, local Selmer conditions, Selmer dimension, or rank upper bound is supplied.",
    ],
    "settings" => settings,
    "factor_base_columns" => length(context.FB.ideals),
    "s_columns" => s_columns,
    "s_prime_ideals_absent_from_factor_base" => missing,
    "principal_relation_count" => length(rows),
    "principal_relation_masks_hex" => [string(mask, base = 16) for mask in masks],
    "relation_generators" => [
      string(generator)
      for generator in vcat(context.R_gen, context.R_rel)
    ],
    "augmented_mod2_rank" => augmented_rank,
    "factor_base_mod2_quotient_dimension" => dimension,
    "elapsed_seconds" => time() - started,
  )
  mkpath(dirname(path))
  temporary = path * ".tmp"
  open(temporary, "w") do io
    TOML.print(io, data; sorted = true)
  end
  mv(temporary, path; force = true)
  return dimension, length(rows), augmented_rank
end

function collect_s_targeted_relations!(
  context,
  s_columns;
  valuation_radius,
  precision,
  enumeration_limit,
)
  s_column_set = Set(s_columns)
  s_ideals = sort(
    [context.FB.ideals[index] for index in s_columns];
    by = ideal -> norm(ideal),
  )
  target_columns = [
    index for index in eachindex(context.FB.ideals)
    if !(index in s_column_set)
  ]
  isempty(s_ideals) && error("targeted S-collection needs visible S ideals")
  field_order = order(context.FB.ideals[1])
  norm_precision = Hecke.nbits(isqrt(abs(discriminant(field_order)))) + 30
  round = 0
  while true
    round += 1
    for (target_offset, target_column) in enumerate(target_columns)
      s_offset = mod1(round + target_offset - 1, length(s_ideals))
      search_ideal = context.FB.ideals[target_column] * s_ideals[s_offset]
      enumerator = Hecke.class_group_small_real_elements_relation_start(
        context,
        search_ideal;
        val = valuation_radius,
        prec = precision,
        limit = enumeration_limit,
      )
      attempts = 0
      while attempts < 8
        element = Hecke.class_group_small_real_elements_relation_next(enumerator)
        effective_norm = abs(
          Hecke.norm_div(element, norm(search_ideal), norm_precision)
        )
        if Hecke.nbits(numerator(effective_norm)) > norm_precision - 10 ||
           enumerator.restart > 2
          break
        end
        attempts += 1
        if Hecke.class_group_add_relation(
          context,
          element,
          effective_norm,
          norm(search_ideal);
          orbit = false,
          integral = true,
        )
          break
        end
      end
    end
  end
end

function main()
  settings = parse_args(ARGS)
  bound = parse(Int, settings["--factor-base-bound"])
  timeout_seconds = parse(Float64, settings["--timeout-seconds"])
  method = parse(Int, settings["--method"])
  valuation_radius = parse(Int, settings["--valuation-radius"])
  precision = parse(Int, settings["--precision"])
  enumeration_limit = parse(Int, settings["--enumeration-limit"])
  random_seed = parse(Int, settings["--random-seed"])
  include_s_in_factor_base = lowercase(settings["--include-s-in-factor-base"])
  include_s_in_factor_base in ("true", "false") ||
    error("--include-s-in-factor-base must be true or false")
  method in (2, 3, 4) || error("--method must be 2, 3, or 4")
  bound >= 2 || error("--factor-base-bound must be at least 2")
  timeout_seconds > 0 || error("--timeout-seconds must be positive")
  valuation_radius >= 0 || error("--valuation-radius must be nonnegative")
  precision >= 32 || error("--precision must be at least 32 bits")
  enumeration_limit >= 0 || error("--enumeration-limit must be nonnegative")
  Random.seed!(random_seed)
  s_primes = [ZZ(value) for value in split(settings["--s-primes"], ",")]

  rational_polynomials, x = polynomial_ring(QQ, "x")
  coefficients = [
    ZZ(value)
    for value in split(settings["--polynomial-ascending"], ",")
  ]
  length(coefficients) == 4 || error("the polynomial must have four coefficients")
  coefficients[4] == 1 || error("the cubic must be monic")
  polynomial = sum(coefficients[index + 1] * x^index for index in 0:3)
  field, _generator = number_field(polynomial, "a")
  order = maximal_order(field)
  started = time()
  emit("$PROTOCOL|stage=class_group_init|status=start|bound=$bound")
  if include_s_in_factor_base == "true"
    standard_factor_base = Hecke.NfFactorBase(order, bound; complete = false)
    ideals = copy(standard_factor_base.ideals)
    for rational_prime in s_primes
      for (prime_ideal, _exponent) in prime_decomposition(order, rational_prime)
        prime_ideal in ideals || push!(ideals, prime_ideal)
      end
    end
    # Initializing directly with the augmented base makes Hecke estimate an
    # automorphism heuristic by enumerating every rational prime up to the
    # largest (possibly enormous) S-prime.  Initialize the field-dependent
    # machinery on the small base, then rebuild the relation module on the
    # exact augmented base.  Only the ordinary-base generator relations are
    # needed because all added columns are killed in the displayed quotient.
    context = Hecke.class_group_init(
      standard_factor_base; add_rels = false, use_aut = false
    )
    context.FB = Hecke.NfFactorBase(order, ideals)
    context.M = Hecke.ModuleCtx_fmpz(length(context.FB.ideals))
    empty!(context.R_gen)
    empty!(context.R_rel)
    empty!(context.RS)
    context.bad_rel = 0
    context.rel_cnt = 0
    context.last = 0
    field_degree = degree(order)
    for ideal in standard_factor_base.ideals
      first_generator = ideal.gen_one
      Hecke.class_group_add_relation(
        context,
        Hecke.nf(order)(first_generator),
        QQ(abs(first_generator)^field_degree),
        ZZ(1);
        orbit = false,
      )
      second_generator = Hecke.nf(order)(ideal.gen_two)
      second_norm = Hecke.norm_div(second_generator, ZZ(1), 600)
      if Hecke.nbits(numerator(second_norm)) < 550
        Hecke.class_group_add_relation(
          context,
          second_generator,
          abs(second_norm),
          ZZ(1);
          orbit = false,
        )
      end
    end
  else
    context = Hecke.class_group_init(
      order, bound; complete = false, use_aut = false
    )
  end
  context.B2 = bound * 1000
  s_columns, missing = s_factor_base_columns(order, context, s_primes)
  dimension, relation_count, augmented_rank = checkpoint!(
    settings["--checkpoint"], "initialized", context, s_columns, missing,
    started, settings,
  )
  emit(
    "$PROTOCOL|stage=class_group_init|status=complete" *
    "|fb=$(length(context.FB.ideals))|S=$(length(s_columns))" *
    "|relations=$relation_count|augmented_rank=$augmented_rank" *
    "|quotient_dimension=$dimension"
  )

  if dimension == 0
    checkpoint!(settings["--checkpoint"], "factor_base_mod2_quotient_closed", context,
                s_columns, missing, started, settings)
    return
  end

  collector = Threads.@spawn begin
    if method == 2
      Hecke.class_group_find_relations2(
        context;
        val = valuation_radius,
        prec = precision,
        limit = enumeration_limit,
      )
    elseif method == 3
      context.expect = Hecke.class_group_expected(context, 100)
      Hecke.class_group_via_lll(context)
    else
      include_s_in_factor_base == "true" ||
        error("--method 4 requires --include-s-in-factor-base true")
      collect_s_targeted_relations!(
        context,
        s_columns;
        valuation_radius = valuation_radius,
        precision = precision,
        enumeration_limit = enumeration_limit,
      )
    end
  end
  previous_dimension = dimension
  previous_relations = relation_count
  monitor_error_reported = false
  while !istaskdone(collector)
    sleep(1)
    current_relations = length(context.R_gen) + length(context.R_rel)
    if current_relations != previous_relations
      try
      dimension, relation_count, augmented_rank = quotient_state(
        context, s_columns
      )
      checkpoint!(settings["--checkpoint"], "collecting", context, s_columns,
                  missing, started, settings)
      emit(
        "$PROTOCOL|stage=relations|status=progress" *
        "|relations=$relation_count|augmented_rank=$augmented_rank" *
        "|quotient_dimension=$dimension|seconds=$(time()-started)"
      )
      previous_dimension = dimension
      previous_relations = relation_count
      if dimension == 0
        checkpoint!(settings["--checkpoint"], "factor_base_mod2_quotient_closed",
                    context, s_columns, missing, started, settings)
        emit("$PROTOCOL|stage=complete|status=factor_base_mod2_quotient_closed")
        exit(0)
      end
      catch exception
        if !monitor_error_reported
          message = sprint(showerror, exception)
          emit("$PROTOCOL|stage=monitor|status=retry|error=$(replace(message, '|' => '/'))")
          monitor_error_reported = true
        end
      end
    end
    if time() - started >= timeout_seconds
      checkpoint!(settings["--checkpoint"], "strict_wall_timeout", context,
                  s_columns, missing, started, settings)
      emit("$PROTOCOL|stage=complete|status=strict_wall_timeout")
      exit(124)
    end
  end

  if istaskfailed(collector)
    try
      fetch(collector)
    catch exception
      checkpoint!(settings["--checkpoint"], "collector_error", context,
                  s_columns, missing, started, settings)
      emit("$PROTOCOL|stage=complete|status=collector_error|error=$(typeof(exception))")
      rethrow()
    end
  end
  checkpoint!(settings["--checkpoint"], "collector_returned", context,
              s_columns, missing, started, settings)
  emit("$PROTOCOL|stage=complete|status=collector_returned")
end

main()
