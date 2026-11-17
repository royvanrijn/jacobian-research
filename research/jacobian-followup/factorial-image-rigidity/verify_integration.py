"""Check this packet against the real ledger; report unrelated global failures."""
from pathlib import Path
import hashlib
import json
import re
import shutil
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
import render_status
import research

ids = {'FIRGAPTHREEPOLAR20261117','FIRGAPTHREEPIVOT20261117','FIRGAPTHREEREC20261117','FIRGAPTHREEREAL20261117','FIRGAPTHREEFIFTH20261117','FIRALLGAPINITIAL20261117','FIRGENERALGAPADJOINT20261117','FIRGENERALGAPRANK20261117','FIRCUBICGAMMASIMPLE20261117','FIRCUBICALLSFC20261117','FIRTILTEDSHIFT20261117','FIRTILTEDPRIMEPOWER20261117','FIRCUBICFOLDS20261117','FIRLOWERONESIMPLE20261117','FIRCUBICPRIMECONF20261117','FIRBESSELUNITTRANSFER20261117','FIRDYADICBINSFC20261117','FIRPADEPARABOLA20261117','FIRUNIFORMADJACENTSFC20261117','FIRPADESEPARATION20261117','FIRADJACENTMONOSFC20261117','FIRUNIFORMTWOSIMPLE20261117','FIRALLUNIVTWOSFC20261117','FIREVENTUALSIMPLE620261117','FIRALLINDEXSFC620261117','FIRANALYTICANNULUS20261117', 'FIRRIGIDITY20261116', 'FIRSFCFAMILY20261116', 'FIRGAMMAADAPTER20261116',
       'FIRFLATAFFINE20261116', 'OP-FACTORIAL-GENERAL20261116',
       'FIREIT20261116','FIRFC2NONFLAT20261116','FIRFC2FULL20261116',
       'FIRRADIAL2RINGS20261116','FIRRADIAL2NORM20261116',
       'FIRCHAINMOMENT20261116','FIRSIMPLEXLINE20261116','FIRPENCILGAMMA20261116',
       'FIRAFFINESIMPLEX20261116','FIRFLATNEXTLINE20261116','FIRPENCILEXEC20261116',
       'FIRTWOANGLINEAR20261116','FIRTWOANGSHAPES20261116',
       'FIRTWOANGGAMMA20261116','FIRTWOANGEXEC20261116',
       'FIRRECIPWINDOW20261116','FIRMARKEDWINDOW20261116','FIRSFCAFFINE20261116',
       'FIRSFCPRIME20261116','FIRWINDOWEXEC20261116',
       'FIRPADICDEFORM20261116','FIRSFCMONOPRIME20261116','FIRPRIMEMARKERS20261116',
       'FIRPRIMEPAIRI20261116','FIRPRIMEWINDOWEXEC20261116',
       'FIRNEWTONLIFT20261116','FIRINDEXEDMIXEDSFC20261116','FIRMULTIPRIMEWINDOW20261116',
       'FIRMIXEDPAIR20261116','FIRMIXEDEXEC20261116','FIRREALROOTPARITY20261116',
       'FIRSLOPECAP20261116','FIRSHIFTEDINDEX20261116','FIRSLOPECOVEREXEC20261116','FIRLARGESQUARE20261116',
       'FIRPRIMEGAP20261116','FIRFOURTEENTH20261116','FIRSMALLKERNEL20261116','FIRGAPKERNELFAMILY20261116',
       'FIREFFECTIVECUTOFF20261116','FIRRESIDUALDEGREES20261116','FIRBOUNDEDMULTIPLIER20261116',
       'FIRDENSITYONE20261116','FIRALLMLOCAL20261116','FIRQUANTDENSITY20261116',
       'FIRCOMMONOFFSETS20261117','FIRSIMULTLOG20261117','FIRSIMULTPOWER20261117',
       'FIRLARGEM20261117','FIRLOGEXPCANDIDATES20261117','FIRALLEXPLOG20261117','FIRALLEXPPOWER20261117',
       'FIRGAMMAMONOWINDOW20261116','FIRDIRMONOWINDOW20261116','FIRMONOWINDOWIDEAL20261116',
       'FIRBLOCKSUMSFC20261116','FIRHYPERWINDOWEXEC20261116',
       'FIRBLOCKRADIAL20261116','FIRMONOMIALTOPFC20261116','FIRCOUPLEDPERIOD20261116','FIRCOUPLEDMARKER20261116',
       'FIRCOUPLEDBLOCKFC20261116','FIRCOUPLEDBLOCKRINGS20261116','FIRCOUPLEDEXEC20261116'}
documents = ('GAP_THREE_POLAR_AND_CONTINUANT.tex','GAP_THREE_WRONSKIAN_RECURRENCE.tex','HIGHER_GAP_ADJOINT_RANK.tex','CUBIC_ADJOINT_STABILITY_SFC.tex','TILTED_CUBIC_PRIME_POWER_AND_FOLD.tex','DYADIC_BINOMIAL_SFC.tex','UNIFORM_PADE_PARABOLA_SFC.tex','PADE_ADJACENT_MONOMIALS.tex','ALL_UNIVARIATE_TWO_TERM_SFC.tex','EFFECTIVE_EXPONENT_SIX_SIMPLICITY.tex','EFFECTIVE_SADDLE_ANNULUS.tex','RIGIDITY_CRITICAL_VALUES.tex','STATUS_AND_LITERATURE.tex','LONG_EIT_AND_FC2_AUDIT.tex',
             'SIMPLEX_PENCILS_AND_AFFINE_LIMITS.tex','TWO_ANGULAR_COORDINATES.tex',
             'SFC_WINDOWS_AND_ROCKS_AUDIT.tex','PRIME_POWER_TWO_TERM_SFC.tex',
             'MIXED_PRIME_RESIDUAL_WINDOWS.tex','HYPERGEOMETRIC_MONOMIAL_WINDOWS.tex',
             'COUPLED_BLOCK_PERIODS.tex','SLOPE_COVERAGE_AND_SHIFTED_WINDOWS.tex','PRIME_GAPS_AND_SMALL_KERNELS.tex',
             'EFFECTIVE_CUTOFFS_AND_RESIDUAL_DEGREES.tex','DENSITY_ONE_SQUAREFREENESS.tex','QUANTITATIVE_SQUAREFREENESS_ALL_EXPONENTS.tex','SIMULTANEOUS_EXPONENTS_AND_POWER_SAVING.tex','ALL_EXPONENTS_AND_LARGE_EXPONENT_CUTOFF.tex')
data = json.loads((ROOT/'MATH_STATUS.json').read_text())
known = {e['id']:e for e in data['entries']}
for id in ids:
    row = known[id]
    assert render_status.REQUIRED_FIELDS <= row.keys()
    assert not row.keys() - (render_status.REQUIRED_FIELDS | render_status.OPTIONAL_FIELDS)
    assert row['kind'] in render_status.KINDS and row['state'] in render_status.STATES
    assert row['proof_type'] in render_status.PROOF_TYPES
    assert all(d in known or d.startswith('external: ') for d in row['dependencies'])
    assert (ROOT/row['canonical_source']).is_file()
    assert all((ROOT/p).is_file() for p in row['software_lock'])
    assert not any(row[k] for k in ('independent_replay','formal_verification','external_review'))
    if row['checker']:
        assert row['artifact_hash'] == 'sha256:'+hashlib.sha256((ROOT/row['checker']).read_bytes()).hexdigest()
    else:
        assert row['kind'] == 'open_problem' and row['state'] == 'open'
        assert row['proof_type'] == 'not-applicable' and row['artifact_hash'] is None

lessons = json.loads((ROOT/'knowledge/lessons.json').read_text())['lessons']
actions = json.loads((ROOT/'knowledge/work_items.json').read_text())['actions']
for row in [next(x for x in lessons if x['id']=='METHOD-FACTORIAL-RADIAL-PURITY'),
            next(x for x in actions if x['id']=='WORK-FACTORIAL-GENERAL')]:
    assert ids <= set(row['claims'])
    assert row['reviewed_claims'] == {id:research.claim_fingerprint(known[id]) for id in row['claims']}
    assert all((ROOT/p).is_file() for p in row['sources'])

ledger = json.loads((HERE/'source_ledger.json').read_text())
for row in ledger['cached_primary_sources']:
    assert hashlib.sha256((HERE/row['file']).read_bytes()).hexdigest() == row['sha256']
for filename in documents:
    text = (HERE/filename).read_text()
    for link in re.findall(r'\\href\{([^}]+)\}', text):
        assert link.startswith('https://') or (HERE/link).is_file(), (filename,link)

receipt = {'date':'2026-11-17', 'packet_integration':'pass', 'new_claims':sorted(ids),
           'source_pins':'pass', 'local_tex_links':'pass', 'latex_syntax':[],
           'general_FC_and_SFC_closed':False, 'current_goal_turn_classification':'progress',
           'previous_goal_turn_classification':'progress',
           'classification_evidence':'The previous saved draft and its 50 exact checks were genuine proof/transcription progress. The current turn completed the inactive-variable extension in the canonical block proof, explicitly controlling the inactive boundary and unbounded tail, and updated the zero-exponent monomial and coupled-family exporters, 539-check receipt and seven existing claim scopes. General FC/SFC remain separate open programme obligations.',
           'inactive_variable_extension_complete':True,
           'completion_scope':'Remove the requirement that every variable occurs in the leading monomial; retain arbitrary complex strictly lower-degree couplings and the proved nonzero block families.',
           'remaining':['general FC in dimension >=3 beyond closed classes','unrestricted SFC and full all-starting-index two-term exponent-six SFC','higher residual information or a simple-root theorem beyond the ODE multiplicity bound','arbitrary phases beyond the closed monomial/pencil/cubic/coupled-product classes','complete bibliography and remaining source audits','remaining EIT V geometric/elliptic and higher odd-degree inputs','coupled angular phases beyond product uniforms, especially arbitrary triangle phases','emptiness of every weighted cubic scale-exception set'],
           'global_navigation':[]}

# The elementary large-exponent cutoff removes the analytic exponent cap entirely.
all_exp=json.loads((HERE/'all_exponent_squarefree_verification.json').read_text())
assert all_exp['status']=='pass' and all_exp['large_exponent_criterion_proved']
assert all_exp['all_positive_exponents_share_one_exceptional_set'] and all_exp['overlap_threshold']=='2^144'
assert all_exp['prime_power_result_already_registered'] and all_exp['exponent_one_explicit']
assert all_exp['published_all_exponent_log_bound'] and all_exp['qualified_all_exponent_power_bound']
assert all_exp['MR_II_source_qualification_retained'] and not all_exp['all_index_SFC_closed']
assert not all_exp['numeric_analytic_starting_threshold_extracted']
for name,expected in all_exp['files'].items():
    assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==expected,name
assert all(known[id]['state']=='proved' for id in ('FIRLARGEM20261117','FIRLOGEXPCANDIDATES20261117','FIRALLEXPLOG20261117','FIRALLEXPPOWER20261117'))
all_exp_audit=json.loads((HERE/'all_exponent_squarefree_audit.json').read_text())
assert all_exp_audit['overlap_threshold_is_not_an_analytic_starting_threshold']
assert all_exp_audit['sparse_warning_and_Section10_exclusion_retained']
assert not all_exp_audit['complete_preprint_validation_claimed']
receipt['all_positive_exponents_refinement']={
    'claims':['FIRLARGEM20261117','FIRLOGEXPCANDIDATES20261117','FIRALLEXPLOG20261117','FIRALLEXPPOWER20261117'],
    'proof':'ALL_EXPONENTS_AND_LARGE_EXPONENT_CUTOFF.tex',
    'proof_sha256':all_exp['files']['ALL_EXPONENTS_AND_LARGE_EXPONENT_CUTOFF.tex'],
    'checks':all_exp['total'],'deterministic_criterion':'2^m>=4*n^2 => squarefree, equality allowed',
    'all_index_prime_power_result':'already registered; m1 explicit r0 case',
    'remaining_bad_pair_conditions':['2^m<4*n^2','at least two distinct prime divisors in m'],
    'overlap_threshold':'2^144; capped and unrestricted bad sets agree beyond it',
    'transfer':'D<=D_all<=D+2^144','exponent_cap_remaining':False,
    'published_route':'D_all(X)<<X/(log X)^(1/50), independently of MR II',
    'qualified_power_route':'D_all(X)<<X^(1-10^-18), dense MR II preprint input retained',
    'all_exponent_reciprocal_sum_converges_using_power_route':True,
    'numeric_analytic_starting_threshold_extracted':False,'small_index_pointwise_audit_claimed':False,
    'all_index_SFC_closed':False,'independent_replay':False,'formal_verification':False,'external_review':False
}

# Simultaneous exponents: independent published-log and source-qualified dense-power routes.
simult=json.loads((HERE/'simultaneous_squarefree_verification.json').read_text())
assert simult['status']=='pass' and simult['published_2016_fallback_proved']
assert simult['preprint_dense_power_application_checked'] and simult['MR_II_preprint_source_qualification_retained']
assert simult['final_epsilon']=={'numerator':1,'denominator':10**18}
assert not any(simult[key] for key in ('complete_preprint_validation_claimed','sparse_Section10_used',
    'Theorem1_7_advertised_exponent_used','analytic_exception_multiplied_by_M','all_index_SFC_closed'))
for name,expected in simult['files'].items():
    assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==expected,name
assert all(known[id]['state']=='proved' for id in ('FIRCOMMONOFFSETS20261117','FIRSIMULTLOG20261117','FIRSIMULTPOWER20261117'))
def simultaneous_dependency_closure(id):
    result=set()
    def collect(current):
        for dep in known[current]['dependencies']:
            if dep in result:continue
            result.add(dep)
            if dep in known:collect(dep)
    collect(id);return result
published_route=simultaneous_dependency_closure('FIRSIMULTLOG20261117')
assert 'FIRSIMULTPOWER20261117' not in published_route and not any('MR II' in dep for dep in published_route)
power_route=simultaneous_dependency_closure('FIRSIMULTPOWER20261117')
assert any('MR II' in dep and 'PREPRINT' in dep for dep in power_route)
all_log_route=simultaneous_dependency_closure('FIRALLEXPLOG20261117')
assert 'FIRSIMULTPOWER20261117' not in all_log_route and not any('MR II' in dep for dep in all_log_route)
all_power_route=simultaneous_dependency_closure('FIRALLEXPPOWER20261117')
assert 'FIRSIMULTPOWER20261117' in all_power_route and any('MR II' in dep and 'PREPRINT' in dep for dep in all_power_route)

simult_audit=json.loads((HERE/'simultaneous_squarefree_audit.json').read_text())
assert simult_audit['warning_source']['publisher_warning_checked_online']
assert not simult_audit['complete_preprint_reconstruction_or_validation_claimed']
assert not simult_audit['corrected_MR_II_version_obtained'] and not simult_audit['published_MR_II_claimed']
receipt['simultaneous_growing_exponents_refinement']={
    'claims':['FIRCOMMONOFFSETS20261117','FIRSIMULTLOG20261117','FIRSIMULTPOWER20261117'],
    'proof':'SIMULTANEOUS_EXPONENTS_AND_POWER_SAVING.tex',
    'proof_sha256':simult['files']['SIMULTANEOUS_EXPONENTS_AND_POWER_SAVING.tex'],
    'checks':simult['total'],'range':'all integer m<=n^(1/16), one common exceptional set',
    'published_2016_fallback':'D(X)<<X/(log X)^(1/50), independent of MR II',
    'additional_dense_preprint_route':'D(X)<<X^(1-10^-18), absolute implied constant',
    'source_qualification':'MR II v1 dense9.2(ii) plus ordinary linear sieve; published sparse warning retained',
    'sparse_Section10_or_Lemma10_3_used':False,'Theorem1_7_advertised_exponent_used':False,
    'complete_preprint_validation_claimed':False,'analytic_exception_multiplied_by_M':False,
    'reciprocal_sum_converges_using_power_route':True,
    'numeric_implied_constant_or_starting_threshold_extracted':False,
    'all_index_SFC_closed':False,'finite_controls_are_not_numerical_density_evidence':True,
    'independent_replay':False,'formal_verification':False,'external_review':False
}

# Quantitative all-fixed-exponent coverage is direct, with controlled growing offsets.
quant=json.loads((HERE/'quantitative_squarefree_verification.json').read_text())
assert quant['status']=='pass' and quant['all_fixed_positive_integer_exponents']
assert quant['explicit_logarithmic_exponent']=='1/50' and quant['growing_offsets_bad_prime_union_controlled']
assert not any(quant[key] for key in ('old_density_one_theorem_used','bootstrap_used','parity_used',
                                    'uniform_in_growing_m_claimed','all_index_SFC_closed',
                                    'numerical_constant_or_starting_threshold_extracted'))
for name,expected in quant['files'].items():
    assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==expected,name
assert all(known[id]['state']=='proved' for id in ('FIRALLMLOCAL20261116','FIRQUANTDENSITY20261116'))
ancestors=set()
def collect_quant_dependencies(id):
    for dep in known[id]['dependencies']:
        if dep in known and dep not in ancestors:
            ancestors.add(dep);collect_quant_dependencies(dep)
collect_quant_dependencies('FIRQUANTDENSITY20261116')
assert not ancestors.intersection({'FIRDENSITYONE20261116','FIRREALROOTPARITY20261116','FIRPRIMEGAP20261116'})
quant_audit=json.loads((HERE/'quantitative_squarefree_audit.json').read_text())
assert quant_audit['H']=='floor(X^(1/8))' and quant_audit['Y']=='sqrt(2mX)' and quant_audit['MR_delta']=='1/4'
assert all(row['imported_not_reproved'] for row in quant_audit['analytic_inputs'])
assert not quant_audit['all_offsets_squarefree_assumed'] and not quant_audit['all_C_smoothness_required']
receipt['quantitative_all_fixed_exponents_refinement']={
    'claims':['FIRALLMLOCAL20261116','FIRQUANTDENSITY20261116'],
    'proof':'QUANTITATIVE_SQUAREFREENESS_ALL_EXPONENTS.tex',
    'proof_sha256':quant['files']['QUANTITATIVE_SQUAREFREENESS_ALL_EXPONENTS.tex'],
    'checks':quant['total'],'rate':quant['quantitative_rate'],'fixed_m':'every positive integer',
    'H':quant_audit['H'],'Y':quant_audit['Y'],'MR_delta':'1/4',
    'bad_prime_progression_union':'O_m(X^(7/8))','MR_power_exception':'X^(799/800)*(log X)^(1/3)',
    'quantitative_prime_gap_exception':'O(X^(2/3)), Filaseta-Rocks equation3.1',
    'old_density_one_theorem_used':False,'bootstrap_used':False,'parity_used':False,
    'numerical_constant_or_starting_threshold_extracted':False,'uniform_in_growing_m_claimed':False,
    'all_index_SFC_closed':False,'finite_controls_are_not_numerical_density_evidence':True,
    'independent_replay':False,'formal_verification':False,'external_review':False
}

# Ordinary natural density one uses a written proof with imported analytic inputs.
density=json.loads((HERE/'density_one_verification.json').read_text())
assert density['status']=='pass' and density['ordinary_natural_density_one_proved']
assert density['fixed_H_then_X_then_H_order_checked'] and density['finite_bad_discriminant_thresholds_respected']
assert density['analytic_uniformity_in_X_dependent_function_checked']
assert not density['all_index_SFC_closed'] and not density['irreducibility_density_one_claimed']
for name,expected in density['files'].items():
    assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==expected,name
assert known['FIRDENSITYONE20261116']['state']=='proved'
density_audit=json.loads((HERE/'density_one_audit.json').read_text())
assert density_audit['density_type']=='ordinary natural density one'
assert density_audit['finite_controls_are_not_numerical_density_evidence']
assert not density_audit['all_C_smoothness_required'] and not density_audit['all_offsets_squarefree_assumed']
assert all(row['published_input_imported_not_independently_reproved'] for row in density_audit['analytic_inputs'])
receipt['ordinary_natural_density_one_refinement']={
    'claim':'FIRDENSITYONE20261116','proof':'DENSITY_ONE_SQUAREFREENESS.tex',
    'proof_sha256':density['files']['DENSITY_ONE_SQUAREFREENESS.tex'],
    'checks':density['total'],'exceptional_count':'o(X)',
    'density_type':'ordinary natural density one',
    'analytic_inputs':'Filaseta-Rocks Theorem1.1/equation3.3; uniform Matomaki-Radziwill Theorem1; classical PNT',
    'source_inputs_imported_not_reproved':True,'X_dependent_function_uniformity_checked':True,
    'limit_order':density_audit['limit_order'],'finite_bad_discriminant_thresholds_respected':True,
    'full_all_C_smoothness_required':False,'fixed_multiplier_or_offset_cutoff':False,
    'remaining':'possible density-zero exceptional set; full all-index SFC',
    'irreducibility_density_one_claimed':False,'all_index_SFC_closed':False,
    'finite_controls_are_not_numerical_density_evidence':True,
    'independent_replay':False,'formal_verification':False,'external_review':False
}

# The explicit cutoff closes effectiveness; full residual degree coverage is open.
effective=json.loads((HERE/'effective_residual_verification.json').read_text())
assert effective['status']=='pass' and effective['effective_fixed_multiplier_cutoff_computed']
assert not effective['all_index_SFC_closed'] and not effective['uniform_nearby_prime_coverage_proved']
for name,expected in effective['files'].items():
    assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==expected,name
assert all(known[id]['state']=='proved' for id in ('FIREFFECTIVECUTOFF20261116','FIRRESIDUALDEGREES20261116'))
assert effective['cutoff_examples']['15']==30770 and effective['cutoff_examples']['100']==4082 and effective['cutoff_examples']['1000']==401
assert effective['C_positive_iff_characteristic_zero_squarefree_for_this_sequence']
assert effective['full_residual_divisor_degree_sets_retained']
assert known['FIRBOUNDEDMULTIPLIER20261116']['state']=='proved'
assert effective['complete_bounded_multiplier_family']=={
    'multipliers':[1,24000],'offsets':[0,12],'prime_condition':'p>6max(a,i)',
    'all_r_ge1':True,'finite_exceptions_remaining':False}
finite=json.loads((HERE/'initial_prime_gaps_500000.json').read_text())
assert finite['limit']==500000 and finite['prime_count']==41538
assert len(finite['primes'])==41538 and finite['maximum_consecutive_prime_gap']==114
assert finite['maximum_gap_endpoints']==[[492113,492227]] and finite['maximum_n_minus_preceding_prime']==113
assert effective['exact_endpoint_comparisons_verified']

receipt['effective_cutoff_and_residual_degree_refinement']={
    'proof':'EFFECTIVE_CUTOFFS_AND_RESIDUAL_DEGREES.tex',
    'proof_sha256':effective['files']['EFFECTIVE_CUTOFFS_AND_RESIDUAL_DEGREES.tex'],
    'checks':effective['total'],'cutoff_formula':'ceil(max(400000,exp(sqrt(a/50)))/(a-2)), a>=15',
    'exact_cutoff_examples':effective['cutoff_examples'],'effective_eventual_cutoff_computed':True,
    'finite_exception_domain':'cutoff criterion alone: prime p>6a and p^r<P0(a); none remain at a<=24000',
    'complete_bounded_multiplier_family':effective['complete_bounded_multiplier_family'],
    'finite_prime_gap_replay':effective['finite_prime_gap_replay'],
    'C_positivity_target':'equivalent to all-index characteristic-zero squarefreeness for this sequence; not proved unboundedly',
    'C_smoothness_target':'stronger modular target; not proved for all offsets',
    'full_residual_divisor_degree_sets_retained':True,
    'uniform_nearby_prime_coverage_proved':False,'eventual_coverage_needs_initial_range_for_all_index':True,
    'fixed_prime_beyond_prefix_example':effective['fixed_prime_offset13_example'],
    'all_index_SFC_closed':False,'source_input':'pinned Dusart2010 Proposition6.8, unconditional imported estimate',
    'independent_replay':False,'formal_verification':False,'external_review':False
}

# Prime-gap bounds and the retained integer kernels have independent exact evidence.
gap=json.loads((HERE/'prime_gap_verification.json').read_text())
assert gap['status']=='pass' and gap['all_kernel_remainders_one']
assert not gap['all_index_SFC_closed'] and not gap['eventual_fixed_multiplier_cutoff_computed']
for name,expected in gap['files'].items():
    assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==expected,name
assert all(known[id]['state']=='proved' for id in ('FIRPRIMEGAP20261116','FIRFOURTEENTH20261116','FIRSMALLKERNEL20261116','FIRGAPKERNELFAMILY20261116'))
kernels=json.loads((HERE/'small_kernel_certificates.json').read_text())
assert [r['offset'] for r in kernels['certificates']]==list(range(2,13))
assert all(r['remainder_after_small_primes']=='1' for r in kernels['certificates'])
receipt['prime_gap_kernel_refinement']={
    'proof':'PRIME_GAPS_AND_SMALL_KERNELS.tex','proof_sha256':gap['files']['PRIME_GAPS_AND_SMALL_KERNELS.tex'],
    'checks':gap['total'],'preceding_prime_bound':'even degree<=2*floor((n-q_n)/4) for n>=30',
    'explicit_bound':'degree<n/14 for n>=30','retained_kernel_offsets':list(range(2,13)),
    'all_integer_kernel_remainders_one':True,'uniform_rectangle':gap['uniform_shifted_rectangle'],
    'effective_eventual_cutoff_computed_in_this_earlier_stage':False,'now_effective_via':'FIREFFECTIVECUTOFF20261116','all_index_SFC_closed':False,
    'source_inputs':'pinned Nagura1952; classical PNT verified at official DLMF27.2.3',
    'independent_replay':False,'formal_verification':False,'external_review':False
}

# New shifted families are proved for all powers; the cap scan stays finite.
slope = json.loads((HERE/'slope_coverage_verification.json').read_text())
assert slope['status']=='pass' and not slope['all_index_exponent_six_closed']
for name,expected in slope['files'].items():
    assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==expected,name
assert all(known[id]['state']=='proved' for id in ('FIRSLOPECAP20261116','FIRSHIFTEDINDEX20261116','FIRSLOPECOVEREXEC20261116','FIRLARGESQUARE20261116'))
probe=json.loads((HERE/'slope_coverage_probe.json').read_text())
assert not probe['all_indices_claimed'] and not probe['residual_factorization_used']
assert [row['terminal_index'] for row in probe['certificates']]==list(range(2,65))
assert probe['survivor_indices']==[]
receipt['slope_and_shifted_refinement']={
    'universal_proof':'SLOPE_COVERAGE_AND_SHIFTED_WINDOWS.tex',
    'proof_sha256':slope['files']['SLOPE_COVERAGE_AND_SHIFTED_WINDOWS.tex'],
    'checks':slope['total'],'slope_cap_shape_one':'d*n;6*n for exponent6',
    'candidate_primes':'p<=2*d*k OR p divides n(n-1)...(n-2k+1)',
    'shifted_polygon':'N=u*p^r+i, p>d*max(u,i); horizontal i plus edge denominator p^r',
    'infinite_family':slope['infinite_shifted_family'],
    'uniform_large_degree_bound':slope['uniform_large_square_degree_bound'],
    'fixed_kernel_discriminant_gcd':'gcd(disc F2,disc F3)=2^10*3^4*5^2',
    'full_square_obstruction':'constant1 plus ordinary factorial nonsquare leading coefficient',
    'historical_one_prime_pattern':'4*137^r+3 retains h=t1 under the older coarse test; stronger prime-gap/kernel stage now covers it',
    'finite_full_slope_scan':'n2..64; no survivors; no all-index extrapolation',
    'all_index_exponent_six_closed':False,'residual_tests_can_use_primes_above_cap':True,
    'independent_replay':False,'external_review':False,'formal_verification':False
}

# The focused SFC refinement has its own analytic proof and exact arithmetic controls.
mixed = json.loads((HERE/'mixed_prime_window_verification.json').read_text())
assert mixed['status']=='pass'
for name, expected in mixed['files'].items():
    assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==expected,name
assert known['FIRREALROOTPARITY20261116']['state']=='proved'
assert 'FIRREALROOTPARITY20261116' in known['FIRMULTIPRIMEWINDOW20261116']['dependencies']
receipt['focused_sfc_refinement'] = {
    'claim':'FIRREALROOTPARITY20261116', 'verification_checks':mixed['total'],
    'equivalence_and_real_root_proof':'MIXED_PRIME_RESIDUAL_WINDOWS.tex, real-place theorem',
    'one_prime_plus_real_place_example':{'n':7,'d':6,'prime':3,'local_degrees':[0,1],'filtered_degrees':[0]},
    'two_prime_even_degree_survivors':{'indices':[21,23],'primes':[2,3],'degrees':[0,2]},
    'all_index_exponent_six_closed':False, 'inactive_variable_limit_already_proved':True,
    'independent_replay':False,'formal_verification':False,'external_review':False
}

# Audit the current extension against its full scope, not just finite identities.
block = json.loads((HERE/'coupled_block_verification.json').read_text())
assert block['status'] == 'pass' and block['checks'] == 539
assert block['inactive_variable_extension_checked'] and not block['general_FC_and_SFC_closed']
for name, expected in block['files'].items():
    assert hashlib.sha256((HERE/name).read_bytes()).hexdigest() == expected, name
draft = json.loads((HERE/'../cross-programme-connections/check_inactive_variable_factorial_extension_20261116.json').read_text())
assert draft['exact_checks'] == 50 and draft['all_exact_checks_pass']
assert all(known[id]['state'] == 'proved' for id in ('FIRBLOCKRADIAL20261116','FIRMONOMIALTOPFC20261116','FIRCOUPLEDBLOCKFC20261116'))
assert 'a_i>=0' in known['FIRMONOMIALTOPFC20261116']['scope']
assert 'ACTIVE degree' in known['FIRBLOCKRADIAL20261116']['scope']
assert 'K(X,Y,0)=0' in known['FIRCOUPLEDBLOCKFC20261116']['scope']
proof_text = (HERE/'COUPLED_BLOCK_PERIODS.tex').read_text()
for formula in (r'Q_2(3,0,0)', r'\prod_{i\in I}(\alpha_i)_{a_i m}', r'm^{-r/2-B}', r'H_{[e]}(dT,Z)'):
    assert formula in proof_text, formula
receipt['inactive_variable_completion_audit'] = {
    'status':'pass',
    'universal_evidence':'COUPLED_BLOCK_PERIODS.tex, Theorem 1 and the single-top-monomial corollary; analytic proof, not a bounded moment scan',
    'proof_sha256':hashlib.sha256((HERE/'COUPLED_BLOCK_PERIODS.tex').read_bytes()).hexdigest(),
    'proof_review':[
        'Exact active-only tilt gives Gamma(d_j*m+A_j,1) radii; angles and inactive Gamma(beta_l,1) laws are unchanged.',
        'On compact positive active radii and bounded Z/m, the relative phase is 1+G(uT,Z/m)/(c*m*product u^d)+O(m^-2), uniformly bounded after mth powers. Fixed Z gives the exponent G(dT,0)/(c*Cd).',
        'Active-degree marker scaling is bounded by a fixed polynomial in original Z; its inactive Gamma expectation is integrable, so conditional concentration and dominated convergence give the marked local limit.',
        'On all scaled radii the absolute kernel converges uniformly on compacts to exp(-sum u-sum v)*product u^d. Its unique maximum is M=exp(-D)*Cd at (d,0); inactive zero boundaries do not require vanishing.',
        'The compact complement has strict gap rho<M and the unbounded tail satisfies K_m<=rho*exp(-W/2). Positive shape weights are locally integrable, even below one; fixed polynomial markers cost only polynomial powers.',
        'The scaled normalization is C*M^m*m^(-r/2-B), so the complete complementary integral is exponentially negligible. Singleton blocks with zero exponents then give the stated nonzero exponential for arbitrary complex lower coefficients.',
        'For block families G(X,0) is unchanged by K(X,Z) with K(X,0)=0. Powers send K to k*F0^(k-1)*K; algebraic-scale eventual nonzero, algebraic-shape complex-scale FC and ordinary rational tail-unit conclusions retain their original hypotheses.'
    ],
    'requirements':[
        {'requirement':'Every single nonconstant leading monomial, allowing absent variables and arbitrary complex lower terms','evidence':'FIRMONOMIALTOPFC20261116; active singleton proof with a_i>=0'},
        {'requirement':'Proposed X^3+Q<=2 limit exp(Q2(3,0,0)/27)','evidence':'Canonical equation (1a); exact cubic exponent controls and retained 50 checks'},
        {'requirement':'Active m scale, fixed inactive Gamma scale, and disappearance of lower inactive couplings','evidence':'Exact tilt, uniform local expansion and full complementary estimate in Theorem 1'},
        {'requirement':'Additional variables with genuine coupling in proved block families','evidence':'Equation (6) permits every homogeneous K of degree D-1 vanishing at inactive0 plus arbitrary R<=D-2; exact power controls'},
        {'requirement':'All positive shapes and polynomial markers remain valid','evidence':'Theorem 1 uses active degree and unchanged inactive moments; subunit shapes and zero-leading/nonzero-moment controls'},
        {'requirement':'Maintained theorem, implementation, registry, source provenance and navigation agree','evidence':'Seven updated registered scopes, source pins, checker hashes, reviewed fingerprints and local links checked here'},
    ],
    'verification_checks':block['checks'], 'retained_draft_checks':50,
    'independent_replay':False, 'external_review':False, 'formal_verification':False,
    'general_FC_and_SFC_closed':False
}

annulus=json.loads((HERE/'saddle_annulus_verification.json').read_text())
assert annulus['status']=='pass' and annulus['all_roots_in_annulus_above_cutoff']
assert annulus['inner_disk_cutoff']==10**6 and annulus['cutoff_is_for_zero_exclusion_only']
assert not any(annulus[k] for k in ('eventual_simplicity_proved','simplicity_cutoff_extracted','finite_prefix_verified','all_index_SFC_closed','independent_replay','formal_verification','external_review'))
for name,digest in annulus['files'].items():assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==digest,('stale saddle-annulus receipt',name)
receipt['effective_saddle_annulus_refinement']={
    'claim':'FIRANALYTICANNULUS20261117','proof':'EFFECTIVE_SADDLE_ANNULUS.tex',
    'checks':annulus['total'],'root_annulus':'1<|z/n^6|<2800 for EVERY n>=10^6',
    'both_saddle_critical_values_excluded':True,
    'cutoff_is_for_zero_exclusion_only':True,
    'eventual_simplicity_proved':False,'finite_prefix_verified':False,
    'next_gate':'global contour contributions and effective differentiated error bounds on all root cancellations in the compact annulus',
    'independent_replay':False,'formal_verification':False,'external_review':False}

simplicity=json.loads((HERE/'exponent_six_simplicity_verification.json').read_text())
assert simplicity['status']=='pass' and simplicity['finite_prefix_complete']
assert simplicity['effective_simplicity_cutoff']==10000 and simplicity['all_index_exponent_six_two_term_SFC']
assert not any(simplicity[k] for k in ('unrestricted_SFC_closed','general_FC_closed','independent_replay','formal_verification','external_review'))
for name,digest in simplicity['files'].items():assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==digest,('stale exponent-six simplicity receipt',name)
receipt['effective_exponent_six_all_index_simplicity']={
    'claims':['FIREVENTUALSIMPLE620261117','FIRALLINDEXSFC620261117'],
    'proof':'EFFECTIVE_EXPONENT_SIX_SIMPLICITY.tex','checks':simplicity['total'],
    'effective_simplicity_cutoff':10000,'complete_finite_prefix':True,
    'ordinary_A_plus_B_Y6_SFC_all_starts':True,'unrestricted_SFC_closed':False,'general_FC_closed':False,
    'independent_replay':False,'formal_verification':False,'external_review':False}
for earlier_name in ('focused_sfc_refinement','slope_and_shifted_refinement'):
    earlier=receipt[earlier_name]
    if 'all_index_exponent_six_closed' in earlier:
        earlier['all_index_exponent_six_closed_in_this_earlier_stage']=earlier.pop('all_index_exponent_six_closed')
    earlier['now_closed_by']='FIRALLINDEXSFC620261117, ordinary Gamma shape1'
receipt['effective_saddle_annulus_refinement']['later_simplicity_cutoff']=10000
receipt['effective_saddle_annulus_refinement']['later_all_index_closure']='FIRALLINDEXSFC620261117'
receipt['remaining']=[item for item in receipt['remaining'] if 'two-term exponent-six' not in item]
receipt['remaining'].append('unrestricted SFC/general FC beyond completed classes; other composite exponents and general supports')

uniform=json.loads((HERE/'uniform_two_term_sfc_verification.json').read_text())
assert uniform['status']=='pass' and uniform['complete_finite_prefix']
assert uniform['effective_uniform_simplicity_cutoff']==10000 and uniform['all_positive_exponents_all_indices_squarefree']
assert uniform['ordinary_exceptional_set_empty'] and uniform['ordinary_A_plus_B_Ym_SFC_all_starts']
assert not any(uniform[k] for k in ('unrestricted_SFC_closed','general_FC_closed','independent_replay','formal_verification','external_review'))
for name,digest in uniform['files'].items():assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==digest,('stale uniform two-term receipt',name)
receipt['all_ordinary_constant_plus_power_families']={
    'claims':['FIRUNIFORMTWOSIMPLE20261117','FIRALLUNIVTWOSFC20261117'],
    'proof':'ALL_UNIVARIATE_TWO_TERM_SFC.tex','checks':uniform['total'],
    'effective_uniform_cutoff':10000,'complete_finite_prefix':True,
    'all_positive_exponents_all_indices_squarefree':True,'ordinary_A_plus_B_Ym_SFC_all_starts':True,
    'ordinary_exceptional_set':'EMPTY, E_m=D_all=0; no MR input in the direct proof',
    'unrestricted_SFC_closed':False,'general_FC_closed':False,
    'independent_replay':False,'formal_verification':False,'external_review':False}
receipt['remaining']=[item for item in receipt['remaining'] if 'other composite exponents' not in item]
receipt['remaining'].append('unrestricted SFC/general FC; genuinely new supports and Gamma parameters beyond ordinary constant-plus-one-variable-power families')
receipt['current_goal_turn_classification']='progress'
receipt['previous_goal_turn_classification']='progress'
receipt['classification_evidence']='The previous goal turn completed ordinary exponent-six simplicity with effective cutoff10000 and a full finite prefix. This goal turn completes a uniform contour argument for6<=m<=sqrt(n), bridges all other exponents by the elementary criterion, and verifies all92935 additional finite pairs. Thus every ordinary A+B*Y^m is closed pointwise, while unrestricted FC/SFC remain open.'
receipt['completion_scope']='Every ordinary constant-plus-one-variable-power family A+B*Y^m, all positive exponents and starting indices; not unrestricted FC/SFC or other Gamma shapes.'

pade=json.loads((HERE/'pade_monomial_separation_verification.json').read_text())
assert pade['status']=='pass' and pade['all_n_a1_through6_proved']
assert pade['unbounded_indexed_strip_family_proved'] and pade['positive_symbolic_certificates_all_n']
assert not any(pade[k] for k in ('all_a_proved','unrestricted_SFC_closed','general_FC_closed','independent_replay','formal_verification','external_review'))
for name,digest in pade['files'].items():assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==digest,('stale Pade separation receipt',name)
receipt['Pade_adjacent_nonconstant_supports']={
    'claims':['FIRPADESEPARATION20261117','FIRADJACENTMONOSFC20261117'],
    'proof':'PADE_ADJACENT_MONOMIALS.tex','checks':pade['total'],
    'ordinary_Xa_A_plus_BX_all_starts_a_values':[1,2,3,4,5,6],
    'indexed_unbounded_family':'a>=2*n-3, terminal n>=2',
    'generic_positive_Schwarz_criterion_proved':True,'necessary_bidiagonal_pencil_proved':True,
    'uniform_all_a_target_proved':False,'unrestricted_SFC_closed':False,'general_FC_closed':False,
    'independent_replay':False,'formal_verification':False,'external_review':False}
receipt['classification_evidence']='Previous goal turn completed every ordinary constant-plus-power family A+B*Y^m. This goal turn proves exact Pade strips, cross determinants and matrix-pencil obstructions, then uses symbolic all-n positive Schwarz certificates to close X^a(A+B*X) at a1..6 and proves the unbounded indexed condition a>=2*n-3. Generic all-a and unrestricted FC/SFC remain open.'
receipt['completion_scope']='Ordinary adjacent nonconstant monomial supports X^a(A+B*X), a1..6 at all starts, plus the indexed strip family; not arbitrary supports or unrestricted FC/SFC.'
receipt['remaining'].append('unbounded a>=7 beyond a>=2n-3 in adjacent nonconstant monomial supports; uniform spectral/shifted criterion unproved')

parabola=json.loads((HERE/'pade_parabola_sfc_verification.json').read_text())
assert parabola['status']=='pass' and parabola['all_a_all_n_adjacent_monomial_SFC']
assert parabola['all_positive_Gamma_shapes'] and parabola['generalized_Pade_coprimality']
assert not any(parabola[k] for k in ('finite_scan_extrapolation','vertical_bound_required','arbitrary_monomial_gaps_closed','unrestricted_SFC_closed','general_FC_closed','independent_replay','formal_verification','external_review'))
for name,digest in parabola['files'].items():assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==digest,('stale uniform parabola receipt',name)
receipt['uniform_adjacent_monomial_SFC']={
    'claims':['FIRPADEPARABOLA20261117','FIRUNIFORMADJACENTSFC20261117'],
    'proof':'UNIFORM_PADE_PARABOLA_SFC.tex','checks':parabola['total'],
    'all_a_all_n_all_positive_Gamma_shapes':True,'finite_prefix_needed':False,
    'stronger_vertical_bound_needed':False,'unrestricted_SFC_closed':False,'general_FC_closed':False,
    'independent_replay':False,'formal_verification':False,'external_review':False}
receipt['Pade_adjacent_nonconstant_supports']['later_uniform_closure']='FIRUNIFORMADJACENTSFC20261117; every a/start/positive shape'
receipt['classification_evidence']='Progress: earlier finite-a Pade strip/matrix work is now extended uniformly by complementary parabola inclusions. A self-contained continued-fraction invariant proves every positive integer parameter shift and every positive Gamma shape. This closes X^a(A+B*X) at all starts/all a, without assuming the stronger vertical bound or scanning more indices. Arbitrary gaps/supports and unrestricted FC/SFC remain open.'
receipt['completion_scope']='ALL X^a(A+B*X), a>=0, all starts and all positive Gamma shapes; plus independent-variable monomial factors. Not arbitrary gaps or unrestricted FC/SFC.'
receipt['remaining']=[item for item in receipt['remaining'] if 'unbounded a>=7' not in item]
receipt['remaining'].append('X^a(A+B*X^d), a>=1,d>=2; genuinely new/multivariate supports beyond completed families')

dyadic=json.loads((HERE/'dyadic_binomial_sfc_verification.json').read_text())
assert dyadic['status']=='pass' and dyadic['dyadic_equal_lower_exponent_and_gap_all_indices']
assert dyadic['integer_Gamma_shapes_through_2q'] and dyadic['arbitrary_rational_unit_sequence_transfer']
assert dyadic['same_variable_scaling_at_consecutive_orders'] and dyadic['exact_resultant_valuation_proved']
assert dyadic['lower_one_multiple_root_equivalence_proved'] and dyadic['cubic_quarter_gap_square_degree_bound_proved']
assert dyadic['conditional_cubic_prime_confinement_proved'] and dyadic['cubic_residual_transfer_proved'] and not dyadic['cubic_all_index_SFC_closed']
assert not any(dyadic[k] for k in ('finite_parameter_extrapolation','other_lower_exponents_or_gaps_proved','all_positive_Gamma_shapes','general_FC_closed','unrestricted_SFC_closed','independent_replay','formal_verification','external_review'))
for name,digest in dyadic['files'].items():assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==digest,('stale dyadic binomial receipt',name)
receipt['dyadic_binomial_SFC']={
    'claims':['FIRBESSELUNITTRANSFER20261117','FIRDYADICBINSFC20261117','FIRLOWERONESIMPLE20261117','FIRCUBICPRIMECONF20261117'],
    'proof':'DYADIC_BINOMIAL_SFC.tex','checks':dyadic['total'],
    'all_indices_q2r_lower_exponent_equal_gap':True,'integer_Gamma_shapes1_to_2q':True,
    'arbitrary_rational_unit_sequence_transfer':True,'exact_resultant_valuation':True,
    'lower_one_equivalence_proved':True,'conditional_cubic_prime_confinement_proved':True,
    'cubic_all_index_SFC_closed':False,'other_lower_exponents_or_gaps_proved':False,'all_positive_Gamma_shapes':False,
    'unrestricted_SFC_closed':False,'general_FC_closed':False,
    'independent_replay':False,'formal_verification':False,'external_review':False}
receipt['classification_evidence']='Previous goal turn proved the uniform adjacent-support family. Current turn proves the new all-index nonadjacent dyadic diagonal, gives an arbitrary-unit Bessel transfer and an exact resultant valuation, and derives the exact lower-one multiple-root bridge plus a conditional quarter-prime-gap confinement/residual transfer for the first remaining cubic. Uniform coverage and simplicity are still open. Universal algebraic proof, no parameter extrapolation. General FC/SFC remain open.'
receipt['completion_scope']='X^q(A+B*X^q),q=2^r,all starts,integer Gamma shapes1..2q; independent-variable monomial factors allowed. Not all lower exponents/gaps or unrestricted FC/SFC.'
receipt['remaining'].append('First remaining ordinary gap X*(A+B*X^2), then other pairs outside completed dyadic diagonal; higher-precision/unit or spectral information required')

tilted=json.loads((HERE/'tilted_cubic_frontier_verification.json').read_text())
assert tilted['status']=='pass' and tilted['uniform_shifted_polygon_and_residuals_proved']
assert tilted['all_prime_power_indices_irreducible'] and tilted['ten_uniform_shifted_index_families_proved']
assert tilted['exact_analytic_fold_data_proved']
assert not any(tilted[k] for k in ('finite_extrapolation_used','uniform_prime_coverage_proved','effective_simplicity_cutoff_proved','all_index_cubic_SFC_closed','general_FC_closed','unrestricted_SFC_closed','independent_replay','formal_verification','external_review'))
for name,digest in tilted['files'].items():assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==digest,('stale tilted cubic frontier receipt',name)
receipt['tilted_cubic_prime_power_and_fold']={
    'claims':['FIRTILTEDSHIFT20261117','FIRTILTEDPRIMEPOWER20261117','FIRCUBICFOLDS20261117'],
    'proof':'TILTED_CUBIC_PRIME_POWER_AND_FOLD.tex','checks':tilted['total'],
    'all_prime_power_INDICES_irreducible':True,'ten_uniform_shifted_index_families':True,
    'exact_fold_data_only':True,'uniform_prime_coverage':False,'effective_simplicity_cutoff':False,
    'cubic_all_starting_indices_closed':False,'general_FC_closed':False,'unrestricted_SFC_closed':False,
    'independent_replay':False,'formal_verification':False,'external_review':False}
receipt['classification_evidence']='Previous goal turn completed the dyadic-gap family and exact lower-one simplicity/next-prime reductions. This turn proves uniform shifted polygons/residuals, irreducibility at EVERY prime-power tilted-cubic index, ten complete shifted index families, and exact saddle-fold geometry. Retained exact counterprofiles disprove convenient fixed-prime and simple exclusion/real-rootedness shortcuts. Uniform cubic coverage/simplicity and unrestricted FC/SFC remain open.'
receipt['completion_scope']='All prime-power terminal indices and ten shifted terminal families for the ordinary cubic, plus exact shift and fold data. Not all cubic starting indices or unrestricted FC/SFC.'
receipt['remaining'].append('Uniform simplicity/prime-residual coverage of T_(n,2) outside completed indexed families, or justified all-root fold analysis with an effective cutoff and complete prefix')

cubic=json.loads((HERE/'cubic_adjoint_stability_verification.json').read_text())
assert cubic['status']=='pass' and cubic['all_index_cubic_SFC_closed'] and cubic['all_positive_Gamma_shapes']
assert not any(cubic[k] for k in ('finite_prefix_required','finite_scan_extrapolation','all_gaps_or_lower_exponents_proved','general_FC_closed','unrestricted_SFC_closed','independent_replay','formal_verification','external_review','novelty_claim'))
for name,digest in cubic['files'].items():assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==digest,('stale cubic adjoint receipt',name)
receipt['whole_cubic_SFC_by_adjoint_stability']={
    'claims':['FIRCUBICGAMMASIMPLE20261117','FIRCUBICALLSFC20261117'],
    'proof':'CUBIC_ADJOINT_STABILITY_SFC.tex','checks':cubic['total'],
    'all_cubic_starting_indices_closed':True,'all_positive_Gamma_shapes':True,
    'simplicity_parameter_range':'ALL n>=1,gamma>0,gamma>=n-2',
    'finite_prefix_needed':False,'classical_input':'GWS half-plane polarization, Wagner Theorem4.1',
    'unrestricted_SFC_closed':False,'general_FC_closed':False,
    'independent_replay':False,'formal_verification':False,'external_review':False}
for key,flag in [('dyadic_binomial_SFC','cubic_all_index_SFC_closed'),('tilted_cubic_prime_power_and_fold','cubic_all_starting_indices_closed')]:
    earlier=receipt[key]
    earlier[flag+'_in_this_earlier_stage']=earlier.pop(flag)
    earlier['later_whole_cubic_closure']='FIRCUBICALLSFC20261117; every index and every positive Gamma shape'
receipt['classification_evidence']='Progress: the preceding turn established prime-power index and shifted cubic classes plus exact folds. This turn proves whole-family cubic SFC by an explicit adjoint and Gamma-convolution half-plane separation, for every starting index and every positive Gamma shape. The proof is uniform, with no cutoff or finite-prefix extrapolation. General FC/unrestricted SFC remain open.'
receipt['completion_scope']='Every X*(A+B*X^2), all starts, all positive Gamma shapes; independent-variable monomial factors allowed. Not arbitrary supports or unrestricted FC/SFC.'
receipt['remaining']=[item for item in receipt['remaining'] if not ('First remaining ordinary gap' in item or 'Uniform simplicity/prime-residual coverage of T_(n,2)' in item)]
receipt['remaining'].append('Genuinely new lower exponents/gaps and arbitrary supports; a valid reduction or new argument toward general FC/unrestricted SFC')

highergap=json.loads((HERE/'higher_gap_adjoint_rank_verification.json').read_text())
assert highergap['status']=='pass' and highergap['uniform_higher_gap_adjoint_proved'] and highergap['exact_repeated_root_rank_equivalence_proved']
assert highergap['full_Wronskian_reconstruction_proved'] and highergap['all_positive_shape_lower_one_bridge_proved']
assert highergap['gap_three_closed_sector_sufficiency_proved'] and highergap['Hurwitz_shortcut_refuted']
assert highergap['all_gap_initial_four_orders_all_positive_shapes_proved']
assert not any(highergap[k] for k in ('gap_three_sector_zero_freeness_proved','higher_gap_all_index_SFC_closed','general_FC_closed','unrestricted_SFC_closed','finite_extrapolation','independent_replay','formal_verification','external_review','novelty_claim'))
for name,digest in highergap['files'].items():assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==digest,('stale higher-gap adjoint receipt',name)
receipt['all_gap_adjoint_and_rank_frontier']={
    'claims':['FIRGENERALGAPADJOINT20261117','FIRGENERALGAPRANK20261117','FIRALLGAPINITIAL20261117'],
    'proof':'HIGHER_GAP_ADJOINT_RANK.tex','checks':highergap['total'],
    'all_gap_adjoint_and_rank_equivalence':True,'full_Wronskian_reconstruction':True,'all_gap_initial_four_orders_all_positive_shapes':True,
    'all_positive_shape_lower_one_bridge':True,'gap_three_closed_sector_suffices':True,
    'uniform_sector_or_interlacing_proved':False,'full_Hurwitz_shortcut_refuted_exactly':True,
    'higher_gap_all_index_SFC_closed':False,'unrestricted_SFC_closed':False,'general_FC_closed':False,
    'independent_replay':False,'formal_verification':False,'external_review':False}
receipt['classification_evidence']='Progress: preceding turn closed the whole cubic support. Current turn proves explicit all-gap polynomial adjoints, their exact repeated-root rank equivalence and Wronskian determinant reconstruction, and extends the lower-one bridge to every positive Gamma shape. A uniform positive-third-central-moment argument also proves the initial orders n1..4 for ALL gaps and positive Gamma shapes. The natural stronger gap3 Hurwitz shortcut is refuted by an exact determinant; a precise cyclic-section/closed-sector uniform gate remains open. Finite interlacing controls are not extrapolated.'
receipt['completion_scope']='ALL gaps d>=2: explicit adjoints and exact algebraic rank characterization; ALL d>=1/shapes: initial orders1..4 simple. No uniform higher-gap all-index simplicity theorem. Whole cubic remains completed; unrestricted FC/SFC remain open.'
receipt['remaining'].append('Gap3 uniform cyclic-section coprimality/interlacing or CLOSED pi/3 sector exclusion; full Hurwitz stability is false')

wr=json.loads((HERE/'gap_three_wronskian_verification.json').read_text())
assert wr['status']=='pass' and wr['uniform_Wronskian_shape_shift_and_recurrence_proved']
assert wr['bivariate_coefficient_positivity_proved'] and wr['no_common_real_section_root_proved']
assert wr['one_real_rooted_section_suffices_proved'] and wr['fifth_order_all_positive_Gamma_shapes_proved']
assert wr['separate_input_interlacing_induction_refuted']
assert not any(wr[k] for k in ('uniform_section_real_rootedness_proved','all_index_gap_three_SFC_closed','general_FC_closed','unrestricted_SFC_closed','finite_extrapolation','independent_replay','formal_verification','external_review','novelty_claim'))
for name,digest in wr['files'].items():assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==digest,('stale gap-three Wronskian receipt',name)
receipt['gap_three_Wronskian_recurrence_and_one_real_section']={
    'claims':['FIRGAPTHREEREC20261117','FIRGAPTHREEREAL20261117','FIRGAPTHREEFIFTH20261117'],
    'proof':'GAP_THREE_WRONSKIAN_RECURRENCE.tex','checks':wr['total'],
    'uniform_recurrence_and_bivariate_positivity':True,'no_common_real_section_roots':True,
    'one_real_rooted_section_suffices':True,'all_positive_shape_initial_block_through5':True,
    'uniform_section_existence_proved':False,'all_index_gap_three_SFC_closed':False,
    'unrestricted_SFC_closed':False,'general_FC_closed':False,
    'independent_replay':False,'formal_verification':False,'external_review':False}
receipt['classification_evidence']='Progress: preceding turn derived all-gap adjoint/rank structure and initial orders1..4. Current turn proves a uniform gap3 Wronskian/Gamma-shift recurrence and coefficient positivity, weakens the uniform gate to one real-rooted section by excluding all common real section roots, and proves terminal5 for ALL positive Gamma shapes by a symbolic positive discriminant. An exact separate-input full-interlacing induction shortcut is refuted. No finite profile extrapolation or unrestricted FC/SFC closure.'
receipt['completion_scope']='Uniform gap3 Wronskian identities, the one-real-section criterion, and initial terminal orders through5 at ALL positive Gamma shapes. Whole cubic remains closed; general FC/unrestricted SFC remain open.'
receipt['remaining']=[item for item in receipt['remaining'] if not item.startswith('Gap3 uniform cyclic-section coprimality/interlacing')]
receipt['remaining'].append('Uniform existence of ONE real-rooted nonzero section or real linear combination at gap3; use actual Gamma-shift compatibility or another valid mechanism. General FC/unrestricted SFC remain open.')

pc=json.loads((HERE/'gap_three_polar_continuant_verification.json').read_text())
assert pc['status']=='pass' and pc['uniform_real_rooted_polar_companion_proved'] and pc['companion_roots_simple_below_minus1_proved']
assert pc['nonzero_pivot_backward_coprimality_criterion_proved'] and pc['zero_pivot_fail_closed']
assert pc['generic_inverse_transfer_to_real_rooted_individual_section_refuted'] and pc['native_prime_squarefree_residual_shortcut_refuted']
assert not any(pc[k] for k in ('positive_weights_required_for_coprimality','uniform_pivot_nonvanishing_proved','uniform_section_existence_proved','all_index_gap_three_SFC_closed','general_FC_closed','unrestricted_SFC_closed','finite_extrapolation','independent_replay','formal_verification','external_review','novelty_claim'))
for name,digest in pc['files'].items():assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==digest,('stale polar/continuant receipt',name)
receipt['gap_three_polar_companion_and_nonzero_pivots']={
    'claims':['FIRGAPTHREEPOLAR20261117','FIRGAPTHREEPIVOT20261117'],
    'proof':'GAP_THREE_POLAR_AND_CONTINUANT.tex','checks':pc['total'],
    'uniform_simple_companion_roots_below_minus1':True,
    'generic_inverse_transfer_to_individual_section_refuted':True,
    'nonzero_pivot_backward_coprimality_criterion':True,'positive_pivots_required':False,'zero_pivots_inconclusive':True,
    'uniform_actual_pivot_nonvanishing':False,'uniform_section_existence':False,
    'all_index_gap_three_SFC_closed':False,'unrestricted_SFC_closed':False,'general_FC_closed':False,
    'independent_replay':False,'formal_verification':False,'external_review':False}
receipt['classification_evidence']='Progress: previous turn proved Wronskian recurrences, the one-real-section criterion and the fifth-order initial block. Current turn proves a uniform simple real-rooted polar companion and exact finite-difference identities, refutes a generic inverse transfer to individual sections, and proves a nonzero-pivot backward certificate that does not require positivity. Exact source-hypothesis and native-prime counterprofiles delimit shortcuts; no uniform actual transfer/pivot theorem or general FC/SFC closure is inferred.'
receipt['completion_scope']='Uniform polar companion at gamma>n-2 and general backward nonzero-pivot coprimality criterion. Finite certificates remain scoped; prior completed support families remain closed. No all-index gap3 or unrestricted FC/SFC closure.'
receipt['remaining'].append('Exploit the actual arithmetic-progression polar factors for a valid section transfer, or prove actual-family backward pivots nonzero; generic companion hyperbolicity alone is insufficient for individual sections.')

if shutil.which('pdflatex'):
    with tempfile.TemporaryDirectory(prefix='factorial-tex-') as build:
        for filename in documents:
            command = ['pdflatex','-draftmode','-halt-on-error','-interaction=batchmode','-output-directory',build,filename]
            for run in range(2):
                result = subprocess.run(command, cwd=HERE, capture_output=True, text=True, timeout=30)
                if result.returncode:
                    log = (Path(build)/Path(filename).with_suffix('.log')).read_text(errors='replace')
                    raise AssertionError(filename+'\n'+log[-5000:])
            log = (Path(build)/Path(filename).with_suffix('.log')).read_text(errors='replace')
            assert 'undefined references' not in log
            receipt['latex_syntax'].append({'file':filename,'status':'pass','mode':'draftmode: no PDF generated'})

for script in ('render_status.py','research.py'):
    command = [sys.executable,str(ROOT/'scripts'/script)]
    if script == 'research.py':
        command.append('render')
    result = subprocess.run(command, cwd=ROOT.parent, capture_output=True, text=True, timeout=30)
    receipt['global_navigation'].append({'script':script,'returncode':result.returncode,
                                         'diagnostic':(result.stderr or result.stdout).strip()[-1600:]})
# Global validation above is reported unchanged. Render only this packet's
# affected navigation through the existing renderer functions; no missing
# artifacts are fabricated and no global validation success is asserted.
import research_work
views = research.render_outputs(data['entries'], lessons, work_data=research_work.load_work(ROOT))
affected_views = ('index/core.md','index/README.md','index/documents.tsv','knowledge/work/geometry.md')
receipt['scoped_navigation'] = []
for name in affected_views:
    path = ROOT/name
    expected = views[path]
    changed = path.read_text() != expected
    if changed:
        path.write_text(expected)
    assert path.read_text() == expected
    receipt['scoped_navigation'].append({'file':name,'status':'pass','regenerated':changed})
expected_status=render_status.render(data)
if (ROOT/'STATUS.md').read_text()!=expected_status:
    (ROOT/'STATUS.md').write_text(expected_status)
receipt['status_view_current'] = (ROOT/'STATUS.md').read_text()==expected_status
(HERE/'integration.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
