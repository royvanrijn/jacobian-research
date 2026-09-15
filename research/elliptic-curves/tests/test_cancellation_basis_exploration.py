"""Arithmetic exposure and timeout regressions for the unfitted successor."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'cas'))
from cancellation_basis_exploration import Exposure, in_box, apply
from cancellation_cloud_policy import CloudExposure


def test_zero_fitted_mass_does_not_remove_unsearched_coordinates():
    base = {'matrix': [1, 0, 0, 1]}
    neighbour = {'matrix': [3, 0, 0, 1]}
    fit = {'anchor_priors': [.2]*3, 'radial_heights': [8000, 32000, 125000],
           'anchor_direction_rates': [.2]*3, 'cloud_prior_exposure': 4}
    # Base H125000 exhausts the empirical CDF, including this neighbour.
    cloud = CloudExposure({'weights': [1.], 'ratios': [[1., 1/3]],
                           'models': [base, neighbour]}, fit, 0)
    cloud.observe_cloud(2, True, 0)
    assert (cloud.exposures <= cloud.covered).all()
    exposure = Exposure()
    exposure.observe([2, 3], base, 125000, True)
    assert not exposure.seen([2, 3], neighbour, 125000)
    witness = exposure.inspect([2, 3], neighbour, base, 125000)
    assert witness['old_height'] > 125000 >= witness['new_height']
    raw = apply(neighbour['matrix'], witness['new_coordinate'])
    assert not in_box(raw, base['matrix'], 125000)


def test_only_completed_exact_boxes_are_deduplicated():
    e = Exposure(); base = {'matrix': [1, 0, 0, 1]}
    equivalent = {'matrix': [0, -2, 2, 0]}
    e.observe([2, 3], base, 125000, False)
    assert not e.seen([2, 3], equivalent, 125000)
    e.observe([2, 3], base, 125000, True)
    assert e.seen([2, 3], equivalent, 125000)
    assert not e.seen([2, -3], equivalent, 125000)
    assert not e.seen([2, 3], equivalent, 125001)


def test_missing_boundary_witness_is_not_an_exclusion():
    e = Exposure(); base = {'matrix': [1, 0, 0, 1]}
    assert e.inspect([2, 3], base, base, 125000) is None
    assert not e.seen([2, 3], base, 125000)


def test_completed_union_includes_projective_infinity():
    assert in_box((1, 0), [1, 0, 0, 1], 1)
    assert in_box((0, 1), [0, 1, 1, 0], 1)
