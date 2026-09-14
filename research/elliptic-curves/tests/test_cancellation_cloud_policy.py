"""Statistical-state invariants, including the timeout non-exclusion rule."""
import math
from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'cas'))
from cancellation_cloud_policy import capped_poisson_mean,CloudExposure


def exposure():
    prepared={'weights':[.4,.6],'ratios':[[1.,2.],[2.,1.]],'models':[{},{}]}
    fit={'anchor_priors':[.2,.1,.05],'anchor_direction_rates':[.2,.1,.05],
        'cloud_prior_exposure':4.,'radial_heights':[1000,8000,32000,125000]}
    return CloudExposure(prepared,fit,0)


def test_poisson_remaining_target_changes_expected_credit():
    for mean in (.001,.2,1.,8.):
        assert capped_poisson_mean(mean,1)==pytest.approx(1-math.exp(-mean))
        assert capped_poisson_mean(mean,2)==pytest.approx(2-(2+mean)*math.exp(-mean))
        assert capped_poisson_mean(mean,1)<=capped_poisson_mean(mean,2)<=min(2.,mean)
    assert capped_poisson_mean(0.,2)==0.
    assert capped_poisson_mean(1000.,2)==2.


def test_positive_cloud_updates_remaining_intensity():
    e=exposure();before=e.shape/e.rate
    e.observe_cloud(0,True,2)
    assert e.shape/e.rate>before
    assert 0<float(e.covered @ e.weights)<1
    assert e.shape==pytest.approx(2.8)


def test_timeout_removes_no_exposure_and_no_intensity():
    e=exposure();before=e.shape,e.rate,e.covered.copy()
    e.observe_cloud(0,False,0)
    assert (e.shape,e.rate)==before[:2]
    assert np.array_equal(e.covered,before[2])
    assert 0 in e.attempted
    with pytest.raises(ArithmeticError):e.observe_cloud(1,False,1)


def test_completed_overlap_cannot_count_as_new_exposure():
    e=exposure();e.observe_cloud(2,True,0);after=e.rate
    e.observe_cloud(0,True,0)
    assert e.rate==after
    assert e.shape/e.rate<.2
