"""Expected independent-cloud yield over the exact V4 residue/box geometry.

Gamma/Poisson updates are scheduling assumptions. Neither their zero scores
nor exhausted numerical exposure are rational-point exclusions.
"""
from math import exp, expm1

import numpy as np
from cancellation_scheduler import AnchorExposure, HEIGHTS


def bucket(index):return 0 if index<4 else 1 if index<16 else 2


def capped_poisson_mean(mean, limit):
    """E[min(Poisson(mean), limit)], with stable small-intensity evaluation."""
    if mean<0 or limit<1:raise ValueError('nonnegative mean and positive limit required')
    if mean==0:return 0.
    probability=exp(-mean);tail=-expm1(-mean);value=tail
    for k in range(1,limit):
        probability*=mean/k
        tail=max(0.,tail-probability)
        value+=tail
    return min(float(limit),value)


class CloudExposure(AnchorExposure):
    def __init__(self,prepared,fit,index):
        super().__init__(prepared,fit,index)
        self.shape=fit['cloud_prior_exposure']*fit['anchor_direction_rates'][bucket(index)]
        self.rate=float(fit['cloud_prior_exposure'])

    def options(self,cost,remaining):
        increments=np.maximum(self.exposures-self.covered,0) @ self.weights
        intensity=self.shape/self.rate
        for job,((model,height),mass) in enumerate(zip(self.jobs,increments)):
            if job in self.attempted or mass<=1e-12:continue
            mean=intensity*float(mass)
            credit=capped_poisson_mean(mean,remaining)
            price=cost.estimate(height,-expm1(-mean))
            yield credit/price,self.index,job,model,height,mean,credit,price

    def observe_cloud(self,job,complete,gains):
        if gains<0 or (not complete and gains):
            raise ArithmeticError('cloud update needs a completed counted box')
        before=float(self.covered @ self.weights)
        super().observe(job,complete)
        if complete:
            self.rate+=max(0.,float(self.covered @ self.weights)-before)
            self.shape+=gains


def choose(anchors,cost,fit,next_index,count,remaining):
    options=[row for anchor in anchors.values() for row in anchor.options(cost,remaining)]
    best=max(options,key=lambda r:(r[0],-r[1],-r[2])) if options else None
    fraction=sum(h<=HEIGHTS[0] for h in fit['radial_heights'])/len(fit['radial_heights'])
    mean=fit['anchor_direction_rates'][bucket(next_index)]*fraction
    score=capped_poisson_mean(mean,remaining)/(cost.estimate(HEIGHTS[0],-expm1(-mean))+.008)
    if next_index<count and (best is None or score>best[0]):
        return {'kind':'prepare','index':next_index,'score':score}
    if best is None:return {'kind':'exhausted'}
    score,index,job,model,height,mean,credit,price=best
    return {'kind':'search','index':index,'job':job,'model':model,'height':height,
        'score':score,'mean_new_directions':mean,'expected_capped_gain':credit,'predicted_cpu':price}
