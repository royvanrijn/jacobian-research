"""Exploratory six-cell standard-chain gaps in the positive cone.

Ritz integrals use scrambled Sobol, and Nyström uses a second Sobol
design with full folding. Neither calculation is directed.
"""
import argparse
from itertools import permutations
import json
import numpy as np
from scipy.integrate import cumulative_simpson
from scipy.linalg import eigh
from scipy.special import ndtr, ndtri, roots_hermite
from scipy.stats import qmc

def facets(samples):
    grid = np.linspace(-8.0, 8.0, 1 + 2**16)
    density = np.exp(-grid**2) / np.sqrt(np.pi) * ndtr(grid)**4
    cumulative = cumulative_simpson(density, x=grid, initial=0)
    mass = cumulative[-1]
    u = np.interp(mass * samples[:, 0], cumulative, grid)
    g = u[:, None] - ndtri(ndtr(u[:, None]) * samples[:, 1:5])
    return mass, g

def row_mass(g, rho, order=32):
    z, weights = roots_hermite(order)
    weights /= np.sqrt(np.pi)
    s = rho / np.sqrt(1-rho*rho) * g
    def integrate(v):
        out = []
        for start in range(0, len(g), 2048):
            m = v[start:start+2048]
            out.append(np.prod(ndtr(m[:, None, :] + z[None, :, None]), axis=2) @ weights)
        return np.concatenate(out)
    result = integrate(s)
    for j in range(4):
        adjacent = s - s[:, j:j+1]/2
        adjacent[:, j] -= s[:, j]
        result += 0.5 * np.exp(-s[:, j]**2/4) * integrate(adjacent)
    return result

def basis(g, degree):
    z = g/(1+g)
    power = [(z**d).sum(axis=1) for d in range(1, 5)]
    exponent = []
    for d in range(1, degree+1):
        for a in range(d+1):
            for b in range((d-a)//2+1):
                for c in range((d-a-2*b)//3+1):
                    rem = d-a-2*b-3*c
                    if rem % 4 == 0:
                        exponent.append((a,b,c,rem//4))
    return np.column_stack([np.prod([p**e for p,e in zip(power, abcd)],axis=0) for abcd in exponent])

def ritz(points, gx, x, rho):
    b = row_mass(gx, rho)
    tau = b.mean()
    joint_mass = 0.0
    first = np.zeros(x.shape[1])
    second = np.zeros((x.shape[1], x.shape[1]))
    e = np.zeros_like(second)
    t = rho / np.sqrt(1-rho*rho) * gx
    conditional = ndtri(points[:, 5:6]) / np.sqrt(2)
    for j in (None, 0):
        if j is None:
            z, factor = t, 1
        else:
            z = t - t[:, j:j+1]/2
            z[:, j] -= t[:, j]
            factor = 2 * np.exp(-t[:, j]**2/4)
        z = z + conditional
        limits = ndtr(z)
        gy = np.sqrt(1-rho*rho) * (z - ndtri(limits * points[:, 6:10]))
        weight = factor * np.prod(limits, axis=1)
        y = basis(gy, 5)
        diff = x-y
        e += (diff.T * weight) @ diff / (2*len(b))
        joint_mass += weight.mean()
        first += weight @ (x+y) / (2*len(b))
        second += ((x.T*weight) @ x +(y.T*weight) @ y) / (2*len(b))
    var = second-np.outer(first,first)/joint_mass
    dims = (6, x.shape[1])
    gaps = [float(eigh(e[:n,:n],var[:n,:n],subset_by_index=(0,0),eigvals_only=True)[0]) for n in dims]
    correction = tau * np.mean(1/b)
    return dict(E_D=float(tau), E_inverse_D=float(np.mean(1/b)), correction=float(correction),
                degree_3_gap=gaps[0], degree_5_gap=gaps[1],
                degree_5_ratio=gaps[1]/((1-rho)*correction),
                joint_row_relative_error=float(joint_mass/tau-1))

def nystrom(sample, rho, dilation):
    p, gaps = facets(sample)
    gaps = dilation * gaps
    a = gaps.sum(axis=1)/6
    x = np.column_stack((a,a,a[:,None]-np.sort(gaps, axis=1)))
    s2 = 1-rho*rho
    prefactor = np.exp(-(rho*rho/(2*s2))*np.sum(x*x,axis=1))
    b = np.zeros((len(x),len(x)))
    for perm in permutations(range(1,6)):
        dot = (x[:,1:] @ x[:,perm].T + np.outer(x[:,0],x[:,0])) * (rho/s2)
        b += (2 if perm[0] == 1 else 1) * np.exp(dot)
    b *= (prefactor[:,None]*prefactor[None,:])/(s2*s2*48)
    b = (b+b.T)/2
    weights = dilation**4*np.exp(-(1-1/dilation**2)*np.sum(x*x,axis=1)/2)
    weights /= weights.sum()
    row = (b @ weights) * p
    scale = np.sqrt(weights/row)
    b *= p * scale[:,None]*scale[None,:]
    values=eigh(b, subset_by_index=(len(x)-3,len(x)-1),eigvals_only=True)
    actual = row_mass(gaps, rho,order=48)
    def weighted(v): return float(weights @ v)
    return dict(gap=float(1-values[-2]), third_eigenvalue=float(values[-3]),
                correction=weighted(row)*weighted(1/row),
                analytic_correction=weighted(actual)*weighted(1/actual),
                analytic_ratio=float(1-values[-2])/((1-rho)*weighted(actual)*weighted(1/actual)),
                row_relative_rms=float(np.sqrt(weighted((row/actual-1)**2))),
                max_diagonal=float(np.diag(b).max()))

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=('ritz','nystrom'),required=True)
    parser.add_argument('--power', type=int, default=15)
    parser.add_argument('--seed',type=int,default=221)
    parser.add_argument('--rho',nargs='+',type=float,default=[.20,.35,.50,.65,.80,.90])
    parser.add_argument('--dilation',type=float,default=1.0)
    args=parser.parse_args()
    points=qmc.Sobol(10 if args.mode=='ritz' else 5,scramble=True,seed=args.seed).random_base2(args.power)
    if args.mode=='ritz':
        mass,g=facets(points)
        features=basis(g,5)
    for rho in args.rho:
        result=ritz(points,g,features,rho) if args.mode=='ritz' else nystrom(points,rho,args.dilation)
        print(json.dumps(dict(rho=rho,mode=args.mode,power=args.power,seed=args.seed,dilation=args.dilation,**result)),flush=True)

if __name__=='__main__':
    main()
