#!/usr/bin/env sage-python
"""Projective all-parameter gate for the frozen67-by25 carrier families.

The source P2 parameter map uses a homogenized chord sextic. The target is
the P1-by-P1 product of a complete branch-quartic pencil and a squared linear
factor. Base points and bad coefficient reductions are deferred, never used
as exclusions. This is a complete finite-projective-image test, not a point
search or a rational-point exhaustion of any halving curve.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import resource
import time

from sage.all import PolynomialRing, QQ, matrix

ROOT=Path(__file__).resolve().parents[2]
NODES=ROOT/"artifacts/generated-results/elkies-k3-r17-one-node-correlated-v1/input.json"
PENCILS=ROOT/"artifacts/generated-results/elkies-k3-r17-correlated-genus-one-v1/pencils.json"
DEFAULT=ROOT/"artifacts/generated-results/elkies-k3-r17-one-node-complete-pairs-v1"
PRIMES=[101,103,107,109,113]


def write_new(path,data):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open("x") as f:json.dump(data,f,indent=2,sort_keys=True);f.write("\n")


def prepare(path):
    source=json.loads(NODES.read_text());targets=json.loads(PENCILS.read_text())
    R=PolynomialRing(QQ,'t');t=R.gen();S=PolynomialRing(R,'L');L=S.gen()
    A=R(source['generic_source']['A']);families=[]
    for trace in source['traces']:
        h,Nx,Ny,M0=[R(trace[k]) for k in ('h','Nx','Ny','M0')]
        M=M0+L*h*h
        numerator=M**4-6*M*M*Nx-8*M*Ny-3*Nx*Nx-4*A*h**4
        channels=[]
        for f in numerator:
            q,rem=f.quo_rem(h**6);assert not rem
            channels.append([str(x) for x in q.list()])
        assert len(channels)==5
        families.append({'index':trace['index'],'word':trace['word'],'channels':channels})
    packet={'schema':'r17-one-node-complete-pairs-input-v1',
            'source_input_sha256':sha256(NODES.read_bytes()).hexdigest(),
            'target_pencils_sha256':sha256(PENCILS.read_bytes()).hexdigest(),
            'source_families':families,'target_matrices':[p['branch_matrix'] for p in targets],
            'limits':{'cpu_seconds':120,'address_space_gib':4,'primes':PRIMES,'pairs':1675},
            'selection':'Exactly the frozen67 norm-six trace nets and25 norm-eight smooth pencils. No new word, desired point, node position or specialization is selected.',
            'identity':'[sum_j c^(4-j)*(a+b*t)^j*q_j(t)] = [ell_r(t)^2 * B_T*Veronese_4(lambda)] in P6, with (a:b:c) in P2 and r,lambda in P1.'}
    write_new(path/'input.json',packet)
    print(json.dumps({'prepared_source_families':67,'target_families':25,'pairs':1675}),flush=True)


def mod_coefficients(data,p):
    out=[]
    for row in data:
        values=[QQ(x) for x in row]
        if any(v.denominator()%p==0 for v in values):return None
        out.append([int(v.numerator()%p)*pow(int(v.denominator()%p),-1,p)%p for v in values])
    return out


def normal(v,p):
    for c in v:
        if c%p:
            inv=pow(c%p,-1,p);return tuple(x*inv%p for x in v)
    return None


def multiply(f,g,p):
    out=[0]*(len(f)+len(g)-1)
    for i,a in enumerate(f):
        for j,b in enumerate(g):out[i+j]=(out[i+j]+a*b)%p
    return out


def target_images(B,p):
    quartics=[]
    for lam in range(p):
        vv=[pow(lam,j,p) for j in range(5)]
        quartics.append([sum(a*b for a,b in zip(row,vv))%p for row in B])
    quartics.append([row[4] for row in B])
    for q in quartics:
        for r in range(p):yield normal(multiply(q,[r*r%p,-2*r%p,1],p),p)
        yield normal(q+[0,0],p)  # r=infinity, the squared homogeneous z factor.


def source_images(channels,p):
    # All affine c=1 parameters, followed by every (a:b:0).
    for a in range(p):
        for b in range(p):
            powers=[[1]]
            for j in range(4):powers.append(multiply(powers[-1],[a,b],p))
            out=[0]*7
            for f,L in zip(channels,powers):
                product=multiply(f,L,p)
                for i,x in enumerate(product):out[i]=(out[i]+x)%p
            yield normal(out,p)
    h2=channels[4]
    for a in range(p):
        L=[1]
        for j in range(4):L=multiply(L,[a,1],p)
        yield normal(multiply(h2,L,p),p)
    yield normal(h2+[0]*(7-len(h2)),p)  # (a:b:c)=(1:0:0).


def run(path):
    resource.setrlimit(resource.RLIMIT_CPU,(120,125));resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3))
    started=time.process_time();packet=json.loads((path/'input.json').read_text())
    remaining={(i,j) for i in range(67) for j in range(25)};certificates=[]
    for p in packet['limits']['primes']:
        before=len(remaining);target_lookup={};good_targets=[]
        for j,rows in enumerate(packet['target_matrices']):
            B=mod_coefficients(rows,p)
            if B is None or not int(matrix(QQ,B).det())%p:continue
            good_targets.append(j)
            for q in target_images(B,p):
                assert q is not None
                target_lookup[q]=target_lookup.get(q,0)|(1<<j)
        excluded=[];source_audits=[]
        for i,row in enumerate(packet['source_families']):
            C=mod_coefficients(row['channels'],p)
            if C is None:
                source_audits.append({'source':i,'status':'DEFERRED_NONINTEGRAL_COEFFICIENTS'});continue
            seen=0;base_points=0;count=0
            for q in source_images(C,p):
                count+=1
                if q is None:base_points+=1
                else:seen|=target_lookup.get(q,0)
            assert count==p*p+p+1
            audit={'source':i,'projective_parameters':count,'base_points':base_points,
                   'intersecting_targets':[j for j in good_targets if (seen>>j)&1]}
            source_audits.append(audit)
            if not base_points:
                for j in good_targets:
                    if (i,j) in remaining and not (seen>>j)&1:
                        remaining.remove((i,j));excluded.append([i,j])
        certificate={'prime':p,'good_targets':good_targets,'source_audits':source_audits,
                     'excluded_pairs':excluded,'remaining_pairs':sorted(map(list,remaining)),
                     'cpu_seconds':time.process_time()-started}
        write_new(path/f'prime-{p}.json',certificate);certificates.append(certificate)
        print(json.dumps({'prime':p,'new_exclusions':before-len(remaining),'remaining':len(remaining),'cpu_seconds':time.process_time()-started}),flush=True)
        if not remaining:break
    result={'schema':'r17-one-node-complete-pair-result-v1',
            'input_sha256':sha256((path/'input.json').read_bytes()).hexdigest(),
            'pairs':1675,'excluded_pairs':1675-len(remaining),'survivors':sorted(map(list,remaining)),
            'prime_certificates':[{'path':f"prime-{r['prime']}.json",'sha256':sha256((path/f"prime-{r['prime']}.json").read_bytes()).hexdigest()} for r in certificates],
            'cpu_seconds':time.process_time()-started,
            'scope':'Exact projective branch-sextic exclusions for the frozen67-by25 family bank, including all rational parameters and repeated-root positions. Other trace classes and higher-pole carriers remain outside this result.'}
    write_new(path/'result.json',result);print(json.dumps({k:result[k] for k in ('pairs','excluded_pairs','survivors','cpu_seconds')}),flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('mode',choices=('prepare','run'));p.add_argument('--output',type=Path,default=DEFAULT)
    a=p.parse_args();prepare(a.output) if a.mode=='prepare' else run(a.output)


if __name__=='__main__':main()
