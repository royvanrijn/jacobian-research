"""General rational polynomial Gaussian coefficient-family safe loci.

Macaulay2 joint holonomic operations, ambient local cohomology and conormal
multiplicity compute the reduced full mixed-safe ideal. Ordinary saturation
then computes the exact normalization ideal over the supplied coefficient
algebra, retaining nilpotents. No numerical universal cutoff is computed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess

import sympy as sp

CORE = Path(__file__).with_name('gaussian_family_core.m2')
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()


def terms(values, count):
    out = {}
    for exponents, coefficient in values:
        if len(exponents) != count or any(type(i) is not int or i < 0 for i in exponents):
            raise ValueError('invalid polynomial exponents')
        if isinstance(coefficient, (float, bool)):
            raise ValueError('exact rational coefficients required')
        c = sp.Rational(coefficient)
        e = tuple(exponents)
        out[e] = out.get(e, sp.S.Zero)+c
    return {e: c for e, c in out.items() if c}


def normalize(data):
    allowed = {'dimension','parameter_count','phase_terms','parameter_ideal','precision_matrix'}
    if set(data)-allowed:
        raise ValueError('unknown family input keys')
    n, r = data['dimension'], data['parameter_count']
    if type(n) is not int or n < 1 or type(r) is not int or r < 1:
        raise ValueError('positive spatial dimension and parameter_count required')
    phase = terms(data['phase_terms'], n+r)
    ideal = [terms(g, r) for g in data.get('parameter_ideal', [])]
    entries = data.get('precision_matrix', [[int(i==j) for j in range(n)] for i in range(n)])
    if len(entries) != n or any(len(row) != n for row in entries):
        raise ValueError('precision matrix has wrong size')
    if any(isinstance(v, (float,bool)) for row in entries for v in row):
        raise ValueError('exact rational precision entries required')
    precision = sp.Matrix([[sp.Rational(v) for v in row] for row in entries])
    if precision != precision.T or not precision.det():
        raise ValueError('precision matrix must be symmetric and nondegenerate')
    return n, r, phase, ideal, precision


def expression(poly, variables):
    return '+'.join('('+str(c)+')*'+'*'.join([v+'^'+str(k) for v,k in zip(variables,e) if k] or ['1'])
                    for e,c in sorted(poly.items())) or '0'


def make_script(data):
    n,r,phase,ideal,precision = normalize(data)
    xs = ['x'+str(i) for i in range(n)]
    cs = ['c'+str(i) for i in range(r)]
    coords = xs+['t']+cs
    ds = ['d'+v for v in coords]
    pairs = ','.join(x+'=>'+d for x,d in zip(coords,ds))
    gens = []
    for i,x in enumerate(xs):
        linear = '+'.join('('+str(precision[i,j])+')*'+y for j,y in enumerate(xs) if precision[i,j])
        gens.append(ds[i]+'+('+linear+')-t*diff('+x+',P)')
    gens.append('dt-P')
    gens.extend('d'+c+'-t*diff('+c+',P)' for c in cs)
    weights = '{'+','.join(['1']*n+['0']*(r+1))+'}'
    base = ','.join(expression(g,cs) for g in ideal) or '0_B'
    return f'''load "{CORE.as_posix()}";
W=QQ[{','.join(coords+ds)},WeylAlgebra=>{{{pairs}}}];
P=sub({expression(phase,xs+cs)},W);
I=ideal({','.join(gens)});
if not isHolonomic(W^1/I) then error "exponential source not holonomic";
M=gQuotient Dintegration(0,I,{weights});
stdio << "SOURCE_PRESENTATION|" << toString presentation M << endl;
B=QQ[{','.join(cs)}];
parameterIdeal=ideal matrix(B,{{{{{base}}}}});
R=gFamilySafeIdeal(M,B,parameterIdeal);
stdio << "RESULT_JSON|" << toJSON R << endl;
print("GAUSSIAN_FAMILY_COMPLETE");exit 0;
'''


def parse_output(output):
    prefix='RESULT_JSON|'
    values=[line[len(prefix):] for line in output.splitlines() if line.startswith(prefix)]
    if len(values)!=1 or 'GAUSSIAN_FAMILY_COMPLETE' not in output:
        raise ValueError('Macaulay2 family calculation did not complete')
    result=json.loads(values[0])
    if result.get('status')!='DECIDED':
        raise ValueError('invalid family result status')
    return result


def run(data, workdir, timeout=30):
    if type(timeout) is not int or timeout < 0:
        raise ValueError('timeout must be nonnegative; zero means uncapped')
    script=make_script(data)
    out=Path(workdir).resolve();out.mkdir(parents=True,exist_ok=True)
    path=out/'family.m2';path.write_text(script)
    try:
        process=subprocess.run(['M2','--no-readline','--silent','--stop',str(path)],capture_output=True,text=True,
                               timeout=timeout or None)
    except subprocess.TimeoutExpired as exc:
        (out/'family.stdout').write_bytes((exc.stdout or b'')+(exc.stderr or b''))
        return dict(status='UNKNOWN',reason='joint-module family CAS timeout',input=data)
    output=process.stdout+process.stderr;(out/'family.stdout').write_text(output)
    if process.returncode:
        raise RuntimeError('Macaulay2 family operation failed; see retained stdout')
    result=parse_output(output)
    result.update(input=data,source_backend='joint exponential connection and FULL unlocalized spatial integration',
                  actual_image_kernel_computed=False,relative_Hom_base_change_assumed=False,
                  ordinary_coefficient_nilpotents_retained=True,
                  coefficient_ring='QQ[c]/parameter_ideal; normalization ideal is its exact I-torsion',
                  universal_numerical_cutoff_computed=False,
                  assurance='classical joint holonomic/local-cohomology/characteristic-multiplicity algorithms executed by Macaulay2; not independently certified finite operator proofs',
                  core_sha256=sha(CORE),cas_files_sha256={p.name:sha(p) for p in [path,out/'family.stdout']})
    return result


def universal_family(dimension, degree):
    """The ENTIRE bounded polynomial phase space, not a support restriction."""
    import itertools
    if type(dimension) is not int or dimension<1 or type(degree) is not int or degree<0:
        raise ValueError('positive dimension and nonnegative degree required')
    monomials=sorted(e for e in itertools.product(range(degree+1),repeat=dimension) if sum(e)<=degree)
    r=len(monomials)
    return dict(dimension=dimension,parameter_count=r,
                phase_terms=[[list(e)+[int(i==j) for j in range(r)],'1'] for i,e in enumerate(monomials)])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input',type=Path)
    parser.add_argument('output',type=Path)
    parser.add_argument('--workdir',type=Path)
    parser.add_argument('--timeout',type=int,default=30)
    args=parser.parse_args()
    result=run(json.loads(args.input.read_text()),args.workdir or args.output.with_suffix('.cas'),args.timeout)
    args.output.write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':
    main()
