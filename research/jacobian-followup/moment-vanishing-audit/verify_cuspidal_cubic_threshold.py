#!/usr/bin/env python3
"""Exact Wick/chart checks for the sharp cuspidal ternary cubic threshold10."""
from pathlib import Path
from math import comb, factorial
from functools import lru_cache
import hashlib
import json
import re
import time
import sympy as S

P = Path(__file__).resolve().parent
ACCEPTED = ['cuspidal_boundary_charts','cuspidal_generic_small_powers_v2']
PARTIAL = ['cuspidal_prefix_geometry','cuspidal_generic_chart',
           'cuspidal_generic_determinant_saturation','cuspidal_generic_small_powers']


def sha(name):
    return hashlib.sha256((P / name).read_bytes()).hexdigest()


def clean(name):
    return '\n'.join(x.strip() for x in (P / name).read_text().splitlines()
                     if x.strip() and not x.lstrip().startswith('//'))


def primitive(poly, variables):
    return S.Poly(poly, *variables, domain=S.QQ).clear_denoms()[1].primitive()[1].as_expr()


def counted_moment(n, variables):
    """Explicit pairing counts, independent of the generating Wick recurrence."""
    terms = {}
    for j in range(n+1):
        p,q,r = 2*j,3*(n-j),j
        for xy in range(min(p,q)+1):
            for xz in range(min(p-xy,r)+1):
                for yz in range(min(q-xy,r-xz)+1):
                    rem = [p-xy-xz,q-xy-yz,r-xz-yz]
                    if any(k < 0 or k % 2 for k in rem):
                        continue
                    loops = [k//2 for k in rem]
                    den = 2**sum(loops)*factorial(xy)*factorial(xz)*factorial(yz)
                    for k in loops:
                        den *= factorial(k)
                    num = comb(n,j)*factorial(p)*factorial(q)*factorial(r)
                    assert num % den == 0
                    key = tuple(loops+[xy,xz,yz])
                    terms[key] = terms.get(key,0)+num//den
    return S.Poly.from_dict(terms, variables, domain=S.QQ)


def main():
    start = time.monotonic()
    for stem in ACCEPTED+PARTIAL+['cuspidal_generic_modular_diagnostic']:
        rec = json.loads((P/(stem+'_run.json')).read_text())
        assert rec['terminal'] and rec['input_sha256'] == sha(stem+'.sing')
        assert rec['log_sha256'] == sha(stem+'_stdout.txt')
        if stem in ACCEPTED or stem == 'cuspidal_generic_modular_diagnostic':
            assert rec['completed'] and not rec['timed_out'] and rec['exit_code'] == 0
            assert '?' not in (P/(stem+'_stdout.txt')).read_text()
        else:
            assert not rec['completed']
            if stem == 'cuspidal_generic_small_powers':
                assert not rec['timed_out'] and '?' in (P/(stem+'_stdout.txt')).read_text()
            else:
                assert rec['timed_out']
    a,b,c,d,e,f = variables = S.symbols('a b c d e f')
    syms = {str(x):x for x in variables}
    data = json.loads((P/'cuspidal_moments_ten_input.json').read_text())
    assert data['even_moment_indices'] == [2,4,6,8,10]
    moments = []; actual_moments = {}
    weights = [(2,0),(0,2),(-4,6),(1,1),(-1,3),(-2,4)]
    for n, pol, scalar in zip(data['even_moment_indices'],data['primitive_moments'],data['nonzero_integer_scalars']):
        counted = counted_moment(n, variables)
        actual_moments[n] = counted.as_expr()
        given = S.Poly(S.sympify(pol,locals=syms),*variables,domain=S.QQ)
        assert counted == given*S.Integer(scalar)
        assert not given.as_expr().has(S.Float)
        for mon in given.monoms():
            assert [sum(k*w[i] for k,w in zip(mon,weights)) for i in range(2)] == [0,3*n]
        moments.append(given.as_expr())
    det = S.det(S.Matrix([[a,d,e],[d,b,f],[e,f,c]]))
    charts = json.loads((P/'cuspidal_covariance_charts_input.json').read_text())
    assert [x['label'] for x in charts] == ['AB11','A0B1D1','A1B0D1','A1B0D0','A0B0D1']
    bsource = (P/'cuspidal_boundary_charts.sing').read_text()
    gsource = (P/'cuspidal_generic_chart.sing').read_text()
    blog = clean('cuspidal_boundary_charts_stdout.txt')
    checked = []
    for row in charts:
        label = row['label']; vs = [syms[x] for x in row['variables']]
        sub = {syms[k]:S.sympify(v) for k,v in row['normalization'].items()}
        eliminated = syms[row['eliminated']]
        solution = S.sympify(row['solution'],locals=syms)
        assert S.expand(moments[0].subs(sub).subs(eliminated,solution)) == 0
        coeff = S.Poly(moments[0].subs(sub),eliminated).coeff_monomial(eliminated)
        assert coeff.is_number and coeff != 0
        expected = [primitive(S.expand(m.subs(sub).subs(eliminated,solution)),vs) for m in moments[1:]]
        actual = [S.sympify(x,locals=syms) for x in row['primitive_remaining_moments']]
        assert expected == actual
        delta = primitive(S.expand(det.subs(sub).subs(eliminated,solution)),vs)
        assert delta == S.sympify(row['determinant_primitive'],locals=syms)
        source = gsource if label == 'AB11' else bsource
        block = source.split('ring r'+label+'=')[1]
        native = block.split('ideal I=')[1].split(';')[0].split(',')
        assert [S.sympify(x.replace('^','**'),locals=syms) for x in native] == actual
        assert S.sympify(block.split('poly delta=')[1].split(';')[0].replace('^','**'),locals=syms) == delta
        if label != 'AB11':
            assert blog.split(label+'_DET_POWER\n')[1].splitlines()[0] == '0'
            assert blog.split(label+'_SATURATED\n')[1].splitlines()[0] == '1'
            # A second small exact calculation checks the boundary identities.
            gb = S.groebner(actual,*vs)
            assert gb.reduce(S.expand(delta**12))[1] == 0
        checked.append({'label':label,'eliminated_variable':str(eliminated),
                        'normalization':row['normalization'],
                        'determinant_power_certified':4 if label == 'AB11' else 12})
    # Bind the completed generic rational basis checkpoint to its accepted reuse.
    raw = (P/'cuspidal_generic_chart_stdout.txt').read_text()
    assert '?' not in raw
    assert raw.split('AB11_BASIS_SIZE\n')[1].splitlines()[0] == '52'
    basis = raw.split('AB11_BASIS\n')[1].split('AB11_DET_POWER')[0].strip()
    replay = (P/'cuspidal_generic_small_powers_v2.sing').read_text()
    assert replay.split('ideal G=')[1].split(';')[0].strip() == basis
    assert 'attrib(G,"isSB",1)' in replay and 'reduce(delta^k,G)' in replay
    assert replay.split('poly delta=')[1].split(';')[0] == gsource.split('poly delta=')[1].split(';')[0]
    log = clean('cuspidal_generic_small_powers_v2_stdout.txt').splitlines()
    tests = [[int(x) for x in log[i+1:i+3]] for i,line in enumerate(log)
             if line == 'DET_POWER_REMAINDER_SIZE']
    assert tests == [[1,11],[2,113],[3,113],[4,0]]
    # The two algebraic branches outside those five charts.
    assert moments[0].subs({a:0,d:0}) == 5*b**3
    assert det.subs({a:0,b:0,d:0}) == 0
    # Elementary reduction of an irreducible repeated-tangent cubic to the cusp.
    X,Y,Z,u,v,w,h = S.symbols('X Y Z u v w h',nonzero=True)
    shift = -w/(3*h)
    z1 = -v-2*w*shift-3*h*shift**2
    z0 = -u-v*shift-w*shift**2-h*shift**3
    F = X**2*Z+u*X**3+v*X**2*Y+w*X*Y**2+h*Y**3
    assert S.expand(F.subs({Y:Y+shift*X,Z:Z+z0*X+z1*Y}, simultaneous=True)-X**2*Z-h*Y**3) == 0
    transform = S.Matrix([X,Y+shift*X,Z+z0*X+z1*Y]).jacobian([X,Y,Z])
    assert S.simplify(transform.det()) == 1
    # Real Gaussian realization and all-order lattice selection of the witness.
    B = S.Matrix([[1/S.sqrt(2),S.I/S.sqrt(2),0],
                  [1/S.sqrt(2),-S.I/S.sqrt(2),0],[0,0,1]])
    assert S.simplify(B*B.T) == S.Matrix([[0,1,0],[1,0,0],[0,0,1]])
    assert B.det() != 0
    witness = []
    for n in range(1,31):
        value = 0
        for j in range(n+1):
            if 2*j == 3*(n-j) and j % 2 == 0:
                value += comb(n,j)*factorial(2*j)*S.factorial2(j-1)
        expected = 0 if n % 10 else comb(n,3*n//5)*factorial(6*n//5)*S.factorial2(3*n//5-1)
        assert value == expected
        witness.append(int(value))
    assert witness[:9] == [0]*9 and witness[9] == 1508855040000
    for n, moment in actual_moments.items():
        assert moment.subs({a:0,b:0,c:1,d:1,e:0,f:0}) == witness[n-1]
    # The symbolic divisibility proof uses gcd(5,3)=1 and parity of j=6r.
    nn,jj,kk = S.symbols('nn jj kk',integer=True)
    assert S.solve(2*jj-3*(nn-jj),jj) == [3*nn/5]
    assert S.expand((3*(10*kk)/5) % 2) == 0
    names = ['cuspidal_moments_input.json','cuspidal_moments_ten_input.json',
             'cuspidal_covariance_charts_input.json','explore_cuspidal_moments.py',
             'explore_cuspidal_moments_ten.py','prepare_cuspidal_covariance_charts.py',
             'cuspidal_cubic_source_audit.json']
    for stem in ACCEPTED+PARTIAL+['cuspidal_generic_modular_diagnostic']:
        names += ['run_'+stem+'.py',stem+'.sing',stem+'_stdout.txt',stem+'_run.json']
    cert = {'claim':'G3CUSP10_20261114','field':'C, characteristic zero',
        'theorem':'Every irreducible homogeneous ternary cubic with a cuspidal projective curve has a nonzero Gaussian pure moment by order10;10 is sharp.',
        'normal_form':'X^2 Z + Y^3 with an arbitrary invertible linear coordinate change',
        'covariance_requirement':'det([[a,d,e],[d,b,f],[e,f,c]])!=0; automatic for an invertible change of three Gaussian coordinates',
        'normalization_weights':{'a':[2,0],'b':[0,2],'c':[-4,6],'d':[1,1],'e':[-1,3],'f':[-2,4]},
        'normalization_charts':checked,
        'generic_completed_checkpoint':{'full_run_completed':False,'basis_phase_completed':True,'basis_size':52,'determinant_power':4,'power_remainder_sizes':tests},
        'witness':{'polynomial':'P=Z^2*T+W^3; Z=(X+iY)/sqrt(2), W=(X-iY)/sqrt(2)',
                   'M1_to_M10':witness[:10],'all_order_formula':'M_n=0 unless n=10r; M_(10r)=binomial(10r,6r)*(12r)!*(6r-1)!!',
                   'M10':witness[9],'singular_point':'[Z:W:T]=[0:0:1], nonisotropic',
                   'global_ternary_cubic_prefix_lower_bound':10},
        'provenance':{name:sha(name) for name in names},
        'boundary':'This closes homogeneous irreducible cuspidal cubics, with no pure-null polynomial in that class. Homogeneous irreducible nodal, nonsplit line-conic and nonhomogeneous singular-leading cubics remain separate. No global arbitrary-cubic upper bound is given. Degenerate Gram matrices in more Gaussian coordinates are not included. Modular determinant saturation is only a diagnostic and is not used to lift the theorem. The all-order witness formula is a written Wick-selection proof, not finite extrapolation.',
        'assurance':'Internal elementary normalization and exact rational CAS/Wick/combinatorial proof; no independent replay, external review or formal verification.'}
    path = P/'cuspidal_cubic_threshold_certificate.json';path.write_text(json.dumps(cert,indent=2)+'\n')
    rec = {'claim':cert['claim'],'completed':True,'cuspidal_pure_threshold':10,
        'global_cubic_prefix_lower_bound':10,'witness_M10':witness[9],
        'checker_sha256':sha(Path(__file__).name),'certificate_sha256':sha(path.name),
        'elapsed_seconds':time.monotonic()-start,'assurance':cert['assurance']}
    (P/'cuspidal_cubic_threshold_verification.json').write_text(json.dumps(rec,indent=2)+'\n')
    print(json.dumps(rec,indent=2))


if __name__ == '__main__':
    main()
