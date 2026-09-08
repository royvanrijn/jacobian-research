#!/usr/bin/env sage-python
"""Standalone exact generic base/model transport; exports only seventeen sections."""
import argparse, hashlib, json
from pathlib import Path
from sage.all import QQ, ZZ, PolynomialRing, EllipticCurve, matrix

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def main(args):
    p=json.loads(args.protocol.read_bytes())
    assert sha(Path(__file__))==p['script_sha256']
    assert sha(args.parent)==p['parent_sha256'] and sha(args.chart)==p['chart_sha256']
    source=json.loads(args.parent.read_bytes()); result=json.loads(args.chart.read_bytes())
    assert result['status']=='PASS'
    assert set(source)=={'a_invariants','basis_weierstrass_coordinates','source_sha256'}
    chart=result['selected']['state']; R=PolynomialRing(QQ,'s'); s=R.gen(); K=R.fraction_field()
    M=matrix(ZZ,2,chart['parameter_matrix']); a,b,c,d=M.list(); assert M.det()
    t=K((a*s+b)/(c*s+d)); u=QQ(chart['weierstrass_u']); h=K(c*s+d)
    def val(v): return K(R(v['numerator'])(t)/R(v['denominator'])(t))
    old=EllipticCurve(K,[val(v) for v in source['a_invariants']])
    A=R(chart['A_coefficients']); B=R(chart['B_coefficients'])
    assert A == -old.c4()/48*h**8/u**4
    assert B == -old.c6()/864*h**12/u**6
    E=EllipticCurve(K,[A,B]); assert E.j_invariant()==old.j_invariant()
    coords=[]
    for raw in source['basis_weierstrass_coordinates']:
        P=old([val(v) for v in raw]); x,y=P.xy()
        X=K((x+old.b2()/12)*h**4/u**2)
        Y=K((y+(old.a1()*x+old.a3())/2)*h**6/u**3)
        # Two original generic sections are rational functions, not polynomials.
        assert X.denominator() and Y.denominator()
        Q=E([X,Y]); assert not Q.is_zero(); coords.append([X,Y])
    assert len(coords)==17
    def serial(v):
        v=K(v); return dict(numerator=list(map(str,v.numerator())),denominator=list(map(str,v.denominator())))
    exported=dict(a_invariants=[serial(v) for v in E.a_invariants()],basis_weierstrass_coordinates=[[serial(v) for v in P] for P in coords],source_sha256=sha(args.chart))
    assert not args.output.exists();args.output.write_text(json.dumps(exported,indent=2)+'\n')
    audit=dict(status='PASS',inputs=dict(parent_sha256=sha(args.parent),chart_sha256=sha(args.chart)),export_sha256=sha(args.output),parameter_matrix=list(map(str,M.list())),weierstrass_u=str(u),section_count=17,coefficient_bits=dict(A=max(abs(ZZ(v)).nbits() for v in A),B=max(abs(ZZ(v)).nbits() for v in B)),boundary='Exact Q(s) isomorphism, j identity, and all seventeen section transports verified independently with Sage; no project arithmetic imports. This is the same fibration and parent. No point search, rank increase, global coefficient optimum, or new parent claimed.')
    assert not args.certificate.exists();args.certificate.write_text(json.dumps(audit,indent=2)+'\n');print('PASS generic chart and17 sections',audit['coefficient_bits'],flush=True)
if __name__=='__main__':
    ap=argparse.ArgumentParser()
    for n in ['protocol','parent','chart','output','certificate']: ap.add_argument('--'+n,type=Path,required=True)
    main(ap.parse_args())
