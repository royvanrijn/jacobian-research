#!/usr/bin/env sage-python
"""Reconstruct the marked MW16 sections over Q(t), then transport exactly.

The reconstruction uses the sanitized template and generic 11952 sections.
No exceptional fibre, point or parameter is an input.
"""
import argparse
from importlib.machinery import SourceFileLoader
from pathlib import Path
import sys
from sage.all import PolynomialRing, QQ, ZZ, Matrix, EllipticCurve

ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT/'elliptic-curves/cas'
sys.path.insert(0,str(CAS))
import certify_compact_r17_candidates as cert
import audit_r17_constant_scaling as scaling
from research_runtime.store import checkpoint
builder = SourceFileLoader('compact_mw16_marking',str(CAS/'build_a1_mw16_target_free_parameter_candidates.sage')).load_module()
helpers = SourceFileLoader('compact_mw16_identity',str(CAS/'reduce_r17_family_base.sage')).load_module()
TEMPLATE = ROOT/'elliptic-curves/data/a1_mw16_family_template_v1.json'
DIRECTORY = ROOT/'artifacts/local/elliptic-curves/compact-mw16-sections-v1'
BASE_DIRECTORY = ROOT/'artifacts/local/elliptic-curves/compact-mw16-base-v2'
R = PolynomialRing(QQ,'z'); z = R.gen(); K = R.fraction_field()


def sources():
    paths = (Path(__file__).resolve(), TEMPLATE, builder.MODEL, builder.CHORD,
             Path(builder.__file__).resolve(), Path(helpers.__file__).resolve(),
             Path(cert.__file__).resolve(), Path(scaling.__file__).resolve())
    return {str(p.relative_to(ROOT)): cert.hashed(p) for p in paths}


def encode(value):
    q = K(value)
    return {'numerator_coefficients_low_to_high': list(map(str,q.numerator().list())),
            'denominator_coefficients_low_to_high': list(map(str,q.denominator().list()))}


def decode(record):
    return K(R(record['numerator_coefficients_low_to_high']))/R(record['denominator_coefficients_low_to_high'])


def coefficients(row):
    template = cert.read(TEMPLATE)
    presentation = next(p for p in template['presentations'] if p['presentation_id']==row['presentation_id'])
    if presentation['fibration_id'] != row['fibration_id']:
        raise ArithmeticError('fibration binding changed')
    A,B = (R(presentation['pencil'][key]) for key in ('A_coefficients_low_to_high','B_coefficients_low_to_high'))
    new_A,new_B = (R(row[key]) for key in ('A_coefficients_low_to_high','B_coefficients_low_to_high'))
    a,b,c,d = map(QQ,row['base_matrix_a_b_c_d']); u = QQ(row['scale_u'])
    if not u or a*d==b*c:
        raise ArithmeticError('singular coordinate map')
    if helpers.homogeneous(A,8,a,b,c,d)!=u**4*new_A or helpers.homogeneous(B,12,a,b,c,d)!=u**6*new_B:
        raise ArithmeticError('coefficient identity failed')
    gram = Matrix(QQ,presentation['generic_height_gram'])
    if gram.nrows()!=16 or not gram.is_positive_definite():
        raise ArithmeticError('generic Gram certificate is not positive definite')
    return presentation,new_A,new_B,K(a*z+b)/(c*z+d),K(c*z+d),u


def reconstruct(row):
    presentation,new_A,new_B,lam,den,u = coefficients(row)
    native = cert.read(builder.MODEL)
    T = PolynomialRing(QQ,'t'); t = T.gen(); L = T.fraction_field()
    oldA,oldB = (T(native['weierstrass_model'][key]) for key in ('A_coefficients_low_to_high','B_coefficients_low_to_high'))
    E = EllipticCurve(L,[oldA,oldB])
    basis = [E(builder.polynomial_from_record(s['X'],T),builder.polynomial_from_record(s['Y'],T)) for s in native['sections']['records']]
    marking = presentation['source_marking']
    def combination(word):
        if len(word)!=len(basis): raise ArithmeticError('marking dimension mismatch')
        return sum((ZZ(n)*P for n,P in zip(word,basis) if n),E(0))
    chord = SourceFileLoader('compact_mw16_chord',str(builder.CHORD)).load_module()
    h,nx,m0,quartic,childA,childB = builder.child_geometry(combination(marking['trace_section_basis_w']),oldA,oldB,T,chord)
    for value,key in ((childA,'A_coefficients_low_to_high'),(childB,'B_coefficients_low_to_high')):
        if list(map(str,value.list())) != presentation['pencil'][key]:
            raise ArithmeticError('reconstructed source pencil differs')
    # Work over Q(z) at the compact base coordinate to avoid carrying huge
    # unreduced polynomial expressions in the old lambda coordinate.
    S = PolynomialRing(K,'v'); v = S.gen()
    fixed_m = S(m0)+lam*S(h)**2
    numerator = fixed_m**2-S(nx)
    sum_x,remainder = numerator.quo_rem(S(h)**2)
    if remainder: raise ArithmeticError('chord x-sum division failed')
    q = S([K(coef(lam)) for coef in quartic.list()])
    words = [marking['new_zero_source_section_basis_coordinates'],*marking['generic_source_section_basis_coordinates']]
    maps = marking['base_maps_lambda_of_old_t']
    if len(words)!=17 or len(maps)!=17: raise ArithmeticError('marking count changed')
    quartic_points=[]
    for i,(word,record) in enumerate(zip(words,maps)):
        P=combination(word)
        num=T(record['numerator_coefficients_low_to_high']); bot=T(record['denominator_coefficients_low_to_high'])
        if num.degree()>1 or bot.degree()>1: raise ArithmeticError('section base map is not Mobius')
        old_t=(num[0]-lam*bot[0])/(lam*bot[1]-num[1])
        if K(num(old_t))/bot(old_t)!=lam: raise ArithmeticError('base inverse failed')
        w=(2*K(P[0](old_t))-sum_x(old_t))/h(old_t)
        if w**2!=q(old_t): raise ArithmeticError('old section misses quartic')
        quartic_points.append((old_t,w))
        print('MW16 QUARTIC SECTION',row['presentation_id'],i,flush=True)
    t0,w0=quartic_points[0]
    if not w0: raise ArithmeticError('zero section is branch point')
    shifted=q(v+t0); ee,dd,cc,bb,aa=[K(shifted[i]) for i in range(5)]
    if ee!=w0**2: raise ArithmeticError('pointed constant term failed')
    a1=dd/w0; a2=cc-dd**2/(4*w0**2); a3=2*w0*bb; a4=-4*w0**2*aa; a6=a2*a4
    b2=a1**2+4*a2; b4=a1*a3+2*a4; b6=a3**2+4*a6
    if 81*(-(b2**2-24*b4)/48)!=K(childA(lam)) or 729*(-(-b2**3+36*b2*b4-216*b6)/864)!=K(childB(lam)):
        raise ArithmeticError('pointed quartic normalization failed')
    sections=[]
    for i,(old_t,w) in enumerate(quartic_points[1:]):
        delta=old_t-t0
        if not delta: raise ArithmeticError('section collides with zero')
        xg=(2*w0*(w+w0)+dd*delta)/delta**2
        yg=(4*w0**2*(w+w0)+2*w0*dd*delta+(2*w0*cc-dd**2/(2*w0))*delta**2)/delta**3
        old_x=9*(xg+b2/12); old_y=27*(yg+(a1*xg+a3)/2)
        if old_y**2!=old_x**3+childA(lam)*old_x+childB(lam): raise ArithmeticError('source short section failed')
        x=old_x*den**4/u**2; y=old_y*den**6/u**3
        if y**2!=x**3+new_A*x+new_B: raise ArithmeticError('compact section equation failed')
        sections.append({'basis_index':i,'X':encode(x),'Y':encode(y)})
        print('MW16 COMPACT SECTION',row['presentation_id'],i,flush=True)
    return sections,presentation['generic_height_gram']


def run(directory,presentation):
    output=directory/(presentation+'.json')
    if output.exists(): raise FileExistsError('preserve section export')
    row=cert.read(BASE_DIRECTORY/(presentation+'.json'))
    if row['status']!='PASS_EXACT_BASE_CHANGE': raise ArithmeticError('coordinate computation incomplete')
    sections,gram=reconstruct(row)
    record={**row,'status':'PASS_EXACT_EQUATIONS_AND_16_SECTION_TRANSPORTS','sections':sections,
            'generic_height_gram':gram,'generic_height_gram_determinant':str(Matrix(QQ,gram).det()),
            'generic_rank_lower_bound':16,'checker_sources':sources(),
            'coordinate_result_sha256':cert.hashed(BASE_DIRECTORY/(presentation+'.json'))}
    checkpoint(output,record)
    print('EXPORTED MW16',presentation,len(sections),flush=True)


def check(path,reconstruct_sections=False):
    data=cert.read(path)
    rows=data.get('families',[data])
    seen=set()
    for row in rows:
        if row['fibration_id'] in seen: raise ArithmeticError('duplicate fibration')
        seen.add(row['fibration_id'])
        if row['checker_sources']!=sources(): raise ArithmeticError('checker or source changed')
        presentation,A,B,base,den,u=coefficients(row)
        sections=row['sections'];gram=presentation['generic_height_gram']
        if len(sections)!=16 or [s['basis_index'] for s in sections]!=list(range(16)):
            raise ArithmeticError('section roster changed')
        for s in sections:
            x,y=decode(s['X']),decode(s['Y'])
            if y**2!=x**3+A*x+B: raise ArithmeticError('section identity failed')
        if gram!=row['generic_height_gram'] or str(Matrix(QQ,gram).det())!=row['generic_height_gram_determinant']:
            raise ArithmeticError('Gram binding changed')
        if reconstruct_sections and reconstruct(row)!=(sections,gram):
            raise ArithmeticError('full marking reconstruction differs')
        print('CHECKED MW16',row['presentation_id'],'full marking' if reconstruct_sections else 'equations and membership',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--directory',type=Path,default=DIRECTORY)
    p.add_argument('--presentation',default='a1-presentation-01')
    p.add_argument('--check',type=Path)
    p.add_argument('--reconstruct',action='store_true')
    a=p.parse_args()
    check(a.check,a.reconstruct) if a.check else run(a.directory,a.presentation)
