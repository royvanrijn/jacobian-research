"""Two-stage Mestre parent/fibre specialization with coherent generic labels.

An outer u selects the parent. T selects its fibre. This adapter preserves
section identities across both parameters; independence is a separate gate.
"""
from fractions import Fraction as Q
from pathlib import Path
import json
from hashlib import sha256
from icarm_curve245_mestre import fermigier_roots
from probe_mestre_fermigier_two_section_local_continuation import reconstructed_second_line,normalized_data
from mestre_root_tuples import SixRootMestreConstruction
from nagao_1994 import primitive_visible_points,quartic_point_to_short_jacobian,short_jacobian_coefficients
ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/'artifacts/generated-results/elliptic-curves/mestre_component_coherent_sections_v1.json'
SOURCE=DATA.with_name('mestre_component_label_audit_v1.json')

def evaluate(coefficients,value):
    out=Q(0)
    for c in reversed(coefficients):out=out*value+Q(c)
    return out

def rational_function(record,value):
    den=evaluate(record['denominator'],value)
    if not den:raise ValueError('section coefficient has a pole at this outer parameter')
    return evaluate(record['numerator'],value)/den

def specialize(outer_u,fibre_T):
    u,t=Q(outer_u),Q(fibre_T);data=json.loads(DATA.read_text())
    if data['source_sha256']!=sha256(SOURCE.read_bytes()).hexdigest():raise ArithmeticError('global section source changed')
    v,c2,m2=reconstructed_second_line(u);c1,m1=normalized_data(u,v)[4:6];source=fermigier_roots(u,v)
    roots=tuple((r-source[0])/(source[1]-source[0]) for r in source)
    if roots!=tuple(rational_function(r,u) for r in data['roots']):raise ArithmeticError('labelled source roots changed')
    construction=SixRootMestreConstruction(roots);visible=primitive_visible_points(construction,t)
    indices=[2*construction.roots.index(roots[j])+s for j in data['fixed_source_root_order'] for s in (0,1)]
    quartic_points=[visible[i] for i in indices]
    for i,(c,m) in enumerate(((c1,m1),(c2,m2))):
        yy=evaluate([rational_function(a,u) for a in data['extra_ordinate_coefficients'][i]],t)/construction.quartic_square_scale
        quartic_points.append((c+m*t,yy))
    points=tuple(quartic_point_to_short_jacobian(construction,t,p) for p in quartic_points)
    model=short_jacobian_coefficients(construction,t)
    if not 4*model[3]**3+27*model[4]**2:raise ValueError('singular fibre')
    if any(y*y!=x*x*x+model[3]*x+model[4] for x,y in points):raise ArithmeticError('section image missed its fibre')
    return model,points
