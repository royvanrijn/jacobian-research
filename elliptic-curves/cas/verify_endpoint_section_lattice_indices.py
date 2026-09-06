#!/usr/bin/env python3
"""Sage-free independent verification of exact known-section lattice indices."""
from fractions import Fraction as F
from pathlib import Path
import math
import certify_compact_r17_candidates as cert
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'

def determinant(a):
    a=[[F(v) for v in row] for row in a];n=len(a);answer=F(1)
    if any(len(row)!=n for row in a):raise ArithmeticError('square determinant required')
    for i in range(n):
        pivot=next((j for j in range(i,n) if a[j][i]),None)
        if pivot is None:return F(0)
        if pivot!=i:a[i],a[pivot]=a[pivot],a[i];answer=-answer
        q=a[i][i];answer*=q
        for j in range(i+1,n):
            factor=a[j][i]/q
            for k in range(i+1,n):a[j][k]-=factor*a[i][k]
    return answer

def main():
    d=cert.read(ART/'endpoint_section_lattice_indices_v1.json');source=cert.read(ART/'endpoint_section_relations_v1.json')
    if d['status']!='PASS' or len(d['rows'])!=21 or source['completed_span_ranks']!=21 or any(cert.hashed(ROOT/n)!=h for n,h in d['sources'].items()):raise ArithmeticError('exact known-section input binding differs')
    count=0
    for r,s in zip(d['rows'],source['rows']):
        rank=s['generic_section_span_rank'];n=len(s['generic_points']);a=r['integer_coordinate_matrix'];u=r['unimodular_transform'];h=r['hermite_matrix'];den=r['coordinate_denominator'];coords=[[F(c,q['denominator']) for c in q['coefficients']] for q in s['relations']]
        if (r['family'],r['endpoint'],r['rank'])!=(s['family'],s['endpoint'],rank) or den!=math.lcm(*(v.denominator for row in coords for v in row)):raise ArithmeticError('coordinate rank/denominator differs')
        if a!=[[int(v*den) for v in row] for row in coords] or len(a)!=n or any(len(row)!=rank for row in a) or len(u)!=n or any(len(row)!=n for row in u):raise ArithmeticError('integer coordinate matrix differs')
        if any(type(v)!=int for matrix in (a,u,h) for row in matrix for v in row):raise ArithmeticError('integral matrix required')
        if abs(determinant(u))!=1 or h!=[[sum(u[i][k]*a[k][j] for k in range(n)) for j in range(rank)] for i in range(n)]:raise ArithmeticError('unimodular lattice transformation differs')
        if any(v for row in h[rank:] for v in row):raise ArithmeticError('lower transformed rows are not zero')
        det=abs(determinant(h[:rank]));index=F(den**rank)/det
        if index.denominator!=1 or index!=r['index_of_certified_subset_in_known_section_lattice']:raise ArithmeticError('exact subgroup index differs')
        for j,P in enumerate(s['independent_points']):
            matches=[i for i,Q in enumerate(s['generic_points']) if Q==P]
            if not matches or coords[matches[0]]!=[int(j==k) for k in range(rank)]:raise ArithmeticError('known subset inclusion differs')
        count+=index>1
    if count!=d['proper_subset_lattices']:raise ArithmeticError('proper known-lattice count differs')
    print('SAGE-FREE EXACT21 SECTION LATTICE CHECKS;',count,'PROPER SUBSETS',flush=True)
if __name__=='__main__':main()
