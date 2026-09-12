"""Cheap commissioning regressions; no MW16 class or point search."""
from fractions import Fraction as F
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from transfer_common import *
import transfer_constructor as constructor
from transfer_lift import quadrics
import run_fresh_constructor_transfer as runner


class Gates(unittest.TestCase):
    def run_gate(self,control,fresh):
        with tempfile.TemporaryDirectory() as tmp:
            folder=Path(tmp);(folder/'control/class').mkdir(parents=True)
            write(folder/'control/class/seal.json',{})
            class Fake(runner.Controller):
                def __init__(self):
                    self.folder=folder;self.p=dict(reference_cap_seconds=7200,roster=[dict(id='control')]+[dict(id=f'fresh-{i:02d}') for i in range(1,9)]);self.calls=[];self.state={}
                def paired(self,case,cap):
                    self.calls.append((case,cap));success=control if case=='control' else case in fresh
                    return dict(transfer_success=success),1000
                def update(self,**kw):self.state.update(kw)
            c=Fake();c.execute();return c

    def test_failed_commission_never_starts_fresh(self):
        c=self.run_gate(False,[]);self.assertEqual(len(c.calls),1);self.assertEqual(c.state['status'],'STOPPED_COMMISSIONING_GATE')

    def test_two_misses_stop_without_replacement(self):
        c=self.run_gate(True,[]);self.assertEqual([x[0] for x in c.calls],['control','fresh-01','fresh-02']);self.assertEqual(c.state['status'],'STOPPED_FIRST_TWO_TRANSFER_GATE')

    def test_fresh_success_allows_fixed_remaining_six(self):
        c=self.run_gate(True,['fresh-02']);self.assertEqual(len(c.calls),9);self.assertTrue(all(cap==2000 for case,cap in c.calls[1:]));self.assertEqual(c.state['status'],'COMPLETED_FROZEN_PILOT')

    def test_kernel_provenance_and_duplicate(self):
        rows=[3,5,6,3];K=kernel(rows);self.assertEqual(len(K),2);self.assertTrue(all(xor(k,rows)==0 for k in K));self.assertEqual(rank(rows),2)


class Arithmetic(unittest.TestCase):
    def test_complete_dyadic_residue_degree_three(self):
        from sage.all import QQ,PolynomialRing,pari
        from research_runtime.local_kummer import LocalSquareclasses
        R=PolynomialRing(QQ,'z');nf=pari.nfinit(R([-1,-1,0,1]));L=LocalSquareclasses(nf,2)
        for P in pari.idealprimedec(nf,2):
            self.assertEqual(int(P[3]),3)
            gens,pairing=constructor.local_generators(nf,P)
            self.assertEqual(rank(list(map(pack,pairing))),5)
            self.assertEqual(rank([pack(L.signature(g)) for g in gens]),5)
            for g in gens:self.assertFalse(any(L.signature(g*g)))

    def test_toy_cold_field_preparation_and_transport(self):
        # Admission is replaced ONLY in this low-rank arithmetic API test.
        # The production entry point always requires an actual rank16 proof.
        with tempfile.TemporaryDirectory() as tmp:
            folder=Path(tmp)
            write(folder/'input.json',dict(id='toy',parameter='1',short_model=['-4','1'],generic_points=[['0','1']]*16))
            with patch.object(constructor,'admission',return_value=(None,{'rank_lower_bound':16,'scope':'TEST_STUB_ONLY'})):
                constructor.prepare(folder)
            d=read(folder/'field.json');R,f,nf=constructor.setup(d)
            self.assertEqual(len(d['generic_classes']),16)
            for g in d['generic_classes']:
                beta=list(map(F,g['beta_ascending']))+[F(0)]
                from blind_constructed_cover_check import norm
                self.assertEqual(norm(beta,list(map(F,d['cubic_ascending']))),F(g['norm']))
            policy=runner.policy();policy['constructor'].update(smooth_bound=97,anchor_prime_bound=97)
            write(folder/'policy.json',policy)
            # Empty toy anchor bank bypasses duplicate toy points only; exercise
            # the actual factor-base, primitive atom and parity implementation.
            with patch.object(constructor,'kernel',return_value=[]):
                collector=constructor.Collector(folder)
            collector.batch([(m,n) for n in range(1,5) for m in range(-4,5) if __import__('math').gcd(m,n)==1],folder/'toy-relations.json')
            self.assertGreater(len(collector.atoms),0)
            classes=dict(columns=[],atoms=[],classes=[dict(factor_labels=[['generic',0],['generic',1]])])
            write(folder/'classes.json',classes)
            write(folder/'classes-verified.json',dict(classes_sha256=sha(folder/'classes.json'),scope='TEST_STUB_ONLY'))
            constructor.compact(folder)
            from transfer_lift import verify_compact
            verify_compact(folder,0)
            self.assertEqual(read(folder/'compact-0-verified.json')['status'],'PASS_EXACT_SQUARE_EQUIVALENCE_PRODUCT_TREE')

    def test_cover_map_sign_and_denominator(self):
        from blind_constructed_cover_check import evaluate
        # E: W^2=Z^3-4Z+1, rational point (0,1), beta=-theta.
        f=[F(1),F(-4),F(0),F(1)];beta=[F(0),F(-1),F(0)]
        _,source=quadrics(beta,f)
        self.assertTrue(all(evaluate(q,[F(1),F(0),F(0),F(1)])==0 for q in source))


if __name__=='__main__':unittest.main()
