"""Small arithmetic, cache-integrity and worker-lifecycle regressions."""

from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction
import json
from pathlib import Path
import os
import signal
import subprocess
import sys
from tempfile import TemporaryDirectory
import time
import unittest

CAS = Path(__file__).resolve().parents[1] / "cas"
sys.path.insert(0, str(CAS))

from research_runtime.arithmetic import ArithmeticContext, CurveModel, TwoTorsionContext
from research_runtime.binary import BinaryBasis, kernel_masks, pack
from research_runtime.store import FactStore, digest
from research_runtime.supervisor import Limits, run, supervise_source


class FactTests(unittest.TestCase):
    def test_concurrent_discovery_once_and_readonly_miss(self):
        with TemporaryDirectory() as directory:
            store = FactStore(directory)
            calls = []
            def builder():
                calls.append(1)
                time.sleep(.02)
                return {"order": 17}
            with ThreadPoolExecutor(4) as pool:
                values = list(pool.map(lambda _: store.discover("order", [1, 2], builder), range(8)))
            self.assertEqual(calls, [1])
            self.assertEqual(values, [{"order": 17}]*8)
            with self.assertRaises(FileNotFoundError):
                store.require("order", [1, 3])
            with self.assertRaises(FileNotFoundError):
                store.require("order", [1, 2], version="2")

    def test_tamper_and_interrupted_publication(self):
        with TemporaryDirectory() as directory:
            store = FactStore(directory)
            def failed():
                raise RuntimeError("interrupted discovery")
            with self.assertRaises(RuntimeError):
                store.discover("field", "a", failed)
            self.assertIsNone(store.get("field", "a"))
            store.discover("field", "a", lambda: [1])
            blob = next((Path(directory)/"objects").glob("*.json"))
            blob.write_text("{}")
            with self.assertRaisesRegex(ValueError, "corrupt"):
                store.discover("field", "a", lambda: [2])

    def test_exact_identity(self):
        self.assertEqual(digest(Fraction(2, 4)), digest(Fraction(1, 2)))
        with self.assertRaises(TypeError):
            digest(.5)
        first = TwoTorsionContext((1, -1, 0, 1), ("theta",))
        second = TwoTorsionContext((1, -1, 0, 1), ("permuted",))
        self.assertNotEqual(first.key, second.key)

    def test_portable_facts_and_conflicts(self):
        with TemporaryDirectory() as first, TemporaryDirectory() as second:
            source, target = FactStore(first), FactStore(second)
            source.discover("order", "field", lambda: {"basis": [1, 2, 3]})
            snapshot = source.snapshot()
            target.import_snapshot(snapshot)
            self.assertEqual(target.require("order", "field"), {"basis": [1, 2, 3]})
            snapshot["facts"][0]["record"]["value"]["basis"][0] = 99
            with self.assertRaises(ValueError):
                target.import_snapshot(snapshot)

    def test_models_and_incomplete_factorization(self):
        model = CurveModel((0, 0, 0, -1, 1))
        self.assertEqual(model.discriminant, -368)
        self.assertEqual(model.two_division_polynomial, ("64", "-16", "0", "1"))
        self.assertTrue(model.contains((0, 1)))
        with self.assertRaises(ValueError):
            ArithmeticContext(model, model, (1,0,0,0), ((2,4),),
                              TwoTorsionContext(model.two_division_polynomial), (0,1,0))


    def test_raw_search_context_defers_arithmetic_and_upgrades_state(self):
        from research_runtime.search_state import raw_state
        from research_runtime.finite_reduction import ReductionCache
        with TemporaryDirectory() as directory:
            store=FactStore(directory);cache=ReductionCache(store)
            state=raw_state([0,0,1,-1,0],[(0,0)],cache=cache,prime_bound=43)
            self.assertEqual(state.rank,1)
            self.assertFalse(state.arithmetic.prepared)
            self.assertEqual(state.arithmetic.record()["minimal_model"],None)
            with self.assertRaises(ValueError):state.arithmetic.require_prepared()
            self.assertTrue(all('order' not in f['record']['key']['namespace'] and 'bnf' not in f['record']['key']['namespace']
                                for f in store.snapshot()['facts']))
            from research_runtime.sage_arithmetic import SageArithmetic
            try:import sage.all
            except ImportError:self.skipTest('arithmetic upgrade requires Sage; raw setup assertions passed')
            prepared=SageArithmetic(store).prepare_context(state.arithmetic,factor_primes=[37],discover=True)
            upgraded=state.with_arithmetic(prepared)
            self.assertEqual(upgraded.reductions,state.reductions)
            self.assertEqual(upgraded.basis,state.basis)
            self.assertTrue(upgraded.arithmetic.prepared)


class BinaryTests(unittest.TestCase):
    def test_incremental_dependencies_and_surviving_combination(self):
        basis = BinaryBasis(3)
        for value in (0b101, 0b011):
            basis, dependency = basis.append(value)
            self.assertIsNone(dependency)
        extended, dependency = basis.append(0b110)
        self.assertEqual(dependency, 0b111)
        self.assertEqual(extended.rank, 2)
        self.assertEqual(basis.column_count, 2)
        self.assertEqual(kernel_masks([[1],[1]]), (3,))
        with self.assertRaises(ValueError):
            pack([2])


class IncrementalStateTests(unittest.TestCase):
    def test_only_new_columns_and_immutable_point_admission(self):
        from research_runtime.finite_reduction import ReductionCache
        from research_runtime.mw_state import MWState,ParityLattice
        model=CurveModel((0,0,0,-1,1));algebra=TwoTorsionContext(model.two_division_polynomial)
        context=ArithmeticContext(model,model,(1,0,0,0),((2,4),(23,1)),algebra,(0,1,0))
        with TemporaryDirectory() as directory:
            cache=ReductionCache(FactStore(directory))
            state=MWState.empty(context,cache=cache,primes=[3,5,7,11],no_two_torsion_prime=3)
            def forbidden_escalation():
                raise AssertionError("independent point requested more primes")
                yield 13
            first=state.adjoin((0,1),cache=cache,extra_primes=forbidden_escalation())
            self.assertEqual((state.rank,first.rank),(0,1))
            self.assertEqual(cache.quotient_builds,4)
            count=cache.point_evaluations
            repeated=first.adjoin((0,1),cache=cache)
            self.assertEqual(repeated.rank,1)
            self.assertEqual(repeated.observations[-1].status,"KNOWN_POINT_UP_TO_SIGN")
            self.assertEqual(cache.point_evaluations,count)
            self.assertEqual(first.kummer_classes,(("0","-1","0"),))
            doubled=first.adjoin(("1/4","-7/8"),cache=cache)
            self.assertEqual(doubled.rank,1)
            self.assertEqual(cache.quotient_builds,4)
            self.assertEqual(cache.point_evaluations,count+4)
            with self.assertRaises(ValueError):first.adjoin((0,2),cache=cache)
            with self.assertRaises(ValueError):
                first.with_geometry((("1",),),height_kind="approximate",parity=ParityLattice(generic_coordinates=((1,0),)))

    def test_cache_restart_checks_retained_signatures(self):
        from research_runtime.finite_reduction import ReductionCache
        with TemporaryDirectory() as directory:
            first=ReductionCache(FactStore(directory))
            expected=first.signature((0,0,0,-1,1),[(0,1)],5)
            second=ReductionCache(FactStore(directory))
            self.assertEqual(second.signature((0,0,0,-1,1),[(0,1)],5),expected)
            self.assertEqual(second.quotient_builds,0)


class SupervisorTests(unittest.TestCase):
    def limits(self, wall=2, rss=256_000_000):
        return Limits(wall, rss, terminate_grace_seconds=.1, poll_seconds=.01)

    @staticmethod
    def live(pid):
        try:
            return Path(f"/proc/{pid}/stat").read_text().rsplit(")", 1)[1].split()[0] != "Z"
        except FileNotFoundError:
            return False

    def test_success_failure_timeout_and_retained_inputs(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            for code, outcome in [("print(1)", "completed"), ("raise ValueError()", "backend_failure"),
                                  ("import time;time.sleep(5)", "strict_wall_timeout")]:
                result = run([sys.executable,"-c",code],limits=self.limits(.3),log_path=root/"log")
                self.assertEqual(result["outcome"], outcome)
            result = supervise_source(sys.executable,
                "from pathlib import Path\nPath(OUTPUT_PATH).write_text('{\"ok\":true}')", {},
                root/"result.json", root/"worker.log", timeout=2, rss_limit_bytes=256_000_000)
            self.assertEqual(result["outcome"], "completed")
            self.assertTrue((root/"worker.log.input.json").exists())
            self.assertTrue((root/"worker.log.worker.py").exists())
            with self.assertRaises(FileExistsError):
                run([sys.executable,"-c","pass"],limits=self.limits(),log_path=root/"log",result_path=root/"result.json")

    def test_descendant_cleaned_after_parent_exits(self):
        with TemporaryDirectory() as directory:
            root=Path(directory)
            code="import subprocess,sys; p=subprocess.Popen([sys.executable,'-c','import time;time.sleep(60)']);print(p.pid,flush=True)"
            result=run([sys.executable,"-c",code],limits=self.limits(),log_path=root/"log")
            pid=int((root/"log").read_text().strip())
            self.assertEqual(result["outcome"],"completed")
            self.assertFalse(self.live(pid))

    def test_rss_includes_children(self):
        with TemporaryDirectory() as directory:
            root=Path(directory)
            code="import subprocess,sys,time;p=subprocess.Popen([sys.executable,'-c','import time;data=bytearray(80_000_000);time.sleep(5)']);print(p.pid,flush=True);time.sleep(5)"
            result=run([sys.executable,"-c",code],limits=self.limits(3,50_000_000),log_path=root/"log")
            self.assertEqual(result["outcome"],"strict_rss_limit")
            self.assertFalse(self.live(int((root/"log").read_text().strip())))

    def test_nested_session_memory_and_cleanup(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            code = "import subprocess,sys,time;p=subprocess.Popen([sys.executable,'-c','import time;data=bytearray(80000000);time.sleep(5)'],start_new_session=True);print(p.pid,flush=True);time.sleep(5)"
            result = run([sys.executable, "-c", code], limits=self.limits(3, 50000000), log_path=root/"log")
            self.assertEqual(result["outcome"], "strict_rss_limit")
            self.assertFalse(self.live(int((root/"log").read_text().strip())))

    def test_watchdog_survives_supervisor_sigkill(self):
        with TemporaryDirectory() as directory:
            root=Path(directory)
            code=("from research_runtime.supervisor import run,Limits;from pathlib import Path;import sys;"
                  f"run([sys.executable,'-c','import os,time;print(os.getpid(),flush=True);time.sleep(60)'],limits=Limits(60,256000000,terminate_grace_seconds=.1),log_path=Path({str(root/'log')!r}),checkpoint_path=Path({str(root/'checkpoint')!r}))")
            parent=subprocess.Popen([sys.executable,"-c",code],env={**os.environ,"PYTHONPATH":str(CAS)},stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
            try:
                deadline=time.monotonic()+5
                while not (root/"checkpoint").exists() or not (root/"log").read_text().strip():
                    if time.monotonic()>deadline:
                        self.fail("worker did not start")
                    time.sleep(.02)
                pid=int((root/"log").read_text().strip())
                parent.kill();parent.wait(timeout=2)
                deadline=time.monotonic()+3
                while self.live(pid) and time.monotonic()<deadline:
                    time.sleep(.02)
                self.assertFalse(self.live(pid))
            finally:
                if parent.poll() is None:
                    parent.kill();parent.wait()
                parent.stderr.close()


try:
    from sage.all import pari
    HAVE_SAGE=True
except ImportError:
    HAVE_SAGE=False


@unittest.skipUnless(HAVE_SAGE,"Sage adapter tests: use sage -python")
class SageArithmeticTests(unittest.TestCase):
    def test_transported_maximal_order_and_bound_polmod_variable(self):
        from unittest.mock import patch
        from research_runtime import pari_context
        from research_runtime.sage_arithmetic import SageArithmetic
        with TemporaryDirectory() as directory:
            adapter = SageArithmetic(FactStore(directory))
            with patch.object(pari_context, "_arithmetic", adapter):
                source = pari("(x-20)^3-(x-20)+1")
                reduced, image = pari_context.prepared_polredabs(source, [23])
                self.assertEqual(image.mod(), reduced)
                self.assertEqual(pari.subst(source, "x", image), 0)
                nf = pari_context.prepared_nf(reduced, [23], discover=False)
                self.assertEqual(nf[2], -23)
                self.assertEqual(nf.nf_get_pol(), reduced)
                self.assertEqual(len(pari_context.prepared_prime_ideals(nf, 2)), 1)

    def test_cached_order_local_characters_and_no_implicit_bnf(self):
        from research_runtime.sage_arithmetic import SageArithmetic
        from research_runtime.local_kummer import LocalSquareclasses
        with TemporaryDirectory() as directory:
            arithmetic=SageArithmetic(FactStore(directory))
            context=arithmetic.prepare([0,0,0,-1,1],factor_primes=[2,23],discover=True)
            nf=arithmetic.nf(context.two_torsion,factor_primes=[2,23],discover=True)
            self.assertEqual(nf[2],-23)
            restarted=SageArithmetic(FactStore(directory))
            self.assertEqual(restarted.prepare([0,0,0,-1,1]),context)
            self.assertEqual(restarted.nf(context.two_torsion)[2],-23)
            self.assertEqual(restarted.scheduling_features(context.two_torsion)["components"][0]["class_group_status"],"UNKNOWN")
            with self.assertRaises(ValueError):
                restarted.bnf(context.two_torsion,discover=True)
            t=pari("Mod(y,y^3-16*y+64)")
            for p in [2,3,5,23]:
                local=LocalSquareclasses(nf,p,arithmetic=arithmetic,context=context.two_torsion)
                for value in [t,t+1,2*t-3,(t+2)**2,pari(2),pari(3),1/(t+1)]:
                    expected=all(pari.nfislocalpower(nf,P,value,2) for P in local.primes)
                    self.assertEqual(local.is_square(value),expected,(p,str(value)))
                classes=[t,t+1,t*(t+1),(t+2)**2]
                basis,coordinates=local.coordinates(classes)
                for value,bits in zip(classes,coordinates):
                    for b,bit in zip(basis,bits):
                        if bit:value/=b
                    self.assertTrue(local.is_square(value))

    def test_split_algebra_and_factor_support_failure(self):
        from research_runtime.sage_arithmetic import SageArithmetic
        with TemporaryDirectory() as directory:
            arithmetic=SageArithmetic(FactStore(directory))
            with self.assertRaisesRegex(ValueError,"cover"):
                arithmetic.prepare([0,0,0,-1,1],factor_primes=[2],discover=True)
            field=arithmetic.field(TwoTorsionContext((0,-1,0,1)),factor_primes=[2],discover=True)
            self.assertEqual([r["degree"] for r in field["components"]],[1,1,1])

    def test_complete_selmer_uses_prepared_field_and_reuses_result(self):
        from research_runtime.sage_arithmetic import SageArithmetic
        with TemporaryDirectory() as directory:
            arithmetic=SageArithmetic(FactStore(directory))
            context=arithmetic.prepare([0,0,0,-1,1],factor_primes=[2,23],discover=True)
            arithmetic.field(context.two_torsion,factor_primes=[2,23],discover=True)
            result=arithmetic.full_selmer(context,requirement="complete-selmer",discover=True)
            self.assertEqual(result["full_selmer_dimension"],1)
            # Independent PARI control on this tiny curve only. Production
            # workers do not invoke this monolithic entry point.
            self.assertEqual(len(pari.ellinit([0,0,0,-1,1]).ellrankinit().ell2cover()),1)
            restarted=SageArithmetic(FactStore(directory))
            def forbidden(*args,**kwargs):raise AssertionError("rediscovered BNF")
            restarted.bnf=forbidden
            self.assertEqual(restarted.full_selmer(context),result)


if __name__ == "__main__":
    unittest.main()
