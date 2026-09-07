import copy
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'cas'))
from research_runtime.store import FactStore, encoded
from research_runtime.memory_store import MemoryFactStore
from research_runtime.finite_reduction import ReductionCache
from research_runtime.search_state import raw_state
from research_runtime.mw_state import MWState


class MemoryFactsTests(unittest.TestCase):
    def test_disk_and_memory_witnesses_are_interchangeable(self):
        with TemporaryDirectory() as directory:
            memory=MemoryFactStore();disk=FactStore(Path(directory))
            a=ReductionCache(memory);b=ReductionCache(disk)
            sa=raw_state([0,0,0,-7,10],[(1,2),(2,2)],cache=a,prime_bound=31)
            sb=raw_state([0,0,0,-7,10],[(1,2),(2,2)],cache=b,prime_bound=31)
            self.assertEqual(sa.record(),sb.record())
            self.assertEqual(memory.snapshot(),disk.snapshot())
            fresh=MemoryFactStore();fresh.import_snapshot(disk.snapshot())
            self.assertEqual(MWState.from_record(sa.record(),cache=ReductionCache(fresh)),sa)

    def test_values_are_immutable_and_tampered_snapshots_fail(self):
        store=MemoryFactStore();value=store.discover('test',{'x':1},lambda:{'rows':[1,2]})
        value['rows'].append(3)
        self.assertEqual(store.require('test',{'x':1}),{'rows':[1,2]})
        bad=copy.deepcopy(store.snapshot());bad['facts'][0]['record']['value']['rows'][0]=7
        with self.assertRaises(ValueError):MemoryFactStore().import_snapshot(bad)
        other=MemoryFactStore();other.discover('test',{'x':1},lambda:{'rows':[7]})
        duplicate=store.snapshot();duplicate['facts']+=other.snapshot()['facts']
        with self.assertRaises(ValueError):MemoryFactStore().import_snapshot(duplicate)
        with self.assertRaises(KeyError):store.require('missing',{})


if __name__=='__main__':unittest.main()
