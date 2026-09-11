"""Scheduling/resume regressions with the existing fake CAS; no arithmetic claim."""
import fcntl
from pathlib import Path
import sys
import tempfile
import threading
import unittest
from unittest.mock import patch

import numpy  # Load before patch.dict(sys.modules); NumPy cannot be reimported.
import test_broad_rank_search as fixtures
from test_broad_rank_search import io, runner
import run_broad_rank_capacity as capacity


class CapacityTests(unittest.TestCase):
    def test_drain_then_resume_increases_concurrency_without_repeating_jobs(self):
        fixture = fixtures.ControllerTests()
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            root, calls, execute = fixture.fixture(folder)
            original = {name: (folder/name).read_bytes() for name in ('plan.json', 'manifest.json')}
            original_guard = runner.guard

            def first(folder, job, command, seconds):
                result = execute(folder, job, command, seconds)
                if 'broad-cases' in job.parts:
                    (folder/'STOP').touch()
                return result

            with patch.object(runner, 'execute', first), patch.dict(sys.modules, fixture.modules()):
                capacity.run_at_capacity(runner, folder, 1)
            self.assertEqual(len(calls), 3)  # score, preflight, one fibre
            self.assertIs(runner.guard, original_guard)
            (folder/'STOP').unlink()
            barrier = threading.Barrier(3, timeout=10)
            mutex = threading.Lock()
            started = active = peak = 0

            def resumed(folder, job, command, seconds):
                nonlocal started, active, peak
                if 'broad-cases' not in job.parts or (job/'seal.json').exists():
                    return execute(folder, job, command, seconds)
                with mutex:
                    started += 1
                    number = started
                    active += 1
                    peak = max(peak, active)
                try:
                    if number <= 3:
                        barrier.wait()
                    return execute(folder, job, command, seconds)
                finally:
                    with mutex:
                        active -= 1

            with patch.object(runner, 'execute', resumed), patch.dict(sys.modules, fixture.modules()):
                capacity.run_at_capacity(runner, folder, 3)
                capacity.run_at_capacity(runner, folder, 3)
            self.assertEqual(peak, 3)
            self.assertEqual(started, 4)
            self.assertEqual(len(calls), 7)
            self.assertTrue((folder/'COMPLETE.json').exists())
            for name, content in original.items():
                self.assertEqual((folder/name).read_bytes(), content)
            for request in root.glob('broad-cases/*/batch-000/request.json'):
                self.assertEqual(io.read(request)['allowance'], 198)

    def test_original_guard_and_controller_lock_remain_enforced(self):
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            root, calls, execute = fixtures.ControllerTests().fixture(folder)
            guard = runner.guard
            with (folder/'controller.lock').open('a') as lock:
                fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
                with self.assertRaises(BlockingIOError):
                    capacity.run_at_capacity(runner, folder, 6)
            self.assertIs(runner.guard, guard)
            (root/'broad-inputs/exclusions.json').write_text('{}')
            with self.assertRaisesRegex(ValueError, 'frozen input/source changed'):
                capacity.run_at_capacity(runner, folder, 6)
            self.assertIs(runner.guard, guard)
            self.assertEqual(calls, [])

    def test_stop_and_capacity_limits(self):
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            root, calls, execute = fixtures.ControllerTests().fixture(folder)
            (folder/'STOP').touch()
            with patch.object(runner, 'execute', execute):
                capacity.run_at_capacity(runner, folder, 6)
            self.assertEqual(calls, [])
            for workers in (0, 9, True, None):
                with self.assertRaises(ValueError):
                    capacity.run_at_capacity(runner, folder, workers)

    def test_revision_binds_capacity_plan_manifest_and_adapter(self):
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            fixtures.ControllerTests().fixture(folder)
            revision = folder/'revision'
            revision.mkdir()
            (revision/'scheduler.py').write_text('retained scheduler')
            policy = {'schema': 'broad-rank.capacity.v1', 'folder': str(folder), 'workers': 6,
                      'plan.json_sha256': io.sha(folder/'plan.json'),
                      'manifest.json_sha256': io.sha(folder/'manifest.json'),
                      'adapter_sha256': io.sha(revision/'scheduler.py')}
            io.write(revision/'capacity.json', policy)
            digest = io.sha(revision/'capacity.json')
            self.assertEqual(capacity.validate_revision(folder, revision, digest)['workers'], 6)
            for path in (folder/'plan.json', folder/'manifest.json',
                         revision/'scheduler.py', revision/'capacity.json'):
                previous = path.read_bytes()
                path.write_text('{}')
                with self.assertRaises(ValueError):
                    capacity.validate_revision(folder, revision, digest)
                path.write_bytes(previous)

    def test_launch_refuses_stop_or_live_jobs(self):
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            root, calls, execute = fixtures.ControllerTests().fixture(folder)
            with patch.object(capacity, 'load_controller', return_value=runner):
                (folder/'STOP').touch()
                with self.assertRaisesRegex(ValueError, 'STOP present'):
                    capacity.launch(folder, 6, 'test')
                (folder/'STOP').unlink()
                import os
                io.write(root/'broad-cases/toy/batch-000/driver-lease.json',
                         {'pid': os.getpid(), 'token': io.token(os.getpid())})
                with self.assertRaisesRegex(ValueError, 'still draining'):
                    capacity.launch(folder, 6, 'test')
            self.assertFalse((folder/'capacity-revisions').exists())


if __name__ == '__main__':
    unittest.main()
