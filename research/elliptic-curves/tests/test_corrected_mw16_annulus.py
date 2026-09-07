"""Independent modular sums check corrected signed, annular and heap ordering."""
import json
import math
import shutil
import struct
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


@unittest.skipUnless(shutil.which('g++'), 'g++ required')
class CorrectedMW16Annulus(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory()
        cls.d = Path(cls.temp.name)
        cls.binary = cls.d/'scanner'
        subprocess.run(['g++', '-O3', '-std=c++17',
                        str(ROOT/'elliptic-curves/cas/newfamily/scan_corrected_mw16_annulus.cpp'),
                        '-o', str(cls.binary)], check=True, capture_output=True, timeout=20)
        cls.local = json.loads((ROOT/'artifacts/generated-results/elliptic-curves/mw16_local_score_corrections_v1.json').read_text())

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def base_caches(self, ties=False):
        paths, tables = [], []
        for index, primes in enumerate(([5, 7], [11, 13])):
            data = b'R17XS001'+struct.pack('<I', 2)
            for p in primes:
                flags = [int(p not in (5, 13)) for _ in range(p+1)]
                units = [0 if ties or not flags[i] else (i % 5-2)*10**6 for i in range(p+1)]
                tables.append((p, units, flags))
                data += struct.pack('<II', p, p+1)+struct.pack('<'+'q'*(p+1), *units)+bytes(flags)
            path = self.d/f'base-{index}.bin'
            path.write_bytes(data+b'ENDXSC01')
            paths.append(path)
        return paths, tables

    def local_cache(self, family, ties=False):
        tables = [t for t in self.local['tables'] if t['family'] == family]
        data = b'MW16LC01'+struct.pack('<I', 2)
        for t in tables:
            rows = t['entries']
            data += struct.pack('<III', t['prime'], t['modulus'], len(rows))
            data += struct.pack('<'+'q'*len(rows), *[0 if ties else r['correction_units'] for r in rows])
            data += bytes(0 if ties else int(r['restored_good']) for r in rows)
        path = self.d/'local.bin'
        path.write_bytes(data+b'ENDLC001')
        return path, tables

    def invoke(self, paths, local, sign, N, M, K, shard, shards, inner):
        return subprocess.run([str(self.binary), *map(str, paths),
                               *map(str, [sign, N, M, K, shard, shards, inner, 4]), str(local)],
                              text=True, capture_output=True, timeout=5)

    def check_frame(self, family, ties=False):
        paths, base = self.base_caches(ties)
        local, tables = self.local_cache(family, ties)
        for sign in (-1, 1):
            for N, M, inner, shard, shards in [(31, 29, 7, 2, 3), (341, 27, 125, 0, 1),
                                               (27, 341, 125, 0, 1), (171, 170, 31, 12, 13)]:
                rows = []
                for d in range(shard+1, M+1, shards):
                    for absolute in range(1, N+1):
                        if max(absolute, d) <= inner or math.gcd(absolute, d) != 1:
                            continue
                        n = sign*absolute
                        score = good = 0
                        for p, units, flags in base:
                            i = n*pow(d, -1, p) % p if d % p else p
                            score += units[i]
                            good += flags[i]
                        if not ties:
                            for t in tables:
                                p, m = t['prime'], t['modulus']
                                if d % p:
                                    i = n*pow(d, -1, m) % m
                                else:
                                    reciprocal = d*pow(n, -1, m) % m
                                    self.assertEqual(reciprocal % p, 0)
                                    i = m+reciprocal//p
                                score += t['entries'][i]['correction_units']
                                good += t['entries'][i]['restored_good']
                        rows.append((n, d, score, good))
                rows.sort(key=lambda r: (-r[2], -r[3], r[1], abs(r[0])))
                for K in (1, 7, 100000):
                    actual = self.invoke(paths, local, sign, N, M, K, shard, shards, inner)
                    self.assertEqual(actual.returncode, 0, actual.stderr)
                    expected = ['CORRECTED_MW16_ANNULUS_V1',
                                f'P {sign} {N} {M} {K} {shard} {shards} {inner} 4']
                    expected += ['C '+' '.join(map(str, r)) for r in rows[:K]]
                    expected += [f'S {len(rows)} {min(K, len(rows))}']
                    self.assertEqual(actual.stdout.splitlines(), expected)

    def test_all_five_families_and_both_local_charts(self):
        for family in sorted({t['family'] for t in self.local['tables']}):
            with self.subTest(family=family):
                self.check_frame(family)

    def test_exact_ties_and_annulus_order(self):
        self.check_frame('a1-fibration-01', True)

    def test_malformed_local_frames_fail_before_output(self):
        paths, _ = self.base_caches()
        local, _ = self.local_cache('a1-fibration-01')
        data = local.read_bytes()
        corruptions = [data[:-1], data+b'x', b'BADMAGIC'+data[8:],
                       data[:16]+struct.pack('<I', 25)+data[20:],
                       data[:24]+struct.pack('<q', 10**14)+data[32:]]
        for damaged in corruptions:
            local.write_bytes(damaged)
            r = self.invoke(paths, local, 1, 31, 29, 7, 0, 1, 7)
            self.assertNotEqual(r.returncode, 0)
            self.assertEqual(r.stdout, '')


if __name__ == '__main__':
    unittest.main()
