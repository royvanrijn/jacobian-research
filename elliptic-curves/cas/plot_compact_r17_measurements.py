#!/usr/bin/env python3
"""Plot the independently checked initial lower bounds; requires matplotlib."""
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
DIRECTORY = ROOT / 'artifacts/generated-results/elliptic-curves'
SOURCE = DIRECTORY / 'compact_r17_initial_measurements_v1.json'


def main():
    data = json.loads(SOURCE.read_text())
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
    fig, axes = plt.subplots(1, 2, figsize=(11, 5.3))
    for cohort, label, colour, marker in (
        ('compact-r17-top64-v1', 'Height 4,096 (53 fresh workers)', '#2166ac', 'o'),
        ('compact-r17-h16384-v1', 'Height 16,384 (58 fresh workers)', '#b35806', '^')):
        rows = [r for r in data['rows'] if r['cohort'] == cohort]
        ranks = [r['rank_certificate']['rank_lower_bound'] for r in rows]
        axes[0].scatter([r['score_units'] / 10**12 for r in rows], ranks, label=label,
                        color=colour, marker=marker, s=36, alpha=.68)
        axes[1].scatter([r['short_coefficient_bits'] for r in rows], ranks,
                        color=colour, marker=marker, s=36, alpha=.68)
    for ax in axes:
        ax.set_ylim(16.5, 28.6); ax.set_yticks(range(17, 29))
        ax.axhline(28, color='#555555', linestyle='--', linewidth=1)
        ax.grid(axis='y', alpha=.16); ax.set_ylabel('Certified rank lower bound after initial 43 charts')
    axes[0].set_xlabel('Full-prime selection score (primes 5–4,093)')
    axes[1].set_xlabel('Largest coefficient numerator/denominator size (bits)')
    axes[0].text(.03, .95, 'Near-record target: ≥28', transform=axes[0].transAxes, va='top', color='#555555')
    axes[0].legend(loc='upper left', bbox_to_anchor=(0, .9), frameon=False, fontsize=9)
    fig.suptitle('A larger parameter sweep did not improve the best fresh initial lower bound', fontsize=13, y=.98)
    fig.text(.05, .035, 'Each mark is one completed fresh worker; known and reused curves are omitted. Lower bounds are not exact ranks.\n'
             'Chart boxes were generally incomplete at the fixed four-second allowance; score/height patterns do not establish causation.', fontsize=9)
    fig.tight_layout(rect=(0, .11, 1, .94))
    paths = []
    for suffix in ('png', 'pdf'):
        path = DIRECTORY / ('compact_r17_initial_measurements_v1.' + suffix)
        if path.exists():
            raise FileExistsError('do not replace a retained research figure')
        fig.savefig(path, dpi=180, metadata={'Creator': 'plot_compact_r17_measurements.py'} if suffix == 'pdf' else {})
        paths.append(path)
    manifest = {'source': str(SOURCE.relative_to(ROOT)), 'source_sha256': hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
                'plot_source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                'matplotlib_version': matplotlib.__version__,
                'outputs': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}
    (DIRECTORY / 'compact_r17_initial_measurements_figure_v1.json').write_text(json.dumps(manifest, indent=2) + '\n')


if __name__ == '__main__':
    main()
