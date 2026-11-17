#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

echo '== exact coordinate and support checks =='
python verify_laurent_reduction.py
python verify_case1_reduction.py
python verify_firstblock_exact.py

echo '== exact Case 2 reconstruction and standard basis =='
python case2_exact_generate.py
Singular -q case2_compact4_exact.sing | tee case2_compact4_exact.out

echo '== exact Case 1 branch reductions =='
python case1_cascade_analysis.py | tee case1_cascade_analysis.out
python case1_cascade_w.py | tee case1_cascade_w.out

echo '== exact Case 1 standard bases =='
Singular -q case1_branch1_after_w_nfmod.sing | tee case1_branch1_after_w_nfmod.out
Singular -q case1_branch2_after_w_nfmod.sing | tee case1_branch2_after_w_nfmod.out

echo 'ALL_EXACT_UNIT'
