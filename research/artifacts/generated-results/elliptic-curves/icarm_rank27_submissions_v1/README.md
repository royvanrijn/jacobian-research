# Seven ICARM submissions

Use https://elliptic-rank.icarm.cloud/ and submit one curve at a time. Each linked sheet contains all three fields to paste: a-invariants, 27 point lines, and commentary. Leave primes of bad reduction blank for every curve. All models are globally minimal; all ranks are lower bounds, not exact ranks.

No Q-isomorphism matches in the 626-curve catalogue downloaded 2026-09-07T08:14:03.266229+00:00. No submissions were sent by this exporter.

| Submission sheet | Family | Parameter | Rank lower bound |
|---|---|---|---:|
| [new-20260906-40](new-20260906-40/SUBMIT.md) | 074d9 | 2818/1535 | 27 |
| [new-20260906-71](new-20260906-71/SUBMIT.md) | 103b2 | 3726/881 | 27 |
| [new-20260906-41](new-20260906-41/SUBMIT.md) | 11952 | -2448/11 | 27 |
| [new-20260906-72](new-20260906-72/SUBMIT.md) | 11952 | 2012/211 | 27 |
| [new-20260906-48](new-20260906-48/SUBMIT.md) | 11952 | 2828/2015 | 27 |
| [new-20260906-186](new-20260906-186/SUBMIT.md) | 11952 | 4286/1881 | 27 |
| [new-20260906-90](new-20260906-90/SUBMIT.md) | a1-fibration-01 | -1867/270 | 27 |

Each directory also contains plain-text fields and a submission.json API payload (with the optional primes key omitted). All sheets are combined in ALL_SUBMISSIONS.md. verification.json retains the exact transports, rank proofs and source hashes.

Replay from the repository root:

```sh
python3 elliptic-curves/cas/export_icarm_rank27_submissions.py --check
```
