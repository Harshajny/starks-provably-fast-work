# Action log

Plain record of every action taken in this project, newest last. "Public" means visible on
provably.fast or GitHub; everything else stayed on this machine.

## 2026-11-09

| # | What | Public? | Result |
|---|------|---------|--------|
| 1 | Downloaded the pvfast client zip and checked its SHA-256 against the published hash | No | Hash matched; unpacked to `pvfast-participant/` |
| 2 | Requested a sign-in link (`pvfast remote login --no-wait`), twice (first code expired) | No (login only) | Approved by you; agent id `agt_b7dde97e…` |
| 3 | Read the claim board, state digests, and threads (read-only) | No | Found 36 cells; 7 open |
| 4 | Ran 3 read-only sub-agents (leaderboard, open gaps, workflows) | No | Leaderboard URL was wrong (my error); gaps and workflows summarised |
| 5 | Cloned the public repo `starknet-innovation/starks.provably.fast` | No | In `starks-provably/` |
| 6 | Ran the starter kit example (`make_instance.py`, `certify.py`) | No | 2,240 of 2,240 challenges certified; toolchain works |
| 7 | Fetched the points ledger from the correct URL (`starks.provably.fast/data/mathematics.json`) | No | 62 people; top: nadai2010 13,678, binharold 8,103, py2755 3,370 |
| 8 | Checked candidate cells live | No | `lift-residual-growth` already finished by others; picked `problem-3a-L1` |
| 9 | Read py2755's living post for problem-3a-L1 | No | Saved as `living_3a.json` |
| 10 | Wrote own code and reproduced the 504 list over F_(p²) | No | 504 members, 1.3 s; see `work3a/README.md` |
| 11 | Claimed `problem-3a-L1` with the client's `claim` command in the rules thread (note: `work3a/claim_note.md`) | **Yes** | Post `bp1_13b554753d6237e40917d346`; promised the odd-fibration census by 2026-11-12. The client also created a local signing key at `~/.config/pvfast/agent_ed25519.json` (private, never uploaded) |
| 12 | Step 3: wrote `binom.c` and ran it once over all subsets for every one- and two-monomial word | No | 20 s; exact maximum 504 (even), 7 (odd), 0 (mixed) over all F_(p²) coefficients |
| 13 | Deduplicated and cross-checked (`dedupe.py`, `crosscheck.py`) | No | Matches step 1 set for set; 497/152/7 re-verified by interpolation. One check-script bug found and fixed (it skipped "universal" sets) |
| 14 | Posted the two-monomial result as a reply to py2755's living post (intent: upper-bound proof), 7 attachments (code, README_post.md, RESULTS.txt). First try failed: max 8 parts per post, so NOTES.md was merged into README_post.md and tally.py left out | **Yes** | Post `bp1_c90aa8fa159a5ee2221aef45` in `bt1_f4a8d9ad544c60349264ef90` |
| 15 | Filed the graph lead (DERIVED_FROM the BabyBear two-monomial census `rgr1_421d61…` and the 504 record `rgr1_60ae30…`) | **Yes** | Proposal `ba71306c-aabb-4d70-aed8-98e7de3e6572`, state PROPOSED (not yet admitted) |
| 16 | Checked the points ledger (generated 2026-11-08 09:02 UTC) | No | We are not on it yet; nothing is accepted. 50 pts ≈ rank 31, 150 ≈ 19, 500 ≈ 13 of 62 |
| 17 | Step 4 (three high monomials): wrote `trinom.c`; test showed ~4.5 GB of output, so rewrote as `trinom2.c` counting distinct members directly | No | Full pass 191 s; 1.55 M distinct-member records |
| 18 | Bounded every word per triple; `linesmax.c` found 72 candidate words; `listcount.c` recounted them exactly with 2 controls | No | Controls 504/497 exact; largest trinomial list 102; all trinomials ≤ 311; so every ≤3-monomial word ≤ 504 |
| 19 | Stopped one slow Python aggregation (`tri_agg.py`, replaced by `tri_stats.py` + `linesmax.c`) | No | No effect on results |
| 20 | Read site updates: ledger (still 8 Nov), recent threads, replies to us | No | py2755 ACCEPTed our post (`bp1_8789b2ae`, independent rerun). Their correction `bp1_5d784ef5`: the ≤3-monomial F_(p²) result already exists (`bp1_cc3e6eb2`, M-0833), so our step-4 draft is NOT posted as new |
| 21 | Rewrote step 4 as an outside review of py2755's `bp1_cc3e6eb2` (M-0833), since that result already existed. First attempt aborted when one blob upload returned no id; confirmed nothing posted, retried once | **Yes** | Review `bp1_0d14e91287c7f17534661b1c`, verdict ACCEPT, with code blobs (`trinom2.c`, `linesmax.c`, `listcount.c`, `RESULTS_step4.txt`); text in `work3a/review_m0833.md` |
| 22 | Four-monomial feasibility probe (`work3a/quadprobe.c`) on 28% of subsets for support t^18, t^22, t^26, t^30 | No | 33 s on 3 threads; ~190 M canonical one-word records per support expected; rank ≤ 2 sets ~1.2 M. Feasible with a two-pass heavy-word filter |
| 23 | Step 5: four monomials, the three mod-4 supports. `quad.c` (two-pass sketch), `concur3.c` (plane/line concurrency), `weigh.c`, `listcount.c` recounts | No | Every four-term word t^r W(t^4) has ≤ 338 members; heavy orbits are exactly 69, 1, 1; class max stays 504 (binomial). ~15 min machine time total |
| 24 | Posted the four-monomial result as a reply to the living post (6 attachments) and filed its graph lead | **Yes** | Post `bp1_92f0a2685f21312926ffcade`; lead proposal `423309ca-56ea-4bfd-932f-b73dce4943e3` (PROPOSED) |
| 25 | Checked leaderboard and replies | No | Ledger still dated 2026-11-08 09:02 UTC (we are not on it yet); no new replies since py2755's ACCEPT; our four-monomial post is the newest post in the thread |
| 26 | Reviewed py2755's problem-5 post `bp1_c741c3aa` (BabyBear, k = 16, A = 18, ≤3 monomials ≤ 715): own code `review_p5/bbbin.c` + `bbdedupe.py` reproduced the full two-monomial class exhaustively. First post attempt failed (short post id); retried once with the full id | **Yes** | Review `bp1_1c7705b0c2c864481e79599d`, verdict PARTIAL: binomial max exactly 715 at x^16+x^30 and x^17+x^31; trinomial census not rerun |
| 27 | Installed the GitHub CLI (`brew install gh`), staged `repo/` (33 files, 236 KB: our code, READMEs, notes, small outputs, post texts; no binaries, raw outputs, others' files or keys) and ran `git init` | No | Waiting for you to log in with `gh auth login` before creating and pushing the repo |
| 28 | 2026-11-10 | Pushed repo/ to public GitHub repo https://github.com/Harshajny/starks-provably-fast-work (33 files, noreply author email) | public |
| 29 | 2026-11-10 | problem-3a-upper (double-points cell): built native-code counter (work3u/native.c), validated vs M-0908/M-0911, three-monomial census: <= 54, so <= 58 for <= 3 monomials. Local only | local |
| 30 | 2026-11-11 | Pushed problem-3a-upper/ to GitHub (commit daf6b23); posted the three-monomial native result as reply to M-0911: post bp1_7e2fc4439123e9f34673a6fa; lead proposal d6d6ffa2-db26-408f-a751-3bd58d382559 | public |
| 31 | 2026-11-11 | Four-monomial native census (quadnat.c, quadgen.c, concur4.c): <= 54, exact max 43, so <= 58 for at most four monomials. Draft work3u/post2_draft.md, not posted | local |

Public actions so far: rows 11, 14, 15, 21, 24, 26, 28, 30.
