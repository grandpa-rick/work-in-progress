# Wake 282: sieve test for Q68-typeA-pair-criterion

Grade: **computed** for everything below. Nothing here is proved. The KL inputs are the same recall-only ones as in PROVE 281
(Lemma S, F1, and P_{a,b} = P_{a^{-1},b^{-1}}).

## Definitions used

The conventions are those of PROVE 281 (`../2026-11-07-q68-s15-separator.md` §3): GW convention C, left cell {Q fixed},
L(w) = Des P(w), and the superstandard cell Q_0. Verdicts are cell-independent by `dualconn.py`.

- **(ii)-element**: the definition is PROVE 281 §3 (ii), verbatim. It is y = (U, Q_0) with Des U ⊋ A (strictly) and μ(y, x) ≠ 0.
  It is not required that y < x. One of the two S_15 (ii)-elements has l(y) = l(x) + 1: it is the star image s_1·x.
- **(i)-element**: s·x ∈ Γ with L(sx) ⊄ A, so μ = 1.
- Lemma A requires μ(y, z) = 1 for every element of the (i) ∪ (ii) set. This is the PROVE 280 Step 2 completeness argument.

## Data reused / computed

- **S_15.** The 48 pairs are the cached `../prove281/certs15_*.json`. No KL was recomputed.
  A re-run of certifyn on 5433 reproduces the cached verdicts (`t15/`).
- **S_14 (control).** All 4 shapes were re-certified (`run14.py`, 2 s).
- **S_16 (new, out-of-sample).**
  - `listn.py 16` gives 227 (λ, A), 21 shapes and 410 candidate pairs (2 min 20 s).
  - certifyn took 17 min wall time on 4 processes.
  - verifyn replays **410/410 verdicts OK**. `verifyn_u.py` is a copy patched only so it does not crash on undecided paths.
  - Result: **112 TWIN, 291 NO, 7 UNDECIDED**. An UNDECIDED pair has a μ(y, z) whose star orbit is exhausted before it reaches a decidable endpoint.
- `sieve.py n` builds `sieve<n>.json`. It records, per pair, every (i) and (ii) element with Des P_y, whether it is ⊋ A,
  whether P_y is an x-half of a candidate pair at (λ, Des P_y), the twin status of that pair with the same P_z,
  μ(y, x), μ(y, z), and l(y).
- `analyze.py n` runs tests 2 and 3 and the negative controls. Its output is in `out_analyze{14,15,16}.txt`.

## Test 1 (S_15): the (ii)-elements, in full

Only the 4 pairs of (5,4,3,3) at A = {2,3,4,6,7,8,10,11,12} and at its mirror {3,4,5,7,8,9,11,12,13} have (ii)-elements.
There are 2 distinct y, one per x. The other 44 pairs have an empty (ii)-set.

| pair | parity | verdict | Des P_y | cand x at Des P_y? | (P_y, P_z) status | μ(y,x) | μ(y,z) | l(y), l(x), l(z) |
|---|---|---|---|---|---|---|---|---|
| A, Pz = 1237b/456f/.. | fail | NO | A ∪ {14} | yes | TWIN | 1 | 0 | 52, 59, 76 |
| A, Pz = 1237f/456b/.. | ok | NO | A ∪ {14} | yes | NO | 1 | 0 | 52, 59, 77 |
| A′, Pz = 12348/.. | ok | NO | A′ ∪ {1} | yes | NO | 1 | 0 | 50, 49, 67 |
| A′, Pz = 1234c/.. | fail | NO | A′ ∪ {1} | yes | TWIN | 1 | 0 | 50, 49, 68 |

In the A′ rows, y = s_1·x is also an (i)-element.

**Structural remark (one line, follows from Lemma P).** μ(y, x) = 1 forces l(y) ≢ l(x) (mod 2). So for a parity-OK (x, z), the pair (P_y, P_z)
always fails parity, and for a parity-failing (x, z) it passes. The parity-failing z at A is the twin partner of y at A′.
In every (λ, A) with a (ii)-element, both z's die: one by parity, the other by y.

## Test 2: is every (ii)-element a candidate x-half at a larger A′?

| n | (ii) instances | distinct y | candidate x-half at (λ, Des P_y) | violating |
|---|---|---|---|---|
| 14 | 0 | 0 | — | — |
| 15 | 4 | 2 | 4 | **0** |
| 16 | 54 | 30 | 54 | **0** |

The claim survives, 58/58 instances, including 54 out of sample.

There is also a stronger observation. **Every decided (ii)-instance has μ(y, z) = 0**: 4/4 at n = 15 and 52/52 at n = 16.
The remaining 7 at n = 16 are undecided.

## Test 3: does a sieve rule predict twin status?

Rules, each stated in one sentence. Free choices are counted beyond the PROVE 281 definitions.

- **R0** (0 choices): twin ⟺ parity OK.
- **R1** (0 choices): twin ⟺ parity OK and the (ii)-set is empty.
- **R2**: twin ⟺ parity OK and every (ii)-y has μ(y, z) = 1. This is Lemma A restricted to (ii). It is not tableau-level.
- **R3** (1 choice, sieve form): twin ⟺ parity OK and no (ii)-y is a candidate x-half at Des P_y.
- **R4** (1 choice): twin ⟺ parity OK and no (ii)-y has (P_y, P_z) non-twin.

On the data R1 ≡ R2 ≡ R3 ≡ R4. R3 ≡ R1 holds by test 2. R4 ≡ R1 holds by the parity remark. R2 ≡ R1 holds because (ii) always kills.

Confusion matrices (TP/FP/FN/TN, TWIN = positive):

| rule | S_14 | S_15 | S_16 (403 decided) |
|---|---|---|---|
| R0 | 4/0/0/0 | 26/2/0/20 (46/48) | 112/105/0/186 (298/403) |
| R1 (= R2 = R3 = R4) | 4/0/0/0 | **26/0/0/22 (48/48)** | **112/86/0/205 (317/403)** |

**The sieve rule fits S_15 exactly and FAILS out of sample at S_16.** It gives 86 false positives and 0 false negatives.
So parity plus no-(ii) is necessary in all the data, but it is not sufficient.

**The 86 false positives are a new failure mode.** These pairs are parity-OK with an empty (ii)-set, and they are killed by an **(i)-element**: an ordinary star image s·x with μ(sx, z) = 0.

- There are 367 killing (i)-instances.
- Of these, 9 have Des ⊋ A, and almost none are candidate x-halves (`fp16.py`).
- The mechanism is mostly immediate Bruhat-incomparability in Q_0 (175 instances at path length 0). The rest are deeper incomparability or F1.
- This failure mode **does not occur at n ≤ 15**. At S_15, the only (i)-kill of a parity-OK pair is s_1·x in the mirror U case, which is also a (ii)-element.
- Two tableau predicates were tried for μ(sx, z). Neither separates (`probe_i.py`):
  - Des P_z ⊆ Des P_sx holds for all 894 cases, both kills and non-kills.
  - Des P_z ⊆ int Des P_sx gives 15 false positives and 527 false negatives at n = 16.

**Shape level, n = 16** (`out_shapes16.txt`). Four shapes carry (★)-violations but have **no twin**. All their pairs are decided NO:

| shape | pairs | NO by parity | NO by (i)-kill |
|---|---|---|---|
| (5,3,3,3,2) | 3 | 2 | 1 |
| (4,4,3,3,2) | 3 | 2 | 1 |
| (5,5,4,2) | 3 | 2 | 1 |
| (5,5,4,1,1) | 3 | 2 | 1 |

So the shape-level hunch "a shape has a twin ⟺ it carries a violation" (9/9 at n = 15) **fails at n = 16: 17/21**. The pairs of (5,5,3,3) and (6,4,3,3) include UNDECIDED ones, but both shapes do have twins.

## Controls

- **S_14.** All 4/4 pairs come out TWIN and every rule reproduces them. The (ii)-set is empty and every (i) has μ = 1. This is trivial agreement, not evidence.
- **Negative controls** (10^4 trials, seed 282). Each one breaks the claimed link; none merely relabels it.

| control | S_15 (rule exact = 48) | S_16 (rule = 317) |
|---|---|---|
| NC1: permute the (ii)-sets among pairs, rule R1 | mean 44.0, P(48) = 0.0008 | mean 297.2, max 317 |
| NC2: shuffle all μ(y, z) values, rule "parity & all (i)+(ii) μ(y,z)=1" (unshuffled = 48 and 403) | mean 25.1, max 31 | mean 291.7, max 298 |
| NC3: shuffle the parity labels, rule R1 | mean 26.7, max 39 | mean 219.6, max 258 |

Agreement drops in every case. NC1 at S_15 is weak: only 2 pairs are parity-OK and carry (ii), which is the small-sample caveat behind the n = 16 failure.

## Verdict

1. **"Every (ii)-element is a candidate x-half at a larger A′": SURVIVES**, 58/58 computed, 54 of them new at n = 16.
   "(ii)-elements always kill" also survives: 56/56 decided.
2. **"Twin ⟺ parity + no (ii)-element" (the top-down sieve): DEAD at n = 16** (86 false positives). The sieve is necessary in the data
   but not sufficient. A third obstruction, star-image incomparability, appears first at n = 16 and is not tied to larger descent sets.
3. **The shape-level hunch is dead at n = 16**: (5,3,3,3,2), (4,4,3,3,2), (5,5,4,2) and (5,5,4,1,1) carry violations but have no twin.
4. Open: the 7 UNDECIDED pairs. In 5 of them the only undecided μ belongs to a (ii)-element. If "(ii) always kills" holds, they are NO.
