# Equal-cofactor polynomial intersections over P={2}

`cofactor_intersections` searches the proposed equation directly with
`ell2 = ell1 + 2`:

```
C1(u1) = 2^ell1*u1^2 + 2*b1*u1 + (b1^2-g1*N)/2^ell1
C2(u2) = 2^ell2*u2^2 + 2*b2*u2 + (b2^2-g2*N)/2^ell2
C1(u1) = C2(u2) = c > 0, c odd
```

The multipliers are held fixed during a run. By default both are the unique g
in {1,3,5,7} with g*N=1 modulo 8. This is a separate executable from the earlier
fixed-cofactor grouping experiment. It discovers a common c per intersection;
after a trivial intersection it continues looking for a useful pair, rather
than locking itself to that trivial pair's c.

## Reproduce

```powershell
cmake --preset vs2022-x64
cmake --build build --config Release --target cofactor_intersections
ctest --test-dir build -C Release --output-on-failure
./build/Release/cofactor_intersections.exe
./build/Release/cofactor_intersections.exe --number 1022117
./build/Release/cofactor_intersections.exe --number 1000036000099 --max-candidates 1000000
```

Options:

```
--number          77
--g1              automatic small multiplier
--g2              automatic small multiplier
--ell-min         3
--ell-max         12
--max-candidates  100000
--keep-searching  optional: continue after a useful match
--output          optional text report path
--self-test       run mathematical invariant checks
```

Require odd N>1, both g_i*N=1 modulo 8, and 3<=ell-min<=ell-max<=256.
The executable exits with code 0 on finding a nontrivial factor, 2 when the
budget finishes without one, and 1 on an error. Reports retain all encountered
intersections, their exact integer certificates, and gcd outcomes. In the
printed `gcd_minus` and `gcd_plus` fields, 0 denotes a trivial gcd (1 or N).

## How the search works

1. For every ell1 in the configured interval, lift all four roots b1 modulo
   M1=2^ell1. Keep only those branches whose quotient C1 is odd. There are two:
   the parity of C1 is constant along a branch because its increment is
   `C1(u+1)-C1(u)=2*A1+M1`, which is even. The other branches have strictly
   larger 2-adic valuation and belong to a different exponent.
2. Initialize each branch at its smallest positive A1=b1+M1*u1 strictly above
   sqrt(g1*N). Merge these progressions in increasing A1 using a small heap.
   This fairly interleaves exponent choices instead of spending the entire
   budget at one ell1. A positive A1 has only one exact valuation, so accepted
   candidates cannot be duplicated across exponents.
3. For each candidate with c=C1(u1), equality forces

   ```
   A2^2 = g2*N + 2^(ell1+2)*c
        = 4*A1^2 + (g2-4*g1)*N.
   ```

   Test this single integer for being a perfect square. This replaces a
   Cartesian search of (u1,u2) with one square-root test per candidate. No
   factorization of N, c, or (g2-4*g1)*N is supplied or computed to find the
   intersections.
4. For a square, set A2 to its positive root, b2=A2 mod 2^ell2 and
   u2=(A2-b2)/2^ell2. The exact equation and odd c automatically imply that b2
   is a lifted root with exact valuation ell2. There is no need to enumerate
   u2 or precompute a table of C2 values.
5. Independently verify both polynomial values and

   ```
   (A2-2*A1)*(A2+2*A1) = (g2-4*g1)*N.
   ```

   Test gcd(A2-2*A1,N) and gcd(A2+2*A1,N). Continue past trivial matches;
   normally stop at the first nontrivial factor.

Only positive odd cofactors and positive square representatives are searched.
Negative cofactors, other exponent gaps, other g choices, and exponents outside
the configured range are not covered. Exhausting the candidate budget does not
prove that no suitable pair exists.

## Results on 2026-11-12

All rows below use ell1=3..12 and the same small multiplier on both sides.
Factors were not inputs to the search. Counts are square-root tests on exact-
valuation candidates, not counts directly comparable to the previous program's
CRT/multiplier candidate count.

| N | g | Candidates | ell1,ell2 | Common c | Result |
|---|---:|---:|---|---:|---|
| 77 | 5 | 5 | 3,5 | 57 | 7 * 11 |
| 1022117 | 5 | 888 | 5,7 | 349587 | 1009 * 1013 |
| 1248133 | 5 | 797 | 4,6 | 656999 | 1009 * 1237 |
| 100160063 | 7 | 14267 | 3,5 | 291021385 | 10007 * 10009 |
| 1000036000099 | 3 | 383591 | 3,5 | 406245874991 | 1000003 * 1000033 |
| 2000009000009 | 1 | 167730 | 5,7 | 33203171875 | 1000003 * 2000003 |
| 1000000016000000063 | 7 | 1000000 (budget) | - | - | no intersection found |

The 40-bit N=1000036000099 run initially found no match within 100000 candidates;
the extended run reached the pair shown above. Both reports are retained.
All reports are in `results/intersection_*.txt`.

For N=77 the search first finds a trivial pair with c=9, then the useful pair:

```
29^2 - 5*77 = 2^3 * 57
47^2 - 5*77 = 2^5 * 57
gcd(47-2*29,77)=11, gcd(47+2*29,77)=7
```

This pair occurs before the original illustrative pair (43,79,c=183) because
both exact-valuation root branches are searched. The original example is also
checked in the tests.

The earlier failed fixed-cofactor search at N=100160063 is now replaced by:

```
55039^2 - 7*100160063 = 2^3 * 291021385
100069^2 - 7*100160063 = 2^5 * 291021385
```

Here the algorithm finds a useful pair with a new cofactor, not five relations
with the old locked cofactor. This is a different target from the previous
experiment and is the reason the comparison is exploratory.

An additional diagnostic with N=77, g1=5, g2=13 finds (A1,A2,c)=(21,35,7).
That cofactor already shares a factor with N, so this is a correctness check
for unequal multipliers rather than strong factoring evidence.

## Interpretation and validation

When g1=g2=g the equation is

```
(2*A1-A2)*(2*A1+A2) = 3*g*N.
```

Thus this implementation is a constrained difference-of-squares search on
3*g*N, with exact 2-adic valuations and c>0. It scans A1 linearly and does not
provide a polynomial-time factoring algorithm. There is no demonstrated
asymptotic improvement over existing factoring methods. Several sample inputs
have close factors and are also easy for ordinary Fermat factorization; no
performance advantage over Fermat is asserted.

The tests compare lifted roots against exhaustive enumeration, reproduce both
useful N=77 examples, reject a trivial congruence, enforce the candidate budget,
and compare all matches in a bounded heap search against an independent plain
odd-integer scan for equal and unequal multipliers. Every recorded match is
also checked against its exact integer identities before it is reported.

## Twenty local OpenMP experiments with separated factors

The executable now supports a deterministic batch of 20 semiprimes. OpenMP
parallelizes complete independent searches with dynamic scheduling; each worker
owns its GMP state and heap. It does not change the per-case candidate order or
budget. The additional `cofactor_openmp_reproducibility` CTest compares serial
and parallel results and verifies the batch input construction.

```powershell
cmake --preset vs2022-x64
cmake --build build --config Release --target cofactor_intersections
./build/Release/cofactor_intersections.exe --batch-separated --max-candidates 5000000 --batch-directory experiments/results/separated_openmp
```

`--threads 0` (the default) uses up to the runtime's available processor/thread
count, capped at the 20 independent jobs. An explicit `--threads 1..20` overrides
it. CMake discovers and links OpenMP for this executable.

The completed local run used **14 actual OpenMP workers**. It tested five
smaller-prime scales from approximately 10^3 through 10^7, each with larger-prime
ratios approximately 3, 10, 30, and 100. Products span 22..54 bits. Both factors
were verified by exact trial division; only their product was passed to search.

| Approximate q/p | Factored within 5 million candidates |
|---:|---:|
| 3 | 5 / 5 |
| 10 | 5 / 5 |
| 30 | 4 / 5 |
| 100 | 3 / 5 |

In total **17/20** cases factored, with **25,964,989** candidates examined and no
reported errors. The three cases exhausting the budget were:

```
100090653535849   = 1000453  * 100045333
3000322878687049  = 10000537 * 300016177
10001166324005807 = 10000583 * 1000058329
```

All 20 input factorizations and all 17 intersection certificates were checked
independently using PowerShell/.NET BigInteger arithmetic. All four CTest suites
passed before the batch. See [full results](results/separated_openmp/summary.md)
and [independent verification](results/separated_openmp/verification.txt).

The elapsed-time records show a discontinuity: fifteen cases took under 1.5
seconds each, while five cases show approximately 7300 seconds. The reports
preserve those observed wall-clock values. The cause was not measured; these
times should not be treated as a clean CPU benchmark or evidence of OpenMP
speedup. Candidate counts and verified certificates remain the reproducible
results. Factor separation and product size both vary in this structured sample,
so the success rates do not isolate the effect of separation alone.
