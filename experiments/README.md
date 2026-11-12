# Fixed-cofactor CRT experiment

The subsequent P={2} polynomial-intersection experiment is documented in
[INTERSECTIONS.md](INTERSECTIONS.md). It searches directly for a pair with
exponents ell and ell+2.

`cofactor_experiment` searches for exact positive relations

```
A_j^2 - g_j N = c * product(p_i ^ e_ji)
```

over a fixed prime set P. It discovers the **first repeated canonical cofactor
c > 1**, then locks c and searches only for additional relations with that same
cofactor. The target defaults to |P|. It does not factor c or add c's factors to P.
All powers of every prime in P are removed before comparing cofactors.

This is a separate, bounded research executable. It does not change the existing
singleton-{2} generator or its random 512-bit demo.

## Build and run (PowerShell)

```powershell
cmake --preset vs2022-x64
cmake --build build --config Release --target cofactor_experiment crt_test
ctest --test-dir build -C Release --output-on-failure
./build/Release/cofactor_experiment.exe --csv build/fixed_cofactor.csv
```

All defaults are reproducible:

```
--number          1022117
--base-size       5
--relations       |P|
--lift-levels     3
--radius          32
--max-candidates  250000
--max-multiplier  1000000000
--csv             (optional output file)
```

`--self-test` runs the mathematical checks. `--help` shows the options. The
executable returns 0 if the relation target was reached, 2 for a bounded search
that did not reach its target, and 1 for an error. Reaching the target does not
imply finding a factor.

## Search performed

1. Choose g0 in {1,3,5,7} with g0*N = 1 (mod 8). Select the first |P| primes
   coprime to g0*N at which g0*N has a square root, including 2. Any factor of N
   noticed while selecting P is reported separately, not used by the search.
2. Vary the forced exponents independently. With three lift levels, the exponent
   of 2 ranges over 3,4,5, and each odd prime's exponent over 1,2,3. The first
   prime's exponent varies fastest. All selected primes are forced in every
   candidate; this is a restricted search, not enumeration of every P-smooth
   right-hand side.
3. Hensel-lift each prime-power root, retaining all four roots at powers of 2
   and both roots at the odd primes. CRT combines all branches into b mod M,
   where M is the product of the forced prime powers.
4. Search u near floor((sqrt(g0*N)-b)/M), in order 0,+1,-1,+2,-2,... relative
   to that center. Only A=b+M*u > 0 with A^2-g0*N > 0 are examined.
5. Evaluate both g=g0 and, when different, the multiplier g=g0+M*k that makes
   the positive raw quotient as small as possible subject to the g bound:

   ```
   F = (A^2-g0*N)/M
   k = min(floor((F-1)/N), floor((max_multiplier-g0)/M))
   A^2-(g0+M*k)*N = M*(F-k*N)
   ```

   This optimizes the raw quotient for a fixed A and branch. It does not
   minimize the canonical cofactor after trial division over every possible k.
6. Remove all prime powers in P, including any beyond the forced exponents.
   Skip c=1 to study a common external cofactor. Deduplicate by
   `(c, min(A mod N, -A mod N), full exponent vector)`. This removes repeated
   identical modular relations and their sign equivalents; it does not claim
   that all remaining relations are independent.
7. Lock onto the first repeated c>1 and collect up to the requested count.
   A lock is never silently replaced by a more successful cofactor. Stop on
   reaching the target, exhausting the profile space, or the candidate budget.

Search sizes are intentionally capped: 2..12 factor-base primes, 1..8 lift
levels, radius at most 10000, and 2..1024 collected relations. CRT branch counts
and exponent-profile counts grow exponentially in the base size. The use of GMP
permits large N, but this is not a claim of efficient large-semiprime factoring.

## Dependency checks

The parity vector contains one coordinate per prime. If c is nonsquare, it also
contains a coordinate equal to 1: an even number of rows must be multiplied for
c to become a square. If c is a square, this extra coordinate is zero.

Thus |P| relations alone do not guarantee a dependency. With nonsquare c, the
general dimension bound guarantees one by |P|+2 rows (possibly sooner), but a
dependency can still give only trivial gcds. Use `--relations 7`, for example,
to collect more than five rows without changing P.

Gaussian elimination records the exact row subsets. Each tested subset is
checked by reconstructing X and Y modulo N and verifying X^2=Y^2 (mod N), then
computing gcd(X-Y,N) and gcd(X+Y,N). All kernel basis vectors are tested. If
needed, combinations of the first twelve are also tested; a failure is not an
exhaustive claim about a larger kernel. Direct gcd(c,N) factors are reported
separately.

Tests compare Hensel and CRT roots against exhaustive small-modulus enumeration,
check the nonsquare-cofactor parity coordinate, distinguish a useful pair from
a merely scaled seed, and verify every certificate in the default search.

## Observations on 2026-11-12

These are small exploratory runs, not a benchmark or complexity result. No
factorization of N was supplied to the search.

| N | P | Locked c | Relations / target | Candidates | Dependency result |
|---|---|---:|---:|---:|---|
| 1022117 | 2,3,7,13,23 | 1343 | 5 / 5 | 26923 | 1009 * 1013 |
| 1248133 | 2,7,17,19,23 | 4245 | 2 / 5 | 1000000 | no dependency |
| 100160063 | 2,5,17,19,37 | 9307 | 2 / 5 | 1000000 | no dependency |
| 1000036000099 | 2,7,11,17,19 | 11023 | 2 / 5 | 1000000 | no dependency |

The first row uses default parameters. The others use `--lift-levels 4
--max-candidates 1000000`; their 250000-candidate, three-level runs also stopped
at two relations with the same respective cofactors. CSV certificates are in
`results/`. Their `forced_exponents` describe M; `actual_exponents` describe the
full relation. Lists use semicolons, matching the `primes` column.

The successful example has exponent vectors, in prime order 2,3,7,13,23:

| A | g | Actual exponents |
|---:|---:|---|
| 555095 | 301397 | 3,1,1,1,1 |
| 1110335 | 1205573 | 3,3,1,1,1 |
| 610433 | 351629 | 5,1,3,1,1 |
| 394529 | 150701 | 6,2,1,1,1 |
| 7611791 | 56511005 | 4,4,3,1,1 |

Every row has cofactor 1343. The parity rank is 2, so these five rows provide
three independent kernel vectors. The first two rows already yield
X=767357 and Y=5562 modulo N, with gcd(X-Y,N)=1009 and gcd(X+Y,N)=1013.
The first repeated cofactor was found at candidate 398; collecting five rows
required 26923 candidates.

The other runs show the central unresolved issue: finding a repeated cofactor
does not guarantee that a fixed-c group will fill. Large permitted multipliers
also contribute to the successful example; it is not a demonstration with
uniformly small g. The CSV files make both limitations inspectable.
