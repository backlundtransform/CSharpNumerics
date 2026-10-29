# Number Theory Roadmap

Add a `CSharpNumerics.Numerics.NumberTheory` sub-section (`Numerics/Numerics/NumberTheory/`) with
integer and modular arithmetic, primality testing, factorisation, continued fractions and finite
fields. Today this functionality is scattered and weak: `NumericExtensions.IsPrime` does O(n) trial
division, `GetPrimeFactors` is plain trial division, and `ModularMultiplyGate` carries its own
private `GCD`. Consolidating it in Numerics gives Quantum (Shor), error correction, transforms and
random generators one shared foundation.

Scope is numerical and educational use. The implementation is **not** constant-time and must be
documented as not suitable for cryptographic purposes.

## Phase 1 — Research & Design
- [ ] Decide the integer types to support: `long` everywhere, `BigInteger` overloads where it matters
      (`System.Numerics.BigInteger` is BCL on all three targets, so no external dependency)
- [ ] Settle overflow strategy for `long` modular multiplication (`UInt128` on net8+/net10,
      `Math.BigMul`/mul-mod fallback on netstandard2.1)
- [ ] Settle API shape: static classes (`ModularArithmetic.Pow(a, e, m)`) plus thin extension methods
      for the common cases (`n.IsPrime()`, `a.Gcd(b)`)
- [ ] Define `ModInt` / `FiniteFieldElement` value types and their operators
- [ ] Plan backward compatibility for `NumericExtensions.IsPrime(int)` and `GetPrimeFactors(int)`
      (keep signatures, delegate to the new implementation)

## Phase 2 — Euclidean Algorithm & Modular Arithmetic
- [ ] `EuclideanAlgorithm`: `Gcd`, `Lcm`, `ExtendedGcd` (returns g, x, y), binary GCD, GCD/LCM of sequences
- [ ] `ModularArithmetic`: `Mod` (always non-negative), `AddMod`, `MulMod` (overflow-safe), `Pow`
      (square-and-multiply), `Inverse` (throws when gcd ≠ 1), `TryInverse`
- [ ] `ModularArithmetic`: Euler's totient φ(n), Carmichael λ(n), multiplicative order, primitive root
- [ ] `ModularArithmetic`: Legendre and Jacobi symbols, modular square root (Tonelli–Shanks)
- [ ] `ChineseRemainderTheorem`: solve for pairwise coprime moduli, plus the general (non-coprime)
      case that reports when no solution exists
- [ ] Replace the private `GCD` in `Physics/Quantum/ModularMultiplyGate.cs` with `EuclideanAlgorithm.Gcd`
- [ ] Unit tests in `NumberTheoryTests` (identities: Bézout, Fermat's little theorem, CRT round-trip)

## Phase 3 — Primality Testing & Factorisation
- [ ] `PrimalityTesting`: deterministic Miller–Rabin for `long` (fixed witness set valid below 2⁶⁴)
- [ ] `PrimalityTesting`: probabilistic Miller–Rabin for `BigInteger` with configurable rounds/seed
- [ ] `PrimeSieve`: sieve of Eratosthenes (`PrimesUpTo(n)`), segmented sieve for ranges
- [ ] `NextPrime` / `PreviousPrime`
- [ ] `IntegerFactorization`: trial division for small factors + Pollard's rho (Brent variant),
      result as `(prime, exponent)` pairs
- [ ] Divisor functions: divisor list, σ(n), τ(n), Möbius μ(n)
- [ ] Make `NumericExtensions.IsPrime` / `GetPrimeFactors` delegate to the new code
      (existing tests in `NumericsTests` must keep passing)
- [ ] Unit tests (Carmichael numbers, strong pseudoprimes, large semiprimes, 2⁶¹−1)

## Phase 4 — Continued Fractions & Shor Integration
- [ ] `ContinuedFractions`: expansion of a rational / double, convergents, best rational approximation
      with bounded denominator
- [ ] Shor classical post-processing: from a measured phase, recover the period r via convergents,
      then factor N from gcd(a^(r/2) ± 1, N)
- [ ] End-to-end test: factor N = 15 and N = 21 using the existing `ModularMultiplyGate` and the
      quantum phase-estimation building blocks
- [ ] Unit tests for continued fractions (√2, π, known convergents)

## Phase 5 — Finite Fields
- [ ] `PrimeField` GF(p) and `ModInt` element type with field operators, inverse and power
- [ ] Polynomial arithmetic over GF(p): add, multiply, divmod, gcd, evaluation
- [ ] Irreducibility test (Rabin) and search for irreducible / primitive polynomials
- [ ] `ExtensionField` GF(pⁿ) with element type; table-based fast path for GF(2⁸)
- [ ] Unit tests (field axioms, multiplicative group order pⁿ−1, AES GF(2⁸) known products)

## Phase 6 — Integration & Docs
- [ ] Number Theoretic Transform (NTT) and exact integer convolution, next to the existing FFT in
      `TransformExtensions`
- [ ] Evaluate a classical Reed–Solomon or BCH code over GF(2⁸) as a follow-up; note it in the
      roadmap backlog if out of scope
- [ ] Update `Numerics/Numerics/Numerics/README.md` with a Number Theory section (example-driven)
- [ ] Add the "not cryptographically secure" note to the README and XML docs

## Notes
- Placement follows the boundary rules: this is pure math with no domain semantics, so it lives in
  Numerics. Physics (Quantum) and Statistics (random generators) consume it.
- `EuclideanAlgorithm` and `ChineseRemainderTheorem` are small. They are kept as separate static
  classes to respect the one public class per file rule, rather than merged into `ModularArithmetic`.
