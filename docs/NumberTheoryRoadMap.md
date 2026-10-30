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
- [x] Decide the integer types to support: `long` everywhere, `BigInteger` overloads where it matters
      (`System.Numerics.BigInteger` is BCL on all three targets, so no external dependency)
- [x] Settle overflow strategy for `long` modular multiplication (`UInt128` on net8+/net10,
      mul-mod fallback on netstandard2.1)
- [x] Settle API shape: static classes (`ModularArithmetic.Pow(a, e, m)`) plus thin extension methods
      for the common cases (`n.IsPrime()`, `a.Gcd(b)`)
- [x] Define `ModInt` / `FiniteFieldElement` value types and their operators
- [x] Plan backward compatibility for `NumericExtensions.IsPrime(int)` and `GetPrimeFactors(int)`
      (keep signatures, delegate to the new implementation)

See [Design decisions](#design-decisions) below.

## Phase 2 — Euclidean Algorithm & Modular Arithmetic
- [x] `EuclideanAlgorithm`: `Gcd`, `Lcm`, `ExtendedGcd` (returns g, x, y), binary GCD, GCD/LCM of sequences
- [x] `ModularArithmetic`: `Mod` (always non-negative), `AddMod`, `MulMod` (overflow-safe), `Pow`
      (square-and-multiply), `Inverse` (throws when gcd ≠ 1), `TryInverse`
- [x] `ModularArithmetic`: Euler's totient φ(n), Carmichael λ(n), multiplicative order, primitive root
- [x] `ModularArithmetic`: Legendre and Jacobi symbols, modular square root (Tonelli–Shanks)
- [x] `ChineseRemainderTheorem`: solve for pairwise coprime moduli, plus the general (non-coprime)
      case that reports when no solution exists
- [x] Replace the private `GCD` in `Physics/Quantum/ModularMultiplyGate.cs` with `EuclideanAlgorithm.Gcd`
- [x] Unit tests in `NumberTheoryTests` (identities: Bézout, Fermat's little theorem, CRT round-trip)

## Phase 3 — Primality Testing & Factorisation
- [ ] `PrimalityTesting`: deterministic Miller–Rabin for `long` (fixed witness set valid below 2⁶⁴)
- [ ] `PrimalityTesting`: probabilistic Miller–Rabin for `BigInteger` with configurable rounds/seed
- [ ] `PrimeSieve`: sieve of Eratosthenes (`PrimesUpTo(n)`), segmented sieve for ranges
- [ ] `NextPrime` / `PreviousPrime`
- [ ] `IntegerFactorization`: trial division for small factors + Pollard's rho (Brent variant),
      result as `(prime, exponent)` pairs
- [ ] Divisor functions: divisor list, σ(n), τ(n), Möbius μ(n)
- [ ] Switch `Totient`, `CarmichaelLambda`, `MultiplicativeOrder` and `PrimitiveRoot` from the interim
      private trial-division `Factor` in `ModularArithmetic` to `IntegerFactorization`
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

## Design decisions

Outcome of Phase 1. Later phases implement against this.

### Location and files

Namespace `CSharpNumerics.Numerics.NumberTheory`, folder `Numerics/Numerics/Numerics/NumberTheory/`.

| File | Kind | Phase |
|------|------|-------|
| `EuclideanAlgorithm.cs` | static class | 2 |
| `ModularArithmetic.cs` | static class | 2 |
| `ChineseRemainderTheorem.cs` | static class | 2 |
| `ModInt.cs` | `readonly struct` | 2 |
| `NumberTheoryExtensions.cs` | static extension class | 2 (grows in 3) |
| `PrimalityTesting.cs`, `PrimeSieve.cs`, `IntegerFactorization.cs` | static classes | 3 |
| `ContinuedFractions.cs` | static class | 4 |
| `PrimeField.cs`, `Polynomial` over GF(p), `ExtensionField.cs`, `ExtensionFieldElement.cs` | classes / struct | 5 |
| `Interfaces/IFiniteField.cs` | interface | 5 |

Tests go in `NumericTest/NumberTheoryTests.cs`.

### Integer types

- **`long` is the primary type.** There are no `int` overloads: `int` converts implicitly to `long`,
  so `ModularArithmetic.Pow(3, 4, 7)` works as is. The exception is the existing
  `NumericExtensions.IsPrime(this int)`, which stays (see backward compatibility).
- **`BigInteger` overloads** for `Gcd`, `ExtendedGcd`, `Lcm`, `Pow`, `Inverse`, CRT and
  `IsProbablePrime`. Where the BCL already has an implementation (`BigInteger.GreatestCommonDivisor`,
  `BigInteger.ModPow`) we wrap it instead of reimplementing it, but normalise the results the same way
  as the `long` versions (see conventions).
- Verified: `BigInteger` and value tuples compile on `netstandard2.1` without extra package references.

### Overflow strategy for `MulMod(a, b, m)` on `long`

The inputs are first reduced to `[0, m)`. Then:

1. **Fast path:** if `m ≤ 3 037 000 499` (⌊√`long.MaxValue`⌋), compute `a * b % m` directly. This
   cannot overflow and covers most practical use.
2. **net8.0 / net10.0 (`#if NET7_0_OR_GREATER`):** `(long)((UInt128)(ulong)a * (ulong)b % (ulong)m)`.
3. **netstandard2.1:** double-and-add on `ulong`, which needs no allocation. It cannot overflow
   because `m ≤ 2⁶³ − 1`, so `x + x < 2⁶⁴`. `Math.BigMul(ulong, ulong, out ulong)` and `UInt128` do
   not exist on netstandard2.1.

Spike results (net10.0, 200 000 random cases and edge cases up to `long.MaxValue`): both strategies
agree with `BigInteger` with no mismatches. Over 2M multiplications near 2⁶³, `UInt128` took about
4 ms, double-and-add about 200 ms, and `BigInteger` about 280 ms with allocations, so double-and-add
is the better fallback.

This is the **first `#if` in the library**. It stays confined to the private `MulMod` helper in
`ModularArithmetic`, and every other method calls that helper.

### Conventions and error handling

- Results of modular operations are always in `[0, m)`, including for negative inputs (`Mod(-3, 7) = 4`).
- `m ≤ 0` throws `ArgumentOutOfRangeException`. `m = 1` is allowed and gives 0. Negative exponents in
  `Pow` mean "inverse first" and throw if no inverse exists.
- `Inverse` throws `ArgumentException` when gcd(a, m) ≠ 1, consistent with `ModularMultiplyGate`.
  `TryInverse(a, m, out inverse)` is the version that does not throw.
- `Gcd` is always ≥ 0, and `Gcd(0, 0) = 0`. Inputs whose absolute value cannot be represented
  (`long.MinValue`) throw `OverflowException`. `Lcm` throws `OverflowException` on overflow (`checked`).
- `ExtendedGcd` returns `(long Gcd, long X, long Y)` with `a·X + b·Y = Gcd`.
- `ChineseRemainderTheorem.Solve(remainders, moduli)` returns `(long Remainder, long Modulus)`. It
  throws if the moduli are not pairwise coprime. `TrySolve` handles the general case, where the
  moduli need not be coprime, and returns `false` when the system is inconsistent.

### API shape

Static classes carry the full API. A small `NumberTheoryExtensions` adds the common cases:

```csharp
using CSharpNumerics.Numerics.NumberTheory;

long g   = EuclideanAlgorithm.Gcd(84, 36);               // 12
var (d, x, y) = EuclideanAlgorithm.ExtendedGcd(240, 46); // 2, -9, 47
long p   = ModularArithmetic.Pow(4, 13, 497);            // 445
long inv = ModularArithmetic.Inverse(3, 11);             // 4
var (r, m) = ChineseRemainderTheorem.Solve(new long[] { 2, 3, 2 }, new long[] { 3, 5, 7 }); // (23, 105)

84L.Gcd(36);          // extensions: Gcd, Lcm, ModPow, ModInverse, IsPrime (long, Phase 3)
```

Extensions are added only for operations that read naturally as `a.Op(b)`. Everything else stays on
the static classes.

### `ModInt` (Phase 2)

`public readonly struct ModInt : IEquatable<ModInt>`

- Fields: `Value` (always in `[0, Modulus)`) and `Modulus`.
- Constructor `ModInt(long value, long modulus)` normalises `value`.
- Operators `+ - * /` and unary `-`. `/` multiplies by the inverse, and throws if there is none.
  Operations that mix moduli throw `ArgumentException`.
- `Pow(long e)`, `Inverse()`, `==`, `!=`, `Equals`, `GetHashCode`, and `ToString()` → `"3 (mod 7)"`.
- No implicit conversions, so a modulus is never chosen by accident. `ModInt` with a non-prime
  modulus is a ring element and `/` works only for units. With a prime modulus it is a GF(p) element.
- No generic math (`INumber<T>`), because it does not exist on netstandard2.1.

### Finite fields (Phase 5)

- `IFiniteField<TElement>` in `Numerics/Interfaces/`: `Characteristic`, `Order`, `Zero`, `One`,
  `Add`, `Subtract`, `Multiply`, `Inverse`, `Pow`, `Elements()`.
- `PrimeField : IFiniteField<ModInt>` wraps a prime p, validated with `PrimalityTesting`.
  `field.Element(3)` creates a `ModInt`.
- `ExtensionField : IFiniteField<ExtensionFieldElement>` for GF(pⁿ), defined by p and an irreducible
  polynomial that is either given or found automatically.
- `ExtensionFieldElement` is a `readonly struct` that holds a reference to its field and an
  immutable coefficient array, with operators `+ - * /`. Mixing fields throws.
- GF(2⁸) gets precomputed log/antilog tables inside `ExtensionField` when p = 2 and n ≤ 16. The
  public API is the same.

### Backward compatibility

- `NumericExtensions.IsPrime(this int)` and `GetPrimeFactors(this int)` keep their signatures and
  namespace (`CSharpNumerics.Numerics`). In Phase 3 they become one-line delegates to
  `PrimalityTesting.IsPrime(long)` and `IntegerFactorization`.
- `GetPrimeFactors` keeps its exact output format: ascending primes listed with repetition
  (78 → `[2, 3, 13]`, 12 → `[2, 2, 3]`). An empty list for n < 2 is the current behaviour and stays.
  The existing tests `TestIsPrime` and `TestGetPrimeFactors` must pass unchanged.
- `IsPrime(this long)` is added in `NumberTheoryExtensions`. It cannot change what `79.IsPrime()` binds
  to: C# does not apply numeric conversions to the receiver of an extension method, so an `int`
  receiver only ever matches the `int` overload. The flip side is that the `long` extensions need an
  `L` suffix or a cast (`84L.Gcd(36)`), which the README must show.

## Notes
- Phase 2: the multiplicative-group functions (`Totient`, `CarmichaelLambda`, `MultiplicativeOrder`,
  `PrimitiveRoot`) factor by trial division until Phase 3. They are slow when the argument has a prime
  factor above about 10¹², so the tests use moduli well below that.
- The test project only targets net10.0, so the netstandard2.1 double-and-add branch of `MulMod`
  is not exercised by the tests. The same algorithm was checked against `BigInteger` in the
  Phase 1 spike.
- Placement follows the boundary rules: this is pure math with no domain semantics, so it lives in
  Numerics. Physics (Quantum) and Statistics (random generators) consume it.
- `EuclideanAlgorithm` and `ChineseRemainderTheorem` are small. They are kept as separate static
  classes to respect the one public class per file rule, rather than merged into `ModularArithmetic`.
