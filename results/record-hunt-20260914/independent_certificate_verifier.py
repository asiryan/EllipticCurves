#!/usr/bin/env python3
"""Exact, standalone rank lower bounds from local 2-Kummer characters.

No Sage, PARI, floating point heights, BSD, GRH, or rank metadata is used.
Input: a JSON object with ``ainvs`` and ``points`` (rational strings).

For the integral general Weierstrass equation, set X=4*x and
V=4*(2*y+a1*x+a3).  Then
    V^2 = f(X) = X^3 + b2*X^2 + 8*b4*X + 16*b6.
At an odd good prime and each root alpha of f in F_p, the standard
2-Kummer character is the quadratic character of X-alpha.  At the
2-torsion point (alpha,0), its value is the character of f'(alpha).
These are homomorphisms E(Q) -> F_2 and kill 2*E(Q).  A matrix of
characters of k points with rank s proves
    rank E(Q) >= s - dim_F2 E(Q)[2].
The number of roots of f modulo a good odd prime gives an upper bound
on dim E(Q)[2], so this checker remains conservative with 2-torsion.
In particular a rootless good reduction proves E(Q)[2]=0.  Full column
rank then proves all input points independent modulo rational torsion.

The rows are Kummer characters across several reductions, NOT a claim
that independence of points in one finite group proves Q-independence.
Failure to reach full rank is inconclusive: points may generate a
subgroup of even index, or additional primes may be needed.

References: J. E. Cremona, "On the computation of Mordell-Weil and
2-Selmer Groups of Elliptic Curves", section 2,
https://johncremona.github.io/papers/filter.pdf
and the leaderboard's implementation (independently coded here):
https://github.com/icarm/elliptic-rank/blob/main/src/verify.ts
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from math import isqrt
from pathlib import Path
import time


def primes_up_to(bound):
    sieve = bytearray(b"\1") * (bound + 1)
    if bound >= 0:
        sieve[0] = 0
    if bound >= 1:
        sieve[1] = 0
    for q in range(2, isqrt(bound) + 1):
        if sieve[q]:
            sieve[q*q:bound+1:q] = b"\0" * ((bound-q*q)//q + 1)
    return [q for q in range(2, bound + 1) if sieve[q]]


def is_prime(n):
    return n >= 2 and all(n % d for d in range(2, isqrt(n) + 1))


def normalized_input(data):
    if len(data["ainvs"]) != 5:
        raise ValueError("Five Weierstrass coefficients are required")
    ainvs = [Fraction(str(a)) for a in data["ainvs"]]
    if any(a.denominator != 1 for a in ainvs):
        raise ValueError("This implementation requires integral ainvs")
    ainvs = [int(a) for a in ainvs]
    points = []
    for xy in data["points"]:
        if len(xy) != 2:
            raise ValueError("Only finite affine points (x,y) are supported")
        points.append(tuple(Fraction(str(z)) for z in xy))
    return ainvs, points


def validate_curve_and_points(data):
    ainvs, points = normalized_input(data)
    a1, a2, a3, a4, a6 = ainvs
    b2 = a1*a1 + 4*a2
    b4 = 2*a4 + a1*a3
    b6 = a3*a3 + 4*a6
    b8 = a1*a1*a6 + 4*a2*a6 - a1*a3*a4 + a2*a3*a3 - a4*a4
    discriminant = -b2*b2*b8 - 8*b4**3 - 27*b6*b6 + 9*b2*b4*b6
    if discriminant == 0:
        raise ValueError("Singular curve")
    for j, (x, y) in enumerate(points):
        if y*y + a1*x*y + a3*y != x**3 + a2*x*x + a4*x + a6:
            raise ValueError(f"Point {j} does not satisfy the equation exactly")
    # Polynomial discriminant is 256*discriminant, so odd good primes
    # are identical on this model.
    return ainvs, points, (b2, 8*b4, 16*b6), discriminant


def fingerprint(ainvs, points):
    normalized = {"ainvs": [str(a) for a in ainvs],
                  "points": [[str(x), str(y)] for x, y in points]}
    payload = json.dumps(normalized, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode()).hexdigest()


def roots_mod(cubic, p):
    a, b, c = (v % p for v in cubic)
    return [x for x in range(p) if ((x + a)*x*x + b*x + c) % p == 0]


def character_row(points, cubic, p, alpha):
    """A bit row of a genuine homomorphism; reject denominator primes."""
    a, b, _ = cubic
    derivative = (3*alpha*alpha + 2*a*alpha + b) % p
    if derivative == 0:
        raise ValueError("The selected root is not simple")
    mask = 0
    for j, (x, _) in enumerate(points):
        X = 4*x
        if X.denominator % p == 0:
            raise ValueError("Prime divides a point denominator")
        residue = (X.numerator * pow(X.denominator, -1, p) - alpha) % p
        if residue == 0:
            residue = derivative
        symbol = pow(residue, (p-1)//2, p)
        if symbol not in (1, p-1):
            raise ValueError("Invalid quadratic character")
        if symbol == p-1:
            mask |= 1 << j
    return mask


def add_row(basis, mask):
    while mask:
        pivot = mask.bit_length() - 1
        if pivot in basis:
            mask ^= basis[pivot]
        else:
            basis[pivot] = mask
            return True
    return False


def bitstring(mask, columns):
    # Leftmost bit refers to input point zero.
    return "".join(str((mask >> j) & 1) for j in range(columns))


def build_certificate(data, max_prime=10000):
    ainvs, points, cubic, disc = validate_curve_and_points(data)
    basis, rows = {}, []
    torsion_bound = 2
    torsion_witness = None
    tested_primes = 0
    largest_prime = None
    n = len(points)
    for p in primes_up_to(max_prime):
        if p == 2 or disc % p == 0:
            continue
        tested_primes += 1
        largest_prime = p
        roots = roots_mod(cubic, p)
        upper_bound = {0: 0, 1: 1, 3: 2}[len(roots)]
        if torsion_witness is None or upper_bound < torsion_bound:
            torsion_bound = upper_bound
            torsion_witness = {"prime": p, "roots": roots,
                               "dimension_upper_bound": torsion_bound}
        if not any(x.denominator % p == 0 for x, _ in points):
            for alpha in roots:
                mask = character_row(points, cubic, p, alpha)
                if add_row(basis, mask):
                    rows.append({"prime": p, "root": alpha,
                                 "bits": bitstring(mask, n)})
        if len(basis) == n and torsion_bound == 0:
            break
    if torsion_witness is None:
        raise ValueError("No good odd prime found; increase --max-prime")
    return {
        "format": "local-2-kummer-rank-certificate-v1",
        "input_sha256": fingerprint(ainvs, points),
        "points_checked": n,
        "matrix_rank": len(basis),
        "torsion_2_dimension_upper_bound": torsion_bound,
        "rank_lower_bound": max(0, len(basis) - torsion_bound),
        "all_points_independent_modulo_torsion": len(basis) == n and torsion_bound == 0,
        "torsion_witness": torsion_witness,
        "independent_rows": rows,
        "good_primes_examined": tested_primes,
        "largest_prime_examined": largest_prime,
        "requested_prime_bound": max_prime,
        "conditional_assumptions": [],
    }


def verify_certificate(data, certificate):
    """Replay all mathematics; do not trust rank, roots, or metadata."""
    ainvs, points, cubic, disc = validate_curve_and_points(data)
    n = len(points)
    if certificate.get("format") != "local-2-kummer-rank-certificate-v1":
        raise ValueError("Unsupported certificate format")
    if fingerprint(ainvs, points) != certificate["input_sha256"]:
        raise ValueError("Certificate belongs to different inputs")
    torsion = certificate["torsion_witness"]
    p = torsion["prime"]
    if not is_prime(p) or p == 2 or disc % p == 0:
        raise ValueError("Invalid prime in torsion witness")
    roots = roots_mod(cubic, p)
    bound = {0: 0, 1: 1, 3: 2}[len(roots)]
    if roots != torsion["roots"] or bound != torsion["dimension_upper_bound"]:
        raise ValueError("Incorrect torsion witness")
    basis = {}
    for row in certificate["independent_rows"]:
        p, alpha = row["prime"], row["root"]
        if not is_prime(p) or p == 2 or disc % p == 0:
            raise ValueError("Invalid character prime")
        a, b, c = cubic
        if not 0 <= alpha < p or (alpha**3+a*alpha**2+b*alpha+c) % p:
            raise ValueError("Not a cubic root modulo the given prime")
        mask = character_row(points, cubic, p, alpha)
        if bitstring(mask, n) != row["bits"]:
            raise ValueError("Incorrect character row")
        if not add_row(basis, mask):
            raise ValueError("Certificate contains a dependent row")
    claims = {
        "points_checked": n,
        "matrix_rank": len(basis),
        "torsion_2_dimension_upper_bound": bound,
        "rank_lower_bound": max(0, len(basis)-bound),
        "all_points_independent_modulo_torsion": len(basis) == n and bound == 0,
        "conditional_assumptions": [],
    }
    for key, value in claims.items():
        if certificate.get(key) != value:
            raise ValueError(f"Incorrect certificate claim: {key}")
    return claims


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="JSON with ainvs and points")
    parser.add_argument("--output", type=Path, help="Write generated certificate JSON")
    parser.add_argument("--verify", type=Path, help="Replay an existing certificate")
    parser.add_argument("--max-prime", type=int, default=10000)
    args = parser.parse_args()
    data = json.loads(args.input.read_text())
    started = time.perf_counter()
    if args.verify:
        cert = json.loads(args.verify.read_text())
    else:
        cert = build_certificate(data, max_prime=args.max_prime)
    result = verify_certificate(data, cert)
    if args.output:
        args.output.write_text(json.dumps(cert, indent=2) + "\n")
    result["elapsed_seconds"] = round(time.perf_counter()-started, 6)
    result["largest_prime_examined"] = cert["largest_prime_examined"]
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
