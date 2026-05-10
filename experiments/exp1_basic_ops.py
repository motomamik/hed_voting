"""
Eksperyment 1 — Reprodukcja Tabeli 4 i 5 z artykułu Yuan et al. (2023).

Zgodnie z sekcją "HED-voting scheme computational complexity":
    - Paillier 1024 bitów (n = pq, p i q po 512 bitach), więc n^2 = 2048 bitów.
    - Operacje wykonywane są w gmpy2 oraz phe (Pallier key generation).
    - Bazową jednostką jest Mul1 — modulo-multiplication w Z*_{n^2}.

Mierzone operacje:
    Mul1 — modulo-multiplication na Z*_{n_1^2}
    Mul2 — modulo-multiplication na Z*_{n_1}
    Pow  — modulo-power on Z*_{n_1^2}
    Inv  — modulo-inverse on Z*_{n_1^2}
    Key  — generacja klucza Paillier (phe.generate_paillier_keypair)
"""
import os, sys, time, json, statistics
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import gmpy2
from gmpy2 import mpz, powmod, invert
import phe

from hed_voting.paillier import _random_prime


# ---------------------------------------------------------------------------
# Parametry zgodne z artykułem
# ---------------------------------------------------------------------------
PAILLIER_BITS = 1024          # n_1 ma 1024 bitów (p, q po 512)
REPEAT        = 200_000       # Liczba powtórzeń dla operacji szybkich
REPEAT_POW    = 2_000         # Liczba powtórzeń dla operacji wolnych
REPEAT_KEY    = 30            # Liczba powtórzeń dla generacji klucza


def gen_modulus():
    p = _random_prime(PAILLIER_BITS // 2)
    q = _random_prime(PAILLIER_BITS // 2)
    n = p * q
    return p, q, n


def measure_mul_n2(n: mpz, repeat: int) -> float:
    """Średni czas pojedynczego mnożenia modulo na Z*_{n^2}."""
    n2 = n * n
    a  = mpz(0xC0FFEE) % n2
    b  = mpz(0xBADBEEF) % n2
    # rozgrzewka
    for _ in range(1000): _ = (a * b) % n2
    t0 = time.perf_counter()
    for _ in range(repeat):
        _ = (a * b) % n2
    return (time.perf_counter() - t0) / repeat


def measure_mul_n(n: mpz, repeat: int) -> float:
    """Średni czas pojedynczego mnożenia modulo na Z*_n."""
    a = mpz(0xC0FFEE) % n
    b = mpz(0xBADBEEF) % n
    for _ in range(1000): _ = (a * b) % n
    t0 = time.perf_counter()
    for _ in range(repeat):
        _ = (a * b) % n
    return (time.perf_counter() - t0) / repeat


def measure_pow_n2(n: mpz, repeat: int) -> float:
    """Średni czas pojedynczego potęgowania modulo na Z*_{n^2}."""
    n2 = n * n
    base = mpz(2) ** 5 + 7
    exp  = (mpz(0xDEADBEEFCAFE) << 1000) | 1   # eksponent ~1024 b
    for _ in range(20): _ = powmod(base, exp, n2)
    t0 = time.perf_counter()
    for _ in range(repeat):
        _ = powmod(base, exp, n2)
    return (time.perf_counter() - t0) / repeat


def measure_inv_n2(n: mpz, repeat: int) -> float:
    """Średni czas pojedynczego odwracania modulo na Z*_{n^2}."""
    n2 = n * n
    a  = mpz(2) ** 5 + 7
    for _ in range(50): _ = invert(a, n2)
    t0 = time.perf_counter()
    for _ in range(repeat):
        _ = invert(a, n2)
    return (time.perf_counter() - t0) / repeat


def measure_paillier_keygen(repeat: int) -> float:
    """Średni czas generacji klucza Paillier (zgodnie z artykułem — biblioteka phe)."""
    # rozgrzewka
    _ = phe.generate_paillier_keypair(n_length=PAILLIER_BITS)
    t0 = time.perf_counter()
    for _ in range(repeat):
        _ = phe.generate_paillier_keypair(n_length=PAILLIER_BITS)
    return (time.perf_counter() - t0) / repeat


def main():
    print("=== Eksperyment 1: pomiar względnego czasu operacji bazowych ===")
    print(f"Klucz Paillier: {PAILLIER_BITS} bitów  (n^2 ≈ {2*PAILLIER_BITS} bitów)")

    p, q, n = gen_modulus()
    print(f"Wygenerowano moduł n o długości {n.bit_length()} bitów.")

    print("\nPomiary (czas pojedynczej operacji w mikrosekundach):")
    t_mul2 = measure_mul_n2(n, REPEAT)
    print(f"  Mul1  (mod n^2)  =  {t_mul2*1e6:.4f} µs   ({REPEAT} powtórzeń)")

    t_mul1 = measure_mul_n(n, REPEAT)
    print(f"  Mul2  (mod n)    =  {t_mul1*1e6:.4f} µs   ({REPEAT} powtórzeń)")

    t_pow  = measure_pow_n2(n, REPEAT_POW)
    print(f"  Pow   (mod n^2)  =  {t_pow*1e6:.4f} µs   ({REPEAT_POW} powtórzeń)")

    t_inv  = measure_inv_n2(n, REPEAT_POW * 5)
    print(f"  Inv   (mod n^2)  =  {t_inv*1e6:.4f} µs   ({REPEAT_POW*5} powtórzeń)")

    t_key  = measure_paillier_keygen(REPEAT_KEY)
    print(f"  Key   (Paillier) =  {t_key*1e3:.4f} ms   ({REPEAT_KEY} powtórzeń)")

    # Względny czas (Mul1 = 1)
    rel = {
        "Mul1": 1.0,
        "Mul2": t_mul1 / t_mul2,
        "Pow":  t_pow  / t_mul2,
        "Inv":  t_inv  / t_mul2,
        "Key":  t_key  / t_mul2,
    }
    print("\nWzględny czas (Mul1 = 1):")
    for k, v in rel.items():
        print(f"  {k:5s} = {v:.4f}")

    # Złożoność algorytmów (Tabela 5 artykułu) — przeliczona na nasz zegar
    complexity = {
        "Setup_HE":     {"formula": "Key",                    "value": rel["Key"]},
        "Setup_CC":     {"formula": "Mul2",                   "value": rel["Mul2"]},
        "User_KeyGen":  {"formula": "Inv",                    "value": rel["Inv"]},
        "Enc_HE":       {"formula": "Pow + Mul1",             "value": rel["Pow"] + rel["Mul1"]},
        "Rsa_Sig":      {"formula": "Pow",                    "value": rel["Pow"]},
        "Enc_CC":       {"formula": "2*Pow",                  "value": 2*rel["Pow"]},
        "Dec_CC":       {"formula": "3*Mul2 + 3*Pow + Inv",   "value": 3*rel["Mul2"] + 3*rel["Pow"] + rel["Inv"]},
        "Rsa_Ver":      {"formula": "Pow",                    "value": rel["Pow"]},
        "Count_CC":     {"formula": "2*Pow + 2*Mul1",         "value": 2*rel["Pow"] + 2*rel["Mul1"]},
        "Dec_HE":       {"formula": "Pow + Mul2",             "value": rel["Pow"] + rel["Mul2"]},
    }
    print("\nZłożoność algorytmów (jednostka = czas Mul1):")
    for k, v in complexity.items():
        print(f"  {k:14s} {v['formula']:24s} = {v['value']:.2f}")

    # Zapis do pliku JSON i CSV
    out = {
        "paillier_bits": PAILLIER_BITS,
        "abs_times_us": {
            "Mul1": t_mul2*1e6, "Mul2": t_mul1*1e6,
            "Pow":  t_pow*1e6,  "Inv":  t_inv*1e6,
            "Key":  t_key*1e6,
        },
        "relative":  rel,
        "complexity": complexity,
    }
    out_dir = os.path.join(os.path.dirname(__file__), "results")
    os.makedirs(out_dir, exist_ok=True)
    with open(os.path.join(out_dir, "exp1_basic_ops.json"), "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nZapisano: {os.path.join(out_dir, 'exp1_basic_ops.json')}")


if __name__ == "__main__":
    main()
