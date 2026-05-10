"""
ElGamal Public-Key Encryption (EPKE)
Implementacja zgodna z Definition 4 oraz sekcjami "Setup_CC", "Enc_CC", "Dec_CC"
z artykułu Yuan et al. (2023).

W artykule:
    - "q1" jest dużą liczbą pierwszą definiującą grupę cykliczną G,
    - g1 jest jej generatorem (oryginalnym pierwiastkiem),
    - x ∈ (1, q1-1) jest kluczem prywatnym,
    - h = g1^x mod q1.

Tutaj realizujemy EPKE w grupie multiplikatywnej Z*_p, gdzie p jest dużą
liczbą pierwszą bezpieczną (p = 2q+1) — najczęściej spotykana realizacja
i identyczna z opisem z artykułu.
"""
from __future__ import annotations
from dataclasses import dataclass
import secrets
import gmpy2
from gmpy2 import mpz, powmod, invert


# Małe liczby pierwsze do sieve trial-division (przyspiesza generację 50-100x)
_SMALL_PRIMES = [
    3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67,
    71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139,
    149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223,
    227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293,
    307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383,
    389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463,
    467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569,
    571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647,
    653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743,
    751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839,
    853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941,
    947, 953, 967, 971, 977, 983, 991, 997
]


# ---------------------------------------------------------------------------
# Bezpieczne liczby pierwsze: p = 2q+1, q pierwsza  (z trial division sieve)
# ---------------------------------------------------------------------------
def _safe_prime(bits: int) -> tuple[mpz, mpz]:
    """
    Generuje safe-prime p = 2q+1, q pierwsza, o p długości `bits`.
    Optymalizacja: trial division przez 168 najmniejszych liczb pierwszych
    eliminuje ~99% kandydatów przed kosztownym testem Millera-Rabina.
    """
    # q ma bits-1 bitów, więc p = 2q+1 ma bits bitów
    qbits = bits - 1
    while True:
        q = mpz(secrets.randbits(qbits) | (1 << (qbits - 1)) | 1)
        # q == 5 mod 6  (warunek konieczny by zarówno q jak i 2q+1 były nieparzyste i niepodzielne przez 3)
        if q % 6 != 5:
            q += (5 - q % 6) % 6
            if q.bit_length() != qbits:
                continue
        # trial division
        if any(q % p == 0 for p in _SMALL_PRIMES):
            continue
        p = 2 * q + 1
        if any(p % pp == 0 for pp in _SMALL_PRIMES):
            continue
        # Test Millera-Rabina
        if not gmpy2.is_prime(q, 5):
            continue
        if not gmpy2.is_prime(p, 25):
            continue
        if not gmpy2.is_prime(q, 25):
            continue
        return p, q


def _generator(p: mpz, q: mpz) -> mpz:
    """Element rzędu 2q (czyli generator Z*_p) w grupie z bezpieczną liczbą pierwszą."""
    while True:
        h = mpz(secrets.randbelow(int(p - 3))) + 2
        if powmod(h, 2, p) != 1 and powmod(h, q, p) != 1:
            return h


# ---------------------------------------------------------------------------
# Klucze i schemat ElGamal
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class ElGamalPublicKey:
    p: mpz
    g: mpz
    h: mpz   # h = g^x mod p


@dataclass(frozen=True)
class ElGamalPrivateKey:
    x: mpz
    pub: ElGamalPublicKey


def setup_cc(key_bits: int = 1024) -> tuple[ElGamalPublicKey, ElGamalPrivateKey]:
    """
    Setup_CC: generacja pary kluczy Counting Center (CCpub, CCpriv).
        CCpub  = (p, g, h)
        CCpriv = x
    """
    p, q = _safe_prime(key_bits)
    g = _generator(p, q)
    x = mpz(secrets.randbelow(int(p - 3))) + 2
    h = powmod(g, x, p)
    pub  = ElGamalPublicKey(p=p, g=g, h=h)
    priv = ElGamalPrivateKey(x=x, pub=pub)
    return pub, priv


def enc_cc(pub: ElGamalPublicKey, m: int) -> tuple[mpz, mpz]:
    """
    Enc_CC: c = (a, b) = (g^y mod p, m * h^y mod p),  y losowe.

    UWAGA: w artykule (sekcja "Enc_CC") wzór jest podany w tej postaci,
    przy czym m jest mapowane na element grupy. Tutaj m jest już liczbą
    całkowitą < p (wartość ciphertextu Paillier rzutowana modulo p).
    """
    p, g, h = pub.p, pub.g, pub.h
    y = mpz(secrets.randbelow(int(p - 3))) + 2
    a = powmod(g, y, p)
    s = powmod(h, y, p)
    b = (mpz(m) * s) % p
    return a, b


def dec_cc(priv: ElGamalPrivateKey, c: tuple[mpz, mpz]) -> mpz:
    """
    Dec_CC: m = b * (a^x)^{-1} mod p
    """
    a, b = c
    p = priv.pub.p
    s = powmod(a, priv.x, p)
    return (mpz(b) * invert(s, p)) % p
