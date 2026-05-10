"""
Paillier Homomorphic Encryption (PHE)
Implementacja zgodna z Definition 3 oraz sekcją "Setup_HE", "Enc_HE", "Count_CC", "Dec_HE"
z artykułu Yuan et al. (2023): "An electronic voting scheme based on homomorphic
encryption and decentralization", PeerJ Computer Science, 9:e1649.

Wszystkie operacje arytmetyczne wykonywane są w gmpy2 (zgodnie z artykułem,
sekcja "Performance Analysis": "Testing with the gmpy2 python module").
"""
from __future__ import annotations
from dataclasses import dataclass
import secrets
import gmpy2
from gmpy2 import mpz, powmod, gcd, invert, lcm


# ---------------------------------------------------------------------------
# Pomocnicze funkcje teorii liczb
# ---------------------------------------------------------------------------
def _random_prime(bits: int) -> mpz:
    """Generuje losową liczbę pierwszą o zadanej liczbie bitów."""
    while True:
        n = mpz(secrets.randbits(bits) | (1 << (bits - 1)) | 1)
        if gmpy2.is_prime(n, 25):
            return n


def _L(x: mpz, n: mpz) -> mpz:
    """Funkcja L(x) = (x - 1) / n  (zgodnie z Definition 3)."""
    return (x - 1) // n


# ---------------------------------------------------------------------------
# Klucze i schemat Paillier
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class PaillierPublicKey:
    n: mpz   # n = p*q
    g: mpz   # generator z Z*_{n^2}

    @property
    def n2(self) -> mpz:
        return self.n * self.n


@dataclass(frozen=True)
class PaillierPrivateKey:
    lam: mpz   # lambda = lcm(p-1, q-1)  -- "k" w artykule
    mu:  mpz   # mu = [L(g^lambda mod n^2)]^{-1} mod n  -- "l" w artykule
    pub: PaillierPublicKey


def setup_he(key_bits: int = 1024) -> tuple[PaillierPublicKey, PaillierPrivateKey]:
    """
    Setup_HE  (Yuan et al., 2023, sekcja "HED-voting solution").

    Generuje parę kluczy Paillier:
        HE_pub  = (n, g)        — klucz publiczny szyfrowania homomorficznego
        HE_priv = (lambda, mu)  — klucz prywatny

    Zgodnie z artykułem (sekcja "HED-Voting Scheme Computational Complexity")
    używamy n o długości 1024 bitów (p, q po 512 bitów).
    """
    half = key_bits // 2
    while True:
        p = _random_prime(half)
        q = _random_prime(half)
        if p != q:
            n = p * q
            # gcd(pq, (p-1)(q-1)) == 1  (warunek poprawności Paillier)
            if gcd(n, (p - 1) * (q - 1)) == 1:
                break

    lam = lcm(p - 1, q - 1)
    n2  = n * n
    g   = n + 1                                  # standardowy, bezpieczny wybór: g = n+1
    mu  = invert(_L(powmod(g, lam, n2), n), n)

    pub  = PaillierPublicKey(n=n, g=g)
    priv = PaillierPrivateKey(lam=lam, mu=mu, pub=pub)
    return pub, priv


def enc_he(pub: PaillierPublicKey, m: int) -> mpz:
    """
    Enc_HE  (Yuan et al., 2023): c'_i = g^m * r^n  mod n^2,
    gdzie r jest losowym elementem Z*_n (gcd(r, n) = 1).
    """
    n, g, n2 = pub.n, pub.g, pub.n2
    while True:
        r = mpz(secrets.randbelow(int(n - 1))) + 1
        if gcd(r, n) == 1:
            break
    m = mpz(m) % n
    return (powmod(g, m, n2) * powmod(r, n, n2)) % n2


def add_he(pub: PaillierPublicKey, c1: mpz, c2: mpz) -> mpz:
    """
    Count_CC  (Yuan et al., 2023): operacja homomorficznej addytywności.
    E(m1) * E(m2) mod n^2  =  E(m1 + m2)
    """
    return (mpz(c1) * mpz(c2)) % pub.n2


def dec_he(priv: PaillierPrivateKey, c: mpz) -> int:
    """
    Dec_HE: m = L(c^lambda mod n^2) * mu  mod n.
    """
    n, n2 = priv.pub.n, priv.pub.n2
    u = powmod(mpz(c), priv.lam, n2)
    return int((_L(u, n) * priv.mu) % n)
