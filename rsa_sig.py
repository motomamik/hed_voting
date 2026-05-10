"""
Podpis RSA (warstwa druga "szyfrowania" w schemacie HED-Voting)

Zgodnie z artykułem Yuan et al. (2023), sekcja "User_KeyGen" oraz "Rsa_Sig":
    upub  = (n2, e)
    upriv = (p2, q2, d)
    s = sig_upriv(c'_i) = (c'_i)^d  mod n2

Weryfikacja (Rsa_Ver):
    c''_i = s^e mod n2
    podpis poprawny  <=>  c''_i ≡ c'_i (mod n2)

Artykuł stosuje "tekstowy" podpis bez funkcji skrótu — realizujemy zgodnie
z opisem (zachowując tym samym tę samą złożoność: 1 modulo-power przy
podpisie, 1 modulo-power przy weryfikacji, jak w Tabeli 4-5 artykułu).
"""
from __future__ import annotations
from dataclasses import dataclass
import secrets
import gmpy2
from gmpy2 import mpz, powmod, gcd, invert


def _random_prime(bits: int) -> mpz:
    while True:
        n = mpz(secrets.randbits(bits) | (1 << (bits - 1)) | 1)
        if gmpy2.is_prime(n, 25):
            return n


@dataclass(frozen=True)
class RSAPublicKey:
    n: mpz
    e: mpz


@dataclass(frozen=True)
class RSAPrivateKey:
    d: mpz
    pub: RSAPublicKey


def user_keygen(key_bits: int = 1024) -> tuple[RSAPublicKey, RSAPrivateKey]:
    """User_KeyGen: generuje (upub, upriv) — klasyczny RSA."""
    half = key_bits // 2
    p = _random_prime(half)
    q = _random_prime(half)
    while p == q:
        q = _random_prime(half)
    n = p * q
    phi = (p - 1) * (q - 1)
    e = mpz(65537)
    if gcd(e, phi) != 1:
        e = mpz(3)
        while gcd(e, phi) != 1:
            e += 2
    d = invert(e, phi)
    pub  = RSAPublicKey(n=n, e=e)
    priv = RSAPrivateKey(d=d, pub=pub)
    return pub, priv


def rsa_sig(priv: RSAPrivateKey, m: int) -> mpz:
    """s = m^d mod n"""
    return powmod(mpz(m) % priv.pub.n, priv.d, priv.pub.n)


def rsa_ver(pub: RSAPublicKey, m: int, s: mpz) -> bool:
    """Sprawdza, czy s^e ≡ m (mod n)."""
    return powmod(mpz(s), pub.e, pub.n) == (mpz(m) % pub.n)
