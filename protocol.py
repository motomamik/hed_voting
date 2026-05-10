"""
HED-Voting — pełny protokół zgodny z artykułem
Yuan, K., Sang, P., Zhang, S., Chen, X., Yang, W., Jia, C. (2023).
"An electronic voting scheme based on homomorphic encryption and decentralization",
PeerJ Computer Science 9:e1649. DOI: 10.7717/peerj-cs.1649.

Kodowanie głosów dla Nc kandydatów (sekcja "Performance Analysis"):
    Każdy głos jest wektorem 0/1 o długości Nc, ale zapisywany jako jedna liczba
    w bazie B = (Nv + 1), gdzie Nv to maks. liczba wyborców. Wówczas
    homomorficzne sumowanie ciphertekstów Paillier sumuje wektory komponentowo,
    a po deszyfracji można odczytać liczbę głosów oddaną na każdego kandydata
    przez kolejne dzielenie modulo. Jest to standardowy trick przy
    Paillier-multi-candidate i odpowiada opisowi z artykułu (rozdział
    "Performance analysis of the voter client": "each voter only needs to cast
    one vote for his or her target candidate ... and does not need to perform
    any operation on the other candidates").
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import List, Tuple
from gmpy2 import mpz

from hed_voting import paillier as P
from hed_voting import elgamal  as E
from hed_voting import rsa_sig  as R


# ---------------------------------------------------------------------------
# Kodowanie głosu  (multi-candidate Paillier)
# ---------------------------------------------------------------------------
def encode_vote(choice: int, n_candidates: int, base: int) -> int:
    """
    choice ∈ {0, ..., n_candidates-1}: indeks wybranego kandydata.
    Zwraca: base^choice  (jedynka na pozycji 'choice' w zapisie pozycyjnym).
    """
    if not 0 <= choice < n_candidates:
        raise ValueError("Nieprawidłowy indeks kandydata.")
    return base ** choice


def decode_tally(total: int, n_candidates: int, base: int) -> List[int]:
    """Z odszyfrowanej sumy wyciąga liczbę głosów na każdego kandydata."""
    counts = []
    for _ in range(n_candidates):
        counts.append(total % base)
        total //= base
    return counts


# ---------------------------------------------------------------------------
# Stany / role
# ---------------------------------------------------------------------------
@dataclass
class HEDSystem:
    """
    Pełna instancja systemu HED-Voting.
        - PHE_pub / PHE_priv  : klucze Paillier (KC -> kandydat A)
        - CC_pub  / CC_priv   : klucze ElGamal (KC -> kandydat B = CC)
        - n_candidates, base  : parametry kodowania głosu
    """
    PHE_pub:  P.PaillierPublicKey
    PHE_priv: P.PaillierPrivateKey
    CC_pub:   E.ElGamalPublicKey
    CC_priv:  E.ElGamalPrivateKey
    n_candidates: int
    base: int


def setup_system(n_candidates: int, n_voters_max: int,
                 paillier_bits: int = 1024,
                 elgamal_bits: int  | None = None) -> HEDSystem:
    """
    Inicjalizacja zgodna z 1. fazą artykułu ("Initialization").
    Klucze są generowane przez Key Center, a klucze prywatne PHE_priv
    i CC_priv są przekazywane dwóm konkurencyjnym kandydatom.

    UWAGA techniczna:
        Ciphertext Paillier należy do Z*_{n^2}, więc ma rozmiar do 2*paillier_bits.
        Aby Enc_CC nie redukował go modulo p (co psuje weryfikację podpisu RSA),
        moduł ElGamal musi być co najmniej tak duży jak n^2. Dlatego domyślnie
        ustawiamy elgamal_bits = 2*paillier_bits + 16.
    """
    if elgamal_bits is None:
        elgamal_bits = 2 * paillier_bits + 16
    if elgamal_bits < 2 * paillier_bits + 2:
        raise ValueError(
            f"elgamal_bits ({elgamal_bits}) musi być >= 2*paillier_bits+2 "
            f"({2*paillier_bits+2}) aby ciphertext Paillier mieścił się w grupie ElGamal."
        )
    php, phs = P.setup_he(paillier_bits)
    egp, egs = E.setup_cc(elgamal_bits)
    base = n_voters_max + 1
    return HEDSystem(
        PHE_pub=php, PHE_priv=phs,
        CC_pub=egp,  CC_priv=egs,
        n_candidates=n_candidates,
        base=base,
    )


# ---------------------------------------------------------------------------
# Strona klienta  (wyborca)
# ---------------------------------------------------------------------------
def voter_cast(sys: HEDSystem,
               voter_priv: R.RSAPrivateKey,
               choice: int) -> Tuple[Tuple[mpz, mpz], mpz]:
    """
    Operacje wyborcy zgodne z sekcją "Performance analysis of the voter client":
        1) Enc_HE   : c'_i = Paillier(m)            -- 1 modulo-power
        2) Rsa_Sig  : s    = c'_i^d mod n2          -- 1 modulo-power
        3) Enc_CC   : C = (g^y, c'_i_mapped * h^y)  -- 2 modulo-powers
       razem 4 mnożenia decentralizowane (Coste) na wyborcę i kandydata.

    Zwraca: (C, s)  -- ciphertext zewnętrzny ElGamal i podpis RSA.
    """
    # 1. Pierwsza warstwa: szyfrowanie homomorficzne Paillier
    m   = encode_vote(choice, sys.n_candidates, sys.base)
    c1  = P.enc_he(sys.PHE_pub, m)

    # 2. Druga warstwa: podpis RSA wyborcy nad c1
    s   = R.rsa_sig(voter_priv, int(c1))

    # 3. Trzecia warstwa: ElGamal kluczem CC, na c1
    #    (zgodnie z opisem ct_i = (c'_i, s); tutaj szyfrujemy c'_i,
    #     a podpis transmitowany jest osobno — co odpowiada artykułowi:
    #     na rys. 2 "Layer-3 encryption" obejmuje (c'_i, s),
    #     ale dla wektora wymiarów liczb i potrzeb deszyfracji wystarczy
    #     przekazać (Enc_CC(c'_i), s) — ekwiwalentnie na potrzeby protokołu.)
    C   = E.enc_cc(sys.CC_pub, int(c1))
    return C, s


# ---------------------------------------------------------------------------
# Strona Counting Center
# ---------------------------------------------------------------------------
def cc_process(sys: HEDSystem,
               voter_pub: R.RSAPublicKey,
               C: Tuple[mpz, mpz],
               s: mpz,
               running_total: mpz | None) -> mpz:
    """
    Operacje CC zgodne z sekcją "Performance analysis of the CC":
        1) Dec_CC   : odzyskuje c'_i z C       -- 1 modulo-power
        2) Rsa_Ver  : sprawdza podpis          -- 1 modulo-power
        3) Count_CC : homomorficzne sumowanie  -- mnożenie wielkich liczb (zaniedbywalne)
       razem 2 mnożenia decentralizowane (Coste) na głos.

    Zwraca: zaktualizowane running_total (Paillier-ciphertext).
    """
    c1_back = int(E.dec_cc(sys.CC_priv, C))
    if not R.rsa_ver(voter_pub, c1_back, s):
        raise ValueError("Podpis RSA wyborcy niepoprawny — głos odrzucony.")
    if running_total is None:
        return mpz(c1_back)
    return P.add_he(sys.PHE_pub, running_total, mpz(c1_back))


# ---------------------------------------------------------------------------
# Faza ogłoszenia wyniku
# ---------------------------------------------------------------------------
def reveal(sys: HEDSystem, total_ct: mpz) -> List[int]:
    """
    Dec_HE: deszyfracja sumy ciphertekstów Paillier kluczem PHE_priv,
    a następnie zdekodowanie liczby głosów per kandydat.
    """
    plain = P.dec_he(sys.PHE_priv, total_ct)
    return decode_tally(plain, sys.n_candidates, sys.base)
