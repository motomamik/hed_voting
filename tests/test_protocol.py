"""
Test poprawności pełnego protokołu HED-Voting:
  - 5 kandydatów,
  - 50 wyborców,
  - oddają głosy zgodnie z zadanym rozkładem,
  - sprawdzamy, czy końcowa lista głosów się zgadza.
"""
import sys, os, random
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from hed_voting import protocol as H
from hed_voting import rsa_sig    as R


def main():
    random.seed(0)
    N_CAND = 5
    N_VOTERS = 50

    # Symulujemy oczekiwany rozkład głosów
    expected = [0] * N_CAND
    choices  = []
    for _ in range(N_VOTERS):
        c = random.randrange(N_CAND)
        choices.append(c)
        expected[c] += 1
    print(f"Oczekiwane wyniki: {expected}")

    print("Inicjalizacja systemu (klucz Paillier 1024 b, ElGamal automatycznie 2*1024+16 b)...")
    sys_ = H.setup_system(n_candidates=N_CAND, n_voters_max=N_VOTERS,
                          paillier_bits=1024)

    print("Generacja kluczy wyborców...")
    voters = [R.user_keygen(1024) for _ in range(N_VOTERS)]

    print("Faza głosowania + zliczania w CC...")
    total = None
    for i, choice in enumerate(choices):
        upub, upriv = voters[i]
        C, s  = H.voter_cast(sys_, upriv, choice)
        total = H.cc_process(sys_, upub, C, s, total)

    print("Ogłoszenie wyników (Dec_HE)...")
    result = H.reveal(sys_, total)
    print(f"Otrzymane wyniki:  {result}")

    assert result == expected, "BŁĄD: wyniki niezgodne z oczekiwanymi!"
    print("OK — protokół działa poprawnie.")

if __name__ == "__main__":
    main()
