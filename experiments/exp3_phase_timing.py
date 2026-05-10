"""
Eksperyment 3 — Czas wykonania poszczególnych faz protokołu HED-Voting.

Mierzymy czas każdego z 10 algorytmów składowych protokołu
(zgodnie z Tabelą 5 z artykułu Yuan et al. 2023):

    Setup_HE      — generacja pary kluczy Paillier
    Setup_CC      — generacja pary kluczy ElGamal
    User_KeyGen   — generacja pary kluczy RSA wyborcy
    Enc_HE        — szyfrowanie homomorficzne głosu
    Rsa_Sig       — podpis RSA wyborcy
    Enc_CC        — szyfrowanie ElGamal
    Dec_CC        — deszyfracja ElGamal w CC
    Rsa_Ver       — weryfikacja podpisu w CC
    Count_CC      — homomorficzne sumowanie w CC
    Dec_HE        — deszyfracja Paillier (ogłoszenie wyniku)

Każdy algorytm uruchamiamy wielokrotnie i raportujemy medianę.
"""
import os, sys, time, json, statistics
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from hed_voting import paillier as P, elgamal as E, rsa_sig as R, protocol as H


REPEAT_FAST = 200    # dla algorytmów < 1 ms
REPEAT_SLOW = 5      # dla algorytmów ~10ms-1s  (Setup_HE, User_KeyGen, Setup_CC)


def time_fn(fn, repeat: int) -> tuple[float, float]:
    """Zwraca (mediana, odch. std) czasu wykonania funkcji w sekundach."""
    times = []
    for _ in range(repeat):
        t0 = time.perf_counter()
        fn()
        times.append(time.perf_counter() - t0)
    return statistics.median(times), statistics.stdev(times) if len(times) > 1 else 0.0


def main():
    out_dir = os.path.join(os.path.dirname(__file__), "results")
    os.makedirs(out_dir, exist_ok=True)

    print("=== Eksperyment 3: czas poszczególnych faz protokołu (1024 b) ===\n")

    # 1. Setup_HE
    print("Generacja kluczy testowych...")
    php, phs = P.setup_he(1024)
    egp, egs = E.setup_cc(2 * 1024 + 16)
    upub, upriv = R.user_keygen(1024)

    # Test plaintext
    test_m = 7
    c1 = P.enc_he(php, test_m)
    sig = R.rsa_sig(upriv, int(c1))
    C   = E.enc_cc(egp, int(c1))

    print("Pomiary (mediana czasu pojedynczej operacji):\n")
    print(f"   {'Algorytm':<14} {'Mediana':>12} {'Std':>10}")
    print("   " + "-"*40)

    # Dobieramy liczbę powtórzeń indywidualnie do każdego algorytmu
    fns = [
        ("Setup_HE",    lambda: P.setup_he(1024),       10),   # ~50-100 ms
        ("Setup_CC",    lambda: E.setup_cc(2*1024+16),  5),    # ~1-3 s (po optymalizacji sieve)
        ("User_KeyGen", lambda: R.user_keygen(1024),    10),   # ~50 ms
        ("Enc_HE",      lambda: P.enc_he(php, test_m),  REPEAT_FAST),
        ("Rsa_Sig",     lambda: R.rsa_sig(upriv, int(c1)), REPEAT_FAST),
        ("Enc_CC",      lambda: E.enc_cc(egp, int(c1)), 50),   # 2050 b -> wolniejsze potęgowanie
        ("Dec_CC",      lambda: E.dec_cc(egs, C),       50),
        ("Rsa_Ver",     lambda: R.rsa_ver(upub, int(c1), sig), REPEAT_FAST),
        ("Count_CC",    lambda: P.add_he(php, c1, c1),  REPEAT_FAST),
        ("Dec_HE",      lambda: P.dec_he(phs, c1),      REPEAT_FAST),
    ]

    rows = []
    for name, fn, rep in fns:
        med, std = time_fn(fn, rep)
        rows.append({"name": name, "median_ms": med*1000, "std_ms": std*1000})
        print(f"   {name:<14} {med*1000:>10.3f} ms {std*1000:>8.3f} ms")

    out = {"paillier_bits": 1024, "rows": rows}
    with open(os.path.join(out_dir, "exp3_phase_timing.json"), "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nZapisano: {os.path.join(out_dir, 'exp3_phase_timing.json')}")

if __name__ == "__main__":
    main()
