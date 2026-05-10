"""
Eksperyment 3 — część szybka: wszystko poza Setup_CC.

Setup_CC mierzymy osobno (skrypt exp3b_setup_cc.py), żeby nie blokować
całego eksperymentu długą generacją safe-prime.
"""
import os, sys, time, json, statistics
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from hed_voting import paillier as P, elgamal as E, rsa_sig as R


def time_fn(fn, repeat: int) -> tuple[float, float]:
    times = []
    for _ in range(repeat):
        t0 = time.perf_counter()
        fn()
        times.append(time.perf_counter() - t0)
    return statistics.median(times), statistics.stdev(times) if len(times) > 1 else 0.0


def main():
    out_dir = os.path.join(os.path.dirname(__file__), "results")
    os.makedirs(out_dir, exist_ok=True)

    print("=== Eksperyment 3a: czas poszczególnych faz protokołu (bez Setup_CC) ===\n")

    print("Generacja kluczy testowych...")
    php, phs    = P.setup_he(1024)
    egp, egs    = E.setup_cc(2 * 1024 + 16)
    upub, upriv = R.user_keygen(1024)

    test_m = 7
    c1 = P.enc_he(php, test_m)
    sig = R.rsa_sig(upriv, int(c1))
    C   = E.enc_cc(egp, int(c1))

    print("Pomiary (mediana czasu pojedynczej operacji):\n")
    print(f"   {'Algorytm':<14} {'Mediana':>14} {'Std':>14} {'#':>5}")
    print("   " + "-"*52)

    fns = [
        ("Setup_HE",    lambda: P.setup_he(1024),               20),
        ("User_KeyGen", lambda: R.user_keygen(1024),            10),
        ("Enc_HE",      lambda: P.enc_he(php, test_m),         300),
        ("Rsa_Sig",     lambda: R.rsa_sig(upriv, int(c1)),     300),
        ("Enc_CC",      lambda: E.enc_cc(egp, int(c1)),         80),
        ("Dec_CC",      lambda: E.dec_cc(egs, C),               80),
        ("Rsa_Ver",     lambda: R.rsa_ver(upub, int(c1), sig), 300),
        ("Count_CC",    lambda: P.add_he(php, c1, c1),       2000),
        ("Dec_HE",      lambda: P.dec_he(phs, c1),             300),
    ]

    rows = []
    for name, fn, rep in fns:
        med, std = time_fn(fn, rep)
        rows.append({"name": name, "median_ms": med*1000, "std_ms": std*1000, "repeat": rep})
        print(f"   {name:<14} {med*1000:>12.4f} ms {std*1000:>12.4f} ms {rep:>5}")

    out = {"paillier_bits": 1024, "elgamal_bits": 2050, "rows": rows}
    with open(os.path.join(out_dir, "exp3a_phase_timing.json"), "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nZapisano: {os.path.join(out_dir, 'exp3a_phase_timing.json')}")


if __name__ == "__main__":
    main()
