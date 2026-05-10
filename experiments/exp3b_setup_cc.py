"""
Eksperyment 3b — pomiar Setup_CC (oddzielnie, bo jest najwolniejszy).

Generacja safe-prime 2050-bitowego z trial-division sieve to ~1-3 s.
"""
import os, sys, time, json, statistics
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from hed_voting import elgamal as E

def main():
    out_dir = os.path.join(os.path.dirname(__file__), "results")
    os.makedirs(out_dir, exist_ok=True)

    print("=== Eksperyment 3b: Setup_CC (safe-prime 2050 b z sieve trial-division) ===\n")

    times = []
    REPEAT = 8
    for i in range(REPEAT):
        t0 = time.perf_counter()
        E.setup_cc(2 * 1024 + 16)
        t = time.perf_counter() - t0
        times.append(t)
        print(f"   próba {i+1}/{REPEAT}: {t*1000:8.1f} ms")

    med = statistics.median(times)
    std = statistics.stdev(times) if len(times) > 1 else 0.0
    print(f"\n   Mediana: {med*1000:.2f} ms   ({statistics.mean(times)*1000:.2f} ± {std*1000:.2f} ms)")

    # Aktualizujemy JSON z exp3a
    fname = os.path.join(out_dir, "exp3a_phase_timing.json")
    if os.path.exists(fname):
        with open(fname) as f:
            data = json.load(f)
        data["rows"].insert(1, {
            "name": "Setup_CC", "median_ms": med*1000,
            "std_ms": std*1000, "repeat": REPEAT,
        })
        out_combined = os.path.join(out_dir, "exp3_phase_timing.json")
        with open(out_combined, "w") as f:
            json.dump(data, f, indent=2)
        print(f"\nZapisano scalony plik: {out_combined}")


if __name__ == "__main__":
    main()
