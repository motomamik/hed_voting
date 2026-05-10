"""
Eksperyment 2 — Reprodukcja Tabeli 3 i Rys. 3 z artykułu Yuan et al. (2023).

Zgodnie z artykułem (sekcja "Performance Analysis"):

    COSTvoter (HED) = 4 * Coste * Nc
    COSTCC    (HED) = 2 * Coste * Nc * Nv
    COSTvoter (HSE) = 9 * Coste * Nc
    COSTCC    (HSE) = 6 * Coste * Nc * Nv  +  4 * Coste

gdzie Coste to czas pojedynczej operacji "decentralizowanego mnożenia"
(w praktyce: potęgowania modulo na Z*_{n_1^2}; w artykule autorzy
mierzą Coste = 0.0068 s na CPU i7-12700H).

Eksperyment składa się z dwóch części:

A) Walidacja end-to-end dla N_v ∈ {100, 200, 500} — uruchamiamy pełną
   symulację HED-Voting z 5 kandydatami, sprawdzamy poprawność wyników
   i mierzymy *rzeczywisty* czas pracy CC. Następnie obliczamy
   eksperymentalne Coste i porównujemy z modelem teoretycznym.

B) Ekstrapolacja teoretyczna dla N_v ∈ {1000, 2000, 4000, 7000, 10000}
   (te same wartości co w artykule), dla 5 kandydatów. Wykres porównawczy
   HED-Voting vs HSE-Voting.

Wynik zapisywany jest do JSON-a i wykresu PNG (zgodnego z Rys. 3).
"""
import os, sys, time, json, random
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import gmpy2
from gmpy2 import mpz, powmod

from hed_voting import protocol as H
from hed_voting import rsa_sig  as R
from hed_voting import paillier as P


# ---------------------------------------------------------------------------
# Parametry zgodne z artykułem
# ---------------------------------------------------------------------------
N_CANDIDATES = 5
N_VOTERS_PAPER = [1000, 2000, 4000, 7000, 10000]   # wartości z Rys. 3
N_VOTERS_SIM   = [100, 200, 500]                   # rzeczywista symulacja


# ---------------------------------------------------------------------------
# Pomiar Coste — czas pojedynczego potęgowania modulo na Z*_{n_1^2}
# ---------------------------------------------------------------------------
def measure_coste(repeat: int = 2000) -> float:
    php, _ = P.setup_he(1024)
    n2 = php.n2
    base = mpz(2) ** 7 + 5
    exp  = (mpz(0xDEADBEEFCAFE) << 1000) | 1
    for _ in range(50): _ = powmod(base, exp, n2)   # rozgrzewka
    t0 = time.perf_counter()
    for _ in range(repeat):
        _ = powmod(base, exp, n2)
    return (time.perf_counter() - t0) / repeat


# ---------------------------------------------------------------------------
# Symulacja "end-to-end" — rzeczywisty pomiar pracy CC dla N_v wyborców
# ---------------------------------------------------------------------------
def simulate_real(n_voters: int, n_candidates: int, seed: int = 42) -> dict:
    """
    Pełna symulacja HED-Voting:
        - Każdy wyborca wysyła n_candidates ciphertextów (zgodnie z modelem
          z artykułu: po jednym ciphertekście dla każdego kandydata,
          z wartością 1 dla wybranego, 0 dla pozostałych).
        - CC dla każdego ciphertextu wykonuje Dec_CC + Rsa_Ver + Count_CC.
    Zwraca słownik z czasami i wynikami.
    """
    random.seed(seed)
    sys_ = H.setup_system(n_candidates=n_candidates,
                          n_voters_max=n_voters,
                          paillier_bits=1024)

    # losujemy preferencje wyborców
    expected = [0] * n_candidates
    choices  = []
    for _ in range(n_voters):
        c = random.randrange(n_candidates)
        choices.append(c)
        expected[c] += 1

    # generujemy klucze RSA wyborców
    voters = [R.user_keygen(1024) for _ in range(n_voters)]

    # FAZA WYBORCY: każdy generuje n_candidates ciphertextów (po 1 na kandydata)
    t0 = time.perf_counter()
    bulletins = []   # list of (upub, [(C, s), (C, s), ...])
    for i in range(n_voters):
        upub, upriv = voters[i]
        choice      = choices[i]
        per_cand = []
        for k in range(n_candidates):
            value = 1 if k == choice else 0
            # Enc_HE
            c1 = P.enc_he(sys_.PHE_pub, value)
            # Rsa_Sig
            sig = R.rsa_sig(upriv, int(c1))
            # Enc_CC
            from hed_voting import elgamal as E
            C = E.enc_cc(sys_.CC_pub, int(c1))
            per_cand.append((C, sig))
        bulletins.append((upub, per_cand))
    t_voters = time.perf_counter() - t0

    # FAZA CC: zliczanie
    from hed_voting import elgamal as E
    totals = [None] * n_candidates
    t0 = time.perf_counter()
    for upub, per_cand in bulletins:
        for k, (C, sig) in enumerate(per_cand):
            c1 = int(E.dec_cc(sys_.CC_priv, C))
            if not R.rsa_ver(upub, c1, sig):
                raise RuntimeError("Niepoprawny podpis!")
            if totals[k] is None:
                totals[k] = mpz(c1)
            else:
                totals[k] = P.add_he(sys_.PHE_pub, totals[k], mpz(c1))
    t_cc = time.perf_counter() - t0

    # Deszyfracja wyników
    results = [P.dec_he(sys_.PHE_priv, t) for t in totals]
    assert results == expected, f"Wyniki niezgodne: {results} vs {expected}"

    return {
        "n_voters":     n_voters,
        "n_candidates": n_candidates,
        "t_voters_s":   t_voters,                 # łączny czas wszystkich wyborców
        "t_voters_per": t_voters / n_voters,      # czas na 1 wyborcę
        "t_cc_s":       t_cc,                     # łączny czas CC
        "t_cc_per":     t_cc / n_voters,          # czas na 1 wyborcę w CC
        "results":      results,
        "expected":     expected,
    }


# ---------------------------------------------------------------------------
# Główny eksperyment
# ---------------------------------------------------------------------------
def main():
    out_dir = os.path.join(os.path.dirname(__file__), "results")
    os.makedirs(out_dir, exist_ok=True)

    print("=== Eksperyment 2: porównanie HED-Voting vs HSE-Voting ===")
    print(f"Liczba kandydatów: {N_CANDIDATES}")

    # ---- A) Pomiar Coste ------------------------------------------------
    print("\n[A] Pomiar Coste (czas modulo-power na Z*_{n^2}, n = 1024 b):")
    coste = measure_coste(repeat=2000)
    print(f"   Coste = {coste*1000:.4f} ms = {coste:.6f} s")
    print(f"   Coste z artykułu (i7-12700H) = 6.8000 ms = 0.006800 s")

    # ---- B) Walidacja end-to-end ---------------------------------------
    print("\n[B] Symulacja end-to-end (sprawdza poprawność i mierzy rzeczywisty czas):")
    real_runs = []
    for nv in N_VOTERS_SIM:
        print(f"   --- N_v = {nv} ---")
        res = simulate_real(nv, N_CANDIDATES)
        real_runs.append(res)
        print(f"      Wyniki:    {res['results']}")
        print(f"      Oczekiwane:{res['expected']}")
        print(f"      Czas CC:   {res['t_cc_s']:.3f} s   "
              f"(per wyborca: {res['t_cc_per']*1000:.2f} ms)")
        # Porównanie z modelem teoretycznym
        model = 2 * coste * N_CANDIDATES * nv
        print(f"      Model teoretyczny CostCC = 2*Coste*Nc*Nv "
              f"= {model:.3f} s  (błąd: {(res['t_cc_s']-model)/model*100:+.1f}%)")

    # Wyznaczamy "Coste empiryczne" z dopasowania t_cc do 2*Coste*Nc*Nv
    # (uśredniamy względem trzech symulowanych N_v)
    cost_emp = sum(r["t_cc_s"] / (2 * N_CANDIDATES * r["n_voters"])
                   for r in real_runs) / len(real_runs)
    print(f"\n   Coste empiryczne (z dopasowania CostCC) = {cost_emp*1000:.4f} ms")

    # ---- C) Ekstrapolacja na N_v z artykułu (Rys. 3) -------------------
    print("\n[C] Wartości teoretyczne dla N_v z artykułu:")
    print(f"   Używamy Coste = {coste*1000:.4f} ms (zmierzone na tej maszynie).")
    print()
    print(f"   {'N_v':>6} | {'HED CostCC [s]':>14} | {'HSE CostCC [s]':>14} | "
          f"{'HED CostVoter [ms]':>18} | {'HSE CostVoter [ms]':>18} | speedup")
    print("   " + "-"*100)

    table = []
    for nv in N_VOTERS_PAPER:
        cc_hed  = 2 * coste * N_CANDIDATES * nv
        cc_hse  = 6 * coste * N_CANDIDATES * nv + 4 * coste
        v_hed   = 4 * coste * N_CANDIDATES * 1000     # w ms (1 wyborca)
        v_hse   = 9 * coste * N_CANDIDATES * 1000
        speedup = cc_hse / cc_hed
        print(f"   {nv:>6} | {cc_hed:>14.3f} | {cc_hse:>14.3f} | "
              f"{v_hed:>18.3f} | {v_hse:>18.3f} | x{speedup:.3f}")
        table.append({
            "n_voters":      nv,
            "cc_hed_s":      cc_hed,
            "cc_hse_s":      cc_hse,
            "voter_hed_ms":  v_hed,
            "voter_hse_ms":  v_hse,
            "speedup":       speedup,
        })

    saving_pct = (1 - 2/6) * 100
    print(f"\n   Względne przyspieszenie HED/HSE w CC: stałe = "
          f"{table[0]['cc_hse_s']/table[0]['cc_hed_s']:.3f}x  "
          f"(={saving_pct:.1f}% mniej czasu — zgadza się z artykułem '~66.7%').")

    # ---- D) Wykres -----------------------------------------------------
    fig, ax = plt.subplots(figsize=(8, 5))
    nvs   = [r["n_voters"] for r in table]
    hed_s = [r["cc_hed_s"] for r in table]
    hse_s = [r["cc_hse_s"] for r in table]

    ax.plot(nvs, hed_s, "o-", color="#2A6FB7", label="HED-Voting (this work)")
    ax.plot(nvs, hse_s, "s-", color="#C0392B", label="HSE-Voting (Fan et al., 2020)")

    # Punkty empiryczne (z naszych symulacji end-to-end)
    sim_nvs = [r["n_voters"]  for r in real_runs]
    sim_t   = [r["t_cc_s"]    for r in real_runs]
    ax.plot(sim_nvs, sim_t, "x", color="#1F6E20", markersize=10,
            label="HED-Voting — pomiar empiryczny")

    ax.set_xlabel("Liczba wyborców  (N_v)")
    ax.set_ylabel("Czas pracy Counting Center  [s]")
    ax.set_title(f"Czas pracy CC w zależności od liczby wyborców  "
                 f"(N_c = {N_CANDIDATES},  Coste ≈ {coste*1000:.2f} ms)")
    ax.grid(True, linestyle=":", alpha=0.7)
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "exp2_cc_time.png"), dpi=130)
    plt.close(fig)
    print(f"\nZapisano wykres: {os.path.join(out_dir, 'exp2_cc_time.png')}")

    # Zapis JSON-a
    out = {
        "n_candidates":  N_CANDIDATES,
        "coste_s":       coste,
        "coste_paper_s": 0.0068,
        "coste_empirical_s": cost_emp,
        "real_runs":     real_runs,
        "extrapolation": table,
    }
    with open(os.path.join(out_dir, "exp2_cc_time.json"), "w") as f:
        json.dump(out, f, indent=2)
    print(f"Zapisano JSON:  {os.path.join(out_dir, 'exp2_cc_time.json')}")

if __name__ == "__main__":
    main()
