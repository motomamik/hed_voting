"""Generuje wykres słupkowy rozkładu czasu faz protokołu HED-Voting."""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

out_dir = os.path.join(os.path.dirname(__file__), "results")
with open(os.path.join(out_dir, "exp3_phase_timing.json")) as f:
    data = json.load(f)

# Ułóżmy w kolejności występowania w protokole
order = ["Setup_HE", "Setup_CC", "User_KeyGen",
         "Enc_HE", "Rsa_Sig", "Enc_CC",
         "Dec_CC", "Rsa_Ver", "Count_CC",
         "Dec_HE"]
rows = {r["name"]: r for r in data["rows"]}
names  = order
medians = [rows[n]["median_ms"] for n in names]

# Klasyfikacja: pomarańczowy = jednorazowa konfiguracja, niebieski = wyborca, zielony = CC
colors = {
    "Setup_HE": "#E67E22", "Setup_CC": "#E67E22", "User_KeyGen": "#E67E22",
    "Enc_HE":   "#2A6FB7", "Rsa_Sig":  "#2A6FB7", "Enc_CC":      "#2A6FB7",
    "Dec_CC":   "#1F6E20", "Rsa_Ver":  "#1F6E20", "Count_CC":    "#1F6E20",
    "Dec_HE":   "#8E44AD",
}

fig, ax = plt.subplots(figsize=(10, 5))
bars = ax.bar(names, medians, color=[colors[n] for n in names])
ax.set_yscale("log")
ax.set_ylabel("Czas wykonania [ms]  (skala log.)")
ax.set_title("Czas wykonania pojedynczej operacji każdego algorytmu HED-Voting\n"
             "(Paillier 1024 b, ElGamal 2050 b safe-prime, RSA 1024 b)")
ax.grid(True, axis="y", linestyle=":", alpha=0.6)

# Podpisy nad słupkami
for bar, m in zip(bars, medians):
    if m < 0.01:
        label = f"{m*1000:.1f} µs"
    elif m < 1:
        label = f"{m:.2f} ms"
    elif m < 1000:
        label = f"{m:.1f} ms"
    else:
        label = f"{m/1000:.1f} s"
    ax.text(bar.get_x() + bar.get_width()/2, m * 1.15, label,
            ha="center", va="bottom", fontsize=8)

# Legenda
from matplotlib.patches import Patch
legend = [
    Patch(color="#E67E22", label="Jednorazowa konfiguracja"),
    Patch(color="#2A6FB7", label="Operacje wyborcy"),
    Patch(color="#1F6E20", label="Operacje CC"),
    Patch(color="#8E44AD", label="Ogłoszenie wyniku"),
]
ax.legend(handles=legend, loc="upper right")

plt.xticks(rotation=20, ha="right")
fig.tight_layout()
fig.savefig(os.path.join(out_dir, "exp3_phase_timing.png"), dpi=130)
plt.close(fig)
print(f"Zapisano: {os.path.join(out_dir, 'exp3_phase_timing.png')}")
