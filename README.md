# HED-Voting — implementacja schematu z Yuan i in. (2023)

Projekt na temat: **„Bezpieczne głosowanie elektroniczne: Implementacja systemu
głosowania, w którym głosy są szyfrowane homomorficznie, a wynik wyborów obliczany
jest bez poznania indywidualnych głosów"** — z wykorzystaniem szyfrowania
homomorficznego Pailliera.

Reprodukcja artykułu:
> Yuan K., Sang P., Zhang S., Chen X., Yang W., Jia C. (2023).
> „An electronic voting scheme based on homomorphic encryption and decentralization",
> *PeerJ Computer Science* 9:e1649. DOI: 10.7717/peerj-cs.1649.

## Wymagania

- Python 3.10+ (testowane na 3.12)
- Zależności:
  ```bash
  pip install gmpy2 phe matplotlib
  ```

## Struktura pakietu

```
hed_voting/
├── __init__.py
├── paillier.py            # Setup_HE, Enc_HE, Count_CC (homomorf.), Dec_HE
├── elgamal.py             # Setup_CC, Enc_CC, Dec_CC  (z optymalizacją sieve)
├── rsa_sig.py             # User_KeyGen, Rsa_Sig, Rsa_Ver
├── protocol.py            # pełny protokół (voter_cast, cc_process, reveal)
│
├── tests/
│   └── test_protocol.py   # test poprawności end-to-end (50 wyborców, 5 kand.)
│
└── experiments/
    ├── exp1_basic_ops.py        # pomiar Mul1/Mul2/Pow/Inv/Key (Tabela 4-5)
    ├── exp2_compare.py          # porównanie HED vs HSE (Tabela 3, Rys. 3)
    ├── exp3a_phase_timing.py    # czas każdego z 10 algorytmów
    ├── exp3b_setup_cc.py        # osobno Setup_CC (długie)
    ├── exp3c_plot.py            # wykres słupkowy faz protokołu
    └── results/                 # JSON-y i wykresy z eksperymentów
```

## Uruchomienie

Z poziomu tego katalogu:

```bash
# Test poprawności pełnego protokołu 
PYTHONPATH=. python3 hed_voting/tests/test_protocol.py

# Eksperyment 1: czas operacji bazowych 
PYTHONPATH=. python3 hed_voting/experiments/exp1_basic_ops.py

# Eksperyment 2: porównanie HED-Voting vs HSE-Voting 
PYTHONPATH=. python3 hed_voting/experiments/exp2_compare.py

# Eksperyment 3a: czas każdej fazy protokołu 
PYTHONPATH=. python3 hed_voting/experiments/exp3a_phase_timing.py

# Eksperyment 3b: Setup_CC 
PYTHONPATH=. python3 hed_voting/experiments/exp3b_setup_cc.py

# Wykres słupkowy faz:
PYTHONPATH=. python3 hed_voting/experiments/exp3c_plot.py
```

Wyniki (JSON-y i pliki PNG) trafiają do `hed_voting/experiments/results/`.

## Kluczowe wyniki

- **Test poprawności**: 100 % zgodność odszyfrowanych wyników z rzeczywistym
  rozkładem głosów (testowano dla 50, 100, 200, 500 wyborców i 5 kandydatów).
- **Przyspieszenie HED-Voting względem HSE-Voting**: dokładnie ×3 (oszczędność
  czasu w CC = 66,7 %), zgodnie z deklaracją artykułu.
- **Złożoność pojedynczych operacji**: Setup_HE ≈ 7 ms, Enc_HE ≈ 1,6 ms,
  Dec_CC ≈ 3,6 ms, Rsa_Ver ≈ 9,5 µs, Count_CC ≈ 2,3 µs, Setup_CC ≈ 50 s
  (jednorazowo).
