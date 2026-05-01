"""Verify the L formula handles ties correctly.

Constructs a sim buffer with deliberate exact ties and a query at the
tied value.  Compares:

  buggy:   (searchsorted_right - 0.5) / M
  fixed:   0.5 * (searchsorted_left + searchsorted_right) / M
  chapter: (#more)/M + 0.5 * (#equal)/M

Also re-runs the continuous Gaussian demo from luck_fit.py to confirm
results are essentially unchanged when there are no ties.
"""

import numpy as np

nlp_sim = np.array([1.0, 1.0, 1.0, 2.0, 2.0, 3.0, 4.0, 5.0])
M = nlp_sim.size

queries = [0.5, 1.0, 2.0, 3.0, 6.0]

def buggy(q):
    idx = np.searchsorted(nlp_sim, q, side="right")
    return (idx - 0.5) / M

def fixed(q):
    lo = np.searchsorted(nlp_sim, q, side="left")
    hi = np.searchsorted(nlp_sim, q, side="right")
    return 0.5 * (lo + hi) / M

def chapter(q):
    n_more = int(np.sum(nlp_sim < q))
    n_equal = int(np.sum(nlp_sim == q))
    return (n_more + 0.5 * n_equal) / M

print(f"sim nlp = {nlp_sim.tolist()}  (M={M})")
print(f"  {'q':>4}  {'buggy':>7}  {'fixed':>7}  {'chapter':>7}  match?")
for q in queries:
    b, f, c = buggy(q), fixed(q), chapter(q)
    print(f"  {q:>4}  {b:>7.4f}  {f:>7.4f}  {c:>7.4f}  "
          f"{'OK' if abs(f - c) < 1e-12 else 'FAIL'}")

print("\n--- continuous case (Gaussian, no ties): re-running luck_fit demo ---\n")
import subprocess, sys, pathlib
here = pathlib.Path(__file__).resolve().parent
result = subprocess.run([sys.executable, str(here / "luck_fit.py")],
                        capture_output=True, text=True)
print(result.stdout)
