# Data

Both datasets ship inside R packages, so nothing is downloaded or redistributed here.

| Dataset | Package | Source | Content |
|---|---|---|---|
| `colon` | `survival` (LGPL ≥ 2) | Moertel CG et al. *N Engl J Med* 1990; Laurie JA et al. *J Clin Oncol* 1989 | 929 patients with resected stage C colon cancer randomised to observation, levamisole, or levamisole + 5-FU. Two records per patient (recurrence, death); this project uses death (`etype == 2`). |
| `smokingcessation` | `netmeta` (GPL ≥ 2) | Hasselblad V. *Med Decis Making* 1998; used in NICE DSU TSD 2 (Dias et al. 2011) | 24 trials, 4 interventions (no contact, self-help, individual counselling, group counselling); number who quit / number randomised per arm. |

Load them with `survival::colon` and `data(smokingcessation, package = "netmeta")`.
