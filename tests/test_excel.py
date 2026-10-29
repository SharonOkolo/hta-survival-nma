import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "excel"))
import build_model as b


def test_excel_formulas_match_r():
    assert len(b.check(b.load_params())) == 3


def test_survival_starts_at_one_and_falls():
    for m in b.MODELS:
        s = [b.surv(m, 1.2, 8.0, t) for t in (0, 1, 5, 20)]
        assert s[0] == 1 and all(x > y for x, y in zip(s, s[1:]))
