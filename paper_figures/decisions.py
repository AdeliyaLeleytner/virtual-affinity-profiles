"""Loaders for the measured decision outcomes behind the manuscript figures.

Reads the tables published with the article. Every figure value comes from
here; none is transcribed by hand.

  data/decisions_paired_uncertainty.csv   Additional file 5
  data/acquisition_curves.csv             Additional file 4
"""
import collections
import csv
import os

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")

DECISIONS = os.path.join(DATA, "decisions_paired_uncertainty.csv")
ACQUISITION = os.path.join(DATA, "acquisition_curves.csv")

HELD_OUT_PROTEIN = "new_proteins_known_molecules"
NEW_MOLECULE = "new_molecules_known_proteins"
EXACT = "exact_four_cell"
CENSORED = "censor_involving_identified_D"

# The eight scores, in the order used in the manuscript tables.
NAME = {
    "raw_ad4": "AutoDock4",
    "raw_gnina_own": "GNINA (own search)",
    "vina": "Vina",
    "gnina": "GNINA (Vina pose)",
    "boltz_affinity": "Boltz-2",
    "boltzina_vina_zero_affinity": "Boltzina (zero)",
    "boltzina_vina_far_affinity": "Boltzina (far)",
    "nesso_affinity": "Nesso-1",
}


def trade_offs(regime, model, stratum=EXACT, metric="pareto_correct_fraction"):
    """{channel: dict} of percentages for one regime and one model.

    Keys: base (the chemistry-and-level baseline), est (this model), d (paired
    difference), lo/hi (2.5th and 97.5th percentiles of 400 chemical-group
    resampling draws). All in percent.
    """
    out = {}
    with open(DECISIONS, encoding="utf-8", errors="replace") as fh:
        for x in csv.DictReader(fh):
            if (x["regime"] == regime and x["model"] == model
                    and x["metric"] == metric
                    and x["scope"] == "native_genes"
                    and x["weighting"] == "actual_assay_pair_equal"
                    and x["measurement_stratum"] == stratum
                    and x["channel"] in NAME and x["channel"] not in out):
                out[x["channel"]] = dict(
                    base=100 * float(x["baseline_estimate"]),
                    est=100 * float(x["estimate"]),
                    d=100 * float(x["arm_minus_baseline"]),
                    lo=100 * float(x["paired_delta_q025"]),
                    hi=100 * float(x["paired_delta_q975"]))
    return out


def growth(dataset, model, metric="pareto_correct_fraction"):
    """{acquisition order: {k: percentage}} for one panel and one model.

    Restricted to the cohort that matches the support of the main decision
    analysis. A second, smaller metadata-matched cohort exists in the same
    table; the two are never averaged.
    """
    out = collections.defaultdict(dict)
    with open(ACQUISITION, encoding="utf-8", errors="replace") as fh:
        for x in csv.DictReader(fh):
            if (x["table_kind"] == "all_growth_preferences"
                    and x["dataset"] == dataset and x["model"] == model
                    and x["weighting"] == "actual_assay_pair_equal"
                    and x["domain"] == "all_source_operational"
                    and x[metric]):
                out[x["order"]][int(x["k"])] = 100 * float(x[metric])
    return out
