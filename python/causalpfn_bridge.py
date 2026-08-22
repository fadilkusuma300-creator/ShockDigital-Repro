#!/usr/bin/env python3
"""CausalPFN command-line bridge used by the R validation workflow.

The R side supplies an already numeric model matrix. Query rows can carry a
`grid_id`; each grid is evaluated separately so standardized-curve uncertainty
can be requested without constructing an excessively large sampling tensor.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from causalpfn import CATEEstimator


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--train", required=True)
    parser.add_argument("--query", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--summary", required=False)
    parser.add_argument("--device", default="auto")
    parser.add_argument("--seed", type=int, default=20260820)
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--ci-samples", type=int, default=1000)
    parser.add_argument("--with-ci", action="store_true")
    return parser.parse_args()


def choose_device(value: str) -> torch.device:
    if value == "auto":
        return torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    return torch.device(value)


def ate_interval(estimator: CATEEstimator, x: np.ndarray, alpha: float, n_samples: int):
    """Return the model's interval for the mean CATE over x.

    Current package versions expose estimate_ate_CI. A fallback to the internal
    helper is retained because some releases contain a return-key mismatch in
    the public wrapper while the underlying interval calculation is available.
    """
    try:
        ci = estimator.estimate_ate_CI(x, alpha=alpha, n_samples=n_samples)
        return float(np.asarray(ci["lower_bound"]).reshape(-1)[0]), float(
            np.asarray(ci["upper_bound"]).reshape(-1)[0]
        )
    except (KeyError, AttributeError):
        ci = estimator._estimate_ate_cate_CI(x, alpha=alpha, n_samples=n_samples)
        return float(np.asarray(ci["ate_lower_bound"]).reshape(-1)[0]), float(
            np.asarray(ci["ate_upper_bound"]).reshape(-1)[0]
        )


def main() -> None:
    args = parse_args()
    np.random.seed(args.seed)
    torch.manual_seed(args.seed)

    train = pd.read_csv(args.train)
    query = pd.read_csv(args.query)
    reserved = {"treatment", "outcome", "grid_id", "analysis_weight"}
    features = [c for c in train.columns if c not in reserved]
    if not features:
        raise ValueError("No feature columns found in training file")
    missing = [c for c in features if c not in query.columns]
    if missing:
        raise ValueError(f"Query file is missing feature columns: {missing}")

    x_train = train[features].to_numpy(dtype=np.float32)
    t_train = train["treatment"].to_numpy(dtype=np.float32)
    y_train = train["outcome"].to_numpy(dtype=np.float32)

    device = choose_device(args.device)
    estimator = CATEEstimator(device=device, verbose=False)
    estimator.fit(x_train, t_train, y_train)

    if "grid_id" in query.columns:
        groups = list(query.groupby("grid_id", sort=True))
    else:
        groups = [(0, query)]

    row_outputs = []
    summaries = []
    for grid_id, group in groups:
        x_query = group[features].to_numpy(dtype=np.float32)
        cate = np.asarray(estimator.estimate_cate(x_query), dtype=float).reshape(-1)
        w = (
            group["analysis_weight"].to_numpy(dtype=float)
            if "analysis_weight" in group.columns
            else np.ones(len(group), dtype=float)
        )
        w = np.where(np.isfinite(w) & (w > 0), w, 0.0)
        if w.sum() <= 0:
            w = np.ones(len(group), dtype=float)
        w = w / w.sum()
        estimate = float(np.sum(w * cate))

        lower = upper = np.nan
        if args.with_ci:
            # R normally supplies a probability-weighted resample with equal
            # query weights, so the native mean-CATE interval targets the same
            # survey-weighted reference distribution up to Monte Carlo error.
            lower, upper = ate_interval(estimator, x_query, args.alpha, args.ci_samples)

        block = pd.DataFrame(
            {
                "grid_id": grid_id,
                "analysis_weight": w,
                "cate": cate,
                "lower": lower,
                "upper": upper,
            }
        )
        row_outputs.append(block)
        summaries.append(
            {
                "grid_id": grid_id,
                "estimate": estimate,
                "lower": lower,
                "upper": upper,
            }
        )

    output = pd.concat(row_outputs, ignore_index=True)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output, index=False)
    if args.summary:
        Path(args.summary).parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(summaries).to_csv(args.summary, index=False)

    metadata = {
        "features": features,
        "n_train": int(len(train)),
        "n_query": int(len(query)),
        "device": str(device),
        "with_ci": bool(args.with_ci),
        "ci_samples": int(args.ci_samples),
    }
    Path(args.output).with_suffix(".json").write_text(
        json.dumps(metadata, indent=2), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
