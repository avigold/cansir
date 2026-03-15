"""
ML model to rank genes by cancer-protective potential.

Uses a semi-supervised approach:
1. Train on labeled genes (known protective vs. known risk)
2. Predict scores for unlabeled genes
3. Rank all genes by predicted protective probability
"""

from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier
from sklearn.metrics import (
    classification_report,
    roc_auc_score,
)
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.preprocessing import StandardScaler


def load_data(path: str = "data/processed/feature_matrix.csv"):
    """Load feature matrix with labels."""
    df = pd.read_csv(path)
    return df


def train_and_predict(
    df: pd.DataFrame,
    feature_cols: list[str],
    n_iterations: int = 3,
) -> pd.DataFrame:
    """
    Semi-supervised learning loop:
    1. Train on labeled data
    2. Predict unlabeled
    3. Add high-confidence predictions to training set
    4. Repeat

    Returns DataFrame with protection_probability for all genes.
    """
    results = df.copy()
    labels = df["label"].copy()

    scaler = StandardScaler()
    X_all = scaler.fit_transform(df[feature_cols].values)

    current_labels = labels.copy()

    for iteration in range(n_iterations):
        # Split labeled vs unlabeled
        labeled_mask = current_labels >= 0
        X_train = X_all[labeled_mask]
        y_train = current_labels[labeled_mask].values

        if len(np.unique(y_train)) < 2:
            print(f"Iteration {iteration}: Not enough classes to train. "
                  "Need both positive and negative examples.")
            break

        n_pos = (y_train == 1).sum()
        n_neg = (y_train == 0).sum()
        print(f"\nIteration {iteration + 1}: "
              f"{n_pos} positive, {n_neg} negative, "
              f"{(~labeled_mask).sum()} unlabeled")

        # Train ensemble
        rf = RandomForestClassifier(
            n_estimators=200,
            max_depth=10,
            class_weight="balanced",
            random_state=42 + iteration,
        )
        gb = GradientBoostingClassifier(
            n_estimators=200,
            max_depth=5,
            learning_rate=0.05,
            random_state=42 + iteration,
        )

        rf.fit(X_train, y_train)
        gb.fit(X_train, y_train)

        # Cross-validation score on labeled data
        if labeled_mask.sum() >= 10:
            n_splits = min(5, min(n_pos, n_neg))
            if n_splits >= 2:
                cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)
                rf_scores = cross_val_score(rf, X_train, y_train, cv=cv, scoring="roc_auc")
                print(f"  RF CV AUC: {rf_scores.mean():.3f} (+/- {rf_scores.std():.3f})")

        # Predict on all data (ensemble average)
        rf_proba = rf.predict_proba(X_all)[:, 1]
        gb_proba = gb.predict_proba(X_all)[:, 1]
        ensemble_proba = (rf_proba + gb_proba) / 2

        # Self-training: add high-confidence predictions
        unlabeled_mask = current_labels == -1
        high_conf_pos = unlabeled_mask & (ensemble_proba > 0.85)
        high_conf_neg = unlabeled_mask & (ensemble_proba < 0.15)

        n_added = high_conf_pos.sum() + high_conf_neg.sum()
        current_labels[high_conf_pos] = 1
        current_labels[high_conf_neg] = 0
        print(f"  Added {high_conf_pos.sum()} positive, "
              f"{high_conf_neg.sum()} negative pseudo-labels")

        if n_added == 0:
            print("  No new pseudo-labels added. Stopping.")
            break

    # Final predictions
    results["protection_probability"] = ensemble_proba
    results["predicted_label"] = (ensemble_proba > 0.5).astype(int)

    # Feature importance (from Random Forest)
    importance = pd.DataFrame({
        "feature": feature_cols,
        "importance": rf.feature_importances_,
    }).sort_values("importance", ascending=False)

    return results, importance


def evaluate_on_known_genes(results: pd.DataFrame):
    """
    Sanity check: do known protective/tumor suppressor genes rank high?
    """
    known_protective = [
        "TP53", "BRCA1", "BRCA2", "APC", "RB1", "PTEN", "MLH1", "MSH2",
        "ATM", "CHEK2", "CDH1", "PALB2", "RAD51", "NBN", "MRE11",
    ]

    found = results[results["gene"].isin(known_protective)]
    if not found.empty:
        print("\n=== Known cancer-related genes in results ===")
        print(found[["gene", "protection_probability", "predicted_label",
                      "label"]].to_string(index=False))


def save_results(
    results: pd.DataFrame,
    importance: pd.DataFrame,
    output_dir: str = "data/processed",
):
    """Save ranked gene list and feature importance."""
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    # Ranked gene list
    ranked = results.sort_values("protection_probability", ascending=False)
    ranked.to_csv(out_path / "ranked_protective_genes.csv", index=False)

    # Feature importance
    importance.to_csv(out_path / "feature_importance.csv", index=False)

    print(f"\nResults saved to {out_path}")
    print(f"\nTop 30 predicted cancer-protective genes:")
    top_cols = ["gene", "protection_probability", "n_cancer_types",
                "mean_or", "in_dna_repair", "in_immune"]
    available = [c for c in top_cols if c in ranked.columns]
    print(ranked[available].head(30).to_string(index=False))

    print(f"\nFeature importance:")
    print(importance.to_string(index=False))


if __name__ == "__main__":
    print("=== Cancer-Protective Gene Classifier ===\n")

    df = load_data()

    # Identify feature columns (exclude metadata and label)
    meta_cols = ["gene", "ensembl_id", "approved_name", "cancer_types",
                 "variants", "label"]
    feature_cols = [c for c in df.columns
                    if c not in meta_cols and df[c].dtype in [np.float64, np.int64]]

    print(f"Features ({len(feature_cols)}): {feature_cols}")

    results, importance = train_and_predict(df, feature_cols)
    evaluate_on_known_genes(results)
    save_results(results, importance)
