"""
=============================================================
  PLANT STRESS CLASSIFIER  —  Deep Learning Pipeline v2
  ─────────────────────────────────────────────────────────
  • Deep MLP (128-128-64) neural network
  • 5-Fold Stratified Cross-Validation
  • Integrated Gradients → gene-level attribution
  • Feature rank table:
      gene | mean_shap | shap_expr_corr | association |
      shap_fc_concordant | rank

    pip install pandas numpy scikit-learn matplotlib seaborn
  ─────────────────────────────────────────────────────────
  Outputs
  -------
    gene_feature_ranks.csv       — full gene rank table
    plant_stress_cv_results.png  — CV + feature plots
    plant_stress_model.pt        — saved model (PyTorch)
          OR
    plant_stress_model.pkl       — saved model (sklearn fallback)
=============================================================
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch
import warnings, os
warnings.filterwarnings("ignore")

from sklearn.model_selection import StratifiedKFold, train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.metrics import (roc_auc_score, confusion_matrix, roc_curve,
                             ConfusionMatrixDisplay, classification_report)
from scipy.stats import pearsonr

# ─────────────────────────────────────────────────────────────
#  CONFIG
# ─────────────────────────────────────────────────────────────
SEED      = 42
N_FOLDS   =10
N_PCA     = 50        # PCA components (~90-95% variance)
IG_STEPS  = 50        # Integrated Gradients interpolation steps
np.random.seed(SEED)

# Detect backend
try:
    import torch
    import torch.nn as nn
    import torch.optim as optim
    from torch.utils.data import DataLoader, TensorDataset
    BACKEND = "pytorch"
    DEVICE  = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    torch.manual_seed(SEED)
    print(f"Backend : PyTorch  |  Device : {DEVICE}")
except ImportError:
    from sklearn.neural_network import MLPClassifier
    import pickle
    BACKEND = "sklearn"
    print("Backend : scikit-learn MLPClassifier (install PyTorch for full ResNet)")
print()

# ─────────────────────────────────────────────────────────────
#  1.  LOAD DATA
# ─────────────────────────────────────────────────────────────
def load_data(csv_path: str):
    df = pd.read_csv(csv_path)
    # Drop ID + near-zero columns + label
    drop = ["sampleID", "Sobic.001G280100", "Sobic.001G279601",
            "Sobic.004G172800", "Sobic.008G015400", "treatment"]
    gene_cols = [c for c in df.columns if c not in drop]
    X  = df[gene_cols].values.astype(np.float32)
    y  = df["treatment"].values.astype(int)
    ids = df["sampleID"].values
    print(f"Samples : {X.shape[0]}  |  Genes : {X.shape[1]}")
    print(f"CK(0) = {(y==0).sum()}  |  Stressed(1) = {(y==1).sum()}\n")
    return X, y, ids, gene_cols


# ─────────────────────────────────────────────────────────────
#  2.  PYTORCH MODEL  (used when PyTorch is available)
# ─────────────────────────────────────────────────────────────
if BACKEND == "pytorch":

    class ResidualBlock(nn.Module):
        def __init__(self, dim, dropout=0.35):
            super().__init__()
            self.block = nn.Sequential(
                nn.Linear(dim, dim), nn.BatchNorm1d(dim),
                nn.ReLU(), nn.Dropout(dropout),
                nn.Linear(dim, dim), nn.BatchNorm1d(dim))
            self.relu = nn.ReLU()
        def forward(self, x):
            return self.relu(self.block(x) + x)

    class ResNet(nn.Module):
        """Deep Residual MLP — robust to vanishing gradients."""
        def __init__(self, in_dim, hidden=128, n_blocks=4, dropout=0.35):
            super().__init__()
            self.stem   = nn.Sequential(
                nn.Linear(in_dim, hidden), nn.BatchNorm1d(hidden), nn.ReLU())
            self.blocks = nn.Sequential(
                *[ResidualBlock(hidden, dropout) for _ in range(n_blocks)])
            self.head   = nn.Linear(hidden, 1)
        def forward(self, x):
            return self.head(self.blocks(self.stem(x))).squeeze(1)

    def build_model():
        return ResNet(in_dim=N_PCA, hidden=128, n_blocks=4, dropout=0.35)

    def train_fold(model, X_tr, y_tr, X_va, y_va,
                   epochs=200, lr=1e-3, patience=25, batch=16):
        model = model.to(DEVICE)
        crit  = nn.BCEWithLogitsLoss()
        opt   = optim.Adam(model.parameters(), lr=lr, weight_decay=1e-4)
        sch   = optim.lr_scheduler.ReduceLROnPlateau(opt, patience=10, factor=0.5)

        def mk(X, y):
            return DataLoader(
                TensorDataset(torch.tensor(X, dtype=torch.float32),
                              torch.tensor(y, dtype=torch.float32)),
                batch_size=batch, shuffle=True,
                drop_last=(len(X) % batch == 1))

        tr_ld = mk(X_tr, y_tr)
        va_ld = DataLoader(
            TensorDataset(torch.tensor(X_va, dtype=torch.float32),
                          torch.tensor(y_va, dtype=torch.float32)),
            batch_size=batch)

        best_val, best_state, wait = 1e9, None, 0
        for epoch in range(1, epochs + 1):
            model.train()
            for xb, yb in tr_ld:
                xb, yb = xb.to(DEVICE), yb.to(DEVICE)
                opt.zero_grad(); crit(model(xb), yb).backward(); opt.step()

            model.eval(); v_loss = 0
            with torch.no_grad():
                for xb, yb in va_ld:
                    xb, yb = xb.to(DEVICE), yb.to(DEVICE)
                    v_loss += crit(model(xb), yb).item() * len(xb)
            v_loss /= len(va_ld.dataset)
            sch.step(v_loss)

            if v_loss < best_val:
                best_val, best_state, wait = v_loss, model.state_dict(), 0
            else:
                wait += 1
                if wait >= patience:
                    break

        model.load_state_dict(best_state)
        return model

    def predict_proba(model, X_p):
        model.eval()
        with torch.no_grad():
            logits = model(torch.tensor(X_p, dtype=torch.float32).to(DEVICE)).cpu().numpy()
        return 1 / (1 + np.exp(-logits))

    def integrated_gradients(model, X_pca, steps=IG_STEPS):
        """PyTorch autograd IG — exact gradients, fast."""
        model.eval()
        X_t    = torch.tensor(X_pca, dtype=torch.float32)
        base_t = torch.zeros_like(X_t)
        alphas = torch.linspace(0, 1, steps).view(-1, 1, 1)
        interp = (base_t.unsqueeze(0) + alphas *
                  (X_t - base_t).unsqueeze(0)).to(DEVICE)
        interp.requires_grad_(True)
        grads = []
        for s in interp:
            g = torch.autograd.grad(model(s).sum(), s)[0]
            grads.append(g.detach().cpu())
        avg_grads = torch.stack(grads).mean(0)
        return ((X_t - base_t) * avg_grads).numpy()

else:
    # ── sklearn fallback ──────────────────────────────────────
    def build_model():
        return MLPClassifier(
            hidden_layer_sizes=(128, 128, 64),
            max_iter=300, random_state=SEED,
            early_stopping=True, validation_fraction=0.15,
            alpha=1e-4)

    def train_fold(model, X_tr, y_tr, X_va, y_va, **kw):
        model.fit(X_tr, y_tr)
        return model

    def predict_proba(model, X_p):
        return model.predict_proba(X_p)[:, 1]

    def integrated_gradients(model, X_pca, steps=IG_STEPS):
        """Numerical finite-difference IG for sklearn models."""
        baseline = np.zeros_like(X_pca)
        ig = np.zeros_like(X_pca)
        eps = 1e-4
        for alpha in np.linspace(0, 1, steps):
            interp = baseline + alpha * (X_pca - baseline)
            pb     = model.predict_proba(interp)[:, 1]
            for d in range(X_pca.shape[1]):
                Xp = interp.copy(); Xp[:, d] += eps
                ig[:, d] += (model.predict_proba(Xp)[:, 1] - pb) / eps
        ig /= steps
        ig *= (X_pca - baseline)
        return ig


# ─────────────────────────────────────────────────────────────
#  3.  FEATURE RANK TABLE
# ─────────────────────────────────────────────────────────────
def build_feature_table(all_ig_gene, all_X_gene, y_all,
                        gene_cols, X_full, y_full):
    """
    Returns DataFrame with columns:
      gene | mean_shap | shap_expr_corr | association |
      shap_fc_concordant | rank
    """
    n = len(gene_cols)
    mean_shap      = np.abs(all_ig_gene).mean(axis=0)
    shap_expr_corr = np.zeros(n)

    for g in range(n):
        ig_g, ex_g = all_ig_gene[:, g], all_X_gene[:, g]
        if ig_g.std() > 1e-9 and ex_g.std() > 1e-9:
            shap_expr_corr[g], _ = pearsonr(ig_g, ex_g)

    # association: positive corr → gene drives "Stressed" prediction
    association = np.where(shap_expr_corr > 0, "Stressed", "CK")

    # fold-change: mean(stressed) - mean(CK) in raw expression space
    fc             = X_full[y_full == 1].mean(0) - X_full[y_full == 0].mean(0)
    mean_ig_signed = all_ig_gene.mean(0)
    concordant     = (np.sign(mean_ig_signed) == np.sign(fc)).astype(int)

    df = pd.DataFrame({
        "gene"              : gene_cols,
        "mean_shap"         : mean_shap,
        "shap_expr_corr"    : shap_expr_corr,
        "association"       : association,
        "shap_fc_concordant": concordant,
    }).sort_values("mean_shap", ascending=False).reset_index(drop=True)
    df["rank"] = df.index + 1
    return df


# ─────────────────────────────────────────────────────────────
#  4.  PLOTTING
# ─────────────────────────────────────────────────────────────
def make_plots(fold_metrics, y_all, all_probs, feature_df):
    fig = plt.figure(figsize=(18, 12))
    fig.suptitle(
        "Plant Stress Classifier — K-Fold CV (Integrated Gradients Feature Ranking)",
        fontsize=13, fontweight="bold")
    gs = gridspec.GridSpec(2, 4, figure=fig, hspace=0.45, wspace=0.38)

    folds = [f"F{m['fold']}" for m in fold_metrics]
    accs  = [m["acc"] for m in fold_metrics]
    aucs  = [m["auc"] for m in fold_metrics]

    # Per-fold accuracy
    ax = fig.add_subplot(gs[0, 0])
    bars = ax.bar(folds, accs, color="#4CAF50", edgecolor="k", lw=0.7)
    ax.axhline(np.mean(accs), ls="--", color="gray",
               label=f"Mean={np.mean(accs):.2f}")
    ax.set_ylim(0, 1.12); ax.set_title("Per-Fold Accuracy"); ax.legend(fontsize=8)
    for b, a in zip(bars, accs):
        ax.text(b.get_x()+b.get_width()/2, b.get_height()+0.01,
                f"{a:.2f}", ha="center", fontsize=8)

    # Per-fold AUC
    ax = fig.add_subplot(gs[0, 1])
    bars = ax.bar(folds, aucs, color="#2196F3", edgecolor="k", lw=0.7)
    ax.axhline(np.mean(aucs), ls="--", color="gray",
               label=f"Mean={np.mean(aucs):.2f}")
    ax.set_ylim(0, 1.12); ax.set_title("Per-Fold ROC-AUC"); ax.legend(fontsize=8)
    for b, a in zip(bars, aucs):
        ax.text(b.get_x()+b.get_width()/2, b.get_height()+0.01,
                f"{a:.2f}", ha="center", fontsize=8)

    # Aggregate confusion matrix
    ax = fig.add_subplot(gs[0, 2])
    cm = confusion_matrix(y_all, all_probs >= 0.5)
    ConfusionMatrixDisplay(cm, display_labels=["CK", "Stressed"]).plot(
        ax=ax, colorbar=False)
    ax.set_title("Aggregate Confusion Matrix")

    # ROC curve
    ax = fig.add_subplot(gs[0, 3])
    fpr, tpr, _ = roc_curve(y_all, all_probs)
    auc = roc_auc_score(y_all, all_probs)
    ax.plot(fpr, tpr, color="#9C27B0", lw=2, label=f"AUC = {auc:.3f}")
    ax.plot([0,1], [0,1], "k--", alpha=0.4)
    ax.set_title("Aggregate ROC Curve"); ax.legend()
    ax.set_xlabel("FPR (1-specificity)"); ax.set_ylabel("TPR (sensitivity)")

    # Top-20 genes horizontal bar
    ax = fig.add_subplot(gs[1, :2])
    top20 = feature_df.head(20)
    colors = ["#FF5722" if a == "Stressed" else "#2196F3"
              for a in top20["association"]]
    ax.barh(top20["gene"][::-1], top20["mean_shap"][::-1],
            color=colors[::-1], edgecolor="k", lw=0.4)
    ax.set_xlabel("Mean |Attribution|  (Integrated Gradients)")
    ax.set_title("Top-20 Most Important Genes")
    ax.legend(handles=[Patch(color="#FF5722", label="→ Stressed"),
                        Patch(color="#2196F3", label="→ CK")], fontsize=8)

    # IG–expression correlation scatter (top 50)
    ax = fig.add_subplot(gs[1, 2])
    top50 = feature_df.head(50)
    sc = ["#FF5722" if a == "Stressed" else "#2196F3"
          for a in top50["association"]]
    ax.scatter(top50["shap_expr_corr"], top50["mean_shap"],
               c=sc, s=35, edgecolors="k", lw=0.3)
    ax.axvline(0, color="gray", ls="--", lw=0.8)
    ax.set_xlabel("IG–Expression Correlation  (r)")
    ax.set_ylabel("Mean |Attribution|")
    ax.set_title("Attribution vs Expression Corr\n(top-50 genes)")
    ax.legend(handles=[Patch(color="#FF5722", label="Stressed"),
                        Patch(color="#2196F3", label="CK")], fontsize=8)

    # Concordance pie
    ax = fig.add_subplot(gs[1, 3])
    nc = feature_df["shap_fc_concordant"].sum()
    nd = len(feature_df) - nc
    ax.pie([nc, nd],
           labels=[f"Concordant\n({nc})", f"Discordant\n({nd})"],
           colors=["#4CAF50", "#FF9800"], autopct="%1.0f%%",
           startangle=90, textprops={"fontsize": 9})
    ax.set_title("FC–Attribution Concordance\n(all genes)")

    plt.savefig("plant_stress_cv_results.png", dpi=150, bbox_inches="tight")
    print("✔  Saved → plant_stress_cv_results.png")


# ─────────────────────────────────────────────────────────────
#  5.  PREDICT NEW SAMPLES
# ─────────────────────────────────────────────────────────────
def predict_new_samples(new_csv: str, model, scaler, pca, gene_cols):
    """
    Predict treatment class for new plant samples.
    new_csv must have a 'sampleID' column + the same gene columns.
    Missing genes are filled with 0.
    """
    df_new = pd.read_csv(new_csv)
    ids    = (df_new["sampleID"].values if "sampleID" in df_new.columns
              else np.arange(len(df_new)))
    X_new  = df_new.reindex(columns=gene_cols, fill_value=0.0).values.astype(np.float32)
    X_s    = scaler.transform(X_new)
    X_p    = pca.transform(X_s).astype(np.float32)
    probs  = predict_proba(model, X_p)
    preds  = (probs >= 0.5).astype(int)
    labels = ["Stressed" if p else "CK (Control)" for p in preds]

    out = pd.DataFrame({
        "sampleID"   : ids,
        "prediction" : labels,
        "P(Stressed)": probs.round(4),
    })
    print(out.to_string(index=False))
    out.to_csv("new_sample_predictions.csv", index=False)
    print("✔  Saved → new_sample_predictions.csv")
    return out


# ─────────────────────────────────────────────────────────────
#  6.  MAIN — K-FOLD CV PIPELINE
# ─────────────────────────────────────────────────────────────
def main():
    CSV_PATH = "ML_Matrix.csv"        # ← change if needed
    X, y, ids, gene_cols = load_data(CSV_PATH)
    n_genes = len(gene_cols)

    skf          = StratifiedKFold(n_splits=N_FOLDS, shuffle=True, random_state=SEED)
    fold_metrics = []
    all_probs    = np.zeros(len(y))
    all_ig_gene  = np.zeros((len(y), n_genes))
    all_X_gene   = np.zeros((len(y), n_genes))

    print("=" * 58)
    print(f"  K-Fold Stratified Cross-Validation  [{BACKEND}]")
    print("=" * 58)

    for fold, (tr_i, va_i) in enumerate(skf.split(X, y), 1):
        print(f"\n── Fold {fold}/{N_FOLDS} ──────────────────────────────────")

        X_tr, X_va = X[tr_i], X[va_i]
        y_tr, y_va = y[tr_i], y[va_i]

        # Pre-process
        scaler = StandardScaler()
        X_tr_s = scaler.fit_transform(X_tr)
        X_va_s = scaler.transform(X_va)

        pca    = PCA(n_components=N_PCA, random_state=SEED)
        X_tr_p = pca.fit_transform(X_tr_s).astype(np.float32)
        X_va_p = pca.transform(X_va_s).astype(np.float32)
        print(f"   PCA variance: {pca.explained_variance_ratio_.cumsum()[-1]*100:.1f}%")

        # Train
        model = build_model()
        model = train_fold(model, X_tr_p, y_tr, X_va_p, y_va)

        # Evaluate
        probs = predict_proba(model, X_va_p)
        acc   = ((probs >= 0.5).astype(int) == y_va).mean()
        auc   = roc_auc_score(y_va, probs)
        fold_metrics.append({"fold": fold, "acc": acc, "auc": auc})
        all_probs[va_i] = probs
        print(f"   Accuracy : {acc:.4f}  |  AUC : {auc:.4f}")
        print(classification_report(y_va, probs >= 0.5,
              target_names=["CK (0)", "Stressed (1)"], digits=3))

        # Integrated Gradients
        print("   Computing Integrated Gradients …", end="", flush=True)
        ig_pca            = integrated_gradients(model, X_va_p, steps=IG_STEPS)
        all_ig_gene[va_i] = ig_pca @ pca.components_   # back-project → gene space
        all_X_gene[va_i]  = X_va                        # raw expression
        print(" done")

    # ── Aggregate report ─────────────────────────────────────
    mean_acc = np.mean([m["acc"] for m in fold_metrics])
    std_acc  = np.std( [m["acc"] for m in fold_metrics])
    mean_auc = np.mean([m["auc"] for m in fold_metrics])
    std_auc  = np.std( [m["auc"] for m in fold_metrics])
    print("\n" + "=" * 58)
    print("  Aggregate Cross-Validation Performance")
    print("=" * 58)
    print(f"  Accuracy : {mean_acc:.4f} ± {std_acc:.4f}")
    print(f"  ROC-AUC  : {mean_auc:.4f} ± {std_auc:.4f}")
    print(classification_report(y, all_probs >= 0.5,
          target_names=["CK (0)", "Stressed (1)"], digits=3))

    # ── Feature ranking ──────────────────────────────────────
    print("── Building Feature Rank Table …")
    feature_df = build_feature_table(all_ig_gene, all_X_gene, y,
                                     gene_cols, X, y)
    feature_df.to_csv("gene_feature_ranks.csv", index=False)
    print(f"✔  Saved → gene_feature_ranks.csv")

    pd.set_option("display.float_format", "{:.8f}".format)
    pd.set_option("display.width", 120)
    print(f"\nTop 20 Ranked Genes:\n")
    print(feature_df.head(20).to_string(index=False))

    # ── Retrain final model on all data ─────────────────────
    print("\n── Retraining final model on full dataset …")
    scaler_f = StandardScaler()
    pca_f    = PCA(n_components=N_PCA, random_state=SEED)
    X_s_f    = scaler_f.fit_transform(X)
    X_p_f    = pca_f.fit_transform(X_s_f).astype(np.float32)
    tr_i, va_i = train_test_split(np.arange(len(y)), test_size=0.15,
                                   stratify=y, random_state=SEED)
    model_f = build_model()
    model_f = train_fold(model_f, X_p_f[tr_i], y[tr_i], X_p_f[va_i], y[va_i])

    if BACKEND == "pytorch":
        torch.save({"model_state": model_f.state_dict(),
                    "scaler": scaler_f, "pca": pca_f,
                    "gene_cols": gene_cols, "n_pca": N_PCA},
                   "plant_stress_model.pt")
        print("✔  Saved → plant_stress_model.pt")
    else:
        import pickle
        with open("plant_stress_model.pkl", "wb") as f:
            pickle.dump({"model": model_f, "scaler": scaler_f,
                         "pca": pca_f, "gene_cols": gene_cols}, f)
        print("✔  Saved → plant_stress_model.pkl")

    # ── Plot ─────────────────────────────────────────────────
    make_plots(fold_metrics, y, all_probs, feature_df)

    # ── Predict new samples (example) ───────────────────────
    # predict_new_samples("new_plants.csv", model_f, scaler_f, pca_f, gene_cols)

    print("\n✅  Pipeline complete.")
    return feature_df


if __name__ == "__main__":
    feature_df = main()
