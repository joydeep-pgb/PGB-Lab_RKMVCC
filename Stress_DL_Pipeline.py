"""
=============================================================
  PLANT STRESS CLASSIFIER  —  Deep Learning Pipeline v4
  ─────────────────────────────────────────────────────────
  Install
  ─────────────────────────────────────────────────────────
    pip install pandas numpy scikit-learn matplotlib scipy
    pip install shap           ← optional but recommended
    pip install torch          ← optional (faster + exact grads)

  ─────────────────────────────────────────────────────────
  Outputs
  -------
    gene_feature_ranks.csv       — gene rank table with SHAP
    plant_stress_cv_results.png  — CV + SHAP feature plots
    plant_stress_model.pt        — saved model (PyTorch)
          OR
    plant_stress_model.pkl       — saved model (sklearn)
=============================================================
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch
import warnings
warnings.filterwarnings("ignore")

from sklearn.model_selection import StratifiedKFold, train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import VarianceThreshold
from sklearn.metrics import (roc_auc_score, confusion_matrix, roc_curve,
                             ConfusionMatrixDisplay, classification_report)
from scipy.stats import pearsonr, wilcoxon

# ─────────────────────────────────────────────────────────────
#  CONFIG
# ─────────────────────────────────────────────────────────────
SEED            = 42
N_FOLDS         = 10
VAR_THRESH      = 0.01
N_BG_SAMPLES    = 100
SHAP_N_SAMPLES  = 100
np.random.seed(SEED)

# ── Detect shap library ───────────────────────────────────────
try:
    import shap as shap_lib
    HAS_SHAP_LIB = True
    print("shap library : available  → will use shap.GradientExplainer")
except ImportError:
    HAS_SHAP_LIB = False
    print("shap library : not found  → using manual GradientSHAP (mathematically equivalent)")
    print("  Install for faster runtime: pip install shap")

# ── Detect PyTorch ────────────────────────────────────────────
try:
    import torch
    import torch.nn as nn
    import torch.optim as optim
    from torch.utils.data import DataLoader, TensorDataset
    BACKEND = "pytorch"
    DEVICE  = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    torch.manual_seed(SEED)
    print(f"Backend      : PyTorch  |  Device : {DEVICE}")
except ImportError:
    from sklearn.neural_network import MLPClassifier
    import pickle
    BACKEND = "sklearn"
    print("Backend      : scikit-learn MLPClassifier")
print()


# ─────────────────────────────────────────────────────────────
#  1.  LOAD DATA
# ─────────────────────────────────────────────────────────────
def load_data(csv_path: str):
    df = pd.read_csv(csv_path)

    # [FIX 4] Data-driven VarianceThreshold — no hardcoded gene drops
    meta_cols      = ["sampleID", "treatment"]
    candidate_cols = [c for c in df.columns if c not in meta_cols]

    vt = VarianceThreshold(threshold=VAR_THRESH)
    vt.fit(df[candidate_cols].values.astype(np.float32))
    gene_cols = [c for c, keep in zip(candidate_cols, vt.get_support()) if keep]
    removed   = len(candidate_cols) - len(gene_cols)
    print(f"VarianceThreshold(var < {VAR_THRESH}): removed {removed} gene(s), "
          f"kept {len(gene_cols)}")

    X   = df[gene_cols].values.astype(np.float32)
    y   = df["treatment"].values.astype(int)
    ids = df["sampleID"].values

    # [FIX 5] Confirm log2-scale: arithmetic mean diff = log2FC
    val_max = float(X.max())
    is_log  = val_max < 30
    print(f"Value range  : 0 – {val_max:.2f}  →  "
          f"{'log2-transformed ✓' if is_log else 'WARNING: possibly raw counts'}")
    print(f"Samples      : {X.shape[0]}  |  Genes : {X.shape[1]}")
    print(f"CK(0) = {(y==0).sum()}  |  Stressed(1) = {(y==1).sum()}\n")
    return X, y, ids, gene_cols


# ─────────────────────────────────────────────────────────────
#  2.  MODEL  (no PCA — operates directly in gene space)
# ─────────────────────────────────────────────────────────────
if BACKEND == "pytorch":

    class ResidualBlock(nn.Module):
        def __init__(self, dim, dropout=0.45):
            super().__init__()
            self.block = nn.Sequential(
                nn.Linear(dim, dim), nn.BatchNorm1d(dim),
                nn.ReLU(), nn.Dropout(dropout),
                nn.Linear(dim, dim), nn.BatchNorm1d(dim))
            self.relu = nn.ReLU()
        def forward(self, x):
            return self.relu(self.block(x) + x)

    class ResNet(nn.Module):
        """
        Residual MLP operating directly in gene space.
        [NEW-3] Input dim = n_genes (397), not N_PCA (50).
        Architecture: 397→64 stem, 2 residual blocks × 64, head 64→1.
        ~42k parameters — appropriate for n=160 samples.
        Dropout raised to 0.45 to compensate for higher input
        dimensionality without PCA compression.
        """
        def __init__(self, in_dim, hidden=64, n_blocks=2, dropout=0.45):
            super().__init__()
            self.stem   = nn.Sequential(
                nn.Linear(in_dim, hidden),
                nn.BatchNorm1d(hidden),
                nn.ReLU())
            self.blocks = nn.Sequential(
                *[ResidualBlock(hidden, dropout) for _ in range(n_blocks)])
            self.head   = nn.Linear(hidden, 1)

        def forward(self, x):
            return self.head(self.blocks(self.stem(x))).squeeze(1)

    def build_model(n_genes):
        # [NEW-3] in_dim = n_genes directly (no PCA)
        return ResNet(in_dim=n_genes, hidden=64, n_blocks=2, dropout=0.45)

    def train_fold(model, X_tr, y_tr, X_va, y_va,
                   epochs=200, lr=1e-3, patience=25, batch=16):
        model = model.to(DEVICE)
        crit  = nn.BCEWithLogitsLoss()
        opt   = optim.Adam(model.parameters(), lr=lr, weight_decay=1e-4)
        sch   = optim.lr_scheduler.ReduceLROnPlateau(opt, patience=10, factor=0.5)

        def mk_loader(X, y, shuffle=True):
            ds = TensorDataset(torch.tensor(X, dtype=torch.float32),
                               torch.tensor(y, dtype=torch.float32))
            return DataLoader(ds, batch_size=batch, shuffle=shuffle,
                              drop_last=(len(X) % batch == 1))

        tr_ld = mk_loader(X_tr, y_tr)
        va_ld = mk_loader(X_va, y_va, shuffle=False)

        best_val, best_state, wait = 1e9, None, 0
        for _ in range(1, epochs + 1):
            model.train()
            for xb, yb in tr_ld:
                xb, yb = xb.to(DEVICE), yb.to(DEVICE)
                opt.zero_grad()
                crit(model(xb), yb).backward()
                opt.step()

            model.eval()
            v_loss = 0.0
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

    def predict_proba(model, X_s):
        model.eval()
        with torch.no_grad():
            logits = model(
                torch.tensor(X_s, dtype=torch.float32).to(DEVICE)
            ).cpu().numpy()
        return 1.0 / (1.0 + np.exp(-logits))

    # ── GradientSHAP (PyTorch) ────────────────────────────────
    def _gradient_shap_torch(model, X_va_s, background):
        """
        [NEW-2] Manual GradientSHAP for PyTorch (used when shap
        library is not installed).

        Algorithm (Smilkov 2017 / Lundberg & Lee 2017):
          For i = 1 … N_SHAP_SAMPLES:
            b      ~ Uniform draw from background set
            alpha  ~ Uniform(0, 1)
            interp = b + alpha * (x − b)
            shap  += autograd_grad(f(interp)) * (x − b)
          shap /= N_SHAP_SAMPLES

        Satisfies the SHAP completeness axiom:
          Σ_j shap_j(x) ≈ f(x) − E[f(background)]

        This is mathematically equivalent to shap.GradientExplainer.
        """
        model.eval()
        X_t  = torch.tensor(X_va_s,  dtype=torch.float32)
        bg_t = torch.tensor(background, dtype=torch.float32)
        n, n_bg = len(X_t), len(bg_t)
        shap_vals = torch.zeros(n, X_t.shape[1])

        for _ in range(SHAP_N_SAMPLES):
            bg_idx = torch.randint(0, n_bg, (n,))
            b      = bg_t[bg_idx]                         # (n, d)
            alpha  = torch.rand(n, 1)                     # (n, 1)
            interp = (b + alpha * (X_t - b)).to(DEVICE)
            interp.requires_grad_(True)
            out    = model(interp)
            grads  = torch.autograd.grad(out.sum(), interp)[0].detach().cpu()
            shap_vals += grads * (X_t - b)

        return (shap_vals / SHAP_N_SAMPLES).numpy()

    def _gradient_shap_lib(model, X_va_s, background):
        """
        [NEW-2] GradientSHAP via shap.GradientExplainer (preferred
        when the shap library is available — same algorithm, optimised
        implementation).

        Root cause of IndexError: shap.GradientExplainer internally does
        `outputs[:, idx]`, which requires a 2D tensor. Our ResNet head ends
        with .squeeze(1), producing shape (n,) — a 1D tensor. Fix: wrap the
        model so it returns shape (n, 1) instead. The wrapper is used only
        for the SHAP explainer; predict_proba continues to use the original.
        """
        class _ModelWrapper(nn.Module):
            """Thin wrapper: unsqueezes output to (n, 1) for shap."""
            def __init__(self, base): super().__init__(); self.base = base
            def forward(self, x):    return self.base(x).unsqueeze(1)

        bg_t     = torch.tensor(background, dtype=torch.float32).to(DEVICE)
        X_t      = torch.tensor(X_va_s,    dtype=torch.float32).to(DEVICE)
        wrapped  = _ModelWrapper(model).to(DEVICE)
        wrapped.eval()
        exp      = shap_lib.GradientExplainer(wrapped, bg_t)
        sv       = exp.shap_values(X_t)
        # GradientExplainer returns list[array] for each output dim,
        # or a single array depending on shap version — handle both.
        # Shape coming in: (n, n_genes, 1) due to the unsqueeze wrapper.
        # Squeeze the trailing dim → (n, n_genes).
        sv_arr = sv[0] if isinstance(sv, list) else sv
        sv_arr = np.array(sv_arr)
        if sv_arr.ndim == 3 and sv_arr.shape[-1] == 1:
            sv_arr = sv_arr.squeeze(-1)
        return sv_arr

    def compute_shap_values(model, X_va_s, X_tr_s, y_tr):
        """
        Dispatch to shap library or manual implementation.
        Background is a stratified sample of the training set so
        SHAP values represent gene contributions relative to the
        real data distribution.
        """
        background = _stratified_background(X_tr_s, y_tr, N_BG_SAMPLES)
        if HAS_SHAP_LIB:
            return _gradient_shap_lib(model, X_va_s, background)
        else:
            return _gradient_shap_torch(model, X_va_s, background)

else:
    # ── sklearn fallback ──────────────────────────────────────
    def build_model(n_genes):
        return MLPClassifier(
            hidden_layer_sizes=(64, 64),
            max_iter=300, random_state=SEED,
            early_stopping=True, validation_fraction=0.15,
            alpha=1e-4)

    def train_fold(model, X_tr, y_tr, X_va, y_va, **kw):
        model.fit(X_tr, y_tr)
        return model

    def predict_proba(model, X_s):
        return model.predict_proba(X_s)[:, 1]

    def _gradient_shap_fd(model, X_va_s, background):
        """
        [NEW-2] Finite-difference GradientSHAP for sklearn models.
        Uses central differences for gradient estimation.
        """
        prob_fn = lambda Z: model.predict_proba(Z)[:, 1]
        n, d    = X_va_s.shape
        n_bg    = len(background)
        eps     = 1e-4
        shap_v  = np.zeros((n, d))

        for _ in range(SHAP_N_SAMPLES):
            b_idx  = np.random.randint(0, n_bg, n)
            b      = background[b_idx]
            alpha  = np.random.rand(n, 1)
            interp = b + alpha * (X_va_s - b)
            for j in range(d):
                Xp = interp.copy(); Xp[:, j] += eps
                Xm = interp.copy(); Xm[:, j] -= eps
                grad_j = (prob_fn(Xp) - prob_fn(Xm)) / (2.0 * eps)
                shap_v[:, j] += grad_j * (X_va_s[:, j] - b[:, j])

        return shap_v / SHAP_N_SAMPLES

    def compute_shap_values(model, X_va_s, X_tr_s, y_tr):
        if HAS_SHAP_LIB:
            background = _stratified_background(X_tr_s, y_tr, N_BG_SAMPLES)
            exp = shap_lib.KernelExplainer(
                lambda Z: model.predict_proba(Z)[:, 1],
                shap_lib.kmeans(background, 10))
            sv = exp.shap_values(X_va_s, nsamples=200)
            return np.array(sv)
        else:
            background = _stratified_background(X_tr_s, y_tr, N_BG_SAMPLES)
            return _gradient_shap_fd(model, X_va_s, background)


# ─────────────────────────────────────────────────────────────
#  HELPER — stratified background sampling
# ─────────────────────────────────────────────────────────────
def _stratified_background(X_tr_s, y_tr, n_total):
    """
    Draw a class-balanced background: n/2 CK + n/2 Stressed samples
    (or as balanced as possible given training class sizes).
    Using a stratified background ensures SHAP values reflect
    contributions relative to the real data distribution.
    """
    n_total = min(n_total, len(X_tr_s))
    n_each  = n_total // 2
    idx0    = np.where(y_tr == 0)[0]
    idx1    = np.where(y_tr == 1)[0]
    rng     = np.random.default_rng(SEED)
    sel0    = rng.choice(idx0, size=min(n_each, len(idx0)), replace=False)
    sel1    = rng.choice(idx1, size=min(n_each, len(idx1)), replace=False)
    return X_tr_s[np.concatenate([sel0, sel1])].astype(np.float32)


# ─────────────────────────────────────────────────────────────
#  3.  FEATURE RANK TABLE
# ─────────────────────────────────────────────────────────────
def _bh_correction(pvals):
    """Benjamini-Hochberg FDR correction."""
    n     = len(pvals)
    order = np.argsort(pvals)
    rank  = np.empty(n, dtype=int)
    rank[order] = np.arange(1, n + 1)
    qvals = np.minimum(1.0, pvals * n / rank)
    for i in range(n - 2, -1, -1):
        qvals[order[i]] = min(qvals[order[i]], qvals[order[i + 1]])
    return qvals


def build_feature_table(all_shap, all_X_gene, y_all, gene_cols, X_full, y_full):
    """
    Gene rank table with columns:
      gene | mean_shap | shap_expr_corr | association |
      shap_fc_concordant | log2fc | pval | qval | rank

    [NEW-2] SHAP values are now direct (no PCA back-projection).
            Column names correctly use 'shap' (v3 used 'ig').
    [FIX 5] log2fc = mean(stressed) − mean(CK) in log2 space.
    [FIX 8] Wilcoxon signed-rank test + BH q-values.
    [FIX 9] No duplicate X/y argument.
    """
    n = len(gene_cols)

    mean_shap      = np.abs(all_shap).mean(axis=0)
    shap_expr_corr = np.zeros(n)

    for g in range(n):
        sv_g, ex_g = all_shap[:, g], all_X_gene[:, g]
        if sv_g.std() > 1e-9 and ex_g.std() > 1e-9:
            shap_expr_corr[g], _ = pearsonr(sv_g, ex_g)

    association = np.where(shap_expr_corr > 0, "Stressed", "CK")

    # [FIX 5] Arithmetic diff = log2FC because data is log2-transformed
    log2fc            = X_full[y_full == 1].mean(0) - X_full[y_full == 0].mean(0)
    mean_shap_signed  = all_shap.mean(0)
    shap_fc_concordant = (
        np.sign(mean_shap_signed) == np.sign(log2fc)).astype(int)

    # [FIX 8] Per-gene Wilcoxon signed-rank test + BH correction
    pvals = np.ones(n)
    for g in range(n):
        sv_g = all_shap[:, g]
        if sv_g.std() > 1e-9:
            try:
                _, pvals[g] = wilcoxon(sv_g, alternative="two-sided")
            except ValueError:
                pvals[g] = 1.0
    qvals = _bh_correction(pvals)

    df = pd.DataFrame({
        "gene"              : gene_cols,
        "mean_shap"         : mean_shap,
        "shap_expr_corr"    : shap_expr_corr,
        "association"       : association,
        "shap_fc_concordant": shap_fc_concordant,
        "log2fc"            : log2fc,
        "pval"              : pvals,
        "qval"              : qvals,
    }).sort_values("mean_shap", ascending=False).reset_index(drop=True)
    df["rank"] = df.index + 1

    sig = (df["qval"] < 0.05).sum()
    print(f"   Significant genes (q < 0.05): {sig} / {n}")
    return df


# ─────────────────────────────────────────────────────────────
#  4.  PLOTTING
# ─────────────────────────────────────────────────────────────
def make_plots(fold_metrics, y_all, all_probs, feature_df):
    fig = plt.figure(figsize=(20, 13))
    fig.suptitle(
        "Plant Stress Classifier — 10-Fold CV  "
        "(GradientSHAP Feature Ranking — Direct Gene Space)",
        fontsize=13, fontweight="bold")
    gs = gridspec.GridSpec(2, 4, figure=fig, hspace=0.48, wspace=0.40)

    folds = [f"F{m['fold']}" for m in fold_metrics]
    accs  = [m["acc"]        for m in fold_metrics]
    aucs  = [m["auc"]        for m in fold_metrics]

    # Per-fold accuracy
    ax = fig.add_subplot(gs[0, 0])
    bars = ax.bar(folds, accs, color="#4CAF50", edgecolor="k", lw=0.7)
    ax.axhline(np.mean(accs), ls="--", color="gray",
               label=f"Mean={np.mean(accs):.2f}")
    ax.set_ylim(0, 1.14)
    ax.set_title("Per-Fold Accuracy")
    ax.legend(fontsize=8)
    for b, a in zip(bars, accs):
        ax.text(b.get_x() + b.get_width() / 2, b.get_height() + 0.01,
                f"{a:.2f}", ha="center", fontsize=7)

    # Per-fold AUC
    ax = fig.add_subplot(gs[0, 1])
    bars = ax.bar(folds, aucs, color="#2196F3", edgecolor="k", lw=0.7)
    ax.axhline(np.mean(aucs), ls="--", color="gray",
               label=f"Mean={np.mean(aucs):.2f}")
    ax.set_ylim(0, 1.14)
    ax.set_title("Per-Fold ROC-AUC")
    ax.legend(fontsize=8)
    for b, a in zip(bars, aucs):
        ax.text(b.get_x() + b.get_width() / 2, b.get_height() + 0.01,
                f"{a:.2f}", ha="center", fontsize=7)

    # Aggregate confusion matrix
    ax = fig.add_subplot(gs[0, 2])
    cm = confusion_matrix(y_all, all_probs >= 0.5)
    ConfusionMatrixDisplay(cm, display_labels=["CK", "Stressed"]).plot(
        ax=ax, colorbar=False)
    ax.set_title("Aggregate Confusion Matrix")

    # ROC curve
    ax   = fig.add_subplot(gs[0, 3])
    fpr, tpr, _ = roc_curve(y_all, all_probs)
    auc  = roc_auc_score(y_all, all_probs)
    ax.plot(fpr, tpr, color="#9C27B0", lw=2, label=f"AUC = {auc:.3f}")
    ax.plot([0, 1], [0, 1], "k--", alpha=0.4)
    ax.set_title("Aggregate ROC Curve")
    ax.legend()
    ax.set_xlabel("FPR (1-specificity)")
    ax.set_ylabel("TPR (sensitivity)")

    # Top-20 SHAP bar (starred if q<0.05)
    ax    = fig.add_subplot(gs[1, :2])
    top20 = feature_df.head(20)
    colors = ["#FF5722" if a == "Stressed" else "#2196F3"
              for a in top20["association"]]
    labels = [f"{g}{'*' if q < 0.05 else ''}"
              for g, q in zip(top20["gene"][::-1], top20["qval"][::-1])]
    ax.barh(labels, top20["mean_shap"][::-1],
            color=colors[::-1], edgecolor="k", lw=0.4)
    ax.set_xlabel("Mean |SHAP|  (GradientSHAP — gene space)")
    ax.set_title("Top-20 Most Important Genes  (* q < 0.05)")
    ax.legend(handles=[Patch(color="#FF5722", label="→ Stressed"),
                        Patch(color="#2196F3", label="→ CK")], fontsize=8)

    # SHAP–expression correlation scatter (top 50)
    ax    = fig.add_subplot(gs[1, 2])
    top50 = feature_df.head(50)
    sc    = ["#FF5722" if a == "Stressed" else "#2196F3"
             for a in top50["association"]]
    ax.scatter(top50["shap_expr_corr"], top50["mean_shap"],
               c=sc, s=35, edgecolors="k", lw=0.3)
    ax.axvline(0, color="gray", ls="--", lw=0.8)
    ax.set_xlabel("SHAP–Expression Correlation  (r)")
    ax.set_ylabel("Mean |SHAP|")
    ax.set_title("SHAP vs Expression Corr\n(top-50 genes, direct gene space)")
    ax.legend(handles=[Patch(color="#FF5722", label="Stressed"),
                        Patch(color="#2196F3", label="CK")], fontsize=8)

    # log2FC–SHAP concordance pie
    ax = fig.add_subplot(gs[1, 3])
    nc = feature_df["shap_fc_concordant"].sum()
    nd = len(feature_df) - nc
    ax.pie([nc, nd],
           labels=[f"Concordant\n({nc})", f"Discordant\n({nd})"],
           colors=["#4CAF50", "#FF9800"], autopct="%1.0f%%",
           startangle=90, textprops={"fontsize": 9})
    ax.set_title("log2FC–SHAP Concordance\n(all genes)")

    plt.savefig("plant_stress_cv_results.png", dpi=150, bbox_inches="tight")
    print("✔  Saved → plant_stress_cv_results.png")


# ─────────────────────────────────────────────────────────────
#  5.  PREDICT NEW SAMPLES
# ─────────────────────────────────────────────────────────────
def predict_new_samples(new_csv: str, model, scaler, gene_cols):
    """
    Predict treatment class for new samples.
    [FIX 6] Missing genes imputed with training-set mean (scaler.mean_).
    [NEW-1] No PCA transform needed — model takes scaled genes directly.
    """
    df_new = pd.read_csv(new_csv)
    ids    = (df_new["sampleID"].values if "sampleID" in df_new.columns
              else np.arange(len(df_new)))

    X_new    = df_new.reindex(columns=gene_cols).values.astype(np.float32)
    nan_mask = np.isnan(X_new)

    # [FIX 6] Mean imputation (training-set means stored in scaler)
    if nan_mask.any():
        n_imp = nan_mask.sum()
        print(f"   Imputing {n_imp} missing value(s) with training-set gene means.")
        for col in range(X_new.shape[1]):
            rows = nan_mask[:, col]
            if rows.any():
                X_new[rows, col] = scaler.mean_[col]

    X_s    = scaler.transform(X_new)         # [NEW-1] no PCA step
    probs  = predict_proba(model, X_s)
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
#  6.  MAIN — 10-FOLD CV PIPELINE
# ─────────────────────────────────────────────────────────────
def main():
    CSV_PATH = "ML_Matrix.csv"        # ← change if needed
    X, y, ids, gene_cols = load_data(CSV_PATH)
    n_genes = len(gene_cols)

    skf          = StratifiedKFold(n_splits=N_FOLDS, shuffle=True, random_state=SEED)
    fold_metrics = []
    all_probs    = np.zeros(len(y))
    all_shap     = np.zeros((len(y), n_genes))   # [NEW-1] gene space directly
    all_X_gene   = np.zeros((len(y), n_genes))

    print("=" * 62)
    print(f"  10-Fold Stratified Cross-Validation  [{BACKEND}]")
    print(f"  SHAP method : {'shap.GradientExplainer' if HAS_SHAP_LIB else 'manual GradientSHAP'}")
    print(f"  Input space : {n_genes} genes (no PCA)")
    print("=" * 62)

    for fold, (tr_i, va_i) in enumerate(skf.split(X, y), 1):
        print(f"\n── Fold {fold}/{N_FOLDS} ──────────────────────────────────────")

        X_tr, X_va = X[tr_i], X[va_i]
        y_tr, y_va = y[tr_i], y[va_i]

        # [NEW-1] Preprocessing: StandardScaler only (no PCA)
        # [FIX 7] Scaler fit on training data only — no leakage
        scaler   = StandardScaler()
        X_tr_s   = scaler.fit_transform(X_tr).astype(np.float32)
        X_va_s   = scaler.transform(X_va).astype(np.float32)

        # Train
        model = build_model(n_genes)
        model = train_fold(model, X_tr_s, y_tr, X_va_s, y_va)

        # Evaluate
        probs = predict_proba(model, X_va_s)
        acc   = ((probs >= 0.5).astype(int) == y_va).mean()
        auc   = roc_auc_score(y_va, probs)
        fold_metrics.append({"fold": fold, "acc": acc, "auc": auc})
        all_probs[va_i] = probs
        print(f"   Accuracy : {acc:.4f}  |  AUC : {auc:.4f}")
        print(classification_report(y_va, probs >= 0.5,
              target_names=["CK (0)", "Stressed (1)"], digits=3))

        # [NEW-2] GradientSHAP — direct gene space, no back-projection
        print("   Computing GradientSHAP …", end="", flush=True)
        all_shap[va_i]   = compute_shap_values(model, X_va_s, X_tr_s, y_tr)
        all_X_gene[va_i] = X_va          # raw log2 expression
        print(" done")

    # ── Aggregate performance ─────────────────────────────────
    mean_acc = np.mean([m["acc"] for m in fold_metrics])
    std_acc  = np.std( [m["acc"] for m in fold_metrics])
    mean_auc = np.mean([m["auc"] for m in fold_metrics])
    std_auc  = np.std( [m["auc"] for m in fold_metrics])
    print("\n" + "=" * 62)
    print("  Aggregate Cross-Validation Performance")
    print("=" * 62)
    print(f"  Accuracy : {mean_acc:.4f} ± {std_acc:.4f}")
    print(f"  ROC-AUC  : {mean_auc:.4f} ± {std_auc:.4f}")
    print(classification_report(y, all_probs >= 0.5,
          target_names=["CK (0)", "Stressed (1)"], digits=3))

    # ── Feature ranking — [FIX 9] no duplicate X/y ───────────
    print("── Building Feature Rank Table …")
    feature_df = build_feature_table(
        all_shap, all_X_gene, y, gene_cols, X, y)
    feature_df.to_csv("gene_feature_ranks.csv", index=False)
    print("✔  Saved → gene_feature_ranks.csv")

    pd.set_option("display.float_format", "{:.6f}".format)
    pd.set_option("display.width", 140)
    print(f"\nTop 20 Ranked Genes:\n")
    print(feature_df.head(20).to_string(index=False))

    # ── Retrain final model on all data — [FIX 7] ────────────
    # Scaler fit only on training split (not full dataset)
    print("\n── Retraining final model on full dataset …")
    tr_i_f, va_i_f = train_test_split(
        np.arange(len(y)), test_size=0.15, stratify=y, random_state=SEED)

    scaler_f  = StandardScaler()
    X_tr_s_f  = scaler_f.fit_transform(X[tr_i_f]).astype(np.float32)  # [FIX 7]
    X_va_s_f  = scaler_f.transform(X[va_i_f]).astype(np.float32)

    model_f   = build_model(n_genes)
    model_f   = train_fold(model_f, X_tr_s_f, y[tr_i_f], X_va_s_f, y[va_i_f])

    if BACKEND == "pytorch":
        torch.save({"model_state": model_f.state_dict(),
                    "scaler"     : scaler_f,
                    "gene_cols"  : gene_cols},      # [NEW-1] no pca key
                   "plant_stress_model.pt")
        print("✔  Saved → plant_stress_model.pt")
    else:
        import pickle
        with open("plant_stress_model.pkl", "wb") as f:
            pickle.dump({"model"    : model_f,
                         "scaler"   : scaler_f,
                         "gene_cols": gene_cols}, f)
        print("✔  Saved → plant_stress_model.pkl")

    # ── Plot ─────────────────────────────────────────────────
    make_plots(fold_metrics, y, all_probs, feature_df)

    # ── Predict new samples (example — uncomment to use) ─────
    # predict_new_samples("new_plants.csv", model_f, scaler_f, gene_cols)

    print("\n✅  Pipeline complete.")
    return feature_df


if __name__ == "__main__":
    feature_df = main()
