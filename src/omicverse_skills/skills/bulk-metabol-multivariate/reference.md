# Multivariate & biomarker — quick commands

## PLS-DA / OPLS-DA on Pareto-scaled data

```python
import omicverse as ov
import matplotlib.pyplot as plt
ov.plot_set()

# Preprocess: PQN → log → Pareto (see preprocessing skill for full chain)
adata = ov.metabol.read_metaboanalyst('human_cachexia.csv', group_col='Muscle loss')
adata = ov.metabol.normalize(adata, method='pqn')
adata = ov.metabol.transform(adata, method='log')
adata = ov.metabol.transform(adata, method='pareto', stash_raw=False)

# 1) PLS-DA baseline
pls = ov.metabol.plsda(adata, group_col='group', n_components=2, scale=False)
print(f'PLS-DA  R²X={pls.r2x:.3f}  R²Y={pls.r2y:.3f}  Q²={pls.q2:.3f}')

# 2) OPLS-DA with one orthogonal component
opls = ov.metabol.opls_da(adata, group_col='group', n_ortho=1, scale=False)
print(f'OPLS-DA R²X={opls.r2x:.3f}  R²Y={opls.r2y:.3f}  Q²={opls.q2:.3f}')
```

## VIP table + score plot + S-plot

```python
# VIP table — sorted descending; columns include 'vip' and 'coef'
vip_df = opls.to_vip_table(adata.var_names)
print(vip_df.head(15))

# VIP bar — colour encodes sign of regression coefficient
fig, ax = ov.metabol.vip_bar(opls, adata.var_names, top_n=15)
ax.set_title(f'OPLS-DA top-15 VIP\n(red = ↑ {opls.group_labels[0]}, blue = ↑ {opls.group_labels[1]})')
plt.tight_layout(); plt.show()

# Score plot — predictive vs orthogonal
fig, ax = plt.subplots(figsize=(4.8, 4))
for grp, marker in zip(opls.group_labels, 'os'):
    mask = (adata.obs['group'] == grp).values
    ax.scatter(opls.scores[mask, 0], opls.x_ortho_scores[mask, 0],
               marker=marker, label=grp, s=36, alpha=0.85,
               edgecolor='black', linewidth=0.5)
ax.axvline(0, c='grey', lw=0.7)
ax.set(xlabel='t[1] — predictive', ylabel='t[o1] — orthogonal',
       title=f'OPLS-DA  R²Y={opls.r2y:.3f}  Q²={opls.q2:.3f}')
ax.legend(frameon=False); plt.tight_layout(); plt.show()

# S-plot — covariance vs correlation against the predictive component
fig, ax = ov.metabol.s_plot(opls, adata, label_top_n=10)
ax.set_title('OPLS-DA S-plot — covariance (x) vs correlation (y)')
plt.tight_layout(); plt.show()
```

## Per-feature AUC with bootstrap CI

```python
# IMPORTANT: roc_feature consumes log-transformed (NOT Pareto-scaled) data.
adata_log = ov.metabol.read_metaboanalyst('human_cachexia.csv', group_col='Muscle loss')
adata_log = ov.metabol.impute(adata_log, method='qrilc', seed=0)
adata_log = ov.metabol.normalize(adata_log, method='pqn')
adata_log = ov.metabol.transform(adata_log, method='log')

auc = ov.metabol.roc_feature(
    adata_log, group_col='group',
    pos_group='cachexic', neg_group='control',
    ci=True, n_bootstrap=500, seed=0,
)
print(auc.head(15))
# Columns: auc, auc_ci_lower, auc_ci_upper, pvalue (vs auc=0.5)
```

## Multi-metabolite biomarker panel (nested CV + permutation null)

```python
panel = ov.metabol.biomarker_panel(
    adata_log,
    group_col='group',
    pos_group='cachexic', neg_group='control',
    features=10,                # int picks top-N inside inner CV; list[str] fixes the panel
    classifier='lr',            # 'lr' (small cohorts), 'rf' (robust), 'svm' (non-linear)
    cv_outer=5, cv_inner=3,
    n_permutations=100,         # use >=100 for any reportable panel; 0 disables
    seed=0,
)

print(f'mean AUC (outer folds): {panel.mean_auc:.3f} ± {panel.std_auc:.3f}')
print(f'permutation p-value   : {panel.permutation_pvalue:.3f}')
print(panel.feature_importance.head())
print(panel.to_frame())          # per-fold AUC table
```

## Held-out ROC curve

```python
from sklearn.metrics import roc_curve, auc as _auc

fpr, tpr, _ = roc_curve(panel.outer_labels, panel.outer_predictions)
fig, ax = plt.subplots(figsize=(5, 5))
ax.plot(fpr, tpr, lw=2, label=f'nested CV AUC = {_auc(fpr, tpr):.3f}')
ax.plot([0, 1], [0, 1], '--', color='gray')
ax.set(xlabel='False positive rate', ylabel='True positive rate', aspect='equal')
ax.legend(loc='lower right')
plt.tight_layout(); plt.show()
```
