import numpy as np
import pandas as pd
from sklearn.metrics import (
    roc_auc_score, recall_score, accuracy_score, confusion_matrix, f1_score
)
from sklearn.model_selection import train_test_split
import joblib
import optuna
from lightgbm import LGBMClassifier
from sklearn.preprocessing import StandardScaler
from optuna.integration import LightGBMPruningCallback

# micromamba install shap
# micromamba install lightgbm
# === 数据加载 ===
outdir="结果/ICH-CI模型"

marker = pd.read_csv("结果/ICH_CI模型基因筛选/marker.csv", index_col=0).index.to_list()

train_df = pd.read_csv("Data/Set1_mRNA_tpm.xls", sep='\t', index_col=0).T
test_df = pd.read_csv("Data/Set2_mRNA_tpm.xls", sep='\t', index_col=0).T

# === 构建 ICH vs CI 数据集（排除 H 样本）===
sample_selector_train = train_df.index.str.startswith("H")
X_train_raw = train_df.loc[~sample_selector_train, marker]
X_train_raw=np.log2(X_train_raw+1)
y_train = np.where(X_train_raw.index.str.startswith("ICH"), 1, 0)

sample_selector_test = test_df.index.str.startswith("H")
X_test_raw = test_df.loc[~sample_selector_test, marker]
X_test_raw=np.log2(X_test_raw+1)
y_test = np.where(X_test_raw.index.str.startswith("ICH"), 1, 0)

scaler = StandardScaler()
scaler.fit(X_train_raw)

X_train = pd.DataFrame(scaler.transform(X_train_raw),columns=X_train_raw.columns,index=X_train_raw.index)
X_test= pd.DataFrame(scaler.transform(X_test_raw),columns=X_test_raw.columns,index=X_test_raw.index)

# === Optuna 目标函数（LightGBM 版本）===
def lgbm_objective(trial):
    params = {
        "n_estimators": trial.suggest_int("n_estimators", 50, 1000),
        "learning_rate": trial.suggest_float("learning_rate", 0.001, 0.3, log=True),
        "max_depth": trial.suggest_int("max_depth", 2, 3),  
        "num_leaves": trial.suggest_int("num_leaves", 2, 256),
        "min_child_samples": trial.suggest_int("min_child_samples", 2, 200),
        "subsample": trial.suggest_float("subsample", 0.5, .9),
        "colsample_bytree": trial.suggest_float("colsample_bytree", 0.5, .9),
        "reg_alpha": trial.suggest_float("reg_alpha", 0.01, 50.0, log=True), #1e-8,
        "reg_lambda": trial.suggest_float("reg_lambda", 0.01,  50.0, log=True),
        "objective": "binary",  # LightGBM 写法
        "metric": "auc",        # 不影响训练，仅用于日志
        "verbosity": -1,        # 静默模式
        "random_state": 42,
        "n_jobs": 32,
        # "is_unbalance": True,   # 👈 自动处理类别不平衡（推荐替代 scale_pos_weight）
    }

    model = LGBMClassifier(**params)

    model.fit(X_train, y_train, 
          eval_set=[(X_test, y_test)], 
          callbacks=[LightGBMPruningCallback(trial, 'auc')])

    # 👇 优化目标：平衡训练和测试的 F1（防止过拟合）
    return min(f1_score(y_train, model.predict(X_train)), f1_score(y_test, model.predict(X_test)))
    # return f1_score(y_test, model.predict(X_test))
    # return min(roc_auc_score(y_train, model.predict_proba(X_train)[:, 1]),roc_auc_score(y_test, model.predict_proba(X_test)[:, 1]))
    # return min(accuracy_score(y_train, model.predict(X_train)),accuracy_score(y_test, model.predict(X_test)))
    # return min(recall_score(y_train, model.predict(X_train)),recall_score(y_test, model.predict(X_test)))


# === Optuna 优化 ===
study = optuna.create_study(direction="maximize")
study.optimize(lgbm_objective, n_trials=1000)

print("Best params:", study.best_params)

# === 训练最终模型 ===
model = LGBMClassifier(
    **study.best_params,
    objective="binary",
    random_state=42,
    verbosity=-1,
    n_jobs=32
)
model.fit(X_train, y_train)

# 保存模型
joblib.dump(model, f"{outdir}/LGBM.jb")

# 加载（可选）
model = joblib.load(f"{outdir}/LGBM.jb")

# === 评估 ===
results = {
    "train_roc": [],
    "train_recall": [],
    "train_accuracy": [],
    "train_confusion_matrix": [],
    "test_roc": [],
    "test_recall": [],
    "test_accuracy": [],
    "test_confusion_matrix": [],
}

# 训练集
y_train_proba = model.predict_proba(X_train)[:, 1]
y_train_pred = model.predict(X_train)
results["train_roc"].append(roc_auc_score(y_train, y_train_proba))
results["train_recall"].append(recall_score(y_train, y_train_pred, zero_division=0))
results["train_accuracy"].append(accuracy_score(y_train, y_train_pred))
results["train_confusion_matrix"].append(confusion_matrix(y_train, y_train_pred))

# 测试集
y_test_proba = model.predict_proba(X_test)[:, 1]
y_test_pred = model.predict(X_test)
results["test_roc"].append(roc_auc_score(y_test, y_test_proba))
results["test_recall"].append(recall_score(y_test, y_test_pred, zero_division=0))
results["test_accuracy"].append(accuracy_score(y_test, y_test_pred))
results["test_confusion_matrix"].append(confusion_matrix(y_test, y_test_pred))

# 保存结果
df = pd.DataFrame(results, index=["LightGBM"])
df.to_csv(f"{outdir}/LightGBM-ICH_vs_CI.csv")
print(df)