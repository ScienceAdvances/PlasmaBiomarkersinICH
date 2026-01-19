import numpy as np
import matplotlib.pyplot as plt
import shap
from sklearn.datasets import load_breast_cancer
from sklearn.model_selection import train_test_split
from xgboost import XGBClassifier
from statsmodels.nonparametric.smoothers_lowess import lowess
from scipy.stats import scoreatpercentile
from scipy.interpolate import interp1d
import matplotlib
import atopos
from pathlib import Path
import shap
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
from scipy.interpolate import interp1d
from sklearn.metrics import r2_score
from matplotlib.colors import Normalize
import matplotlib.cm as cm
from matplotlib.lines import Line2D
import pathlib
from atopos.pl import Palette
matplotlib.rcParams['font.sans-serif'] = ['Arial', 'Microsoft YaHei'] # 指定默认字体
matplotlib.rcParams['axes.unicode_minus'] = False # 解决保存图像是负号’-'显示为方块的问题

# 设置字体渲染方式（关键！）
matplotlib.rcParams['svg.fonttype'] = 'none'  # SVG保存时保留文本为可编辑字体
matplotlib.rcParams['pdf.fonttype'] = 42      # PDF保存时使用TrueType字体

# ------------------------------------------------------
# 2. 计算 SHAP 值
# ------------------------------------------------------
# 初始化 SHAP 解释器（树模型专用）
explainer = shap.TreeExplainer(model)
# 计算测试集的 SHAP 值（也可使用训练集，根据需求选择）
shap_values = explainer.shap_values(X_test)  # 输出形状: (n_samples, n_features)

# ------------------------------------------------------
# 3. 选择一个特征进行分析（以第0个特征为例，可替换为其他特征）
# ------------------------------------------------------
feature_names = X_test.columns.tolist()  # 特征名称列表
feature_names=marker
for i in feature_names:
    plot_shap_dependence(shap_values, feature_name=i, zone_name="")
    plt.close('all')
explainer_values = explainer(X_test)
outdir="结果/ICH-CI模型"
Path(f"{outdir}/瀑布图/").mkdir(exist_ok=True,parents=True)
for i in range(X_test.shape[0]):
    shap.plots.waterfall(explainer_values[i],max_display=10)
    plt.tight_layout()
    plt.savefig(f"{outdir}/瀑布图/waterfall_{X_test.index[i]}.pdf", dpi=300, bbox_inches='tight')
    plt.close('all')


shap.summary_plot(shap_values, X_test)
plt.savefig(f'{outdir}/shap_summary.pdf', dpi=300)
plt.close('all')
shap.summary_plot(shap_values, X_test, plot_type="bar")
plt.savefig(f'{outdir}/shap_bar.pdf', dpi=300)
plt.close('all')

colors=["#a2d2e7", "#ffc17f", "#cf9f88", "#6fb3a8", "#b3e19b", "#50aa4b",
    "#ff9d9f", "#f36569", "#3581b7", "#cdb6da", "#704ba3", "#9a7fbd", "#dba9a8", "#e43030", "#e99b78", "#ff8831",
    "#e75b58", "#e49ac3", "#ab3181", "#20452e", "#bb9369", "#8c529a", "#e4d1dc", "#52a65e", "#f0b971", "#f2b09e",
    "#d4e6a1", "#55c0f2", "#496d87"]
def ROC(y_true:dict, y_hat:dict,xlabel=None,ylabel=None, outdir=None,filename='',palette=colors):
    """
    plot multiple ROC curves in one figure
    y_true: dict
    y_hat: dict
    outdir: output directory
    filename: filename
    """
    from sklearn.metrics import RocCurveDisplay
    colors = dict(zip(y_hat,palette))
    for i in y_hat.keys():
        RocCurveDisplay.from_predictions(y_true[i], y_hat[i], name = i, color= colors[i], ax=plt.gca(),linewidth=2)
    plt.plot([0,1],[0,1], linestyle="dashed",color = "grey");
    if ylabel:
        plt.ylabel(ylabel)
    else:
        plt.ylabel("True Positive Rate");
    if xlabel:
        plt.xlabel(xlabel);
    else:
        plt.xlabel("False Positive Rate");
    p = plt.gcf()
    p.set_size_inches(6, 6)
    plt.close()
    if outdir:
        p.savefig(pathlib.Path(outdir).joinpath(f"{filename}.pdf"))
    return p


ROC(
    y_true = {'Train':y_train, 'Test': y_test },
    y_hat = {'Train':y_train_proba, 'Test': y_test_proba },
    xlabel='1-Specificity',
    ylabel='Sensitivity',
    outdir=outdir,
    filename='ICH-CI',
)
    y_true: Dict[str, np.ndarray],
    y_hat: Dict[str, np.ndarray],
    xlabel: str,
    ylabel: str,
    outdir: Union[pathlib.Path, str] = pathlib.Path.cwd(),
    filename: str = "CalibrationCurve",
    palette: List[str] = Palette.set2

calibration_curve(
    y_true = {'Train':y_train, 'Test': y_test },
    y_hat = {'Train':y_train_proba, 'Test': y_test_proba },
    xlabel='Predicted probability',
    ylabel='Observed probability',
    outdir=outdir,
    palette = colors,
    filename='calibration_curve',
)

precision_recall_curve(
    y_true = {'Train':y_train, 'Test': y_test },
    y_hat = {'Train':y_train_proba, 'Test': y_test_proba },
    xlabel='Recall',
    ylabel='Precision',
    outdir=outdir,
    palette = colors,
    filename='PrecisionRecallCurve',
)


ROC(
y_hat= {'Fold':y_train_proba,  'ADAM17':X_train['ADAM17'].values,'TMPRSS5':X_train['TMPRSS5'].values,'PLAU':X_train['PLAU'].values,'ADAMTS13':X_train['ADAMTS13'].values},
y_true = {'Fold':y_train,  'ADAM17':y_train,'TMPRSS5':np.where(y_train==0,1,0),'PLAU':y_train,'ADAMTS13':np.where(y_train==0,1,0) },
xlabel='1-Specificity',
ylabel='Sensitivity',
outdir=outdir,
filename='Train',
)


plt.close('all')
ROC(
y_hat= {'Fold':y_test_proba,  'ADAM17':X_test['ADAM17'].values,'TMPRSS5':X_test['TMPRSS5'].values,'PLAU':X_test['PLAU'].values,'ADAMTS13':X_test['ADAMTS13'].values},
y_true = {'Fold':y_test,  'ADAM17':y_test,'TMPRSS5':np.where(y_test==0,1,0),'PLAU':y_test,'ADAMTS13':np.where(y_test==0,1,0) },
xlabel='1-Specificity',
ylabel='Sensitivity',
outdir=outdir,
filename='Test',
)

