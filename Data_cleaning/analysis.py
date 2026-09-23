#%%
# ライブラリのインポート
import os
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
import matplotlib.pyplot as plt
import seaborn as sns
from linearmodels import PanelOLS


# %%
# --- 設定 ---
# データの置き場所。環境変数 YIELD_DATA_DIR で上書きできる
DATA_DIR = Path(os.environ.get("YIELD_DATA_DIR", "."))

# %%
# ファイル読み込みと基本統計量
df = pd.read_csv(DATA_DIR / "model" / "kanto_consol.csv")
df.describe()

#%%
# 欠損値と型
df.info()

#%%
# 回帰に使う説明変数
features = ['heat', 'GSR', 'APCP', 'TMP_mea', 'paddy_ratio']

#%%
# 相関を見る
df[features + ['yield']].corr()

#%%
# VIFを見る（説明変数のみで計算。目的変数(yield)や識別子(year, city_id)は対象外）
corr_features = df[features].corr()
vif = np.diag(np.linalg.inv(corr_features.values))
df_vif = pd.DataFrame({'VIF': vif}, index=corr_features.columns)
df_vif.sort_values('VIF')

#%%
# yieldおよび説明変数に欠損値がある行を削除
# （GSR/TMP_meaに日射量データ欠測由来のNaNが一部の市町村・年で発生しているため）
df_analysis = df.dropna(subset=['yield'] + features).copy()
df_analysis.info()

#%%
# 係数プロット用の関数
def plot_coefficients(model, title, filename):
    coef_df = model.params.drop('const').reset_index()
    coef_df.columns = ['Feature', 'Coefficient']
    plt.figure(figsize=(10, 6))
    sns.barplot(x='Coefficient', y='Feature', data=coef_df.sort_values('Coefficient'))
    plt.title(title)
    plt.axvline(0, color='black', linestyle='--')
    plt.savefig(filename)
    plt.show()


year_range = f"{df_analysis['year'].min()}-{df_analysis['year'].max()}"

# %%
# 重回帰分析（ベースモデル）
X = sm.add_constant(df_analysis[features])
y = df_analysis['yield']

model = sm.OLS(y, X).fit()
print(model.summary())
plot_coefficients(model, f'Impact of Environmental Factors on Rice Yield ({year_range})', 'coefficients_base.png')

#%%
# 高温と収量の関係
sns.regplot(x='heat', y='yield', data=df_analysis, scatter_kws={'alpha': 0.5}, line_kws={'color': 'red'})
plt.title('Relationship between Heat Days and Rice Yield')
plt.xlabel('Number of days >= 34C (July-Aug)')
plt.ylabel('Yield')
plt.grid(True)
plt.savefig('heat_yield_regression.png')
plt.show()

# %%
# 交差項(heat × paddy_ratio)を追加したモデル
df_analysis['heat_paddy'] = df_analysis['heat'] * df_analysis['paddy_ratio']
features_interaction = features + ['heat_paddy']

X = sm.add_constant(df_analysis[features_interaction])
y = df_analysis['yield']

model_interaction = sm.OLS(y, X).fit()
print(model_interaction.summary())
plot_coefficients(
    model_interaction,
    f'Impact on Rice Yield with heat x paddy_ratio interaction ({year_range})',
    'coefficients_interaction.png'
)

# %%
# 固定効果モデル（市町村固定効果・年固定効果）
df_panel = df_analysis.set_index(['city_id', 'year'])

X = sm.add_constant(df_panel[features_interaction])
y = df_panel['yield']

model_fe = PanelOLS(y, X, entity_effects=True, time_effects=True)
res_fe = model_fe.fit()
print(res_fe.summary)
plot_coefficients(
    res_fe,
    f'Fixed Effects: Impact on Rice Yield ({year_range})',
    'coefficients_fe.png'
)

# %%
# TMP_meaの2乗項を追加したモデル（気温と収量の非線形性を確認）
df_panel['tmp_mea_sq'] = df_panel['TMP_mea'] ** 2
features_quad = features_interaction + ['tmp_mea_sq']

X = sm.add_constant(df_panel[features_quad])
y = df_panel['yield']

model_quad = sm.OLS(y, X).fit()
print(model_quad.summary())

model_fe_quad = PanelOLS(y, X, entity_effects=True, time_effects=True)
res_fe_quad = model_fe_quad.fit()
print(res_fe_quad.summary)
plot_coefficients(
    res_fe_quad,
    f'Fixed Effects with TMP_mea^2: Impact on Rice Yield ({year_range})',
    'coefficients_fe_quad.png'
)

# %%
# paddy_ratioの分布
plt.figure(figsize=(10, 6))
sns.histplot(df_analysis['paddy_ratio'], bins=30, kde=True)
plt.title(f'Distribution of paddy_ratio ({year_range})')
plt.xlabel('paddy_ratio')
plt.savefig('paddy_ratio_distribution.png')
plt.show()

#%%
# paddy_ratioの時系列推移（都県別の年平均）
paddy_by_year_pref = (
    df_analysis.groupby(['year', 'prefecture'])['paddy_ratio']
    .mean()
    .reset_index()
)
plt.figure(figsize=(10, 6))
sns.lineplot(x='year', y='paddy_ratio', hue='prefecture', data=paddy_by_year_pref, marker='o')
plt.title(f'paddy_ratio over time by prefecture ({year_range})')
plt.xlabel('Year')
plt.ylabel('paddy_ratio (mean)')
plt.grid(True)
plt.savefig('paddy_ratio_trend.png')
plt.show()

#%%
# paddy_ratioと収量の関係
sns.regplot(x='paddy_ratio', y='yield', data=df_analysis, scatter_kws={'alpha': 0.5}, line_kws={'color': 'red'})
plt.title('Relationship between paddy_ratio and Rice Yield')
plt.xlabel('paddy_ratio')
plt.ylabel('Yield')
plt.grid(True)
plt.savefig('paddy_ratio_yield_regression.png')
plt.show()

#%%
# paddy_ratioとheatの関係（交差項の符号が逆転する要因の確認用）
sns.regplot(x='paddy_ratio', y='heat', data=df_analysis, scatter_kws={'alpha': 0.5}, line_kws={'color': 'red'})
plt.title('Relationship between paddy_ratio and Heat Days')
plt.xlabel('paddy_ratio')
plt.ylabel('Number of days >= 34C (July-Aug)')
plt.grid(True)
plt.savefig('paddy_ratio_heat_regression.png')
plt.show()

"""
interactive fixed effect
group structure panel
functional analysis

"""
#%%
# レポート作成
import sweetviz as sv

report = sv.analyze(df_panel)
report.show_html('report.html')
# %%
