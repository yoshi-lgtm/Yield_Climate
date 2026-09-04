# Yield_Climate

関東地方の水稲収量と気象条件の関係を、市町村×年のパネルデータで分析する。

## 環境変数

**動かす前に、自分の環境に合わせて設定すること。**

| 変数 | 用途 | 既定値 |
| --- | --- | --- |
| `YIELD_DATA_DIR` | データの置き場所（この下に `yield/` `results/` などができる） | `.`（カレントディレクトリ） |
| `ESTAT_APPID` | e-Stat API の appId。[e-Stat](https://www.e-stat.go.jp/api/) で取得する | なし（未設定だとエラー） |

```bash
# bash
export YIELD_DATA_DIR=/path/to/data
export ESTAT_APPID=xxxxxxxx
```

```powershell
# PowerShell
$env:YIELD_DATA_DIR = "D:\data"
$env:ESTAT_APPID    = "xxxxxxxx"
```

**API キーはコードに直接書かない。** 環境変数か `.env`（gitignore 済み）で渡す。

## データ

大容量の生データ（`*.csv` `*.gpkg` `*.nc` など）と生成物は git 管理外。
気象データは [メッシュ農業気象データ](https://amu.rd.naro.go.jp/) から `AMD_Tools4.py` で取得する
（このモジュール自体も配布条件のため管理外）。

## 実行順序

`Data_cleaning/` の番号順に実行する。

| スクリプト | 入力 | 出力 |
| --- | --- | --- |
| `01read_and_write.py` | メッシュ農業気象データ（API） | `MetData_Output/` |
| `02merge_land_use_data_year.py` | `LandUseData/` | `merged_data.gpkg` |
| `03calculate_all_def_debuged.py` | `ShapeFile/`, `merged_data.gpkg` | `avg_temp/`, `paddy/`, `total_puddyratio/` |
| `04yield_data_kanto_debuged.py` | e-Stat API | `yield/kanto_yield.csv` |
| `05make_df_kanto_debuged.py` | 上記すべて | `results/kanto.csv` |
| `analysis.py` | `results/kanto.csv` | 図（PNG） |

## 分析

- `research/model/model.R` — fixest による固定効果モデル（TWFE）
- `research/model/fe_model.ipynb` — Python 側での検討
