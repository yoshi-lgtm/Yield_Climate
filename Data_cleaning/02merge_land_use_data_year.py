
"""
1kmメッシュの土地利用データを結合し、GPKG形式で保存するスクリプト
- 入力: 複数のZIPファイル（各ZIP内にシェープファイルが含まれる）
- 出力: 結合された土地利用データをGPKG形式で保存
- 座標系変換: JGD1941 (EPSG:4301) から JGD2000 (EPSG:4612)
"""
# %%
# Import necessary libraries
import pandas as pd
import geopandas as gpd
import zipfile
import os
import glob
import shutil
import re
import tempfile
import time

# %%
# Define paths
input_dir = "../LandUseData/"
output_file = "merged_data.gpkg"
# OneDrive配下だと同期プロセスがファイルをロックし rmtree が PermissionError になることがあるため、
# システムの一時フォルダ（OneDrive管理外）を使う
temp_dir = os.path.join(tempfile.gettempdir(), "landuse_extract")


def safe_rmtree(path, retries=5, delay=0.5):
    """OneDrive等による一時的なファイルロックを考慮し、リトライ付きで削除する"""
    for attempt in range(retries):
        try:
            shutil.rmtree(path)
            return
        except PermissionError:
            if attempt == retries - 1:
                raise
            time.sleep(delay)

# %%
# Get list of zip files
zip_files = glob.glob(os.path.join(input_dir,"*.zip"))


# %%
# Process each zip file
gdf_list = []
for zip_path in zip_files:
  file_base = os.path.basename(zip_path)
  print(f"処理中: {file_base}...")

  # ファイル名から数字（下2桁）を抽出
  match = re.search(r'-(\d{2})_', file_base) # 末尾が「数字2桁.zip」の場合
  if match:
      short_year = int(match.group(1))
        # 下2桁から西暦4桁を復元 (70以上なら1900年代、それ以外は2000年代と判定)
      full_year = 1900 + short_year if short_year >= 70 else 2000 + short_year
  else:
      full_year = 0 # 判別できない場合

  if os.path.exists(temp_dir):
    safe_rmtree(temp_dir)
  
  with zipfile.ZipFile(zip_path, "r") as zip_ref:
    zip_ref.extractall(temp_dir)
    
  shp_files = glob.glob(os.path.join(temp_dir, "**", "*.shp"), recursive=True)
  target_shp = shp_files[0]

  try:
    gdf = gpd.read_file(target_shp, encoding='cp932')
  except:          
    gdf = gpd.read_file(target_shp)  

  gdf = gdf.set_crs(epsg=4301, allow_override=True)  # Set to JGD1941（旧日本測地系）

  # 3. 読み込んだデータに年度列を追加
  gdf['year'] = full_year

  gdf_list.append(gdf)
  print(f"年度 {full_year} として処理完了")
  gdf.head

if os.path.exists(temp_dir):
    shutil.rmtree(temp_dir)

# %%
print("全データの結合中...")
merged_gdf = gpd.GeoDataFrame(pd.concat(gdf_list, ignore_index=True),
                              crs=gdf_list[0].crs)
print("結合完了")

# %%
merged_gdf.to_crs(epsg=4612, inplace=True)  # Convert to JGD2000

# %%
# Save the merged GeoDataFrame
print(f"ファイルに書き出し中: {output_file}")
if output_file.endswith(".gpkg"):
    merged_gdf.to_file(output_file, driver="GPKG")
else:
    merged_gdf.to_file(output_file, encoding='cp932')
print("書き出し完了")

# %%
# カラム一覧を表示
print(gdf.columns.tolist())

# %%
