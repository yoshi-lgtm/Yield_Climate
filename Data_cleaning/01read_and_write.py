"""
このプログラムでは、メッシュ農業気象データのサイトからAMD_Tools4モジュールを使い
栃木県域の気象データをダウンロードし、データをnetCDFに変換して、「MetData_Output」
に保存する。
"""


# %%
import os
import numpy as np
import Yield_Climate.AMD_Tools4 as amd
import xarray as xr
import rioxarray
import pandas as pd

# %%

# --- 設定 ---
itsu = ['1980-01-01', '2025-12-31'] # 全期間
doko = [34.8, 37.3, 138.3, 141.0] # 領域（関東7都県の本土をカバー。東京都の伊豆・小笠原諸島は対象外）
#facs = ["TMP_mea", "TMP_max", "TMP_min", "APCP", "GSR"] # 気象要素(APCP: 降水量, GSR: 全天日射量)
facs = [ "APCP", "GSR"]
# 保存先のルートディレクトリ
output_root = "MetData_Output"

# 開始年と終了年を数値で取得
start_year = int(itsu[0].split("-")[0])
end_year   = int(itsu[1].split("-")[0])

# %%
# --- ループ処理 ---
# 1. 気象要素 (facs) のループ
for fac in facs:
    print(f"=== 要素: {fac} の処理を開始します ===")
    
    # 保存用ディレクトリの作成 (例: MetData_Output/TMP_mea/)
    save_dir = os.path.join(output_root, fac)
    os.makedirs(save_dir, exist_ok=True)

    # 2. 年 (Year) のループ
    for year in range(start_year, end_year + 1):
        
        # その年の期間設定 (例: ['1980-01-01', '1980-12-31'])
        current_itsu = [f"{year}-01-01", f"{year}-12-31"]
        
        print(f"  > {year}年 のデータを取得中...")

        try:
            # (A) データ取得
            # amd.GetMetDataX は単一要素・指定期間で呼ぶ想定
            da = amd.GetMetDataX(fac, current_itsu, doko)

            # データが空でないか簡易チェック (必要に応じて)
            if da is None or da.size == 0:
                print(f"    [Warning] {year}年のデータがありません。スキップします。")
                continue

            # (B) GIS情報の設定 (rioxarray)
            # 空間次元の定義 (経度=x, 緯度=y)
            da.rio.set_spatial_dims(x_dim="lon", y_dim="lat", inplace=True) #daに変更を加える
            da = da.rename({"lon": "x", "lat": "y"})
            
            # CRS定義 (JGD2000: EPSG:4612)
            da.rio.write_crs("EPSG:4612", inplace=True) #daに変更を加える

            # (C) 書き出し
            # ファイル名: 1980.nc
            filename = f"data_{year}.nc"
            output_path = os.path.join(save_dir, filename)
            
            #netCDF保存
            da.to_netcdf(output_path)
            print(f"    保存完了: {output_path}")

            # (D) メモリ解放
            # 次のループに行く前にオブジェクトを削除してメモリを空ける
            del da

        except Exception as e:
            print(f"    [Error] {year}年の処理中にエラーが発生しました: {e}")

print("\n全ての処理が完了しました。")



# %%
