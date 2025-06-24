# This code is to plot the figure 6

import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors

df=gpd.read_file('data/taxi_zones/taxi_zones.shp')

if df.crs != "EPSG:4326":
    df = df.to_crs(epsg=4326)

df['centroid'] = df.geometry.centroid
df['longitude'] = df.centroid.x
df['latitude'] = df.centroid.y

df1=df[df.borough=='Manhattan']

# read the dis data in full time
dis_our=pd.read_csv("Dis_our.csv")
dis_dynamic=pd.read_csv("Dis_dynamic.csv")
dis_structure=pd.read_csv("Dis_structure.csv")
df1=df1.assign(dis0=0)
df1=df1.assign(dis1=0)
df1=df1.assign(dis2=0)
for i in range(len(df1)):
  df1.iloc[i,df1.columns.get_loc('dis0')]=dis_our.iloc[i,-1]
  df1.iloc[i,df1.columns.get_loc('dis1')]=dis_dynamic.iloc[i,-1]
  df1.iloc[i,df1.columns.get_loc('dis2')]=dis_structure.iloc[i,-1]


fig, axes = plt.subplots(1, 3, figsize=(18, 6))
cmap = cm.YlOrRd
names = ['Our', 'No dynamic', 'No structure']
for idx, col in enumerate(['dis0', 'dis1', 'dis2']):
    vmax = df1[col].quantile(0.95)
    df1.plot(
        column=col,
        vmax=vmax,
        legend=True,
        ax=axes[idx],
        cmap=cmap,
        edgecolor='black',
        legend_kwds={'shrink': 0.3}
    )
    axes[idx].set_title(names[idx], fontsize=16, fontweight='bold')
    axes[idx].set_axis_off()
plt.tight_layout()
plt.savefig('appendix_full_time.png', dpi=300, bbox_inches='tight')
plt.show()
