# This code is to plot the manhattan map

import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import matplotlib.patches as mpatches
import pandas as pd

df=gpd.read_file('data/taxi_zones/taxi_zones.shp')

if df.crs != "EPSG:4326":
    df = df.to_crs(epsg=4326)

df['centroid'] = df.geometry.centroid
df['longitude'] = df.centroid.x
df['latitude'] = df.centroid.y

df1=df[df.borough=='Manhattan']

dis=pd.read_csv("data/mahh.csv")

df1=df1.assign(dis0=0)
for i in range(len(df1)):
  df1.iloc[i,-1]=dis.iloc[i,-1]

colors = ['#fdae61', '#2c7bb6', '#DDA0DD']  
cmap = ListedColormap(colors)
df1['dis0'] = df1['dis0'].astype('category')
fig, ax = plt.subplots(1, 1, figsize=(10, 8))
df1.plot(column='dis0',
         cmap=cmap,
         ax=ax,
         edgecolor='black',
         legend=False)  

category_labels = ['Downtown', 'Midtown', 'Uptown']  
handles = [mpatches.Patch(color=colors[i], label=category_labels[i]) for i in range(3)]
ax.legend(handles=handles, loc='upper left', bbox_to_anchor=(0.02, 0.85), fontsize=12)
ax.set_axis_off()
plt.savefig('map.png', dpi=300, bbox_inches='tight')
plt.show()