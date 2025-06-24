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

dis=pd.read_csv("Dis.csv")

df1=df1.assign(dis0=0)
for i in range(len(df1)):
  df1.iloc[i,-1]=dis.iloc[i,-1]

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
fig, ax = plt.subplots(1, 1, figsize=(10, 8))
cmap = cm.YlOrRd
df1.plot(column='dis0',vmax=2*1e12,legend=True,ax=ax,cmap=cmap,edgecolor='black'
,legend_kwds={'label': "dis",'shrink': 0.3})
ax.set_axis_off()
plt.savefig('Figure 6.png', dpi=300, bbox_inches='tight')
plt.show()
