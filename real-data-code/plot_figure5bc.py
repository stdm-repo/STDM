# This code is to plot the figure 5b and 5c

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

dis=pd.read_csv("Dis2.csv")

df1=df1.assign(dis0=0)
for i in range(len(df1)):
  df1.iloc[i,-1]=dis.iloc[i,0]
df1=df1.assign(dis1=0)
for i in range(len(df1)):
  df1.iloc[i,-1]=dis.iloc[i,1]

custom_vmin = min(df1['dis0'].min(),df1['dis1'].min())
custom_vmax = max(df1['dis0'].max(),df1['dis1'].max())

fig, ax = plt.subplots(1, 1, figsize=(10, 8))
ax.axis('off')
df1.plot(column='dis0',legend=False,ax=ax,cmap = 'BuPu',edgecolor="Black",vmin=custom_vmin,vmax= custom_vmax)
plt.savefig('real_4.png', dpi=300, bbox_inches='tight')
plt.show()

fig, ax = plt.subplots(1, 1, figsize=(10, 8))
ax.axis('off')
df1.plot(column='dis1',legend=True,ax=ax,cmap = 'BuPu',edgecolor="Black",vmin=custom_vmin,vmax= custom_vmax,
         legend_kwds={'shrink': 0.3})
plt.savefig('real_5.png', dpi=300, bbox_inches='tight')
plt.show()