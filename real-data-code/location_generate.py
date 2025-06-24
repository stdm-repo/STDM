# This code is to generate location ID for manhattan 

import geopandas as gpd

df=gpd.read_file('data/taxi_zones/taxi_zones.shp')

if df.crs != "EPSG:4326":
    df = df.to_crs(epsg=4326)

df['centroid'] = df.geometry.centroid

df['longitude'] = df.centroid.x
df['latitude'] = df.centroid.y

df1=df[df.borough=='Manhattan']
df2=df1.loc[:,['LocationID','longitude','latitude']]

df2.to_csv('data/location.csv',index=False, encoding='utf-8')