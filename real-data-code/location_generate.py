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

# Change the three LocationIDs with value 103 to 103, 104, 105 in order
count_103 = 0
for i in range(len(df2)):
    if df2.iloc[i]['LocationID'] == 103:
        count_103 += 1
        if count_103 == 1:
            df2.iloc[i, df2.columns.get_loc('LocationID')] = 103
        elif count_103 == 2:
            df2.iloc[i, df2.columns.get_loc('LocationID')] = 104
        elif count_103 == 3:
            df2.iloc[i, df2.columns.get_loc('LocationID')] = 105

df2.to_csv('data/location.csv',index=False, encoding='utf-8')
