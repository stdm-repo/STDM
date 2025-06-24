library(arrow)
library(dplyr)

##########################
id_list=read.csv("data/location.csv")
df2019=c()
for(i in 1:12){
  if(i<10){
    df1=read_parquet(paste("data/yellow_tripdata_2019-0",i,".parquet",sep=""))
  }else{
    df1=read_parquet(paste("data/yellow_tripdata_2019-",i,".parquet",sep=""))
  }
  df1_new=df1 %>%
    group_by(PULocationID) %>%
    summarise(count = n(),
              avg_customer = mean(passenger_count, na.rm = TRUE),
              avg_mile = mean(trip_distance, na.rm = TRUE))
  df1_new=as.data.frame(df1_new)
  if(i<10){
    df1_new$time=paste("20190",i,sep="")
  }else{
    df1_new$time=paste("2019",i,sep="")
  }
  df1_new=merge(id_list,df1_new,by.x="LocationID",by.y="PULocationID")
  df2019=rbind(df2019,df1_new)
  print(i)
}

df2018=c()
for(i in 1:12){
  if(i<10){
    df1=read_parquet(paste("data/yellow_tripdata_2018-0",i,".parquet",sep=""))
  }else{
    df1=read_parquet(paste("data/yellow_tripdata_2018-",i,".parquet",sep=""))
  }
  df1_new=df1 %>%
    group_by(PULocationID) %>%
    summarise(count = n(),
              avg_customer = mean(passenger_count, na.rm = TRUE),
              avg_mile = mean(trip_distance, na.rm = TRUE))
  df1_new=as.data.frame(df1_new)
  if(i<10){
    df1_new$time=paste("20180",i,sep="")
  }else{
    df1_new$time=paste("2018",i,sep="")
  }
  df1_new=merge(id_list,df1_new,by.x="LocationID",by.y="PULocationID")
  df2018=rbind(df2018,df1_new)
  print(i)
}

df2017=c()
for(i in 1:12){
  if(i<10){
    df1=read_parquet(paste("data/yellow_tripdata_2017-0",i,".parquet",sep=""))
  }else{
    df1=read_parquet(paste("data/yellow_tripdata_2017-",i,".parquet",sep=""))
  }
  df1_new=df1 %>%
    group_by(PULocationID) %>%
    summarise(count = n(),
              avg_customer = mean(passenger_count, na.rm = TRUE),
              avg_mile = mean(trip_distance, na.rm = TRUE))
  df1_new=as.data.frame(df1_new)
  if(i<10){
    df1_new$time=paste("20170",i,sep="")
  }else{
    df1_new$time=paste("2017",i,sep="")
  }
  df1_new=merge(id_list,df1_new,by.x="LocationID",by.y="PULocationID")
  df2017=rbind(df2017,df1_new)
  print(i)
}

df2016=c()
for(i in 1:12){
  if(i<10){
    df1=read_parquet(paste("data/yellow_tripdata_2016-0",i,".parquet",sep=""))
  }else{
    df1=read_parquet(paste("data/yellow_tripdata_2016-",i,".parquet",sep=""))
  }
  df1_new=df1 %>%
    group_by(PULocationID) %>%
    summarise(count = n(),
              avg_customer = mean(passenger_count, na.rm = TRUE),
              avg_mile = mean(trip_distance, na.rm = TRUE))
  df1_new=as.data.frame(df1_new)
  if(i<10){
    df1_new$time=paste("20160",i,sep="")
  }else{
    df1_new$time=paste("2016",i,sep="")
  }
  df1_new=merge(id_list,df1_new,by.x="LocationID",by.y="PULocationID")
  df2016=rbind(df2016,df1_new)
  print(i)
}

df2015=c()
for(i in 1:12){
  if(i<10){
    df1=read_parquet(paste("data/yellow_tripdata_2015-0",i,".parquet",sep=""))
  }else{
    df1=read_parquet(paste("data/yellow_tripdata_2015-",i,".parquet",sep=""))
  }
  df1_new=df1 %>%
    group_by(PULocationID) %>%
    summarise(count = n(),
              avg_customer = mean(passenger_count, na.rm = TRUE),
              avg_mile = mean(trip_distance, na.rm = TRUE))
  df1_new=as.data.frame(df1_new)
  if(i<10){
    df1_new$time=paste("20150",i,sep="")
  }else{
    df1_new$time=paste("2015",i,sep="")
  }
  df1_new=merge(id_list,df1_new,by.x="LocationID",by.y="PULocationID")
  df2015=rbind(df2015,df1_new)
  print(i)
}

df2014=c()
for(i in 1:12){
  if(i<10){
    df1=read_parquet(paste("data/yellow_tripdata_2014-0",i,".parquet",sep=""))
  }else{
    df1=read_parquet(paste("data/yellow_tripdata_2014-",i,".parquet",sep=""))
  }
  df1_new=df1 %>%
    group_by(PULocationID) %>%
    summarise(count = n(),
              avg_customer = mean(passenger_count, na.rm = TRUE),
              avg_mile = mean(trip_distance, na.rm = TRUE))
  df1_new=as.data.frame(df1_new)
  if(i<10){
    df1_new$time=paste("20140",i,sep="")
  }else{
    df1_new$time=paste("2014",i,sep="")
  }
  df1_new=merge(id_list,df1_new,by.x="LocationID",by.y="PULocationID")
  df2014=rbind(df2014,df1_new)
  print(i)
}

df2013=c()
for(i in 1:12){
  if(i<10){
    df1=read_parquet(paste("data/yellow_tripdata_2013-0",i,".parquet",sep=""))
  }else{
    df1=read_parquet(paste("data/yellow_tripdata_2013-",i,".parquet",sep=""))
  }
  df1_new=df1 %>%
    group_by(PULocationID) %>%
    summarise(count = n(),
              avg_customer = mean(passenger_count, na.rm = TRUE),
              avg_mile = mean(trip_distance, na.rm = TRUE))
  df1_new=as.data.frame(df1_new)
  if(i<10){
    df1_new$time=paste("20130",i,sep="")
  }else{
    df1_new$time=paste("2013",i,sep="")
  }
  df1_new=merge(id_list,df1_new,by.x="LocationID",by.y="PULocationID")
  df2013=rbind(df2013,df1_new)
  print(i)
}

df2012=c()
for(i in 1:12){
  if(i<10){
    df1=read_parquet(paste("data/yellow_tripdata_2012-0",i,".parquet",sep=""))
  }else{
    df1=read_parquet(paste("data/yellow_tripdata_2012-",i,".parquet",sep=""))
  }
  df1_new=df1 %>%
    group_by(PULocationID) %>%
    summarise(count = n(),
              avg_customer = mean(passenger_count, na.rm = TRUE),
              avg_mile = mean(trip_distance, na.rm = TRUE))
  df1_new=as.data.frame(df1_new)
  if(i<10){
    df1_new$time=paste("20120",i,sep="")
  }else{
    df1_new$time=paste("2012",i,sep="")
  }
  df1_new=merge(id_list,df1_new,by.x="LocationID",by.y="PULocationID")
  df2012=rbind(df2012,df1_new)
  print(i)
}

df2011=c()
for(i in 1:12){
  if(i<10){
    df1=read_parquet(paste("data/yellow_tripdata_2011-0",i,".parquet",sep=""))
  }else{
    df1=read_parquet(paste("data/yellow_tripdata_2011-",i,".parquet",sep=""))
  }
  df1_new=df1 %>%
    group_by(PULocationID) %>%
    summarise(count = n(),
              avg_customer = mean(passenger_count, na.rm = TRUE),
              avg_mile = mean(trip_distance, na.rm = TRUE))
  df1_new=as.data.frame(df1_new)
  if(i<10){
    df1_new$time=paste("20110",i,sep="")
  }else{
    df1_new$time=paste("2011",i,sep="")
  }
  df1_new=merge(id_list,df1_new,by.x="LocationID",by.y="PULocationID")
  df2011=rbind(df2011,df1_new)
  print(i)
}

data=rbind(df2011,df2012,df2013,df2014,df2015,df2016,df2017,df2018,df2019)
write.csv(data, "data/combined_data.csv", row.names = FALSE)

