# data prep ----
for (plant_code in plants){
  # reading ----
  gfs_data=fread(paste0(data_path, "/", plant_code, "/", plant_code, '_GFS.csv'))
  arpege_data=fread(paste0(data_path, "/", plant_code, "/", plant_code, '_ARPEGE.csv'))
  production_data=fread(paste0(data_path, "/", plant_code, "/", plant_code,'_production.csv'))
  production_data$corrected_production = ifelse(production_data$corrected_production < 0, 0, production_data$corrected_production)
  gfs_data[,c('typ'):=NULL]
  arpege_data[,c('typ'):=NULL]
  production_data[,c('production', "dt", "max_production", "reliable", "cap_usage"):=NULL]
  
  # take common interval & merge production ----
  
  gfs_data2=gfs_data
  arpege_data2=arpege_data
  
  min = max(min(gfs_data$date), min(arpege_data$date))
  max = min(max(gfs_data$date), max(arpege_data$date))
  
  gfs_data=gfs_data[date>=min & date<=max]
  arpege_data=arpege_data[date>=min & date<=max]
  
  wide_data_gfs <- dcast(gfs_data[variable %in% c('VGRD','UGRD')],date+hour+lat+lon+level~variable,value.var='value')
  wide_data_gfs[,ws:=sqrt(VGRD^2+UGRD^2)]
  wide_data_gfs <- merge(wide_data_gfs,production_data,by=c('date','hour')) %>% drop_na()
  
  wide_data_arpege <- dcast(arpege_data[variable %in% c('VGRD','UGRD')],date+hour+lat+lon+level~variable,value.var='value')
  wide_data_arpege[,ws:=sqrt(VGRD^2+UGRD^2)]
  wide_data_arpege <- merge(wide_data_arpege,production_data,by=c('date','hour')) %>% drop_na()
  
  # all available time ----
  wide_data_gfs2 <- dcast(gfs_data2[variable %in% c('VGRD','UGRD')],date+hour+lat+lon+level~variable,value.var='value')
  wide_data_gfs2[,ws:=sqrt(VGRD^2+UGRD^2)]
  wide_data_gfs2 <- merge(wide_data_gfs2,production_data,by=c('date','hour')) %>% drop_na()
  
  wide_data_arpege2 <- dcast(arpege_data2[variable %in% c('VGRD','UGRD')],date+hour+lat+lon+level~variable,value.var='value')
  wide_data_arpege2[,ws:=sqrt(VGRD^2+UGRD^2)]
  wide_data_arpege2 <- merge(wide_data_arpege2,production_data,by=c('date','hour')) %>% drop_na()
  
  # save_rm ----
  save_rm_object(gfs_data, output_path=output_path, code=plant_code)
  save_rm_object(wide_data_gfs, output_path=output_path, code=plant_code)
  save_rm_object(arpege_data, output_path=output_path, code=plant_code)
  save_rm_object(wide_data_arpege, output_path=output_path, code=plant_code)
  save_rm_object(gfs_data2, output_path=output_path, code=plant_code)
  save_rm_object(wide_data_gfs2, output_path=output_path, code=plant_code)
  save_rm_object(arpege_data2, output_path=output_path, code=plant_code)
  save_rm_object(wide_data_arpege2, output_path=output_path, code=plant_code)
}
