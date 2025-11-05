# AHOG.O data ----
nbins = c(6, 9, 12, 18, 24, 36)
nwp_set = c("gfs", "arpege")
for (plant_code in plants){
  for (nwp in nwp_set){
    qr = readRDS(file = paste0(output_path, "/", plant_code, "/", "data", "/", paste0("data.", nwp,".topws.qr2.1"))) %>% select(-topws)
    
    HOG_data_ws <- readRDS(file = paste0(output_path, "/", plant_code, "/", "data", "/", "wide_data_",nwp,"_ws")) %>% mutate(angle=atan2(mean_v, mean_u)* 180/pi) %>% select(date, hour, corrected_production, ws, angle)
    tohog.qr = merge(qr, HOG_data_ws, by=c("date", "hour", "corrected_production"))
    
    for (nb in nbins){
      name = paste0("data.", nwp,".agghog", as.character(nb))
      
      data.agghog.1 <- bins_circular0(HOG_data_ws, nbin=nb)
      data.agghog.1 <- setNames(data.agghog.1, c(names(data.agghog.1[,c(1,2,3)]), paste0("bin_", names(data.agghog.1[,-c(1,2,3)]))))
      poly1 <- setNames(data.agghog.1[,-c(1,2,3)]^2, c(paste0(names(data.agghog.1[,-c(1,2,3)]),'_2')))
      data.agghog.2 <- as.data.frame(cbind(data.agghog.1, poly1))
      poly2 <- setNames(data.agghog.1[,-c(1,2,3)]^3, c(paste0(names(data.agghog.1[,-c(1,2,3)]),'_3')))
      data.agghog.3 <- as.data.frame(cbind(data.agghog.2, poly2))
      
      for(p in c(3)){
        namewanted = paste0(name, ".qr2.", p)
        objname = paste0("data.agghog.", p)
        save_rm_wname(namewanted, objname, output_path=output_path, code=plant_code)
      }
    }
    
  }
  
  
  for (nb in nbins){
    data.all.agghog.3 <- merge(readRDS(file = paste0(output_path, "/", plant_code, "/", "data", "/", paste0("data.gfs.agghog", as.character(nb), ".qr2.3"))), readRDS(file = paste0(output_path, "/", plant_code, "/", "data", "/", paste0("data.arpege.agghog", as.character(nb), ".qr2.3"))), by=c('date','hour', "corrected_production"), suffixes=c(".gfs", ".arpege"))
    save_rm_wname(paste0("data.all.agghog", as.character(nb), ".qr2.3") , "data.all.agghog.3", output_path=output_path, code=plant_code)
  }
  
}

# AHOG.O val utilization model glm ----
nbins = c(6, 9, 12, 18, 24, 36)
nwp = "all"
hog_data_names = paste0("data.",nwp,".agghog", nbins, ".qr2.3")
for (plant_code in plants){
  for (dataname in hog_data_names){
    
    data=readRDS(file=paste0(output_path, "/", plant_code, "/", "data", "/", dataname)) %>% arrange(date, hour)
    colnames(data) = make.names(names(data))
    capacity = max(data$corrected_production)
    data$corrected_production = data$corrected_production / capacity
    train = test_train(data)[2][[1]]
    registerDoMC(cores = 5)
    set.seed(1)
    mymodel <- cv.glmnet(x=as.matrix(train[,-c(1,2,3)]), y=train$corrected_production, alpha=1, nfolds=5, family=quasibinomial(), parallel=TRUE)
    name_model = paste0("model.glm.cv.", str_sub(dataname,6,-3), ".r.3")
    saveRDS(mymodel, file = paste0(output_path, "/", plant_code, "/", "models", "/", name_model))
    
    
  }
}

# AHOG.O val results ----
nbins = c(6, 9, 12, 18, 24, 36)
nwp = "all"
hog_data_names = paste0("data.",nwp,".agghog", nbins, ".qr2.3")
outputs.sensitivity.train = data.frame()
outputs.sensitivity.test = data.frame()
for (plant_code in plants){
  train_outputs = c()
  test_outputs = c()
  for (dataname in hog_data_names){
    data=readRDS(file=paste0(output_path, "/", plant_code, "/", "data", "/", dataname)) %>% arrange(date, hour)
    colnames(data) = make.names(names(data))
    capacity = max(data$corrected_production)
    data$corrected_production = data$corrected_production / capacity
    train = test_train(data)[2][[1]]
    test = test_train(data)[1][[1]]
    pred_test = test[,c(1:3)]
    pred_train = train[,c(1:3)]
    
    name_model = paste0("model.glm.cv.", str_sub(dataname,6,-3), ".r.3")
    mymodel=readRDS(file=paste0(output_path, "/", plant_code, "/", "models", "/", name_model))
    pred_test$pred = predict(mymodel, type="response", as.matrix(test[,-c(1,2,3)]))
    pred_test$pred = ifelse(pred_test$pred > capacity, capacity, pred_test$pred)
    pred_test$pred = ifelse(pred_test$pred < 0, 0, pred_test$pred)
    pred_train$pred = predict(mymodel, type="response", as.matrix(train[,-c(1,2,3)]))
    pred_train$pred = ifelse(pred_train$pred > capacity, capacity, pred_train$pred)
    pred_train$pred = ifelse(pred_train$pred < 0, 0, pred_train$pred)
    
    result_test = wmape(pred_test$pred, pred_test$corrected_production)*100
    result_train = wmape(pred_train$pred, pred_train$corrected_production)*100
    train_outputs = append(train_outputs, result_train)
    test_outputs = append(test_outputs, result_test)
  }
  outputs.sensitivity.train = rbind(outputs.sensitivity.train, train_outputs)
  outputs.sensitivity.test = rbind(outputs.sensitivity.test, test_outputs)
}
colnames(outputs.sensitivity.train) = paste0("bin", nbins)
colnames(outputs.sensitivity.test) = paste0("bin", nbins)
outputs.sensitivity.train$min_bin  <- apply(outputs.sensitivity.train, 1, function(x) colnames(outputs.sensitivity.train)[which.min(x)])
outputs.sensitivity.train$min_bin = substring(outputs.sensitivity.train$min_bin, 4)
fwrite(outputs.sensitivity.train, file=paste0(output_path, "/", "outputs.",nwp,".AHog.sensitivity.qr2.r.train.cv", ".csv"))
fwrite(outputs.sensitivity.test, file=paste0(output_path, "/", "outputs.",nwp,".AHog.sensitivity.qr2.r.test.cv", ".csv"))

nbins = c(6, 9, 12, 18, 24, 36)

bin_sensitivity = read.csv(paste0(output_path, "/","outputs.all.AHog.sensitivity.qr2.r.train.cv.csv")) %>% select(-min_bin)
params = c()
for (r in 1:47){
  train = as.numeric(bin_sensitivity[r,])
  dif = diff(train, lag=1)
  nbin=6
  for(c in 1:5){
    if (dif[c]<0){
      nbin = nbins[c+1]
    }else{
      break
    }
  }
  params = c(params, nbin)
}
bin_sensitivity$min_bin = params
fwrite(bin_sensitivity, file=paste0(output_path, "/", "outputs.","all",".AHog.sensitivity.qr2.r.train.cv", ".csv"))


# AHOG.O results ----

nwp = "all"
bin_sensitivity = read.csv(paste0(output_path, "/","outputs.all.AHog.sensitivity.qr2.r.train.cv.csv"))
hogbins = bin_sensitivity$min_bin


for (i in 1:47){
  plant_code = plants[i]
  nbin = hogbins[i]
  
  dataname = paste0("data.",nwp,".agghog", nbin, ".qr2.3")
  data=readRDS(file=paste0(output_path, "/", plant_code, "/", "data", "/", dataname)) %>% arrange(date, hour)
  colnames(data) = make.names(names(data))
  capacity = max(data$corrected_production)
  data$corrected_production = data$corrected_production / capacity
  train = test_train(data)[2][[1]]
  test = test_train(data)[1][[1]]
  pred_test = test[,c(1:3)]
  pred_train = train[,c(1:3)]
  
  name_model = paste0("model.glm.cv.", str_sub(dataname,6,-3), ".r.3")
  mymodel=readRDS(file=paste0(output_path, "/", plant_code, "/", "models", "/", name_model))
  pred_test$pred = predict(mymodel, type="response", as.matrix(test[,-c(1,2,3)]))
  pred_test$pred = ifelse(pred_test$pred > capacity, capacity, pred_test$pred)
  pred_test$pred = ifelse(pred_test$pred < 0, 0, pred_test$pred)
  pred_train$pred = predict(mymodel, type="response", as.matrix(train[,-c(1,2,3)]))
  pred_train$pred = ifelse(pred_train$pred > capacity, capacity, pred_train$pred)
  pred_train$pred = ifelse(pred_train$pred < 0, 0, pred_train$pred)
  
  fwrite(pred_test, file=paste0(output_path, "/", plant_code, "/", "preds", "/", "test.all.AHog", ".qr2.r.glm"))
  fwrite(pred_train, file=paste0(output_path, "/", plant_code, "/", "preds", "/", "train.all.AHog", ".qr2.r.glm"))
  
}




# LHOG.O server1 data ----
nbins = c(6, 9, 12, 18, 24, 36)

nwp_set = c("gfs", "arpege")
for (plant_code in plants){
  for (nwp in nwp_set){
    ws.qr.1 = readRDS(file = paste0(output_path, "/", plant_code, "/", "data", "/", paste0("data.", nwp,".topws.qr2.1"))) %>% select(-topws)
    
    HOG_data_loc = readRDS(file = paste0(output_path, "/", plant_code, "/", "data", "/", "wide_data_",nwp)) %>% mutate(angle=atan2(VGRD, UGRD)* 180/pi)
    HOG_data_loc$loc = paste0(HOG_data_loc$lat, "-", HOG_data_loc$lon)
    HOG_data_loc = HOG_data_loc %>% select(date, hour, corrected_production, loc, ws, angle)
    tohog.qr.1 = merge(ws.qr.1, HOG_data_loc, by=c("date", "hour", "corrected_production"))
    
    
    for (nb in nbins){
      dat <- bins_circular00(tohog.qr.1, nbin=nb)
      hog.qr.1 = reshape(data.table(dat), idvar = c("date", "hour", "corrected_production"), timevar = "loc", direction="wide") # no fill=0
      hog.qr.1 <- setNames(hog.qr.1, c(names(hog.qr.1[,c(1,2,3)]), paste0("bin_", names(hog.qr.1[,-c(1,2,3)]))))
      poly1 <- setNames(hog.qr.1[,-c(1,2, 3)]^2, c(paste0(names(hog.qr.1[,-c(1,2,3)]),'_2'))  )
      hog.qr.2 <- as.data.frame(cbind(hog.qr.1, poly1))
      poly2 <- setNames(hog.qr.1[,-c(1,2, 3)]^3, c(paste0(names(hog.qr.1[,-c(1,2,3)]),'_3'))  )
      hog.qr.3 <- as.data.frame(cbind(hog.qr.2, poly2))
      
      save_rm_wname(paste0("data.", nwp,".lochog", as.character(nb), ".qr2.3"), "hog.qr.3", output_path=output_path, code=plant_code)
      
    }
    
  }
  
  
  for (nb in nbins){
    data.all.lochog.3 <- merge(readRDS(file = paste0(output_path, "/", plant_code, "/", "data", "/", paste0("data.gfs.lochog", as.character(nb), ".qr2.3"))), readRDS(file = paste0(output_path, "/", plant_code, "/", "data", "/", paste0("data.arpege.lochog", as.character(nb), ".qr2.3"))), by=c('date','hour', "corrected_production"), suffixes=c(".gfs", ".arpege"))
    save_rm_wname(paste0("data.all.lochog", as.character(nb), ".qr2.3") , "data.all.lochog.3", output_path=output_path, code=plant_code)
  }
  
}

# LHOG.O val utilization model glm ----
nbins = c(6, 9, 12, 18, 24, 36)
nwp = "all"
hog_data_names = paste0("data.",nwp,".lochog", nbins, ".qr2.3")
library(doParallel)
doParallel::registerDoParallel(cores = 5)
for (plant_code in plants[15:47]){
  for (dataname in hog_data_names){
    data=readRDS(file=paste0(output_path, "/", plant_code, "/", "data", "/", dataname)) %>% arrange(date, hour)
    colnames(data) = make.names(names(data))
    capacity = max(data$corrected_production)
    data$corrected_production = data$corrected_production / capacity
    train = test_train(data)[2][[1]]
    set.seed(1)
    mymodel <- cv.glmnet(x=as.matrix(train[,-c(1,2,3)]), y=train$corrected_production, alpha=1, nfolds=5, family=quasibinomial(), parallel=TRUE)
    name_model = paste0("model.glm.cv.", str_sub(dataname,6,-3), ".r.3")
    saveRDS(mymodel, file = paste0(output_path, "/", plant_code, "/", "models", "/", name_model))
  }
}
# LHOG.O val results ----

nbins = c(6, 9, 12, 18, 24, 36)
nwp = "all"
hog_data_names = paste0("data.",nwp,".lochog", nbins, ".qr2.3")
outputs.sensitivity.train = data.frame()
outputs.sensitivity.test = data.frame()
for (plant_code in plants){
  train_outputs = c()
  test_outputs = c()
  for (dataname in hog_data_names){
    data=readRDS(file=paste0(output_path, "/", plant_code, "/", "data", "/", dataname)) %>% arrange(date, hour)
    colnames(data) = make.names(names(data))
    capacity = max(data$corrected_production)
    data$corrected_production = data$corrected_production / capacity
    train = test_train(data)[2][[1]]
    test = test_train(data)[1][[1]]
    pred_test = test[,c(1:3)]
    pred_train = train[,c(1:3)]
    
    name_model = paste0("model.glm.cv.", str_sub(dataname,6,-3), ".r.3")
    mymodel=readRDS(file=paste0(output_path, "/", plant_code, "/", "models", "/", name_model))
    pred_test$pred = predict(mymodel, type="response", as.matrix(test[,-c(1,2,3)]))
    pred_test$pred = ifelse(pred_test$pred > capacity, capacity, pred_test$pred)
    pred_test$pred = ifelse(pred_test$pred < 0, 0, pred_test$pred)
    pred_train$pred = predict(mymodel, type="response", as.matrix(train[,-c(1,2,3)]))
    pred_train$pred = ifelse(pred_train$pred > capacity, capacity, pred_train$pred)
    pred_train$pred = ifelse(pred_train$pred < 0, 0, pred_train$pred)
    
    result_test = wmape(pred_test$pred, pred_test$corrected_production)*100
    result_train = wmape(pred_train$pred, pred_train$corrected_production)*100
    train_outputs = append(train_outputs, result_train)
    test_outputs = append(test_outputs, result_test)
  }
  outputs.sensitivity.train = rbind(outputs.sensitivity.train, train_outputs)
  outputs.sensitivity.test = rbind(outputs.sensitivity.test, test_outputs)
}
colnames(outputs.sensitivity.train) = paste0("bin", nbins)
colnames(outputs.sensitivity.test) = paste0("bin", nbins)
outputs.sensitivity.train$min_bin  <- apply(outputs.sensitivity.train, 1, function(x) colnames(outputs.sensitivity.train)[which.min(x)])
outputs.sensitivity.train$min_bin = substring(outputs.sensitivity.train$min_bin, 4)
fwrite(outputs.sensitivity.train, file=paste0(output_path, "/", "outputs.",nwp,".LHog.sensitivity.qr2.r.train.cv", ".csv"))
fwrite(outputs.sensitivity.test, file=paste0(output_path, "/", "outputs.",nwp,".LHog.sensitivity.qr2.r.test.cv", ".csv"))

bin_sensitivity = read.csv(paste0(output_path, "/","outputs.all.LHog.sensitivity.qr2.r.train.cv.csv")) %>% select(-min_bin)
params = c()
for (r in 1:47){
  train = as.numeric(bin_sensitivity[r,])
  dif = diff(train, lag=1)
  nbin=6
  for(c in 1:5){
    if (dif[c]<0){
      nbin = nbins[c+1]
    }else{
      break
    }
  }
  params = c(params, nbin)
}
bin_sensitivity$min_bin = params
fwrite(bin_sensitivity, file=paste0(output_path, "/", "outputs.","all",".LHog.sensitivity.qr2.r.train.cv", ".csv"))


# LHOG.O results ----
nwp = "all"
bin_sensitivity = read.csv(paste0(output_path, "/","outputs.all.LHog.sensitivity.qr2.r.train.cv.csv"))
hogbins = bin_sensitivity$min_bin
for (i in 1:47){
  plant_code = plants[i]
  nbin = hogbins[i]
  dataname = paste0("data.",nwp,".lochog", nbin, ".qr2.3")
  data=readRDS(file=paste0(output_path, "/", plant_code, "/", "data", "/", dataname)) %>% arrange(date, hour)
  colnames(data) = make.names(names(data))
  capacity = max(data$corrected_production)
  data$corrected_production = data$corrected_production / capacity
  train = test_train(data)[2][[1]]
  test = test_train(data)[1][[1]]
  pred_test = test[,c(1:3)]
  pred_train = train[,c(1:3)]
  
  name_model = paste0("model.glm.cv.", str_sub(dataname,6,-3), ".r.3")
  mymodel=readRDS(file=paste0(output_path, "/", plant_code, "/", "models", "/", name_model))
  pred_test$pred = predict(mymodel, type="response", as.matrix(test[,-c(1,2,3)]))
  pred_test$pred = ifelse(pred_test$pred > capacity, capacity, pred_test$pred)
  pred_test$pred = ifelse(pred_test$pred < 0, 0, pred_test$pred)
  pred_train$pred = predict(mymodel, type="response", as.matrix(train[,-c(1,2,3)]))
  pred_train$pred = ifelse(pred_train$pred > capacity, capacity, pred_train$pred)
  pred_train$pred = ifelse(pred_train$pred < 0, 0, pred_train$pred)
  
  fwrite(pred_test, file=paste0(output_path, "/0predictions/preds/", paste0("test.preds_",plant_code,"_LHog.qr2.r.csv")))
  fwrite(pred_train, file=paste0(output_path, "/0predictions/preds/", paste0("train.preds_",plant_code,"_LHog.qr2.r.csv")))
  
}
# AwB.O data wrt val ----

bin_sensitivity = read.csv(paste0(output_path, "/","outputs.all.AHog.sensitivity.qr2.r.train.cv.csv"))
hogbins = bin_sensitivity$min_bin

nwp_set = c("gfs", "arpege")
for (i in 1:length(plants)){
  plant_code = plants[i]
  nbin = hogbins[i]
  
  for (nwp in nwp_set){
    
    qr = readRDS(file = paste0(output_path, "/", plant_code, "/", "data", "/", paste0("data.", nwp,".topws.qr2.1"))) %>% select(-topws)
    
    
    wide.aggwsangle <- readRDS(file = paste0(output_path, "/", plant_code, "/", "data", "/", "wide_data_",nwp,"_ws")) %>% mutate(angle=atan2(mean_v, mean_u)* 180/pi) %>% select(date, hour, corrected_production, ws, angle)
    wide.aggwsangle = merge(qr, wide.aggwsangle, by=c("date", "hour", "corrected_production"))
    
    small_number = 0.01
    max_ang = 180
    min_ang = -180
    
    brk_points=seq(min_ang-small_number,max_ang+small_number,length=nbin+1)
    wide.aggwsangle=wide.aggwsangle[order(angle)]
    wide.aggwsangle[,bin_id:=cut(angle,brk_points)]
    wide.aggwsangle[,bin_number:=as.numeric(bin_id)]
    
    data.binaggws.1 = dcast(data.table(wide.aggwsangle),date+hour+corrected_production~bin_number,value.var="ws",fill=0)
    
    poly1 <- setNames(data.binaggws.1[,-c(1,2,3)]^2, c(paste0(names(data.binaggws.1[,-c(1,2,3)]),'_2')))
    data.binaggws.2 = as.data.frame(cbind(data.binaggws.1, poly1))
    poly2 <- setNames(data.binaggws.1[,-c(1,2,3)]^3, c(paste0(names(data.binaggws.1[,-c(1,2,3)]),'_3')))
    data.binaggws.3 = as.data.frame(cbind(data.binaggws.2, poly2))
    
    save_rm_wname(paste0("data.",nwp,".binaggws.qr2.3"), paste0("data.binaggws.3"), folder="data", output_path=output_path, code=plant_code)
    
  }
  
  data.all.binaggws.qr2.3 <- merge(readRDS(file = paste0(output_path, "/", plant_code, "/", "data", "/", "data.gfs.binaggws.qr2.3")), readRDS(file = paste0(output_path, "/", plant_code, "/", "data", "/", "data.arpege.binaggws.qr2.3")), by=c('date','hour', "corrected_production"), suffixes=c(".gfs", ".arpege"))
  save_rm_name("data.all.binaggws.qr2.3", output_path=output_path, code=plant_code)
  
  
}



# AwB.O utilization model glm wrt val ----

for (i in 1:length(plants)){
  plant_code = plants[i]
  
  dataname = "data.all.binaggws.qr2.3"
  data=readRDS(file=paste0(output_path, "/", plant_code, "/", "data", "/", dataname)) %>% arrange(date, hour)
  colnames(data) = make.names(names(data))
  capacity = max(data$corrected_production)
  data$corrected_production = data$corrected_production / capacity
  train = test_train(data)[2][[1]]
  registerDoParallel(cores = 5)
  set.seed(1)
  mymodel <- cv.glmnet(x=as.matrix(train[,-c(1,2,3)]), y=train$corrected_production, alpha=1, nfolds=5, family=quasibinomial(), parallel=TRUE)
  name_model = paste0("model.glm.cv.", str_sub(dataname,6,-3), ".r.3")
  saveRDS(mymodel, file = paste0(output_path, "/", plant_code, "/", "models", "/", name_model))
  
  
}


nwp_set = c("all")
for (nwp in nwp_set){
  outputs = data.frame()
  for (i in 1:length(plants)){
    plant_code=plants[i]
    dataname = paste0("data.", "all",".binaggws", ".qr2.3")
    data=readRDS(file=paste0(output_path, "/", plant_code, "/", "data", "/", dataname)) %>% arrange(date, hour)
    colnames(data) = make.names(names(data))
    capacity = max(data$corrected_production)
    data$corrected_production = data$corrected_production / capacity
    test = test_train(data)[1][[1]]
    pred_test = test[,c(1:3)]
    
    name_model = paste0("model.glm.cv.", str_sub(dataname,6,-3), ".r.3")
    mymodel=readRDS(file=paste0(output_path, "/", plant_code, "/", "models", "/", name_model))
    pred_test$pred = predict(mymodel, type="response", as.matrix(test[,-c(1,2,3)]))
    fwrite(pred_test, file=paste0(output_path, "/", plant_code, "/", "preds", "/", "test.all.BAw.qr2.r", ".glm"))
    
  }
  
}

# LwB.O data wrt val ----

bin_sensitivity = read.csv(paste0(output_path, "/","outputs.all.LHog.sensitivity.qr2.r.train.cv.csv"))
hogbins = bin_sensitivity$min_bin


nwp_set = c("gfs", "arpege")
for (i in 1:length(plants)){
  plant_code = plants[i]
  nbin = hogbins[i]
  
  for (nwp in nwp_set){
    
    ws.qr.1=readRDS(file = paste0(output_path, "/", plant_code, "/", "data", "/", "data.",nwp,".aggws.qr2.1"))
    wide.locwsangle <- readRDS(file = paste0(output_path, "/", plant_code, "/", "data", "/", "wide_data_", nwp)) %>% mutate(angle=atan2(VGRD, UGRD)* 180/pi)
    wide.locwsangle$loc = paste0(wide.locwsangle$lat, "-", wide.locwsangle$lon)
    wide.locwsangle = wide.locwsangle %>% select(-lat, -lon, -level, -UGRD, -VGRD)
    
    colnames(ws.qr.1) = c(colnames(ws.qr.1)[1:3], "aw")
    tohog.qr.1 = merge(ws.qr.1, wide.locwsangle, by=c("date", "hour", "corrected_production")) %>% select(-aw)
    
    small_number = 0.01
    max_ang = 180
    min_ang = -180
    
    brk_points=seq(min_ang-small_number,max_ang+small_number,length=nbin+1)
    tohog.qr.1=tohog.qr.1[order(angle)]
    tohog.qr.1[,bin_id:=cut(angle,brk_points)]
    tohog.qr.1[,bin_number:=as.numeric(bin_id)]
    
    data.binlocws.1 = dcast(data.table(tohog.qr.1), date+hour+corrected_production~bin_number+loc,value.var="ws",fill=0)
    
    poly1 <- setNames(data.binlocws.1[,-c(1,2,3)]^2, c(paste0(names(data.binlocws.1[,-c(1,2,3)]),'_2')))
    data.binlocws.2 = as.data.frame(cbind(data.binlocws.1, poly1))
    poly2 <- setNames(data.binlocws.1[,-c(1,2,3)]^3, c(paste0(names(data.binlocws.1[,-c(1,2,3)]),'_3')))
    data.binlocws.3 = as.data.frame(cbind(data.binlocws.2, poly2))
    
    save_rm_wname(paste0("data.", nwp,".binlocws.qr2.3"), "data.binlocws.3", output_path=output_path, code=plant_code)
  }
  
  data.all.binlocws.3 <- merge(readRDS(file = paste0(output_path, "/", plant_code, "/", "data", "/", paste0("data.gfs.binlocws", ".qr2.3"))), readRDS(file = paste0(output_path, "/", plant_code, "/", "data", "/", paste0("data.arpege.binlocws", ".qr2.3"))), by=c('date','hour', "corrected_production"), suffixes=c(".gfs", ".arpege"))
  save_rm_wname(paste0("data.all.binlocws", ".qr2.3") , "data.all.binlocws.3", output_path=output_path, code=plant_code)
  
  
}



# LwB.O utilization model glm wrt val ----

bin_sensitivity = read.csv(paste0(output_path, "/","outputs.all.LHog.sensitivity.qr2.r.train.cv.csv"))
hogbins = bin_sensitivity$min_bin

library(doParallel)
doParallel::registerDoParallel(cores = 5)
for (i in 1:length(plants)){
  plant_code = plants[i]
  dataname = "data.all.binlocws.qr2.3"
  data=readRDS(file=paste0(output_path, "/", plant_code, "/", "data", "/", dataname)) %>% arrange(date, hour)
  colnames(data) = make.names(names(data))
  capacity = max(data$corrected_production)
  data$corrected_production = data$corrected_production / capacity
  train = test_train(data)[2][[1]]
  set.seed(1)
  mymodel <- cv.glmnet(x=as.matrix(train[,-c(1,2,3)]), y=train$corrected_production, alpha=1, nfolds=5, family=quasibinomial(), parallel=TRUE)
  name_model = paste0("model.glm.cv.", str_sub(dataname,6,-3), ".r.3")
  saveRDS(mymodel, file = paste0(output_path, "/", plant_code, "/", "models", "/", name_model))
  
  
}

nwp_set = c("all")
for (nwp in nwp_set){
  outputs = data.frame()
  for (i in 1:length(plants)){
    plant_code=plants[i]
    dataname = paste0("data.", "all",".binlocws", ".qr2.3")
    data=readRDS(file=paste0(output_path, "/", plant_code, "/", "data", "/", dataname)) %>% arrange(date, hour)
    colnames(data) = make.names(names(data))
    capacity = max(data$corrected_production)
    data$corrected_production = data$corrected_production / capacity
    test = test_train(data)[1][[1]]
    pred_test = test[,c(1:3)]
    
    name_model = paste0("model.glm.cv.", str_sub(dataname,6,-3), ".r.3")
    mymodel=readRDS(file=paste0(output_path, "/", plant_code, "/", "models", "/", name_model))
    pred_test$pred = predict(mymodel, type="response", as.matrix(test[,-c(1,2,3)]))
    fwrite(pred_test, file=paste0(output_path, "/0predictions/preds/", paste0("test.preds_",plant_code,"_BLw.qr2.r.csv")))
    
  }
  
}


# LuvHOG.O.U data ----
tree="locuv"
nwp_set = c("gfs", "arpege")
depth_set = c(3, 4, 5, 6)
nbin_set = c(3, 6, 9)

for (plant_code in plants){
  for (nwp in nwp_set){
    wide.locuv <- readRDS(file = paste0(output_path, "/", plant_code, "/", "data", "/", "wide_data_",nwp))
    wide.locuv$loc = paste0(wide.locuv$lat, "-", wide.locuv$lon)
    wide.locuv = wide.locuv %>% select(-lat, -lon, -level)
    locuv0 = dcast(data.table(wide.locuv), date+hour+corrected_production~loc,value.var=c("UGRD", "VGRD")) # no fill=0
    
    ws.qr.1 = readRDS(file = paste0(output_path, "/", plant_code, "/", "data", "/", paste0("data.", nwp,".topws.qr2.1"))) %>% select(-topws)
    locuv = merge(ws.qr.1, locuv0, by=c("date", "hour", "corrected_production"))
    
    
    trees.locuv = make_trees_qr2(data=nwp, tree="locuv", depth=depth_set)
    save_rm_element_names(trees.locuv, "tree", output_path=output_path, code=plant_code)
    
    save_rm_wname(paste0("trees.", nwp, ".locuv.qr2"), "trees.locuv", folder="tree", output_path=output_path, code=plant_code)
    
    #save_tree_plots(data=nwp, tree="locuv", output_path=output_path, code=plant_code)
    
    put_leaves_qr2(data=nwp, tree="locuv", output_path=output_path, code=plant_code)
    
    save_rm_wname(paste0(nwp, ".locuv.qr2"), "locuv", output_path=output_path, code=plant_code)
  }
  
}

for (plant_code in plants){
  for (nwp in nwp_set){
    for(d in depth_set){
      for(n in nbin_set){
        make_adj_hog_bins_poly_qr2(data=nwp, tree=tree, depth=d, numbin=n, output_path=output_path, code=plant_code)
      }
    }
  }
  
  for(depth in depth_set){
    for(numbin in nbin_set){
      datanamestocombine = paste0("data.", nwp_set, ".", "tree.", tree, ".d", depth, ".b", numbin, ".qr2.3")
      data.treehog.3 <- merge(readRDS(file = paste0(output_path, "/", plant_code, "/", "data", "/", datanamestocombine[1])), readRDS(file = paste0(output_path, "/", plant_code, "/", "data", "/", datanamestocombine[2])), by=c('date','hour', "corrected_production"), suffixes=c(".gfs", ".arpege"))
      save_rm_wname(paste0("data.", "all", ".", "tree.", tree, ".d", depth, ".b", numbin, ".qr2.3"), "data.treehog.3",output_path=output_path, code=plant_code)
      
    }
  }
}


# LuvHOG.O.U val utilization model glm ----
nwp = "all"
doParallel::registerDoParallel(5)
for (plant_code in plants[1:47]){
  for(depth in depth_set){
    for(numbin in nbin_set){
      dataname = paste0("data.", nwp, ".tree.", tree, ".d", depth, ".b", numbin, ".qr2.3")
      data=readRDS(file=paste0(output_path, "/", plant_code, "/", "data", "/", dataname)) %>% arrange(date, hour)
      colnames(data) = make.names(names(data))
      capacity = max(data$corrected_production)
      data$corrected_production = data$corrected_production / capacity
      train = test_train(data)[2][[1]]
      #registerDoMC(cores = 5)
      set.seed(1)
      mymodel <- cv.glmnet(x=as.matrix(train[,-c(1,2,3)]), y=train$corrected_production, alpha=1, nfolds=5, family=quasibinomial(), parallel=TRUE, trace.it=1)
      name_model = paste0("model.glm.cv.", str_sub(dataname,6,-3), ".r.3")
      saveRDS(mymodel, file = paste0(output_path, "/", plant_code, "/", "models", "/", name_model))
    }
  }
}
# LuvHOG.O.U val results ----
depth_set = c(3, 4, 5, 6)
nbin_set = c(3, 6, 9)
tree="locuv"
for (plant_code in plants[1:47]){
  train_outputs = data.frame()
  test_outputs = data.frame()
  for(depth in depth_set){
    for(numbin in nbin_set){
      dataname = paste0("data.", nwp, ".tree.", tree, ".d", depth, ".b", numbin, ".qr2.3")
      data=readRDS(file=paste0(output_path, "/", plant_code, "/", "data", "/", dataname)) %>% arrange(date, hour)
      colnames(data) = make.names(names(data))
      capacity = max(data$corrected_production)
      data$corrected_production = data$corrected_production / capacity
      train = test_train(data)[2][[1]]
      test = test_train(data)[1][[1]]
      pred_test = test[,c(1:3)]
      pred_train = train[,c(1:3)]
      
      name_model = paste0("model.glm.cv.", str_sub(dataname,6,-3), ".r.3")
      mymodel=readRDS(file=paste0(output_path, "/", plant_code, "/", "models", "/", name_model))
      pred_test$pred = predict(mymodel, type="response", as.matrix(test[,-c(1,2,3)]))
      pred_test$pred = ifelse(pred_test$pred > capacity, capacity, pred_test$pred)
      pred_test$pred = ifelse(pred_test$pred < 0, 0, pred_test$pred)
      pred_train$pred = predict(mymodel, type="response", as.matrix(train[,-c(1,2,3)]))
      pred_train$pred = ifelse(pred_train$pred > capacity, capacity, pred_train$pred)
      pred_train$pred = ifelse(pred_train$pred < 0, 0, pred_train$pred)
      
      result_test_ = wmape(pred_test$pred, pred_test$corrected_production)*100
      result_train_ = wmape(pred_train$pred, pred_train$corrected_production)*100
      
      result_train = c(plant_code, depth, numbin, result_train_)
      result_test = c(plant_code, depth, numbin, result_test_)
      train_outputs = rbind(train_outputs, result_train)
      test_outputs = rbind(test_outputs, result_test)
    }
  }
  colnames(train_outputs) = c("Plant", "Depth", "Number of Bins", "Train")
  colnames(test_outputs) = c("Plant", "Depth", "Number of Bins", "Test")
  fwrite(train_outputs, file=paste0(output_path, "/", plant_code, "/", "outputs.",nwp,".tree",tree,".qr2.train.glm.cv", ".csv"))
  fwrite(test_outputs, file=paste0(output_path, "/", plant_code, "/", "outputs.",nwp,".tree",tree,".qr2.test.glm.cv", ".csv"))
}


# LuvHOG.O.U results ----
tree = "locuv"
depth_nbin = read.csv(file=paste0(output_path, "/", "outputs.all.tree",tree,".sensitivity.qr2.r.train", ".glm.cv.csv"))

for (i in 1:47){
  plant_code=plants[i]
  dataname = paste0("data.", nwp, ".tree.", tree, ".d", depth_nbin$Depth[i], ".b", depth_nbin$Number.of.Bins[i], ".qr2.3")
  data=readRDS(file=paste0(output_path, "/", plant_code, "/", "data", "/", dataname)) %>% arrange(date, hour)
  colnames(data) = make.names(names(data))
  capacity = max(data$corrected_production)
  data$corrected_production = data$corrected_production / capacity
  train = test_train(data)[2][[1]]
  
  name_model = paste0("model.glm.cv.", str_sub(dataname,6,-3), ".r.3")
  mymodel <- readRDS(file=paste0(output_path, "/", plant_code, "/", "models", "/", name_model))
  
  test = test_train(data)[1][[1]]
  pred_test = test[,c(1:3)]
  pred_test$pred = predict(mymodel, type="response", as.matrix(test[,-c(1,2,3)]))
  
  train = test_train(data)[2][[1]]
  pred_train = train[,c(1:3)]
  pred_train$pred = predict(mymodel, type="response", as.matrix(train[,-c(1,2,3)]))
  
  fwrite(pred_test, file=paste0(output_path, "/0predictions/preds/", paste0("test.preds_",plant_code,"_LuvHOG.qr2.r.csv")))
  fwrite(pred_train, file=paste0(output_path, "/0predictions/preds/", paste0("train.preds_",plant_code,"_LuvHOG.qr2.r.csv")))
  
}

# LwdHOG.O.U data ----
tree="locwsangle"
nwp_set = c("gfs", "arpege")
depth_set = c(3, 4, 5, 6)
nbin_set = c(3, 6, 9)

for (plant_code in plants){
  for (nwp in nwp_set){
    locwsangle0 = readRDS(file = paste0(output_path, "/", plant_code, "/", "data", "/", nwp, ".locwsangle"))
    ws.qr.1 = readRDS(file = paste0(output_path, "/", plant_code, "/", "data", "/", paste0("data.", nwp,".topws.qr2.1"))) %>% select(-topws)
    locwsangle = merge(ws.qr.1, locwsangle0, by=c("date", "hour", "corrected_production"))
    
    trees.locwsangle = make_trees_qr2(data=nwp, tree="locwsangle", depth=depth_set)
    save_rm_element_names(trees.locwsangle, "tree", output_path=output_path, code=plant_code)
    
    save_rm_wname(paste0("trees.", nwp, ".locwsangle.qr2"), "trees.locwsangle", folder="tree", output_path=output_path, code=plant_code)
    
    #save_tree_plots(data=nwp, tree="locwsangle", output_path=output_path, code=plant_code)
    
    put_leaves_qr2(data=nwp, tree="locwsangle", output_path=output_path, code=plant_code)
    
    save_rm_wname(paste0(nwp, ".locwsangle.qr2"), "locwsangle", output_path=output_path, code=plant_code)
  }
  
}


for (plant_code in plants){
  for (nwp in nwp_set){
    for(d in depth_set){
      for(n in nbin_set){
        make_adj_hog_bins_poly_qr2(data=nwp, tree=tree, depth=d, numbin=n, output_path=output_path, code=plant_code)
      }
    }
  }
  
  for(depth in depth_set){
    for(numbin in nbin_set){
      datanamestocombine = paste0("data.", nwp_set, ".", "tree.", tree, ".d", depth, ".b", numbin, ".qr2.3")
      data.treehog.3 <- merge(readRDS(file = paste0(output_path, "/", plant_code, "/", "data", "/", datanamestocombine[1])), readRDS(file = paste0(output_path, "/", plant_code, "/", "data", "/", datanamestocombine[2])), by=c('date','hour', "corrected_production"), suffixes=c(".gfs", ".arpege"))
      save_rm_wname(paste0("data.", "all", ".", "tree.", tree, ".d", depth, ".b", numbin, ".qr2.3"), "data.treehog.3",output_path=output_path, code=plant_code)
      
    }
  }
}


# LwdHOG.O.U val utilization model glm ----
nwp = "all"
doParallel::registerDoParallel(5)
for (plant_code in plants[1:47]){
  for(depth in depth_set){
    for(numbin in nbin_set){
      dataname = paste0("data.", nwp, ".tree.", tree, ".d", depth, ".b", numbin, ".qr2.3")
      data=readRDS(file=paste0(output_path, "/", plant_code, "/", "data", "/", dataname)) %>% arrange(date, hour)
      colnames(data) = make.names(names(data))
      capacity = max(data$corrected_production)
      data$corrected_production = data$corrected_production / capacity
      train = test_train(data)[2][[1]]
      #registerDoMC(cores = 5)
      set.seed(1)
      mymodel <- cv.glmnet(x=as.matrix(train[,-c(1,2,3)]), y=train$corrected_production, alpha=1, nfolds=5, family=quasibinomial(), parallel=TRUE, trace.it=1)
      name_model = paste0("model.glm.cv.", str_sub(dataname,6,-3), ".r.3")
      saveRDS(mymodel, file = paste0(output_path, "/", plant_code, "/", "models", "/", name_model))
    }
  }
}

# LwdHOG.O.U val results ----
depth_set = c(3, 4, 5, 6)
nbin_set = c(3, 6, 9)
tree="locwsangle"
for (plant_code in plants[1:47]){
  train_outputs = data.frame()
  test_outputs = data.frame()
  for(depth in depth_set){
    for(numbin in nbin_set){
      dataname = paste0("data.", nwp, ".tree.", tree, ".d", depth, ".b", numbin, ".qr2.3")
      data=readRDS(file=paste0(output_path, "/", plant_code, "/", "data", "/", dataname)) %>% arrange(date, hour)
      colnames(data) = make.names(names(data))
      capacity = max(data$corrected_production)
      data$corrected_production = data$corrected_production / capacity
      train = test_train(data)[2][[1]]
      test = test_train(data)[1][[1]]
      pred_test = test[,c(1:3)]
      pred_train = train[,c(1:3)]
      
      name_model = paste0("model.glm.cv.", str_sub(dataname,6,-3), ".r.3")
      mymodel=readRDS(file=paste0(output_path, "/", plant_code, "/", "models", "/", name_model))
      pred_test$pred = predict(mymodel, type="response", as.matrix(test[,-c(1,2,3)]))
      pred_test$pred = ifelse(pred_test$pred > capacity, capacity, pred_test$pred)
      pred_test$pred = ifelse(pred_test$pred < 0, 0, pred_test$pred)
      pred_train$pred = predict(mymodel, type="response", as.matrix(train[,-c(1,2,3)]))
      pred_train$pred = ifelse(pred_train$pred > capacity, capacity, pred_train$pred)
      pred_train$pred = ifelse(pred_train$pred < 0, 0, pred_train$pred)
      
      result_test_ = wmape(pred_test$pred, pred_test$corrected_production)*100
      result_train_ = wmape(pred_train$pred, pred_train$corrected_production)*100
      
      result_train = c(plant_code, depth, numbin, result_train_)
      result_test = c(plant_code, depth, numbin, result_test_)
      train_outputs = rbind(train_outputs, result_train)
      test_outputs = rbind(test_outputs, result_test)
    }
  }
  colnames(train_outputs) = c("Plant", "Depth", "Number of Bins", "Train")
  colnames(test_outputs) = c("Plant", "Depth", "Number of Bins", "Test")
  fwrite(train_outputs, file=paste0(output_path, "/", plant_code, "/", "outputs.",nwp,".tree",tree,".qr2.train.glm.cv", ".csv"))
  fwrite(test_outputs, file=paste0(output_path, "/", plant_code, "/", "outputs.",nwp,".tree",tree,".qr2.test.glm.cv", ".csv"))
}


outputs.sensitivity.train = data.frame()
for (i in 1:47){
  plant_code = plants[i]
  depth_nbin = read.csv(file=paste0(output_path, "/",plant_code, "/", "outputs.all.tree",tree,".qr2.train", ".glm.cv.csv"))
  depth_nbin$complexity = 2^(depth_nbin$Depth) * depth_nbin$Number.of.Bins
  depth_nbin = depth_nbin[with(depth_nbin, order(complexity, Number.of.Bins)), ]  
  depth_nbin$dif = c(diff(depth_nbin$Train, lag=1),0)
  depth_nbin_test = read.csv(file=paste0(output_path, "/",plant_code, "/", "outputs.all.tree",tree,".qr2.test", ".glm.cv.csv")) %>% as.data.table()
  
  d=3
  n=3
  er = depth_nbin$Train[1]
  
  for(c in 1:(length(depth_nbin$dif)-1)){
    if (depth_nbin$dif[c]<0){
      d = depth_nbin$Depth[c+1]
      n = depth_nbin$Number.of.Bins[c+1]
      er = depth_nbin$Train[c+1]
    }else{
      break
    }
  }
  params_ = c(plant_code, d, n, er, depth_nbin_test[Depth==d & Number.of.Bins==n,]$Test)
  outputs.sensitivity.train = rbind(outputs.sensitivity.train, params_)
}
colnames(outputs.sensitivity.train)=c("plant", "Depth", "Number.of.Bins", "Train", "Test")
fwrite(outputs.sensitivity.train, file=paste0(output_path, "/", "outputs.",nwp,".tree",tree,".sensitivity.qr2.r.train.glm.cv", ".csv"))





# LwdHOG.O.U results ----

tree = "locwsangle"
depth_nbin = read.csv(file=paste0(output_path, "/", "outputs.all.tree",tree,".sensitivity.qr2.r.train", ".glm.cv.csv"))

for (i in 1:47){
  plant_code=plants[i]
  dataname = paste0("data.", nwp, ".tree.", tree, ".d", depth_nbin$Depth[i], ".b", depth_nbin$Number.of.Bins[i], ".qr2.3")
  data=readRDS(file=paste0(output_path, "/", plant_code, "/", "data", "/", dataname)) %>% arrange(date, hour)
  colnames(data) = make.names(names(data))
  capacity = max(data$corrected_production)
  data$corrected_production = data$corrected_production / capacity
  train = test_train(data)[2][[1]]
  
  name_model = paste0("model.glm.cv.", str_sub(dataname,6,-3), ".r.3")
  mymodel <- readRDS(file=paste0(output_path, "/", plant_code, "/", "models", "/", name_model))
  
  test = test_train(data)[1][[1]]
  pred_test = test[,c(1:3)]
  pred_test$pred = predict(mymodel, type="response", as.matrix(test[,-c(1,2,3)]))
  
  train = test_train(data)[2][[1]]
  pred_train = train[,c(1:3)]
  pred_train$pred = predict(mymodel, type="response", as.matrix(train[,-c(1,2,3)]))
  
  fwrite(pred_test, file=paste0(output_path, "/0predictions/preds/", paste0("test.preds_",plant_code,"_LwdHOG.qr2.r.csv")))
  fwrite(pred_train, file=paste0(output_path, "/0predictions/preds/", paste0("train.preds_",plant_code,"_LwdHOG.qr2.r.csv")))
  
}






# combined Luv + LHOG RF ----

bin_sensitivity = read.csv(paste0(output_path, "/","outputs.all.LHog.sensitivity.qr2.r.train.cv.csv"))
hogbins = bin_sensitivity$min_bin
for (i in 1:47){
  plant_code = plants[i]
  nb=hogbins[i]
  
  data.all.lochog.1 <- merge(readRDS(file = paste0(output_path, "/", plant_code, "/", "data", "/", paste0("data.gfs.lochog", as.character(nb), ".qr2.1"))), readRDS(file = paste0(output_path, "/", plant_code, "/", "data", "/", paste0("data.arpege.lochog", as.character(nb), ".qr2.1"))), by=c('date','hour', "corrected_production"), suffixes=c(".gfs", ".arpege"))
  data.all.luvlochog.1 <- merge(data.all.lochog.1,   readRDS(file = paste0(output_path, "/", plant_code, "/", "data", "/", paste0("data.all.locuv", ".qr2.1"))) , by=c('date','hour', "corrected_production"), suffixes=c(".lhog", ".luv"))
  
  save_rm_wname(paste0("data.all.luvlochog.qr2.1") , "data.all.luvlochog.1", output_path=output_path, code=plant_code)
  
}

nwp_set = c("all")
names = paste0("data.all.luvlochog.qr2.1")
for (plant_code in plants){
  for (name in names){
    model_caret_ranger_cv(name, output_path=output_path, code=plant_code)
  }
}

nwp_set = c("all")
for (nwp in nwp_set){
  outputs = data.frame()
  tuning = data.frame()
  name = paste0("data.all.luvlochog.qr2.1")
  for (plant_code in plants){
    output = c(plant_code)
    tune = c(plant_code)
    tune_ = data.frame(matrix(ncol = 0, nrow = 1))
    
    data=readRDS(file=paste0(output_path, "/", plant_code, "/", "data", "/", name)) %>% arrange(date, hour)
    colnames(data) = make.names(names(data))
    test = test_train(data)[1][[1]]
    pred_test = test[,c(1:3)]
    
    name_model = paste0("model.caret.rf.cv.", substring(name, 6))
    mymodel=readRDS(file=paste0(output_path, "/", plant_code, "/", "models", "/", name_model))
    
    pred_test$pred = predict(mymodel, test)
    
    output = append(output, modelPerformance_caret_ranger(name, model=mymodel, output_path=output_path, code=plant_code)[1])
    outputs = rbind(outputs, output)
    
    
    tune = append(tune, cbind(as.character(list(unique(mymodel[["results"]][["mtry"]]))), list(mymodel[["bestTune"]][["mtry"]]), list(mymodel[["bestTune"]][["splitrule"]])))
    tuning = rbind(tuning,tune)
    
    fwrite(pred_test, file=paste0(output_path, "/0predictions/preds/", paste0("test.preds_",plant_code,"_Luv.LHOG.rf.csv")))
    
  }
  colnames(outputs) = c("plant", paste0("luvlhog.qr2", ".glm.cv"))
  fwrite(outputs, file=paste0(output_path, "/", paste0("outputs.",nwp,".LuvLHOG.qr2.rf.cv", ".csv")))
  
  tuning[] <- lapply(tuning, as.character)
  colnames(tuning) = c("plant", paste0(rep(c("luvlhog"),3), c(".mtry", ".best_mtry", ".best_splitrule")))
  fwrite(tuning, file=paste0(output_path, "/", paste0("tuning.",nwp,".LuvLHOG.qr2.rf.cv", ".csv")))
}




# combined Lhog + HogLwd glm ----

hogbins_loc = read.csv(paste0(output_path, "/","outputs.all.LHog.sensitivity.qr2.r.train.cv.csv"))$min_bin
depth_nbin_Lwd = read.csv(file=paste0(output_path, "/", "outputs.all.tree","locwsangle",".sensitivity.qr2.r.train", ".glm.cv.csv"))
tree = "locwsangle"
nwp="all"
newname = paste0("data.all.lochoglocwsangle.qr2.3")
for (i in 1:length(plants)){
  plant_code = plants[i]
  
  names = c(paste0("data.",nwp,".lochog", hogbins_loc[i], ".qr2.3"), paste0("data.", nwp, ".tree.", tree, ".d", depth_nbin_Lwd$Depth[i], ".b", depth_nbin_Lwd$Number.of.Bins[i], ".qr2.3"))
  data3 = merge(readRDS(file = paste0(output_path, "/", plant_code, "/", "data", "/", names[1])), readRDS(file = paste0(output_path, "/", plant_code, "/", "data", "/", names[2])), by=c('date','hour', "corrected_production"), suffixes=c(".Lhog", ".HogLwd"))
  
  save_rm_wname(newname, "data3", output_path=output_path, code=plant_code)
}

library(doParallel)
doParallel::registerDoParallel(cores = 5)
for (i in 1:length(plants)){
  plant_code = plants[i]
  dataname = "data.all.lochoglocwsangle.qr2.3"
  data=readRDS(file=paste0(output_path, "/", plant_code, "/", "data", "/", dataname)) %>% arrange(date, hour)
  colnames(data) = make.names(names(data))
  capacity = max(data$corrected_production)
  data$corrected_production = data$corrected_production / capacity
  train = test_train(data)[2][[1]]
  set.seed(1)
  mymodel <- cv.glmnet(x=as.matrix(train[,-c(1,2,3)]), y=train$corrected_production, alpha=1, nfolds=5, family=quasibinomial(), parallel=TRUE)
  name_model = paste0("model.glm.cv.", str_sub(dataname,6,-3), ".r.3")
  saveRDS(mymodel, file = paste0(output_path, "/", plant_code, "/", "models", "/", name_model))
  
  
}

nwp_set = c("all")
for (nwp in nwp_set){
  outputs = data.frame()
  for (i in 1:length(plants)){
    plant_code=plants[i]
    dataname = paste0("data.all.lochoglocwsangle.qr2.3")
    data=readRDS(file=paste0(output_path, "/", plant_code, "/", "data", "/", dataname)) %>% arrange(date, hour)
    colnames(data) = make.names(names(data))
    capacity = max(data$corrected_production)
    data$corrected_production = data$corrected_production / capacity
    test = test_train(data)[1][[1]]
    pred_test = test[,c(1:3)]
    
    name_model = paste0("model.glm.cv.", str_sub(dataname,6,-3), ".r.3")
    mymodel=readRDS(file=paste0(output_path, "/", plant_code, "/", "models", "/", name_model))
    pred_test$pred = predict(mymodel, type="response", as.matrix(test[,-c(1,2,3)]))
    fwrite(pred_test, file=paste0(output_path, "/0predictions/preds/", paste0("test.preds_",plant_code,"_LHOG.LwdHOG.qr2.r.csv")))
    
  }
  
}



# combined Lw + Lhog + HogLwd glm ----

hogbins_loc = read.csv(paste0(output_path, "/","outputs.all.LHog.sensitivity.qr2.r.train.cv.csv"))$min_bin
nwp_set = c("gfs", "arpege")
for (i in 1:47){
  plant_code = plants[i]
  names=c("data.all.locws.qr2.3", "data.all.lochoglocwsangle.qr2.3")
  data3 = merge(readRDS(file = paste0(output_path, "/", plant_code, "/", "data", "/", names[1])), readRDS(file = paste0(output_path, "/", plant_code, "/", "data", "/", names[2])), by=c('date','hour', "corrected_production"), suffixes=c(".Lw", ".comb"))
  save_rm_wname("data.all.lwlochoglocwsangle.qr2.3", "data3", output_path=output_path, code=plant_code)
  
}



library(doParallel)
doParallel::registerDoParallel(cores = 5)
for (i in 1:length(plants)){
  plant_code = plants[i]
  dataname = "data.all.lwlochoglocwsangle.qr2.3"
  data=readRDS(file=paste0(output_path, "/", plant_code, "/", "data", "/", dataname)) %>% arrange(date, hour)
  colnames(data) = make.names(names(data))
  capacity = max(data$corrected_production)
  data$corrected_production = data$corrected_production / capacity
  train = test_train(data)[2][[1]]
  set.seed(1)
  mymodel <- cv.glmnet(x=as.matrix(train[,-c(1,2,3)]), y=train$corrected_production, alpha=1, nfolds=5, family=quasibinomial(), parallel=TRUE)
  name_model = paste0("model.glm.cv.", str_sub(dataname,6,-3), ".r.3")
  saveRDS(mymodel, file = paste0(output_path, "/", plant_code, "/", "models", "/", name_model))
  
  
}

nwp_set = c("all")
for (nwp in nwp_set){
  outputs = data.frame()
  for (i in 1:length(plants)){
    plant_code=plants[i]
    dataname = paste0("data.all.lwlochoglocwsangle.qr2.3")
    data=readRDS(file=paste0(output_path, "/", plant_code, "/", "data", "/", dataname)) %>% arrange(date, hour)
    colnames(data) = make.names(names(data))
    capacity = max(data$corrected_production)
    data$corrected_production = data$corrected_production / capacity
    test = test_train(data)[1][[1]]
    pred_test = test[,c(1:3)]
    
    name_model = paste0("model.glm.cv.", str_sub(dataname,6,-3), ".r.3")
    mymodel=readRDS(file=paste0(output_path, "/", plant_code, "/", "models", "/", name_model))
    pred_test$pred = predict(mymodel, type="response", as.matrix(test[,-c(1,2,3)]))
    fwrite(pred_test, file=paste0(output_path, "/0predictions/preds/", paste0("test.preds_",plant_code,"_Lw.LHOG.LwdHOG.qr2.r.csv")))
    
  }
  
}




# ensemble Luv.RF + LHOG ----
library(quadprog)

nwp="all"
outputs = data.frame()
for (plant_code in plants){
  output = c(plant_code)
  train_LHog = read.csv(file=paste0(output_path, "/0predictions/preds/", "train.preds_", plant_code, "_LHOG.qr2.r.csv")) %>% select(date, hour, corrected_production, pred) %>% dplyr::rename(actual=corrected_production, LHOG=pred)
  test_LHog = read.csv(file=paste0(output_path, "/0predictions/preds/", "test.preds_", plant_code, "_LHOG.qr2.r.csv")) %>% select(date, hour, corrected_production, pred) %>% dplyr::rename(actual=corrected_production, LHOG=pred)
  train_luv = read.csv(paste0(output_path, "/0predictions/preds/", paste0("train.",plant_code,".Luv.csv"))) %>% select(date, hour, corrected_production, pred) %>% dplyr::rename(Luv=pred)
  test_luv = read.csv(paste0(output_path, "/0predictions/preds/", paste0("test.",plant_code,".Luv.csv"))) %>% select(date, hour, corrected_production, pred) %>% dplyr::rename(Luv=pred)
  capacity = max(max(train_luv$corrected_production),max(test_luv$corrected_production))
  train_luv$Luv = train_luv$Luv / capacity
  test_luv$Luv = test_luv$Luv / capacity
  train_luv=train_luv %>% select(-corrected_production)
  test_luv=test_luv %>% select(-corrected_production)
  
  train = merge(train_LHog, train_luv, by=c("date", "hour"))
  test = merge(test_LHog, test_luv, by=c("date", "hour"))
  
  
  #data = rbind(train, test)
  
  X <- as.matrix(train[,4:5])
  y <- as.matrix(train$actual)
  
  d <- t(X) %*% y
  D <- t(X) %*% X
  A <- cbind(rep(1,2),diag(2))#constraint LHS
  b0 <- c(1,numeric(2)) # constraint RHS
  set.seed(1)
  sc <- norm(D,"2")
  soln = solve.QP(D/sc,d/sc,A,b0,meq = 1, factorized=FALSE )
  
  w1 <- soln$solution  # Your model weights
  w1
  
  test$pred = as.matrix(test[,4:5]) %*% w1
  
  pred_test = test %>% select(date, hour, actual, pred)
 
  wmape_ = wmape(pred_test$pred, pred_test$actual)*100
  
  output = append(output, wmape_)
  outputs = rbind(outputs, output)
  fwrite(pred_test, file=paste0(output_path, "/0predictions/preds/", paste0("preds.",plant_code,".E.Luv.LHOG.csv")))
  
  
}
colnames(outputs) = c("plant","E.Luv.LHOG")
fwrite(outputs, file=paste0(output_path, "/", paste0("outputs.E.Luv.LHOG.csv")))



# ensemble Lhog + HogLwd ----

nwp="all"

outputs = data.frame()
for (plant_code in plants){
  output = c(plant_code)
  train_LHog = read.csv(file=paste0(output_path, "/0predictions/preds/", "train.preds_", plant_code, "_LHOG.qr2.r.csv")) %>% select(date, hour, corrected_production, pred) %>% dplyr::rename(actual=corrected_production, LHOG=pred)
  test_LHog = read.csv(file=paste0(output_path, "/0predictions/preds/", "test.preds_", plant_code, "_LHOG.qr2.r.csv")) %>% select(date, hour, corrected_production, pred) %>% dplyr::rename(actual=corrected_production, LHOG=pred)
  train_lwdhog = read.csv(file=paste0(output_path, "/0predictions/preds/", "train.preds_", plant_code, "_LwdHOG.qr2.r.csv")) %>% select(date, hour, corrected_production, pred) %>% dplyr::rename(actual=corrected_production, LwdHOG=pred)
  test_lwdhog = read.csv(file=paste0(output_path, "/0predictions/preds/", "test.preds_", plant_code, "_LwdHOG.qr2.r.csv")) %>% select(date, hour, corrected_production, pred) %>% dplyr::rename(actual=corrected_production, LwdHOG=pred)
  
  train = merge(train_LHog, train_lwdhog, by=c("date", "hour", "actual"))
  test = merge(test_LHog, test_lwdhog, by=c("date", "hour", "actual"))
  
  #data = rbind(train, test)
  
  X <- as.matrix(train[,4:5])
  y <- as.matrix(train$actual)
  
  d <- t(X) %*% y
  D <- t(X) %*% X
  A <- cbind(rep(1,2),diag(2))#constraint LHS
  b0 <- c(1,numeric(2)) # constraint RHS
  set.seed(1)
  sc <- norm(D,"2")
  soln = solve.QP(D/sc,d/sc,A,b0,meq = 1, factorized=FALSE )
  
  #soln <- solve.QP(D,d,A,b0,meq = 1)
  w1 <- soln$solution  # Your model weights
  w1
  
  test$pred = as.matrix(test[,4:5]) %*% w1
  
  pred_test = test %>% select(date, hour, actual, pred)
  
  #capacity = max(train$actual)
  #pred = ifelse(pred > capacity, capacity, pred)
  #wmape = wmape(ifelse(pred < 0, 0, pred), test$actual)*100
  wmape_ = wmape(pred_test$pred, pred_test$actual)*100
  
  output = append(output, wmape_)
  outputs = rbind(outputs, output)
  fwrite(pred_test, file=paste0(output_path, "/0predictions/preds/", paste0("preds.",plant_code,".E.LHOG.LwdHOG.csv")))
  
  
}
colnames(outputs) = c("plant","E.LHOG.LwdHOG")
fwrite(outputs, file=paste0(output_path, "/", paste0("outputs.E.LHOG.LwdHOG.csv")))

# ensemble Lw + Lhog + HogLwd ----

outputs = data.frame()
for (plant_code in plants){
  output = c(plant_code)
  train_lw = read.csv(file=paste0(output_path, "/", plant_code, "/", "preds", "/train.", "all", "Lw.QR.R", ".glm")) %>% select(date, hour, corrected_production, Lw.QR.R) %>% dplyr::rename(actual=corrected_production, Lw=Lw.QR.R)
  test_lw = read.csv(file=paste0(output_path, "/", plant_code, "/", "preds", "/", "all", "Lw.QR.R", ".glm")) %>% select(date, hour, corrected_production, Lw.QR.R) %>% dplyr::rename(actual=corrected_production, Lw=Lw.QR.R)
  train_LHog = read.csv(file=paste0(output_path, "/0predictions/preds/", "train.preds_", plant_code, "_LHOG.qr2.r.csv")) %>% select(date, hour, corrected_production, pred) %>% dplyr::rename(actual=corrected_production, LHOG=pred)
  test_LHog = read.csv(file=paste0(output_path, "/0predictions/preds/", "test.preds_", plant_code, "_LHOG.qr2.r.csv")) %>% select(date, hour, corrected_production, pred) %>% dplyr::rename(actual=corrected_production, LHOG=pred)
  train_lwdhog = read.csv(file=paste0(output_path, "/0predictions/preds/", "train.preds_", plant_code, "_LwdHOG.qr2.r.csv")) %>% select(date, hour, corrected_production, pred) %>% dplyr::rename(actual=corrected_production, LwdHOG=pred)
  test_lwdhog = read.csv(file=paste0(output_path, "/0predictions/preds/", "test.preds_", plant_code, "_LwdHOG.qr2.r.csv")) %>% select(date, hour, corrected_production, pred) %>% dplyr::rename(actual=corrected_production, LwdHOG=pred)
  
  train = merge(train_lw, train_LHog, by=c("date", "hour", "actual"))
  train = merge(train, train_lwdhog, by=c("date", "hour", "actual"))
  test = merge(test_lw, test_LHog, by=c("date", "hour", "actual"))
  test = merge(test, test_lwdhog, by=c("date", "hour", "actual"))
  
  #data = rbind(train, test)
  
  X <- as.matrix(train[,4:6])
  y <- as.matrix(train$actual)
  
  d <- t(X) %*% y
  D <- t(X) %*% X
  A <- cbind(rep(1,3),diag(3))#constraint LHS
  b0 <- c(1,numeric(3)) # constraint RHS
  set.seed(1)
  sc <- norm(D,"2")
  soln = solve.QP(D/sc,d/sc,A,b0,meq = 1, factorized=FALSE )
  
  w1 <- soln$solution  # Your model weights
  w1
  
  test$pred = as.matrix(test[,4:6]) %*% w1
  
  pred_test = test %>% select(date, hour, actual, pred)
  
  wmape_ = wmape(pred_test$pred, pred_test$actual)*100
  
  output = append(output, wmape_)
  outputs = rbind(outputs, output)
  fwrite(pred_test, file=paste0(output_path, "/0predictions/preds/", paste0("preds.",plant_code,".E.Lw.LHOG.LwdHOG.csv")))
  
  
}
colnames(outputs) = c("plant","E.Lw.Lhog.HogLwd")
fwrite(outputs, file=paste0(output_path, "/", paste0("outputs.E.Lw.LHOG.LwdHOG.csv")))



