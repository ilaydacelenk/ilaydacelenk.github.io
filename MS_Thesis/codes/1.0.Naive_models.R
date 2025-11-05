# Aw data ----
nwp_set = c("gfs", "arpege")
for (plant_code in plants){
  for (nwp in nwp_set){
    wide_ws <- readRDS(file = paste0(output_path, "/", plant_code, "/", "data", "/", "wide_data_", nwp)) %>% group_by(date, hour) %>% select(-level) %>% summarise(mean_u=mean(UGRD), mean_v=mean(VGRD), corrected_production=mean(corrected_production), meanws=mean(ws)) %>% data.table() %>% mutate(ws=sqrt(mean_v^2+mean_u^2))
    # ws: sqrt(mean_u^2 + mean_v^2)
    # meanws: mean(ws) where ws values are available for each location
    data.aggws.1 <- wide_ws %>% select(date, hour, corrected_production, ws)
    data.aggws.2 <- data.aggws.1 %>% mutate(ws2=ws^2)
    data.aggws.3 <- data.aggws.2 %>% mutate(ws3=ws^3)
    
    save_rm_wname(paste0("wide_data_",nwp,"_ws"), paste0("wide_ws"), folder="data", output_path=output_path, code=plant_code)
    for(p in c(1, 3)){
      save_rm_wname(paste0("data.",nwp,".aggws.", p), paste0("data.aggws.", p), folder="data", output_path=output_path, code=plant_code)
    }
    
  }
  
  data.all.aggws.1 <- merge(readRDS(file = paste0(output_path, "/", plant_code, "/", "data", "/", "data.gfs.aggws.1")), readRDS(file = paste0(output_path, "/", plant_code, "/", "data", "/", "data.arpege.aggws.1")), by=c('date','hour', "corrected_production"), suffixes=c(".gfs", ".arpege"))
  data.all.aggws.3 <- merge(readRDS(file = paste0(output_path, "/", plant_code, "/", "data", "/", "data.gfs.aggws.3")), readRDS(file = paste0(output_path, "/", plant_code, "/", "data", "/", "data.arpege.aggws.3")), by=c('date','hour', "corrected_production"), suffixes=c(".gfs", ".arpege"))
  save_rm_name("data.all.aggws.1", output_path=output_path, code=plant_code)
  save_rm_name("data.all.aggws.3", output_path=output_path, code=plant_code)
  
}
# Aw power glm ----
nwp_set = c("gfs", "arpege", "all")
names = paste0("data.", nwp_set, ".aggws.3")
for (plant_code in plants){
  for (dataname in names){
    data=readRDS(file=paste0(output_path, "/", plant_code, "/", "data", "/", dataname)) %>% arrange(date, hour)
    colnames(data) = make.names(names(data))
    train = test_train(data)[2][[1]]
    registerDoMC(cores = 5)
    mymodel <- cv.glmnet(x=as.matrix(train[,-c(1,2,3)]), y=train$corrected_production, nfolds=5, alpha=1, family="gaussian")
    name_model = paste0("model.glm.cv.", substring(dataname, 6))
    saveRDS(mymodel, file = paste0(output_path, "/", plant_code, "/", "models", "/", name_model))
  }
}
for (nwp in nwp_set){
  for (plant_code in plants){
    dataname = paste0("data.", nwp, ".aggws", ".3")
    data=readRDS(file=paste0(output_path, "/", plant_code, "/", "data", "/", dataname)) %>% arrange(date, hour)
    colnames(data) = make.names(names(data))
    capacity = max(data$corrected_production)
    train = test_train(data)[2][[1]]
    test = test_train(data)[1][[1]]
    pred_test = test[,c(1:3)]
    
    name_model = paste0("model.glm.cv.", substring(dataname, 6))
    mymodel=readRDS(file=paste0(output_path, "/", plant_code, "/", "models", "/", name_model))
    pred_test$Aw = predict(mymodel, type="response", as.matrix(test[,-c(1,2,3)]))
    pred_test$Aw = ifelse(pred_test$Aw > capacity, capacity, pred_test$Aw)
    pred_test$Aw = ifelse(pred_test$Aw < 0, 0, pred_test$Aw)
    fwrite(pred_test, file=paste0(output_path, "/", plant_code, "/", "preds", "/", nwp, "Aw", ".glm"))
    
  }
  
}
# Aw utilization model glm ----
nwp_set = c("gfs", "arpege", "all")
names = paste0("data.", nwp_set, ".aggws.3")
for (plant_code in plants){
  for (dataname in names){
    data=readRDS(file=paste0(output_path, "/", plant_code, "/", "data", "/", dataname)) %>% arrange(date, hour)
    colnames(data) = make.names(names(data))
    capacity = max(data$corrected_production)
    data$corrected_production = data$corrected_production / capacity
    train = test_train(data)[2][[1]]
    #registerDoParallel(5)
    registerDoMC(cores = 5)
    mymodel <- cv.glmnet(x=as.matrix(train[,-c(1,2,3)]), y=train$corrected_production, alpha=1, nfolds=5, family=quasibinomial(), parallel=TRUE)
    name_model = paste0("model.glm.cv.", str_sub(dataname,6,-3), ".r.3")
    saveRDS(mymodel, file = paste0(output_path, "/", plant_code, "/", "models", "/", name_model))
  }
}

for (nwp in nwp_set){
  outputs = data.frame()
  for (plant_code in plants){
    dataname = paste0("data.", nwp, ".aggws", ".3")
    data=readRDS(file=paste0(output_path, "/", plant_code, "/", "data", "/", dataname)) %>% arrange(date, hour)
    colnames(data) = make.names(names(data))
    capacity = max(data$corrected_production)
    data$corrected_production = data$corrected_production / capacity
    test = test_train(data)[1][[1]]
    pred_test = test[,c(1:3)]
    
    name_model = paste0("model.glm.cv.", str_sub(dataname,6,-3), ".r.3")
    mymodel=readRDS(file=paste0(output_path, "/", plant_code, "/", "models", "/", name_model))
    pred_test$Aw.R = predict(mymodel, type="response", as.matrix(test[,-c(1,2,3)]))
    fwrite(pred_test, file=paste0(output_path, "/", plant_code, "/", "preds", "/", nwp, "Aw.R", ".glm"))
    
  }
  
}


# top cor loc ----
top_locs = data.frame()
nwp_set = c("gfs", "arpege")
for (plant_code in plants){
  top_loc = c(plant_code)
  for (nwp in nwp_set){
    wide <- readRDS(file = paste0(output_path, "/", plant_code, "/", "data", "/", "wide_data_", nwp))
    wide$loc = paste0(wide$lat, "-", wide$lon)
    cor_wide = wide[,list(correl=cor(ws,corrected_production)),by=list(loc)] %>% arrange(desc(correl))
    top = cor_wide$loc[1]
    top_loc = c(top_loc, top)
    wide_top = wide[loc == top, ]
    save_rm_wname(paste0("wide_data_toploc_", nwp), "wide_top", output_path=output_path, code=plant_code)
  }
  top_locs = rbind(top_locs, top_loc)
}

colnames(top_locs) = c("plant" , nwp_set)
fwrite(top_locs, file=paste0(output_path, "/", paste0("top_cor_locs", ".csv")))

# outlier elimination ----
nwp_set = c("gfs", "arpege")
for (plant_code in plants){
  for (nwp in nwp_set){
    wide_top <- readRDS(file = paste0(output_path, "/", plant_code, "/", "data", "/", "wide_data_toploc_", nwp)) %>% select(date, hour, corrected_production, ws) %>% rename(topws = ws) 
    
    qr01 = rq(corrected_production ~ 0 + I(topws) + I(topws^2) + I(topws^3), tau=0.01, data=wide_top)
    qr99 = rq(corrected_production ~ 0 + I(topws) + I(topws^2) + I(topws^3), tau=0.99, data=wide_top)
    
    data_pred01 = predict(qr01, wide_top)
    data_pred99 = predict(qr99, wide_top)
    topws.qr.1 = wide_top[corrected_production > data_pred01 & corrected_production < data_pred99]
    save_rm_wname(paste0("data.", nwp,".topws.qr2.1"), "topws.qr.1",output_path=output_path, code=plant_code)
    
    
    
  }
}


# Aw.O data ----
nwp_set = c("gfs", "arpege")
for (plant_code in plants){
  for (nwp in nwp_set){
    qr = readRDS(file = paste0(output_path, "/", plant_code, "/", "data", "/", paste0("data.", nwp,".topws.qr2.1")))
    
    #aggws.1
    nonqr = readRDS(file = paste0(output_path, "/", plant_code, "/", "data", "/", paste0("data.",nwp,".aggws", ".1")))
    aggws.qr.1 = merge(qr, nonqr, by=c("date", "hour", "corrected_production")) %>% select(-topws)
    #aggws.3
    nonqr = readRDS(file = paste0(output_path, "/", plant_code, "/", "data", "/", paste0("data.",nwp,".aggws", ".3")))
    aggws.qr2.3 = merge(qr, nonqr, by=c("date", "hour", "corrected_production")) %>% select(-topws)
    
    save_rm_wname(paste0("data.", nwp,".aggws.qr2.1"), "aggws.qr.1",output_path=output_path, code=plant_code)
    save_rm_wname(paste0("data.", nwp,".aggws.qr2.3"), "aggws.qr2.3",output_path=output_path, code=plant_code)
    
  }
  
  data.all.aggws.qr2.1 <- merge(readRDS(file = paste0(output_path, "/", plant_code, "/", "data", "/", "data.gfs.aggws.qr2.1")), readRDS(file = paste0(output_path, "/", plant_code, "/", "data", "/", "data.arpege.aggws.qr2.1")), by=c('date','hour', "corrected_production"), suffixes=c(".gfs", ".arpege"))
  data.all.aggws.qr2.3 <- merge(readRDS(file = paste0(output_path, "/", plant_code, "/", "data", "/", "data.gfs.aggws.qr2.3")), readRDS(file = paste0(output_path, "/", plant_code, "/", "data", "/", "data.arpege.aggws.qr2.3")), by=c('date','hour', "corrected_production"), suffixes=c(".gfs", ".arpege"))
  save_rm_name("data.all.aggws.qr2.1", output_path=output_path, code=plant_code)
  save_rm_name("data.all.aggws.qr2.3", output_path=output_path, code=plant_code)
}

# Aw.O power model glm ----
nwp_set = c("gfs", "arpege", "all")
names = paste0("data.", nwp_set, ".aggws.qr2.3")
for (plant_code in plants){
  for (dataname in names){
    data=readRDS(file=paste0(output_path, "/", plant_code, "/", "data", "/", dataname)) %>% arrange(date, hour)
    colnames(data) = make.names(names(data))
    train = test_train(data)[2][[1]]
    registerDoMC(cores = 5)
    mymodel <- cv.glmnet(x=as.matrix(train[,-c(1,2,3)]), y=train$corrected_production, alpha=1, nfolds=5, family="gaussian")
    name_model = paste0("model.glm.cv.", substring(dataname, 6))
    saveRDS(mymodel, file = paste0(output_path, "/", plant_code, "/", "models", "/", name_model))
  }
}

for (nwp in nwp_set){
  outputs = data.frame()
  for (plant_code in plants){
    dataname = paste0("data.", nwp, ".aggws", ".qr2.3")
    data=readRDS(file=paste0(output_path, "/", plant_code, "/", "data", "/", dataname)) %>% arrange(date, hour)
    colnames(data) = make.names(names(data))
    capacity = max(data$corrected_production)
    train = test_train(data)[2][[1]]
    test = test_train(data)[1][[1]]
    pred_test = test[,c(1:3)]
    
    name_model = paste0("model.glm.cv.", substring(dataname, 6))
    mymodel=readRDS(file=paste0(output_path, "/", plant_code, "/", "models", "/", name_model))
    pred_test$Aw.QR = predict(mymodel, type="response", as.matrix(test[,-c(1,2,3)]))
    pred_test$Aw.QR = ifelse(pred_test$Aw.QR > capacity, capacity, pred_test$Aw.QR)
    pred_test$Aw.QR = ifelse(pred_test$Aw.QR < 0, 0, pred_test$Aw.QR)
    fwrite(pred_test, file=paste0(output_path, "/", plant_code, "/", "preds", "/", nwp, "Aw.QR", ".glm"))
    
  }
  
}

# Aw.O utilization model glm ----
nwp_set = c("gfs", "arpege", "all")
names = paste0("data.", nwp_set, ".aggws.qr2.3")
for (plant_code in plants){
  for (dataname in names){
    data=readRDS(file=paste0(output_path, "/", plant_code, "/", "data", "/", dataname)) %>% arrange(date, hour)
    colnames(data) = make.names(names(data))
    capacity = max(data$corrected_production)
    data$corrected_production = data$corrected_production / capacity
    train = test_train(data)[2][[1]]
    #registerDoParallel(5)
    registerDoMC(cores = 5)
    mymodel <- cv.glmnet(x=as.matrix(train[,-c(1,2,3)]), y=train$corrected_production, alpha=1, nfolds=5, family=quasibinomial(), parallel=TRUE)
    name_model = paste0("model.glm.cv.", str_sub(dataname,6,-3), ".r.3")
    saveRDS(mymodel, file = paste0(output_path, "/", plant_code, "/", "models", "/", name_model))
  }
}

for (nwp in nwp_set){
  outputs = data.frame()
  for (plant_code in plants){
    dataname = paste0("data.", nwp, ".aggws", ".qr2.3")
    data=readRDS(file=paste0(output_path, "/", plant_code, "/", "data", "/", dataname)) %>% arrange(date, hour)
    colnames(data) = make.names(names(data))
    capacity = max(data$corrected_production)
    data$corrected_production = data$corrected_production / capacity
    test = test_train(data)[1][[1]]
    pred_test = test[,c(1:3)]
    
    name_model = paste0("model.glm.cv.", str_sub(dataname,6,-3), ".r.3")
    mymodel=readRDS(file=paste0(output_path, "/", plant_code, "/", "models", "/", name_model))
    pred_test$Aw.QR.R = predict(mymodel, type="response", as.matrix(test[,-c(1,2,3)]))
    fwrite(pred_test, file=paste0(output_path, "/", plant_code, "/", "preds", "/", nwp,"Aw.QR.R", ".glm"))
    
  }
  
}
# Aw.O power model rf ----

nwp_set = c("all")
names = paste0("data.", nwp_set, ".aggws.qr2.1")
for (plant_code in plants){
  for (name in names){
    model_caret_ranger_cv(name, output_path=output_path, code=plant_code)
  }
}

nwp_set = c("all")
for (nwp in nwp_set){
  outputs = data.frame()
  tuning = data.frame()
  name = paste0("data.", nwp, ".aggws.qr2", ".1")
  for (plant_code in plants){
    output = c(plant_code)
    tune = c(plant_code)
    tune_ = data.frame(matrix(ncol = 0, nrow = 1))
    
    name_model = paste0("model.caret.rf.cv.", substring(name, 6))
    mymodel=readRDS(file=paste0(output_path, "/", plant_code, "/", "models", "/", name_model))
    
    output = append(output, modelPerformance_caret_ranger(name, model=mymodel, output_path=output_path, code=plant_code)[1])
    outputs = rbind(outputs, output)
    
    
    tune = append(tune, cbind(as.character(list(unique(mymodel[["results"]][["mtry"]]))), list(mymodel[["bestTune"]][["mtry"]]), list(mymodel[["bestTune"]][["splitrule"]])))
    tuning = rbind(tuning,tune)
    
    
    saveRDS(list(modelPerformance_caret_ranger(name, model=mymodel, output_path=output_path, code=plant_code)[-c(1, 2)]), file=paste0(output_path, "/", plant_code, "/", "preds", "/", paste0("preds.rf.cv.", substring(name, 6))))
  }
  colnames(outputs) = c("plant", paste0("aggws.qr2", ".glm.cv"))
  fwrite(outputs, file=paste0(output_path, "/", paste0("outputs.",nwp,".aggws.qr2.rf.cv", ".csv")))
  
  tuning[] <- lapply(tuning, as.character)
  colnames(tuning) = c("plant", paste0(rep(c("aggws"),3), c(".mtry", ".best_mtry", ".best_splitrule")))
  fwrite(tuning, file=paste0(output_path, "/", paste0("tuning.",nwp,".aggws.qr2.rf.cv", ".csv")))
}

# Lw.O data ----
nwp_set = c("gfs", "arpege")
for (plant_code in plants){
  for (nwp in nwp_set){
    qr = readRDS(file = paste0(output_path, "/", plant_code, "/", "data", "/", paste0("data.", nwp,".topws.qr2.1")))
    
    wide.locws <- readRDS(file = paste0(output_path, "/", plant_code, "/", "data", "/", "wide_data_", nwp)) %>% select(-level)
    wide.locws$loc = paste0(wide.locws$lat, "-", wide.locws$lon)
    wide.locws = wide.locws %>% select(-lat, -lon, -UGRD, -VGRD)
    data.locws.1 = dcast(data.table(wide.locws), date+hour+corrected_production~loc,value.var="ws") # fill=0 yapmadık
    
    poly1 <- setNames(data.locws.1[,-c(1,2, 3)]^2, c(paste0(names(data.locws.1[,-c(1,2, 3)]),'_2'))  )
    data.locws.2 <- as.data.frame(cbind(data.locws.1, poly1))
    poly2 <- setNames(data.locws.1[,-c(1,2, 3)]^3, c(paste0(names(data.locws.1[,-c(1,2, 3)]),'_3'))  )
    data.locws.3 <- as.data.frame(cbind(data.locws.2, poly2))
    
    locws.qr2.1 = merge(qr, data.locws.2, by=c("date", "hour", "corrected_production")) %>% select(-topws)
    locws.qr2.3 = merge(qr, data.locws.3, by=c("date", "hour", "corrected_production")) %>% select(-topws)
    
    save_rm_wname(paste0("data.", nwp,".locws.qr2.1"), "locws.qr2.1",output_path=output_path, code=plant_code)
    save_rm_wname(paste0("data.", nwp,".locws.qr2.3"), "locws.qr2.3",output_path=output_path, code=plant_code)
    
    
  }
  
  data.all.locws.qr2.1 <- merge(readRDS(file = paste0(output_path, "/", plant_code, "/", "data", "/", "data.gfs.locws.qr2.1")), readRDS(file = paste0(output_path, "/", plant_code, "/", "data", "/", "data.arpege.locws.qr2.1")), by=c('date','hour', "corrected_production"), suffixes=c(".gfs", ".arpege"))
  save_rm_name("data.all.locws.qr2.1", output_path=output_path, code=plant_code)
  
  data.all.locws.qr2.3 <- merge(readRDS(file = paste0(output_path, "/", plant_code, "/", "data", "/", "data.gfs.locws.qr2.3")), readRDS(file = paste0(output_path, "/", plant_code, "/", "data", "/", "data.arpege.locws.qr2.3")), by=c('date','hour', "corrected_production"), suffixes=c(".gfs", ".arpege"))
  save_rm_name("data.all.locws.qr2.3", output_path=output_path, code=plant_code)
  
  
}

# Lw.O utilization model glm ----
nwp_set = c("all")
names = paste0("data.", nwp_set, ".locws.qr2.3")
for (plant_code in plants){
  for (dataname in names){
    data=readRDS(file=paste0(output_path, "/", plant_code, "/", "data", "/", dataname)) %>% arrange(date, hour)
    colnames(data) = make.names(names(data))
    capacity = max(data$corrected_production)
    data$corrected_production = data$corrected_production / capacity
    train = test_train(data)[2][[1]]
    #registerDoParallel(5)
    registerDoMC(cores = 5)
    mymodel <- cv.glmnet(x=as.matrix(train[,-c(1,2,3)]), y=train$corrected_production, alpha=1, nfolds=5, family=quasibinomial(), parallel=TRUE)
    name_model = paste0("model.glm.cv.", str_sub(dataname,6,-3), ".r.3")
    saveRDS(mymodel, file = paste0(output_path, "/", plant_code, "/", "models", "/", name_model))
  }
}

for (nwp in nwp_set){
  outputs = data.frame()
  for (plant_code in plants){
    dataname = paste0("data.", nwp, ".locws", ".qr2.3")
    data=readRDS(file=paste0(output_path, "/", plant_code, "/", "data", "/", dataname)) %>% arrange(date, hour)
    colnames(data) = make.names(names(data))
    capacity = max(data$corrected_production)
    data$corrected_production = data$corrected_production / capacity
    test = test_train(data)[1][[1]]
    pred_test = test[,c(1:3)]
    
    name_model = paste0("model.glm.cv.", str_sub(dataname,6,-3), ".r.3")
    mymodel=readRDS(file=paste0(output_path, "/", plant_code, "/", "models", "/", name_model))
    pred_test$Lw.QR.R = predict(mymodel, type="response", as.matrix(test[,-c(1,2,3)]))
    fwrite(pred_test, file=paste0(output_path, "/", plant_code, "/", "preds", "/", nwp,"Lw.QR.R", ".glm"))
    
  }
  
}


# Lw.O power model rf ----
nwp_set = c("all")
names = paste0("data.", nwp_set, ".locws.qr2.1")
for (plant_code in plants){
  for (name in names){
    model_caret_ranger_cv(name, output_path=output_path, code=plant_code)
  }
}

nwp_set = c("all")
for (nwp in nwp_set){
  outputs = data.frame()
  tuning = data.frame()
  name = paste0("data.", nwp, ".locws.qr2", ".1")
  for (plant_code in plants){
    output = c(plant_code)
    tune = c(plant_code)
    tune_ = data.frame(matrix(ncol = 0, nrow = 1))
    
    name_model = paste0("model.caret.rf.cv.", substring(name, 6))
    mymodel=readRDS(file=paste0(output_path, "/", plant_code, "/", "models", "/", name_model))
    
    output = append(output, modelPerformance_caret_ranger(name, model=mymodel, output_path=output_path, code=plant_code)[1])
    outputs = rbind(outputs, output)
    
    
    tune = append(tune, cbind(as.character(list(unique(mymodel[["results"]][["mtry"]]))), list(mymodel[["bestTune"]][["mtry"]]), list(mymodel[["bestTune"]][["splitrule"]])))
    tuning = rbind(tuning,tune)
    
    
    saveRDS(list(modelPerformance_caret_ranger(name, model=mymodel, output_path=output_path, code=plant_code)[-c(1, 2)]), file=paste0(output_path, "/", plant_code, "/", "preds", "/", paste0("preds.rf.cv.", substring(name, 6))))
  }
  colnames(outputs) = c("plant", paste0("locws.qr2", ".glm.cv"))
  fwrite(outputs, file=paste0(output_path, "/", paste0("outputs.",nwp,".locws.qr2.rf.cv", ".csv")))
  
  tuning[] <- lapply(tuning, as.character)
  colnames(tuning) = c("plant", paste0(rep(c("locws"),3), c(".mtry", ".best_mtry", ".best_splitrule")))
  fwrite(tuning, file=paste0(output_path, "/", paste0("tuning.",nwp,".locws.qr2.rf.cv", ".csv")))
}


# Auv.O data ----
nwp_set = c("gfs", "arpege")
for (plant_code in plants){
  for (nwp in nwp_set){
    qr = readRDS(file = paste0(output_path, "/", plant_code, "/", "data", "/", paste0("data.", nwp,".topws.qr2.1")))
    nonqr <- readRDS(file = paste0(output_path, "/", plant_code, "/", "data", "/", "wide_data_",nwp,"_ws")) %>% select(date, hour, corrected_production, mean_u, mean_v)
    
    #agguv.1
    agguv.qr.1 = merge(qr, nonqr, by=c("date", "hour", "corrected_production")) %>% select(-topws)
    
    save_rm_wname(paste0("data.", nwp,".agguv.qr2.1"), "agguv.qr.1",output_path=output_path, code=plant_code)
    
  }
  
  data.all.agguv.qr2.1 <- merge(readRDS(file = paste0(output_path, "/", plant_code, "/", "data", "/", "data.gfs.agguv.qr2.1")), readRDS(file = paste0(output_path, "/", plant_code, "/", "data", "/", "data.arpege.agguv.qr2.1")), by=c('date','hour', "corrected_production"), suffixes=c(".gfs", ".arpege"))
  save_rm_name("data.all.agguv.qr2.1", output_path=output_path, code=plant_code)
}


# Auv.O power model rf ----

nwp_set = c("all")
names = paste0("data.", nwp_set, ".agguv.qr2.1")
for (plant_code in plants){
  for (name in names){
    model_caret_ranger_cv(name, output_path=output_path, code=plant_code)
  }
}

nwp_set = c("all")
for (nwp in nwp_set){
  outputs = data.frame()
  tuning = data.frame()
  name = paste0("data.", nwp, ".agguv.qr2", ".1")
  for (plant_code in plants){
    output = c(plant_code)
    tune = c(plant_code)
    tune_ = data.frame(matrix(ncol = 0, nrow = 1))
    
    name_model = paste0("model.caret.rf.cv.", substring(name, 6))
    mymodel=readRDS(file=paste0(output_path, "/", plant_code, "/", "models", "/", name_model))
    
    output = append(output, modelPerformance_caret_ranger(name, model=mymodel, output_path=output_path, code=plant_code)[1])
    outputs = rbind(outputs, output)
    
    
    tune = append(tune, cbind(as.character(list(unique(mymodel[["results"]][["mtry"]]))), list(mymodel[["bestTune"]][["mtry"]]), list(mymodel[["bestTune"]][["splitrule"]])))
    tuning = rbind(tuning,tune)
    
    
    saveRDS(list(modelPerformance_caret_ranger(name, model=mymodel, output_path=output_path, code=plant_code)[-c(1, 2)]), file=paste0(output_path, "/", plant_code, "/", "preds", "/", paste0("preds.rf.cv.", substring(name, 6))))
  }
  colnames(outputs) = c("plant", paste0("agguv.qr2", ".glm.cv"))
  fwrite(outputs, file=paste0(output_path, "/", paste0("outputs.",nwp,".agguv.qr2.rf.cv", ".csv")))
  
  tuning[] <- lapply(tuning, as.character)
  colnames(tuning) = c("plant", paste0(rep(c("agguv"),3), c(".mtry", ".best_mtry", ".best_splitrule")))
  fwrite(tuning, file=paste0(output_path, "/", paste0("tuning.",nwp,".agguv.qr2.rf.cv", ".csv")))
}


# Luv.O data ----
nwp_set = c("gfs", "arpege")
for (plant_code in plants){
  for (nwp in nwp_set){
    qr = readRDS(file = paste0(output_path, "/", plant_code, "/", "data", "/", paste0("data.", nwp,".topws.qr2.1")))
    
    wide.locuv <- readRDS(file = paste0(output_path, "/", plant_code, "/", "data", "/", "wide_data_", nwp))
    wide.locuv$loc = paste0(wide.locuv$lat, "-", wide.locuv$lon)
    wide.locuv = wide.locuv %>% select(-lat, -lon, -ws, -level)
    nonqr = dcast(data.table(wide.locuv), date+hour+corrected_production~loc,value.var=c("UGRD", "VGRD")) # fill=0 yapmadık
    
    #agguv.1
    locuv.qr.1 = merge(qr, nonqr, by=c("date", "hour", "corrected_production")) %>% select(-topws)
    
    save_rm_wname(paste0("data.", nwp,".locuv.qr2.1"), "locuv.qr.1",output_path=output_path, code=plant_code)
    
  }
  
  data.all.locuv.qr2.1 <- merge(readRDS(file = paste0(output_path, "/", plant_code, "/", "data", "/", "data.gfs.locuv.qr2.1")), readRDS(file = paste0(output_path, "/", plant_code, "/", "data", "/", "data.arpege.locuv.qr2.1")), by=c('date','hour', "corrected_production"), suffixes=c(".gfs", ".arpege"))
  save_rm_name("data.all.locuv.qr2.1", output_path=output_path, code=plant_code)
}

# Luv.O power model rf ----
nwp_set = c("all")
names = paste0("data.", nwp_set, ".locuv.qr2.1")
for (plant_code in plants){
  for (name in names){
    model_caret_ranger_cv(name, output_path=output_path, code=plant_code)
  }
}

nwp_set = c("all")
for (nwp in nwp_set){
  outputs = data.frame()
  tuning = data.frame()
  name = paste0("data.", nwp, ".locuv.qr2", ".1")
  for (plant_code in plants){
    output = c(plant_code)
    tune = c(plant_code)
    tune_ = data.frame(matrix(ncol = 0, nrow = 1))
    
    name_model = paste0("model.caret.rf.cv.", substring(name, 6))
    mymodel=readRDS(file=paste0(output_path, "/", plant_code, "/", "models", "/", name_model))
    
    output = append(output, modelPerformance_caret_ranger(name, model=mymodel, output_path=output_path, code=plant_code)[1])
    outputs = rbind(outputs, output)
    
    
    tune = append(tune, cbind(as.character(list(unique(mymodel[["results"]][["mtry"]]))), list(mymodel[["bestTune"]][["mtry"]]), list(mymodel[["bestTune"]][["splitrule"]])))
    tuning = rbind(tuning,tune)
    
    
    saveRDS(list(modelPerformance_caret_ranger(name, model=mymodel, output_path=output_path, code=plant_code)[-c(1, 2)]), file=paste0(output_path, "/", plant_code, "/", "preds", "/", paste0("preds.rf.cv.", substring(name, 6))))
  }
  colnames(outputs) = c("plant", paste0("locuv.qr2", ".glm.cv"))
  fwrite(outputs, file=paste0(output_path, "/", paste0("outputs.",nwp,".locuv.qr2.rf.cv", ".csv")))
  
  tuning[] <- lapply(tuning, as.character)
  colnames(tuning) = c("plant", paste0(rep(c("locuv"),3), c(".mtry", ".best_mtry", ".best_splitrule")))
  fwrite(tuning, file=paste0(output_path, "/", paste0("tuning.",nwp,".locuv.qr2.rf.cv", ".csv")))
}

# Lwd.O data ----
nwp_set = c("gfs", "arpege")
for (plant_code in plants){
  for (nwp in nwp_set){
    qr = readRDS(file = paste0(output_path, "/", plant_code, "/", "data", "/", paste0("data.", nwp,".topws.qr2.1")))
    
    wide.locwsangle <- readRDS(file = paste0(output_path, "/", plant_code, "/", "data", "/", "wide_data_", nwp)) %>% mutate(angle=atan2(VGRD, UGRD)* 180/pi)
    wide.locwsangle$loc = paste0(wide.locwsangle$lat, "-", wide.locwsangle$lon)
    wide.locwsangle = wide.locwsangle %>% select(-lat, -lon, -level)
    nonqr = dcast(data.table(wide.locwsangle), date+hour+corrected_production~loc,value.var=c("ws", "angle")) # do not fill=0
    
    #locwd.1
    locwsangle.qr2.1 = merge(qr, nonqr, by=c("date", "hour", "corrected_production")) %>% select(-topws)
    
    save_rm_wname(paste0("data.",nwp,".locwsangle.qr2.1"), "locwsangle.qr2.1", output_path=output_path, code=plant_code)
    
  }
  
  data.all.locwsangle.qr2.1 <- merge(readRDS(file = paste0(output_path, "/", plant_code, "/", "data", "/", "data.gfs.locwsangle.qr2.1")), readRDS(file = paste0(output_path, "/", plant_code, "/", "data", "/", "data.arpege.locwsangle.qr2.1")), by=c('date','hour', "corrected_production"), suffixes=c(".gfs", ".arpege"))
  save_rm_name("data.all.locwsangle.qr2.1", output_path=output_path, code=plant_code)
}





# Lwd.O power model rf ----


nwp_set = c("all")
names = paste0("data.", nwp_set, ".locwsangle.qr2.1")
for (plant_code in plants){
  for (name in names){
    model_caret_ranger_cv(name, output_path=output_path, code=plant_code)
  }
}

nwp_set = c("all")
for (nwp in nwp_set){
  outputs = data.frame()
  tuning = data.frame()
  name = paste0("data.", nwp, ".locwsangle.qr2", ".1")
  for (plant_code in plants){
    output = c(plant_code)
    tune = c(plant_code)
    tune_ = data.frame(matrix(ncol = 0, nrow = 1))
    
    name_model = paste0("model.caret.rf.cv.", substring(name, 6))
    mymodel=readRDS(file=paste0(output_path, "/", plant_code, "/", "models", "/", name_model))
    
    output = append(output, modelPerformance_caret_ranger(name, model=mymodel, output_path=output_path, code=plant_code)[1])
    outputs = rbind(outputs, output)
    
    
    tune = append(tune, cbind(as.character(list(unique(mymodel[["results"]][["mtry"]]))), list(mymodel[["bestTune"]][["mtry"]]), list(mymodel[["bestTune"]][["splitrule"]])))
    tuning = rbind(tuning,tune)
    
    
    saveRDS(list(modelPerformance_caret_ranger(name, model=mymodel, output_path=output_path, code=plant_code)[-c(1, 2)]), file=paste0(output_path, "/", plant_code, "/", "preds", "/", paste0("preds.rf.cv.", substring(name, 6))))
  }
  colnames(outputs) = c("plant", paste0("locwsangle.qr2", ".glm.cv"))
  fwrite(outputs, file=paste0(output_path, "/", paste0("outputs.",nwp,".locwsangle.qr2.rf.cv", ".csv")))
  
  tuning[] <- lapply(tuning, as.character)
  colnames(tuning) = c("plant", paste0(rep(c("locwsangle"),3), c(".mtry", ".best_mtry", ".best_splitrule")))
  fwrite(tuning, file=paste0(output_path, "/", paste0("tuning.",nwp,".locwsangle.qr2.rf.cv", ".csv")))
}







