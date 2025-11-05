# 0.1 prep fncs ----
rm(list = ls(all.names = TRUE))
gc()

save_rm_object <- function(object, folder="data", output_path, code){
  saveRDS(object, file = paste0(output_path, "/", code, "/", folder, "/", deparse(substitute(object))))
  rm(list = as.character(substitute(object)), pos = 1, inherits = TRUE)
}

save_rm_object_overall <- function(object, folder="data", output_path){
  saveRDS(object, file = paste0(output_path, "/", deparse(substitute(object))))
  rm(list = as.character(substitute(object)), pos = 1, inherits = TRUE)
}

save_rm_name <- function(object_name, folder="data", output_path, code){
  saveRDS(get(object_name), file = paste0(output_path, "/", code, "/", folder, "/", object_name))
  rm(list = object_name, pos = 1, inherits = TRUE)
}

save_rm_element_names <- function(object_names, folder="data", output_path, code){
  # folder = "data", "tree", "tree_plot", "plot"
  for(name in object_names){
    saveRDS(get(name), file = paste0(output_path, "/", code, "/", folder, "/", name))
    rm(list = name, pos = 1, inherits = TRUE)
  }
}

save_rm_wname <- function(name, object_name, folder="data", output_path, code){
  saveRDS(get(object_name), file = paste0(output_path, "/", code, "/", folder, "/", name))
  rm(list = object_name, pos = 1, inherits = TRUE)
}


# 0.2 Std HOG fncs ----
# not loc based
bins_circular0 <- function(data, nbin=36){
  # for all directions between -180,180
  tmp=data.table(data)
  small_number = 0.01
  max_ang = 180
  min_ang = -180
  
  brk_points=seq(min_ang-small_number,max_ang+small_number,length=nbin+1)
  tmp=tmp[order(angle)]
  tmp[,bin_id:=cut(angle,brk_points)]
  tmp[,bin_number:=as.numeric(bin_id)]
  
  #tmp[,min_angle_bin:=min(angle),by=list(bin_id)]
  #tmp[,max_angle_bin:=max(angle),by=list(bin_id)]
  tmp[,min_angle_bin:=-180+(360/nbin)*(bin_number-1),by=list(bin_id)]
  tmp[,max_angle_bin:=-180+(360/nbin)*(bin_number),by=list(bin_id)]
  
  tmp[,comp2:=(angle-min_angle_bin)/(max_angle_bin-min_angle_bin)]
  tmp[,comp1:=1-comp2]
  
  tmp[,ws_comp1:=comp1*ws]
  tmp[,ws_comp2:=comp2*ws]
  
  tmp[,neighbor_bin:=bin_number+1]
  tmp[neighbor_bin==(nbin+1),neighbor_bin:=1]
  
  tmp1 = dcast(data.table(tmp),date+hour+corrected_production~bin_number,value.var="ws_comp1",fill=0)
  tmp2 = dcast(data.table(tmp),date+hour+corrected_production~neighbor_bin,value.var="ws_comp2",fill=0)
  tmp_all = rbind(tmp1, tmp2, fill=T)
  tmp_all[is.na(tmp_all)] <- 0
  
  tmp_sum = aggregate(. ~ date+hour+corrected_production, data = tmp_all, FUN=sum)
  
  return(tmp_sum)
}
# loc based
bins_circular00 <- function(data, nbin=36){
  # for all directions between -180,180
  tmp=data.table(data)
  small_number = 0.01
  max_ang = 180
  min_ang = -180
  
  brk_points=seq(min_ang-small_number,max_ang+small_number,length=nbin+1)
  tmp=tmp[order(angle)]
  tmp[,bin_id:=cut(angle,brk_points)]
  tmp[,bin_number:=as.numeric(bin_id)]
  
  #tmp[,min_angle_bin:=min(angle),by=list(bin_id)]
  #tmp[,max_angle_bin:=max(angle),by=list(bin_id)]
  tmp[,min_angle_bin:=-180+(360/nbin)*(bin_number-1),by=list(bin_id)]
  tmp[,max_angle_bin:=-180+(360/nbin)*(bin_number),by=list(bin_id)]
  
  tmp[,comp2:=(angle-min_angle_bin)/(max_angle_bin-min_angle_bin)]
  tmp[,comp1:=1-comp2]
  
  tmp[,ws_comp1:=comp1*ws]
  tmp[,ws_comp2:=comp2*ws]
  
  tmp[,neighbor_bin:=bin_number+1]
  tmp[neighbor_bin==(nbin+1),neighbor_bin:=1]
  
  tmp1 = dcast(data.table(tmp),date+hour+loc+corrected_production~bin_number,value.var="ws_comp1",fill=0)
  tmp2 = dcast(data.table(tmp),date+hour+loc+corrected_production~neighbor_bin,value.var="ws_comp2",fill=0)
  tmp_all = rbind(tmp1, tmp2, fill=T)
  tmp_all[is.na(tmp_all)] <- 0
  
  tmp_sum = aggregate(. ~ date+hour+loc+corrected_production, data = tmp_all, FUN=sum)
  
  return(tmp_sum)
}

# SupervisedHOG fncs ----
bins_circular <- function(data, nbin=5){
  # only for angle=360
  tmp=data.table(data)
  small_number=0.01
  max_ang=180
  min_ang=-180
  
  brk_points=seq(min_ang-small_number,max_ang+small_number,length=nbin+1)
  tmp=tmp[order(angle)]
  tmp[,bin_id:=cut(angle,brk_points)]
  tmp[,bin_number:=as.numeric(bin_id)]
  
  #tmp[,min_angle_bin:=min(angle),by=list(bin_id)]
  #tmp[,max_angle_bin:=max(angle),by=list(bin_id)]
  tmp[,min_angle_bin:=-180+(360/nbin)*(bin_number-1),by=list(bin_id)]
  tmp[,max_angle_bin:=-180+(360/nbin)*(bin_number),by=list(bin_id)]
  
  tmp[,comp2:=(angle-min_angle_bin)/(max_angle_bin-min_angle_bin)]
  tmp[,comp1:=1-comp2]
  
  tmp[,ws_comp1:=comp1*ws]
  tmp[,ws_comp2:=comp2*ws]
  
  tmp[,neighbor_bin:=bin_number+1]
  tmp[neighbor_bin==(nbin+1),neighbor_bin:=1]
  
  tmp1 = dcast(data.table(tmp),date+hour+corrected_production+leaves~bin_number,value.var="ws_comp1",fill=0)
  tmp2 = dcast(data.table(tmp),date+hour+corrected_production+leaves~neighbor_bin,value.var="ws_comp2",fill=0)
  tmp_all = rbind(tmp1, tmp2, fill=T)
  tmp_all[is.na(tmp_all)] <- 0
  return(tmp_all)
}
bins_circular_mean <- function(data){
  # only for angle<60
  tmp=data.table(data)
  tmp_all = tmp %>% select(date, hour, corrected_production, leaves, ws)
  return(tmp_all)
}

bins_noncircular <- function(data, nbin=6, minbound, maxbound, anglerange){
  # if angle != 360
  
  tmpp3=data.table(data)
  
  step = anglerange / (nbin-1)
  small_number=0.01
  brk_points=seq(minbound-small_number,maxbound+small_number,length=nbin)
  tmpp3=as.data.table(tmpp3)
  tmpp3=tmpp3[order(angle)]
  tmpp3[,bin_id:=cut(angle,brk_points)]
  tmpp3[,bin_number:=as.numeric(bin_id)]
  
  tmpp3[,min_angle_bin:=minbound+step*(bin_number-1),by=list(bin_id)]
  tmpp3[,max_angle_bin:=minbound+step*(bin_number),by=list(bin_id)]
  
  tmpp3[,comp2:=(angle-min_angle_bin)/(max_angle_bin-min_angle_bin)]
  tmpp3[,comp1:=1-comp2]
  
  tmpp3[,ws_comp1:=comp1*ws]
  tmpp3[,ws_comp2:=comp2*ws]
  tmpp3[,neighbor_bin:=bin_number+1]
  
  tmpp4 = dcast(data.table(tmpp3),date+hour+corrected_production+leaves~bin_number,value.var="ws_comp1",fill=0)
  tmpp5 = dcast(data.table(tmpp3),date+hour+corrected_production+leaves~neighbor_bin,value.var="ws_comp2",fill=0)
  tmp_all = rbind(tmpp4, tmpp5, fill=T)
  tmp_all[is.na(tmp_all)] <- 0
  return(tmp_all)
}



# SupervisedHOG fncs qr2 ----

make_trees_qr2 = function(data="gfs", tree="locuv", depth=c(3, 4, 5, 6, 7, 8, 9)){
  trees = NULL
  temp_dat_all = get(paste0(tree))
  temp_dat = test_train(temp_dat_all, cut_date=365)[[2]] #train
  for (d in depth){
    dat = paste0("tree.", data, ".", tree, ".depth", as.character(d), ".qr2")
    if (tree=="agguv"){
      assign(dat,rpart(corrected_production ~ mean_u + mean_v, temp_dat, control=rpart.control(maxdepth=d, minsplit=200, minbucket=100, cp=0)), pos=1)
      trees = c(trees,dat)
    }else if(tree=="locuv"){
      assign(dat,rpart(corrected_production ~ ., temp_dat[,-c(1,2)], control=rpart.control(maxdepth=d, minsplit=200, minbucket=100, cp=0)), pos=1)
      trees = c(trees,dat)
    }else if(tree=="aggwsangle"){
      assign(dat,rpart(corrected_production ~ ws + angle, temp_dat, control=rpart.control(maxdepth=d, minsplit=200, minbucket=100, cp=0)), pos=1)
      trees = c(trees,dat) 
    }else if(tree=="locwsangle"){
      assign(dat,rpart(corrected_production ~ ., temp_dat[,-c(1,2)], control=rpart.control(maxdepth=d, minsplit=200, minbucket=100, cp=0)), pos=1)
      trees = c(trees,dat) 
    }
    
  }
  
  return(trees)
}

put_leaves_qr2 = function(data="gfs", tree="locuv", output_path, code=plant_code){
  meanuv.wsangle <- mutate(readRDS(file = paste0(output_path, "/", code, "/", "data", "/", "data.", data, ".agguv.qr2.1")), angle=atan2(mean_v, mean_u)* 180/pi, ws=sqrt(mean_v^2+mean_u^2))
  temp_dat = get(paste0(tree)) #get(paste0("wide.", data, ".", tree))
  trees = readRDS(file = paste0(output_path, "/", code, "/", "tree", "/", "trees.", data, ".", tree, ".qr2"))
  if (tree %in% c("locuv", "locwsangle")){
    for (tree_name in trees){
      temp_tree = readRDS(file = paste0(output_path, "/", code, "/", "tree", "/", tree_name))
      pred = data.frame(rpart.predict.leaves(temp_tree, temp_dat))
      pred = setNames(pred, c('leaves'))
      new_data_name = paste0("data.", data, ".", "tree.", tree, ".d", str_sub(tree_name, start= -5, end=-5), ".qr2")
      temp_data2 = as.data.frame(cbind(temp_dat, pred))
      temp_data2$leaves = as.character(temp_data2$leaves)
      temp_data3 <- merge(temp_data2, get(paste0("meanuv.wsangle")),by=c('date','hour', "corrected_production"))
      assign(new_data_name, temp_data3, pos=1)
      saveRDS(get(new_data_name), file = paste0(output_path, "/", code, "/", "data", "/", new_data_name))
      rm(list = new_data_name, pos = 1, inherits = TRUE)
    }
  } else {
    for (tree_name in trees){
      temp_tree = readRDS(file = paste0(output_path, "/", code, "/", "tree", "/", tree_name))
      pred = data.frame(rpart.predict.leaves(temp_tree, temp_dat))
      pred = setNames(pred, c('leaves'))
      new_data_name = paste0("data.", data, ".", "tree.", tree, ".d", str_sub(tree_name, start= -5, end=-5), ".qr2")
      temp_data2 = as.data.frame(cbind(temp_dat, pred))
      temp_data2$leaves = as.character(temp_data2$leaves)
      assign(new_data_name, temp_data2, pos=1)
      saveRDS(get(new_data_name), file = paste0(output_path, "/", code, "/", "data", "/", new_data_name))
      rm(list = new_data_name, pos = 1, inherits = TRUE)
    }
  }
}

adjusted_hog_data_qr2 <- function(depth=3, numbin=5, data="gfs", tree="agguv", output_path, code=plant_code){
  name = paste0("data.", data, ".tree.", tree, ".d", as.character(depth), ".qr2") # get(name)
  tmp <- readRDS(file = paste0(output_path, "/", code, "/", "data", "/", name))
  #if (tree == "locuv.2"){
  #  removed = 1:((ncol(tmp)-4)/2) + (ncol(tmp)-4)/3 + 3
  #  tmp = tmp[,-removed]
  #}else if (tree == "locuv.3"){
  #  removed = 1:((ncol(tmp)-4)/3) + (ncol(tmp)-4)/3 + 3
  #  tmp = tmp[,-removed]
  #}
  nbin=numbin
  
  nodes <- unique(tmp$leaves)
  
  #ggplot(tmp, aes(x=mean_u, y=mean_v, col=leaves)) + geom_point()
  
  #tmp %>% group_by(leaves) %>% summarise(min = min(angle, na.rm= TRUE), max = max(angle, na.rm= TRUE))
  
  dat = c()
  for (node in nodes){
    small_number=0.01
    tmpp = filter(tmp, leaves == node) %>% data.table()
    
    hpts <- chull(tmpp$mean_u, tmpp$mean_v)
    hpts <- c(hpts)
    pts = tmpp[hpts,] %>% select(mean_u, mean_v, angle) %>% as.data.frame()
    rownames(pts) = c(as.character(hpts))
    diff.matrix = outer(pts$angle, pts$angle, '-')
    colnames(diff.matrix) = c(as.character(hpts))
    diff.matrix = cbind(c(as.character(hpts)), diff.matrix)
    colnames(diff.matrix) = c("index", c(as.character(hpts)))
    a = pivot_longer(data.frame(diff.matrix), !index, names_to = "other", values_to = "anglediff")
    a$anglediff = as.double(a$anglediff)
    index1=unlist(a[which.max(a$anglediff),][1])
    index2=substring(as.character(unlist(a[which.max(a$anglediff),][2])), 2)
    min_angle_info =pts[row.names(pts) == index1,]
    max_angle_info = pts[row.names(pts) == index2,]
    middle_angle_info = arrange(pts, angle)[round(nrow(pts)/2),]
    b = pts[c(as.character(index1), as.character(index2)),]
    
    tmpp2=tmpp
    b2=b
    
    
    if (b$angle[1] >0 && b$angle[2] >0){ # 1
      max_ang=max(b$angle)
      min_ang=min(b$angle)
      brk_points=seq(min_ang-small_number,max_ang+small_number,length=nbin)
      angle.to.be.binned = max_ang-min_ang # max(a$anglediff)
      if(pts$angle >= min_ang && pts$angle <= max_ang){ # 1a
        # since middle_angle_info is on the convex hull
        # acute angle no adjustment
      }else{ # 1b
        # other side
        angle.to.be.binned = 360 - angle.to.be.binned
        tmpp2$angle = ifelse(tmpp$angle < min_ang, tmpp$angle + 360, tmpp$angle)
        min_ang = max(b$angle)
        max_ang = min(b$angle) + 360
        brk_points=seq(min_ang-small_number,max_ang+small_number,length=nbin)
        
      }
    }else if(b$angle[1] <0 && b$angle[2] <0){ #2
      max_ang=max(b$angle)
      min_ang=min(b$angle)
      brk_points=seq(min_ang-small_number,max_ang+small_number,length=nbin)
      angle.to.be.binned = max_ang-min_ang  #max(a$anglediff)
      if(pts$angle >= min_ang && pts$angle <= max_ang){ #2a
        # since middle_angle_info is on the convex hull
        # acute angle no adjustment
      }else{ #2b
        # other side
        angle.to.be.binned = 360 - angle.to.be.binned
        tmpp2$angle = -tmpp$angle
        min_ang = -max(b$angle)
        max_ang = -min(b$angle)
        tmpp2$angle = ifelse(tmpp2$angle < min_ang, tmpp2$angle + 360, tmpp2$angle)
        min_ang = -min(b$angle)
        max_ang = -max(b$angle) + 360
        brk_points=seq(min_ang-small_number,max_ang+small_number,length=nbin)
      }
    }else if(b$angle[1] * b$angle[2] < 0){ # one positive one negative 4 cases, each with 2 sub-cases
      max_ang=max(b$angle)
      min_ang=min(b$angle)
      brk_points=seq(min_ang-small_number,max_ang+small_number,length=nbin)
      angle.to.be.binned = max_ang+abs(min_ang) # max(a$anglediff)
      u1=b$mean_u[1]
      v1=b$mean_v[1]
      u2=b$mean_u[2]
      v2=b$mean_v[2]
      
      if(u1<=0 && u2<=0 && v1*v2<=0){ #3 then 2nd & 3rd quadrant
        if(min(tmpp$angle) <= -175 && max(tmpp$angle) >= 175 && sum(pts$angle[pts$angle>0]<max_ang)<5 && sum(pts$angle[pts$angle<0]>min_ang)<5){ #3a
          # acute angle passes -180
          max_ang = max(pts$angle[pts$angle<0]) + 360
          min_ang = min(pts$angle[pts$angle>0])
          tmpp2$angle = ifelse(tmpp$angle < 0, tmpp$angle + 360, tmpp$angle)
          brk_points=seq(min_ang-small_number,max_ang+small_number,length=nbin)
          angle.to.be.binned = max_ang-min_ang # max(a$anglediff)
        }else if(min(tmpp$angle) <= -175 && max(tmpp$angle) >= 175 && sum(pts$angle[pts$angle>0]<max_ang)>5 && sum(pts$angle[pts$angle<0]>min_ang)>5){# covers 360
          angle.to.be.binned=360
        }else{# 3b look outer pts
          min_ang = min(tmpp$angle)
          max_ang = max(tmpp$angle)
          angle.to.be.binned = max_ang-min_ang
          brk_points=seq(min_ang-small_number,max_ang+small_number,length=nbin)
        }
      }
      if(u1>=0 && u2>=0 && v1*v2<=0){ #4 then 1st & 4th quadrant
        if(sum(pts$angle<= min_ang) >=1 || sum(pts$angle>= max_ang) >=1){ #4a
          #passes -180
          tmpp2$angle = ifelse(tmpp$angle < 0, tmpp$angle + 360, tmpp$angle)
          max_ang = min_ang + 360
          min_ang = max(b$angle)
          brk_points=seq(min_ang-small_number,max_ang+small_number,length=nbin)
          angle.to.be.binned = max_ang-min_ang # max(a$anglediff)
        }# 4b else no adjustment
      }
      if(u1*u2<=0 && v1*v2<=0 && u1*v1<=0){ #5 then 2nd & 4th quadrant
        if(min(tmpp$angle) <= -175 && max(tmpp$angle) >= 175){ #5a
          #passes -180
          tmpp2$angle = ifelse(tmpp$angle < 0, tmpp$angle + 360, tmpp$angle)
          max_ang = min_ang + 360
          min_ang = max(b$angle)
          brk_points=seq(min_ang-small_number,max_ang+small_number,length=nbin)
          angle.to.be.binned = max_ang-min_ang # max(a$anglediff)
        }else{# 5b look outer
          min_ang = min(tmpp$angle)
          max_ang = max(tmpp$angle)
          angle.to.be.binned = max_ang-min_ang
          brk_points=seq(min_ang-small_number,max_ang+small_number,length=nbin)
        }
      }
      if(u1*u2<=0 && v1*v2<=0 && u1*v1>=0){ #6 then 1st & 3rd quadrant
        if(min(tmpp$angle) <= -175 && max(tmpp$angle) >= 175){ #6b
          #passes -180
          tmpp2$angle = ifelse(tmpp$angle < 0, tmpp$angle + 360, tmpp$angle)
          max_ang = min_ang + 360
          min_ang = max(b$angle)
          brk_points=seq(min_ang-small_number,max_ang+small_number,length=nbin)
          angle.to.be.binned = max_ang-min_ang # max(a$anglediff)
        }else{# 6a look outer
          min_ang = min(tmpp$angle)
          max_ang = max(tmpp$angle)
          angle.to.be.binned = max_ang-min_ang
          brk_points=seq(min_ang-small_number,max_ang+small_number,length=nbin)
        }
      }
    }
    
    if (angle.to.be.binned==360){
      tmp_node = bins_circular(tmpp, nbin=numbin)
    } else if(angle.to.be.binned<60){
      tmp_node = bins_circular_mean(tmpp)
    } else {
      tmp_node = bins_noncircular(tmpp2, nbin=numbin, minbound=min_ang, maxbound=max_ang, anglerange=angle.to.be.binned)
    }
    
    tmp_dat <- aggregate(. ~ date+hour+corrected_production+leaves, data = tmp_node, FUN=sum)
    tmp_dat <- setNames(tmp_dat, c(names(tmp_dat[,c(1,2,3,4)]), paste0("node", node ,"_bin", names(tmp_dat[,-c(1,2,3,4)]))))
    
    dat = bind_rows(dat, tmp_dat)
  }
  
  dat[is.na(dat)] <- 0
  
  return(select(dat, -leaves))
}

make_adj_hog_bins_poly_qr2 = function(data="gfs", tree="agguv", depth=3, numbin=5, output_path, code=plant_code){
  dat = adjusted_hog_data_qr2(depth, numbin, data, tree, output_path, code)
  poly1 <- setNames(dat[,-c(1,2,3)]^2, c(paste0(names(dat[,-c(1,2, 3)]),'_2'))  )
  dat2 <- as.data.frame(cbind(dat, poly1))
  poly2 <- setNames(dat[,-c(1,2,3)]^3, c(paste0(names(dat[,-c(1,2, 3)]),'_3'))  )
  dat3 <- as.data.frame(cbind(dat2, poly2))
  name3 = paste0("data.", data, ".", "tree.", tree, ".d", depth, ".b", numbin, ".qr2.3")
  assign(name3, dat3, pos=1)
  saveRDS(get(name3), file = paste0(output_path, "/", code, "/", "data", "/", name3))
  rm(list = name3, pos = 1, inherits = TRUE)
}
save_tree_plots_qr2 = function(data="gfs", tree="locuv", output_path, code=plant_code){
  trees = readRDS(file = paste0(output_path, "/", code, "/", "tree", "/", "trees.", data, ".", tree, ".qr2"))
  for (tree_name in trees){
    png(file=paste0(output_path, "/", code, "/", "tree_plot", "/", tree_name, ".png"), width=720, height=720)
    temp = readRDS(file = paste0(output_path, "/", code, "/", "tree", "/", tree_name))
    fancyRpartPlot(temp)
    dev.off()
  }
}

# 0.3 model fncs general ----

wmape = function(pred, actual){
  return(sum(abs(pred-actual))/sum(actual))
}

test_train <- function(data, cut_date=365){
  max_date = as.Date(max(data$date),format='%Y-%B-%d')
  test_data <- data %>% filter(date>max_date-cut_date) # last 365 days
  train_data <- data %>% filter(date<=max_date-cut_date) # remove last 365 days
  test_data[is.na(test_data)] <- 0
  train_data[is.na(train_data)] <- 0
  return(list(test_data, train_data))
}

wmapeSummary <- function (data,
                          lev = NULL,
                          model = NULL) {
  c(WMAPE=100*sum(abs(data$pred-data$obs))/sum(data$obs))
}



# 0.4 Random Forest Model: ranger from caret ----

model_caret_ranger <- function(dataname, output_path, code=plant_code){
  #registerDoParallel(cores=3)
  
  data=readRDS(file=paste0(output_path, "/", code, "/", "data", "/", dataname)) %>% arrange(date, hour)
  colnames(data) = make.names(names(data))
  
  test = test_train(data)[1][[1]]
  train = test_train(data)[2][[1]]
  
  initial_window_size = floor(ceiling(nrow(train)*0.6)/24) * 24
  
  TimeControl <- trainControl(method = "timeslice",
                              initialWindow = initial_window_size,
                              horizon = 720,
                              skip = 720,
                              fixedWindow = FALSE,
                              summaryFunction = wmapeSummary
                              #verboseIter = TRUE,
                              #allowParallel = TRUE
  )
  
  if (ncol(train)-3==1){ # otherwise /6 results in 0
    m = c(1)
  }else if(ncol(train)-3<=5){ # otherwise /6 results in 0
    m = c(1, 2)
  }else{
    m = c(floor((ncol(train)-3)/6), floor((ncol(train)-3)/3))
  }
  
  rf_grid=expand.grid(mtry = m, # Number of variables to possibly split at each node
                      min.node.size=c(5), # Minimum size of terminal nodes. If large, then smaller trees are grown, so less time
                      splitrule=c("variance", "extratrees") # For regression "variance", "extratrees", "maxstat" or "beta"
  )
  
  set.seed(1)
  model_time = caret::train(x = data.frame(train[,-c(1,2,3)]), 
                            y = train$corrected_production, 
                            method = "ranger", 
                            trControl = TimeControl, 
                            #num.trees = 200, # for simple ranger 200 was used before
                            tuneGrid = rf_grid,
                            metric='WMAPE',
                            maximize=FALSE
  )
  
  name_model = paste0("model.caret.rf.", substring(dataname, 6))
  assign(name_model, model_time, pos=1)
  saveRDS(get(name_model), file = paste0(output_path, "/", code, "/", folder="models", "/", name_model))
  rm(list = name_model, pos = 1, inherits = TRUE)
}

modelPerformance_caret_ranger <- function(dataname, model, output_path, code=plant_code){
  data=readRDS(file=paste0(output_path, "/", code, "/", "data", "/", dataname)) %>% arrange(date, hour)
  colnames(data) = make.names(names(data))
  test = test_train(data)[1][[1]]
  train = test_train(data)[2][[1]]
  
  pred_test = predict(model,test)
  pred_train = predict(model,train)
  capacity = max(train$corrected_production)
  pred_test = ifelse(pred_test > capacity, capacity, pred_test)
  wmape_test = wmape(ifelse(pred_test < 0, 0, pred_test), test$corrected_production)*100
  wmape_train = wmape(as.double(ifelse(pred_train < 0, 0, pred_train)), train$corrected_production)*100
  wmape_test
  wmape_train
  return(c(wmape_test, wmape_train, pred_test))
}

model_caret_ranger_cv <- function(dataname, output_path, code=plant_code){
  #registerDoParallel(cores=3)
  
  data=readRDS(file=paste0(output_path, "/", code, "/", "data", "/", dataname)) %>% arrange(date, hour)
  colnames(data) = make.names(names(data))
  
  test = test_train(data)[1][[1]]
  train = test_train(data)[2][[1]]
  
  initial_window_size = floor(ceiling(nrow(train)*0.6)/24) * 24
  
  TrainControl <- trainControl(method = "cv",
                               number = 5,
                               summaryFunction = wmapeSummary
  )
  
  if (ncol(train)-3==1){ # otherwise /6 results in 0
    m = c(1)
  }else if(ncol(train)-3<=5){ # otherwise /6 results in 0
    m = c(1, 2)
  }else{
    m = c(floor((ncol(train)-3)/6), floor((ncol(train)-3)/3))
  }
  
  rf_grid=expand.grid(mtry = m, # Number of variables to possibly split at each node
                      min.node.size=c(5), # Minimum size of terminal nodes. If large, then smaller trees are grown, so less time
                      splitrule=c("variance", "extratrees") # For regression "variance", "extratrees", "maxstat" or "beta"
  )
  
  set.seed(1)
  model_time = caret::train(x = data.frame(train[,-c(1,2,3)]), 
                            y = train$corrected_production, 
                            method = "ranger", 
                            trControl = TrainControl, 
                            #num.trees = 200, # for simple ranger 200 was used before
                            tuneGrid = rf_grid,
                            metric='WMAPE',
                            maximize=FALSE
  )
  
  name_model = paste0("model.caret.rf.cv.", substring(dataname, 6))
  assign(name_model, model_time, pos=1)
  saveRDS(get(name_model), file = paste0(output_path, "/", code, "/", folder="models", "/", name_model))
  rm(list = name_model, pos = 1, inherits = TRUE)
}


# libraries and initialization ----
require(data.table)
require(gratia)
require(lubridate)
require(dplyr)
require(tidyr)
require(zoo)
library(caret)
library(stringr)
library(dplyr)
library(doParallel)
library(ranger)
library(reshape)
library(ggplot2)
library(viridis)
library(tidyr)
library(stringr)
require(rpart)
require(rpart.plot)
require(rattle)
require(treeClust)
require(tsutils)
library(EnvNJ)
library(quantreg)
library(doParallel)
library(glmnet)
output_path = "~/Desktop/tez2.nosync/Outputs"
data_path = "~/Desktop/tez2.nosync/Data"
plants = c("40W000000000573X")

output_path = "/Volumes/Elements1/TEZ_final/Outputs"


