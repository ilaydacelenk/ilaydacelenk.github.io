wd = "~/Desktop/tez2.nosync"
data_path = "~/Desktop/tez2.nosync/Data"
dir.create("~/Desktop/tez2.nosync/Outputs")
output_path = "~/Desktop/tez2.nosync/Outputs"

plants = list.files(data_path)
dir.create(paste0(output_path,"/", "0predictions"))

for (p in plants){
  dir.create(paste0(output_path,"/", p))
  dir.create(paste0(output_path,"/", p, "/", "data"))
  dir.create(paste0(output_path,"/", p, "/", "models"))
  dir.create(paste0(output_path,"/", p, "/", "plot"))
  dir.create(paste0(output_path,"/", p, "/", "preds"))
  dir.create(paste0(output_path,"/", p, "/", "tree"))
  dir.create(paste0(output_path,"/", p, "/", "tree_plot"))
}






