library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(png)
library(rjson)

# Save tissue_lowres_image.png file
writePNG(GetImage(hj3, mode='raw'), 'pngsuka.png')

# Save scale factors to json scalefactors_json.json
#Radius(hj3@images[["slice1"]])
json_coordinates <- toJSON(ScaleFactors(hj3@images[["slice1"]]))
write(json_coordinates, file="scalefactors_json.json")

# Save tissue coordinates to tissue_positions_list.csv
write.csv(hj3@images[["slice1"]]@coordinates, 'tissue_positions_list.csv')

Sys.glob("rds_data/*.Rds")
list.files('rds_data')

dir.create('spatial_data')

output_dir <- "spatial_data"
data_dir <- "rds_data"

for (file in list.files(data_dir)) {
  file_path <- paste(data_dir, file, sep='/')
  output_dir_path <- paste(output_dir, tools::file_path_sans_ext(file), sep='/')
  # create dir for object
  dir.create(output_dir_path)
  # read Rds object
  seurat_object <- readRDS(file_path)
  # Save tissue png file
  writePNG(GetImage(seurat_object, mode='raw'), paste(output_dir_path, 'tissue_lowres_image.png', sep='/'))
  # Save scale factors
  json_coordinates <- toJSON(ScaleFactors(seurat_object@images[["slice1"]]))
  write(json_coordinates, file=paste(output_dir_path, "scalefactors_json.json", sep='/'))
  # Save tissue coordinates
  write.csv(seurat_object@images[["slice1"]]@coordinates, paste(output_dir_path, 'tissue_positions_list.csv', sep='/'))
  # Convert Seurat to h5ad
  file_name <- paste(output_dir_path, '/', tools::file_path_sans_ext(file),'.h5Seurat', sep='')
  SaveH5Seurat(seurat_object, filename = file_name)
  Convert(file_name, dest = "h5ad")
}

