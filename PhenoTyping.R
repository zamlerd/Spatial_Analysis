install.packages("devtools")
devtools::install_github("akoyabio/phenoptr")
devtools::install_github("akoyabio/rtree")
devtools::install_github("akoyabio/tiff")
devtools::install_github("akoyabio/phenoptrReports")

usethis::browse_github_pat()

usethis::edit_r_environ()
phenoptrReports::spatial_map_viewer('~/Documents/Thesis/IMPx/IMPx4/Vectra/Dan-Images/
                                    Phenotyping/Consolidated_data.txt','~/Documents/Thesis/
                                    IMPx/IMPx4/Vectra/Dan-Images/Myeloid MSI')

phenoptrReports::spatial_map_viewer('~/Documents/Thesis/IMPx/IMPx4/Vectra/Dan-Images/
                                    Phenotyping/T/Consolidated_data.txt', '~/Documents/
                                    Thesis/IMPx/IMPx4/Vectra/Dan-Images/T MSI')
