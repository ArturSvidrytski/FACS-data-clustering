# Data import -------------------------------------------------------------
library(dplyr)
library(plotly)
library(rgl)
library(shiny)
path <- 'D:/1_Work/PROJECTS/Side_projects/RPr_Cell_Survival_Analysis/filtered_data'
filt.exper.paths <- dir(path = path, full.names = T)
files <- dir(path = filt.exper.paths[1], full.names = T)

annex.data <- read.csv(file = files[1])

# First visualization -----------------------------------------------------


library(plotly)
library(rgl)
library(shiny)

tmp <- plot_ly(annex.data
                 , x=~RED.HLog
               , y=~GRN.HLog
               , z=~YEL.HLog
               , size = 1)
htmlwidgets::saveWidget(as_widget(tmp), file = 'graph.html')

plot3d(annex.data$RED.HLog, annex.data$GRN.HLog, annex.data$YEL.HLog, )
api_create(tmp, filename="scatter3d-basic")





plot( annex.data$GRN.HLin,  annex.data$RED.HLin )
plot( annex.data$GRN.HLog,  annex.data$RED.HLog )
plot( annex.data$FSC.HLin,  annex.data$SSC.HLin )
plot( annex.data$FSC.HLog,  annex.data$SSC.HLog )
plot( annex.data$GRN.HLog,  annex.data$YEL.HLog )
plot( annex.data$RED.HLog,  annex.data$YEL.HLog )


# Creation of Tibble for all data -----------------------------------------
create.tib.cl.info <- function(){
tib.cluster.info <- tibble::tibble(   cl_info_id = integer(0)
                                      , exper_name = character(0)
                                      , well_name = character(0)
                                      , well_number = character(0)
                                      , cl_method = character(0)
                                      , cl_metric = character(0)
                                      , cl_amount = integer(0)
                                      , cl_def_approach = character(0)
                                      , dim_amount = integer(0)
                                      , cell_amount = integer(0)
                                      , cl_cell_amount = integer(0)
                                      , cl_group_id = integer(0)
                                      , cl_position = character(0)
                                      , quality_Q1_mean = numeric(0)
                                      , quality_Q2_mean = numeric(0)
                                      , quality_Q3 = numeric(0)
                                      , quality_Q1_med = numeric(0)
                                      , quality_Q2_med = numeric(0)
                                      , cl_radius_mean = numeric(0)
                                      , cl_radius_med = numeric(0)
                                      , cl_variance_mean = numeric(0)
                                      , cl_variance_med = numeric(0)
                                      , cent_GRN_mean = numeric(0)
                                      , cent_RED_mean = numeric(0)
                                      , cent_YEL_mean = numeric(0)
                                      , var_GRN_mean = numeric(0)
                                      , var_RED_mean = numeric(0)
                                      , var_YEL_mean = numeric(0)
                                      , cent_GRN_med = numeric(0)
                                      , cent_RED_med = numeric(0)
                                      , cent_YEL_med = numeric(0)
                                      , var_GRN_med = numeric(0)
                                      , var_RED_med = numeric(0)
                                      , var_YEL_med = numeric(0)
                                      
)
return(tib.cluster.info)
}
# k-means -----------------------------------------------------------------
library(tidyverse)
tib.cluster.info <- create.tib.cl.info()
cl.group_id <- 0

for( exper.num in 1:length(filt.exper.paths) ){
  cur.exper <- filt.exper.paths[exper.num]
  print('-----------------------------------------')
  print( paste0( 'Experiment ', exper.num, ' of ', length(filt.exper.paths)  ) )
  print( cur.exper )
  print('-----------------------------------------')
  files <- dir(path = cur.exper, full.names = T)
  
  for( file.num in 1:length(files) ){
    if( ( ( exper.num == 13 ) & (file.num == 19) ) |
        (( exper.num == 15 ) & (file.num == 13 ) ) ){
      next
    }
    cur.file <- files[file.num]
    annex.data.ini <- read.csv(file = cur.file)
    print( paste0('File No.: ', file.num, ' of ', length(files) )) 
    
    for( dim.amount in c(2,3) ){
      if( dim.amount == 2 ){
        annex.data <- annex.data.ini %>%
          select( GRN.HLog, RED.HLog ) %>%
          filter( GRN.HLog>0 & RED.HLog>0 )
      } else if( dim.amount == 3 ) {
        annex.data <- annex.data.ini %>%
          select( GRN.HLog, RED.HLog, YEL.HLog ) %>%
          filter( GRN.HLog>0 & RED.HLog>0 & YEL.HLog>0 )
      }
      
      for( cur.cl.method in 'k-means' ){
        for( cur.d.method in 'euclidian' ){
          for( cl.amount in c(2,3) ){
            for( cl.def.approach in c('no_ajustment', 'm2') ){
              # avoid pointless data
              if( cl.amount == 3 & cl.def.approach == 'm2' ){
                next
              }
              
              cl.res <- kmeans(x = annex.data
                               , centers = cl.amount
                               , iter.max = 1000
                               , nstart = 25)
              
              labels <- cl.res$cluster
              
              if( cl.def.approach == 'm2' ){
                labels[annex.data$GRN.HLog > 2] <- max(labels)+1
              }
              
              tib.postproc.cl.dat <- retrieve.cl.info(data = annex.data
                                                      , labels = labels)
              
              cur.file.basename <- basename(cur.file)
              
              cl.group_id <- cl.group_id + 1
              
              for( i in 1:nrow(tib.postproc.cl.dat) ){
                cur.line <- tib.postproc.cl.dat[i, ]
                tib.cluster.info[nrow(tib.cluster.info)+1,] <-
                  tibble(  nrow(tib.cluster.info)+1
                           , basename(cur.exper)
                           , substr(  cur.file.basename
                                      , start = nchar(cur.file.basename)-6
                                      , stop = nchar(cur.file.basename)-6)
                           , substr(cur.file.basename
                                    , start = nchar(cur.file.basename)-5
                                    , stop = nchar(cur.file.basename)-4)
                           , cur.cl.method
                           , cur.d.method
                           , cl.amount
                           , cl.def.approach
                           , dim.amount
                           , nrow(annex.data)
                           , cur.line$cl_cell_amount
                           , cl.group_id
                           , NA
                           , cur.line$quality_Q1_mean
                           , cur.line$quality_Q2_mean
                           , cur.line$quality_Q3
                           , cur.line$quality_Q1_med
                           , cur.line$quality_Q2_med
                           , cur.line$cl_radius_mean
                           , cur.line$cl_radius_med
                           , cur.line$cl_variance_mean
                           , cur.line$cl_variance_med
                           , cur.line$cent_GRN_mean
                           , cur.line$cent_RED_mean
                           , cur.line$cent_YEL_mean
                           , cur.line$var_GRN_mean
                           , cur.line$var_RED_mean
                           , cur.line$var_YEL_mean
                           , cur.line$cent_GRN_med
                           , cur.line$cent_RED_med
                           , cur.line$cent_YEL_med
                           , cur.line$var_GRN_med
                           , cur.line$var_RED_med
                           , cur.line$var_YEL_med
                  )
              }
              
              
              # do.call(
              #   paste
              #   , as.list( 
              #     c(as.character(tib.cluster.info[nrow(tib.cluster.info), 1:8])
              #       , sep = '_' ) ) )
              
              
              figname <- paste0( 'Q1mean:'
                                 , round(tib.cluster.info$quality_Q1_mean[nrow(tib.cluster.info)], 7)
                                 , ' / Q1med:'
                                 , round(tib.cluster.info$quality_Q1_med[nrow(tib.cluster.info)], 7)
                                 , ' / Q2mean:'
                                 , round(tib.cluster.info$quality_Q2_mean[nrow(tib.cluster.info)], 7)
                                 , ' / Q2med:'
                                 , round(tib.cluster.info$quality_Q2_med[nrow(tib.cluster.info)], 7)
                                 , ' / Q3:'
                                 , round(tib.cluster.info$quality_Q3[nrow(tib.cluster.info)], 5)
              )
              
              # filename <- paste0( 'figures/'
              #                     , 'clamount'
              #                     , cl.amount
              #                     , '_'
              #                     , cur.cl.method
              #                     , '_'
              #                     , cur.d.method
              #                     , '_'
              #                     , sprintf('%03d', file.num)
              #                     , '.tif')
              
              filename <- 
                paste0( 
                  'figures/'
                  , do.call(
                    paste
                    , as.list( 
                      c(as.character(tib.cluster.info[nrow(tib.cluster.info), c(2,3,4,7,8,5,6)])
                        , sep = '_' ) ) )
                  , paste0( '_', dim.amount, 'D' )
                )
              
              if( dim.amount == 2 ){
                save.figure( annex.data
                             , labels
                             , figname
                             , file = paste0(filename, '.png') )
              } else if( dim.amount == 3 ){
                # my.view.zoom<-par3d()$zoom
                # my.view.userMatrix<-par3d()$userMatrix
                # my.view.windowRect<-par3d()$windowRect
                open3d(zoom = my.view.zoom
                       , userMatrix = my.view.userMatrix
                       , windowRect=my.view.windowRect)
                plot3d(annex.data, col = labels, main = figname  )
                rgl.snapshot( paste0(filename, '.png'), fmt="png", top=TRUE)
                
                rgl.close()
                
                # print projection
                save.figure( annex.data[, c(1,2)]
                             , labels
                             , figname
                             , file = paste0(filename, '_projection.png') )
              }
              
            }
          }
        }
      }
    }
  }
}



# k-medoids ---------------------------------------------------------------
#library('kmed')
library('cluster')
#library('factoextra')

d.methods <- c( 'euclidian' )

#fviz_nbclust(annex.data[ , c(3,5) ], pam, method = "silhouette")





for( file.num in 1:length(files) ){
  #for( file.num in c(14,15,16,20,22,24,25) ){
  print( paste0('File No.: ', file.num, ' of ', length(files) )) 
  annex.data <- read.csv(file = files[file.num])
  
  for( cur.cl.method in c(F) ){
    for( cur.d.method in d.methods ){
      
      cl.res <- pam(x = annex.data[, c(3,5) ]
                    , k = 2
                    , stand = cur.cl.method
                    , pamonce = F )
      
      tiff(filename = paste0( 'figures/'
                              , sprintf('%03d', file.num)
                              , '_'
                              , cur.cl.method
                              , '_'
                              , cur.d.method
                              , '_'
                              , '_'
                              
                              , '.tif'), width = 2048, height = 1024  )
      
      plot( annex.data[, c(3,5) ]
            , col = cl.res$clustering
            , main = paste0( 'FileNo:'
                             , file.num
                             , ' / '
                             , cur.cl.method
                             , ' / '
                             , cur.d.method
                             , ' / ') )
      dev.off()
      
    }
  }
}



# GMMs --------------------------------------------------------------------
library('ClusterR')

cl.res <- GMM( data=annex.data[, c(3,5)]
               , gaussian_comps = 2
               , dist_mode = 'eucl_dist',
)

cl.res <- predict_GMM(data=annex.data[, c(3,5)]
                      , cl.res$centroids
                      , cl.res$covariance_matrices
                      , cl.res$weights)  

table( cl.res$cluster_labels )

plot( annex.data[, c(3,5) ]
      , col = cl.res$cluster_labels + 1)



# Hierarchcal clustering ------------------------------------------------

library(parallel)
#tib.cluster.info <- tib.cluster.info[0,]

cl.methods <- c( 'ward.D2'
                 , 'complete'
                 , 'average'
                 , 'mcquitty'
                 , 'median'
                 , 'centroid' )

# cl.methods <- c( 'ward.D2'
#                  , 'average' )

d.methods <- c( 'euclidian'
                , 'maximum'
                , 'manhattan'
                #, 'canberra'
)

# d.methods <- c( 'canberra' )
stopCluster(cl)
ncores <- detectCores()-1

cl <- makeCluster(ncores)
clusterExport(cl, c('create.tib.cl.info'
                    , 'cl.methods'
                    , 'filt.exper.paths'
                    , 'd.methods'
                    , 'save.figure'
                    , 'retrieve.cl.info'
                    , 'my.view.zoom'
                    , 'my.view.windowRect'
                    , 'my.view.userMatrix'))

clusterEvalQ(cl, library( 'tidyverse' ) )
clusterEvalQ(cl, library( 'dplyr' ) )
clusterEvalQ(cl, library( 'plotly' ) )
clusterEvalQ(cl, library( 'rgl' ) )
clusterEvalQ(cl, library( 'shiny' ) )

unlink('parallel_exper_log.txt')
unlink('parallel_file_log.txt')

lst.cl.info <- parLapply( cl, 1:(length(filt.exper.paths)-0), function(exper.num){
  cur.exper <- filt.exper.paths[exper.num]
  write( paste0( 'Experiment ', exper.num, ' of ', length(filt.exper.paths)  ), file = 'parallel_exper_log.txt', append = T)
  print('-----------------------------------------')
  print( paste0( 'Experiment ', exper.num, ' of ', length(filt.exper.paths)  ) )
  print( cur.exper )
  print('-----------------------------------------')
  files <- dir(path = cur.exper, full.names = T)
  
  tib.cluster.info <- create.tib.cl.info()
  cl.group_id <- 0

  for( file.num in 1:length(files) ){
    if( ( ( exper.num == 13 ) & (file.num == 19) ) |
         (( exper.num == 15 ) & (file.num == 13 ) ) ){
      next
    }
    
    cur.file <- files[file.num]
    annex.data.ini <- read.csv(file = cur.file)
    print( paste0('File No.: ', file.num, ' of ', length(files) )) 
    write( paste0( 'Experiment ', exper.num, ' of ', length(filt.exper.paths), ' | ', 'File No.: ', file.num, ' of ', length(files)  ), file = 'parallel_file_log.txt', append = T)
    
    for( dim.amount in c(2,3) ){
      if( dim.amount == 2 ){
        annex.data <- annex.data.ini %>%
          select( GRN.HLog, RED.HLog ) %>%
          filter( GRN.HLog>0 & RED.HLog>0 )
      } else if( dim.amount == 3 ) {
        annex.data <- annex.data.ini %>%
          select( GRN.HLog, RED.HLog, YEL.HLog ) %>%
          filter( GRN.HLog>0 & RED.HLog>0 & YEL.HLog>0 )
      }
      
      for( cur.cl.method in cl.methods ){
        for( cur.d.method in d.methods ){
          for( cl.amount in c(2,3) ){
            for( cl.def.approach in c('no_ajustment', 'm2') ){
              # avoid pointless data
              if( cl.amount == 3 & cl.def.approach == 'm2' ){
                next
              }
              
              cl.res <- hclust(dist(x = annex.data, method = cur.d.method)
                               , method = cur.cl.method)
              labels <- cutree(cl.res, k = cl.amount)
              cl.group_id <- cl.group_id + 1
              if( cl.def.approach == 'm2' ){
                labels[annex.data$GRN.HLog > 2] <- max(labels)+1
                cl.amount <- 3
              }
              
              tib.postproc.cl.dat <- retrieve.cl.info(data = annex.data
                                                      , labels = labels)
              
              cur.file.basename <- basename(cur.file)
              
              for( i in 1:nrow(tib.postproc.cl.dat) ){
                cur.line <- tib.postproc.cl.dat[i, ]
                tib.cluster.info[nrow(tib.cluster.info)+1,] <-
                  tibble(  nrow(tib.cluster.info)+1
                           , basename(cur.exper)
                           , substr(  cur.file.basename
                                      , start = nchar(cur.file.basename)-6
                                      , stop = nchar(cur.file.basename)-6)
                           , substr(cur.file.basename
                                    , start = nchar(cur.file.basename)-5
                                    , stop = nchar(cur.file.basename)-4)
                           , cur.cl.method
                           , cur.d.method
                           , cl.amount
                           , cl.def.approach
                           , dim.amount
                           , nrow(annex.data)
                           , cur.line$cl_cell_amount
                           , cl.group_id
                           , NA
                           , cur.line$quality_Q1_mean
                           , cur.line$quality_Q2_mean
                           , cur.line$quality_Q3
                           , cur.line$quality_Q1_med
                           , cur.line$quality_Q2_med
                           , cur.line$cl_radius_mean
                           , cur.line$cl_radius_med
                           , cur.line$cl_variance_mean
                           , cur.line$cl_variance_med
                           , cur.line$cent_GRN_mean
                           , cur.line$cent_RED_mean
                           , cur.line$cent_YEL_mean
                           , cur.line$var_GRN_mean
                           , cur.line$var_RED_mean
                           , cur.line$var_YEL_mean
                           , cur.line$cent_GRN_med
                           , cur.line$cent_RED_med
                           , cur.line$cent_YEL_med
                           , cur.line$var_GRN_med
                           , cur.line$var_RED_med
                           , cur.line$var_YEL_med
                  )
              }
              

              # do.call(
              #   paste
              #   , as.list( 
              #     c(as.character(tib.cluster.info[nrow(tib.cluster.info), 1:8])
              #       , sep = '_' ) ) )
              
              
              figname <- paste0( 'Q1mean:'
                , round(tib.cluster.info$quality_Q1_mean[nrow(tib.cluster.info)], 7)
                , ' / Q1med:'
                , round(tib.cluster.info$quality_Q1_med[nrow(tib.cluster.info)], 7)
                , ' / Q2mean:'
                , round(tib.cluster.info$quality_Q2_mean[nrow(tib.cluster.info)], 7)
                , ' / Q2med:'
                , round(tib.cluster.info$quality_Q2_med[nrow(tib.cluster.info)], 7)
                , ' / Q3:'
                , round(tib.cluster.info$quality_Q3[nrow(tib.cluster.info)], 5)
                )

              # filename <- paste0( 'figures/'
              #                     , 'clamount'
              #                     , cl.amount
              #                     , '_'
              #                     , cur.cl.method
              #                     , '_'
              #                     , cur.d.method
              #                     , '_'
              #                     , sprintf('%03d', file.num)
              #                     , '.tif')
              
              filename <- 
                paste0( 
                  'figures/'
                  , do.call(
                    paste
                    , as.list( 
                      c(as.character(tib.cluster.info[nrow(tib.cluster.info), c(2,3,4,7,8,5,6)])
                        , sep = '_' ) ) )
                  , paste0( '_', dim.amount, 'D' )
                )
              
              if( dim.amount == 2 ){
                save.figure( annex.data
                             , labels
                             , figname
                             , file = paste0(filename, '.png') )
              } else if( dim.amount == 3 ){
                # my.view.zoom<-par3d()$zoom
                # my.view.userMatrix<-par3d()$userMatrix
                # my.view.windowRect<-par3d()$windowRect
                open3d(zoom = my.view.zoom
                       , userMatrix = my.view.userMatrix
                       , windowRect=my.view.windowRect)
                plot3d(annex.data, col = labels, main = figname  )
                rgl.snapshot( paste0(filename, '.png'), fmt="png", top=TRUE)
                
                rgl.close()
                
                
                # print projection
                save.figure( annex.data[, c(1,2)]
                             , labels
                             , figname
                             , file = paste0(filename, '_projection.png') )
                
              }
              
            }
          }
        }
      }
    }
  }
  return(tib.cluster.info)
})
stopCluster(cl)


# Spectral clustering: FCD ------------------------------------------------

library('fcd')

cl.res <- spectral.clustering(A = annex.data[, c(3,5)], K = 2 )

# Spectral clustering: kernlab --------------------------------------------

library('kernlab')

system.time({cl.res <- kernlab::specc(as.matrix(annex.data[1:1100, c(3,5)]), centers = 2)})
plot( annex.data$FSC.HLin[1:1100]
      ,  annex.data$SSC.HLin[1:1100]
      , col = cl.res@.Data)


# Test --------------------------------------------------------------------


# 
# 
# cars$time <- cars$dist/cars$speed
# 
# ui <- fluidPage(
#   hr("how do we get the plot inside this app window rather than in a popup?"),
#   
#   rglwidgetOutput("plot",  width = 800, height = 600)
# )
# 
# server <- (function(input, output) {
#   
#   output$plot <- renderRglwidget({
#     rgl.open(useNULL=T)
#     scatter3d(x=cars$speed, y=cars$dist, z=cars$time, surface=FALSE, ellipsoid = TRUE)
#     rglwidget()
#   })
# })   
# shinyApp(ui = ui, server = server)
# 
