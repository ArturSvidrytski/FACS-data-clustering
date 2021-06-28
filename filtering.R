# Data import -------------------------------------------------------------

path <- 'D:/1_Work/PROJECTS/Side_projects/FACS/Raw_Data'
experiment.paths <- dir( path, full.names = T)
files <- dir(path = path, full.names = T )
annex.data <- read.csv(file = files[10])
annex.data <- read.csv(file = files[15])

# First visualization -----------------------------------------------------



plot( annex.data$GRN.HLin,  annex.data$RED.HLin )
plot( annex.data$GRN.HLog,  annex.data$RED.HLog )
plot( annex.data$FSC.HLin,  annex.data$SSC.HLin )



# Hierarchcal clustering ------------------------------------------------


cl.methods <- c( #'ward.D2'
                  'average' )

for( i in 1:length(experiment.paths) ){
  
  cur.exper.path <- experiment.paths[i]
  
  files <- dir(path = cur.exper.path, full.names = T)
  
  cur.fig.path <- paste(sep = '/'
                        , 'figures'
                        , paste0( sprintf('%03d', i)
                                  , '_'
                                  , basename(cur.exper.path) ) )
  
  if( dir.exists(cur.fig.path) ){
    unlink( cur.fig.path, recursive = T )
  }
  
  dir.create(cur.fig.path)
  
  for( file.num in 1:length(files) ){
    #for( file.num in c(14,15,16,20,22,24,25) ){
    print( paste0('File No.: ', file.num, ' of ', length(files) )) 
    annex.data <- read.csv(file = files[file.num])
    
    for( cur.cl.method in cl.methods ){
      cl.res <- hclust(dist(x = annex.data[, c(1,2)], method = 'canberra')
                       , method = cur.cl.method)
      
      tiff(filename = paste0( cur.fig.path
                              , '/'
                              , sprintf('%03d', file.num)
                              , '_'
                              , cur.cl.method
                              , '_'
                              , 'canberra'
                              , '_'
                              , '_'
                              
                              , '.tif'), width = 2048, height = 1024  )
      plot( annex.data$FSC.HLin
            ,  annex.data$SSC.HLin
            , col = cutree(cl.res, k = 2) + 1
            , main = paste0( 'FileNo:'
                             , file.num
                             , ' / '
                             , cur.cl.method
                             , ' / '
                             , 'canberra'
                             , ' / ') )
      dev.off()
    }
  }
}
  


# Filtering ---------------------------------------------------------------
library(tidyverse)

filtered.data.path <- 'filtered_data'
cl.methods <- c( #'ward.D2'
  'average' )

for( i in 1:length(experiment.paths) ){
  cur.exper.path <- experiment.paths[i]
  files <- dir(path = cur.exper.path, full.names = T)
  cur.filt.data.path <- paste(sep='/'
                              , filtered.data.path
                              , paste0( sprintf('%03d', i)
                                        , '_'
                                        , basename(cur.exper.path) ))
  if( dir.exists(cur.filt.data.path) ){
    unlink( cur.filt.data.path, recursive = T )
  }
  
  dir.create(cur.filt.data.path)
  
  for( file.num in 1:length(files) ){
    #for( file.num in c(14,15,16,20,22,24,25) ){
    print( paste0('File No.: ', file.num, ' of ', length(files) )) 
    cur.file <- files[file.num]
    annex.data <- read.csv(file = cur.file)
    
    for( cur.cl.method in cl.methods ){
      cl.res <- hclust(dist(x = annex.data[, c(1,2)], method = 'canberra')
                       , method = cur.cl.method)
      
      annex.data$labels <- cutree(cl.res, k = 2)
      
      to.filter <- as.integer(
        annex.data %>%
          group_by(labels) %>%
          summarise( avgFSC = mean(FSC.HLin)
                     , avgSSC = mean(SSC.HLin) ) %>%
          mutate( sum = avgFSC+avgSSC) %>%
          top_n(n=-1, wt = sum) %>%
          select(labels) )
      
      annex.data.filt <- annex.data %>%
        filter( labels != to.filter ) %>%
        select(-labels)
      
      write.csv(x = annex.data.filt
                , file = paste(sep='/'
                               , cur.filt.data.path
                               , basename(cur.file))
                , row.names = F)
      

    }
  }
}




# Test --------------------------------------------------------------------



data(iris)
iris <- as.matrix(iris[,1:4])

## find suitable eps parameter using a k-NN plot for k = dim + 1
## Look for the knee!
kNNdistplot(iris, k = 5)
abline(h=.5, col = "red", lty=2)