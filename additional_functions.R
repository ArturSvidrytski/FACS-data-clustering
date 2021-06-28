save.figure <- function(data, labels, figname, file, wdth = 2048, hght = 1024  ){
  png(filename = file, width = wdth, height = hght  )
  plot( data
  , pch = 20
  , col = labels
  , main = figname )
dev.off()
}

retrieve.cl.info <- function( data, labels ){
  
  library(tibble)
  
  cl.lbls <- unique(labels)
  cl.amount <- length(cl.lbls)
  
  tib.cur.cl.inf <- tibble(
      cl_cell_amount = integer(cl.amount)
    , quality_Q1_mean = numeric(cl.amount)
    , quality_Q2_mean = numeric(cl.amount)
    , quality_Q3 = numeric(cl.amount)
    , quality_Q1_med = numeric(cl.amount)
    , quality_Q2_med = numeric(cl.amount)
    , cl_radius_mean = numeric(cl.amount)
    , cl_radius_med = numeric(cl.amount)
    , cl_variance_mean = numeric(cl.amount)
    , cl_variance_med = numeric(cl.amount)
    , cent_GRN_mean = numeric(cl.amount)
    , cent_RED_mean = numeric(cl.amount)
    , cent_YEL_mean = numeric(cl.amount)
    , var_GRN_mean = numeric(cl.amount)
    , var_RED_mean = numeric(cl.amount)
    , var_YEL_mean = numeric(cl.amount)
    , cent_GRN_med = numeric(cl.amount)
    , cent_RED_med = numeric(cl.amount)
    , cent_YEL_med = numeric(cl.amount)
    , var_GRN_med = numeric(cl.amount)
    , var_RED_med = numeric(cl.amount)
    , var_YEL_med = numeric(cl.amount)
  )

  quality_Q1_mean <- 0
  quality_Q2_mean <- 0
  quality_Q3 <- 0
  quality_Q1_med <- 0
  quality_Q2_med <- 0
  
  for( cl.idx in 1:cl.amount ){
    is.lbl <- labels == cl.lbls[cl.idx]
    cur.data <- data[is.lbl, ]
    
    cl_cell_amount <- sum( is.lbl )
    
    cent_mean <- colMeans( cur.data )
    cent_med <- apply( cur.data, MARGIN = 2, FUN = median )
    var_mean <- apply( cur.data, MARGIN = 2, FUN = var )
    var_med <-  
      apply( cur.data - cent_med, MARGIN = 2, FUN = function(x){
        sum( x^2 )/length(x)
      })

    cmean.dist <- dist( rbind( cent_mean, cur.data ) )
    cl_radius_mean <- max( cmean.dist[1:nrow(cur.data)] )
    cl_variance_mean <- mean( cmean.dist[1:nrow(cur.data)]^2 )
    quality_Q1_mean <- quality_Q1_mean + sum( cmean.dist[1:nrow(cur.data)] )
    quality_Q2_mean <- quality_Q2_mean + sum( cmean.dist[1:nrow(cur.data)]^2 )
    
    cmed.dist <- dist( rbind( cent_med, cur.data ) )
    cl_radius_med <- max( cmed.dist[1:nrow(cur.data)] )
    cl_variance_med <- mean( cmed.dist[1:nrow(cur.data)]^2 )
    quality_Q1_med <- quality_Q1_med + sum( cmed.dist[1:nrow(cur.data)] )
    quality_Q2_med <- quality_Q2_med + sum( cmed.dist[1:nrow(cur.data)]^2 )
    
    quality_Q3 <- quality_Q3 + sum(dist(cur.data))
    
    
    tib.cur.cl.inf[ cl.idx , c(  'cl_cell_amount'
                               , 'cl_radius_mean'
                               , 'cl_radius_med'
                               , 'cl_variance_mean'
                               , 'cl_variance_med' 
                               , 'cent_GRN_mean'
                               , 'cent_RED_mean'
                               , 'cent_YEL_mean' 
                               , 'var_GRN_mean'
                               , 'var_RED_mean'
                               , 'var_YEL_mean'
                               , 'cent_GRN_med'
                               , 'cent_RED_med'
                               , 'cent_YEL_med'
                               , 'var_GRN_med'
                               , 'var_RED_med'
                               , 'var_YEL_med') ] <- 
      c(   cl_cell_amount
         , cl_radius_mean
         , cl_radius_med
         , cl_variance_mean
         , cl_variance_med
         , cent_mean[1:3]
         , var_mean[1:3]
         , cent_med[1:3]
         , var_med[1:3]
         )
  }
  
  quality_Q1_mean <- quality_Q1_mean / nrow( data )
  quality_Q2_mean <- quality_Q2_mean / nrow( data )
  quality_Q3 <- quality_Q3 / nrow( data )
  quality_Q1_med <- quality_Q1_med / nrow( data )
  quality_Q2_med <- quality_Q2_med / nrow( data )
  
  tib.cur.cl.inf$quality_Q1_mean  <- quality_Q1_mean
  tib.cur.cl.inf$quality_Q2_mean <- quality_Q2_mean
  tib.cur.cl.inf$quality_Q3 <- quality_Q3
  tib.cur.cl.inf$quality_Q1_med <- quality_Q1_med 
  tib.cur.cl.inf$quality_Q2_med <- quality_Q2_med 
  
  return( tib.cur.cl.inf )
}


