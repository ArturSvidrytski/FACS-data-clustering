#library(psych)
#library(PerformanceAnalytics)  
library(xlsx)
library(tidyverse)
library(reshape2)


# Combine results ---------------------------------------------------------


for( i in 1:length(lst.cl.info) ){
  tmp.cl.info <- lst.cl.info[[i]]
  print(nrow(tmp.cl.info))
  tmp.cl.info$cl_group_id <- tmp.cl.info$cl_group_id + max(tib.cluster.info$cl_group_id)
  tib.cluster.info <- union_all(tib.cluster.info , tmp.cl.info)
}
rm( tmp.cl.info )

tib.cluster.info <- arrange(tib.cluster.info
                            , exper_name
                            , well_name
                            , well_number
                            , cl_group_id)
tib.cluster.info$cl_info_id <- 1:nrow(tib.cluster.info)

#spoiled tib.cluster.info / same data in bkp
#tib.cluster.info.bkp <- tib.cluster.info

# Add descriptive data ----------------------------------------------------
tib.cluster.info <- tib.cluster.info.bkp

experiment.desc <- read.xlsx( file = 'Experiment_descriptions.xlsx', sheetIndex = 1, header = T  )
experiment.desc <- experiment.desc[,-1]
experiment.desc <-  melt( experiment.desc, id = c(1:5 ) )
names(experiment.desc)[c(6,7)] <- c( 'well_fullname', 'substance' )
experiment.desc <- experiment.desc[!is.na(experiment.desc$substance),]
experiment.desc$well_fullname <- as.character(experiment.desc$well_fullname)
experiment.desc$NOTE <- as.character(experiment.desc$NOTE)
experiment.desc$date <- as.character(experiment.desc$date)
experiment.desc$exper_name <- as.character(experiment.desc$exper_name)
experiment.desc <- as.tibble(experiment.desc)

tib.cluster.info$well_fullname <- 
  paste0(tib.cluster.info$well_name
         , as.integer(tib.cluster.info$well_number) )

tib.cluster.info <- tib.cluster.info %>%
  right_join(experiment.desc, by = c('exper_name', 'well_fullname')  )

tib.cluster.info <- rename( tib.cluster.info, cl_cell_state = cl_position)


# Cleaning ----------------------------------------------------------------
library(tidyverse)

tib.cluster.info <- tib.cluster.info %>%
  filter(   !is.na(cent_GRN_mean) 
            , !is.na(cent_RED_mean) 
            , !is.na(cent_GRN_med ) 
            , !is.na(cent_RED_med )  )

# Define dead and alive groups of cells -----------------------------------
library(tidyverse)

tib.cluster.info$cl_cell_state <- NA

tib.cluster.info <-
  union_all(
    tib.cluster.info %>%
      filter( is.na(cent_YEL_mean) ) %>%
      mutate(cl_dist_from_zero = sqrt(cent_GRN_mean^2 + 
                                        cent_RED_mean^2) )
    ,
    tib.cluster.info %>%
      filter( !is.na(cent_YEL_mean) ) %>%
      mutate(cl_dist_from_zero = sqrt(cent_GRN_mean^2 + 
                                        cent_RED_mean^2 + 
                                        cent_YEL_mean^2) )
  )

length(unique(tib.cluster.info$cl_group_id))

tib.cl_info_id.alive <- tib.cluster.info %>%
  group_by(cl_group_id) %>%
  summarise( min_cl = min(cl_dist_from_zero)) %>%
  inner_join( tib.cluster.info
              , by = c('cl_group_id'='cl_group_id'
                       , 'min_cl'='cl_dist_from_zero')  ) %>%
  select( cl_info_id )

vect.cl_info_id.alive <- pull(tib.cl_info_id.alive, cl_info_id )
rm(tib.cl_info_id.alive)

tib.cluster.info$cl_cell_state[tib.cluster.info$cl_info_id %in%
                                 vect.cl_info_id.alive ] <- 'alive'

tib.cluster.info$cl_cell_state[ is.na(tib.cluster.info$cl_cell_state) ] <- 'dead'

tib.cluster.info$cl_dist_from_zero <- NULL

tib.cluster.info <- tib.cluster.info %>%
  group_by(cl_group_id, cl_cell_state) %>%
  summarise( 
    cl_info_id = min( cl_info_id )
    , exper_name = min( exper_name )
    , well_name = min( well_name )
    , well_number = min( well_number )
    , cl_method = min( cl_method )
    , cl_metric = min( cl_metric )
    , cl_amount = min( cl_amount )
    , cl_def_approach = min( cl_def_approach )
    , dim_amount = min( dim_amount )
    , cell_amount = min( cell_amount )
    , cl_cell_amount = sum( cl_cell_amount )
    , quality_Q1_mean = min( quality_Q1_mean )
    , quality_Q2_mean = min( quality_Q2_mean )
    , quality_Q3 = min( quality_Q3 )
    , quality_Q1_med = min( quality_Q1_med )
    , quality_Q2_med = min( quality_Q2_med )
    , cl_radius_mean = sum( cl_radius_mean )
    , cl_radius_med = sum( cl_radius_med )
    , cl_variance_mean = sum( sqrt(cl_variance_mean) )^2
    , cl_variance_med = sum( sqrt(cl_variance_med) )^2
    , cent_GRN_mean = mean( cent_GRN_mean )
    , cent_RED_mean = mean( cent_RED_mean )
    , cent_YEL_mean = mean( cent_YEL_mean )
    , var_GRN_mean = sum( sqrt(var_GRN_mean) )^2
    , var_RED_mean = sum( sqrt(var_RED_mean) )^2
    , var_YEL_mean = sum( sqrt(var_YEL_mean) )^2
    , cent_GRN_med = mean( cent_GRN_med )
    , cent_RED_med = mean( cent_RED_med )
    , cent_YEL_med = mean( cent_YEL_med )
    , var_GRN_med = sum( sqrt(var_GRN_med) )^2
    , var_RED_med = sum( sqrt(var_RED_med) )^2
    , var_YEL_med = sum( sqrt(var_YEL_med) )^2
    , well_fullname = min( well_fullname )
    , date = min( date )
    , pretr_time = min( pretr_time )
    , tr_time = min( tr_time )
    , NOTE = min( NOTE )
    , substance = min( substance )
  )

tib.cluster.info$cl_info_id <- 1:nrow(tib.cluster.info)

# Prepare tibble for analysis ---------------------------------------------

tib.cell.surv <- tib.cluster.info %>%
  filter( cl_cell_state == 'alive' ) %>%
  select(
      cell_surv_id = cl_info_id
    , exper_name
    , well_name
    , well_number
    , cl_method
    , cl_metric
    , cl_amount
    , cl_def_approach
    , dim_amount
    , cell_amount
    , alive_cells_amount = cl_cell_amount
    , cl_group_id
    , quality_Q1_mean
    , quality_Q2_mean
    , quality_Q3
    , quality_Q1_med
    , quality_Q2_med
    , well_fullname
    , date
    , pretr_time
    , tr_time
    , NOTE
    , substance
  ) %>%
  inner_join(
    tib.cluster.info %>%
      filter( cl_cell_state == 'dead' ) %>%
      select(cl_group_id, dead_cells_amount = cl_cell_amount)
    , by = 'cl_group_id'
  ) %>%
  mutate(   alive_cells_perc = alive_cells_amount/cell_amount
          , dead_cells_perc = dead_cells_amount/cell_amount)

tib.cell.surv <- ungroup(tib.cell.surv)

#tib.cell.surv$NOTE[is.na(tib.cell.surv$NOTE)] <- 'replicable'

#tib.cell.surv <- tib.cell.surv.bkp
tib.cell.surv <-
  tib.cell.surv %>%
  filter(NOTE=='replicable')



sort(unique(tib.cell.surv$substance))
substance.order <- 
  c(1,39,35,40,36,37,38,41,33,34,2,14 # 13
    ,20,18,21,17,19,22,16,15,23,27,31,28,29,30 # 27
    ,32,25,26,24,9,10,12,11,13,3,7,6,8,5,4)
sort(unique(tib.cell.surv$substance))[substance.order]

tib.cell.surv$substance <-
  factor(tib.cell.surv$substance
         , levels = sort(unique(tib.cell.surv$substance))[substance.order]
         , ordered = T )
