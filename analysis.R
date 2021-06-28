library(ggplot2)
library(dplyr)
library(PerformanceAnalytics)

hist( tib.cluster.info$cent_GRN_mean - tib.cluster.info$cent_GRN_med, breaks = 30 )
hist( tib.cluster.info$cent_YEL_mean - tib.cluster.info$cent_YEL_med, breaks = 30 )
hist( tib.cluster.info$cent_RED_mean - tib.cluster.info$cent_RED_med, breaks = 30 )

pairs( tib.cluster.info[, c(12,13,14)] )
system.time({pairs.panels(tib.cluster.info[1:4e4, c(12,13,14)], scale = TRUE, ellipses = T, lm = T, density = T )})
system.time( {chart.Correlation(tib.cluster.info[1:4e4, c(12,13,14)], histogram = TRUE) })



unique(tib.cell.surv$substance)

ungroup(tib.cell.surv)
tib.cell.surv %>%
  dplyr::filter( substance == '0'
                 , tr_time == 16
                 , NOTE == 'not_replicable') %>%
  summarise( count = length( tr_time ))




# plot --------------------------------------------------------------------

sort(unique(tib.cell.surv$substance))

ggplot( tib.cell.surv %>%
          dplyr::filter( substance %in% c('g5'
                                        , 'g5+p5'
                                        , 'g5+p10'
                                        , 'g5+p100'
                                        , 'g5+p150'
                                        , 'g5+p250'
                                        , 'g5+p50'
                                        , 'g5+p500'
                                        )
                      #   , tr_time %in% c(16)
                         , cl_def_approach %in% c('no_ajustment', 'm2')
                         , !(cl_method %in% c('median', 'mcquitty', 'centroid', 'average', 'complete') )
                         , cl_amount %in% c(2, 3)
                         , dim_amount %in% c(2,3)
                         )
        , aes( as.factor(substance) , dead_cells_perc)
          ) +
  geom_violin(fill = "grey80"
              , colour = "black"
              , draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_jitter( aes( colour = cl_method )) +
  ylim(0,1)

  

# plot 3 ------------------------------------------------------------------

ggplot( tib.cell.surv %>%
          dplyr::filter( 1==1
          #  substance %in% c('g8')
          #   , tr_time %in% c(16)
          , pretr_time == 0
          , cl_def_approach %in% c('!no_ajustment', 'm2')
          , !(cl_method %in% c('median', 'mcquitty', 'centroid', 'average', 'complete') )
          , cl_amount %in% c(2, 3)
          , dim_amount %in% c(2,3)
          )
        , aes( substance , dead_cells_perc)
) +
  geom_violin(fill = "grey80"
              , colour = "black"
              , draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_jitter( size = 1.5, aes( shape = tr_time, color = exper_name )) +
  #scale_color_gradient2(midpoint = 13, low="blue", mid = 'white', high="red") +
  ylim(0,1) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  stat_summary(fun.y=median, geom="line", aes(group=1))  + 
  stat_summary(fun.y=median, geom="point")

ggsave(filename = 'overview.png', device = 'png', width = 28, height = 14)
  
  
# plot 2 ------------------------------------------------------------------

chart.Correlation(tib.cell.surv %>%
                    select( quality_Q1_mean
                            , quality_Q2_mean
                            , quality_Q3
                            , dead_cells_perc
                            , dead_cells_amount)
                    , histogram = TRUE)


chart.Correlation(cbind(tib.cell.surv %>%
                    filter( cl_amount == 2
                            , cl_def_approach %in% c('no_ajustment', 'm2')
                            , !(cl_method %in% c('median', 'mcquitty', 'centroid', 'average', 'complete') )
                            ) %>%
                    select( dead_cells_amount)
                  , tib.cell.surv %>%
                    filter( cl_amount == 3
                            , cl_def_approach %in% c('no_ajustment', 'm2')
                            , !(cl_method %in% c('median', 'mcquitty', 'centroid', 'average', 'complete') )
                            ) %>%
                    select( dead_cells_amount) )
                  , histogram = TRUE)


# dsfa --------------------------------------------------------------------
