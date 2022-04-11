#assessment of cell lines

assess_and_compare_files <- function(lineage_str,file_paths){
  #designed to take as input (paths to) new csv files and
  #plot some diagnostics compared to existing cell lines
  library(dplyr)
  library(lazychromosomes)
  library(zoo)
  library(ggplot2)
  library(patchwork)
  
  assess_laziness_for_cell <- function(job_str){
    dat <- process_jobset(job_str,K=Inf,max_missing=0.95) %>%
      mutate(filename=job_str)
    in_anaphase_df <- dat %>% group_by(Frame,SisterPairID) %>%
      mutate(kk_dist = sqrt(diff(Position_1)^2+diff(Position_2)^2+diff(Position_3)^2)) %>% 
      group_by(Frame) %>%
      summarise(kk_dist=max(kk_dist,na.rm=T)) %>% 
      filter(kk_dist>1.8) #more than 1.8um separation
    t_ana_frame <- min(in_anaphase_df$Frame) + 15 #very rough estimation to give something to work with
    laziness_df <- annotate_anaphase_onset_for_cell(dat,method="manual",t_ana_frame=t_ana_frame) %>%
      get_laziness() 
    return(laziness_df)
  }
  
  assess_files <- function(files){
    files %>% 
      purrr::map_dfr(assess_laziness_for_cell,.id="cell") %>%
      mutate(filename=files[as.integer(cell)])
  }
  
  
  compute_directional_correlation_statistic <- function(job_str){
    #see armond et al 2015, supp material section 2.5
    #load data  
    dat <- process_jobset(job_str,K=Inf,max_missing=0.25) %>%
      mutate(filename=job_str)
    in_anaphase_df <- dat %>% group_by(Frame,SisterPairID) %>%
      mutate(kk_dist = sqrt(diff(Position_1)^2+diff(Position_2)^2+diff(Position_3)^2)) %>% 
      group_by(Frame) %>%
      summarise(kk_dist=max(kk_dist,na.rm=T)) %>% 
      filter(kk_dist>1.8) #more than 1.8um separation
    K_max <- min(in_anaphase_df$Frame) - 15
    dat <- process_jobset(job_str,K=K_max,max_missing=0.25)
    dat %>% 
      group_by(SisterPairID,Frame) %>%
      summarise(centre_position = mean(Position_1),
                y=mean(Position_2),
                z=mean(Position_3)) %>%
      group_by(SisterPairID) %>%
      mutate(change = c(NA,diff(centre_position)),
             sign_of_change = sign(change),
             dsign_of_change = abs(c(NA,diff(sign_of_change)))/2) %>% 
      summarise(k = sum(dsign_of_change,na.rm=T),
                n = max(Frame),
                q = pbinom(k,size=n,prob=0.5),
                directional_corr = 1-2*min(q,1-q),
                y = mean(y,na.rm=T),
                z = mean(z,na.rm=T))
  }
  
  assess_bad_oscillators <- function(job_str,oscillator_cutoff=0.4,K=Inf){
    stopifnot(file.exists(job_str))
    dat <- process_jobset(job_str,K=K,max_missing=0.25,plot_opt=0)
    
    #can input data that includes anaphase but only want to fit to metaphase data
    #use kk distance as a crude approximation to indicate where metaphase is to fit to this
    in_anaphase_df <- dat %>% group_by(Frame,SisterPairID) %>%
      mutate(kk_dist = sqrt(diff(Position_1)^2+diff(Position_2)^2+diff(Position_3)^2)) %>% 
      group_by(Frame) %>%
      summarise(kk_dist=max(kk_dist,na.rm=T)) %>% 
      filter(kk_dist>1.8) #more than 1.8um separation
    K_max <- min(in_anaphase_df$Frame) - 15 #go 15 frames (~30s) earlier than this
    if (K_max<20){
      return(tibble())
    }
    dat <- process_jobset(job_str,K=K_max,max_missing=0.25,plot_opt=0)
    
    pairIDs <- dat$SisterPairID %>% unique()
    
    #assess oscillators
    bad_oscillators <- dat %>%
      group_by(SisterPairID,SisterID) %>%
      summarise(av_x = mean(Position_1,na.rm=T),
                DAP = mean(abs(Position_1-av_x),na.rm=T),
                poor_oscillator = DAP < oscillator_cutoff)
    return(bad_oscillators)
  }
  
  #################
  
  
  
  MC212_files <- list.files(here::here("data/MC212/"),full.names = TRUE, recursive = TRUE, pattern = ".csv")
  
  MC191_files <- list.files(here::here('data/MC191/'),
                            full.names = TRUE,recursive = TRUE,pattern=".csv")
  MC194_files <- list.files("data/MC194/",
                            recursive = TRUE,full.names=TRUE,pattern="ome.csv")
  
  tracking_files <- list(input = file_paths,
                         MC212=MC212_files,
                         MC191=MC191_files,
                         MC194=MC194_files)
  names(tracking_files) <- c(lineage_str,"MC212","MC191","MC194")
  mdf <- tracking_files %>%
    purrr::map_dfr(assess_files,.id="lineage") 
  
  # mdf %>% 
  #   group_by(cell,SisterPairID,lineage) %>% 
  #   summarise(laziness=max(laziness,na.rm=T),
  #             opposite_laziness=max(opposite_laziness,na.rm=T)) %>% 
  #   mutate(is_lazy = (laziness>3.0) & (opposite_laziness>3.0)) %>%
  #   group_by(lineage,cell) %>% 
  #   summarise(total=n(),
  #             num_lazy = sum(is_lazy)) %>%
  #   mutate(any_lazy = num_lazy>0) %>%
  #   group_by(lineage) %>%
  #   summarise(total=n(),
  #             num_lazy_cells = sum(any_lazy),
  #             proportion=num_lazy_cells/total)
  
  
  lazy_plt <- mdf %>% 
    group_by(cell,SisterPairID,lineage) %>% 
    summarise(laziness=max(laziness,na.rm=T),
              opposite_laziness=max(opposite_laziness,na.rm=T)) %>% 
    filter(opposite_laziness>2.0) %>% 
    ggplot(aes(laziness,color=lineage)) + 
    stat_ecdf(geom="step") + 
    theme_bw() + 
    coord_cartesian(xlim = c(1.8, 5),
                    ylim=c(0.75,1)) + 
    scale_color_manual(values=RColorBrewer::brewer.pal(n=4,"PiYG")) + 
    labs(x="Maximum laziness",y="Empirical CDF",color="Lineage") + 
    theme(legend.position="bottom")
  lazy_plt 
  ggsave(here::here("plots/laziness_comparison.eps"),width=210,height=100,units="mm")
  
  directional_correlation_df <- tracking_files %>%
    purrr::map_dfr(function(x) {purrr::map_dfr(x,compute_directional_correlation_statistic,.id="cell")},.id="lineage") 
  
  p1 <- ggplot(directional_correlation_df, aes(directional_corr,color=lineage)) +
    stat_ecdf(geom = "step") + 
    theme_bw() + 
    scale_color_manual(values=RColorBrewer::brewer.pal(n=4,"PiYG")) +
    labs(x="Directional correlation",
         y="Empirical CDF",
         color="Lineage") + 
    theme(legend.position="bottom")
  
  p2 <- ggplot(directional_correlation_df, aes(y,z,color=directional_corr)) + 
    geom_point(size=0.5) + 
    theme_bw() + 
    facet_wrap(.~lineage) + 
    theme(legend.position="bottom") + 
    labs(x="y position (um)",
         y="z position (um)",
         color="Directional\ncorrelation ")
  p1 | p2 
  ggsave(here::here("plots/directional_correlation_comparison.eps"),
         width=210,height=100,units="mm")
  
  # bad_oscillators_MC212_df <- purrr::map(MC212_files, assess_bad_oscillators) %>% 
  #   bind_rows(.id="cell")
  # bad_oscillators_MC191_df <- purrr::map(MC191_files, assess_bad_oscillators) %>% 
  #   bind_rows(.id="cell")
  # bad_oscillators_MC194_df <- purrr::map(MC194_files, assess_bad_oscillators) %>% 
  #   bind_rows(.id="cell")
  bad_oscillators_df <- purrr::map_dfr(tracking_files,function(x) purrr::map(x, assess_bad_oscillators) %>% 
                                         bind_rows(.id="cell"), .id="lineage")
  
  plt <- ggplot(bad_oscillators_df, aes(lineage,DAP)) + 
    geom_jitter(alpha=0.3) + 
    geom_violin(draw_quantiles=0.5) + 
    theme_bw() + 
    labs(x="Cell lineage",y="DAP (um)")
  plt
  ggsave(here::here("plots/deviation_from_average_position.eps"),
         width=6,height=4)
  
  # q1 <- inner_join(directional_correlation_df,bad_oscillators_df,
  #            by=c("cell","SisterPairID","lineage")) %>% 
  #   ggplot(aes(directional_corr,DAP,color=n)) + 
  #   geom_point() + 
  #   geom_smooth(method="lm") + 
  #   theme_bw() + facet_wrap(.~lineage)
  
  # q1 <- inner_join(directional_correlation_df,bad_oscillators_df,
  #                  by=c("cell","SisterPairID","lineage")) %>% 
  #   ggplot(aes(n,color=lineage)) + geom_density() + theme_bw() 
  # q1prime <- inner_join(directional_correlation_df,bad_oscillators_df,
  #                       by=c("cell","SisterPairID","lineage")) %>% 
  #   ggplot(aes(n,DAP,color=directional_corr)) + 
  #   geom_point() + theme_bw() + facet_wrap(.~lineage) 
  # q1 <- q1 / q1prime
  
  # q2 <- inner_join(directional_correlation_df,bad_oscillators_df,
  #            by=c("cell","SisterPairID","lineage")) %>%
  # ggplot(aes(y,z,color=DAP)) + 
  #   geom_point(size=0.5) + 
  #   theme_dark() + 
  #   facet_wrap(.~lineage) + 
  #   theme(legend.position="bottom") + 
  #   labs(x="y position (um)",
  #        y="z position (um)",
  #        color="DAP (um)") + 
  #   scale_colour_gradient2(
  #     low = RColorBrewer::brewer.pal(n=4,"PiYG")[1],
  #     mid = "white",
  #     high = RColorBrewer::brewer.pal(n=4,"PiYG")[4],
  #     midpoint = 0.4
  #   )
  
  #lazy_plt / {p1 | p2} / plt / {q1 | q2}
  wrap_plots(lazy_plt,p1,p2,plt)
}


