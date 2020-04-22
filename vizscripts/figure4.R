library("tidyverse")
library("ggplot2")
library("cowplot")
library("gridExtra")
library("RColorBrewer")

threshold = 0

makecdf_KGE <- function(){
  exp_forTidy = read.csv("./stats/stats/dischargeAssim_cal10_log_stats.csv")
  ref_forTidy = read.csv("./stats/stats/discharge_cal10_baseline_stats.csv")
  # meanflow = read.csv("./average_discharge_gauge.csv")
  names(exp_forTidy) = c("r", "RMSE", "NRMSE", "RRMSE", "RBIAS", "NSE", "KGE", "hrrid")
  exp_forTidy$case = "assimilated"
  # exp_forTidy = left_join(exp_forTidy, meanflow, by="hrrid")
  names(ref_forTidy) = c("r", "RMSE", "NRMSE", "RRMSE", "RBIAS", "NSE", "KGE", "hrrid")
  ref_forTidy$case = "baseline"
  # ref_forTidy = left_join(ref_forTidy, meanflow, by="hrrid")
  
  ref_forTidyb = ref_forTidy %>% subset(NSE < threshold) #%>% subset(NSE > -10000)
  exp_forTidyb = exp_forTidy %>% filter(hrrid %in% ref_forTidyb$hrrid) #%>% subset(NSE > -10000)
  ref_forTidyb = ref_forTidyb %>% filter(hrrid %in% exp_forTidyb$hrrid)
  all_forTidyb = rbind(exp_forTidyb, ref_forTidyb)
  ggtidyb = all_forTidyb %>% select(hrrid, case, KGE) %>% pivot_longer(c(-hrrid, -case), values_to="values", names_to="stats")
  ggtidyb$case <- factor(ggtidyb$case, levels = c("assimilated", "baseline"), ordered = TRUE)
  
  # boxplot
  # ggplot(ggtidyb, aes(x=stats, y=values, fill=case)) + geom_boxplot() + coord_cartesian(ylim=c(-10, 10))
    # ggplot(ggtidy, aes(x=stats, y=values, fill=case)) + geom_boxplot() + ylim(-5, 5)
  
  cdf = ggplot(ggtidyb, aes(values)) + stat_ecdf(geom="step", aes(linetype=case), size=0.75, color="#F1C40F") +
    coord_cartesian(xlim=c(-5, 1)) + theme_light() + xlab("values") + ylab("density") +
    ggtitle(expression("CDF [KGE]: baseline NSE"<~"0 [N=225]")) +
    theme(axis.text.x=element_text(size=18),
          axis.text.y=element_text(size=18),
          axis.title=element_text(size=20),
          plot.title=element_text(size=15),
          legend.position=c(0.15,0.825),
          legend.text=element_text(size=15),
          legend.title=element_blank())
  return(cdf)
}

makecdf_NRMSE <- function(){
  exp_forTidy = read.csv("./stats/stats/dischargeAssim_cal10_log_stats.csv")
  ref_forTidy = read.csv("./stats/stats/discharge_cal10_baseline_stats.csv")
  # meanflow = read.csv("./average_discharge_gauge.csv")
  names(exp_forTidy) = c("r", "RMSE", "NRMSE", "RRMSE", "RBIAS", "NSE", "KGE", "hrrid")
  exp_forTidy$case = "assimilated"
  # exp_forTidy = left_join(exp_forTidy, meanflow, by="hrrid")
  names(ref_forTidy) = c("r", "RMSE", "NRMSE", "RRMSE", "RBIAS", "NSE", "KGE", "hrrid")
  ref_forTidy$case = "baseline"
  # ref_forTidy = left_join(ref_forTidy, meanflow, by="hrrid")
  
  ref_forTidyb = ref_forTidy %>% subset(NSE < threshold) #%>% subset(NSE > -10000)
  exp_forTidyb = exp_forTidy %>% filter(hrrid %in% ref_forTidyb$hrrid) #%>% subset(NSE > -10000)
  ref_forTidyb = ref_forTidyb %>% filter(hrrid %in% exp_forTidyb$hrrid)
  all_forTidyb = rbind(exp_forTidyb, ref_forTidyb)
  ggtidyb = all_forTidyb %>% select(hrrid, case, NRMSE) %>% pivot_longer(c(-hrrid, -case), values_to="values", names_to="stats")
  ggtidyb$case <- factor(ggtidyb$case, levels = c("assimilated", "baseline"), ordered = TRUE)
  
  # boxplot
  # ggplot(ggtidyb, aes(x=stats, y=values, fill=case)) + geom_boxplot() + coord_cartesian(ylim=c(-10, 10))
  # ggplot(ggtidy, aes(x=stats, y=values, fill=case)) + geom_boxplot() + ylim(-5, 5)
  
  cdf = ggplot(ggtidyb, aes(values)) + stat_ecdf(geom="step", aes(linetype=case), size=0.75, color="#2E4053") +
    coord_cartesian(xlim=c(0, 5)) + theme_light() + xlab("values") + ylab("density") +
    ggtitle(expression("CDF [NRMSE]: baseline NSE"<~"0 [N=225]")) +
    theme(axis.text.x=element_text(size=18),
          axis.text.y=element_text(size=18),
          axis.title=element_text(size=20),
          plot.title=element_text(size=15),
          legend.position=c(0.15,0.825),
          legend.text=element_text(size=15),
          legend.title=element_blank())
  return(cdf)
}

makediffhist_kge <- function(nsethreshold=0, flag="below"){
  exp = read.csv("./stats/stats/dischargeAssim_cal10_log_stats.csv")
  ref = read.csv("./stats/stats/discharge_cal10_baseline_stats.csv")
  # meanflow = read.csv("./average_discharge_gauge.csv")
  all = left_join(exp, ref, by="hrrid")
  
  if (flag == "up") {
    all_sel_below = all %>% subset(nse_baseline >= nsethreshold) #%>% subset(nse_assimilated > -10000)
    hist_below = ggplot(all_sel_below, aes(x=(kge_assimilated-kge_baseline))) + 
      geom_histogram(bins=30, fill="#F1C40F", color="black") + xlim(-2.5,5) +
      theme_light() + xlab("improvement in KGE") + ylab("count") + 
      labs(title=expression("CDF [NRMSE]: baseline NSE">=~"0 [N=178]")) +
      theme(axis.text.x=element_text(size=18),
            axis.text.y=element_text(size=18),
            axis.title=element_text(size=20),
            plot.title=element_text(size=15),
            legend.text=element_text(size=15),
            legend.title=element_text(size=18))
    print(summary(all_sel_below$kge_assimilated-all_sel_below$kge_baseline))
  } else {
    all_sel_below = all %>% subset(nse_baseline < nsethreshold) #%>% subset(nse_assimilated > -10000)
    hist_below = ggplot(all_sel_below, aes(x=(kge_assimilated-kge_baseline))) + 
      geom_histogram(bins=30, fill="#F1C40F", color="black") + xlim(-2.5,5) +
      theme_light() + xlab("improvement in KGE") + ylab("count") + 
      labs(title=expression("CDF [NRMSE]: baseline NSE"<~"0 [N=225]")) +
      theme(axis.text.x=element_text(size=18),
            axis.text.y=element_text(size=18),
            axis.title=element_text(size=20),
            plot.title=element_text(size=15),
            legend.text=element_text(size=15),
            legend.title=element_text(size=18))
    print(summary(all_sel_below$kge_assimilated-all_sel_below$kge_baseline))
  }
  
  all_sel_below$impr = (all_sel_below$kge_assimilated-all_sel_below$kge_baseline)
  print(nrow(all_sel_below))
  print(nrow(all_sel_below[all_sel_below$impr > 0, ]))
  return(hist_below)
}

makediffhist_nrmse <- function(nsethreshold=0, flag="below"){
  exp = read.csv("./stats/stats/dischargeAssim_cal10_log_stats.csv")
  ref = read.csv("./stats/stats/discharge_cal10_baseline_stats.csv")
  # meanflow = read.csv("./average_discharge_gauge.csv")
  all = left_join(exp, ref, by="hrrid")
  if (flag == "up") {
    all_sel_below = all %>% subset(nse_baseline >= threshold) #%>% subset(nse_assimilated > -10000)
    hist_below = ggplot(all_sel_below, aes(x=-1*(nrmse_assimilated-nrmse_baseline))) +
      geom_histogram(bins=30, fill="#2E4053", color="black") + xlim(-2.5,5) +
      theme_light() + xlab("improvement in NRMSE") + ylab("count") +
      labs(title=expression("CDF [NRMSE]: baseline NSE">=~"0 [N=178]")) +
      theme(axis.text.x=element_text(size=18),
            axis.text.y=element_text(size=18),
            axis.title=element_text(size=20),
            plot.title=element_text(size=15),
            legend.text=element_text(size=15),
            legend.title=element_text(size=18))
    print(summary(-1*(all_sel_below$nrmse_assimilated-all_sel_below$nrmse_baseline)))
  } else {
    all_sel_below = all %>% subset(nse_baseline < nsethreshold) #%>% subset(nse_assimilated > -10000)
    hist_below = ggplot(all_sel_below, aes(x=-1*(nrmse_assimilated-nrmse_baseline))) +
      geom_histogram(bins=30, fill="#2E4053", color="black") + xlim(-2.5,5) +
      theme_light() + xlab("improvement in NRMSE") + ylab("count") +
      labs(title=expression("CDF [NRMSE]: baseline NSE"<~"0 [N=225]")) +
      theme(axis.text.x=element_text(size=18),
            axis.text.y=element_text(size=18),
            axis.title=element_text(size=20),
            plot.title=element_text(size=15),
            legend.text=element_text(size=15),
            legend.title=element_text(size=18))
    print(summary(-1*(all_sel_below$nrmse_assimilated-all_sel_below$nrmse_baseline)))
  }
  
  all_sel_below$impr = -1*(all_sel_below$nrmse_assimilated-all_sel_below$nrmse_baseline)
  print(summary(all_sel_below))
  print(nrow(all_sel_below[all_sel_below$impr > 0, ]))
  return(hist_below)
}

cdf_kge = makecdf_KGE()
cdf_kge
cdf_nrmse = makecdf_NRMSE()
cdf_nrmse
hist_kge_b = makediffhist_kge(threshold, flag="below")
hist_kge_b
hist_kge_u = makediffhist_kge(threshold, flag="up")
hist_kge_u
hist_nrmse_b = makediffhist_nrmse(threshold, flag="below")
hist_nrmse_b
hist_nrmse_u = makediffhist_nrmse(threshold, flag="up")
hist_nrmse_u
# outhist_below = plot_grid(hist_kge_b, hist_nrmse_b, hist_kge_u, hist_nrmse_u, ncol=2, labels=c("b", "c", "d", "e"))
# outhist_below
# outhist = plot_grid(cdf, outhist_below, nrow=2, labels=c("a", NA))
outhist = plot_grid(cdf_kge,
                    cdf_nrmse,
                    hist_kge_b,
                    hist_nrmse_b,
                    hist_kge_u,
                    hist_nrmse_u,
                    labels=c("a", "b", "c", "d", "e", "f"), ncol=2, nrow=3)
outhistl <- plot_grid(outhist) + draw_grob(legend, x=-0.09, y=0.89, width=0.25, height=0.25)
outhistl
save_plot("final/kge_nrmse_improvement.png", outhistl, base_height=12, base_width=15, dpi=300)
