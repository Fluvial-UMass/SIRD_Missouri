library("tidyverse")
library("ggplot2")
library("cowplot")
library("gridExtra")
library("RColorBrewer")


lm_eqn <- function(df){
  df = df %>% subset(nse_baseline > -5) %>% subset(nse_assimilated-nse_baseline > -2.5) %>%
    subset(nse_assimilated-nse_baseline < 2.5)
  y=df$nse_baseline
  x=df$nse_assimilated - df$nse_baseline
  data = data.frame(x=x, y=y)
  m <- lm(y ~ x, data)
  b = format(unname(coef(m)[2]), digits = 3)
  r2 = format(summary(m)$r.squared, digits = 2)
  print(r2)
  return(c(b, r2))
}

threshold = 0

# exp = read.csv("./meritvic_syears_assim_cal10_stats.csv")
exp = read.csv("./stats/stats/dischargeAssim_cal10_log_stats.csv")
# ref = read.csv("./tmpDischarge_10_stats.csv")
ref = read.csv("./stats/stats/discharge_cal10_baseline_stats.csv")
meanflow = read.csv("./average_discharge_gauge.csv")
all = left_join(exp, ref, by="hrrid")
eqn = lm_eqn(all)
rsq = as.expression(substitute(italic(r)^2~"="~r2, 
                            list(r2 = as.numeric(eqn[2]))))
scatter = ggplot(all) + geom_smooth(method="lm", aes(y=(nse_assimilated-nse_baseline), x=nse_baseline)) +
  geom_rect(data=data.frame(ymin=-Inf, ymax=Inf, xmin=-Inf, xmax=threshold),
            aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="#CB4335", alpha=0.75) +
  geom_rect(data=data.frame(ymin=-Inf, ymax=Inf, xmin=threshold, xmax=Inf),
            aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="#2E86C1", alpha=0.75) +
  labs(title="N=403") +
  geom_point(aes(y=(nse_assimilated-nse_baseline), x=nse_baseline)) + ylim(-1, 2.5) + xlim(-5, 1) +
  ylab("improvement in NSE") + xlab("baseline NSE") + theme(axis.title=element_text(face="bold")) + theme_light() +
  theme(axis.text.x=element_text(size=18),
        axis.text.y=element_text(size=18),
        axis.title=element_text(size=20),
        plot.title=element_text(size=15),
        legend.text=element_text(size=15),
        legend.title=element_text(size=18)) +
  geom_text(x = 0.5, y = 2.25, label = paste("Slope: ", eqn[1]), parse = TRUE, size=6) +
  geom_text(x = 0.5, y = 2.00, label = rsq, parse = TRUE, size=6)
# add R2 send figures 
scatter 

all_sel_below = all %>% subset(nse_baseline < threshold)
hist_below = ggplot(all_sel_below, aes(x=(nse_assimilated-nse_baseline))) + geom_histogram(bins=30, fill="#CB4335") + xlim(-1,5) +
  xlab(expression("improvement in NSE (baseline NSE"<~"0)")) + ylab("count") + theme(axis.title=element_text(face="bold")) + theme_light() +
  labs(title=expression("baseline NSE"<~"0 [N=225]")) +
  theme(axis.text.x=element_text(size=18),
        axis.text.y=element_text(size=18),
        axis.title=element_text(size=20),
        plot.title=element_text(size=15),
        legend.text=element_text(size=15),
        legend.title=element_text(size=18))
hist_below
all_sel_above = all %>% subset(nse_baseline >= threshold)
hist_above = ggplot(all_sel_above, aes(x=(nse_assimilated-nse_baseline))) + geom_histogram(bins=30, fill="#2E86C1") + xlim(-1,5) +
  xlab(expression("improvement in NSE (baseline NSE">=~"0)")) + ylab("count") + theme(axis.title=element_text(face="bold")) + theme_light() +
  labs(title=expression("improvement in NSE">=~"0 [N=178]")) +
  theme(axis.text.x=element_text(size=18),
        axis.text.y=element_text(size=18),
        axis.title=element_text(size=20),
        plot.title=element_text(size=15),
        legend.text=element_text(size=15),
        legend.title=element_text(size=18))
hist_above

exp_ = exp %>% select(hrrid, nse_assimilated) # %>% subset(nse_assimilated > -10000)
colnames(exp_) = c("hrrid", "nse")
exp_$case <- "assimilated"
ref_ = ref %>% select(hrrid, nse_baseline)
ref_ = ref_ %>% filter(hrrid %in% exp_$hrrid)
colnames(ref_) = c("hrrid", "nse")
ref_$case <- "baseline"
cdfdf = rbind(exp_, ref_)
cdf = ggplot(cdfdf, aes(nse)) + stat_ecdf(geom="step", aes(color=case, linetype=case), size=1.50) +
  coord_cartesian(xlim=c(-2.5, 1)) + theme_light() + theme(legend.position=c(0.15,0.825), legend.title = element_blank()) +
  xlab("NSE") + ylab("density") +
  ggtitle("CDF: baseline NSE (N=403)") +
  theme(axis.text.x=element_text(size=18),
        axis.text.y=element_text(size=18),
        axis.title=element_text(size=20),
        plot.title=element_text(size=15),
        legend.text=element_text(size=15))
cdf = cdf + scale_color_grey(start=0.15, end=0.75)
cdf


exp_forTidy = read.csv("./stats/stats/dischargeAssim_cal10_logdelta_stats.csv")
ref_forTidy = read.csv("./stats/stats/discharge_cal10_baseline_stats.csv")
meanflow = read.csv("./average_discharge_gauge.csv")
names(exp_forTidy) = c("r", "RMSE", "NRMSE", "RRMSE", "RBIAS", "NSE", "KGE", "hrrid")
exp_forTidy$case = "assimilated"
# exp_forTidy = left_join(exp_forTidy, meanflow, by="hrrid")
names(ref_forTidy) = c("r", "RMSE", "NRMSE", "RRMSE", "RBIAS", "NSE", "KGE", "hrrid")
ref_forTidy$case = "baseline"
# ref_forTidy = left_join(ref_forTidy, meanflow, by="hrrid")

ref_forTidyb = ref_forTidy %>% subset(NSE < threshold)# %>% subset(NSE > -10000)
exp_forTidyb = exp_forTidy %>% filter(hrrid %in% ref_forTidyb$hrrid)# %>% subset(NSE > -10000)
ref_forTidyb = ref_forTidyb %>% filter(hrrid %in% exp_forTidyb$hrrid)
all_forTidyb = rbind(exp_forTidyb, ref_forTidyb)
ggtidyb = all_forTidyb %>% select(hrrid, case, NSE) %>% pivot_longer(c(-hrrid, -case), values_to="values", names_to="stats")
ggtidyb$case <- factor(ggtidyb$case, levels = c("assimilated", "baseline"), ordered = TRUE)

cdfi = ggplot(ggtidyb, aes(values)) + stat_ecdf(geom="step", aes(linetype=case), size=1.50, color="#CB4335") +
  coord_cartesian(xlim=c(-10, 1)) + theme_light() + xlab("NSE") + ylab("density") + ggtitle(expression("CDF: baseline NSE"<~"0 [N=225]")) +
  theme(axis.text.x=element_text(size=18),
        axis.text.y=element_text(size=18),
        axis.title=element_text(size=20),
        plot.title=element_text(size=15),
        legend.position=c(0.15,0.825),
        legend.title = element_blank(),
        legend.text=element_text(size=15))
cdfi

out = plot_grid(scatter,
                hist_below,
                cdf,
                hist_above,
                cdfi, 
                labels=c("a","b","c","d", "e"), ncol=2, nrow=3)
legend <- get_legend(cdf +
                     theme(legend.title = element_blank(), 
                           legend.text = element_text(size = 18)) +
                     guides(shape = guide_legend(override.aes = list(size = 18))))
outl <- plot_grid(out)
# outl
save_plot(paste("./final/improvement.png"), outl, ncol=2, base_height=15, base_width =7, dpi=300)

