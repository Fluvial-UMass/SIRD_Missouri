library("tidyverse")
library("ggplot2")
library("cowplot")
library("gridExtra")
library("RColorBrewer")

assignCase <- function(df, name) {
  df$case = name
  return(df)
}

assignSubCase <- function(df, name) {
  df$subCase = name
  return(df)
}

threshold = 0

# read data
assim_uc = read.csv("./stats/stats/dischargeAssim_uncal_log_stats.csv")
assim_uc = assignCase(assim_uc, "caseA")
assim_uc = assignSubCase(assim_uc, "assimilated")

assim_03 = read.csv("./stats/stats/dischargeAssim_cal03_log_stats.csv") #%>% select(-rbias)
assim_03 = assignCase(assim_03, "caseB")
assim_03 = assignSubCase(assim_03, "assimilated")

assim_10 = read.csv("./stats/stats/dischargeAssim_cal10_log_stats.csv") #%>% select(-rbias)
names(assim_10) = c("r", "rmse", "nrmse", "rrmse", "rbias", "nse", "kge", "hrrid")
assim_10 = assignCase(assim_10, "caseC")
assim_10 = assignSubCase(assim_10, "assimilated")

insert_uc = read.csv("./stats/stats/dischargeInsert_summarized_uncal_insert_stats.csv")
insert_uc = assignCase(insert_uc, "caseA")
insert_uc = assignSubCase(insert_uc, "inserted")

insert_03 = read.csv("./stats/stats/dischargeInsert_summarized_cal03_insert_stats.csv")
insert_03 = assignCase(insert_03, "caseB")
insert_03 = assignSubCase(insert_03, "inserted")

insert_10 = read.csv("./stats/stats/dischargeInsert_summarized_cal10_insert_stats.csv")
insert_10 = assignCase(insert_10, "caseC")
insert_10 = assignSubCase(insert_10, "inserted")

ref_uc = read.csv("./stats/stats/discharge_uncal_baseline_stats.csv")
ref_uc = assignCase(ref_uc, "caseA")
ref_uc = assignSubCase(ref_uc, "baseline")

ref_03 = read.csv("./stats/stats/discharge_cal10_baseline_stats.csv")
ref_03 = assignCase(ref_03, "caseB")
ref_03 = assignSubCase(ref_03, "baseline")

ref_10 = read.csv("./stats/stats/discharge_cal10_baseline_stats.csv")
names(ref_10) = c("r", "rmse", "nrmse", "rrmse", "rbias", "nse", "kge", "hrrid")
ref_10 = assignCase(ref_10, "caseC")
ref_10 = assignSubCase(ref_10, "baseline")


# filter nse < target
ref_10 = ref_10 %>% subset(nse < threshold)

assim_uc = assim_uc %>% filter(hrrid %in% ref_10$hrrid)
assim_03 = assim_03 %>% filter(hrrid %in% ref_10$hrrid)
assim_10 = assim_10 %>% filter(hrrid %in% ref_10$hrrid)

insert_uc = insert_uc %>% filter(hrrid %in% ref_10$hrrid)
insert_03 = insert_03 %>% filter(hrrid %in% ref_10$hrrid)
insert_10 = insert_10 %>% filter(hrrid %in% ref_10$hrrid)

ref_uc = ref_uc %>% filter(hrrid %in% ref_10$hrrid)
ref_03 = ref_03 %>% filter(hrrid %in% ref_10$hrrid)


# combine
all = rbind(assim_uc, assim_03, assim_10, insert_uc, insert_03, insert_10, ref_uc, ref_03, ref_10)
all_nse = all %>% select(hrrid, nse, case, subCase)
tidyall_nse = all_nse %>% pivot_longer(-c(hrrid, case, subCase), names_to="nse", values_to="values")
assim_nse = rbind(assim_uc, assim_03, assim_10) %>% select(hrrid, nse, case, subCase)
tidyassim_nse = assim_nse %>% pivot_longer(-c(hrrid, case, subCase), names_to="nse", values_to="values")
# boxplot
ggplot(tidyall_nse, aes(x=subCase, y=values, fill=case)) + geom_boxplot() + coord_cartesian(ylim=c(-20, 1))
ggplot(tidyassim_nse, aes(x=case, y=values, fill=case)) + geom_boxplot() + coord_cartesian(ylim=c(-5, 1))
# cdf
ggplot(tidyall_nse, aes(values)) + stat_ecdf(geom="step", aes(color=case, linetype=subCase), size=0.75) +
  coord_cartesian(xlim=c(-10, 1)) + theme_light() + xlab("NSE") + ylab("density")

# histgram
takeDifferenceNse <- function(df_exp, df_ref, case, subCase) {
  df_exp_s = df_exp[order(df_exp$hrrid), ] %>% select(nse)
  df_ref_s = df_ref[order(df_exp$hrrid), ] %>% select(nse)
  diff = (df_exp_s - df_ref_s)
  diff$case = case
  diff$subCase = subCase
  return(diff)
}

caseA_a=takeDifferenceNse(assim_uc, ref_uc, "caseA", "assimilated")
caseB_a=takeDifferenceNse(assim_03, ref_03, "caseB", "assimilated")
caseC_a=takeDifferenceNse(assim_10, ref_10, "caseC", "assimilated")
caseA_i=takeDifferenceNse(insert_uc, ref_uc, "caseA", "inserted")
caseB_i=takeDifferenceNse(insert_03, ref_03, "caseB", "inserted")
caseC_i=takeDifferenceNse(insert_10, ref_10, "caseC", "inserted")
histtidydf = rbind(caseA_a, caseB_a, caseC_a, caseA_i, caseB_i, caseC_i) %>%
  pivot_longer(-c(case, subCase), names_to="nse", values_to="values")
histtidydf$subCase =  factor(histtidydf$subCase, levels = c("inserted", "assimilated"), ordered = TRUE)
box = ggplot(histtidydf, aes(x=case, y=values, fill=subCase)) + geom_boxplot() + coord_cartesian(ylim=c(-2.5, 5)) +
  xlab("experiment case") + ylab("NSE change") +
  labs(title=expression("baseline NSE"<~"0 [N=225]")) +
  theme_light() + scale_fill_discrete(name="sub case") +
  theme(axis.text.x=element_text(size=18),
        axis.text.y=element_text(size=18),
        axis.title=element_text(size=20),
        plot.title=element_text(size=15),
        legend.text=element_text(size=15),
        legend.title=element_blank())
#
# nrow(caseA_a %>% subset(nse < -2.5)) 7
# nrow(caseB_a %>% subset(nse < -2.5)) 11
# nrow(caseC_a %>% subset(nse < -2.5)) 9
# nrow(caseA_a %>% subset(nse > 5)) 23
# nrow(caseB_a %>% subset(nse > 5)) 24
# nrow(caseC_a %>% subset(nse > 5)) 24

# nrow(caseA_i %>% subset(nse < -2.5)) 23
# nrow(caseB_i %>% subset(nse < -2.5)) 0
# nrow(caseC_i %>% subset(nse < -2.5)) 0
# nrow(caseA_i %>% subset(nse > 5)) 0
# nrow(caseB_i %>% subset(nse > 5)) 1
# nrow(caseC_i %>% subset(nse > 5)) 0

nrow(caseA_a %>% subset(nse < 0)) # 22
nrow(caseB_a %>% subset(nse < 0)) # 27
nrow(caseC_a %>% subset(nse < 0)) # 22
nrow(caseA_i %>% subset(nse < 0)) # 186
nrow(caseB_i %>% subset(nse < 0)) # 48
nrow(caseC_i %>% subset(nse < 0)) # 111
box

caseA_a=takeDifferenceNse(assim_uc, ref_uc, "caseA", "assimilated")
caseB_a=takeDifferenceNse(assim_03, ref_03, "caseB", "assimilated")
caseC_a=takeDifferenceNse(assim_10, ref_10, "caseC", "assimilated")
densetidydf = rbind(caseA_a, caseB_a, caseC_a) %>%
  pivot_longer(-c(case, subCase), names_to="nse", values_to="values")
density = ggplot(densetidydf%>%filter(subCase == "assimilated"), aes(values, color=case)) + geom_density(position="identity", fill="black", alpha=0.15, size=1) + xlim(-5, 10) +
  xlab("NSE change (assimilated)") + ylab("density") +
  labs(title=expression("baseline NSE"<~"0 [N=225]")) +
  theme_light() + scale_fill_brewer(palette="Paired") +
  theme(axis.text.x=element_text(size=18),
        axis.text.y=element_text(size=18),
        axis.title=element_text(size=20),
        plot.title=element_text(size=15),
        legend.text=element_text(size=15),
        legend.title=element_blank())
density

out = plot_grid(box, density, nrow=2, labels=c("a", "b"))
out
save_plot("./caseboxplot/assim_insert_density.png", out, nrow=2, base_height =5, base_width = 12, dpi=300)
