library(mountgrass)
library(ggplot2)
library(tidyr)
library(readr)


# load data 
setwd(dir = "C:/Users/Clemt/Documents/PhD/PhD/")

showtime_g <- read_rds("Resources/2017-12-17_showtime2.rds")

summary(showtime_g)

showtime <- showtime_g %>% group_by(rain, plasticity, tau, age, par_v) %>% transmute(
  RGR = RGR,
  BM_ag = BM_ag,
  BM_bg = BM_bg,
  BM_tot = BM_ag + BM_bg+ stem,
  stem = stem,
  PAS = 1 / (1+1/as_s),
  PAR = 1 / (1+1/as_r),
  RMF = BM_bg / (BM_ag + BM_bg),
  NPP = NPP,
  GPP = NPP + respiration,
  respiration = - respiration,
  BMprev = (BM_ag + BM_bg + stem) / (1+ RGR) ,
  turnover = - (BM_ag + BM_bg + stem - BMprev + NPP),
  to_rate = turnover / BM_tot,
  growth = BM_tot - BMprev
)

showtime$plasticity <- factor(showtime$plasticity, levels = c("none", "fixed-equilibrium", "fixed", "plastic"))

showtime %>% melt(id.vars = c("age","rain", "plasticity", "tau", "par_v")) %>% ggplot(aes(age, value)) + 
  geom_line(aes(group = interaction(rain, tau, par_v), colour = rain, linetype = as.factor(tau))) + facet_grid(variable~plasticity, scales = "free")

#☺plot BM
par <- unique(showtime$par_v)[2]
showtime2  <- showtime %>% filter(par_v == par, rain == 8, tau == 0) %>% mutate(active_s = BM_ag * PAS,
                                             structure_s = BM_ag * (1-PAS),
                                             stem_s = stem * (1-RMF),
                                             stem_r = - stem * RMF,
                                             structure_r = - BM_bg * (1 - PAR),
                                             active_r = - BM_bg * PAR)
plot <- showtime2 %>%
  melt(id.vars = c("par_v", "age", "rain", "plasticity", "tau"), measure.vars = c("active_s", "structure_s", "stem_s", "stem_r", "structure_r", "active_r") ) %>%
  mutate(variable = factor(variable, levels = c("active_s", "structure_s",  "stem_s", "active_r", "structure_r","stem_r"))) %>%
  arrange(variable) %>%
  ggplot(aes(age, value)) + geom_area(aes(fill = variable)) + facet_wrap(~ plasticity ) +
  theme_mountgrass() +
  scale_fill_manual(values = c("#DA4426", "#F37820", "#FBD475", "#065938",  "#178E5B","#9CB76C"))

path <- paste0(c(getwd(), "/Thesis/2_PP/Figures/"), collapse = "")

ggsave(paste0(c(path, "BM_var_algo.pdf"), collapse = ""), plot, width = 135, heigh = 110, units = "mm")

plot2 <- showtime %>% filter(par_v == par, rain == 8)%>%
  ggplot(aes(age, BM_bg + BM_ag + stem)) + geom_line(aes(group = interaction(plasticity, tau), linetype = as.factor(tau), colour = as.factor(plasticity)), size = 1)+
  # scale_y_log10()+
  facet_wrap(~plasticity)+
  theme_mountgrass() + scale_colour_mountgrass()
plot2


plot3 <- showtime %>% filter(par_v == par, rain == 8)%>%
  melt(measure.vars = c("PAR", "PAS", "RMF")) %>%
  ggplot(aes(age, value)) + geom_line(aes(group = interaction(variable, tau), linetype = as.factor(tau), colour = variable), size = 1)+
  # scale_y_log10()+
  facet_grid(variable~plasticity, scale = "free")+
  theme_mountgrass() + scale_colour_mountgrass()
plot3

showtime %>% filter(tau == 0) %>% melt(measure.vars = c("RMF", "PAR", "PAS")) %>% ggplot(aes(age, value)) + facet_grid(variable~plasticity, scales = "free") +
  geom_line(aes(colour = as.factor(rain), group = interaction(par_v, tau, plasticity ,rain), linetype = as.factor(tau)), alpha = 0.2) +
  theme_mountgrass() + scale_colour_mountgrass()



showtime_med <- showtime %>% group_by(rain, plasticity, tau, age) %>% summarise_all(median)

showtime_med %>% ggplot(aes(age, BM_tot)) + facet_grid(rain~plasticity) + geom_line(aes(colour = as.factor(tau)))+
  theme_mountgrass() + scale_colour_mountgrass()
showtime_med %>% ggplot(aes(age, PAS)) + facet_grid(rain~plasticity) + geom_line(aes(colour = as.factor(tau)))+
  theme_mountgrass() + scale_colour_mountgrass()


# About phenotypic variations:

# Does tau (0 or 1) lead to more variaitons ? (probably)



# is range and sd varying with the algorithm ? (and rain ?)

# plotboxplot sd (or range) function of algo (facet x), tau (x) (ignore non plastic) for the 3 rain (facet)

showtime %>% filter(plasticity != "none") %>% group_by(rain, plasticity, tau, par_v) %>% summarise(sdRMF = sd(RMF), rangeRMF = (max(RMF) - min(RMF))) %>%
  ggplot(aes(plasticity, sdRMF)) + geom_boxplot(aes(group = interaction(tau, plasticity, rain), colour = plasticity)) + geom_jitter(aes(colour = plasticity), alpha = 0.2) +
  facet_grid(rain~tau)+
  theme_mountgrass() + scale_colour_mountgrass()


showtime %>% filter(plasticity != "none", tau == 0) %>% group_by(rain, plasticity, tau, par_v) %>% summarise(sdRMF = sd(RMF), rangeRMF = (max(RMF) - min(RMF))) %>%
  melt(measure.vars = c("sdRMF", "rangeRMF")) %>%
  ggplot(aes(plasticity, value)) + geom_boxplot(aes(group = interaction(tau, plasticity, rain), colour = plasticity)) + geom_jitter(aes(colour = plasticity), alpha = 0.2) +
  facet_wrap(~variable, scales = "free_y")+
  theme_mountgrass() + scale_colour_mountgrass()

showtime %>% filter(plasticity != "none", tau == 1) %>% group_by(rain, plasticity, tau, par_v) %>% summarise(sdRMF = sd(RMF), rangeRMF = (max(RMF) - min(RMF))) %>%
  group_by(plasticity) %>% summarise(sdMRF = median(sdRMF), rangeRMF = median(rangeRMF)) %>% print()

# Effect of pl on perf:
showtime %>% filter(tau == 0, age == 100) %>% group_by(par_v, rain) %>% mutate(BM_tot_rel = BM_tot/min(BM_tot)) %>%
  ggplot(aes(plasticity, BM_tot_rel)) + geom_boxplot(aes(colour = plasticity))+
  facet_grid(.~rain)+
  theme_mountgrass() + scale_colour_mountgrass()


# Effect of pl on perf:
showtime %>% filter(tau == 0, age == 100) %>% group_by(par_v, rain) %>% dcast(par_v+rain~plasticity, value.var = "BM_tot") %>% 
  mutate(diff = (plastic - none)/none, diff2 = (fixed - none)/none) %>% summary()
