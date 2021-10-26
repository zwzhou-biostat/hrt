
library(rstan)
library(dplyr)
library(ggplot2)
library(tidyr)
library(matrixStats)
library(gridExtra)
library(parallel)
library(data.table)
library(gridExtra)

source(".\\code\\adj_inc_EpiEstim.R")

MA_cases_count <- readRDS(".\\data\\MA_all_counties_cases_count.rds")
avg_mobility_matrix_list <- readRDS(".\\data\\avg_mobility_matrix_ma_all_counties_list.rds")
geoid_link <- fread(".\\data\\all-geocodes-v2017.csv")

p <- sqrt(length(avg_mobility_matrix_list[[1]]))

geoid_link_mass <- geoid_link[geoid_link$`State Code (FIPS)` == 25, ]
ma_fips <- seq(25001, 25027, 2)

geoid_link_mass_unique <- geoid_link_mass[!duplicated(geoid_link_mass$`County Code (FIPS)`), ] %>% as.data.frame()
county_names <- geoid_link_mass_unique$`Area Name (including legal/statistical area description)`[-1]
county_names <- gsub(" County", "", county_names)

geoids <- ma_fips

case_count <- data.frame(date = rep(seq(as.Date("2020-07-01"), as.Date("2021-01-28"), by = "day"), each = p), geoid = rep(geoids, length(seq(as.Date("2020-07-01"), as.Date("2021-01-28"), by = "day"))))

MA_cases_count$cdc_report_dt <- as.Date(MA_cases_count$cdc_report_dt)

case_count <- dplyr::left_join(case_count, MA_cases_count, by = c("date" = "cdc_report_dt", "geoid" = "county_fips_code"))

case_count$n[is.na(case_count$n)] <- 0

case_count$county <- factor(case_count$geoid, levels = geoids, labels = county_names)

MA_cases_count_wide <- case_count %>% select(., date, geoid, n) %>% spread(., key = geoid, value = n)

MA_cases_count_wide_sub <- MA_cases_count_wide[MA_cases_count_wide$date >= "2020-07-01", ]

case_matrix <- as.matrix(MA_cases_count_wide_sub[, -1])

#### INPUT
start_day <- 12 ## day for start estimating
window_length <- 7 ## smoothing window for estimaiton
si_shape <- 3.45 ## shape parameter for serical interval assuming gamma distribution
si_rate <- 0.66 ## rate parameter for serical interval assuming gamma distribution
si_t_max <- 14 ## maximum number of days with non-zero probability for serial interval
P <- lapply(avg_mobility_matrix_list[1:length(seq(as.Date("2020-07-01"), as.Date("2021-01-28"), by = "day"))], function(x){matrix(x, p)})

#### OUTPUT
## the direct output from estimate_R_adj_inc keeps all outputs the same as EpiEsim, but in a list format, with each element in the list correspond to one region

## obtain result from main function
r_scale_result <- estimate_R_adj_inc(case_matrix, P, si_shape, si_rate, si_t_max, start_day, window_length)


geoid_link <- fread("/rprojectnb2/hrtgrp/data/all-geocodes-v2017.csv")
geoid_link_mass <- geoid_link[geoid_link$`State Code (FIPS)` == 25, ]
ma_fips <- seq(25001, 25027, 2)

geoid_link_mass_unique <- geoid_link_mass[!duplicated(geoid_link_mass$`County Code (FIPS)`), ] %>% as.data.frame()
county_names <- geoid_link_mass_unique$`Area Name (including legal/statistical area description)`[-1]
county_names <- gsub(" County", "", county_names)




dat_regions <- do.call("rbind", lapply(1:ncol(case_matrix), function(x){
  dat <- as.data.frame(r_scale_result[[x]]$R[, c(3, 5, 11)])
  dat$region <- county_names[x]
  dat$day <- 1:nrow(dat)
  dat
}))

library(tidyr)
dat_plot <- group_by(dat_regions, region, day) %>% summarise(r = mean(`Mean(R)`), rl = mean(`Quantile.0.025(R)`), ru = mean(`Quantile.0.975(R)`))


included_days <- seq(as.Date("2020-07-01"), as.Date("2021-01-28"), by = "day")

dat_plot$day <- rep(included_days[start_day:(dim(case_matrix)[1] - window_length)], p)


g2 <- ggplot(dat_plot, aes(x = day)) +
  # geom_smooth(data,mapping=aes(y=Rt,x=x, color = County),size=1.2)+lims(y = c(0, 5)) +
  geom_line(dat_plot,mapping=aes(y=r,x=day),size=0.6)+#lims(y = c(0, 3)) +
  geom_ribbon(dat_plot,mapping=aes(x=day,ymin=rl,ymax=ru),alpha=0.4) +
  geom_hline(yintercept = 1,color = "black", size=0.6,alpha=0.5)  +
  theme_classic() +  xlab("Days") + ylab('Reproduction Number')  +
  facet_wrap( ~ region, nrow = 2, scales = "free")  # + scale_y_continuous(trans = "log10")


dat_regions_R_est <- spread(dat_regions[, c(1, 4, 5)], region, `Mean(R)`)[, -1]
dat_regions_R_lb_est <- spread(dat_regions[, c(2, 4, 5)], region, `Quantile.0.025(R)`)[, -1]
dat_regions_R_ub_est <- spread(dat_regions[, c(3, 4, 5)], region, `Quantile.0.975(R)`)[, -1]
dat_regions_inc <- as.matrix(case_matrix)


estimatated_incidence <- sapply((1:nrow(dat_regions_R_est))+6, function(t){
  rev(w[1:min(length(w), t)])%*%dat_regions_inc[c(max(t-length(w) + 1, 1): t), ,drop = F]%*%diag(dat_regions_R_est[t, ])%*%P[[t]]
})


estimatated_incidence_lb <- sapply((1:nrow(dat_regions_R_est))+6, function(t){
  rev(w[1:min(length(w), t)])%*%dat_regions_inc[c(max(t-length(w) + 1, 1): t), ,drop = F]%*%diag(dat_regions_R_lb_est[t, ])%*%P[[t]]
})


estimatated_incidence_ub <- sapply((1:nrow(dat_regions_R_est))+6, function(t){
  rev(w[1:min(length(w), t)])%*%dat_regions_inc[c(max(t-length(w) + 1, 1): t), ,drop = F]%*%diag(dat_regions_R_ub_est[t, ])%*%P[[t]]
})


included_days <- seq(as.Date("2020-07-01"), as.Date("2021-01-28"), by = "day")

estimatated_incidence <- as.data.frame(t(estimatated_incidence))


names(estimatated_incidence) <- county_names

estimatated_incidence_long <- gather(estimatated_incidence)
estimatated_incidence_long$day <- rep(included_days[start_day:(dim(case_matrix)[1] - window_length)], p)

estimatated_incidence_long$yl <- gather(as.data.frame(t(estimatated_incidence_lb)))$value

estimatated_incidence_long$yh <- gather(as.data.frame(t(estimatated_incidence_ub)))$value

g1 <- ggplot(estimatated_incidence_long, aes(x = day)) +
  # geom_smooth(data,mapping=aes(y=y,x=x, color = County),size=1.2) +
  geom_line(estimatated_incidence_long,mapping=aes(y=value,x=day),size=0.6) +
  # geom_line(data = data, aes(y = y_real, color = County)) +
  geom_ribbon(estimatated_incidence_long,mapping=aes(x=day,ymin=yl,ymax=yh),alpha=0.3) +
  facet_wrap( ~ key, nrow = 2, scales = "free")  +
  theme_classic() +  xlab("Days") + ylab('Infections')


jpeg(".\\output\\ma_all_counties_scale_inc.jpg", res = 300, w = 5000, h = 2500)
grid.arrange(g1,g2)
dev.off()
