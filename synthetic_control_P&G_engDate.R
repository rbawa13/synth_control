library(dplyr)
library(data.table)
library(reshape2)
library(CausalImpact)
options(scipen = 9999)


###
setwd("~/Downloads/")

file1 <- fread('2024-08-21 2_25pm.csv', header = T, sep = ",")
file1 <- as.data.frame(file1)
file1$ENTRY_TIMESTAMP <- as.Date(file1$ENTRY_TIMESTAMP, format = "%Y-%m-%d")
str(file1)
file1 <- na.omit(file1)
file1$unit <- "One"
file1$unit <- as.factor(as.character(file1$unit))


####
file2 <- fread('Engagement_date_data1.csv', header = T, sep = ",")
str(file2)
file2$DAY <- as.Date(file2$DAY, format = "%Y-%m-%d")
file2$unit <- "One"
file2$unit <- as.factor(as.character(file2$unit))

####functions
select_cols <- function(filename, date_col, unit_col, chosen_metric)
{
  file <- dplyr::select(filename, date_col, unit_col, chosen_metric)
  return(file)
}

exp_temp_plot <- function(filename, date_col, unit_col, chosen_metric)
{
  file <- select_cols(filename, date_col, unit_col, chosen_metric)
  meltdf <- melt(file,id=c(date_col, unit_col), variable.name = chosen_metric)
  names(meltdf) <- c("Date", "Unit", "Variable", "Value")
  return(meltdf)
}

exp_plot <- function(filename, date_col, unit_col, chosen_metric, intvn_date)
{
  file <- exp_temp_plot(filename, date_col, unit_col, chosen_metric)
  p1 <- ggplot(file,aes_string(x="Date",y="Value",colour="Unit",group="Unit"))
  p1 <- p1 + geom_line(size = 0.9, alpha = 0.7) + scale_color_brewer(palette = "Set1")
  p1 <- p1 + geom_vline(xintercept=as.numeric(as.Date(intvn_date)), linetype=2)
  p1 <- p1 + theme_classic() + theme(legend.position="bottom")
  return(p1)
}

exp_plotly <- function(filename, date_col, unit_col, chosen_metric, intvn_date)
{
  p1 <- exp_plot(filename, date_col, unit_col, chosen_metric, intvn_date)
  p1 <- ggplotly(p1)
  p1 <- p1 
  return(p1)
}


CI_format <- function(filename, date_col, unit_col, chosen_metric, treatment_unit)
{
  file <- select_cols(filename, date_col, unit_col, chosen_metric)
  wide_data <-  dcast(file, paste(paste(date_col, collapse = " + "), "~", paste(unit_col, collapse = " + ")), 
                      value.var = chosen_metric, fun.aggregate = mean)
  wide_data[is.na(wide_data)] <- 0
  time <- wide_data[,1]
  wide_data <- wide_data[,-1]
  test1 <- as.zoo(xts(wide_data, order.by=time))
  names(test1) <- chosen_metric
  return(test1)
}

synth_control2 <- function(input_file, trt_start_date, data_type="daily", seasons = 52)
{
  
  #input_file = input file from CI_format function
  #trt_start_date = Date when the treatment started
  #data_type = What time level the data is (weekly, daily, hourly or monthly)
  
  set.seed(1234)
  time_df <- time(input_file)
  time_df <- sort(time_df)
  min_time <- min(time_df)
  max_time <- max(time_df)
  trt_start <- as.Date(trt_start_date, format = "%Y-%m-%d")
  pre_end <- trt_start - 1
  
  pre_period = c(min_time, pre_end)
  post_period = c(trt_start, max_time)
  
  unit_names =  names(input_file)
  unit_names1 = gsub(" ", "_", unit_names, fixed = TRUE)
  
  names(input_file) = unit_names1
  
  if(data_type == "weekly")
    model = CausalImpact(input_file, pre_period, post_period, 
                         model.args = list(niter = 1000, nseasons = 7, season.duration = 7, dynamic.regression = F))
  else if(data_type == "week of the year")
    model = CausalImpact(input_file, pre_period, post_period, 
                         model.args = list(niter = 1000, nseasons = 52))
  else if(data_type == "daily" | data_type == "day of the week")
    model = CausalImpact(input_file, pre_period, post_period, 
                         model.args = list(niter = 1000, nseasons = 7, dynamic.regression = F))
  else if(data_type == "day of the year")
    model = CausalImpact(input_file, pre_period, post_period, 
                         model.args = list(niter = 1000, nseasons = 200))
  else if(data_type == "monthly")
    model = CausalImpact(input_file, pre_period, post_period, 
                         model.args = list(niter = 1000, nseasons = 12))
  else if(data_type == "quarterly")
    model = CausalImpact(input_file, pre_period, post_period, 
                         model.args = list(niter = 1000, nseasons = 4))
  else if(data_type == "custom" | data_type == "multiple seasons")
    model = CausalImpact(input_file, pre_period, post_period, 
                         model.args = list(niter = 1000, nseasons = seasons))
  
  summary_df <- data.frame(model$summary)
  abs_diff <- round(summary_df$AbsEffect[1], 4)
  abs_diff_U <- round(summary_df$AbsEffect.upper[1], 4)
  abs_diff_L <- round(summary_df$AbsEffect.lower[1], 4)
  
  cum_diff <- round(summary_df$AbsEffect[2], 4)
  cum_diff_U <- round(summary_df$AbsEffect.upper[2], 4)
  cum_diff_L <- round(summary_df$AbsEffect.lower[2], 4)
  
  pval <- round(summary_df$p[1], 4)
  print(paste("Absolute effect of the treatment is: ", abs_diff, sep = ""))
  print(paste("Credible intervals on the absolute effect of the treatment are: ", abs_diff_L, "  -  ", abs_diff_U, sep = ""))
  
  print(paste("Cumulative effect of the treatment is: ", cum_diff, sep = ""))
  print(paste("Credible intervals on the cumulative effect of the treatment are: ", 
              cum_diff_L, "  -  ", cum_diff_U, sep = ""))
  
  print(paste("The p-value on the treatment effect is: ", pval, sep = ""))
  
  
  plot(model)
  return(model)
  
}

mean_diff <- function(pre, post) {
  mean(post) - mean(pre)
}

# Bootstrap function
bootstrap <- function(pre, post, n_bootstrap = 1000, alpha = 0.05) {
  bootstrap_diffs <- numeric(n_bootstrap)
  
  for (i in 1:n_bootstrap) {
    # Resample with replacement
    resample_pre <- sample(pre, size = length(pre), replace = TRUE)
    resample_post <- sample(post, size = length(post), replace = TRUE)
    
    # Calculate the difference in means for this resample
    bootstrap_diffs[i] <- mean_diff(resample_pre, resample_post)
  }
  
  # Mean of the bootstrap differences
  mean_diff_bootstrap <- mean(bootstrap_diffs)
  
  # Confidence interval
  lower_bound <- quantile(bootstrap_diffs, alpha / 2)
  upper_bound <- quantile(bootstrap_diffs, 1 - alpha / 2)
  
  # P-value calculation using two-tailed test
  observed_diff <- mean_diff(pre, post)
  p_value <- mean(abs(bootstrap_diffs) >= abs(observed_diff))
  
  list(mean_diff_bootstrap = mean_diff_bootstrap, 
       confidence_interval = c(lower_bound, upper_bound), 
       p_value = p_value, 
       bootstrap_diffs = bootstrap_diffs)
}

####exploratory plots
tmp_plot1 <- exp_plotly(file2, date_col = 'DAY', unit_col = 'unit', 
                      chosen_metric = 'CNT_STARTED_ENG_DATE', intvn_date = '2024-07-15')
tmp_plot1

#####
sc_input1 <- CI_format(filename = file2, date_col = 'DAY', unit_col = 'unit', 
                       chosen_metric = 'CNT_STARTED_ENG_DATE', treatment_unit = "One")


###find the right seasonality for the dataset
pgram <- spec.pgram(sc_input1$CNT_STARTED_ENG_DATE, log="no", plot = T, spans = c(2,2))
pgram_df <- data.frame(Freq = pgram$freq, 
                       Spectrum = pgram$spec)
pgram_df <- pgram_df[order(-pgram_df$Spectrum),]
top_pgram <- head(pgram_df, 2)
top_pgram
pgram_max_order1 = round(1/top_pgram$Freq[1], 0)
pgram_max_order2 = round(1/top_pgram$Freq[2], 0)
print(pgram_max_order1) ### major seasonality1
print(pgram_max_order2) ### major seasonality2

seasons = min(pgram_max_order1, pgram_max_order2)

###run synthetic control
sc_res1 <- synth_control2(input_file = sc_input1, trt_start_date = "2024-07-15", 
                          data_type = "custom", seasons = seasons)
plot(sc_res1)
summary(sc_res1)

###### D7 Sessions
plot_d7 <- exp_plotly(file2, date_col = 'DAY', unit_col = 'unit', 
                        chosen_metric = 'AVG_D7_SESSIONS', intvn_date = '2024-07-15')
plot_d7

#####
sc_input_d7 <- CI_format(filename = file2, date_col = 'DAY', unit_col = 'unit', 
                       chosen_metric = 'AVG_D7_SESSIONS', treatment_unit = "One")

sc_res_d7 <- synth_control2(input_file = sc_input_d7, trt_start_date = "2024-07-15", 
                          data_type = "custom", seasons = seasons)
plot(sc_res_d7)
summary(sc_res_d7)

#####PCT_COMPLETED_ENG_STEP
plot_eng_step <- exp_plotly(file2, date_col = 'DAY', unit_col = 'unit', 
                      chosen_metric = 'PCT_COMPLETED_ENG_STEP', intvn_date = '2024-07-15')
plot_eng_step

#####
sc_input_eng_step <- CI_format(filename = file2, date_col = 'DAY', unit_col = 'unit', 
                         chosen_metric = 'PCT_COMPLETED_ENG_STEP', treatment_unit = "One")

sc_res_eng_step <- synth_control2(input_file = sc_input_eng_step, trt_start_date = "2024-07-15", 
                            data_type = "custom", seasons = seasons)
plot(sc_res_eng_step)
summary(sc_res_eng_step)

#####CNT_PROVIDED_ENG_DATE
plot_pr_eng_date <- exp_plotly(file2, date_col = 'DAY', unit_col = 'unit', 
                            chosen_metric = 'CNT_PROVIDED_ENG_DATE', intvn_date = '2024-07-15')
plot_pr_eng_date

#####
sc_input_pr_eng_date <- CI_format(filename = file2, date_col = 'DAY', unit_col = 'unit', 
                               chosen_metric = 'CNT_PROVIDED_ENG_DATE', treatment_unit = "One")

sc_res_pr_eng_date <- synth_control2(input_file = sc_input_pr_eng_date, trt_start_date = "2024-07-15", 
                                  data_type = "custom", seasons = seasons)
plot(sc_res_pr_eng_date)
summary(sc_res_pr_eng_date)

#####PCT_PROVIDED_ENG_DATE
plot_pct_eng_date <- exp_plotly(file2, date_col = 'DAY', unit_col = 'unit', 
                               chosen_metric = 'PCT_PROVIDED_ENG_DATE', intvn_date = '2024-07-15')
plot_pct_eng_date

#####
sc_input_pct_eng_date <- CI_format(filename = file2, date_col = 'DAY', unit_col = 'unit', 
                                  chosen_metric = 'PCT_PROVIDED_ENG_DATE', treatment_unit = "One")

sc_res_pct_eng_date <- synth_control2(input_file = sc_input_pct_eng_date, trt_start_date = "2024-07-15", 
                                     data_type = "custom", seasons = seasons)
plot(sc_res_pct_eng_date)
summary(sc_res_pct_eng_date)

#####"PCT_COMPLETED_ONBOARDING" 
plot_pct_comp_onboard <- exp_plotly(file2, date_col = 'DAY', unit_col = 'unit', 
                                chosen_metric = 'PCT_COMPLETED_ONBOARDING', intvn_date = '2024-07-15')
plot_pct_comp_onboard

#####
sc_input_pct_comp_onboard <- CI_format(filename = file2, date_col = 'DAY', unit_col = 'unit', 
                                   chosen_metric = 'PCT_COMPLETED_ONBOARDING', treatment_unit = "One")

sc_res_pct_comp_onboard <- synth_control2(input_file = sc_input_pct_comp_onboard, trt_start_date = "2024-07-15", 
                                      data_type = "custom", seasons = seasons)
plot(sc_res_pct_comp_onboard)
summary(sc_res_pct_comp_onboard)

######"PCT_VIEWED_GC_STEP"       
plot_pct_view_gc <- exp_plotly(file2, date_col = 'DAY', unit_col = 'unit', 
                                    chosen_metric = 'PCT_VIEWED_GC_STEP', intvn_date = '2024-07-15')
plot_pct_view_gc

#####
sc_input_pct_view_gc <- CI_format(filename = file2, date_col = 'DAY', unit_col = 'unit', 
                                       chosen_metric = 'PCT_VIEWED_GC_STEP', treatment_unit = "One")

sc_res_pct_view_gc <- synth_control2(input_file = sc_input_pct_view_gc, trt_start_date = "2024-07-15", 
                                          data_type = "custom", seasons = seasons)
plot(sc_res_pct_view_gc)
summary(sc_res_pct_view_gc)


#######"CNT_COMPLETED_ONBOARDING"
plot_cnt_onboard <- exp_plotly(file2, date_col = 'DAY', unit_col = 'unit', 
                               chosen_metric = 'CNT_COMPLETED_ONBOARDING', intvn_date = '2024-07-15')
plot_cnt_onboard

#####
sc_input_cnt_onboard <- CI_format(filename = file2, date_col = 'DAY', unit_col = 'unit', 
                                  chosen_metric = 'CNT_COMPLETED_ONBOARDING', treatment_unit = "One")

sc_res_cnt_onboard <- synth_control2(input_file = sc_input_cnt_onboard, trt_start_date = "2024-07-15", 
                                     data_type = "custom", seasons = seasons)
plot(sc_res_cnt_onboard)
summary(sc_res_cnt_onboard)

######ONBOARDING TIME
####outlier removal
Q <- quantile(file1$AVG_ONBOARDING_TIME, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(file1$AVG_ONBOARDING_TIME)
up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range

file1_1 <- subset(file1, file1$AVG_ONBOARDING_TIME > (Q[1] - 1.5*iqr) & file1$AVG_ONBOARDING_TIME < (Q[2]+1.5*iqr))

pre_data <- file1_1 %>% filter(ENTRY_TIMESTAMP < '2024-07-15')%>% select(AVG_ONBOARDING_TIME)
pre_data1 <- na.omit(pre_data$AVG_ONBOARDING_TIME)
post_data <- file1_1 %>% filter(ENTRY_TIMESTAMP >= '2024-07-15')%>% select(AVG_ONBOARDING_TIME)
post_data1 <- na.omit(post_data$AVG_ONBOARDING_TIME)

####plotdata
pre_data$unit <- "pre"
post_data$unit <- "post"
####
plot_data <- data.frame(rbind(pre_data, post_data))
str(plot_data)
plot_data$unit <- as.factor(plot_data$unit)

p<-ggplot(plot_data, aes(x=AVG_ONBOARDING_TIME, color=unit)) +
  geom_density()
p


####
result <- bootstrap(pre_data1, post_data1)

# Output the results
cat("Bootstrap Mean Difference:", result$mean_diff_bootstrap, "\n")
cat("95% Confidence Interval:", result$confidence_interval, "\n")
cat("P-Value:", result$p_value, "\n")

# Plot the distribution of bootstrap differences
hist(result$bootstrap_diffs, breaks = 50, main = "Bootstrap Distribution of Mean Differences",
     xlab = "Mean Difference", col = "lightblue", border = "white")

# Add mean difference line
abline(v = result$mean_diff_bootstrap, col = "red", lwd = 2, lty = 2)

# Add confidence interval lines
abline(v = result$confidence_interval[1], col = "darkgreen", lwd = 2, lty = 2)
abline(v = result$confidence_interval[2], col = "darkgreen", lwd = 2, lty = 2)

# Add legend
legend("topright", legend = c("Mean Difference", "95% Confidence Interval"), 
       col = c("red", "darkgreen"), lty = 2, lwd = 2)



#####
file4 <- file1_1 %>%
  group_by(ENTRY_TIMESTAMP) %>%
  dplyr::summarize(Mean_metric = mean(AVG_ONBOARDING_TIME, na.rm=TRUE))

file4$unit <- "One"
file4$unit <- as.factor(as.character(file4$unit))

####
plot_onboard_time <- exp_plotly(file4, date_col = 'ENTRY_TIMESTAMP', unit_col = 'unit', 
                               chosen_metric = 'Mean_metric', intvn_date = '2024-07-15')
plot_onboard_time

#####
sc_input_onboard_time <- CI_format(filename = file4, date_col = 'ENTRY_TIMESTAMP', unit_col = 'unit', 
                                  chosen_metric = 'Mean_metric', treatment_unit = "One")

sc_res_onboard_time <- synth_control2(input_file = sc_input_onboard_time, trt_start_date = "2024-07-15", 
                                     data_type = "custom", seasons = seasons)
plot(sc_res_onboard_time)
summary(sc_res_onboard_time)




