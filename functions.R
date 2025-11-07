## LI Xingyao (s2502246)
# Place your function definitions that may be needed in the report.Rmd, including function documentation.
# You can also include any needed library() calls here

data(ghcnd_stations, package = "StatCompLab")
data(ghcnd_values, package = "StatCompLab")
library(tidyverse)

# Define a function to divide seasons into summer and winter
add_season <- function(as) {
  as$Season <- case_when(
    as$Month %in% c(1, 2, 3, 10, 11, 12) ~ "Winter",
    as$Month %in% c(4, 5, 6, 7, 8, 9) ~ "Summer",
  )
  as
}
 
# Add season to ghcnd
ghcnd <- left_join(ghcnd_values, ghcnd_stations, by = "ID")%>%
  add_season()
head(ghcnd)

# Seasonal variability for temperature
ghcnd %>%
  filter(Element %in% c("TMIN", "TMAX")) %>%
  filter(Year == "2002") %>%
  group_by(Name, Element, Month, Season) %>%
  ggplot(aes(DecYear, Value, colour = Element)) +
  geom_point() +
  facet_wrap(~ Name)

# Seasonal variability for precipitation
ghcnd %>%
  filter(Element %in% c("PRCP")) %>%
  filter(Year == "2002") %>%
  group_by(Name, Element, Month, Season) %>%
  ggplot(aes(DecYear, Value, colour = Element)) +
  geom_point() +
  facet_wrap(~ Name)

## 2
# Add summer column to ghcnd_values
summer_months <- c(4, 5, 6, 7, 8, 9)
ghcnd_values <- ghcnd_values %>%
  mutate(Summer = ifelse(Month %in% summer_months, TRUE, FALSE)) %>%
  # The Name column has been added here for future use
  left_join(ghcnd_stations %>% select(ID, Name), by = "ID")


# Calculate the average summer PRCP for each station
summer_av <- ghcnd_values %>%
  filter(Element == "PRCP") %>%
  filter(Summer == TRUE) %>%
  group_by(ID) %>% 
  summarise(Value = mean(Value))

# Calculate the average winter PRCP for each station
winter_av <- ghcnd_values %>%
  filter(Element == "PRCP") %>%
  filter(Summer == FALSE) %>%
  group_by(ID) %>% 
  summarise(Value = mean(Value))

# Compute observed (true) test statistic 
obs_stat <- abs(summer_av$Value - winter_av$Value)


## Below is the computation of the Monte Carlo permutation test 
p_value_CI <- function(data, alpha=0.05, n=1000){
  count <- rep(0, 8)
  ci_lower <- vector("numeric", 8)
  ci_upper <- vector("numeric", 8)
  se <- rep(0, 8)
  
  for (i in 1:n){
    # Generate a random permutation of summer indicator variable
    summer_perm <- sample(data$Summer, length(data$Summer), replace = FALSE)
    # Apply the permutation to the data frame
    ghcnd_perm <- data %>%
      mutate(Summer = summer_perm) %>%
      group_by(ID)
    # Calculate the null distribution 
    null_dist <- ghcnd_perm %>%
      summarise(winter_avg = mean(Value[!Summer], na.rm = TRUE),
                summer_avg = mean(Value[Summer], na.rm = TRUE), .groups="drop") %>%
      # The difference in means between summer and winter temperatures
      mutate(diff = abs(summer_avg - winter_avg))
    count <- count + (null_dist$diff >= obs_stat)
  }
  p_values <- count/n
  
  # Compute the standard error
  se <- sqrt(p_values*(1-p_values)/n)
  z_alpha <- qnorm(1-alpha/2)
  # Compute confidence interval
  for (j in 1:8){
    if(p_values[j] > 0) {
      ci_lower[j] <- p_values[j] - z_alpha*se[j]
      ci_upper[j] <- p_values[j] + z_alpha*se[j]
    } else {
      ci_lower[j] <- 0
      ci_upper[j] <- 1 - 0.025^(1/n)
    }
  }
  # Return a table with Name, p_value, and confidence interval
  result <- data.frame(Name = unique(ghcnd_values$Name), 
                       p_value = p_values, 
                       CI_lower = ci_lower, 
                       CI_upper = ci_upper)
  return(result)
}

  
# Get the p_value and confidence interval for each station
p_values <- p_value_CI(data = ghcnd_values, alpha = 0.05, n = 1000)
p_values




## 3
# Define a function to get Value_sqrt_avg
ghcnd_values_sqrt <- ghcnd_values %>%
  filter(Element == "PRCP") %>%
  group_by(ID, Year, Month) %>%
  summarise(Value_sqrt_avg = sqrt(mean(Value)), DecYear = mean(DecYear), .groups="drop") %>%
  left_join(ghcnd_stations %>% select(ID, Name, Longitude, Latitude, Elevation), by = "ID")
  

# M0 model
M0 <- lm(Value_sqrt_avg ~ Longitude + Latitude + Elevation + DecYear, data = ghcnd_values_sqrt)

# M1 model
M1 <- lm(Value_sqrt_avg ~ Longitude + Latitude + Elevation + DecYear +
           cos(2*pi*DecYear) + sin(2*pi*DecYear), data = ghcnd_values_sqrt)

# M2 model
M2 <- lm(Value_sqrt_avg ~ Longitude + Latitude + Elevation + DecYear +
           cos(2*pi*DecYear) + sin(2*pi*DecYear) +
           cos(4*pi*DecYear) + sin(4*pi*DecYear), data = ghcnd_values_sqrt)

# M3 model
M3 <- lm(Value_sqrt_avg ~ Longitude + Latitude + Elevation + DecYear +
           cos(2*pi*DecYear) + sin(2*pi*DecYear) +
           cos(4*pi*DecYear) + sin(4*pi*DecYear) +
           cos(6*pi*DecYear) + sin(6*pi*DecYear), data = ghcnd_values_sqrt)

# M4 model
M4 <- lm(Value_sqrt_avg ~ Longitude + Latitude + Elevation + DecYear +
           cos(2*pi*DecYear) + sin(2*pi*DecYear) +
           cos(4*pi*DecYear) + sin(4*pi*DecYear) +
           cos(6*pi*DecYear) + sin(6*pi*DecYear) +
           cos(8*pi*DecYear) + sin(8*pi*DecYear), data = ghcnd_values_sqrt)
# Make the summary
summary(M0)
summary(M1)
summary(M2)
summary(M3)
summary(M4)
# Print the results of the coefficients' estimation
coef(M0)
coef(M1)
coef(M2)
coef(M3)
coef(M4)


##4
library(tidyverse)
library(StatCompLab)
# Define a function to fit the models
fit_models <- function(K, data){
  data <- data %>%
    mutate(Month = factor(Month, levels = 1:12, labels = month.abb))
  
  # Create the formula 
  if (K == 0) {
    f <- as.formula("Value_sqrt_avg ~ Longitude + Latitude + Elevation + DecYear")
  } else if (K > 0 & K <= 4) {
    f <- as.formula(paste0("Value_sqrt_avg ~ Longitude + Latitude + Elevation + DecYear + ",
                           paste0("cos(2*pi*", 1:K, "*DecYear) + sin(2*pi*", 1:K, "*DecYear)", collapse = " + ")))
  } else {
    stop("Invalid value of K")
  }
  
  # Fit the model
  fit <- lm(f, data = data)
  return(fit)
}

# Cross Validation for stations
cv_station <- function(K, score_type){
  # create an empty data frame to store the results
  ave_score <- data.frame(Name = integer(),
                          Year = integer(),
                          Month = integer(),
                          pred_score = numeric())
  for (i in unique(ghcnd_values_sqrt$Name)){
    prcp_station <- subset(ghcnd_values_sqrt, Name == i)
    stats <- subset(ghcnd_values_sqrt, Name != i)
    # Fit the model
    fit <- fit_models(K, data = stats)
    # Subset data for the year and month
    for (j in unique(prcp_station$Year)){
      prcp_year <- prcp_station %>% filter(Year == j)
      for (k in unique(prcp_year$Month)){
        prcp_data <- prcp_year %>% filter(Month == k)
        # Predict precipitation using the model fit
        pred <- predict(fit, newdata = prcp_data, se.fit = TRUE)
        pred_score <- score_prediction(fit, pred, prcp_data, K, score_type)
        # Check if pred_score is not NA
        if(!is.na(pred_score)){    
          ave_score <- rbind(ave_score, data.frame(Name = i,
                                                   Year = j,
                                                   Month = k,
                                                   pred_score = pred_score))
        }
      }
    }
  }
  # Calculate the average score for each station and summarise
  Score <- ave_score %>%
    group_by(Name, Month) %>%
    summarise(Name=Name, Month=Month, pred_score = mean(pred_score, na.rm = TRUE), .groups = "drop") %>% 
    distinct()
  return(Score)
}

# Define a function to calculate the SE and DS scores 
score_prediction <- function(fit, pred, prcp_data, K, score_type){
  # Calculate residual variance
  residual_variance <- sum(fit$residuals^2)/fit$df.residual
  # Calculate standard deviation of prediction
  sd_pred <- sqrt(pred$se.fit^2 + residual_variance)
  # Calculate SE score
  if (score_type=="se"){
    score <- proper_score(type = "se", obs = prcp_data$Value_sqrt_avg,
                          mean = pred$fit)
  }
  # Calculate DS score
  else if (score_type=="ds"){
    score <- proper_score(type = "ds", obs = prcp_data$Value_sqrt_avg,
                          mean = pred$fit, sd = sd_pred)
  }
  return(score)
}

# Run cv_station function to create ave_score object
score_M0_se <- cv_station(K = 0, score_type = "se")
score_M1_se <- cv_station(K = 1, score_type = "se")
score_M2_se <- cv_station(K = 2, score_type = "se")
score_M3_se <- cv_station(K = 3, score_type = "se")
score_M4_se <- cv_station(K = 4, score_type = "se")

# Calculate the average pred_score for each station across all months
mean_score_M0 <- score_M0_se %>%
  group_by(Name) %>%
  summarize(mean_pred_score = mean(pred_score, na.rm = TRUE)) %>%
  arrange(mean_pred_score)

mean_score_M1 <- score_M1_se %>%
  group_by(Name) %>%
  summarize(mean_pred_score = mean(pred_score, na.rm = TRUE)) %>%
  arrange(mean_pred_score)

mean_score_M2 <- score_M2_se %>%
  group_by(Name) %>%
  summarize(mean_pred_score = mean(pred_score, na.rm = TRUE)) %>%
  arrange(mean_pred_score)

mean_score_M3 <- score_M3_se %>%
  group_by(Name) %>%
  summarize(mean_pred_score = mean(pred_score, na.rm = TRUE)) %>%
  arrange(mean_pred_score)

mean_score_M4 <- score_M4_se %>%
  group_by(Name) %>%
  summarize(mean_pred_score = mean(pred_score, na.rm = TRUE)) %>%
  arrange(mean_pred_score)

df_se <- data.frame(
  Name = mean_score_M0$Name,
  scores_M0 = mean_score_M0$mean_pred_score,
  scores_M1 = mean_score_M1$mean_pred_score,
  scores_M2 = mean_score_M2$mean_pred_score,
  scores_M3 = mean_score_M3$mean_pred_score,
  scores_M4 = mean_score_M4$mean_pred_score
)

library(reshape2)
# reshape the data into long format
df_long_se <- reshape2::melt(df_se, id.vars = "Name", variable.name = "Model", value.name = "MeanScore")

library(ggplot2)
# create a bar chart for se scores
ggplot(df_long_se, aes(x = Name, y = MeanScore, fill = Model)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "", y = "Mean Score", fill = "Model") +
  ggtitle("Mean Scores by Name and Model for SE") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


# Run cv_station function to create ave_score objects for each month
score_M0_ds <- cv_station(K = 0, score_type = "ds")
score_M1_ds <- cv_station(K = 1, score_type = "ds")
score_M2_ds <- cv_station(K = 2, score_type = "ds")
score_M3_ds <- cv_station(K = 3, score_type = "ds")
score_M4_ds <- cv_station(K = 4, score_type = "ds")

# Calculate the average pred_score for each station across all months
mean_score1_M0 <- score_M0_ds%>%
  group_by(Name) %>%
  summarize(mean_pred_score = mean(pred_score, na.rm = TRUE)) %>%
  arrange(mean_pred_score)

mean_score1_M1 <- score_M1_ds %>%
  group_by(Name) %>%
  summarize(mean_pred_score = mean(pred_score, na.rm = TRUE)) %>%
  arrange(mean_pred_score)

mean_score1_M2 <- score_M2_ds %>%
  group_by(Name) %>%
  summarize(mean_pred_score = mean(pred_score, na.rm = TRUE)) %>%
  arrange(mean_pred_score)

mean_score1_M3 <- score_M3_ds %>%
  group_by(Name) %>%
  summarize(mean_pred_score = mean(pred_score, na.rm = TRUE)) %>%
  arrange(mean_pred_score)

mean_score1_M4 <- score_M4_ds %>%
  group_by(Name) %>%
  summarize(mean_pred_score = mean(pred_score, na.rm = TRUE)) %>%
  arrange(mean_pred_score)


df_ds <- data.frame(
  Name = mean_score1_M0$Name,
  scores_M0 = mean_score1_M0$mean_pred_score,
  scores_M1 = mean_score1_M1$mean_pred_score,
  scores_M2 = mean_score1_M2$mean_pred_score,
  scores_M3 = mean_score1_M3$mean_pred_score,
  scores_M4 = mean_score1_M4$mean_pred_score
)

# reshape the data into long format
df_long_ds <- reshape2::melt(df_ds, id.vars = "Name", variable.name = "Model", value.name = "MeanScore")

# create a bar chart for ds scores
ggplot(df_long_ds, aes(x = Name, y = MeanScore, fill = Model)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "", y = "Mean Score", fill = "Model") +
  ggtitle("Mean Scores by Name and Model for DS") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

cv_station(K = 0, score_type = "se")
cv_station(K = 1, score_type = "se")
cv_station(K = 2, score_type = "se")
cv_station(K = 3, score_type = "se")
cv_station(K = 4, score_type = "se")
cv_station(K = 0, score_type = "ds")
cv_station(K = 1, score_type = "ds")
cv_station(K = 2, score_type = "ds")
cv_station(K = 3, score_type = "ds")
cv_station(K = 4, score_type = "ds")
