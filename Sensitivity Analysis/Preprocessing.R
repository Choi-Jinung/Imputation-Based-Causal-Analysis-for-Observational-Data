getwd()
setwd(getwd())
rm(list=ls())
library(readxl)
library(dplyr)

# Import 'indiegogo' and 'Crowdfunding' data
indiegogo_foodtravel <- read_excel("indiegogo_foodtravel.xlsx")
Crowdfunding <- read_excel("Crowdfunding proposal- F&B Jan 2024_ HJOO_Combined.xlsx")

# Merge the two data by 'project_id'
raw.data <- merge(indiegogo_foodtravel, Crowdfunding, by="project_id")

# Make 'duration' variable by counting the dates
raw.data$open_date <- as.Date(raw.data$open_date, format="%Y-%m-%d")
raw.data$close_date <- as.Date(raw.data$close_date, format="%Y-%m-%d")
raw.data$duration <- as.numeric(raw.data$close_date - raw.data$open_date)
raw.data <- raw.data %>% select(fund_USD, duration, visuals, videos, word_count, Classification)

# Log-transformation scailing
scailing <- function(x) log(x+1) 
raw.data <- raw.data %>% mutate(fund_USD=scailing(fund_USD), duration=scailing(duration), 
                                visuals=scailing(visuals), videos=scailing(videos), word_count=scailing(word_count))
raw.data <- raw.data %>% filter(fund_USD>0)
data <- raw.data %>% mutate(y=fund_USD, x1=duration, x2=visuals, x3=videos, x4=word_count, t=Classification)

# t=='user' -> treatment group, t=='self' -> control group
data <- data %>% mutate(t=case_when(t=="user" ~ 1,
                                    t=="self" ~ 0,
                                    TRUE ~ NA)) %>% mutate(t=as.factor(t))
data <- data %>% filter(is.na(t)==F)
data <- data %>% mutate(x1_2=x1^2, x2_2=x2^2, x3_2=x3^2, x4_2=x4^2,
                        x12=x1*x2, x13=x1*x3, x14=x1*x4, x23=x2*x3, x24=x2*x4, x34=x3*x4)
data <- data %>% select(y, x1, x2, x3, x4, x1_2, x2_2, x3_2, x4_2, x12, x13, x14, x23, x24, x34, t)

# Control group data
data_t0 <- data %>% filter(t==0) %>% rename(y0=y); data_t0$t <- NULL
write.csv(data_t0, "data_t0.csv")

# Treatment group data
data_t1 <- data %>% filter(t==1) %>% rename(y1=y); data_t1$t <- NULL
write.csv(data_t1, "data_t1.csv")
