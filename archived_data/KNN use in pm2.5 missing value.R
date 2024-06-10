library(DMwR2)
library(magrittr)
library(dplyr) # re column names
library(BBmisc)
CA_pm2.5_4.1 <- read.csv("D:\\研究所Meeting\\COVID19\\PM2.5 in CA\\4_1_2022.csv")
CA_pm2.5_4.2 <- read.csv("D:\\研究所Meeting\\COVID19\\PM2.5 in CA\\4_2_2022.csv")
CA_pm2.5_4.3 <- read.csv("D:\\研究所Meeting\\COVID19\\PM2.5 in CA\\4_3_2022.csv")
CA_pm2.5_4.4 <- read.csv("D:\\研究所Meeting\\COVID19\\PM2.5 in CA\\4_4_2022.csv")
CA_pm2.5_4.5 <- read.csv("D:\\研究所Meeting\\COVID19\\PM2.5 in CA\\4_5_2022.csv")
CA_pm2.5_4.6 <- read.csv("D:\\研究所Meeting\\COVID19\\PM2.5 in CA\\4_6_2022.csv")
CA_pm2.5_4.7 <- read.csv("D:\\研究所Meeting\\COVID19\\PM2.5 in CA\\4_7_2022.csv")
CA_pm2.5_4.8 <- read.csv("D:\\研究所Meeting\\COVID19\\PM2.5 in CA\\4_8_2022.csv")
CA_pm2.5_4.9 <- read.csv("D:\\研究所Meeting\\COVID19\\PM2.5 in CA\\4_9_2022.csv")
CA_pm2.5_4.10 <- read.csv("D:\\研究所Meeting\\COVID19\\PM2.5 in CA\\4_10_2022.csv")
CA_pm2.5_4.11 <- read.csv("D:\\研究所Meeting\\COVID19\\PM2.5 in CA\\4_11_2022.csv")
CA_pm2.5_4.12 <- read.csv("D:\\研究所Meeting\\COVID19\\PM2.5 in CA\\4_12_2022.csv")
CA_pm2.5_4.13 <- read.csv("D:\\研究所Meeting\\COVID19\\PM2.5 in CA\\4_13_2022.csv")
CA_pm2.5_4.14 <- read.csv("D:\\研究所Meeting\\COVID19\\PM2.5 in CA\\4_14_2022.csv")
CA_pm2.5_4.15 <- read.csv("D:\\研究所Meeting\\COVID19\\PM2.5 in CA\\4_15_2022.csv")
CA_pm2.5_4.16 <- read.csv("D:\\研究所Meeting\\COVID19\\PM2.5 in CA\\4_16_2022.csv")
CA_pm2.5_4.1out <- cbind(CA_pm2.5_4.1[,1],knnImputation(CA_pm2.5_4.1[,-1]))  %>%rename(county = "CA_pm2.5_4.1[, 1]")
CA_pm2.5_4.2out <- cbind(CA_pm2.5_4.2[,1],knnImputation(CA_pm2.5_4.2[,-1]))  %>%rename(county = "CA_pm2.5_4.2[, 1]")
CA_pm2.5_4.3out <- cbind(CA_pm2.5_4.3[,1],knnImputation(CA_pm2.5_4.3[,-1]))  %>%rename(county = "CA_pm2.5_4.3[, 1]")
CA_pm2.5_4.4out <- cbind(CA_pm2.5_4.4[,1],knnImputation(CA_pm2.5_4.4[,-1]))  %>%rename(county = "CA_pm2.5_4.4[, 1]")
CA_pm2.5_4.5out <- cbind(CA_pm2.5_4.5[,1],knnImputation(CA_pm2.5_4.5[,-1]))  %>%rename(county = "CA_pm2.5_4.5[, 1]")
CA_pm2.5_4.6out <- cbind(CA_pm2.5_4.6[,1],knnImputation(CA_pm2.5_4.6[,-1]))  %>%rename(county = "CA_pm2.5_4.6[, 1]")
CA_pm2.5_4.7out <- cbind(CA_pm2.5_4.7[,1],knnImputation(CA_pm2.5_4.7[,-1]))  %>%rename(county = "CA_pm2.5_4.7[, 1]")
CA_pm2.5_4.8out <- cbind(CA_pm2.5_4.8[,1],knnImputation(CA_pm2.5_4.8[,-1]))  %>%rename(county = "CA_pm2.5_4.8[, 1]")
CA_pm2.5_4.9out <- cbind(CA_pm2.5_4.9[,1],knnImputation(CA_pm2.5_4.9[,-1]))  %>%rename(county = "CA_pm2.5_4.9[, 1]")
CA_pm2.5_4.10out <- cbind(CA_pm2.5_4.10[,1],knnImputation(CA_pm2.5_4.10[,-1]))  %>%rename(county = "CA_pm2.5_4.10[, 1]")
CA_pm2.5_4.11out <- cbind(CA_pm2.5_4.11[,1],knnImputation(CA_pm2.5_4.11[,-1]))  %>%rename(county = "CA_pm2.5_4.11[, 1]")
CA_pm2.5_4.12out <- cbind(CA_pm2.5_4.12[,1],knnImputation(CA_pm2.5_4.12[,-1]))  %>%rename(county = "CA_pm2.5_4.12[, 1]")
CA_pm2.5_4.13out <- cbind(CA_pm2.5_4.13[,1],knnImputation(CA_pm2.5_4.13[,-1]))  %>%rename(county = "CA_pm2.5_4.13[, 1]")
CA_pm2.5_4.14out <- cbind(CA_pm2.5_4.14[,1],knnImputation(CA_pm2.5_4.14[,-1]))  %>%rename(county = "CA_pm2.5_4.14[, 1]")
CA_pm2.5_4.15out <- cbind(CA_pm2.5_4.15[,1],knnImputation(CA_pm2.5_4.15[,-1]))  %>%rename(county = "CA_pm2.5_4.15[, 1]")
CA_pm2.5_4.16out <- cbind(CA_pm2.5_4.16[,1],knnImputation(CA_pm2.5_4.16[,-1]))  %>%rename(county = "CA_pm2.5_4.16[, 1]")

#cbind all pm2.5 data in one data as one day repeat data 
#normaliz data 
CA_pm2.5_total_out <- cbind(CA_pm2.5_4.1out,CA_pm2.5_4.2out[,4],CA_pm2.5_4.3out[,4],CA_pm2.5_4.4out[,4],CA_pm2.5_4.5out[,4],
                            CA_pm2.5_4.6out[,4],CA_pm2.5_4.7out[,4],CA_pm2.5_4.8out[,4],CA_pm2.5_4.9out[,4],CA_pm2.5_4.10out[,4],
                            CA_pm2.5_4.11out[,4],CA_pm2.5_4.12out[,4],CA_pm2.5_4.13out[,4],CA_pm2.5_4.14out[,4],CA_pm2.5_4.15out[,4],
                            CA_pm2.5_4.16out[,4]) 
write.csv(CA_pm2.5_total_out,file="D:\\研究所Meeting\\COVID19\\PM2.5 in CA\\CA_pm2.5_total_out.csv",row.names = FALSE)

read.csv("D:\\研究所Meeting\\COVID19\\PM2.5 in CA\\CA_pm2.5_total_out.csv")
CA_pm2.5_total_out$several_days <- scale(as.numeric(apply(CA_pm2.5_total_out[,4:19], 1, mean)))
CA_pm2.5_total_KNNout <- CA_pm2.5_total_out[,c(1:3,20)]
