library("modelr")
library(smotefamily)
library(imbalance)
CA_case <- read.csv("D:\\研究所Meeting\\COVID19\\CA government data\\2022_4_18 CA case.csv")[-59,]
#Confirmed rate(單日確診率)
#CA_case$Class <-rep(as.character(c('positive', 'negative')), times=c(20, 38))
CA_case$Class <-rep(as.character(c('positive', 'negative')), times=29)
#CA_case$Class <-rep(as.character(c('Yes')), times=c(58))

CA_case$Confirmed_rate <- as.numeric(CA_case$reported_cases/CA_case$population)
CA_case <- as.data.frame(CA_case[,c(3,4,16,17)])

#new_data <- mwmote(data=CA_case, numInstances = 500,cclustering = 3,classAttr = "Class")

new_data <- pdfos(data=CA_case, numInstances = 500,classAttr = "Class")
#new_data <- neater(data=CA_case, numInstances = 100,k=3,classAttr = "Class")


plot(CA_case[,1:2],xlim=c(-125,-114),ylim=c(32,42))
plot(new_data[,1:2],xlim=c(-125,-114),ylim=c(32,42))
plot(new_data[,1:2])


summary(CA_case[,1:2])
summary(new_data[,1:2])
