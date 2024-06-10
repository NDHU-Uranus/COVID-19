library(autoFRK)
library(BBmisc)
library(caret)
library(dplyr)
library(dismo)
library(fields)
library(geoR)
library(GGally)
library(ggplot2)
library(gstat)
library(lubridate)
library(mapdata)
library(maps)
library(maptools)
library(mvtnorm)
library(MASS)
library(numDeriv)
library(rgdal)
library(stats)
library(stringr)
library(scales)
library(splines)
library(tidyverse)
library(tmap)
library(viridis)



###------------Worldwide old data only to 2020/12/31-------------------###
# FROM: https://www.ecdc.europa.eu (This site is no longer updated)
data <- read.csv("https://opendata.ecdc.europa.eu/covid19/casedistribution/csv", na.strings = "", fileEncoding = "UTF-8-BOM" , stringsAsFactors = F)
head(data)
#Date format modification
data$date_reported <- dmy(data$dateRep)
#add death rate
data$death_rate <- percent(data$deaths/data$cases)
data[is.na(data)] <- 0
# total cases(Confirmed) worldwide to date
sum(data$cases)
# total deaths worldwide to date
sum(data$deaths)

# total cases and max single day by country
total_cases<-data %>% 
  group_by(countriesAndTerritories) %>% 
  summarise(cases_sum = sum(cases), cases_max = max(cases)) %>% 
  arrange(desc(cases_sum)) %>% as.data.frame()

#ggplot(data = total_cases ) +geom_bar( aes(x = countriesAndTerritories,fill=cases_sum)) +coord_flip()

# total deaths and max single day by country
total_deaths<-data %>% 
  group_by(countriesAndTerritories) %>% 
  summarise(deaths_sum = sum(deaths), deaths_max = max(deaths)) %>% 
  arrange(desc(deaths_sum)) %>% as.data.frame()

#plot cases in us
us <- data[data$countriesAndTerritories == 'United_States_of_America',]
head(us)
US_cases <- ggplot(us, 
                   aes(date_reported, as.numeric(cases))) +
  geom_col(fill = 'blue', alpha = 0.6) + 
  theme_minimal(base_size = 14) +
  xlab(NULL) + ylab(NULL) + 
  scale_x_date(date_labels = "%Y/%m/%d")
US_cases + labs(title="Daily COVID-19 Cases in US")

#plot death in us
US_deaths <- ggplot(us, 
                    aes(date_reported, as.numeric(deaths))) +
  geom_col(fill = 'purple', alpha = 0.6) + 
  theme_minimal(base_size = 14) +
  xlab(NULL) + ylab(NULL) + 
  scale_x_date(date_labels = "%Y/%m/%d")
US_deaths + labs(title="Daily COVID-19 Deaths in US")


#plot death rate
us[is.na(us)] <- 0
ggplot(us, aes(x = date_reported, y = death_rate,group = 1)) + geom_line() + 
  geom_point(size = 4,shape = 20, colour = "darkred")


### UK weekly Confirmed
#plot cases in uk
uk <- data[data$countriesAndTerritories == 'United_Kingdom',]
head(uk)
UK_cases <- ggplot(uk, 
                   aes(date_reported, as.numeric(cases))) +
  geom_col(fill = 'blue', alpha = 0.6) + 
  theme_minimal(base_size = 14) +
  xlab(NULL) + ylab(NULL) + 
  scale_x_date(date_labels = "%Y/%m/%d")
UK_cases + labs(title="weekly COVID-19 Cases in UK")

###UK weekly death
#plot cases in uk
uk <- data[data$countriesAndTerritories == 'United_Kingdom',]
head(uk)
UK_cases <- ggplot(uk, 
                   aes(date_reported, as.numeric(deaths))) +
  geom_col(fill = 'blue', alpha = 0.6) + 
  theme_minimal(base_size = 14) +
  xlab(NULL) + ylab(NULL) + 
  scale_x_date(date_labels = "%Y/%m/%d")
UK_cases + labs(title="weekly COVID-19 Deaths in UK")


####----------Worldwide new data to 2021/12/28------------------#
Worldwide_data <- read.csv("D:\\研究所Meeting\\COVID19\\csse_covid_19_data\\World_covid_19_(14days)week_reports\\owid-covid-data(modify).csv")
Worldwide_data[is.na(Worldwide_data)] <- 0
#Specify date information(2021/12/29)
#let date from "/" to "-"
Worldwide_data$date<-ymd(Worldwide_data$date)
Worldw_28 <- Worldwide_data[Worldwide_data$date=="2021-12-28",]

# total cases single day by country
Ww28total_cases<-Worldw_28 %>% 
  group_by(location) %>% 
  summarise(total_cases = total_cases) %>% 
  arrange(desc(total_cases)) %>% as.data.frame()

# total deaths single day by country
Ww28total_deaths<-Worldw_28 %>% 
  group_by(location) %>% 
  summarise(total_deaths = total_deaths) %>% 
  arrange(desc(total_deaths)) %>% as.data.frame()

# total cases(Confirmed) worldwide to date
Ww28total_cases$total_cases[1]
# total deaths worldwide to date
Ww28total_deaths$total_deaths[1]


#boxplot cases in us 2021/12/28
wwus <- Worldwide_data[Worldwide_data$location == 'United States',]
head(wwus)
wwUS_cases <- ggplot(wwus,
                   aes(date, as.numeric(total_cases))) +
  geom_col(fill = 'blue', alpha = 0.5) + 
  theme_minimal(base_size = 9) +
  xlab(NULL) + ylab(NULL)+
  # Only draw the area at the specified time
  scale_x_continuous(limits = c(2021/01/01, 2021/12/28)) + 
  scale_x_date(date_labels = "%Y/%m/%d")
wwUS_cases + labs(title="Daily COVID-19 Cases in US")

#boxplot death in us 2022/12/28
wwUS_deaths <- ggplot(wwus, 
                    aes(date, as.numeric(total_deaths))) +
  geom_col(fill = 'purple', alpha = 0.5) + 
  theme_minimal(base_size = 9) +
  xlab(NULL) + ylab(NULL) + 
  scale_x_date(date_labels = "%Y/%m/%d")
wwUS_deaths + labs(title="Daily COVID-19 Deaths in US")


#------compair education & poverty & income with different counties--------------------#
califor_poverty <- read.csv("D:\\研究所Meeting\\COVID19\\Poverty\\california_P(modify).csv")
califor_education <- read.csv("D:\\研究所Meeting\\COVID19\\Education\\california_E(modify).csv")
califor_unemployment <- read.csv("D:\\研究所Meeting\\COVID19\\Unemployment\\california_U(modify).csv")
califor_population <- read.csv("D:\\研究所Meeting\\COVID19\\populations\\California counties population(2020).csv")
califor_population$population_rate <- califor_population[,2]/sum(califor_population[,2])
califor_housunit <- read.csv("D:\\研究所Meeting\\COVID19\\housing units\\housing_units_-_single_multi_mobile_-_by_county_2018(California).csv")
califor_location <- read.csv("D:\\研究所Meeting\\COVID19\\ca location\\county location.csv")
califor_location_modi <- read.csv("D:\\研究所Meeting\\COVID19\\ca location\\county location(modify).csv")

# califor hospital beds(Occupation rate)
califor_bed <- read.csv("D:\\研究所Meeting\\COVID19\\Hospital beds\\Hospitals(california).csv")
califor_bed$bed_rate <- califor_bed[,2]/sum(califor_bed[,2])

# califor confirm(2021/2/27)
califor_conf <- read.csv("D:\\研究所Meeting\\COVID19\\csse_covid_19_data\\California data\\statewide_cases(2021_2_27modify).csv")
califor_confirm <- califor_conf[c(-3,-4,-5)]
califor_confirm$confirm_rate <- califor_confirm[,2]/sum(califor_confirm[,2])
#califor_confirm_rate <-percent(califor_confirm[,2]/sum(califor_confirm[,2]))

# texas_poverty <- read.csv("D:\\研究所Meeting\\COVID19\\Poverty\\Texas_P(modify).csv")
# texas_education <- read.csv("D:\\研究所Meeting\\COVID19\\Education\\Texas_E(modify).csv")
# texas_unemployment <- read.csv("D:\\研究所Meeting\\COVID19\\Unemployment\\Texas_U(modify).csv")
# texas_bed <- read.csv("D:\\研究所Meeting\\COVID19\\Poverty\\california(modify).csv")
# texas_population  <- read.csv("D:\\研究所Meeting\\COVID19\\")

# merge california data
fact1 <- merge(califor_poverty,califor_education,
                            by.x = 'subregion',
                            by.y = 'Name')

fact2 <- merge(fact1,califor_unemployment,
                      by.x = 'subregion',
                      by.y = 'Name')

fact3 <- merge(fact2,califor_population,
                      by.x = 'subregion',
                      by.y = 'County')


fact4 <- merge(fact3,califor_housunit,
               by.x = 'subregion',
               by.y = 'County')

fact5 <- merge(fact4,califor_bed,
               by.x = 'subregion',
               by.y = 'COUNTY')

fact6 <- merge(fact5,califor_confirm,
               by.x = 'subregion',
               by.y = 'county')
#merge death cases
#califor_conf$county <- tolower(califor_conf$county)
fact7 <- merge(fact6,califor_conf[c(-2,-4,-5)],
               by.x = 'subregion',
               by.y = 'county')

# ffact <- merge(texas_poverty,texas_education,
#                by.x = 'Name',
#                by.y = 'Name')
# 
# texas_fact <- merge(ffact,texas_unemployment,
#                by.x = 'Name',
#                by.y = 'Name')

# total confirm in california
total_ca <- sum(fact6$confirmed)
total_ca

#-----------------california mapping with data --------------------------#
states <- map_data("state")
ca_df <- states %>%
  filter(region == "california")

head(ca_df)

# california county map
counties <- map_data("county")
ca_county <- counties %>%
  filter(region == "california")

head(ca_county)

# merge with county location 
fact5$subregion <- tolower(fact5$subregion)
califor_factor <- full_join(ca_county, fact5, by = "subregion") 
head(califor_factor)
#add confirmed
fact6$subregion <- tolower(fact6$subregion)
ca_fact <- full_join(ca_county, fact6, by = "subregion")
head(ca_fact)

###eliminate los angeles 
ca_fac <- ca_fact[-c(845:895),]

# californic base map
ca_base <- ggplot(data = ca_df, mapping = aes(x = long, y = lat, group = group)) + 
  coord_quickmap() + 
  geom_polygon(color = "black", fill = "gray") 
ca_base + theme_void()

ca_basemap<-ca_base + theme_void() + 
  geom_polygon(data = ca_county, fill = NA, color = "white") +
  geom_polygon(color = "black", fill = NA)  # get the state border back on top
ca_basemap
#-------------------------about  poverty------------------------#
#map
califor_poverty_map <- ca_base + 
  geom_polygon(data = califor_factor, aes(fill = poverty.percent.), color = "white") +
  geom_polygon(color = 'pink', fill = NA) +#frame color
  labs(title = 'califor poverty') + theme(plot.title = element_text(hjust = 0.5))
califor_poverty_map
#normalization (value between 0 to 1)
califor_poverty$poverty_normal <- unlist(predict(preProcess(califor_poverty[2], method=c("range")), califor_poverty[2]))


#bar plot
ggplot(califor_poverty,aes(subregion,poverty.percent.))+ 
  geom_col()+ geom_bar(stat = "identity",fill = 'pink') +
  scale_fill_brewer(palette = 'Blues') +
  labs(title = 'califor poverty') + theme(plot.title = element_text(hjust = 0.5)) + 
  geom_text(aes(label = str_c(round(poverty.percent.,1),''), y = poverty.percent.), color = 'steelblue') +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) + labs(x = 'county',y = 'poverty')
#order bar
ggplot(califor_poverty,aes(reorder(subregion , -poverty.percent.), poverty.percent.))+ 
  geom_col()+ geom_bar(stat = "identity",fill = 'pink') +
  scale_fill_brewer(palette = 'Blues') +
  labs(title = 'califor poverty') + theme(plot.title = element_text(hjust = 0.5)) + 
  geom_text(aes(label = str_c(round(poverty.percent.,1),''), y = poverty.percent.), color = 'steelblue') +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) + labs(x = 'county',y = 'poverty')



#------------------------about confirmed--------------------------# 
#changing to lower case
califor_location$county <- tolower(califor_location$county)
califor_confirm$county <- tolower(califor_confirm$county)
cali_confir_longlat <- merge(califor_location,califor_confirm,
                             by.x = 'county',
                             by.y = 'county')
head(cali_confir_longlat)
#dot map
ca_basemap+
  geom_point(aes(x = long, y = lat, group = NULL),data = cali_confir_longlat)


#map 
califor_confirmed_map <- ca_base + 
  geom_polygon(data = ca_fact, aes(fill = confirmed), color = "white") +
  geom_polygon(color = 'pink', fill = NA) +
  scale_fill_gradientn(colours = plasma(7))  +
  labs(title = 'califor confirmed case') + theme(plot.title = element_text(hjust = 0.5))
#zoom in
#+ xlim(-123, -121.0) + ylim(36, 38)
califor_confirmed_map

#map (elimate los angele)
califor_confirmed_map <- ca_base + 
  geom_polygon(data = ca_fac, aes(fill = confirmed), color = "white") +
  geom_polygon(color = 'pink', fill = NA) +
  scale_fill_gradientn(colours = plasma(7))  +
  labs(title = 'califor confirmed case') + theme(plot.title = element_text(hjust = 0.5))
califor_confirmed_map

###eliminate los angeles(no los)
califor_confirnolos<- califor_confirm[-19,]
# #standardization
califor_confirnolos[2]<-apply(califor_confirnolos[2], 2, scale) #1=row, 2=column
head(califor_confirnolos)
#normalization (value between 0 to 1)
califor_confirm$confirm_normal <- unlist(predict(preProcess(califor_confirm[2], method=c("range")), califor_confirm[2]))
califor_confirnolos$confirm_normal <- unlist(predict(preProcess(califor_confirnolos[2], method=c("range")), califor_confirnolos[2]))

#bar plot
ggplot(califor_confirm,aes(county,confirmed))+ 
  geom_col()+ geom_bar(stat = "identity",fill = 'pink') +
  scale_fill_brewer(palette = 'Blues') +
  labs(title = 'califor confirmed case') + theme(plot.title = element_text(hjust = 0.5)) +
  geom_text(aes(label = str_c(round(confirmed,1),''), y = confirmed), color = 'steelblue') +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) + labs(x = 'county',y = 'confirmed')
#order bar
ggplot(califor_confirm,aes(reorder(county , -confirmed), confirmed))+ 
  geom_col()+ geom_bar(stat = "identity",fill = 'pink') +
  scale_fill_brewer(palette = 'Blues') +
  
  labs(title = 'califor confirmed case') + theme(plot.title = element_text(hjust = 0.5))+
  geom_text(aes(label = str_c(round(confirmed,1),''), y = confirmed), color = 'steelblue') +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) + labs(x = 'county',y = 'confirmed')

##bar plot after normalization
#bar 
ggplot(califor_confirnolos,aes(county,confirm_normal))+ 
  geom_col()+ geom_bar(stat = "identity",fill = 'pink') +
  scale_fill_brewer(palette = 'Blues') +
  labs(title = 'califor confirmed case') + theme(plot.title = element_text(hjust = 0.5)) +
  geom_text(aes(label = str_c(round(confirm_normal,1),''), y = confirm_normal), color = 'steelblue') +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) + labs(x = 'county',y = 'confirmed')
#order bar
ggplot(califor_confirnolos,aes(reorder(county , -confirm_normal), confirm_normal))+ 
  geom_col()+ geom_bar(stat = "identity",fill = 'pink') +
  scale_fill_brewer(palette = 'Blues') +
  labs(title = 'califor confirmed case') + theme(plot.title = element_text(hjust = 0.5)) +
  geom_text(aes(label = str_c(round(confirm_normal,1),''), y = confirm_normal), color = 'steelblue') +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) + labs(x = 'county',y = 'confirmed')

#---------------------------about education--------------------------------# 
#map
califor_education_map <- ca_base + 
  geom_polygon(data = califor_factor, aes(fill = education.2015.2019.), color = "white") +
  geom_polygon(color = "black", fill =  NA) +
  scale_fill_gradientn(colours = terrain.colors(7)) +
  labs(title = 'califor education') + theme(plot.title = element_text(hjust = 0.5))
  #+theme_void()
califor_education_map
#normalization (value between 0 to 1)
califor_education$education_normal <- unlist(predict(preProcess(califor_education[2], method=c("range")), califor_education[2]))

#bar
ggplot(califor_education,aes(Name,education.2015.2019.))+ 
  geom_col()+ geom_bar(stat = "identity",fill = 'pink') +
  scale_fill_brewer(palette = 'Blues') +
  labs(title = 'califor education') + theme(plot.title = element_text(hjust = 0.5)) +
  geom_text(aes(label = str_c(round(education.2015.2019.,1),''), y = education.2015.2019.), color = 'steelblue') +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) + labs(x = 'county',y = 'education')
#order bar
ggplot(califor_education,aes(reorder(Name , -education.2015.2019.), education.2015.2019.))+ 
  geom_col()+ geom_bar(stat = "identity",fill = 'pink') +
  scale_fill_brewer(palette = 'Blues') +
  labs(title = 'califor education') + theme(plot.title = element_text(hjust = 0.5)) +
  geom_text(aes(label = str_c(round(education.2015.2019.,1),''), y = education.2015.2019.), color = 'steelblue') +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) + labs(x = 'county',y = 'education')
#------------------------about unemployment------------------------------# 
#map
califor_unemployment_map <- ca_base + 
  geom_polygon(data = califor_factor, aes(fill = Unemployment.2019.), color = "white") +
  geom_polygon(color = "black", fill = NA) +
  scale_fill_gradientn(colours = viridis(7)) +
  labs(title = 'califor unemployment') + theme(plot.title = element_text(hjust = 0.5)) 
  #+theme_void()
califor_unemployment_map
#normalization (value between 0 to 1)
califor_unemployment$unemployment_normal <- unlist(predict(preProcess(califor_unemployment[2], method=c("range")), califor_unemployment[2]))

#bar
ggplot(califor_unemployment,aes(Name,Unemployment.2019.))+ 
  geom_col()+ geom_bar(stat = "identity",fill = 'pink') +
  scale_fill_brewer(palette = 'Blues') +
  labs(title = 'califor Unemployment') + theme(plot.title = element_text(hjust = 0.5)) +
  geom_text(aes(label = str_c(round(Unemployment.2019.,1),''), y = Unemployment.2019.), color = 'steelblue') +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) + labs(x = 'county',y = 'Unemployment')
#order bar
ggplot(califor_unemployment,aes(reorder(Name , -Unemployment.2019.), Unemployment.2019.))+ 
  geom_col()+ geom_bar(stat = "identity",fill = 'pink') +
  scale_fill_brewer(palette = 'Blues') +
  labs(title = 'califor Unemployment') + theme(plot.title = element_text(hjust = 0.5)) +
  geom_text(aes(label = str_c(round(Unemployment.2019.,1),''), y = Unemployment.2019.), color = 'steelblue') +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) + labs(x = 'county',y = 'Unemployment')



#------------------------------------about hospital beds------------------------------------# 
#map
califor_bed_map <- ca_base + 
  geom_polygon(data = califor_factor, aes(fill = BEDS), color = "white") +
  geom_polygon(color = "black", fill = NA) +
  scale_fill_gradientn(colours = viridis(10)) + 
  labs(title = 'califor hospital beds') + theme(plot.title = element_text(hjust = 0.5))
  #theme_void()
califor_bed_map

###eliminate los angeles(no los)
califor_bednolos<- califor_bed[-19,]
#standardization
califor_bednolos[2]<-apply(califor_bednolos[2], 2, scale) #1=row, 2=column
head(califor_bednolos)
#normalization (value between 0 to 1)
califor_bed$bed_normal<-unlist(predict(preProcess(califor_bed[2], method=c("range")), califor_bed[2]))
califor_bednolos$bed_normal<-unlist(predict(preProcess(califor_bednolos[2], method=c("range")), califor_bednolos[2]))


#bar plot
ggplot(califor_bed,aes(COUNTY,BEDS))+ 
  geom_col()+ geom_bar(stat = "identity",fill = 'pink') +
  scale_fill_brewer(palette = 'Blues') +
  labs(title = 'califor hospital beds') + theme(plot.title = element_text(hjust = 0.5)) +
  geom_text(aes(label = str_c(round(BEDS,1),''), y = BEDS), color = 'steelblue') +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) + labs(x = 'county',y = 'hospital beds')
#order bar
ggplot(califor_bed,aes(reorder(COUNTY , -BEDS), BEDS))+ 
  geom_col()+ geom_bar(stat = "identity",fill = 'pink') +
  scale_fill_brewer(palette = 'Blues') +
  
  labs(title = 'califor hospital beds') + theme(plot.title = element_text(hjust = 0.5))+
  geom_text(aes(label = str_c(round(BEDS,1),''), y = BEDS), color = 'steelblue') +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) + labs(x = 'county',y = 'hospital beds')

##bar plot after eliminate los angeles 
#bar 
ggplot(califor_bednolos,aes(COUNTY,BEDS))+ 
  geom_col()+ geom_bar(stat = "identity",fill = 'pink') +
  scale_fill_brewer(palette = 'Blues') +
  labs(title = 'califor hospital beds') + theme(plot.title = element_text(hjust = 0.5)) +
  geom_text(aes(label = str_c(round(BEDS,1),''), y = BEDS), color = 'steelblue') +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) + labs(x = 'county',y = 'hospital beds')
#order bar
ggplot(califor_bednolos,aes(reorder(COUNTY , -BEDS), BEDS))+ 
  geom_col()+ geom_bar(stat = "identity",fill = 'pink') +
  scale_fill_brewer(palette = 'Blues') +
  labs(title = 'califor hospital beds') + theme(plot.title = element_text(hjust = 0.5)) +
  geom_text(aes(label = str_c(round(BEDS,1),''), y = BEDS), color = 'steelblue') +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) + labs(x = 'county',y = 'hospital beds')


#----------------------------about population---------------------------# 
#map
califor_population_map <- ca_base + 
  geom_polygon(data = califor_factor, aes(fill = population), color = "white") +
  geom_polygon(color = "black", fill = NA) +
  scale_fill_gradientn(colours = inferno(7))
  #theme_void()
califor_population_map

#ca housunit_single map
califor_housunit_map <- ca_base + 
  geom_polygon(data = califor_factor, aes(fill = poverty.percent.), color = "white") +
  geom_polygon(color = "black", fill = NA) +
  scale_fill_gradientn(colours = terrain.colors(7))
  #theme_void()
califor_housunit_map

#----------------------------about incident rate---------------------------# 
#map
califor_incident_rate_map <- ca_base + 
  geom_polygon(data = califor_factor, aes(fill = population), color = "white") +
  geom_polygon(color = "black", fill = NA) +
  scale_fill_gradientn(colours = inferno(7))
#theme_void()
califor_incident_rate_map



#-----------------------------multiple bar comparisons------------------------------#
install.packages("devtools")
library(devtools)
install_github("kassambara/easyGgplot2")
library(ggplot2)
library(easyGgplot2)
#confirm-poverty
#changing to lower case
califor_confirm$county <- tolower(califor_confirm$county)
#merge data
conf_pove <- merge(califor_confirm[c(1,4)],califor_poverty[c(1,3)],
               by.x = 'county',
               by.y = 'subregion')
#Convert format
conpov <- gather(conf_pove,'confirm_normal','poverty_normal',key = "group",value = "value")

ggplot2.barplot(data=conpov, xName='county', yName="value",
                groupName='group', groupColors=c('#999999','#E69F00'),
                position=position_dodge(),
                #background and line colors
                backgroundColor="white", color="black", 
                xtitle="country", ytitle="value(normalization)", 
                mainTitle="comparison",
                removePanelGrid=TRUE,removePanelBorder=TRUE,
                axisLine=c(0.5, "solid", "black"), xtickLabelRotation=45
) 

#------------------------------creat model--------------------------------------#

# remove family data & rate
ca_modmodel <- fact7[c(-1,-5,-7,-8,-9,-10,-12,-14)]
# see correlation 
ggpairs(ca_modmodel, mapping=aes(color='yellow'))
#normalization data
normali_model<-predict(preProcess(ca_modmodel, method=c("range")), ca_modmodel)
#normalization data (no los)
normali_nolosmodel<-predict(preProcess(ca_modmodel[-19,], method=c("range")), ca_modmodel[-19,])

#--------draw the lm line on ggpairs--------#
lowerFn <- function(data, mapping, method = "lm", ...) {
  p <- ggplot(data = data, mapping = mapping) +
    geom_point(colour = "blue") +
    geom_smooth(method = method, color = "red", ...)
  p
}
#all data
ggpairs(normali_model, lower = list(continuous = wrap(lowerFn, method = "lm")),
        diag = list(continuous = wrap("barDiag", colour = "pink")),
        upper = list(continuous = wrap("cor", size = 5)))

#no los angeles data
ggpairs(normali_nolosmodel, lower = list(continuous = wrap(lowerFn, method = "lm")),
        diag = list(continuous = wrap("barDiag", colour = "pink")),
        upper = list(continuous = wrap("cor", size = 5)))


#---------------------------------model----------------------------------#
# y = confirm
cali_mod1 <- glm.nb(confirm_rate ~ ., data = ca_modmodel)
summary(cali_mod1)
with(cali_mod1, cbind(res.deviance = deviance, df = df.residual,
               covid1231_to_sub_glm = pchisq(deviance, df.residual, lower.tail=FALSE)))
par(mfrow=c(2,2))
plot(cali_mod1)

cali_mod2 <- glm(confirm_rate ~ ., data = normali_nolosmodel)
summary(cali_mod2)
with(cali_mod2, cbind(res.deviance = deviance, df = df.residual,
               covid1231_to_sub_glm = pchisq(deviance, df.residual, lower.tail=FALSE)))
par(mfrow=c(2,2))
plot(cali_mod2)

cali_mod3 <- lm(confirm_rate ~ ., data = ca_modmodel)
summary(cali_mod3)
par(mfrow=c(2,2))
plot(cali_mod3)

#normalize model
cali_mod4 <- lm(confirm_rate ~ ., data = normali_model)
summary(cali_mod4)
par(mfrow=c(2,2))
plot(cali_mod4)

# y = death
cali_mod5 <- glm(deaths ~ ., data = normali_nolosmodel)
summary(cali_mod5)
with(cali_mod5, cbind(res.deviance = deviance, df = df.residual,
                      covid1231_to_sub_glm = pchisq(deviance, df.residual, lower.tail=FALSE)))
par(mfrow=c(2,2))
plot(cali_mod5)

#---------------Huber M estimator-------------------#
X=as.matrix(ca_modmodel[-6])   #parameter
Y=as.matrix(ca_modmodel[6]) #confirm data
#normalization
XX=as.matrix(normali_model[-6])   #parameter
YY=as.matrix(normali_model[6]) #confirm data
#normalization(no los angle)
XXX=as.matrix(normali_nolosmodel[c(-6)])   #parameter
YYY=as.matrix(normali_nolosmodel[6]) #confirm data

#--------------Multi-Resolution Thin-plate Spline Basis Functions--------#
# about "mrts" function I put it below the code
### this is example in mrts
originalPar <- par(no.readonly = TRUE)
knot <- seq(0, 1, l = 25)
b <- mrts(knot, 25)
x0 <- seq(0, 1, l = 58)
# prediction
bx <- predict(b, x0)
par(mfrow = c(5, 5), mar = c(0, 0, 0, 0))
for (i in 1:25) {
  plot(bx[, i], type = "l", axes = FALSE)
  box()
}
par(originalPar)

# use covid-19 location in mrts
# m by d matrix (d<=3) for m locations of d-dimensional 
b <- as.matrix(mrts(knot = coor, 25))
Xx<-seq(from = min(coor$long),to = max(coor$long),length.out =22 )
Yy<-seq(from = min(coor$lat),to = max(coor$lat),length.out =22 )
gridpo <- expand.grid(x = Xx, y = Yy)

#eleminate grid point outside map 
grid_poin <- gridpo[-c(1:17,23:38,45:58,67:78,89:97,111:119,133:140
                       ,155:161,221:224,243:246,265:267,287:289,309,331,353,
                       177:182,199:204,176,196:198,217:220,237:242,258:264,279:286,
                       299:308,320:330,341:352,363:374,385:396,407:418,429:440,451:462,473:484),]

#phi s0 TPS basis function
b_grid <- as.matrix(mrts(knot = grid_poin, 25))

XXXX=as.matrix(bx)   #parameter
YYYY=as.matrix(ca_modmodel[6]) #confirm data

# orthogonal bspline function
#De:dimension , L:number of basis lower&upper:lower or upper data value , crep_input:x(explain variable)
B_orth <- function(De = k, L, lower , upper , crep_input = NULL, approx = FALSE){
  
  knot.vector <- seq(lower, upper, length.out = De + L + 1 - 2 *De)
  ikn = knot.vector[-c(1 ,length(knot.vector))]
  bkn = knot.vector[ c(1 ,length(knot.vector))]
  
  if(approx){
    # 是否用數值積分
    if(is.null(crep_input)){
      # [0 1]等切割10e6份
      nrep <- 1e6
      crep <- 1:nrep
      crep <- crep/nrep 
    }
    else if(length(crep_input) < 1e6){
      # 小於10e6份則用U(0,1)補上去
      nrep <- 1e6
      crep <- sort(c(crep_input, runif(1e6-length(crep_input) ,min=lower ,max=upper)))
    }
    else{
      nrep <- length(crep_input)
      crep <- crep_input
    }
  }else{
    # 不用數值積分
    if(is.null(crep_input)){
      nrep <- 1e4
      crep <- 1:nrep
      crep <- crep/nrep
    }else{
      nrep <- length(crep_input)
      crep <- crep_input
    }
  }
  ##### Compute a0, bs0, bs from Toshio Hondo
  
  bs0 <- bs(crep, knots = ikn ,Boundary.knots = bkn, degree = De, intercept = T)
  LL <- dim(bs0)[2]# this is the dimension of the basis
  # total knot number + De -1
  
  # new spline basis for matrix computation ( nrep * LL )
  bs1 <- matrix(0, nrow = nrep, ncol = LL)
  
  # coefficient matrix LL * LL
  a0 <- diag(LL)
  
  # initialization of bs1
  sqLL <- 1/sqrt(LL)
  bs1[,1] <- sqLL* sqrt(1 /(upper - lower) ) 
  bs1[,2] <- sqrt(24 / (2* (upper - lower)^3) ) *sqLL*(crep - (upper + lower)/2 )
  # 利用bs0矩陣的 2至4行(即頭尾沒有用到)，作為bs1的 3至5行
  for (ii in 3:LL)
    bs1[,ii]<- bs0[,(ii-1)]
  
  # computation of a0
  for (jj in 3:LL){ #for (jj in 3:3)
    
    ee <- as.numeric(bs1[,jj])
    
    for (ii in 1:(jj - 1))
      a0[jj,ii] <- -LL * ((upper - lower) * mean(ee*bs1[,ii]))
    
    # compute B_jj(z) by Gram-Schmit
    for (ii in 1:(jj - 1))
      ee <- ee + a0[jj,ii] *as.numeric(bs1[,ii])
    
    nee <- sqrt((upper - lower) * mean(ee**2)) # norm of B_jj(z)
    
    bs1[,jj] <- ee /nee *sqLL # normalize B_jj(z)
    a0[jj,] <- a0[jj,] /nee *sqLL
  }
  
  if(is.null(crep_input)||approx == FALSE){
    b_function <- bs1
    z <- crep
  }else{
    idx <- sapply(1:length(crep_input), function(i) which(crep == crep_input[i])[1])
    b_function <- bs1[unlist(idx),]
    z <- crep_input
  }
  return(list(b_function = b_function, z = z))
}

bspline <- function(De,L,data){ #m=number of variable
  A <- {}
  #A[[1]] <- B_orth(De, L, lower = min (data), upper = max (data), crep_input = unlist(data), approx = FALSE)
  #A[[1]] <- A[[1]]$b_function
  for(i in 1:ncol(data)){
    A[[i]] <- B_orth(De, L, lower = min (data[,i]), upper = max (data[,i]), crep_input = unlist(data[,i]), approx = FALSE)
    A[[i]] <- A[[i]]$b_function[,-1] #delete the constant trems
  }
  basis <- matrix(unlist(A),nrow = nrow(data))
  return(basis)
  }

XXXXX<-bspline(3,5,ca_modmodel[-6]) #based function
YYYYY=as.matrix(ca_modmodel[6]) #confirm data
#unlist

#------Huber----------#

#Standardization
stan_e<-function(ei){
  sigma<-median(abs(ei-median(ei)))/0.6745
  e_stan<-ei/sigma   
  return(e_stan)
}

new_residu<-function(X,Y,beta){
  tt<- X %*% beta
  new_res <- as.matrix(Y-tt)
  return(new_res)
}

# huber weight function(second differential)
hu_weifun <-function(residual.wei, k ){
  temp_2<-c() 
  # j=1
  # k=1.345
  for (j in 1:length(residual.wei)) {
    if(abs(residual.wei[j]) <= k){
      temp_2[j] <- 1
    } 
    else {
      temp_2[j] <- k/abs(residual.wei[j])
    }
  }
  return(temp_2)
}

# weighted least squares
WLS<-function(X,Y,W){
  w_b = solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% Y
  return(w_b)
}

#iterative reweighted least squares
IRLS<-function(X,Y,max.iter,conv.eps){
  beta_init <- solve(t(X) %*% X) %*% t(X) %*% Y #initialize beta
  beta_prev <- beta_init               #beta_{t-1} (for comparisons)
  mean_fun<- X %*% beta_init
  e_i.prev<- as.matrix(Y-mean_fun)
  
  for(ii in 1:max.iter){
    standar_ei <- stan_e(new_residu(X,Y,beta_prev)) 
    W <-  diag(hu_weifun(residual.wei=standar_ei, k = 1.345))
    new_beta <- WLS(X=X,Y=Y,W=W)
    if( sum((beta_prev-new_beta)^2) <= conv.eps){
      break
    } else {
      e_i.prev <- standar_ei
      beta_prev <- new_beta
    }
    if(ii%%5 == 0) cat("iterative: ", ii, "\n")
  }
  return(list(new_beta,e_i.prev)) 
}

#par(mfrow=c(2,2))
# outcome(original data)
outcome1 <- IRLS(X=X,Y=Y,max.iter=300,conv.eps=1e-10)
outcome1[[1]]
irls.residu1 <-outcome1[[2]]
qqnorm(irls.residu1)
qqline(irls.residu1, col = "steelblue", lwd = 2)

fitt_valu1 <- X %*% outcome1[[1]]
plot(x = fitt_valu1 , y = irls.residu1,main="Residual vs fitted(original)")

# outcome(normalization)
outcome2 <- IRLS(X=XX,Y=YY,max.iter=300,conv.eps=1e-10)
outcome2[[1]]
irls.residu2 <-outcome2[[2]]
qqnorm(irls.residu2)
qqline(irls.residu2, col = "steelblue", lwd = 2)

fitt_valu2 <- X %*% outcome2[[1]]
plot(x = fitt_valu2 , y = irls.residu2,main="Residual vs fitted(normalization)")

# outcome(normalization no los angle)
outcome3 <- IRLS(X=XXX,Y=YYY,max.iter=300,conv.eps=1e-10)
outcome3[[1]]
irls.residu3 <-outcome3[[2]]
qqnorm(irls.residu3)
qqline(irls.residu3, col = "steelblue", lwd = 2)

fitt_valu3 <- X %*% outcome[[1]]
plot(x = fitt_valu3 , y = irls.residu3,main="Residual vs fitted")

#### TPS basis function
# b :K=25 TPS basis(not contain grid point)
outcome4 <- IRLS(X=b,Y=YYYY,max.iter=300,conv.eps=1e-10)
outcome4[[1]]
irls.residu4 <-outcome4[[2]]
qqnorm(irls.residu4)
qqline(irls.residu4, col = "steelblue", lwd = 2)

fitt_valu4 <- X %*% outcome4[[1]]
plot(x = fitt_valu4 , y = irls.residu4,main="Residual vs fitted(original)")

# orthogonal b-spline 
outcome5 <- IRLS(X=XXXXX,Y=YYYYY,max.iter=300,conv.eps=1e-10)
outcome5[[1]]
irls.residu5 <-outcome5[[2]]
qqnorm(irls.residu5)
qqline(irls.residu5, col = "steelblue", lwd = 2)

fitt_valu5 <- XXXXX %*% outcome5[[1]]
plot(x = fitt_valu5 , y = irls.residu5,main="Residual vs fitted(original)")

# = = = = = = = = = = = = = = = generate data (use in variogram) = = = = = = = = = = = = = =#

muldata <- cali_confir_longlat[5]
coor <-  cali_confir_longlat[c(2,3)]


# = = = = = =  My semivariogram = = = = = = = = = #

# = = = = = = data = = = = = #
data_variance <-c()
for (i in 1:length(muldata$confirm_rate)){
  for (j in 1:length(muldata$confirm_rate)){
    #if(i==j){next}
    #else{
      vari <- (1/2)*(muldata[,c("confirm_rate")][i]-muldata[,c("confirm_rate")][j])^2
      data_variance <- c(data_variance,vari)
     # }
  }
  if(i%%5 == 0) cat("running: ", i, "\n")
}

# = = = = = = = = = coordinate = =  = = = = = = #

data_coordinate <-c()
for (i in 1:length(coor[,c("long")])){
  for (j in 1:length(coor[,c("long")])){
   # if(i==j){next}
    #else{
      disten<-sqrt((coor[,c("long")][i]-coor[,c("long")][j])^2+(coor[,c("lat")][i]-coor[,c("lat")][j])^2) 
      data_coordinate <- c(data_coordinate ,disten)
     # }
  }
  if(i%%5 == 0) cat("running: ", i, "\n")
}

# = = = = = = plot the variogram = = = = = = = #
variogram <-as.data.frame(cbind(data_variance,data_coordinate))
plot(x = data_coordinate,y = data_variance,xlab = "u",ylab = "V(u)")# ,xlim = c(0,12)

# = = = = = = = = my sample variogram = = = = = = = #
aa <-variogram[order(variogram$data_coordinate,decreasing = F),]
s = 0.6                  #bin_width
l=floor(min(aa$data_coordinate))   #lower
u=floor(max(aa$data_coordinate))   #upper
y=aa$data_coordinate
t=aa$data_variance
n=(u-l)/s
gg <-rep(0,n)
tt <-rep(0,n)
#sd_tt <-rep(NA,n) #standard deviation of the values in each bin
num_k <-rep(NA,n) # number of pairs in each bin

for (i in 1:n) {
  # i=16
  k <-which((y>s*(i-1)) & y<=s*i)
  num_k[i] <- length(k)/2
  #sd_tt[i] <- sd(t[k])
  tt[i] <- mean(t[k])
}

for (i in 1:n) {
  lllll <- (i-0.5)*s
  gg[i] <- lllll
}

j <-as.data.frame(cbind(gg,tt))
plot(x = j$gg,y = j$tt,xlab = "u",ylab = "V(u)",main=("myself"))
lines.variomodel(tt, cov.pars = c(2.5,12.5), cov.model = "exp", kap = 1.5, nug = 0.2)
#cov.pars : a vector or matrix with the values for the partial sill (sigma square) and range (phi) parameters.
#nugget : a scalar with the value of the nugget (tau square) parameter.
par(mfrow=c(1,2))

points(x=j$gg,y=j$tt,col='red')

# = = = = model fit (use (weighted)Least Square)= = = = #

# uvec : a vector with values used to define the variogram binning
variogram()
vario.b <- variog(coords = coor, data = muldata,uvec=seq(0,12.25762, l=15)) #, max.dist=1)
# variogram cloud
vario.c <- variog(coords = coor, data = muldata, op="cloud")
#binned variogram and stores the cloud
vario.bc <- variog(coords = coor, data = muldata, bin.cloud=TRUE)

par(mfrow=c(2,2))
plot(vario.b, main="binned variogram") 
plot(vario.c, main="variogram cloud")
plot(vario.bc, bin.cloud=TRUE, main="clouds for binned variogram")  

ini.vals <- expand.grid(seq(0.005,0.1,l=20), seq(3,12,l=50)) # sigma^2 (partial sill) and phi (range parameter)
ols <- variofit(vario.b,fix.nug=TRUE, wei="equal", ini=ini.vals)
summary(ols)
lines(ols,lty=3,col=3,lwd=2)
# weight have "npairs" ,"cressie" ,"equal"
wls <- variofit(vario.b, ini=ini.vals, fix.nug=F,wei="cressie")
summary(wls)
lines(wls, lty=2,col=4,lwd=2)


# = = = = = = = = kriging model(by package) = = = = = = = = = = = = #

#convert this basic data frame into a spatial points data frame
coordinates(cali_confir_longlat) <- ~ long + lat
plot(cali_confir_longlat)

#establish the Coordinate Reference System (CRS).
proj4string(cali_confir_longlat) <- CRS('+proj=longlat +datum=NAD83')

#Then reproject to a more appropriate CRS, such as Teale Albers(亚尔勃斯投影). 
#Note the units=km, which is needed to fit the variogram.
TA <- CRS("+proj=aea +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 +x_0=0 +y_0=-4000000 +datum=NAD83 +units=km +ellps=GRS80")
aq <- spTransform(cali_confir_longlat, TA)

# create a template raster to interpolate to. 
# interpolate across California, so bring in the file counties.shp.
boroughs <- readOGR(dsn ="D:\\研究所Meeting\\COVID19\\R code\\CA_Counties( outline boundary)", layer = "CA_Counties_TIGER2016") 
boroughoutline <- fortify(boroughs, region = "COUNTYFP")

# Reproject to have the same CRS as TA.
ca <- spTransform(boroughs, TA)

# plot the points on CA 
plot(ca, border='gray')
points(aq, cex=.5, col='red')
#size of raster(點陣圖)
r <- raster(ca)
res(r) <- 10  # 10 km if your CRS's units are in km
#Go from raster to a SpatialGrid object
g <- as(r, 'SpatialGrid')
#variogram cloud
vcloud <- variogram(confirm_rate~1, locations=aq, width=20, cloud = TRUE)
plot(vcloud)
#sample variogram
sam_vari <- variogram(confirm_rate~1, locations=aq, width=20)
plot(sam_vari)
# fit model and give the original point
expfit_samviri <- fit.variogram(sam_vari, model = vgm(psill = 0.008, model = "Exp", range = 150, nugget = 0))
expfit_samviri 

plot(variogramLine(expfit_samviri, 400), type='l', ylim=c(0,0.005), col='blue', main = 'Exponential variogram model')
points(sam_vari[,2:3], pch=20, col='red')

ordikrig <- gstat(formula = confirm_rate~1, locations = aq, model=expfit_samviri)
# predict or interpolate for our grid g
pred_ordikrig <- predict(ordikrig, g)
# Convert kriged surface to a raster object for clipping(轉換為點陣圖)
rastpred_ordikrig <- raster(pred_ordikrig)
rastpred_ordikrig <- mask(rastpred_ordikrig, ca)

krigi_ca<-tm_shape(rastpred_ordikrig) + 
  tm_raster(n=10, palette="RdBu", auto.palette.mapping=FALSE, 
            title=" kriging in covid-19") +
  tm_legend(legend.outside=TRUE)
krigi_ca

#run 5-fold cross-validation to estimate the test prediction error
RMSE <- function(observed, predicted) {
  sqrt(mean((predicted - observed)^2, na.rm=TRUE))
}

set.seed(1234)
kf <- kfold(nrow(aq))
rmseok <- rep(NA, 5)
for (k in 1:5) {
  test <- aq[kf == k, ]
  train <- aq[kf != k, ]
  gscv <- gstat(formula = confirm_rate~1, locations = train, model=fve.o)
  p <- predict(gscv, newdata = test, debug.level=0)$var1.pred
  rmseok[k] <- RMSE(test$confirm_rate, p)
}
# 5-fold root mean squared error(RMSE)
mean(rmseok)

#------------another way to draw kriging ----------------#
## 1. Create a grid from the values in your points dataframe
## first get the range in data
x.range <- as.integer(range(cali_confir_longlat@coords[,1]))
y.range <- as.integer(range(cali_confir_longlat@coords[,2]))
##2. Create a grid with a slightly larger extent
plot(cali_confir_longlat)
#use the locator to click 4 points beyond the extent of the plot
#and use those to set your x and y extents
locator(4)
x.range <- as.integer(c(-123,-115))
y.range <- as.integer(c(33,41))
## now expand your range to a grid with spacing that you'd like to use in your interpolation
#here we will use 200m grid cells:
grd <- expand.grid(x=seq(from=x.range[1], to=x.range[2], by=200), y=seq(from=y.range[1], to=y.range[2], by=200))

## convert grid to SpatialPixel class
coordinates(grd) <- ~ x+y
gridded(grd) <- TRUE

## test it out - this is a good way of checking that your sample points are all well within your grid. If they are not, try some different values in you r x and y ranges:
plot(grd, cex=1.5)
points(cali_confir_longlat, pch=1, col='red', cex=1)
title("Interpolation Grid and Sample Points")


#sigma(sill), phi(range),tau(nugget)
model.variog<-vgm(psill=0.01, model="Exp", nugget=0, range=4.65)
krig<-krige(formula=confirm_rate ~ 1, locations=cali_confir_longlat, newdata=grd, model=model.variog)

krig.output=as.data.frame(krig)
names(krig.output)[1:3]<-c("long","lat","var1.pred")

plot<-ggplot(data=krig.output,aes(x=long,y=lat))#start with the base-plot and add the Kriged data to it
layer1<-c(geom_tile(data=krig.output,aes(fill=var1.pred)))#then create a tile layer and fill with predicted
layer2<-c(geom_path(data=boroughoutline,aes(long, lat, group=group),colour = "grey40", size=1))#then create an outline
plot+layer1+layer2+scale_fill_gradient(low="#FEEBE2", high="#7A0177")+coord_equal()


#-------------in real data pm2.5 kriging part-----------------------#
# mean function residu
residu<-as.matrix(irls.residu4)
wls$cov.pars[1]
wls$cov.pars[2]
wls$lambda
wls$nugget
grid_poin
TPS_beta_hat <- outcome4[[1]]
# bx:grid point into TPS
phi_s0 <- t(t(TPS_beta_hat) %*% t(b_grid))
  
c_s <-c()
for (i in 1:length(coor[,c("long")])){
  for (j in 1:length(grid_poin[,c("x")])){
    # if(i==j){next}
    #else{
    disten<-sqrt((coor[,c("long")][i]-grid_poin[,c("x")][j])^2+(coor[,c("lat")][i]-grid_poin[,c("y")][j])^2) 
    c_s <- c(c_s ,disten)
    # }
  }
  if(i%%5 == 0) cat("running: ", i, "\n")
}
c_s0 <- matrix(c_s,nrow = nrow(grid_poin))

zigma_theta <-as.matrix(exp(-dist(coor, method = "euclidean", diag = F, upper = FALSE, p = 2))/wls$lambda)
# Xx<-seq(from = min(coor$long),to = max(coor$long),length.out =20 )
# Yy<-seq(from = min(coor$lat),to = max(coor$lat),length.out =20 )
# grid_poin <- expand.grid(x = Xx, y = Yy)
# plot(grid_poin)

y_hat =  phi_s0 + c_s0  %*% solve(zigma_theta+(ols$nugget+ols$cov.pars[1])*diag(58)) %*% residu

# grid point & y_hat
krig_data <- cbind(grid_poin,y_hat)

#plot in map 
#dot map
ca_basemap+
  geom_point(aes(x = x, y = y, group = NULL),data = krig_data)


#map 
krining_map <- ca_base + 
  geom_polygon(data = krig_data , aes(fill = confirm_rate), color = "white") +
  geom_polygon(color = 'pink', fill = NA) +
  scale_fill_gradientn(colours = plasma(7))  +
  labs(title = 'califor confirmed rate') + theme(plot.title = element_text(hjust = 0.5))
#zoom in
#+ xlim(-123, -121.0) + ylim(36, 38)
krining_map


#####----------------------model simulation-----------------------------#####
##[generate data way1]
# multivariate normal (This is covariance tern)
library(mvtnorm)
mul<-function(n) {
  #random generate coords
  #n=100
  #x <- runif(n,max = 10,min = 0)
  #y <- runif(n,max = 10,min = 0)
  # x <- rnorm(n,mean = 0,sd = 2)
  # y <- rnorm(n,mean = 0,sd = 2)
  #generate grid data
   x <- seq(from = 0, to = 10, length = n)
   y <- seq(from = 0, to = 10, length = n)
  xxx<-expand.grid(x = x, y = y)
  #k=1
  #ele_cov = Exp.cov(x1 = xxx[1:20,],x2 = xxx[31:50,], theta = 1.5)#this is same as the next line
  #rdist:computes the pairwise distances between observations in one matrix and returns a 'dist' object
  ele_cov = stationary.cov(x1 = xxx,x2 =NULL, theta = 1, Distance = "rdist", Covariance="Exponential")
  asd =as.vector(t(rmvnorm(n/n,sigma = ele_cov))) 
  return(list(x=xxx,data=asd,ele_cov)) 
}

multi_data<-mul(30)

covar<- multi_data[[2]]
multicoor <- multi_data[[1]] 
grid_po <-multicoor


# (linear mean tern)
simean<-function(x,y){
  # z1=3*x+4*y
  z2 = 11*x^1/2+1.8*y^1/2
  #z3=6*x+0.2*y^3
  # sim_mean<-cbind(z1,z2,z3)
  #  return(sim_mean)
  return(z2)
}
# 

simumean <- simean(multicoor$x,multicoor$y)
simumean

muldata <- simumean + covar
muldata
summary(muldata)
# ordered value
multidata <- muldata
#random rreange value
#set.seed(123)
multidata <- as.data.frame(sample(muldata)) 

#[generate data way2]
#generate grid points(Xx & Yy coordinates)
Xx<-seq(from = 0,to = 1,length.out =15 )
Yy<-seq(from = 0,to = 1,length.out =15 )
grid_po <- expand.grid(x = Xx, y = Yy)
#creat value on observation
Zz<-runif(min = 0,max = 10,n =length(Xx)*length(Yy))

#plot points
#cex:Modify the size of drawing symbols
#lwd:Modify the wide of point
plot(grid_po,cex = 2 ,main = 'Regular points')
#merge coordinates & value
gridpo<-cbind(grid_po,multidata)
colnames(gridpo)[3] <- "value"
# color plot
ggplot(gridpo, aes(x, y))+
  geom_point(aes(color = value)) +
  scale_color_viridis(option = "D")

# heat map
library(plotly)
# random value
heap <- plot_ly(x =gridpo$x ,y =gridpo$y ,z =gridpo$value ,type = "contour") 
heap

#random draw the grid point(with sample)
# This random sample size will equal to line 1248(identity matrix:diag)
nn <- 100 # number of pick points
rangridpo<-sample(1:nrow(gridpo), nn)
#simudata:simulation data(from random draw)
simudata<-gridpo[rangridpo,]
plot(grid_po,cex = 2 ,lwd = 2,main = 'Regular points with Random selected')
# random pick the points and marked in red 
points(simudata[1:2],cex = 2,lwd = 2, col = "red")
# other points are not be select will marked in blue 
simu_grid_poin<-gridpo[-rangridpo,]
points(simu_grid_poin[1:2],cex = 2,lwd = 2, col = "blue")


# use simulation location in mrts to generate basis
bbaa <- 5 # number of basis
simubas <- as.matrix(mrts(knot = simudata[1:2], bbaa))
# plot the 25 basis
par(mfrow = c(5, 5), mar = c(0, 0, 0, 0))
for (i in 1:25) {
  plot(simubas[, i], type = "l", axes = FALSE)
  box()
}

# use basis in Huber M estimate
# TPS basis function
# b :K=25 TPS basis(not contain grid point)
#oc:outcome
oc<- IRLS(X=simubas,Y=simudata[,3],max.iter=100,conv.eps=1e-5)
oc[[1]]
irls.simresidu <-oc[[2]]
qqnorm(irls.simresidu,main="Normal Q-Q plot with mean")
qqline(irls.simresidu, col = "steelblue", lwd = 2)

fitt_simvalu <- simubas %*% oc[[1]]
plot(x = fitt_simvalu , y = irls.simresidu,main="Residual vs fitted(original)")

# = = = = variogram fit (use (weighted)Least Square)= = = = #

# uvec : a vector with values used to define the variogram binning
#variogram()
simuvario.b <- variog(coords = simudata[1:2], data = simudata[,3],uvec=seq(0,10, l=30)) #, max.dist=1)
# variogram cloud
simuvario.c <- variog(coords = simudata[1:2], data = simudata[,3], op="cloud")
#binned variogram and stores the cloud
simuvario.bc <- variog(coords = simudata[1:2], data = simudata[,3], bin.cloud=TRUE)

par(mfrow=c(2,2))
plot(simuvario.b, main="binned variogram") 
plot(simuvario.c, main="variogram cloud")
plot(simuvario.bc, bin.cloud=TRUE, main="clouds for binned variogram")  

simuini.vals <- expand.grid(seq(0.1,600,l=15), seq(0.1,10,l=15)) # sigma^2 (partial sill) and phi (range parameter)
simuols <- variofit(simuvario.b,fix.nug=T, wei="equal",cov.model="exp", ini=simuini.vals)
summary(simuols)
lines(simuols,lty=3,col=3,lwd=2)
# weight have "npairs" ,"cressie" ,"equal"
# I use the first outcome(cov.pars) in ols to be the initial value in wls
simuini.valss <- expand.grid(seq(simuols$cov.pars[1],10,l=15), seq(simuols$cov.pars[2],1,l=15)) # sigma^2 (partial sill) and phi (range parameter)
simuwls <- variofit(simuvario.b, ini=simuini.valss, fix.nug=T,wei="cressie",cov.model="exp")
summary(simuwls)
lines(simuwls, lty=2,col=4,lwd=2)
# choose weight type ols or wls
simuweight <- simuols
simuweight <- simuwls

###-------use in kriging----------###
# mean function residu
# τ^2 + σ^2, corresponds to the variance of the observation process Y
simresidu<-as.matrix(irls.simresidu)
simuweight$cov.pars[1] #sill
simuweight$cov.pars[2] #range
simuweight$kappa # exponential lamda
simuweight$nugget
simu_grid_poin # other unknown points ,see line 1165
TPS_simu_beta_hat <- oc[[1]]
# simu_b_grid:grid point(simu_grid_poin) into TPS basis
# about 'bbaa' see line 1200
simu_b_grid <- as.matrix(mrts(knot = simu_grid_poin[1:2], bbaa))
simu_phi_s0 <- t(t(TPS_simu_beta_hat) %*% t(simu_b_grid))

simu_c_s <-c()
for (i in 1:length(simudata[,c("x")])){
  for (j in 1:length(simu_grid_poin[,c("x")])){
    # if(i==j){next}
    #else{
    simu_disten<-sqrt((simudata[,c("x")][i]-simu_grid_poin[,c("x")][j])^2+(simudata[,c("y")][i]-simu_grid_poin[,c("y")][j])^2) 
    simu_c_s <- c(simu_c_s ,simu_disten)
    # }
  }
  if(i%%5 == 0) cat("running: ", i, "\n")
}
simu_c_s0 <- simuweight$cov.pars[1]*exp(-matrix(simu_c_s,nrow = nrow(simu_grid_poin))/simuweight$cov.pars[2])

simu_zigma_theta <- simuweight$cov.pars[1]*as.matrix(exp(-dist(simudata[-3], method = "euclidean", diag = F, upper = FALSE, p = 2))/simuweight$cov.pars[2])
# Xx<-seq(from = min(coor$long),to = max(coor$long),length.out =20 )
# Yy<-seq(from = min(coor$lat),to = max(coor$lat),length.out =20 )
# grid_poin <- expand.grid(x = Xx, y = Yy)
# plot(grid_poin)

# About diag see line 1158
simu_y_hat =  simu_phi_s0 + simu_c_s0  %*% solve(simu_zigma_theta + (simuweight$nugget+simuweight$cov.pars[1])*diag(nn)) %*% simresidu

# grid point & y_hat
simu_krig_data <- cbind(simu_grid_poin[1:2],simu_y_hat)

# check MSE & resudual
library(Metrics)
simu_MSE = mse(simu_grid_poin$value, simu_krig_data$simu_y_hat)
simu_MSE
simu_residual <- simu_grid_poin$value - simu_krig_data$simu_y_hat
qqnorm(simu_residual,main="Normal Q-Q plot with model fit")
qqline(simu_residual, col = "steelblue", lwd = 2)

# Compare y,y hat,residual
simu_comptable <-cbind(simu_grid_poin$value,simu_krig_data$simu_y_hat,simu_residual)
head(simu_comptable,10)

# predict heat map vs real heat map
preheap <- plot_ly(x =simu_krig_data$x ,y =simu_krig_data$y ,z =simu_krig_data$simu_y_hat ,type = "contour" ) %>% colorbar(title = "Estimate(left)")
preheap
reaheap <- plot_ly(x =simu_grid_poin$x ,y =simu_grid_poin$y ,z =simu_grid_poin$value ,type = "contour") %>% colorbar(title = "Real(right)")
reaheap
fig <- subplot(preheap,reaheap)
fig
# plot data in persp(three-dimensional space)
library(graphics)
persp(x = sort(grid_po$x) , y = sort(grid_po$y), z = value, theta=45, phi=35, r=3, expand=0.6, axes=T,
      ticktype="detailed", xlab="x", ylab="y", zlab="value")
filled.contour(x=0:20, y=0:20, z, color.palette=gray.colors)

#plot in map 
#dot map
ca_basemap+
  geom_point(aes(x = x, y = y, group = NULL),data = simu_krig_data)


#map 
krining_map <- ca_base + 
  geom_polygon(data = krig_data , aes(fill = confirm_rate), color = "white") +
  geom_polygon(color = 'pink', fill = NA) +
  scale_fill_gradientn(colours = plasma(7))  +
  labs(title = 'califor confirmed rate') + theme(plot.title = element_text(hjust = 0.5))
#zoom in
#+ xlim(-123, -121.0) + ylim(36, 38)
krining_map

#------------------------CA temperature ---------------------------#
ca_temp <- read.csv("D:\\研究所Meeting\\COVID19\\temperature in CA\\ca county(temperature123total).csv")

#----------mrts struction--------------------#
mrts<-function (knot, k, x = NULL, maxknot = 5000) 
{
  # to make sure your R is in 64-bit
  is64bit <- length(grep("64-bit", sessionInfo()$platform)) > 
    0
  if ((!is64bit) & (max(NROW(x), NROW(knot)) > 20000)) {
    stop("Use 64-bit version of R for such a volume of data!")
  }
  # number of columns 
  if (NCOL(knot) == 1) {
    xobs <- as.matrix(as.double(as.matrix(knot)))
  }
  else {
    xobs <- apply(knot, 2, as.double)
  }
  #uniquecombs:find the unique rows in a matrix
  Xu <- uniquecombs(cbind(xobs))
  if (is.null(x) & length(Xu) != length(xobs)) {
    x <- xobs
  }
  colnames(Xu) <- NULL
  n <- n.Xu <- NROW(Xu)
  ndims <- NCOL(Xu)
  if (k < (ndims + 1)) {
    stop("k-1 can not be smaller than the number of dimensions!")
  }
  if (maxknot < n) {
    bmax <- maxknot
    Xu <- subKnot(Xu, bmax)
    Xu <- as.matrix(Xu)
    if (is.null(x)) {
      x <- knot
    }
    n <- NROW(Xu)
    n.Xu <- n
  }
  xobs_diag <- diag(sqrt(n/(n - 1))/apply(xobs, 2, sd), ndims)
  if (!is.null(x)) {
    if (NCOL(x) == 1) {
      x <- as.matrix(as.double(as.matrix(x)))
    }
    else {
      x <- as.matrix(array(as.double(as.matrix(x)), dim(x)))
    }
    if (k - ndims - 1 > 0) {
      result <- mrtsrcpp_predict0(Xu, xobs_diag, x, k - 
                                    ndims - 1)
    }
    else {
      X2 <- scale(Xu, scale = FALSE)
      shift <- colMeans(Xu)
      nconst <- sqrt(diag(t(X2) %*% X2))
      X2 <- cbind(1, t((t(x) - shift)/nconst) * sqrt(n))
      result <- list(X = X2[, 1:k])
      x <- NULL
    }
  }
  else {
    if (k - ndims - 1 > 0) {
      result <- mrtsrcpp(Xu, xobs_diag, k - ndims - 1)
    }
    else {
      X2 <- scale(Xu, scale = FALSE)
      shift <- colMeans(Xu)
      nconst <- sqrt(diag(t(X2) %*% X2))
      X2 <- cbind(1, t((t(Xu) - shift)/nconst) * sqrt(n))
      result <- list(X = X2[, 1:k])
    }
  }
  obj <- result$X
  if (is.null(result$nconst)) {
    X2 <- scale(Xu, scale = FALSE)
    result$nconst <- sqrt(diag(t(X2) %*% X2))
  }
  attr(obj, "UZ") <- result$UZ
  attr(obj, "Xu") <- Xu
  attr(obj, "nconst") <- result$nconst
  attr(obj, "BBBH") <- result$BBBH
  attr(obj, "class") <- c("matrix", "mrts")
  class(obj) <- "mrts"
  if (is.null(x)) {
    return(obj)
  }
  else {
    shift <- colMeans(attr(obj, "Xu"))
    X2 <- sweep(cbind(x), 2, shift, "-")
    X2 <- cbind(1, sweep(X2, 2, attr(obj, "nconst"), 
                         "/"))
    if (k - ndims - 1 > 0) {
      obj0 <- as.matrix(cbind(X2, result$X1))
    }
    else {
      obj0 <- as.matrix(X2)
    }
    dimnames(obj) <- NULL
    aname <- names(attributes(obj))
    attributes(obj0) <- c(attributes(obj0), attributes(obj)[setdiff(aname, 
                                                                    c("dim", "dimnames"))])
    return(obj0)
  }
}
