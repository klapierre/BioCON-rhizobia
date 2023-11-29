setwd('C:\\Users\\kjkomatsu\\OneDrive - UNCG\\NSF BioCON rhizobia\\data\\BioCON data')

##### Read in climate data ####
climate <- read.csv('e080_Daily climate summary_2023.csv') %>% 
  #create year, month, day columns
  mutate(date=as.Date(as.character(Date), format='%m/%d/%Y'),
         year=as.numeric(format(date, '%Y')),
         month=as.factor(format(date, '%B')),
         day=as.numeric(format(date, '%d'))) %>% 
  #convert precip from in to cm
  mutate(precip_cm=Precip.inches.*2.54) %>% 
  #take average of min and max temps
  mutate(avg_temp=(MinTemp.degF.+MaxTemp.degF.)/2)


#### Growing season climate summaries ####
climateGrowingSeason <- climate %>% 
  filter(month %in% c('May', 'June', 'July', 'August')) %>% 
  group_by(year) %>% 
  summarise(growing_precip_cm=sum(precip_cm),
            growing_temp_min=mean(MinTemp.degF.),
            growing_temp_max=mean(MaxTemp.degF.),
            growing_temp=mean(avg_temp)) %>% 
  ungroup()


#### Growing year climate summaries ####
climateGrowingYear <- climate %>% 
  mutate(growing_year=ifelse(month %in% c('September', 'October', 'November', 'December'), year+1,
                             year)) %>% 
  filter(growing_year!=2024) %>% 
  group_by(growing_year) %>% 
  summarise(annual_precip_cm=sum(precip_cm, na.rm=T),
            avg_temp_min=mean(MinTemp.degF.),
            avg_temp_max=mean(MaxTemp.degF.),
            avg_temp=mean(avg_temp)) %>% 
  ungroup() %>% 
  rename(year=growing_year)


#### Checking normality ####
shapiro.test(climateGrowingSeason$growing_precip_cm)
qqnorm(climateGrowingSeason$growing_precip_cm)

shapiro.test(climateGrowingYear$annual_precip_cm)
qqnorm(climateGrowingYear$annual_precip_cm)

#merge all climate data together
climateSummary <- left_join(climateGrowingSeason, climateGrowingYear)

#cleanup
rm(list=setdiff(ls(), "climateSummary"))