setwd('C:\\Users\\kjkomatsu\\OneDrive - UNCG\\NSF BioCON rhizobia\\data\\BioCON data')

#climate data
climate <- read.csv('e080_Daily climate summary.csv')

#get year,  month, day columns
climate$date <- as.Date(as.character(climate$date), format='%m/%d/%Y') #make date column into a date
climate$year <- as.numeric(format(climate$date, '%Y')) #year as numeric
climate$month <- as.factor(format(climate$date, '%B')) #month in text format
climate$day <- as.numeric(format(climate$date, '%d')) #day as numeric

#convert precip from in to cm
climate$precip_cm <- climate$precip_in*2.54

#take average of min and max temps
climate$avg_temp <- (climate$temp_min+climate$temp_max)/2

#get calendar year annual precip and mean max and min temp
climateAnnual <- ddply(climate, c('year'), summarise,
                       tot_precip_cm=sum(precip_cm),
                       avg_temp_min=mean(temp_min),
                       avg_temp_max=mean(temp_max),
                       avg_temp=mean(avg_temp))

#get growing season metrics
climateDailyGrowingSeason <- subset(climate, subset=(gs=='growing'))
climateGrowingSeason <- ddply(climateDailyGrowingSeason, c('year'), summarise,
                              gs_precip_cm=sum(precip_cm),
                              gs_temp_min=mean(temp_min),
                              gs_temp_max=mean(temp_max),
                              gs_temp=mean(avg_temp))

#get growing year annual precip and mean max and min temp
climateGrowingYear <- ddply(climate, c('gs_year'), summarise,
                            gy_precip_cm=sum(precip_cm),
                            gy_temp_min=mean(temp_min),
                            gy_temp_max=mean(temp_max),
                            gy_temp=mean(avg_temp))
names(climateGrowingYear)[names(climateGrowingYear)=='gs_year'] <- 'year'

#checking normality
shapiro.test(climateAnnual$tot_precip_cm)
qqnorm(climateAnnual$tot_precip_cm)

shapiro.test(climateGrowingSeason$gs_precip_cm)
qqnorm(climateGrowingSeason$gs_precip_cm)

shapiro.test(climateGrowingYear$gy_precip_cm)
qqnorm(climateGrowingYear$gy_precip_cm)

#merge all climate data together
climateInitial <- merge(climateAnnual, climateGrowingSeason)
climate <- merge(climateInitial, climateGrowingYear)

#cleanup
rm(list=c('climateAnnual', 'climateDailyGrowingSeason', 'climateGrowingSeason', 'climateGrowingYear', 'climateInitial'))
