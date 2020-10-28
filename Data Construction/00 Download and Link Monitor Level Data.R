# install.packages("devtools")
library(devtools)

# Download package tarball from CRAN archive
# url <- "https://cran.r-project.org/src/contrib/Archive/zipcode/zipcode_1.0.tar.gz"
# pkgFile <- "zipcode_1.0.tar.gz"
# download.file(url = url, destfile = pkgFile)
# 
# # Install package
# install.packages(pkgs=pkgFile, type="source", repos=NULL)
# 
# # Delete package tarball
# unlink(pkgFile)
# 
# devtools::install_github("czigler/arepa")
rm(list=ls())
library(zipcode)
library(arepa)

datdir_in = "Other Data Files/"

## ---------------------------------------------------------------------------------------- ##
##   Step 1: Download AQS data, and link zip codes to monitoring locations                  ##
## ---------------------------------------------------------------------------------------- ##
years = 1990:2013
parameter_code = 88101 # for PM 2.5
within_km = 9.656 # 6 miles  ### 4.828 #3 miles #####   
observation_percent = 67 

##--   Get annual PM10 AQS data from EPA
#-- Downloaded on 7/6/2015 - don't download again unless you're prepared to absorb changes to the raw data (which should be minor)
# get_AQS_data_annual(1990:2013)  ## get more years than you'll use for the main analysis.

aqsdat = load_annual_average(years)  ## only load the data that you'll use for the main analysis
PM = subset_monitors(MONITORS = aqsdat, parameter_code, 
                     observation_percent = observation_percent, 
                     monitor_info_file = paste(datdir_in, "/monitor_list.csv", sep = ""))
length(unique(PM$Monitor)) ## 1625 Monitors

## Get zip codes
ZIP = get_zip_codes()

## --- Construct two spatial linkages
## -- Spatial Linkage #1 (non unique)
#  -- Link each monitor to all zipcodes within the specified radius
PMzip_link <- spatial_link_index(PM, "Latitude", "Longitude", "Monitor",
                                 ZIP, "Latitude.zip", "Longitude.zip", "zip",
                                 within = within_km, closest = FALSE)
length(unique(PMzip_link$Monitor)) #1581 Monitors with a zip code located within within_km 
length(unique(PMzip_link$zip)) # 10918 Zip codes located within within_km of a monitor

## -- Spatial Linkage #2 (uniqe)
## ----- Spatial index to join monitors and Zip codes, and assign every 
## ----- zip code to its closest monitor
PMzip_link_unique <- spatial_link_index(PM, "Latitude", "Longitude", "Monitor",
                                 ZIP, "Latitude.zip", "Longitude.zip", "zip",
                                 within = within_km, closest = TRUE)
length(unique(PMzip_link_unique$Monitor)) #1522 Monitors (1522<1581 because some monitors are only within within_km of zip codes that are all closer to another monitor)
length(unique(PMzip_link_unique$zip)) #10918 zip codes


## ----------------------------------------------------------------------------------- ##
##   Step 2: Get zipcode level Census data, aggregate to the Monitor Level             ##
## ----------------------------------------------------------------------------------- ##

### 2000 Census
censusdat = datdir_in
d_census = fread(paste(censusdat, "subset_census_2000_with_zip_and_percentage.csv", sep=""))[, V1 := NULL]
##-----  Variables selected are determined in census_2000.R script and already subsetted here
subset_census = d_census

#### ---  Link census data to PM monitors using spatial linkage #1 (non-unique) 
PMzip_link$ZIP = as.character(PMzip_link$zip)

inPM = unique(PMzip_link$ZIP)
incens = unique(subset_census$ZIP)

inPMnotcens = inPM[!(inPM %in% incens)]
length(inPMnotcens) ## 740 (year 2000) zip codes in PM data not in census data

PMzip_link_census <- merge(PMzip_link, subset_census, by = "ZIP", all.x = TRUE)
length(unique(PMzip_link_census$Monitor)) #1581 Monitors
length(unique(PMzip_link_census$ZIP)) #10918 Zip Codes

### -- Aggregate Census Data to the Monitor Level
CENSUS_Monitor = PMzip_link_census[, list(TotPop = sum(TotPop, na.rm=TRUE), 
                                           PctUrban = sum(Urban, na.rm = TRUE)/sum(TotPop, na.rm = TRUE),
                                           PctInUAs = sum(InUAs, na.rm = TRUE)/sum(TotPop, na.rm = TRUE),
                                           PctInUCs = sum(InUCs, na.rm = TRUE)/sum(TotPop, na.rm = TRUE),
                                           PctRural = sum(Rural, na.rm = TRUE)/sum(TotPop, na.rm = TRUE),
                                           PctWhite = sum(White1, na.rm = TRUE)/sum(TotPop, na.rm = TRUE),
                                           PctBlack = sum(Black1, na.rm = TRUE)/sum(TotPop, na.rm = TRUE),
                                           PctHisp = sum(HispPop, na.rm = TRUE)/sum(TotPop, na.rm = TRUE),
                                           PctHighSchool = sum(HighSchool, na.rm = TRUE)/sum(Over25, na.rm = TRUE),
                                           MedianHHInc = mean(MedianHHInc, na.rm = TRUE),
                                           PctPoor = sum(Poor, na.rm = TRUE)/sum(PovUniverse, na.rm = TRUE),
                                           PctFemale = sum(Female, na.rm = TRUE)/sum(TotPop, na.rm = TRUE),
                                           TotHUs = sum(TotHUs, na.rm = TRUE),
                                           PctUrbanHUs = sum(UrbanHUs, na.rm = TRUE)/sum(TotHUs, na.rm = TRUE),
                                           PctRuralHUs = sum(RuralHUs, na.rm = TRUE)/sum(TotHUs, na.rm = TRUE),
                                           PctOccupied = sum(Occupied, na.rm = TRUE)/sum(TotHUs, na.rm = TRUE),
                                           PctMovedIn5 = sum(MovedIn5, na.rm = TRUE)/sum(Occupied, na.rm = TRUE),
                                           MedianHValue = mean(MedianHValue, na.rm = TRUE),
                                           PopPerSQM = sum(TotPop, na.rm = TRUE)/sum(LandSQMI, na.rm = TRUE)), by = "Monitor"]


## places with missing TotPop and TotHUs will be given values of 0 with the above code.  Set them back to missing.
CENSUS_Monitor$TotPop[CENSUS_Monitor$TotPop==0] = NA
CENSUS_Monitor$TotHUs[CENSUS_Monitor$TotHUs==0] = NA
colMeans(is.na(CENSUS_Monitor)) ## 0.005 (year 2000) missing census data
length(unique(PMzip_link_census$Monitor))  #1581 Monitors


## -------------------------------------------------------------------------------------- ##
##   Step 3: Get Weather data from ASOS and merge to PM monitors                          ##
## -------------------------------------------------------------------------------------- ##
asos = fread(paste(datdir_in, 'Monitor_ASOS.csv', sep = ""))
setkeyv(asos, c("Year", "station"))
asos = unique(asos)
length(unique(asos$station)) #1782 stations

## - Some basic data cleaning here because there are some seriously extreme values (many orders of magnitude too high) - ##
# avrelh is relative humidity as a %
asos$avrelh[asos$avrelh>100] = NA
# avdewf is average dew point temperature, in Fahrenheit
asos$avdewf[asos$avdewf>200] = NA
# avtmpf is average temperature, in Fahrenheit, and looks fine

dim(unique(asos[, c("Year", "lon", "lat"), with = FALSE]))
with(asos, table(table(Year, station)))

## For every monitoring location, find all asos monitors that are within 50km
## Note: this is NOT a unique linkage
PMasos_link = spatial_link_index(PM, "Latitude", "Longitude", "Monitor",
                                 asos, "lat", "lon", "station",
                                 within = 150, closest = FALSE)
length(unique(PMasos_link$station)) #1717 stations
length(unique(PMasos_link$Monitor)) #1581 monitors 

asos_withmonitor = merge(asos, PMasos_link, by = "station", allow.cartesian = TRUE)
length(unique(asos_withmonitor$station)) #1717 stations
length(unique(asos_withmonitor$Monitor)) #1581 monitors

## -- Aggregate to the monitor level
asos_Monitor = asos_withmonitor[, list(avdewf = mean(avdewf, na.rm = TRUE), ## average dew point temp, in Fahrenheit
                                       avtmpf = mean(avtmpf, na.rm = TRUE), ## average temperature, in Fahrenheit
                                       avrelh = mean(avrelh, na.rm = TRUE), ## average relative humidity, in %
                                       nmonitors = length(station)),
                                by = c("Monitor", "Year")]
length(unique(asos_Monitor$Monitor)) #1581 monitors


## ----------------------------------------------------------------------------------- ##
##   Step 4: Get zipcode level Medicare data, aggregate to the Monitor Level           ##
## ----------------------------------------------------------------------------------- ##

##----- Read in simulated annual zip code level Medicare data
BIGMED <- fread(paste(datdir_in, "FakeMedicare.csv", sep = ""))
## -- Notes on the Denominator Variables in Medicare:
#   - Person_year_FFS: use this for the hospitalization outcomes, which are only observed for FFS beneficiaries
#   - Total_den_FFS: use this as the denominator for the _rate variables (e.g., to get the number of Females calculate Female_rate*Total_den_FFS)
#   - Total_den_for_death_MA_FFS: use this as the denominator for Total_death_FFS_MA, since death is measured for all FFS and managed care beneficiaries


## -----  Make a separate data set for each condition
conditions = names(with(BIGMED, table(Condition)))
condition_list = list()
for (i in 1:length(conditions)){
  condition_list[[i]] = subset(BIGMED, Condition == conditions[i])
  condition_list[[i]] = condition_list[[i]][, c("zip", "year", "Total_admission_FFS"), with = FALSE]
  setnames(condition_list[[i]], "Total_admission_FFS", conditions[i])
}

## ----- Now make a data set for all_cause death and covariates, which all get their own fields (and are repeated for each value of 'Condition').  
## ----- Do this by taking only unique zip/year combinations
setkeyv(BIGMED, c("zip", "year"))
BIGMED_unique_zip_year = unique(BIGMED)
BIGMED_unique_zip_year = BIGMED_unique_zip_year[, !(names(BIGMED_unique_zip_year) %in% c("Condition", "Total_admission_FFS")), with = FALSE]

length(unique(BIGMED_unique_zip_year$zip)) ## 38475 unique zip codes across all years
uzips = rep(NA, length(condition_list))
for (i in 1:length(condition_list))
  uzips[i] = length(unique(condition_list[[i]]$zip))
max(uzips) ## maximum number of unique zips in a condition is 38298, so each condition list is a subset of the alldat_unique_zip_year data set

MEDICARE = merge(BIGMED_unique_zip_year, condition_list[[1]], by = c("zip", "year"), all.x = TRUE)
for (i in 2:length(condition_list))
  MEDICARE = merge(MEDICARE, condition_list[[i]], by = c("zip", "year"), all.x = TRUE)

length(unique(MEDICARE$zip)) #38457 unique zip codes

MEDICARE$Year = MEDICARE$year

## -- Unscramble the zip code (zip code is scrambled for data storage purposes)
# temp_zip<-MEDICARE$zipcode_R
# temp_zip<-formatC(temp_zip, width = 5, format = "d", flag = "0")
# zip<-unlist(lapply(temp_zip, function(x) as.numeric(paste(rev(strsplit(as.character(x),"")[[1]]),collapse=""))))
MEDICARE$zip = as.numeric(MEDICARE$zip)
MEDICARE$zip<-formatC(MEDICARE$zip, width = 5, format = "d", flag = "0")
#MEDICARE = cbind(MEDICARE, zip)
# MEDICARE[order(ZIP), c("zipcode_R", "ZIP"), with = FALSE]

## - Restrict only to zip codes that have been linked to monitors
## - NOTE: For Medicare data, make sure to use PM-Zip linkage #2 (unique)!
MEDICARE_subset = subset(MEDICARE, zip %in% PMzip_link_unique$zip)
inMED = unique(MEDICARE_subset$zip)
length(inMED) # 8376 zip codes
inPM = unique(PMzip_link_unique$zip)
length(inPM) # 10918 zip codes
notinMED = inPM[!(inPM %in% inMED)]
length(unique(notinMED)) #2542 zip codes that were linked to a PM monitor but are not in the Medicare data


## -- Get the linked monitor for each zip code that appears in the Medicare data
MEDICARE_withmonitor = merge(MEDICARE_subset, PMzip_link_unique, by = "zip", all.x = TRUE)
length(unique(MEDICARE_withmonitor$zip)) #8376 zip codes
length(unique(MEDICARE_withmonitor$Monitor)) #1510 monitors


##----- Aggragate Medicare Data to the Monitor Level
## Also need to take rate variables in Medicare and transform back to counts
## so that they can be aggragated.

MEDICARE_Monitor <- MEDICARE_withmonitor[, list(mean_age = mean(mean_age, na.rm = TRUE),
                            Female_rate = sum(Female_rate*Total_den_FFS)/sum(Total_den_FFS),
                            White_rate = sum(White_rate*Total_den_FFS)/sum(Total_den_FFS),
                            Black_rate = sum(Black_rate*Total_den_FFS)/sum(Total_den_FFS),
                            Person_year_FFS = sum(Person_year_FFS, na.rm = TRUE),
                            Total_den_FFS = sum(Total_den_FFS, na.rm = TRUE),
                            ALL_CVD = sum(ALL_CVD, na.rm = TRUE),
                            AMI = sum(AMI, na.rm = TRUE),
                            COPD = sum(COPD, na.rm = TRUE),
                            CV_stroke = sum(CV_stroke, na.rm = TRUE),
                            HF = sum(HF, na.rm = TRUE),
                            HRD = sum(HRD, na.rm = TRUE),
                            IHD = sum(IHD, na.rm = TRUE),
                            PVD = sum(PVD, na.rm = TRUE),
                            RTI = sum(RTI, na.rm = TRUE),
                            Tot_den_for_death_MA_FFS = sum(Tot_den_for_death_MA_FFS, na.rm = TRUE),
                            total_death_FFS_MA = sum(total_death_FFS_MA, na.rm = TRUE),
                            nzip = length(unique(ZIP))),
                     by = c("Year", "Monitor")]
length(unique(MEDICARE_Monitor$Monitor)) #1510 Monitors with Medicare information

## ----------------------------------------------------- ##
##   Step 5: Merge All Data at the Monitor Level         ##
## ----------------------------------------------------- ##

length(unique(PM$Monitor))## 1625 monitors in total
unique(PM$Year) ## 1999 - 2013

length(unique(CENSUS_Monitor$Monitor)) ## 1581 have census data available
table(CENSUS_Monitor$Monitor %in% PM$Monitor) ## a subset of the 1625

length(unique(asos_Monitor$Monitor)) ## 1581 have asos data available
table(asos_Monitor$Monitor %in% PM$Monitor) ## a subset of the 1625
unique(asos$Year) ## 1997 - 2014

length(unique(MEDICARE_Monitor$Monitor)) ## 1510 have Medicare data available
table(MEDICARE_Monitor$Monitor %in% PM$Monitor) ## a subset of the 1625
unique(MEDICARE$Year) ## 2001 - 2013

## - First merge (with cartesian) all of the data sources that are time varying
## - the idea is to get any Monitor/Year combination that appears in any of the data sources
PM_MEDICARE = merge(PM, MEDICARE_Monitor, by = c("Monitor", "Year"), all.x = TRUE, all.y = TRUE, allow.cartesian = TRUE)
length(unique(PM_MEDICARE$Monitor)) #1625 Monitors
unique(PM_MEDICARE$Year) ## 1997 - 2013

PM_MEDICARE_ASOS = merge(PM_MEDICARE, asos_Monitor, by = c("Monitor", "Year"), all.x = TRUE, all.y = TRUE, allow.cartesian = TRUE)
length(unique(PM_MEDICARE_ASOS$Monitor)) #1625 Monitors
unique(PM_MEDICARE_ASOS$Year) ## 1997 - 2014

## - Now merge census, which is only for 1 year of data
PM_MEDICARE_ASOS_CENSUS = merge(PM_MEDICARE_ASOS, CENSUS_Monitor, by = "Monitor", all.x = TRUE)
length(unique(PM_MEDICARE_ASOS_CENSUS$Monitor)) #1625 Monitors
unique(PM_MEDICARE_ASOS_CENSUS$Year) ## 1997 - 2014
with(PM_MEDICARE_ASOS_CENSUS, table(is.na(TotPop), Year))

PM_merged = PM_MEDICARE_ASOS_CENSUS

save(PM_merged, file = "PM_merged.Rda")

