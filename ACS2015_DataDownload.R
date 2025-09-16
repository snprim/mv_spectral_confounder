
#############################################
##### Author: K. Lloyd Hill #################
#############################################


#American Community Survey -- US ZCTAs

#this script downloads a selection of demographic variables from the US census's
#American Community Survey (ACS). ACS data downloaded for all ZIP code tabulated
#areas (ZCTA) in the US for the year 2015.

#set working directory
#setwd()

#load packages
library(tidycensus)
library(tidyverse)

#if needed, load your census API key
#census_api_key("YOUR KEY HERE", install = T)

#import list of ACS variables to download
myvars <- read.csv(file = "data/inputs/acs5_variables1.csv") %>% 
          select(VAR_ID, VAR_LABEL)

#extract list of variable codes
varlist <- myvars$VAR_ID

#download ACS data with tidycensus's 'get_acs()' function
acsdata <-  get_acs(geography = "zcta", 
            dataset = "acs5", 
            variables = varlist, 
            year = 2015)

#add variable labels to ACS results (eg, "TOTAL POPULATION" instead of "B01003_001")
acsdata <-  acsdata %>% 
            left_join(., myvars, by = c("variable"="VAR_ID"))

#reshape data from long to wide giving each ZCTA a row and each ACS variable a column
acsdata <-  acsdata %>% 
            select(GEOID, VAR_LABEL, estimate) %>% 
            pivot_wider(names_from = VAR_LABEL, values_from = estimate)

#create new variables by combining ACS results. For example, to estimate total
#population over 65, add together population estimates for 65-66, 67-69, etc...
#for both males and females.
acsdata <-  acsdata %>%
            mutate(
              POP_65over = MALE_65to66 + MALE_67to69 + MALE_70to74 + 
                MALE_75to79 + MALE_80to84 + MALE_85over + FEMALE_65to66 + 
                FEMALE_67to69 + FEMALE_70to74 + FEMALE_75to79 + FEMALE_80to84 + 
                FEMALE_85over,
              POP_85over = MALE_85over + FEMALE_85over,
              BLACK_65over = BLACK_65to74M + BLACK_75to84M + BLACK_85overM +
                BLACK_65to74F + BLACK_75to84F + BLACK_85overF,
              WHITE_65over = WHITE_65to74M + WHITE_75to84M + WHITE_85overM +
                WHITE_65to74F + WHITE_75to84F + WHITE_85overF,
              AMIND_65over = AMIND_65to74M + AMIND_75to84M + AMIND_85overM + 
                AMIND_65to74F + AMIND_75to84F + AMIND_85overF,
              ASIAN_65over = ASIAN_65to74M + ASIAN_75to84M + ASIAN_85overM +
                ASIAN_65to74F + ASIAN_75to84F + ASIAN_85overF,
              PACISL_65over = PACISL_65to74M + PACISL_75to84M + PACISL_85overM +
                PACISL_65to74F + PACISL_75to84F + PACISL_85overF,
              OTHER_65over = OTHER_65to74M + OTHER_75to84M + OTHER_85overM + 
                OTHER_65to74F + OTHER_75to84F + OTHER_85overF,
              MIXED_65over = MIXED_65to74M + MIXED_75to84M + MIXED_85overM + 
                MIXED_65to74F + MIXED_75to84F + MIXED_85overF,
              HISPLAT_65over = HISPLAT_65to74M + HISPLAT_75to84M + 
                HISPLAT_85overM + HISPLAT_65to74F + HISPLAT_75to84F + 
                HISPLAT_85overF,
              EDU_HSorHigher = HIGHSCH_M + SOMECOL1_M + SOMECOL2_M + 
                ASSOCIATES_M + BACHELORS_M + MASTERS_M + PROFESS_M + PHD_M + 
                HIGHSCH_F + SOMECOL1_F + SOMECOL2_F + ASSOCIATES_F + 
                BACHELORS_F + MASTERS_F + PROFESS_F + PHD_F,
              EDU_BACHorHigher = BACHELORS_M + MASTERS_M + PROFESS_M + PHD_M + 
                BACHELORS_F + MASTERS_F + PROFESS_F + PHD_F,
              EDU_HSorHigher65over = HIGHSCH_65overM + SOMECOLL_65overM + 
                ASSOCIATES_65overM + BACHELORS_65overM + GRADUATES_65overM + 
                HIGHSCH_65overF + SOMECOLL_65overF + ASSOCIATES_65overF + 
                BACHELORS_65overF + GRADUATES_65overF,
              EDU_BACHorHigher65over = BACHELORS_65overM + GRADUATES_65overM + 
                BACHELORS_65overF + GRADUATES_65overF,
              POVERTY_65over = POVERTY_65to74M + POVERTY_75overM + 
                POVERTY_65to74F + POVERTY_75overF) %>%
        #remove the extra columns
          select(GEOID, TOTAL, MALE, FEMALE, POP_65over, POP_85over, BLACK_POP,
            BLACK_65over, WHITE_POP, WHITE_65over, AMIND_POP, AMIND_65over, 
            ASIAN_POP, ASIAN_65over, PACISL_POP, PACISL_65over, OTHER_POP, 
            OTHER_65over, MIXED_POP, MIXED_65over, HISPLAT_POP, HISPLAT_65over, 
            EDU_HSorHigher, EDU_BACHorHigher, EDU_HSorHigher65over, 
            EDU_BACHorHigher65over, POVERTY_POP, POVERTY_65over) %>%
        #also create a "non-white" estimate by subtracting total population by
        #white population
          mutate(NONWHITE_POP = TOTAL - WHITE_POP, 
            NONWHITE_65over = POP_65over - WHITE_65over)

#ACS data are avaiable as total population estimates, however we are more
#interested in ZCTA-level percentages. Calculate population perecentages from
#ACS estimates
acsdata <-  acsdata %>% 
            mutate(across(3:30, ~ ./TOTAL * 100, .names = '{col}_prct')) %>%
          #also create a variable to describe the 85 and over population as a 
          #percentage of the population over 65
            mutate(POP_85over_prct65over = POP_85over/POP_65over*100)

#rename and rearrange columns as needed
acsdata <-  acsdata %>%
            rename("ZCTA" = GEOID) %>%
            select(ZCTA, TOTAL, MALE, MALE_prct, FEMALE, FEMALE_prct, 
              POP_65over, POP_65over_prct, POP_85over, POP_85over_prct, POP_85over_prct65over, 
              BLACK_POP, BLACK_POP_prct, BLACK_65over, BLACK_65over_prct, 
              WHITE_POP, WHITE_POP_prct, WHITE_65over, WHITE_65over_prct,
              AMIND_POP, AMIND_POP_prct, AMIND_65over, AMIND_65over_prct, 
              ASIAN_POP, ASIAN_POP_prct, ASIAN_65over, ASIAN_65over_prct,
              PACISL_POP, PACISL_POP_prct, PACISL_65over, PACISL_65over_prct, 
              OTHER_POP, OTHER_POP_prct, OTHER_65over, OTHER_65over_prct,
              MIXED_POP, MIXED_POP_prct, MIXED_65over, MIXED_65over_prct,
              NONWHITE_POP, NONWHITE_POP_prct, NONWHITE_65over, NONWHITE_65over_prct,
              HISPLAT_POP, HISPLAT_POP_prct, HISPLAT_65over, HISPLAT_65over_prct,
              EDU_HSorHigher, EDU_HSorHigher_prct, EDU_HSorHigher65over, EDU_HSorHigher65over_prct,
              EDU_BACHorHigher, EDU_BACHorHigher_prct, EDU_BACHorHigher65over, EDU_BACHorHigher65over_prct,
              POVERTY_POP_prct, POVERTY_65over_prct)


#EJ SCREEN ----

#EPA's environmental justice screening tool makes use of several ACS variables to
#characterize social and economic vulnerability of communities. EJscreen is 
#published at the census tract and block group but not at ZCTA. This script
#calculates the ACS-based variables from the EJScreen but for all US ZCTAs

#import list of ACS variables to download
myvars2 <- read.csv(file = "data/inputs/acs5_variables2.csv")

#extract list of variable codes
varlist2 <- myvars2$VAR_ID

#download ACS data with tidycensus's 'get_acs()' function
acsdata2 <-  get_acs(geography = "zcta", 
                    dataset = "acs5", 
                    variables = varlist2, 
                    year = 2015)

#reshape data from long to wide giving each ZCTA a row and each ACS variable a column
acsdata2 <-   acsdata2 %>% 
              select(GEOID, variable, estimate) %>% 
              pivot_wider(names_from = variable, values_from = estimate)

#calculate metrics
acsdata2 <- acsdata2 %>%
            rename("ZCTA"=GEOID) %>%
            mutate(POC_prct = (B03002_001 - B03002_003)/B03002_001*100,
              LOWINC_prct = (C17002_001 - C17002_008)/C17002_001*100,
              UNEMP_prct = B23025_005/B23025_003*100,
              LNGISO_prct = (B16002_004 + B16002_007 + B16002_010 + B16002_013)/
                B16002_001*100,
              NOHSCHL_prct = (B15002_003 + B15002_004 + B15002_005 + B15002_006 + 
                B15002_007 + B15002_008 + B15002_009 + B15002_010 + B15002_020 + 
                B15002_021 + B15002_022 + B15002_023 + B15002_024 + B15002_025 + 
                B15002_026 + B15002_027)/B15002_001*100) %>%
            select(ZCTA, POC_prct, LOWINC_prct, UNEMP_prct, LNGISO_prct, 
              NOHSCHL_prct)

acsdata3 <- left_join(acsdata, acsdata2, by = "ZCTA")

acsdata3[is.na(acsdata3)] <- 0

#export data
write.csv(acsdata3, file = "data/outputs/US_ZCTA_ACSvars_2015.csv", row.names = F)
