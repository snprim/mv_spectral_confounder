library(findSVI)
library(sf)
library(tidyverse)
library(spdep)


# set the correct path
setwd("")

# health outcome
# data not available to public

health <- read_csv("health_data.csv")

match <- c("HYPERT", "CHRNKIDN", "HYPERL", "CHF", "DIABETES")

conditions <- health %>% dplyr::select(contains(match))

rates <- conditions %>% dplyr::select(contains(c("IncRateYear")))
health2 <- cbind(health[,1:4], rates)

health3 <- health2 %>% dplyr::mutate(ZCTA = paste(str_pad(ZCTA, width = 5, side = "left", pad = 0)))

health4 <- health3[health3$BASF_YR_NUM == 2017 | health3$BASF_YR_NUM == 2018 | health3$BASF_YR_NUM == 2019,]

health_avg <- health4 %>% na.omit() %>% group_by(ZCTA) %>%
  dplyr::summarise(HYPERT_avg_IncRateYear = mean(HYPERT_IncRateYear),
                   CHRNKIDN_avg_IncRateYear = mean(CHRNKIDN_IncRateYear),
                   HYPERL_avg_IncRateYear = mean(HYPERL_IncRateYear),
                   CHF_avg_IncRateYear = mean(CHF_IncRateYear),
                   DIABETES_avg_IncRateYear = mean(DIABETES_IncRateYear),
                   GroupCount_avg = mean(GroupCount))
health_avg2 <- health_avg[health_avg$GroupCount_avg>100,]
summary(health_avg2$CHF_avg_IncRateYear)
summary(health_avg2$DIABETES_avg_IncRateYear)
health_avg3 <- health_avg2[health_avg2$CHF_avg_IncRateYear>0,]
health_avg4 <- health_avg3[health_avg3$DIABETES_avg_IncRateYear>0,]


# gather SVI data 

svi <- find_svi(year = 2018, geography = "zcta", state = 'US')
svi2 <- svi %>% na.omit()

# SES 

ses <- read_csv("data/outputs/US_ZCTA_ACSvars_2015.csv")
ses2 <- ses[ses$TOTAL>0,] %>% dplyr::select(c("ZCTA", "TOTAL")) %>% na.omit()

# crosswalk between ZCTA, state, region
# for 48 states 

state_zcta_region2 <- read_csv("data/inputs/zcta_state_region.csv")

# combine svi, health, state_zcta_region2, ses

svi_region <- left_join(svi2, state_zcta_region2, by = c("GEOID" = "ZCTA"))
svi_region_health <- full_join(svi_region, health_avg4, by = c("GEOID" = "ZCTA")) %>% na.omit()

svi_region_health_ses <- left_join(svi_region_health, ses2, by = c("GEOID" = "ZCTA")) %>% na.omit()

svi_region_health_ses2 <- svi_region_health_ses %>% dplyr::select(-c("year", "state", "FIPS", "STATE_NAME"))

data_south <- svi_region_health_ses2[svi_region_health_ses2$Region=="South",] 

# shape file 

shape_zcta <- read_sf("tl_2023_us_zcta520/tl_2023_us_zcta520.shp")

data_south_shape <- right_join(shape_zcta, data_south, by = c("ZCTA5CE20" = "GEOID")) %>% na.omit()

zcta_order <- readRDS("data/outputs/zcta_order.rds")

# find neighborhood matrix 
# these files are provided in local_files
# but make sure the ZCTAs are ordered the same way as outputs/zcta_order above

adj.mat_south <- nb2mat(poly2nb(data_south_shape), style = "B", zero.policy = T)
M_south <- diag(rowSums(adj.mat_south))
Qmat_south   <- M_south - adj.mat_south
eig_south <- eigen(Qmat_south)
Gamma_south <- eig_south$vectors
W_south <- eig_south$values

# urbanicity 

urban <- read_csv("data/inputs/nhgis0015_ds172_2010_zcta.csv")
urban2 <- urban %>% dplyr::select(c("NAME", "H7W001", "H7W002", "H7W003", "H7W004", "H7W005", "H7W006"))
urban3 <- urban2 %>% dplyr::rename(Total=H7W001, 
                                   Urban=H7W002,
                                   UrbanArea=H7W003,
                                   UrbanCluster=H7W004,
                                   Rural=H7W005,
                                   NotDefined=H7W006)
urban3$ZCTA <- str_sub(urban3$NAME,7,11)
urban3$Urbanicity <- ifelse(urban3$Urban/urban3$Total>0.5, 1, 0)
urban4 <- urban3 %>% dplyr::select(c("ZCTA", "Urbanicity"))

# final full data 
data_south_shape_urban <- left_join(data_south_shape, urban4, by = c("ZCTA5CE20" = "ZCTA"))

# create Y, X, Z for analysis

data_south_shape_urban2 <- st_drop_geometry(data_south_shape_urban)

exposures <- data_south_shape_urban2 %>% 
  dplyr::select(c("RPL_theme1", "RPL_theme2", "RPL_theme3", "RPL_theme4")) %>% 
  dplyr::rename(theme1 = RPL_theme1,
                theme2 = RPL_theme2,
                theme3 = RPL_theme3,
                theme4 = RPL_theme4)

outcomes <- data_south_shape_urban2 %>% dplyr::select(contains("IncRateYear")) %>% log(.)

covariates <- cbind(log(data_south_shape_urban2$TOTAL), data_south_shape_urban2$Urbanicity)

Y <- as.matrix(outcomes)
X <- as.matrix(exposures)
Z <- as.matrix(covariates)

Y_standardized <- scale(Y)
Z_standardized <- scale(Z)

Y_star <- t(Gamma_south) %*% Y_standardized
X_star <- t(Gamma_south) %*% X
Z_star <- t(Gamma_south) %*% cbind(1, Z_standardized)
