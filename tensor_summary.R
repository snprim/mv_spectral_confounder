library(tidyverse)
library(rTensor)

# these scripts contain necessary functions 
source("functions_analysis.R")
source("cross_validation.R")


# visualize tensor

fit_svi_K10_L10_stB <- readRDS("fit_svi_K10_L10_stB.rds")
E <- 4
R <- 5
Test <- fit_svi_K10_L10_stB

# calculate gamma from Tl, Te, Tr
n <- nrow(Y_star)
Test_gamma <- array(0, c(L,E,R,9000))
for (i in 1:9000){
  Test_gamma[,,,i] <- find_tol_gamma(Test$Tl[,,i], Test$Te[,,i], Test$Tr[,,i]) 
}

mean_gamma <- rowMeans(Test_gamma, dims = 3)

mean_gamma <- as.tensor(mean_gamma)
set.seed(18)
rank <- 2
cp_decomp <- cp(mean_gamma, num_components = rank, max_iter = 5000, tol = 1e-12)

tensor_margin <- c("Basis functions", "Exposures", "Outcomes")
tensor_dim <- c(L, E, R)

# rearrange cp decomposition data 
# multiply L and R data by (-1) for easier interpretation

tensor_L <- list()
for (k in 1:rank){
  tensor_L[[k]] <- data.frame(K = k, 
                              component = "Basis functions",
                              coef = (-1)* cp_decomp$U[[1]][,k],
                              name = 1:10)
  
}

tensor_L_all <- rbind(tensor_L[[1]],
                      tensor_L[[2]])

tensor_E <- list()
for (k in 1:rank){
  tensor_E[[k]] <- data.frame(K = k, 
                              component = "Exposures",
                              coef = cp_decomp$U[[2]][,k],
                              name = c("Theme 1", "Theme 2", "Theme 3", "Theme 4"))
  
}

tensor_E_all <- rbind(tensor_E[[1]],
                      tensor_E[[2]])

tensor_R <- list()
for (k in 1:rank){
  tensor_R[[k]] <- data.frame(K = k, 
                              component = "Outcomes",
                              coef = (-1)* cp_decomp$U[[3]][,k],
                              name = c("HPT", "CKD", "HPL", "CHF", "DBT"))
}

tensor_R_all <- rbind(tensor_R[[1]],
                      tensor_R[[2]])

tensor_data <- rbind(tensor_L_all, tensor_E_all, tensor_R_all)

L_level = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")
all_levels = c(L_level, exposure_names, c("HPT", "CKD", "HPL", "CHF", "DBT"))

ggplot(data = tensor_data, aes(x = factor(name, levels = all_levels), group = factor(K), color = factor(K), y = coef)) + 
  geom_line() + 
  facet_wrap(~factor(component, levels = tensor_margin), scales = "free") + 
  labs(color = "Rank", x = "", y = "Tensor coefficient") + 
  theme(plot.title = element_text(size=10, face="bold.italic", hjust = 0.5),
        axis.text.x = element_text(size=7)) + 
  ylim(c(-0.3, 0.7)) +
  scale_color_manual(values = c("1" = "#35B779FF",
                                "2" = "#440154FF"))