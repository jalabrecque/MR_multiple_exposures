library(here)
setwd(here())

library(magrittr)
setwd("V:/HomeDir/495055(J. Labrecque)/Analyses/MR_time/MR_multiple_exposures")
files <- list.files(pattern = "_MR_")
lapply(files, source)

param.grid <- expand.grid(time0 = c(0, 0.5), B0_A1 = c(0, 0.5), A0_B1 = c(0, 0.5))

iter <- 1000

## Network MR


## A is null
t0 <- proc.time()
nw_MR_anull <- lapply(1:nrow(param.grid), function(x) {
  res_mid <- replicate(iter, {
    sim_list <- network_MR_simulation(n=5000, 
                                      A0_Y = 0,
                                      B0_Y = param.grid$time0[x],
                                      B0_A1 = param.grid$B0_A1[x],
                                      A0_B1 = param.grid$A0_B1[x],
                                      A1_Y = 0)
    network_MR_analysis(sim_list)
  }) %>% rowMeans()
}) %>% do.call(rbind, .) %>% as.data.frame %>% round(2)
proc.time() - t0


## B is null
t0 <- proc.time()
nw_MR_bnull <- lapply(1:nrow(param.grid), function(x) {
  res_mid <- replicate(iter, {
    sim_list <- network_MR_simulation(n=5000, 
                                      A0_Y = param.grid$time0[x],
                                      B0_Y = 0,
                                      B0_A1 = param.grid$B0_A1[x],
                                      A0_B1 = param.grid$A0_B1[x],
                                      B1_Y = 0)
    network_MR_analysis(sim_list)
  }) %>% rowMeans()
}) %>% do.call(rbind, .) %>% as.data.frame %>% round(2)
proc.time() - t0



## MVMR


## Effect of A is null, B is a mediator
t0 <- proc.time()
MV_MR_anull_medi <- lapply(1:nrow(param.grid), function(x) {
  res_mid <- replicate(iter, {
    sim_list <- multivariable_MR_simulation(n=5000, 
                                            B0_Y = param.grid$time0[x],
                                            B0_A1 = param.grid$B0_A1[x],
                                            A0_B1 = param.grid$A0_B1[x],
                                            A0_Y = 0,
                                            A1_Y = 0)$medi
    multivariable_MR_analysis(sim_list)
  }) %>% rowMeans()
}) %>% do.call(rbind, .) %>% as.data.frame %>% round(2)
MV_MR_anull_medi <- cbind(param.grid, MV_MR_anull_medi)
MV_MR_anull_medi[order(MV_MR_anull_medi$time0, MV_MR_anull_medi$A0_B1, MV_MR_anull_medi$B0_A1),]
proc.time() - t0


## Effect of A is null, B is confounder
t0 <- proc.time()
MV_MR_anull_conf <- lapply(1:nrow(param.grid), function(x) {
  res_mid <- replicate(iter, {
    sim_list <- multivariable_MR_simulation(n=5000, 
                                            B0_Y = param.grid$time0[x],
                                            B0_A1 = param.grid$B0_A1[x],
                                            A0_B1 = param.grid$A0_B1[x],
                                            A0_Y = 0,
                                            A1_Y = 0)$conf
    multivariable_MR_analysis(sim_list)
  }) %>% rowMeans()
}) %>% do.call(rbind, .) %>% as.data.frame %>% round(2)
MV_MR_anull_conf <- cbind(param.grid, MV_MR_anull_conf)
MV_MR_anull_conf[order(MV_MR_anull_conf$time0, MV_MR_anull_conf$A0_B1, MV_MR_anull_conf$B0_A1),]
proc.time() - t0

H0_test <- list(nw_MR_anull = nw_MR_anull,
                nw_MR_bnull = nw_MR_bnull,
                MV_MR_anull_medi = MV_MR_anull_medi,
                MV_MR_anull_conf = MV_MR_anull_conf)


save(H0_test, file = "V:/HomeDir/495055(J. Labrecque)/Analyses/MR_time/output/H0_testing.Rdata")


## FACTORIAL 


## No intx, anull

t0 <- proc.time()
fac_MR_anull <- lapply(1:nrow(param.grid), function(x) {
  res_mid <- replicate(iter, {
    sim_list <- factorial_MR_simulation(n = 30000,
                                        A0_Y = param.grid$time0[x],
                                        B0_Y = param.grid$time0[x],
                                        A0_B1 = param.grid$A0_B1[x],
                                        B0_A1 = param.grid$B0_A1[x],
                                        intx_0=0, intx_1 = 0)$medi
    factorial_MR_analysis(sim_list)
  }) %>% rowMeans()
}) %>% do.call(rbind, .) %>% as.data.frame %>% round(2)
fac_nointx_MR_anull <- cbind(param.grid, fac_MR_anull)
proc.time() - t0
