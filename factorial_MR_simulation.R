


#' Factorial MR longitudinal data simulation
#'
#' @param n Number of Observations
#' @param alpha_U Strength of effect of unmeasured confounder on A, B and Y
#' @param B0_A0 Parameter of causal effect of variable before underscore on variable after (in some scenarios the direction is reversed)
#' @param A0_A1 Parameter of causal effect of variable before underscore on variable after
#' @param B0_B1 Parameter of causal effect of variable before underscore on variable after
#' @param A1_Y Parameter of causal effect of variable before underscore on variable after
#' @param B1_Y Parameter of causal effect of variable before underscore on variable after
#' @param A0_B1 Parameter of causal effect of variable before underscore on variable after
#' @param B0_A1 Parameter of causal effect of variable before underscore on variable after
#' @param A0_Y Parameter of causal effect of variable before underscore on variable after
#' @param B0_Y Parameter of causal effect of variable before underscore on variable after
#' @param intx_0 Interaction between A0 and B0 in their effect on Y
#' @param intx_1 Interaction between A1 and B1 in their effect on Y
#'
#' @return A list of four lists each of which is composed of dataframe of data generated under a different causal structure (confounding, collider, pleiotropy, mediation) and vector of parameters from the model
#' @export
#'
#' @examples factorial_MR_simulation()
#' 
#' @author Jeremy A Labrecque
#' @references Labrecque JA, Swanson SA. Mendelian randomization with multiple exposures: The importance of thinking about time. 2019.
#' @references  Sanderson E, Davey Smith G, Windmeijer F, et al. An examination of multivariable Mendelian randomization in the single-sample and two-sample summary data settings. International Journal of Epidemiology. 2018;1-25.
#' 
factorial_MR_simulation <- function(n = 30000, alpha_U = 0.3, 
                                    B0_A0 = 1, B1_A1 = 1,
                                    A0_A1 = 1, B0_B1 = 1,
                                    A1_Y = 1, B1_Y = 1,
                                    A0_B1 = 0, B0_A1 = 0,
                                    A0_Y = 0, B0_Y = 0,
                                    intx_0 = 0, intx_1 = 0) {
  library(magrittr)
  
  # Create umeasured confounder
  U <- rnorm(n, 0, 1)
  
  # Create 30 SNPs with different relationships with A and B
  G <- lapply(1:30,function(x) {
    return(rbinom(n = n, size = 2, prob = 0.3))
  }) %>% do.call(cbind,.)
  
  # Parameters for the effects of the SNPs
  alpha_GA <- c(seq(0.05,0.5,length.out = 10),
                seq(0.5,0.05,length.out = 10),
                rep(0,10))
  
  alpha_GB <- c(rep(0,10),
                seq(0.05,0.5,length.out = 10),
                seq(0.05,0.5,length.out = 10))
  
  # CONFOUNDING SCENARIO-------------------------------------------------------
  # Data generation of A, B and Y
  B0 <- (G %*% alpha_GB) + alpha_U*U
  A0 <- (G %*% alpha_GA) + B0_A0*B0 + alpha_U*U 
  B1 <- B0_B1*B0 + A0_B1*A0 + alpha_U*U
  A1 <- A0_A1*A0 + B1_A1*B1 + B0_A1*A0 + alpha_U*U 
  Y <- A0_Y*A0 + B0_Y*B0 + intx_0*A0*B0 + A1_Y*A1 + B1_Y*B1 + intx_1*A1*B1 + alpha_U*U
  
  ds_conf <- data.frame(U, G, A0, B0, A1, B1, Y)
  params_conf <- c(A0_B1 = A0_B1, B0_A1 = B0_A1, A0_Y = A0_Y, B0_Y = B0_Y,
                   A1_Y = A1_Y, B1_Y = B1_Y, intx_0 = intx_0, intx_1 = intx_1)
  conf <- list(data = ds_conf,
               params = params_conf)
  
  # COLLIDER SCENARIO-----------------------------------------------------------
  # Data generation of A, B and Y
  A0 <- (G %*% alpha_GA) + alpha_U*U
  B0 <- (G %*% alpha_GB) + B0_A0*A0 + alpha_U*U
  A1 <- A0_A1*A0 + B0_A1*B0 + alpha_U*U
  Y  <- A0_Y*A0 + B0_Y*B0 + A1_Y*A1 + intx_0*A0*B0 + alpha_U*U
  B1 <- B0_B1*B0 + A0_B1*A0 + B1_A1*A1 + B1_Y*Y + alpha_U*U 
  
  ds_coll <- data.frame(U, G, A0, B1, A1, Y)
  params_coll <- c(A0_B1 = A0_B1, B0_A1 = B0_A1, A0_Y = A0_Y, B0_Y = B0_Y,
                   A1_Y = A1_Y, B1_Y = 0, intx_0 = intx_0, intx_1 = 0)
  coll <- list(data = ds_coll,
               params = params_coll)
  
  # PLEIOTROPY SCENARIO---------------------------------------------------------
  # Data generation of A, B and Y
  A0 <- (G %*% alpha_GA) + alpha_U*U  
  B0 <- (G %*% alpha_GB) + alpha_U*U
  A1 <- A0_A1*A0 + B0_A1*B0  + alpha_U*U
  B1 <- B0_B1*B0 + A0_B1*A0  + alpha_U*U
  Y <- A0_Y*A0 + B0_Y*B0 + intx_0*A0*B0 + A1_Y*A1 + B1_Y*B1 + intx_1*A1*B1  + alpha_U*U
  
  ds_plei <- data.frame(U, G, A0, B0, A1, B1, Y)
  params_plei <- c(A0_B1 = A0_B1, B0_A1 = B0_A1, A0_Y = A0_Y, B0_Y = B0_Y,
                   A1_Y = A1_Y, B1_Y = B1_Y, intx_0 = intx_0, intx_1 = intx_1)
  plei <- list(data = ds_plei,
               params = params_plei)
  
  # It's worth noting that adding effects of A on B or B on A means that this is no longer a true "pleiotropy" scenario as described in Sanderson 2018
  
  
  # MEDIATION SCENARIO----------------------------------------------------------
  # Data generation of A, B and Y
  A0 <- (G %*% alpha_GA) + alpha_U*U  
  B0 <- (G %*% alpha_GB) + B0_A0*A0 + alpha_U*U
  A1 <- A0_A1*A0 + B0_A1*B0  + alpha_U*U
  B1 <- B0_B1*B0 + A0_B1*A0 + B1_A1*A1 + alpha_U*U
  Y <- A0_Y*A0 + B0_Y*B0 + intx_0*A0*B0 + A1_Y*A1 + B1_Y*B1 + intx_1*A1*B1  + alpha_U*U
  
  ds_medi <- data.frame(U, G, A0, B0, A1, B1, Y)
  params_medi <- c(A0_B1 = A0_B1, B0_A1 = B0_A1, A0_Y = A0_Y, B0_Y = B0_Y,
                   A1_Y = A1_Y, B1_Y = B1_Y, intx_0 = intx_0, intx_1 = intx_1)
  medi <- list(data = ds_medi,
               params = params_medi)
  
  
  # Output
  
  return(list(conf = conf,
              coll = coll,
              plei = plei,
              medi = medi))
  
}





