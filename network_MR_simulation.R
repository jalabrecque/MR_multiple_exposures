

#' Data generation for network MR simulation
#'
#' @param n Number of observations
#' @param A0_B1 
#' @param B0_A1 
#' @param A0_Y 
#' @param B0_Y 
#' @param A0_A1 
#' @param B0_B1 
#' @param A0_B0 
#' @param A1_B1 
#' @param A1_Y 
#' @param B1_Y 
#' @param alpha_Ga Effect of gene Ga on A0
#' @param alpha_Gb Effect of gene Gb on B0
#'
#' @return A dataset and vector of parameters
#' @export
#'
#' @examples
#' network_MR_simulation()
#' 
#' @author Jeremy A Labrecque
#' @references Labrecque JA, Swanson SA. Mendelian randomization with multiple exposures: The importance of thinking about time. 2019.
#' @references Burgess S, Daniel RM, Butterworth AS, et al. Network Mendelian randomization: Using genetic variants as instrumental variables to investigate mediation in causal pathways. International Journal of Epidemiology [electronic article]. 2015;44(2):484-495. (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4469795/pdf/dyu176.pdf)

network_MR_simulation <- function(n = 30000, 
                                  A0_B1 = 0, B0_A1 = 0, A0_Y = 0, B0_Y = 0,
                                  A0_A1 = 1, B0_B1 = 1,
                                  A0_B0 = 1, A1_B1 = NULL,
                                  A1_Y = 1, B1_Y = 1,
                                  alpha_Ga = 0.6, alpha_Gb = 1.2) {
  library(magrittr)
  
  # If A1_B1 is not specified, assign the same value as A0_B0
  if (is.null(A1_B1)) A1_B1 <- A0_B0
  
  # Create 3 unmeasured confounders: u1, u2, u3
  U <- lapply(1:3,function(x) {
    return(rnorm(n = n,mean =  0,sd = 1))
  }) %>% do.call(cbind,.) %>% as.data.frame %>% setNames(c("u1", "u2", "u3"))
  
  # Create 5 error terms
  E <- lapply(1:5,function(x) {
    return(rnorm(n = n,mean =  0,sd = 1))
  }) %>% do.call(cbind,.) %>% as.data.frame %>% setNames(c("ea0", "eb0", "ea1", "eb1", "ey"))
  
  # Create two genes
  Ga <- rbinom(n, 2, 0.3)
  Gb <- rbinom(n, 2, 0.3)
  
  # Data generation where A_t causes B_t (i.e. B mediates the A->Y relationship)
  A0 <- Ga*alpha_Ga + 1*U$u1 + 1*U$u2 + E$ea0
  B0 <- Gb*alpha_Gb + A0_B0*A0 + 1*U$u1 + 1*U$u3 + E$eb0
  A1 <- A0_A1*A0 + B0_A1*B0  + 1*U$u1 + 1*U$u2 + E$ea1
  B1 <- B0_B1*B0 + A0_B1*A0 + A1_B1*A1 + 1*U$u1 + 1*U$u3 + E$eb1
  Y <- A0_Y*A0 + B0_Y*B0 + A1_Y*A1 + B1_Y*B1  + 1*U$u1 + 1*U$u2 + 1*U$u3 + E$ey
  
  # Collect into a data frame, create vector of parameters
  ds <- data.frame(U, E, Ga, Gb, A0, B0, A1, B1, Y)
  params <- c(A0_B1 = A0_B1, B0_A1 = B0_A1, A0_Y = A0_Y, B0_Y = B0_Y,
              A1_Y = A1_Y, B1_Y = B1_Y, A0_A1 = A0_A1, B0_B1 = B0_B1,
              A0_B0 = A0_B0, A1_B1 = A1_B1)
  sim <- list(data = ds,
              params = params)
  
  
  # Output
  return(sim)
  
}





