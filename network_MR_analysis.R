#' Estimation for network MR simulation
#'
#' @param sim_list Output from network_MR_simulation function
#'
#' @return A vector of the true and estimated total, direct and indirect effect and the bias of each
#' @export
#'
#' @examples
#' network_MR_analysis(network_MR_simulation())
#' 
#' @author Jeremy A Labrecque
#' @references Labrecque JA, Swanson SA. Mendelian randomization with multiple exposures: The importance of thinking about time. 2019.
#' @references Burgess S, Daniel RM, Butterworth AS, et al. Network Mendelian randomization: Using genetic variants as instrumental variables to investigate mediation in causal pathways. International Journal of Epidemiology. 2015;44(2):484-495. (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4469795/pdf/dyu176.pdf)
#' 
#' 
#' 
network_MR_analysis <- function(sim_list) { 
  
  ds <- sim_list$data
  params <- sim_list$params
  
  # Calculate the true values of total, direct and indirect effects
  true_TE <- (params["A0_Y"] + params["A1_Y"] + params["A0_B0"]*params["B0_Y"] + 
    params["A0_B1"]*params["B1_Y"] + 
    params["A0_B0"]*params["B0_B1"]*params["B1_Y"] + 
    params["A1_B1"]*params["B1_Y"]) %>% unname
  true_DE <- (params["A0_Y"] + params["A1_Y"]) %>% unname
  true_IE <- (true_TE - true_DE) %>% unname
  
  # Estimate the total, indirect and direct effects
  TE <- (lm(Y ~ Ga, data  = ds)$coef["Ga"]/lm(A1 ~ Ga, data = ds)$coef["Ga"]) %>% unname
  IE <- (lm(B1 ~ Ga, data  = ds)$coef["Ga"]/lm(A1 ~ Ga, data = ds)$coef["Ga"]*
    lm(Y ~ Gb, data  = ds)$coef["Gb"]/lm(B1 ~ Gb, data = ds)$coef["Gb"]) %>% unname
  DE <- (TE - IE) %>% unname

  return(c("TRUE_TE" = true_TE,
           "TRUE_DE" = true_DE,
           "TRUE_IE" = true_IE,
           "EST_TE" = TE,
           "EST_DE" = DE,
           "EST_IE" = IE,
           "BIAS_TE" = TE - true_TE,
           "BIAS_DE" = DE - true_DE,
           "BIAS_IE" = IE - true_IE,
           params[1:4]))
  
}
