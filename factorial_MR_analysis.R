#' Factorial MR analysis
#'
#' @param sim_list A list of data.frame and a vector of parameters
#'
#' @return A vector of the true and estimated main effects of A and B and their interaction and the bias of each
#' @export
#'
#' @examples factorial_MR_analysis(factorial_MR_simulation()$conf)
#' 
#' #' @author Jeremy A Labrecque
#' @references Labrecque JA, Swanson SA. Mendelian randomization with multiple exposures: The importance of thinking about time. 2019.
#' @references   Rees JMB, Foley C, Burgess S. Factorial Mendelian randomization: using genetic variants to assess interactions. International Journal of Epidemiology. 2019.

factorial_MR_analysis <- function(sim_list) {
  
  ds <- sim_list$data
  params <- sim_list$params
  
  # First stage
  S1 <- lm(formula = paste0("cbind(A1,B1) ~ ", paste0("X",1:30,collapse = " + ")), data = ds)
  pred <- predict(S1)
  pred <- cbind(pred,pred[,1]*pred[,2])
  
  # Second stage
  out <- lm(ds$Y ~ pred)$coef[-1]
  names(out) <- c("est_A", "est_B","intx")
  
  true <- c(true_A=sum(params[c("A0_Y", "A1_Y")]),
               true_B=sum(params[c("B0_Y", "B1_Y")]),
               true_intx=sum(params[c("intx_0","intx_1")]))
  
  bias <- out - true 
  names(bias) <- paste0("bias_", c("A", "B","intx"))
  
  return(c(true, out, bias))
  
}