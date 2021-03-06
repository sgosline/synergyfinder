#' HSA synergy score based on highest single agent (HSA) model
#'
#' A function to calculate HSA synergy score based on HSA model
#'
#' @param response.mat a dose-response matrix with concentrations as row names and column names
#' @param correction a parameter to specify if baseline correction is used or not. Defaults to TRUE.
#' @param correction.fixed a parameter to specify which parameters of L.4 or LL.4 functions are fixed and at what value they are fixed for 
#' baseline correction. NAs for parameter that are not fixed. More details in function \code{\link[drc]{drm}}.
#' @param nan.handle a parameter to specify if L.4 function or LL.4 function is used when fitting with LL.4 produces
#' NaNs.
#' @return A matrix of HSA synergy scores for all the dose pairs for a drug combination. For a
#' does pair with at least one zero concentration, 0 is used as the synergy score.
#' @author Liye He \email{liye.he@helsinki.fi}
#' @examples
#' data("mathews_screening_data")
#' data <- ReshapeData(mathews_screening_data)
#' delta.score <- HSA(data$dose.response.mats[[1]])
HSA <- function(response.mat, correction = TRUE, correction.fixed = c(NA, NA, NA, NA), nan.handle = c("LL4", "L4")) {
  if(correction) {
    # correct the response data
    response.mat <- BaselineCorrectionSD(response.mat, correction.fixed, nan.handle)$corrected.mat
  }
  drug1.response <- response.mat[, 1]
  drug2.response <- response.mat[1, ]
  # reference matrix
  ref.mat <- response.mat
  for (i in 2:nrow(response.mat)) {
    for (j in 2:ncol(response.mat)) {
      ref.mat[i, j] <- ifelse(drug1.response[i] > drug2.response[j],
                              drug1.response[i], drug2.response[j])
    }
  }
  # synergy matrix
  syn.mat <- response.mat - ref.mat
  syn.mat
}
