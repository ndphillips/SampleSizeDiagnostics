#' Calculate sample size for diagnostic studies using Buderer (1996)'s method and
#'  based on code from https://github.com/sg6667/SampleSizeDiagnostics
#'
#' @param sn numeric. Sensitivity.
#' @param sp numeric. Specificity.
#' @param p numeric. Prevalence.
#' @param w numeric. Max clinically acceptable confidence interval width
#' @param ci numeric. Confidence interval.
#'
#' @return A list with the sample size and the parameters used.
#'
#' @details
#'  Abstract of Buderer (1996): Careful consideration of statistical issues related to the choice of a sample size is critical for achieving
#' meaningful results in research studies designed to evaluate diagnostic tests. When assessing the ability of a
#' diagnostic test to screen for disease, the parameters sensitivity, specificity, and predictive values are of interest.
#' Study sample size requirements can be calculated based on a clinically acceptable degree of precision. the
#' hypothesized values of sensitivity and specificity, and the estimated prevalence of disease in the target population. The simple methods and tables in this paper guide the researcher when deciding how many subjects
#' to sample in a study designed to estimate both the sensitivity and the specificity of a diagnostic test, given a
#' specified precision and estimated disease prevalence.
#'
#' @references Buderer, N. M. F. (1996). Statistical methodology: I. Incorporating the prevalence of disease into the sample size calculation for sensitivity and specificity. Academic Emergency Medicine, 3(9), 895-900.
#'
#' @examples
#'
#' ss_diagnostics(sn = .9, sp = .85, p = .20, w = .10, ci = .95)
#' ss_diagnostics(sn = .8, sp = .9, p = .10, w = .10, ci = .95)
#'
ss_diagnostics <- function(sn = .9,
                           sp = .85,
                           p = .20,
                           w = .10,
                           ci = .95) {
  z <- qnorm(ci)

  a_c <- (z^2) * sn * (1 - sn) / (w^2)
  n1 <- a_c / p

  n1_int <- floor(n1)
  if (n1 != n1_int) {
    n1 <- n1_int + 1
  }

  b_d <- (z^2) * sp * (1 - sp) / (w^2)
  n2 <- b_d / (1 - p)
  n2_int <- floor(n2)
  if (n2 != n2_int) {
    n2 <- n2_int + 1
  }

  n <- max(c(n1, n2))

  out <- list(
    n = n,
    params = list(
      sn = sn,
      sp = sp,
      p = p,
      ci = ci
    )
  )

  return(out)
}
