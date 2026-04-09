#' INPRES School Construction Data (Duflo 2001)
#'
#' Microdata from the 1995 Indonesian Intercensal Survey (SUPAS), used by
#' Duflo (2001) to estimate the returns to schooling from the INPRES
#' school construction program, and reanalyzed by De Chaisemartin and
#' D'Haultfoeuille (2018) as a fuzzy DID application.
#'
#' @format A list with the following components:
#' \describe{
#'   \item{Y}{Numeric vector. Log monthly wages (lwage).}
#'   \item{D}{Binary vector. Completed primary school (yeduc >= 6).}
#'   \item{G}{Binary vector. Treatment district (education increased
#'     between cohorts, CD'H Section 5.2 classification).}
#'   \item{Ti}{Binary vector. Young cohort (born 1968-1972, exposed to INPRES).}
#'   \item{X}{Numeric matrix (N x 197). Covariate dictionary: quartile-binned
#'     district characteristics plus pairwise interactions.}
#'   \item{cluster}{Numeric vector. District identifier for clustered SEs
#'     (290 districts).}
#' }
#'
#' @details
#' The sample includes both sexes with positive wages, born 1957-1962 (old
#' cohort, unexposed) or 1968-1972 (young cohort, exposed to INPRES). Group
#' classification follows CD'H (2018) Section 5.2: districts where the education
#' distribution changed significantly between cohorts (chi-squared test, p <= 0.5)
#' and mean schooling increased are assigned G=1. Districts with no significant
#' change are G=0.
#'
#' @source Duflo, E. (2001). Schooling and Labor Market Consequences of School
#'   Construction in Indonesia: Evidence from an Unusual Policy Experiment.
#'   *American Economic Review*, 91(4), 795-813.
#'
#' @references
#' De Chaisemartin, C. and D'Haultfoeuille, X. (2018). Fuzzy
#' Differences-in-Differences. *Review of Economic Studies*, 85(2), 999-1028.
#'
#' @examples
#' \donttest{
#' data(duflo)
#' # Fuzzy DID with DML
#' fit <- didml(duflo$Y, duflo$D, duflo$G, duflo$Ti, duflo$X,
#'              iv = TRUE, estimator = "wald", method = "lasso", K = 5,
#'              cluster = duflo$cluster)
#' print(fit)
#' }
"duflo"
