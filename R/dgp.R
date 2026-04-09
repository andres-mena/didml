#' Simulate Fuzzy DID Data (Mena 2026 Monte Carlo DGP)
#'
#' Generates data from the simulation DGP used in Mena (2026), Section 5.
#' The design satisfies conditional parallel trends (Assumption 4) but
#' violates unconditional parallel trends. Treatment assignment is fuzzy
#' with compliance depending on covariates. Covariates include 20 signal
#' variables and 80 noise variables correlated with the signal.
#'
#' @param n Integer, sample size.
#' @param p Integer, number of signal covariates (default 20).
#' @param p_aug Integer, total covariates including noise (default 100).
#'   Noise covariates \eqn{X_{21}, \ldots, X_{100}} are correlated
#'   (\eqn{\rho = 0.3}) with signal covariates \eqn{X_1, \ldots, X_5}.
#' @param tau Numeric, true treatment effect (LATE, default 1).
#' @param rho Numeric, Toeplitz correlation for signal covariates
#'   (default 0.5).
#' @param seed Integer, random seed (default \code{NULL}).
#'
#' @return A list with components:
#'   \describe{
#'     \item{Y}{Numeric vector. Observed outcome.}
#'     \item{D}{Binary vector. Treatment take-up (fuzzy).}
#'     \item{G}{Binary vector. Group assignment.}
#'     \item{Ti}{Binary vector. Time period.}
#'     \item{X}{Numeric matrix (\code{n} x \code{p_aug}). Covariates.}
#'     \item{Y0}{Numeric vector. Potential outcome under control.}
#'     \item{Y1}{Numeric vector. Potential outcome under treatment.}
#'     \item{tau}{The true LATE used to generate the data.}
#'   }
#'
#' @details
#' The DGP is:
#' \itemize{
#'   \item \strong{Covariates:} \eqn{X \sim N(0, \Sigma)},
#'     \eqn{\Sigma_{jk} = \rho^{|j-k|} \cdot 2.25}, \eqn{p_0 = 20}.
#'   \item \strong{Group:} \eqn{G \sim Bernoulli(\Lambda(0.6 X_1 + 0.3 X_3))}.
#'   \item \strong{Time:} \eqn{T \sim Bernoulli(0.5)}.
#'   \item \strong{Outcome:}
#'     \eqn{g_0(X) = \Lambda(X_1) + 0.3 X_3 + 0.3 X_1 X_3 + 0.2 X_2^2 + 0.15 X_4},
#'     \eqn{h(X) = 0.4 X_1 + 0.2 X_1^2 + 0.15 X_3},
#'     \eqn{Y(0) = 0.3 G + h(X) T + g_0(X) + \epsilon},
#'     \eqn{Y(1) = Y(0) + \tau + \zeta}.
#'   \item \strong{Treatment:}
#'     \eqn{V = G + 0.4 X_1 + 0.2 X_3 + \nu},
#'     \eqn{D(t) = 1\{V > v_{gt}\}} with
#'     \eqn{v_{11} = 0}, \eqn{v_{10} = v_{01} = v_{00} = 2}.
#'   \item \strong{Noise covariates:} \eqn{X_j = 0.3 X_{parent} +
#'     \sqrt{1 - 0.09} \epsilon_j} for \eqn{j = 21, \ldots, p_{aug}},
#'     where parent cycles through \eqn{X_1, \ldots, X_5}.
#' }
#'
#' Conditional parallel trends holds because \eqn{h(X)} does not depend on
#' \eqn{G}. Unconditional parallel trends fails because \eqn{G} and
#' \eqn{h(X)} share \eqn{X_1, X_3}.
#'
#' @examples
#' # Generate data and estimate
#' d <- dgp_didml(n = 1000, seed = 42)
#' fit <- didml(d$Y, d$D, d$G, d$Ti, d$X,
#'              iv = TRUE, method = "lasso", K = 5, seed = 1)
#' print(fit)
#' cat("True LATE:", d$tau, "\n")
#'
#' # Sharp DID version (set D = G * Ti)
#' d$D <- as.integer(d$G * d$Ti)
#' fit_sharp <- didml(d$Y, d$D, d$G, d$Ti, d$X,
#'                    iv = FALSE, dml = TRUE, method = "ols", K = 5)
#'
#' @references
#' Mena, A. (2026). Double Debiased Machine Learning for
#' Difference-in-Differences under Imperfect Compliance. Working Paper,
#' Brown University. Section 5 (Simulations).
#'
#' @export
dgp_didml <- function(n = 1000L, p = 20L, p_aug = 100L, tau = 1,
                       rho = 0.5, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- as.integer(n)
  p <- as.integer(p)

  # --- Signal covariates X1..X_p ---
  Sigma <- stats::toeplitz(rho^(0:(p - 1))) * 2.25
  X <- MASS::mvrnorm(n, rep(0, p), Sigma)
  colnames(X) <- paste0("X", seq_len(p))

  # --- Group and time ---
  Ti <- stats::rbinom(n, 1, 0.5)
  pG <- stats::plogis(0.6 * X[, 1] + 0.3 * X[, 3])
  G <- stats::rbinom(n, 1, pG)

  # --- Outcome functions ---
  g0 <- stats::plogis(X[, 1]) + 0.3 * X[, 3] + 0.3 * X[, 1] * X[, 3] +
    0.2 * X[, 2]^2 + 0.15 * X[, 4]
  h_x <- 0.4 * X[, 1] + 0.2 * X[, 1]^2 + 0.15 * X[, 3]

  # --- Treatment (fuzzy compliance) ---
  V <- G + 0.4 * X[, 1] + 0.2 * X[, 3] + stats::rnorm(n)
  # Thresholds: (G=1,T=1) -> 0, all others -> 2
  D1 <- Ti * as.integer(V > 0) + (1 - Ti) * as.integer(V > 2)
  D0 <- as.integer(V > 2)
  D <- G * D1 + (1 - G) * D0

  # --- Potential outcomes ---
  Y0 <- 0.3 * G + h_x * Ti + g0 + stats::rnorm(n)
  Y1 <- Y0 + tau + stats::rnorm(n, 0, 0.5)
  Y <- D * Y1 + (1 - D) * Y0

  # --- Augment with noise covariates ---
  if (p_aug > p) {
    p_new <- p_aug - p
    alpha_corr <- 0.3
    parents <- ((seq_len(p_new) - 1L) %% 5L) + 1L
    sd_x <- sqrt(2.25)
    X_noise <- matrix(NA_real_, n, p_new)
    for (j in seq_len(p_new)) {
      eps <- stats::rnorm(n, 0, sd_x)
      X_noise[, j] <- alpha_corr * X[, parents[j]] +
        sqrt(1 - alpha_corr^2) * eps
    }
    colnames(X_noise) <- paste0("X", (p + 1):(p + p_new))
    X <- cbind(X, X_noise)
  }

  list(Y = Y, D = D, G = G, Ti = Ti, X = X,
       Y0 = Y0, Y1 = Y1, tau = tau)
}
