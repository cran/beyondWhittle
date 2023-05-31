#' Time series model X_t=e_t, E[e_t]=0
#' @keywords internal
psd_dummy_model <- function() {
  theta_dim <- 0
  excludeBoundary <- T
  get_noise <- function(data, theta, ...) {
    # mean centered version
    if (is.matrix(data)) {
      apply(data, 2, center)
    } else {
      center(data)
    }
  }
  get_data <- function(noise, theta, ...) {
    noise
  }
  propose_next_theta <- function(data, f, previous_theta, ...) {
    # dummy
    numeric(0)
  }
  initialize_theta <- function(data, ...) {
    # dummy
    numeric(0)
  }
  lprior_theta <- function(theta, ...) {
    # dummy
    0
  }
  model_params <- list(theta_dim=theta_dim,
                       get_noise=get_noise,
                       get_data=get_data,
                       propose_next_theta=propose_next_theta,
                       lprior_theta=lprior_theta,
                       initialize_theta=initialize_theta,
                       excludeBoundary=excludeBoundary)
  return(model_params)
}
