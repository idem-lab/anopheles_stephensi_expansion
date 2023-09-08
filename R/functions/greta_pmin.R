#' .. content for \description{} (no empty lines) ..
#'
#' equivalent to base::pmin(a, b), but for greta arrays and differentiable, and
#' works inside the greta.dynamics transition function
#'
#' @title
#' @param a
#' @param b
#' @return
#' @author Nick Golding
#' @export
greta_pmin <- function(a, b) {
  mask <- a > b
  # the following line is an utterly ridiculous way of computing 1-x, which
  # avoids defining a new data greta array inside the greta.dynamics transition
  # function, due to a bug in the shapes of tensors
  anti_mask <- exp(log1p(-mask))
  a * anti_mask + b * mask
}
