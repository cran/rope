#' A simulated data set for a scale-free network of 200 nodes
#'
#' The data set contains 175 observations for each node, the true network
#' structure dat was used to generate data and edge presence counts from glasso
#' over 100 bootstraps.
#'
#' @format A list containing:
#' \describe{
#'   \item{x}{A matrix of 175 observations (rows) for 200 variabels (columns)}
#'   \item{g}{The generating network structure (as a vector)}
#'   \item{B}{100, the number of bootstraps used when counting edge presence}
#'   \item{lambda}{The range of penalization used for glasso (the first 9
#'     generate U-shaped histograms)}
#'   \item{W}{A matrix of length(lambda) rows and 200*199/2 columns containing
#'     presence counts for each edge and each level of penalization}
#'   \item{Wlist}{A list of length(lamdba) containing matrices of size 200 by
#'     200, the data in W but in an alternative format}
#'   \item{gmatrix}{A 200 by 200 matrix, the data in g but in an alternative
#'     format}
#' }
'scalefree'
