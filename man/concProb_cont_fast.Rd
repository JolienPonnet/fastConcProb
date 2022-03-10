\name{concProb_cont_fast}

\alias{concProb_cont_fast}

\title{Fast concordance probability calculation in continuous setting.}

\description{Fast concordance probability calculation in continuous setting.}

\usage{concProb_cont_fast(obs, pred, w = NULL, nu = 0, tiesPred = FALSE, 
tiesObs = FALSE, verb = FALSE)}

\arguments{
  \item{obs}{Vector containing the observations. This is an obligatory argument,
  without default value.}
  \item{pred}{Vector containing the predictions. This is an obligatory argument,
  without default value.}
  \item{w}{Vector containing the weights for each observation-prediction pair. 
  Default weights are 1.}
  \item{nu}{Numeric vector of length 1 referring to the value that 
  the observations have to differ at least from one another in order to be 
  taken into account by the algorithm. Default value is 0.}
  \item{tiesPred}{Logical indicating whether the ties in the predictions should 
  be taken into account while calculating the concordance probability. 
  Default value is FALSE.}
  \item{tiesObs}{Logical indicating whether the ties in the observations should 
  be taken into account while calculating the concordance probability. 
  Default value is FALSE.}
  \item{verb}{Logical indicating whether some intermediate output should be 
  printed. Default value is FALSE.}
}

\value{
  \code{concProb_cont_fast} returns a list that contains the following items:
    \item{concProb}{the concordance probability, obtained following the supplied
    values for the arguments}
    \item{time}{the computing time of this function}
    \item{nb.conc}{the number of concordant pairs}
    \item{nb.disc}{the number of discordant pairs}
    \item{nb.ties.pred}{the number of ties in the predictions}
    \item{nb.ties.obs}{the number of ties in the observations}
    \item{nb.ties.pred.obs}{the number of ties in both predictions and observations}
}

\author{Jolien Ponnet, Jakob Raymaekers, Robin Van Oirbeek and Tim Verdonck}

\examples{
  library(fastConcProb)

  nSamples = 5000
  set.seed(123)
  obs   <- rnorm(nSamples)
  pred <- obs + rnorm(nSamples)

  f_cont_concProb <- concProb_cont_fast(obs, pred)
}
