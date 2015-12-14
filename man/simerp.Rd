\name{simerp}
\alias{simerp}
\docType{data}
\title{ERP dataset from a simulation study presented in detail in Sheu, Perthame, Lee and Causeur (2016)}
\description{
The dataset contains 20 ERP curves (time frames every 4 ms from 0 to 1000 ms) and a variable y for behavioral score.
}
\usage{data(simerp)}
\format{
  A sample of 20 simulated ERP curves, each curve lasting for one thousand milliseconds (ms) with one record per 4ms. The last column gives a behavioral score for the directed forgetting
  experiment described in Sheu, Perthame, Lee and Causeur (2016). 
}
\details{
The dataset is simulated under assumption of a multivariate regression model with normal residuals. The covariate, y, is a recognition score, as observed in a directed
forgetting experiment described in Sheu, Perthame, Lee and Causeur (2016). Only the slope coefficients corresponding to time points in the interval [450 ms ; 550 ms] are true non-nulls, with a peak at 500 ms. 
The residual covariance of the model mimics the sample residual covariance of the ERP datasets collected in the directed forgetting experiment.
}
\source{
Sheu, C.-F., Perthame, E., Lee Y.-S. and Causeur, D. 2016. Accounting for time dependence in large-scale multiple testing of event-related potential data. To appear in Annals of Applied Statistics.
}
\examples{
data(simerp)
erpplot(simerp[,1:251],frames=seq(0,1001,4),xlab="Time (ms)",
   ylab="ERP",main="Simulated ERP")
}
\keyword{datasets}
