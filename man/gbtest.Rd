\name{gbtest}
\alias{gbtest}
\title{The Guthrie-Buchwald procedure for significance analysis of ERP data}
\description{Monte-Carlo implementation of the Guthrie-Buchwald procedure (see Guthrie and Buchwald, 1991) which accounts for the
auto-correlation among test statistics to control erroneous detections of short intervals.}
\usage{
gbtest(dta, design, design0 = NULL, graphthresh = 0.05, nbsamples = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dta}{Data frame containing the ERP curves: each column corresponds to a time frame and each row to a curve.}
  \item{design}{Design matrix of the full model for the relationship between the ERP and the experimental variables. Typically the output of the function
  model.matrix}
  \item{design0}{Design matrix of the null model. Typically a submodel of the full model, obtained by removing columns from design. Default is 
  NULL, corresponding to the model with no covariates.}
  \item{graphthresh}{Graphical threshold (see Guthrie and Buchwald, 1991). Default is 0.05. As the FDR control level, the smaller is the graphical threshold, the more conservative is the
  procedure.}
  \item{nbsamples}{Number of samples in the Monte-Carlo method to estimate the residual covariance. Default is 1000.}
}
\details{
The Guthrie-Buchwald method starts from a preliminary estimation of r, the lag-1 autocorrelation, among test statistics.
Then, the null distribution of the lengths of the intervals $I_alpha = t : pvalue_t <= alpha $, where alpha is the
so-called graphical threshold parameter of the method, is obtained using simulations of p-values $p_t$ associated to
auto-regressive t-test process of order 1 with mean 0 and auto-correlation r. Such an interval $I_alpha$ is declared
significant if its length exceeds the $(1-alpha)-$quantile of the null distribution. Note that the former method is
designed to control erroneous detections of short significant intervals but not to control any type-I error rate.
}
\value{
\item{nbsignifintervals}{Number of significant intervals.}
\item{intervals}{List of length nbsignifintervals which components give the indices of each significant intervals.}
\item{significant}{Indices of the time points for which the test is positive.}
\item{signal}{Estimated signal: a pxT matrix, where p is the difference between the numbers of parameters in the full and null models and T the number of frames.}
\item{rho}{Estimated lag-1 auto-correlation.}
\item{r2}{R-squared values for each of the T linear models.}
}
\references{
Guthrie, D. and Buchwald, J.S. (1991). Significance testing of difference potentials. Psychophysiology, 28, 240-244.

Sheu, C.-F., Perthame, E., Lee Y.-S. and Causeur, D. (2016). Accounting for time dependence in large-scale multiple testing of event-related potentials data. To appear in Annals of Applied Statistics.
}
\author{
David Causeur - david.causeur@agrocampus-ouest.fr and Mei-Chen Chu (National Cheng-Kung University, Tainan, Taiwan)
}
\seealso{
\code{\link{erpavetest}}, \code{\link{erpfatest}}, \code{\link{erptest}}
}
\examples{
require(mnormt)
data(erpcz)
data(simerp)

# Paired comparison of ERP curves

tests = gbtest(erpcz[,1:251],design=model.matrix(~Subject+Instruction,data=erpcz),
   design0=model.matrix(~Subject,data=erpcz),nbsamples=500)

frames = seq(0,1001,4)
plot(frames,tests$signal,type="l",xlab="Time (ms)",
  ylab="Difference ERP curves")
points(frames[tests$significant],rep(0,length(tests$significant)),
   pch=16,col="blue")
title("Paired comparison at electrode CZ")

# Tests for significance of correlations

tests = gbtest(simerp[,1:251],design=model.matrix(~y,data=simerp),nbsamples=500)
plot(frames,sign(tests$signal)*sqrt(tests$r2),type="l",xlab="Time (ms)",
   ylab="Correlation",ylim=c(-1,1))
points(frames[tests$significant],rep(-1,length(tests$significant)),
   pch=16,col="blue")
title("Simulation")
}
\keyword{ ERP }
\keyword{ Guthrie-Buchwald procedure }
\keyword{ Multiple testing }
