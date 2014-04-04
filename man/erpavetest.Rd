\name{erpavetest}
\alias{erpavetest}
\title{
Significance testing of averaged ERPs. The entire ERP recording time is first partitioned into a pre-determined number
of equal intervals. The averaged ERPs for each time intervals are the input for analysis.
}
\description{
The function first calculates averaged ERP values within a predetermined number of equally-spaced intervals then tests
for significance of the relationship between averaged ERPs and covariates in a linear model framework.
}
\usage{
erpavetest(dta, design, design0 = NULL, nintervals = 10, method = "none", alpha = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dta}{Data frame containing the ERP curves: each column corresponds to a time frame and each row to a curve.}
  \item{design}{Design matrix of the full model for the relationship between the ERP and the experimental variables.
           Typically the output of the function model.matrix}
  \item{design0}{Design matrix of the null model. Typically a submodel of the full model, obtained by removing columns from design.
           Default is NULL, corresponding to the model with no covariates.}
  \item{nintervals}{Number of intervals in the partition of the whole interval of observation. Default is 10.}
  \item{method}{FDR- or FWER- controlling multiple testing procedures as available in the function p.adjust. Default is "none".}
  \item{alpha}{The FDR or FWER control level. Default is 0.05}
}
\value{
\item{pval}{p-values of the tests.}
\item{correctedpval}{Corrected p-values, for the multiplicity of tests. Depends on the multiple testing method (see function p.adjust).}
\item{significant}{Indices of the time points for which the test is positive.}
\item{segments}{Factor giving the membership of timepoints to each interval in the partition.}
\item{breaks}{Breakpoints of the partition.}
\item{test}{F-statistics.}
\item{df1}{Residual degrees of freedom for the full model.}
\item{df0}{Residual degrees of freedom for the null model.}
\item{signal}{Estimated signal: a pxT matrix, where p is the difference between the numbers of parameters in the full and null models and T the number of frames.}
\item{r2}{R-squared values for each of the T linear models.}
}
\author{
David Causeur - david.causeur@agrocampus-ouest.fr
}
\seealso{
\code{\link{erptest}}, \code{\link{erpfatest}}, \code{\link{gbtest}}, \code{\link[stats]{p.adjust}}, \code{\link[fdrtool]{pval.estimate.eta0}}
}
\examples{

require(mnormt)
require(fdrtool)

data(erpcz)
data(simerp)

# Paired t-tests for the comparison of ERP curves between two groups

tests = erpavetest(erpcz[,1:251],design=model.matrix(~Subject+Instruction,data=erpcz),
   design0=model.matrix(~Subject,data=erpcz))

frames = seq(0,1001,4)
plot(frames,tests$signal,type="l",xlab="Time (ms)",ylab="Difference ERP curves")
points(frames[tests$significant],rep(0,length(tests$significant)),pch=16,col="blue")
abline(v=frames[tests$breaks],lty=2,col="darkgray")
title("Paired comparison at electrode CZ")

# Tests for significance of correlations

tests = erpavetest(simerp[,1:251],design=model.matrix(~y,data=simerp))

plot(frames,sign(tests$signal)*sqrt(tests$r2),type="l",
   xlab="Time (ms)",ylab="Correlation",ylim=c(-1,1))
points(frames[tests$significant],rep(-1,length(tests$significant)),
   pch=16,col="blue")
abline(v=frames[tests$breaks],lty=2,col="darkgray")
title("Simulation")
}
\keyword{ERP data}
\keyword{FDR}
\keyword{Multiple testing}
