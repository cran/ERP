\name{erptest}
\alias{erptest}
\title{FDR- and FWER-controlling Multiple testing of ERP data}
\description{
Classical FDR- and FWER-controlling multiple testing procedures for ERP data in a linear model framework.
}
\usage{
erptest(dta, design, design0 = NULL, method = "BH", alpha = 0.05, pi0 = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dta}{Data frame containing the ERP curves: each column corresponds to a time frame and each row to a curve.}
  \item{design}{Design matrix of the full model for the relationship between the ERP and the experimental variables. Typically the output of the function
  model.matrix}
  \item{design0}{Design matrix of the null model. Typically a submodel of the full model, obtained by removing columns from design. Default is 
  NULL, corresponding to the model with no covariates.}
  \item{method}{FDR- or FWER- controlling multiple testing procedures as available in the function p.adjust. Default is "BH", 
  for the Benjamini-Hochberg procedure (see Benjamini and Hochberg, 1995).}
  \item{alpha}{The FDR or FWER control level. Default is 0.05}
  \item{pi0}{An estimate of the proportion of true null hypotheses, which can be plugged into an FDR controlling multiple testing procedure to improve
  its efficiency. Default is 1, corresponding to the classical FDR controlling procedures. If NULL, the proportion is estimated using 
  the function pval.estimate.eta0 of package fdrtool, with the method proposed by Storey and Tibshirani (2003).}
}
\details{The multiple testing corrections of the p-values are described in the help documentation of the function p.adjust du package stats.
}
\value{
\item{pval}{p-values of the tests.}
\item{correctedpval}{Corrected p-values, for the multiplicity of tests. Depends on the multiple testing method (see function p.adjust).}
\item{significant}{Indices of the time points for which the test is positive.}
\item{pi0}{Value for pi0: if the input argument pi0 is NULL, the output is the estimated proportion of true null hypotheses using the method 
by Storey and Tibshirani (2003).}
\item{test}{F-statistics.}
\item{df1}{Residual degrees of freedom for the full model.}
\item{df0}{Residual degrees of freedom for the null model.}
\item{signal}{Estimated signal: a pxT matrix, where p is the difference between the numbers of parameters in the full and null models and T the number of frames.}
\item{r2}{R-squared values for each of the T linear models.}
}
\references{
Benjamini, Y., and Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the Royal Statistical Society Series B, 289-300.

Benjamini, Y., and Yekutieli, D. (2001). The control of the false discovery rate in multiple testing under dependency. Annals of Statistics 29, 1165-1188.

Holm, S. (1979). A simple sequentially rejective multiple test procedure. Scandinavian Journal of Statistics, 6, 65-70.

Hommel, G. (1988). A stagewise rejective multiple test procedure based on a modified Bonferroni test. Biometrika, 75, 383-386.

Hochberg, Y. (1988). A sharper Bonferroni procedure for multiple tests of significance. Biometrika, 75, 800-803.

Shaffer, J. P. (1995). Multiple hypothesis testing. Annual Review of Psychology, 46, 561-576.

Sarkar, S. (1998). Some probability inequalities for ordered MTP2 random variables: a proof of Simes conjecture. Annals of Statistics, 26, 494-504.

Sarkar, S., and Chang, C. K. (1997). Simes' method for multiple hypothesis testing with positively dependent test statistics. Journal of the American Statistical Association, 92, 1601-1608.

Sheu, C.-F., Perthame, E., Lee Y.-S. and Causeur, D. (2016). Accounting for time dependence in large-scale multiple testing of event-related potential data. To appear in Annals of Applied Statistics.

Wright, S. P. (1992). Adjusted P-values for simultaneous inference. Biometrics, 48, 1005-1013.
}
\author{
David Causeur - david.causeur@agrocampus-ouest.fr
}
\seealso{
\code{\link{erpavetest}}, \code{\link{erpfatest}}, \code{\link{gbtest}}, \code{\link[stats]{p.adjust}}, \code{\link[fdrtool]{pval.estimate.eta0}}
}
\examples{
require(mnormt)
require(fdrtool)

data(erpcz)
data(simerp)

# Paired comparison of ERP curves

tests = erptest(erpcz[,1:251],design=model.matrix(~Subject+Instruction,data=erpcz),
   design0=model.matrix(~Subject,data=erpcz))

frames = seq(0,1001,4)
plot(frames,tests$signal,type="l",xlab="Time (ms)",
   ylab="Difference ERP curves")
points(frames[tests$significant],rep(0,length(tests$significant)),
   pch=16,col="blue")
title("Paired comparison at electrode CZ")

# Independent two-group comparison of ERP curves

tests = erptest(erpcz[,1:251],design=model.matrix(~Instruction,data=erpcz))

frames = seq(0,1001,4)
plot(frames,tests$signal,type="l",xlab="Time (ms)",
   ylab="Difference ERP curves")
points(frames[tests$significant],rep(0,length(tests$significant)),
   pch=16,col="blue")
title("Independent comparison at electrode CZ")

# Tests for significance of correlations

tests = erptest(simerp[,1:251],design=model.matrix(~y,data=simerp))
plot(frames,sign(tests$signal)*sqrt(tests$r2),type="l",xlab="Time (ms)",
   ylab="Correlation",ylim=c(-1,1))
points(frames[tests$significant],rep(-1,length(tests$significant)),
   pch=16,col="blue")
title("Simulation")

}
\keyword{ERP data}
\keyword{FDR}
\keyword{Multiple testing}
