\name{erpplot}
\alias{erpplot}
\title{
Plot of ERP curves.
}
\description{
Wrapper for matplot (package graphics) to display ERP curves.
}
\usage{
erpplot(dta, frames = NULL, ylim = NULL, ...)
}
\arguments{
  \item{dta}{Data frame containing the ERP curves: each column corresponds to a time frame and each row to a curve.}
  \item{frames}{Sequence of time frames. Default is NULL, in which case frames is just the sequence of intergers between one and the total number of frames.}
  \item{ylim}{Limits for the y-axis. Default is NULL (set into the function).}
  \item{\dots}{Graphical parameters (see par) and any further arguments of plot, typically plot.default, may also be supplied as arguments to this function. 
  Hence, the high-level graphics control arguments described under par and the arguments to title may be supplied to this function.}
}
\value{The function generates a plot, but does not return any other numerical outputs.}
\author{
David Causeur (Agrocampus, Rennes, France), Ching-Fan Sheu (National Cheng-Kung University, Tainan, Taiwan) and Mei-Chen Chu (National Cheng-Kung University, Tainan, Taiwan) 
}
\seealso{\code{\link{matplot}}}
\examples{
data(erpcz)
erpplot(erpcz[,1:251],frames=seq(0,1001,4),xlab="Time (ms)",
   ylab="ERP",main="ERP at electrode CZ")

data(simerp)
erpplot(simerp[,1:251],frames=seq(0,1001,4),xlab="Time (ms)",
   ylab="ERP",main="Simulated ERP")
}
