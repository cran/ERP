\name{erpcz}
\alias{erpcz}
\docType{data}
\title{ERP dataset from a directed forgetting experiment as presented in detail in Sheu, Perthame, Lee and Causeur (2016)}
\description{
The experiment was conducted using an experimental paradigm similar to that described in Lee, Lee and Fawcett (2013). The objective of the study was 
to examine the electrophysiological process associated with directed forgetting. 
}
\usage{data(erpcz)}
\format{
  40 ERP curves obtained at electrode CZ for 20 subjects in two experimental conditions (TBR and TBF), each curve lasting for one thousand milliseconds (ms) 
  with one record per 4ms. The last two columns, Subject and Instruction, provide the subject identifier and the group membership respectively.
}
\details{
In the study phase of the experiment, 20 participants were instructed 
(with a `+' or `X' cue, respectively) to either remember or forget a stimulus word that was presented on a computer screen. ERPs were recorded throughout 
each of 90 trials - half for to-be-remembered (TBR) and half for to-be-forgotten (TBF).  
ERPs from nine electrode positions - 3 each from frontal, central and posterior regions - during the study phase were analyzed. The curent dataset provides
the ERP curves at electrode CZ (central on the miline). The corresponding experimental covariates (subject and instruction) are provided in the last two columns.
}
\source{
Sheu, C.-F., Perthame, E., Lee Y.-S. and Causeur, D. 2016. Accounting for time dependence in large-scale multiple testing of event-related potential data. To appear in Annals of Applied Statistics.
}
\references{
Lee, Y-S., Lee, H.M., Fawcett, J.-M. (2013). Intentional forgetting reduces color-naming interference: evidence from item-method directed forgetting. 
Journal of Experimental Psychology. Learning, Memory and Cognition, 39 (1), 220-236.
}
\examples{
data(erpcz)
erpplot(erpcz[,1:251],frames=seq(0,1001,4),xlab="Time (ms)",
   ylab="ERP",main="ERP at electrode CZ")
}
\keyword{datasets}
