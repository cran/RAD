\name{MDMnMod}
\alias{MDMnMod}
\title{MDMnMod}
\description{MDMnMod}
\usage{
MDMnMod(MDMn.form, data, ID, dist="MDMn", scale.covar=FALSE, est.var=TRUE, calc.resid=TRUE, trace=TRUE)
}
\arguments{
\item{ MDMn.form}{an object of class "formula" (or one that can be coerced to that class):a symbolic description of the model to be fitted. Must include both the response and covariates}
\item{ data}{a data frame containing the variables in the formula. Must contain a column with the ordered (most to least) abundances for each site. Covariates must be repeated for each row belonging to a site.}
\item{ ID}{Vector of site identifiers for all observations. For every abundance (row) that comes from a site, the same identifier must be used. The length of ID must equal nrows(data)}
\item{ dist}{one of either multinomial "Mn" or Dirichlet Multinomial "DMN" or Modified Dirichlet Multinomial "MDMn" (default)}
\item{ scale.covar}{Should the model matrix be scaled (TRUE/FALSE)? Useful if models do not converge. Setting this to TRUE will save the column means and sd for later prediction.}
\item{ est.var}{estimate the variance-covariance matrix using newton-raphson}
\item{ calc.resid}{should residuals be calculated?}
\item{ trace}{print model trace}
}
\details{Fits a selected distribution to the vector of relative abundances. See Foster and Dunstan 2009 for details.}
\value{
\item{ coef}{coefficents}
\item{ vcov}{variance-covariance matrix}
\item{ logl}{log-likelihood}
\item{ AIC}{AIC}
\item{ residuals}{model residuals, calculated using PIT}
\item{ fitted}{matrix of fitted values (nij and pij)}
\item{ mean.X}{column means for model matrix if scale.covar is TRUE}
\item{ sd.X}{column sds for model matrix if scale.covar is TRUE}
\item{ formula}{model formula}
}
\section{Authors}{Piers Dunstan and Scott Foster}
\examples{
head(n.data)
## format for the data object
n.data$N.scale <- n.data$N/n.data$area
n.data$S.scale <- n.data$S/n.data$area
nij.form <- nij~1 + N.scale + S.scale  + depth + O2_AV + temp_AV
model.nij <- MDMnMod(nij.form, data=n.data, ID=n.data$i, dist="MDMn", scale.covar=TRUE, est.var=FALSE, calc.resid=TRUE, trace=TRUE)
plot(model.nij$fitted$nij,model.nij$residuals)
plot(log(model.nij$fitted$nij),model.nij$residuals)
}
\keyword{misc}
