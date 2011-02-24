\name{negBinMod}
\alias{negBinMod}
\title{negBinMod}
\description{negBinMod fits a negative bionomial model to data. Use to fit Total Abundance (Ni)}
\usage{
negBinMod(NB.form, data, est.var=TRUE, scale.covar=FALSE , trace=TRUE)
}
\arguments{
\item{ NB.form}{an object of class "formula" (or one that can be coerced to that class):a symbolic description of the model to be fitted. An offset in the form offset(log(area sampled)) should be included in the model formula to adjust for differences in area sampled between sites.}
\item{ data}{a data frame containing the all variables in the formula including the offset and the response. Each row should be a different site.}
\item{ scale.covar}{Should the model matrix be scaled (TRUE/FALSE)? Useful if models do not converge. Setting this to TRUE will save the column means and sd for later prediction.}
\item{ est.var}{estimate the variance-covariance matrix using newton-raphson}
\item{ trace}{print model trace}
}
\details{negBinMod used glm.nb (MASS) to get initial parameters & then uses nlminb to fit model. Allows calculation of se for overdispersion parameter. Output is otherwise identical to glm.nb}
\value{
\item{ coef}{coefficents}
\item{ vcov}{variance-covariance matrix}
\item{ logl}{log-likelihood}
\item{ AIC}{AIC}
\item{ residuals}{model residuals, calculated using PIT}
\item{ fitted}{fitted values}
\item{ mean.X}{column means for model matrix if scale.covar is TRUE}
\item{ sd.X}{column sds for model matrix if scale.covar is TRUE}
\item{ formula}{model formula}
}
\section{Authors}{Piers Dunstan and Scott Foster}
\examples{
N.form <- N~1  + depth + O2_AV + temp_AV + offset(log(area))
model.N <- negBinMod(N.form,data=fish.data,est.var=TRUE,scale.covar =FALSE)
plot(model.N$fitted,model.N$residuals)
abline(h=0)
}
\keyword{misc}
