\name{truncMod}
\alias{truncMod}
\title{truncMod}
\description{truncMod fits a truncated negative bionomial model to data. Use to fit Species Richness (S\_i)}
\usage{
truncMod(trunc.form, trunc.pts, data, dist = 'NB', scale.covar=FALSE, est.var = TRUE, trace = TRUE)
}
\arguments{
\item{ trunc.form}{an object of class "formula" (or one that can be coerced to that class):a symbolic description of the model to be fitted. An offset in the form offset(log(area sampled)) should be included in the model formula to adjust for differences in area sampled between sites.}
\item{ trunc.pts}{A vector of truncation points of model. Must be of length equal to the numbers of rows(sites) of data. Set to N for fitting truncated species richness}
\item{ data}{a data frame containing the variables in the formula including the offset and the response. Each row should be a different site.}
\item{ dist}{one of either \code{NB} (negative binomial) or \code{poission}}
\item{ scale.covar}{Should the model matrix be scaled (TRUE/FALSE)? Useful if models do not converge. Setting this to TRUE will save the column means and sd for later prediction.}
\item{ est.var}{estimate the variance-covariance matrix using newton-raphson}
\item{ trace}{print model trace}
}
\details{truncMod used glm.nb (MASS) or glm to get initial parameters for untruncated likelihoood & then uses nlminb to fit truncated model. Allows calculation of se for overdispersion parameter.}
\value{
\item{ coef}{coefficents}
\item{ vcov}{variance-covariance matrix}
\item{ logl}{log-likelihood}
\item{ AIC}{AIC}
\item{ residuals}{model residuals, calculated using PIT}
\item{ fitted}{fitted values}
\item{ sds}{dfdf}
\item{ sums}{model sums}
\item{ mean.X}{column means for model matrix if scale.covar is TRUE}
\item{ sd.X}{column sds for model matrix if scale.covar is TRUE}
\item{ formula}{model formula}
}
\section{Authors}{Piers Dunstan and Scott Foster}
\examples{
fish.data$N.scale <- fish.data$N/fish.data$area
S.form <- S~1 + N.scale + depth + O2_AV + temp_AV + offset(log(area))
model.S <- truncMod(S.form,fish.data$N,data=fish.data,est.var=TRUE,scale.covar=TRUE)
plot(model.S$fitted,model.S$residuals)
abline(h=0)
}
\keyword{misc}
