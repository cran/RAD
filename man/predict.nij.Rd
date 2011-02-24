\name{predict.nij
}
\alias{predict.nij
}
\title{predict.nij
}
\description{Predict values for the abundances of species 1 to S and the derivative of the RAD curve at species rank 1.
}
\usage{predict.nij(object,new.obs,N=NA,S=NA,...)
}
\arguments{
\item{ object}{a fitted model from negBinMod, of class nij.}
\item{ new.obs}{new covariates used to predict nij. Must contain the same names as in the formula for negBinMod}
\item{ N}{N is the predicted values returned from predict.N}
\item{ S}{S is the predicted values returned from predict.S}
\item{ ... }{additional arguments to be passed to predict}
}
\details{
Uses mvtnorm to obtain new parameters for the glm using calculated variance-covariance matrix, assuming multivariate normal for parameters.

new.obs will be automatically scaled using details in model.nij if scale.covar==TRUE
}
\value{
\item{ deriv.eta}{derivative of the RAD curve at species rank 1}
\item{ nij}{vector of predicted abundances}
}
\section{Authors}{Piers Dunstan and Scott Foster
}
\examples{
N.form <- N~1  + depth + O2_AV + temp_AV + offset(log(area))
model.N <- negBinMod(N.form,data=fish.data,est.var=TRUE,scale.covar =FALSE)
fish.data$N.scale <- fish.data$N/fish.data$area
S.form <- S~1 + N.scale + depth + O2_AV + temp_AV + offset(log(area))
model.S <- truncMod(S.form,fish.data$N,data=fish.data,est.var=TRUE,scale.covar=TRUE)
n.data$N.scale <- n.data$N/n.data$area
n.data$S.scale <- n.data$S/n.data$area
nij.form <- nij~1 + N.scale + S.scale  + depth + O2_AV + temp_AV
model.nij <- MDMnMod(nij.form, data=n.data, ID=n.data$i, dist="MDMn", scale.covar=TRUE, est.var=TRUE, calc.resid=TRUE, trace=TRUE)
new.obs <- as.data.frame(t(mean(fish.data)))
N.pred <- rep(NA,100)
S.pred <- rep(NA,100)
eta.pred <- rep(NA,100)
for(n in 1:100)\{
  N.pred[n] <- predict.N(model.N,new.obs,offset=10)$N
  # insert estimate for N into data frame for it is used, not the original value
  new.obs$N.scale <- N.pred[n]/10
  S.pred[n] <- predict.S(model.S,new.obs,N.pred[n],offset=10)$S
  new.obs$S.scale <- S.pred[n]/10
  eta.pred[n] <- unlist(predict.nij(model.nij,new.obs,N.pred[n],S.pred[n])$deriv.eta)
\}
hist(eta.pred)
}
\keyword{misc}
