\name{predict.S
}
\alias{predict.S
}
\title{predict.S
}
\description{Predict values for S and E(S) given a model for S and new covariates.
}
\usage{predict.S(object,new.obs,N=NA, offset=1,dist="NB",...)
}
\arguments{
\item{ object}{a fitted model from truncMod, of class S}
\item{ new.obs}{new covariates used to predict N. Must contain the same names as in the formula for truncMod}
\item{ N}{if S is conditional on N, then N is the predicted values returned from predict.N}
\item{ offset}{offset for sampling}
\item{ dist}{distribution used for S, either negative binomial ("NB") or poisson ("poisson")}
\item{ ... }{additional arguments to be passed to predict}
}
\details{
Uses mvtnorm to obtain new parameters for the glm using calculated variance-covariance matrix, assuming multivariate normal for parameters.

new.obs will be automatically scaled using details in model.S if scale.covar==TRUE
}
\value{
\item{ S}{predicted S}
\item{ expect.S}{expected S}
}
\section{Authors}{Piers Dunstan and Scott Foster
}
\examples{
N.form <- N~1  + depth + O2_AV + temp_AV + offset(log(area))
model.N <- negBinMod(N.form,data=fish.data,est.var=TRUE,scale.covar =FALSE)
fish.data$N.scale <- fish.data$N/fish.data$area
S.form <- S~1 + N.scale + depth + O2_AV + temp_AV + offset(log(area))
model.S <- truncMod(S.form,fish.data$N,data=fish.data,est.var=TRUE,scale.covar=TRUE)
new.obs <- as.data.frame(t(mean(fish.data)))
N.pred <- rep(NA,100)
S.pred <- rep(NA,100)
for(n in 1:100)\{
  N.pred[n] <- predict.N(model.N,new.obs,offset=10)$N
  # insert estimate for N into data frame for it is used, not the original value
  new.obs$N.scale <- N.pred[n]/10
  S.pred[n] <- predict.S(model.S,new.obs,N.pred[n],offset=10)$S
\}
hist(S.pred)
}
\keyword{misc}
