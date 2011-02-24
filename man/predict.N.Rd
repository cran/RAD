\name{predict.N
}
\alias{predict.N
}
\title{predict.N
}
\description{Predict values for N and E(N) given a model for N and new covariates.
}
\usage{predict.N(object,new.obs,offset=1,...)
}
\arguments{
\item{ object}{a fitted model from negBinMod, of class N.}
\item{ new.obs}{new covariates used to predict N. Must contain the same names as in the formula for negBinMod}
\item{ offset}{Predict sample for an area of this dimension}
\item{ ... }{additional arguments to be passed to predict}
}
\details{
Uses mvtnorm to obtain new parameters for the glm using calculated variance-covariance matrix, assuming multivariate normal for parameters.

new.obs will be automatically scaled using details in model.N if scale.covar==TRUE
}
\value{
\item{ N}{predicted N}
\item{ expect.N}{expected N}
}
\section{Authors}{Piers Dunstan and Scott Foster
}
\examples{
N.form <- N~1  + depth + O2_AV + temp_AV + offset(log(area))
model.N <- negBinMod(N.form,data=fish.data,est.var=TRUE,scale.covar =FALSE)
N.pred <- rep(NA,100)
for(n in 1:100) N.pred[n] <- predict.N(model.N,fish.data[5,],offset=10)$N
hist(N.pred)
}
\keyword{misc}
