\name{RoMLE.error}
\alias{SpatialRoMLE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Robust Maximum Likelihood Estimation for Spatial Error Model
}
\description{
Provides robust maximum likelihood estimation for spatial error model for symmetric weight matrix
}
\usage{
RoMLE.error(initial.beta,  initial.s2, initial.lambda, W, y, x, phi.function=NULL, converge.v=NULL, iter=NULL, print.values=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{initial.beta}{
  initial value of coefficients
  }
  \item{initial.s2}{
  initial value of sigma square
  }
  \item{initial.lambda}{
  initial value of lambda
  }
  \item{W}{
  a spatial weight matrix
  }
  \item{y}{
  dependent variable
  }
  \item{x}{
  independent variables
  }
  \item{phi.function}{
  a robust m-estimator function
  }
  \item{converge.v}{
  converge value for fisher scoring algorithm
  }
  \item{iter}{
  iteration number for fisher scoring algorithm
  }
  \item{print.values}{
  printing estimated values for each step until converge
  }
}
\details{
  \preformatted{phi.function} should be set as 1 for Cauchy, 2 for Welsch, 3 for Insha and 4 for Logistic
  \preformatted{converge.v} can be set as 1e-04
  \preformatted{iter} iter is set by users (e.g. 100)
  \preformatted{print.values} should be set TRUE or FALSE
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{coefficients }{coefficients values}
  \item{lambda }{lambda value}
  \item{s2 }{s2 value}
%% ...
}
\references{
Yildirim, V. and Kantar, Y.M. (2019). Robust estimation of spatial error model. in Journal of Statistical Computation and Simulation (in minor revision).
}
\author{
Vural Yildirim <vurall_yildirim@hotmail.com>

Yeliz Mert Kantar
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
#spdep library can be used for wieght matrix from listw
#require(spdep)
#W <- as(listw, "CsparseMatrix")

#Yildirim, V. and Kantar, Y.M. (2019). Robust estimation of spatial error model. in Journal of Statistical Computation and Simulation (in minor revision).

data(TRlistwmat)
data(TRQW)
data(TRunemploymentdata)
data(TRunemploymentcoefs)

y <- TRunemploymentdata$unemployment
x <- TRunemploymentdata$urbanization

#initial values was taken from MLE
initial.beta <- TRunemploymentcoefs[1:2,2]
initial.lambda <- TRunemploymentcoefs[3,2]
initial.s2 <- TRunemploymentcoefs[4,2]

RoMLE.error(initial.beta, initial.s2, initial.lambda, W=TRQW, y, x, phi.function=3, converge.v=0.0001, iter=100, print.values=TRUE)


#Yildirim, V , Mert Kantar, Y . (2019). Spatial Statistical Analysis of Participants in The Individual Pension System of Turkey. Eskisehir Teknik Universitesi Bilim Ve Teknoloji Dergisi B - Teorik Bilimler , 7 (2) , 184-194 . Doi: 10.20290/Estubtdb.518706

data(TRQW)
data(besdata)
data(bescoefs)

y <- besdata[,2]
x <- besdata[,3:9]

#initial values was taken from MLE
initial.beta <- bescoefs[1:8,2]
initial.lambda <- bescoefs[9,2]
initial.s2 <- bescoefs[10,2]

RoMLE.error(initial.beta, initial.s2, initial.lambda, W=TRQW, y, x, phi.function=3, converge.v=0.0001, iter=100, print.values=TRUE)

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{  }% use one of  RShowDoc("KEYWORDS")
\keyword{ }% use one of  RShowDoc("KEYWORDS")
\keyword{  }% __ONLY ONE__ keyword per line
\keyword{  }% __ONLY ONE__ keyword per line
