\name{stackPlot}
\alias{stackPlot}
\title{Stack Plot}
\description{
  Stack plot of an object.  This function attempts to mimic a function
  called \code{stack.plot} in S+WAVELETS.  It is mostly a hacked version
  of \code{plot.ts} in \R.
}
\usage{stackPlot(x, plot.type = c("multiple", "single"), panel = lines, 
    log = "", col = par("col"), bg = NA, pch = par("pch"), cex = par("cex"),
    lty = par("lty"), lwd = par("lwd"), ann = par("ann"), xlab = "Time",
    main = NULL, oma = c(6, 0, 5, 0), layout = NULL,
    same.scale = 1:dim(x)[2], ...)
}
\arguments{
  \item{x}{\code{ts} object}
  \item{layout}{Doublet defining the dimension of the panel.  If not
    specified, the dimensions are chosen automatically.}
  \item{same.scale}{Vector the same length as the number of series to be
    plotted.  If not specified, all panels will have unique axes.}
  \item{plot.type,panel,log,col,bg,pch,cex,lty,lwd,ann,xlab,main,oma,...}{See \code{plot.ts}.}
}
%\value{}
\details{
  Produces a set of plots, one for each element (column) of \code{x}.
}
%\references{}
%\seealso{}
\author{Brandon Whitcher}
%\examples{}
\keyword{hplot}
