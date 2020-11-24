\name{snapshot}
\alias{snapshot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{

snapshot
}
\description{
Snapshot can be used to take a snapshot of the current PUMA 3D view and save it as a PNG image.
}
\usage{
snapshot(output3DFileName)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{output3DFileName}{
The path for the PNG image to be saved. It should be provided as a character vector: "C:/Datasets/puma3d.png"
}
}

\references{
\url{https://github.com/drsayedmetwaly/PUMA}
}
\author{
Sayed Metwaly, \email{sayed.metwaly@abdn.ac.uk}
}

\seealso{
\code{\link{PUMA}}
}
\examples{
## Save current PUMA 3D view in as a PNG image

snapshot("C:/PUMA/picture.png")

}

