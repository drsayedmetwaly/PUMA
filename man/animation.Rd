\name{animation}
\alias{animation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
animation
}
\description{
Makes 3D spining rotations of the PUMA visualization and can save then as animated GIF files. Animation can be stopped by pressing the ESC key. If the PUMA visualization window is closed, it will not be possible to make animations or take snapshots.
}
\usage{
animation(spin_axisXYZ = c(0, 0, 1), spin_duration = 30, spin_fps = 20, spin_startTime = 0, output3DFileName = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{spin_axisXYZ}{
The axis/axes of spining. By default the Z axis is used.
}
  \item{spin_duration}{
The duration of the spining animation. By default it is 30 seconds.
}
  \item{spin_fps}{
Animation frames per seconds. By default it is 20 frames per second.
}
  \item{spin_startTime}{
Spin start time. By default it is 0.
}
  \item{output3DFileName}{
The name of the output GIF animation. It must not contain spaces, path nor the file extension (.gif). If no output file name is entered, the animation will be done on screen only and will not be saved.
}
}

\references{
\url{https://github.com/drsayedmetwaly/PUMA}
}
\author{
Sayed Metwaly, \email{sayed.metwaly@abdn.ac.uk}
}
\note{
Installation of ImageMagick program is required to save PUMA animations as animated GIF files.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{PUMA}}, \code{\link{snapshot}}
}
\examples{
## Shows a 3D spining animation of the current PUMA visualization and save it in an animated GIF file.

animation("PUMA3D")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
