\name{loadcsv}
\alias{loadcsv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
loadcsv(CSV_File, idCol = 1, classCol = 2, excludeColNo = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{CSV_File}{
%%     ~~Describe \code{CSV_File} here~~
}
  \item{idCol}{
%%     ~~Describe \code{idCol} here~~
}
  \item{classCol}{
%%     ~~Describe \code{classCol} here~~
}
  \item{excludeColNo}{
%%     ~~Describe \code{excludeColNo} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (CSV_File, idCol = 1, classCol = 2, excludeColNo = NULL) 
{
    Load_CSV <- read.csv(CSV_File, sep = ",", check.names = TRUE, 
        stringsAsFactors = FALSE, comment.char = "")
    if (is.null(excludeColNo)) {
        csv.data <- Load_CSV[, c(idCol, classCol, c(1:ncol(Load_CSV))[-c(idCol, 
            classCol)])]
    }
    else {
        csv.data <- Load_CSV[, c(idCol, classCol, c(1:ncol(Load_CSV))[-c(idCol, 
            classCol, excludeColNo)])]
    }
    rownames(csv.data) <- csv.data[, 1]
    dat <- csv.data[, c(3:ncol(csv.data))]
    dat$Class <- factor(csv.data[, 2])
    return(dat)
  }
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
