\name{loadcsv}
\alias{loadcsv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{

loadcsv

}
\description{
Formats data sets in CSV files into PUMA data format.
}
\usage{
loadcsv(CSV_File, idCol = 1, classCol = 2, excludeColNo = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{CSV_File}{
The path for the CSV file written as a character vector: "C:/Datasets/mydataset.csv"
}
  \item{idCol}{
The number of the column that contains sample IDs. By default, loadcsv assumes ID column is the first column in the CSV file.
}
  \item{classCol}{
The number of the column that contains classes. By default, loadcsv assumes the column that contain classes is the second column in the CSV file.
}
  \item{excludeColNo}{
The number of the columns to be excluded (optional).
}
}
\details{
PUMA data sets have specific format with observations (samples) in rows and variables (metabolites)in columns. The row names are the sample IDs and the last column in the data set should be named "Class" and should contain the classes. All variable (metabolite) columns must have numeric values.
}
\references{
\url{https://github.com/drsayedmetwaly/PUMA}
}
\author{
Sayed Metwaly, \email{sayed.metwaly@@abdn.ac.uk}
}


\seealso{
\code{\link{PUMA}}
}
\examples{
## Import a data set from CSV file

library(PUMA)

dat<-loadcsv("C:/PUMA/Dataset1.csv")

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
