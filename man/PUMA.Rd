\name{Probabilistic Universal Model Approximator}

\alias{PUMA}

\title{Probabilistic Universal Model Approximator (PUMA)}

\description{
PUMA is a function for 3D visualization of the decision boundary of classification algorithms.
}

\usage{
PUMA(inputData, learner, task, measures = mmce, override_defaulat_ploting_PCs = NULL, cv = 10L, pointsize = 2, grp.cols = c("darkblue", "green", "darkred"), LearnerName = toupper(getLearnerShortName(learner)), contours3D, contours3Dcolors, engine3D = "rgl", material3D = "default", showTitles = TRUE, spin_axisXYZ = c(0, 0, 1), spin_duration = 30, spin_fps = 20, spin_startTime = 0, output3DFileName)
}

\arguments{
  \item{inputData}{
Data set with sample ID as row names, variables (metabolites) as columns and a "Class" column containing groups at the last column of the data set. Data sets in CSV files can be imported using the built-in function \code{\link{loadcsv}}.
}
  \item{learner}{
MLR learner. This parameter defines the model (classification algorithm) of interest.
}
  \item{task}{
MLR task. This parameter defines the data set and classes to be used by the classification model.
}
  \item{measures}{
MLR measures. This parameter defines the method used for assessing performance of the classification algorithm.
}
  \item{override_defaulat_ploting_PCs}{
By default PUMA uses principle components PC1, PC2 and PC3 for plotting. If other set of PCs are required for plotting instead they can be entered here.
}
  \item{cv}{
The number of cross validations (max 10).
}
  \item{pointsize}{
The size of points (scores) in the 3D PCA plot.
}
  \item{grp.cols}{
Groups colors
}
  \item{LearnerName}{
By default LearnerName is the classification algorithm name. This is used in the plot title.
}
  \item{contours3D}{
The required numbers of decision boundaries. By default one decision boundary at 0.5 probability of class assignment is used.
}
  \item{contours3Dcolors}{
Color of the 3D contours used to demarcate the decision boundary
}
  \item{engine3D}{
By default "rgl" engine is used. Other alternatives include "standard" and "grid".
}
  \item{material3D}{
This argument is required only for "standard" and "grid" engines 3D engines. Possible values are the character strings "dull", "shiny", "metal", and "default"
}
  \item{showTitles}{
By default showTitles is TRUE and this shows the title and subtitle in the 3D plot.
}
  \item{spin_axisXYZ}{
Defines the axes used for rotating the final PUMA 3D plot. Formula: c(X,Y,Z) for each axis 0=No, 1=Yes.
}
  \item{spin_duration}{
If duration = 0 the 3D output will be a still .png picture
}
  \item{spin_fps}{
The number of (3D animation) frames/seconds to be used for building the animation movie
}
  \item{spin_startTime}{
Start time of spinning of the 3D plot
}
  \item{output3DFileName}{
The file to be used for saving 3D animation
}
}

\details{
Classification models are initially fitted using the versatile MLR library. PCA transformation of the data set is performed and a 3D scatter plot of the PCA scores is plotted using the principal components (PCs) of choice. PUMA then calculates the coordinates of all the points present in the 3D PCA space (PCA 3D grid). Matrix manipulation is carried out to reverse all PCA 3D grid points from the principal components (PCs) domain back to the original variables (metabolites) scale. PUMA then uses the trained classification model to predict the class assignment of the reversed PCA 3D grid points and estimate the probability of class assignment at each point. Once completed, this step has effectively generated a probabilistic map containing the probability of class assignment at each point of the entire PCA 3D space. PUMA finally overlays a 3D contour to join points in the 3D PCA space that have 0.5 probability of class assignment and as such demarcate the decision boundary of the classification algorithm.
}

\references{
\url{https://github.com/drsayedmetwaly/PUMA}
}

\author{
Sayed Metwaly, \email{sayed.metwaly@abdn.ac.uk}
}

\note{
PUMA depends on the following packages devtools, mlr, checkmate, Hmisc, rgl, misc3d and sinkr. The built-in function installdep() can be used to install all the required dependency packages before running PUMA for the first time. Installation of ImageMagick program is required to save PUMA animations as animated GIF files.
The built-in function \code{\link{loadcsv}} can import data sets from CSV files and prepare them for direct use in PUMA. PUMA 3D plots are interactive and can be manipulated using the mouse. The built-in function \code{\link{snapshot}} can be used to take a snapshot of the current view. The function \code{\link{animation}} makes 3D spining rotations of the PUMA visualization and can save them as animated GIF files. Installation of ImageMagick program is required to save PUMA animations as animated GIF files. If the PUMA visualization window is closed, it will not be possible to make animations or take snapshots.

}

\seealso{
\code{\link{loadcsv}}, \code{\link{snapshot}}, \code{\link{animation}}
}

\examples{

## Example 1: Fit a PLS-DA model and plot it using PUMA

library(PUMA)
library(mlr)

trainTask <- makeClassifTask(data = PUMA_Example_Dataset1, target = "Class")

PUMA(learner = "classif.plsdaCaret", measures=mmce, task = trainTask, inputData = PUMA_Example_Dataset1, cv=7, pointsize=2, grp.cols=c("orangered","blueviolet"), LearnerName="Partial Least Squares Discriminant Analysis (PLS-DA)", contours3D=1, contours3Dcolors="green", engine3D="rgl", material3D="default", showTitles=TRUE, spin_axisXYZ=c(0, 0, 1), spin_duration=30, spin_fps=20, spin_startTime=0, output3DFileName="PUMA3D")



## Example 2: Fit an RF model and plot it using PUMA

library(PUMA)
library(mlr)

trainTask <- makeClassifTask(data = PUMA_Example_Dataset2, target = "Class")

PUMA(learner = "classif.randomForest", measures=mmce, task = trainTask, inputData = PUMA_Example_Dataset2, cv=7, pointsize=2, grp.cols=c("orangered","blueviolet"), LearnerName="Random Forest (RF)", contours3D=1, contours3Dcolors="green", engine3D="rgl", material3D="default", showTitles=TRUE, spin_axisXYZ=c(0, 0, 1), spin_duration=30, spin_fps=20, spin_startTime=0, output3DFileName="PUMA3D")



## Example 3: Fit an SVM (polynomial kernel) model and plot it using PUMA

library(PUMA)
library(mlr)

trainTask <- makeClassifTask(data = PUMA_Example_Dataset1, target = "Class")

PUMA(learner = makeLearner("classif.svm", kernel = "polynomial"), measures=mmce, task = trainTask, inputData = PUMA_Example_Dataset1, cv=7, pointsize=2, grp.cols=c("orangered","blueviolet"), LearnerName="Support Vector Machine (SVM) with a polynomial kernel", contours3D=1, contours3Dcolors="green", engine3D="rgl", material3D="default", showTitles=TRUE, spin_axisXYZ=c(0, 0, 1), spin_duration=30, spin_fps=20, spin_startTime=0, output3DFileName="PUMA3D")

}


\keyword{ PUMA }
\keyword{ PCA reversal }
\concept{ PCA Probabilistic mapping }


