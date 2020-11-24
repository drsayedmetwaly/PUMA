\name{PUMA}

\alias{PUMA}

\title{Probabilistic Universal Model Approximator (PUMA)}

\description{
A function for 3D visualization of decision boundary for classification algorithms. Classification models are built using MLR. PUMA plots a 3D scatter plot of the PCA scores. PUMA then acquires each possible point on the selected PCA coordinates (PCA grid) and reverse it to the original variables scale. The classification model is then utilized to predict the class of each PCA grid point and generates the probability for such assignment. A 3D contor is then projected using grid points with 0.5 probability to demarcate the decision boundary.

\usage{
PUMA(inputData, learner, task, measures = mmce, override_defaulat_ploting_PCs = NULL, cv = 10L, pointsize = 2, grp.cols = c("darkblue", "green", "darkred"), LearnerName = toupper(getLearnerShortName(learner)), contours3D, contours3Dcolors, engine3D = "rgl", material3D = "default", showTitles = TRUE, spin_axisXYZ = c(0, 0, 1), spin_duration = 30, spin_fps = 20, spin_startTime = 0, output3DFileName)
}

\arguments{
  \item{inputData}{
Data set with sample ID as row names, variables (metabolites) as columns and class column is the last column in the data set.
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
By default PUMA uses principle components PC1, PC2 and PC3 for plotting. If other PCs are required for plotting instead they can be entered here.
}
  \item{cv}{
The number of cross validations (max 10).
}
  \item{pointsize}{
The size of points (scores) in the 3D PCA plot.
}
  \item{grp.cols}{
Groups colours
}
  \item{LearnerName}{
By default LearnerName is the classification algorithm name. This is used in the plot title.
}
  \item{contours3D}{
The required numbers of decision boundaries. By default one decision boundary at 0.5 probability of class assignment is used.
}
  \item{contours3Dcolors}{
Colour of the 3D contours used to demarcate the decision boundary
}
  \item{engine3D}{
By default "rgl" engine is used. Other alternatives include "standard" and "grid".
}
  \item{material3D}{
This argument is required only for "standard" and "grid" engines 3D 3ngines. Possible values are the character strings "dull", "shiny", "metal", and "default"
}
  \item{showTitles}{
By default showTitles is ture and this shows the title and subtile in the 3D plot.
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

\details{
The built-in function \code{\link{loadcsv}} can import data sets from csv files and prepare them for use in PUMA. PUMA 3D plots are interactive and can be manipulated using the mouse. The function \code{\link{snapshot}} can be used to take a snapshot of the current view.
}

\references{
\url{https://github.com/drsayedmetwaly/PUMA}
}

\author{
Sayed Metwaly, \email{sayed.metwaly@@abdn.ac.uk}
}

\note{
PUMA depends on the following packages devtools, mlr, checkmate, Hmisc, rgl, misc3d and sinkr. The built-in function \code{\link{installdep}} can be used to install all the required packages before running PUMA for the first time. ImageMagick installation is required to save animations as animated GIF files
}

\seealso{
The built-in function \code{\link{loadcsv}} can import data sets from csv files and prepare them for use in PUMA. PUMA 3D plots are interactive and can be manipulated using the mouse. The function \code{\link{snapshot}} can be used to take a snapshot of the current view.
}

\examples{

## Example 1: Fit a PLS-DA model and plot it using PUMA

## Load libraries
library(PUMA)
library(mlr)

## Create an mlr task
trainTask <- makeClassifTask(data = PUMA_Example_Dataset1, target = "Class")

## Fit a PLS-DA model, plot it, do animation and save the output as an animated GIF file
PUMA(learner = "classif.plsdaCaret", measures=mmce, task = trainTask, inputData = PUMA_Example_Dataset1, cv=7, pointsize=2, grp.cols=c("orangered","blueviolet"), LearnerName="Partial Least Squares Discriminant Analysis (PLS-DA)", contours3D=1, contours3Dcolors="green", engine3D="rgl", material3D="default", showTitles=TRUE, spin_axisXYZ=c(0, 0, 1), spin_duration=30, spin_fps=20, spin_startTime=0, output3DFileName="PUMA3D")



## Example 2: Fit an RF model and plot it using PUMA

## Load libraries
library(PUMA)
library(mlr)

## Create an mlr task
trainTask <- makeClassifTask(data = PUMA_Example_Dataset2, target = "Class")

## Fit an RF model, plot it, do animation and save the output as an animated GIF file
PUMA(learner = "classif.randomForest", measures=mmce, task = trainTask, inputData = PUMA_Example_Dataset1, cv=7, pointsize=2, grp.cols=c("orangered","blueviolet"), LearnerName="Random Forest (RF)", contours3D=1, contours3Dcolors="green", engine3D="rgl", material3D="default", showTitles=TRUE, spin_axisXYZ=c(0, 0, 1), spin_duration=30, spin_fps=20, spin_startTime=0, output3DFileName="PUMA3D")



## Example 3: Fit an SVM (polynomial kernel) model and plot it using PUMA

## Load libraries
library(PUMA)
library(mlr)

## Create an mlr task
trainTask <- makeClassifTask(data = PUMA_Example_Dataset1, target = "Class")

## Fit an SVM (polynomial kernel) model, plot it, do animation and save the output as an animated GIF file
PUMA(learner = makeLearner("classif.svm", kernel = "polynomial"), measures=mmce, task = trainTask, inputData = PUMA_Example_Dataset1, cv=7, pointsize=2, grp.cols=c("orangered","blueviolet"), LearnerName="Support Vector Machine (SVM) with a polynomial kernel", contours3D=1, contours3Dcolors="green", engine3D="rgl", material3D="default", showTitles=TRUE, spin_axisXYZ=c(0, 0, 1), spin_duration=30, spin_fps=20, spin_startTime=0, output3DFileName="PUMA3D")


}


\keyword{ Probabilistic Universal Model Approximator (PUMA) }
\keyword{ PCA reversal }
\concept{ PCA Probabilistic mapping }


