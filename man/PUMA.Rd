\name{Probabilistic Universal Model Approximator}

\alias{PUMA}

\title{Probabilistic Universal Model Approximator (PUMA)}

\description{
PUMA is a function for 2D/3D visualization of the decision boundary of classification algorithms.
}

\usage{
PUMA(inputData, learner, task, measures=mmce, override_default_ploting_PCs=NULL,
cv = 7L, plot_dimensions=3, pointsize = 2, grp.cols = c("darkblue", "green", "darkred"),
LearnerName=toupper(getLearnerShortName(learner)), contours3D, contours3Dcolors,
engine3D="rgl", material3D="default", showTitles3D=TRUE, spin3D_axisXYZ = c(0, 0, 1),
spin3D_duration = 30, spin3D_fps=20, spin3D_startTime = 0, output3DFileName = NULL,
prob.alpha2D = TRUE, err.mark2D = "cv", err.col2D = "black", err.size2D = 4,
err.shape2D=4)
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
By default PUMA uses principle components "PC1", "PC2" and "PC3" for plotting. If other set of PCs are required for plotting instead they can be entered here.

if "maxd2" is entered, PUMA will objectively measure class separation in different combinations of PCs by means of Mahalanobis distance and will automatically detremine the otpimal PCs (have the maximal class separation) to use for plotting.
}
  \item{cv}{
The number of cross validations (max 10).
}
  \item{plot_dimensions}{
The dimensions of the output plot. This can be 2 (for a 2D plot) or 3 (for a 3D plot).
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
  \item{showTitles3D}{
By default showTitles3D is TRUE and this shows the title and subtitle in the 3D plot.
}
  \item{spin3D_axisXYZ}{
Defines the axes used for rotating the final PUMA 3D plot. Formula: c(X,Y,Z) for each axis 0=No, 1=Yes.
}
  \item{spin3D_duration}{
If duration = 0 the 3D output will be a still .png picture
}
  \item{spin3D_fps}{
The number of (3D animation) frames/seconds to be used for building the animation movie
}
  \item{spin3D_startTime}{
Start time of spinning of the 3D plot
}
  \item{output3DFileName}{
The file to be used for saving 3D animation
}


  \item{prob.alpha2D}{
Change background transparency according to class probability in 2D plots
}

  \item{err.mark2D}{
Controls the plotting of error marks. It can be "cv" (default value, this mark samples misclassified in cross-validation), "train" (to mark samples misclassified in model training) or "none" (shows no error marks).
}

  \item{err.col2D}{
Controls the colour of error marks (default value is "black").
}

  \item{err.size2D}{
Controls the size of error marks (default value is 4).
}

  \item{err.shape2D}{
Controls the shape of error marks (default value is 4 which is X).
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

## Example 1: Fit a PLS-DA model and plot it in 2D using PUMA

library(PUMA)
library(mlr)

trainTask <- makeClassifTask(data = PUMA_Example_Dataset1, target = "Class")

PUMA(learner = "classif.plsdaCaret", measures=mmce, task = trainTask,
     inputData = PUMA_Example_Dataset1,
     override_default_ploting_PCs = "maxd2", cv=7, plot_dimensions=2, pointsize=2,
     grp.cols=c("orangered","blueviolet"),
     LearnerName="Partial Least Squares Discriminant Analysis (PLS-DA)",
     contours3D=1, contours3Dcolors="green", engine3D="rgl", material3D="default",
     showTitles3D=TRUE, spin3D_axisXYZ=c(0, 0, 1), spin3D_duration=30, spin3D_fps=20,
     spin3D_startTime=0, output3DFileName="PUMA3D",
     prob.alpha2D = TRUE,
     err.mark2D = "cv",
     err.col2D = "black", err.size2D=4, err.shape2D=4)


## Example 2: Fit an RF model and plot it in 3D using PUMA

library(PUMA)
library(mlr)

trainTask <- makeClassifTask(data = PUMA_Example_Dataset2, target = "Class")

PUMA(learner = "classif.randomForest", measures=mmce, task = trainTask,
     inputData = PUMA_Example_Dataset2,
     override_default_ploting_PCs = "maxd2", cv=7, plot_dimensions=3, pointsize=2,
     grp.cols=c("orangered","blueviolet"),
     LearnerName="Partial Least Squares Discriminant Analysis (PLS-DA)",
     contours3D=1, contours3Dcolors="green", engine3D="rgl", material3D="default",
     showTitles3D=TRUE, spin3D_axisXYZ=c(0, 0, 1), spin3D_duration=30, spin3D_fps=20,
     spin3D_startTime=0, output3DFileName="PUMA3D",
     prob.alpha2D = TRUE,
     err.mark2D = "cv",
     err.col2D = "black", err.size2D=4, err.shape2D=4)


## Example 3: Fit an SVM (polynomial kernel) model and plot it in 3D using PUMA

library(PUMA)
library(mlr)

trainTask <- makeClassifTask(data = PUMA_Example_Dataset1, target = "Class")

PUMA(learner = makeLearner("classif.svm", kernel = "polynomial"), measures=mmce,
task = trainTask, inputData = PUMA_Example_Dataset2,
     override_default_ploting_PCs = "maxd2", cv=7, plot_dimensions=3, pointsize=2,
     grp.cols=c("orangered","blueviolet"),
     LearnerName="Partial Least Squares Discriminant Analysis (PLS-DA)",
     contours3D=1, contours3Dcolors="green", engine3D="rgl", material3D="default",
     showTitles3D=TRUE, spin3D_axisXYZ=c(0, 0, 1), spin3D_duration=30, spin3D_fps=20,
     spin3D_startTime=0, output3DFileName="PUMA3D",
     prob.alpha2D = TRUE,
     err.mark2D = "cv",
     err.col2D = "black", err.size2D=4, err.shape2D=4)


}


\keyword{ PUMA }
\keyword{ PCA reversal }
\concept{ PCA Probabilistic mapping }


