\name{Probabilistic Universal Model Approximator}

\alias{PUMA}

\title{Probabilistic Universal Model Approximator (PUMA)}

\description{
PUMA is a function for 2D/3D visualization of the decision boundary of classification algorithms.
}

\usage{
PUMA(inputData, targetColumn="Class",

                   model,
                   model_crossValid = 7L,
                   model_perfMetric = mmce,
                   model_name=toupper(getLearnerShortName(model)),

                   plot_dimensions = 3,
                   plot_PCs = "default",
                   plot_pointSize = 2,
                   plot_groupColors = c("darkblue", "green", "darkred"),
                   plot_dBoundary.level = 0.5,
                   plot_dBoundary.color = "black",
                   plot_dBoundary.size = 1,
                   plot_dBoundary.alpha = 0.5,

                   option3D_engine = "rgl",
                   option3D_material = "default",
                   option3D_showTitles = TRUE,
                   option3D_spin.axisXYZ = c(0, 0, 1),
                   option3D_spin.duration = 30,
                   option3D_spin.fps = 20,
                   option3D_spin.startTime = 0,
                   option3D_outputFileName = NULL,

                   option2D_bgProbAlpha = TRUE,
                   option2D_err.mark = "cv",
                   option2D_err.color = "black",
                   option2D_err.pointSize = 4,
                   option2D_err.pointShape = 4)
}

\arguments{
  \item{inputData}{
Data set with sample ID as row names, variables (metabolites) as columns and a "Class" column containing groups at the last column of the data set. Data sets in CSV files can be imported using the built-in function \code{\link{loadcsv}}.
}

  \item{targetColumn}{
This argument identifies the name of data set column that contains the classes.
}

  \item{model}{
MLR learner. This parameter defines the model (classification algorithm) of interest.
}

  \item{model_crossValid}{
The number of cross validations (max 10).
}

  \item{model_perfMetric}{
MLR measures. This parameter defines the method used for assessing performance of the classification algorithm.
}

  \item{model_name}{
By default model_name is the classification algorithm name. This parameter sets the title of the model to show in the final plot.
}

  \item{plot_dimensions}{
The dimensions of the output plot. This can be 2 (for a 2D plot) or 3 (for a 3D plot).
}

  \item{plot_PCs}{
By default PUMA uses principle components "PC1", "PC2" and "PC3" for plotting. If the user requires other set of PCs for plotting instead they can be entered here.

if "maxd2" is entered, PUMA will objectively measure class separation in different combinations of PCs by means of Mahalanobis distance and will automatically detremine the otpimal PCs (have the maximal class separation) to use for plotting.
}

  \item{plot_pointSize}{
The size of observation points (scores) in the 2D/3D PUMA plot.
}

  \item{plot_groupColors}{
Groups colors
}


  \item{decisionBoundary.level}{
By default one decision boundary at 0.5 probability of class assignment is used. Note: The probability level of more than one contour can be entered e.g decisionBoundary_Level = c(0.3, 0.5, 0.8)
}
  \item{decisionBoundary.color}{
Color of the contours used to demarcate the decision boundary and other probability levels.
}
  \item{decisionBoundary.size}{
Size of the contours used to demarcate the decision boundary and other probability levels.
}
  \item{decisionBoundary.alpha}{
transparency of the contours used to demarcate the decision boundary and other probability levels.
}
  \item{option3D_engine}{
By default "rgl" engine is used. Other alternatives include "standard" and "grid".
}
  \item{option3D_material}{
This argument is required only for "standard" and "grid" engines 3D engines. Possible values are the character strings "dull", "shiny", "metal", and "default"
}
  \item{option3D_showTitles}{
By default option3D_showTitles is TRUE and this shows the title and subtitle in the 3D plot.
}
  \item{option3D_spin.axisXYZ}{
Defines the axes used for rotating the final PUMA 3D plot. Formula: c(X,Y,Z) for each axis 0=No, 1=Yes.
}
  \item{option3D_spin.duration}{
If duration = 0 the 3D output will be a still .png picture
}
  \item{option3D_spin.fps}{
The number of (3D animation) frames/seconds to be used for building the animation movie
}
  \item{option3D_spin.startTime}{
Start time of spinning of the 3D plot
}
  \item{option3D_outputFileName}{
The file to be used for saving 3D animation
}


  \item{option2D_bgProbAlpha}{
Change background transparency according to class probability in 2D plots
}

  \item{option2D_err.mark}{
Controls the plotting of error marks. It can be "cv" (default value, this mark samples misclassified in cross-validation), "train" (to mark samples misclassified in model training) or "none" (shows no error marks).
}

  \item{option2D_err.color}{
Controls the colour of error marks (default value is "black").
}

  \item{option2D_err.pointSize}{
Controls the size of error marks (default value is 4).
}

  \item{option2D_err.pointShape}{
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

## Example 1: Fit a PLS-DA model on mice data set and plot it in 2D using PUMA

## Load libraries
library(PUMA)

## Fit PLS-DA model, plot it in 2D
PUMA(
     inputData = PUMA::mice, targetColumn="Class",

     model = "classif.plsdaCaret",
     model_perfMetric=mmce, model_crossValid=7,
     model_name="Partial Least Squares Discriminant Analysis (PLS-DA)",

     plot_dimensions=2, plot_pointSize=2, plot_PCs = "maxd2",
     plot_groupColors=c("orangered","blueviolet"),
     plot_dBoundary.level=0.5,
     plot_dBoundary.color="black",
     plot_dBoundary.size=1,
     plot_dBoundary.alpha=1,

     option3D_engine="rgl", option3D_material="default", option3D_showTitles=TRUE,
     option3D_spin.axisXYZ=c(0, 0, 1), option3D_spin.duration=30,
     option3D_spin.fps=20, option3D_spin.startTime=0,
     option3D_outputFileName="PUMA3D",

     option2D_bgProbAlpha = TRUE,
     option2D_err.mark = "cv", option2D_err.color = "black",
     option2D_err.pointSize=4, option2D_err.pointShape=4
     )



## Example 2: Fit an RF model on mice data set and plot it in 3D using PUMA

## Load libraries
library(PUMA)

## Fit a RF model, plot it, do animation and save the output as an animated GIF file
PUMA(
     inputData = PUMA::mice, targetColumn="Class",

     model = "classif.randomForest",
     model_perfMetric=mmce, model_crossValid=7,
     model_name="Random Forest (RF)",

     plot_dimensions=3, plot_pointSize=2, plot_PCs = "maxd2",
     plot_groupColors=c("orangered","blueviolet"),
     plot_dBoundary.level=0.5,
     plot_dBoundary.color="black",
     plot_dBoundary.size=1,
     plot_dBoundary.alpha=0.5,

     option3D_engine="rgl", option3D_material="default", option3D_showTitles=TRUE,
     option3D_spin.axisXYZ=c(0, 0, 1), option3D_spin.duration=30,
     option3D_spin.fps=20, option3D_spin.startTime=0,
     option3D_outputFileName="PUMA3D",

     option2D_bgProbAlpha = TRUE,
     option2D_err.mark = "cv", option2D_err.color = "black",
     option2D_err.pointSize=4, option2D_err.pointShape=4
     )



## Example 3: Fit an SVM (polynomial kernel) model on mice data set and plot it using PUMA

## Load libraries
library(PUMA)

## Fit a 3D SVM (polynomial kernel) model, plot decision boundary (at .50 probabilty) and other contours at 0.025 and 0.975 probabilities, do animation and save the output as an animated GIF file
PUMA(
     inputData = PUMA::mice, targetColumn="Class",
     model = makeLearner("classif.svm", kernel = "polynomial"),
     model_perfMetric=mmce, model_crossValid=7,
     model_name="Support Vector Machine (SVM) with a polynomial kernel",

     plot_dimensions=3, plot_pointSize=2, plot_PCs = "maxd2",
     plot_groupColors=c("orangered","blueviolet"),
     plot_dBoundary.level=c(0.025,0.5,0.975),
     plot_dBoundary.color=c("blueviolet","black","orangered"),
     plot_dBoundary.size=c(1,1,1),
     plot_dBoundary.alpha=c(1,1,1),

     option3D_engine="rgl", option3D_material="default", option3D_showTitles=TRUE,
     option3D_spin.axisXYZ=c(0, 0, 1), option3D_spin.duration=30,
     option3D_spin.fps=20, option3D_spin.startTime=0,
     option3D_outputFileName="PUMA3D",

     option2D_bgProbAlpha = TRUE,
     option2D_err.mark = "cv", option2D_err.color = "black",
     option2D_err.pointSize=4, option2D_err.pointShape=4
     )


}

\Datasets{
PUMA package contains 4 normalized data sets:
\item{arcene}{
This is a normalized version of the publically available Arcene mass-spectrometry data set which contains 9961 variables and compares the mass spectrometry data of cancer (1, n=88) versus normal subjects (2, n=112). Noteworthy, the Arcene dataset was one of the datasets used in Neural Information Processing Systems (NIPS) 2003 Feature Selection Challenge and contains 3000 'probes' - added distractor features (of random orders and patterns) that are having no predictive power (pure noise).

For more details about Arcene data set see \url{ https://archive.ics.uci.edu/ml/datasets/Arcene }
}

\item{coffee}{
This is a normalized verion of the Coffee data set which contains nuclear magnetic resonance (NMR) profiles of Colombian (n=15) vs non-Colombian (n=19) coffee bean samples using 1175 normalized spectral bins.

For more details about Coffee data set see \url{ https://doi.org/10.1155/2016/8564584 }
}

\item{diabetes}{
This is a normalized verion of The Qatar Metabolomics Study on Diabetes (QMDiab)publicly available data set containing untargeted plasma metabolomics (LC/MS+, LC/MS-, and GC/MS) from 179 diabetic patients and 177 non-diabetic normal controls.

For more details about QMDiab data set see \url{ https://doi.org/10.6084/m9.figshare.5904022.v5 }
}

\item{mice}{
This is a normalized version of the Mice data set which compares the difference in 28 proteins between control mice (c-SC-s: control mice, not stimulated to learn, injected with saline. n=135) and mice with trisomy 21 i.e. Down syndrome (t-CS-m: trisomy mice, stimulated to learn, injected with memantine. n=135).

For more details about Mice data set see \url{ https://doi.org/10.6084/m9.figshare.1421985.v1 }
}

}

\keyword{ PUMA }
\keyword{ PCA reversal }
\concept{ PCA Probabilistic mapping }


