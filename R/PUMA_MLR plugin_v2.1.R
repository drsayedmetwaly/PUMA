#-------------------------------------------------------------------------------------
#Probabilistic Universal Model Approximator (PUMA)
#
#A function for 3D visualization of decision boundary for classification algorithms


  PUMA = function( inputData, learner, task, measures=mmce,
                   override_defaulat_ploting_PCs=NULL,
                   cv = 10L,
                   pointsize = 2,
                   grp.cols = c("darkblue", "green", "darkred"),
                   LearnerName=toupper(getLearnerShortName(learner)),
                   contours3D, contours3Dcolors, engine3D="rgl",
                   material3D="default",
                   showTitles=TRUE,
                   spin_axisXYZ = c(0, 0, 1),
                   spin_duration = 30,
                   spin_fps=20,
                   spin_startTime = 0,
                   output3DFileName

                   ) {

    library(mlr)
    library(checkmate)
    library(Hmisc)
    library(sinkr)
    library(rgl)
    library(misc3d)

    learner = checkLearner(learner)
    assert(
      checkClass(task, "ClassifTask"),
    )
    td = getTaskDesc(task)
    inputData_ClassCol<-which( colnames(inputData)=="Class")

    features = getTaskFeatureNames(task)
    f.no = length(features)
    if (f.no==0) stopf("No features (metabolites) were input. Exiting...")
    if (f.no<3) stopf("Too few features/metabolites. PUMA algorithm is designed for 3D visualization of multivariate models (3 or more features/metabolites). Exiting...")

    taskdim=3

      PCA_Class<-inputData$Class
      PCA_DataSet<-inputData[,-c(inputData_ClassCol)]

      pcaTransform <- prcomp(PCA_DataSet, scale = TRUE, rank. =f.no)

      pcaScores <- as.data.frame(pcaTransform$x)
      All.PCs<-colnames(pcaScores)

      pcaScores$Class<-as.factor(PCA_Class)

      if (is.null(override_defaulat_ploting_PCs)){
        PCs <- c("PC1","PC2","PC3")
      }else{
        PCs <-override_defaulat_ploting_PCs
      }

      unusedPCs<-setdiff(All.PCs,PCs)

    cv = asCount(cv)

    assertNumber(pointsize, lower = 1)

    task = subsetTask(task, features = features)
    learner = setHyperPars(learner)

    target = td$target
    data = getTaskData(task)
    y = getTaskTargets(task)

    if (td$type == "classif" && hasLearnerProperties(learner, "prob"))
      learner = setPredictType(learner, "prob")
    mod = train(learner, task)
    pred.train = predict(mod, task)
    yhat = pred.train$data$response
    perf.train = performance(pred.train, task = task, measures = measures)
    if (cv > 0L) {
      cv = crossval(learner, task, iters = cv, measures = measures, show.info = FALSE)
      perf.cv = cv$aggr
      pred.cv = cv$pred
    } else {
      perf.cv = NA_real_
    }

      x1n = PCs[1L]
      x1 = pcaScores[, x1n]

      x2n = PCs[2L]
      x2 = pcaScores[, x2n]

      x3n = PCs[3L]
      x3 = pcaScores[, x3n]

      gridsize=50

      xs <- seq(min(x1), max(x1), length.out = gridsize)
      ys <- seq(min(x2), max(x2), length.out = gridsize)
      zs <- seq(min(x3), max(x3), length.out = gridsize)

      grid = expand.grid(xs,ys,zs)
      colnames(grid) = PCs

      grid2=grid
      for (i in 1:length(unusedPCs)){
        grid2$newCol<-0
        colnames(grid2)[ncol(grid2)]<-unusedPCs[i]
      }
      grid2<-as.matrix(grid2)
      reversedMat<-t(t(grid2 %*% t(pcaTransform$rotation)) * pcaTransform$scale+ pcaTransform$center)

      reversedDF<-as.data.frame(reversedMat)
      colnames(reversedDF)<-features

    pred.grid = predict(mod, newdata = reversedDF)

    grid[, target] = pred.grid$data$response

        if (hasLearnerProperties(learner, "prob")) {
          grid$value = pred.grid$data[,1]
        }else{
          grid$value = as.numeric(as.factor(pred.grid$data[,1]))-1
          contours3D=1
        }

          vran <- range(grid$value)
          levs <- seq(vran[1], vran[2], length.out=contours3D+2)[-c(1, contours3D+2)]

          if (missing(contours3Dcolors)){
            levcols <- val2col(levs, jetPal(contours3D), zlim = vran)
          }else{
            levcols <- contours3Dcolors[contours3D]
          }

          if (missing(grp.cols)){
            pointsColor <- val2col(as.numeric(as.factor(data$Class))-1,
                                   col=jetPal(length(unique(data$Class))),
                                   zlim=as.numeric(as.factor(data$Class))-1)
          }else{
            pointsColor <- grp.cols[1:length(unique(data$Class))]
          }

          fun <- function(x,y,z){return(grid$value)}

          xlab = colnames(grid)[1]
          ylab = colnames(grid)[2]
          zlab = colnames(grid)[3]

          title = sprintf("%s: %s", LearnerName, paramValueToString(learner$par.set, learner$par.vals))
          subtitle = sprintf("Train: %s; CV: %s", paste(names(perf.train),": ",perf.train, sep=""),
                          paste(names(perf.cv),": ",perf.cv, sep=""))


          while (rgl.cur() > 0) { rgl.close() }

          open3d()

          par3d(windowRect = c(20, 30, 800, 800))

          with(pcaScores, spheres3d(x=pcaScores[,xlab], y=pcaScores[,ylab], z=pcaScores[,zlab],
                               col=ifelse(pcaScores$Class==td$class.levels[1], pointsColor[1],
                                           ifelse(pcaScores$Class==td$class.levels[2], pointsColor[2],
                                                   pointsColor[3])),
                               radius = pointsize/10))

          contour3d(fun, level = levs,
                    x=xs, y=ys, z=zs,
                    color=levcols, material=material3D,
                    engine=engine3D, add=TRUE, alpha=0.5
          )

          axes3d(edges = "bbox", labels = TRUE, tick = TRUE,
                 box = TRUE, expand = 1.03)

          if (showTitles==TRUE){
            title3d(title, subtitle, xlab, ylab, zlab)
          } else {
            title3d("", "", xlab, ylab, zlab)
          }

          legend3d("topright", legend = paste(td$class.levels), pch = 19,
                   col = pointsColor, cex=1, inset=c(0.02))

          if (!missing(output3DFileName)){
            timeStamp <-  strftime(Sys.time(),"%Y-%m-%d_%H.%M.%S")
            FileNameSuffix = sprintf("_%s_%s", toupper(getLearnerShortName(learner)),timeStamp)


            if (spin_duration>0){
              play3d(spin3d(axis = spin_axisXYZ), duration = spin_duration,
                      fps=spin_fps,movie = paste(output3DFileName,FileNameSuffix,sep=""), frames = output3DFileName,
                      convert = TRUE, clean = TRUE, verbose = FALSE,top = TRUE,
                      type = "gif", startTime = spin_startTime, dir = getwd())
            }else{
              snapshot3d(paste(output3DFileName,FileNameSuffix,".png",sep=""))
            }
          }





  }

#End of PUMA function
#--------------------------------------------------------------------------------------

  installdep<-function(){
    install.packages(c("mlr","checkmate","Hmisc","rgl","misc3d"))
    library(devtools)
    devtools::install_github("marchtaylor/sinkr", force = TRUE)

  }

#--------------------------------------------------------------------------------------
  loadcsv<-function(CSV_File, idCol=1, classCol=2, excludeColNo=NULL)
  {
    Load_CSV <-read.csv(CSV_File, sep=",", check.names=TRUE, stringsAsFactors=FALSE, comment.char="")

    #rearrange data columns/exclude columns
    if (is.null(excludeColNo)){
      csv.data<-Load_CSV[,c(idCol,classCol, c(1:ncol(Load_CSV))[-c(idCol,classCol)])]
    } else {
      csv.data<-Load_CSV[,c(idCol,classCol, c(1:ncol(Load_CSV))[-c(idCol,classCol,excludeColNo)])]
    }

    #Sets the rownames of the data to sample names
    rownames(csv.data) <- csv.data[, 1]

    #Removes the 1st column (sample ID) and 2nd column (classes) from the dataset
    dat <- csv.data[,c(3:ncol(csv.data))]
    #Creates a new column containing Class data at the end of the dataset
    dat$Class<-factor(csv.data[, 2])

    return(dat)
  }

#--------------------------------------------------------------------------------------
  snapshot<-function(output3DFileName){
    timeStamp <-  strftime(Sys.time(),"%Y-%m-%d_%H.%M.%S")
    FileNameSuffix = sprintf("_%s",timeStamp)

    snapshot3d(paste(output3DFileName,FileNameSuffix,".png",sep=""))
  }

#--------------------------------------------------------------------------------------
  animation<-function(spin_axisXYZ=c(0, 0, 1), spin_duration=30, spin_fps=20, spin_startTime=0, output3DFileName=NULL){
    timeStamp <-  strftime(Sys.time(),"%Y-%m-%d_%H.%M.%S")
    FileNameSuffix = sprintf("_%s",timeStamp)

    if (!is.null(output3DFileName)){
      movie3d(spin3d(axis = spin_axisXYZ), duration = spin_duration,
              fps=spin_fps,movie = paste(output3DFileName,FileNameSuffix,sep=""),
              frames = output3DFileName,
              convert = TRUE, clean = TRUE, verbose = FALSE,top = TRUE,
              type = "gif", startTime = spin_startTime, dir = getwd())
    } else {
      play3d(spin3d(axis = spin_axisXYZ), duration = spin_duration,
             startTime = spin_startTime)

    }
  }
