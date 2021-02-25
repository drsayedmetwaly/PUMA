#-------------------------------------------------------------------------------------
#Probabilistic Universal Model Approximator (PUMA)
#
#A function for 2D/3D visualization of decision boundary for classification algorithms


  PUMA = function( inputData, learner, task, measures=mmce,
                   override_default_ploting_PCs=NULL,
                   cv = 7L,
                   plot_dimensions=3,
                   pointsize = 2,
                   grp.cols = c("darkblue", "green", "darkred"),
                   LearnerName=toupper(getLearnerShortName(learner)),
                   contours3D, contours3Dcolors, engine3D="rgl",
                   material3D="default",
                   showTitles3D=TRUE,
                   spin3D_axisXYZ = c(0, 0, 1),
                   spin3D_duration = 30,
                   spin3D_fps=20,
                   spin3D_startTime = 0,
                   output3DFileName = NULL,
                   prob.alpha2D = TRUE,
                   err.mark2D = "cv",
                   err.col2D = "black", err.size2D = 4, err.shape2D=4


                   ) {

    library(mlr)
    library(checkmate)
    library(Hmisc)
    library(BBmisc)
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
    if (f.no<3) stopf("Too few features/metabolites. PUMA algorithm is designed for 2D/3D visualization of multivariate models (3 or more features/metabolites). Exiting...")

    assertNumber(plot_dimensions,
                 na.ok = FALSE,
                 lower = 2,
                 upper = 3,
                 finite = FALSE,
                 null.ok = FALSE,
                 add = NULL)

    taskdim=plot_dimensions

      PCA_Class<-inputData$Class
      PCA_DataSet<-inputData[,-c(inputData_ClassCol)]

      pcaTransform <- prcomp(PCA_DataSet, scale = TRUE, rank. =f.no)

      pcaScores <- as.data.frame(pcaTransform$x)
      All.PCs<-colnames(pcaScores)

      pcaScores$Class<-as.factor(PCA_Class)

      if (is.null(override_default_ploting_PCs)){
        PCs <- c("PC1","PC2","PC3")[1:taskdim]
      }else{
        override_default_ploting_PCs<-toupper(override_default_ploting_PCs)
        if (override_default_ploting_PCs=="MAXD2"){

          S1 = pcaScores[which(pcaScores$Class==td$class.levels[1]),]
          S2 = pcaScores[which(pcaScores$Class!=td$class.levels[1]),]

          library(arrangements)
          if (length(All.PCs) > 1000){combNo=1000}else{combNo=length(All.PCs)}
          perMat<-combinations(combNo, taskdim)
          Dist<-vector(length = nrow(perMat), mode = "numeric")
          for (i in 1:nrow(perMat)){
            tryCatch(
              {
                Dist[i]<-D.sq(S1[,perMat[i,]],S2[,perMat[i,]])[['D.sq']]
              },
              error = function(e) {}

            )

          }
          perMat<-cbind(perMat, Dist)
          perMat<-perMat[order(perMat[,"Dist"],decreasing=TRUE),]
          PCs <-colnames(pcaScores[,perMat[1,c(1:taskdim)]])

        }else {
          PCs <- override_default_ploting_PCs
        }

      }


        if ((length(unique(PCs %nin% All.PCs))==1 && unique(PCs %nin% All.PCs)==TRUE) || length(unique(PCs %nin% All.PCs))>1
            #PCs %nin% All.PCs
            ) {
          stopf("Principal component [%s] does not exist.",
                  paste(PCs[PCs %nin% All.PCs], collapse=", ") )
              }

      unusedPCs<-setdiff(All.PCs,PCs)

    cv = asCount(cv)

    assertNumber(pointsize, lower = 1)
    assertFlag(prob.alpha2D)
    assertChoice(err.mark2D, choices = c("train", "cv", "none"))
    assertString(err.col2D)
    assertNumber(err.size2D, lower = 1)


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


    if (taskdim==2L){

      x1n = PCs[1L]
      x1 = pcaScores[, x1n]

      x2n = PCs[2L]
      x2 = pcaScores[, x2n]

      gridsize=100

      xs <- seq(min(x1), max(x1), length.out = gridsize)
      ys <- seq(min(x2), max(x2), length.out = gridsize)

      grid = expand.grid(xs,ys)
      colnames(grid) = PCs

      grid2=grid
      for (i in 1:length(unusedPCs)){
        grid2$newCol<-0
        colnames(grid2)[ncol(grid2)]<-unusedPCs[i]
      }
      grid2<-as.matrix(grid2[,All.PCs])
      reversedMat<-t(t(grid2 %*% t(pcaTransform$rotation)) * pcaTransform$scale+ pcaTransform$center)

      reversedDF<-as.data.frame(reversedMat)
      colnames(reversedDF)<-features

      pred.grid = predict(mod, newdata = reversedDF)

      grid[, target] = pred.grid$data$response

      if (hasLearnerProperties(learner, "prob")) {
        grid$value = pred.grid$data[,1]
      }else{
        grid$value = as.numeric(as.factor(pred.grid$data[,1]))-1
      }


      if (missing(grp.cols)){
        pointsColor <- val2col(as.numeric(as.factor(data$Class))-1,
                               col=jetPal(length(unique(data$Class))),
                               zlim=as.numeric(as.factor(data$Class))-1)
      }else{
        pointsColor <- grp.cols[1:length(unique(data$Class))]
      }

      backgColor <- sapply(pointsColor, lighten, USE.NAMES=FALSE)


      # Define error points
      pcaScores$.err = if (err.mark2D == "train")
          y != yhat
        else if (err.mark2D == "cv")
          y != pred.cv$data[order(pred.cv$data$id), "response"]


          # Initialize plot object p
          p = ggplot(grid, aes_string(x = x1n, y = x2n))

          # Define background transparency
          if (hasLearnerProperties(learner, "prob") && prob.alpha2D) {
            prob = apply(getPredictionProbabilities(pred.grid, cl = td$class.levels), 1, max)
            grid$.prob.pred.class = prob
            p = p + geom_tile(data = grid, mapping = aes_string(fill = target, alpha = ".prob.pred.class"),
                              show.legend = TRUE)
            p = p + scale_alpha(limits = range(grid$.prob.pred.class))
          } else {
            p = p + geom_tile(mapping = aes_string(fill = target))
          }

          # # Plot correctly predicted points
          # p = p + geom_point(data = subset(pcaScores, !pcaScores$.err),
          #                    mapping = aes_string(x = x1n, y = x2n, color  = target),
          #                    size = pointsize)
          #
          # # Plot error points
          # p = p + geom_point(data = subset(pcaScores, pcaScores$.err),
          #                    mapping = aes_string(x = x1n, y = x2n, color  = target ),
          #                    size = err.size2D, show.legend = FALSE)

          # Plot all points
          p = p + geom_point(data = pcaScores,
                             mapping = aes_string(x = x1n, y = x2n, color  = target),
                             size = pointsize)
          p = p + scale_color_manual(values=pointsColor)


          # Mark error points
          if (err.mark2D != "none" && any(pcaScores$.err)) {
            # p = p + geom_point(data = subset(pcaScores, pcaScores$.err),
            #                    mapping = aes_string(x = x1n, y = x2n),
            #                    size = err.size2D + 1.5, show.legend = FALSE)
            # p = p + geom_point(data = subset(pcaScores, pcaScores$.err),
            #                    mapping = aes_string(x = x1n, y = x2n),
            #                    size = err.size2D + 1, col = err.col2D, show.legend = FALSE)
            p = p + geom_point(data = subset(pcaScores, pcaScores$.err),
                               mapping = aes_string(x = x1n, y = x2n),
                               size = err.size2D, shape=err.shape2D, col = err.col2D)
          }


          p  = p + guides(alpha = FALSE)


          p = p + scale_fill_manual(values = backgColor)


          title = sprintf("%s: %s", LearnerName, paramValueToString(learner$par.set, learner$par.vals))
          title = sprintf("%s\nTrain: %s; CV: %s", title, paste(names(perf.train),": ",perf.train, sep=""),
                          paste(names(perf.cv),": ",perf.cv, sep=""))
          p = p + ggtitle(title)


          return(p)

      }





    if (taskdim==3L){
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
      grid2<-as.matrix(grid2[,All.PCs])
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

          if (showTitles3D==TRUE){
            title3d(title, subtitle, xlab, ylab, zlab)
          } else {
            title3d("", "", xlab, ylab, zlab)
          }

          legend3d("topright", legend = paste(td$class.levels), pch = 19,
                   col = pointsColor, cex=1, inset=c(0.02))


            if (spin3D_duration==0){
              snapshot(output3DFileName=output3DFileName)

            }else{
              animation(spin3D_axisXYZ=spin3D_axisXYZ,
                        spin3D_duration=spin3D_duration,
                        spin3D_fps=spin3D_fps,
                        spin3D_startTime=spin3D_startTime,
                        output3DFileName=output3DFileName)
            }


          }

  }

#End of PUMA function
#--------------------------------------------------------------------------------------

  installdep<-function(){
    install.packages(c("mlr","checkmate","Hmisc","rgl","misc3d", "caret", "plsdepot",
                       "pls", "randomForest", "mda", "bst", "adabag", "glmnet",
                       "deepnet", "SwarmSVM", "gbm", "MASS", "nnet", "deepnet",
                       "xgboost", "kernlab", "survival","RWeka", "rpart","fpc",
                       "FNN", "h2o", "gbm", "klaR", "modeltools", "ICSNP", "BBmisc"))
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

    if(!is.null(output3DFileName)){
    snapshot3d(paste(output3DFileName,FileNameSuffix,".png",sep=""))
    }
  }

#--------------------------------------------------------------------------------------
  animation<-function(spin3D_axisXYZ=c(0, 0, 1), spin3D_duration=30, spin3D_fps=20, spin3D_startTime=0, output3DFileName=NULL){
    timeStamp <-  strftime(Sys.time(),"%Y-%m-%d_%H.%M.%S")
    FileNameSuffix = sprintf("_%s",timeStamp)

    if (!is.null(output3DFileName)){
      movie3d(spin3d(axis = spin3D_axisXYZ), duration = spin3D_duration,
              fps=spin3D_fps,movie = paste(output3DFileName,FileNameSuffix,sep=""),
              frames = output3DFileName,
              convert = TRUE, clean = TRUE, verbose = FALSE,top = TRUE,
              type = "gif", startTime = spin3D_startTime, dir = getwd())
    } else {
      play3d(spin3d(axis = spin3D_axisXYZ), duration = spin3D_duration,
             startTime = spin3D_startTime)

    }
  }

#--------------------------------------------------------------------------------------
# Mahalanobis Distance calculation Function from https://stackoverflow.com/a/34708113/5731401
# Custom function to calculate Mahalanobis distance between 2 groups centroids
  D.sq <- function (g1, g2) {
    dbar <- as.vector(colMeans(g1) - colMeans(g2))
    S1 <- cov(g1)
    S2 <- cov(g2)
    n1 <- nrow(g1)
    n2 <- nrow(g2)
    V <- as.matrix((1/(n1 + n2 - 2)) * (((n1 - 1) * S1) + ((n2 - 1) * S2)))
    D.sq <- t(dbar) %*% solve(V) %*% dbar
    res <- list()
    res$D.sq <- D.sq
    res$V <- V
    res
  }

#--------------------------------------------------------------------------------------
  darken <- function(color, factor=1.4){
    if (factor < 1) stop("factor needs to be > 1.0")
    col <- col2rgb(color)
    col <- col/factor
    col <- rgb(t(col), maxColorValue=255)
    col
  }

#--------------------------------------------------------------------------------------
  lighten <- function(color, factor = 0.5) {
    if ((factor > 1) | (factor < 0)) stop("factor needs to be within [0,1]")
    col <- col2rgb(color)
    col <- col + (255 - col)*factor
    col <- rgb(t(col), maxColorValue=255)
    col
  }



