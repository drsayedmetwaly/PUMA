#-------------------------------------------------------------------------------------
#Probabilistic Universal Model Approximator (PUMA)
#
#A function for 2D/3D visualization of decision boundary for classification algorithms


  PUMA = function( inputData, targetColumn="Class",

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
                   option2D_err.pointShape = 4


                   ) {

    library(mlr)
    library(checkmate)
    library(Hmisc)
    library(BBmisc)
    library(rgl)
    library(misc3d)


    ## Create an mlr task
    task <- makeClassifTask(data = inputData, target = targetColumn)#"Class")


    model = checkLearner(model)
    assert(
      checkClass(task, "ClassifTask"),
    )
    td = getTaskDesc(task)
    inputData_ClassCol<-which( colnames(inputData)==targetColumn)#"Class")

    features = getTaskFeatureNames(task)
    f.no = length(features)
    if (f.no==0) stopf("No features/variables were input. Exiting...")
    if (f.no<3) stopf("Too few features/variables. PUMA algorithm is designed for 2D/3D visualization of multivariate models (3 or more features/variables). Exiting...")

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

      plot_PCs<-toupper(plot_PCs)
      if (plot_PCs=="DEFAULT"){
        PCs <- c("PC1","PC2","PC3")[1:taskdim]
      }else{
        if (plot_PCs=="MAXD2"){

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
          PCs <- plot_PCs
        }

      }


        if ((length(unique(PCs %nin% All.PCs))==1 && unique(PCs %nin% All.PCs)==TRUE) || length(unique(PCs %nin% All.PCs))>1
            #PCs %nin% All.PCs
            ) {
          stopf("Principal component [%s] does not exist.",
                  paste(PCs[PCs %nin% All.PCs], collapse=", ") )
              }

      unusedPCs<-setdiff(All.PCs,PCs)

    model_crossValid = asCount(model_crossValid)

    assertNumber(plot_pointSize, lower = 1)
    assertFlag(option2D_bgProbAlpha)
    assertChoice(option2D_err.mark, choices = c("train", "cv", "none"))
    assertString(option2D_err.color)
    assertNumber(option2D_err.pointSize, lower = 1)


    task = subsetTask(task, features = features)
    model = setHyperPars(model)

    target = td$target
    data = getTaskData(task)
    y = getTaskTargets(task)

    if (td$type == "classif" && hasLearnerProperties(model, "prob"))
      model = setPredictType(model, "prob")
    mod = train(model, task)
    pred.train = predict(mod, task)
    yhat = pred.train$data$response
    perf.train = performance(pred.train, task = task, measures = model_perfMetric)
    if (model_crossValid > 0L) {
      cv = crossval(model, task, iters = model_crossValid, measures = model_perfMetric, show.info = FALSE)
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

      # if (hasLearnerProperties(model, "prob")) {
      #   grid$value = pred.grid$data[,1]
      # }else{
      #   grid$value = as.numeric(as.factor(pred.grid$data[,1]))-1
      # }



      contours2D = length(plot_dBoundary.level)
      if (sum(plot_dBoundary.level)==0){contours2D=1}

      if (hasLearnerProperties(model, "prob")) {
        grid$value = pred.grid$data[,1]
      }else{
        grid$value = as.numeric(as.factor(pred.grid$data[,1]))-1
        contours2D=1
      }

      vran <- range(grid$value)

      if (mean(plot_dBoundary.level, na.rm=TRUE)>1){#number of contours entered
        levs <- seq(vran[1], vran[2], length.out=contours2D+2)[-c(1, contours2D+2)]
      } else {# probability levels entered
        levs <- plot_dBoundary.level[1:contours2D]
      }





      if (missing(plot_groupColors)){
        pointsColor <- val2col(as.numeric(as.factor(data$Class))-1,
                               col=jetPal(length(unique(data$Class))),
                               zlim=as.numeric(as.factor(data$Class))-1)
      }else{
        pointsColor <- plot_groupColors[1:length(unique(data$Class))]
      }

      backgColor <- sapply(pointsColor, lighten, USE.NAMES=FALSE)


      # Define error points
      # pcaScores$.err = if (option2D_err.mark == "train")
      #     y != yhat
      #   else if (option2D_err.mark == "cv")
      #     y != pred.cv$data[order(pred.cv$data$id), "response"]
      if (option2D_err.mark == "train"){ pcaScores$.err <- (y != yhat)}
      if (option2D_err.mark == "cv") {
        pcaScores$.err <- (y != pred.cv$data[order(pred.cv$data$id), "response"])
      }

          # Initialize plot object p
          p = ggplot(grid, aes_string(x = x1n, y = x2n))

          # Plot background +/- transparency
          if (hasLearnerProperties(model, "prob") && option2D_bgProbAlpha) {
            prob = apply(getPredictionProbabilities(pred.grid, cl = td$class.levels), 1, max)
            grid$.prob.pred.class = prob
            p = p + geom_tile(data = grid, mapping = aes_string(fill = target, alpha = ".prob.pred.class"),
                              show.legend = TRUE)
            p = p + scale_alpha(limits = range(grid$.prob.pred.class))
          } else {
            backgColor <- sapply(backgColor, lighten, USE.NAMES=FALSE)
            p = p + geom_tile(mapping = aes_string(fill = target))
          }

          # Plot decision boundary
          if (sum(plot_dBoundary.level)!=0){
            for (i in 1:contours2D){
              grid$.contour<-(grid$value<levs[i]+0.005 & grid$value>levs[i]-0.005) #levs[i])

              tryCatch(
                {
                  p = p + geom_smooth(data = subset(grid, grid$.contour),
                                      mapping = aes_string(x = x1n, y = x2n),
                                      #method = "gam",
                                      size = plot_dBoundary.size[i],
                                      color = plot_dBoundary.color[i],
                                      alpha = plot_dBoundary.alpha[i],
                                      se = FALSE)
                },
                error = function(e) {}

              )

            }

          }

          # # Plot correctly predicted points
          # p = p + geom_point(data = subset(pcaScores, !pcaScores$.err),
          #                    mapping = aes_string(x = x1n, y = x2n, color  = target),
          #                    size = plot_pointSize)
          #
          # # Plot error points
          # p = p + geom_point(data = subset(pcaScores, pcaScores$.err),
          #                    mapping = aes_string(x = x1n, y = x2n, color  = target ),
          #                    size = option2D_err.pointSize, show.legend = FALSE)

          # Plot all points
          p = p + geom_point(data = pcaScores,
                             mapping = aes_string(x = x1n, y = x2n, color  = target),
                             size = plot_pointSize)
          p = p + scale_color_manual(values=pointsColor)

          # Mark error points
          if (option2D_err.mark != "none" && any(pcaScores$.err)) {
            # p = p + geom_point(data = subset(pcaScores, pcaScores$.err),
            #                    mapping = aes_string(x = x1n, y = x2n),
            #                    size = option2D_err.pointSize + 1.5, show.legend = FALSE)
            # p = p + geom_point(data = subset(pcaScores, pcaScores$.err),
            #                    mapping = aes_string(x = x1n, y = x2n),
            #                    size = option2D_err.pointSize + 1, col = option2D_err.color, show.legend = FALSE)
            p = p + geom_point(data = subset(pcaScores, pcaScores$.err),
                               mapping = aes_string(x = x1n, y = x2n),
                               size = option2D_err.pointSize, shape=option2D_err.pointShape, col = option2D_err.color)
          }


          p  = p + guides(alpha = FALSE)


          p = p + scale_fill_manual(values = backgColor)


          title = sprintf("%s: %s", model_name, paramValueToString(model$par.set, model$par.vals))
          title = sprintf("%s\nTrain: %s; model_crossValid: %s", title, paste(names(perf.train),": ",perf.train, sep=""),
                          paste(names(perf.cv),": ",perf.cv, sep=""))
          p = p + ggtitle(title)


          #return(p)
          suppressMessages(print(p))

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

    contours3D = length(plot_dBoundary.level)


        if (hasLearnerProperties(model, "prob")) {
          grid$value = pred.grid$data[,1]
        }else{
          grid$value = as.numeric(as.factor(pred.grid$data[,1]))-1
          contours3D=1
        }

    vran <- range(grid$value)

    if (mean(plot_dBoundary.level, na.rm=TRUE)>1){#number of contours entered
         levs <- seq(vran[1], vran[2], length.out=contours3D+2)[-c(1, contours3D+2)]
    } else {# probability levels entered
         levs <- plot_dBoundary.level[1:contours3D]
    }

          if (missing(plot_dBoundary.color)){
            levcols <- val2col(levs, jetPal(contours3D), zlim = vran)
          }else{
            levcols <- plot_dBoundary.color[1:contours3D]
          }

          if (missing(plot_groupColors)){
            pointsColor <- val2col(as.numeric(as.factor(data$Class))-1,
                                   col=jetPal(length(unique(data$Class))),
                                   zlim=as.numeric(as.factor(data$Class))-1)
          }else{
            pointsColor <- plot_groupColors[1:length(unique(data$Class))]
          }

          fun <- function(x,y,z){return(grid$value)}

          xlab = colnames(grid)[1]
          ylab = colnames(grid)[2]
          zlab = colnames(grid)[3]

          title = sprintf("%s: %s", model_name, paramValueToString(model$par.set, model$par.vals))
          subtitle = sprintf("Train: %s; model_crossValid: %s", paste(names(perf.train),": ",perf.train, sep=""),
                          paste(names(perf.cv),": ",perf.cv, sep=""))


          while (rgl.cur() > 0) { rgl.close() }

          open3d()

          par3d(windowRect = c(20, 30, 800, 800))

          with(pcaScores, spheres3d(x=pcaScores[,xlab], y=pcaScores[,ylab], z=pcaScores[,zlab],
                               col=ifelse(pcaScores$Class==td$class.levels[1], pointsColor[1],
                                           ifelse(pcaScores$Class==td$class.levels[2], pointsColor[2],
                                                   pointsColor[3])),
                               radius = plot_pointSize/10))

          # Plot decision boundary
          if (sum(plot_dBoundary.level)!=0){
            for (i in 1:contours3D){
              tryCatch(
                {
                  contour3d(fun, level = levs[i],
                            x=xs, y=ys, z=zs,
                            color=levcols[i], material=option3D_material,
                            engine=option3D_engine, add=TRUE, alpha=plot_dBoundary.alpha[i]
                  )
                },
                error = function(e) { message(paste("Could not plot decision boundary at level [",levs[i],"]. ", e))
                  }
              )

            }

          }




          # contour3d(fun, level = levs,
          #           x=xs, y=ys, z=zs,
          #           color=levcols, material=option3D_material,
          #           engine=option3D_engine, add=TRUE, alpha=plot_dBoundary.alpha
          # )

          axes3d(edges = "bbox", labels = TRUE, tick = TRUE,
                 box = TRUE, expand = 1.03)

          if (option3D_showTitles==TRUE){
            title3d(title, subtitle, xlab, ylab, zlab)
          } else {
            title3d("", "", xlab, ylab, zlab)
          }

          legend3d("topright", legend = paste(td$class.levels), pch = 19,
                   col = pointsColor, cex=1, inset=c(0.02))


            if (option3D_spin.duration==0){
              snapshot(option3D_outputFileName=option3D_outputFileName)

            }else{
              animation(option3D_spin.axisXYZ=option3D_spin.axisXYZ,
                        option3D_spin.duration=option3D_spin.duration,
                        option3D_spin.fps=option3D_spin.fps,
                        option3D_spin.startTime=option3D_spin.startTime,
                        option3D_outputFileName=option3D_outputFileName)
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

  }

#-------------------------------------------------------------------------------------
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
  snapshot<-function(option3D_outputFileName){
    timeStamp <-  strftime(Sys.time(),"%Y-%m-%d_%H.%M.%S")
    FileNameSuffix = sprintf("_%s",timeStamp)

    if(!is.null(option3D_outputFileName)){
    snapshot3d(paste(option3D_outputFileName,FileNameSuffix,".png",sep=""))
    }
  }

#--------------------------------------------------------------------------------------
  animation<-function(option3D_spin.axisXYZ=c(0, 0, 1), option3D_spin.duration=30, option3D_spin.fps=20, option3D_spin.startTime=0, option3D_outputFileName=NULL){
    timeStamp <-  strftime(Sys.time(),"%Y-%m-%d_%H.%M.%S")
    FileNameSuffix = sprintf("_%s",timeStamp)

    if (!is.null(option3D_outputFileName)){
      movie3d(spin3d(axis = option3D_spin.axisXYZ), duration = option3D_spin.duration,
              fps=option3D_spin.fps,movie = paste(option3D_outputFileName,FileNameSuffix,sep=""),
              frames = option3D_outputFileName,
              convert = TRUE, clean = TRUE, verbose = FALSE,top = TRUE,
              type = "gif", startTime = option3D_spin.startTime, dir = getwd())
    } else {
      play3d(spin3d(axis = option3D_spin.axisXYZ), duration = option3D_spin.duration,
             startTime = option3D_spin.startTime)

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

#--------------------------------------------------------------------------------------
# val2col: Convert values to color levels
# In marchtaylor/sinkr: Collection of functions with emphasis in multivariate data analysis

val2col <-  function (z, zlim, col = heat.colors(12), breaks)
  {
    if (!missing(breaks)) {
      if (length(breaks) != (length(col) + 1)) {
        stop("must have one more break than color")
      }
    }
    if (missing(breaks) & !missing(zlim)) {
      breaks <- seq(zlim[1], zlim[2], length.out = (length(col) +
                                                      1))
    }
    if (missing(breaks) & missing(zlim)) {
      zlim <- range(z, na.rm = TRUE)
      breaks <- seq(zlim[1], zlim[2], length.out = (length(col) +
                                                      1))
    }
    CUT <- cut(z, breaks = breaks, include.lowest = TRUE)
    colorlevels <- col[match(CUT, levels(CUT))]
    return(colorlevels)
  }

