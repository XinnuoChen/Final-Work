############FUNCIONS###########

# Real data generation
twovar.data.generation <- function(N, mean, sd, varname=c("height","weight", "group"),groupseed=123){
  # mean1 dataframe de mitjana para cada grup i cada variable. columna-variable, fila-grup
  # groupnumber es nombre de cluster dins d'una base de dades
  groupnumber <- nrow(mean)
  n <- N/groupnumber # nombre d'obervacions per cluster
  
  # Comprovació
  if(nrow(mean) != nrow(sd)) {
    stop("The number of parameters does not match the number of groups.")
  }
  # comprovar que el nombre de parametres es correcte d'acord als grups que es vol crear
  
  # Generació de les dades
  set.seed(groupseed)
  var1 <- c()
  var2 <- c()
  for(i in 1:groupnumber){ #per cada base grup
    var1.group <- rnorm(n, mean = mean[i,1], sd = sd[i,1]) #fila grup, col variable
    var2.group <- rnorm(n, mean = mean[i,2], sd = sd[i,2])
    
    var1 <- c(var1, var1.group)
    var2 <- c(var2, var2.group)
  }
  
  group <- rep(1:groupnumber, each=n)
  data <- data.frame(var1, var2, group)
  colnames(data) <- varname
  return(data)
}


normal.data.generation <- function(n, mu, sigma = NULL, seed = NULL)
{
  if(!is.null(seed)) set.seed(seed)
  
  # Si no passes cap matriu, fem una de correlacions AR(1) (rho = 0.4)
  if(is.null(sigma)){
    rho   <- 0.4
    idx   <- 0:7
    sigma <- outer(idx, idx, function(i, j) rho^abs(i - j))  # AR(1)
  }
  
  stopifnot(length(mu) == 8,
            all(dim(sigma) == c(8, 8)),
            isSymmetric(sigma),
            all(eigen(sigma, only.values = TRUE)$values > 0))  # ha de ser definida positiva
  
  mat <- rmvnorm(n, mean = mu, sigma = sigma)
  as.data.frame(mat)
}


# Synthetic data generation with cp & minbucket
synt.data.generation <- function(real.data, cp.values, minbucket.values, m){
  
  specks.value <- data.frame(cp = numeric(0), minbucket = numeric(0), SPECKS = numeric(0))
  # sds.default <- list() # list on guardo els objectes syn
  synt.data <- list() # list on guardo totes les dades sintetiques generades
  
  for(i in 1:length(cp.values)){
    for(j in 1:length(minbucket.values)){
      syn.obj <- syn(real.data,  # Generació amb CART
                     method = rep("cart", ncol(real.data)),
                     cart.control = list(cp = cp.values[i], minbucket = minbucket.values[j]),
                     m = m,
                     print.flag = FALSE)
      # syn.obj <- as.list(syn.obj)
      # sds.default <- c(sds.default, syn.obj)
      synt.data <- c(synt.data, syn.obj$syn)
      
      for(k in 1:m){
        names(synt.data[[k]]) <- names(real.data)
        
        data <- rbind(
          cbind(real.data, label = 0),
          cbind(synt.data[[k]], label = 1)
        )
        
        # Calcular SPECKS para la combinación actual de cp y minbucket
        specks.val <- specks(data = data)$Specks
        specks.value <- rbind(specks.value, data.frame(minbucket = minbucket.values[j], cp = cp.values[i], SPECKS = specks.val))
      }
    }
  }  
  results <- list("Specks" = specks.value, "Syn data" = synt.data)
  return(results)
}


# K decision & clustering 
k.decision <- function(data, k.real = 0, synthetic = TRUE){
  if(synthetic == FALSE){
    kopt.obj <- NbClust(data = data, diss = NULL, distance = "euclidean", method = "kmeans", index = "all", alphaBeale = 0.1)
    kopt <- as.numeric(names(table(as.vector(kopt.obj$Best.nc[1,])))[which.max(table(as.vector(kopt.obj$Best.nc[1,])))])
    sdata <- scale(data)
    clust <- kmeans(sdata, kopt, nstart = 25)
    
    clustering.result <- list("K optim" = kopt, "Clustering" = clust)
    return(clustering.result)
  } else {
    kopt.obj <- NbClust(data = data, diss = NULL, distance = "euclidean", method = "kmeans", index = "all", alphaBeale = 0.1)
    kopt <- as.numeric(names(table(as.vector(kopt.obj$Best.nc[1,])))[which.max(table(as.vector(kopt.obj$Best.nc[1,])))])
    
    if(kopt == k.real){
      sdata <- scale(data)
      clust <- kmeans(sdata, kopt, nstart = 25) #resultado de clustering
      
      clustering.result <- list("K optim" = kopt, "Clustering" = clust)
      return(clustering.result)
    }
  }
}

# SPECKS CALCULATION
specks <- function(data, sample=FALSE, minbucket=5, cp = 1e-3){   #### With pvalue
  
  if(sample==TRUE){
    data[,3] <- sample(data[,3])
  }
  
  cart_model <- rpart(
    formula = label ~ .,             # label com a resposta, la resta com a predictors
    data = data,
    method = "class",                # arbre de classificació
    control = rpart.control(mincriterion = 0, 
                            minbucket = minbucket, 
                            cp = cp))
  
  pred.prob <- predict(cart_model, type = "prob")[, 2]
  
  label.vector <- data$label
  
  prop.score.real  <- pred.prob[label.vector == 0]
  prop.score.synt  <- pred.prob[label.vector == 1]
  
  ecdf.real  <- ecdf(prop.score.real)
  ecdf.synth <- ecdf(prop.score.synt)
  
  all.scores <- sort(unique(c(prop.score.real, prop.score.synt)))
  
  specks.value <- max(abs(ecdf.real(all.scores) - ecdf.synth(all.scores)))
  p.value <- suppressWarnings(ks.test(prop.score.real, prop.score.synt)$p)
  
  return(data.frame(
    "Specks" = specks.value, 
    "P-value" = p.value))
}

# Label switching function
label.switching <- function(km.real, km.synt){
  if(!is.null(km.synt)){
    centroids <- rbind(km.real$centers, km.synt$centers)
    dist.global <- as.matrix(dist(centroids, method = "euclidean"))
    n.real <- nrow(km.real$centers) #numero de cluster en real
    n.synt <- nrow(km.synt$centers) #numero de cluster en synt
    dist <- dist.global[1:n.real, (n.real + 1) : (n.real + n.synt)]
    label.ass <- solve_LSAP(dist, maximum = FALSE) 
    
    # Reetiquetació els clústers sintètics
    km.synt$recluster <- km.real$cluster # nova variable
    
    for (i in 1:length(label.ass)) {
      km.synt$recluster[km.synt$cluster == label.ass[i]] <- i
    }
  }
  return(km.synt)
}


# Metrics function
cluster.metrics <- function(km.real, km.synt, k.real, k.synt, real.data, synt.data){
  # Nombre de cluster
  clust.num <- k.synt/k.real
  
  # Nombre d'observacions dins de clúster -> accuracy
  conf.mat <- confusionMatrix(data = as.factor(km.synt$recluster), reference = as.factor(km.real$cluster))
  accuracy <- conf.mat$overall["Accuracy"]
  
  # Distribució d'observacions en clúster -> Índex de Gini
  gini.real <- Gini(km.real$size)
  gini.synt <- Gini(km.synt$size)
  gini <- gini.synt/gini.real
  
  # Distancia d'observacions a cluster -> Coeficient de silhouette (entre -1 i 1)
  sil.real <- mean(silhouette(km.real$cluster, dist(real.data))[,3])
  sil.synt <- mean(silhouette(km.synt$recluster, dist(synt.data))[,3])
  sil <- sil.synt/sil.real
  
  # Distancia mitjana entre centroides
  lreal.data <- real.data
  lreal.data$label <- as.factor(km.real$cluster)
  lsynt.data <- synt.data
  lsynt.data$label <- as.factor(km.synt$recluster)
  
  real.centroid <- aggregate(cbind(height, weight) ~ label,
                             data = lreal.data,
                             FUN = base::mean)
  synt.centroid <- aggregate(cbind(height, weight) ~ label,
                             data = lsynt.data,
                             FUN = base::mean)
  real.centroid <- real.centroid[order(real.centroid$label), ]
  synt.centroid <- synt.centroid[order(synt.centroid$label), ]
  distances <- sqrt(
    (real.centroid$height - synt.centroid$height)^2 + 
      (real.centroid$weight - synt.centroid$weight)^2
  )
  mean.distance <- mean(distances)
  
  # Variancia mitjana entre centroides
  real.var <- aggregate(cbind(height, weight) ~ label,
                        data = lreal.data,
                        FUN = var)
  synt.var <- aggregate(cbind(height, weight) ~ label,
                        data = lsynt.data,
                        FUN = var)
  real.var <- real.var[order(real.var$label), ]
  synt.var <- synt.var[order(synt.var$label), ]
  real.var <- real.var[,-1]
  synt.var <- synt.var[,-1]
  
  dif.var <- abs(real.var - synt.var)
  mean.var <- mean(unlist(dif.var))
  
  result <- list("Clus.num" = clust.num, "Accuracy" = accuracy, "Gini" = gini, "Sil" = sil, "Mean distance" = mean.distance, "Mean variance" = mean.var)
  return(result)
}

# plot.metrics
plot.metrics <- function(df, specks.col = "SPECKS") {
  num_cols <- sapply(df, is.numeric)
  num_cols[specks.col] <- FALSE
  
  long.df <- pivot_longer(df,
                          cols      = names(num_cols)[num_cols],
                          names_to  = "Metric",
                          values_to = "Value")
  
  ggplot(long.df,
         aes(x = .data[[specks.col]], y = Value)) +   # tidy eval
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = FALSE, colour = "red", formula = y ~ x) +
    facet_wrap(~ Metric, scales = "free_y") +
    labs(x = specks.col, y = NULL,
         title = "Relació de cada mètrica amb SPECKS") +
    theme_minimal()
}
