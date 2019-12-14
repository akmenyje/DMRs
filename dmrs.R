### data: data frame of signal per position + column of coordinates (start)
### name1, name2: groups (column names of data)
### ff: vector of fetal fraction in samples
### bandwidth: smoothing window
### cutoff: percentage of highest CG values to include in candidate regions. CG values are defined as the difference of smoothed signal between two groups
### cg: minimum number of CGs per region
### permut_n: number of permutations
### design: model matrix

library(dplyr)
library(rlist)

perm.reg <- function(data, cutoff = 0.25, bandwidth, cg = 5, name1 = "XX", name2 = "XY", ff, permut_n = 5, design = design){
  
  x <- which(colnames(data) %in% c(name1, name2))
  pos <- data$start
  data <- data.frame(data[,x], pos)
  
  regions <- function(design = design){
    
    
    colnames(design) <- c("intercept", "coeff")
    colnames(data) <- c(design[, "coeff"], "start")
    
    
    #### smoothing and difference 
    
    group1 <- data[,colnames(data) == "0"]
    group2 <- data[,colnames(data) == "1"]
    
    group1 <- group1 %>% rowMeans
    group2 <- group2 %>% rowMeans
    
    smooth_1 <- unlist(ksmooth(pos, group1, x.points = pos, kernel = "normal", bandwidth = bandwidth)$y)
    smooth_2 <- unlist(ksmooth(pos, group2, x.points = pos, kernel = "normal", bandwidth = bandwidth)$y)
    
    diffs <- smooth_1 - smooth_2
    cutoff <- quantile(abs(diffs), cutoff)
    
    ##### regions
    
    chr <- rep(1, length(pos))
    
    cluster <- bumphunter::clusterMaker(chr = chr, pos = pos, maxGap = 100, assumeSorted = TRUE)
    index <- bumphunter::getSegments(f = cluster, x = diffs, assumeSorted = TRUE, cutoff = cutoff)
    
    index <- index[seq_len(2)]
    
    index1 <- index[[1]][which(lengths(index[[1]]) >= cg)]
    index2 <- index[[2]][which(lengths(index[[2]]) >= cg)]
    
    
    indexes <- c(index1, index2)
    
    message(length(index1)," ", "female "," ", " ", "candidates")
    message(length(index2)," ", "male "," ", " ", "candidates")
    
    #indexes <- sample(indexes, 100)
    
    ### GLS
    
    index <- seq_along(indexes)
    
    l <- ncol(data) - 1
    
    gls.cov <- function(L, pos, correlation = corAR1(form = ~1 |sample), design, data, f, l){
      
      Ll <- as.vector(sapply(pos[L], function(x) rep(x, l)))
      
      ff <- rep(f, length(L))
      
      library(dplyr)
      library(nlme)
      
      colnames(data) <- c(1:l, "start")
      signal <- data %>% filter(start %in% unique(Ll)) %>% subset(select = -(l+1)) %>% apply(1, as.vector) %>% as.vector
      
      
      dat <- data.frame(pos = Ll, gender = rep(design[,2], length(L)), sample = rep(1:l, length(L)), data = signal, ff = ff)
      model <- formula(data ~ gender + factor(pos) + ff)
      
      fit <- summary(gls(model, data = dat, correlation = correlation))
      
      stat <- fit$tTable[2, 3]
      value <- fit$tTable[2, 1]
      f <- c(stat, value)
    }
    
    re <- do.call("rbind", lapply(indexes, gls.cov, pos = pos, design = design, data = data, f = ff, l  =  l))
    
    colnames(re) <- c("stat", "value")
    re <- as.data.frame(re)
    indexRanges <- as.data.frame(IRanges::IRanges(unlist(lapply(indexes, min)), 
                                         unlist(lapply(indexes, max))))
    
    coord <- as.data.frame(IRanges::IRanges(unlist(lapply(index, function(x) pos[min(indexes[[x]])])),
                                   unlist(lapply(index, function(x) pos[max(indexes[[x]])])))) 
    
    re <- as.data.frame(cbind(re, CpG = indexRanges$width, start = coord$start, end = coord$end, width = coord$width))
  }
  
  ##### candidate regions
  re <- regions(design = design)
  
  ### permutations
  sampleSize <- min(length(which(design[,2] == 1)), length(which(design[,2] == 0)))
  
  perms <- combn(seq(1, ncol(data)-1), sampleSize)
  
  ### Leaving only those permutation in which: 
  ### 1. samples of different groups are mixed 
  ### 2. groups have imbalanced number of group1 and group2 samples  
  
  coef_uniq <- NULL
  coef_equal <- NULL
  for (i in 1: ncol(perms)){
    p <- perms[,i]
    coef_uniq[i] <- length(unique(design[p, 2]))
    coef_equal <- c(coef_equal, sum(design[p, 2]))
  }
  
  ### permutations with no mixing between groups
  rmv <- which(coef_uniq != 2)
  ### permutations unbalanced mixing between groups
  rmv2 <- which(coef_equal > round(sampleSize/2))
  
  perms <- perms[,-c(rmv, rmv2)]
  
  samp <- sample(1:ncol(perms), permut_n)
  perms <- perms[, samp]
  
  #### Performing permutations. Null distribution of statistics
  
  fstat <- NULL
  coord <- NULL
  
  for(i in 1:ncol(perms)){
    ### design matrix
    designr <- design
    reorder <- perms[,i]
    designr[,2] <- 0
    designr[reorder, 2] <- 1
    ### 
    fre <- regions(design = designr)
    fstat <- c(fstat, fre$stat)
    ### pool of regions' coordinates for the enrichment analysis
    coord <- rbind(coord, cbind(start = fre$start, end = fre$end))
  }
  
  QU <- quantile(fstat, 0.975)
  QD <- quantile(fstat, 0.025)
  
  ### candidates
  
  rmv <- which(re$stat <= QU & re$stat >= QD)
  res <- re[-rmv,]
  
  
  if (length(res$stat) == 0) {
    message("no candidates regions found")
    res <- NULL
  } else {
    rstat <- res$stat
    l <- length(rstat)
    p <- rep(NA, l)
    
    ### p value
    for (i in 1:l) {
      if(rstat[i] >= QU) {
        p[i] <- length(fstat[fstat >= rstat[i]])/length(fstat)
      } else {
        p[i] <- length(fstat[fstat <= rstat[i]])/ length(fstat)
      }
    }
    
    res$pval <- p
    ### adjusted p value
    res$qval <- p.adjust(res$pval, method = "BH")
  }
  res <- res[order(abs(res$stat), decreasing = TRUE),]
  res <- list(res, coord)
}

# glmnet
# training regions

train.dmrs <- function(data, ff, bandwidth){
  
  library(rlist)
  library(dplyr)
  
  dat <- data %>% subset(select = -(i+2))
  
  names <- unique(colnames(data)[-c(1,2)])
  name1 <- names[1]
  name2 <- names[2]
  names1 <- grep(name1, colnames(data))
  names2 <- grep(name2, colnames(data))
  
  nm <- NULL
  nm[names1] <- name1
  nm[names2] <- name2
  nm[c(1,2)] <- c("chr", "start")
  colnames(dat) <- nm[-(i+2)]
  
  design <- model.matrix(~colnames(dat)[3:ncol(dat)])
  
  
  ff <- ff[-i]
  
  re <- perm.reg(data = dat, cutoff = 0.25, cg = 5, bandwidth = bandwidth, name1 = name1, name2 = name2, ff = ff, design = design)[[1]]
  re
}

# glmnet regression and predictions 
# data: nelogoritmuotas signalas
# res: train.dmrs results (list)

predict.lasso <- function(alpha, res, data){
  
  names <- unique(colnames(data)[-c(1,2)])
  name1 <- names[1]
  name2 <- names[2]
  pos <- data[[2]]
  
  response <- NULL
  l <- 1:length(res)
  for(i in l){
    train_dat <- data %>% subset(select = -(i+2))
    names <- NULL
    names[grep(name1, colnames(train_dat))] <- name1
    names[grep(name2, colnames(train_dat))] <- name2
    names[1:2] <- c("chr", "start")
    colnames(train_dat) <- names
    
    train_reg <- res[[i]]
    
    st <- which(pos %in% train_reg$start)
    end <- which(pos %in% train_reg$end)
    
    x.train <- NULL
    for (x in 1:length(st)){
      index <- st[x]:end[x]
      dat <- train_dat[index, -c(1,2)] 
      x.train <- rbind(x.train, apply(dat, 2, sum))
    }
    
    x.train <- t(x.train)
    
    y.train <- colnames(train_dat)[-c(1,2)]
    y.train <- ifelse(y.train == name1, 1, 0)
    
    model <- cv.glmnet(x.train, y.train, alpha = alpha, family = "binomial")
    
    #### test
    test_dat <- data %>% subset(select = (i+2))
    
    x.test <- NULL
    for (y in 1:length(st)){
      index <- st[y]:end[y]
      x.test <- rbind(x.test, sum(dat))
    }
    
    x.test <- t(x.test) 
    response[i] <- predict(model, s = model$lambda.min, newx=x.test, type = "response")
    
  }
  
  resp <- data.frame(cbind(response, colnames(data)[-c(1,2)]))
  colnames(resp) <- c("response", "real_class")
  resp
}
