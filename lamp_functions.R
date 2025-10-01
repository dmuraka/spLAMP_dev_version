
################################################################
####################### Install packages #######################
################################################################

pkgs       <- c("FNN", "fields", "nloptr", "dbscan", "ranger")
to_install <- setdiff(pkgs, rownames(installed.packages()))
if(length(to_install) > 0){
  install.packages(to_install)
}

library(FNN)
library(fields)
library(nloptr)
library(dbscan)
library(ranger)


################################################################
####################### Internal functions #####################
################################################################

### Initialization
initial_fun <- function(x, y, coords, x_sel=NULL, train_rat, func){
  n              <- length(y)
  if(train_rat<1){
    set.seed(4321)
    id_train     <- sample(n, round(n*train_rat))
    set.seed(1234)
  } else {
    id_train     <- 1:n
  }
  
  if(is.null(x_sel)) x_sel <- apply(x,2,sd)!=0
  xname          <- c("Intercept", names(data.frame(x[,x_sel])))
  x              <- as.matrix(cbind(1,x[,x_sel]))
  coords         <- as.matrix(coords)
  nx             <- ncol(x)
  xx_inv         <- solve(t(x)%*%x)
  beta_int       <- xx_inv %*% t(x)%*%y
  pred           <- x%*%beta_int
  resid          <- y - pred
  sig2           <- sum(resid^2)/(n-nx)
  beta_int_vcov  <- sig2*xx_inv
  beta           <- matrix(beta_int  , nrow = n, ncol = nx, byrow = TRUE)
  beta_v         <- matrix(diag(beta_int_vcov), nrow = n, ncol = nx, byrow = TRUE)
  return(list(xx_inv=xx_inv, beta_int=beta_int, x=x, id_train=id_train,
              beta=beta, beta_v=beta_v, pred=pred, resid=resid, n=n, nx=nx,
              x_sel=x_sel,xname=xname,coords=coords))
}

### Kernel
kfun        <- function(dist, band, kernel){
  if(kernel=="gau"){
    wei   <- exp(-dist^2/band^2)
  } else if(kernel=="exp"){
    wei   <- exp(-dist/band)
  }
  return(wei)
}

### Minor coefficients tuning
bopt_core   <- function(par, bands, BBBB, beta_int, is_vc, 
                       nx, x, y, n_bid, id_train=NULL) {
  xbeta    <- matrix(0, nrow = nrow(x), ncol = nx)
  for(j in 1:nx) {
    if(is_vc[j]) {
      w    <- exp(-par[j] / bands)
      w    <- w / w[1]
      bbb  <- Reduce(`+`, lapply(1:n_bid, function(i) w[i] * BBBB[[i]][, j]))
      xbeta[, j] <- x[, j] * (beta_int[j, 1] + bbb)
    } else {
      xbeta[, j] <- x[, j] * (beta_int[j, 1])
    }
  }
  resid      <- y - rowSums( xbeta[, !is_vc, drop = FALSE] )
  xbeta_tt   <- xbeta[ , is_vc, drop = FALSE]
  vpar       <- solve(crossprod(xbeta_tt), crossprod(xbeta_tt, resid))
  if(!is.null(id_train)){
    xbeta_test <- xbeta[-id_train, is_vc, drop = FALSE]
    sse        <- sum((resid[-id_train] - as.vector(xbeta_test %*% vpar))^2)
  } else {
    sse        <- sum((resid - as.vector(xbeta_tt %*% vpar))^2)
  }
  return( list(sse = sse, vpar = vpar ) )
}

### Single-scale spatial process modeling via local model aggregation (gPoE)
lwr         <- function(coords, resid, x, band, b_old, i_hat, vc,ridge, coords_old=NULL,
                        kernel,hetero,id_train,y,beta=NULL,coords0,x0, 
                        sel_id=NULL,func="lamp"){
  
  n            <- nrow(coords)
  nx           <- ncol(x)
  if(kernel=="gau"){
    threshold  <- sqrt(-log(0.05))*band ## Gaussian kernel:0.01 is the kernel value at the threshold range
  } else if(kernel=="exp"){
    threshold  <- -log(0.05)*band
  }
  
  if(is.null(sel_id)){
    area         <- (max(coords[,1])-min(coords[,1]))^2 + (max(coords[,2])-min(coords[,2]))^2
    n_knot       <- round(1.5*area/band^2)
    if( n_knot < n ){
      set.seed(4321)
      coords_k_tmp<- kmeans(coords,n_knot,iter.max=5)$centers
      set.seed(1234)
      sel_id     <- get.knnx(coords, coords_k_tmp, 1)$nn.index
      coords_cent<- coords[sel_id,]
      sel_list   <- 1:nrow(coords_cent)
    } else {
      n_knot     <- n
      coords_cent<- coords
      sel_list   <- 1:n_knot
      sel_id     <- NA
    }
  } else if(is.na(sel_id[1])){
    n_knot       <- n
    coords_cent  <- coords
    sel_list     <- 1:n_knot
  } else {
    n_knot       <- length(sel_id)
    coords_cent  <- coords[sel_id,]
    sel_list     <- 1:n_knot
  }
  
  ################# Prior coefficient variance
  B_var          <- matrix(Inf,nrow=n_knot,ncol=nx)
  if( !is.null(b_old) & ridge==TRUE ){
    if(!is.null(coords_old)&hetero==TRUE){
      n_knot2    <- min(n_knot-1,nrow(coords_old)-1,20)
      knnx_0     <- get.knnx(coords_old,coords_cent,n_knot2)
      for(i in 1:nx){
        B_var0   <- matrix(b_old[,i][knnx_0$nn.index],ncol=n_knot2)
        wei_var0 <- exp(-knnx_0$nn.dist^2/band^2)
        wei_var  <- wei_var0/rowSums(wei_var0)
        B_var[,i]<- rowSums(wei_var * B_var0^2)
      }
    } else {
      for(i in 1:nx) B_var[,i]<- mean(b_old[,i]^2)
    }
  }
  
  ################# gPoE-based coefficient evaluation
  b_all        <- matrix(0,n,nx) # coefficients         : (wei/V[y])*b
  bv_inv_all   <- matrix(0,n,nx) # coefficient variances
  pv_inv_all   <- matrix(0,n,nx) # predictive variances : wei/v[y]
  b_old        <- matrix(0,length(sel_list),nx) # reset b_old
  if(!is.null(coords0)){
    n0         <- nrow(coords0)
    b_all0     <- matrix(0,n0,nx) # coefficients         : (wei/V[y])*b
    bv_inv_all0<- matrix(0,n0,nx) # coefficient variances
    pv_inv_all0<- matrix(0,n0,nx) # predictive variances : wei/v[y]
  }
  
  ### Estimate individual model and aggregate them (gPoE)
  id_train_flag<- logical(n)
  id_train_flag[id_train] <- TRUE
  query        <- coords_cent[sel_list,,drop=FALSE]
  dbnn         <- frNN(x = coords, query = query, eps = threshold, sort = FALSE) 
  if (!is.null(coords0)) {
    dbnn0      <- frNN(x = coords0, query = query, eps = threshold, sort = FALSE)
  }
  
  for(sel in sel_list){
    samp         <- dbnn$id[[sel]]
    flag         <- id_train_flag[samp]
    samp_hv      <- samp[flag]
    if( length(samp_hv) <= 5 ) next
    dist         <- dbnn$dist[[sel]]
    wei          <- kfun(dist=dist, band=band, kernel = kernel)
    wei_hv       <- wei[flag]
    wx_sel       <- wei_hv*x[samp_hv, , drop = FALSE]
    wxy_sel_csum <- crossprod(wx_sel, wei_hv*resid[samp_hv])
    wxxw_sel_csum<- colSums(wx_sel * wx_sel)
    b_sel0       <- as.vector(wxy_sel_csum)/ wxxw_sel_csum
    b_old[sel, ] <- b_sel0
    
    samp0        <- wei0 <- NULL
    if(!is.null(coords0)){
      samp0      <- dbnn0$id[[sel]]
      dist0      <- dbnn0$dist[[sel]]
      if( length(samp0)>0 ){
        wei0     <- kfun(dist=dist0, band=band, kernel = kernel)
      }
    }

    x_samp          <- x[samp, , drop = FALSE]
    wei_sq          <- wei^2
    for(j in (1:nx)[vc] ){
      resid_sub     <- resid[samp] - x_samp[,j]*b_sel0[j]
      sigma         <- sum( (wei*resid_sub)^2 ) /(length(samp)-1)
      lambda        <- sigma / B_var[sel,j]
      wxxw_lambda   <- wxxw_sel_csum[j] + lambda
      b_sel         <- wxy_sel_csum[j]/wxxw_lambda
      bv_sel        <- sigma/wxxw_lambda
      
      pv_sel        <- x[samp,j]^2/wxxw_sel_csum[j]*sigma + sigma/wei_sq # Var[yhat] sigma- x[samp,j]^2/wxxw_sel_csum[j]*sigma# V[y]=x%*%[sigma*solve(xwwx)]%*%t(x)+sigma
      wei2_pv_sel   <- wei_sq/pv_sel                           # wei^2/V[y]
      b_all[samp,j] <- b_all[samp,j] + wei2_pv_sel*b_sel       # (wei^2/V[y])*b
      bv_inv_all[samp,j]<- bv_inv_all[samp,j] + wei2_pv_sel/bv_sel#wei^2/V[b]
      pv_inv_all[samp,j]<- pv_inv_all[samp,j] + wei2_pv_sel    # wei^2/V[y]

      if( !is.null(coords0) & length(samp0)>0){
        pv_sel0     <- x0[samp0,j]^2/wxxw_sel_csum[j]*sigma + sigma/wei0^2#sigma- x[samp,j]^2/wxxw_sel_csum[j]*sigma# V[y]=x%*%[sigma*solve(xwwx)]%*%t(x)+sigma
        wei2_pv_sel0<- wei0^2/pv_sel0
        b_all0[samp0,j] <- b_all0[samp0,j] + wei2_pv_sel0*b_sel        # (wei^2/V[y])*b
        bv_inv_all0[samp0,j]<- bv_inv_all0[samp0,j] + wei2_pv_sel0/bv_sel#wei^2/V[b]
        pv_inv_all0[samp0,j]<- pv_inv_all0[samp0,j] + wei2_pv_sel0    # wei^2/V[y]
      }
    }
  }
  
  ################# selection of vc through CV
  run            <- FALSE
  if( func == "lamp_hv" ){
    b_all_sel    <- b_all[-id_train,vc,drop=FALSE]/pv_inv_all[-id_train,vc,drop=FALSE]
    b_all_sel[is.nan(b_all_sel)]<-0
    pred_hv      <- rowSums(x[-id_train,vc,drop=FALSE]*b_all_sel)
    resid_hv     <- resid[-id_train]-pred_hv
    sse_hv0      <- sum((resid[-id_train])^2)
    sse_hv       <- sum((resid_hv)^2)
    run          <- ifelse(sse_hv<sse_hv0, TRUE, FALSE)
    if(run){
      beta_new     <- beta[-id_train,  ,drop=FALSE]
      beta_new[,vc]<- beta[-id_train,vc,drop=FALSE] + b_all_sel# coefficient (resol: 1:R)
    } else {
      sse_hv     <- sse_hv0
    }
  } else {
    sse_hv      <- NA
    run         <- TRUE
  }
  
  if( run ){
    bv_all          <- bv_inv_all
    b_all[,vc]      <- b_all[,vc]/pv_inv_all[,vc]
    b_all[,-vc]         <- 0
    b_all[is.nan(b_all)]<- 0
    
    bv_inv_all[, vc]<- bv_inv_all[, vc]/pv_inv_all[, vc]
    bv_all[, vc]    <- 1/bv_inv_all[, vc]
    bv_all[, -vc]   <- NA
    bv_all[is.nan(bv_all)]<-Inf
    
    pred            <- rowSums(x*b_all)
    if( !is.null(coords0) ){
      bv_all0          <- bv_inv_all0
      b_all0[,vc]      <- b_all0[,vc]/pv_inv_all0[,vc]
      b_all0[,-vc]     <- 0
      b_all0[is.nan(b_all0)]<-0
      b_all0[is.na(b_all0)] <-0
      
      bv_inv_all0[, vc]<- bv_inv_all0[, vc]/pv_inv_all0[, vc]
      bv_all0[, vc]    <- 1/bv_inv_all0[, vc]
      bv_all0[, -vc]   <- NA
      bv_all0[is.nan(bv_all0)]<-Inf
      pred0       <- rowSums(x0*b_all0)

    } else {
      b_all0 <-bv_all0<-pred0<-NULL
    }
    
    return(list(beta=b_all, beta_v=bv_all, pred=pred, sel_id=sel_id,
                coords_cent=coords_cent,B_var=B_var,
                beta0=b_all0,beta0_v=bv_all0, pred0=pred0, 
                b_old=b_old, run=run,sse_hv=sse_hv,vc_sel=vc))#sse_hv0=sse_hv0, 
  } else {
    return(list(run=FALSE,B_var=B_var))
  }
}


################################################################
######################### Main functions #######################
################################################################

###########################################################
### Holdout validation optimizing the number of resolutions 
#
# lamp_hv(y, x, coords, train_rat=0.75,alpha=0.9, ridge=TRUE, kernel="exp", rf=FALSE)
#
### input
# y        : Data vector
# x        : Matrix of covariates
# coords   : Matrix of 2-dimensional point coordinates
# train_rat: Ratio of training samples
# alpha    : Decay ratio of the bandwidth values over resolutions
# ridge    : If TRUE, ridge parameter is used to regularize each local model
# kernel   : Kernel function for local model weighting 
#            ("exp": exponential kernel, "gau": Gaussian kernel)
# rf       : If TRUE, random forest is additionally trained to capture
#            non-linear patterns and/or higher-order interactions
#
### output
# param    : Internal parameters
# sse_hv   : Sum of squares error for the validation samples

lamp_hv     <- function(y, x, coords, train_rat=0.75,
                        alpha=0.9, ridge=TRUE, kernel="exp", rf=FALSE){
  init           <- initial_fun(x=x,y=y,coords=coords,train_rat=train_rat,
                                x_sel=NULL,func="lamp_hv")
  xx_inv         <- init$xx_inv
  beta_int       <- init$beta_int
  beta           <- init$beta
  coords         <- init$coords
  pred           <- init$pred
  resid          <- init$resid
  x              <- init$x
  x_sel          <- init$x_sel
  xname          <- init$xname
  n              <- init$n
  nx             <- init$nx
  id_train       <- init$id_train
  vc             <- 1
  hetero         <- FALSE
  
  BBB            <- list(NULL)
  max_d          <- sqrt(diff(range(coords[,1]))^2+diff(range(coords[,2]))^2)/3
  Bands          <- max_d*alpha^(1:100)
  accept_num     <- 5
  
  ##################### main loop for feature extraction
  print("--- SSE: Learning multi-resolution spatial processes ---")
  coords_old     <- NULL
  sel_id_list   <- list(NULL)
  b_old          <- NULL
  bands          <- NULL
  SSE            <- NULL
  count          <- 0
  VCmat          <- NULL
  for(i in 1:length(Bands)){
    band         <- Bands[i]
    lmod         <- lwr(coords=coords, resid=resid, x=x, band=band, coords_old=coords_old,
                        b_old=b_old,vc=vc, id_train=id_train,ridge=ridge, kernel=kernel,
                        beta=beta,y=y, coords0=NULL, x0=NULL, sel_id=NULL,
                        hetero=hetero,func="lamp_hv")
    run          <- lmod$run
    if(run==TRUE){
      bands      <- c(bands, band)
      b_old      <- lmod$b_old
      coords_old <- lmod$coords_cent
      
      beta_add   <- lmod$beta
      pred_add   <- lmod$pred
      beta       <- beta + beta_add
      pred       <- pred + pred_add
      resid      <- y - pred
      SSE        <- c(SSE ,lmod$sse_hv)
      vc_sel     <- lmod$vc_sel
      vcmat      <- rep(0,nx);vcmat[vc_sel]<-1
      VCmat      <- rbind(VCmat,vcmat)
      
      beta_int_add  <- xx_inv %*% t(x)%*%resid
      pred0_add  <- x%*%beta_int_add
      beta       <- sweep(beta, 2, beta_int_add, "+")
      pred       <- pred  + pred0_add
      resid      <- resid - pred0_add
      
      beta_add_m <- colMeans(beta_add)
      BBB[[i]]   <- sweep(beta_add, 2, beta_add_m, "-") # centered process
      sel_id_list[[i]]<- lmod$sel_id
      beta_int   <- beta_int + beta_int_add + beta_add_m# de-centered coefficients
      count      <- 0
    } else {
      VCmat      <-rbind(VCmat,rep(0,nx))
      SSE        <-c(SSE, SSE[length(SSE)])
      
      if(i>5) count      <- count + 1
      if(count==accept_num) break
    }
    
    print_add<-ifelse(i<10,"  "," ")
    print( paste0( formatC(SSE[length(SSE)], digits = 7, format = "g", flag = "#"),
                   " (Resolution",print_add, i,")"), quote = FALSE )
  }
  
  print("--- SSE: Coefficient adjustment ---")
  bid        <- which(sapply(BBB, length) > 0)
  n_bid      <- length(bid)
  is_vc      <- (1:nx) %in% vc
  if(n_bid>0){
    BBBB     <- BBB[bid]
    bopt_obj <- (function(bands, BBBB, beta_int, is_vc, nx, 
                          x, y, n_bid, id_train) {
      function(par) {
        out <- try(bopt_core(par, bands = bands, BBBB = BBBB,
                             beta_int = beta_int, is_vc = is_vc, nx = nx,
                             x = x, y = y, n_bid = n_bid, id_train=id_train),
                   silent = TRUE)
        if (inherits(out, "try-error") || !is.finite(out$sse)) {
          return(.Machine$double.xmax)
        }
        out$sse
      }
    })(bands, BBBB, beta_int, is_vc, nx, x, y, n_bid, id_train)
    
    v_opt0   <- nloptr(x0 = rep(0, sum(is_vc)),eval_f = bopt_obj,
                       opts = list(algorithm = "NLOPT_LN_BOBYQA", maxeval = 500))
    v_test   <- bopt_core(v_opt0$solution, bands=bands, BBBB=BBBB, 
                          beta_int=beta_int, is_vc=is_vc, nx=nx, 
                          x=x, y=y, n_bid=n_bid,id_train=id_train)
    if(v_test$sse< SSE[length(SSE)]){
      vpar   <- c(v_test$vpar[1], v_opt0$solution)
    } else {
      vpar   <- c(1, 0)
      print("no adjustment")
    }
    
  } else {
    vpar       <- c(NA,NA)
    message("Warning: No residual spatial process was detected.")
  }
  
  xbeta      <- matrix(0,nrow=n,ncol=nx)
  for(j in 1:nx){
    if(is_vc[j]&!is.na(vpar[j])){
      w_0     <- exp(-vpar[sum(is_vc)+j]/bands)
      w       <- vpar[j]* w_0/w_0[1]#vpar[j]^2 fixed
      w[w<0]  <- 0
      b       <- Reduce("+", lapply(1:n_bid, function(i) w[i]*BBB[bid][[i]][, j]))
      xbeta[,j]<- x[,j]*(b + beta_int[j,1])
    } else {
      xbeta[,j]<- x[,j]*(beta_int[j,1])
    }
  }
  
  pred       <- rowSums(xbeta)
  resid      <- y-pred
  sse_hv     <- sum( resid[-id_train]^2 )
  SSE        <- c(SSE,sse_hv)
  print(sse_hv)
  
  rpar       <- data.frame(mtry=NA, min.node.size=NA)
  if(rf==TRUE){
    print("--- SSE: Additional learning (ranfom forest) ---")
    rf_data    <- data.frame(resid=resid,x[,-1],coords)
    mtry_all   <- c( round( (nx+1)/5), round( (nx+1)/3), round( (nx+1)/2))
    param_grid <- expand.grid(mtry = unique(mtry_all),min.node.size = c(1, 5, 10))
    param_grid <- rbind(data.frame(mtry=NA,min.node.size=NA),param_grid)
    
    for (i in 2:nrow(param_grid)){
      params    <- param_grid[i, ]
      rf_mod    <- ranger(formula = resid ~ .,data = rf_data[id_train,],
                          classification = FALSE,probability = FALSE, 
                          verbose = FALSE, mtry = params$mtry, num.trees = 500,
                          min.node.size = params$min.node.size, seed=123)
      resid_rf  <- resid[-id_train] - predict(rf_mod, data=rf_data[-id_train,])$predictions
      sse_rf    <- sum( resid_rf^2 )
      if(sse_rf <  sse_hv){
        sse_hv  <- sse_rf
        rpar    <- params
      }
    }
    SSE        <- c(SSE,sse_hv)
    print(sse_hv)
  }
  
  ##################### summary
  param      <- list(bands=bands,vpar=vpar,rpar=rpar,alpha=alpha,ridge=ridge,
                     id_train_hv=id_train,vc=vc,x_sel=x_sel,SSE=SSE, 
                     sel_id_list=sel_id_list,
                     VCmat=VCmat,hetero=hetero,kernel=kernel,rf=rf)
  
  return(list(param=param,sse_hv=sse_hv))
}

######################################################
### Local aggregate multiscale process (LAMP) modeling 
#
# lamp(y, x, coords, x0=NULL, coords0=NULL, mod_hv)
#
### input
# y        : Data vector
# x        : Matrix of covariates
# coords   : Matrix of 2-dimensional point coordinates
# x0       : Matrix of covariates at prediction sites
# coords0  : Matrix of 2-dimensional point coordinates at prediction sites
# mod_hv   : Output object of the lamp_hv function
#
### output
# beta     : Regression coefficients and their standard errors
# pred     : Predictive values (sample sites; vector)
# pred0    : Predictive values (prediction sites; vector)
# BBB      : Mean of the spatial process in each resolution (sample sites; list) 
# BBBv     : Variance of the spatial process in each resolution (sample sites; list)
# BBB0     : Mean of the spatial process in each resolution (prediction sites; list) 
# BBB0v    : Variance of the spatial process in each resolution (prediction sites; list) 
# param    : Internal parameters

lamp        <- function(y, x, coords, x0=NULL, coords0=NULL, mod_hv){
  
  bands          <- mod_hv$param$bands
  vpar           <- mod_hv$param$vpar
  rpar           <- mod_hv$param$rpar
  sel_id_list    <- mod_hv$param$sel_id_list
  alpha          <- mod_hv$param$alpha
  ridge          <- mod_hv$param$ridge
  vc             <- mod_hv$param$vc
  x_sel          <- mod_hv$param$x_sel
  SSE            <- mod_hv$param$SSE
  VCmat          <- mod_hv$param$VCmat
  hetero         <- mod_hv$param$hetero
  kernel         <- mod_hv$param$kernel
  rf             <- mod_hv$param$rf
  
  init           <- initial_fun(x=x,y=y,coords=coords,x_sel=x_sel,func="lamp",train_rat=1)
  xx_inv         <- init$xx_inv
  beta_int       <- init$beta_int
  beta           <- init$beta
  coords         <- init$coords
  pred           <- init$pred
  resid          <- init$resid
  x              <- init$x
  x_sel          <- init$x_sel
  xname          <- init$xname
  n              <- init$n
  nx             <- init$nx
  id_train       <- init$id_train
  
  if(!is.null(coords0)){
    n0           <- nrow(coords0)
    x0           <- as.matrix(cbind(1,x0[,x_sel]))
    pred0        <- x0 %*% beta_int
  } else {
    n0   <- x0   <- NA
    pred0        <- NULL
  }
  
  ##################### main loop for feature extraction
  b_old          <- NULL
  BBB   <- BBBv  <- list(NULL)
  BBB0  <- BBB0v <- list(NULL)
  if(!is.null(bands)){
    for(i in 1:length(bands)){
      vc           <- which(VCmat[i,]==1)
      lmod         <- lwr(coords=coords, resid=resid, x=x, band=bands[i], hetero=hetero,
                          b_old=b_old, vc=vc, id_train=id_train,ridge=ridge,kernel=kernel,
                          x0=x0, coords0=coords0, sel_id=sel_id_list[[i]], func="lamp")
      beta_add     <- lmod$beta
      beta_v_add   <- lmod$beta_v
      beta_v_add[is.infinite(beta_v_add)]<-0
      pred_add     <- lmod$pred
      b_old        <- lmod$b_old
      pred         <- pred + pred_add
      resid        <- y - pred
      beta_int_add <- xx_inv %*% t(x)%*%resid
      pred_int_add <- x%*%beta_int_add
      pred         <- pred  + pred_int_add
      resid        <- resid - pred_int_add
      
      beta_add_m   <- colMeans(beta_add)
      BBB[[i]]     <- sweep(beta_add, 2, beta_add_m, "-")
      BBBv[[i]]    <- beta_v_add
      beta_int     <- beta_int + beta_int_add + beta_add_m
      
      if(!is.null(coords0)){
        beta0_add     <- lmod$beta0
        beta0_v_add   <- lmod$beta0_v
        beta0_v_add[is.infinite(beta0_v_add)]<-0#tentative
        pred0_add     <- lmod$pred0
        pred0         <- pred0 + pred0_add
        pred0_int_add <- x0 %*% beta_int_add
        pred0         <- pred0 + pred0_int_add
        
        BBB0[[i]]     <- sweep(beta0_add, 2, beta_add_m, "-")
        BBB0v[[i]]    <- beta0_v_add
      }
      print(i)
    }
  } else {
    message("Warning: No residual spatial process was modeled")
  }
  
  pred_pre       <- rowSums(x*beta)
  
  ##################### tuning
  beta           <- matrix(beta_int[,1], nrow = n, ncol = nx, byrow = TRUE)
  if(!is.null(coords0)){
    beta0        <- matrix(beta_int[,1], nrow=n0,ncol=nx, byrow=TRUE)
  }
  
  n_band_x       <- apply(VCmat,2,function(x) sum(x==1))
  n_bid          <- length(bands)
  vpar_coef      <- bopt_core(vpar[2], bands=bands, BBBB=BBB, 
                        beta_int=beta_int, is_vc=ifelse(n_band_x>0,1,0), nx=nx, 
                        x=x, y=y, n_bid=n_bid,id_train=NULL)$vpar[1]
  for(j in 1:nx){
    if(n_band_x[j]>0){
      w_0        <- exp(-vpar[sum(n_band_x>0) + j]/bands)
      w          <- vpar_coef*w_0/w_0[1]#vpar[j]
      w[w<0]     <-0
      b          <- Reduce("+", lapply(1:n_band_x[j], function(i) w[i]*BBB[[i]][, j]))
      beta[,j]   <- beta[,j] + b
      if(!is.null(coords0)){
        b0       <- Reduce("+", lapply(1:n_band_x[j], function(i) w[i]*BBB0[[i]][, j]))
        beta0[,j]<- beta0[,j] + b0
      }
    }
  }
  
  pred           <- rowSums(x*beta)
  beta           <- as.data.frame(beta)
  names(beta)    <- xname
  rf_mod         <- rf_xname<-NULL
  if(rf==TRUE&!is.na(rpar[1])){
    rf_data     <- data.frame(resid=y-pred,x[,-1],coords)
    names(rf_data)[-1]<-c(xname[-1],"px","py")
    rf_xname    <- names(rf_data)[-1]
    rf_mod      <- ranger(formula = resid ~ .,data = rf_data,quantreg = TRUE,# keep.inbag = TRUE,
                         classification = FALSE,probability = FALSE, 
                         verbose = FALSE, mtry = rpar$mtry, num.trees = 500,
                         min.node.size = rpar$min.node.size, seed=123)
    pred        <- pred + rf_mod$predictions
  }
  
  if(!is.null(coords0)){
    pred0       <- rowSums(x0*beta0)
    if(rf==TRUE&!is.na(rpar[1])){
      rf_data0  <- data.frame(x0[,-1],coords0)
      names(rf_data0)<-rf_xname
      pred0     <- pred0 + predict(rf_mod, data=rf_data0)$predictions
    }
    beta0       <- as.data.frame(beta0)
    names(beta0)<- xname
  }
  
  ######### coefficients
  sig_pre       <- sum( (y - pred_pre)^2)/(n-nx)
  v_diag        <- Reduce("+",BBBv)[,1] + sig_pre
  beta_int_se   <- sqrt(diag(solve(crossprod(x, 1/v_diag * x))))
  beta_int_summ <- data.frame(coef=beta_int,coef_se=beta_int_se) 
  
  ######### summary parameters
  param2<-list(n=n,n0=n0,nx=nx,y=y,x=x,x0=x0,VCmat=VCmat,bands=bands,
               coords=coords,coords0=coords0,vpar=vpar,vc=mod_hv$param$vc, 
               xx_inv=xx_inv, rf_mod=rf_mod, pred_pre=pred_pre,rf_xname=rf_xname)
  
  return(list(beta=beta_int_summ,pred=pred,pred0=pred0,
              BBB=BBB,BBBv=BBBv, BBB0=BBB0, BBB0v=BBB0v, param=param2))
}

######################################################
### Simulation evaluating standard deviation and percentiles of the prediction
#
# lamp_sim(mod, nsim=200, pred_sample=FALSE, pred0_sample=FALSE,
#          probs=c(0.025,0.5,0.975))
#
### input
# mod         : Output object of the lamp function   
# nsim        : Number of Monte Carlo simulations
# pred_sample : If TRUE, simulated predictive values at sample sites are returned
# pred0_sample: If TRUE, simulated predictive values at prediction sites are returned
# probs       : Prediction percentiles. By default, 
#               2.5, 50.0, and 97.5 percentiles are evaluated
#
### output
# pred        : Matrix of predictive mean and standard deviation (sample sites)
# pred0       : Matrix of predictive mean and standard deviation (prediction sites)
# pred_qt     : Matrix of prediction percentiles at probs (sample sites)
# pred0_qt    : Matrix of prediction percentiles at probs (prediction sites) 
# Pred_sim    : Matrix of nsim simulated predictive values (sample sites)
# Pred0_sim   : Matrix of nsim simulated predictive values (prediction sites)

lamp_sim    <- function(mod, nsim=200, pred_sample=FALSE, pred0_sample=FALSE,
                        probs=c(0.025,0.5,0.975)){
  
  n            <- mod$param$n
  n0           <- mod$param$n0
  nx           <- mod$param$nx
  y            <- mod$param$y
  x            <- mod$param$x
  x0           <- mod$param$x0
  coords       <- mod$param$coords
  coords0      <- mod$param$coords0
  pred         <- mod$pred
  pred0        <- mod$pred0
  resid        <- y - pred
  sig          <- sum(resid^2)/(n-nx)
  bands        <- mod$param$bands
  beta_int     <- mod$beta$coef
  beta_int_se  <- mod$beta$coef_se
  vpar         <- mod$param$vpar
  vc           <- mod$param$vc
  rf_mod       <- mod$param$rf_mod
  rf_xname     <- mod$param$rf_xname
  n_band_x     <- apply(mod$param$VCmat,2,function(x) sum(x==1))
  
  ######### Functions
  sim_lamp_fun <- function(BBB, BBBv, x, n, nx,nsim,n_band_x, beta_int, beta_int_se,
                         vpar,vc,bands,sig){
    BBB_sim      <- BBB
    Pred_sim     <- matrix(0, nrow = n , ncol = nsim, byrow = TRUE)
    for(i in 1:nsim){
      for(j in 1:nx){    ####### simulating spatial process
        if(n_band_x[j]>0){
          for(l in 1:n_band_x[j]){
            BBB_sim[[l]][,j]<-rnorm(n=n,mean=BBB[[l]][, j],sd=sqrt(BBBv[[l]][,j]))
          }
        }
      }
      
      beta_sim       <- matrix(0, nrow = n , ncol = nx, byrow = TRUE)
      for(j in 1:nx){
        beta_sim[,j] <- rnorm(n=1,mean=beta_int[j],sd=beta_int_se[j])
        if(n_band_x[j]>0){
          w_0       <- exp(-vpar[length(vc) + j]/bands)
          w         <- vpar[j]*w_0/w_0[1]
          w[w<0]    <- 0
          bsim        <- Reduce("+", lapply(1:n_band_x[j], function(l) w[l]*BBB_sim[[l]][,j]))
          beta_sim[,j]<- beta_sim[,j] + bsim
        }
      }
      
      ######### simulation for the predictive values
      Pred_sim[,i]   <- rowSums(x*beta_sim)  + rnorm(n=n ,sd=sqrt(sig))
    }
    return(list(Pred_sim=Pred_sim, beta_sim=beta_sim, BBB_sim=BBB_sim))
  }
  
  sample_from_qrf <- function(rf_qmat, n, n_draw = 100) {
    U      <- matrix(runif(n * n_draw), nrow = n, ncol = n_draw)
    draws  <- matrix(NA_real_, nrow = n, ncol = n_draw)
    for (i in 1:n) {
      qi   <- rf_qmat[i, ]
      if (is.unsorted(qi)) qi <- sort(qi)
      draws[i, ] <- approx(x = qs, y = qi, xout = U[i, ], ties = "ordered")$y
    }
    draws
  }
  
  pred_summary<- function(probs, n, pred, Pred_sim){
    pred_qt   <- NULL
    if( length( probs ) > 0 ){
      pred_qt <- matrix(NA,nrow=n,ncol=length(probs))
      for(pid in 1:length(probs)){
        pred_qt[,pid]<- apply(Pred_sim,1,function(x) quantile(x,probs=probs[pid]))
      }
      pred_qt        <- data.frame(pred_qt)
      names(pred_qt) <- paste0("q",probs)
    }
    pred_ms  <- data.frame(pred=pred,pred_se=apply(Pred_sim ,1,sd) ) 
    return( list(pred_ms=pred_ms,pred_qt=pred_qt) )
  }
  
  ######### Simulation for sapmle sites
  Pred_sim <- sim_lamp_fun(BBB=mod$BBB, BBBv=mod$BBBv, x=x, n=n, nx=nx, nsim=nsim, 
                        n_band_x=n_band_x,beta_int=beta_int, beta_int_se=beta_int_se,
                        vpar=vpar,vc=vc,bands=bands,sig=sig)$Pred_sim
  if(!is.null(rf_mod)){
    qs      <- seq(0, 1, length.out = 201)
    rf_dat  <- data.frame(x[,-1],coords)
    names(rf_dat)<-rf_xname
    rf_qmat <- predict(rf_mod, data = rf_dat, type = "quantiles", quantiles = qs)$predictions
    Pred_sim<- Pred_sim + sample_from_qrf(rf_qmat, n, n_draw = nsim)
  }
  pred_summ  <- pred_summary(probs=probs, n=n, pred=pred  , Pred_sim=Pred_sim)
  pred_ms    <- pred_summ$pred_ms
  pred_qt    <- pred_summ$pred_qt
  
  ######### Simulation for prediction sites
  pred0_ms   <- pred0_qt <- Pred0_sim <- NULL
  if( !is.na(n0) ){
    Pred0_sim  <- sim_lamp_fun(BBB=mod$BBB0, BBBv=mod$BBB0v, x=x0, n=n0, nx=nx, nsim=nsim, 
                              n_band_x=n_band_x,beta_int=beta_int, beta_int_se=beta_int_se,
                              vpar=vpar,vc=vc,bands=bands,sig=sig)$Pred_sim
    if(!is.null(rf_mod)){
      rf_dat0  <- data.frame(x0[,-1],coords0)
      names(rf_dat0)<- rf_xname
      rf_qmat0 <- predict(rf_mod, data = rf_dat0, type = "quantiles", quantiles = qs)$predictions
      Pred0_sim<- Pred0_sim + sample_from_qrf(rf_qmat0, n0, n_draw = nsim)
    }
    pred0_summ <- pred_summary(probs=probs, n=n0, pred=pred0, Pred_sim=Pred0_sim)
    pred0_ms   <- pred0_summ$pred_ms
    pred0_qt   <- pred0_summ$pred_qt
  }
  
  return(list(pred=pred_ms,pred0=pred0_ms,pred_qt=pred_qt,pred0_qt=pred0_qt,
              Pred_sim=Pred_sim,Pred0_sim=Pred0_sim))
}

