# clear environment
rm(list = ls())

# load libraries
library(foreach)
library(doParallel)
library(matlib)  # for Ginv
library(pracma) # for geo_median

######################################
########## Source Functions ##########
######################################


trace <- function(A){ sum(diag(A))}


# pilot_data: matrix containing n rows of pilot data, p predictors, and response
# candidates: N by p matrix of candidate points for future design. Does not include intercept.
# initial_inds: list of n0 < N indices of initial points that correspond to rows in candidates. Ideally equally spaced in candidate region. 
# formula: formula for a GLM. Must use variable names from pilot_data
# family: one of "binomial" or "poisson"
# B: Number of bootstrap iterations, default = 100.
# tries up to n_reboot times (5)
# returns a list of 4 elements:
# [[1]] opt weights for untrimmed BAG
# [[2]] list of betahats (candidate parameter vectors)
# [[3]] opt weights for BAG trimmed at 0.05
# [[4]] opt weights for BAG trimmed at 0.025
# [[5]] opt weights for BAG with geo_median
get_bag_design_weights <- function(pilot_data, candidates, initial_inds, formula, family, B = 100,
                                     n_reboot = 10, max_iter = 100){
  n <- nrow(pilot_data)
  N <- nrow(candidates)
  weights <- matrix(nrow = B, ncol = N)
  betahats <- matrix(nrow = B, ncol = ncol(pilot_data)) 
  # ncol is ok because pilot_data has dimension p+1 due to Y
  for(i in 1:B){
    print(paste("Bootstrap Sample", i, sep = " "))
    found_sample <- FALSE
    reboot_counter <- 0
    while(!found_sample){
      sample_inds <- sample(1:n, n, replace = TRUE)
      samplei <- pilot_data[sample_inds,]
      modeli <- glm(formula, family, data = samplei)
      betahat <- modeli$coefficients
      p <- length(betahat)
      if(all(!is.na(betahat)) & modeli$converged){
        # found_sample <- TRUE
        betahats[i,] <- betahat
        print(betahat)
        tryCatch({
          out1 <- get_opt_design(initial_inds, candidates, betahats[i,], family, max_iter = max_iter)
          locally_opt_inds <- out1[[1]]
          locally_opt_weights <- numeric(N)
          locally_opt_weights[locally_opt_inds] <- out1[[2]]
          weights[i,] <- locally_opt_weights
        }, error=function(e){cat("Local Design Not Found for BAG:",conditionMessage(e), "\n")})
        
        # print("finding opt design")
        # out1 <- get_opt_design(initial_inds, candidates, betahats[i,], family)
        # print("found opt deisgn")
        # locally_opt_inds <- out1[[1]]
        # locally_opt_weights <- numeric(N)
        # locally_opt_weights[locally_opt_inds] <- out1[[2]]
        # weights[i,] <- locally_opt_weights
        
        if(!is.na(weights[i,1])){
          found_sample <- TRUE
        }
        else{
          reboot_counter <- reboot_counter + 1
          # print(paste("Reboot Attempt", reboot_counter, sep = " "))
          if(reboot_counter > n_reboot){
            # print("Bootstrap sample discarded")
            weights[i,] <- NA  # this will be ignored in the optimal weights since na.rm = TRUE
            break  # move on to the next iteration
          }
        }
        
      }
      
    }
    
    
    
    
    
    
  }
  
  # opt_weights <- round(colMeans(weights, na.rm = TRUE),2)
  opt_weights <- colMeans(weights, na.rm = TRUE)
  
  # try "cleaning up" some of these weights
  # if the weight is less than a tolerance, set it to 0 and re-normalize
  opt_weights2 <- opt_weights
  opt_weights3 <- opt_weights
  opt_weights2[opt_weights < 0.02] <- 0
  opt_weights3[opt_weights < 0.01] <- 0
  opt_weights2 <- opt_weights2/sum(opt_weights2)
  opt_weights3 <- opt_weights3/sum(opt_weights3)
  
  opt_weights4 <- geo_median(weights)$p
  opt_weights4 <- opt_weights4/sum(opt_weights4)
  opt_weights5 <- opt_weights4
  opt_weights6 <- opt_weights5
  opt_weights5[opt_weights4 < 0.02] <- 0
  opt_weights6[opt_weights4 < 0.01] <- 0
  opt_weights5 <- opt_weights5/sum(opt_weights5)
  opt_weights6 <- opt_weights6/sum(opt_weights6)
  
  return(list(opt_weights, betahats, opt_weights2, opt_weights3, opt_weights4,
              opt_weights5, opt_weights6))
  
}



# initial_inds: a vector of indices for n0 initial points (taken from candidates)
# candidates: an N times p matrix of candidate points (typically a regular grid)
# beta: vector of length p+1 of unknown regression coefficients
# family: one of "poisson" or "binomial"
# returns a list of two elements:
## [[1]] indicies for support points (from candidates)
## [[2]] weights for the corresponding support points
get_opt_design <- function(initial_inds, candidates, beta, family, tol = 1e-6,
                           max_iter = 200){
  t <- 1
  n0 <- length(initial_inds)
  old_weights <- rep(1/n0, n0)
  N <- nrow(candidates)
  current_inds <- initial_inds
  p <- length(beta)
  
  keep_iterating <- TRUE
  final_weights <- NA
  final_points <- NA
  
  # Xcandidates <- cbind(1,candidates)  # add 1 for intercept
  # etacands <- Xcandidates%*%beta
  
  while(keep_iterating){
    
    # print(paste("old weights: ", paste(round(old_weights,3), collapse = ",")))
    # print(paste("old points: ", paste(current_inds, collapse = ",")))
    
    newton1 <- newton_weights2(current_inds = current_inds, old_weights = old_weights,
                               beta = beta, candidates = candidates, family = family, max_iter = max_iter, tol = tol)
    new_inds <- newton1[[1]]
    new_weights <- newton1[[2]]
    new_points <- candidates[new_inds,]
    
    
    # print(paste("new weights: ", paste(round(new_weights,3), collapse = ",")))
    # print(paste("new points: ", paste(new_inds, collapse = ",")))
    
    
    I0 <- get_info_mat(new_points, beta, new_weights, family)
    I0inv <- Ginv(I0)
    sigma_inv <- I0
    # part <- I0inv%*%Ginv(I0inv)%*%I0inv
    n <- length(new_inds)
    dvec <- rep(-Inf, N)
    # diff1 <- trace(I0%*%part)
    remaining_cands <- setdiff(1:N, new_inds)
    for(i in remaining_cands){
      
      info_xi <- get_info_mat(candidates[i,], beta = beta, old_weights = 1, family = family)
      dvec[i] <- trace( sigma_inv%*%I0inv%*%(info_xi - I0)%*%I0inv)
      # dvec[i] <- trace(info_xi%*%part) - diff1
      # dvec[i] <- trace(info_xi%*%part) - p
      
      
    }
    
    max_ind <- which.max(dvec)
    # print(paste("Delta = ", dvec[max_ind], sep = ""))
    if( dvec[max_ind] <= tol){
      # this means the current design is optimal
      final_weights <- new_weights
      final_inds <- new_inds
      keep_iterating <- FALSE
    }
    else{
      # add index of best candidate to current inds
      current_inds <- c(new_inds, max_ind)
      # print(paste("Added point: ", max_ind))
      # update the weights
      old_weights <- c(new_weights, 0)
      # update counter
      t <- t+1
      # print(t)
      if(t > max_iter){
        # print(paste("Optimal design not found in ", max_iter, " iterations.", sep = ""))
        # stop(paste("Optimal design not found in ", max_iter, " iterations.", sep = ""))
        # return(list(NA, NA))
        warning(paste("Optimal design not found in ", max_iter, " iterations, returning best design found so far.", sep = ""))
        keep_iterating <- FALSE
        return(list(new_inds, new_weights))
      }
    }
    
    
  }
  
  
  return(list(final_inds, final_weights))
  
}



newton_weights2 <- function(current_inds, old_weights, beta, candidates, family, tol = 1e-6,
                            max_iter = 100){
  
  
  n <- length(current_inds)
  if(n > 1){
    iter <- 1
    weights <- newton_weights(current_inds = current_inds, old_weights = old_weights, beta=beta, 
                              candidates = candidates, family= family, tol=tol)
    inds <- current_inds
    while(min(weights) < tol & length(inds) > 1  & iter < max_iter){
      
      iter <- iter + 1
      # remove support point with smallest weight
      smallest_ind <- which.min(weights)
      if(weights[smallest_ind] < tol){
        # print(paste("Removing point: ", inds[smallest_ind]))
        inds <- inds[-smallest_ind]  # remove the smallest support point
        weights <- weights[-smallest_ind] # remove the corresponding weight for the smallest point
      }
      if(length(inds) > 1){
        weights <- newton_weights(current_inds = inds, old_weights = weights, beta=beta, 
                                  candidates = candidates, family= family, tol=tol, max_iter = max_iter)
      }
      
    }
  }
  else{
    weights <- 1
    inds <- current_inds
  }
  if(length(inds) == 1){
    weights <- 1
  }
  
  
  return( list(inds, weights) )
  
  
}


# helper function to find information matrix
get_info_mat <- function(support_points, beta, old_weights, family){
  if(is.vector(support_points)){
    X_full <- matrix(c(1,support_points), nrow = 1)
    X_full_weighted <- X_full # if there is only one support point, then weight = 1
  }
  else{
    X_full <- cbind(1, support_points)
    X_full <- as.matrix(X_full)
    X_full_weighted <- sweep(X_full, MARGIN = 1, old_weights, '*')
  }
  eta <- as.vector(X_full%*%beta)
  vs <- NA
  if(family == "binomial"){
    vs <- exp(eta)/((1+exp(eta))^2)
  }
  
  if(family == "poisson"){
    vs <- exp(eta)
  }
  V <- diag(vs)
  if(length(vs) == 1){
    I0 <- vs*t(X_full)%*%X_full_weighted
  }
  else{
    I0 <- t(X_full)%*%V%*%X_full_weighted
  }
  return(I0)
}


# newton_weights
# returns a list of new support point indices followed by new weights
newton_weights <- function(current_inds,old_weights, beta, candidates, family,
                           tol = 1e-6, max_iter = 100){
  new_weights <- NA
  new_inds <- current_inds
  alpha <- 1 # step size
  diff <- Inf 
  iter <- 1
  indic <- TRUE
  
  while((diff > tol)  & (iter < max_iter) & indic  ){
    support_points <- candidates[new_inds,]
    n <- length(new_inds)
    m <- length(new_inds)  # for now, n = m
    first_derivative <- numeric(m-1)
    second_derivative <- matrix(NA, nrow = m-1, ncol = m-1)
    
    I0 <- get_info_mat(support_points, beta, old_weights, family)
    I0inv <- Ginv(I0)
    sigma_inv <- I0
    
    # find first derivatives
    Ixm <- get_info_mat(support_points[m,], beta, old_weights, family)
    for(i in 1:(m-1)){
      Ixi <- get_info_mat(support_points[i,], beta, old_weights, family)
      # first_derivative[i] <- sum(diag((n*(Ixi - Ixm))%*%I0inv ))
      first_derivative[i] <- trace( -sigma_inv%*%I0inv%*%(n*(Ixi - Ixm))%*%I0inv)
    }
    
    # find second derivatives
    for(i in 1:(m-1)){
      
      Ixi <- get_info_mat(support_points[i,], beta, old_weights, family)
      Ii <- n*(Ixi - Ixm)
      for(j in 1:(m-1)){
        
        Ixj <- get_info_mat(support_points[j,], beta, old_weights, family)
        Ij <- n*(Ixj - Ixm)
        mat1 <- I0inv%*%Ij%*%I0inv%*%Ii%*%I0inv + I0inv%*%Ii%*%I0inv%*%Ij%*%I0inv
        mat2 <- I0inv%*%Ij%*%I0inv
        mat3 <- I0inv%*%Ii%*%I0inv
        tmp_mat <- sigma_inv%*%mat1 - sigma_inv%*%mat2%*%sigma_inv%*%mat3
        second_derivative[i,j] <- trace(tmp_mat)
        
      }
      
    }
    
    second_deriv_inv <- Ginv(second_derivative)
    new_weights[1:(m-1)] <- old_weights[1:(m-1)] - alpha*t(second_deriv_inv%*%first_derivative)
    
    new_weights[m] <- 1-sum(new_weights[1:(m-1)])
    if( min(new_weights[1:(m-1)]) < 0 | (sum(new_weights[1:(m-1)]) > 1)){
      if(alpha > 0.00001){  alpha <- alpha/2  }
      else{ 
        indic <- FALSE
      }
    }
    else{
      diff <- sum((new_weights[1:(m-1)] - old_weights[1:(m-1)])^2)
      iter <- iter + 1
    }
    old_weights <- new_weights
    
    
    
    
    
    
    
  }
  
  return(new_weights)
  
}


get_logD_criterion <- function(candidates, weights, beta, family){
  X <- cbind(1, candidates)
  W <- diag(weights)
  
  eta <- as.vector(X%*%beta)
  if(family == "poisson"){
    vs <- exp(eta)
  }
  if(family == "binomial"){
    vs <- exp(eta)/(1+exp(eta))
  }
  V <- diag(vs)
  if(nrow(V) == 1 || nrow(W) == 1){
    M <- t(X)%*%X*W[1,1]*V[1,1]
  } else{
    M <- t(X)%*%V%*%W%*%X
  }
  return(log(det(M)))
}


# searches for a maximin design using a pre-defined grid over betas
get_maximin_design <- function(initial_inds, betagrid, candidates, family,
                               max_iter = 100){
  N <- nrow(candidates)
  n_beta <- nrow(betagrid)
  p <- ncol(betagrid)
  weights <- matrix(nrow = n_beta, ncol = N)
  opt_design_inds <- list()
  log_deffs <- numeric(n_beta)
  for(i in 1:n_beta){
    print(i)
    tryCatch({
      out1 <- get_opt_design(initial_inds, candidates, betagrid[i,], family, max_iter = max_iter)
      locally_opt_inds <- out1[[1]]
      opt_design_inds[[i]] <- locally_opt_inds
      locally_opt_weights <- numeric(N)
      locally_opt_weights[locally_opt_inds] <- out1[[2]]
      weights[i,] <- locally_opt_weights
    }, error=function(e){cat("Local Design Not Found for Maximin:",conditionMessage(e), "\n")})
    
  }
  
  worst_designs <- matrix(nrow = n_beta, ncol = N)
  worst_Deff_diffs <- numeric(n_beta)
  for(i in 1:n_beta){
    
    log_deffii <- get_logD_criterion(candidates, weights[i,],
                                     betagrid[i,],family=family)
    
    # for each set of weights, what is the worst beta? (lowest relative D-efficiency)
    Deff_diffs <- numeric(n_beta)
    for(j in 1:n_beta){
      
      log_deffij <- get_logD_criterion(candidates, weights[j,],
                                       betagrid[i,],family=family)
      Deff_diffs[j] <- log_deffij - log_deffii
      # print(log_deffij)
      # if(is.na(log_deffij) || is.na(log_deffii)){
      #   
      # } else{
      #   Deff_diffs[j] <- log_deffij - log_deffii
      # }
      
    }
    worst_case <- which.min(Deff_diffs)
    if(length(worst_case) == 0){
      worst_Deff_diffs[i] <- -Inf # don't choose this design
      
    } else{
      worst_Deff_diffs[i] <- Deff_diffs[worst_case]
      worst_designs[i,] <- weights[worst_case,]
    }
    
    
  }
  
  # choose the best of the worst designs
  best_worse_case <- which.max(worst_Deff_diffs)
  best_worse_design <- worst_designs[best_worse_case,]
  
  
  return(best_worse_design)
  
}


# this function finds optimal designs for each value in betagrid
# and then uses kmeans clustering to compute the final design
# k is the number of support points for the final cluster design
get_cluster_design <- function(initial_inds, betagrid, candidates, family, k = 5,
                               max_iter = 100){
  
  N <- nrow(candidates)
  n_beta <- nrow(betagrid)
  p <- ncol(betagrid)
  ncol_cands <- ncol(candidates)
  weights <- matrix(nrow = n_beta, ncol = N)
  design_pts <- NULL
  opt_design_inds <- list()
  log_deffs <- numeric(n_beta)
  for(i in 1:n_beta){
    print(i)
    tryCatch({
      out1 <- get_opt_design(initial_inds, candidates, betagrid[i,], family, max_iter = max_iter)
      locally_opt_inds <- out1[[1]]
      opt_design_inds[[i]] <- locally_opt_inds
      locally_opt_weights <- numeric(N)
      locally_opt_weights[locally_opt_inds] <- out1[[2]]
      weights[i,] <- locally_opt_weights
      for(j in 1:N){
        if(ceiling(N*locally_opt_weights[j]) > 0){
          nreps <- ceiling(N*locally_opt_weights[j])
          new_row <- candidates[j,]
          new_pts <- matrix(nrow = nreps, ncol = ncol_cands, rep(new_row, nreps), byrow = TRUE )
          design_pts <- rbind(design_pts, new_pts)
        }
      }
    }, error=function(e){cat("Local Design Not Found:",conditionMessage(e), "\n")})
    
  }
  
  n_unique = nrow(unique(design_pts))
  if(k > n_unique){
    cluster.out1 <- kmeans(design_pts, centers = n_unique)
    final_weights <- rep(1/n_unique, n_unique)
  } else{
    cluster.out1 <- kmeans(design_pts, centers = k)
    final_weights <- rep(1/k, k)
  }
  final_design_pts <- cluster.out1$centers
  # final_weights <- rep(1/k, k)
  
  # also return the design from the BAG_cluster method
  # mean_weights <- colMeans(weights, na.rm = TRUE)
  # weighted_candidates <- candidates*mean_weights
  # cluster.out2 <- kmeans(weighted_candidates, centers = k)
  # final_weights2 <- rep(1/k,k)
  # final_design_pts2 <- cluster.out2$centers
  
  return(list(final_design_pts, final_weights))
  
}




######################################
###### Simulation Parameters #########
######################################
set.seed(1234) # set random seed
n <- 50 # pilot data size 
# n <- 200
sigma <- 0.5 # standard deviation of pilot data
# sigma <- 2
p <- 2  # number of predictors
B <- 50 # number of bootstrap resamples
# family = "binomial" 
family = "poisson"
n_sim <- 30 # number of pilot datasets
n_iter <- 25 # number of times designs are found for each pilot dataset
n_beta_eval <- 100 # number of betas used for efficiency evaluation
# robust_scenario <- 1 # this is 1 or 2
robust_scenario <- 2

n.cores <- 30
my.cluster <- parallel::makeCluster(n.cores)
registerDoParallel(my.cluster)

# generate candidate design points
xlist <- list()
for(i in 1:p){
  xlist[[i]] <- seq(from = -1, to = 1, by = 0.1)
  
}

candidates <- expand.grid(xlist)
xnames <- paste("x",1:p, sep = "")
colnames(candidates) <- xnames
candidates <- as.matrix(candidates)
N = nrow(candidates)

if(robust_scenario == 1){
    eval_beta_grid <- matrix(rnorm(n_beta_eval*(p+1), mean = c(-3,-4,-5),sd = 0.5), nrow = n_beta_eval , ncol = p+1)
}
if(robust_scenario == 2){
    eval_beta_grid <- matrix(rnorm(n_beta_eval*(p+1), mean = c(1,3,1),sd = 0.5), nrow = n_beta_eval , ncol = p+1)
}




simout <- foreach(ind = 1:n_sim, .packages = c("matlib","pracma"), .combine = 'rbind', .errorhandling = "remove" ) %dopar% {

    # generate true beta
    true_beta <- runif(p+1, min = 0, max = 1)
    init_inds <- round(seq(from = 1, to = nrow(candidates), length.out = length(true_beta)+1))


    # generate pilot data of size n
    keep_sampling <- TRUE

    while(keep_sampling){
        X_pilot <- matrix(nrow = n, ncol = p)
        Y <- numeric(n)
        for(ell in 1:p){
            X_pilot[,ell] <- rnorm(n, mean = 0, sd = sigma)
        }
  
        eta <- cbind(1,X_pilot)%*%true_beta  # 1 is for the intercept
        if(family == "binomial"){
                Y <- rbinom(n, size = 1, prob = exp(eta)/(1+exp(eta)))
        }
        else{
                Y <- rpois(n, lambda = exp(eta))
        }
    pilot_data <- cbind(X_pilot, Y)
    xnames <- paste("x",1:p, sep = "")
    colnames(pilot_data) <- c(xnames,"y")
    pilot_data <- as.data.frame(pilot_data)
    formula <- paste("y ~ ", paste(xnames,collapse="+"), sep = "")
    local_model <-  glm(formula, data = pilot_data, family = family)
    if(local_model$converged){

        # get locally optimal design
        # if there is no locally optimal design, resample
        betahat <- local_model$coefficients
        local_design <- get_opt_design(initial_inds = init_inds, candidates = candidates, beta = betahat, family = family)
        locally_opt_inds <- local_design[[1]]
        locally_opt_weights <- numeric(N)
        locally_opt_weights[locally_opt_inds] <- local_design[[2]]
        n_local <- length(locally_opt_inds)
    
        if(!is.na(locally_opt_weights[1])){
            keep_sampling <- FALSE
        }

        }
  
    }


    avg_local_Ds <- numeric(n_iter)
    avg_BAG_Ds <- numeric(n_iter)
    avg_BAG_trimmed0.05Ds <- numeric(n_iter)
    avg_BAG_trimmed0.025Ds <- numeric(n_iter)
    avg_BAG_geo_medians <- numeric(n_iter)
    avg_BAG_geo_median0.05Ds <- numeric(n_iter)
    avg_BAG_geo_median0.025Ds <- numeric(n_iter)
    avg_cluster_Ds <- numeric(n_iter)
    avg_maximin_Ds <- numeric(n_iter)


    n_BAGs <- numeric(n_iter)
    n_BAG_trimmed0.05Ds <- numeric(n_iter)
    n_BAG_trimmed0.025Ds <- numeric(n_iter)
    n_BAG_geo_medians <- numeric(n_iter)
    n_BAG_geo_median0.05Ds <- numeric(n_iter)
    n_BAG_geo_median0.025Ds <- numeric(n_iter)
    n_clusters <- numeric(n_iter)
    n_maximins <- numeric(n_iter)


    relD_BAG_true <- numeric(n_iter)
    relD_BAG0.05_true <- numeric(n_iter)
    relD_BAG0.025_true <- numeric(n_iter)
    relD_BAGgeo_true <- numeric(n_iter)
    relD_BAGgeo0.05_true <- numeric(n_iter)
    relD_BAGgeo0.025_true <- numeric(n_iter)
    relD_cluster_true <- numeric(n_iter)
    relD_maximin_true <- numeric(n_iter)


    for(iter in 1:n_iter){
  
        # get BAG design (and candidates)
        out1 <- get_bag_design_weights(pilot_data, candidates, init_inds, formula, family = family, B)
        BAG_weights <- out1[[1]]
        betagrid <- out1[[2]]
        BAG_weights_trimmed0.05 <- out1[[3]]
        BAG_weights_trimmed0.025 <- out1[[4]]
        BAG_weights_geo_median <- out1[[5]]
        BAG_weights_geo_median0.05 <- out1[[6]]
        BAG_weights_geo_median0.025 <- out1[[7]]
        n_BAGs[iter] <- sum(BAG_weights > 0)
        n_BAG_trimmed0.05Ds[iter] <- sum(BAG_weights_trimmed0.05 > 0)
        n_BAG_trimmed0.025Ds[iter] <- sum(BAG_weights_trimmed0.025 > 0)
        n_BAG_geo_medians[iter] <- sum(BAG_weights_geo_median > 0)
        n_BAG_geo_median0.05Ds[iter] <- sum(BAG_weights_geo_median0.05 > 0)
        n_BAG_geo_median0.025Ds[iter] <- sum(BAG_weights_geo_median0.025 > 0)
  
        # get cluster design
        out_cluster <- get_cluster_design(init_inds, betagrid, candidates, family = family, k = 5)
        n_clusters[iter] <- sum(out_cluster[[2]] > 0)
  
        # get maximin design
        maximin_weights <- get_maximin_design(init_inds, betagrid, candidates, family = family, max_iter = 100)
        n_maximins[iter] <- sum(maximin_weights > 0)
  
  
        relD_BAG <- numeric(n_beta_eval)
        relD_BAGtrimmed0.05 <- numeric(n_beta_eval)
        relD_BAGtrimmed0.025 <- numeric(n_beta_eval)
        relD_BAGgeo_median <- numeric(n_beta_eval)
        relD_BAGgeo_median0.05 <- numeric(n_beta_eval)
        relD_BAGgeo_median0.025 <- numeric(n_beta_eval)
        relD_cluster <- numeric(n_beta_eval)
        relD_maximin <- numeric(n_beta_eval)
        logD_local <- numeric(n_beta_eval)
  
  
        # store the relative D-efficiencies for each design to the local design
        for(i in 1:n_beta_eval){
            logD_local[i] <- get_logD_criterion(candidates, locally_opt_weights, beta = eval_beta_grid[i,], family = family)
            relD_BAG[i] <- exp(get_logD_criterion(candidates, BAG_weights, beta = eval_beta_grid[i,], family = family)- logD_local[i])
            relD_BAGtrimmed0.05[i] <- exp(get_logD_criterion(candidates, BAG_weights_trimmed0.05, beta = eval_beta_grid[i,], family = family) - logD_local[i])
            relD_BAGtrimmed0.025[i] <- exp(get_logD_criterion(candidates, BAG_weights_trimmed0.025, beta = eval_beta_grid[i,], family = family) - logD_local[i])
            relD_BAGgeo_median[i] <- exp(get_logD_criterion(candidates, BAG_weights_geo_median, beta = eval_beta_grid[i,], family = family) - logD_local[i])
            relD_BAGgeo_median0.05[i] <- exp(get_logD_criterion(candidates, BAG_weights_geo_median0.05, beta = eval_beta_grid[i,], family = family) - logD_local[i])
            relD_BAGgeo_median0.025[i] <- exp(get_logD_criterion(candidates, BAG_weights_geo_median0.025, beta = eval_beta_grid[i,], family = family) - logD_local[i])
            relD_cluster[i] <- exp(get_logD_criterion(out_cluster[[1]], out_cluster[[2]], beta = eval_beta_grid[i,], family = family) - logD_local[i])
            relD_maximin[i] <- exp(get_logD_criterion(candidates, maximin_weights, beta = eval_beta_grid[i,], family = family) - logD_local[i])
        }
        
        # set infinite values to NA so that they will be removed with na.rm = TRUE
        relD_cluster[is.infinite(relD_cluster)] <- NA
        relD_maximin[is.infinite(relD_maximin)] <- NA
        relD_BAG[is.infinite(relD_BAG)] <- NA
        relD_BAGtrimmed0.05[is.infinite(relD_BAGtrimmed0.05)] <- NA
        relD_BAGtrimmed0.025[is.infinite(relD_BAGtrimmed0.025)] <- NA
        relD_BAGgeo_median[is.infinite(relD_BAGgeo_median)] <- NA
        relD_BAGgeo_median0.05[is.infinite(relD_BAGgeo_median0.05)] <- NA
        relD_BAGgeo_median0.025[is.infinite(relD_BAGgeo_median0.025)] <- NA
  
        avg_cluster_Ds[iter] <- mean(relD_cluster, na.rm = TRUE)
        avg_maximin_Ds[iter] <- mean(relD_maximin, na.rm = TRUE)
        avg_BAG_Ds[iter] <- mean(relD_BAG, na.rm = TRUE)
        avg_BAG_trimmed0.05Ds[iter] <- mean(relD_BAGtrimmed0.05, na.rm = TRUE)
        avg_BAG_trimmed0.025Ds[iter] <- mean(relD_BAGtrimmed0.025, na.rm = TRUE)
        avg_BAG_geo_medians[iter] <- mean(relD_BAGgeo_median, na.rm = TRUE)
        avg_BAG_geo_median0.05Ds[iter] <- mean(relD_BAGgeo_median0.05, na.rm = TRUE)
        avg_BAG_geo_median0.025Ds[iter] <- mean(relD_BAGgeo_median0.025, na.rm = TRUE)
  
        }


    avg_BAG_D <- mean(avg_BAG_Ds, na.rm = TRUE)
    avg_BAG_trimmed0.05_D <- mean(avg_BAG_trimmed0.05Ds, na.rm = TRUE)
    avg_BAG_trimmed0.025_D <- mean(avg_BAG_trimmed0.025Ds, na.rm = TRUE)
    avg_BAG_geo_median_D <- mean(avg_BAG_geo_medians, na.rm = TRUE)
    avg_BAG_geo_median0.05_D <- mean(avg_BAG_geo_median0.05Ds, na.rm = TRUE)
    avg_BAG_geo_median0.025_D <- mean(avg_BAG_geo_median0.025Ds, na.rm = TRUE)
    avg_cluster_D <- mean(avg_cluster_Ds, na.rm = TRUE)
    avg_maximin_D <- mean(avg_maximin_Ds, na.rm = TRUE)
    avg_n_BAG <- mean(n_BAGs)
    avg_n_BAG_trimmed0.05 <- mean(n_BAG_trimmed0.05Ds)
    avg_n_BAG_trimmed0.025 <- mean(n_BAG_trimmed0.025Ds)
    avg_n_BAG_geo_median <- mean(n_BAG_geo_medians)
    avg_n_BAG_geo_median0.05 <- mean(n_BAG_geo_median0.05Ds)
    avg_n_BAG_geo_median0.025 <- mean(n_BAG_geo_median0.025Ds)
    avg_n_cluster <- mean(n_clusters)
    avg_n_maximin <- mean(n_maximins)
    c(avg_BAG_D, avg_BAG_trimmed0.05_D, avg_BAG_trimmed0.025_D, avg_BAG_geo_median_D, avg_BAG_geo_median0.05_D, avg_BAG_geo_median0.025_D,
            avg_cluster_D, avg_maximin_D, n_local, avg_n_BAG, avg_n_BAG_trimmed0.05,
            avg_n_BAG_trimmed0.025, avg_n_BAG_geo_median, avg_n_BAG_geo_median0.05, avg_n_BAG_geo_median0.025, avg_n_cluster, avg_n_maximin)
}


namestr <- paste("RobustBAGsim_withmedian_n", n, "sigma", sigma, family, "p", p, "robust",robust_scenario, ".csv", sep = "")
# colnames(simout) = c("local_D", "BAG_D", "cluster_D", "maximin_D", "n_local", "n_BAG", "n_cluster", "n_maximin")
colnames(simout) <- c("BAG_D", "BAG_trimmed0.02_D", "BAG_trimmed0.01_D", "BAG_geo_median_D", "BAG_geo_median0.02_D", "BAG_geo_median0.01D", "cluster_D", "maximin_D", "n_local", 
                   "n_BAG", "n_BAG_trimmed0.02", "n_BAG_trimmed0.01", "n_BAG_geo_median", "n_BAG_geo_median0.02","n_BAG_geo_median0.01", "n_cluster", "n_maximin")
write.csv(simout, namestr)



parallel::stopCluster(my.cluster)







