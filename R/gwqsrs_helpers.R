

# function to select variables
rs.select_vars <- function(data, na.action, formula, mix_name, ...){
  allvars = all.vars(formula)
  other_vars <- c(...)
  data$wqs = 0
  data <- data[, c(allvars, mix_name, other_vars)]
  if(missing(na.action)) na.action <- na.omit
  dtf <- gWQS:::na_action(data, na.action)
  return(dtf)
}


# function to sample from q_name
rs.sample_f = function(i, q_name, n_vars){

  if(!is.numeric(n_vars) & !is.null(n_vars)) stop("n_vars must be either numeric or NULL")
  if(is.numeric(n_vars) & n_vars >= length(q_name)) stop("n_vars must be lower than the number of elements in mix_name")
  slctd_vars = sample(q_name, n_vars, replace=FALSE)

  return(slctd_vars)
}

# function that call optim.f to estimate parameters for each bootstrap sample
rs.estimate_param = function(slctd_vars, dtf, Q, b1_pos, b1_constr, family, zilink, zero_infl, formula, ff, weights, stratified, control){

  # param = gWQS:::optim.f(dtf, Q[,slctd_vars], b1_pos, b1_constr, family, zilink, zero_infl, formula, ff, weights, stratified, control, lp=0, ln=0)
  param = gWQS:::optim.f(dtf, Q[,slctd_vars], b1_pos, b1_constr, family, zilink, zero_infl, formula, ff, weights, control, lp=0, ln=0)

  return(param)
}


# function to be passed to future_lapply to fit the models
rs.model.fit_f = function(i, param, slctd_vars, dtf, Q, b1_pos, family, zilink, formula, ff, weights, stratified, zero_infl){

  # b_fit <- gWQS:::model.fit(param[[i]]$par_opt, dtf, Q[,slctd_vars[[i]]], family, zilink, formula, ff, weights, stratified, b1_pos, zero_infl)
  b_fit <- gWQS:::model.fit(param[[i]]$par_opt, dtf, Q[,slctd_vars[[i]]], family, zilink, formula, ff, weights, b1_pos, zero_infl)

  return(b_fit)
}

set_par_names <- function(i, slctd_vars, param, q_name, family){

  temp <- param[[i]]$par_opt
  if(family$family == "multinomial"){
    param[[i]]$par_opt <- matrix(NA, length(q_name), dim(temp)[2])
    param[[i]]$par_opt[which(q_name %in% slctd_vars[[i]]),] <- temp
    rownames(param[[i]]$par_opt) <- q_name
  }
  else{
    param[[i]]$par_opt <- rep(NA, length(q_name))
    names(param[[i]]$par_opt) <- q_name
    param[[i]]$par_opt[slctd_vars[[i]]] <- temp
  }

  return(param[[i]])
}

# function that calls the optimization function and the function to fit the model for each bootstrap sample
rs.par.modl.est <- function(dtf, Q, formula, ff, weights, rs, n_vars, b1_pos, b1_constr, family, zilink, zero_infl, plan_strategy, seed, stratified, control){

  if (family$family %in% c("gaussian", "quasipoisson")) ts = "t"
  else if (family$family %in% c("binomial", "poisson", "multinomial", "negbin")) ts = "z"

  if(!is.numeric(rs)) stop("'rs' must be a number")
  if(is.null(seed)) seed = TRUE

  plan(plan_strategy)
  slctd_vars = future_lapply(X=1:rs, FUN = rs.sample_f, q_name = colnames(Q), n_vars = n_vars, future.seed = seed)

  plan(plan_strategy)
  param = lapply(X = slctd_vars, FUN = rs.estimate_param, dtf = dtf, Q = Q, b1_pos = b1_pos,
                        b1_constr = b1_constr, family = family, zilink = zilink, zero_infl = zero_infl,
                        formula = formula, ff = ff, weights = weights, stratified = stratified, control = control)
  # param = future_lapply(X = slctd_vars, FUN = rs.estimate_param, dtf = dtf, Q = Q, b1_pos = b1_pos,
  #                       b1_constr = b1_constr, family = family, zilink = zilink, zero_infl = zero_infl,
  #                       formula = formula, ff = ff, weights = weights, stratified = stratified, control = control, future.seed = FALSE)

  plan(plan_strategy)
  b_fit = future_lapply(X = 1:rs, FUN = rs.model.fit_f, param = param, slctd_vars = slctd_vars, dtf = dtf, Q = Q,
                        b1_pos = b1_pos, family = family, zilink = zilink, formula = formula, ff = ff,
                        weights = weights, stratified = stratified, zero_infl = zero_infl, future.seed = FALSE)

  plan(plan_strategy)
  # param <- lapply(X = 1:rs, FUN = set_par_names, slctd_vars, param, q_name = colnames(Q), family = family)
  param <- future_lapply(X = 1:rs, FUN = set_par_names, slctd_vars, param, q_name = colnames(Q), family = family,
                         future.seed = FALSE)

  conv <- c(sapply(param, function(i) i$conv))
  nfuneval <- c(sapply(param, function(i) i$nfuneval))
  if(family$family == "multinomial"){
    n_levels <- dim(param[[1]]$par_opt)[2]+1
    wqs_site <- which(grepl("^wqs_", rownames(b_fit[[1]]$m_f$sum_stat)))
    wght_matrix <- lapply(1:(n_levels-1), function(j) do.call("rbind", lapply(param, function(i) i$par_opt[,j])))
    b1 <- lapply(wqs_site, function(j) sapply(b_fit, function(i) i$m_f$nlm_out$estimate[j]))
    se <- lapply(wqs_site, function(j) sapply(b_fit, function(i) i$m_f$sum_stat$Standard_Error[j]))
    stat <- lapply(wqs_site, function(j) sapply(b_fit, function(i) i$m_f$sum_stat$stat[j]))
    p_val <- lapply(wqs_site, function(j) sapply(b_fit, function(i) i$m_f$sum_stat$p_value[j]))
  }
  else{
    wght_matrix <- do.call("rbind", lapply(param, function(i) i$par_opt))
    if(zero_infl){
      b1_count <- sapply(b_fit, function(i) i$m_f$coefficients$count["wqs"])
      se_count <- sapply(b_fit, function(i) summary(i$m_f)$coefficients$count["wqs", "Std. Error"])
      stat_count <- sapply(b_fit, function(i) summary(i$m_f)$coefficients$count["wqs", paste0(ts, " value")])
      p_val_count <- sapply(b_fit, function(i) summary(i$m_f)$coefficients$count["wqs", gsub("x", ts, "Pr(>|x|)")])
      if("wqs" %in% names(b_fit[[1]]$m_f$coefficients$zero)){
        b1_zero <- sapply(b_fit, function(i) i$m_f$coefficients$zero["wqs"])
        se_zero <- sapply(b_fit, function(i) summary(i$m_f)$coefficients$zero["wqs", "Std. Error"])
        stat_zero <- sapply(b_fit, function(i) summary(i$m_f)$coefficients$zero["wqs", paste0(ts, " value")])
        p_val_zero <- sapply(b_fit, function(i) summary(i$m_f)$coefficients$zero["wqs", gsub("x", ts, "Pr(>|x|)")])
      }
      else b1_zero <- se_zero <- stat_zero <- p_val_zero <- NULL
    }
    else{
      b1 <- sapply(b_fit, function(i) summary(i$m_f)$coefficients["wqs", "Estimate"])
      se <- sapply(b_fit, function(i) summary(i$m_f)$coefficients["wqs", "Std. Error"])
      # stat <- sapply(b_fit, function(i){
      #   chisq <- wald.test(coef(i$m_f), vcov(i$m_f), Terms = which(names(i$m_f$coefficients)=="wqs"))
      #   chisq$result$chi2[1]
      # })
      stat <- sapply(b_fit, function(i) summary(i$m_f)$coefficients["wqs", paste0(ts, " value")])
      p_val <- sapply(b_fit, function(i) summary(i$m_f)$coefficients["wqs", gsub("x", ts, "Pr(>|x|)")])
    }
    n_levels <- 1
  }

  n_non_conv = sum(conv == 2)
  if(n_non_conv == 0 & control$trace) message(paste0("The optimization function always converged\n"))
  else if(n_non_conv == rs) stop("The optimization function never converged\n")
  else if(control$trace) message(paste0("The optimization function did not converge ", n_non_conv, " time/times\n"))

  # estimate mean weight for each component (exclude weights from iterations with failed convergence)
  if (family$family == "multinomial"){
    bres <- Map(cbind, wght_matrix, b1, se, stat, p_val)
    bres <- lapply(bres, as.data.frame)
    bres <- lapply(bres, setNames, c(colnames(Q), "b1", "Std_Error", "stat", "p_val"))
    strata_names <- gsub("wqs_", "", rownames(b_fit[[1]]$m_f$sum_stat)[wqs_site])
    names(bres) <- strata_names
  }
  else {
    if(zero_infl){
      if(is.null(b1_zero)){
        bres <- as.data.frame(cbind(wght_matrix, b1_count, se_count, stat_count, p_val_count))
        names(bres) <- c(colnames(Q), "b1_count", "Std_Error_count", "stat_count", "p_val_count")
      }
      else{
        bres <- as.data.frame(cbind(wght_matrix, b1_count, se_count, stat_count, p_val_count, b1_zero, se_zero, stat_zero, p_val_zero))
        names(bres) <- c(colnames(Q), "b1_count", "Std_Error_count", "stat_count", "p_val_count", "b1_zero", "Std_Error_zero", "stat_zero", "p_val_zero")
      }
    }
    else{
      bres <- as.data.frame(cbind(wght_matrix, b1, se, stat, p_val))
      names(bres) <- c(colnames(Q), "b1", "Std_Error", "stat", "p_val")
    }
    strata_names <- NULL
  }

  par_model_out <- list(bres, conv, slctd_vars, nfuneval, n_levels, strata_names)
  names(par_model_out) <- c("bres", "conv", "slctd_vars", "nfuneval", "n_levels", "strata_names")

  return(par_model_out)
}


# function to estimate mean weights for each component
rs.mean_weight_f = function(mix_name, bres, conv, b1_pos, family, n_levels, strata_names, zero_infl){

  if (family$family == "multinomial"){
    mean_weight <- lapply(1:(n_levels-1), function(i){
      # bres[[i]][mix_name][is.na(bres[[i]][mix_name])] = 0
      if(b1_pos[i]) w_t = apply(bres[[i]][bres[[i]]$b1 > 0 & conv!=2, mix_name], 2, weighted.mean, abs(bres[[i]][bres[[i]]$b1 > 0 & conv!=2, "stat"]), na.rm = T)
      else if(!b1_pos[i]) w_t = apply(bres[[i]][bres[[i]]$b1 < 0 & conv!=2, mix_name], 2, weighted.mean, abs(bres[[i]][bres[[i]]$b1 < 0 & conv!=2, "stat"]), na.rm = T)
      if (all(is.nan(w_t)))
        stop(paste0("There are no ", ifelse(b1_pos[i], "positive", "negative"), " b1 in the bootstrapped models for ", strata_names[i]))
      if(any(is.na(w_t))) w_t[is.na(w_t)] <- 0
      w_t <- w_t/sum(w_t)
      return(w_t)
    })
    mean_weight <- list.cbind(mean_weight)
  }
  else{
    # bres[mix_name][is.na(bres[mix_name])] <- 0
    if(zero_infl){
      if(b1_pos) mean_weight = apply(bres[bres$b1_count > 0 & conv!=2, mix_name], 2, weighted.mean, abs(bres[bres$b1_count > 0 & conv!=2, "stat_count"]), na.rm = T)
      else mean_weight = apply(bres[bres$b1_count < 0 & conv!=2, mix_name], 2, weighted.mean, abs(bres[bres$b1_count < 0 & conv!=2, "stat_count"]), na.rm = T)
    }
    else{
      if(b1_pos) mean_weight = apply(bres[bres$b1 > 0 & conv!=2, mix_name], 2, weighted.mean, abs(bres[bres$b1 > 0 & conv!=2, "stat"]), na.rm = T)
      else mean_weight = apply(bres[bres$b1 < 0 & conv!=2, mix_name], 2, weighted.mean, abs(bres[bres$b1 < 0 & conv!=2, "stat"]), na.rm = T)
    }
    if(all(is.nan(mean_weight)))
      stop("There are no ", ifelse(b1_pos, "positive", "negative"), " b1 in the bootstrapped models")
    if(any(is.na(mean_weight))) mean_weight[is.na(mean_weight)] <- 0
    mean_weight <- mean_weight/sum(mean_weight)
  }

  return(mean_weight)
}
