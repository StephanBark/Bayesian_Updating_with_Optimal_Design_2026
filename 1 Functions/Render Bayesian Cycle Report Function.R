# Report tools for fit_general_multiple_cycles objects -------------------------
#
# This script implements a report workflow for the object returned by
# fit_general_multiple_cycles(). The main user-facing entry
# point is render_fit_general_multiple_cycles_report().
#
# Top-level functions:
#   1) build_fit_general_multiple_cycles_report_spec()
#   2) render_bayesian_cycle_report_core()
#   3) render_fit_general_multiple_cycles_report()
#

build_fit_general_multiple_cycles_report_spec <- function(
    fit_object,
    reference_values = NULL,
    residual_trace_max = 12,
    matrix_entry_max = 15,
    pairs_max = 6,
    digits = 4
) {
  `%||%` <- function(x, y) {
    if (is.null(x) || (length(x) == 1 && is.na(x))) y else x
  }
  
  required_pkgs <- c("dplyr", "purrr", "rstan")
  missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(missing_pkgs) > 0) {
    stop(
      "The following packages are required but not installed: ",
      paste(missing_pkgs, collapse = ", "),
      call. = FALSE
    )
  }
  
  if (!inherits(fit_object, "fit_general_multiple_cycles")) {
    warning(
      "fit_object does not inherit from class 'fit_general_multiple_cycles'. ",
      "The report builder will still try to use the expected object structure."
    )
  }
  
  if (is.null(fit_object$cycles) || length(fit_object$cycles) == 0) {
    stop("fit_object$cycles is empty. Nothing to report.", call. = FALSE)
  }
  
  if (is.null(reference_values)) {
    reference_values <- numeric(0)
  }
  if (!is.numeric(reference_values)) {
    stop("reference_values must be a named numeric vector.", call. = FALSE)
  }
  if (length(reference_values) > 0 && is.null(names(reference_values))) {
    stop("reference_values must be named.", call. = FALSE)
  }
  
  sanitize_effect_name <- function(effect_name) {
    gsub(":", "_", effect_name)
  }
  
  resolve_digits <- function(d) {
    if (is.null(d) || length(d) == 0 || is.na(d[1]) || !is.finite(d[1])) return(4L)
    as.integer(d[1])
  }
  
  format_num <- function(x, digits = NULL) {
    digits <- resolve_digits(digits)
    if (length(x) == 0 || all(is.na(x))) return(NA_character_)
    if (inherits(x, "difftime")) x <- as.numeric(x)
    if (is.logical(x)) return(as.character(x))
    if (is.numeric(x)) {
      rounded <- round(x, digits)
      if (isTRUE(all.equal(rounded, as.integer(rounded), tolerance = 1e-12))) {
        return(format(as.integer(rounded), trim = TRUE, scientific = FALSE))
      }
      return(format(rounded, nsmall = digits, trim = TRUE, scientific = FALSE))
    }
    as.character(x)
  }
  
  flatten_named_object <- function(x, component_name = "") {
    rows <- list()
    
    rec <- function(obj, path) {
      if (is.null(obj)) {
        rows[[length(rows) + 1]] <<- data.frame(
          component = path,
          field = NA_character_,
          value = "NULL",
          stringsAsFactors = FALSE
        )
        return(invisible(NULL))
      }
      
      if (is.atomic(obj) && is.null(dim(obj)) && length(obj) <= 10) {
        rows[[length(rows) + 1]] <<- data.frame(
          component = path,
          field = NA_character_,
          value = paste(obj, collapse = ", "),
          stringsAsFactors = FALSE
        )
        return(invisible(NULL))
      }
      
      if (is.matrix(obj) || is.array(obj)) {
        rows[[length(rows) + 1]] <<- data.frame(
          component = path,
          field = paste0(class(obj)[1], " ", paste(dim(obj), collapse = "x")),
          value = paste(capture.output(print(round(obj, resolve_digits(digits)))), collapse = " "),
          stringsAsFactors = FALSE
        )
        return(invisible(NULL))
      }
      
      if (is.list(obj)) {
        nms <- names(obj)
        if (is.null(nms)) nms <- paste0("[[", seq_along(obj), "]]")
        if (length(obj) == 0) {
          rows[[length(rows) + 1]] <<- data.frame(
            component = path,
            field = NA_character_,
            value = "empty list",
            stringsAsFactors = FALSE
          )
          return(invisible(NULL))
        }
        for (i in seq_along(obj)) {
          next_path <- if (nzchar(path)) paste0(path, " / ", nms[i]) else nms[i]
          rec(obj[[i]], next_path)
        }
        return(invisible(NULL))
      }
      
      rows[[length(rows) + 1]] <<- data.frame(
        component = path,
        field = class(obj)[1],
        value = paste(capture.output(print(obj)), collapse = " "),
        stringsAsFactors = FALSE
      )
    }
    
    rec(x, component_name)
    dplyr::bind_rows(rows)
  }
  
  effect_table <- function(effect_list, effect_role = c("fixed", "random")) {
    effect_role <- match.arg(effect_role)
    if (is.null(effect_list) || length(effect_list) == 0) {
      return(data.frame(
        effect = character(0),
        role = character(0),
        type = character(0),
        variables = character(0),
        stringsAsFactors = FALSE
      ))
    }
    
    nms <- names(effect_list)
    if (is.null(nms)) {
      nms <- vapply(effect_list, function(x) {
        if (is.character(x) && length(x) == 1) x else paste(x$vars %||% character(0), collapse = ":")
      }, character(1))
    }
    
    dplyr::bind_rows(lapply(seq_along(effect_list), function(i) {
      eff <- effect_list[[i]]
      if (is.character(eff) && length(eff) == 1) {
        data.frame(
          effect = nms[i],
          role = effect_role,
          type = ifelse(identical(eff, "beta_0"), "intercept", "character"),
          variables = eff,
          stringsAsFactors = FALSE
        )
      } else {
        eff_vars <- eff$vars %||% character(0)
        eff_type <- eff$type %||% NA_character_
        if (length(eff_vars) == 1 && identical(eff_vars[[1]], "beta_0")) {
          eff_type <- "intercept"
        }
        data.frame(
          effect = nms[i],
          role = effect_role,
          type = eff_type,
          variables = paste(eff_vars, collapse = ", "),
          stringsAsFactors = FALSE
        )
      }
    }))
  }
  
  cycle_label <- function(cycle_obj, i) {
    yrs <- cycle_obj$data_summary$years_included %||% NULL
    if (is.null(yrs) || length(yrs) == 0) {
      paste0("Window ", i)
    } else {
      paste0("Window ", i)
    }
  }
  
  extract_upper_tri_names <- function(arr_name, arr, max_entries = matrix_entry_max) {
    if (is.null(arr) || is.null(dim(arr)) || length(dim(arr)) != 3) return(character(0))
    p <- dim(arr)[2]
    idx <- which(upper.tri(matrix(TRUE, p, p), diag = TRUE), arr.ind = TRUE)
    if (nrow(idx) > max_entries) idx <- idx[seq_len(max_entries), , drop = FALSE]
    apply(idx, 1, function(z) sprintf("%s[%d,%d]", arr_name, z[1], z[2]))
  }
  
  extract_vector_names <- function(param_name, mat, max_entries = residual_trace_max) {
    if (is.null(mat) || is.null(dim(mat)) || length(dim(mat)) != 2) return(character(0))
    keep <- seq_len(min(ncol(mat), max_entries))
    sprintf("%s[%d]", param_name, keep)
  }
  
  get_call_value <- function(call_obj, name) {
    if (is.null(call_obj)) return(NA_character_)
    call_list <- as.list(call_obj)
    if (!(name %in% names(call_list))) return(NA_character_)
    paste(deparse(call_list[[name]]), collapse = "")
  }
  
  collect_model_variables <- function(fit_object) {
    vars <- character(0)
    
    fixed_effects <- fit_object$model_spec$fixed_effects %||% list()
    random_effects <- fit_object$model_spec$random_effects %||% list()
    
    for (eff in fixed_effects) {
      if (is.list(eff)) {
        vars <- c(vars, eff$vars %||% character(0))
      } else if (is.character(eff) && length(eff) == 1 && !identical(eff, "beta_0")) {
        vars <- c(vars, eff)
      }
    }
    
    for (eff in random_effects) {
      if (is.list(eff)) {
        vars <- c(vars, eff$vars %||% character(0))
      } else if (is.character(eff)) {
        vars <- c(vars, eff)
      }
    }
    
    vars <- unique(vars)
    vars <- setdiff(vars, c("yield", "beta_0"))
    unique(c(vars, "environment"))
  }
  
  extract_stan_data_list <- function(cycle_obj) {
    out <- tryCatch(cycle_obj$stan_fit@stan_args[[1]]$data, error = function(e) NULL)
    if (is.list(out)) return(out)
    NULL
  }
  
  extract_parameter_shape <- function(cycle_obj, param_name) {
    pm <- cycle_obj$posterior_means %||% list()
    ps <- cycle_obj$posterior_samples %||% list()
    
    if (param_name %in% names(pm)) {
      x <- pm[[param_name]]
      if (!is.null(dim(x))) return(as.integer(dim(x)))
      return(as.integer(length(x)))
    }
    
    if (param_name %in% names(ps)) {
      x <- ps[[param_name]]
      d <- dim(x)
      if (is.null(d)) return(as.integer(length(x)))
      if (length(d) <= 1) return(as.integer(d))
      return(as.integer(d[-1]))
    }
    
    NULL
  }
  
  resolve_cycle_variable_levels <- function(cycle_obj, fit_object, parsed_priors) {
    vars_to_report <- collect_model_variables(fit_object)
    known_levels <- stats::setNames(as.list(rep(NA_integer_, length(vars_to_report))), vars_to_report)
    
    stan_data_list <- extract_stan_data_list(cycle_obj)
    if (is.list(stan_data_list)) {
      for (var in vars_to_report) {
        n_key <- paste0("n_", var, "s")
        if (n_key %in% names(stan_data_list)) {
          known_levels[[var]] <- as.integer(stan_data_list[[n_key]][1])
        }
      }
    }
    
    if ("environment" %in% names(known_levels) && is.na(known_levels[["environment"]])) {
      known_levels[["environment"]] <- as.integer(cycle_obj$data_summary$n_environments %||% NA_integer_)
    }
    
    if ("year" %in% names(known_levels) && is.na(known_levels[["year"]])) {
      yrs <- cycle_obj$data_summary$years_included %||% NULL
      if (!is.null(yrs)) known_levels[["year"]] <- as.integer(length(unique(yrs)))
    }
    
    fixed_effects <- fit_object$model_spec$fixed_effects %||% list()
    for (nm in names(fixed_effects)) {
      eff <- fixed_effects[[nm]]
      if (!is.list(eff)) next
      eff_vars <- eff$vars %||% character(0)
      if (length(eff_vars) == 1 && !identical(eff_vars[[1]], "beta_0")) {
        beta_name <- paste0("beta_", eff_vars[[1]])
        shape <- extract_parameter_shape(cycle_obj, beta_name)
        if (!is.null(shape) && length(shape) == 1 && is.na(known_levels[[eff_vars[[1]]]])) {
          known_levels[[eff_vars[[1]]]] <- as.integer(shape[1] + 1L)
        }
      }
    }
    
    random_effects <- fit_object$model_spec$random_effects %||% list()
    for (nm in names(random_effects)) {
      eff <- random_effects[[nm]]
      if (!is.list(eff)) next
      eff_vars <- eff$vars %||% character(0)
      u_name <- paste0("u_", sanitize_effect_name(nm))
      shape <- extract_parameter_shape(cycle_obj, u_name)
      if (is.null(shape)) next
      
      if (identical(eff$type %||% NA_character_, "simple") && length(eff_vars) == 1) {
        if (is.na(known_levels[[eff_vars[[1]]]])) {
          known_levels[[eff_vars[[1]]]] <- as.integer(shape[1])
        }
      } else if (length(eff_vars) == 2 && length(shape) == 2) {
        if (is.na(known_levels[[eff_vars[[1]]]])) known_levels[[eff_vars[[1]]]] <- as.integer(shape[1])
        if (is.na(known_levels[[eff_vars[[2]]]])) known_levels[[eff_vars[[2]]]] <- as.integer(shape[2])
      }
    }
    
    changed <- TRUE
    while (isTRUE(changed)) {
      changed <- FALSE
      for (nm in names(random_effects)) {
        eff <- random_effects[[nm]]
        if (!is.list(eff)) next
        eff_vars <- eff$vars %||% character(0)
        if (length(eff_vars) < 3) next
        
        u_name <- paste0("u_", sanitize_effect_name(nm))
        shape <- extract_parameter_shape(cycle_obj, u_name)
        if (is.null(shape) || length(shape) != 2) next
        
        other_vars <- eff_vars[-length(eff_vars)]
        last_var <- eff_vars[length(eff_vars)]
        row_product <- as.integer(shape[1])
        col_size <- as.integer(shape[2])
        
        if (is.na(known_levels[[last_var]])) {
          known_levels[[last_var]] <- col_size
          changed <- TRUE
        }
        
        other_known <- vapply(other_vars, function(v) known_levels[[v]] %||% NA_integer_, integer(1))
        if (sum(is.na(other_known)) == 1) {
          missing_var <- other_vars[is.na(other_known)]
          denom <- prod(other_known[!is.na(other_known)])
          if (is.finite(denom) && denom > 0 && row_product %% denom == 0) {
            solved <- as.integer(row_product / denom)
            if (is.na(known_levels[[missing_var]]) && solved > 0) {
              known_levels[[missing_var]] <- solved
              changed <- TRUE
            }
          }
        }
      }
    }
    
    known_levels
  }
  
  summary_lines <- tryCatch(capture.output(summary(fit_object)), error = function(e) paste("summary() failed:", e$message))
  print_lines <- tryCatch(capture.output(print(fit_object)), error = function(e) paste("print() failed:", e$message))
  
  parsed_priors <- fit_object$parameters$initial_priors %||% list()
  first_post <- fit_object$cycles[[1]]$posterior_samples
  
  pretty_prior_type <- function(x) {
    x <- x %||% NA_character_
    x <- as.character(x)[1]
    gsub("_", " ", x)
  }
  
  format_matrix_text <- function(x) {
    if (is.null(x)) return("NA")
    paste(capture.output(print(round(x, resolve_digits(digits)))), collapse = "\n")
  }
  
  wrap_pre <- function(x) {
    paste0("<pre class='cell-pre'>", x, "</pre>")
  }
  
  build_prior_table <- function(prior_obj, prior_template) {
    if (is.null(prior_template) || length(prior_template) == 0) {
      return(data.frame(
        Component = character(0),
        Type = character(0),
        Value = character(0),
        stringsAsFactors = FALSE
      ))
    }
    
    rows <- lapply(names(prior_template), function(effect_name) {
      template <- prior_template[[effect_name]]
      prior_type <- template$type %||% NA_character_
      type_label <- pretty_prior_type(prior_type)
      stan_effect_name <- sanitize_effect_name(effect_name)
      
      value_text <- "NA"
      
      structured_effect <- (
        is.list(prior_obj) &&
          effect_name %in% names(prior_obj) &&
          is.list(prior_obj[[effect_name]]) &&
          !is.null(prior_obj[[effect_name]]$type)
      )
      
      if (structured_effect) {
        obj <- prior_obj[[effect_name]]
        
        if (identical(prior_type, "inverse_gamma")) {
          value_text <- paste(
            c(
              paste0("alpha = ", format_num(obj$alpha %||% NA_real_, digits)),
              paste0("beta = ", format_num(obj$beta %||% NA_real_, digits))
            ),
            collapse = "\n"
          )
        } else if (identical(prior_type, "inverse_wishart")) {
          value_lines <- c(
            if (!is.null(obj$blocks)) paste0("blocks = ", paste(obj$blocks, collapse = ", ")) else NULL,
            paste0("nu = ", format_num(obj$nu %||% NA_real_, digits)),
            "psi =",
            format_matrix_text(obj$psi %||% NULL)
          )
          value_text <- paste(value_lines, collapse = "\n")
        }
      } else {
        if (identical(prior_type, "inverse_gamma")) {
          if (identical(effect_name, "residual")) {
            alpha_val <- prior_obj[["alpha_var_resid_env_mean"]] %||% prior_obj[["alpha_var_resid_env"]] %||% NA_real_
            beta_val <- prior_obj[["beta_var_resid_env_mean"]] %||% prior_obj[["beta_var_resid_env"]] %||% NA_real_
          } else {
            alpha_val <- prior_obj[[paste0("alpha_var_", stan_effect_name)]] %||% NA_real_
            beta_val <- prior_obj[[paste0("beta_var_", stan_effect_name)]] %||% NA_real_
          }
          
          value_text <- paste(
            c(
              paste0("alpha = ", format_num(alpha_val, digits)),
              paste0("beta = ", format_num(beta_val, digits))
            ),
            collapse = "\n"
          )
        } else if (identical(prior_type, "inverse_wishart")) {
          nu_val <- prior_obj[[paste0("nu_", stan_effect_name)]] %||% NA_real_
          psi_val <- prior_obj[[paste0("psi_", stan_effect_name)]] %||% NULL
          
          value_lines <- c(
            if (!is.null(template$blocks)) paste0("blocks = ", paste(template$blocks, collapse = ", ")) else NULL,
            paste0("nu = ", format_num(nu_val, digits)),
            "psi =",
            format_matrix_text(psi_val)
          )
          value_text <- paste(value_lines, collapse = "\n")
        }
      }
      
      data.frame(
        Component = effect_name,
        Type = type_label,
        Value = wrap_pre(value_text),
        stringsAsFactors = FALSE
      )
    })
    
    dplyr::bind_rows(rows)
  }
  
  prior_names <- names(parsed_priors) %||% character(0)
  non_residual_names <- setdiff(prior_names, "residual")
  
  inverse_gamma_effects <- non_residual_names[
    vapply(parsed_priors[non_residual_names], function(x) identical(x$type %||% NA_character_, "inverse_gamma"), logical(1))
  ]
  inverse_wishart_effects <- non_residual_names[
    vapply(parsed_priors[non_residual_names], function(x) identical(x$type %||% NA_character_, "inverse_wishart"), logical(1))
  ]
  
  scalar_hist_params <- paste0("var_", vapply(inverse_gamma_effects, sanitize_effect_name, character(1)))
  scalar_hist_params <- scalar_hist_params[scalar_hist_params %in% names(first_post)]
  scalar_hist_params <- unique(scalar_hist_params)
  
  scalar_trace_params <- scalar_hist_params[scalar_hist_params != "var_resid_env_mean"]
  
  covariance_blocks <- stats::setNames(
    lapply(inverse_wishart_effects, function(effect_name) {
      arr_name <- paste0("Sigma_", sanitize_effect_name(effect_name))
      arr <- first_post[[arr_name]]
      list(
        effect_name = effect_name,
        param_name = arr_name,
        entries = extract_upper_tri_names(arr_name, arr, max_entries = matrix_entry_max)
      )
    }),
    inverse_wishart_effects
  )
  covariance_blocks <- covariance_blocks[vapply(covariance_blocks, function(x) length(x$entries) > 0, logical(1))]
  
  residual_trace_source <- NULL
  residual_trace_params <- character(0)
  if ("var_resid_env" %in% names(first_post) && is.matrix(first_post$var_resid_env)) {
    residual_trace_source <- "var_resid_env"
    residual_trace_params <- extract_vector_names("var_resid_env", first_post$var_resid_env, max_entries = residual_trace_max)
  } else if ("sigma_resid_env" %in% names(first_post) && is.matrix(first_post$sigma_resid_env)) {
    residual_trace_source <- "sigma_resid_env"
    residual_trace_params <- extract_vector_names("sigma_resid_env", first_post$sigma_resid_env, max_entries = residual_trace_max)
  }
  
  is_diagonal_entry <- function(entry_name) {
    idx <- sub("^.*\\[([0-9]+),([0-9]+)\\]$", "\\1,\\2", entry_name)
    parts <- strsplit(idx, ",", fixed = TRUE)[[1]]
    length(parts) == 2 && identical(parts[1], parts[2])
  }
  
  pair_params <- unique(c(
    scalar_trace_params,
    unlist(lapply(covariance_blocks, function(x) {
      diag_entries <- x$entries[vapply(x$entries, is_diagonal_entry, logical(1))]
      if (length(diag_entries) == 0) {
        diag_entries <- x$entries
      }
      diag_entries
    }), use.names = FALSE)
  ))
  if (length(pair_params) > pairs_max) pair_params <- pair_params[seq_len(pairs_max)]
  
  call_obj <- fit_object$session_info$call %||% NULL
  
  global_settings <- data.frame(
    setting = c("chains", "thin", "seed", "init", "n_cycles_fitted", "total_observations"),
    value = c(
      get_call_value(call_obj, "chains"),
      get_call_value(call_obj, "thin"),
      fit_object$session_info$seed %||% get_call_value(call_obj, "seed"),
      get_call_value(call_obj, "init"),
      fit_object$model_spec$n_cycles_fitted %||% length(fit_object$cycles),
      fit_object$model_spec$total_observations %||% sum(vapply(fit_object$cycles, function(x) x$data_summary$n_obs %||% 0, numeric(1)))
    ),
    stringsAsFactors = FALSE
  )
  
  cycle_mcmc_settings <- fit_object$model_spec$mcmc_settings %||% data.frame(message = "MCMC settings not available.")
  if (is.data.frame(cycle_mcmc_settings)) {
    cycle_mcmc_settings[] <- lapply(cycle_mcmc_settings, as.character)
  }
  
  model_effects_df <- dplyr::bind_rows(
    effect_table(fit_object$model_spec$fixed_effects %||% list(), "fixed"),
    effect_table(fit_object$model_spec$random_effects %||% list(), "random")
  )
  
  cycle_data_summary_df <- dplyr::bind_rows(lapply(seq_along(fit_object$cycles), function(i) {
    cyc <- fit_object$cycles[[i]]
    var_levels <- resolve_cycle_variable_levels(cyc, fit_object, parsed_priors)
    out <- data.frame(
      window = i,
      n_obs = cyc$data_summary$n_obs %||% NA_real_,
      stringsAsFactors = FALSE
    )
    if (length(var_levels) > 0) {
      for (nm in names(var_levels)) {
        out[[paste0("n_levels_", nm)]] <- as.integer(var_levels[[nm]])
      }
    }
    out
  }))
  
  cycle_timing_df <- dplyr::bind_rows(lapply(seq_along(fit_object$cycles), function(i) {
    cyc <- fit_object$cycles[[i]]
    cycle_time_raw <- cyc$cycle_time %||% NA_real_
    if (inherits(cycle_time_raw, "difftime")) cycle_time_raw <- as.numeric(cycle_time_raw, units = "mins")
    prior_time_raw <- cyc$post_parm_time %||% NA_real_
    if (inherits(prior_time_raw, "difftime")) prior_time_raw <- as.numeric(prior_time_raw, units = "mins")
    data.frame(
      window = i,
      sampling_time_min = cycle_time_raw,
      prior_update_time_sec = prior_time_raw * 60,
      stringsAsFactors = FALSE
    )
  }))
  
  total_computing_time_min <- sum(
    c(cycle_timing_df$sampling_time_min, cycle_timing_df$prior_update_time_sec / 60),
    na.rm = TRUE
  )
  
  cycle_summary_df <- dplyr::bind_rows(lapply(seq_along(fit_object$cycles), function(i) {
    cyc <- fit_object$cycles[[i]]
    conv <- fit_object$diagnostics$stan_convergence[[i]] %||% list()
    yrs <- cyc$data_summary$years_included %||% NULL
    data.frame(
      window = i,
      years = if (is.null(yrs)) NA_real_ else length(unique(yrs)),
      converged = conv$stan_converged %||% NA,
      rhat_max = format_num(conv$rhat_max %||% NA_real_),
      ess_min = format_num(conv$ess_min %||% NA_real_),
      divergent = conv$divergent_transitions %||% NA_real_,
      max_treedepth_hits = cyc$diagnostics$max_treedepth_hits %||% NA_real_,
      ebfmi_min = format_num(cyc$diagnostics$ebfmi_min %||% NA_real_),
      warnings = paste(c(conv$our_warnings, conv$stan_warnings), collapse = " | "),
      stringsAsFactors = FALSE
    )
  }))
  
  cycle_prior_tables <- stats::setNames(
    lapply(seq_along(fit_object$cycles), function(i) {
      build_prior_table(
        prior_obj = fit_object$cycles[[i]]$priors_used %||% list(),
        prior_template = parsed_priors
      )
    }),
    paste0("Window_", seq_along(fit_object$cycles))
  )
  
  final_prior_result <- fit_object$parameters$final_posteriors$result %||% fit_object$parameters$final_posteriors %||% list()
  final_posterior_update_df <- build_prior_table(
    prior_obj = final_prior_result,
    prior_template = parsed_priors
  )
  
  scalar_evolution_df <- dplyr::bind_rows(lapply(seq_along(fit_object$cycles), function(i) {
    cyc <- fit_object$cycles[[i]]
    dplyr::bind_rows(lapply(scalar_hist_params, function(param) {
      x <- cyc$posterior_samples[[param]]
      if (is.null(x)) return(NULL)
      q <- stats::quantile(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm = TRUE)
      data.frame(
        window = i,
        label = cycle_label(cyc, i),
        parameter = param,
        q2.5 = unname(q[1]),
        q25 = unname(q[2]),
        q50 = unname(q[3]),
        q75 = unname(q[4]),
        q97.5 = unname(q[5]),
        stringsAsFactors = FALSE
      )
    }))
  }))
  
  scalar_x_limits <- stats::setNames(vector("list", length(scalar_hist_params)), scalar_hist_params)
  for (param in scalar_hist_params) {
    pooled <- unlist(lapply(fit_object$cycles, function(cyc) {
      cyc$posterior_samples[[param]] %||% numeric(0)
    }), use.names = FALSE)
    pooled <- pooled[is.finite(pooled)]
    if (length(pooled) == 0) next
    lower <- if (all(pooled >= 0)) 0 else stats::quantile(pooled, 0.005, na.rm = TRUE)
    upper <- stats::quantile(pooled, 0.995, na.rm = TRUE)
    ref <- if (param %in% names(reference_values)) reference_values[[param]] else NA_real_
    if (!is.na(ref)) upper <- max(upper, ref, na.rm = TRUE)
    scalar_x_limits[[param]] <- c(lower, upper * 1.05)
  }
  
  covariance_x_limits <- list()
  for (block_name in names(covariance_blocks)) {
    entries <- covariance_blocks[[block_name]]$entries
    block_param <- covariance_blocks[[block_name]]$param_name
    covariance_x_limits[[block_name]] <- stats::setNames(vector("list", length(entries)), entries)
    for (entry in entries) {
      pooled <- unlist(lapply(fit_object$cycles, function(cyc) {
        arr <- cyc$posterior_samples[[block_param]]
        if (is.null(arr)) return(numeric(0))
        idx <- gsub(paste0("^", block_param, "\\[(.*)\\]$"), "\\1", entry)
        ij <- strsplit(idx, ",")[[1]]
        arr[, as.integer(ij[1]), as.integer(ij[2])]
      }), use.names = FALSE)
      pooled <- pooled[is.finite(pooled)]
      if (length(pooled) == 0) next
      lower <- stats::quantile(pooled, 0.005, na.rm = TRUE)
      upper <- stats::quantile(pooled, 0.995, na.rm = TRUE)
      ref <- if (entry %in% names(reference_values)) reference_values[[entry]] else NA_real_
      if (!is.na(ref)) {
        lower <- min(lower, ref)
        upper <- max(upper, ref)
      }
      covariance_x_limits[[block_name]][[entry]] <- c(lower, upper * 1.05)
    }
  }
  
  # Unify x-axis range within each covariance block so that every entry in
  # an effect-specific grid shares a single common x-axis range (Chapter 8).
  for (block_name in names(covariance_x_limits)) {
    entry_lims <- covariance_x_limits[[block_name]]
    all_lowers <- vapply(entry_lims, function(lim) if (!is.null(lim)) lim[1] else NA_real_, numeric(1))
    all_uppers <- vapply(entry_lims, function(lim) if (!is.null(lim)) lim[2] else NA_real_, numeric(1))
    if (any(is.finite(all_lowers)) && any(is.finite(all_uppers))) {
      unified <- c(min(all_lowers, na.rm = TRUE), max(all_uppers, na.rm = TRUE))
      for (entry in names(entry_lims)) {
        if (!is.null(entry_lims[[entry]])) {
          covariance_x_limits[[block_name]][[entry]] <- unified
        }
      }
    }
  }
  
  compute_hist_density_max <- function(sample_data, n_bins = 16) {
    sample_data <- as.numeric(sample_data)
    sample_data <- sample_data[is.finite(sample_data)]
    if (length(sample_data) < 2) return(0)
    if (min(sample_data) == max(sample_data)) return(0)
    brks <- seq(min(sample_data), max(sample_data), length.out = n_bins)
    h <- tryCatch(graphics::hist(sample_data, breaks = brks, plot = FALSE), error = function(e) NULL)
    if (is.null(h)) return(0)
    max(h$density, na.rm = TRUE)
  }
  
  compute_hist_count_max <- function(sample_data, n_bins = 16) {
    sample_data <- as.numeric(sample_data)
    sample_data <- sample_data[is.finite(sample_data)]
    if (length(sample_data) < 2) return(0)
    if (min(sample_data) == max(sample_data)) return(0)
    brks <- seq(min(sample_data), max(sample_data), length.out = n_bins)
    h <- tryCatch(graphics::hist(sample_data, breaks = brks, plot = FALSE), error = function(e) NULL)
    if (is.null(h)) return(0)
    max(h$counts, na.rm = TRUE)
  }
  
  scalar_y_limits <- stats::setNames(vector("list", length(scalar_hist_params)), scalar_hist_params)
  for (param in scalar_hist_params) {
    max_density <- 0
    for (cyc in fit_object$cycles) {
      draws <- cyc$posterior_samples[[param]]
      if (is.null(draws)) next
      d <- compute_hist_density_max(draws)
      if (d > max_density) max_density <- d
    }
    if (max_density > 0) {
      scalar_y_limits[[param]] <- c(0, max_density * 1.05)
    }
  }
  
  covariance_y_limits <- list()
  for (block_name in names(covariance_blocks)) {
    entries <- covariance_blocks[[block_name]]$entries
    block_param <- covariance_blocks[[block_name]]$param_name
    covariance_y_limits[[block_name]] <- stats::setNames(vector("list", length(entries)), entries)
    for (entry in entries) {
      max_count <- 0
      idx_str <- gsub(paste0("^", block_param, "\\[(.*)\\]$"), "\\1", entry)
      ij <- as.integer(strsplit(idx_str, ",")[[1]])
      for (cyc in fit_object$cycles) {
        arr <- cyc$posterior_samples[[block_param]]
        if (is.null(arr)) next
        draws <- arr[, ij[1], ij[2]]
        ct <- compute_hist_count_max(draws)
        if (ct > max_count) max_count <- ct
      }
      if (max_count > 0) {
        covariance_y_limits[[block_name]][[entry]] <- c(0, max_count * 1.05)
      }
    }
  }
  
  list(
    fit_object = fit_object,
    reference_values = reference_values,
    digits = digits,
    print_lines = print_lines,
    summary_lines = summary_lines,
    model_effects_df = model_effects_df,
    cycle_data_summary_df = cycle_data_summary_df,
    cycle_prior_tables = cycle_prior_tables,
    final_posterior_update_df = final_posterior_update_df,
    global_settings = global_settings,
    cycle_mcmc_settings = cycle_mcmc_settings,
    cycle_summary_df = cycle_summary_df,
    scalar_hist_params = scalar_hist_params,
    scalar_trace_params = scalar_trace_params,
    scalar_evolution_df = scalar_evolution_df,
    scalar_x_limits = scalar_x_limits,
    covariance_blocks = covariance_blocks,
    covariance_x_limits = covariance_x_limits,
    scalar_y_limits = scalar_y_limits,
    covariance_y_limits = covariance_y_limits,
    residual_trace_source = residual_trace_source,
    residual_trace_params = residual_trace_params,
    pair_params = pair_params,
    cycle_timing_df = cycle_timing_df,
    total_computing_time_min = total_computing_time_min,
    cycle_labels = stats::setNames(vapply(seq_along(fit_object$cycles), function(i) cycle_label(fit_object$cycles[[i]], i), character(1)), seq_along(fit_object$cycles))
  )
}

render_bayesian_cycle_report_core <- function(
    report_spec,
    output_file = "bayesian_multiple_cycles_report.html",
    report_title = "Bayesian mixed-model multi-window summary",
    base_family = "serif",
    image_width = 1800,
    image_height = 1100,
    image_res = 170,
    autocorr_lag_max = 30
) {
  `%||%` <- function(x, y) {
    if (is.null(x) || (length(x) == 1 && is.na(x))) y else x
  }
  
  required_pkgs <- c(
    "ggplot2", "dplyr", "htmltools", "knitr", "patchwork",
    "rstan", "base64enc"
  )
  missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(missing_pkgs) > 0) {
    stop(
      "The following packages are required but not installed: ",
      paste(missing_pkgs, collapse = ", "),
      call. = FALSE
    )
  }
  
  resolve_digits <- function(d) {
    if (is.null(d) || length(d) == 0 || is.na(d[1]) || !is.finite(d[1])) return(4L)
    as.integer(d[1])
  }
  
  fit_object <- report_spec$fit_object
  reference_values <- report_spec$reference_values %||% numeric(0)
  digits <- resolve_digits(report_spec$digits)
  
  format_num <- function(x, digits = NULL) {
    digits <- resolve_digits(digits)
    if (length(x) == 0 || all(is.na(x))) return(NA_character_)
    if (inherits(x, "difftime")) x <- as.numeric(x)
    if (is.logical(x)) return(as.character(x))
    if (is.numeric(x)) {
      rounded <- round(x, digits)
      if (isTRUE(all.equal(rounded, as.integer(rounded), tolerance = 1e-12))) {
        return(format(as.integer(rounded), trim = TRUE, scientific = FALSE))
      }
      return(format(rounded, nsmall = digits, trim = TRUE, scientific = FALSE))
    }
    as.character(x)
  }
  
  prettify_name <- function(x) {
    x <- gsub("^var_", "", x)
    x <- gsub("^Sigma_", "", x)
    x <- gsub("^sigma_", "", x)
    x <- gsub("_", " × ", x)
    x <- gsub("\\[", " [", x)
    x
  }
  
  as_table_tag <- function(df, caption = NULL, escape = TRUE) {
    if (is.null(df) || nrow(df) == 0) {
      return(htmltools::tags$p("No rows available."))
    }
    htmltools::HTML(
      knitr::kable(
        df,
        format = "html",
        caption = caption,
        table.attr = 'class="report-table"',
        escape = escape
      )
    )
  }
  
  gg_to_img_tag <- function(plot_obj, width = image_width, height = image_height, res = image_res) {
    tf <- tempfile(fileext = ".png")
    grDevices::png(tf, width = width, height = height, res = res)
    withCallingHandlers(
      print(plot_obj),
      warning = function(w) {
        if (grepl("does not name a graphical parameter", conditionMessage(w)))
          invokeRestart("muffleWarning")
      }
    )
    grDevices::dev.off()
    uri <- paste0("data:image/png;base64,", base64enc::base64encode(tf))
    htmltools::tags$img(src = uri, class = "report-img")
  }
  
  expr_to_img_tag <- function(expr, width = image_width, height = image_height, res = image_res) {
    tf <- tempfile(fileext = ".png")
    grDevices::png(tf, width = width, height = height, res = res)
    eval(expr)
    grDevices::dev.off()
    uri <- paste0("data:image/png;base64,", base64enc::base64encode(tf))
    htmltools::tags$img(src = uri, class = "report-img")
  }
  
  fit_inverse_gamma_equal_prob <- function(sample_data) {
    sample_data <- as.numeric(sample_data)
    sample_data <- sample_data[is.finite(sample_data) & sample_data > 0]
    if (length(sample_data) < 30 || stats::sd(sample_data) == 0) return(NULL)
    
    n_bins <- max(10, as.integer(length(sample_data) / 100))
    breaks <- unique(stats::quantile(sample_data, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE))
    if (length(breaks) < 5) return(NULL)
    
    mids <- tryCatch((graphics::hist(sample_data, breaks = breaks, plot = FALSE))$mids, error = function(e) NULL)
    if (is.null(mids) || length(mids) < 5) return(NULL)
    
    neg_loglik <- function(theta, data) {
      alpha <- exp(theta[1])
      beta <- exp(theta[2])
      n <- length(data)
      term1 <- n * (alpha * log(beta) - lgamma(alpha))
      term2 <- -(alpha + 1) * sum(log(data))
      term3 <- -beta * sum(1 / data)
      -(term1 + term2 + term3)
    }
    
    starts <- lapply(seq(0.1, 2.0, by = 0.1), function(z) c(log(z), log(z)))
    lower_bounds <- c(log(1e-5), log(1e-5))
    upper_bounds <- c(log(1.5e2), log(1.5e2))
    
    best_fit <- NULL
    best_value <- Inf
    for (start in starts) {
      fit <- tryCatch(
        stats::optim(
          par = start,
          fn = neg_loglik,
          data = mids,
          method = "L-BFGS-B",
          lower = lower_bounds,
          upper = upper_bounds
        ),
        error = function(e) NULL
      )
      if (!is.null(fit) && is.finite(fit$value) && fit$value < best_value) {
        best_fit <- fit
        best_value <- fit$value
      }
    }
    
    if (is.null(best_fit)) return(NULL)
    list(alpha = exp(best_fit$par[1]), beta = exp(best_fit$par[2]))
  }
  
  dinvgamma_stable <- function(x, alpha, beta) {
    ifelse(
      x <= 0,
      0,
      exp(alpha * log(beta) - lgamma(alpha) - (alpha + 1) * log(x) - beta / x)
    )
  }
  
  cycle_label <- function(i) report_spec$cycle_labels[[as.character(i)]] %||% paste0("Window ", i)
  
  get_scalar_draws <- function(cycle_obj, param_name) {
    cycle_obj$posterior_samples[[param_name]]
  }
  
  get_cov_entry_draws <- function(cycle_obj, array_name, entry_name) {
    arr <- cycle_obj$posterior_samples[[array_name]]
    if (is.null(arr)) return(NULL)
    idx <- gsub(paste0("^", array_name, "\\[(.*)\\]$"), "\\1", entry_name)
    ij <- as.integer(strsplit(idx, ",")[[1]])
    arr[, ij[1], ij[2]]
  }
  
  build_hist_plot <- function(sample_data, param_name, ref_value = NA_real_, x_limits = NULL, y_limits = NULL, ig_alpha = NULL, ig_beta = NULL) {
    df <- data.frame(value = as.numeric(sample_data))
    df <- df[is.finite(df$value), , drop = FALSE]
    if (nrow(df) == 0) {
      return(
        ggplot2::ggplot() +
          ggplot2::annotate("text", x = 1, y = 1, label = paste("No draws for", param_name)) +
          ggplot2::theme_void()
      )
    }
    
    hist_breaks <- seq(min(df$value), max(df$value), length.out = 15 + 1)
    
    p <- ggplot2::ggplot(df, ggplot2::aes(x = value)) +
      ggplot2::geom_histogram(
        ggplot2::aes(y = after_stat(density)),
        breaks = hist_breaks,
        fill = "grey80",
        color = "white"
      ) +
      ggplot2::labs(title = prettify_name(param_name), x = NULL, y = NULL) +
      ggplot2::theme_bw(base_family = base_family) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 12),
        axis.text = ggplot2::element_text(size = 9)
      )
    
    if (!is.null(ig_alpha) && !is.null(ig_beta) && is.finite(ig_alpha) && is.finite(ig_beta)) {
      x_grid <- seq(1e-8, min(10, mean(df$value) * 5), length.out = 2000)
      pdf_df <- data.frame(x = x_grid, pdf = dinvgamma_stable(x_grid, ig_alpha, ig_beta))
      p <- p + ggplot2::geom_line(
        data = pdf_df,
        ggplot2::aes(x = x, y = pdf),
        color = "blue",
        linewidth = 0.7,
        inherit.aes = FALSE
      )
    }
    
    if (!is.na(ref_value)) {
      p <- p + ggplot2::geom_vline(xintercept = ref_value, color = "red", linewidth = 0.8)
    }
    has_x <- !is.null(x_limits) && length(x_limits) == 2 && all(is.finite(x_limits))
    has_y <- !is.null(y_limits) && length(y_limits) == 2 && all(is.finite(y_limits))
    if (has_x || has_y) {
      p <- p + ggplot2::coord_cartesian(
        xlim = if (has_x) x_limits else NULL,
        ylim = if (has_y) y_limits else NULL
      )
    }
    
    p
  }
  
  build_cov_hist_plot <- function(sample_data, param_name, ref_value = NA_real_, x_limits = NULL, y_limits = NULL) {
    df <- data.frame(value = as.numeric(sample_data))
    df <- df[is.finite(df$value), , drop = FALSE]
    if (nrow(df) == 0) {
      return(
        ggplot2::ggplot() +
          ggplot2::annotate("text", x = 1, y = 1, label = paste("No draws for", param_name)) +
          ggplot2::theme_void()
      )
    }
    
    hist_breaks <- seq(min(df$value), max(df$value), length.out = 15 + 1)
    
    p <- ggplot2::ggplot(df, ggplot2::aes(x = value)) +
      ggplot2::geom_histogram(
        breaks = hist_breaks,
        fill = "grey80",
        color = "white"
      ) +
      ggplot2::labs(title = prettify_name(param_name), x = NULL, y = NULL) +
      ggplot2::theme_bw(base_family = base_family) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 12),
        axis.text = ggplot2::element_text(size = 9)
      )
    
    if (!is.na(ref_value)) {
      p <- p + ggplot2::geom_vline(xintercept = ref_value, color = "red", linewidth = 0.8)
    }
    has_x <- !is.null(x_limits) && length(x_limits) == 2 && all(is.finite(x_limits))
    has_y <- !is.null(y_limits) && length(y_limits) == 2 && all(is.finite(y_limits))
    if (has_x || has_y) {
      p <- p + ggplot2::coord_cartesian(
        xlim = if (has_x) x_limits else NULL,
        ylim = if (has_y) y_limits else NULL
      )
    }
    
    p
  }
  
  build_trace_df <- function(stan_fit, pars) {
    pars <- unique(pars)
    pars <- pars[!is.na(pars) & nzchar(pars)]
    if (length(pars) == 0) return(NULL)
    
    arr <- tryCatch(rstan::extract(stan_fit, pars = pars, permuted = FALSE, inc_warmup = FALSE), error = function(e) NULL)
    if (is.null(arr) || is.null(dim(arr))) return(NULL)
    if (length(dim(arr)) == 2) {
      arr <- array(arr, dim = c(dim(arr)[1], dim(arr)[2], 1), dimnames = list(NULL, NULL, pars))
    }
    
    pnames <- dimnames(arr)[[3]] %||% pars
    data.frame(
      Iteration = rep(seq_len(dim(arr)[1]), times = dim(arr)[2] * dim(arr)[3]),
      Chain = factor(rep(rep(seq_len(dim(arr)[2]), each = dim(arr)[1]), times = dim(arr)[3])),
      Parameter = rep(pnames, each = dim(arr)[1] * dim(arr)[2]),
      Value = as.vector(arr),
      stringsAsFactors = FALSE
    )
  }
  
  build_trace_plot <- function(stan_fit, pars, title, ncol = 3) {
    trace_df <- build_trace_df(stan_fit, pars)
    if (is.null(trace_df) || nrow(trace_df) == 0) return(NULL)
    
    ggplot2::ggplot(trace_df, ggplot2::aes(x = Iteration, y = Value, group = interaction(Chain, Parameter), color = Chain)) +
      ggplot2::geom_line(alpha = 0.7, linewidth = 0.25) +
      ggplot2::facet_wrap(~Parameter, scales = "free_y", ncol = ncol) +
      ggplot2::theme_bw(base_family = base_family) +
      ggplot2::theme(
        legend.position = "bottom",
        plot.title = ggplot2::element_text(face = "bold"),
        strip.text = ggplot2::element_text(size = 9)
      ) +
      ggplot2::labs(title = title, x = "Iteration", y = "Draw value")
  }
  
  build_acf_df <- function(stan_fit, pars, lag_max = autocorr_lag_max) {
    pars <- unique(pars)
    pars <- pars[!is.na(pars) & nzchar(pars)]
    if (length(pars) == 0) return(NULL)
    
    arr <- tryCatch(rstan::extract(stan_fit, pars = pars, permuted = FALSE, inc_warmup = FALSE), error = function(e) NULL)
    if (is.null(arr) || is.null(dim(arr))) return(NULL)
    if (length(dim(arr)) == 2) {
      arr <- array(arr, dim = c(dim(arr)[1], dim(arr)[2], 1), dimnames = list(NULL, NULL, pars))
    }
    
    pnames <- dimnames(arr)[[3]] %||% pars
    out <- list()
    idx <- 1L
    for (p in seq_along(pnames)) {
      for (ch in seq_len(dim(arr)[2])) {
        ac <- stats::acf(arr[, ch, p], plot = FALSE, lag.max = lag_max, na.action = na.pass)$acf
        if (length(ac) > 1) {
          out[[idx]] <- data.frame(
            Lag = seq_len(length(ac) - 1),
            ACF = ac[-1],
            Chain = factor(ch),
            Parameter = pnames[p],
            stringsAsFactors = FALSE
          )
          idx <- idx + 1L
        }
      }
    }
    if (length(out) == 0) return(NULL)
    dplyr::bind_rows(out)
  }
  
  build_acf_plot <- function(stan_fit, pars, title, lag_max = autocorr_lag_max, ncol = 3) {
    acf_df <- build_acf_df(stan_fit, pars, lag_max = lag_max)
    if (is.null(acf_df) || nrow(acf_df) == 0) return(NULL)
    
    ggplot2::ggplot(acf_df, ggplot2::aes(x = Lag, y = ACF, fill = Chain)) +
      ggplot2::geom_col(position = "dodge") +
      ggplot2::facet_wrap(~Parameter, scales = "free_y", ncol = ncol) +
      ggplot2::theme_bw(base_family = base_family) +
      ggplot2::theme(
        legend.position = "bottom",
        plot.title = ggplot2::element_text(face = "bold"),
        strip.text = ggplot2::element_text(size = 9)
      ) +
      ggplot2::labs(title = title, x = "Lag", y = "Autocorrelation")
  }
  
  build_pairs_tag <- function(stan_fit, pars, title) {
    pars <- unique(pars)
    pars <- pars[!is.na(pars) & nzchar(pars)]
    if (length(pars) < 2) {
      return(htmltools::tags$p("Not enough parameters available for a pairs plot."))
    }
    
    tf <- tempfile(fileext = ".jpeg")
    grDevices::jpeg(tf, width = 1200, height = 1200, res = 150, quality = 85)
    graphics::par(family = base_family)
    try(graphics::pairs(stan_fit, pars = pars), silent = TRUE)
    grDevices::dev.off()
    uri <- paste0("data:image/jpeg;base64,", base64enc::base64encode(tf))
    
    htmltools::tagList(
      htmltools::tags$h4(title),
      htmltools::tags$img(src = uri, class = "report-img")
    )
  }
  
  build_hmc_table <- function(stan_fit, max_treedepth_val = NA_real_) {
    sampler_params <- tryCatch(rstan::get_sampler_params(stan_fit, inc_warmup = FALSE), error = function(e) NULL)
    if (is.null(sampler_params)) return(data.frame(message = "Sampler parameters could not be extracted."))
    
    out <- lapply(seq_along(sampler_params), function(ch) {
      sp <- sampler_params[[ch]]
      energy <- sp[, "energy__"]
      bfmi <- if (stats::var(energy) > 0) mean(diff(energy)^2) / stats::var(energy) else NA_real_
      data.frame(
        Chain = ch,
        Divergent = sum(sp[, "divergent__"] > 0),
        Max_Treedepth_Hits = if (!is.na(max_treedepth_val) && "treedepth__" %in% colnames(sp)) sum(sp[, "treedepth__"] >= max_treedepth_val) else NA_integer_,
        Mean_Accept_Stat = mean(sp[, "accept_stat__"]),
        Mean_Stepsize = mean(sp[, "stepsize__"]),
        Mean_Leapfrog = mean(sp[, "n_leapfrog__"]),
        E_BFMI = bfmi,
        stringsAsFactors = FALSE
      )
    })
    dplyr::bind_rows(out)
  }
  
  build_gelman_geweke <- function(stan_fit, pars) {
    if (!requireNamespace("coda", quietly = TRUE)) {
      return(list(
        gelman = data.frame(message = "Package 'coda' is not installed; Gelman-Rubin diagnostics were skipped."),
        geweke = data.frame(message = "Package 'coda' is not installed; Geweke diagnostics were skipped.")
      ))
    }
    
    gelman_rows <- list()
    geweke_rows <- list()
    
    for (param in pars) {
      mcmc_list <- tryCatch(rstan::As.mcmc.list(stan_fit, pars = param), error = function(e) NULL)
      if (is.null(mcmc_list)) next
      
      gel <- tryCatch(coda::gelman.diag(mcmc_list)$psrf, error = function(e) NULL)
      if (!is.null(gel)) {
        gelman_rows[[length(gelman_rows) + 1]] <- data.frame(
          Parameter = param,
          Point_Estimate = gel[1],
          Upper_CI = gel[2],
          stringsAsFactors = FALSE
        )
      }
      
      gw <- tryCatch(coda::geweke.diag(mcmc_list), error = function(e) NULL)
      if (!is.null(gw)) {
        zvals <- unlist(lapply(gw, function(x) x$z), use.names = FALSE)
        geweke_rows[[length(geweke_rows) + 1]] <- data.frame(
          Parameter = param,
          Chain = seq_along(zvals),
          Geweke_Z = zvals,
          stringsAsFactors = FALSE
        )
      }
    }
    
    list(
      gelman = if (length(gelman_rows)) dplyr::bind_rows(gelman_rows) else data.frame(message = "No Gelman-Rubin diagnostics could be computed."),
      geweke = if (length(geweke_rows)) dplyr::bind_rows(geweke_rows) else data.frame(message = "No Geweke diagnostics could be computed.")
    )
  }
  
  summarise_stan_parameters <- function(stan_fit, pars) {
    pars <- unique(pars)
    pars <- pars[!is.na(pars) & nzchar(pars)]
    if (length(pars) == 0) return(data.frame())
    
    sm <- tryCatch(rstan::summary(stan_fit, pars = pars)$summary, error = function(e) NULL)
    if (is.null(sm)) return(data.frame())
    
    df <- data.frame(Parameter = rownames(sm), sm, row.names = NULL, check.names = FALSE)
    keep <- intersect(c("Parameter", "mean", "sd", "2.5%", "25%", "50%", "75%", "97.5%", "n_eff", "Rhat"), names(df))
    df[, keep, drop = FALSE]
  }
  
  overall_diag <- fit_object$diagnostics %||% list()
  
  sampling_time_plot_tag <- NULL
  prior_update_time_plot_tag <- NULL
  if (!is.null(report_spec$cycle_timing_df) && nrow(report_spec$cycle_timing_df) > 0) {
    timing_df <- report_spec$cycle_timing_df
    timing_df$window <- factor(timing_df$window)
    
    sampling_plot <- ggplot2::ggplot(timing_df, ggplot2::aes(x = window, y = sampling_time_min)) +
      ggplot2::geom_col(fill = "#486581", width = 0.6) +
      ggplot2::theme_bw(base_family = base_family) +
      ggplot2::labs(
        title = "Sampling time per window",
        x = "Window",
        y = "Time (minutes)"
      ) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold"),
        axis.text = ggplot2::element_text(size = 10)
      )
    sampling_time_plot_tag <- gg_to_img_tag(sampling_plot)
    
    prior_update_plot <- ggplot2::ggplot(timing_df, ggplot2::aes(x = window, y = prior_update_time_sec)) +
      ggplot2::geom_col(fill = "#334E68", width = 0.6) +
      ggplot2::theme_bw(base_family = base_family) +
      ggplot2::labs(
        title = "Prior update time per window",
        x = "Window",
        y = "Time (seconds)"
      ) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold"),
        axis.text = ggplot2::element_text(size = 10)
      )
    prior_update_time_plot_tag <- gg_to_img_tag(prior_update_plot)
  }
  
  evolution_plot_tag <- NULL
  if (!is.null(report_spec$scalar_evolution_df) && nrow(report_spec$scalar_evolution_df) > 0) {
    evolution_plot <- ggplot2::ggplot(report_spec$scalar_evolution_df, ggplot2::aes(x = window, y = q50)) +
      ggplot2::geom_linerange(ggplot2::aes(ymin = q2.5, ymax = q97.5), linewidth = 0.35) +
      ggplot2::geom_linerange(ggplot2::aes(ymin = q25, ymax = q75), linewidth = 1.00) +
      ggplot2::geom_point(size = 1.4) +
      ggplot2::facet_wrap(~parameter, scales = "free_y", ncol = 3) +
      ggplot2::theme_bw(base_family = base_family) +
      ggplot2::labs(
        title = "Posterior evolution across windows",
        x = "Window",
        y = "Posterior median with 50% and 95% intervals"
      ) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold"),
        strip.text = ggplot2::element_text(size = 9)
      )
    evolution_plot_tag <- gg_to_img_tag(evolution_plot)
  }
  
  get_evolution_ig_params <- function(cycle_index, param_name) {
    sanitize_effect_name <- function(x) gsub(":", "_", x)
    stan_name <- sub("^var_", "", param_name)
    alpha_key <- paste0("alpha_var_", stan_name)
    beta_key <- paste0("beta_var_", stan_name)
    n_cycles <- length(fit_object$cycles)
    evol <- fit_object$parameters$evolution
    final <- fit_object$parameters$final_posteriors$result %||% fit_object$parameters$final_posteriors %||% list()
    if (cycle_index <= length(evol)) {
      src <- evol[[cycle_index]]
    } else if (cycle_index == n_cycles) {
      src <- final
    } else {
      return(list(alpha = NULL, beta = NULL))
    }
    if (is.null(src)) return(list(alpha = NULL, beta = NULL))
    alpha_val <- src[[alpha_key]]
    beta_val <- src[[beta_key]]
    if (is.null(alpha_val) || is.null(beta_val)) return(list(alpha = NULL, beta = NULL))
    list(alpha = alpha_val, beta = beta_val)
  }
  
  scalar_sections <- lapply(seq_along(fit_object$cycles), function(i) {
    cyc <- fit_object$cycles[[i]]
    if (length(report_spec$scalar_hist_params) == 0) {
      return(htmltools::tagList(
        htmltools::tags$h3(cycle_label(i)),
        htmltools::tags$p("No scalar inverse-gamma variance components were detected.")
      ))
    }
    
    plots <- lapply(report_spec$scalar_hist_params, function(param) {
      ref <- if (param %in% names(reference_values)) reference_values[[param]] else NA_real_
      ig <- get_evolution_ig_params(i, param)
      build_hist_plot(
        sample_data = get_scalar_draws(cyc, param),
        param_name = param,
        ref_value = ref,
        x_limits = report_spec$scalar_x_limits[[param]],
        y_limits = report_spec$scalar_y_limits[[param]],
        ig_alpha = ig$alpha,
        ig_beta = ig$beta
      )
    })
    
    grid_ncol <- if (length(plots) <= 2) length(plots) else 2L
    combined <- patchwork::wrap_plots(plots, ncol = grid_ncol) + patchwork::plot_annotation(title = cycle_label(i))
    
    htmltools::tagList(
      htmltools::tags$h3(cycle_label(i)),
      gg_to_img_tag(combined, height = 1400)
    )
  })
  
  covariance_window_sections <- lapply(seq_along(fit_object$cycles), function(i) {
    cyc <- fit_object$cycles[[i]]
    block_tags <- lapply(names(report_spec$covariance_blocks), function(block_name) {
      block_info <- report_spec$covariance_blocks[[block_name]]
      plots <- lapply(block_info$entries, function(entry_nm) {
        ref <- if (entry_nm %in% names(reference_values)) reference_values[[entry_nm]] else NA_real_
        build_cov_hist_plot(
          sample_data = get_cov_entry_draws(cyc, block_info$param_name, entry_nm),
          param_name = entry_nm,
          ref_value = ref,
          x_limits = report_spec$covariance_x_limits[[block_name]][[entry_nm]],
          y_limits = report_spec$covariance_y_limits[[block_name]][[entry_nm]]
        )
      })
      
      combined <- patchwork::wrap_plots(plots, ncol = 2) +
        patchwork::plot_annotation(title = paste0(cycle_label(i), " \u2014 ", block_info$param_name))
      
      htmltools::tagList(
        htmltools::tags$h4(paste0(cycle_label(i), " \u2014 ", block_info$param_name)),
        gg_to_img_tag(combined, height = 1700)
      )
    })
    htmltools::tagList(block_tags)
  })
  
  build_per_chain_trace_plot <- function(stan_fit, pars, chain_id, title) {
    pars <- unique(pars)
    pars <- pars[!is.na(pars) & nzchar(pars)]
    if (length(pars) == 0) return(NULL)
    arr <- tryCatch(rstan::extract(stan_fit, pars = pars, permuted = FALSE, inc_warmup = FALSE), error = function(e) NULL)
    if (is.null(arr) || is.null(dim(arr))) return(NULL)
    if (length(dim(arr)) == 2) {
      arr <- array(arr, dim = c(dim(arr)[1], dim(arr)[2], 1), dimnames = list(NULL, NULL, pars))
    }
    n_chains <- dim(arr)[2]
    if (chain_id > n_chains) return(NULL)
    pnames <- dimnames(arr)[[3]] %||% pars
    trace_df <- data.frame(
      Iteration = rep(seq_len(dim(arr)[1]), times = length(pnames)),
      Parameter = rep(pnames, each = dim(arr)[1]),
      Value = as.vector(arr[, chain_id, ]),
      stringsAsFactors = FALSE
    )
    if (nrow(trace_df) == 0) return(NULL)
    ggplot2::ggplot(trace_df, ggplot2::aes(x = Iteration, y = Value)) +
      ggplot2::geom_line(alpha = 0.7, linewidth = 0.25, color = "#333333") +
      ggplot2::facet_wrap(~Parameter, scales = "free_y", ncol = 3) +
      ggplot2::theme_bw(base_family = base_family) +
      ggplot2::theme(
        legend.position = "none",
        plot.title = ggplot2::element_text(face = "bold"),
        strip.text = ggplot2::element_text(size = 9)
      ) +
      ggplot2::labs(title = title, x = "Iteration", y = "Draw value")
  }
  
  build_per_chain_acf_plot <- function(stan_fit, pars, chain_id, title, lag_max = autocorr_lag_max) {
    pars <- unique(pars)
    pars <- pars[!is.na(pars) & nzchar(pars)]
    if (length(pars) == 0) return(NULL)
    arr <- tryCatch(rstan::extract(stan_fit, pars = pars, permuted = FALSE, inc_warmup = FALSE), error = function(e) NULL)
    if (is.null(arr) || is.null(dim(arr))) return(NULL)
    if (length(dim(arr)) == 2) {
      arr <- array(arr, dim = c(dim(arr)[1], dim(arr)[2], 1), dimnames = list(NULL, NULL, pars))
    }
    n_chains <- dim(arr)[2]
    if (chain_id > n_chains) return(NULL)
    pnames <- dimnames(arr)[[3]] %||% pars
    out <- list()
    idx <- 1L
    for (p in seq_along(pnames)) {
      ac <- stats::acf(arr[, chain_id, p], plot = FALSE, lag.max = lag_max, na.action = na.pass)$acf
      if (length(ac) > 1) {
        out[[idx]] <- data.frame(
          Lag = seq_len(length(ac) - 1),
          ACF = ac[-1],
          Parameter = pnames[p],
          stringsAsFactors = FALSE
        )
        idx <- idx + 1L
      }
    }
    if (length(out) == 0) return(NULL)
    acf_df <- dplyr::bind_rows(out)
    ggplot2::ggplot(acf_df, ggplot2::aes(x = Lag, y = ACF)) +
      ggplot2::geom_col(fill = "#666666") +
      ggplot2::facet_wrap(~Parameter, scales = "free_y", ncol = 3) +
      ggplot2::theme_bw(base_family = base_family) +
      ggplot2::theme(
        legend.position = "none",
        plot.title = ggplot2::element_text(face = "bold"),
        strip.text = ggplot2::element_text(size = 9)
      ) +
      ggplot2::labs(title = title, x = "Lag", y = "Autocorrelation")
  }
  
  collect_all_stan_params <- function(stan_fit) {
    sm <- tryCatch(rstan::summary(stan_fit)$summary, error = function(e) NULL)
    if (is.null(sm)) return(character(0))
    rownames(sm)
  }
  
  categorize_parameters <- function(all_pars, fit_object) {
    random_effects <- fit_object$model_spec$random_effects %||% list()
    fixed_effects <- fit_object$model_spec$fixed_effects %||% list()
    sanitize_effect_name <- function(x) gsub(":", "_", x)
    
    ig_var_pars <- character(0)
    fixed_pars <- character(0)
    random_pars <- character(0)
    resid_pars <- character(0)
    other_pars <- character(0)
    
    for (p in all_pars) {
      if (grepl("^lp__$", p)) next
      if (grepl("^var_resid_env\\[", p) || grepl("^sigma_resid_env\\[", p) || grepl("^var_resid_env$", p) || grepl("^sigma_resid_env$", p) || grepl("^var_resid_env_mean", p)) {
        resid_pars <- c(resid_pars, p)
      } else if (grepl("^var_", p)) {
        ig_var_pars <- c(ig_var_pars, p)
      } else if (grepl("^Sigma_", p)) {
        ig_var_pars <- c(ig_var_pars, p)
      } else if (grepl("^beta_", p)) {
        fixed_pars <- c(fixed_pars, p)
      } else if (grepl("^u_", p)) {
        random_pars <- c(random_pars, p)
      } else {
        other_pars <- c(other_pars, p)
      }
    }
    
    list(
      random_effect_variances = ig_var_pars,
      fixed_effects = fixed_pars,
      random_effects = random_pars,
      residual_error_variances = resid_pars,
      other = other_pars
    )
  }
  
  diagnostics_sections <- lapply(seq_along(fit_object$cycles), function(i) {
    cyc <- fit_object$cycles[[i]]
    stan_fit <- cyc$stan_fit
    
    all_pars <- collect_all_stan_params(stan_fit)
    par_cats <- categorize_parameters(all_pars, fit_object)
    
    # Drop redundant sigma_resid_env (standard-deviation) parameters;
    # var_resid_env (variance) already carries the same information.
    resid_pars_filtered <- par_cats$residual_error_variances[
      !grepl("^sigma_resid_env", par_cats$residual_error_variances)
    ]
    
    summary_rv <- summarise_stan_parameters(stan_fit, par_cats$random_effect_variances)
    summary_fe <- summarise_stan_parameters(stan_fit, par_cats$fixed_effects)
    summary_resid <- summarise_stan_parameters(stan_fit, resid_pars_filtered)
    
    geweke_rows <- list()
    geweke_pars <- unique(c(par_cats$random_effect_variances, par_cats$fixed_effects, resid_pars_filtered))
    for (param in geweke_pars) {
      mcmc_list <- tryCatch(rstan::As.mcmc.list(stan_fit, pars = param), error = function(e) NULL)
      if (is.null(mcmc_list)) next
      gw <- tryCatch(coda::geweke.diag(mcmc_list), error = function(e) NULL)
      if (!is.null(gw)) {
        zvals <- unlist(lapply(gw, function(x) x$z), use.names = FALSE)
        geweke_rows[[length(geweke_rows) + 1]] <- data.frame(
          Parameter = param,
          Chain = seq_along(zvals),
          Geweke_Z = zvals,
          stringsAsFactors = FALSE
        )
      }
    }
    geweke_df <- if (length(geweke_rows)) dplyr::bind_rows(geweke_rows) else data.frame(message = "No Geweke diagnostics could be computed.")
    
    hmc_table <- build_hmc_table(stan_fit, max_treedepth_val = cyc$mcmc_settings$max_treedepth %||% NA_real_)
    
    n_chains <- tryCatch(stan_fit@sim$chains, error = function(e) 4L)
    if (is.null(n_chains)) n_chains <- 4L
    
    chain_trace_tags <- lapply(seq_len(n_chains), function(ch) {
      group_tags <- htmltools::tagList()
      
      if (length(report_spec$scalar_trace_params) > 0) {
        tp_scalar <- build_per_chain_trace_plot(stan_fit, report_spec$scalar_trace_params,
                                                ch, paste0(cycle_label(i), " \u2014 chain ", ch, " trace: scalar random effect variances"))
        ap_scalar <- build_per_chain_acf_plot(stan_fit, report_spec$scalar_trace_params,
                                              ch, paste0(cycle_label(i), " \u2014 chain ", ch, " ACF: scalar random effect variances"),
                                              lag_max = autocorr_lag_max)
        group_tags <- htmltools::tagList(
          group_tags,
          if (!is.null(tp_scalar)) gg_to_img_tag(tp_scalar, height = 1300) else htmltools::tags$p("No scalar variance trace data available."),
          if (!is.null(ap_scalar)) gg_to_img_tag(ap_scalar, height = 1300) else htmltools::tags$p("No scalar variance ACF data available.")
        )
      }
      
      for (block_name in names(report_spec$covariance_blocks)) {
        block_info <- report_spec$covariance_blocks[[block_name]]
        if (length(block_info$entries) > 0) {
          tp_cov <- build_per_chain_trace_plot(stan_fit, block_info$entries,
                                               ch, paste0(cycle_label(i), " \u2014 chain ", ch, " trace: ", block_info$param_name))
          ap_cov <- build_per_chain_acf_plot(stan_fit, block_info$entries,
                                             ch, paste0(cycle_label(i), " \u2014 chain ", ch, " ACF: ", block_info$param_name),
                                             lag_max = autocorr_lag_max)
          group_tags <- htmltools::tagList(
            group_tags,
            if (!is.null(tp_cov)) gg_to_img_tag(tp_cov, height = 1300) else htmltools::tags$p(paste0("No trace data for ", block_info$param_name, ".")),
            if (!is.null(ap_cov)) gg_to_img_tag(ap_cov, height = 1300) else htmltools::tags$p(paste0("No ACF data for ", block_info$param_name, "."))
          )
        }
      }
      
      if (length(report_spec$residual_trace_params) > 0) {
        resid_trace_filtered <- report_spec$residual_trace_params[
          !grepl("^sigma_resid_env", report_spec$residual_trace_params)
        ]
        if (length(resid_trace_filtered) > 0) {
          tp_resid <- build_per_chain_trace_plot(stan_fit, resid_trace_filtered,
                                                 ch, paste0(cycle_label(i), " \u2014 chain ", ch, " trace: var_resid_env"))
          ap_resid <- build_per_chain_acf_plot(stan_fit, resid_trace_filtered,
                                               ch, paste0(cycle_label(i), " \u2014 chain ", ch, " ACF: var_resid_env"),
                                               lag_max = autocorr_lag_max)
          group_tags <- htmltools::tagList(
            group_tags,
            if (!is.null(tp_resid)) gg_to_img_tag(tp_resid, height = 1300) else htmltools::tags$p("No residual trace data available."),
            if (!is.null(ap_resid)) gg_to_img_tag(ap_resid, height = 1300) else htmltools::tags$p("No residual ACF data available.")
          )
        }
      }
      
      htmltools::tags$details(
        htmltools::tags$summary(paste0("Chain ", ch)),
        group_tags
      )
    })
    
    pairs_tags <- htmltools::tagList()
    
    if (length(report_spec$scalar_trace_params) >= 2) {
      pairs_tags <- htmltools::tagList(
        pairs_tags,
        build_pairs_tag(stan_fit, report_spec$scalar_trace_params, paste0(cycle_label(i), " \u2014 pairs plot: scalar random effect variances"))
      )
    }
    
    for (block_name in names(report_spec$covariance_blocks)) {
      block_info <- report_spec$covariance_blocks[[block_name]]
      entries <- block_info$entries
      diag_entries <- entries[vapply(entries, function(e) {
        idx <- sub("^.*\\[([0-9]+),([0-9]+)\\]$", "\\1,\\2", e)
        parts <- strsplit(idx, ",", fixed = TRUE)[[1]]
        length(parts) == 2 && identical(parts[1], parts[2])
      }, logical(1))]
      first_row_entries <- entries[grepl("\\[1,", entries)]
      
      if (length(diag_entries) >= 2) {
        pairs_tags <- htmltools::tagList(
          pairs_tags,
          build_pairs_tag(stan_fit, diag_entries, paste0(cycle_label(i), " \u2014 pairs plot: ", block_info$param_name, " diagonal"))
        )
      }
      if (length(first_row_entries) >= 2) {
        pairs_tags <- htmltools::tagList(
          pairs_tags,
          build_pairs_tag(stan_fit, first_row_entries, paste0(cycle_label(i), " \u2014 pairs plot: ", block_info$param_name, " first row"))
        )
      }
    }
    
    htmltools::tags$details(
      htmltools::tags$summary(class = "window-summary", cycle_label(i)),
      htmltools::tags$div(class = "note-box",
                          htmltools::tags$p(
                            htmltools::tags$strong("Stored window diagnostics: "),
                            paste0(
                              "converged = ", fit_object$diagnostics$stan_convergence[[i]]$stan_converged %||% NA,
                              "; rhat_max = ", format_num(fit_object$diagnostics$stan_convergence[[i]]$rhat_max %||% NA_real_),
                              "; ess_min = ", format_num(fit_object$diagnostics$stan_convergence[[i]]$ess_min %||% NA_real_),
                              "; divergent_transitions = ", fit_object$diagnostics$stan_convergence[[i]]$divergent_transitions %||% NA_real_,
                              "; max_treedepth_hits = ", cyc$diagnostics$max_treedepth_hits %||% NA_real_,
                              "; ebfmi_min = ", format_num(cyc$diagnostics$ebfmi_min %||% NA_real_)
                            )
                          )
      ),
      htmltools::tags$h4("Stan parameter summary \u2014 random effect variances"),
      as_table_tag(summary_rv),
      htmltools::tags$h4("Stan parameter summary \u2014 fixed effects"),
      as_table_tag(summary_fe),
      htmltools::tags$h4("Stan parameter summary \u2014 residual error variances"),
      as_table_tag(summary_resid),
      htmltools::tags$h4("HMC diagnostics"),
      as_table_tag(hmc_table),
      htmltools::tags$details(
        htmltools::tags$summary("Geweke diagnostics"),
        as_table_tag(geweke_df)
      ),
      htmltools::tags$h4("Trace and autocorrelation plots by chain"),
      chain_trace_tags,
      htmltools::tags$h4("Pairs plots"),
      pairs_tags
    )
  })
  
  stan_code_tag <- NULL
  if (!is.null(fit_object$model_spec$stan_model)) {
    stan_code_tag <- htmltools::tags$details(
      htmltools::tags$summary("Show generated Stan code from the fit_general_multiple_cycles R function pipeline."),
      htmltools::tags$pre(fit_object$model_spec$stan_model)
    )
  }
  
  session_df <- data.frame(
    field = c("timestamp", "seed", "r_version", "rstan_version", "asreml_version"),
    value = c(
      as.character(fit_object$session_info$timestamp %||% NA),
      fit_object$session_info$seed %||% NA,
      fit_object$session_info$r_version %||% NA,
      fit_object$session_info$package_versions$rstan %||% NA,
      fit_object$session_info$package_versions$asreml %||% NA
    ),
    stringsAsFactors = FALSE
  )
  
  variable_catalog_df <- data.frame(
    Variable = c(
      "var_<effect>",
      "Sigma_<effect>",
      "sigma_resid_env / var_resid_env",
      "beta_<variable>",
      "u_<effect>",
      "alpha (inverse-gamma)",
      "beta (inverse-gamma)",
      "nu (inverse-Wishart)",
      "psi (inverse-Wishart)",
      "n_eff",
      "Rhat",
      "Geweke Z",
      "E-BFMI",
      "divergent transitions",
      "max treedepth hits"
    ),
    Description = c(
      "Scalar variance component with an inverse-gamma prior. One per random effect that is modelled independently (e.g. var_year, var_trial).",
      "Covariance matrix sampled from an inverse-Wishart prior. One per block of correlated random effects (e.g. Sigma_genotype_year). Upper-triangular entries include variances on the diagonal and covariances off-diagonal.",
      "Environment-specific residual error variances (or standard deviations). A separate variance is estimated for each environment level to allow heteroscedastic residuals.",
      "Fixed-effect regression coefficients. beta_0 is the intercept; beta_<variable> are coefficients for categorical or continuous covariates using reference coding.",
      "Best linear unbiased predictor (BLUP) draws for the random effect named <effect>. These are the posterior samples of individual random-effect levels.",
      "Shape parameter of the inverse-gamma prior or posterior. Larger alpha concentrates mass away from zero, implying a more informative prior on the variance.",
      "Rate (scale) parameter of the inverse-gamma distribution. Together with alpha, determines the prior expected variance: E[sigma^2] = beta / (alpha - 1) for alpha > 1.",
      "Degrees-of-freedom parameter of the inverse-Wishart distribution. Must exceed p + 1 (p = matrix dimension). Larger nu yields a tighter prior around psi / (nu - p - 1).",
      "Scale matrix of the inverse-Wishart distribution. Its structure encodes the expected covariance pattern. The prior expectation is psi / (nu - p - 1).",
      "Effective sample size from HMC. Values much smaller than the nominal post-warmup draws indicate high autocorrelation; n_eff > 100 per parameter is a common minimum threshold.",
      "Potential scale reduction factor (Gelman-Rubin statistic). Values close to 1.00 indicate chain convergence; Rhat > 1.01 is usually considered a warning.",
      "Geweke convergence diagnostic comparing early and late portions of each chain. Approximate standard-normal under convergence; |Z| > 2 suggests the chain has not converged.",
      "Energy Bayesian Fraction of Missing Information. Values below 0.2 indicate that the HMC sampler may have difficulty exploring the posterior efficiently.",
      "Number of HMC transitions that diverged. Non-zero counts indicate the sampler encountered regions of high curvature; model reparametrisation may be needed.",
      "Number of transitions that hit the maximum tree depth. Frequent hits mean the sampler is capping exploration steps and may need a larger max_treedepth."
    ),
    stringsAsFactors = FALSE
  )
  
  n_windows <- length(fit_object$cycles)
  
  report_doc <- htmltools::tags$html(
    htmltools::tags$head(
      htmltools::tags$meta(charset = "utf-8"),
      htmltools::tags$title(report_title),
      htmltools::tags$style(htmltools::HTML(
        "
        body {
          font-family: Georgia, 'Times New Roman', serif;
          margin: 32px auto;
          max-width: 1260px;
          line-height: 1.5;
          color: #1f1f1f;
          padding: 0 16px 48px 16px;
          background: #ffffff;
        }
        h1, h2, h3, h4 { color: #102a43; }
        h1 { border-bottom: 3px solid #d9e2ec; padding-bottom: 12px; }
        h2 { margin-top: 34px; border-bottom: 1px solid #e6edf3; padding-bottom: 6px; }
        .subtitle { color: #486581; margin-top: -10px; }
        .note-box {
          background: #f7fafc;
          border-left: 5px solid #486581;
          padding: 10px 14px;
          margin: 14px 0;
        }
        .report-img {
          width: 100%;
          height: auto;
          border: 1px solid #d9e2ec;
          margin: 12px 0 20px 0;
          box-shadow: 0 1px 4px rgba(0,0,0,0.06);
        }
        table.report-table {
          border-collapse: collapse;
          width: 100%;
          margin: 12px 0 22px 0;
          font-size: 0.94rem;
        }
        table.report-table th, table.report-table td {
          border: 1px solid #d9e2ec;
          padding: 8px 10px;
          vertical-align: top;
        }
        table.report-table th {
          background: #f0f4f8;
          text-align: left;
        }
        code, pre {
          background: #f7fafc;
          padding: 2px 4px;
          border-radius: 3px;
        }
        pre {
          padding: 12px;
          overflow-x: auto;
          white-space: pre-wrap;
        }
        pre.cell-pre {
          margin: 0;
          padding: 0;
          background: transparent;
          border-radius: 0;
          white-space: pre-wrap;
          font-family: 'Courier New', Courier, monospace;
        }
        details { margin: 12px 0 18px 0; }
        summary { cursor: pointer; font-weight: 600; }
        summary.chapter-summary { font-size: 1.3em; color: #102a43; border-bottom: 1px solid #e6edf3; padding-bottom: 6px; margin-top: 34px; }
        summary.window-summary { font-size: 1.1em; color: #243b53; border-bottom: 1px solid #e6edf3; padding-bottom: 4px; margin-top: 18px; }
        .small { color: #617d98; font-size: 0.92rem; }
        .carousel-container { position: relative; }
        .carousel-slide { display: none; }
        .carousel-slide.active { display: block; }
        .carousel-nav { display: flex; justify-content: center; align-items: center; gap: 12px; margin: 10px 0; }
        .carousel-btn { background: #486581; color: #fff; border: none; padding: 6px 16px; border-radius: 4px; cursor: pointer; font-size: 1rem; }
        .carousel-btn:hover { background: #334E68; }
        .carousel-label { font-weight: 600; min-width: 120px; text-align: center; }
        "
      )),
      htmltools::tags$script(htmltools::HTML("
        function initCarousel(containerId) {
          var container = document.getElementById(containerId);
          if (!container) return;
          var slides = container.querySelectorAll('.carousel-slide');
          var label = container.querySelector('.carousel-label');
          var idx = 0;
          function show(n) {
            if (n < 0) n = 0;
            if (n >= slides.length) n = slides.length - 1;
            for (var j = 0; j < slides.length; j++) slides[j].classList.remove('active');
            slides[n].classList.add('active');
            idx = n;
            if (label) label.textContent = 'Window ' + (idx + 1) + ' of ' + slides.length;
          }
          container.querySelector('.carousel-prev').addEventListener('click', function() { show(idx - 1); });
          container.querySelector('.carousel-next').addEventListener('click', function() { show(idx + 1); });
          show(0);
        }
        document.addEventListener('DOMContentLoaded', function() {
          var carousels = document.querySelectorAll('.carousel-container');
          carousels.forEach(function(c) { initCarousel(c.id); });
        });
      "))
    ),
    htmltools::tags$body(
      htmltools::tags$h1(report_title),
      htmltools::tags$p(
        class = "subtitle",
        paste0(
          "Report generated on ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
          " from a fit_general_multiple_cycles R function pipeline output object with ", n_windows, " fitted window(s) calculated at ,",
          format(cycle_fits_US_zones_model$session_info$timestamp, "%Y-%m-%d %H:%M:%S"), "."
        )
      ),
      
      htmltools::tags$h2("1. Model overview"),
      htmltools::tags$div(
        class = "note-box",
        htmltools::tags$p(htmltools::tags$strong("Overall converged flag: "), as.character(overall_diag$overall_converged %||% NA)),
        htmltools::tags$p(htmltools::tags$strong("Total divergent transitions: "), overall_diag$total_divergent %||% NA_real_),
        htmltools::tags$p(htmltools::tags$strong("Max Rhat across windows: "), format_num(overall_diag$max_rhat %||% NA_real_)),
        htmltools::tags$p(htmltools::tags$strong("Min ESS across windows: "), format_num(overall_diag$min_ess %||% NA_real_)),
        htmltools::tags$p(htmltools::tags$strong("Total observations: "), fit_object$model_spec$total_observations %||% sum(vapply(fit_object$cycles, function(x) x$data_summary$n_obs %||% 0, numeric(1)))),
        htmltools::tags$p(htmltools::tags$strong("Total computing time (min): "), format_num(report_spec$total_computing_time_min))
      ),
      stan_code_tag,
      
      htmltools::tags$details(
        htmltools::tags$summary(class = "chapter-summary", "2. Data summary"),
        htmltools::tags$h3("Specified mixed-model effects"),
        as_table_tag(report_spec$model_effects_df),
        htmltools::tags$h3("Observations and variable levels by window"),
        as_table_tag(report_spec$cycle_data_summary_df)
      ),
      
      htmltools::tags$details(
        htmltools::tags$summary(class = "chapter-summary", "3. MCMC hyperparameters and session information"),
        htmltools::tags$h3("Global information"),
        as_table_tag(report_spec$global_settings),
        htmltools::tags$h3("Window-specific MCMC settings"),
        as_table_tag(report_spec$cycle_mcmc_settings),
        htmltools::tags$h3("Session information"),
        as_table_tag(session_df)
      ),
      
      htmltools::tags$details(
        htmltools::tags$summary(class = "chapter-summary", "4. Prior specification and updating"),
        htmltools::tagList(
          lapply(seq_along(report_spec$cycle_prior_tables), function(i) {
            htmltools::tagList(
              htmltools::tags$h3(paste0("Priors used in ", cycle_label(i))),
              as_table_tag(report_spec$cycle_prior_tables[[i]], escape = FALSE)
            )
          })
        ),
        htmltools::tags$h3("Priors after final window update"),
        as_table_tag(report_spec$final_posterior_update_df, escape = FALSE)
      ),
      
      htmltools::tags$details(
        htmltools::tags$summary(class = "chapter-summary", "5. Window convergence summary"),
        as_table_tag(report_spec$cycle_summary_df),
        sampling_time_plot_tag %||% htmltools::tags$p("No sampling time data available."),
        prior_update_time_plot_tag %||% htmltools::tags$p("No prior update time data available.")
      ),
      
      htmltools::tags$details(
        htmltools::tags$summary(class = "chapter-summary", "6. Posterior evolution across windows"),
        htmltools::tags$p(
          class = "small",
          "For each scalar inverse-gamma variance component, the thick interval is the central 50 % posterior interval and the thin interval is the central 95 % interval."
        ),
        evolution_plot_tag %||% htmltools::tags$p("No scalar variance parameters were available for an evolution plot.")
      ),
      
      htmltools::tags$details(
        htmltools::tags$summary(class = "chapter-summary", "7. Posterior distribution plots \u2014 inverse-gamma variance components"),
        htmltools::tags$p(
          class = "small",
          "Histograms of posterior draws for var_<effect> parameters. The blue curve is the inverse-gamma density with the parametrisation obtained during the Bayesian updating pipeline. Named reference REML estimated values are overlaid as red vertical lines. On the y-axis is the density function of the inverse gamma, on the x-axis is the values of the MCMC drawn variance components. Note that in all histograms the breaks are fixed to 15. Wider bars hint to a wider spread of the values plotted along the x-axis."
        ),
        htmltools::tags$div(
          id = "carousel-scalar", class = "carousel-container",
          htmltools::tags$div(class = "carousel-nav",
                              htmltools::tags$button(class = "carousel-btn carousel-prev", "\u25C0 Previous"),
                              htmltools::tags$span(class = "carousel-label", ""),
                              htmltools::tags$button(class = "carousel-btn carousel-next", "Next \u25B6")),
          htmltools::tagList(lapply(seq_along(scalar_sections), function(k) {
            htmltools::tags$div(class = if (k == 1) "carousel-slide active" else "carousel-slide", scalar_sections[[k]])
          }))
        )
      ),
      
      htmltools::tags$details(
        htmltools::tags$summary(class = "chapter-summary", "8. Posterior distribution plots \u2014 inverse-Wishart covariance blocks"),
        htmltools::tags$p(
          class = "small",
          "Histograms of upper-triangular entries of Sigma_<effect> posterior arrays associated with inverse-Wishart priors. Named reference REML estimated values are overlaid as red vertical lines. On the y-axis is the posterior sample count, on the x-axis is the values of the MCMC drawn variance components. Note that in all histograms the breaks are fixed to 15. Wider bars hint to a wider spread of the values plotted along the x-axis."
        ),
        if (length(covariance_window_sections) > 0) {
          htmltools::tags$div(
            id = "carousel-cov", class = "carousel-container",
            htmltools::tags$div(class = "carousel-nav",
                                htmltools::tags$button(class = "carousel-btn carousel-prev", "\u25C0 Previous"),
                                htmltools::tags$span(class = "carousel-label", ""),
                                htmltools::tags$button(class = "carousel-btn carousel-next", "Next \u25B6")),
            htmltools::tagList(lapply(seq_along(covariance_window_sections), function(k) {
              htmltools::tags$div(class = if (k == 1) "carousel-slide active" else "carousel-slide", covariance_window_sections[[k]])
            }))
          )
        } else {
          htmltools::tags$p("No inverse-Wishart covariance blocks were detected.")
        }
      ),
      
      htmltools::tags$details(
        htmltools::tags$summary(class = "chapter-summary", "9. Advanced convergence diagnostics by window"),
        htmltools::tags$p(
          class = "small",
          "Each window includes Stan parameter summaries (grouped by random-effect variances, fixed effects, and residual error variances), HMC diagnostics, Geweke diagnostics, per-chain trace and autocorrelation plots (grouped by scalar variance components, covariance blocks, and residual components), and pairs plots."
        ),
        diagnostics_sections
      ),
      
      htmltools::tags$details(
        htmltools::tags$summary(class = "chapter-summary", "10. Variable catalog"),
        htmltools::tags$p(
          class = "small",
          "Explanation of key variables, column names, and diagnostics used throughout this report."
        ),
        as_table_tag(variable_catalog_df)
      )
    )
  )
  
  htmltools::save_html(report_doc, file = output_file, background = "white")
  invisible(normalizePath(output_file, mustWork = FALSE))
}

render_fit_general_multiple_cycles_report <- function(
    fit_object,
    output_file = "bayesian_multiple_cycles_report.html",
    reference_values = NULL,
    report_title = "Bayesian mixed-model multi-window summary",
    residual_trace_max = 12,
    autocorr_lag_max = 30,
    matrix_entry_max = 15,
    pairs_max = 6,
    digits = 4,
    base_family = "serif",
    image_width = 1800,
    image_height = 1100,
    image_res = 170
) {
  report_spec <- build_fit_general_multiple_cycles_report_spec(
    fit_object = fit_object,
    reference_values = reference_values,
    residual_trace_max = residual_trace_max,
    matrix_entry_max = matrix_entry_max,
    pairs_max = pairs_max,
    digits = digits
  )
  
  render_bayesian_cycle_report_core(
    report_spec = report_spec,
    output_file = output_file,
    report_title = report_title,
    base_family = base_family,
    image_width = image_width,
    image_height = image_height,
    image_res = image_res,
    autocorr_lag_max = autocorr_lag_max
  )
}
