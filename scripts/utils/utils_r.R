if (!require("pacman")) install.packages("pacman")
pkgs<-c("data.table", "terra","MASS"
)
pacman::p_load(char=pkgs,install=T)
rm(pkgs)

pred_gam_newd <- function(mod, varx = "x", vary = NA, var_cat = NA, 
                          extra_x = c(0.1, 0.1), extra_y = c(0.1, 0.1), 
                          n_gridlines = 25, too_far_val = 0.2, se=FALSE,
                          exclude = NULL, type = "response", newdata_check = TRUE) {
  
  data_mod <- data.table(mod$model)
  
  if (is.na(vary)) {
    # 1D case
    x_seq <- seq(
      min(data_mod[[varx]]) - diff(range(data_mod[[varx]])) * extra_x[1],
      max(data_mod[[varx]]) + diff(range(data_mod[[varx]])) * extra_x[2],
      length.out = n_gridlines
    )
    
    if (!is.na(var_cat[1])) {
      # Expand grid over x and all combinations of categorical variables
      cat_levels <- lapply(var_cat, \(v) unique(data_mod[[v]]))
      vars <- do.call(expand.grid, c(list(x_seq), cat_levels))
      setnames(vars, c(varx, var_cat))
    } else {
      vars <- data.table(x_seq)
      setnames(vars, varx)
    }
    
  } else {
    # 2D case
    x_seq <- seq(
      min(data_mod[[varx]]) - diff(range(data_mod[[varx]])) * extra_x[1],
      max(data_mod[[varx]]) + diff(range(data_mod[[varx]])) * extra_x[2],
      length.out = n_gridlines
    )
    y_seq <- seq(
      min(data_mod[[vary]]) - diff(range(data_mod[[vary]])) * extra_y[1],
      max(data_mod[[vary]]) + diff(range(data_mod[[vary]])) * extra_y[2],
      length.out = n_gridlines
    )
    
    if (!is.na(var_cat[1])) {
      cat_levels <- lapply(var_cat, \(v) unique(data_mod[[v]]))
      vars <- do.call(expand.grid, c(list(x_seq, y_seq), cat_levels))
      setnames(vars, c(varx, vary, var_cat))
    } else {
      vars <- expand.grid(x_seq, y_seq)
      setnames(vars, c(varx, vary))
    }
  }
  
  data_new <- data.table(vars)
  
  pred = predict.gam(mod, data_new, type = type, exclude = exclude, newdata.guaranteed = !newdata_check, se.fit=se)
  if(se==TRUE){
    fit=pred$fit
    se=pred$se.fit
    pred = data.table(fit = as.numeric(fit),se  = as.numeric(se)
    )
  } else {
    fit=pred
    pred = data.table(fit = as.numeric(fit))
  }

  if (type == "terms") {
    pred_names <- colnames(pred)
    inv_link <- mod$family$linkinv
    data_new[, (pred_names) := pred[, lapply(.SD, inv_link), ]]
    if (!is.na(var_cat)) {
      coalesce_terms(data_new, var_cat, pred_names)
    }
  } else {
    data_new = cbind(data_new, pred)
  }
  
  if (!is.na(vary)) {
    # Only use `exclude.too.far()` when there's a 2D surface
    if (!is.na(var_cat)) {
      for (i in unique(data_mod[[var_cat]])) {
        data_new[get(var_cat) == i, too_far := exclude.too.far(
          get(varx), get(vary),
          data_mod[get(var_cat) == i, get(varx)],
          data_mod[get(var_cat) == i, get(vary)],
          too_far_val
        )]
      }
    } else {
      data_new[, too_far := exclude.too.far(
        get(varx), get(vary),
        data_mod[[varx]], data_mod[[vary]],
        too_far_val
      )]
    }
  } else {
    # In 1D, too_far doesn't apply
    data_new[, too_far := FALSE]
  }
  
  return(data_new)
}


kde_2d <- function(x, y, levels = c(0.25, 0.5, 0.75), n_grid = 100, 
                   bw="auto", bw_adjust = 1, expand = 0.05) {
  valid_idx <- !is.na(x) & !is.na(y)
  x <- x[valid_idx]
  y <- y[valid_idx]
  
  # Expanded ranges
  x_range <- range(x)
  y_range <- range(y)
  x_margin <- diff(x_range) * expand
  y_margin <- diff(y_range) * expand
  xlim <- c(x_range[1] - x_margin, x_range[2] + x_margin)
  ylim <- c(y_range[1] - y_margin, y_range[2] + y_margin)
  
  ## kernel density estimation
  if (isTRUE(identical(bw, "auto"))) {
    h <- MASS::bandwidth.nrd(c(x, y)) * bw_adjust
    kde_result <- kde2d(x, y, n = n_grid, h = h, lims = c(xlim, ylim))
  } else if (is.na(bw) || is.null(bw)) {
    kde_result <- kde2d(x, y, n = n_grid, lims = c(xlim, ylim))
  } else {
    h <- bw
    kde_result <- kde2d(x, y, n = n_grid, lims = c(xlim, ylim))
  }

  density_values <- kde_result$z
  sorted_density <- sort(c(density_values), decreasing = TRUE)
  cumulative_density <- cumsum(sorted_density) / sum(sorted_density)
  
  level_thresholds <- sapply(levels, function(level) {
    sorted_density[which.min(abs(cumulative_density - level))]
  })
  
  contour_list <- lapply(seq_along(level_thresholds), function(i) {
    threshold <- level_thresholds[i]
    contour_segments <- contourLines(kde_result$x, kde_result$y, kde_result$z, levels = threshold)
    
    rbindlist(lapply(seq_along(contour_segments), function(j) {
      segment <- contour_segments[[j]]
      data.table(x = segment$x, y = segment$y, level = levels[i], group = j)
    }), use.names = TRUE, fill = TRUE)
  })
  
  contour_dt <- rbindlist(contour_list, use.names = TRUE, fill = TRUE, idcol = "level_idx")
  return(contour_dt)
}


kde_grid_2d <- function(
    x, y, n_grid = 100, bw_adjust = 1, expand = 0.1
) {
  ok <- !is.na(x) & !is.na(y)
  x <- x[ok]; y <- y[ok]
  
  # bandwidth (scalar, then diagonal 2x2 for ks::kde when needed)
  h  <- MASS::bandwidth.nrd(c(x, y)) * bw_adjust
  H  <- diag(rep(h, 2))
  
  if (is.numeric(n_grid) && length(n_grid) == 1) {
    xr <- range(x); yr <- range(y)
    xm <- diff(xr) * expand; ym <- diff(yr) * expand
    xlim <- c(xr[1] - xm, xr[2] + xm)
    ylim <- c(yr[1] - ym, yr[2] + ym)
    
    k <- MASS::kde2d(x, y, n = n_grid, h = h, lims = c(xlim, ylim))
    
    grid <- data.table(expand.grid(x = k$x, y = k$y))
    grid[, density := as.vector(k$z)]
    dx <- mean(diff(k$x)); dy <- mean(diff(k$y))
    grid[, pmass := density * dx * dy]
    return(grid[])
  }
  
  # Accept list(x=..., y=...) or list(seq_x, seq_y)
  if (is.list(n_grid)) {
    x_seq <- if (!is.null(n_grid$x)) n_grid$x else n_grid[[1]]
    y_seq <- if (!is.null(n_grid$y)) n_grid$y else n_grid[[2]]
    
    kd <- ks::kde(x = cbind(x, y), H = H, eval.points = list(x_seq, y_seq))
    
    grid <- CJ(x = x_seq, y = y_seq)
    grid[, density := as.vector(kd$estimate)]
    dx <- mean(diff(x_seq)); dy <- mean(diff(y_seq))
    grid[, pmass := density * dx * dy]
    return(grid[])
  }
  
  stop("n_grid must be either a single number or a list of two numeric vectors (x- and y-sequences).")
}

cosine_sim <- function(a, b) {
  sim <- sum(a * b) / (sqrt(sum(a^2)) * sqrt(sum(b^2)))
  (sim + 1) / 2
}

mean_cosine_sim <- function(mat) {
  if (nrow(mat) < 2) return(NA_real_)
  norm_mat <- mat / sqrt(rowSums(mat^2))
  sim_matrix <- norm_mat %*% t(norm_mat)
  mean_sim = mean((sim_matrix[lower.tri(sim_matrix)] + 1) / 2)
}

mahalanobis_dist <- function(row, mean_vec, cov_mat) {
  delta <- as.numeric(row) - mean_vec
  sqrt(as.numeric(t(delta) %*% solve(cov_mat) %*% delta))
}

mirror_plot <- function(p) {
  g <- ggplotGrob(p)
  # find all panel grobs
  panels <- grep("^panel", g$layout$name)
  # reverse x scale inside each panel viewport
  for (i in panels) {
    g$grobs[[i]] <- editGrob(g$grobs[[i]],
                             vp = viewport(xscale = c(1, 0)))  # mirror horizontally
  }
  ggplotify::as.ggplot(g)
}
