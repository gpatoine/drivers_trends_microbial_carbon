# date: 2022-06-15
# author: Guillaume Patoine <guillaume.patoine@idiv.de>
# This script is meant to be sourced. It loads all the packages needed (which
# should first be installed) and functions used in the analysis. In some cases, 
# Roxygen comments are used to document the functions.


# ************* tmst.R *************
# create timestamp
# taken from package gpatoine/gptools

#' Timestamp
#'
#' Squishes, especially useful for file names
#'
#' @param ext character extension
#' @param time logical Should time (HMS) also be included?
#' @param prefix character to be added before, defaults to "_c"
#'
#' @return character squished timestamp
#' @export
tmst <- function(ext = NULL, time = T, prefix = "_c") {
  
  if (!is.null(ext)) {
    
    if(!stringi::stri_sub(ext,1,1) == ".") {
      ext <- paste0(".", ext)
    }
    
  }
  
  if (time) {
    paste0(prefix, format(Sys.time(), "%Y-%m-%d_%H%M%S"), ext)
  } else {
    paste0(prefix, format(Sys.time(), "%Y-%m-%d"), ext)
  }
}

#' Last timestamped
#'
#' Read last timestamped (RDS) file
#' Default is to load file, but can return only name
#'
#' @param fold folder path
#' @param pattern regex pattern passed to list.files
#' @param load logical wanna load to file or just get it's name. use FALSE if the file is not rds format
#' @param prev int previous version before last
#'
#' @return R object read from RDS file
#' @export
last_tmst <- function(fold, pattern = "", load = TRUE, prev = 0) {
  files <- list.files(fold, pattern = pattern, full.names = TRUE)
  file <- files %>% sort(TRUE) %>% .[1 + prev]
  
  if (load & tools::file_ext(file) == "rds") {
    message("Reading ", basename(file))
    readRDS(file)
    
  }  else file
  
}


# Mahalanobis distance ----------------------------------------------------

#' Mahalanobis distance
#'
#' @param dat data.frame
#' @param world world grid
#' @param vars variables used
#' @param xy columns names
#'
#' @return
mahadist <- function (dat, world, vars, xy = c("X", "Y")){
  #check names in both df, stop if not
  
  stopifnot(all(c(vars %in% names(dat), vars %in% names(world))))
  
  dat_sub <- dat %>% dplyr::select(all_of(vars)) %>%
    tidyr::drop_na()
  
  
  dat_diff_nrow <- nrow(dat) - nrow(dat_sub)
  if (dat_diff_nrow > 0) {
    message(dat_diff_nrow, " dat entries with NA values removed")
  }
  
  world_sub <- world %>% dplyr::select(all_of(xy), all_of(vars)) %>%
    tidyr::drop_na()
  
  wrld_diff_nrow <- nrow(world) - nrow(world_sub)
  if (wrld_diff_nrow > 0) {
    message(wrld_diff_nrow, " dat entries with NA values removed")
  }
  
  world_calc <- world_sub %>% dplyr::select(-c(all_of(xy)))
  
  mu <- colMeans(dat_sub) #vector of means
  sigma <- cov(dat_sub) #covariance matrix
  limit97 <- qchisq(.975, df = length(dat_sub))
  limit50 <- qchisq(.5, df = length(dat_sub))
  
  mahaDist <- mahalanobis(world_calc, mu, sigma)
  
  world_calc <- world_calc %>%
    mutate(mahaDistance = mahaDist,
           mahatype = case_when(
             is.na(mahaDist) ~ NA_character_,
             mahaDistance < limit50 ~ "ok",
             mahaDistance < limit97 ~ "chisq > 0.5",
             TRUE ~ "chisq > 0.975"
           ))
  
  world_calc <- world_calc %>%
    mutate(mahatype = factor(world_calc$mahatype,
                             levels = c("ok", "chisq > 0.5", "chisq > 0.975")))
  
  outliers <- length(which(world_calc$mahaDistance > limit97))
  # world_calc %>% filter(mahatype == "chisq > 0.975") %>% nrow #check
  
  message(paste0(outliers, " outliers at 97.5% limit (", round(outliers/nrow(world_calc)*100, 2), "%)"))
  
  world_calc <- bind_cols(world_sub %>% dplyr::select(all_of(xy)), world_calc)
  
  return(world_calc)
  
}



# plot raster layer with ggplot ------------------------------------------------

gp_gplot <- function(x, maxpixels = 5e+4, filt_val = NULL) {
  
  x <- raster::sampleRegular(x, maxpixels, asRaster = TRUE)
  
  coords <- xyFromCell(x, seq_len(ncell(x)))
  dat <- utils::stack(as.data.frame(values(x)))
  names(dat) <- c('value', 'variable')
  dat <- cbind(coords, dat)
  # dat$value %>% unique
  
  if (!is.null(filt_val)) {
    dat <- dat %>% filter(value == filt_val)
  }
  
  ggplot2::ggplot(data=dat, ggplot2::aes(x = x, y = y))+ #, ...
    geom_raster(aes(fill = value))+
    scale_fill_viridis_c(na.value = NA)+
    coord_fixed()
}



# ********** df_slope_calculation.R **********

#' df_slope_lm
#'
#' Calculate slope using lm function. Slower but provides p-values.
#'
#' @param df data
#' @param min_val minimum number of datapoints needed to calculate slope. Otherwise NAs returned.
#'
#' @return vector with slope, p-value and number of entries
#' @export
#'
#' @examples
#' pred <- pred %>% mutate(slope_calc = map(data, ~ df_slope_ma(.x)))
df_slope_lm <- function(df, min_val = 10){

  nent <- nrow(df)  #!is.na(df$stock)

  if (nent < min_val){
    slope <- pval <- NA

  } else { #lm

    mod <- lm(stock ~ year, data = df)
    slope <- coef(mod)[2]
    pval <- summary(mod)$coefficients %>% .[length(.)]

  }

  c(slope = slope, pval = pval, n = nent) #use list if n should be integer

}


#' df_slope_ma
#'
#' Calculate slope using matrix algebra. Faster but no p-values.
#'
#' @param df data
#' @param min_val minimum number of datapoints needed to calculate slope. Otherwise NAs returned.
#'
#' @return vector with slope, p-value and number of entries
#' @export
#'
#' @examples
#' pred <- pred %>% mutate(slope_calc = map(data, ~ df_slope_ma(.x)))
df_slope_ma <- function(df, min_val = 10){

  nent <- nrow(df)  #!is.na(df$pred)

  if (nent == 22) { #fast OLS
    X <- cbind(1, 1:nrow(df))
    invXtX <- solve(t(X) %*% X) %*% t(X)
    slope <- (invXtX %*% df$pred)[2]

  } else  if (nent < min_val){
    slope <- NA

  } else { #longer GLS

    Xf <- model.matrix(~ df$year)
    idx   <- !is.na(df$pred)
    Xsvd  <- svd(Xf[idx,])
    Xplus <- tcrossprod(Xsvd$v %*% diag(Xsvd$d^(-2)) %*% t(Xsvd$v), Xf[idx, ])
    # list(coefs=(Xplus %*% y[idx]))
    slope <- (Xplus %*% df$pred[idx])[2]

  }

  c(slope, nent)

}



# ********** ggplot_helpers.R **********

#' draw_key_polygon3
#'
#' To adjust legend with spaces. The function is called internally by ggplot when key_glyph = "polygon3" is specified in the geom call.
#' To use with factorial raster map (binned/cut).
#'
#' reference:
#' https://stackoverflow.com/questions/11366964/is-there-a-way-to-change-the-spacing-between-legend-items-in-ggplot2
#' https://github.com/tidyverse/ggplot2/issues/3180
#'
#' @param data data
#' @param params params
#' @param size size
#'
#' @return dunno
#' @export
#'
#' @examples geom_raster(key_glyph = "polygon3")
draw_key_polygon3 <- function(data, params, size) {
  lwd <- min(data$size, min(size) / 4)

  grid::rectGrob(
    width = grid::unit(0.8, "npc"),
    height = grid::unit(0.8, "npc"),
    gp = grid::gpar(
      col = NA,
      fill = alpha(data$fill, data$alpha),
      lty = data$linetype,
      lwd = lwd * .pt,
      linejoin = "mitre"
    ))
}



# ********** load_dataset.R **********

#' Load processed dataset
#'
#' Should load the most recent cleaned up dataset.
#'
#' @return tibble
#' @export
glc_proc_data <- function() {

  # readRDS(here("rawdata", "glc_proc_data.rds"))
  
  readxl::read_xlsx(here("rawdata/glc_cmic_data.xlsx"), sheet = "data") %>% 
    mutate(land_cover = factor(land_cover, levels = c("C", "FB", "FC", "FT", "G", "S")))

}


#' glc_rf_model
#'
#' Load most recent random forest model
#'
#' @return caret train object
#' @export
glc_rf_model <- function() {
  readRDS(here("derived", "06-1-rf_model.rds"))
}


#' glc_fullmak_df
#'
#' Load fullmask
#'
#' @return tibble
#' @export
glc_fullmask_df <- function() {
  readRDS(here("derived", "06-3-fullmask_df.rds"))
}


#' glc_predicted_dataset
#'
#' All predictions, all years, with stocks calculated
#'
#' @return dataset
#' @export
glc_predicted_dataset <- function() {
  readRDS(here("derived", "07-cmic_stocks_dataset_all_years.rds"))
  
}



# ********** misc.R **********

#' glc_land_cover_classes
#'
#' @return tibble with land_cover values used
#' @export
glc_land_cover_classes <- function(){
  tibble::tribble(
     ~code, ~value,              ~class,
       "C",    10L,          "cropland",
      "FT",    50L,   "tropical forest",
      "FB",    60L,  "broadleaf forest",
      "FC",    70L, "coniferous forest",
       "S",   120L,         "shrubland",
       "G",   130L,         "grassland",
       "B",   150L,      "sparse, bare"
     )
}


#' LC numeric to character
#'
#' @param x numeric vector of land cover codes
#'
#' @return character vector as factor
#' @export
glc_LC_num2chr <- function(x) {
  tibble(value = x) %>%
    left_join(glc_land_cover_classes()[1:2], by = "value") %>%
    select(land_cover = code) %>%
    mutate(land_cover = as.factor(land_cover)) %>%
    pull
  
}


#' LC character to numeric
#'
#' @param x character vector of land cover values
#'
#' @return numeric vector
#' @export
glc_LC_chr2num <- function(x) {
  tibble(code = x) %>%
    left_join(glc_land_cover_classes()[1:2], by = "code") %>%
    select(land_cover = value) %>%
    # mutate(land_cover = as.factor(land_cover)) %>%
    pull
  
}


#' makeQuantiles
#'
#' convenience function to make quantiles with cut
#'
#' @param x vectorof values
#' @param probs sections
#'
#' @return
#' @export
makeQuantiles <- function(x, probs = seq(0, 1, by = 0.125)) {
  cut(x, breaks = quantile(x, probs = probs), include.lowest = TRUE)
}



# ********** pIDfromXY.R **********

#' pIDfromXY
#'
#' Returns the cell number (point ID) from resampled raster layer based on XY coordinates.
#'
#'
#' @param x numeric
#' @param y numeric
#' @param csiz numeric Default to 0.05
#'
#' @return numeric
#' @export
#'
#' @examples df %>% mutate(pid = glc_pIDfromXY(x,y))
glc_pIDfromXY <- function(x, y, csiz = 0.05) {

  assertthat::assert_that(all(between(x, -180, 180)),
                          all(between(y, -90, 90)))

  x <- round(x, 4)
  y <- round(y, 4)

  nccol <- 2*180/csiz
  ncrow <- 2*90/csiz

  cellcol <- (x + 180) / csiz + 0.5
  cellrow <- (90 - y) / csiz + 0.5

  pID <- (cellrow-1) * nccol + cellcol
  pID %>% round %>% as.integer

}


#' XY from pID
#'
#' @param pid point ID number, equivalent to raster cell number
#' @param csiz cell size. Default to 0.05
#'
#' @return matrix of x and y values
#' @export
#'
#' @examples df %>% bind_cols(glc_XYfrompID(pid))
glc_XYfrompID <- function(pid, csiz = 0.05) {
  pid <- pid %>% as.numeric %>% round %>% as.integer

  nccol <- 2*180/csiz
  ncrow <- 2*90/csiz

  cellrow <- pid %/% nccol + 1
  cellcol <- pid %% nccol

  y <- 90 - (cellrow - 0.5) * csiz
  x <- -180 + (cellcol - 0.5) * csiz

  tibble(x = x, y = y)

}


#' X from pID
#'
#' @inheritParams glc_XYfrompID
#'
#' @return point ID
#' @export
#'
#' @examples df %>% mutate(x = glc_XfrompID(pid),
#' y = glc_YfrompID(pid))
glc_XfrompID <- function(pid, csiz = 0.05) {
  XYfrompID(pid, csiz)[,1]
  
}


#' Y from pID
#'
#' @inheritParams glc_XYfrompID
#'
#' @return point ID
#' @export
#'
#' @examples df %>% mutate(x = glc_XfrompID(pid),
#' y = glc_YfrompID(pid))
glc_YfrompID <- function(pid, csiz = 0.05) {
  XYfrompID(pid, csiz)[,2]
  
}



# ********** retrieve_data_layers.R **********

#' All layers
#'
#' List all layers
#'
#' @return character vector
#' @export
glc_layers <- function() {
  c(static = c("elev", "clay", "nitrogen", "phh2o", "sand", "soc"),
    dynamic = c("land_cover", "tmean", "prec", "ndvi"))
  
}



#' glc_dynrange
#'
#' Prints temporal range of dynamic layers
#'
#' @return tibble
#' @export
glc_dynrange <- function() {
  tibble::tribble(
           ~layer, ~from,  ~to,
     "land_cover",  1992, 2015,
          "tmean",  1979, 2013,
           "prec",  1979, 2013,
           "ndvi",  1981, 2016
     )
}


#' glc_get_raster
#'
#' Find the path and reads raster file.
#' Available layers are "elev", "phh2o", "soc", "sand", "clay", "nitrogen",
#' "tmean", "prec", "ndvi", "land_cover", "latitude".
#'
#' @param layer character One of "elev", "phh2o", "soc", "sand", "clay", "nitrogen",
#' "tmean", "prec", "ndvi", "land_cover", or "latitude".
#' @param year integer Year value from 1992 to 2013.
#' @param depth character One of the three SoilGrids depth range: "0-5", "5-15" or "15-30".
#'
#' @return RasterLayer
#' @export
glc_get_raster <- function(layer, year = NULL, depth = "5-15") {

  slgrds <- c("phh2o", "soc", "sand", "clay", "nitrogen")


  # static ------------------------------------------------------------------

  if (layer == "elev") {
    rpath <- here("geodata/static/WorldClim/wc2.1_30s_elev/wc2.1_30s_elev.tif")
    ras <- raster::raster(rpath)

  } else if (layer %in% slgrds) {
    rpath <- here("geodata/static/SoilGrids/tif_files", paste0(layer, "_", depth, "cm_mean.tif"))
    ras <- raster::raster(rpath)


    # dynamic -----------------------------------------------------------------

  } else if (layer == "tmean") {
    rpath <- here("geodata/dynamic/CHELSA/year_tmean",
                       paste0("CHELSA_", layer, "_", year, "_mean.tif"))
    ras <- raster::raster(rpath)

  } else if (layer == "prec") {
    rpath <- list.files(here("geodata/dynamic/CHELSA/year_prec"),
                        paste0("prec_", year, "_sum"),
                        full.names = TRUE)
    ras <- raster::raster(rpath)

  } else if (layer == "ndvi") {
    rpath <- list.files(here("geodata/dynamic/NDVI/yearly/yearly_mean/geotiff"),
                        paste0("NDVI_", year, "_mean_c"),
                        full.names = TRUE)
    ras <- raster::raster(rpath)

  } else if (layer == "land_cover") {
    # rpath <- here("geodata/dynamic/ESA/ESACCI-LC-L4-LCCS-Map-300m-P1Y-1992_2015-v2.0.7.tif") # need to use reclassified
    rpath <- list.files(here("geodata/dynamic/ESA/reclass"),
                        paste0("ESA_reclassified_", year, "_c20"),
                        full.names = TRUE)

    ras <- raster::raster(rpath)

  } else {
    stop("Not a valid layer")

  }

  ras

}



#' glc_get_resamp
#'
#' Get resampled global raster
#'
#' @inheritParams glc_get_raster
#'
#' @return RasterLayer
#' @export
glc_get_resamp <- function(layer, year = NULL, depth = "5-15") {

  slgrds <- c("phh2o", "soc", "sand", "clay", "nitrogen")

  chls <- c("tmean", "prec")

  resamp_path <- function(...) here("geodata/resampled_0p05deg", ...)


  # static ------------------------------------------------------------------

  if (layer == "elev") {
    rpath <- resamp_path("static/elev") %>% last_tmst("elev_c", load = F)

  } else if (layer == "latitude") {
    rpath <- resamp_path("static/latitude") %>% last_tmst("latitude_from_ndvi_c", load = F)

  } else if (layer %in% slgrds) {
    rpath <- resamp_path("static", layer) %>% last_tmst(paste0(layer, "_", depth, "cm_mean_resamp_c"), load = F)


    # dynamic -----------------------------------------------------------------

  } else if (layer %in% chls) {
    rpath <- resamp_path("dynamic", layer) %>% last_tmst(paste0(layer, "_", year, "_c"), load = F)

  } else if (layer == "ndvi") {
    rpath <- resamp_path("dynamic", layer) %>% last_tmst(paste0("NDVI_", year, "_mean_c"), load = F)

  } else if (layer == "land_cover") {
    rpath <- resamp_path("dynamic", layer) %>% last_tmst(paste0(layer, "_", year, "_c"), load = F)

  } else {
    stop("Not a valid layer")

  }

  stopifnot(file.exists(rpath))
  message("Reading ", rpath)

  raster::raster(rpath)

}


# https://rdrr.io/cran/biscale/src/R/bi_legend.R
bi_legend <- function(pal, dim = 3, xlab, ylab, size = 10){
  
  # global binding
  bi_class = bi_fill = x = y = NULL
  
  # check parameters
  if (missing(pal) == TRUE){
    stop("A palette must be specified for the 'pal' argument.")
  }
  
  if ("bi_pal_custom" %in% class(pal) == TRUE) {
    
    if (dim == 2 & length(pal) != 4){
      stop("There is a mismatch between the length of your custom palette object and the given dimensions.")
    } else if (dim == 3 & length(pal) != 9){
      stop("There is a mismatch between the length of your custom palette object and the given dimensions.")
    }
    
  } else if ("bi_pal_custom" %in% class(pal) == FALSE){
    
    if (pal %in% c("Brown", "DkBlue", "DkCyan", "DkViolet", "GrPink") == FALSE){
      stop("The given palette is not one of the allowed options for bivariate mapping. Please choose one of: 'Brown', 'DkBlue', 'DkCyan', 'DkViolet', or 'GrPink'.")
    }
    
  }
  
  if (is.numeric(dim) == FALSE){
    stop("The 'dim' argument only accepts the numeric values '2' or '3'.")
  }
  
  if (dim != 2 & dim != 3){
    stop("The 'dim' argument only accepts the numeric values '2' or '3'.")
  }
  
  if (missing(xlab) == TRUE){
    xlab <- "x var "
  }
  
  if (is.character(xlab) == FALSE){
    stop("The 'xlab' argument must be a character string.")
  }
  
  if (missing(ylab) == TRUE){
    ylab <- "y var "
  }
  
  if (is.character(ylab) == FALSE){
    stop("The 'ylab' argument must be a character string.")
  }
  
  if (is.numeric(size) == FALSE){
    stop("The 'size' argument must be a numeric value.")
  }
  
  # nse
  xQN <- rlang::quo_name(rlang::enquo(xlab))
  yQN <- rlang::quo_name(rlang::enquo(ylab))
  
  # obtain palette
  if ("bi_pal_custom" %in% class(pal) == TRUE) {
    
    x <- pal
    
  } else if ("bi_pal_custom" %in% class(pal) == FALSE){
    
    if (pal == "DkViolet"){
      x <- pal_dkviolet(n = dim)
    } else if (pal == "GrPink"){
      x <- pal_grpink(n = dim)
    } else if (pal == "DkBlue"){
      x <- pal_dkblue(n = dim)
    } else if (pal == "DkCyan"){
      x <- pal_dkcyan(n = dim)
    } else if (pal == "Brown"){
      x <- pal_brown(n = dim)
    }
    
  }
  
  # create tibble for plotting
  x <- dplyr::tibble(
    bi_class = names(x),
    bi_fill = x
  )
  
  # reformat
  leg <- tidyr::separate(x, bi_class, into = c("x", "y"), sep = "-")
  leg <- dplyr::mutate(leg, x = as.integer(x), y = as.integer(y))
  
  # create ggplot2 legend object
  legend <- ggplot2::ggplot() +
    ggplot2::geom_tile(data = leg, mapping = ggplot2::aes(x = x, y = y, fill = bi_fill)) +
    ggplot2::scale_fill_identity() +
    ggplot2::labs(x = substitute(paste(xQN, ""%->%"")), y = substitute(paste(yQN))) +
    bi_theme() +
    ggplot2::theme(axis.title = ggplot2::element_text(size = size)) +
    ggplot2::coord_fixed()
  
  # return output
  return(legend)
  
}
