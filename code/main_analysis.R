# date: 2022-06-15
# author: Guillaume Patoine <guillaume.patoine@idiv.de>
# description: This is the script related to the publication "Drivers and trends 
# of global soil microbial carbon over two decades".

# The project uses the following folder structure
#
# project_root
# ├─ code 
# ├─ derived 
# ├─ geodata
# ├─ output
# │   ├─ figures
# │   └─ tables
# └─ rawdata

# NOTE some of the functions require additional datasets that are not provided
# in this repository, dur to space limitation. For example, `glc_get_resamp()`
# is a convenience function to load a raster layer, based on the variable name,
# year and soil depth. Code sections based on that function won't work, unless
# these layers are available.


# Load packages -----------------------------------------------------------

library(raster)
library(tools)
library(sf)
library(tidyverse)
library(here)
library(doParallel)
library(foreach)
library(caret)
library(CAST)
library(cowplot)
library(grid)
library(gridExtra)
library(ggExtra)
library(forcats)
library(magick)
library(biscale)
library(gt)

# The packages below need to be installed, but do no need to be loaded.
# The needed functions are called directly using ::
# library(readxl)
# library(fasterize)


# source functions from other script --------------------------
source(here::here("code/functions.R"))

# figure dimensions and default theme
gg_width = 11
gg_height = 5.7
ggplot2::theme_set(ggplot2::theme_bw())


# main dataset -------------------------------------------------

cmic0 <- glc_proc_data()


# ************************************************************
# ---------------- Random forest model (06-1) ----------------
# ************************************************************

predictors <- glc_layers() %>% unname
cmic <- as.data.frame(cmic0) # needed

set.seed(202)
model <- train(x = cmic[,predictors], 
               y = cmic$Cmic,
               method = "rf",
               importance = TRUE,
               tuneGrid = expand.grid(mtry = c(2:4)), # length(predictors) or 2:6
               trControl = trainControl(method = "cv", 
                                        number = 20,
                                        p = 0.75,
                                        savePredictions = TRUE))

model # most often mtry = 2

saveRDS(model, here("derived", "06-1-rf_model.rds"))

# RMSE + R2
model$results %>% as_tibble %>% filter(mtry == model$bestTune %>% unlist) %>% select(RMSE, Rsquared)


# variable importance
varImp(model) %>% plot
varImp(model, scale = FALSE) %>% plot


## variability assessment models -----------------------------------------------
set.seed(202)
models <- map(1:100, ~ {print(.x); train(x = cmic[,predictors], 
                                         y = cmic$Cmic,
                                         method = "rf",
                                         importance = TRUE,
                                         tuneGrid = expand.grid(mtry = c(2:4)),
                                         trControl = trainControl(method = "cv", 
                                                                  number = 20,
                                                                  p = 0.75))}
)

vimp_df <- map_dfr(models, ~ varImp(.x, scale = F)$importance %>% as_tibble(rownames = "variable"))
vimp_df <- vimp_df %>% rename(importance = Overall)

# update names
ren_tib = tibble(variable = c("elev", "clay", "nitrogen", "phh2o", "sand", "soc", "land_cover", 
                              "tmean", "prec", "ndvi"),
                 recode = c("Elevation", "Clay", "Nitrogen", "pH", "Sand", "SOC",
                            "Land cover", "Temperature", "Precipitation", "NDVI"))

vimp_df <- vimp_df %>% left_join(ren_tib, by = "variable") %>% select(-variable, variable = recode)

# plot variable importance
ggplot(vimp_df, aes(importance, y = reorder(variable, importance, FUN = stats::median),
                    fill = reorder(variable, importance, FUN = stats::median))) +
  geom_boxplot() +
  # scale_fill_viridis_d(direction = -1, alpha = 0.8) +
  scale_fill_grey(start = 0.35, end = 0.9) +
  ylab("Variable") +
  xlab("Variable importance") + # Mean decrease in accuracy when permuted OR Mean error increase when permuted
  theme_bw() +
  theme(legend.position = "none")

ggsave(here("output/figures", "06-1-rf_model_varImp.png"), width = 7, height = 6)



# ************************************************************
# ------- Environmental coverage analysis (06-2,3) -----------
# ************************************************************

## Mahalanobis distance --------------------------------------

### world grid -------------------------------------------------------------

# no LC
vars <- glc_layers() %>% .[!. %in% c("land_cover")] %>% unname

# stack all layers
st <- stack(map(vars, ~glc_get_resamp(.x, 2013, "5-15")))
names(st) <- vars

# ~2 min
grid <- as.data.frame(st, xy = TRUE, na.rm = TRUE)
head(grid)

# maha cmic
maha_cmic <- mahadist(dat = cmic, world = grid, vars = vars, xy = c("x", "y"))
mahamask <- maha_cmic %>% filter(mahatype == "chisq > 0.975") %>% select(x, y)


### Figure supp ------------------------------------------------------------

# adjust legend fct 

# https://stackoverflow.com/questions/11366964/is-there-a-way-to-change-the-spacing-between-legend-items-in-ggplot2
# https://github.com/tidyverse/ggplot2/issues/3180

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



# figure
cmic <- mahadist(dat = cmic, world = grid, vars = vars, xy = c("x", "y"))

ggplot(cmic, aes(x, y, fill = mahatype))+
  borders(fill = "grey90", colour = NA, ylim = c(-60, 90))+
  geom_raster(key_glyph = "polygon3")+
  coord_fixed(xlim = c(-180, 180), expand = FALSE)+
  scale_fill_manual(values = c("<0.50" = "#0c1079ff",
                               "0.50-0.75" = "#1b93abff",
                               "0.75-0.975" = "#38df09ff",
                               ">0.975 (Outliers)" = "#fea77fff"),
                    name = "Quantile distribution")+
  theme_void()+
  theme(legend.position = c(0.1, 0.18),
        legend.title = element_text(face = "bold"),
        # legend.key = element_rect(size = 3, fill = "white", colour = NA), # didn't work
        legend.key.width = unit(1.2, "cm")
        # legend.key.height = unit(1, "cm")
  ) +
  annotate("text", -Inf, Inf, hjust = -0.5, vjust = 1.5, label = "A", size = 16/.pt,
           fontface = "bold")

ggsave(here("output/figures", "05-maha_quant_suppFig.png"),
       width = gg_width, height = gg_height)



## Area of applicability ----------------------------------------------------

rmse <- function(pred,obs){sqrt( mean((pred - obs)^2, na.rm = TRUE) )}

# cmic

cmic # as.data.frame needed to avoid issues with tibble in CAST

### RF model ----------------------------------------------------------------

mod <- model

predictors <- mod$trainingData %>% names %>% .[-length(.)]

# prediction layers as stack
raspred <- stack(map(predictors, ~glc_get_resamp(.x, 2013, "5-15")))
# cbind(predictors, names(raspred) %>% word(sep = "_"))
names(raspred) <- predictors


# ~2min
preddf <- as.data.frame(raspred, xy = TRUE, na.rm = TRUE)
all(names(preddf %>% select(-c(x, y))) %in% names(mod$trainingData))

# fix land_cover col
preddf <- preddf %>% 
  mutate(land_cover = glc_LC_num2chr(land_cover)) %>% 
  filter(land_cover %in% levels(mod$trainingData$land_cover)) %>% 
  mutate(land_cover = factor(land_cover)) #remove unused levels

# predict from stack.as.df (~1min)
prediction <- predict(mod, preddf)


### aoa calculation ---------------------------------------------------------

cl <- makeCluster(10) #8-10
registerDoParallel(cl)

# with variable weighting: (~5-10 min)
AOA <- aoa(preddf, model = mod, returnTrainDI = TRUE, cl = cl)
# AOA$AOA %>% table
# attributes(AOA)$aoa_stats

saveRDS(AOA, here("derived", "06-2-AOA_object.rds"))

stopCluster(cl)


## put together ------------------------------------------------------------

preds <- mod$pred[mod$pred$mtry==mod$bestTune$mtry,]

absError <- abs(preds$pred[order(preds$rowIndex)]-preds$obs[order(preds$rowIndex)])

preddf <- preddf %>% mutate(DI = AOA$DI,
                            AOA = AOA$AOA,
                            DI_nw = AOA_noWeights$DI,
                            AOA_nw = AOA_noWeights$AOA)

### figure suppl ------------------------------------------------------------

thr_stats <- attributes(AOA)$aoa_stats$threshold_stats
thresh <- attributes(AOA)$aoa_stats$threshold

pfig <- preddf %>% mutate(di_quant = factor(case_when(
  DI > thresh ~ ">0.95 (Outliers)",
  DI > thresh*(0.75/0.95) ~ "0.75-0.95",
  DI > thresh*(0.5/0.95) ~ "0.50-0.75",
  TRUE ~ "<0.50"), 
  levels = c("<0.50", "0.50-0.75", "0.75-0.95", ">0.95 (Outliers)"))
  )


ggplot(pfig, aes(x, y, fill = di_quant))+
  borders(fill = "grey90", colour = NA, ylim = c(-60, 90))+
  geom_raster(key_glyph = "polygon3")+
  coord_fixed(xlim = c(-180, 180), expand = FALSE)+
  scale_fill_manual(values = c("<0.50" = "#0c1079ff",
                               "0.50-0.75" = "#1b93abff",
                               "0.75-0.95" = "#38df09ff",
                               ">0.95 (Outliers)" = "#fea77fff"),
                    name = "Quantile distribution")+
  theme_void()+
  theme(legend.position = c(0.1, 0.18),
        legend.title = element_text(face = "bold"),
        legend.key.width = unit(1.2, "cm")
  )+
  annotate("text", -Inf, Inf, hjust = -0.5, vjust = 1.5, label = "B", size = 16/.pt,
           fontface = "bold")

ggsave(here("output/figures", "06-2-AOA_quant_suppFig.png"),
       width = gg_width, height = gg_height)


## Create mask

preddf

aoamask <- preddf %>% filter(AOA == 0) %>% select(x, y)
aoamask$aoam <- 1
aoamask <- aoamask %>% mutate(pid = glc_pIDfromXY(x, y)) %>% 
  select(-c(x,y))
sum(!aoa_df$AOA) == nrow(aoamask)


# mahalanobis mask
mahamask
mahamask$maham <- 1
mahamask <- mahamask %>% mutate(pid = glc_pIDfromXY(x, y)) %>% 
  select(-c(x,y))

head(aoamask)
head(mahamask)


# combine masks
fullmask <- full_join(aoamask, mahamask, by = "pid")

fullmask <- fullmask %>% select(pid, maha_mask = maham, aoa_mask = aoam) %>% 
  mutate(mask = 1) # as.numeric(maha_mask | aoa_mask)

saveRDS(fullmask, here("derived", "06-3-fullmask_df.rds"))

head(fullmask)


### plot --------------------------------------------------------------------

preddf <- preddf %>% mutate(pid = glc_pIDfromXY(x, y)) %>% 
  left_join(fullmask, by = "pid")

preddf %>% names

sum(preddf$maha_mask == 1, na.rm = TRUE)
sum(preddf$aoa_mask == 1, na.rm = TRUE)

# Compare both
preddf <- preddf %>% mutate(Category = case_when(
  maha_mask == 1 & aoa_mask == 1 ~ "Outlier with both methods",
  maha_mask == 1 ~ "Mahalanobis outlier",
  aoa_mask == 1 ~ "Area of applicability outlier",
  TRUE ~ "Predicted"
) %>% factor(levels = c("Mahalanobis outlier", "Area of applicability outlier",
                        "Outlier with both methods", "Predicted", "Deficient data")))

preddf %>% count(Category)

cmic <- glc_proc_data()


# need to separate legend
(p_mask <- ggplot(preddf, aes(x, y, fill = Category))+
    borders(size = 0.3, fill = "grey90", ylim = c(-60, 90), colour = NA,
            show.legend = TRUE)+
    geom_raster(alpha = 0.7)+
    # borders(size = 0.3, fill = NA, ylim = c(-60, 90))+
    scale_fill_manual(values = c("Mahalanobis outlier" = "#2a7fff",
                                 "Area of applicability outlier" = "#ff4a4a",
                                 "Outlier with both methods" = "#cc78df",
                                 "Predicted" = "#8ef284",
                                 "Deficient data" = "grey90"),
                      drop = FALSE)+
    geom_point(data = cmic, aes(x = longitude, y = latitude), 
               size = 1, inherit.aes = F, shape = 1, alpha = 0.5)+
    theme_void()+
    theme(
      # legend.position = c(0.1, 0.2),
      legend.title = element_text(face = "bold"))+
    coord_fixed())

p_mask_noleg <- p_mask +
  theme(legend.position = "none")

leg_pt <- get_legend(ggplot(mtcars %>% mutate(col = "Sampling site"), aes(wt, mpg, color = col)) +
                       geom_point(shape = 1)+
                       scale_color_manual(values = "black", name = NULL)+
                       theme(legend.justification = c(-0.011,1)))

ggdraw(leg_pt)

leg1 <- get_legend(p_mask +
                     scale_fill_manual(breaks = c("Predicted"),
                                       values = c("Mahalanobis outlier" = "#2a7fff",
                                                  "Area of applicability outlier" = "#ff4a4a",
                                                  "Outlier with both methods" = "#cc78df",
                                                  "Predicted" = "#8ef284",
                                                  "Deficient data" = "grey90"),
                                       labels = "Within thresholds of both methods",
                                       drop = TRUE)+
                     labs(fill = "Predicted")+
                     theme(legend.justification = c(0,1)))
ggdraw(leg1)

leg2 <- get_legend(p_mask +
                     scale_fill_manual(breaks = c("Mahalanobis outlier", "Area of applicability outlier",
                                                  "Outlier with both methods", "Deficient data"),
                                       values = c("Mahalanobis outlier" = "#2a7fff",
                                                  "Area of applicability outlier" = "#ff4a4a",
                                                  "Outlier with both methods" = "#cc78df",
                                                  "Predicted" = "#8ef284",
                                                  "Deficient data" = "grey90"),
                                       # labels = "Within thresholds of both methods",
                                       drop = FALSE)+
                     labs(fill = "Excluded locations")+
                     theme(legend.justification = c(0,1)))

ggdraw(leg2)

legs_both <- plot_grid(leg1, leg2,
                       ncol = 1, rel_heights = c(1.2, 2.5))

lc_tags <- tibble(lc_code = c("B", "C", "FB", "FC", "FT", "G", "S", "W"),
                  lc_full = c("Bare", "Cropland", "Broadleaf forest", "Coniferous forest",
                              "Tropical forest", "Grassland", "Shrubland", "Wetlands"))

lc_count <- cmic %>% count(land_cover) %>% mutate(land_cover = factor(land_cover)) %>%
  left_join(lc_tags, by = c("land_cover" = "lc_code"))

lc_count %>% pull(lc_full)

fct_ord <- lc_count %>% arrange(desc(n)) %>% pull(lc_full)
lc_count <- lc_count %>% mutate(lc_full = factor(lc_full, levels = fct_ord))

# removed colors for clarity
fill_val <- c("Cropland" = "#ffe699",
              "Tropical forest" = "#e2efda",
              "Broadleaf forest" = "#a9d08e",
              "Coniferous forest" = "#548235",
              "Shrubland" = "#f8cbad",
              "Grassland" = "#f4b084")

(lcb1 <- ggplot(lc_count, aes(n, lc_full, fill = lc_full)) +
    geom_col(orientation = "y", show.legend = FALSE, fill = "grey20") +
    scale_fill_manual(values = fill_val)+
    xlab("Number of samples")+
    theme_classic() +
    theme(axis.title.y = element_blank(),
          plot.background = element_blank(),
          panel.background = element_blank()))

(finalPlot <- ggdraw(p_mask_noleg) +
    draw_plot(lcb1, 0.58, 0.1, 0.24, 0.24)+
    draw_plot(leg_pt, 0.055, 0.42, 0.2, 0.06, hjust = 0)+
    draw_plot(legs_both, 0.06, 0.12, 0.2, 0.3, hjust = 0, halign = -1))


ggsave(plot = finalPlot,
       here("output/figures", "06-3-Figure_1_masks.png"),
       width = gg_width, height = gg_height)



# ************************************************************
# ----------------- Prediction SD (06-7) ---------------------
# ************************************************************

# Prediction intervals for 2013 with 100 runs

# data used
cmic <- glc_proc_data() %>% as.data.frame
predictors
lc_levs <- levels(cmic$land_cover)

# global layers
st <- stack(map(predictors, ~glc_get_resamp(.x, 2013, "5-15")))
names(st) <- predictors

gridyear <- as.data.frame(st, xy = TRUE, na.rm = TRUE)
gridyear <- gridyear %>% mutate(pid = glc_pIDfromXY(x, y)) # similar to cellFromXY()

# LC to text
gridyear <- gridyear %>% 
  mutate(land_cover = glc_LC_num2chr(land_cover)) %>% 
  filter(land_cover %in% lc_levs) %>%
  mutate(land_cover = factor(land_cover))

# mask
fullmask <- fullmask %>% select(pid, mask)

gridyear <- gridyear %>% left_join(fullmask, by = "pid") %>% 
  filter(is.na(mask)) %>% select(-mask)


# fct to run RF
make_pred_2013 <- function() {
  
  cat("\nstarting")
  
  model <- train(x = cmic[,predictors], 
                 y = cmic$Cmic,
                 method = "rf",
                 importance = TRUE,
                 tuneGrid = expand.grid(mtry = c(2:4)),
                 trControl = trainControl(method = "cv", 
                                          number = 20,
                                          p = 0.75,
                                          savePredictions = FALSE))
  
  cat("\n    training done")
  
  # predict (100 sec)
  gridpred <- gridyear %>% mutate(pred = predict(model, gridyear))
  
  cat("\n        predict done")
  
  stmst <- tmst(ext = paste0("_", stringi::stri_rand_strings(1,4)))
  
  # save raster
  df_ras <- gridpred %>% select(x, y, pred)
  p_ras <- rasterFromXYZ(df_ras)
  crs(p_ras) <- crs(st)
  p_ras2 <- extend(p_ras, st)
  
  writeRaster(p_ras2, here("derived", "06-7-RF_runs",
                                paste0("One_run_2013", stmst, ".tif")))
  
  cat("\n            full run done")
  
}

# Requires a multi-core cluster
cl <- makeCluster(6, outfile = "")
registerDoParallel(cl)

foreach(ii = 1:100, .packages = c("caret")) %dopar% {
  
  source(here::here("code/glc_functions.R"))
  make_pred_2013()

}

stopCluster(cl)


## calculate stocks --------------------------------------------------------------

st <- stack(list.files(here("derived/06-7-RF_runs"), 
                       pattern = "\\.tif$", full.names = TRUE))

# calculate stocks
# stock = pred * bulk * (1000-cofr) * 12.01 / 10^8
bulk <- raster(here("geodata/resampled_0p05deg/static/bdod",
                         "bdod_5-15cm_mean_resamp.tif"))
cofr <- raster(here("geodata/resampled_0p05deg/static/cfvo",
                         "cfvo_5-15cm_mean_resamp.tif"))

bulk_cofr <- bulk * (1000-cofr) * 12.01 / 10^8

# stock is in tonnes/ha (weight/area)
st_stock <- st * bulk_cofr


## SD and plot -------------------------------------------------------------

pred_sd <- calc(st_stock, sd, na.rm = TRUE, progress = "text",
                filename = here("derived", "06-7-pred_sd_2013_100runs.tif"))

# is sd and mean value strongly correlated? yea, a bit
pred_mean <- calc(st_stock, mean, na.rm = TRUE, progress = "text")

st2 <- stack(pred_mean, pred_sd)
st_df <- as.data.frame(st2, xy = TRUE, na.rm = TRUE)
st_df <- st_df %>% rename(val = 3, sd = 4)


## Finish plot SD ----------------------------------------------------------

# better: SD divided by mean
rel_sd <- pred_sd / pred_mean

gp_gplot(rel_sd, 5e+5)+ #Inf
  borders(fill = "grey90", colour = NA, ylim = c(-60, 90))+
  coord_fixed(xlim = c(-180, 180), ylim = c(-65, NA), expand = FALSE)+
  geom_raster(aes(fill = value))+
  scale_fill_distiller(na.value = NA, palette = "YlOrRd", direction = 1,
                       name = "Coefficient of variation") +
  theme_void()+
  theme(legend.position = c(0.14, 0.22),
        legend.title = element_text(face = "bold"),
        # legend.key = element_rect(size = 3, fill = "white", colour = NA), # didn't work
        legend.key.width = unit(0.8, "cm")
        # legend.key.height = unit(1, "cm")
  ) +
  annotate("text", -Inf, Inf, hjust = -0.5, vjust = 1.5, label = "A", size = 16/.pt,
           fontface = "bold")

ggsave(here("output/figures", "06-7-prediction_cofVar_suppFig.png"),
       width = gg_width, height = gg_height)



# ************************************************************
# ---------------- Temporal predictions (06-5) ---------------
# ************************************************************

# create prediction files
predictors <- model$trainingData %>% names %>% .[-length(.)]

lc_levs <- levels(model$trainingData$land_cover)

# Predictions all years ---------------------------------------------------
# predict, save df

cl <- makeCluster(6, outfile = "")
registerDoParallel(cl)

# 16 min with 6 cores, without saving raster
foreach (iyear = 1992:2013, .packages = c("caret")) %dopar% {
  
  source(here::here("code/functions.R"))
  
  # global layers
  st <- stack(map(predictors, ~glc_get_resamp(.x, iyear, "5-15")))
  
  # rename
  names(st) <- predictors
  
  # as.data.frame (20 sec)
  gridyear <- as.data.frame(st, xy = TRUE, na.rm = TRUE)
  gridyear <- gridyear %>% mutate(pid = glc_pIDfromXY(x, y)) # pid = rownames(gridyear) or cellFromXY()
  
  # land_cover to text
  gridyear <- gridyear %>% 
    mutate(land_cover = glc_LC_num2chr(land_cover)) %>% 
    filter(land_cover %in% lc_levs) %>%
    mutate(land_cover = factor(land_cover))
  
  # mask
  fullmask <- glc_fullmask_df() %>% select(pid, mask)
  
  gridyear <- gridyear %>% left_join(fullmask, by = "pid") %>% 
    filter(is.na(mask)) %>% select(-mask)
  
  
  # predict (100 sec)
  gridyear$pred <- predict(model, gridyear)
  
  saveRDS(gridyear, here("derived", "06-prediction_years",
                              paste0("06-predictions_cmic_", iyear, ".rds")))
  
}

stopCluster(cl)



# ************************************************************
# ------ Partial dependence and response curves (07-1 + 11-1) ------
# ************************************************************
#purpose: Set all other variables constant (e.g. 10%, mean, 90%)
# Plot predictions for each variable over the range + buffer can check for
# non-linear relationships and if the shape of the regression is maintained
# across other variable.

cmic_clean <- glc_proc_data()

# RF model ----------------------------------------------------------------

model
predictors <- model$trainingData %>% names %>% .[-length(.)]

# prediction variables ----------------------------------------------------
# remove masks

rasts <- stack(map(predictors, ~glc_get_resamp(.x, 2013, "5-15")))
names(rasts) <- predictors

# ~ 2 min
rdf <- as.data.frame(rasts, xy = TRUE, na.rm = TRUE) %>% 
  mutate(pid = glc_pIDfromXY(x, y))


# mask
fullmask <- glc_fullmask_df() %>% select(pid, mask)

rdf <- rdf %>% left_join(fullmask, by = "pid") %>% filter(is.na(mask))

# check
ggplot(rdf, aes(x, y, fill = as.factor(land_cover)))+
  borders()+
  geom_raster()+
  scale_fill_viridis_d(direction = -1)+
  coord_fixed()

# land_cover to text
rdf <- rdf %>% mutate(land_cover = glc_LC_num2chr(land_cover))

rdf$land_cover %>% unique

# to plot, all
predictors # all but land_cover
toplot <- predictors %>% .[!. == "land_cover"]

# toplot <- c("elev", "clay", "nitrogen", "phh2o", "sand", "soc",
#             "tmean", "prec", "ndvi")

## Using boxplots --------------------------------------------

dat <- rdf

# predictors %>% dput
# Elevation, Clay content (%), Nitrogen content, pH, Sand content (%), Soil organic carbon, Land cover type, Mean temperature, Precipitation, NDVI

ren_axis <- function(var) {
  
  tib <- tibble(var = c("elev", "clay", "nitrogen", "phh2o", "sand", "soc", "land_cover", 
                        "tmean", "prec", "ndvi"),
                label = c("Elevation (m)", "Clay content (%)","Nitrogen content", "Soil pH", 
                          "Sand content (%)", "Soil organic carbon", "Land cover type", 
                          "Mean temperature (°C)", "Precipitation", "NDVI"))
  
  lab <- tib$label[tib$var == var]
  
  if (length(lab) == 1) lab else "error"
  
}

# plot per factor boxplot
plprfac_box <- function(dat, model, xvar, npts = 150, probs = c(0.25, 0.75), facet = FALSE) {
  # dat = rdf; xvar = "tmean"
  
  #***
  # npts = 150; probs = c(0.25, 0.75); dat = rdf; xvar = "tmean"
  #***
  
  
  #random forest predictors
  train_data <- model$trainingData
  predictors <- train_data %>% names %>% .[-length(.)]
  
  
  # only use predictor columns, remove NA
  dat <- dat %>% select(all_of(predictors)) %>% drop_na
  
  # same classes
  stopifnot(
    isTRUE(all.equal(map_chr(train_data, class) %>% .[-length(.)], 
                     map_chr(dat, class))))
  
  
  # factor (character treated as factor)
  non_num <- which(!map_lgl(dat, is.numeric))
  
  fac_predictor <- names(dat)[non_num]
  num_predictors <- names(dat)[-non_num]
  levs <- unique(dat %>% pull(non_num))
  datnum <- dat %>% select(where(is.numeric))
  
  
  quantdf <- function(x, probs = c(0.25, 0.75)){
    
    map_dfc(x, ~ quantile(.x, na.rm = T, probs = probs))
    
  }
  
  
  #set quantiles based on factors
  # `$`(dat, "land_cover") #could also be used
  quant_dfs <- map(levs, ~ dat %>% filter(!!sym(fac_predictor) == .x) %>%
                     select(all_of(num_predictors)) %>% 
                     quantdf(probs = probs)) %>% set_names(levs)
  
  
  #*** create df with points to predict
  # pred_cst <- num_predictors[-which(num_predictors == xvar)]
  # npred <- length(num_predictors)
  # nquant <- length(probs)
  # graphs <- nquant^(npred-1)
  # nrows <- graphs * npts
  
  
  # deprec: xvect should fit LC
  # xvect <- seq(min(dat[,xvar], na.rm = T), max(dat[,xvar], na.rm = T), length.out = npts) 
  
  
  ranges <- map(levs, ~ dat %>% filter(!!sym(fac_predictor) == .x) %>% pull({{xvar}}) %>% range)
  
  # TODO scale npts to vector length (just for efficiency)
  # don't fox boxplots
  xvects <- map(ranges, ~ seq(.x[1], .x[2], length.out = npts)) %>% set_names(levs)
  
  
  # test
  # dftopred <- cross_df(c(quant_dfs[[1]] %>% select(-{{xvar}}), list(xvar = xvect))) #
  
  # expand_grid can be used instead to repeat whole df
  dftopred <- imap_dfr(quant_dfs, ~ cross_df(c(.x %>% select(-{{xvar}}), set_names(xvects[.y], xvar), land_cover = .y)))
  
  
  # checks
  # dftopred %>% slice_sample(n = 200) %>% View
  # map_dbl(dftopred, ~ length(unique(.x)))
  
  
  # need to ID single lines
  
  code <- dftopred %>% select(-{{xvar}}) %>% mutate(across(where(is.numeric), ~as.numeric(as.factor(.x)))) %>% 
    unite(code, everything())
  dftopred <- dftopred %>% bind_cols(code)
  
  # dftopred <- dftopred %>% mutate(across(c(where(is.numeric), -{{xvar}}),
  #                                        ~as.numeric(as.factor(.x)))) %>% 
  #   unite(code, everything(), -{{xvar}}, remove = FALSE)
  
  
  #predict
  dftopred <- dftopred %>% mutate(land_cover = factor(land_cover, levels = levels(model$trainingData$land_cover)))
  dftopred <- dftopred %>% mutate(predicted = predict(model, dftopred))
  
  
  # boxplots
  dftopred %>% names
  frange <- dftopred %>% pull({{xvar}}) %>% range
  
  br = seq(frange[1]-0.0001, frange[2], length.out = 10)
  labs = rowMeans(cbind(br[-1], br[-length(br)])) %>% #map2_dbl(br[-1], br[-length(br)], ~mean(c(.x, .y)))
    signif(4)
  
  df2 <- dftopred %>% 
    mutate(gr = cut(dftopred[,xvar, drop = TRUE], breaks= seq(frange[1]-0.0001, frange[2], length.out = 10),
                    labels = labs))
  
  # rem FC ***
  # df2 <- df2 %>% filter(!land_cover == "FC")
  
  # Use full names
  lc_labs <- glc_land_cover_classes()[c(1,3)] %>% mutate(class = str_to_sentence(class))
  df2 <- df2 %>% left_join(lc_labs, by = c("land_cover" = "code")) %>% 
    mutate(land_cover = class)
  
  col_val <- c("Cropland" = "#ffe699",
               "Tropical forest" = "#e2efda",
               "Broadleaf forest" = "#a9d08e",
               "Coniferous forest" = "#548235",
               "Shrubland" = "#f8cbad",
               "Grassland" = "#f4b084")
  
  lc_lvls <- c("Cropland", "Tropical forest", "Broadleaf forest", "Coniferous forest", "Shrubland", "Grassland")
  
  df2 <- df2 %>% mutate(land_cover = factor(land_cover,
                                            levels = lc_lvls))
  
  # TODO Transform x axis values to proper unit, esp. tmean, pH, etc.
  # df2 <- df2 %>% mutate()
  
  # df2$gr %>% unique %>% as.numeric
  
  p <- ggplot(df2, aes_string("gr", "predicted"))+
    # geom_jitter(alpha = 0.02, position = position_dodge()) + #violin
    geom_point(data = df2 %>% slice_sample(prop = 0.1),
               aes(color = land_cover),  
               size = 0.5, alpha = 0.08, position = position_jitter(),
               show.legend = FALSE)+
    geom_boxplot(aes(fill = land_cover), outlier.shape = NA, color = "grey10")+
    labs(x = ren_axis(xvar),
         y = NULL,
         fill = "Land cover type")+
    scale_fill_manual(values = col_val)+
    scale_color_manual(values = col_val)+
    theme_classic()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1.4, vjust = 1.4),
          legend.title = element_text(face = "bold"),
          panel.grid.major.y = element_line(color = "grey92", size = 0.5),
          panel.grid.minor.y = element_line(color = "grey92", size = 0.5))
  
  if(xvar == "tmean") {
    p + scale_x_discrete(labels = function(x) format(round(as.numeric(as.character(x))/10 -273, 1)))
  } else if (xvar %in% c("clay", "sand", "phh2o")) {
    p + scale_x_discrete(labels = function(x) format(round(as.numeric(as.character(x))/10, 1)))
  } else if (xvar == "elev") {
    p + scale_x_discrete(labels = function(x) format(round(as.numeric(as.character(x)), 1)))
    
  } else {
    p
  }
  
}

# examples
plprfac_box(dat = rdf, model, xvar = "tmean")
plprfac_box(rdf, model, xvar = "prec")

all_plots <- map(toplot, ~ plprfac_box(rdf, model, .x))

leg <- get_legend(all_plots[[1]]+
                    guides(colour = guide_legend(override.aes = list(alpha = 1))))
# ggdraw(leg)

no_leg <- map(all_plots, ~ .x +
                theme(legend.position = "none"))

# 5x2 display 
all_grid <- plot_grid(plot_grid(plotlist = no_leg[1:8], ncol = 2),
                      plot_grid(no_leg[[9]], leg, ncol = 2),
                      ncol = 1, rel_heights=c(4, 1))

# for common y axis
y.grob <- textGrob("Predicted Soil Microbial Carbon", 
                   gp=gpar(fontsize=11),
                   rot=90)

main_final <- plot_grid(y.grob, all_grid, rel_widths = c(0.05, 0.95))

ggsave(plot = main_final,
         here("output/figures", "07-1-plot_preds_wLC_boxplot.png"),
         width = gg_width, height = gg_height*2.5)


## Partial dependance plots -----------------------------------------------

# prep --------------------------------------------------------------------

mod1 <- glc_rf_model()

range_df <- glc_predicted_dataset()

lc_lvls <- model$trainingData %>% pull(land_cover) %>% levels


lc_labs <- glc_land_cover_classes()[c(1,3)] %>% mutate(class = str_to_sentence(class)) %>% pull(2,1)

col_val <- c("C" = "#ffe699",
             "FT" = "#e2efda",
             "FB" = "#a9d08e",
             "FC" = "#548235",
             "S" = "#f8cbad",
             "G" = "#f4b084")


# fct ---------------------------------------------------------------------

ren_axis <- function(var) {
  
  tib <- tibble(var = c("elev", "clay", "nitrogen", "phh2o", "sand", "soc", "land_cover", 
                        "tmean", "prec", "ndvi"),
                label = c("Elevation (m)", "Clay content (%)","Nitrogen content", "Soil pH", 
                          "Sand content (%)", "Soil organic carbon", "Land cover type", 
                          "Mean temperature (°C)", "Precipitation", "NDVI"))
  
  lab <- tib$label[tib$var == var]
  
  if (length(lab) == 1) lab else "error"
  
}


# do work -----------------------------------------------------------------

# function to make plot for variable, one line per LC

# default range from training data
# use median values from training data
glc_pdp <- function(xvar, model = mod1, range_df = NULL, npts = 100) {
  
  # xvar = "tmean"
  
  train_data0 <- model$trainingData
  train_data <- train_data0 %>% .[-length(.)]
  predictors_num <- train_data %>% select(-land_cover) %>% names 
  
  tr_sp <- train_data %>%
    group_split(land_cover) %>%
    set_names(levels(train_data$land_cover))
  
  # create dataset for each land cover
  lvl_lc <- levels(train_data$land_cover)
  
  mk_pr_df <- function(dat) {
    # dat <- tr_sp[[1]]
    
    lc1 <- dat$land_cover %>% unique
    stopifnot(length(lc1) == 1)
    
    # xvals
    if (is.null(range_df)) {
      xrange <- range(dat[xvar])
      
    } else {
      xrange <- range_df %>% 
        filter(land_cover == lc1) %>%
        pull(all_of(xvar)) %>% 
        range(na.rm = TRUE)
      
    }
    
    xvals <- seq(from = xrange[1], to = xrange[2], length.out = npts)
    
    # medians
    meds <- dat %>% select(-c(land_cover, all_of(xvar))) %>% map_dbl(median, na.rm = TRUE)
    
    prdf <- tibble(xvals = xvals,
                   land_cover = lc1) %>% 
      bind_cols(meds %>% bind_rows) #meds %>% as.list %>% as_tibble # or tibble(!!!meds)
    
    names(prdf)[1] <- xvar
    
    prdf
    
  }
  
  to_predict <- map_dfr(tr_sp, mk_pr_df)
  
  preds <- to_predict %>% mutate(predicted = predict(model, to_predict))
  
  p <- ggplot(preds, aes_string(xvar, "predicted", color = "land_cover"))+
    geom_point(data = train_data0, aes(y = .outcome), alpha = 0.4, fill = NA)+
    geom_line(size = 1)+
    ylim(0, NA)+ #270 for curves only
    labs(x = ren_axis(xvar),
         y = NULL,
         color = "Land-cover type")+
    scale_color_manual(values = col_val,
                       labels = lc_labs)
  
  if (!exists("leg")) leg <<- get_legend(p)
  # ggdraw(leg)
  
  p <- p + theme(legend.position = "none")
  
  if(xvar == "tmean") {
    p + scale_x_continuous(labels = function(x) format(round(as.numeric(as.character(x))/10 -273, 1)))
  } else if (xvar %in% c("clay", "sand", "phh2o")) {
    p + scale_x_continuous(labels = function(x) format(round(as.numeric(as.character(x))/10, 1)))
  }
  
  (pbot <- p +
      ylim(0, 270)+
      theme(plot.margin = margin(-0.5,5.5,5.5,5.5),
            panel.border = element_blank())+
      annotate(geom = 'segment', y = Inf, yend = -Inf, color = 'black', x = -Inf, xend = -Inf, size = 1)+
      annotate(geom = 'segment', y = Inf, yend = Inf, linetype = "dashed", color = 'black', x = -Inf, xend = Inf, size = 1)+
      annotate(geom = 'segment', y = Inf, yend = -Inf, color = 'black', x = Inf, xend = Inf, size = 0.5)+
      annotate(geom = 'segment', y = -Inf, yend = -Inf, color = 'black', x = -Inf, xend = Inf, size = 0.5)
  )
  
  (ptop <- p +
      ylim(271, 1000)+
      xlab(NULL)+
      theme(axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            plot.margin = margin(5.5,5.5,-2,5.5),
            panel.border = element_blank())+
      annotate(geom = 'segment', y = Inf, yend = -Inf, color = 'black', x = -Inf, xend = -Inf, size = 1)+
      annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+
      annotate(geom = 'segment', y = Inf, yend = -Inf, color = 'black', x = Inf, xend = Inf, size = 0.5))
  # annotate(geom = 'segment', y = -Inf, yend = -Inf, linetype = "dashed", color = 'black', x = -Inf, xend = Inf, size = 1))
  
  plot_grid(ptop, pbot, nrow = 2, rel_heights = c(1, 3), align = "v")
  
}

plots <- map(predictors_num, glc_pdp)

# 5x2 display 
all_grid <- plot_grid(plot_grid(plotlist = plots[1:8], ncol = 2),
                      plot_grid(plots[[9]], leg, ncol = 2),
                      ncol = 1, rel_heights=c(4, 1))

# for common y axis
y.grob <- textGrob("Predicted Soil Microbial Carbon", 
                   gp=gpar(fontsize=11),
                   rot=90)

main_final <- plot_grid(y.grob, all_grid, rel_widths = c(0.05, 0.95))

ggsave(plot = main_final,
         here("output/figures", "11-1-partial_dependence_plots.png"),
         width = gg_width, height = gg_height*2.5)

         
# ************************************************************
# -------------- Stock calculation (07-2) -----------------
# ************************************************************
#purpose: process predicted values, calculate cmic stocks, save dataset

# load data ----------------------------------------------

# load all with specified year
pred_files <- list.files(here("derived", "06-prediction_years"), pattern = ".rds", full.names = TRUE)

# load all predictions (~ 2min)
preds <- map_dfr(pred_files, readRDS, .id = "year")

# fix year
preds$year <- as.numeric(preds$year) + 1991

# calculate cstocks
bulk <- raster(here("geodata/resampled_0p05deg/static/bdod", "bdod_5-15cm_mean_resamp.tif"))

cofr <- raster(here("geodata/resampled_0p05deg/static/cfvo", "cfvo_5-15cm_mean_resamp_c.*tif"))

rarea <- area(cofr)

st <- stack(bulk, cofr, rarea)
names(st) <- c("bulk", "cofr", "area")

# ~ 2 min
bulk_df <- as.data.frame(st, xy = TRUE, na.rm = TRUE)
bulk_df <- bulk_df %>% mutate(pid = glc_pIDfromXY(x, y)) %>% 
  select(-c(x, y))

preds <- preds %>% left_join(bulk_df, by = "pid")


# stock is in tonnes/ha (weight/area)
# cell_stock is the total weight for the cell in tonnes
preds <- preds %>% mutate(stock = pred * bulk * (1000-cofr) * 12.01 / 10^8,
                          cell_stock = stock * area * 100) 

# remove missing stocks
preds <- preds %>% drop_na(stock)

# *** save ***
preds %>% names %>% writeLines
preds <- preds %>% select(-c(elev, clay, nitrogen, phh2o, sand, soc, bulk, cofr, area))
saveRDS(preds, here("derived", "07-cmic_stocks_dataset_all_years.rds"))

# load with glc_predicted_dataset()



# ************************************************************
# ----------------- Slope calculation (07-3) ------------------
# ************************************************************
#purpose: calculate and plot predicted cmic slope per grid cell.
# using df_slope_lm()

preds <- glc_predicted_dataset()


# slope per pixel ---------------------------------------------------------

# nest dataframes for each pixel and calculate slope
pred2 <- preds %>% select(pid, year, stock) %>% nest(data = c(year, stock))
pred2 <- pred2 %>% bind_cols(glc_XYfrompID(.$pid))

# example with one point
dat1 <- pred2$data[[1]]
dat1$stock[4] <- NA
df_slope_lm(df = dat1[-4,])
ggplot(dat1, aes(year, stock)) +
  geom_smooth(method = "lm") +
  geom_point()


# run for all (~ 45 min). Much faster because not saving/writing full lm models
pred2 <- pred2 %>% 
  mutate(slope_calc = imap(data, ~ {if (.y %% 50000 == 0) cat(.y, "done\n");
    df_slope_lm(.x)} )) 

# pred2 <- pred2 %>% slice(which(!misind)) # what was this line for?

slo_df <- pred2 %>% mutate(slope = map_dbl(slope_calc, 1),
                           pval = map_dbl(slope_calc, 2),
                           nona = map_dbl(slope_calc, 3)) %>% 
  select(-c(data, slope_calc))

slo_df <- slo_df %>% filter(!is.na(slope))

saveRDS(slo_df, here("derived", "07-3-slope_per_pix_LM.rds"))


# figure ------------------------------------------------------------------
# *** Load from here for figure

slo_df <- readRDS(here("derived", "07-3-slope_per_pix_LM.rds"))
slo_df <- slo_df %>% filter(!is.na(slope))

# main figures -------------------------------------------------------------

# without masking
slo_df <- slo_df %>% mutate(slope_cut = cut(slope, breaks = quantile(slope, probs = seq(0, 1, length.out = 6)),
                                            include.lowest = TRUE))
# levels(slo_df$slope_cut) <- rev(levels(slo_df$slope_cut))

# manual cuts
breaks <- quantile(slo_df$slope, probs = seq(0, 1, length.out = 4))
mean(abs(breaks[2:3]))
breaks %>% dput
breaks <- c(`0%` = -0.222891884675364, `33.33333%` = -0.0007, 
            `66.66667%` = 0.0007, `100%` = 0.201048117709902
)

slo_df <- slo_df %>% mutate(slope_cut = cut(slope, breaks = breaks,
                                            include.lowest = TRUE))

ggplot(slo_df, aes(x, y, fill = slope_cut))+
  borders(fill = "grey90", colour = NA, ylim = c(-60, 90))+
  geom_raster()+ # key_glyph = "polygon3"
  coord_fixed()+
  scale_fill_manual(values = c("#dc7656", "#8c57db", "#57d0db"),
                    labels = function(x) signif(x, 3))+
  guides(fill = guide_bins(title = "Rate of change (t/ha)",
                           axis = FALSE, show.limits = T, title.position = "top"))+ #, title.hjust = 0
  # labs(fill = "Rate of change")+
  theme_void()+
  theme(legend.title = element_text(face = "bold"),
        legend.position = c(0.61, 0),
        legend.direction = "horizontal",
        legend.key.width = unit(1.6, "cm")) +
  annotate("text", -Inf, Inf, hjust = -0.5, vjust = 1.5, label = "B", size = 16/.pt,
           fontface = "bold")

ggsave(here("output/figures", "07-3-stock_slope_cut.png"),
       width = gg_width, height = gg_height)


## relative ----------------------------------------------------------------

pred_mean <- preds %>% group_by(pid) %>% 
  summarise(mstock = mean(stock, na.rm = TRUE))

slo_df <- slo_df %>% left_join(pred_mean, by = "pid")
slo_df <- slo_df %>% mutate(slo_rel = slope/mstock * 100) #percentage change per year


# manual cuts
breaks <- quantile(slo_df$slo_rel, probs = seq(0, 1, length.out = 4)) #6
breaks <- c(`0%` = -13.4415961717228, `33.33333%` = -0.13, 
            `66.66667%` = 0.13, `100%` = 19.3287866611536)

slo_df <- slo_df %>% mutate(slo_rel_cut = cut(slo_rel, breaks = breaks,
                                              include.lowest = TRUE))

# slo_df %>% count(slo_rel_cut)
# slo_df$slo_rel %>% sort(F) %>% head(20)
# slo_df$slo_rel %>% sort(T) %>% head(20)

(p_slope_rel <- ggplot(slo_df, aes(x, y, fill = slo_rel_cut))+
    borders(fill = "grey90", colour = NA, ylim = c(-60, 90))+
    geom_raster()+
    coord_fixed(xlim = c(-180, 180), expand = FALSE)+
    # scale_fill_manual(values = c("red", "pink", "yellow", "lightblue", "blue"))+
    scale_fill_manual(values = c("#dc7656", "#8c57db", "#57d0db"))+
    theme_void()+
    guides(fill = guide_bins(title = "Rate of change for\n1992-2013 (% per year)",
                             axis = FALSE, show.limits = T, title.position = "top"))+
    theme(legend.position = c(0.63, 0),
          legend.direction = "horizontal",
          legend.title = element_text(face = "bold"),
          # legend.title.align = 0,
          legend.key.width = unit(1.2, "cm"))#+
  # annotate("text", -Inf, Inf, hjust = -0.5, vjust = 1.5, label = "C", size = 16/.pt,
  #          fontface = "bold")
)




# plot pval with cuts
slo_df <- slo_df %>% mutate(pval_cut = cut(pval, breaks = c(0, 0.05, 0.1, 0.5, 1), include.lowest = TRUE))
slo_df$pval_cut %>% levels

ggplot(slo_df, aes(x, y, fill = pval_cut))+
  borders(fill = "grey90", colour = NA, ylim = c(-60, 90))+
  geom_raster()+
  coord_fixed(xlim = c(-180, 180), expand = FALSE)+
  scale_fill_manual(values = c("blue", "lightblue", "pink", "red"))+
  guides(fill = guide_bins(title = "P-value of the temporal trend",
                           axis = FALSE, show.limits = T, title.position = "top"))+
  # labs(fill = "p-value")+
  theme_void()+
  theme(legend.title = element_text(face = "bold"),
        legend.position = c(0.61, 0),
        legend.direction = "horizontal",
        legend.key.width = unit(1.2, "cm")) +
  annotate("text", -Inf, Inf, hjust = -0.5, vjust = 1.5, label = "C", size = 16/.pt,
           fontface = "bold")

ggsave(here("output/figures", "07-3-slope_pvalue_supp.png"),
       width = gg_width, height = gg_height)



# ************************************************************
# ----------------- GIF Animation (07-4) -------------------
# ************************************************************
#purpose: prediction maps all years, saved as animated GIF, Cmic stocks

# make figures for GIF ----------------------------------------------
# NOTE masked already and stocks calculated

# predictions
preds <- glc_predicted_dataset()

# quantiles and figures ---------------------------------------------------

# quantiles
preds <- preds %>%
  mutate(CmicStock = cut(stock, breaks = quantile(stock, probs = seq(0, 1, length.out = 7)),
                         include.lowest = TRUE))

preds$CmicStock %>% unique %>% as.character %>% paste0(., '" = "', ., '",\n"') %>% cat(sep = "")

preds <- preds %>%
  mutate(CmicStock = fct_recode(CmicStock,
                                "(1.5,7.0]" = "(1.55,6.96]",
                                "(0.67,1.5]" = "(0.672,1.55]",
                                "(0.50,0.67]" = "(0.492,0.672]",
                                "(0.38,0.50]" = "(0.385,0.492]",
                                "(0.30,0.38]" = "(0.313,0.385]",
                                "[0.06,0.30]" = "[0.0598,0.313]"
  ))

# rename

time_sca <- tibble(x = c(-150, -150, -50), y = c(77, 79, 79),
                   xend = c(-50, -150, -50), yend = c(77, 75, 75))

# save a png for each year
for (iyear in 1992:2013) {
  # iyear <- 1992
  
  cat("*** ", iyear, " ***")
  
  pred_iyear <- preds %>% filter(year == iyear)
  
  # NOTE warning about uneven intervals due to float rounding, can ignore
  p <- ggplot(pred_iyear, aes(x, y, fill = CmicStock)) + 
    borders(ylim = c(-60, 90), fill = "grey90", colour = NA)+
    geom_raster()+ #key_glyph = "polygon3"
    coord_fixed(xlim = c(-180, 180), expand = FALSE)+
    scale_fill_viridis_d()+ # guide = guide_legend(reverse = TRUE)
    theme_void()+
    guides(fill = guide_bins(title = "Microbial Carbon Stock (t/ha)",
                             axis = FALSE, show.limits = T, title.position = "top"))+
    theme(legend.position = c(0.63, 0),
          legend.direction = "horizontal",
          legend.title = element_text(face = "bold"),
          # legend.title.align = 0,
          legend.key.width = unit(1.2, "cm")) +
    # moving year bar
    annotate("rect", xmin = -158, xmax = -42, ymin = 65, ymax = 84,
             fill = "white", alpha = 0.3, color = "grey")+
    geom_segment(data = time_sca, aes(x, y, xend = xend, yend = yend),
                 inherit.aes = FALSE)+
    annotate(geom = "text", x = -150, y = 70, label = 1992)+
    annotate(geom = "text", x = -50, y = 70, label = 2013)+
    annotate("point", x = -150 + (-50 - -150) * (iyear - 1992) / (2013 - 1992),
             y = 77, color = "black", fill = "grey40", shape = 21, size = 5)
  
  ggsave(plot = p,
         here("output/anim", "07-pred_bins_gif", paste0("pred_bins_", iyear, ".png")),
         width = gg_width, height = gg_height)
  
}


# 2013 for main text ------------------------------------------------------

pred_2013 <- preds %>% filter(year == 2013)

p2013 <- ggplot(pred_2013, aes(x, y, fill = CmicStock)) + 
  borders(ylim = c(-60, 90), fill = "grey90", colour = NA)+
  geom_raster()+ #key_glyph = "polygon3"
  coord_fixed(xlim = c(-180, 180), expand = FALSE)+
  scale_fill_viridis_d(direction = -1, begin = 0.05)+ # guide = guide_legend(reverse = TRUE)
  theme_void()+
  guides(fill = guide_bins(title = "Microbial Carbon Stock (t/ha)",
                           axis = FALSE, show.limits = T, title.position = "top"))+
  theme(legend.position = c(0.63, 0),
        legend.direction = "horizontal",
        legend.title = element_text(face = "bold"),
        # legend.title.align = 0,
        legend.key.width = unit(1.2, "cm"))#+
# annotate("text", -Inf, Inf, hjust = -0.5, vjust = 1.5, label = "B", size = 16/.pt,
#          fontface = "bold")


# add var importance for main figure
vimp_df <- readRDS(here("derived", "06-1-vimp_df.rds"))

p_vimp <- ggplot(vimp_df, aes(importance, y = reorder(variable, importance, FUN = stats::median),
                              fill = reorder(variable, importance, FUN = stats::median))) +
  geom_boxplot() +
  # scale_fill_viridis_d(direction = -1, alpha = 0.8) +
  scale_fill_grey(start = 0.35, end = 0.9) +
  # ylab("Variable") +
  ylab(NULL) +
  xlab("Variable importance") +
  theme_bw() +
  theme(legend.position = "none",
        panel.border = element_blank(),
        axis.line = element_line(colour = "black", 
                                 size = rel(1)),
        plot.margin = unit(c(20, 5.5, 5.5, 5.5), "pt")
  )

p_vimpl <- ggdraw(p_vimp) + 
  draw_figure_label("A", "top.left", fontface = "bold", size = 16)

fplot <- ggdraw(p2013)+
  draw_plot(p_vimp, 0.01, 0.04, 0.25, 0.42)
fplot

ggsave(plot = fplot,
         here("output/figures", paste0("07-4-predictions_2013_main.png")),
         width = gg_width, height = gg_height)


# Full figure 2 -----------------------------------------------------------
# make figure C in script 07-3-lm_slope

# full_fig2 <- plot_grid(fplot, p_slope_rel, ncol = 1)

p2013l <- ggdraw(p2013) + 
  draw_figure_label("B", "top.left", fontface = "bold", size = 16)

p_slope_rell <- ggdraw(p_slope_rel) + 
  draw_figure_label("C", "top.left", fontface = "bold", size = 16)

full_fig2 <- plot_grid(
  plot_grid(p_vimpl, p2013l, nrow = 1, rel_widths = c(0.3, 0.7)),
  p_slope_rell,
  ncol = 1,
  rel_heights = c(0.7, 1)
)

ggsave(plot = full_fig2,
         here("output/figures", "07-4-Figure_2_fullABC.png"),
         width = gg_width, height = gg_height*1.7)


# make GIF with magick ----------------------------------------------------
# otherwise make GIF with GIMP

## list file names and read in
imgs <- list.files(here("output/anim/07-pred_bins_gif"), 
                   pattern = ".png",
                   full.names = TRUE)
img_list0 <- map(imgs, image_read)
image_info(img_list0[[1]])

img_list <- map(img_list0, ~image_scale(.x, "1200"))

## join the images together
img_joined <- image_join(img_list)

## animate at 2 frames per second
img_animated <- image_animate(img_joined, fps = 2)

## view animated image
# img_animated

## save to disk
image_write(image = img_animated,
            path = here("output/anim", "07-pred_bins_magick.gif"))


# ************************************************************
# ------------------ Bivariate map (08-0) -------------------
# ************************************************************
#purpose: bivariate map, mostly for proof of concept
# use rate of change (%) and Cmic val

# Cmic vals
preds <- glc_predicted_dataset()
pred_mean <- preds %>% group_by(pid) %>% 
  summarise(mstock = mean(stock, na.rm = TRUE))

# % slope
slo_df <-readRDS(here("derived", "07-3-slope_per_pix_LM.rds"))
slo_df <- slo_df %>% left_join(pred_mean, by = "pid")
df <- slo_df %>% mutate(slo_rel = slope/mstock * 100) #percentage change per year
df <- df %>% drop_na(mstock, slo_rel)


# make bivariate map ------------------------------------------------------
# https://timogrossenbacher.ch/2019/04/bivariate-maps-with-ggplot2-and-sf/
# biscale: https://cran.r-project.org/web/packages/biscale/vignettes/biscale.html
# 
# additional references:
# continuous palettes: https://www.datalorax.com/post/creating-bivariate-color-palettes/
# https://rpubs.com/ayushbipinpatel/593942
# concepts: https://www.joshuastevens.net/cartography/make-a-bivariate-choropleth-map/
# 
# Only used for vector graphics, so will need to adjust

br_cmst <- quantile(df$mstock, probs = seq(0, 1, length.out = 4)) #6
br_rslo <- quantile(df$slo_rel, probs = seq(0, 1, length.out = 4)) #6
mean(abs(br_rslo[2:3]))
br_rslo[2:3] <- c(-0.13, 0.13)

bim_df <- df %>% mutate(cmst_class = case_when(mstock < br_cmst[2] ~ 1,
                                               mstock < br_cmst[3] ~ 2,
                                               TRUE ~ 3),
                        rslo_class = case_when(slo_rel  < br_rslo[2] ~ 1,
                                               slo_rel < br_rslo[3] ~ 2,
                                               TRUE ~ 3),
                        bi_class = paste0(cmst_class, "-", rslo_class))


bidf <- df %>% bi_class(x = mstock, y = slo_rel, style = "quantile", dim = 3)

summ_bi <- . %>% group_by(bi_class) %>% summarise(mincmic = min(mstock),
                                                  maxcmic = max(mstock),
                                                  minslo = min(slo_rel),
                                                  maxslo = max(slo_rel),
                                                  n = n())

bidf %>% summ_bi
bim_df %>% summ_bi

# heatmap of the legend
bim_red <- bim_df %>% slice_sample(n = 200000)

pheat <- ggplot(bim_red, aes(mstock, slo_rel, color = bi_class))+ 
  geom_point(size = 0.8, alpha = 0.1, show.legend = FALSE)+
  scale_color_manual(values = c("#f3ccc1", "#d4c0f2", "#c0edf2",
                                "#dc7656", "#8c57db", "#57d0db",
                                "#933b1f", "#4e1f93", "#1f8993"))+
  labs(x = "Microbial carbon stock (t/ha)",
       y = "Relative rate of change (% per year)") #+
# annotate("text", -Inf, Inf, hjust = -0.5, vjust = 1.5, label = "B", size = 16/.pt,
#          fontface = "bold")

ggsave(here("output/figures", "08-bivariate_legend_heatmap.png"),
         width = gg_width*.6, height = gg_height*.7)


# full map
man_pal <- bi_pal_manual("#f3ccc1", "#d4c0f2", "#c0edf2",
                         "#dc7656", "#8c57db", "#57d0db",
                         "#933b1f", "#4e1f93", "#1f8993", preview = F)

scales::show_col(c("#dc7656", "#8c57db", "#57d0db"))

(p <- ggplot(bim_df, aes(x, y, fill = bi_class))+
    borders(ylim = c(-60, 90), fill = "grey90", colour = NA)+
    geom_raster(show.legend = FALSE)+
    coord_fixed(xlim = c(-180, 180), expand = FALSE)+
    bi_scale_fill(pal = man_pal, dim = 3) +
    theme_void()) #bi_theme()


legend <- bi_legend(pal = man_pal,
                    dim = 3,
                    xlab = "Microbial Carbon Stock",
                    ylab = "Increase\n\n\n\n\nDecrease",
                    size = 8) +
  theme(axis.title.y.left = element_text(angle = 0, vjust = 0.15))
legend

finalPlot <- ggdraw(p) +
  draw_plot(legend, 0.01, 0.2, 0.28, 0.28)
finalPlot +
  annotate("text", -Inf, Inf, hjust = -0.5, vjust = 1.5, label = "A", size = 16/.pt,
           fontface = "bold")

ggsave(here("output/figures", "08-bivariate_plot.png"),
         width = gg_width, height = gg_height)



# ************************************************************
# -------------------- Regional analysis ---------------------
# ************************************************************
# regional
# 09-0 to 09-4

# IPBES regions taken from https://zenodo.org/record/3928281#.Yhzup-iZOUk
ipbes <- st_read(here("geodata/ipbes_regions/ipbes_regions_subregions_shape_1.1/IPBES_Regions_Subregions2.shp"))

ipbes$subfct <- factor(ipbes$Sub_Region)
levels(ipbes$subfct) %>% dput
levs <- c("Antarctica", "Caribbean", "Central Africa", "Central and Western Europe", 
          "Central Asia", "East Africa and adjacent islands", "Eastern Europe", 
          "Mesoamerica", "North-East Asia", "North Africa", "North America", 
          "Oceania", "South-East Asia", "South America", "South Asia", 
          "Southern Africa", "West Africa", "Western Asia")

match_tbl <- tibble(name = levels(ipbes$subfct),
                    code = 1:length(name))


# template
rtemp <- glc_get_resamp("ndvi", 2013)

# fun = "last" gives same result, polygons don't ovelap
raspol <- fasterize::fasterize(ipbes, rtemp, field = "subfct", fun = "first")

# s1 <- spplot(raspol, maxpixels = 500000)
writeRaster(raspol, here("geodata", "09-ipbes_subregions_alphaOrder.tif"))




#purpose: temporal regional slope analysis

# load data ----------------------------------------------
preds <- glc_predicted_dataset()
preds %>% names %>% writeLines

preds <- preds %>% 
  select(pid, x, y, year, tmean, prec, ndvi, land_cover, pred, stock, cell_stock)

# NOTE masked areas already removed

# Match IPBES region
ipras <- raster(here("geodata", "09-ipbes_subregions_alphaOrder.tif"))
names(ipras) <- "ipbco"

levs <- c("Antarctica", "Caribbean", "Central Africa", "Central and Western Europe", 
          "Central Asia", "East Africa and adjacent islands", "Eastern Europe", 
          "Mesoamerica", "North-East Asia", "North Africa", "North America", 
          "Oceania", "South-East Asia", "South America", "South Asia", 
          "Southern Africa", "West Africa", "Western Asia")

match_tbl <- tibble(ipbsr = levs,
                    code = 1:length(ipbsr))

ip_df <- as.data.frame(ipras, xy = TRUE, na.rm = TRUE)
ip_df <- ip_df %>% 
  left_join(match_tbl, by = c("ipbco" = "code")) %>% 
  mutate(pid = glc_pIDfromXY(x, y)) %>% 
  select(-c(x,y))

preds <- preds %>% left_join(ip_df, by = "pid")
preds <- preds  %>% drop_na(ipbsr)


# nest per cell --------------------------------------------------------

prn <- preds %>% select(pid, year, cell_stock, ipbsr, prec, tmean, ndvi, land_cover) %>%
  group_by(pid, ipbsr) %>%
  nest %>% ungroup

# remove if missing temporal datapoints
nrow(prn) #2658500
prn <- prn %>% filter(map_lgl(data, ~ nrow(.x) == 22))
nrow(prn) #2596897

jreg <- prn %>% unnest(cols = c(data))


# global temporal trend ---------------------------------------------------

# NOTE cell_stock is in tonnes
glob_sum <- jreg %>% group_by(year) %>% 
  summarise(cstock = sum(cell_stock))

# mean value cstock
mcs <- glob_sum$cstock %>% mean #4.34 Pg Cmic

glob_sum <- glob_sum %>% mutate(change = cstock - cstock[1],
                                perc_ch = change / cstock[1] * 100,
                                ysca = year - 1992)

glob_sum <- glob_sum %>% mutate(mod_fix = "full")


# Calculation global cmstock changes ---------------------------------------

# how much C loss in 21 years?
# NOTE units in tonnes
mody <- lm(cstock ~ ysca, glob_sum)

# absolute change per year
abs_ch_p_year <- mody$coefficients["ysca"]

#over studied period
abs_ch_p_year * 21

# relative change per year
rel_ch_p_year <- mody$coefficients["ysca"] / mcs * 100

# global change over studied period, in percent
# (mody$coefficients["ysca"] * 21) / mcs * 100
rel_ch_p_year * 21



ggplot(glob_sum, aes(ysca, perc_ch)) +
  geom_smooth(method = "lm")

pval <- if (pval < 0.001) "< 0.001" else paste("=", round(pval, 3))

tpc <- paste0(round((2013-1992)*mod$coefficients[1], 2), " %\np.val ", pval)

(p <- ggplot(glob_sum, aes(ysca, perc_ch))+
    geom_line()+
    # geom_abline(slope = mod$coefficients[1], intercept = 0)+
    geom_smooth(method = "lm", formula = y ~ 0 + x)+
    scale_x_continuous(labels = function(x) x + 1992)+
    ggtitle("Global percentage change in microbial carbon stocks")+
    labs(x = "Year", y = "Percentage change")+
    annotate("text", x = 20, y = 3, label = tpc))

ggsave(here("output/figures", "09-global_cmicst_trend.png"),
         width = gg_width*.6, height = gg_height*.8)

gl_tib <- tibble(ipbsr = "Global",
                 slope = coef(mod)[2],
                 pval = coef(summary(mod))[2,4], 
                 rsq = summary(mod)$r.squared)

# regional trends ---------------------------------------------------------

reg_sum <- jreg %>% group_by(ipbsr, year) %>% 
  summarise(cstock = sum(cell_stock),
            prec = mean(prec))

reg_sum_n <- reg_sum %>% group_by(ipbsr) %>% nest

# fixed scale with zero
(p <- ggplot(reg_sum, aes(year, cstock))+
    geom_line(size = 1)+
    geom_smooth(method = "lm", size = 0.5)+
    theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=1))+
    ylim(0, NA)+
    facet_wrap(~ipbsr))

ggsaveme(here("output/figures", "09-regional_trends_ipbsub.png"),
         width = gg_width*1.5, height = gg_height*1.5)

names(reg_sum)
reg_sum <- reg_sum %>% group_by(ipbsr) %>% mutate(change = cstock - first(cstock),
                                                  perc_ch = change / first(cstock) * 100,
                                                  ysca = year - 1992)

reg_sum <- reg_sum %>% mutate(mod_fix = "full")

ggplot(reg_sum, aes(ysca, perc_ch))+
  geom_line(size = 0.5)+
  geom_smooth(method = "lm", formula = y ~ 0 + x)+
  scale_x_continuous(labels = function(x) x + 1992)+
  ggtitle("Global percentage change in microbial carbon stocks")+
  labs(x = "Year", y = "Percentage change")+
  facet_wrap(~ipbsr)


fct <- function(x) {
  
  mod <- lm(perc_ch ~ year, data = reg_sum %>% filter(ipbsr == x))
  slope <- coef(mod)[2]
  pval <- coef(summary(mod))[2,4] 
  rsq <- summary(mod)$r.squared
  c(slope, pval, rsq)
  
}

mod_par <- lm(perc_ch ~ year, data = reg_sum)

fct(x = unique(reg_sum$ipbsr)[1])

vals <- map(unique(reg_sum$ipbsr), fct)

tib <- tibble(ipbsr = unique(reg_sum$ipbsr),
              slope = map_dbl(vals, 1),
              pval = map_dbl(vals, 2),
              rsq = map_dbl(vals, 3))

tib %>% arrange(pval)

tib <- bind_rows(gl_tib, tib)
tib <- tib %>% mutate(ipbsr = factor(ipbsr, levels = levels(reg_gl$ipbsr)))
levels(tib$ipbsr)


# add global
reg_names <- reg_sum$ipbsr %>% unique
reg_gl <- bind_rows(glob_sum %>% mutate(ipbsr = "Global"), reg_sum)
reg_gl <- reg_gl %>% mutate(ipbsr = factor(ipbsr, levels = c("Global", reg_names)))
levels(reg_gl$ipbsr)

# with bars instead
p_bars <- ggplot(reg_gl, aes(ysca, perc_ch))+
  # geom_rect(data = tib, aes(fill = slope), xmin = -Inf,xmax = Inf,
  #           ymin = -Inf, ymax = Inf, alpha = 0.3, inherit.aes = FALSE) +
  # scale_fill_gradient2()+
  geom_col(size = 0.5, aes(fill = perc_ch > 0))+
  scale_fill_manual(guide = FALSE, breaks = c(TRUE, FALSE), values=c("lightblue", "pink"))+
  geom_smooth(method = "lm", color = "black")+ #, formula = y ~ 0 + x
  scale_x_continuous(labels = function(x) x + 1992)+
  # ggtitle("Regional percentage changes in microbial carbon stocks")+
  labs(x = "Year", y = "Percentage change")+
  # coord_cartesian(ylim=c(-60, 60))+
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=1))+
  geom_text(data = tib, aes(label = paste0("slope = ", format(round(slope, 2), digits = 2),
                                           "\np.val = ", format(round(pval, 3), digits = 3),
                                           "\nR² = ", format(round(rsq, 2), digits = 2))),
            x = Inf, y = Inf, vjust = 1.1, hjust = 1.1, inherit.aes = FALSE)+
  facet_wrap(~ ipbsr)
# annotate("text", -Inf, Inf, hjust = -0.5, vjust = 1.5, label = "A", size = 16/.pt,
#          fontface = "bold")

p_bars

ggdraw(p_bars)+
  draw_plot_label("A")

ggsave(here("output/figures", "09-regio_trends_stand_bars.png"),
         width = gg_width*1.5, height = gg_height*2)



# ************************************************************
# ------------ Dynamic variables trends (09-2) -------------
# ************************************************************

# trends in other variables
# fixed y, scale for prec, but leave as is for temp and ndvi

library(glcpck)
library(cowplot)

jreg <- readRDS(here("derived", "09-jreg_data_b4summ.rds"))
names(jreg)


# add area
cofr <- raster(here("geodata/resampled_0p05deg/static/cfvo"), "cfvo_5-15cm_mean_resamp.tif")
rarea <- area(cofr)
names(rarea) <- "area"

area_df <- as.data.frame(rarea, xy = TRUE, na.rm = TRUE)
area_df <- area_df %>% mutate(pid = glc_pIDfromXY(x, y)) %>% 
  select(-c(x, y))

jreg <- jreg %>% left_join(area_df, by = "pid")


# weigh per area ----------------------------------------------------------
# prec is in Kg/m2, area in km2

jreg <- jreg %>% mutate(prec_Mt = prec * area * 1e+6 / 1e+9,  #in megatonnes
                        tmean_w = tmean * area,
                        ndvi_w = ndvi * area)

reg_sum <- jreg %>% drop_na(ipbsr) %>% 
  group_by(ipbsr, year) %>%
  summarise(cstock = sum(cell_stock),
            prec = sum(prec_Mt),
            tmean_ws = sum(tmean_w),
            # ndvi_mean = mean(ndvi),
            ndvi_ws = sum(ndvi_w),
            area = sum(area)) %>% 
  mutate(tmean = tmean_ws / area,
         ndvi = ndvi_ws / area)

names(reg_sum)  
head(jreg)

reg_sum %>% filter(year == 1992, ipbsr == "Eastern Europe") %>% pull(ndvi) %>% mean



# prec --------------------------------------------------------------------
# as percentage

precdf <- reg_sum %>% group_by(ipbsr) %>% mutate(change = prec - first(prec),
                                                 perc_ch = change / first(prec) * 100,
                                                 ysca = year - 1992)


# ggplot(precdf, aes(ysca, prec, color = ipbsr))+
#   geom_line(size = 0.5)+
#   # geom_smooth(method = "lm", formula = y ~ 0 + x)+
#   scale_x_continuous(labels = function(x) x + 1992)+
#   ggtitle("Regional percentage change in precipitation")+
#   labs(x = "Year", y = "Percentage change")

ggplot(precdf, aes(ysca, perc_ch))+
  geom_line(size = 0.5)+
  geom_smooth(method = "lm", formula = y ~ 0 + x)+
  scale_x_continuous(labels = function(x) x + 1992)+
  ggtitle("Regional percentage change in precipitation")+
  labs(x = "Year", y = "Percentage change")+
  facet_wrap(~ipbsr)


fct <- function(x) {
  
  mod <- lm(perc_ch ~ year, data = precdf %>% filter(ipbsr == x))
  slope <- coef(mod)[2]
  pval <- coef(summary(mod))[2,4] 
  rsq <- summary(mod)$r.squared
  list(x, slope, pval, rsq) %>% set_names("ipbsr", "slope", "pval", "rsq")
  
}

lm(perc_ch ~ 0 + ysca, data = precdf %>% filter(ipbsr == "South Asia"))

# mod_par <- 

fct(x = unique(precdf$ipbsr)[1])

tib <- map_dfr(unique(precdf$ipbsr), fct)

# tib <- tib %>% mutate(col = )

(p1 <- ggplot(precdf, aes(ysca, perc_ch))+
    geom_rect(data = tib, aes(fill = slope), xmin = -Inf,xmax = Inf,
              ymin = -Inf, ymax = Inf, alpha = 0.3, inherit.aes = FALSE) +
    scale_fill_gradient2()+
    geom_line(size = 0.5)+
    geom_smooth(method = "lm", formula = y ~ x, color = "black", size = 0.7, linetype = "dashed")+
    scale_x_continuous(labels = function(x) x + 1992)+
    # ggtitle("Regional percentage changes in precipitation")+
    labs(x = "Year", y = "Percentage change")+
    # coord_cartesian(ylim=c(-100, 100))+
    theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=1))+
    geom_text(data = tib, aes(label = paste0("slope = ", round(slope, 2),
                                             "\np.val = ", round(pval, 3),
                                             "\nR2 = ", round(rsq, 3))),
              x = Inf, y = Inf, vjust = 1.1, hjust = 1.1, inherit.aes = FALSE,
              size = 3.3)+
    facet_wrap(~ipbsr))

ggdraw(p1)+
  annotate("text", -Inf, Inf, hjust = -0.5, vjust = 1.5, label = "B", size = 16/.pt,
           fontface = "bold")

ggsave(here("output/figures", "09-regio_trends_dynvar_prec.png"),
         width = gg_width, height = gg_height*1.3)



# tmean -------------------------------------------------------------------
# plot change

tmeandf <- reg_sum %>% group_by(ipbsr) %>% mutate(tmean = tmean / 10,
                                                  change = tmean - first(tmean),
                                                  # perc_ch = change / first(tmean) * 100,
                                                  ysca = year - 1992)


# ggplot(tmeandf, aes(ysca, tmean, color = ipbsr))+
#   geom_line(size = 0.5)+
#   geom_smooth(method = "lm", formula = y ~ x)+
#   scale_x_continuous(labels = function(x) x + 1992)+
#   ggtitle("Regional percentage change in precipitation")+
#   labs(x = "Year", y = "Percentage change")

ggplot(tmeandf, aes(ysca, change))+
  geom_line(size = 0.5)+
  geom_smooth(method = "lm", formula = y ~ 0 + x)+
  scale_x_continuous(labels = function(x) x + 1992)+
  ggtitle("Regional change in temperature")+
  labs(x = "Year", y = "Temperature change (in C)")+
  facet_wrap(~ipbsr)


fct <- function(x) {
  
  mod <- lm(change ~ year, data = tmeandf %>% filter(ipbsr == x))
  slope <- coef(mod)[2]
  pval <- coef(summary(mod))[2,4] 
  rsq <- summary(mod)$r.squared
  list(x, slope, pval, rsq) %>% set_names("ipbsr", "slope", "pval", "rsq")
  
}

fct(x = unique(tmeandf$ipbsr)[17])

tib <- map_dfr(unique(tmeandf$ipbsr), fct)

p2 <- ggplot(tmeandf, aes(ysca, change))+
  geom_rect(data = tib, aes(fill = slope), xmin = -Inf,xmax = Inf,
            ymin = -Inf, ymax = Inf, alpha = 0.3, inherit.aes = FALSE) +
  scale_fill_gradient2()+
  geom_line(size = 0.5)+
  geom_smooth(method = "lm", formula = y ~ x, color = "black", linetype = "dashed")+
  scale_x_continuous(labels = function(x) x + 1992)+
  # ggtitle("Regional change in temperature")+
  labs(x = "Year", y = "Temperature change (in C)")+
  # coord_cartesian(ylim=c(-100, 100))+
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=1))+
  geom_text(data = tib, aes(label = paste0("slope = ", round(slope, 3),
                                           "\np.val = ", round(pval, 3),
                                           "\nR2 = ", round(rsq, 3))),
            x = Inf, y = Inf, vjust = 1.1, hjust = 1.1, inherit.aes = FALSE,
            size = 3.3)+
  facet_wrap(~ipbsr)

ggdraw(p2)+
  annotate("text", -Inf, Inf, hjust = -0.5, vjust = 1.5, label = "A", size = 16/.pt,
           fontface = "bold")

ggsave(here("output/figures", "09-regio_trends_dynvar_tmean.png"),
         width = gg_width, height = gg_height*1.3)




# ndvi --------------------------------------------------------------------
# also plot untransformed change

ndvidf <- reg_sum %>% group_by(ipbsr) %>% mutate(change = ndvi - first(ndvi),
                                                 # perc_ch = change / first(ndvi) * 100,
                                                 ysca = year - 1992)

ggplot(ndvidf, aes(ysca, ndvi, color = ipbsr))+
  geom_line(size = 0.5)+
  geom_smooth(method = "lm", formula = y ~ x, se = F)+
  scale_x_continuous(labels = function(x) x + 1992)+
  ggtitle("Regional percentage change in precipitation")+
  labs(x = "Year", y = "Percentage change")

ggplot(ndvidf, aes(ysca, change, color = ipbsr))+
  geom_line(size = 0.5)+
  geom_smooth(method = "lm", formula = y ~ 0 + x, se = F)+
  scale_x_continuous(labels = function(x) x + 1992)+
  ggtitle("Regional percentage change in precipitation")+
  labs(x = "Year", y = "Percentage change")

ggplot(ndvidf, aes(ysca, change))+
  geom_line(size = 0.5)+
  geom_smooth(method = "lm", formula = y ~ x)+
  scale_x_continuous(labels = function(x) x + 1992)+
  ggtitle("Regional change in NDVI")+
  labs(x = "Year", y = "NDVI change")+
  facet_wrap(~ipbsr)


fct <- function(x) {
  
  mod <- lm(change ~ year, data = ndvidf %>% filter(ipbsr == x))
  slope <- coef(mod)[2]
  pval <- coef(summary(mod))[2,4] 
  rsq <- summary(mod)$r.squared
  list(x, slope, pval, rsq) %>% set_names("ipbsr", "slope", "pval", "rsq")
  
}

fct(x = unique(ndvidf$ipbsr)[1])

tib <- map_dfr(unique(ndvidf$ipbsr), fct)

p3 <- ggplot(ndvidf, aes(ysca, change))+
  geom_rect(data = tib, aes(fill = slope), xmin = -Inf,xmax = Inf,
            ymin = -Inf, ymax = Inf, alpha = 0.3, inherit.aes = FALSE) +
  scale_fill_gradient2()+
  geom_line(size = 0.5)+
  geom_smooth(method = "lm", formula = y ~ x, color = "black", linetype = "dashed")+
  scale_x_continuous(labels = function(x) x + 1992)+
  ggtitle("Regional change in NDVI")+
  labs(x = "Year", y = "NDVI change")+
  # coord_cartesian(ylim=c(-100, 100))+
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=1))+
  geom_text(data = tib, aes(label = paste0("slope = ", round(slope, 4),
                                           "\np.val = ", round(pval, 3),
                                           "\nR2 = ", round(rsq, 3))),
            x = Inf, y = Inf, vjust = 1.1, hjust = 1.1, inherit.aes = FALSE,
            size = 3.3)+
  facet_wrap(~ipbsr)

ggdraw(p3)+
  annotate("text", -Inf, Inf, hjust = -0.5, vjust = 1.5, label = "C", size = 16/.pt,
           fontface = "bold")

ggsave(here("output/figures", "09-regio_trends_dynvar_ndvi.png"),
         width = gg_width, height = gg_height*1.3)



# ************************************************************
# ------- Fixed global change drivers (09-3 + 09-4) -----------
# ************************************************************

# run predictions with variables fixed one by one
# based on 06-0-RF_and_predictions_v2

# cmic --------------------------------------------------------------------

cmic <- glc_proc_data()

# dplyr::select(new_ID1, longitude, latitude, Reference,
#               Cmic, tmean, prec, clay, nitrogen, phh2o, sand, 
#               soc, ndvi, elev, land_cover)


# RF model ----------------------------------------------------------------


model <- glc_rf_model()
predictors <- model$trainingData %>% names %>% .[-length(.)]
lc_lvls <- levels(model$trainingData$land_cover)


# Predictions all years ---------------------------------------------------

# non-fixed dynamic predictors
fclim_dp <- c("tmean", "prec")
fLC_dp <- c("ndvi", "land_cover")


# par setup
cl <- makeCluster(5, outfile = "")
registerDoParallel(cl)

# predict, save rds, (raster), png
# takes ~ 50min
foreach (iyear = 1992:2013, .packages = "glcpck") %dopar% {
  # iyear <- 1992
  
  for (tofix in c("clim", "LC")){
    # tofix <- "clim"
    
    fixvars <- if (tofix == "clim") fclim_dp else if (tofix == "LC") fLC_dp
    
    # either use iyear or starting year 1992
    st <- stack(c(map(predictors[!predictors %in% fixvars], ~glc_get_resamp(.x, iyear, "5-15")),
                  map(fixvars, ~glc_get_resamp(.x, 1992, "5-15"))))
    
    # rename
    names(st) <- c(predictors[!predictors %in% fixvars], fixvars)
    
    # as.data.frame (~ 30 sec, ~ 2min with na.rm = T)
    gridyear <- as.data.frame(st, xy = TRUE, na.rm = T) %>%
      mutate(pid = glc_pIDfromXY(x,y))
    
    # land_cover to text
    gridyear <- gridyear %>% 
      mutate(land_cover = glc_LC_num2chr(land_cover)) %>% 
      filter(land_cover %in% lc_lvls) %>% 
      mutate(land_cover = factor(land_cover))
    
    # mask
    fullmask <- glc_fullmask_df() %>% select(pid, mask)
    
    gridyear <- gridyear %>% left_join(fullmask, by = "pid") %>% 
      filter(is.na(mask)) %>% select(-mask)
    
    # predict2013 (100 sec)
    gridyear$pred <- predict(model, gridyear)
    
    saveRDS(gridyear, here("derived/09-pred_fixed_vars",
                           paste0("09-fixed_", tofix, "_", iyear, ".rds")))
    
    ggplot(gridyear, aes(x, y, fill = pred)) +
      borders()+
      geom_raster()+
      coord_fixed()+
      scale_fill_viridis_c()+
      annotate(geom = "label", x = -168, y = -62, label = iyear)
    
    ggsave(here("derived/09-pred_fixed_vars", paste0("09-fixed_", tofix, "_", iyear, ".png")),
           width = gg_width, height = gg_height)
    
  }
  
}

stopCluster(cl)


# compare all three models in their predictions ---------------------------

# Global

list.files(here("derived/09-pred_fixed_vars/savepoint"))

globs <- list.files(here("derived/09-pred_fixed_vars/savepoint"),
                    "glob_sum", full.names = TRUE)

glob_sum <- map_dfr(globs, readRDS)

ggplot(glob_sum, aes(year, cstock, color = mod_fix))+
  geom_line()+
  geom_smooth(method = "lm")

ggplot(glob_sum, aes(year, change, color = mod_fix, fill = mod_fix))+
  geom_line()+
  geom_smooth(method = "lm", alpha = 0.3)


(p <- ggplot(glob_sum, aes(ysca, perc_ch, color = mod_fix))+
    geom_line()+
    # geom_abline(slope = mod$coefficients[1], intercept = 0)+
    geom_smooth(method = "lm", formula = y ~ 0 + x)+
    scale_x_continuous(labels = function(x) x + 1992)+
    # ggtitle("Global percentage change in microbial carbon stocks")+
    labs(x = "Year", y = "Percentage change")
  # annotate("text", x = Inf, y = Inf, hjust=1.1, vjust=1.2, label = tpc)
)

ggsave(file.path("derived/09-pred_fixed_vars/figures",
                 paste0("global_fixed_", tofix, ".png")),
       width = gg_width*.6, height = gg_height*.8)


# regional ----------------------------------------------------------------

regs <- list.files(here("derived/09-pred_fixed_vars/savepoint"),
                   "reg_sum", full.names = TRUE)

reg_suma <- map_dfr(regs, readRDS)

reg_suma <- reg_suma %>% mutate(mod_fix = factor(mod_fix, levels = c("full", "clim", "LC")))

# library(forcats)

reg_suma <- reg_suma %>% mutate(mod_fix = fct_recode(mod_fix, "All" = "full",
                                                     "Only land cover" = "clim",
                                                     "Only climate" = "LC"))

p_smooth <- ggplot(reg_suma, aes(ysca, perc_ch, color = mod_fix, fill = mod_fix))+
  geom_line(size = 0.5)+
  geom_smooth(method = "lm", size = 0.5, 
              formula = y ~ 0 + x,
              alpha = 0.3)+
  scale_x_continuous(labels = function(x) x + 1992)+
  # ggtitle("Regional percentage changes in microbial carbon stocks")+
  labs(x = "Year", y = "Percentage change", color = "Free predictors", fill = "Free predictors")+
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=1))+
  facet_wrap(~ipbsr)


ggdraw(p_smooth)+
  draw_plot_label("B")

ggsave(here("derived/09-pred_fixed_vars/figures",
                 paste0("compare_all_smooth.png")),
       width = gg_width*1.2, height = gg_height*1.4)


# plot only slope

fct <- function(dat) {
  # dat <- reg_n$data[[1]]
  
  mod <- lm(perc_ch ~ year, data = dat)
  slope <- coef(mod)[2]
  stdErr <- summary(mod)$coef %>% as_tibble %>% pull(`Std. Error`) %>% .[2]
  pval <- coef(summary(mod))[2,4] 
  rsq <- summary(mod)$r.squared
  list(slope, stdErr, pval, rsq) %>% set_names(c("slope", "stdErr", "pval", "rsq"))
  
}

glob_df <- glob_sum %>% group_by(mod_fix) %>% nest %>% mutate(ipbsr = "Global")

glob_df <- glob_df %>% 
  mutate(mod_fix = factor(mod_fix)) %>% 
  mutate(mod_fix = fct_recode(mod_fix, "All" = "full",
                              "Only land cover" = "clim",
                              "Only climate" = "LC"))

reg_n <- reg_suma %>% group_by(ipbsr, mod_fix) %>% nest %>% ungroup

reg_n <- reg_n %>% bind_rows(glob_df)

reg_n <- reg_n %>% bind_cols(map_dfr(.$data, fct))

# reg_n <- reg_n %>% mutate(mod_fix = factor(mod_fix, levels = c("full", "clim", "LC")))

reg_n <- reg_n %>% mutate(cross0 = as.numeric(abs(slope) < 1.96*stdErr))
reg_n <- reg_n %>% mutate(signi = as.numeric(pval < 0.05))

# dada <- reg_n %>% filter(ipbsr == "Global\n", mod_fix == "Only climate")

reg_n %>% count(cross0, signi)
dada <- reg_n %>% filter(cross0 == signi) %>% .$data %>% .[[1]]

fct(dat = dada)
0.1456444 - 1.96*0.07211675
qplot(data = dada, year, cstock)
qplot(data = dada, year, perc_ch)

# Reorder regions per continent
main_order <- c("Global\n", "Africa\n", "Americas\n", "Asia-Pacific\n", "Europe-\nCentral Asia")

reg_order <- c("Global",
               "Central Africa", "East Africa and\nadjacent islands", "North Africa", "Southern Africa", "West Africa",
               "Caribbean", "Mesoamerica", "North America", "South America",
               "North-East Asia", "Oceania", "South-East Asia", "South Asia", "Western Asia",
               "Central and\nWestern Europe", "Central Asia", "Eastern Europe")

reg_n <- reg_n %>% 
  mutate(ipbsr = if_else(ipbsr == "East Africa and adjacent islands", "East Africa and\nadjacent islands", ipbsr)) %>%
  mutate(ipbsr = if_else(ipbsr == "Central and Western Europe", "Central and\nWestern Europe", ipbsr))

reg_n <- reg_n %>% mutate(main_region = case_when(
  ipbsr == "Global" ~ "Global\n",
  ipbsr %in% c("Central Africa", "East Africa and\nadjacent islands", "North Africa", "Southern Africa", "West Africa") ~ "Africa\n",
  ipbsr %in% c("Caribbean", "Mesoamerica", "North America", "South America") ~ "Americas\n",
  ipbsr %in% c("North-East Asia", "Oceania", "South-East Asia", "South Asia", "Western Asia") ~ "Asia-Pacific\n",
  ipbsr %in% c("Central and\nWestern Europe", "Central Asia", "Eastern Europe") ~ "Europe-\nCentral Asia") %>% 
    factor(levels = main_order)
)

reg_n <- reg_n %>% mutate(ipbsr = factor(ipbsr, levels = reg_order))

ggplot(reg_n, aes(mod_fix, slope, fill = mod_fix))+
  geom_col()+
  facet_wrap(~ipbsr)

ggplot(reg_n, aes(slope, mod_fix, color = mod_fix, alpha = cross0))+
  geom_vline(aes(xintercept = 0), color = "grey", linetype = 2)+
  geom_point()+ #shape = 15
  geom_errorbarh(aes(xmin = slope - stdErr, xmax = slope + stdErr), height = 0)+
  scale_y_discrete(limits = rev)+
  scale_alpha(range = c(0.5, 1))+
  facet_wrap(~ipbsr)+
  theme(legend.position = "none")

# All together
ggplot(reg_n, aes(slope, ipbsr, color = mod_fix, alpha = -signi))+
  geom_vline(aes(xintercept = 0), color = "grey", linetype = 2)+
  geom_point(position = position_nudge(y = rep(c(-0.2, 0.2, 0), each = 17)), size = 1.7)+ #shape = 15
  geom_errorbarh(aes(xmin = slope - 1.96*stdErr, xmax = slope + 1.96*stdErr), height = 0,
                 position = position_nudge(y = rep(c(-0.2, 0.2, 0), each = 17)),
                 size = 1)+
  scale_y_discrete(limits = rev)+
  scale_alpha(range = c(0.5, 1), guide = "none")+
  scale_color_discrete(name = "Free predictors",
                       labels = c(full = "All", clim = "Only land cover", LC = "Only climate"))+
  theme(panel.grid.major.y = element_line(size = 9, color = "#DDDDDD77"))+
  ylab("Europe-Central Asia |       Asia-Pacific       |       Americas      |       Africa     |   ")


# with facets
p1 <- ggplot(reg_n, aes(slope, ipbsr, color = mod_fix, alpha = -signi))+
  geom_vline(aes(xintercept = 0), color = "grey", linetype = 2)+
  geom_point(position = position_nudge(y = rep(c(-0.2, 0.2, 0), each = 17)), size = 1.7)+ #shape = 15
  geom_errorbarh(aes(xmin = slope - 1.96*stdErr, xmax = slope + 1.96*stdErr), height = 0,
                 position = position_nudge(y = rep(c(-0.2, 0.2, 0), each = 17)),
                 size = 1)+
  scale_y_discrete(limits = rev)+
  scale_alpha(range = c(0.5, 1), guide = "none")+
  scale_color_discrete(name = "Free predictors",
                       labels = c(full = "All", clim = "Only land cover", LC = "Only climate"))+
  theme(panel.grid.major.y = element_line(size = 9, color = "#DDDDDD77"))+
  facet_wrap(~main_region, scales = "free_y", ncol = 1, strip = "left")

leg <- p1 %>% get_legend

# split dfs
# reg_dfs <- split(reg_n, reg_n$main_region)

reg_n <- reg_n %>% mutate(mod_fix = factor(mod_fix, levels = c("All", "Only land cover", "Only climate")))

mk_facet <- function(region) {
  
  df_sub <- reg_n %>% filter(main_region == region)
  
  p <- ggplot(df_sub, aes(slope, ipbsr, color = mod_fix, alpha = -cross0, shape = mod_fix))+
    geom_vline(aes(xintercept = 0), color = "grey", linetype = 2)+
    geom_point(position = position_nudge(y = rep(c(-0.2, 0.2, 0), each = nrow(df_sub)/3)), size = 1.7)+ #shape = 15
    geom_errorbarh(aes(xmin = slope - 1.96*stdErr, xmax = slope + 1.96*stdErr), height = 0,
                   position = position_nudge(y = rep(c(-0.2, 0.2, 0), each = nrow(df_sub)/3)),
                   size = 1)+
    scale_y_discrete(limits = rev)+
    scale_alpha(range = c(0.5, 1), guide = "none")+
    scale_shape_discrete(name = "Free predictors",
                         labels = c(full = "All", clim = "Only land cover", LC = "Only climate"))+
    scale_color_brewer(name = "Free predictors",
                       labels = c(full = "All", clim = "Only land cover", LC = "Only climate"),
                       palette = "Dark2")+
    theme(panel.grid.major.y = element_line(size = 9, color = "#DDDDDD77"))+
    # facet_wrap(~main_region, scales = "free_y", ncol = 1, strip = "left")+ylab(NULL)+
    ylab(if(region == "Global\n") "" else region)+
    coord_cartesian(x = c(-0.73, 0.55))+
    xlab(NULL)+
    annotate(geom = 'segment', y = Inf, yend = -Inf, color = 'black', x = -Inf, xend = -Inf, size = 1)+
    annotate(geom = 'segment', y = Inf, yend = -Inf, color = 'black', x = Inf, xend = Inf, size = 0.5)+
    theme(legend.position = "none",
          panel.border = element_blank(),
          plot.margin = margin(0,0.1,-2.4,0.1),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 13)
    )
  # p
  
  if (region == "Global\n") {
    p + annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1)+
      theme(axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            plot.margin = margin(0.1,0.1,-2.45,0.1)
      )
  } else if (region == "Europe-\nCentral Asia") {
    p + 
      xlab("Rate of change (% per year)")+
      annotate(geom = 'segment', y = -Inf, yend = -Inf, color = 'black', x = -Inf, xend = Inf, size = 0.5)+
      theme(plot.margin = margin(0.1,0.1,0,0.1))
  } else {
    p + theme(axis.ticks.x = element_blank(),
              axis.text.x = element_blank())
  }
}

mk_facet(region = "Global\n")

leg <- (mk_facet(region = "Africa\n") +
          theme(legend.position = "bottom")+
          guides(shape = guide_legend(label.position = "top",
                                      title.position = "left", title.vjust = 0.5))) %>% 
  get_legend()


all_plots <- map(main_order, mk_facet)
all_plots[[3]]
all_plots[[5]]

for_plot <- plot_grid(leg, plot_grid(plotlist = all_plots, ncol = 1, rel_heights = c(1.3,5,4,5,4.3), align = "v"),
                      nrow = 2, rel_heights = c(0.08, 0.92)) #+

for_plot

ggsave(plot = for_plot,
       here("derived/09-pred_fixed_vars/figures",
                 paste0("compare_all_dots.png")),
       width = gg_width*0.5, height = gg_height*1)


# make table --------------------------------------------------------------

library(gt)

reg_n$data[[1]] %>% qplot(data = ., ysca, perc_ch)


# *slope is in relative terms
reg_n <- reg_n %>% mutate(cmst = map_dbl(data, ~ mean(.x$cstock)), #.x$cstock[1]    #cmic stock
                          abs_rate = slope * cmst * .01,
                          abs_tot = abs_rate * 21,
                          rel_tot = slope * 21)

reg_n %>% filter(ipbsr == "Global")

reg4tab <- reg_n %>% mutate(ci_95 = stdErr * 1.96) %>% 
  arrange(mod_fix, main_region) %>% 
  select(ipbsr, mod_fix, slope, ci_95,
         abs_rate, abs_tot, rel_tot, rsq)

gtab <- reg4tab %>% group_by(mod_fix) %>% gt

gtab_fmt <- gtab %>% 
  tab_spanner(
    label = html("Absolute (UNIT)"),
    columns = c(abs_rate, abs_tot)
  ) %>% 
  tab_spanner(
    label = html("Relative (% per year)"),
    columns = c(slope, ci_95, rel_tot)
  ) %>% 
  fmt_number(slope, decimals = 3) %>% 
  fmt_number(ci_95, decimals = 3) %>% 
  fmt_number(rsq, decimals = 3) %>% 
  fmt_number(abs_rate, decimals = 3) %>% 
  fmt_number(abs_tot, decimals = 3) %>% 
  cols_label(
    ipbsr = "IPBES region",
    mod_fix = "Free predictors",
    slope = "Rate of change",
    ci_95 = "95% CI",
    rsq = "R²",
    abs_rate = "Rate of change",
    abs_tot = "Total over 22 y",
    rel_tot = "Total over 22 y"
  ) %>% 
  tab_style(
    style = list(
      cell_text(weight = "bold")
    ),
    locations = cells_body(
      columns = c(slope, ci_95),
      rows = abs(slope) > ci_95
    )
  )

gtab_fmt

gtsave(gtab_fmt, here("output/tables", "09-4_fig3_table_reg_n.png"))



# ************************************************************
# ------------ Proportion cells predicted (09-5) ------------
# ************************************************************
#purpose: calculate percentage of cells calculated

iyear <- 2013

model <- glc_rf_model()
predictors <- model$trainingData %>% names %>% .[-length(.)]

# global layers
st <- stack(map(predictors, ~glc_get_resamp(.x, 2013, "5-15")))
rarea <- area(st)
st <- stack(st, rarea)
names(st) <- c(predictors, "area")
gridyear <- as.data.frame(st, xy = TRUE, na.rm = TRUE)
gridyear <- gridyear %>% mutate(pid = glc_pIDfromXY(x, y)) # pid = rownames(gridyear) or cellFromXY()
map_dbl(gridyear, ~sum(is.na(.x)))

# mask
fullmask <- glc_fullmask_df() %>% select(pid, mask)

# preds
p2013 <- readRDS(here("derived/06-prediction_years", "06-predictions_cmic_2013.rds")) %>% 
  select(pid, pred)

dat <- gridyear %>%
  left_join(fullmask, by = "pid") %>% 
  left_join(p2013, by = "pid")


# IPBES region ------------------------------------------------------------
dat <- dat %>% left_join(ip_df, by = "pid") #see code above for ip_df

# calc percentage ---------------------------------------------------------
# dat has all locations where predictor data was available

# percentage pixel coverage globally = 47.7 %
(dat$pred %>% {sum(!is.na(.))}) / nrow(dat) * 100

# perc area = 50.2 %
dat <- dat %>% mutate(predicted = !is.na(pred))
(dat %>% mutate(pred_a = predicted * area) %>% pull %>% sum) /
  sum(dat$area) * 100

ggplot(dat, aes(x, y, fill = pred))+ # elev
  geom_raster()+
  scale_fill_continuous(na.value = "grey80")+
  coord_fixed()


# perc covered from pixels with data. Nb. pixels
dat_sum <- dat %>% mutate(is_pred = !is.na(pred)) %>% 
  group_by(ipbsr) %>% 
  summarise(n = n(),
            preds = sum(!is.na(pred))) %>% 
  mutate(perc = preds / n * 100) %>% 
  arrange(perc)

dat_sum

# calculate area
dat_sumA <- dat %>% mutate(is_pred = !is.na(pred)) %>% 
  group_by(ipbsr, is_pred) %>% 
  summarise(n = n(),
            area = sum(area)) %>% 
  ungroup

dat_sA2 <- dat_sumA %>% select(-n) %>% 
  pivot_wider(names_from = is_pred, values_from = area) %>% 
  mutate(tot = `FALSE` + `TRUE`,
         perc_area = `TRUE` / tot *100) %>% 
  arrange(perc_area)


# combine both
dat_comb <- dat_sum %>% select(ipbsr, perc_n = perc, n_pix = n, pred_pix = preds) %>% 
  left_join(dat_sA2 %>% select(ipbsr, perc_area, tot_area = tot, pred_area = `TRUE`), by = "ipbsr")

qplot(data = dat_comb, perc_n, perc_area)+
  geom_label(aes(label = ipbsr))


# some discrepencies, esp. for northern regions
with(dat_comb, cbind(perc_n/perc_area, ipbsr)) %>% as_tibble %>% arrange(V1)


# mean Cstock per region --------------------------------------------------
# some repeated code from above, could be cleaned up

slo_df <- readRDS(here("derived/07-3-slope_per_pix_LM.rds"))
preds <- glc_predicted_dataset()

slo_df %>% wrnam
preds %>% wrnam

slo_df <- slo_df %>% filter(nona == 22)

preds <- preds %>% filter(pid %in% slo_df$pid)

# predictions nested
prenest <- preds %>% select(pid, x, y, year, stock, cell_stock) %>% group_by(pid) %>%  nest() %>% 
  mutate(nrow = map_dbl(data, ~ nrow(.x)))

prenest %>% pull(nrow) %>% table

prenest <- prenest %>% mutate(mcellst = map_dbl(data, ~ mean(.x$cell_stock, na.rm = TRUE))) %>% 
  select(-data)


# match ipbes subregion
ipras <- raster(here("geodata/09-ipbes_subregions_alphaOrder.tif"))
names(ipras) <- "ipbco"

levs <- c("Antarctica", "Caribbean", "Central Africa", "Central and Western Europe", 
          "Central Asia", "East Africa and adjacent islands", "Eastern Europe", 
          "Mesoamerica", "North-East Asia", "North Africa", "North America", 
          "Oceania", "South-East Asia", "South America", "South Asia", 
          "Southern Africa", "West Africa", "Western Asia")

match_tbl <- tibble(ipbsr = levs,
                    code = 1:length(ipbsr))

ip_df <- as.data.frame(ipras, xy = TRUE, na.rm = TRUE)
ip_df <- ip_df %>% 
  left_join(match_tbl, by = c("ipbco" = "code")) %>% 
  mutate(pid = glc_pIDfromXY(x, y)) %>% 
  select(-c(x,y))

prenest <- prenest %>% left_join(ip_df, by = "pid")
prenest <- prenest %>% filter(!is.na(ipbsr))


reg_summ <- prenest %>% group_by(ipbsr) %>% 
  summarise(cstock = sum(mcellst))

dat_comb <- dat_comb %>% left_join(reg_summ, by = "ipbsr")


# make table --------------------------------------------------------------
cat(names(dat_comb), sep = ' = "",\n')

dat_comb$pred_pix %>% sum

tab <- dat_comb %>% drop_na(ipbsr) %>%
  arrange(ipbsr) %>% 
  select(ipbsr, n_pix, pred_pix, perc_n, tot_area, pred_area, perc_area, cstock) %>% 
  gt(rowname_col = "ipbsr") %>% 
  tab_stubhead(label = "Region") %>% 
  tab_spanner(
    label = "Pixels",
    columns = c(n_pix, pred_pix, perc_n)
  ) %>%
  tab_spanner(
    label = html("Area (km<sup>2</sup>)"),
    columns = c(tot_area, pred_area, perc_area)
  ) %>% 
  cols_label(
    perc_n = "%",
    n_pix = "total",
    pred_pix = "predicted",
    perc_area = "%",
    tot_area = "total",
    pred_area ="predicted",
    cstock = html("Microbial carbon<br>stock (t)")
  ) %>% 
  fmt_number(
    columns = c(perc_n, perc_area),
    decimals = 1,
    use_seps = FALSE
  ) %>% 
  fmt_number(
    columns = c(pred_pix, n_pix, tot_area, pred_area, cstock),
    decimals = 0,
    use_seps = TRUE
  )

tab
gtsave(tab, here("output/tables", "09-5-perc_cells_pred.png"))
