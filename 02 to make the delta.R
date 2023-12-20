

# Fabio Castro - Llanos 
# Alliance Bioversity - CIAT 
# December 2023

# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(terra, fs, sf, tidyverse, cptcity, RColorBrewer, rgeos, gtools, glue, rnaturalearthdata, rnaturalearth)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Functions to use --------------------------------------------------------

# Load data ---------------------------------------------------------------
fles <- dir_ls('./data/tif/cmip6/30s', regexp = '.tif$')
fles <- as.character(fles)

# Future data
prec.ftre <- grep('prec', fles, value = T)
prec.ftre <- rast(prec.ftre)

tmax.ftre <- grep('tmax', fles, value = T)
tmax.ftre <- rast(tmax.ftre)

tmin.ftre <- grep('tmin', fles, value = T)
tmin.ftre <- rast(tmin.ftre)

# Baseline data
prec.bsln <- geodata::worldclim_country(country = 'COL', var = 'prec', path = 'tmpr')
tmin.bsln <- geodata::worldclim_country(country = 'COL', var = 'tmin', path = 'tmpr')
tmax.bsln <- geodata::worldclim_country(country = 'COL', var = 'tmax', path = 'tmpr')

# Vector data 
wrld <- ne_countries(returnclass = 'sf', scale = 50)
limt <- wrld[wrld$sov_a3 == 'COL',]
limt <- vect(limt)

# To extract by mask  -----------------------------------------------------

# Baseline
prec.bsln <- terra::crop(prec.bsln, limt) %>% terra::mask(., limt)
tmin.bsln <- terra::crop(tmin.bsln, limt) %>% terra::mask(., limt)
tmax.bsln <- tmax.bsln %>% terra::crop(., limt) %>% terra::mask(., limt)

# Future
prec.ftre <- terra::crop(prec.ftre, limt) %>% terra::mask(., limt)
tmin.ftre <- terra::crop(tmin.ftre, limt) %>% terra::mask(., limt)
tmax.ftre <- terra::mask(tmax.ftre, limt) %>% terra::mask(., limt)

# To calculate the delta --------------------------------------------------

prec.ftre <- sum(prec.ftre)
prec.bsln <- sum(prec.bsln)

tmin.ftre <- mean(tmin.ftre)
tmin.bsln <- mean(tmin.bsln)

tmax.ftre <- mean(tmax.ftre)
tmax.bsln <- mean(tmax.bsln)

prec.dfrn <- prec.ftre - prec.bsln
prec.dfrn <- prec.dfrn / prec.bsln * 100

tmin.dfrn <- tmin.ftre - tmin.bsln
tmax.dfrn <- tmax.ftre - tmax.bsln

# To make the stack -------------------------------------------------------
prec.stck <- c(prec.bsln, prec.ftre)
tmin.stck <- c(tmin.bsln, tmin.ftre)
tmax.stck <- c(tmax.bsln, tmax.ftre)

prec.tble <- terra::as.data.frame(prec.stck, xy = T) %>% setNames(c('x', 'y', 'bsl', 'ftr')) %>% as_tibble() 
tmin.tble <- terra::as.data.frame(tmin.stck, xy = T) %>% setNames(c('x', 'y', 'bsl', 'ftr')) %>% as_tibble() 
tmax.tble <- terra::as.data.frame(tmax.stck, xy = T) %>% setNames(c('x', 'y', 'bsl', 'ftr')) %>% as_tibble() 

prec.tble <- mutate(prec.tble, dfr = (ftr - bsl) / bsl * 100)
tmin.tble <- mutate(tmin.tble, dfr = (ftr - bsl))
tmax.tble <- mutate(tmax.tble, dfr = ftr - bsl)

# To make the map  --------------------------------------------------------

# Minimum temperature -----------------------------------------------------
tmin.prds <- tmin.tble %>% 
  dplyr::select(x, y, bsl, ftr) %>% 
  gather(period, value, -x, -y) %>% 
  mutate(period = ifelse(period == 'bsl', 'Línea base', 'Futuro (2050)'), 
         period = factor(period, levels = c('Línea base', 'Futuro (2050)'))) 

gtmin <- ggplot() + 
  geom_tile(data = tmin.prds, aes(x = x, y = y, fill = value)) + 
  facet_wrap(~period) +
  scale_fill_gradientn(colors = brewer.pal(n = 9, name = 'YlOrRd'), 
                       labels = seq(0, 30, 5), 
                       breaks = seq(0, 30, 5),
                       name = 'Temperatura (°C)') +
  geom_sf(data = wrld, fill = NA, col = 'grey30') +
  geom_sf(data = st_as_sf(limt), fill = NA, col = 'grey40') +
  coord_sf(xlim = ext(limt)[1:2], ylim = ext(limt)[3:4]) +
  labs(x = 'Lon', y = 'Lat', fill = 'Temperatura (°C)', caption = 'Datos de línea base: Worldclim v2.1 (Hijmans et al., 2014)\nDatos de futuro: CMIP6 - SSP 370 Ensemble de 15 GCMs') +
  theme_light() + 
  theme(legend.position = 'bottom',
        axis.text.y = element_text(angle = 90, hjust = 0.5),
        axis.title = element_text(face = 'bold'),
        plot.caption = element_text(hjust = 0.5),
        legend.key.width = unit(3, 'line')) + 
  guides(fill = guide_legend( 
    direction = 'horizontal',
    keyheight = unit(1.15, units = "mm"),
    keywidth = unit(15, units = "mm"),
    title.position = 'top',
    title.hjust = 0.5,
    label.hjust = .5,
    nrow = 1,
    byrow = T,
    reverse = F,
    label.position = "bottom"
  )) +
  annotation_scale(location =  "bl", width_hint = 0.5, text_family = 'Segoe UI', text_col = 'grey60', bar_cols = c('grey60', 'grey99'), line_width = 0.2) +
  annotation_north_arrow(location = "tr", which_north = "true", 
                         pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"), 
                         style = north_arrow_fancy_orienteering(text_family = 'Segoe UI', text_col = 'grey40', line_col = 'grey60', fill = c('grey60', 'grey99'))) 

ggsave(plot = gtmin, filename = glue('./png/maps/tmin_periods.png'), units = 'in', width = 9, height = 7, dpi = 300)

tmin.dfrn <- tmin.tble %>% 
  dplyr::select(x, y, dfr) %>% 
  gather(var, value, -x, -y) %>% 
  mutate(var = ifelse(var == 'dfr', 'Delta', NA)) 

gtmin.dfrn <- ggplot() + 
  geom_tile(data = tmin.dfrn, aes(x = x, y = y, fill = value)) + 
  facet_wrap(~var) +
  scale_fill_gradientn(colors = brewer.pal(n = 9, name = 'YlOrRd'), 
                       # labels = seq(1.5, 2.5, 2), 
                       # breaks = seq(1.5, 2.5, 2),
                       name = 'Temperatura (°C)') +
  geom_sf(data = wrld, fill = NA, col = 'grey30') +
  geom_sf(data = st_as_sf(limt), fill = NA, col = 'grey40') +
  coord_sf(xlim = ext(limt)[1:2], ylim = ext(limt)[3:4]) +
  labs(x = 'Lon', y = 'Lat', fill = 'Temperatura (°C)') +
  theme_light() + 
  theme(legend.position = 'bottom',
        axis.text.y = element_text(angle = 90, hjust = 0.5),
        axis.title = element_text(face = 'bold'),
        plot.caption = element_text(hjust = 0.5),
        legend.key.width = unit(3, 'line')) + 
  guides(fill = guide_legend( 
    direction = 'horizontal',
    keyheight = unit(1.15, units = "mm"),
    keywidth = unit(15, units = "mm"),
    title.position = 'top',
    title.hjust = 0.5,
    label.hjust = .5,
    nrow = 1,
    byrow = T,
    reverse = F,
    label.position = "bottom"
  )) +
  annotation_scale(location =  "bl", width_hint = 0.5, text_family = 'Segoe UI', text_col = 'grey60', bar_cols = c('grey60', 'grey99'), line_width = 0.2) +
  annotation_north_arrow(location = "tr", which_north = "true", 
                         pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"), 
                         style = north_arrow_fancy_orienteering(text_family = 'Segoe UI', text_col = 'grey40', line_col = 'grey60', fill = c('grey60', 'grey99'))) 
gtmin.dfrn

ggsave(plot = gtmin.dfrn, filename = './png/maps/tmin_delta.png', units = 'in', width = 6, height = 9, dpi = 300)

gtmin.all <- ggpubr::ggarrange(gtmin, gtmin.dfrn, ncol = 2, nrow = 1)
ggsave(plot = gtmin.all, filename = './png/maps/tmin_vals-delta.png', units = 'in', width = 13, height = 7, dpi = 300)

# Maximum temperature -----------------------------------------------------
tmax.prds <- tmax.tble %>% 
  dplyr::select(x, y, bsl, ftr) %>% 
  gather(period, value, -x, -y) %>% 
  mutate(period = ifelse(period == 'bsl', 'Línea base', 'Futuro (2050)'), 
         period = factor(period, levels = c('Línea base', 'Futuro (2050)'))) 

gtmax <- ggplot() + 
  geom_tile(data = tmax.prds, aes(x = x, y = y, fill = value)) + 
  facet_wrap(~period) +
  scale_fill_gradientn(colors = brewer.pal(n = 9, name = 'YlOrRd'), 
                       labels = seq(0, 30, 5), 
                       breaks = seq(0, 30, 5),
                       name = 'Temperatura (°C)') +
  geom_sf(data = wrld, fill = NA, col = 'grey30') +
  geom_sf(data = st_as_sf(limt), fill = NA, col = 'grey40') +
  coord_sf(xlim = ext(limt)[1:2], ylim = ext(limt)[3:4]) +
  labs(x = 'Lon', y = 'Lat', fill = 'Temperatura (°C)', caption = 'Datos de línea base: Worldclim v2.1 (Hijmans et al., 2014)\nDatos de futuro: CMIP6 - SSP 370 Ensemble de 15 GCMs') +
  theme_light() + 
  theme(legend.position = 'bottom',
        axis.text.y = element_text(angle = 90, hjust = 0.5),
        axis.title = element_text(face = 'bold'),
        plot.caption = element_text(hjust = 0.5),
        legend.key.width = unit(3, 'line')) + 
  guides(fill = guide_legend( 
    direction = 'horizontal',
    keyheight = unit(1.15, units = "mm"),
    keywidth = unit(15, units = "mm"),
    title.position = 'top',
    title.hjust = 0.5,
    label.hjust = .5,
    nrow = 1,
    byrow = T,
    reverse = F,
    label.position = "bottom"
  )) +
  annotation_scale(location =  "bl", width_hint = 0.5, text_family = 'Segoe UI', text_col = 'grey60', bar_cols = c('grey60', 'grey99'), line_width = 0.2) +
  annotation_north_arrow(location = "tr", which_north = "true", 
                         pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"), 
                         style = north_arrow_fancy_orienteering(text_family = 'Segoe UI', text_col = 'grey40', line_col = 'grey60', fill = c('grey60', 'grey99'))) 

ggsave(plot = gtmax, filename = glue('./png/maps/tmax_periods.png'), units = 'in', width = 9, height = 7, dpi = 300)

tmax.dfrn <- tmax.tble %>% 
  dplyr::select(x, y, dfr) %>% 
  gather(var, value, -x, -y) %>% 
  mutate(var = ifelse(var == 'dfr', 'Delta', NA)) 

gtmax.dfrn <- ggplot() + 
  geom_tile(data = tmax.dfrn, aes(x = x, y = y, fill = value)) + 
  facet_wrap(~var) +
  scale_fill_gradientn(colors = brewer.pal(n = 9, name = 'YlOrRd'), 
                       # labels = seq(1.5, 2.5, 2), 
                       # breaks = seq(1.5, 2.5, 2),
                       name = 'Temperatura (°C)') +
  geom_sf(data = wrld, fill = NA, col = 'grey30') +
  geom_sf(data = st_as_sf(limt), fill = NA, col = 'grey40') +
  coord_sf(xlim = ext(limt)[1:2], ylim = ext(limt)[3:4]) +
  labs(x = 'Lon', y = 'Lat', fill = 'Temperatura (°C)') +
  theme_light() + 
  theme(legend.position = 'bottom',
        axis.text.y = element_text(angle = 90, hjust = 0.5),
        axis.title = element_text(face = 'bold'),
        plot.caption = element_text(hjust = 0.5),
        legend.key.width = unit(3, 'line')) + 
  guides(fill = guide_legend( 
    direction = 'horizontal',
    keyheight = unit(1.15, units = "mm"),
    keywidth = unit(15, units = "mm"),
    title.position = 'top',
    title.hjust = 0.5,
    label.hjust = .5,
    nrow = 1,
    byrow = T,
    reverse = F,
    label.position = "bottom"
  )) +
  annotation_scale(location =  "bl", width_hint = 0.5, text_family = 'Segoe UI', text_col = 'grey60', bar_cols = c('grey60', 'grey99'), line_width = 0.2) +
  annotation_north_arrow(location = "tr", which_north = "true", 
                         pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"), 
                         style = north_arrow_fancy_orienteering(text_family = 'Segoe UI', text_col = 'grey40', line_col = 'grey60', fill = c('grey60', 'grey99'))) 
gtmax.dfrn

ggsave(plot = gtmax.dfrn, filename = './png/maps/tmax_delta.png', units = 'in', width = 6, height = 9, dpi = 300)

gtmax.all <- ggpubr::ggarrange(gtmax, gtmax.dfrn, ncol = 2, nrow = 1)
ggsave(plot = gtmax.all, filename = './png/maps/tmax_vals-delta.png', units = 'in', width = 13, height = 7, dpi = 300)

# Average temperature -----------------------------------------------------

tmin.tble
tmax.tble
colnames(tmin.tble) <- c('x', 'y', 'bsl_tmin', 'ftr_tmin', 'dfr_tmin')
colnames(tmax.tble) <- c('x', 'y', 'bsl_tmax', 'ftr_tmax', 'dfr_tmax')

tasm <- inner_join(tmin.tble, tmax.tble)
tasm <- mutate(tasm, bsl_tavg = (bsl_tmin + bsl_tmax) / 2, ftr_tavg = (ftr_tmin + ftr_tmax) / 2)
tasm <- mutate(tasm, dfr_tavg = ftr_tavg - bsl_tavg)

tasm.prds <- tasm %>% 
  dplyr::select(x, y, bsl_tavg, ftr_tavg) %>% 
  gather(period, value, -x, -y) %>% 
  mutate(period = ifelse(period == 'bsl_tavg', 'Línea base', 'Futuro (2050)'), 
         period = factor(period, levels = c('Línea base', 'Futuro (2050)'))) 

gtasm <- ggplot() + 
  geom_tile(data = tasm.prds, aes(x = x, y = y, fill = value)) + 
  facet_wrap(~period) +
  scale_fill_gradientn(colors = brewer.pal(n = 9, name = 'YlOrRd'), 
                       labels = seq(0, 30, 5), 
                       breaks = seq(0, 30, 5),
                       name = 'Temperatura (°C)') +
  geom_sf(data = wrld, fill = NA, col = 'grey30') +
  geom_sf(data = st_as_sf(limt), fill = NA, col = 'grey40') +
  coord_sf(xlim = ext(limt)[1:2], ylim = ext(limt)[3:4]) +
  labs(x = 'Lon', y = 'Lat', fill = 'Temperatura (°C)', caption = 'Datos de línea base: Worldclim v2.1 (Hijmans et al., 2014)\nDatos de futuro: CMIP6 - SSP 370 Ensemble de 15 GCMs') +
  theme_light() + 
  theme(legend.position = 'bottom',
        axis.text.y = element_text(angle = 90, hjust = 0.5),
        axis.title = element_text(face = 'bold'),
        plot.caption = element_text(hjust = 0.5),
        legend.key.width = unit(3, 'line')) + 
  guides(fill = guide_legend( 
    direction = 'horizontal',
    keyheight = unit(1.15, units = "mm"),
    keywidth = unit(15, units = "mm"),
    title.position = 'top',
    title.hjust = 0.5,
    label.hjust = .5,
    nrow = 1,
    byrow = T,
    reverse = F,
    label.position = "bottom"
  )) +
  annotation_scale(location =  "bl", width_hint = 0.5, text_family = 'Segoe UI', text_col = 'grey60', bar_cols = c('grey60', 'grey99'), line_width = 0.2) +
  annotation_north_arrow(location = "tr", which_north = "true", 
                         pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"), 
                         style = north_arrow_fancy_orienteering(text_family = 'Segoe UI', text_col = 'grey40', line_col = 'grey60', fill = c('grey60', 'grey99'))) 


ggsave(plot = gtasm, filename = glue('./png/maps/tavg_periods.png'), units = 'in', width = 9, height = 7, dpi = 300)

tavg.dfrn <- tasm %>% 
  dplyr::select(x, y, dfr_tavg) %>% 
  gather(var, value, -x, -y) %>% 
  mutate(var = ifelse(var == 'dfr_tavg', 'Delta', NA)) 

gtasm.dfrn <- ggplot() + 
  geom_tile(data = tavg.dfrn, aes(x = x, y = y, fill = value)) + 
  facet_wrap(~var) +
  scale_fill_gradientn(colors = brewer.pal(n = 9, name = 'YlOrRd'), 
                       # labels = seq(1.5, 2.5, 2), 
                       # breaks = seq(1.5, 2.5, 2),
                       name = 'Temperatura (°C)') +
  geom_sf(data = wrld, fill = NA, col = 'grey30') +
  geom_sf(data = st_as_sf(limt), fill = NA, col = 'grey40') +
  coord_sf(xlim = ext(limt)[1:2], ylim = ext(limt)[3:4]) +
  labs(x = 'Lon', y = 'Lat', fill = 'Temperatura (°C)') +
  theme_light() + 
  theme(legend.position = 'bottom',
        axis.text.y = element_text(angle = 90, hjust = 0.5),
        axis.title = element_text(face = 'bold'),
        plot.caption = element_text(hjust = 0.5),
        legend.key.width = unit(3, 'line')) + 
  guides(fill = guide_legend( 
    direction = 'horizontal',
    keyheight = unit(1.15, units = "mm"),
    keywidth = unit(15, units = "mm"),
    title.position = 'top',
    title.hjust = 0.5,
    label.hjust = .5,
    nrow = 1,
    byrow = T,
    reverse = F,
    label.position = "bottom"
  )) +
  annotation_scale(location =  "bl", width_hint = 0.5, text_family = 'Segoe UI', text_col = 'grey60', bar_cols = c('grey60', 'grey99'), line_width = 0.2) +
  annotation_north_arrow(location = "tr", which_north = "true", 
                         pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"), 
                         style = north_arrow_fancy_orienteering(text_family = 'Segoe UI', text_col = 'grey40', line_col = 'grey60', fill = c('grey60', 'grey99'))) 
gtasm.dfrn

ggsave(plot = gtasm.dfrn, filename = './png/maps/tavg_delta.png', units = 'in', width = 6, height = 9, dpi = 300)

gtasm.all <- ggpubr::ggarrange(gtasm, gtasm.dfrn, ncol = 2, nrow = 1)
ggsave(plot = gtasm.all, filename = './png/maps/tavg_vals-delta.png', units = 'in', width = 13, height = 7, dpi = 300)

# Precipitation -----------------------------------------------------------

prec.prds <- prec.tble %>% 
  dplyr::select(x, y, bsl, ftr) %>% 
  gather(period, value, -x, -y) %>% 
  mutate(period = ifelse(period == 'bsl', 'Línea base', 'Futuro (2050)'), 
         period = factor(period, levels = c('Línea base', 'Futuro (2050)'))) 

gprec <- ggplot() + 
  geom_tile(data = prec.prds, aes(x = x, y = y, fill = value)) + 
  facet_wrap(~period) +
  scale_fill_gradientn(colors = brewer.pal(n = 9, name = 'BrBG'), 
                       labels = seq(200, 8800, 800), 
                       breaks = seq(200, 8800, 800),
                       name = 'Precipitación (mm)') +
  geom_sf(data = wrld, fill = NA, col = 'grey30') +
  geom_sf(data = st_as_sf(limt), fill = NA, col = 'grey40') +
  coord_sf(xlim = ext(limt)[1:2], ylim = ext(limt)[3:4]) +
  labs(x = 'Lon', y = 'Lat', fill = 'Precipitación (mm)', caption = 'Datos de línea base: Worldclim v2.1 (Hijmans et al., 2014)\nDatos de futuro: CMIP6 - SSP 370 Ensemble de 15 GCMs') +
  theme_light() + 
  theme(legend.position = 'bottom',
        axis.text.y = element_text(angle = 90, hjust = 0.5),
        axis.title = element_text(face = 'bold'),
        plot.caption = element_text(hjust = 0.5),
        legend.key.width = unit(3, 'line')) + 
  guides(fill = guide_legend( 
    direction = 'horizontal',
    keyheight = unit(1.15, units = "mm"),
    keywidth = unit(15, units = "mm"),
    title.position = 'top',
    title.hjust = 0.5,
    label.hjust = .5,
    nrow = 1,
    byrow = T,
    reverse = F,
    label.position = "bottom"
  )) +
  annotation_scale(location =  "bl", width_hint = 0.5, text_family = 'Segoe UI', text_col = 'grey60', bar_cols = c('grey60', 'grey99'), line_width = 0.2) +
  annotation_north_arrow(location = "tr", which_north = "true", 
                         pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"), 
                         style = north_arrow_fancy_orienteering(text_family = 'Segoe UI', text_col = 'grey40', line_col = 'grey60', fill = c('grey60', 'grey99'))) 

ggsave(plot = gprec, filename = glue('./png/maps/prec_periods.png'), units = 'in', width = 9, height = 7, dpi = 300)

prec.dfrn <- prec.tble %>% 
  dplyr::select(x, y, dfr) %>% 
  gather(var, value, -x, -y) %>% 
  mutate(var = ifelse(var == 'dfr', 'Delta', NA)) 

gprec.dfrn <- ggplot() + 
  geom_tile(data = prec.dfrn, aes(x = x, y = y, fill = value)) + 
  facet_wrap(~var) +
  scale_fill_gradientn(colors = brewer.pal(n = 9, name = 'BrBG'), 
                       labels = seq(-10, 8, 4),
                       breaks = seq(-10, 8, 4),
                       name = 'Precipitación (%)') +
  geom_sf(data = wrld, fill = NA, col = 'grey30') +
  geom_sf(data = st_as_sf(limt), fill = NA, col = 'grey40') +
  coord_sf(xlim = ext(limt)[1:2], ylim = ext(limt)[3:4]) +
  labs(x = 'Lon', y = 'Lat', fill = 'Precipitación (%)') +
  theme_light() + 
  theme(legend.position = 'bottom',
        axis.text.y = element_text(angle = 90, hjust = 0.5),
        axis.title = element_text(face = 'bold'),
        plot.caption = element_text(hjust = 0.5),
        legend.key.width = unit(3, 'line')) + 
  guides(fill = guide_legend( 
    direction = 'horizontal',
    keyheight = unit(1.15, units = "mm"),
    keywidth = unit(15, units = "mm"),
    title.position = 'top',
    title.hjust = 0.5,
    label.hjust = .5,
    nrow = 1,
    byrow = T,
    reverse = F,
    label.position = "bottom"
  )) +
  annotation_scale(location =  "bl", width_hint = 0.5, text_family = 'Segoe UI', text_col = 'grey60', bar_cols = c('grey60', 'grey99'), line_width = 0.2) +
  annotation_north_arrow(location = "tr", which_north = "true", 
                         pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"), 
                         style = north_arrow_fancy_orienteering(text_family = 'Segoe UI', text_col = 'grey40', line_col = 'grey60', fill = c('grey60', 'grey99'))) 
gprec.dfrn

ggsave(plot = gprec.dfrn, filename = './png/maps/prec_delta.png', units = 'in', width = 6, height = 9, dpi = 300)

gprec.all <- ggpubr::ggarrange(gprec, gprec.dfrn, ncol = 2, nrow = 1)
ggsave(plot = gprec.all, filename = './png/maps/send/prec_vals-delta.png', units = 'in', width = 13, height = 7, dpi = 300)



