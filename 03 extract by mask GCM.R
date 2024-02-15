



# Fabio Castro - Llanos 
# Alliance Bioversity - CIAT 
# January 2024

# Load libraries ----------------------------------------------------------
require(pacman)
pacman::p_load(terra, fs, sf, tidyverse, RColorBrewer, rgeos, gtools, glue, rnaturalearthdata, rnaturalearth)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Load data ---------------------------------------------------------------
wrld <- ne_countries(scale = 50, returnclass = 'sf')
col  <- wrld[wrld$sov_a3 == 'COL',]

# Raster data
root <- '//ALLIANCEDFS.ALLIANCE.CGIAR.ORG/data_cluster19/GLOBAL/climate/Agroclimas/data/ipcc_6ar_wcl_downscaled'
ssps <- dir_ls(root)
dirs <- dir_ls(ssps[3]) %>% grep('2050s', ., value = T) %>% as.character() %>% dir_ls() 
dirs


# Future ------------------------------------------------------------------

## Temperature -------------------------------------------------------------
fles.tasm <- as.character(dir_ls(dir_ls(as.character(grep('ec_earth3', dirs, value = T)))))
fles.tasm <- grep('_tm', fles.tasm, value = T)
rstr.tmin <- rast(grep('_tmin', fles.tasm, value = T))
rstr.tmin <- mask(crop(rstr.tmin, col), col)
rstr.tmin <- mean(rstr.tmin)
rstr.tmax <- rast(grep('_tmax', fles.tasm, value = T))
rstr.tmax <- mask(crop(rstr.tmax, col), col)
rstr.tmax <- mean(rstr.tmax)

### To calculate average temperature 
rstr.tavg <- (rstr.tmax + rstr.tmin) / 2
tavg.ftre <- rstr.tavg

## Precipitation -----------------------------------------------------------
fles.prec <- as.character(dir_ls(dir_ls(as.character(grep('ec_earth3', dirs, value = T)))))
fles.prec <- grep('_prec_', fles.prec, value = T)
rstr.prec <- rast(grep('_prec_', fles.prec, value = T))
rstr.prec <- mask(crop(rstr.prec, col), col)
rstr.prec <- sum(rstr.prec)
prec.ftre <- rstr.prec

# Baseline ----------------------------------------------------------------
prec.bsln <- geodata::worldclim_country(country = 'COL', var = 'prec', path = 'tmpr')
tmin.bsln <- geodata::worldclim_country(country = 'COL', var = 'tmin', path = 'tmpr')
tmax.bsln <- geodata::worldclim_country(country = 'COL', var = 'tmax', path = 'tmpr')

## Extract by mask baseline ---------
prec.bsln <- terra::crop(prec.bsln, col)
prec.bsln <- terra::mask(prec.bsln, col)
tmin.bsln <- terra::crop(tmin.bsln, col)
tmin.bsln <- terra::mask(tmin.bsln, col)
tmax.bsln <- terra::crop(tmax.bsln, col)
tmax.bsln <- terra::mask(tmax.bsln, col)

prec.bsln <- sum(prec.bsln)
tmax.bsln <- mean(tmax.bsln)
tmin.bsln <- mean(tmin.bsln)
tavg.bsln <- (tmax.bsln + tmin.bsln) / 2

# To calculate the delta --------------------------------------------------
prec.dfrn <- ((prec.ftre - prec.bsln) /prec.bsln) * 100
tavg.dfrn <- tavg.ftre - tavg.bsln

# To make the stack -------------------------------------------------------
stck <- c(
  prec.bsln, prec.ftre, prec.dfrn, tavg.bsln, tavg.ftre, tavg.dfrn
)

names(stck) <- c(
  'prec.bsln', 'prec.ftre', 'prec.dfrn', 'tavg.bsln', 'tavg.ftre', 'tavg.dfrn'
)

## Raster to table 
tble <- stck %>% 
  terra::as.data.frame(
    xy = T
  ) %>% 
  as_tibble() %>% 
  mutate(
    gid = 1:nrow(.)
  ) %>% 
  gather(period, value, -c(gid, x, y)) %>% 
  mutate(
    period = factor(
      period, levels = c('prec.bsln', 'prec.ftre', 'prec.dfrn', 'tavg.bsln', 'tavg.ftre', 'tavg.dfrn')
    )
  )

# To make the map ---------------------------------------------------------

## Precipitation ---------------------------------------------------------
tble.prec <- tble %>% 
  filter(
    period %in% c('prec.bsln', 'prec.ftre', 'prec.dfrn')
  ) %>% 
  mutate(
    period = ifelse(
      period == 'prec.bsln',
      'Línea base', 
      ifelse(
        period == 'prec.ftre',
        'Futuro', 
        'Diferencia'
      )
    ),
    period = factor(
      period, 
      levels = c(
        'Línea base', 'Futuro', 'Diferencia'
      )
    )
  )

# Periods
gprec.prds <- ggplot() + 
  geom_tile(data = tble.prec %>% 
              filter(period %in% c('Línea base', 'Futuro')), 
            aes(
              x = x, 
              y = y, 
              fill = value)
            ) + 
  facet_wrap(~period) +
  scale_fill_gradientn(colors = brewer.pal(n = 9, name = 'BrBG'), 
                       labels = seq(200, 8800, 800), 
                       breaks = seq(200, 8800, 800),
                       name = 'Precipitación (mm)') +
  geom_sf(data = wrld, fill = NA, col = 'grey30') +
  geom_sf(data = st_as_sf(col), fill = NA, col = 'grey40') +
  coord_sf(xlim = ext(col)[1:2], ylim = ext(col)[3:4]) +
  labs(x = 'Lon', y = 'Lat', fill = 'Precipitación (mm)', caption = 'Datos de línea base: Worldclim v2.1 (Hijmans et al., 2014)\nDatos de futuro: CMIP6 - SSP 370 GCM: Ec-Earth3-Veg-Lr') +
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

gprec.prds
ggsave(plot = gprec.prds, filename = './png/maps/v2/prec_prds.png', units = 'in', width = 11, height = 9, dpi = 300)

# Difference
prec.dfrn <- tble.prec %>% 
  filter(period == 'Diferencia') %>% 
  mutate(period = ifelse(period == 'Diferencia', 'Delta', NA)) 

gprec.dfrn <- ggplot() + 
  geom_tile(data = prec.dfrn, aes(x = x, y = y, fill = value)) + 
  facet_wrap(~period) +
  scale_fill_gradientn(colors = brewer.pal(n = 9, name = 'BrBG'), 
                       labels = seq(-10, 16, 4),
                       breaks = seq(-10, 16, 4),
                       name = 'Precipitación (%)') +
  geom_sf(data = wrld, fill = NA, col = 'grey30') +
  geom_sf(data = st_as_sf(col), fill = NA, col = 'grey40') +
  coord_sf(xlim = ext(col)[1:2], ylim = ext(col)[3:4]) +
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

gprec.all <- ggpubr::ggarrange(gprec.prds, gprec.dfrn, ncol = 2, nrow = 1)
ggsave(plot = gprec.all, filename = './png/maps/v2/prec_vals-delta_ssp370.png', units = 'in', width = 13, height = 7, dpi = 300)

## Temperature -----------------------------------------------------------
tble.tasm <- tble %>% 
  filter(
    period %in% c('tavg.bsln', 'tavg.ftre', 'tavg.dfrn')
  ) %>% 
  mutate(
    period = ifelse(
      period == 'tavg.bsln',
      'Línea base', 
      ifelse(
        period == 'tavg.ftre',
        'Futuro', 
        'Diferencia'
      )
    ),
    period = factor(
      period, 
      levels = c(
        'Línea base', 'Futuro', 'Diferencia'
      )
    )
  )

gtasm.prds <- ggplot() + 
  geom_tile(data = tble.tasm %>% 
              filter(period %in% c('Línea base', 'Futuro')), 
            aes(
              x = x, 
              y = y, 
              fill = value)
  ) + 
  facet_wrap(~period) +
  scale_fill_gradientn(colors = brewer.pal(n = 9, name = 'YlOrRd'), 
                       labels = seq(-3, 32, 6), 
                       breaks = seq(-3, 32, 6),
                       name = 'Temperatura (°C)') +
  geom_sf(data = wrld, fill = NA, col = 'grey30') +
  geom_sf(data = st_as_sf(col), fill = NA, col = 'grey40') +
  coord_sf(xlim = ext(col)[1:2], ylim = ext(col)[3:4]) +
  labs(x = 'Lon', y = 'Lat', fill = 'Temperatura (°C)', caption = 'Datos de línea base: Worldclim v2.1 (Hijmans et al., 2014)\nDatos de futuro: CMIP6 - SSP 370 GCM: Ec-Earth3-Veg-Lr') +
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

ggsave(plot = gtasm.prds, filename = './png/maps/v2/tasm_prds.png', units = 'in', width = 11, height = 9, dpi = 300)

# Difference
tasm.dfrn <- tble.tasm %>% 
  filter(period == 'Diferencia') %>% 
  mutate(period = ifelse(period == 'Diferencia', 'Delta', NA)) 

gtasm.dfrn <- ggplot() + 
  geom_tile(data = tasm.dfrn,
            aes(
              x = x, 
              y = y, 
              fill = value)
  ) + 
  facet_wrap(~period) +
  scale_fill_gradientn(colors = brewer.pal(n = 9, name = 'YlOrRd'), 
                       name = 'Temperatura (°C)') +
  geom_sf(data = wrld, fill = NA, col = 'grey30') +
  geom_sf(data = st_as_sf(col), fill = NA, col = 'grey40') +
  coord_sf(xlim = ext(col)[1:2], ylim = ext(col)[3:4]) +
  labs(x = 'Lon', y = 'Lat', fill = 'Temperatura (°C)', caption = '') +
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

gtasm.all <- ggpubr::ggarrange(gtasm.prds, gtasm.dfrn, ncol = 2, nrow = 1)
ggsave(plot = gtasm.all, filename = './png/maps/v2/tasm_vals-delta_ssp370.png', units = 'in', width = 13, height = 9, dpi = 300)

terra::writeRaster(x = stck, filename = './png/maps/v2/stack.tif')
