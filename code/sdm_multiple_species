---
title: "Species Distribution Models"
subtitle: "Bats distribution analysis"
output:
  html_document: 
    toc: yes
    toc_float: yes
    highlight: zenburn
    theme: cosmo
    css: style.css
    keep_md: yes
    df_print: kable
editor_options:
  chunk_output_type: console
  markdown:
    wrap: 72
---

::: metadata
Version 0.0\
Created: 20230805 patrick\
Updated: 20230915 patrick\
:::

# Setup

```{r setup, echo=FALSE, message=FALSE}
library(biomod2)   #now biomod2 uses terra
library(terra)
#library(ENMTools)  #to perform raster correlation
library(usdm)
#library(rgdal)     
library(raster)    #to read tif files
library(car)
library(sf) 

# define here where maxent.jar is (needed for MAXENT models)
#path_to_maxent.jar <- "/srv/scratch/maxent/maxent.jar"

```
With this script I'll demonstrate how to work on subset of data, for example separating a data set based on the species it contains and how to produce
analysis on them, both for present and future conditions.

# Prepare data

## Environmental Data

```{r}
EGV_Alps <- file.path("/","srv","gis","species_distribution","BASE","EGV Alps","EGV")

# in this case we have point data already preprocessed as a RDS file 
bats.data <- file.path("/","srv","gis","species_distribution","data", "PNGP_bats_allData.RData")

EGV <- raster::stack(EGV_Alps)
```

VIF has already been carried on the previous script on Capra Ibex Data.
It lead to the removal of several factors.

```{r}
EGV <- subset(EGV, c("BIO_02", "BIO_03", "BIO_05", "BIO_06", "BIO_07", "BIO_13", "BIO_14", "Roughness", "Urban.areas", "Seminatural.vegetation", "Woodlands", "Natural.grasslands", "Bare.rocks", "Water.bodies", "Elevation", "Sparsely.vegetated.areas", "Aspect.N", "Aspect.E", "Aspect.S", "Aspect.W", "Slope.level","Slope.moderate", "Slope.strong", "Slope.steep", "Solar.Radiation"))
```

## Process Species presence data

```{r}
bats <- readRDS(bats.data) # a dataframe with lat-lon points (epsg:4326)
# remove missing coordinates
bats <- bats[!is.na(bats$longitude) & !is.na(bats$latitude),] 
bats <- st_as_sf(x=bats, coords=c('longitude', 'latitude'), crs='epsg:4326')
# reproject into 3035, as EGV
bats <- st_transform(bats, "epsg:3035") # same projection as the EGV raster
```

```{r}
# get bounding box
extent <- st_bbox(bats)
# calculate width and height
width <- extent$xmax - extent$xmin
height <- extent$ymax - extent$ymin
# which one is greater?
bufsize <- 0.1 * ifelse(width > height, width, height)
# round to nearest order of magnitude
bufsize <- round(bufsize , -floor(log10(bufsize)) )
# size up a new, wider, bounding box
bufmatrix <- matrix(c(-bufsize, -bufsize, bufsize, bufsize), ncol=2) 
names(bufmatrix) <- names(extent)
AOI_extent <- extent + bufmatrix
# crop, note we're forcing to have a SpatRaster (i.e. "terra" object)
EGV_cropped <- as(terra::crop(EGV, AOI_extent), "SpatRaster")
```

### Species Distribution modeling with `biomod2`

We have multiple species here, and different sampling techniques: just
have a glimpse with a crosstabulation

```{r}
xtabs(~species+method, bats)

```

#### Rearrange Data

First, better rename species using standard binomials, note that some
"species" are actually *groups* of species that cannot be further
identified using bioacoustics (i.e. "statico" and "transetto" methods).

Create a GENSPE -\> latin name lookup table and applies it

```{r}
genspe2name <- function(aGENSPE) {
  # lookup table here
  lookup <- data.frame(genspe=c("BARBAR", "EPTNIL", "EPTSER", "EPTSER/VESMUR", "HYPSAV", "HYPSAV/PIPSP", "MYOCRY", "MYODAU", "MYOMYO/MYOBLY", "MYOMYS", "MYOSP", "NYCLEI", "NYCLEI/EPTSER", "NYCNOC", "NYCSP", "PIPKUH", "PIPKUH/PIPNAT", "PIPPIP", "PIPPYG", "PLEAUR", "PLEMAC", "PLESP", "TADTEN", "TADTEN/NYCLAS", "VESMUR", "VESMUR/NYCLEI", "VESMUR/NYCSP"), name=c("Barbastella barbastellus", "Eptesicus nilssonii", "Eptesicus serotinus", "E. serotinus vel V. murinus", "Hypsugo savii", "H. savii vel P. pipistrellus", "Myotis crypticus", "Myotis daubentonii", "M. myotis vel blythii", "Myotis mystacinus", "Myotis sp.", "Nycatlus leisleri", "N. leisleri vel E. serotinus", "Nyctalus noctula", "Nyctalus sp.", "Pipistrelus kuhlii", "P. kuhlii vel nathusii", "Pipistrellus pipistrellus", "Pipistrellus pygmaeus", "Plecotus auritus", "Plecotus macrobullaris", "Plecotus sp.", "Tadarida teniotis", "T. Teniotis vel N. lasiopterus", "Vespertilio murinus", "V. murinus vel N. leisleri", "V. murinus vel Nyctalus"))
  return(lookup[match(aGENSPE, lookup$genspe),'name'])
}

# add a nre species.name column...
bats$species.name <- genspe2name(bats$species)

unique(bats$species.name)                 #28 different species.


# redo the crosstabulation (more readable)
xtabs(~species.name+method, bats)
```

Some thoughts on our distribution data modeling presence at small scales (here we have a ~50 x 40 km study area!) is rather difficult, and modeling highly elusive species as bats, further complicates things.

Moreover: the three methods used here do not have similar yields for all the species (although one can't say that they're selective by species). A realistic approach could be like this:
   - exclude species having too small a sample size 
   - create models based on the most affordable data (i.e. captures only)
   - create a further batch of models using the other two methods, that have a greater yield but less precision
   

Here the facilities coming from biomod2 come handy, since we have to repeat many times the modeling step, because we have different species, but we want to set up and configure once and for all the modeling engine, since we have on the other hand a single study area and a single set of EGVs.

First, we consider only "meaningful" species, that us we exclude "species groups" that exist just because the bioacoustics methods (automatic classification of ultrasounds) allow that level of discrimination, and we also discard species with too few (say, less than 10) presence records.
 
# Species to keep:
#   cattura statico transetto what to do
# Barbastella barbastellus             0      25         6 this species is almost impossible to
#                                                          catch, make a comparison between static
#                                                          and cinematic ultrasound sampling
# E. serotinus vel V. murinus          0      42         0 exclude from analysis
# Eptesicus nilssonii                  2     108         2 just build a model with static data
# Eptesicus serotinus                  4     294        37 build all three models
# H. savii vel P. pipistrellus         0       1         0 exclude from analysis
# Hypsugo savii                       65    1207        80 build all three models
# M. myotis vel blythii                0       6         0 exclude from analysis
# Myotis crypticus                     1       0         0 exclude from analysis
# Myotis daubentonii                   2       0         0 exclude from analysis
# Myotis mystacinus                   13       0         0 just build a model with capture data
# Myotis sp.                           0    1260        39 exclude from analysis
# N. leisleri vel E. serotinus         0       1         0 exclude from analysis
# Nycatlus leisleri                    0       5         0 exclude from analysis
# Nyctalus noctula                     0       1         0 exclude from analysis
# Nyctalus sp.                         0       0         2 exclude from analysis
# P. kuhlii vel nathusii               0     811        76 exclude from analysis
# Pipistrellus pipistrellus           13    8336       466 build all three models
# Pipistrellus pygmaeus                0     555         0 just build a model with static data
# Pipistrelus kuhlii                   2       0         1 exclude from analysis
# Plecotus auritus                     3       1         0 exclude from analysis
# Plecotus macrobullaris               2       0         0 exclude from analysis
# Plecotus sp.                         0     182         5 exclude from analysis
# T. Teniotis vel N. lasiopterus       0       3         0 exclude from analysis
# Tadarida teniotis                    0      21         0 just build a model with static data
# V. murinus vel N. leisleri           0       3         0 exclude from analysis 
# V. murinus vel Nyctalus              0     103        45 exclude from analysis
# Vespertilio murinus                  1       1         0 exclude from analysis
# 
# Final summary of species to model:
#   
# Barbastella barbastellus             statico, transetto
# Eptesicus nilssonii                  statico
# Eptesicus serotinus                  statico, transetto, cattura
# Hypsugo savii                        statico, transetto, cattura
# Myotis mystacinus                    cattura
# Pipistrellus pipistrellus            statico, transetto, cattura
# Pipistrellus pygmaeus                statico
# Tadarida teniotis                    statico

```{r}
species.to.keep <- c("BARBAR", "EPTNIL", "EPTSER", "HYPSAV", "MYOSP", "PIPPIP", "PIPPYG", "PIPKUH/PIPNAT","TADTEN")
bats <- bats %>% dplyr::filter(species %in% species.to.keep)

# redo the crosstabulation (once again)
xtabs(~species.name+method, bats)


unique(bats$species.name)            #We are left with only 8 species.
```

Now we need to extract and prepare the 15 data sets...

```{r}
batsBB_s <- bats[bats$species.name =="Barbastella barbastellus" & bats$method=="statico",]
batsBB_t <- bats[bats$species.name =="Barbastella barbastellus" & bats$method=="transetto",]
batsEN_s <- bats[bats$species.name =="Eptesicus nilssonii" & bats$method=="statico",]
batsES_s <- bats[bats$species.name =="Eptesicus serotinus" & bats$method=="statico",]  
batsES_t <- bats[bats$species.name =="Eptesicus serotinus" & bats$method=="transetto",]  
batsES_c <- bats[bats$species.name =="Eptesicus serotinus" & bats$method=="cattura",]  
batsHS_s <- bats[bats$species.name =="Hypsugo savii" & bats$method=="statico",]  
batsHS_t <- bats[bats$species.name =="Hypsugo savii" & bats$method=="transetto",]  
batsHS_c <- bats[bats$species.name =="Hypsugo savii" & bats$method=="cattura",]  
batsMS_c <- bats[bats$species.name =="Myotis sp." & bats$method=="statico",]  
batsPPi_s <- bats[bats$species.name =="Pipistrellus pipistrellus" & bats$method=="statico",]  
batsPPi_t <- bats[bats$species.name =="Pipistrellus pipistrellus" & bats$method=="transetto",]
batsPPi_c <- bats[bats$species.name =="Pipistrellus pipistrellus" & bats$method=="cattura",]  
batsPPy_s <- bats[bats$species.name =="Pipistrellus pygmaeus" & bats$method=="statico",]
batsTT_s <- bats[bats$species.name =="Tadarida teniotis" & bats$method=="statico",] 

batsPPi_2 <- batsPPi[seq(1, nrow(batsPPi), 2), ]
```

Could we just make a dataset for each specie, and then one for each
method? Observing whether some species or some methods are more adequate
for building sdms, while mantaining larger data samples.

I produced all of them for now, we can see later which one should we
use.

```{r}
batsBB <- bats[bats$species.name =="Barbastella barbastellus",]
batsEN <- bats[bats$species.name =="Eptesicus nilssonii",]
batsES <- bats[bats$species.name =="Eptesicus serotinus",]  
batsHS <- bats[bats$species.name =="Hypsugo savii",]  
batsMS <- bats[bats$species.name =="Myotis sp.",] 
batsPPi <- bats[bats$species.name =="Pipistrellus pipistrellus",]  
batsPPy <- bats[bats$species.name =="Pipistrellus pygmaeus",]  
batsTT <- bats[bats$species.name =="Tadarida teniotis",]  

bats_s <- bats[bats$method=="statico",] 
bats_t <- bats[bats$method=="transetto",] 
bats_c <- bats[bats$method=="cattura",]
```

Now that we have our data it is time to start applying biodiversity
models.

```{r}
extent <- st_bbox(batsMS)
# calculate width and height
width <- extent$xmax - extent$xmin
height <- extent$ymax - extent$ymin
# which one is greater?
bufsize <- 0.1 * ifelse(width > height, width, height)
# round to nearest order of magnitude
bufsize <- round(bufsize , -floor(log10(bufsize)) )
# size up a new, wider, bounding box
bufmatrix <- matrix(c(-bufsize, -bufsize, bufsize, bufsize), ncol=2) 
names(bufmatrix) <- names(extent)
AOI_extent <- extent + bufmatrix
# crop, note we're forcing to have a SpatRaster (i.e. "terra" object)
EGV_cropped <- as(terra::crop(EGV, AOI_extent), "SpatRaster")
```


# Biodiversity Models

## Format Data

First, we need to specify all the important factors, variables and
information to use.

We also format the data using the appropriate functions.

Some sparse configuration things to keep in mind, and transfer in our
`biomod2` configuration: - the machine we are working on
(`tayra.uagra.net`) has `r parallel::detectCores()` CPUs, so we use a
fair amount of them - having just presence points, we'll make `biomod2`
randomly generate 10% more "pseudo-absence" points

::: comment
The `biomod2` package offers some handy methods to run several models
*unattended*.

This means that one have to **always** rely on `biomod2` toolchain, even
for trivial things such as preparing presence points and EGVs.

This is necessary so that `biomod2` can *automaGically* arrange things
in case some modeling technique has some particular data formatting
needs: as an user, you want having just one dataset and let R fiddle
with variations.

Also note that is here in preparing data for `biomod2` that things can
go south and, when called into action, `biomod2` fails and spits out
error messages. You've been warned.
:::

For this example I'll work on "Myotis sp.", a specific species of bats.
```{r}
beginCluster()

batsdata <- BIOMOD_FormatingData(resp.var = rep(1, nrow(st_coordinates(batsMS))), # we just repeat "1" for as many coordinates we have
                                 resp.xy = st_coordinates(batsMS), # and here are the coordinates of those "1"
                                 expl.var = EGV_cropped, # the rasters we prepared
                                 resp.name = "Bats",
                                 # now deal with pseudorandom pseudoabsences
                                 PA.nb.rep = 1, # number of pseudoabsence sets that will be created
                                 PA.nb.absences = ceiling(1.1 * nrow(st_coordinates(batsMS))), # number of pseudoabsence points to generate
                                 PA.strategy = "random", # how to generate pseudo absences
                                 filter.raster = TRUE, # we want automatic filtering, no more than a single point per raster cell
                                 na.rm = TRUE)         # please eliminate NA.)
```

Define modeling options

```{r set options}
# other modeling angine specific options here
myBiomodOptions <- BIOMOD_ModelingOptions(
  MAXENT = list(path_to_maxent.jar = path_to_maxent.jar, memory_allocated = 4096)
)

myBiomodOptions <- BIOMOD_ModelingOptions(
  GLM = list(type = 'quadratic', interaction.level = 1),
  GBM = list(n.trees = 500),
  GAM = list(algo = 'GAM_mgcv'),
  RF = list(mtry=5)
)

```

## Calculate Models

```{r run biomod2, message=TRUE, warning=TRUE, include=FALSE}
gc()
batsmodels <- BIOMOD_Modeling(bm.format = batsdata, 
                              nb.cpu = floor(0.8 * parallel::detectCores()), # use 80%  
                              bm.options = myBiomodOptions,
                              CV.strategy = "random",
                              models = c("RF", "GLM", "GBM", "MAXNET"),
                              CV.nb.rep = 3,
                              metric.eval = c("KAPPA",'TSS','ROC'),
                              CV.perc = 0.75,
                              var.import = 2,
                              modeling.id = "demo_bats_PPi")

batsmodels
```

### Evaluate models

Now, check how the different models performed. The `get_evaluations()`
function produces a dataframe with any possible performance indicator...

```{r evaluate}
#get models evaluations
batsevals <- get_evaluations(batsmodels) 
```

Specific performance indicators can be "extracted" simply indicating
dimension names, just list them with `dimnames()`, or (easier) have a
look at the `batsevals` dataframe. Also, model scores can be plotted...

```{r}
bm_PlotEvalMean(bm.out = batsmodels, dataset = 'calibration')

#for validation data set as well...
bm_PlotEvalMean(bm.out = batsmodels, dataset = 'validation')

#or grouped by performance indicator
bm_PlotEvalBoxplot(bm.out = batsmodels, dataset = 'validation', group.by = c('algo','algo'))
```

### Variable importance

To have a look at relative variable importance, just call
`get_variables_importance()` to have another longish dataframe.

```{r}
get_variables_importance(batsmodels)

bm_PlotVarImpBoxplot(bm.out = batsmodels, group.by = c('expl.var', 'algo', 'algo'))

```

### Response curves

Further detail can be obtained looking at the per-variable response
curves, which represent How much of the response variable can be
explained by each explanatory variable in percentage; the fixed metric
chosen is median. We are considering all models, but we can also limit
to some of them by setting a specific algo, PA or run.

```{r responses, message=TRUE, warning=TRUE, include=FALSE}
mods <- get_built_models(batsmodels, run = 'RUN1')

eval_strip <-
bm_PlotResponseCurves(
bm.out = batsmodels,
models.chosen = get_built_models(batsmodels, run="RUN1"),
show.variables= get_formal_data(batsmodels,'expl.var.names'),
do.bivariate = FALSE,
fixed.var.metric = 'median',
legend = FALSE)
#display_title = FALSE,
```

### Model Projections

```{r}
batsproj <- BIOMOD_Projection(bm.mod = batsmodels,
proj.name = 'Current',
new.env = stack(EGV_cropped),
models.chosen = get_built_models(batsmodels, run = "RUN1"),
metric.binary = 'all',
metric.filter = 'all',
build.clamping.mask = TRUE)

plot(batsproj)
```

## Ensemble Models

```{r}
bats_ens_models <-
BIOMOD_EnsembleModeling(
bm.mod = batsmodels,
em.by = 'all',
metric.select = 'TSS',
metric.select.thresh = 0.7z,
metric.eval = c('TSS','ROC'),
em.algo = c('EMmean',"EMcv","EMwmean", "EMca"),
var.import = 1)
get_evaluations(bats_ens_models)
get_variables_importance(bats_ens_models)
```

### Ensemble Models Variable Importance

```{r}
bm_PlotVarImpBoxplot(bm.out = bats_ens_models, group.by = c('expl.var', 'algo', 'algo'))
```

### Ensemble Model Projection

```{r}
# Project ensemble models (building single projections)
batsEns_proj <- BIOMOD_EnsembleForecasting(bm.em = bats_ens_models,
proj.name = 'Current EM',
new.env = stack(EGV_cropped),
models.chosen = 'all',
metric.binary = 'all',
metric.filter = 'all')
batsEns_proj
plot(batsEns_proj)
```

# Future conditions projections

## Load Data

First, we find the directory containing the simulated scenario data.

We then upload the 2021-2040 and 2041-2060 version of future climate
simulated scenarios for bioclimatic variables.

```{r}

#prec_30_126_path <- file.path("/","srv","archivio","Cartografia","WORLDCLIM_v2","future", "EC-Earth3-Veg", "ssp126", "geodata.ucdavis.edu", "cmip6", "30s","EC-Earth3-Veg","ssp126", "wc2.1_30s_prec_EC-Earth3-Veg_ssp126_2021-2040.tif")

#prec_50_126_path <- file.path("/","srv","archivio","Cartografia","WORLDCLIM_v2","future", "EC-Earth3-Veg", "ssp126", "geodata.ucdavis.edu", "cmip6", "30s","EC-Earth3-Veg","ssp126", "wc2.1_30s_prec_EC-Earth3-Veg_ssp126_2041-2060.tif")


bioc_30_126_path <- file.path("/","srv","archivio","Cartografia","WORLDCLIM_v2","future", "EC-Earth3-Veg", "ssp126", "geodata.ucdavis.edu", "cmip6", "30s","EC-Earth3-Veg","ssp126", "wc2.1_30s_bioc_EC-Earth3-Veg_ssp126_2021-2040.tif")

bioc_50_126_path <- file.path("/","srv","archivio","Cartografia","WORLDCLIM_v2","future", "2041-2060", "wc2.1_30s_bioc_EC-Earth3-Veg_ssp126_2041-2060.tif")

bioc_70_126_path <- file.path("/","srv","archivio","Cartografia","WORLDCLIM_v2","future", "EC-Earth3-Veg", "ssp126", "geodata.ucdavis.edu", "cmip6", "30s","EC-Earth3-Veg","ssp126", "wc2.1_30s_bioc_EC-Earth3-Veg_ssp126_2061-2080.tif")

bioc_30_370_path <- file.path("/","srv","archivio","Cartografia","WORLDCLIM_v2","future", "EC-Earth3-Veg", "ssp370", "geodata.ucdavis.edu", "cmip6", "30s","EC-Earth3-Veg","ssp370", "wc2.1_30s_bioc_EC-Earth3-Veg_ssp370_2021-2040.tif")

bioc_70_370_path <- file.path("/","srv","archivio","Cartografia","WORLDCLIM_v2","future", "EC-Earth3-Veg", "ssp370", "geodata.ucdavis.edu", "cmip6", "30s","EC-Earth3-Veg","ssp370", "wc2.1_30s_bioc_EC-Earth3-Veg_ssp370_2061-2080.tif")

```


The wc2.1_30s_bioc_EC-Earth3-Veg_ssp370_2041-2060 is composed of empty data; we use the next one instead.


### Extent and Reference System Matching


First we change the reference system of the World Maps (EPSG:4326) to
the one used for the models (EPSG:3035), then we crop the area to make
it correspond with our species records bounding box (just like we did
before).

```{r}
extent <- st_bbox(bats)
# calculate width and height
width <- extent$xmax - extent$xmin
height <- extent$ymax - extent$ymin
# which one is greater?
bufsize <- 0.1 * ifelse(width > height, width, height)
# round to nearest order of magnitude
bufsize <- round(bufsize , -floor(log10(bufsize)) )
# size up a new, wider, bounding box
bufmatrix <- matrix(c(-bufsize, -bufsize, bufsize, bufsize), ncol=2) 
names(bufmatrix) <- names(extent)
AOI_extent <- extent + bufmatrix


bats_WG <- bats %>%  st_as_sf(coords = c("lon", "lat")) %>% 
  st_bbox() %>% 
  st_as_sfc()

bats_WG <- st_transform(bats_WG, crs="epsg:4326")

extent2 <- st_bbox(bats_WG)
width <- extent2$xmax - extent2$xmin
height <- extent2$ymax - extent2$ymin
# which one is greater?
bufsize <- 0.1 * ifelse(width > height, width, height)
# round to nearest order of magnitude
bufsize <- round(bufsize , -floor(log10(bufsize)) )
# size up a new, wider, bounding box
bufmatrix2 <- matrix(c(-bufsize, -bufsize, bufsize, bufsize), ncol=2) 
names(bufmatrix2) <- names(extent2)
AOI_extent2 <- extent2 + bufmatrix2
```

I was able to make it work. The only thing left is to change them back
to EPSG:3035 (not achieved yet) and crop them to the same size of our
environmental data (they are very similar, but due to conversion they
have a few decimals of difference).

### Create new Environment Datasets

```{r}
epsg <- "epsg:3035"

#bioc_30_126 <- raster::stack(bioc_30_126_path)
#bioc_30_126 <- as(terra::crop(bioc_30_126, AOI_extent2),"SpatRaster")
#bioc_30_126 <- terra::project(bioc_30_126, epsg, method="bilinear")
#bioc_30_126 <- terra::project(bioc_30_126, EGV_cropped)

bioc_50_126 <- raster::stack(bioc_50_126_path)
bioc_50_126 <- as(terra::crop(bioc_50_126, AOI_extent2),"SpatRaster")
bioc_50_126 <- terra::project(bioc_50_126, epsg, method="bilinear")
bioc_50_126 <- terra::project(bioc_50_126, EGV_cropped)

bioc_70_126 <- raster::stack(bioc_70_126_path)
bioc_70_126 <- as(terra::crop(bioc_70_126, AOI_extent2),"SpatRaster")
bioc_70_126 <- terra::project(bioc_70_126, epsg, method="bilinear")
bioc_70_126 <- terra::project(bioc_70_126, EGV_cropped)

bioc_30_370 <- raster::stack(bioc_30_370_path)
bioc_30_370 <- as(terra::crop(bioc_30_370, AOI_extent2),"SpatRaster")
bioc_30_370 <- terra::project(bioc_30_370, epsg, method="bilinear")
bioc_30_370 <- terra::project(bioc_30_370, EGV_cropped)

bioc_70_370 <- raster::stack(bioc_70_370_path)
bioc_70_370 <- as(terra::crop(bioc_70_370, AOI_extent2),"SpatRaster")
bioc_70_370 <- terra::project(bioc_70_370, epsg, method="bilinear")
bioc_70_370 <- terra::project(bioc_70_370, EGV_cropped)
```

we plot the variables and species together just to make sure they are in
the same reference system.

```{r}
plot(bioc_70_126$wc2.1_30s_bioc_EC.Earth3.Veg_ssp126_2061.2080_2)
plot(bats$geometry, add=TRUE)

```

It worked, but the names need to be changed.

```{r}
#names(bioc_30_126) <- c("BIO_01", "BIO_02", "BIO_03", "BIO_04", "BIO_05", "BIO_06", "BIO_07", "BIO_08", "BIO_09", "BIO_10", "BIO_11", "BIO_12","BIO_13","BIO_14",  "BIO_15", "BIO_16", "BIO_17", "BIO_18", "BIO_19")

names(bioc_30_370) <- c("BIO_01", "BIO_02", "BIO_03", "BIO_04", "BIO_05", "BIO_06", "BIO_07", "BIO_08", "BIO_09", "BIO_10", "BIO_11", "BIO_12","BIO_13","BIO_14",  "BIO_15", "BIO_16", "BIO_17", "BIO_18", "BIO_19")

names(bioc_50_126) <- c("BIO_01", "BIO_02", "BIO_03", "BIO_04", "BIO_05", "BIO_06", "BIO_07", "BIO_08", "BIO_09", "BIO_10", "BIO_11", "BIO_12","BIO_13","BIO_14",  "BIO_15", "BIO_16", "BIO_17", "BIO_18", "BIO_19")  

names(bioc_70_126) <- c("BIO_01", "BIO_02", "BIO_03", "BIO_04", "BIO_05", "BIO_06", "BIO_07", "BIO_08", "BIO_09", "BIO_10", "BIO_11", "BIO_12","BIO_13","BIO_14",  "BIO_15", "BIO_16", "BIO_17", "BIO_18", "BIO_19")

names(bioc_70_370) <- c("BIO_01", "BIO_02", "BIO_03", "BIO_04", "BIO_05", "BIO_06", "BIO_07", "BIO_08", "BIO_09", "BIO_10", "BIO_11", "BIO_12","BIO_13","BIO_14",  "BIO_15", "BIO_16", "BIO_17", "BIO_18", "BIO_19")

```

Now we join the simulated bioclimatic variables with the environmental
variables we already had, and project the models on these new future
conditions to estimate the expected changes in species distribution.

```{r}
EGV_126_50 <- EGV_cropped
EGV_126_70 <- EGV_cropped
EGV_370_30 <- EGV_cropped
EGV_370_70 <- EGV_cropped


#EGV_126_30 <- dropLayer(EGV_126_30, c(1:6))

EGV_126_50[["BIO_02"]] <- bioc_50_126$BIO_02 
EGV_126_50[["BIO_03"]] <- bioc_50_126$BIO_03
EGV_126_50[["BIO_06"]] <- bioc_50_126$BIO_06
EGV_126_50[["BIO_13"]] <- bioc_50_126$BIO_13
EGV_126_50[["BIO_14"]] <- bioc_50_126$BIO_14

EGV_126_70[["BIO_02"]] <- bioc_70_126$BIO_02
EGV_126_70[["BIO_03"]] <- bioc_70_126$BIO_03
EGV_126_70[["BIO_06"]] <- bioc_70_126$BIO_06
EGV_126_70[["BIO_13"]] <- bioc_70_126$BIO_13
EGV_126_70[["BIO_14"]] <- bioc_70_126$BIO_14

EGV_370_30[["BIO_02"]] <- bioc_30_370$BIO_02
EGV_370_30[["BIO_03"]] <- bioc_30_370$BIO_03
EGV_370_30[["BIO_06"]] <- bioc_30_370$BIO_06
EGV_370_30[["BIO_13"]] <- bioc_30_370$BIO_13
EGV_370_30[["BIO_14"]] <- bioc_30_370$BIO_14

EGV_370_70[["BIO_02"]] <- bioc_30_370$BIO_02
EGV_370_70[["BIO_03"]] <- bioc_30_370$BIO_03
EGV_370_70[["BIO_06"]] <- bioc_30_370$BIO_06
EGV_370_70[["BIO_13"]] <- bioc_30_370$BIO_13
EGV_370_70[["BIO_14"]] <- bioc_30_370$BIO_14

```

## SSP126-2041-2060

#### Future Projection

```{r}
bats_126_50_proj <- BIOMOD_Projection(bm.mod = batsmodels,
                                  proj.name = 'Future126',
                                  new.env = EGV_126_50,
                                  models.chosen = get_built_models(batsmodels, run = "RUN1"),
                                  metric.binary = 'all',
                                  metric.filter = 'all',
                                  build.clamping.mask = TRUE)
bats_126_50_proj
plot(bats_126_50_proj)
```

#### Ensemble Model Projection

```{r}
bats_126_50_ens_proj <- BIOMOD_EnsembleForecasting(bm.em = bats_ens_models,
                                             proj.name = 'Future126EM',
                                             new.env = EGV_126_50,
                                             models.chosen = 'all',
                                             metric.binary = 'all',
                                             metric.filter = 'all')
bats_126_50_ens_proj
plot(bats_126_50_ens_proj)
```

## SSP126-2061-2080

#### Future Projection

```{r}
bats_126_70_proj <- BIOMOD_Projection(bm.mod = batsmodels,
                                  proj.name = 'Future126',
                                  new.env = EGV_126_70,
                                  models.chosen = get_built_models(batsmodels, run = "RUN1"),
                                  metric.binary = 'all',
                                  metric.filter = 'all',
                                  build.clamping.mask = TRUE)
bats_126_70_proj
plot(bats_126_70_proj)
```

#### Ensemble Model Projection

```{r}
bats_126_70_ens_proj <- BIOMOD_EnsembleForecasting(bm.em = bats_ens_models,
                                             proj.name = 'Future126EM',
                                             new.env = EGV_126_70,
                                             models.chosen = 'all',
                                             metric.binary = 'all',
                                             metric.filter = 'all')
bats_126_70_ens_proj
plot(bats_126_70_ens_proj)
```

## SSP370-2021-2040

#### Future Projection

```{r}
bats_370_30_proj <- BIOMOD_Projection(bm.mod = batsmodels,
                                  proj.name = 'Future30',
                                  new.env = EGV_370_30,
                                  models.chosen = get_built_models(batsmodels, run = "RUN1"),
                                  metric.binary = 'all',
                                  metric.filter = 'all',
                                  build.clamping.mask = TRUE)
bats_370_30_proj
plot(bats_370_30_proj)
```

#### Ensemble Model Projection

```{r}
bats_370_30_ens_proj <- BIOMOD_EnsembleForecasting(bm.em = bats_ens_models,
                                             proj.name = 'Future30EM',
                                             new.env = EGV_370_30,
                                             models.chosen = 'all',
                                             metric.binary = 'all',
                                             metric.filter = 'all')
bats_370_30_ens_proj
plot(bats_370_30_ens_proj)
```

## SSP370-2061-2080

#### Future Projection

```{r}
bats_370_70_proj <- BIOMOD_Projection(bm.mod = batsmodels,
                                  proj.name = 'Future30',
                                  new.env = EGV_370_70,
                                  models.chosen = get_built_models(batsmodels, run = "RUN1"),
                                  metric.binary = 'all',
                                  metric.filter = 'all',
                                  build.clamping.mask = TRUE)
bats_370_70_proj
plot(bats_370_70_proj)
```

#### Ensemble Model Projection

```{r}
bats_370_70_ens_proj <- BIOMOD_EnsembleForecasting(bm.em = bats_ens_models,
                                             proj.name = 'Future30EM',
                                             new.env = EGV_370_70,
                                             models.chosen = 'all',
                                             metric.binary = 'all',
                                             metric.filter = 'all')
bats_370_70_ens_proj
plot(bats_370_70_ens_proj)
```
