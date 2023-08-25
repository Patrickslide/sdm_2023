---
title: "Species Distribution Models"
subtitle: "_Capra ibex_ distribution analysis"
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
Created: 20230426 patrick\
Updated: 20230815 patrick\
:::

# Setup

Set up is a bit awkward, since all the R world is moving away from some
geospatial packages (namely `raster`, `sp` and friends).

Try loading the esential minimum of packages, and note that `biomod2`
since versione 2.0 is capable to use `terra`.

```{r setup, echo=FALSE, message=FALSE}
library(biomod2)   #now biomod2 uses terra
library(terra)
library(ENMTools)  #to perform raster correlation
library(usdm)
#library(rgdal)     
library(raster)    #to read tif files
library(car)
library(sf) 
library(tidyverse)
# define here where maxent.jar is (needed for MAXENT models)
#path_to_maxent.jar <- "/srv/scratch/maxent/maxent.jar"

```

::: comment
Always name the first chunk `setup`, and make it quiet and silent.
You're supposed to just load packages in the `setup`chunk.
:::

# Pepare data

## Data paths

Define where all the relevant data sets are.

To derive ecogeographical variables (EGVs), we will rely on a pre-made
raster dataset covering the whole Alps Convention area, at 100 m spatial
resolution, including variables derived from the `EU-DEM 25` digital
elevation model, CORINE Land Cover dataset, and WORLDCLIM dataset.

```{r set paths}
EGV_Alps <- file.path("/","srv","gis","species_distribution","BASE","EGV Alps","EGV")
ibex.DSN <- file.path("/","srv","gis","species_distribution","BASE","ibex")
ibex.Layer <- "punti_2018"
```

## Load Environmental Data

Load the EGV bundle data. Note that being `r EGV_Alps` on "old"
`raster::stack` object, we have to use `raster::stack`.

Furthermore, some layers aren't necessary for the project, so some will
be dropped them to save storage.

::: comment
Note that in this case we're using just a few rasters from `EGV_Alps`,
this is just for testing, once all will work we'll use more layers.
:::

:::comment (Patrick) For the second attempt, i included all BIO except
8, 9, 18, 19. :::

```{r load data}
EGV <- raster::stack(EGV_Alps)
# keep variables, by name (never use numbers to address array elements...) 
```

To check whether all is working, let's plot one of the layers.

```{r plot}
plot(EGV$Elevation)
```

This encompasses the whole Alpine space, possibly too much, we'll have
to do some work to restrict the study area.

## Load species presence data

Now it load the ibex data. We should have a total of 773 observations,
with 14 fields each.

We also need to turn the points reference system in `epsg:3035`, since
EGV are in that Coordinate Reference System.

```{r load ibex data}
ibex <- st_read(ibex.DSN, ibex.Layer)
ibex <- st_transform(ibex, crs="epsg:3035") # same projection as the EGV raster
```

## Visualization

Now, both the raster `EGV` and the point `ibex` should match when
plotted: usually real men don't plot, they trust their `epsg` codes, but
oh well.

```{r plot all}
plot(EGV$Elevation)
plot(ibex, add=TRUE)
```

# Resize the Area of Interest

As seen in the plot above, we need to resize the area covered by `EGV`:
this is because it is a nonsense to calculate a spatial distribution of
an Alpine species in the Po plain... and also because the smaller the
AOI is, the faster calculations are.

Problem is that we've to devise how big the AOI is from points (i.e.
from species distribution). Usually, we get the *extent* of the presence
points we have and enlarge it by a 10%: the first step (getting the
extent) is trivial, but having the correct size (that "plus 10%") to
make a buffer is anything but trivial.

## Extent matching

Plotting (see above) the `EGV` and `ibex` data together we see that the
area occupied by the points is *smaller*. Best if we crop the raster to
only consider the area containing ibex observations.

Usually this is made taking into account not just the *bounding box*
(i.e. the area exactly occupied by the points), but something slightly
bigger: the standard procedure goes like this: - extract the bounding
box (the *extent* in GIS- peak) from the points - calculate the major
dimension, get 10% of that - use that size to make a buffer to the
bounding box - get the bounding box from the buffered bounding box -
crop the \`EGV raster

```{r crop}
# get bounding box
extent <- st_bbox(ibex)
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

::: comment
Note that the procedure here must be repeated when we'll decide how many
layers to keep from the original `EGV_Alps`, since here we're cropping
just the EGV subset we created before.
:::

## Plot, just to be safe

Although it is time consuming (but not now, since we've heawily reduced
the `EGV` dataset), let's have a plot.

```{r plot cropped}
plot(EGV_cropped$Elevation)
plot(ibex, add=TRUE)
```

So far, so good. We now have presence points (`ibex`) and (a few)
ecogeographical variables (`EGV_cropped`) to use to build a model.

## Factorial Analysis

Of all the raster layers available, some may be highly correlated and
including them would just be redundant. Through factorial analysis I
seek to determine which are the least important, including only for this
step even the ones I'd normally ignore.

Removing unnecessary layers will improve our models and reduce
collinearity.

### Pearson correlation Matrix

Function taken from this link:
<https://rdrr.io/github/danlwarren/ENMTools/src/R/raster.cor.plot.R>

```{r}
raster.cor.plot <- function(env, method = "pearson"){
  n.layers <- length(names(env))

  cor.mat <- as.data.frame(matrix(nrow = n.layers, ncol = n.layers))
  colnames(cor.mat) <- names(env)
  rownames(cor.mat) <- names(env)

  for(i in 1:n.layers){
    for(j in i:n.layers){
      cor.mat[i,j] <- abs(raster.cor(env[[i]], env[[j]], method = method))
      cor.mat[j,i] <- cor.mat[i,j]
    }
  }
  colnames(cor.mat) <- names(env)
  row.names(cor.mat) <- names(env)

  # Plot variables in MDS space of correlation
  d <- as.dist(1 - cor.mat)
  mds.cor <- as.data.frame(cmdscale(d))
  colnames(mds.cor) <- c("X", "Y")

  cor.mds.plot <- ggplot(mds.cor, aes_string(x = "X", y = "Y"))  +
    geom_text(aes(label = rownames(mds.cor))) + theme_bw()

  # Make heatmap
  melted.cor.mat <- reshape2::melt(as.matrix(cor.mat))
  colnames(melted.cor.mat) <- c("Var1", "Var2", "value")

  cor.heatmap <- ggplot(data = melted.cor.mat, aes_string(x="Var1", y="Var2", fill="value")) +
    geom_tile() + scale_fill_viridis_c() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

  output <- list(cor.mds.plot = cor.mds.plot,
                 cor.heatmap = cor.heatmap)

  return(output)

}
```

First, I compute a correlation matrix for the first 15 variables of EGV
and solar radiation level.

```{r}
EGV1 <- subset(EGV_cropped, c("Elevation", "Aspect.flat", "Aspect.N", "Aspect.NE", "Aspect.E", "Aspect.SE", "Aspect.S", "Aspect.SW", "Aspect.W", "Aspect.NW", "Slope.level","Slope.moderate", "Slope.strong", "Slope.steep", "Solar.Radiation"))

raster.cor.plot(EGV1)

```

Then, we do the same for all environmental type variables. This wasn't
successfull due to the high number of NA's in some layers, which are
also the primary cause of halts in the models I tried to produce before.
I excluded all the problematic ones, being left with the one here.

::: Update: The next three sections are unnecessary.:::

```{r}
EGV2 <- subset(EGV_cropped, c("Roughness", "Urban.areas", "Seminatural.vegetation", "Woodlands", "Natural.grasslands", "Moors.and.heathland", "Bare.rocks", "Sparsely.vegetated.areas","Burnt.areas", "Marshes.and.bogs", "Water.bodies"))

raster.cor.plot(EGV2)

```

Lastly, we repeat the process for all WorldClim's variables.

```{r}
EGV3 <- subset(EGV_cropped, c("BIO_01", "BIO_02", "BIO_03", "BIO_04", "BIO_05", "BIO_06", "BIO_07", "BIO_08", "BIO_09", "BIO_10", "BIO_11", "BIO_12","BIO_13","BIO_14", "BIO_15", "BIO_16", "BIO_17", "BIO_18", "BIO_19"))

raster.cor.plot(EGV3)
```

Based on the results of the correlation matrix, BIO_8, BIO_9, BIO_10,
BIO_11, and BIO_16, BIO_17, BIO_18, BIO_19 could be omitted, as BIO_04
to BIO_11 and BIO_12 to BIO_19 have a very high pearson coefficient of
correlation between each other and only the two most important variables
of each group will be kept.

### VIF

We performed VIF (Variance Inflation Factor) among all possible factor,
to test inter-factor collinearity and exclude some of them which could
lead to redundancy and estimation errors.

Source: <https://rdrr.io/rforge/usdm/man/vif.html>

```{r}
#v_1 <- vif(EGV_cropped)
v_2 <- vifcor(EGV_cropped, th=0.9, method='pearson')
v_3 <- vifstep(EGV_cropped, th=10, method='pearson')

```

```{r}
plotVIF <- function(aPack) {
    # cast the RasterStack iunto a dataframe without NA
    d <- na.omit(as.data.frame(aPack))
    ct <- cor(d)                        #Correlation
    dm <- abs(as.dist(ct))
    dm[is.na(dm)] <- 0
    dc <- hclust(1-dm)
    dmplot <- plot(dc, xlab="", ylab="Correlation distance (1-r)", main="EGV correlation", sub="")
    print(dmplot)
    invisible(dmplot)
}

plotVIF(EGV_cropped)
```

Based on the results of both tests I removed some of the layers.

```{r}
EGV_cropped <- subset(EGV_cropped, c("BIO_02", "BIO_03","BIO_06", "BIO_13", "BIO_15", "Roughness", "Woodlands", "Natural.grasslands", "Bare.rocks", "Sparsely.vegetated.areas", "Marshes.and.bogs", "Water.bodies", "Elevation", "Aspect.flat", "Aspect.N", "Aspect.E", "Aspect.S", "Aspect.W", "Slope.level", "Slope.moderate", "Slope.strong", "Slope.steep", "Solar.Radiation"))
```

# Species Distribution modeling with `biomod2`

Now that we have our data it's time to start applying biodiversity
models.

Please refer to `biomod2` vignettes, available at
[<https://cran.r-project.org/web/packages/biomod2/vignettes/>].

### Format Data

First, we need to specify all the important factors, variables and
information to use.

We also format the data using the appropriate functions.

Some sparse configuration things to keep in mind, and transfer in our
`biomod2` configuration: - the machine we're working on
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

```{r prepare data}
#Begin the cluster.
beginCluster()
# pass to biomod2 presence data and EGVs


ibexdata <- BIOMOD_FormatingData(resp.var = rep(1, nrow(st_coordinates(ibex))), # we just repeat "1" for as many coordinates we have
                                 resp.xy = st_coordinates(ibex), # and here are the coordinates of those "1"
                                 expl.var = EGV_cropped, # the rasters we prepared
                                 resp.name = "Capra.ibex",
                                 # now deal with pseudorandom pseudoabsences
                                 PA.nb.rep = 1, # number of pseudoabsence sets that will be created
                                 PA.nb.absences = ceiling(1.25 * nrow(st_coordinates(ibex))), # number of pseudoabsence points to generate
                                 PA.strategy = "random", # how to generate pseudo absences
                                 # other options
                                 filter.raster = TRUE, # we want automatic filtering, no more than a single point per raster cell
                                 na.rm = TRUE # please eliminate NA.
                                )

```

## Define modeling options

As well, we also prepare a `myBiomodOptions` data structure, copying all
defaults from `BIOMOD_ModelingOptions()` and changing/adding just the
parameters we're interested in.

```{r set options}
# set default modeling options, just copy them over
myBiomodOptions <- BIOMOD_ModelingOptions()

# add maxent.jar path
myBiomodOptions <- BIOMOD_ModelingOptions(
  MAXENT = list(path_to_maxent.jar = path_to_maxent.jar, memory_allocated = 4096)
)

# other modeling angine specific options here
myBiomodOptions <- BIOMOD_ModelingOptions(
  GLM = list(type = 'quadratic', interaction.level = 1),
  GBM = list(n.trees = 500),
  GAM = list(algo = 'GAM_mgcv'),
  RF = list(mtry=5)
)

```

## Calculate models

Finally, call `BIOMOD_Modeling` to calculate the models. :::comment Note
that since `BIOMOD_Modeling()` output can be rather long, console output
in this code chunk has been "silenced". :::

```{r run biomod2, message=TRUE, warning=TRUE, include=FALSE}

gc()                           #To free unused memory

ibexmodels <- BIOMOD_Modeling(bm.format = ibexdata, 
                              nb.cpu = floor(0.8 * parallel::detectCores()), # use 80% cores
                              bm.options = myBiomodOptions,
                              # change the following parameters if need be...
                              models = c("RF", "GLM", "GBM", "GAM",
                                         "MAXNET"),
                              CV.nb.rep = 2,
                              metric.eval = c("KAPPA",'TSS','ROC'),
                              CV.perc = 0.8,
                              var.import = 2,
                              modeling.id = "demo_1"
                              )
```

How model engines "behaved" can be checked "calling by name" the model
output, so that it "prints" itself:

```{r check}
ibexmodels
```

## Evaluate models

Now, check how the different models performed. The `get_evaluations()`
function produces a dataframe with any possible performance indicator...

```{r evaluate}
#get models evaluations
ibexevals <- get_evaluations(ibexmodels) 
```

Specific performance indicators can be extracted simply indicating
dimension names (just list them with `dimnames()`), or (easier) have a
look at the `ibexevals` dataframe.

Also, model scores can be plotted...

```{r plotscores}
# scores for calibration data set
bm_PlotEvalMean(bm.out = ibexmodels, dataset = 'calibration')
# and for validation data set as well...
bm_PlotEvalMean(bm.out = ibexmodels, dataset = 'validation')
# or grouped by performance indicator:
bm_PlotEvalBoxplot(bm.out = ibexmodels, group.by = c('algo', 'algo'))
```

### Variable importance

To have a look at relative variable importance, just call
`get_variables_importance()` to have another longish dataframe, like
this:

```{r importance0}
ibex_varImp <- get_variables_importance(ibexmodels)
apply(ibex_varImp, c(1,2), mean)

```

Better plot some importance scores, the function is
`bm_PlotVarImpBoxplot()`:

```{r importance1}
# overall importance
bm_PlotVarImpBoxplot(bm.out = ibexmodels, group.by = c('expl.var', 'algo', 'algo'))
 
# grouped by replicate:
#bm_PlotVarImpBoxplot(bm.out = ibexmodels, group.by = c('expl.var', 'algo', 'run'))
```

### Response curves

Further detail can be had looking at the per-variable response curves:
How much of the response variable can be explained by each explanatory
variable in percentage; the fixed metric chosen is median. We are
considering all models, but we can also limit to some of them by setting
a specific algo, PA or run.

```{r responses, message=TRUE, warning=TRUE, include=FALSE}
eval_strip <-
bm_PlotResponseCurves(
bm.out = ibexmodels,
models.chosen = get_built_models(ibexmodels, run = "allRun"),
show.variables= get_formal_data(ibexmodels,'expl.var.names'),
do.bivariate = FALSE,
fixed.var.metric = 'median',
legend = FALSE)
#display_title = FALSE,
```

#### Model Projections

```{r}
ibexproj <- BIOMOD_Projection(bm.mod = ibexmodels,
proj.name = 'Current',
new.env = stack(EGV_cropped),
models.chosen = get_built_models(ibexmodels, run = "allRun"),
metric.binary = 'all',
metric.filter = 'all',
build.clamping.mask = TRUE)

plot(ibexproj)
```

### Ensemble Models

```{r}
ibex_ens_models <-
BIOMOD_EnsembleModeling(
bm.mod = ibexmodels,
em.by = 'all',
metric.select = 'TSS',
metric.select.thresh = 0.7,
metric.eval = c('TSS','ROC'),
em.algo = c('EMmean',"EMcv","EMwmean", "EMca"),
var.import = 1)
get_evaluations(ibex_ensemble_models)
get_variables_importance(ibex_ensemble_models)
```

#### Ensemble Models Variable Importance

```{r}
bm_PlotVarImpBoxplot(bm.out = ibex_ens_models, group.by = c('expl.var', 'algo', 'algo'))
```

#### Ensemble Models Projections

```{r}
ibexEns_proj <- BIOMOD_EnsembleForecasting(bm.em = ibex_ens_models,
proj.name = 'Current EM',
new.env = stack(EGV_cropped),
models.chosen = 'all',
metric.binary = 'all',
metric.filter = 'all')
ibexEns_proj
plot(ibexEns_proj)
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


bioc_30_126 <- raster::stack(bioc_30_126_path)
bioc_30_126 <- as(terra::crop(bioc_30_126, AOI_extent2),"SpatRaster")
bioc_30_126 <- terra::project(bioc_30_126, epsg, method="bilinear")
bioc_30_126 <- terra::project(bioc_30_126, EGV_cropped)

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
plot(bioc_70_126$BIO_02)
plot(bats$geometry, add=TRUE)

```

It worked, but the names need to be changed.

```{r}
names(bioc_30_126) <- c("BIO_01", "BIO_02", "BIO_03", "BIO_04", "BIO_05", "BIO_06", "BIO_07", "BIO_08", "BIO_09", "BIO_10", "BIO_11", "BIO_12","BIO_13","BIO_14",  "BIO_15", "BIO_16", "BIO_17", "BIO_18", "BIO_19")

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
EGV_126_50[["BIO_15"]] <- bioc_50_126$BIO_15

EGV_126_70[["BIO_02"]] <- bioc_70_126$BIO_02
EGV_126_70[["BIO_03"]] <- bioc_70_126$BIO_03
EGV_126_70[["BIO_06"]] <- bioc_70_126$BIO_06
EGV_126_70[["BIO_13"]] <- bioc_70_126$BIO_13
EGV_126_70[["BIO_15"]] <- bioc_70_126$BIO_15

EGV_370_30[["BIO_02"]] <- bioc_30_370$BIO_02
EGV_370_30[["BIO_03"]] <- bioc_30_370$BIO_03
EGV_370_30[["BIO_06"]] <- bioc_30_370$BIO_06
EGV_370_30[["BIO_13"]] <- bioc_30_370$BIO_13
EGV_370_30[["BIO_15"]] <- bioc_30_370$BIO_15

EGV_370_70[["BIO_02"]] <- bioc_30_370$BIO_02
EGV_370_70[["BIO_03"]] <- bioc_30_370$BIO_03
EGV_370_70[["BIO_06"]] <- bioc_30_370$BIO_06
EGV_370_70[["BIO_13"]] <- bioc_30_370$BIO_13
EGV_370_70[["BIO_15"]] <- bioc_30_370$BIO_15

```

## SSP126-2041-2060

#### Future Projection

```{r}
ibex_126_50_proj <- BIOMOD_Projection(bm.mod = ibexmodels,
                                  proj.name = 'Future126',
                                  new.env = EGV_126_50,
                                  models.chosen = get_built_models(ibexmodels, run = "allRun"),
                                  metric.binary = 'all',
                                  metric.filter = 'all',
                                  build.clamping.mask = TRUE)
ibex_126_50_proj
plot(ibex_126_50_proj)
```

#### Ensemble Model Projection

```{r}

ibex_126_50_ens_proj <- BIOMOD_EnsembleForecasting(bm.em = ibex_ens_models,
                                             proj.name = 'Future126EM',
                                             new.env = EGV_126_50,
                                             models.chosen = 'all',
                                             metric.binary = 'all',
                                             metric.filter = 'all')
ibex_126_50_ens_proj
plot(ibex_126_50_ens_proj)
```

## SSP370-2021-2040

#### Future Projection

```{r}
ibex_370_30_proj <- BIOMOD_Projection(bm.mod = ibexmodels,
                                  proj.name = 'Future30',
                                  new.env = EGV_370_30,
                                  models.chosen = get_built_models(ibexmodels, run = "allRun"),
                                  metric.binary = 'all',
                                  metric.filter = 'all',
                                  build.clamping.mask = TRUE)
ibex_370_30_proj
plot(ibex_370_30_proj)
```

#### Ensemble Model Projection

```{r}
ibex_370_30_ens_proj <- BIOMOD_EnsembleForecasting(bm.em = ibex_ens_models,
                                             proj.name = 'Future30',
                                             new.env = EGV_370_30,
                                             models.chosen = 'all',
                                             metric.binary = 'all',
                                             metric.filter = 'all')
ibex_370_30_ens_proj
plot(ibex_370_30_ens_proj)
```

## SSP370-2061-2080

#### Future Projection

```{r}
ibex_370_70_proj <- BIOMOD_Projection(bm.mod = ibexmodels,
                                  proj.name = 'Future30',
                                  new.env = EGV_370_70,
                                  models.chosen = get_built_models(ibexmodels, run = "RUN1"),
                                  metric.binary = 'all',
                                  metric.filter = 'all',
                                  build.clamping.mask = TRUE)
ibex_370_70_proj
plot(ibex_370_70_proj)
```

#### Ensemble Model Projection

```{r}
ibex_370_70_ens_proj <- BIOMOD_EnsembleForecasting(bm.em = ibex_ens_models,
                                             proj.name = 'Future30EM',
                                             new.env = EGV_370_70,
                                             models.chosen = 'all',
                                             metric.binary = 'all',
                                             metric.filter = 'all')
ibex_370_70_ens_proj
plot(ibex_370_70_ens_proj)
```

