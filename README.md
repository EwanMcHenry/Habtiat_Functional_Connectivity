



# Script explainer

Dependant on `U.utilities` package

``` r
devtools::install_github("EwanMcHenry/U.utilities")
```

`~/00 Master - Connectivity metric computation01.R` - Runs majority of relevant scripts, covering:
- Configuration of model and data directories
- Computation of metric
- Presentation of maps etc.

# Data spec

Most data is loaded by `~/subparts of calculation/sub02- load data01.R`

LCM data^[[UKCEH Land Cover Maps \| UK Centre for Ecology & Hydrology](https://www.ceh.ac.uk/data/ukceh-land-cover-maps)] - main basis of model, 
- Defines focal habitat extent, extent of negative impact from neighbouring land use and cost of dispersal.
- Directories in `lcm.directs` list object, from `U.utilities` package, which provides a list of years and corresponding lcm maps for GB and Northern Ireland. 

AWI^[[Ancient Woodland Inventory (Scotland) - data.gov.uk](https://www.data.gov.uk/dataset/c2f57ed9-5601-4864-af5f-a6e73e977f54/ancient-woodland-inventory-scotland1); [Ancient Woodland (England)](https://naturalengland-defra.opendata.arcgis.com/datasets/ancient-woodland-england); [Ancient Woodland Inventory 2021 \| DataMapWales](https://datamap.gov.wales/layers/inspire-nrw:NRW_ANCIENT_WOODLAND_INVENTORY_2021); [Back on the Map - Ancient Tree Inventory](https://ati.woodlandtrust.org.uk/back-on-the-map) ] 
- Used to define higher quality woodland 
- Directories in `~/subparts of calculation/sub01 - configuration.R`

NWSS^[[Native Woodland Survey Data Explorer \| Scottish Forestry](https://www.forestry.gov.scot/native-woodland-survey-scotland-data-explorer)] 
- Used to define native conifer woodland

`dispers.costs` - permeability and negative edge impact extent data
- Located in `~\\Data\\hab costs and edge effects Eycott 2011.csv`
- Largely from [[Eycott2011FillingEvidenceGapsExpertOpinionUse|Eycott et al. 2011]]^[Eycott, A. E., Marzano, M., & Watts, K. (2011). Filling evidence gaps with expert opinion: The use of Delphi analysis in least-cost modelling of functional connectivity. _Landscape and Urban Planning_, _103_(3–4), 400–409. [https://doi.org/10.1016/j.landurbplan.2011.08.014](https://doi.org/10.1016/j.landurbplan.2011.08.014)]
- Small amount of curation within `~/subparts of calculation/sub01 - configuration.R`
# outputs

## Connectivity metric object

Includes:
- all.hexgrids[[1]]$hexgrid
- combo.hexgrid

 "grid_id" - ID of countries hexgrid (over square crop of countries object)
 "hex.grid0" - geometry of hexes
 "hex.ha" - area of hexes in hectares
 "lcm.ncells" - N of cells from CEH LCM in hex (fractions possible)
 "bl.ncells" - N of broadleaf cells from CEH LCM in hex (fractions possible)
 "landnotcoastal.ncells" - N of land cells not coastal from CEH LCM in hex (fractions possible)
 "n.clumps" - N habitat patches in hex
 "tot.patch.ha" - ha of habitat in hex
 "tot.aw.patch.ha" - ha of habitat within ancient woodland inventory area in hex
 "tot.edge.patch.ha" - ha of habitat within edge habitat in hex  
 "tot.awedge.patch.ha" - ha of habitat within ancient woodland and path-edge in hex

* standardised - adjusted for proportion of hex with "land"
  "hex.standardised.leastcost.eca" - primary connectivity metric - hex equivalent connected area (standardised for land area), using least cost distance between patches, and quality weighting ancient woodland and edge habitat, standardised to 1 in landscape

not standardised:
  "hex.standardised.euclid.eca" - hex equivalent connected area (standardised for land area), using euclidean distance between patches, and quality weighting ancient woodland and edge habitat, standardised to 1 in landscape
  "hex.leastcost.eca" - total hex equivalent connected area, using least cost distance between patches, and quality weighting ancient woodland and edge habitat
  "hex.euclid.eca" - hex equivalent connected area, using euclidean distance between patches, and quality weighting ancient woodland and edge habitat
  
* scaled - average cell in landscape has cost of 1
  "hex.standardised.scaled.leastcost.eca" - hex equivalent connected area (standardised for land area), using least cost distance between patches, cell cost scaled so that average cost in landscape = 1, and quality weighting ancient woodland and edge habitat
  "hex.scaled.leastcost.eca" - hex equivalent connected area, using least cost distance between patches with cell cost scaled so that average cost in landscape = 1, and quality weighting ancient woodland and edge habitat, 
  "mean.scaled.ecolog.cost.not.sea" - mean scaled ecological cost within hex
  "median.scaled.ecolog.cost.not.sea"   - median scaled ecological cost within hex
 
 
 
 
 
 