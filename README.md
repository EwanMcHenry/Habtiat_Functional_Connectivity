# Habtiat_Functional_Connectivity

Master - Connectivity metric... Runs the relevent scripts within Subparts calculation


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
  "hex.standardised.euclid.eca" - hex equivalent connected area (standardised for land area), using euclidean distance between patches, and quality weighting ancient woodland and edge habitat, standardised to 1 in landscape
  "hex.leastcost.eca" - total hex equivalent connected area, using least cost distance between patches, and quality weighting ancient woodland and edge habitat
  "hex.euclid.eca" - hex equivalent connected area, using euclidean distance between patches, and quality weighting ancient woodland and edge habitat
  
* scaled - average in landscape has cost of 1
  "hex.standardised.scaled.leastcost.eca" - hex equivalent connected area (standardised for land area), using least cost distance between patches, cell cost scaled so that average cost in landscape = 1, and quality weighting ancient woodland and edge habitat
  "hex.scaled.leastcost.eca" - hex equivalent connected area, using least cost distance between patches with cell cost scaled so that average cost in landscape = 1, and quality weighting ancient woodland and edge habitat, 
  "mean.scaled.ecolog.cost.not.sea" - mean scaled ecological cost within hex
  "median.scaled.ecolog.cost.not.sea"   - median scaled ecological cost within hex
 
 
 
 
 
 