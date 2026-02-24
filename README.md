# Community-level functional stability in East African herbivores amid major environmental change

This is the repository with data and code for the manuscript "Community-level functional stability in East African herbivores amid major environmental change". 

## Datasets

### Genus occurrences

Genus occurrences are provided in the genus_number.csv file, which includes taxonomic information along with geographic data and stratigraphic units. Ages are given by the member or formation, or as point dates.
The data draws from an extensive specimen level list published by Fortelius et al. (2016) with taxonomic updates as well as from the NOW Database of Fossil Mammals. Modern occurrences were retrieved from IUCN species ditribution maps via Phylacine.

### Genus traits

Genus functional traits are given in the genus_traits.xlsx file, which lists dental functional traits along with body mass (estimates). Traits from literature are indicated and given with a reference in the corresponding sheet.

### Pedogenic carbonate stable isotopes

Stable isotopes data from pedogenic carbonates are given in the file pedogenic_carbonate_isotopes.xlsx. This dataset comprises C and O stable isotopes (VPDB) from literature for the relevant regions and ages.
Sources are given in the reference sheet.

### Dental enamel stabel isotopes from large herbivorous mammals

Stable isotopes data from large mammal herbivore dental enamel are given in the file large_mammal_isotopes.xlsx. This dataset comprises C (VPDB) from literature for the relevant  regions,  ages, and taxa, reflecting dietary ecology of the respective communities.
Sources are given in the reference sheet.

## Code

The R script (R version 4.3.1) used in our analyses is proved in the file Traits-EnvironProxies.R. This script provides the functions for the temporal trait analysis, for pairwise correlations of individual functional traits with environmental and dietary proxies, and peforms the canonical correlation analysis for multidimensional trait and proxy combinations.
