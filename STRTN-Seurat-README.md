# Visualization the results using Seurat (R package)

STRT-N library results using Seurat package (v3.0.2) ([Butler et al., 2018](https://www.nature.com/articles/nbt.4096)) can be analysed for data reduction and visualization the expression values and quality check values by STRTN-Seurat.sh.

## Dependencies
- [stringr](https://stringr.tidyverse.org/)
- [dplyr](https://dplyr.tidyverse.org/)
- [ggplot2](https://ggplot2.tidyverse.org/)
- [cowplot](https://www.rdocumentation.org/packages/cowplot/versions/1.1.1)
- [ggbeeswarm](https://github.com/eclarke/ggbeeswarm)
- [forcats](https://forcats.tidyverse.org/)
- [Seurat](https://satijalab.org/seurat/)

## Requirements
- `Example-BarcodesStages.txt` : Sample explanation with developmental stages (in `src` directory).
   #### Example
     |     |     |
     | :-: | :-: |
     | 1 | oocyte | 
     | 2 | oocyte | 
     | 3 | oocyte | 
     | 4 | oocyte | 
     | 5 | oocyte | 
     | 6 | oocyte | 
     | 7 | oocyte | 
     | 8 | oocyte | 
     | 9 | 2cell | 
- Quality check report for all samples, `-QC.txt` (in `out` directory).
 
## Example usage
For general users:
```
./STRTN-Seurat.sh -w /mnt/c/Users/gamyaz/STRTN-Pipeline
```
For CSC users:
```
sbatch -A <project_id> ./STRTN-Seurat-CSC.sh -w /scratch/<project_id>
```

## Parameters
- __Mandatory__

   | Name | Description |
   | :--- | :--- |
   | `-w, --working` | /PATH/to/the working directory. |
   
- __Optional__

   | Name&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|Description|
   | :--- | :--- |
   | `-h, --help`| | Show usage.|
   | `-v, --version`| | Show version.|
     
Note: In the case for 48 mouse samples, outliers and control sample were excluded. To represent %80 of data, `npcs` argument as 33 for PCA and `dims` argument as 4 for UMAP were used.

## Outputs
Outputs are provided in `out` directory.

- __`OUTPUT`-QC-BeeswarmPlots.pdf__ <br>
Visualization quality check values for each developmental stage using BeeswarmPlots.

- __Rplots.pdf__ <br>
Elbow, JackStraw, PCA, UMAP and violin plots.
