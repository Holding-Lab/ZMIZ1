## Quick Start

Below is a quick example that will generate the figures.

```R
# Install Package
library(devtools)
install_github("andrewholding/ZMIZ1")
# You may also need to install packages from BioConductor e.g. Vulcan

# Load package
library(ZMIZ1)

figure1b_qPlexRIME()
figure2a_PLA()
figure3c_ChIPqPCR()
figure3e_E2F2Timecourse()
figure4_cellGrowth()
figure5c_VIPER()
figure5e_ki67()

figureS2() #PLA Controls
figureS9_siZMIZ1_ESR1_ZMIZ1_reads()
figureS11() #GSEA Part 1
figureS12_E2_gene_set_volcano()
figureS13mcf7()  #GSEA Part 2a
figureS13t47d()  #GSEA Part 2b
figureS14_Survival_METABRIC()
figureS14_Survival_TCGA()
figureS15_METABRIC()
figureS15_TCGA()

```
