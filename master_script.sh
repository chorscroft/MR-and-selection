
## Run a basic simulation

cd /Users/ba20701/Documents/MR-and-selection/Simulations/basic
slim basic_sim.txt
open -a Rstudio basic.R

## Run a basic simulation with two phenotypes

cd /Users/ba20701/Documents/MR-and-selection/Simulations/two_pheno
slim two_pheno_sim.txt
open -a Rstudio two_pheno.R

## Run a simulation using number of children as a proxy for fitness
## Basic simulation with one phenotype

cd /Users/ba20701/Documents/MR-and-selection/Simulations/children
slim children.txt
Rscript createdatfile.R $PWD --vanilla
open -a Rstudio MRinR_children.R