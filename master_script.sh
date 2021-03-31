
## Run a basic simulation

cd /Users/ba20701/Documents/MR-and-selection/Simulations/basic
slim basic_sim.txt > o.txt
open -a Rstudio basic.R

## Run a basic simulation with two phenotypes

cd /Users/ba20701/Documents/MR-and-selection/Simulations/two_pheno
slim two_pheno_sim.txt > o.txt
open -a Rstudio two_pheno.R

## Run a simulation using number of children as a proxy for fitness
## Basic simulation with one phenotype

cd /Users/ba20701/Documents/MR-and-selection/Simulations/children
slim children.txt > o.txt
Rscript createdatfile.R $PWD
open -a Rstudio MRinR_children.R

## Run a children simulation with a large number of individuals

cd /Users/ba20701/Documents/MR-and-selection/Simulations/children/children_large_n
slim children_large_n.txt
Rscript /Users/ba20701/Documents/MR-and-selection/Simulations/children/createdatfile.R $PWD

## Run the children simulation a large number of times 

cd /Users/ba20701/Documents/MR-and-selection/Simulations/children/multisim
for i in {1..1000}
do
	mkdir sim_$i
	cp ../children.txt sim_$i/children.txt
	cd sim_$i
	slim children.txt > o.txt
	Rscript /Users/ba20701/Documents/MR-and-selection/Simulations/children/createdatfile.R $PWD
	cd ..
done

## Run the children simulation with seven genes affecting two phenotypes:
## three for one, three for the other, one for both

cd /Users/ba20701/Documents/MR-and-selection/Simulations/children/two_pheno_7_gene
slim children_2p7g.txt > o.txt
Rscript /Users/ba20701/Documents/MR-and-selection/Simulations/children/createdatfile.R $PWD

## Run the children simulation with seven genes affecting two phenotypes for 100,000 individuals:
## three for one, three for the other, one for both

cd /Users/ba20701/Documents/MR-and-selection/Simulations/children/two_pheno_7_gene_large
slim children_2p7g_large.txt > o.txt
Rscript /Users/ba20701/Documents/MR-and-selection/Simulations/children/createdatfile.R $PWD




## Run the children simulation with two phenotypes caused by one gene

cd /Users/ba20701/Documents/MR-and-selection/Simulations/children/two_pheno_one_gene
slim children_2p1g.txt > o.txt
Rscript /Users/ba20701/Documents/MR-and-selection/Simulations/children/createdatfile.R $PWD

## multi generational sim
for i in {1..20}
do
	cd /Users/ba20701/Documents/MR-and-selection/Simulations/children/multigen/selStop_gen_$i
	slim children_multigen_$i.txt > o.txt
	for g in {0..19}
	do
		Rscript /Users/ba20701/Documents/MR-and-selection/Simulations/children/createdatfile.R $PWD output_$g.txt output_$((g+1)).txt dat_$g.txt
	done
done



#### IN PROGRESS



## Run the children simulation with one phenotype and 100 genes
cd /Users/ba20701/Documents/MR-and-selection/Simulations/children/one_pheno_100_gene
slim children_1p100g.txt > o.txt
Rscript /Users/ba20701/Documents/MR-and-selection/Simulations/children/createdatfile.R $PWD





