mr_sel_path=$(cat config.json  | python -c 'import json,sys;obj=json.load(sys.stdin);print obj["mr_sel_path"]')

## Run a basic simulation

cd $mr_sel_path/Simulations/basic
slim basic_sim.txt > o.txt
open -a Rstudio basic.R

## Run a basic simulation with two phenotypes

cd $mr_sel_path/Simulations/two_pheno
slim two_pheno_sim.txt > o.txt
open -a Rstudio two_pheno.R

## Run a simulation using number of children as a proxy for fitness
## Basic simulation with one phenotype

cd $mr_sel_path/Simulations/children
slim children.txt > o.txt
Rscript createdatfile.R $PWD
open -a Rstudio MRinR_children.R

## Run a children simulation with a large number of individuals

cd $mr_sel_path/Simulations/children/children_large_n
slim children_large_n.txt
Rscript $mr_sel_path/Simulations/children/createdatfile.R $PWD

## Run the children simulation a large number of times 

cd $mr_sel_path/Simulations/children/multisim
for i in {1..1000}
do
	mkdir sim_$i
	cp ../children.txt sim_$i/children.txt
	cd sim_$i
	slim children.txt > o.txt
	Rscript $mr_sel_path/Simulations/children/createdatfile.R $PWD
	cd ..
done

## Run the children simulation with seven genes affecting two phenotypes:
## three for one, three for the other, one for both

cd $mr_sel_path/Simulations/children/two_pheno_7_gene
slim children_2p7g.txt > o.txt
Rscript $mr_sel_path/Simulations/children/createdatfile.R $PWD

## Run the children simulation with seven genes affecting two phenotypes for 100,000 individuals:
## three for one, three for the other, one for both

cd $mr_sel_path/Simulations/children/two_pheno_7_gene_large
slim children_2p7g_large.txt > o.txt
Rscript $mr_sel_path/Simulations/children/createdatfile.R $PWD

## Run the children simulation with two phenotypes caused by one gene

cd $mr_sel_path/Simulations/children/two_pheno_one_gene
slim children_2p1g.txt > o.txt
Rscript $mr_sel_path/Simulations/children/createdatfile.R $PWD

## multi generational sim
## Selection stops after 10 generations
## run sim 100 times 
for i in {1..100}
do
	mkdir $mr_sel_path/Simulations/children/multigen/selStop_sim_$i
	cp children_multigen.txt $mr_sel_path/Simulations/children/multigen/selStop_sim_$i/
	cd $mr_sel_path/Simulations/children/multigen/selStop_sim_$i
	slim children_multigen.txt > o.txt
	for g in {0..19}
	do
		Rscript $mr_sel_path/Simulations/children/createdatfile.R $PWD output_$g.txt output_$((g+1)).txt dat_$g.txt
	done
	cd ..
done

## Get grandchildren information for each of the simulations for each generation
for i in {1..100}
do
	cd $mr_sel_path/Simulations/children/multigen/selStop_sim_$i
	for g in {0..18}
	do
		Rscript $mr_sel_path/Simulations/children/createdatfile_gc.R $PWD dat_$g.txt dat_$((g+1)).txt dat_gc_$g.txt
	done
	cd ..
done



#### IN PROGRESS



## Run the children simulation with one phenotype and 100 genes
cd $mr_sel_path/Simulations/children/one_pheno_100_gene
slim children_1p100g.txt > o.txt
Rscript $mr_sel_path/Simulations/children/createdatfile.R $PWD





