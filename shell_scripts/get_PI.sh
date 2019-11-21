cd /storage/home/users/pjt6/newton/final_genome2/lindley_newton_SNPS/


values="100 500 1000 2000 3000 5000 7000 10000 15000"


for v in ${values}
do
	echo "Running PI ${v}"
	cmd="vcftools --vcf Newton_Lindley_SNPS_vs_newton1.0.vcf 
    --window-pi ${v} 
    --out Newton_vs_lindley_PI.${v}.sitepi  " 
	echo ${cmd}
	eval ${cmd}
	wait
	done
	


