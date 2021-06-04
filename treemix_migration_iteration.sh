for m in {1..10}
	do
	for i in {1..3}
	do
		for k in {1..5}
			do
			treemix -i ~/analyses/snp_to_vcfs/treemix_files/afrpapeurd3.tmix.gz -o ~/analyses/snp_to_vcfs/treemix_outputs/migration_selection/migr_select.${i}.${m}.${k} -m ${m} -noss -root Yoruba -k ${k}
		done
	done
done

