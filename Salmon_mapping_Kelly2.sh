for fn in  /mnt/s/AG/AG-Scholz-NGS/Daten/RNASeq_Kelly2_P3302/*_R1*;
	do
	samp=$(echo "`basename ${fn}`" | cut -c16-31);
	R1=$fn;
	R2=$(echo "$R1" | sed 's/R1/R2/');
	echo "Processing Sample: $samp";
	test -f $R1 && echo "--> File: $R1 exists"
	test -f $R2 && echo "--> File: $R2 exists"
	
	salmon quant -i /mnt/s/AG/AG-Scholz-NGS/Daten/Salmon/index/human_ensh38_index -l A \
	-1 $R1 \
	-2 $R2 \
	-p 30 --validateMappings --gcBias -o /mnt/s/AG/AG-Scholz-NGS/Daten/Simon/'RNA-Seq Kelly2 P3302'/quants/${samp}_quant
	
done
