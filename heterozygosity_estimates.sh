#!/bin/bash

#change these paths to your own files
vcf_complete="path/to/your/VCF/with/autosomes"
samples="path/to/a/list/of/samples.txt"
output_dir="path/to/output/directory"
windows="path/to/windows.txt" #you can find the script to produce this file in the same repository
output_summary_file="path/to/output/summary_file"

for sample in $(cat $samples);
do
(

echo "Processing sample $sample"

#remove pre-existing output

output_file="${output_dir}/het_${sample}_autosome.txt"
if [ -f "$output_file" ]; then
    rm "$output_file"
fi
touch "$output_file"


echo "Checking whether single-sample vcf already exists"
if [ ! -f "${output_dir}/singleQ0_${sample}.vcf.gz.tbi" ]; then

	echo "Creating single-sample vcf"
	bcftools view -s "$sample" -m1 -M2 -V indels -i 'QUAL>=30 & FMT/DP>=2 & FMT/DP<=300' "$vcf_complete" -Oz --threads 8 -o "${output_dir}/single_${sample}.vcf.gz"
	vcf_tmp="${output_dir}/singleQ0_${sample}.vcf.gz"

	echo "Creating index for single-sample vcf"
	tabix -p vcf "$vcf_tmp" || { echo "Error with tabix when indexing $vcf_tmp"; exit 1; }

else
    echo "Indexed single-sample vcf already exists, proceeding."
    vcf_tmp="${output_dir}/singleQ0_${sample}.vcf.gz"
fi

echo "Starting to calculate heterozygosity by 100-kbp windows in sample $sample"

while read -r region; do
    {
        chrom_short=$(echo "$region" | awk -F':' '{print $1}')
        
        echo "Processing sample $sample in window $region from chromosome $chrom_short"

        #Get number of callable sites
        callable=$(tabix -h "$vcf_tmp" "$region" | grep -v '#' | wc -l)

        if [ "$callable" -gt 0 ]; then
            echo "Window $region is callable"
            value=${region#*:}
            start=$(echo "$value" | cut -d- -f1)
            tabix -h "$vcf_tmp" "$region" | \
            vcfhetcount | tail -n 1 | \
            awk -v l="$callable" -v chrom="$chrom_short" -v start="$start" '{print chrom"\t"start+50000"\t"$1"\t"l"\t"$1/l}' | \
            tee -a "$output_file"
        fi
    }
done < "$windows"

) &
done
wait

echo "Heterozygosity estimates have been calculated in independent files for each sample, for windows of 100 kbp"


#to summarize your results in mean and median genome-wide estimates per sample:

cd ${output_dir}
echo "Producing a summary file with your results..."

echo -e "Sample\tMedian\tMean" > "$output_summary_file"

for file in het_*_autosome.txt; do
    sample=$(echo "$file" | cut -d'_' -f2)
    sample=$(echo "$sample" | sed 's/down$//')

    awk 'NR>1 {print $5}' "$file" | awk '
    {
        values[NR] = $1;
        sum += $1;
    }
    END {
        n = NR;
        if (n == 0) {
            median = "NA";  # Maneja archivos vacíos después del header
            mean = "NA";
        } else {
            asort(values);
            if (n % 2 == 1) {
                median = values[(n+1)/2];
            } else {
                median = (values[n/2] + values[n/2+1]) / 2;
            }
            mean = sum / n;
        }
        print median, mean;
    }' | awk -v s="$sample" '{print s, $1, $2}' OFS='\t' >> "$output_summary_file"
done

"All finished. Your results are writen on ${output_summary_file}"
