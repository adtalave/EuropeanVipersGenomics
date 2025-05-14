# From raw reads to a VCF file
Here you can follow the pipeline to get a raw VCF from raw Illumina reads.

## Filtering reads with Fastp

We use fastp v0.24.0 to filter out bases below PhredScore 30, trim poly-X and -G tails, overrepresented sequences, and adapters, allowing for Pair-End adapter detection.
Sequencing library names are changed to original sample names on filtered fastq files.

```bash
#!/bin/bash

library=$(sed -n "${SLURM_ARRAY_TASK_ID}p" /list/of/library/IDs.txt)
sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" /list/of/sample/IDs.txt)
project_dir="/path/to/directory/containing/raw/data"
scratch_dir="/path/to/output/directory"

echo "Processing library: $library, now renamed as $sample"
fastp  \
    -i "$project_dir/${library}_1.fq.gz"  \
    -I "$project_dir/${library}_2.fq.gz"  \
    -o "$scratch_dir/${sample}_1_fastp.fastq.gz"  \
    -O "$scratch_dir/${sample}_2_fastp.fastq.gz"  \
	--detect_adapter_for_pe  \
	--correction  \
	--trim_poly_g  \
	--trim_poly_x  \
	--overrepresentation_analysis  \
	--thread 4  \
	--qualified_quality_phred 30  \
	--adapter_sequence AGATCGGAAGAG \
	--adapter_sequence_r2 AGATCGGAAGAG
			
echo "Filtering reads and renaming libraries finished"
```

### Quality Control of filtered reads with fastQC and multiQC


```bash
#!/bin/bash
input_dir="/path/to/input/dir"
output_dir="/path/to/output/dir"
sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" /list/of/sample/IDs.txt)

echo "Performing Quality Control of sample $sample with FastQC"
fastqc -t 2 -o "$output_dir" "$input_dir/${sample}_1_fastp.fastq.gz" "$input_dir/${sample}_2_fastp.fastq.gz"

echo "Finishing fastQC"
```

All fastQC files can be merged for a thorough visualization with multiQC.
```bash
multiqc "$output_dir" -o "$output_dir"
```

## Mapping and duplicate removal
Filtered reads can now be mapped on our reference genome (.fasta or .fna) with bwa mem algorithm.


```bash
#!/bin/bash
samples="/list/of/sample/IDs.txt"
reference="/path/to/reference/genome"
output_dir="/path/to/output/dir"
input_dir="/path/to/input/dir"
i=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $samples)

echo "Mapping sample $i"
bwa mem -t 85 -M $reference "$input_dir/${i}_1_fastp.fastq.gz" \
"$input_dir/${i}_2_fastp.fastq.gz" -R "@RG\tID:$i\tSM:$i" \
| samtools sort -@ 15 -o "$output_dir/$i.bam"

echo "Indexing sample $i"
samtools index "$output_dir/$i.bam"
echo "Sample $i already processed"

echo "Job finished"
```
Mapped reads are filtered by mapping quality (30), and we discard secondary and supplementary alignments (-F 2308) with samtools.
Optical and PCR duplicates are identified and removed with picard MarkDuplicates.

```bash
#!/bin/bash
#SBATCH --job-name=job_name
#SBATCH --output=job_name%A.out
#SBATCH --error=job_name%A.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --time=1:00:00
#SBATCH --array=1-344

module load miniforge
module load samtools
module load java-openjdk/17.0.11+9
module load picard/3.1.1
source activate crea1_env

picard="/path/to/picard.jar"
samples="/list/of/samples.txt"

input_dir="/patch/to/input/directory"
output_dir="/patch/to/output/directory"
output_dir2="/patch/to/second/output/directory"
i=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $samples)

echo "Filtering BAM from sample $i"
samtools view -@ 20 -q 30 -F 2308 -b "$input_dir/${i}.bam" -o "$output_dir/${i}_30.bam"

echo "Deduplicating with picard"
java -jar $picard MarkDuplicates \
	-INPUT "$output_dir/${i}_30.bam" \
	-OUTPUT "$output_dir2/${i}_mkdup.bam" \
	-METRICS_FILE "$output_dir2/${i}.txt" \
	-OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
	-ASSUME_SORT_ORDER coordinate \
	-REMOVE_DUPLICATES True \
	-CREATE_INDEX True \
	-TMP_DIR /tmp

mkdir -p "$output_dir2/coverages"
echo "Calculating final coverage for $i"
samtools depth -@ 20 -a "$output_dir2/${i}_mkdup.bam" | awk '{sum+=$3} END { print "Average = ",sum/NR}'| tee -a "$output_dir2/coverages/${i}_coverage.txt"

echo "Job finished." 
```

For high-coverage samples; you might want to downsample the BAM files for a better comparison with samtools view before calling SNPs.

```bash
samtools view -b -s ratio_to_keep/1 originalbam.bam > downsampled.bam
```

# SNP-calling with GATK
## Haplotype calling

We'll use GATK to call haplotypes for our data. To increase HaplotypeCaller's efficiency, we can parallelize GATK among all samples (if you have enough resources).

```bash
#!/bin/bash
samples="list(of/samples.txt"
gatk="path/to/gatk.jar"
chromosomes="/list/of/chromosomes"
reference="/path/to/reference/genome"
output_dir="/path/to/output/directory"

for sample in $(cat $samples);
do
for i in $(cat $chromosomes)
do
(java -jar $gatk HaplotypeCaller --reference $reference --input "$sample"_mkdup.bam -O "$sample"_"$i".g.vcf -ERC BP_RESOLUTION -L "$i") &
done
wait
done

## CombineGVCFs
Now we combine all files for each chromosome. So after this step, we should have one VCF per chromosome. Some of those files might be heavy.

```bash
#!/bin/bash

reference="/path/to/reference/genome"
chromosomes="/list/of/chromosomes"
out="/path/to/output/directory"
input="path/to/folder/with/haplotype/VCFs"
gatk="path/to/gatk.jar"

for i in $(cat "$chromosomes"); do
    echo "Processing chromosome $i..."

    files=$(ls "$input"/*"$i".g.vcf)

    if [ -z "$files" ]; then
        echo "No files found for chromosome $i"
        continue
    fi

    variant_args=""
    for file in $files; do
        variant_args+=" --variant $file"
    done

    "$gatk" --java-options CombineGVCFs \
    --reference "$reference" \
    $variant_args \
    --convert-to-base-pair-resolution \
    --output combined_"$i".g.vcf.gz
done
```



## Genotype calling

Then, we can use GATK's GenotypeGVCFs tool to call genotypes. 

```bash

#!/bin/bash
reference="path/to/reference/genome"
chromosomes="list/of/chromosomes"
out="path/to/output/directory"
input="path/to/folder/with/haplotype/VCFs"
gatk="path/to/gatk.jar"

for i in $(cat $chromosomes)
do
$gatk --java-options GenotypeGVCFs --reference $reference \
  --variant combined_"$i".g.vcf.gz \
  --include-non-variant-sites \
  --output "$out"/"$i"_dataset.vcf.gz
bcftools index "$out"/"$i"_dataset.vcf.gz
done
```

Now we should have a VCF file per each chromosome and its index. Final step to get a VCF with all positions is concatenating the chromosomes:

```
bcftools concat *_dataset.vcf.gz -o final_file.vcf.gz
```
