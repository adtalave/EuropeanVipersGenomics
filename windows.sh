#This script creates a file listing 100-kbp windows in the format needed to estimate heterozygosity following the script provided within this repository.

genome="path/to/your/reference/genome"

samtools faidx "$genome".fasta
cut -f1,2 "$genome".fai > sizes.genome

#You can keep the ones with a total size higher than 10Mb (or your desired threshold) to only keep chromosomes
awk '$2 >= 10000000' sizes.genome > complete_sizes.genome

#By modifying this line you can filter out the sexual chromosomes to calculate autosomal windows only.
grep -v "NAME_OF_SEX_CHROM" complete_sizes.genome | grep -v "OZ078259.1" > autosome_sizes.genome

#By changing the flag -w you can select a different window length.
bedtools makewindows -g autosome_sizes.genome -w 100000 > windows.txt

#We add +1 to the beginning of every window:
awk '{ $2 = $2 + 1; print }' windows.txt > modified_windows.txt

#We re-format into: "CHR:start-end"
awk ' { print $1 ":" $2 "-" $NF } ' modified_windows.txt > complete_windows.txt 

#remove intermediate files;
rm modified_windows.txt
rm windows.txt
rm sizes.genome

echo "Your list of 100-kbp windows are written in the file complete_windows.txt now."
