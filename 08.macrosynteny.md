# Macrosynteny

From the reference genome .fasta and .gff files form each species to compare, we generate a species.cds file with the software [AGAT](https://agat.readthedocs.io/en/latest/tools/agat_sp_extract_sequences.html).
Then we generate a species.bed with jcvi toolkit.

```bash
agat_sp_extract_sequences.pl --clean_final_stop --clean_internal_stop --fasta genome.fasta --gff genome.gff --output species.cds --type cds

#for TSEBRA-annotated genomes, you might need to reformat:
#sed 's/	transcript	/	mRNA	/g'  species.gff > species_new.gff

python3 -m jcvi.formats.gff bed --type=mRNA --key=ID genome_checked.gff -o species.bed
```
Once we have species.cds and species.bed for each species (e.g. species1.cds, species1.bed, species2.cds, and species2.bed), we looked for orthologs.

#Con los species.cds y species.bed podemos buscar ortólogos. Se deberían llamar igual las parejas de archivos (i.e. species1.cds & species1.bed)
```bash
python3 -m jcvi.compara.catalog ortholog --no_strip_names --notex 'species1' 'species2'
python3 -m jcvi.compara.synteny screen --minspan=30 --simple species1.species2.anchors species1.species2.anchors.new #we simplify anchor blocks
```

To plot the karyotypes, we need a list of the chromosomes of each species to be included, i.e., the seqids file:

```
CHR1,CHR2,CHR3,CHR4,CHRZ
CHR1,CHR2,CHR3,CHR4, CHR5,CHR6,CHRZ
```
And a layout file like this one:



```

# y, xstart, xend, rotation, color, label, va,  bed
 .8,     .1,    .8,       0, #000000, species1, top, species1.bed
 .6,     .1,    .8,       0, #000000, species2, top, species2.bed

# edges
e, 0, 1, species1.species2.anchors.simple
```

```
python -m jcvi.graphics.karyotype --chrstyle=roundrect --basepair --outfile=macrosynteny_unstranded.png --figsize=10x8 --dpi=500 --format=png seqids layout
```
This command produces the karyotype plot, but it's likely that you want certain chromosomes to be reversed, which can be done with samtools faidx.
This example calculates the reverse sequenced for CHR1 and CHR2:

```
samtools faidx -o species.revComp --reverse-complement --mark-strand no species.fasta CHR1 CHR2
```
We remove these from the reference fasta,
```
echo -e "CHR1\nCHR2" > species_revchrom12.ids
seqkit grep --invert-match -f species_revchrom.ids -o species.normal species.fasta
```
And then, combined reversed and original chromosomes into a new reference.
```
cat species.normal species.revComp | seqkit sort -N -o species_stranded.fa
```
Last, you can reannotate this new .fasta (species_stranded.fa) using the original .fasta and .gff files of each species with GeMoMa and re-do macrosynteny with MCScan.
´´´
java -Xmx<memory>G -jar /path/to/GeMoMa/GeMoMa-1.9.jar CLI GeMoMaPipeline threads=<threads> outdir=output/folder GeMoMa.Score=ReAlign AnnotationFinalizer.r=NO o=true t=species_stranded.fa i=species a=genome.gff g=genome.fasta
´´´
