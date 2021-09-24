mkdir 20210921_BleachedPairs_TagSeq

nano /data/putnamlab/KITT/hputnam/20210921_BleachedPairs_TagSeq/BaseSpaceDownload.sh

#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/KITT/hputnam/20210921_BleachedPairs_TagSeq

module load IlluminaUtils/2.11-GCCcore-9.3.0-Python-3.8.2

bs download project -i 293898618 -o /data/putnamlab/KITT/hputnam/20210921_BleachedPairs_TagSeq

sbatch /data/putnamlab/KITT/hputnam/20210921_BleachedPairs_TagSeq/BaseSpaceDownload.sh

Run ID: 
Project ID: 293898618

# QC raw files

### Raw Sequence QC
```
nano /data/putnamlab/KITT/hputnam/20210921_BleachedPairs_TagSeq/scripts/qc.sh
```

```
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/KITT/hputnam/20210921_BleachedPairs_TagSeq/

module load FastQC/0.11.9-Java-11 
module load MultiQC/1.9-intel-2020a-Python-3.8.2

#make qc output folder
mkdir raw_qc/

#run fastqc on raw data
fastqc Raw_Data/*.fastq.gz -o raw_qc/

#Compile MultiQC report from FastQC files
multiqc ./raw_qc
mv multiqc_report.html raw_qc/raw_qc_multiqc_report.html
mv multiqc_data raw_qc/raw_multiqc_data

echo "Initial QC of Seq data complete." $(date)
```

```
sbatch /data/putnamlab/KITT/hputnam/20210921_BleachedPairs_TagSeq/scripts/qc.sh
```

### Assuming they are all run on the same sequencer, we are searching for the read count with @A00842:
```
zgrep -c "@A00842:" ELS*
ELS10_S154_L002_R1_001.fastq.gz:4029092
ELS11_S155_L002_R1_001.fastq.gz:9424605
ELS12_S156_L002_R1_001.fastq.gz:6464861
ELS13_S157_L002_R1_001.fastq.gz:5316658
ELS14_S158_L002_R1_001.fastq.gz:5180350
ELS15_S159_L002_R1_001.fastq.gz:5085566
ELS16_S160_L002_R1_001.fastq.gz:7602907
ELS17_S161_L002_R1_001.fastq.gz:4660744
ELS18_S162_L002_R1_001.fastq.gz:4132829
ELS19_S163_L002_R1_001.fastq.gz:5747439
ELS1_S145_L002_R1_001.fastq.gz:4704830
ELS20_S164_L002_R1_001.fastq.gz:5683815
ELS21_S165_L002_R1_001.fastq.gz:2855452
ELS22_S166_L002_R1_001.fastq.gz:12176183
ELS23_S167_L002_R1_001.fastq.gz:5687433
ELS24_S168_L002_R1_001.fastq.gz:7192473
ELS25_S169_L002_R1_001.fastq.gz:4637094
ELS26_S170_L002_R1_001.fastq.gz:8664050
ELS27_S171_L002_R1_001.fastq.gz:6515034
ELS28_S172_L002_R1_001.fastq.gz:7296565
ELS29_S173_L002_R1_001.fastq.gz:6410050
ELS2_S146_L002_R1_001.fastq.gz:4300211
ELS30_S174_L002_R1_001.fastq.gz:6836136
ELS31_S175_L002_R1_001.fastq.gz:4751888
ELS32_S176_L002_R1_001.fastq.gz:5036467
ELS33_S177_L002_R1_001.fastq.gz:3866534
ELS34_S178_L002_R1_001.fastq.gz:7716555
ELS35_S179_L002_R1_001.fastq.gz:4331840
ELS36_S180_L002_R1_001.fastq.gz:7551260
ELS37_S181_L002_R1_001.fastq.gz:4672725
ELS38_S182_L002_R1_001.fastq.gz:5195458
ELS39_S183_L002_R1_001.fastq.gz:8614131
ELS3_S147_L002_R1_001.fastq.gz:3484849
ELS40_S184_L002_R1_001.fastq.gz:1529935
ELS4_S148_L002_R1_001.fastq.gz:4689271
ELS5_S149_L002_R1_001.fastq.gz:5368385
ELS6_S150_L002_R1_001.fastq.gz:5481764
ELS7_S151_L002_R1_001.fastq.gz:5482118
ELS8_S152_L002_R1_001.fastq.gz:5348723
ELS9_S153_L002_R1_001.fastq.gz:6583875
```

### Sequence trimming
```
nano /data/putnamlab/KITT/hputnam/20210921_BleachedPairs_TagSeq/scripts/trim.sh
```

```
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/KITT/hputnam/20210921_BleachedPairs_TagSeq/Raw_Data/

# load modules needed
module load fastp/0.19.7-foss-2018b

# Make an array of sequences to trim
array1=($(ls *.fastq.gz))

# fastp loop; trim the Read 1 TruSeq adapter sequence; trim poly x default 10 (to trim polyA)
for i in ${array1[@]}; do
        fastp --in1 ${i} --out1 clean.${i} --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --trim_poly_x 6 -q 30 -y -Y 50
done
echo "Read trimming of adapters complete." $(date)

```

```
sbatch /data/putnamlab/KITT/hputnam/20210921_BleachedPairs_TagSeq/scripts/trim.sh
```

mv Raw_Data/clean.*.gz  Clean_Data/

### Trimmed Sequence QC
```
nano /data/putnamlab/KITT/hputnam/20210921_BleachedPairs_TagSeq/scripts/trimmed_qc.sh
```

```
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/KITT/hputnam/20210921_BleachedPairs_TagSeq/Clean_Data

module load FastQC/0.11.9-Java-11 
module load MultiQC/1.9-intel-2020a-Python-3.8.2

#make qc output folder
mkdir clean_qc/

#run fastqc on data
fastqc clean*.fastq.gz -o clean_qc/

#Compile MultiQC report from FastQC files
multiqc ./clean_qc
mv multiqc_report.html clean_qc/clean_qc_multiqc_report.html
mv multiqc_data clean_qc/clean_multiqc_data

echo "Trimmed QC of Seq data complete." $(date)
```

```
sbatch /data/putnamlab/KITT/hputnam/20210921_BleachedPairs_TagSeq/scripts/trimmed_qc.sh
```



#Align to reference genome

```
mkdir refs
cd refs
wget http://cyanophora.rutgers.edu/montipora/Mcap.genome_assembly.fa.gz
wget http://cyanophora.rutgers.edu/montipora/Mcap.GFFannotation.gff
#unzip reference genome
gunzip refs/Mcap.genome_assembly.fa.gz
```

```
nano /data/putnamlab/KITT/hputnam/20210921_BleachedPairs_TagSeq/scripts/genome_index.sh
```

```
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=200GB
#SBATCH --account=putnamlab
#SBATCH -p putnamlab
#SBATCH -D /data/putnamlab/KITT/hputnam/20210921_BleachedPairs_TagSeq/

# load modules needed
module load HISAT2/2.2.1-foss-2019b

# index the reference genome for Montipora capitata output index to working directory
hisat2-build -f /data/putnamlab/KITT/hputnam/20210921_BleachedPairs_TagSeq/refs/Mcap.genome_assembly.fa ./Mcapitata_ref # called the reference genome (scaffolds)
echo "Referece genome indexed." $(date)

```

```
sbatch /data/putnamlab/KITT/hputnam/20210921_BleachedPairs_TagSeq/scripts/genome_index.sh
```

```
nano /data/putnamlab/KITT/hputnam/20210921_BleachedPairs_TagSeq/scripts/align.sh
```

```
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=200GB
#SBATCH --account=putnamlab
#SBATCH -p putnamlab
#SBATCH -D /data/putnamlab/KITT/hputnam/20210921_BleachedPairs_TagSeq/Clean_Data/

# load modules needed
module load HISAT2/2.2.1-foss-2019b
module load SAMtools/1.9-foss-2018b

# This script exports alignments as bam files
# sorts the bam file because Stringtie takes a sorted file for input (--dta)
# removes the sam file because it is no longer needed
array=($(ls clean*)) # call the clean sequences - make an array to align
for i in ${array[@]}; do
        sample_name=`echo $i| awk -F [.] '{print $2}'`
	hisat2 -p 8 --dta -x ../Mcapitata_ref -U ${i} -S ${sample_name}.sam
        samtools sort -@ 8 -o ${sample_name}.bam ${sample_name}.sam
    		echo "${i} bam-ified!"
        rm ${sample_name}.sam
done
```
```
sbatch /data/putnamlab/KITT/hputnam/20210921_BleachedPairs_TagSeq/scripts/align.sh
```

### Alignment Info
slurm-88028.out

Sample.ID|0 times | 1 times | >1 times | overall
---|---|---|---|---|
ELS10|28.93|58.02 |13.05 | 71.07|
ELS11|24.73|61.68 |13.59 | 75.27|
ELS12|26.94|59.10 |13.96 |73.06 |
ELS13|24.54|61.61 |13.85 |75.46 |
ELS14|27.89|59.39 |12.72 |72.11 |
ELS15|26.08|59.98 |13.94 |73.92 |
ELS16|30.09|56.97 |12.94 |69.91 |
ELS17|24.45|62.04 |13.51 |75.55 |
ELS18|23.04|63.55 |13.41 |76.96 |
ELS19|28.42|57.82 |13.76 |71.58 |
ELS1|33.67|53.82 |12.51 |66.33 |
ELS20|26.01|61.01 |12.99 |73.99 |
ELS21|22.84|62.25 |14.92 |77.16 |
ELS22|29.90|57.69 |12.41 |70.10 |
ELS23|22.30|62.54 |15.16 |77.70 |
ELS24|23.52|62.86 |13.62 |76.48 |
ELS25|23.43|62.08 |14.49 |76.57 |
ELS26|20.50|65.45 |14.05 |79.50 |
ELS27|22.82|63.58 |13.60 |77.18 |
ELS28|27.15|60.19 |12.66 |72.85 |
ELS29|26.27|58.99 |14.74 |73.73 |
ELS2|25.18|58.89 |15.93 |74.82 |
ELS30|26.39|59.06 |14.56 |73.61 |
ELS31|25.48|61.62 |12.89 |74.52 |
ELS32|20.99|64.90 |14.12 |79.01 |
ELS33|25.35|60.19 |14.46 |74.65 |
ELS34|21.42|64.27 |14.31 |78.58 |
ELS35|22.81|62.53 |14.66 |77.19 |
ELS36|23.56|60.98 |15.47 |76.44 |
ELS37|19.26|65.56 |15.18 |80.74 |
ELS38|24.95|61.27 |13.78 |75.05 |
ELS39|33.53|54.33 |12.14 |66.47 |
ELS3|32.69|54.26 |13.05 |67.31 |
ELS40|21.99|62.46 |15.55 |78.01 |
ELS4|24.36|61.22 |14.41 |75.64 |
ELS5|32.33|54.80 |12.87 |67.67 |
ELS6|31.54|54.29 |14.27 |68.55 |
ELS7|23.68|62.39 |13.93 |76.32 |
ELS8|26.03|60.28 |13.69 |73.87 |
ELS9|22.95|62.93 |14.11 |77.05 |



# Quantify expression

```
nano /data/putnamlab/KITT/hputnam/20210921_BleachedPairs_TagSeq/scripts/quant_stringtie.sh
```

```
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=250GB
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/KITT/hputnam/20210921_BleachedPairs_TagSeq/Clean_Data/

#load packages 
module load StringTie/2.1.4-GCC-9.3.0

# loop through files to quantify
sh -c 'for i in "ELS1" "ELS2" "ELS3" "ELS4" "ELS5" "ELS6" "ELS7" "ELS8" "ELS9" "ELS10" "ELS11" "ELS12" "ELS13" "ELS14" "ELS15" "ELS16" "ELS17" "ELS18" "ELS19" "ELS20" "ELS21" "ELS22" "ELS23" "ELS24" "ELS25" "ELS26" "ELS27" "ELS28" "ELS29" "ELS30" "ELS31" "ELS32" "ELS33" "ELS34" "ELS35" "ELS36" "ELS37" "ELS38" "ELS39" "ELS40"
do
stringtie ${i}.bam -p 8 -e -G ../refs/Mcap.GFFannotation.fixed.gff -o ${i}.gtf 
echo "${i} gtf complete!"
done'

echo "gtf generation complete!"
```

```
sbatch /data/putnamlab/KITT/hputnam/20210921_BleachedPairs_TagSeq/scripts/quant_stringtie.sh
```

# Assess the performance of the assembly

*Gffcompare is a tool that can compare, merge, annotate and estimate accuracy of GFF/GTF files when compared with a reference annotation*

Using the StringTie merge mode, merge the assembly-generated GTF files to assess how well the predicted transcripts track to the reference annotation file. This step requires the TXT file,  (mergelist.txt). This file lists all of the file names to be merged. *Make sure ```gtf_mergelist.txt``` is in the StringTie program directory*.

++StringTie Arguments Used++:  
- --merge - Distinct from the assembly usage mode used above, in the merge mode, StringTie takes as input a list of GTF/GFF files and merges/assembles these transcripts into a non-redundant set of transcripts.
- -p - Specify number of processers
- -G - Specify reference annotation file. With this option, StringTie assembles the transfrags from the input GTF files with the reference sequences
- -o - Name of output file
- <mergelist.txt> - File listing all filenames to be merged. Include full path.

```
nano /data/putnamlab/KITT/hputnam/20210921_BleachedPairs_TagSeq/scripts/merge_stringtie.sh
```

```
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=250GB
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/KITT/hputnam/20210921_BleachedPairs_TagSeq/Clean_Data/

#load packages 
module load StringTie/2.1.4-GCC-9.3.0
module load GffCompare/0.12.1-GCCcore-8.3.0

#make gtf_list.txt file 
ls *.gtf > gtf_mergelist.txt

stringtie --merge -p 8 -G ../refs/Mcap.GFFannotation.fixed.gff -o Mcap_stringtie_merged.gtf gtf_mergelist.txt

gffcompare -r ../refs/Mcap.GFFannotation.fixed.gff -o Mcap_gffCompared Mcap_stringtie_merged.gtf
```

```
sbatch /data/putnamlab/KITT/hputnam/20210921_BleachedPairs_TagSeq/scripts/merge_stringtie.sh
```

# Generate counts matrix
```
nano /data/putnamlab/KITT/hputnam/20210921_BleachedPairs_TagSeq/scripts/prepDE.sh
```

```
#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --mem=100GB
#SBATCH --account=putnamlab
#SBATCH -D /data/putnamlab/KITT/hputnam/20210921_BleachedPairs_TagSeq/Clean_Data/

#load packages 
module load Python/2.7.15-foss-2018b #Python 

#make gtf list text file 
for filename in ELS*.gtf; do echo $filename $PWD/$filename; done > list_path_GTF.txt

#Compile the gene count matrix 
python ../scripts/prepDE.py -g Mcap_Pairs_gene_count_matrix.csv -i list_path_GTF.txt

echo "Gene count matrix compiled." $(date)
```
```
sbatch /data/putnamlab/KITT/hputnam/20210921_BleachedPairs_TagSeq/scripts/prepDE.sh
```