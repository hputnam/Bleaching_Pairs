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

### Sequence QC
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

