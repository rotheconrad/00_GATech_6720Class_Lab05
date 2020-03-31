# 00_GATech_6720Class_Lab05

This is an assignment developed as Bioinformatics Lab 05 for Class 6720 - Environemntal Microbial Genomics at GA Tech.

The objective of this lab is to identify the most abundant population for which a metagenome assembled genome (MAG) can be recovered. To do this, we will assembled an Illumina sequenced metagenomic sample from Pensacola beach sands collected during the Macodomonis oil spill. We will start with our metagenome in interleaved fasta format. The paired reads in this file have already been quality controlled for us. Our task will be to assembled the reads into contigs using the IDBA-UD assembler, cluster the contigs into bins (MAGs) using MaxBin2 and find the most abundant MAG, and then identify the nearest taxonomic assignment as well as the completeness, contamination, and quality scores for each MAG using the Microbial Genomes Atlas (MiGA). This lab has instructions optimized for the MicrobiomeOS virtual machine (VM) running Ubuntu. The virtual machine can be downloaded [here](http://enve-omics.ce.gatech.edu/microbiomeos/).

This lab will require XXgb of disk space. I recommend to create and mount a virtual harddisk drive to the VM following these instructions:

Add a second hard disk to your VM to increase disk space
Takes 5 min or less.

0. Make sure your VM is powered off.
1. Open the Oracle VM VirtualBox Manager
2. Select the MicrobiomeOS and click Settings
3. Navigate to the Storage tab
4. Select Controller: SATA inside the Storage Devices box on the left and click the diskdrive icon with the green "+" (Not the CD icon)
5. Click Create
6. Select VHD (Virtual Hard Disk) and click continue
7. Select Fixed size and click continue
8. Select desired size (Need at least XXgb for lab 5) and click Create
9. Select the new drive we just created from beneath the Not Attached arrow and click Choose (It should be something like MicrobiomeOS_1.vhd)
10. We have now created and added a second virtual harddrive to our machine. Click Ok.
11. Start the MicrobeOS virtual machine and let it load.
12. From the left bar at the top select search machine and search for the word "Disks". You should find an application called Disks. Open this application.
13. In the Disks application you will see two disks and cd/dvd drive in the left column. The top disk should be the main disk and the second disk (or third etc) will be the new disk we added. Select the new disk. In the right frame at the bottom, the Device should read /dev/sdb (sdb depends on the number of disks on the system. sda is the main disk, sdb is the second, sdc is third and etc.)
14. Underneath Volumes select the gears icon and click format partition. Choose a name for the new disk such as Lab5, for Type select the Ext4 file system, and click format
15. If it asks you if you are sure, click format again. If it asks for a password, the password for the microbiomeOS is microbiome.
15. You should now see the play symbol next to the gear icon. You must now click the play play button. This will mount the partition.
16. Below the play button you will see a Contents label showing the path where your new disk is Mounted at. It should be something like /midea/microbiome/Lab5. This is the new folder to install and store files for Lab5.

## Step 00: Required tools :: Conda with Python 3.6+, IDBA-UD, and MaxBin2.

### Conda with Python 3.6+

From the [mini conda website](https://docs.conda.io/en/latest/miniconda.html), select the appropriate installer and download it. For the VM choose Linux installers Python 3.7 (most modern computers are 64-bit) and save the file. Open a terminal window and navigate to the Miniconda installer (try the Downloads folder). Once you locate the file, you can run it like so (change the file name if yours is different):

```bash
bash Miniconda3-latest-Linux-x86_64.sh
```

Follow the instructions and do a default installation. Answer yes where appropriate.

You can learn more about Conda environments [here](https://docs.conda.io/projects/conda/en/latest/user-guide/index.html).

### IDBA and MaxBin2

Utilizing the Conda system makes installing programs and their depencies much easier. For instance, looking at the [MaxBin2 README.txt file](https://sourceforge.net/projects/maxbin2/files/) you'll notice that installation requires some prerequisites and auxiliary software packages. But, thanks to the [BioConda Project](https://bioconda.github.io/) most commonly used bioinformatics software (along with all dependencies) can be quickly installed. You'll also notice that IDBA-UD is an auxilliary package for MaxBin2. This means that with one conda install command for MaxBin2 we get IDBA-UD as well.

```bash
# Create a conda environment for lab 5
conda create -n EnveomicsLab5
# Activate the lab 5 environment
conda activate EnveomicsLab5
# install MaxBin2
conda install -c bioconda maxbin2
```