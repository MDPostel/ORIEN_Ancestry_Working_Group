{\rtf1\ansi\ansicpg1252\cocoartf2636
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs24 \cf0 #!/usr/bin/bash\
#SBATCH --partition=oneweek\
#SBATCH --mem=16GB\
#SBATCH --cpus-per-task=8\
#SBATCH --time=7-00:00:00\
#SBATCH --mail-user=mpostel@usc.edu\
#SBATCH --mail-type=FAIL\
#SBATCH --output=slurm-%j.out\
#SBATCH --job-name=eigen\
\
\
echo "\
genotypename:    /project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/subset/AIMs/isec.bed\
snpname:         /project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/subset/AIMs/isec.bim\
indivname:       /project/davidwcr_264/Projects/Mackenzie/ORIEN/orien_ancestry/subset/AIMs/isec.fam\
outputformat:    EIGENSTRAT\
genotypeoutname: Intermember_ancestry_eigen.geno\
snpoutname:      Intermember_ancestry_eigen.snp\
indivoutname:   Intermember_ancestry_eigen.ind\
familynames:     NO " > parfile\
\
convertf -p parfile\
\
}