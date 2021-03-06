# TSDfinder workflow
TSDfinder is developed on the basis of [SPLITREADER](https://github.com/LeanQ/SPLITREADER) which is a bioinformatic pipeline dedicated to the discovery of non-reference TE insertions with Target Site Duplications (TSDs) through the use of Illumina short genome sequence reads and a reference genome.

### 1. Clone the repository
```
git clone https://github.com/bioxfu/TSDfinder
```

### 2. Initiate the project
```
source init.sh
```

### 3. Download SRA data
```
bash download_sra.sh
```

### 4. Create *config.yaml* and *Snakefile* based on the examples

### 5. Dry run the workflow to check any mistakes
```
./dry_run.sh
```

### 6. Start the workflow
```
./run.sh
```

### 7. Check the workflow progress in *nohup.out* 

### 8. Remove the temporary files
```
./clean.sh
```

