# BenchMapL
Snakemake pipeline for compare multiple long read mappers on viruses dataset 

*<h2>Authors </h2>* 

Students of the Bioinformatics Master of the University of Bordeaux 

Thomas BAUDEAU – thomas.baudeau@etu.u-bordeaux.fr 

<h2>Install </h2>*

**Installation**

Clone the repository:

<pre>
git clone https://github.com/ThomasBaudeau/BenchMapL <br> 
cd BenchMapL
</pre>


**Requirement**

Snakemake workflows : https://snakemake.readthedocs.io/en/stable/getting_started/installation.html <br> 
Conda : https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html <br> 


## How to start with Visugraph :

The structue of the workflow is build as following:

    ├── .gitignore
    ├── README.md
    ├── LICENSE.md
    ├── workflow
    │   ├── benchmarks
    |   │   ├── benchresult.txt
    |   │   └── ...
    |   ├── data
    |   |   ├── sample
    |   |   │   ├── species_length_error-rates.fasta
    |   |   │   └── ...
    |   |   └── ref_species.fasta
    │   ├── envs
    |   │   ├── tool1.yaml
    |   │   └── tool2.yaml
    │   ├── scripts
    |   │   ├── script1.py
    |   │   └── script2.R
    │   ├── helps_tool
    |   │   ├── tools_helpdoc.txt
    |   │   └── ...
    │   ├── plots
    |   │   ├── plot1.pdf
    |   │   └── plot2.pdf
    |   ├── resu
    |   │   ├── resu.bam
    |   │   └── ...
    |   └── Snakefile
    ├── workflow2
    |   ├── data
    |   |   ├── supl
    |   |   │   ├── model.model
    |   |   │   └── ...
    |   |   └── ref_species.fasta
    │   ├── envs
    |   │   ├── tool1.yaml
    |   │   └── tool2.yaml
    │   ├── scripts
    |   │   ├── script1.py
    |   │   └── script2.R
    │   ├── result
    |   │   ├── nano
    |   │   └── pbsim
    |   └── Snakefile
    ├── config
        ├── config.yaml
        └── config2.yaml



 #### BenchMapL Usage:

 ##### Workflow2 : Data Generation Part

