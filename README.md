MetAssemble
===========

Content

1. Overview
2. Dependencies
3. Installation
4. Usage

1. Overview
===========

MetAssemble is a pipeline that runs several metagenomic assembly strategies
combining Velvet, Meta-Velvet, Minimus2, Ray and Bambus2 on Illumina paired end
reads. The pipeline was originally developed to validate the performance of the
individual strategies, but can be used to perform the assembly strategies
without validation as well. The pipeline is written in GNU make and not very
user friendly for the average user, but if you are familiar with GNU make you
shouldn't have too many troubles getting it to run. The only other metagenomics
assembly pipeline that I am aware of is [metAMOS], which seems to be an effort
towards a more user-friendly approach if you are looking for that. A reason for
using MetAssemble instead is because it allows one to schedule parts of the
assembly pipeline with sbatch. Different steps in the assembly pipeline require
different resources. Velvet for instance runs on only one node, whereas Ray
runs over multiple. MetAssemble allows you to specify resource usage per rule
with [gnu-make-job-scheduler]. Furthermore GNU make makes sure intermediate
output files don't have to be recomputed in case of an error.

[metAMOS]: https://github.com/treangen/metAMOS/
[gnu-make-job-scheduler]: https://github.com/inodb/gnu-make-job-scheduler

2. Dependencies
===============
Dependencies need to be installed by oneself. There is no automated way to do this at the moment.
One can however check if the dependencies are met by running

    bash test/dependencies/test_dependencies.sh

Do note that it is not necessary to install all programs if you only want to do
a subset of the assemblies that MetAssemble covers. MetAssemble requires the
following programs to perform all different assemblies:

Supported input:

- Illumina fastq CASAVA v1.8 paired end reads

Running the MetAssemble pipeline (scripts/Makefile) requires

- GNU make (tested on v3.81)

The Makefile features four steps of the metagenomic assembly pipeline:

1. Read processing.
    - Quality trimming
        * [sickle](https://github.com/najoshi/sickle) 

2. Assembling contigs
    - De Bruijn graph
        * [Velvet 1.2.01](http://www.ebi.ac.uk/~zerbino/velvet/)
        * [Meta-Velvet 1.2.01](http://metavelvet.dna.bio.keio.ac.jp/)
        * [Ray](http://denovoassembler.sourceforge.net/)

3. Merging contigs
    - With cd-hit and minimus2. See [Angus](http://ged.msu.edu/angus/metag-assembly-2011/velvet-multik.html).
        * [cd-hit-est](http://weizhong-lab.ucsd.edu/cd-hit/) NOTE: Increase MAX_SEQ when compiling. I put it at 3000000.
        * [Minimus2 from Amos](http://sourceforge.net/apps/mediawiki/amos/index.php?title=Minimus2)
        * [MUMmer 3.23](http://sourceforge.net/projects/mummer/files/)
    - Cut up contigs and merge with Newbler RunAssembly 2.6
        * scripts/process-reads/cut-up-fasta.py requires [Biopython](http://biopython.org/wiki/Main_Page)
        * Newbler RunAssembly 2.6 (COMMERCIAL)

4. Scaffolding
    - Construct linkage information by mapping reads to contigs
        * [BWA](http://bio-bwa.sourceforge.net/)
    - Scaffold contigs
        * [Bambus2](http://sourceforge.net/apps/mediawiki/amos/index.php?title=Bambus2)


3. Installation
===============
After installing all the dependencies the scripts should work as is. You can do
a test run with `cd test && make test`, which downloads a small set from the
HMP project and runs a subset of all different assembly strategies in the
MetAssemble pipeline.

4. Usage
========
See the example in test/. There is a test/Makefile and a test/Makefile-sbatch
which set some input paramaters and then include scripts/Makefile and
scripts/Makefile-sbatch respectively. Hopefully that is clear enough to help
you understand how to run your own subset of the available assembly strategies.
If you want to change the resource usage per rule, change Makefile-sbatch
accordingly. In the future I might add automatic computation of the resource
usage. For assembly this is unfortunately still a problem, since it depends on
the complexity of your sample and not just the filesize. The specified resource
usage is for a library of ~1M and a mixed community of 60 bacteria and
archaeae.

To see which assemblies have been created:

    make echoexisting

All assemblies, created or not:

    make echoall

To create all:

    make all

Only show commands:

    make -n all

Only make velvet:

    make velvet

Schedule rules with sbatch:

    make -f Makefile-sbatch all

For more rules check in the scripts/parameters.mk file.
