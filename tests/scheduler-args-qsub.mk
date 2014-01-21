# Scheduler arguments for test lib
#
# COMMENTS WITH #
# CAUTION: DON'T WRITE COMMENTS AFTER = SIGN, THEY ARE ASSIGNED TO THE VARIABLE

# Choose the scheduler, gnu-make-job-scheduler only supports qsub and sbatch
SCHEDULER=qsub

# standard options
# The mkdir at the end is for the output file of slurm in case the directory
# doesn't exist. 
SCHEDULER_STD_OPT=-V -N $@ -d $(shell pwd) -e $(shell pwd)/$@-pbs.err -o $(shell pwd)/$@-pbs.out $(shell mkdir -p $(@D))

## Individual options for each rule
#SCHEDULER_GUNZIP_OPT=-t 1-00:00:00 -p core
#SCHEDULER_QTRIM_OPT=-t 06:00:00 -p core
#
## velvet
#SCHEDULER_VELVETH_OPT=-t 01-00:00:00 -p node
## Scheduler velvetg options, also used for metavelvetg
#SCHEDULER_VELVETG_OPT=-t 02:00:00 -p node
#SCHEDULER_METAVELVETG_OPT=$(SCHEDULER_VELVETG_OPT)
#
## Minimus2
#SCHEDULER_MINIMUS2_VELVET_OPT=-t 01-00:00:00 -p node
#SCHEDULER_MINIMUS2_METAVELVET_OPT=$(SCHEDULER_MINIMUS2_VELVET_OPT)
#SCHEDULER_MINIMUS2_RAY_OPT=$(SCHEDULER_MINIMUS2_VELVET_OPT)
#
## Newbler
#SCHEDULER_NEWBLER_VELVET_OPT=-t 01-00:00:00 -p node
#SCHEDULER_NEWBLER_METAVELVET_OPT=$(SCHEDULER_NEWBLER_VELVET_OPT)
#SCHEDULER_NEWBLER_RAY_OPT=-$(SCHEDULER_NEWBLER_VELVET_OPT)

# ray
SCHEDULER_RAY_OPT=-l procs=8,walltime=01:00:00 
# should be short, only moves the Scaffolds.fasta file to the proper output dir
SCHEDULER_RAY_MV_SCAF_OPT=-l procs=1,walltime=00:01:00

# bambus2
#SCHEDULER_BAMBUS2_OPT=-t 01-00:00:00 -p node
