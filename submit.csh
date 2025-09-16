  #!/bin/tcsh
  #BSUB -W 1440                  ## time required by minutes, the job get killed if goes beyond 
  #BSUB -n 10                    ## number of cores required
  #BSUB -R span[hosts=1]         ## `hosts=1` all cores in one node. Or `ptile=1` number of cores per node 
  #BSUB -R "rusage[mem=20GB]"    ## memories required
  ##BSUB -x                      ## for exclusive use of a node, use only if memory intensive
  #BSUB -J simulation[1-4]         ## job name AND job array
  #BSUB -o out.%J
  #BSUB -e err.%J
  
  /usr/local/usrapps/bjreich/sprim/AAA/bin/Rscript simulation.R $LSB_JOBINDEX
