### Indicators for running specific tables or figures
### Change the indicator to FALSE to skip an experiment.
run_Figure1_Table1_TableB1=TRUE
run_Table2=TRUE
run_Table3_Figure2_Figure3_TableB2=TRUE

######################################
#### Figure 1, Table 1, Table B.1 ####
######################################
if(run_Figure1_Table1_TableB1==TRUE){
  source("Part 1/Figure 1 & Table 1 & Table B1.R")
}

#################
#### Table 2 ####
#################
if(run_Table2==TRUE){
  source("Part 2/Table 2.R")
}

################################################
#### Table 3, Figure 2, Figure 3, Table B.2 ####
################################################
if(run_Table3_Figure2_Figure3_TableB2==TRUE){
  source("Part 3/Table 3 & Figure 2 & Figure 3 & Table B2.R")
}