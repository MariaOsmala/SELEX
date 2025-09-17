#!/bin/bash
#SBATCH --job-name=kmers
#SBATCH --output=kmers_outs/kmers_%A_%a.out
#SBATCH --error=kmers_errs/kmers_%A_%a.err
#SBATCH --account=project_2007567
#SBATCH --partition=longrun
#SBATCH --ntasks=1
#SBATCH --time=7-00:00:00
#SBATCH --mem-per-cpu=20G #
#SBATCH --cpus-per-task=1
#SBATCH --gres=nvme:30
#SBATCH --array=1-10 #174 #1-3933/3636

#sacct --format JobID%-20,State -j 28590553 | grep FAILED

#29406957_1-10



# 29343847_1-300


# 29345862_301-600

# 29345938_601-900
# 29359880_901-1000 

# 29347707_1-10 1001-1010
# 29347902_11-300 1011-1300
# 29347983_301-400 1301-1400
# 29348008_401-700 1401-1700
# 29348079 1701-2000
# 29355999 2001-2200
# 29356094 2200-2500
# 29356233 2501-2800
#29357154_2801-3000
# 29357891_3001-3200
# 29357964_3201-3400
# 29358009_3401-636




#Failed cases

# 29360086, COMPLETED
#21,63,180,294,297,300,302,306,384,386,409,411,619,658,685,686,706,730,738,819,830,832,841,842,869,892,893

#29360117, this is still running
#1324

#29360117,29360086

# sacct \
#   --format=JobID,Submit,Cluster,State,Elapsed,ElapsedRaw,ExitCode,User,Group,QOS,AllocCPUS,NNodes,NTasks,TotalCPU,ReqMem,MaxRSS,TIMELIMIT,ElapsedRaw,CPUTime,TotalCPU \
#   --parsable2 --units=K \
#   -j 29343847,29345862,29345938,29347707,29359880,29347902,29347983,29348008,29348079,29355999,29356094,29356233,29357154,29357891,29357964,29358009 |
# awk -F'|' '
# function jobbase(id) { sub(/\..*$/,"",id); return id }
# function worse(s1,s2) {
#     order["FAILED"]=4; order["CANCELLED"]=3; order["TIMEOUT"]=3
#     order["COMPLETED"]=2; order["RUNNING"]=1; order["PENDING"]=0
#     return (order[s1]>order[s2]?s1:s2)
# }
# NR==1 { header=$0; next }
# {
#   base=jobbase($1)
# 
#   # Strings: keep first if unset
#   if (!(base in submit))   submit[base]=$2
#   if (!(base in cluster))  cluster[base]=$3
#   state[base]   = (state[base]=="" ? $4 : worse(state[base], $4))
#   if (!(base in exitcode)) exitcode[base]=$6
#   if (!(base in user))     user[base]=$7
#   if (!(base in group))    group[base]=$8
#   if (!(base in qos))      qos[base]=$9
#   if (!(base in reqmem))   reqmem[base]=$15
#   if (!(base in timelimit)) timelimit[base]=$17
# 
#   # Numeric: take max
#   if($5+0  > elapsed_hms[base]) elapsed_hms[base]=$5
#   if($6+0  > elapsedraw[base])  elapsedraw[base]=$6
#   if($13+0 > ntasks[base])      ntasks[base]=$13
#   if($14+0 > totalcpu[base])    totalcpu[base]=$14
#   if($16+0 > maxrss[base])      maxrss[base]=$16
#   if($18+0 > elapsedraw[base])  elapsedraw[base]=$18
#   if($19+0 > cputime[base])     cputime[base]=$19
#   if($20+0 > totalcpu2[base])   totalcpu2[base]=$20
# 
#   if($10+0 > alloccpus[base])   alloccpus[base]=$10
#   if($11+0 > nnodes[base])      nnodes[base]=$11
# }
# END {
#   OFS="\t"
#   print "JobID","Submit","Cluster","State","ExitCode","User","Group","QOS","AllocCPUS","NNodes","NTasks","TotalCPU","ReqMem","MaxRSS","TimeLimit","ElapsedRaw","CPUTime"
#   for (j in state)
#     print j,submit[j],cluster[j],state[j],exitcode[j],user[j],group[j],qos[j],alloccpus[j],nnodes[j],ntasks[j],totalcpu[j],reqmem[j],maxrss[j],timelimit[j],elapsedraw[j],cputime[j]
# }' > jobs_summary.tsv


# Account             AdminComment        AllocCPUS           AllocNodes
# AllocTRES           AssocID             AveCPU              AveCPUFreq
# AveDiskRead         AveDiskWrite        AvePages            AveRSS
# AveVMSize           BlockID             Cluster             Comment
# Constraints         ConsumedEnergy      ConsumedEnergyRaw   Container
# CPUTime             CPUTimeRAW          DBIndex             DerivedExitCode
# Elapsed             ElapsedRaw          Eligible            End
# ExitCode            Extra               FailedNode          Flags
# GID                 Group               JobID               JobIDRaw
# JobName             Layout              Licenses            MaxDiskRead
# MaxDiskReadNode     MaxDiskReadTask     MaxDiskWrite        MaxDiskWriteNode
# MaxDiskWriteTask    MaxPages            MaxPagesNode        MaxPagesTask
# MaxRSS              MaxRSSNode          MaxRSSTask          MaxVMSize
# MaxVMSizeNode       MaxVMSizeTask       McsLabel            MinCPU
# MinCPUNode          MinCPUTask          NCPUS               NNodes
# NodeList            NTasks              Partition           Planned
# PlannedCPU          PlannedCPURAW       Priority            QOS
# QOSRAW              Reason              ReqCPUFreq          ReqCPUFreqGov
# ReqCPUFreqMax       ReqCPUFreqMin       ReqCPUS             ReqMem
# ReqNodes            ReqTRES             Reservation         ReservationId
# Start               State               StdErr              StdIn
# StdOut              Submit              SubmitLine          Suspended
# SystemComment       SystemCPU           Timelimit           TimelimitRaw
# TotalCPU            TRESUsageInAve      TRESUsageInMax      TRESUsageInMaxNode
# TRESUsageInMaxTask  TRESUsageInMin      TRESUsageInMinNode  TRESUsageInMinTask
# TRESUsageInTot      TRESUsageOutAve     TRESUsageOutMax     TRESUsageOutMaxNode
# TRESUsageOutMaxTask TRESUsageOutMin     TRESUsageOutMinNode TRESUsageOutMinTask
# TRESUsageOutTot     UID                 User                UserCPU
# WCKey               WCKeyID             WorkDir






#Request local storage using the --gres flag in the job submission:
#--gres=nvme:<local_storage_space_per_node> #SBATCH --gres=nvme:50
#The amount of space is given in GB (check maximum sizes from the list above). 
#For example, to request 100 GB of storage, use option --gres=nvme:100. 
#The local storage reservation is on a per node basis.

#Use the environment variable $LOCAL_SCRATCH in your batch job scripts to access the local storage on each node.

# Load r-env
module load r-env/442

# Clean up .Renviron file in home directory
if test -f ~/.Renviron; then
    sed -i '/TMPDIR/d' ~/.Renviron
fi

# Specify a temp folder path
echo "TMPDIR=/scratch/project_2013895/tmp///" >> ~/.Renviron


jid="${SLURM_ARRAY_JOB_ID:+${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}}"
trap 'seff "${jid:-$SLURM_JOBID}"' EXIT

# Run the R script
srun apptainer_wrapper exec Rscript --no-save ../code/generate_kmers_from_seed_better.R $SLURM_ARRAY_TASK_ID 0


