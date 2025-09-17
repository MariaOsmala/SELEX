#!/bin/bash
#SBATCH --job-name=Ham12_kmers
#SBATCH --output=Ham12_kmers_outs/Ham12_kmers_%A_%a.out
#SBATCH --error=Ham12_kmers_errs/Ham12_kmers_%A_%a.err
#SBATCH --account=project_2006203
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --time=3-00:00:00
#SBATCH --mem-per-cpu=100G #
#SBATCH --cpus-per-task=1
#SBATCH --gres=nvme:30
#SBATCH --array=601-700 #1-3636

#sacct --format JobID%-20,State -j 28590553 | grep FAILED

#29364750_1-10
#29367071_[11-100]
#29368948_101-300
#29400536_301-600
#29400548_601-900
#29400561_901-1000
#29400986_1-200 + 1000 29400986_[161-200]   CANCELLED+ 
#29407077_201-400 +1000 29407077_[306-400]   CANCELLED
#29414563_401-500 +1000
#29414564_501-600 +1000

#29400150
#"8,17,21,65,67,186,200,264,272,294,297,300"

sacct \
  --format=JobID,Submit,Cluster,State,Elapsed,ElapsedRaw,ExitCode,User,Group,QOS,AllocCPUS,NNodes,NTasks,TotalCPU,ReqMem,MaxRSS,TIMELIMIT,ElapsedRaw,CPUTime,TotalCPU \
  --parsable2 --units=K \
  -j 29364750,29367071,29368948,29400536,29400548,29400561,29400986,29407077 |
awk -F'|' '
function jobbase(id) { sub(/\..*$/,"",id); return id }
function worse(s1,s2) {
    order["FAILED"]=4; order["CANCELLED"]=3; order["TIMEOUT"]=3
    order["COMPLETED"]=2; order["RUNNING"]=1; order["PENDING"]=0
    return (order[s1]>order[s2]?s1:s2)
}
NR==1 { header=$0; next }
{
  base=jobbase($1)

  # Strings: keep first if unset
  if (!(base in submit))   submit[base]=$2
  if (!(base in cluster))  cluster[base]=$3
  state[base]   = (state[base]=="" ? $4 : worse(state[base], $4))
  if (!(base in exitcode)) exitcode[base]=$6
  if (!(base in user))     user[base]=$7
  if (!(base in group))    group[base]=$8
  if (!(base in qos))      qos[base]=$9
  if (!(base in reqmem))   reqmem[base]=$15
  if (!(base in timelimit)) timelimit[base]=$17

  # Numeric: take max
  if($5+0  > elapsed_hms[base]) elapsed_hms[base]=$5
  if($6+0  > elapsedraw[base])  elapsedraw[base]=$6
  if($13+0 > ntasks[base])      ntasks[base]=$13
  if($14+0 > totalcpu[base])    totalcpu[base]=$14
  if($16+0 > maxrss[base])      maxrss[base]=$16
  if($18+0 > elapsedraw[base])  elapsedraw[base]=$18
  if($19+0 > cputime[base])     cputime[base]=$19
  if($20+0 > totalcpu2[base])   totalcpu2[base]=$20

  if($10+0 > alloccpus[base])   alloccpus[base]=$10
  if($11+0 > nnodes[base])      nnodes[base]=$11
}
END {
  OFS="\t"
  print "JobID","Submit","Cluster","State","ExitCode","User","Group","QOS","AllocCPUS","NNodes","NTasks","TotalCPU","ReqMem","MaxRSS","TimeLimit","ElapsedRaw","CPUTime"
  for (j in state)
    print j,submit[j],cluster[j],state[j],exitcode[j],user[j],group[j],qos[j],alloccpus[j],nnodes[j],ntasks[j],totalcpu[j],reqmem[j],maxrss[j],timelimit[j],elapsedraw[j],cputime[j]
}' > kmers_and_Hamming12.tsv


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
srun apptainer_wrapper exec Rscript --no-save ../code/generate_kmers_and_Hamming_1_and_2_from_seed.R $SLURM_ARRAY_TASK_ID 1000


