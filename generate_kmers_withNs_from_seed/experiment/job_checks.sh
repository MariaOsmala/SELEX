sacct \
  --format=JobID,Submit,Cluster,State,Elapsed,ElapsedRaw,ExitCode,User,Group,QOS,AllocCPUS,NNodes,NTasks,TotalCPU,ReqMem,MaxRSS,TIMELIMIT,ElapsedRaw,CPUTime,TotalCPU \
  --parsable2 --units=K \
  -j 30362167,30362189,30362272,30362467,30362880  |
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
}' > jobs_summary.tsv


sacct -j 30362167,30362189,30362272,30362467,30362880 -P --units=K \
  -o JobID,Submit,Cluster,State,Elapsed,ElapsedRaw,ExitCode,User,Group,QOS,AllocCPUS,NNodes,NTasks,TotalCPU,ReqMem,MaxRSS,TIMELIMIT,ElapsedRaw,CPUTime,TotalCPU \
| awk -F'|' 'NR==1 || $1 ~ /\.0$/' > jobs_summary_steps0.csv