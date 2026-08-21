library(readr)
library(tidyverse)
library(igraph)
library(RCy3)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggseqlogo)
library(gridExtra)
library(grid)
library(magick)

rm(list=ls())

#This is done, cytospace should be up to date 1.7.2026

TFBS_path="/Users/osmalama/projects/TFBS/"

metadata=read_delim(paste0("/Users/osmalama/projects/TFBS/Data/SELEX-motif-collection/Supplementary_Table_1_Submission_July2026.tsv"), 
                    delim="\t")

metadata %>% select(type) %>% table(useNA = "always") 
metadata %>% filter(representative=="yes") %>% select(type) %>% table(useNA = "always") #these are only representatives,

#ssDNA binding are counted as monomers
#composite       dimeric     monomeric       spacing ssDNA binding    tetrameric      trimeric       unknown          <NA> 
#  485           186           212           308            18             9             3            11             0 

metadata %>% filter(representative=="yes") %>% nrow() #1232

metadata %>% select(included_in_S2A) %>% table()

metadata %>% filter(experiment!="CAP-SELEX") %>% select(TF1_family) %>% table(useNA="always") %>% as.data.frame()

metadata %>% filter(experiment!="CAP-SELEX"& is.na(TF1_family)) %>% pull(ID)
#THHEX_TCGAG40NCATT_YT_NCAATTNNNNNNNNNNAATTGN_1_3 Homeodomain

metadata %>% filter(experiment!="CAP-SELEX"& TF1_family=="Znf") %>% pull(ID)
#"XPA_HT-SELEX_TCACGC40NACT_KX_NCACCTCACAN_1_4" Znf_XPA

#These are sorted_Lambert2018_families
tmp=metadata %>% filter(experiment=="CAP-SELEX") %>% select(TF1_family, TF2_family) %>%
rowwise() %>%
  mutate(family = paste0(sort(c(TF1_family, TF2_family)), collapse = "_"))

tmp %>% select(family) %>% table(useNA="always") %>% as.data.frame()
      
cytoscapePing()
cytoscapeVersionInfo()
#apiVersion cytoscapeVersion 
#"v1"         "3.10.4" 

setwd("~/projects/TFBS/code/motif-networks")

main_network_sid=getNetworkSuid()

network_names=getNetworkList(getSUID=FALSE)
network_SUIDs=getNetworkList(getSUID=TRUE)



# 1) Get node IDs (the default key column is usually "name")
nodes <- getTableColumns(table = "node", columns = "name", network = main_network_sid)

# 2) Make your new values (must align to nodes$name)
nodes = nodes %>% left_join(metadata %>% select(ID, type), by=c("name"="ID"))
nodes$type[is.na(nodes$type)]=""


# 4) Load/merge into Cytoscape node table 
loadTableData(nodes, data.key.column = "name", table = "node", table.key.column = "name", network = main_network_sid)

#Check that the network is correct: 
nodes <- getTableColumns(table = "node", network = main_network_sid)
nodes <- nodes %>% filter(node_type=="motif")

(nodes %>% select(representative)) %>% table() #ok

nodes %>% filter(representative=="YES") %>% select(type) %>% table(useNA="always")

nodes %>% filter(representative=="YES") %>% select(TF1_copies) %>% table(useNA="always")
nodes %>% select(Lambert2018_families) %>% table(useNA="always") %>% as.data.frame()

# Lambert2018_families Freq
# 1                     AP-2   17
# 2                   BED ZF    2
# 3                     bHLH  144
# 4                     bZIP  144
# 5                  C2H2 ZF  262
# 6         C2H2 ZF; AT hook    4
# 7                  CCCH ZF    2
# 8                    CENPB    1
# 9                      CSD    2
# 10        CUT; Homeodomain   20
# 11                      DM    8
# 12                     E2F   20
# 13                    EBF1    1
# 14                     Ets   97
# 15            Ets; AT hook    8
# 16                Forkhead   87
# 17                    GATA   18
# 18                     GCM    8
# 19              Grainyhead   15
# 20                     HMG   14
# 21                 HMG/Sox   74
# 22             Homeodomain  566 !
# 23 Homeodomain; Paired box   37
# 24        Homeodomain; POU   76
# 25                     HSF   24
# 26                     IRF   31
# 27                MADS box   10
# 28                    MEIS    5
# 29                Myb/SANT    9
# 30                     NRF    1
# 31        Nuclear receptor  170
# 32                     p53    5
# 33              Paired box    9
# 34                     POU    4
# 35                Prospero    2
# 36                     Rel   27
# 37                     RFX   20
# 38                     RRM    1
# 39                    Runt    9
# 40                    SAND    8
# 41                    SMAD    9
# 42                   T-box   53
# 43                     TEA    7
# 44                    TFAP    3
# 45                 Znf_XPA    1
# 46                    <NA>    0

# Check it:
#getTableColumns("node", c("name", "score"), network = main_network_sid)


#Then for TF_pair network

cytoscapePing()
cytoscapeVersionInfo()
#apiVersion cytoscapeVersion 
#"v1"         "3.10.4" 

setwd("~/projects/TFBS/code/motif-networks")

main_network_sid=getNetworkSuid()

network_names=getNetworkList(getSUID=FALSE)
network_SUIDs=getNetworkList(getSUID=TRUE)


# 1) Get node IDs (the default key column is usually "name")
nodes <- getTableColumns(table = "node", columns = "name", network = main_network_sid)

# 2) Make your new values (must align to nodes$name)

nodes = nodes %>% left_join(metadata %>% select(ID, type), by=c("name"="ID"))

nodes$type[is.na(nodes$type)]=""

#rename column 

#nodes <- nodes %>% 
#  rename(type = type_review)


# 4) Load/merge into Cytoscape node table (creates "score" column if missing)
loadTableData(nodes, data.key.column = "name", table = "node", table.key.column = "name", network = main_network_sid)

#Check nodes

#Check that the network is correct: 
nodes <- getTableColumns(table = "node", network = main_network_sid)
nodes <- nodes %>% filter(node_type=="motif")

(nodes %>% select(representative)) %>% table() #ok
nodes  %>% select(type) %>% table(useNA="always") #ok
# composite   spacing      <NA> 
#   1441       457         0 
nodes %>% filter(representative=="YES") %>% select(type) %>% table(useNA="always")
#composite   spacing      <NA> 
#  485       308         0 

nodes %>% select(sorted_Lambert2018_families) %>% table(useNA="always") %>% as.data.frame()
sorted_Lambert2018_families Freq
# 1                                  AP-2_bHLH    6
# 2                      AP-2_CUT; Homeodomain    5
# 3                                   AP-2_E2F    1
# 4                                   AP-2_Ets    2
# 5                              AP-2_Forkhead    2
# 6                           AP-2_Homeodomain    5
# 7                                 AP-2_T-box    1
# 8                                  bHLH_bHLH   14
# 9                                  bHLH_bZIP    2
# 10                     bHLH_CUT; Homeodomain   12
# 11                                  bHLH_E2F    4
# 12                                  bHLH_Ets   93
# 13                             bHLH_Forkhead   36
# 14                                  bHLH_GCM   15
# 15                          bHLH_Homeodomain   59
# 16                     bHLH_Homeodomain; POU    1
# 17                                  bHLH_IRF    2
# 18                             bHLH_Myb/SANT    8
# 19                                  bHLH_RFX    6
# 20                                 bHLH_Runt    3
# 21                                 bHLH_SAND    1
# 22                                bHLH_T-box   14
# 23                                  bHLH_TEA   29
# 24                                 bZIP_bZIP   21
# 25                     bZIP_CUT; Homeodomain    2
# 26                                  bZIP_Ets   61
# 27                             bZIP_Forkhead   11
# 28                                  bZIP_GCM    4
# 29                        bZIP_GntR-type HTH    1
# 30                              bZIP_HMG/Sox    8
# 31                          bZIP_Homeodomain   17
# 32                     bZIP_Homeodomain; POU    1
# 33                                  bZIP_IRF    3
# 34                             bZIP_Prospero    1
# 35                                  bZIP_Rel    5
# 36                                  bZIP_RFX    2
# 37                                bZIP_T-box   23
# 38                                  bZIP_TEA   12
# 39                               C2H2 ZF_Ets   17
# 40                          C2H2 ZF_Forkhead    1
# 41                           C2H2 ZF_HMG/Sox    2
# 42                       C2H2 ZF_Homeodomain    1
# 43                               C2H2 ZF_RFX    3
# 44                              C2H2 ZF_Runt    1
# 45                             C2H2 ZF_T-box    4
# 46                               C2H2 ZF_TEA    3
# 47                      CUT; Homeodomain_E2F    2
# 48                      CUT; Homeodomain_Ets   14
# 49                 CUT; Homeodomain_Forkhead    7
# 50                     CUT; Homeodomain_GATA    4
# 51                      CUT; Homeodomain_GCM    2
# 52                  CUT; Homeodomain_HMG/Sox    3
# 53              CUT; Homeodomain_Homeodomain   27
# 54         CUT; Homeodomain_Homeodomain; POU    4
# 55                      CUT; Homeodomain_IRF    1
# 56                 CUT; Homeodomain_MADS box    1
# 57                 CUT; Homeodomain_Myb/SANT    1
# 58                      CUT; Homeodomain_Rel    2
# 59                      CUT; Homeodomain_RFX    2
# 60                    CUT; Homeodomain_T-box    4
# 61                      CUT; Homeodomain_TEA    1
# 62                                   E2F_Ets    3
# 63                              E2F_Forkhead    3
# 64                                   E2F_GCM    1
# 65                           E2F_Homeodomain    1
# 66                                 E2F_T-box    6
# 67                                   Ets_Ets    7
# 68                              Ets_Forkhead   97
# 69                                  Ets_GATA    5
# 70                                   Ets_GCM   23
# 71                               Ets_HMG/Sox    6
# 72                           Ets_Homeodomain  223
# 73               Ets_Homeodomain; Paired box    3
# 74                      Ets_Homeodomain; POU    7
# 75                                   Ets_IRF   16
# 76                              Ets_Myb/SANT    1
# 77                            Ets_Paired box    6
# 78                                   Ets_RFX    6
# 79                                  Ets_Runt    3
# 80                                 Ets_T-box   50
# 81                                   Ets_TEA   29
# 82                             Forkhead_GATA    3
# 83                              Forkhead_GCM    8
# 84                          Forkhead_HMG/Sox    2
# 85                      Forkhead_Homeodomain   33
# 86          Forkhead_Homeodomain; Paired box    1
# 87                 Forkhead_Homeodomain; POU   13
# 88                              Forkhead_HSF    1
# 89                              Forkhead_IRF    1
# 90                         Forkhead_MADS box    1
# 91                 Forkhead_Nuclear receptor    1
# 92                       Forkhead_Paired box    2
# 93                             Forkhead_Runt    2
# 94                            Forkhead_T-box    9
# 95                              Forkhead_TEA    7
# 96                              GATA_HMG/Sox    1
# 97                          GATA_Homeodomain    1
# 98              GATA_Homeodomain; Paired box    1
# 99                                  GATA_IRF    1
# 100                          GATA_Paired box    1
# 101                               GATA_T-box    8
# 102                                 GATA_TEA    5
# 103                              GCM_HMG/Sox    5
# 104                          GCM_Homeodomain   18
# 105                                GCM_T-box    5
# 106                        GntR-type HTH_TEA    1
# 107                      HMG/Sox_Homeodomain   26
# 108          HMG/Sox_Homeodomain; Paired box    1
# 109                 HMG/Sox_Homeodomain; POU    9
# 110                            HMG/Sox_T-box    8
# 111                              HMG/Sox_TEA    4
# 112                  Homeodomain_Homeodomain   93
# 113      Homeodomain_Homeodomain; Paired box   31
# 114             Homeodomain_Homeodomain; POU   15
# 115                          Homeodomain_IRF   10
# 116                     Homeodomain_MADS box    7
# 117                     Homeodomain_Myb/SANT    2
# 118             Homeodomain_Nuclear receptor    3
# 119                   Homeodomain_Paired box   33
# 120                     Homeodomain_Prospero    4
# 121                          Homeodomain_Rel   10
# 122                          Homeodomain_RFX    5
# 123                         Homeodomain_Runt    6
# 124                         Homeodomain_SAND    1
# 125                        Homeodomain_T-box  260
# 126                          Homeodomain_TEA   69
# 127 Homeodomain; Paired box_Homeodomain; POU    1
# 128         Homeodomain; Paired box_Myb/SANT    1
# 129             Homeodomain; Paired box_Runt    1
# 130            Homeodomain; Paired box_T-box    8
# 131              Homeodomain; Paired box_TEA    3
# 132                Homeodomain; POU_MADS box    1
# 133                     Homeodomain; POU_Rel    1
# 134                    Homeodomain; POU_Runt    1
# 135                   Homeodomain; POU_T-box    9
# 136                     Homeodomain; POU_TEA    2
# 137                                IRF_T-box    5
# 138                                  IRF_TEA    1
# 139                             MADS box_RFX    1
# 140                           MADS box_T-box    1
# 141                             MADS box_TEA    1
# 142                           Myb/SANT_T-box    3
# 143                   Nuclear receptor_T-box    4
# 144                     Nuclear receptor_TEA    2
# 145                         Paired box_T-box    9
# 146                           Paired box_TEA    3
# 147                            Prospero_Runt    1
# 148                                RFX_T-box    1
# 149                                  RFX_TEA    2
# 150                               Runt_T-box    3
# 151                                T-box_TEA    9
# 152                                     <NA>    0