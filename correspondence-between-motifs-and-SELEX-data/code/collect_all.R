library("readr")
library("dplyr")
library("tidyverse")
library("stringr")

rm(list=ls())
source("find_correspondence_between_motifs_and_SELEX_data_Jolma2013.R") #OK
source("find_correspondence_between_motifs_and_SELEX_data_Morgunova2015.R") #OK
source("find_correspondence_between_motifs_and_SELEX_data_Nitta2015.R") #OK
source("find_correspondence_between_motifs_and_SELEX_data_Yin2017.R") #OK
source("find_correspondence_between_motifs_and_SELEX_data_Jolma2015.R") #OK
source("find_correspondence_between_motifs_and_SELEX_data_Xie2025.R") #OK

# Jolma2013 ---------------------------------------------------------------


Jolma2013=read_delim("/projappl/project_2013895/motif_metadata/metadata_Jolma2013.tsv")

Jolma2013$study="Jolma2013"

Jolma2013%>% filter(is.na(Allas_SELEX_filename)) %>% pull(ID)
Jolma2013%>% filter(is.na(Allas_SELEX_background_filename)) %>% pull(ID)
# "EHF_HT-SELEX_TCTTGA20NGTG_AG_NACCCGGAAGTA_2_3"  
# "ELF3_HT-SELEX_TGACCT20NCCA_AG_NACCCGGAAGTAN_2_4"
# "ELF4_HT-SELEX_TGACTC20NTCA_AG_AACCCGGAAGTR_2_3" 


# Yin2017 -----------------------------------------------------------------

Yin2017=read_delim("/projappl/project_2013895/motif_metadata/metadata_Yin2017.tsv")

Yin2017$study="Yin2017"

Yin2017=Yin2017 %>% filter(experiment=="HT-SELEX")

Yin2017%>% filter(is.na(Allas_SELEX_filename)) %>% pull(ID)
Yin2017%>% filter(is.na(Allas_SELEX_background_filename)) %>% pull(ID)
#"ZSCAN5A_HT-SELEX_TCGCCC40NCAT_KR_NYGTCCCYCCCCAAANMN_2_2"  "ZSCAN31_HT-SELEX_TGGAGA40NCCA_KV_GCATAACKGCCCTGCKKCN_2_4"

# "NANOG_HT-SELEX_40NATTAATAATAA_KAL_TTAATKGN_1_3b0"               "HOXB6_HT-SELEX_40NTTAATATATAT_KAK_YTAATTRY_1_3b0"              
# "UNCX_HT-SELEX_40NATTTATTATAT_KAI_KTAATTAN_1_3b0"                "TCFL5_HT-SELEX_40NTTTAATTTATT_KAL_NCACGYGCAN_1_3b0"            
# "MYCN_HT-SELEX_40NAATTAATAAAT_KAL_NKCACGTGGN_1_3b0"              "NR2C2_HT-SELEX_40NATTAATTATTA_KAI_NRGGTCAN_1_3b0"              
# "ZSCAN5A_HT-SELEX_TCGCCC40NCAT_KR_NYGTCCCYCCCCAAANMN_2_2"        "ZSCAN31_HT-SELEX_TGGAGA40NCCA_KV_GCATAACKGCCCTGCKKCN_2_4"      
# "ZSCAN29_HT-SELEX_40NATTAAATTTAT_KAK_MMGYGTAGMCGKCTACACNN_2_3b0" "KLF17_HT-SELEX_40NTAATTAATTTA_KAM_NMCCACGCWCCCMYY_2_3b0"       
# "ZNF460_HT-SELEX_40NTTTATTTATAA_KAK_NAACGCCCCCCGN_1_3b0"        

# Morgunova2015 ---------------------------------------------------------------
 
Morgunova2015=read_delim("/projappl/project_2013895/motif_metadata/metadata_Morgunova2015.tsv")

Morgunova2015$study="Morgunova2015"

#Background missing, clone is DBD

# Nitta2015 ---------------------------------------------------------------

Nitta2015=read_delim("/projappl/project_2013895/motif_metadata/metadata_Nitta2015.tsv")

Nitta2015$study="Nitta2015"

Nitta2015%>% filter(is.na(Allas_SELEX_filename)) %>% pull(ID)
Nitta2015%>% filter(is.na(Allas_SELEX_background_filename)) %>% pull(ID)
#"SREBF2_HT-SELEX_TGAGAT20NGA_Y_NRTCACGCCAYN_1_3"      "ONECUT2_HT-SELEX_TGGGCG30NCGT_AH_NNCGATCRATAWNN_1_2" "ONECUT1_HT-SELEX_TAGCTC20NTCT_Y_CRATCRATAWN_1_3"  


# Jolma2015 ---------------------------------------------------------------

# cycle missing: BACH1_HT-SELEX_TTCCCC20NCCC_AL_ATGACTCAT_1_NA, this was removed

Jolma2015=read_delim("/projappl/project_2013895/motif_metadata/metadata_Jolma2015.tsv")

Jolma2015$study="Jolma2015"

Jolma2015%>% filter(is.na(Allas_SELEX_filename)) %>% pull(ID)
Jolma2015%>% filter(is.na(Allas_SELEX_background_filename)) %>% pull(ID)

# [1] "CEBPG_ATF4_CAP-SELEX_TGCGTC40NTTA_AAB_NNATGAYGCAAT_1_3b0" "BACH1_HT-SELEX_TTCCCC20NCCC_AL_ATGACTCAT_1_NA"            "FOS_HT-SELEX_TGAACT40NAAG_KR_NGATGACGTCATCR_2_4"         
# [4] "FOXA1_HT-SELEX_TTCTAA40NAAT_KN_TRNGTAAACA_1_3b1"          "IRF2_HT-SELEX_TTGCCC40NCTC_AAF_NAANCGAAASYR_1_3"          "NR1D2_HT-SELEX_TGAATT40NTAA_KR_TRGGTYASTAGGTCA_2_3"      
# [7] "RORB_HT-SELEX_TTCGGG40NGAG_KS_AANTAGGTCAGTAGGTCA_2_4"     "RORB_HT-SELEX_TTCGGG40NGAG_KS_AWNTAGGTCATGACCTANWT_2_4"   "SOX17_HT-SELEX_TATGCT40NACT_KO_ACCGAACAAT_1_4b2"


# Xie2025 -----------------------------------------------------------------

Xie2025=read_delim("/projappl/project_2013895/motif_metadata/metadata_Xie2025.tsv")

Xie2025$study="Xie2025"

Xie2025%>% filter(is.na(Allas_SELEX_filename)) %>% pull(ID)
Xie2025%>% filter(is.na(Allas_SELEX_background_filename)) %>% pull(ID)

# "VSX2_TBX6_TTTCGG40NAGA_YAAII_NAGGTGTTAATTN_1_3"   "ARNT_TGGTCG40NTAA_YT_NGTCACGTGACN_1_3"           
# "ESR2_TTCATA40NTCA_YT_NGGTCANNNNNNNNTGACCN_1_3"    "FOXK2_TTGAAC40NTAGC_YJI_NRTAAAYAAAYAN_1_3"       
# "FOXM1_TGACTG40NATG_YT_NGTAAACAAYRN_1_3"           "FOXM1_TGACTG40NATG_YT_NRYAAAYAAACAN_1_3"         
# "FOXM1_TGACTG40NATG_YT_NYGCATCAYAACRN_1_3"         "FOXP2_TCGATT40NTTG_YT_NRCGTAAACAAN_1_3"          
# "HELT_TAAGCC40NAGT_YT_NTCACGTGAC_1_3"              "OTP_TAGATG40NTAT_YPIII_NTAATTRNNNNTAATTRN_1_3"   
# "SOX11_TAGCCC40NGTA_YT_NACAATNNNNATTGTN_1_3"       "TCF23_TGTTTA40NACG_YT_NGCCATTTGGTN_1_3"          
# "THHEX_TCGAG40NCATT_YT_NCAATTNNNNNNNNNNAATTGN_1_3"

# Combine all -------------------------------------------------------------

Jolma2013$cycle=as.character(Jolma2013$cycle)
Jolma2015$cycle=as.character(Jolma2015$cycle)
Morgunova2015$cycle=as.character(Morgunova2015$cycle) #THis is not a SELEX motif, do not add
Nitta2015$cycle=as.character(Nitta2015$cycle)
Yin2017$cycle=as.character(Yin2017$cycle)
Xie2025$cycle=as.character(Xie2025$cycle)

SELEX_data_for_motifs=dplyr::bind_rows(Jolma2013, Jolma2015, Nitta2015, Yin2017, Xie2025) #3635

#For how many there are both SELEX signal and background signal

SELEX_data_for_motifs %>% filter(!is.na(CSC_SELEX_filename) & !is.na(CSC_SELEX_background_filename)) %>% nrow() #3596

tmp=SELEX_data_for_motifs %>% filter(!is.na(CSC_SELEX_filename) & !is.na(CSC_SELEX_background_filename))
tmp$study=NULL
#write_delim(tmp, "/projappl/project_2013895/motif_metadata/motifs_with_SELEX_signal_and_background.tsv", delim="\t")
#write_delim(tmp, "/projappl/project_2013895/motif_metadata/motifs_with_SELEX_signal_and_background_06102025.tsv", delim="\t") #3594
write_delim(tmp, "/projappl/project_2013895/motif_metadata/motifs_with_SELEX_signal_and_background_19032026.tsv", delim="\t") #3596


# Add the info to the metadata --------------------------------------------

metadata=read_delim("/projappl/project_2013895/motif_metadata/metadata_final.tsv")

setdiff(names(metadata), names(tmp))
setdiff(names(tmp), names(metadata))

metadata_combined=metadata %>% left_join(SELEX_data_for_motifs %>% select(ID, cycle_background, unique_background, 
                                                        CSC_SELEX_filename,
                                                        CSC_SELEX_background_filename,
                                                        Allas_SELEX_filename,
                                                        Allas_SELEX_background_filename, 
                                                        run_accession_signal,
                                                        study_accession_signal,               
                                                        secondary_study_accession_signal,
                                                        sample_accession_signal,
                                                        secondary_sample_accession_signal,
                                                        experiment_accession_signal,          
                                                        submission_accession_signal,
                                                        read_count_signal,
                                                        base_count_signal,
                                                        fastq_ftp_signal,                     
                                                        submitted_ftp_signal, sra_ftp_signal,
                                                        run_accession_background,
                                                        study_accession_background,           
                                                        secondary_study_accession_background,
                                                        sample_accession_background,
                                                        secondary_sample_accession_background,
                                                        experiment_accession_background,
                                                        submission_accession_background,
                                                        read_count_background,base_count_background,
                                                        fastq_ftp_background,                 
                                                        submitted_ftp_background,
                                                        sra_ftp_background
                                                        ), by="ID")

metadata_combined <- metadata_combined %>%
  relocate(cycle_background, .after = cycle)

metadata_combined <- metadata_combined %>%
  relocate(unique_background, .after = cycle_background)

metadata_combined$clone[metadata_combined$study=="Morgunova"]="DBD"

#write_delim(metadata_combined, "/projappl/project_2013895/motif_metadata/metadata_final_with_SELEX_data.tsv", delim="\t")
#write_delim(metadata_combined, "/projappl/project_2013895/motif_metadata/metadata_final_with_SELEX_data_06102025.tsv", delim="\t")
write_delim(metadata_combined, "/projappl/project_2013895/motif_metadata/metadata_final_with_SELEX_data_19032026.tsv", delim="\t")


missing=metadata_combined %>% filter(experiment!="Methyl-HT-SELEX") %>% filter(is.na(CSC_SELEX_filename)| is.na(CSC_SELEX_background_filename)) %>% 
  select(ID, study, experiment, CSC_SELEX_filename, CSC_SELEX_background_filename)

#write_delim(missing, "/projappl/project_2013895/motif_metadata/motifs_with_missing_SELEX_data.tsv", delim="\t")
#write_delim(missing, "/projappl/project_2013895/motif_metadata/motifs_with_missing_SELEX_data_06102025.tsv", delim="\t") # 41 missing
write_delim(missing, "/projappl/project_2013895/motif_metadata/motifs_with_missing_SELEX_data_19032026.tsv", delim="\t") #39 missing


