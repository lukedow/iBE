## ------------- MULTIPLE FILE IMPORT -------------------- ##
## ----Open multiple text files into separate dataframes----
# Set working directory to location of .txt files
setwd("/path/to/files")

library(stringi)

# Create list of file names
filenames <- gsub("\\.tab$","", list.files(pattern="\\.tab$"))

# Run import for loop
for(i in filenames)
  assign(i, read.table(paste(i, ".tab", sep=""), header = TRUE))

# Rename columns in all using for loop
for(df in filenames)
  assign(df, setNames(get(df),  c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","DP","QD","AD")))

# Separate AD value into reference (REF) and variant (VAR) allele counts
for(df in filenames)
  assign(df, (separate((get(df)), AD, c("WT","VAR"), ",")))

## ----Convert and calculate VAF----
#Convert counts to numeric values and calculate VAF
RA690_SW_D8_filtered_variants$WT <- RA690_SW_D8_filtered_variants$WT %>% as.numeric()
RA690_SW_D8_filtered_variants$VAR <- RA690_SW_D8_filtered_variants$VAR %>% as.numeric()
RA690_SW_D8_filtered_variants$VAF <- RA690_SW_D8_filtered_variants$VAR/(RA690_SW_D8_filtered_variants$WT + RA690_SW_D8_filtered_variants$VAR)

RA690_ON_D8_filtered_variants$WT <- RA690_ON_D8_filtered_variants$WT %>% as.numeric()
RA690_ON_D8_filtered_variants$VAR <- RA690_ON_D8_filtered_variants$VAR %>% as.numeric()
RA690_ON_D8_filtered_variants$VAF <- RA690_ON_D8_filtered_variants$VAR/(RA690_ON_D8_filtered_variants$WT + RA690_ON_D8_filtered_variants$VAR)

RA690_OFF_D4_filtered_variants$WT <- RA690_OFF_D4_filtered_variants$WT %>% as.numeric()
RA690_OFF_D4_filtered_variants$VAR <- RA690_OFF_D4_filtered_variants$VAR %>% as.numeric()
RA690_OFF_D4_filtered_variants$VAF <- RA690_OFF_D4_filtered_variants$VAR/(RA690_OFF_D4_filtered_variants$WT + RA690_OFF_D4_filtered_variants$VAR)

# Create CHROM:POS ID
RA690_SW_D8_filtered_variants$POS_ID <- as.factor(paste(RA690_SW_D8_filtered_variants$CHROM, RA690_SW_D8_filtered_variants$POS, sep = ":"))
RA690_ON_D8_filtered_variants$POS_ID <- as.factor(paste(RA690_ON_D8_filtered_variants$CHROM, RA690_ON_D8_filtered_variants$POS, sep = ":"))
RA690_OFF_D4_filtered_variants$POS_ID <- as.factor(paste(RA690_OFF_D4_filtered_variants$CHROM, RA690_OFF_D4_filtered_variants$POS, sep = ":"))

# Subset for read depth
RA690_SW_D8_filtered_variants <- RA690_SW_D8_filtered_variants[RA690_SW_D8_filtered_variants$DP > 20,]
RA690_ON_D8_filtered_variants <- RA690_ON_D8_filtered_variants[RA690_ON_D8_filtered_variants$DP > 20,]
RA690_OFF_D4_filtered_variants <- RA690_OFF_D4_filtered_variants[RA690_OFF_D4_filtered_variants$DP > 20,]

# Define length of REF and ALT alleles
RA690_SW_D8_filtered_variants$REF_length <- stri_length(RA690_SW_D8_filtered_variants$REF)
RA690_SW_D8_filtered_variants$ALT_length <- stri_length(RA690_SW_D8_filtered_variants$ALT)

RA690_ON_D8_filtered_variants$REF_length <- stri_length(RA690_ON_D8_filtered_variants$REF)
RA690_ON_D8_filtered_variants$ALT_length <- stri_length(RA690_ON_D8_filtered_variants$ALT)

RA690_OFF_D4_filtered_variants$REF_length <- stri_length(RA690_OFF_D4_filtered_variants$REF)
RA690_OFF_D4_filtered_variants$ALT_length <- stri_length(RA690_OFF_D4_filtered_variants$ALT)

# Specify edit type (based ON 1:1 substitution OR "other = indel")
RA690_SW_D8_filtered_variants$MUT_type <- ifelse(RA690_SW_D8_filtered_variants$REF_length == "1" 
                                                 & RA690_SW_D8_filtered_variants$ALT_length == "1", 
                                                 "SNP","indel")

RA690_ON_D8_filtered_variants$MUT_type <- ifelse(RA690_ON_D8_filtered_variants$REF_length == "1" 
                                                 & RA690_ON_D8_filtered_variants$ALT_length == "1", 
                                                 "SNP","indel")

RA690_OFF_D4_filtered_variants$MUT_type <- ifelse(RA690_OFF_D4_filtered_variants$REF_length == "1" 
                                                 & RA690_OFF_D4_filtered_variants$ALT_length == "1", 
                                                 "SNP","indel")

# Specify edit detail
RA690_SW_D8_filtered_variants$MUT_detail <- ifelse(RA690_SW_D8_filtered_variants$MUT_type == "SNP" &
                                                   RA690_SW_D8_filtered_variants$REF == "C" 
                                                 & RA690_SW_D8_filtered_variants$ALT == "T" | 
                                                   RA690_SW_D8_filtered_variants$REF == "G" 
                                                 & RA690_SW_D8_filtered_variants$ALT == "A", 
                                                 "C>T",ifelse(RA690_SW_D8_filtered_variants$MUT_type == "SNP" &
                                                                RA690_SW_D8_filtered_variants$REF == "C" 
                                                              & RA690_SW_D8_filtered_variants$ALT == "G" | 
                                                                RA690_SW_D8_filtered_variants$REF == "G" 
                                                              & RA690_SW_D8_filtered_variants$ALT == "C",
                                                              "C>G", ifelse(RA690_SW_D8_filtered_variants$MUT_type == "SNP" &
                                                                              RA690_SW_D8_filtered_variants$REF == "C" 
                                                                            & RA690_SW_D8_filtered_variants$ALT == "A" | 
                                                                              RA690_SW_D8_filtered_variants$REF == "G" 
                                                                            & RA690_SW_D8_filtered_variants$ALT == "T",
                                                                            "C>A",ifelse(RA690_SW_D8_filtered_variants$MUT_type == "SNP" &
                                                                                           RA690_SW_D8_filtered_variants$REF == "A" 
                                                                                         & RA690_SW_D8_filtered_variants$ALT == "T" | 
                                                                                           RA690_SW_D8_filtered_variants$REF == "T" 
                                                                                         & RA690_SW_D8_filtered_variants$ALT == "A",
                                                                                         "A>T",ifelse(RA690_SW_D8_filtered_variants$MUT_type == "SNP" &
                                                                                                        RA690_SW_D8_filtered_variants$REF == "A" 
                                                                                                      & RA690_SW_D8_filtered_variants$ALT == "C" | 
                                                                                                        RA690_SW_D8_filtered_variants$REF == "T" 
                                                                                                      & RA690_SW_D8_filtered_variants$ALT == "G",
                                                                                                      "A>C",ifelse(RA690_SW_D8_filtered_variants$MUT_type == "SNP" &
                                                                                                                        RA690_SW_D8_filtered_variants$REF == "A" 
                                                                                                                      & RA690_SW_D8_filtered_variants$ALT == "G" | 
                                                                                                                        RA690_SW_D8_filtered_variants$REF == "T" 
                                                                                                                      & RA690_SW_D8_filtered_variants$ALT == "C",
                                                                                                                   "A>G",""))))))

RA690_ON_D8_filtered_variants$MUT_detail <- ifelse(RA690_ON_D8_filtered_variants$MUT_type == "SNP" &
                                                     RA690_ON_D8_filtered_variants$REF == "C" 
                                                   & RA690_ON_D8_filtered_variants$ALT == "T" | 
                                                     RA690_ON_D8_filtered_variants$REF == "G" 
                                                   & RA690_ON_D8_filtered_variants$ALT == "A", 
                                                   "C>T",ifelse(RA690_ON_D8_filtered_variants$MUT_type == "SNP" &
                                                                  RA690_ON_D8_filtered_variants$REF == "C" 
                                                                & RA690_ON_D8_filtered_variants$ALT == "G" | 
                                                                  RA690_ON_D8_filtered_variants$REF == "G" 
                                                                & RA690_ON_D8_filtered_variants$ALT == "C",
                                                                "C>G", ifelse(RA690_ON_D8_filtered_variants$MUT_type == "SNP" &
                                                                                RA690_ON_D8_filtered_variants$REF == "C" 
                                                                              & RA690_ON_D8_filtered_variants$ALT == "A" | 
                                                                                RA690_ON_D8_filtered_variants$REF == "G" 
                                                                              & RA690_ON_D8_filtered_variants$ALT == "T",
                                                                              "C>A",ifelse(RA690_ON_D8_filtered_variants$MUT_type == "SNP" &
                                                                                             RA690_ON_D8_filtered_variants$REF == "A" 
                                                                                           & RA690_ON_D8_filtered_variants$ALT == "T" | 
                                                                                             RA690_ON_D8_filtered_variants$REF == "T" 
                                                                                           & RA690_ON_D8_filtered_variants$ALT == "A",
                                                                                           "A>T",ifelse(RA690_ON_D8_filtered_variants$MUT_type == "SNP" &
                                                                                                          RA690_ON_D8_filtered_variants$REF == "A" 
                                                                                                        & RA690_ON_D8_filtered_variants$ALT == "C" | 
                                                                                                          RA690_ON_D8_filtered_variants$REF == "T" 
                                                                                                        & RA690_ON_D8_filtered_variants$ALT == "G",
                                                                                                        "A>C",ifelse(RA690_ON_D8_filtered_variants$MUT_type == "SNP" &
                                                                                                                       RA690_ON_D8_filtered_variants$REF == "A" 
                                                                                                                     & RA690_ON_D8_filtered_variants$ALT == "G" | 
                                                                                                                       RA690_ON_D8_filtered_variants$REF == "T" 
                                                                                                                     & RA690_ON_D8_filtered_variants$ALT == "C",
                                                                                                                     "A>G",""))))))


RA690_OFF_D4_filtered_variants$MUT_detail <- ifelse(RA690_OFF_D4_filtered_variants$MUT_type == "SNP" &
                                                     RA690_OFF_D4_filtered_variants$REF == "C" 
                                                   & RA690_OFF_D4_filtered_variants$ALT == "T" | 
                                                     RA690_OFF_D4_filtered_variants$REF == "G" 
                                                   & RA690_OFF_D4_filtered_variants$ALT == "A", 
                                                   "C>T",ifelse(RA690_OFF_D4_filtered_variants$MUT_type == "SNP" &
                                                                  RA690_OFF_D4_filtered_variants$REF == "C" 
                                                                & RA690_OFF_D4_filtered_variants$ALT == "G" | 
                                                                  RA690_OFF_D4_filtered_variants$REF == "G" 
                                                                & RA690_OFF_D4_filtered_variants$ALT == "C",
                                                                "C>G", ifelse(RA690_OFF_D4_filtered_variants$MUT_type == "SNP" &
                                                                                RA690_OFF_D4_filtered_variants$REF == "C" 
                                                                              & RA690_OFF_D4_filtered_variants$ALT == "A" | 
                                                                                RA690_OFF_D4_filtered_variants$REF == "G" 
                                                                              & RA690_OFF_D4_filtered_variants$ALT == "T",
                                                                              "C>A",ifelse(RA690_OFF_D4_filtered_variants$MUT_type == "SNP" &
                                                                                             RA690_OFF_D4_filtered_variants$REF == "A" 
                                                                                           & RA690_OFF_D4_filtered_variants$ALT == "T" | 
                                                                                             RA690_OFF_D4_filtered_variants$REF == "T" 
                                                                                           & RA690_OFF_D4_filtered_variants$ALT == "A",
                                                                                           "A>T",ifelse(RA690_OFF_D4_filtered_variants$MUT_type == "SNP" &
                                                                                                          RA690_OFF_D4_filtered_variants$REF == "A" 
                                                                                                        & RA690_OFF_D4_filtered_variants$ALT == "C" | 
                                                                                                          RA690_OFF_D4_filtered_variants$REF == "T" 
                                                                                                        & RA690_OFF_D4_filtered_variants$ALT == "G",
                                                                                                        "A>C",ifelse(RA690_OFF_D4_filtered_variants$MUT_type == "SNP" &
                                                                                                                       RA690_OFF_D4_filtered_variants$REF == "A" 
                                                                                                                     & RA690_OFF_D4_filtered_variants$ALT == "G" | 
                                                                                                                       RA690_OFF_D4_filtered_variants$REF == "T" 
                                                                                                                     & RA690_OFF_D4_filtered_variants$ALT == "C",
                                                                                                                     "A>G",""))))))


# Label transitions and transversion
RA690_SW_D8_filtered_variants$MUT_TSTV <- ifelse(RA690_SW_D8_filtered_variants$MUT_type == "SNP" &
                                                  RA690_SW_D8_filtered_variants$MUT_detail == "C>T" |
                                                   RA690_SW_D8_filtered_variants$MUT_detail == "A>G", 
                                                 "transition",ifelse(RA690_SW_D8_filtered_variants$MUT_type == "indel","","transversion"))
RA690_SW_D8_filtered_variants$condition <- "SW"

RA690_ON_D8_filtered_variants$MUT_TSTV <- ifelse(RA690_ON_D8_filtered_variants$MUT_type == "SNP" &
                                                   RA690_ON_D8_filtered_variants$MUT_detail == "C>T" |
                                                   RA690_ON_D8_filtered_variants$MUT_detail == "A>G", 
                                                 "transition",ifelse(RA690_ON_D8_filtered_variants$MUT_type == "indel","","transversion"))
RA690_ON_D8_filtered_variants$condition <- "ON"

RA690_OFF_D4_filtered_variants$MUT_TSTV <- ifelse(RA690_OFF_D4_filtered_variants$MUT_type == "SNP" &
                                                   RA690_OFF_D4_filtered_variants$MUT_detail == "C>T" |
                                                   RA690_OFF_D4_filtered_variants$MUT_detail == "A>G", 
                                                 "transition",ifelse(RA690_OFF_D4_filtered_variants$MUT_type == "indel","","transversion"))
RA690_OFF_D4_filtered_variants$condition <- "OFF"

##---- rbind sample conditions ----
RA698_stack <- rbind(RA698_OFF_D4_filtered_variants,
                    RA698_ON_D8_filtered_variants,
                    RA698_SW_D8_filtered_variants)

RA691_stack <- rbind(RA691_OFF_D4_filtered_variants,
                     RA691_ON_D8_filtered_variants,
                     RA691_SW_D8_filtered_variants)

RA690_stack <- rbind(RA690_OFF_D4_filtered_variants,
                     RA690_ON_D8_filtered_variants,
                     RA690_SW_D8_filtered_variants)

# Eliminate values appearing more than 3 times
RA698_stack_thrice <- RA698_stack[RA698_stack$POS_ID %in% names(which(table(RA698_stack$POS_ID) <=3 )), ]
RA691_stack_thrice <- RA691_stack[RA691_stack$POS_ID %in% names(which(table(RA691_stack$POS_ID) <=3 )), ]
RA690_stack_thrice <- RA690_stack[RA690_stack$POS_ID %in% names(which(table(RA690_stack$POS_ID) <=3 )), ]

# Eliminate values appearing more than 2 times
RA698_stack_twice <- RA698_stack[RA698_stack$POS_ID %in% names(which(table(RA698_stack$POS_ID) <=2 )), ]
RA691_stack_twice <- RA691_stack[RA691_stack$POS_ID %in% names(which(table(RA691_stack$POS_ID) <=2 )), ]
RA690_stack_twice <- RA690_stack[RA690_stack$POS_ID %in% names(which(table(RA690_stack$POS_ID) <=2 )), ]

# Eliminate values appearing more than once
RA698_stack_unique <- RA698_stack[RA698_stack$POS_ID %in% names(which(table(RA698_stack$POS_ID) <=1 )), ]
RA691_stack_unique <- RA691_stack[RA691_stack$POS_ID %in% names(which(table(RA691_stack$POS_ID) <=1 )), ]
RA690_stack_unique <- RA690_stack[RA690_stack$POS_ID %in% names(which(table(RA690_stack$POS_ID) <=1 )), ]

# Eliminate values appearing more than 3 times
RA698_stack_common <- RA698_stack[RA698_stack$POS_ID %in% names(which(table(RA698_stack$POS_ID) == 3 )), ]
RA691_stack_common <- RA691_stack[RA691_stack$POS_ID %in% names(which(table(RA691_stack$POS_ID) == 3 )), ]
RA690_stack_common <- RA690_stack[RA690_stack$POS_ID %in% names(which(table(RA690_stack$POS_ID) == 3 )), ]

# Count number edits in each category
sum(RA698_stack_thrice$condition == "SW" & RA698_stack_thrice$DP > 20 & RA698_stack_thrice$MUT_detail == "C>T")
sum(RA691_stack_thrice$condition == "SW" & RA691_stack_thrice$DP > 20 & RA691_stack_thrice$MUT_detail == "C>T")
sum(RA690_stack_thrice$condition == "SW" & RA690_stack_thrice$DP > 20 & RA690_stack_thrice$MUT_detail == "C>T")

sum(RA698_stack_thrice$condition == "SW" & RA698_stack_thrice$DP > 20)

# plot unique variants
RA690_stack_unique[RA690_stack_unique$QD > 2 & RA690_stack_unique$DP >20 & RA690_stack_unique$MUT_detail == "C>T",] %>%
  mutate(condition = fct_relevel(condition, levels=c("OFF","ON","SW"))) %>% 
  ggplot() +
  geom_jitter( aes(x = condition, y = VAF, color = DP), size = 2, alpha = 0.6, width = 0.3, height = 0) +
  scale_color_viridis(option = "D")+
  ggtitle("RA690") +
  xlab("Treatment") +
  ylab("Fraction C-to-U editing") +
  theme_bw(base_size = 12)+
  facet_grid(cols = vars(MUT_detail))
# plot 'twice' variants
RA690_stack_twice[RA690_stack_twice$QD > 2 & RA690_stack_twice$DP >20 & RA690_stack_twice$MUT_detail == "C>T",] %>%
  mutate(condition = fct_relevel(condition, levels=c("OFF","ON","SW"))) %>% 
  ggplot() +
  geom_jitter( aes(x = condition, y = VAF, color = DP), size = 2, alpha = 0.6, width = 0.3, height = 0) +
  scale_color_viridis(option = "D")+
  ggtitle("RA690") +
  xlab("Treatment") +
  ylab("Fraction C-to-U editing") +
  theme_bw(base_size = 12)+
  facet_grid(cols = vars(MUT_detail))
# plot 'thrice' variants
RA690_stack_thrice[RA690_stack_twice$QD > 2 & RA690_stack_thrice$DP >20 & RA690_stack_thrice$MUT_detail == "C>T",] %>%
  mutate(condition = fct_relevel(condition, levels=c("OFF","ON","SW"))) %>% 
  ggplot() +
  geom_jitter( aes(x = condition, y = VAF, color = DP), size = 2, alpha = 0.6, width = 0.3, height = 0) +
  scale_color_viridis(option = "D")+
  ggtitle("RA690") +
  xlab("Treatment") +
  ylab("Fraction C-to-U editing") +
  theme_bw(base_size = 12)+
  facet_grid(cols = vars(MUT_detail))
# plot only common variants for all editing types
RA698_stack_common[RA698_stack_common$QD > 2 & RA698_stack_common$DP >20,] %>%
  mutate(condition = fct_relevel(condition, levels=c("OFF","ON","SW"))) %>% 
  ggplot() +
  geom_jitter( aes(x = condition, y = VAF, color = MUT_detail), size = 2, alpha = 1, width = 0.3, height = 0) +
  scale_color_viridis_d(option = "D")+
  ggtitle("RA698 | common") +
  xlab("Treatment") +
  ylab("Fraction editing") +
  theme_bw(base_size = 12)+
  facet_grid(cols = vars(MUT_detail))

# plot common excluded for all editing type
RA698_stack_twice[RA698_stack_twice$QD > 2 & RA698_stack_twice$DP >15 & RA698_stack_twice$MUT_detail == "C>T",] %>%
  mutate(condition = fct_relevel(condition, levels=c("OFF","ON","SW"))) %>% 
  ggplot() +
  geom_jitter( aes(x = condition, y = VAF, color = DP), size = 2, alpha = 1, width = 0.3, height = 0) +
  scale_color_viridis(option = "D")+
  ggtitle("RA698 | common excluded") +
  xlab("Treatment") +
  ylab("Fraction editing") +
  theme_bw(base_size = 12)+
  facet_grid(cols = vars(MUT_detail))

RA690_stack_twice[RA690_stack_twice$QD > 2 & RA690_stack_twice$DP >20,] %>%
  mutate(condition = fct_relevel(condition, levels=c("OFF","ON","SW"))) %>% 
  ggplot() +
  geom_jitter( aes(x = condition, y = VAF, color = MUT_detail), size = 2, alpha = 1, width = 0.3, height = 0) +
  scale_color_viridis_d()+
  ggtitle("RA690 | common variants excluded") +
  xlab("Condition") +
  ylab("Fraction editing") +
  theme_bw(base_size = 10)+
  facet_grid(cols = vars(MUT_type))

RA691_stack_twice[RA691_stack_twice$QD > 2 & RA691_stack_twice$DP >20,] %>%
  mutate(condition = fct_relevel(condition, levels=c("OFF","ON","SW"))) %>% 
  ggplot() +
  geom_jitter( aes(x = condition, y = VAF, color = MUT_detail), size = 2, alpha = 1, width = 0.3, height = 0) +
  scale_color_viridis_d()+
  ggtitle("RA691 | common variants excluded") +
  xlab("Condition") +
  ylab("Fraction editing") +
  theme_bw(base_size = 10)+
  facet_grid(cols = vars(MUT_type))

RA698_stack_twice[RA698_stack_twice$QD > 2 & RA698_stack_twice$DP >20,] %>%
  mutate(condition = fct_relevel(condition, levels=c("OFF","ON","SW"))) %>% 
  ggplot() +
  geom_jitter( aes(x = condition, y = VAF, color = MUT_detail), size = 2, alpha = 1, width = 0.3, height = 0) +
  scale_color_viridis_d()+
  ggtitle("RA698 | common variants excluded") +
  xlab("Condition") +
  ylab("Fraction editing") +
  theme_bw(base_size = 10)+
  facet_grid(cols = vars(MUT_type))

##---- cbind sample conditions ----
RA698_merge <- list(RA698_OFF_D4_filtered_variants,
                    RA698_ON_D8_filtered_variants,
                    RA698_SW_D8_filtered_variants) %>% reduce(full_join, by = "POS_ID")

# need to relabel some of the missing REF/ALT information
RA698_merge$REF <- ifelse(!is.na(RA698_merge$REF.x), paste(RA698_merge$REF.x),
                          ifelse(!is.na(RA698_merge$REF.y), paste(RA698_merge$REF.y),
                                        paste(RA698_merge$REF.y.y)))
RA698_merge$ALT <- ifelse(!is.na(RA698_merge$ALT.x), paste(RA698_merge$ALT.x),
                          ifelse(!is.na(RA698_merge$ALT.y), paste(RA698_merge$ALT.y),
                                        paste(RA698_merge$ALT.y.y)))

RA698_merge_VAF <- RA698_merge[,-c(1:11,14:30,32:40,43:48,50:55)]
RA698_merge_VAF <- RA698_merge_VAF[,c(2,4,5,1,3,6)]
colnames(RA698_merge_VAF) <- c("POS_ID","REF","ALT","OFF","ON_D8","SW_D8")

# Define length of REF and ALT alleles
RA698_merge_VAF$REF_length <- stri_length(RA698_merge_VAF$REF)
RA698_merge_VAF$ALT_length <- stri_length(RA698_merge_VAF$ALT)

# label mutation type
RA698_merge_VAF$MUT_type <- ifelse(RA698_merge_VAF$REF_length == "1" 
                                                 & RA698_merge_VAF$ALT_length == "1", 
                                                 "SNP","indel")

# label variant detail
RA698_merge_VAF$MUT_detail <- ifelse(RA698_merge_VAF$MUT_type == "SNP" &
                                       RA698_merge_VAF$REF == "C" 
                                     & RA698_merge_VAF$ALT == "T" | 
                                       RA698_merge_VAF$REF == "G" 
                                     & RA698_merge_VAF$ALT == "A", 
                                     "C>T",ifelse(RA698_merge_VAF$MUT_type == "SNP" &
                                                    RA698_merge_VAF$REF == "C" 
                                                  & RA698_merge_VAF$ALT == "G" | 
                                                    RA698_merge_VAF$REF == "G" 
                                                  & RA698_merge_VAF$ALT == "C",
                                                  "C>G", ifelse(RA698_merge_VAF$MUT_type == "SNP" &
                                                                  RA698_merge_VAF$REF == "C" 
                                                                & RA698_merge_VAF$ALT == "A" | 
                                                                  RA698_merge_VAF$REF == "G" 
                                                                & RA698_merge_VAF$ALT == "T",
                                                                "C>A",ifelse(RA698_merge_VAF$MUT_type == "SNP" &
                                                                               RA698_merge_VAF$REF == "A" 
                                                                             & RA698_merge_VAF$ALT == "T" | 
                                                                               RA698_merge_VAF$REF == "T" 
                                                                             & RA698_merge_VAF$ALT == "A",
                                                                             "A>T",ifelse(RA698_merge_VAF$MUT_type == "SNP" &
                                                                                            RA698_merge_VAF$REF == "A" 
                                                                                          & RA698_merge_VAF$ALT == "C" | 
                                                                                            RA698_merge_VAF$REF == "T" 
                                                                                          & RA698_merge_VAF$ALT == "G",
                                                                                          "A>C",ifelse(RA698_merge_VAF$MUT_type == "SNP" &
                                                                                                         RA698_merge_VAF$REF == "A" 
                                                                                                       & RA698_merge_VAF$ALT == "G" | 
                                                                                                         RA698_merge_VAF$REF == "T" 
                                                                                                       & RA698_merge_VAF$ALT == "C",
                                                                                                       "A>G",""))))))
# label transitions and transversions
RA698_merge_VAF$MUT_TSTV <- ifelse(RA698_merge_VAF$MUT_type == "SNP" &
                                     RA698_merge_VAF$MUT_detail == "C>T" |
                                     RA698_merge_VAF$MUT_detail == "A>G", 
                                                 "transition",ifelse(RA698_merge_VAF$MUT_type == "indel","","transversion"))

# Define common SNPs between sample groups
RA698_merge_VAF$status <- ifelse(!is.na(RA698_merge_VAF$OFF) & !is.na(RA698_merge_VAF$ON_D8) & !is.na(RA698_merge_VAF$SW_D8), "common",
                                 ifelse(is.na(RA698_merge_VAF$OFF) & !is.na(RA698_merge_VAF$ON_D8) & is.na(RA698_merge_VAF$SW_D8),"dox_only",
                                        ifelse(is.na(RA698_merge_VAF$OFF) & is.na(RA698_merge_VAF$ON_D8) & !is.na(RA698_merge_VAF$SW_D8),"SW_only",
                                               ifelse(is.na(RA698_merge_VAF$OFF) & !is.na(RA698_merge_VAF$ON_D8) & !is.na(RA698_merge_VAF$SW_D8),"Dox_SW",
                                                      ""))))
                                                                                      


RA698_merge_VAF[RA698_merge_VAF$status == "common",] %>% ggplot() + 
  geom_point(aes(x = OFF, y = ON_D8, color = MUT_type)) + 
  ggtitle("RA698 common variants - ON DOX vs OFF")+facet_grid(cols = vars(MUT_type))

# Comparing VAF for variants common among all conditions
RA698_merge_VAF[RA698_merge_VAF$status == "common",] %>% 
  ggscatter(x = "OFF", y = "ON_D8", 
          add = "reg.line", conf.int = TRUE, alpha = 1, color = "MUT_type", 
          cor.coef = TRUE, cor.method = "pearson",
          title = "RA698 organoids: common variants | ON vs OFF",
          xlab = "NO DOX", 
          ylab = "ON DOX")+
  theme_bw()+
  xlim(0,1)+
  ylim(0,1)

# Comparing VAF for variants common among all conditions
RA698_merge_VAF %>% 
  ggscatter(x = "OFF", y = "ON_D8", 
            add = "reg.line", conf.int = TRUE, alpha = 1, color = "MUT_type", 
            cor.coef = TRUE, cor.method = "pearson",
            title = "RA698 organoids: all variants | ON vs OFF",
            xlab = "NO DOX", 
            ylab = "ON DOX")+
  theme_bw()+
  xlim(0,1)+
  ylim(0,1)

RA690_VAF_unique <- RA690_merge_VAF %>% filter(RA690_merge_VAF$RA690_OFF_D4_VAF == "NA" | RA690_merge_VAF$RA690_ON_D8_VAF == "NA" | RA690_merge_VAF$RA690_ON_D8_VAF == "NA")


# plot common excluded for all editing type
RA698_stack_twice[RA698_stack_twice$QD > 2 & RA698_stack_twice$DP >15 & RA698_stack_twice$MUT_detail == "C>T",] %>% 
  ggplot(aes(x = condition, y = VAF, color = DP)) +
  geom_jitter(size = 2, alpha = 1, width = 0.3, height = 0) +
  scale_color_viridis(option = "D")+
  ggtitle("RA698 | common excluded") +
  xlab("Treatment") +
  ylab("Fraction editing") +
  theme_bw(base_size = 12)+
  stat_compare_means(comparisons = list(c("OFF","ON"),c("OFF","SW"),c("ON","SW")))
