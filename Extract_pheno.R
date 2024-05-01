library(dplyr)
library(tidyr)

tb=read.csv("/share/ScratchGeneral/zheqia/proj_DR/mocha/cache/DR2010-2011_ANZRAGP2_SampleSheet.csv", sep=",", skip=10, na.strings=c(""), h=T) # 2380
tb$gtc=paste0(tb$SentrixBarcode_A,"_", tb$SentrixPosition_A)
tb=tb[,c("Sample_ID", "Sample_Group", "Gender", "gtc")]
df=read.csv("/share/ScratchGeneral/zheqia/proj_DR/DR_phenotype.csv", sep=",", na.strings=c(""), h=T)

# check numbers
dim(df)                        # 3782
dim(tb)                        # 2380
length(unique(df$ID))          # 3781
length(unique(tb$Sample_ID))   # 2371

names(tb)[1]="ID"
# check duplicates
tb[tb$ID %in% tb[which(duplicated(tb[,"ID", drop=F])),]$ID,]
# some of the duplicates are belong to different Sample Groups, but seems to have re-genotyped idats data
tb=tb[!duplicated(tb),]        # 2378

length(intersect(tb$ID, df$ID))               # 1705
length(unique(intersect(tb$ID, df$ID)))       # 1705

# merge dataframes
out=inner_join(tb, df, by="ID")               # 1712
out[out$ID %in% out[which(duplicated(out[,"ID", drop=F])),]$ID,]

# integrate mLOX, mLOY, autosomal mCA
auto=read.table("/share/ScratchGeneral/zheqia/proj_DR/mocha/mosaic_phenos/DR_from_gtc_biallelic.auto.lines_autosomal_mCAs", h=T)
mlox=read.table("/share/ScratchGeneral/zheqia/proj_DR/mocha/mosaic_phenos/DR_from_gtc_biallelic.X_loss.lines_mLOX", h=F)
mloy=read.table("/share/ScratchGeneral/zheqia/proj_DR/mocha/mosaic_phenos/DR_from_gtc_biallelic.Y_loss.lines_mLOY", h=F)
rmv=read.table("/share/ScratchGeneral/zheqia/proj_DR/mocha/mosaic_phenos/DR_from_gtc_biallelic.remove.lines_samples_to_be_excluded_from_assoc_analyses", h=F)
names(auto)="gtc"
names(mlox)="gtc"
names(mloy)="gtc"
names(rmv)="gtc"
auto$auto_mCA=1
mlox$mLOX=1
mloy$mLOY=1
rmv$remove=1

out=left_join(out, auto)
out=left_join(out, mlox)
out=left_join(out, mloy)
out=left_join(out, rmv)
out=out %>% mutate(auto_mCA=replace_na(auto_mCA, 0)) %>%
        mutate(mLOX=replace_na(mLOX, 0)) %>%
        mutate(mLOY=replace_na(mLOY, 0)) %>%
        mutate(remove=replace_na(remove, 0))
out = filter(out, remove !=1)

out = out %>%
  add_count(ID) %>%
  filter((n == 1) | (n > 1 & Sample_Group == "DR_2011")) %>%
  select(-n)

write.table(out, "DR_pheno_mCA.tsv", row.names=F, col.names=T, quote=F, sep='\t')
out=out[,c("ID", "Sample_Group", "Gender", "gtc", "Sex", "Deceased", "Ethnicity", "EverSmoked", "ExSmoker", "YearsofSmoking", "Noofcigarettesday", "bmi", "vtdr", "t1dm", "t2dm", "diabtype", "age_at_recruitment", "auto_mCA", "mLOX", "mLOY", "remove")]
write.table(out, "DR_pheno_mCA_selected.tsv", row.names=F, col.names=T, quote=F, sep='\t')