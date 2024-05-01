library(dplyr)

# Load data
# This data is located at /share/ScratchGeneral/zheqia/proj_DR
pheno=read.csv("DR_pheno_mCA_selected.tsv", h=T, sep="\t")
identical(as.character(apply(pheno[,"Gender", drop=F], 2, 
                             function(x) substr(x,1,1))), pheno$Sex) # TRUE

# Concordance between EverSmoked and ExSmoker
table(pheno$EverSmoked,pheno$ExSmoker)
# 347 NeverSmoked were not ExSmokers
# 668 EverSmoked were ExSmokers
# 429 unkown Eversmoked were unkown for ExSmoker
# 159 EverSmoked were unkown for ExSmoker (this is OK, just treat them as EverSmoked)
# 24 NeverSmoked were unkown for ExSmoker (this is OK, just treat them as EverSmoked)
# 83 EverSmoked were not ExSmokers (this is a mistake)
filter(pheno, EverSmoked == "Yes" & ExSmoker == "No") # seems like these are heavy smokers...
# Aftering checking, it is reasonable to use EverSmoked as the smoking status

pheno=pheno %>% mutate(vtdr = case_when(
  is.na(vtdr) ~ NA_real_,             
  vtdr == "No VTDR" ~ 0,              
  vtdr == "VTDR" ~ 1)) %>% 
  mutate(Smoking = case_when(
    is.na(EverSmoked) ~ NA_real_,     
    EverSmoked == "No" ~ 0,
    EverSmoked == "Yes" ~ 1)) 
pheno$t1dm=as.factor(pheno$t1dm)
pheno$t2dm=as.factor(pheno$t2dm)
pheno=pheno %>% 
  mutate(t1dm = case_when(
    is.na(t1dm) ~ 0,             
    t1dm == "T1DM" ~ 1,
    TRUE ~ as.numeric(NA)
  )) %>%
  mutate(t2dm = case_when(
    is.na(t2dm) ~ 0,              
    t2dm == "T2DM" ~ 1,
    TRUE ~ as.numeric(NA)
  ))
pheno$t1dm=as.factor(pheno$t1dm)
pheno$t2dm=as.factor(pheno$t2dm)
pheno$vtdr=as.factor(pheno$vtdr)
pheno$Sex=as.factor(pheno$Sex)
pheno$Smoking=as.factor(pheno$Smoking)
pheno$diabtype=as.factor(pheno$diabtype)
pheno$age_at_recruitment=as.numeric(pheno$age_at_recruitment)

table(pheno$diabtype, pheno$vtdr)
mat=matrix(c(159,239,338,894), nrow = 2, byrow = TRUE,
          dimnames = list(c("Type 1 Diabetes", "Type 2 Diabetes"),
                          c("VTDR", "No VTDR")))
chisq.test(mat)

table(pheno$auto_mCA, pheno$vtdr)
mat=matrix(c(1069,480,67,18), nrow = 2, byrow = TRUE,
           dimnames = list(c("No", "Yes"),
                           c("No VTDR", "VTDR")))
chisq.test(mat)

table(pheno$mLOX, pheno$vtdr)
mat=matrix(c(1103,489,33,9), nrow = 2, byrow = TRUE,
           dimnames = list(c("No", "Yes"),
                           c("No VTDR", "VTDR")))
chisq.test(mat)

table(pheno$mLOY, pheno$vtdr)
mat=matrix(c(1032,466,104,32), nrow = 2, byrow = TRUE,
           dimnames = list(c("No", "Yes"),
                           c("No VTDR", "VTDR")))
chisq.test(mat)

# check missing data
library(Amelia)
missmap(pheno, main = "Missing values vs observed")

model=glm(vtdr ~ Sex + age_at_recruitment + auto_mCA + Smoking, data = pheno, family=binomial(link='logit'))
summary(model)
anova(model, test="Chisq")
coef_summary=summary(model)$coefficients
or=exp(coef_summary[, "Estimate"])
ci_lower=exp(coef_summary[, "Estimate"] - 1.96 * coef_summary[, "Std. Error"])
ci_upper=exp(coef_summary[, "Estimate"] + 1.96 * coef_summary[, "Std. Error"])
results=data.frame(
  OR = or,
  LowerCI = ci_lower,
  UpperCI = ci_upper
)
round(results, 2)

female=filter(pheno, Sex=="F")
model=glm(vtdr ~ age_at_recruitment + mLOX + Smoking, data = female, family=binomial(link='logit'))
summary(model)
anova(model, test="Chisq")
coef_summary=summary(model)$coefficients
or=exp(coef_summary[, "Estimate"])
ci_lower=exp(coef_summary[, "Estimate"] - 1.96 * coef_summary[, "Std. Error"])
ci_upper=exp(coef_summary[, "Estimate"] + 1.96 * coef_summary[, "Std. Error"])
results=data.frame(
  OR = or,
  LowerCI = ci_lower,
  UpperCI = ci_upper
)
round(results, 2)

male=filter(pheno, Sex=="M")
model=glm(vtdr ~ age_at_recruitment + mLOY + Smoking, data = male, family=binomial(link='logit'))
summary(model)
anova(model, test="Chisq")
coef_summary=summary(model)$coefficients
or=exp(coef_summary[, "Estimate"])
ci_lower=exp(coef_summary[, "Estimate"] - 1.96 * coef_summary[, "Std. Error"])
ci_upper=exp(coef_summary[, "Estimate"] + 1.96 * coef_summary[, "Std. Error"])
results=data.frame(
  OR = or,
  LowerCI = ci_lower,
  UpperCI = ci_upper
)
round(results, 2)

# Plotting
library(ggplot2)
library(ggpattern)

# df = data.frame(table(pheno$diabtype, pheno$vtdr))
# names(df) = c("Diabetes_Type", "VTDR", "Count")
# 
# ggplot(df, aes(x = Diabetes_Type, y = Count, fill = VTDR)) +
#   geom_col(width=0.45, alpha=0.8, position = "stack") +
#   scale_fill_manual(values = c("lightblue", "darkblue")) +
#   labs(x = "Diabetes", y = "", fill = "VTDR") +
#   coord_flip() +
#   theme_minimal()

# DM
tb = pheno
tb$diabtype = factor(tb$diabtype, levels = c("T2DM", "T1DM"), labels = c(0, 1))
tb$vtdr = factor(tb$vtdr)

tb = tb %>% 
  group_by(diabtype, vtdr) %>% 
  summarise(Count = n(), .groups = 'drop')
tb = filter(tb, !is.na(vtdr) & !is.na(diabtype))
tb$vtdr = factor(tb$vtdr, levels = c(0, 1), labels = c("No VTDR", "VTDR"))
tb$event="Types of\n    Diabetes"
names(tb)[1]="Status"
df=tb

# auto mCAs
tb = pheno
tb$auto_mCA = factor(tb$auto_mCA, levels = c(0, 1))
tb$vtdr = factor(tb$vtdr)

tb = tb %>% 
  group_by(auto_mCA, vtdr) %>% 
  summarise(Count = n(), .groups = 'drop')
tb = filter(tb, !is.na(vtdr))
tb$vtdr = factor(tb$vtdr, levels = c(0, 1), labels = c("No VTDR", "VTDR"))
tb$event="Autosomal\n    mCAs"
names(tb)[1]="Status"
df=rbind(df, tb)

# mLOX
tb = pheno
tb$mLOX = factor(tb$mLOX, levels = c(0, 1))
tb$vtdr = factor(tb$vtdr)

tb = tb %>% 
  group_by(mLOX, vtdr) %>% 
  summarise(Count = n(), .groups = 'drop')
tb = filter(tb, !is.na(vtdr))
tb$vtdr = factor(tb$vtdr, levels = c(0, 1), labels = c("No VTDR", "VTDR"))
tb$event="mLOX"
names(tb)[1]="Status"
df=rbind(df, tb)

# mLOY
tb = pheno
tb$mLOY = factor(tb$mLOY, levels = c(0, 1))
tb$vtdr = factor(tb$vtdr)

tb = tb %>% 
  group_by(mLOY, vtdr) %>% 
  summarise(Count = n(), .groups = 'drop')
tb = filter(tb, !is.na(vtdr))
tb$vtdr = factor(tb$vtdr, levels = c(0, 1), labels = c("No VTDR", "VTDR"))
tb$event="mLOY"
names(tb)[1]="Status"
df=rbind(df, tb)
df$event=factor(df$event, levels=c("mLOY", "mLOX", "Autosomal\n    mCAs", "Types of\n    Diabetes"))

# Add prevalence
tmp = df %>%
  group_by(event, Status) %>%
  summarise(TotalCount = sum(Count), .groups = 'drop')

df = df %>% left_join(tmp, by = c("event", "Status"))
df = df %>%
  mutate(Prevalence = if_else(vtdr == "VTDR", Count / TotalCount * 100, NA_real_),
         PrevalenceLabel = if_else(!is.na(Prevalence), sprintf("%.1f%%", Prevalence), ""))

df$event = factor(df$event, levels = c("mLOY", "mLOX", "Autosomal\n    mCAs", "Types of\n    Diabetes"))
df$Significance <- ""
df$Significance[df$event == "Types of\n    Diabetes"] <- "*"

ggplot(df, aes(x = event, y=Count, fill = Status, pattern = vtdr)) +
  geom_bar_pattern(width = 0.5, alpha=0.6,
                   pattern_angle = 45, pattern_density = 0.01, 
                   pattern_spacing = 0.01, pattern_key_scale_factor = 0.6,
                   position = "stack", stat = "identity") +
  scale_fill_manual(values = c("0" = "lightblue", "1" = "darkblue")) +
  scale_pattern_manual(values = c("No VTDR" = "none", "VTDR" = "stripe")) +
  labs(x = "", y = "", fill = "Category", pattern = "VTDR") +
  theme_minimal()  + 
  coord_flip() +
  geom_text(aes(label = PrevalenceLabel, vjust=5, hjust=-0.1), size=3,
                y=c(0, 0, 0, 450, 0, -50, 0, 200, 0, -50, 0, 200, 0, -50, 0, 200)) +
  geom_text(aes(label = Significance, y = 1650), vjust = 0.8, hjust = 0, 
            color = "black", fontface = "bold") +
  guides(pattern = guide_legend(override.aes = list(fill = "lightgrey"), order=2),
         fill = guide_legend(override.aes = list(pattern = "none"), order=1)) 
