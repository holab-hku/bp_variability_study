setwd("~/OneDrive - connect.hku.hk/Metagenomics/Hypertension")
setwd("C:/Users/GordonQ/OneDrive - connect.hku.hk/Metagenomics/Hypertension")
setwd("C:/Users/gordo/OneDrive - The University of Hong Kong - Connect/Metagenomics/Hypertension")

library(mixOmics)
library(ggplot2)
library(reshape2)
library(vegan)
library(ggpubr)
library(rstatix)
library(grid)
library(gridExtra)
library(gplots)
library(factoextra)

###Load Data & pre-processing###
norm_n_filter <- function(dat){
  dat <- sweep(dat,2,colSums(dat),`/`)
  filter.feat <- names(which(rowSums(dat*100 > 0.1) > .1*ncol(dat)))
  dat <- dat[rownames(dat) %in% filter.feat,]
  dat <- sweep(dat,2,colSums(dat),`/`)
  dat <- as.data.frame(t(dat))
  return(dat)
}
{
  #Meta data
  primary.dat <- read.csv("Data/Clinical/300 Cohort final_checked 2022_02_24/Primary_data.csv",header = T,row.names = 1)
  colnames(primary.dat) <- gsub("\\.\\..*", "", colnames(primary.dat))
  rownames(primary.dat) <- toupper(rownames(primary.dat))
  primary.dat["R0446","Hypertension.Staging"] <- 1
  primary.dat <- primary.dat[primary.dat$Anti.hypertensive.Agents==0,]
  colnames(primary.dat) <- gsub("day\\.time","Daytime",colnames(primary.dat))
  colnames(primary.dat) <- gsub("night\\.time","Nighttime",colnames(primary.dat))
  colnames(primary.dat) <- gsub("24h","24_hour",colnames(primary.dat))
  colnames(primary.dat) <- gsub("Coefficient\\.of\\.Variation","CoV",colnames(primary.dat))
  colnames(primary.dat) <- gsub("Standard\\.Deviation","SD",colnames(primary.dat))
  
  
  secondary.dat <- read.csv("Data/Clinical/300 Cohort final_checked 2022_02_24/Secondary_data.csv", header = T,row.names = 1)
  colnames(secondary.dat) <- gsub("\\.\\..*", "", colnames(secondary.dat))
  rownames(secondary.dat) <- toupper(rownames(secondary.dat))
  
  blood.stool.dat <- read.csv("Data/Clinical/GM cohort_raw data_no HTN meds_SCFA outiers removed_PV_05May22.csv", header = T,row.names = 1)
  colnames(blood.stool.dat) <- gsub("\\.\\..*", "", colnames(blood.stool.dat))
  rownames(blood.stool.dat) <- toupper(rownames(blood.stool.dat))
  blood.stool.dat <- blood.stool.dat[,71:105]
  blood.stool.dat[] <- lapply(blood.stool.dat, function(x) as.numeric(as.character(x)))
  immune_markers <- colnames(blood.stool.dat)[23:35]
  blood.stool.dat$Total.SCFAs <- blood.stool.dat$Total.SCFAs.1
  blood.stool.dat = subset(blood.stool.dat, select = -c(Total.SCFAs.1))
  
  #Dietary Data
  diet.dat <- read.csv("Data/Clinical/HEI_20220422/score.csv", header = T, row.names = 1)
  colnames(diet.dat) <- gsub("\\.\\..*", "", colnames(diet.dat))
  rownames(diet.dat) <- toupper(rownames(diet.dat))
    
  #Extra Dietary Data
  food.groups.dat <- read.csv("Data/Clinical/HEI_20220422/food_groups.csv", header = T, row.names = 1)
  colnames(food.groups.dat) <- gsub("\\.\\..*", "", colnames(food.groups.dat))
  rownames(food.groups.dat) <- toupper(rownames(food.groups.dat))
  
  #Diet group Data
  diet.groups.dat <- read.csv("Data/Clinical/HEI_20220422/diet_groups.csv", header = T, row.names = 1)
  colnames(diet.groups.dat) <- gsub("\\.\\..*", "", colnames(diet.groups.dat))
  rownames(diet.groups.dat) <- toupper(rownames(diet.groups.dat))
  diet.groups.dat$med_original_unhealthy <- ifelse(diet.groups.dat$med_original >= 4, 0, 1)
  diet.groups.dat$DASH_mullen_unhealthy <- ifelse(diet.groups.dat$DASH_mullen >= 4.5, 0, 1)
  
  #Dietary Selenium Data
  selenium <- read.csv("Data/Clinical/HEI_20220422/selenium.csv", header = T, row.names = 1)
  rownames(selenium) <- toupper(rownames(selenium))
  
  #Dietary Macro/Micro - nutrients data
  nutrients.dat <- read.csv("Data/Clinical/HEI_20220422/18012022_diet_241 cohort.csv", header = T, row.names = 1)
  nutri_feats <- c("Carbohydrate_available_g","Total_fat_g","Protein_g","Retinol_µg","Beta_carotene_µg","Vitamin_C_mg" ,"Vitamin_E_mg","Potassium_mg","Magnesium_mg","Calcium_mg",
                   "Phosphorus_mg","Iron_mg","Zinc_mg","Manganese_µg","Iodine_µg")
  nutrients.dat <- nutrients.dat[,nutri_feats]
  rownames(nutrients.dat) <- toupper(rownames(nutrients.dat))
  
  #Sleep data
  sleep_data <- read.csv("Data/Clinical/241_sleep_data.csv", row.names = 1)
  colnames(sleep_data) <- gsub("Latency.*","Latency",colnames(sleep_data))
  
  #MBPS data
  mbps_data <- read.csv("Data/Clinical/241_MBPS_data.csv", row.names = 1)
  
  #Merge all tables
  pat_order <- rownames(primary.dat)
  meta.data <- cbind(primary.dat,blood.stool.dat[pat_order,],secondary.dat[pat_order,],diet.dat[pat_order,],food.groups.dat[pat_order,],diet.groups.dat[pat_order,],selenium[pat_order,,drop=F], nutrients.dat[pat_order,],sleep_data[pat_order,],mbps_data[pat_order,])
  meta.data <- meta.data[,unlist(lapply(meta.data, is.numeric))]
  
  #Metagenomic data
  read.data <- read.table("Data/300/metaphlan3_res.tsv", sep = "\t", header = T, row.names = 1)
  read.data <- read.data[grepl("k__Bacteria",rownames(read.data)),]
  read.data <- read.data[,-1]
  colnames(read.data) <- substr(colnames(read.data), 1, 5); colnames(read.data) <- toupper(colnames(read.data))
  common <- intersect(colnames(read.data),rownames(primary.dat))
  read.data <- read.data[,common]
  
  #Species data
  species.data <- read.data[grepl("s__",rownames(read.data)),]
  rownames(species.data) <- gsub(".*s__","",rownames(species.data))
  species.data <- norm_n_filter(species.data)
  species.data <- species.data[rownames(meta.data),]
  #Genus data
  genus.data <- read.data[grepl("g__[^_]+$",rownames(read.data)),]
  rownames(genus.data) <- gsub(".*\\|g__","",rownames(genus.data))
  genus.data <- norm_n_filter(genus.data)
  genus.data <- genus.data[rownames(meta.data),]
  #Phylum data
  phylum.data <- read.data[grepl("p__[^_]+$",rownames(read.data)),]
  rownames(phylum.data) <- gsub(".*\\|p__","",rownames(phylum.data))
  phylum.data <- norm_n_filter(phylum.data)
  phylum.data <- phylum.data[rownames(meta.data),]
  ###Arcsin-sqrt transformation for compositional data###
  transformed.species.data <- asin(sqrt(species.data))
  transformed.genus.data <- asin(sqrt(genus.data))
  transformed.phylum.data <- asin(sqrt(phylum.data))
  
  #Pathway data
  path.data <- read.csv("Data/300/humann3_pathwayabundances.tsv", sep = "\t", header = T, row.names = 1)
  path.data <- read.csv("Data/300/humann3_ko_unstratified.tsv", sep = "\t", header = T, row.names = 1)
  colnames(path.data) <- substr(colnames(path.data), 1, 5); colnames(path.data) <- toupper(colnames(path.data))
  colnames(path.data) <- gsub("R(\\d*)\\.","R0\\1",colnames(path.data))
  common <- intersect(colnames(path.data),rownames(primary.dat))
  path.data <- path.data[,common]
  # Mapping ko ids
  ko_mapping <- read.csv("Data/full_KO_table.csv")
  grouped.data <- path.data
  grouped.data$group <- ko_mapping[match(rownames(grouped.data),ko_mapping$KO.ID.D),"KO.Desc.C"]
  grouped.data <- aggregate(. ~ group, data = grouped.data, FUN = sum)
  rownames(grouped.data) <- grouped.data$group; grouped.data <- grouped.data[,-1]
  grouped.data <- norm_n_filter(grouped.data)
  grouped.data <- grouped.data[rownames(grouped.data),]
  transformed.grouped.data <- asin(sqrt(grouped.data))
  colnames(transformed.grouped.data) <- make.names(colnames(transformed.grouped.data))
  
  meta.data$Sex <- ifelse(meta.data$Sex == 0, "Men", 'Women')
  male <- rownames(meta.data[meta.data$Sex=="Men",])
  female <- rownames(meta.data[meta.data$Sex=="Women",])
  meta.data$Hypertension.Staging <- as.factor(meta.data$Hypertension.Staging)
  meta.data$Smoking <- as.factor(meta.data$Smoking)
  meta.data$Menopause_code <- as.factor(meta.data$Menopause_code)
  meta.data$med_original_cat <- as.factor(meta.data$med_original_unhealthy)
  meta.data$DASH_mullen_cat <- as.factor(meta.data$DASH_mullen_unhealthy)
  meta.data$Dipping.classification <- as.factor(meta.data$Dipping.classification)
  colnames(meta.data) <- gsub("MSBPS","MBPS", colnames(meta.data))
  
  #meta.data$Dipping.classification.binary <- as.factor(ifelse(meta.data$Dipping.classification==2,0,1))
  
  merged.dat <- cbind(meta.data,transformed.species.data[rownames(meta.data),],transformed.grouped.data[rownames(meta.data),])
}

###Plotting functions###
stacked_bar <- function(dat,feat,max_feat=7,title=""){
  if(ncol(dat)>=max_feat){
    other_col <- colnames(dat)[which(!colnames(dat) %in% names(sort(colSums(dat),decreasing = T)[1:max_feat-1]))] 
    dat$Other <- rowSums(dat[,other_col])
    dat <- dat[,-which(colnames(dat) %in% other_col)]
  }
  dat.m <- reshape2::melt(as.matrix(dat))
  colnames(dat.m) <- c("Patient","Feat","Abund")
  dat.m$Lab <- as.factor(meta.data[match(dat.m$Patient,rownames(meta.data)),feat])
  dat.m$Sex <- as.factor(meta.data[match(dat.m$Patient,rownames(meta.data)),"Sex"])
  all.facet <- dat.m
  all.facet$Sex <- "All"
  dat.m <- rbind(dat.m,all.facet)
  dat.m <- dat.m[!is.na(dat.m$Lab),]
  print(ggplot(dat.m, aes(fill=Feat, y=Abund, x=Lab)) + 
          geom_bar(position="fill", stat="identity") + xlab(feat) + ylab("Relative Abundance") + scale_fill_brewer(palette="Dark2") + 
          ggtitle(title) + facet_wrap(~Sex))
}
my_boxplot <- function(dat,var,label, ylab="", title="", hide_ns = T){
  dat <- data.frame(y = dat[,var], lab = dat[,label], Sex = dat$Sex)
  all.facet <- dat
  all.facet$Sex <- "All"
  dat <- rbind(dat,all.facet)
  dat <- dat[!is.na(dat$lab),]
  stat.test <- try(dat %>% group_by(Sex) %>% pairwise_wilcox_test(y ~ lab, p.adjust.method = "BH") %>% add_y_position(), silent = T)
  if(is(stat.test, "try-error")){return(0)}
  if(sum(grepl("\\*",stat.test$p.adj.signif))>0){
    print(ggboxplot(dat,"lab", "y", color = "lab", add = "jitter") + stat_pvalue_manual(stat.test, label = "p.adj", hide.ns = hide_ns, label.size = 6) + 
          xlab(label) + ylab(ylab) + ggtitle(title) + facet_wrap(~Sex) + 
          theme(axis.title = element_text(size = 14, face="bold"), axis.text = element_text(size = 12, face="bold"), strip.text = element_text(size=12, face="bold"),
                legend.position = "none"))
  }
}
my_scatterplot <- function(dat,var,label, ylab="", title="", pos="bottom"){
  dat <- data.frame(y = dat[,var], lab = dat[,label], Sex = dat$Sex)
  all.facet <- dat
  all.facet$Sex <- "All"
  dat <- rbind(dat,all.facet)
  dat <- dat[!is.na(dat$lab),]
  print(ggscatter(dat,"lab", "y", add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), conf.int = T) + facet_wrap(~Sex) + 
          xlab(label) + ylab(ylab) + ggtitle(title) + stat_cor(method = "pearson", label.x.npc = "left", label.y.npc =pos, col = "red"))
}

lm_model <- function(dat,LHS,RHS){
  f <- formula(paste(LHS,"~",paste(RHS,collapse = '+')))
  return(lm(f,data = dat))
}
glm_model <- function(dat,LHS,RHS){
  f <- formula(paste(LHS,"~",paste(RHS,collapse = '+')))
  return(glm(f,data = dat, family = 'binomial'))
}
permanova_model <- function(dat,LHS,RHS){
  for(var in RHS){
    dat <- dat[!is.na(dat[,var]),]
  }
  LHS <- LHS[rownames(dat),]
  f <- formula(paste("LHS ~",paste(RHS,collapse = '+')))
  return(adonis2(f,data = dat, permutations = 1000, method = 'bray'))
}
anova_model <- function(dat,LHS,RHS){
  f <- formula(paste(LHS,"~",paste(RHS,collapse = '+')))
  return(aov(f,data = dat))
}

format_res <- function(dat,filt=F,filt_thres=0.05){
  dat$pval <- as.numeric(dat$pval)
  if(filt){
    filt.feat <- dat[dat$pval<=filt_thres,"feature"]
    dat <- dat[dat[,"feature"]%in%filt.feat,]
  }
  dat$sym <- ifelse(dat$pval>=0.05, "", ifelse(dat$pval<=0.01,ifelse(dat$pval<=0.001,"***", "**"), "*"))
  dat$log <- -log10(dat$pval)*as.numeric(dat$direction)
  dat$feature <- gsub("\\."," ",dat$feature)
  dat$feature <- gsub("X","",dat$feature)
  dat$feature <- gsub("_","-",dat$feature)
  dat$gender <- factor(dat$gender, levels = c("All","Women","Men"))
  return(dat)
}
format_res_pnova <- function(dat,filt=F,filt_thres=0.05){
  dat$pval <- as.numeric(dat$pval)
  if(filt){
    filt.feat <- dat[dat$pval<=filt_thres,"feature"]
    dat <- dat[dat[,"feature"]%in%filt.feat,]
  }
  dat$sym <- ifelse(dat$pval>=0.05, "", ifelse(dat$pval<=0.01,ifelse(dat$pval<=0.001,"***", "**"), "*"))
  dat$log <- -log10(dat$pval)
  dat$feature <- gsub("\\."," ",dat$feature)
  dat$feature <- gsub("X","",dat$feature)
  dat$feature <- gsub("_","-",dat$feature)
  dat$gender <- factor(dat$gender, levels = c("All","Women","Men"))
  return(dat)
}
### ENV Variables
SBP_cov <- c("X24_hour.SBP.CoV","Daytime.SBP.CoV","Nighttime.SBP.CoV")
SBP_surge <- c("Sleep.through.MBPS","Prewaking.MBPS")
DBP_cov <- c("X24_hour.DBP.CoV","Daytime.DBP.CoV","Nighttime.DBP.CoV")
SBP_sd <- c("X24_hour.SBP.SD", "Daytime.SBP.SD","Nighttime.SBP.SD")
DBP_sd <- c("X24_hour.DBP.SD", "Daytime.DBP.SD","Nighttime.DBP.SD")
secondary_categorical <- c("Dipping.classification")
acids_variables <- c(grep("Blood.SCFA",colnames(blood.stool.dat),value = T),"Total.SCFAs")
diet_variables_continuous <- c(names(which(colSums(is.na(diet.dat))!=nrow(diet.dat))),"med_fish","med_original","DASH_mullen","Selenium_ug")
diet_variables_categorical <- c("med_original_unhealthy","DASH_mullen_unhealthy")
food_group_variables <- names(which(colSums(is.na(food.groups.dat))!=nrow(food.groups.dat)))
vascular_dmg_variables <- c('Mean.baPWV','FMD','Mean.CCA.IMT','Mean.ABI')
macro_micro_nutrients <-  c("Energy_kcal","Carbohydrate_available_g","Total_fat_g","Protein_g","fiber_g","Retinol_µg","Beta_carotene_µg","Vitamin_C_mg" ,"Vitamin_E_mg","Potassium_mg","Magnesium_mg","Calcium_mg",
                            "Phosphorus_mg","Iron_mg","Zinc_mg","Manganese_µg","Iodine_µg","diet_sodium_mg")
stool_scfa_variables <- grep("^Stool",colnames(blood.stool.dat), value = T)

probiotic_species <- c("Bacteroides_dorei", "Bacteroides_stercoris",
"Bacteroides_plebeius",
"Bacteroides_coprocola",
"Bifidobacterium_pseudocatenulatum",
"Roseburia_intestinalis",
"Fusicatenibacter_saccharivorans",
"Ruthenibacterium_lactatiformans",
"Firmicutes_bacterium_CAG_83",
"Faecalibacterium_prausnitzii")

m2_covar <- c('Age','BMI')
m3_covar <- c(m2_covar,'Smoking','Sodium.Intake.based.on.spot.urine','Fasting.Glucose','Triglyceride','LDL.Cholesterol','Serum.HDL.Cholesterol', 'Latency')
m3_covar_sbp<- c(m2_covar,'Smoking','Sodium.Intake.based.on.spot.urine','Fasting.Glucose','Triglyceride','LDL.Cholesterol','Serum.HDL.Cholesterol', 'Latency','X24_hour.Mean.SBP')
m3_covar_dbp <- c(m2_covar,'Smoking','Sodium.Intake.based.on.spot.urine','Fasting.Glucose','Triglyceride','LDL.Cholesterol','Serum.HDL.Cholesterol', 'Latency','X24_hour.Mean.DBP')
m4_covar <- c(m3_covar,'Fatty.liver.by.CAP.score')
m4_covar_sbp <- c(m3_covar_sbp,'Fatty.liver.by.CAP.score')
m4_covar_dbp <- c(m3_covar_dbp,'Fatty.liver.by.CAP.score')

# F/B Ratio
{
merged.dat$f.b.ratio <- phylum.data$Firmicutes/phylum.data$Bacteroidetes
cor.matrix.full <-  data.frame(gender = character(), model = character(), feature = character(), pval = character(), direction = character(),coefficient = character(), se = character())
  
cor.matrix <- data.frame(gender = character(), model = character(), feature = character(), pval = character(), direction = character(),coefficient = character(), se = character())
for(feature in SBP_cov){
  m1 <- lm_model(merged.dat,'f.b.ratio',feature)
  m1.f <- lm_model(merged.dat[female,],'f.b.ratio',feature)
  m1.m <- lm_model(merged.dat[male,],'f.b.ratio',feature)
  m2 <- lm_model(merged.dat,'f.b.ratio',c(feature,m2_covar,'Sex'))
  m2.f <- lm_model(merged.dat[female,],'f.b.ratio',c(feature,m2_covar))
  m2.m <- lm_model(merged.dat[male,],'f.b.ratio',c(feature,m2_covar))
  m3 <- lm_model(merged.dat,'f.b.ratio',c(feature,m3_covar_sbp,'Sex'))
  m3.f <- lm_model(merged.dat[female,],'f.b.ratio',c(feature,m3_covar_sbp,'Menopause_code'))
  m3.m <- lm_model(merged.dat[male,],'f.b.ratio',c(feature,m3_covar_sbp))
  m4 <- lm_model(merged.dat,'f.b.ratio',c(feature,m4_covar_sbp,'Sex'))
  m4.f <- lm_model(merged.dat[female,],'f.b.ratio',c(feature,m4_covar_sbp,'Menopause_code'))
  m4.m <- lm_model(merged.dat[male,],'f.b.ratio',c(feature,m4_covar_sbp))
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m1",feature,summary(m1)$coefficients[2,4],sign(summary(m1)$coefficients[2,1]),summary(m1)$coefficients[2,1],summary(m1)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m1",feature,summary(m1.f)$coefficients[2,4],sign(summary(m1.f)$coefficients[2,1]),summary(m1.f)$coefficients[2,1],summary(m1.f)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m1",feature,summary(m1.m)$coefficients[2,4],sign(summary(m1.m)$coefficients[2,1]),summary(m1.m)$coefficients[2,1],summary(m1.m)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m2",feature,summary(m2)$coefficients[2,4],sign(summary(m2)$coefficients[2,1]),summary(m2)$coefficients[2,1],summary(m2)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m2",feature,summary(m2.f)$coefficients[2,4],sign(summary(m2.f)$coefficients[2,1]),summary(m2.f)$coefficients[2,1],summary(m2.f)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m2",feature,summary(m2.m)$coefficients[2,4],sign(summary(m2.m)$coefficients[2,1]),summary(m2.m)$coefficients[2,1],summary(m2.m)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3",feature,summary(m3)$coefficients[2,4],sign(summary(m3)$coefficients[2,1]),summary(m3)$coefficients[2,1],summary(m3)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3",feature,summary(m3.f)$coefficients[2,4],sign(summary(m3.f)$coefficients[2,1]),summary(m3.f)$coefficients[2,1],summary(m3.f)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3",feature,summary(m3.m)$coefficients[2,4],sign(summary(m3.m)$coefficients[2,1]),summary(m3.m)$coefficients[2,1],summary(m3.m)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m4",feature,summary(m4)$coefficients[2,4],sign(summary(m4)$coefficients[2,1]),summary(m4)$coefficients[2,1],summary(m4)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m4",feature,summary(m4.f)$coefficients[2,4],sign(summary(m4.f)$coefficients[2,1]),summary(m4.f)$coefficients[2,1],summary(m4.f)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m4",feature,summary(m4.m)$coefficients[2,4],sign(summary(m4.m)$coefficients[2,1]),summary(m4.m)$coefficients[2,1],summary(m4.m)$coefficients[2,2])
}
cor.matrix <- format_res(cor.matrix)
cor.matrix$feature <- factor(cor.matrix$feature, levels = c("Nighttime SBP CoV","Daytime SBP CoV","24-hour SBP CoV"))
pdf("Figures/300_bp_variability/fb_ratio_SBP_CoV.pdf", width = 8, height = 5)
print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
        coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender) + 
        theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"), 
              legend.text = element_text(size = 12, face="bold"), legend.title = element_text(size = 12, face="bold")) + ggtitle("F/B Ratio"))
dev.off()
cor.matrix.full <- rbind(cor.matrix.full,cor.matrix)

cor.matrix <- data.frame(gender = character(), model = character(), feature = character(), pval = character(), direction = character(),coefficient = character(), se = character())
for(feature in SBP_sd){
  m1 <- lm_model(merged.dat,'f.b.ratio',feature)
  m1.f <- lm_model(merged.dat[female,],'f.b.ratio',feature)
  m1.m <- lm_model(merged.dat[male,],'f.b.ratio',feature)
  m2 <- lm_model(merged.dat,'f.b.ratio',c(feature,m2_covar,'Sex'))
  m2.f <- lm_model(merged.dat[female,],'f.b.ratio',c(feature,m2_covar))
  m2.m <- lm_model(merged.dat[male,],'f.b.ratio',c(feature,m2_covar))
  m3 <- lm_model(merged.dat,'f.b.ratio',c(feature,m3_covar_sbp,'Sex'))
  m3.f <- lm_model(merged.dat[female,],'f.b.ratio',c(feature,m3_covar_sbp,'Menopause_code'))
  m3.m <- lm_model(merged.dat[male,],'f.b.ratio',c(feature,m3_covar_sbp))
  m4 <- lm_model(merged.dat,'f.b.ratio',c(feature,m4_covar_sbp,'Sex'))
  m4.f <- lm_model(merged.dat[female,],'f.b.ratio',c(feature,m4_covar_sbp,'Menopause_code'))
  m4.m <- lm_model(merged.dat[male,],'f.b.ratio',c(feature,m4_covar_sbp))
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m1",feature,summary(m1)$coefficients[2,4],sign(summary(m1)$coefficients[2,1]),summary(m1)$coefficients[2,1],summary(m1)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m1",feature,summary(m1.f)$coefficients[2,4],sign(summary(m1.f)$coefficients[2,1]),summary(m1.f)$coefficients[2,1],summary(m1.f)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m1",feature,summary(m1.m)$coefficients[2,4],sign(summary(m1.m)$coefficients[2,1]),summary(m1.m)$coefficients[2,1],summary(m1.m)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m2",feature,summary(m2)$coefficients[2,4],sign(summary(m2)$coefficients[2,1]),summary(m2)$coefficients[2,1],summary(m2)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m2",feature,summary(m2.f)$coefficients[2,4],sign(summary(m2.f)$coefficients[2,1]),summary(m2.f)$coefficients[2,1],summary(m2.f)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m2",feature,summary(m2.m)$coefficients[2,4],sign(summary(m2.m)$coefficients[2,1]),summary(m2.m)$coefficients[2,1],summary(m2.m)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3",feature,summary(m3)$coefficients[2,4],sign(summary(m3)$coefficients[2,1]),summary(m3)$coefficients[2,1],summary(m3)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3",feature,summary(m3.f)$coefficients[2,4],sign(summary(m3.f)$coefficients[2,1]),summary(m3.f)$coefficients[2,1],summary(m3.f)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3",feature,summary(m3.m)$coefficients[2,4],sign(summary(m3.m)$coefficients[2,1]),summary(m3.m)$coefficients[2,1],summary(m3.m)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m4",feature,summary(m4)$coefficients[2,4],sign(summary(m4)$coefficients[2,1]),summary(m4)$coefficients[2,1],summary(m4)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m4",feature,summary(m4.f)$coefficients[2,4],sign(summary(m4.f)$coefficients[2,1]),summary(m4.f)$coefficients[2,1],summary(m4.f)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m4",feature,summary(m4.m)$coefficients[2,4],sign(summary(m4.m)$coefficients[2,1]),summary(m4.m)$coefficients[2,1],summary(m4.m)$coefficients[2,2])
}
cor.matrix <- format_res(cor.matrix)
cor.matrix$feature <- factor(cor.matrix$feature, levels = c("Nighttime SBP SD","Daytime SBP SD","24-hour SBP SD"))
pdf("Figures/300_bp_variability/fb_ratio_SBP_SD.pdf", width = 8, height = 5)
print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
        coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender) + 
        theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"), 
              legend.text = element_text(size = 12, face="bold"), legend.title = element_text(size = 12, face="bold")) + ggtitle("F/B Ratio"))
dev.off()
cor.matrix.full <- rbind(cor.matrix.full,cor.matrix)

cor.matrix <- data.frame(gender = character(), model = character(), feature = character(), pval = character(), direction = character(),coefficient = character(), se = character())
for(feature in DBP_cov){
  m1 <- lm_model(merged.dat,'f.b.ratio',feature)
  m1.f <- lm_model(merged.dat[female,],'f.b.ratio',feature)
  m1.m <- lm_model(merged.dat[male,],'f.b.ratio',feature)
  m2 <- lm_model(merged.dat,'f.b.ratio',c(feature,m2_covar,'Sex'))
  m2.f <- lm_model(merged.dat[female,],'f.b.ratio',c(feature,m2_covar))
  m2.m <- lm_model(merged.dat[male,],'f.b.ratio',c(feature,m2_covar))
  m3 <- lm_model(merged.dat,'f.b.ratio',c(feature,m3_covar_dbp,'Sex'))
  m3.f <- lm_model(merged.dat[female,],'f.b.ratio',c(feature,m3_covar_dbp,'Menopause_code'))
  m3.m <- lm_model(merged.dat[male,],'f.b.ratio',c(feature,m3_covar_dbp))
  m4 <- lm_model(merged.dat,'f.b.ratio',c(feature,m4_covar_dbp,'Sex'))
  m4.f <- lm_model(merged.dat[female,],'f.b.ratio',c(feature,m4_covar_dbp,'Menopause_code'))
  m4.m <- lm_model(merged.dat[male,],'f.b.ratio',c(feature,m4_covar_dbp))
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m1",feature,summary(m1)$coefficients[2,4],sign(summary(m1)$coefficients[2,1]),summary(m1)$coefficients[2,1],summary(m1)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m1",feature,summary(m1.f)$coefficients[2,4],sign(summary(m1.f)$coefficients[2,1]),summary(m1.f)$coefficients[2,1],summary(m1.f)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m1",feature,summary(m1.m)$coefficients[2,4],sign(summary(m1.m)$coefficients[2,1]),summary(m1.m)$coefficients[2,1],summary(m1.m)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m2",feature,summary(m2)$coefficients[2,4],sign(summary(m2)$coefficients[2,1]),summary(m2)$coefficients[2,1],summary(m2)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m2",feature,summary(m2.f)$coefficients[2,4],sign(summary(m2.f)$coefficients[2,1]),summary(m2.f)$coefficients[2,1],summary(m2.f)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m2",feature,summary(m2.m)$coefficients[2,4],sign(summary(m2.m)$coefficients[2,1]),summary(m2.m)$coefficients[2,1],summary(m2.m)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3",feature,summary(m3)$coefficients[2,4],sign(summary(m3)$coefficients[2,1]),summary(m3)$coefficients[2,1],summary(m3)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3",feature,summary(m3.f)$coefficients[2,4],sign(summary(m3.f)$coefficients[2,1]),summary(m3.f)$coefficients[2,1],summary(m3.f)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3",feature,summary(m3.m)$coefficients[2,4],sign(summary(m3.m)$coefficients[2,1]),summary(m3.m)$coefficients[2,1],summary(m3.m)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m4",feature,summary(m4)$coefficients[2,4],sign(summary(m4)$coefficients[2,1]),summary(m4)$coefficients[2,1],summary(m4)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m4",feature,summary(m4.f)$coefficients[2,4],sign(summary(m4.f)$coefficients[2,1]),summary(m4.f)$coefficients[2,1],summary(m4.f)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m4",feature,summary(m4.m)$coefficients[2,4],sign(summary(m4.m)$coefficients[2,1]),summary(m4.m)$coefficients[2,1],summary(m4.m)$coefficients[2,2])
}
cor.matrix <- format_res(cor.matrix)
cor.matrix$feature <- factor(cor.matrix$feature, levels = c("Nighttime DBP CoV","Daytime DBP CoV","24-hour DBP CoV"))
pdf("Figures/300_bp_variability/fb_ratio_DBP_CoV.pdf", width = 8, height = 5)
print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
  coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender) + 
    theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"), 
          legend.text = element_text(size = 12, face="bold"), legend.title = element_text(size = 12, face="bold")) + ggtitle("F/B Ratio"))
dev.off()
cor.matrix.full <- rbind(cor.matrix.full,cor.matrix)

cor.matrix <- data.frame(gender = character(), model = character(), feature = character(), pval = character(), direction = character(),coefficient = character(), se = character())
for(feature in DBP_sd){
  m1 <- lm_model(merged.dat,'f.b.ratio',feature)
  m1.f <- lm_model(merged.dat[female,],'f.b.ratio',feature)
  m1.m <- lm_model(merged.dat[male,],'f.b.ratio',feature)
  m2 <- lm_model(merged.dat,'f.b.ratio',c(feature,m2_covar,'Sex'))
  m2.f <- lm_model(merged.dat[female,],'f.b.ratio',c(feature,m2_covar))
  m2.m <- lm_model(merged.dat[male,],'f.b.ratio',c(feature,m2_covar))
  m3 <- lm_model(merged.dat,'f.b.ratio',c(feature,m3_covar_dbp,'Sex'))
  m3.f <- lm_model(merged.dat[female,],'f.b.ratio',c(feature,m3_covar_dbp,'Menopause_code'))
  m3.m <- lm_model(merged.dat[male,],'f.b.ratio',c(feature,m3_covar_dbp))
  m4 <- lm_model(merged.dat,'f.b.ratio',c(feature,m4_covar_dbp,'Sex'))
  m4.f <- lm_model(merged.dat[female,],'f.b.ratio',c(feature,m4_covar_dbp,'Menopause_code'))
  m4.m <- lm_model(merged.dat[male,],'f.b.ratio',c(feature,m4_covar_dbp))
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m1",feature,summary(m1)$coefficients[2,4],sign(summary(m1)$coefficients[2,1]),summary(m1)$coefficients[2,1],summary(m1)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m1",feature,summary(m1.f)$coefficients[2,4],sign(summary(m1.f)$coefficients[2,1]),summary(m1.f)$coefficients[2,1],summary(m1.f)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m1",feature,summary(m1.m)$coefficients[2,4],sign(summary(m1.m)$coefficients[2,1]),summary(m1.m)$coefficients[2,1],summary(m1.m)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m2",feature,summary(m2)$coefficients[2,4],sign(summary(m2)$coefficients[2,1]),summary(m2)$coefficients[2,1],summary(m2)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m2",feature,summary(m2.f)$coefficients[2,4],sign(summary(m2.f)$coefficients[2,1]),summary(m2.f)$coefficients[2,1],summary(m2.f)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m2",feature,summary(m2.m)$coefficients[2,4],sign(summary(m2.m)$coefficients[2,1]),summary(m2.m)$coefficients[2,1],summary(m2.m)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3",feature,summary(m3)$coefficients[2,4],sign(summary(m3)$coefficients[2,1]),summary(m3)$coefficients[2,1],summary(m3)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3",feature,summary(m3.f)$coefficients[2,4],sign(summary(m3.f)$coefficients[2,1]),summary(m3.f)$coefficients[2,1],summary(m3.f)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3",feature,summary(m3.m)$coefficients[2,4],sign(summary(m3.m)$coefficients[2,1]),summary(m3.m)$coefficients[2,1],summary(m3.m)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m4",feature,summary(m4)$coefficients[2,4],sign(summary(m4)$coefficients[2,1]),summary(m4)$coefficients[2,1],summary(m4)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m4",feature,summary(m4.f)$coefficients[2,4],sign(summary(m4.f)$coefficients[2,1]),summary(m4.f)$coefficients[2,1],summary(m4.f)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m4",feature,summary(m4.m)$coefficients[2,4],sign(summary(m4.m)$coefficients[2,1]),summary(m4.m)$coefficients[2,1],summary(m4.m)$coefficients[2,2])
}
cor.matrix <- format_res(cor.matrix)
cor.matrix$feature <- factor(cor.matrix$feature, levels = c("Nighttime DBP SD","Daytime DBP SD","24-hour DBP SD"))
pdf("Figures/300_bp_variability/fb_ratio_DBP_SD.pdf", width = 8, height = 5)
print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
        coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender) + 
        theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"), 
              legend.text = element_text(size = 12, face="bold"), legend.title = element_text(size = 12, face="bold")) + ggtitle("F/B Ratio"))
dev.off()
cor.matrix.full <- rbind(cor.matrix.full,cor.matrix)

cor.matrix <- data.frame(gender = character(), model = character(), feature = character(), pval = character(), direction = character(),coefficient = character(), se = character())
for(feature in SBP_surge){
  m1 <- lm_model(merged.dat,'f.b.ratio',feature)
  m1.f <- lm_model(merged.dat[female,],'f.b.ratio',feature)
  m1.m <- lm_model(merged.dat[male,],'f.b.ratio',feature)
  m2 <- lm_model(merged.dat,'f.b.ratio',c(feature,m2_covar,'Sex'))
  m2.f <- lm_model(merged.dat[female,],'f.b.ratio',c(feature,m2_covar))
  m2.m <- lm_model(merged.dat[male,],'f.b.ratio',c(feature,m2_covar))
  m3 <- lm_model(merged.dat,'f.b.ratio',c(feature,m3_covar_sbp,'Sex'))
  m3.f <- lm_model(merged.dat[female,],'f.b.ratio',c(feature,m3_covar_sbp,'Menopause_code'))
  m3.m <- lm_model(merged.dat[male,],'f.b.ratio',c(feature,m3_covar_sbp))
  m4 <- lm_model(merged.dat,'f.b.ratio',c(feature,m4_covar_sbp,'Sex'))
  m4.f <- lm_model(merged.dat[female,],'f.b.ratio',c(feature,m4_covar_sbp,'Menopause_code'))
  m4.m <- lm_model(merged.dat[male,],'f.b.ratio',c(feature,m4_covar_sbp))
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m1",feature,summary(m1)$coefficients[2,4],sign(summary(m1)$coefficients[2,1]),summary(m1)$coefficients[2,1],summary(m1)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m1",feature,summary(m1.f)$coefficients[2,4],sign(summary(m1.f)$coefficients[2,1]),summary(m1.f)$coefficients[2,1],summary(m1.f)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m1",feature,summary(m1.m)$coefficients[2,4],sign(summary(m1.m)$coefficients[2,1]),summary(m1.m)$coefficients[2,1],summary(m1.m)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m2",feature,summary(m2)$coefficients[2,4],sign(summary(m2)$coefficients[2,1]),summary(m2)$coefficients[2,1],summary(m2)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m2",feature,summary(m2.f)$coefficients[2,4],sign(summary(m2.f)$coefficients[2,1]),summary(m2.f)$coefficients[2,1],summary(m2.f)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m2",feature,summary(m2.m)$coefficients[2,4],sign(summary(m2.m)$coefficients[2,1]),summary(m2.m)$coefficients[2,1],summary(m2.m)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3",feature,summary(m3)$coefficients[2,4],sign(summary(m3)$coefficients[2,1]),summary(m3)$coefficients[2,1],summary(m3)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3",feature,summary(m3.f)$coefficients[2,4],sign(summary(m3.f)$coefficients[2,1]),summary(m3.f)$coefficients[2,1],summary(m3.f)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3",feature,summary(m3.m)$coefficients[2,4],sign(summary(m3.m)$coefficients[2,1]),summary(m3.m)$coefficients[2,1],summary(m3.m)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m4",feature,summary(m4)$coefficients[2,4],sign(summary(m4)$coefficients[2,1]),summary(m4)$coefficients[2,1],summary(m4)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m4",feature,summary(m4.f)$coefficients[2,4],sign(summary(m4.f)$coefficients[2,1]),summary(m4.f)$coefficients[2,1],summary(m4.f)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m4",feature,summary(m4.m)$coefficients[2,4],sign(summary(m4.m)$coefficients[2,1]),summary(m4.m)$coefficients[2,1],summary(m4.m)$coefficients[2,2])
}
cor.matrix <- format_res(cor.matrix)
pdf("Figures/300_bp_variability/fb_ratio_SBP_surge.pdf", width = 8, height = 5)
print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(name = "-log10", midpoint=0, low="blue", mid="white",high="red", space ="Lab") +
        coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender) + 
        theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"), 
              legend.text = element_text(size = 12, face="bold"), legend.title = element_text(size = 12, face="bold")) + ggtitle("F/B Ratio"))
dev.off()
cor.matrix.full <- rbind(cor.matrix.full,cor.matrix)
write.csv(cor.matrix.full[,c('gender','model','feature','pval','coefficient','se')], "Figures/300_bp_variability/fb_ratio_vs_BP_variability.csv",row.names = F)
}

# Alpha Diversity analysis
{
merged.dat$diversity <- vegan::diversity(transformed.species.data, index = 'shannon')
cor.matrix.full <-  data.frame(gender = character(), model = character(), feature = character(), pval = character(), direction = character(),coefficient = character(), se = character())
  
cor.matrix <-  data.frame(gender = character(), model = character(), feature = character(), pval = character(), direction = character(),coefficient = character(), se = character())
for(feature in SBP_cov){
  m1 <- lm_model(merged.dat,'diversity',feature)
  m1.f <- lm_model(merged.dat[female,],'diversity',feature)
  m1.m <- lm_model(merged.dat[male,],'diversity',feature)
  m2 <- lm_model(merged.dat,'diversity',c(feature,m2_covar,'Sex'))
  m2.f <- lm_model(merged.dat[female,],'diversity',c(feature,m2_covar))
  m2.m <- lm_model(merged.dat[male,],'diversity',c(feature,m2_covar))
  m3 <- lm_model(merged.dat,'diversity',c(feature,m3_covar_sbp,'Sex'))
  m3.f <- lm_model(merged.dat[female,],'diversity',c(feature,m3_covar_sbp,'Menopause_code'))
  m3.m <- lm_model(merged.dat[male,],'diversity',c(feature,m3_covar_sbp))
  m4 <- lm_model(merged.dat,'diversity',c(feature,m4_covar_sbp,'Sex'))
  m4.f <- lm_model(merged.dat[female,],'diversity',c(feature,m4_covar_sbp,'Menopause_code'))
  m4.m <- lm_model(merged.dat[male,],'diversity',c(feature,m4_covar_sbp))
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m1",feature,summary(m1)$coefficients[2,4],sign(summary(m1)$coefficients[2,1]),summary(m1)$coefficients[2,1],summary(m1)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m1",feature,summary(m1.f)$coefficients[2,4],sign(summary(m1.f)$coefficients[2,1]),summary(m1.f)$coefficients[2,1],summary(m1.f)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m1",feature,summary(m1.m)$coefficients[2,4],sign(summary(m1.m)$coefficients[2,1]),summary(m1.m)$coefficients[2,1],summary(m1.m)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m2",feature,summary(m2)$coefficients[2,4],sign(summary(m2)$coefficients[2,1]),summary(m2)$coefficients[2,1],summary(m2)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m2",feature,summary(m2.f)$coefficients[2,4],sign(summary(m2.f)$coefficients[2,1]),summary(m2.f)$coefficients[2,1],summary(m2.f)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m2",feature,summary(m2.m)$coefficients[2,4],sign(summary(m2.m)$coefficients[2,1]),summary(m2.m)$coefficients[2,1],summary(m2.m)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3",feature,summary(m3)$coefficients[2,4],sign(summary(m3)$coefficients[2,1]),summary(m3)$coefficients[2,1],summary(m3)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3",feature,summary(m3.f)$coefficients[2,4],sign(summary(m3.f)$coefficients[2,1]),summary(m3.f)$coefficients[2,1],summary(m3.f)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3",feature,summary(m3.m)$coefficients[2,4],sign(summary(m3.m)$coefficients[2,1]),summary(m3.m)$coefficients[2,1],summary(m3.m)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m4",feature,summary(m4)$coefficients[2,4],sign(summary(m4)$coefficients[2,1]),summary(m4)$coefficients[2,1],summary(m4)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m4",feature,summary(m4.f)$coefficients[2,4],sign(summary(m4.f)$coefficients[2,1]),summary(m4.f)$coefficients[2,1],summary(m4.f)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m4",feature,summary(m4.m)$coefficients[2,4],sign(summary(m4.m)$coefficients[2,1]),summary(m4.m)$coefficients[2,1],summary(m4.m)$coefficients[2,2])
}
cor.matrix <- format_res(cor.matrix)
cor.matrix$feature <- factor(cor.matrix$feature, levels = c("Nighttime SBP CoV","Daytime SBP CoV","24-hour SBP CoV"))
pdf("Figures/300_bp_variability/shannon_diversity_SBP_CoV.pdf", width = 8, height = 5)
print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(name = "-log10", midpoint=0, low="blue", mid="white",high="red", space ="Lab", limits=c(-2.2,1)) +
        coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender) + 
        theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"), 
              legend.text = element_text(size = 12, face="bold"), legend.title = element_text(size = 12, face="bold")) + ggtitle("Shannon Diversity"))
dev.off()
cor.matrix.full <- rbind(cor.matrix.full,cor.matrix)

cor.matrix <-  data.frame(gender = character(), model = character(), feature = character(), pval = character(), direction = character(),coefficient = character(), se = character())
for(feature in SBP_sd){
  m1 <- lm_model(merged.dat,'diversity',feature)
  m1.f <- lm_model(merged.dat[female,],'diversity',feature)
  m1.m <- lm_model(merged.dat[male,],'diversity',feature)
  m2 <- lm_model(merged.dat,'diversity',c(feature,m2_covar,'Sex'))
  m2.f <- lm_model(merged.dat[female,],'diversity',c(feature,m2_covar))
  m2.m <- lm_model(merged.dat[male,],'diversity',c(feature,m2_covar))
  m3 <- lm_model(merged.dat,'diversity',c(feature,m3_covar_sbp,'Sex'))
  m3.f <- lm_model(merged.dat[female,],'diversity',c(feature,m3_covar_sbp,'Menopause_code'))
  m3.m <- lm_model(merged.dat[male,],'diversity',c(feature,m3_covar_sbp))
  m4 <- lm_model(merged.dat,'diversity',c(feature,m4_covar_sbp,'Sex'))
  m4.f <- lm_model(merged.dat[female,],'diversity',c(feature,m4_covar_sbp,'Menopause_code'))
  m4.m <- lm_model(merged.dat[male,],'diversity',c(feature,m4_covar_sbp))
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m1",feature,summary(m1)$coefficients[2,4],sign(summary(m1)$coefficients[2,1]),summary(m1)$coefficients[2,1],summary(m1)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m1",feature,summary(m1.f)$coefficients[2,4],sign(summary(m1.f)$coefficients[2,1]),summary(m1.f)$coefficients[2,1],summary(m1.f)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m1",feature,summary(m1.m)$coefficients[2,4],sign(summary(m1.m)$coefficients[2,1]),summary(m1.m)$coefficients[2,1],summary(m1.m)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m2",feature,summary(m2)$coefficients[2,4],sign(summary(m2)$coefficients[2,1]),summary(m2)$coefficients[2,1],summary(m2)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m2",feature,summary(m2.f)$coefficients[2,4],sign(summary(m2.f)$coefficients[2,1]),summary(m2.f)$coefficients[2,1],summary(m2.f)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m2",feature,summary(m2.m)$coefficients[2,4],sign(summary(m2.m)$coefficients[2,1]),summary(m2.m)$coefficients[2,1],summary(m2.m)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3",feature,summary(m3)$coefficients[2,4],sign(summary(m3)$coefficients[2,1]),summary(m3)$coefficients[2,1],summary(m3)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3",feature,summary(m3.f)$coefficients[2,4],sign(summary(m3.f)$coefficients[2,1]),summary(m3.f)$coefficients[2,1],summary(m3.f)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3",feature,summary(m3.m)$coefficients[2,4],sign(summary(m3.m)$coefficients[2,1]),summary(m3.m)$coefficients[2,1],summary(m3.m)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m4",feature,summary(m4)$coefficients[2,4],sign(summary(m4)$coefficients[2,1]),summary(m4)$coefficients[2,1],summary(m4)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m4",feature,summary(m4.f)$coefficients[2,4],sign(summary(m4.f)$coefficients[2,1]),summary(m4.f)$coefficients[2,1],summary(m4.f)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m4",feature,summary(m4.m)$coefficients[2,4],sign(summary(m4.m)$coefficients[2,1]),summary(m4.m)$coefficients[2,1],summary(m4.m)$coefficients[2,2])
}
cor.matrix <- format_res(cor.matrix)
cor.matrix$feature <- factor(cor.matrix$feature, levels = c("Nighttime SBP SD","Daytime SBP SD","24-hour SBP SD"))
pdf("Figures/300_bp_variability/shannon_diversity_SBP_SD.pdf", width = 8, height = 5)
print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(name = "-log10", midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
        coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender) + 
        theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"), 
              legend.text = element_text(size = 12, face="bold"), legend.title = element_text(size = 12, face="bold")) + ggtitle("Shannon Diversity"))
dev.off()
cor.matrix.full <- rbind(cor.matrix.full,cor.matrix)

cor.matrix <-  data.frame(gender = character(), model = character(), feature = character(), pval = character(), direction = character(),coefficient = character(), se = character())
for(feature in DBP_cov){
  m1 <- lm_model(merged.dat,'diversity',feature)
  m1.f <- lm_model(merged.dat[female,],'diversity',feature)
  m1.m <- lm_model(merged.dat[male,],'diversity',feature)
  m2 <- lm_model(merged.dat,'diversity',c(feature,m2_covar,'Sex'))
  m2.f <- lm_model(merged.dat[female,],'diversity',c(feature,m2_covar))
  m2.m <- lm_model(merged.dat[male,],'diversity',c(feature,m2_covar))
  m3 <- lm_model(merged.dat,'diversity',c(feature,m3_covar_dbp,'Sex'))
  m3.f <- lm_model(merged.dat[female,],'diversity',c(feature,m3_covar_dbp,'Menopause_code'))
  m3.m <- lm_model(merged.dat[male,],'diversity',c(feature,m3_covar_dbp))
  m4 <- lm_model(merged.dat,'diversity',c(feature,m4_covar_dbp,'Sex'))
  m4.f <- lm_model(merged.dat[female,],'diversity',c(feature,m4_covar_dbp,'Menopause_code'))
  m4.m <- lm_model(merged.dat[male,],'diversity',c(feature,m4_covar_dbp))
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m1",feature,summary(m1)$coefficients[2,4],sign(summary(m1)$coefficients[2,1]),summary(m1)$coefficients[2,1],summary(m1)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m1",feature,summary(m1.f)$coefficients[2,4],sign(summary(m1.f)$coefficients[2,1]),summary(m1.f)$coefficients[2,1],summary(m1.f)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m1",feature,summary(m1.m)$coefficients[2,4],sign(summary(m1.m)$coefficients[2,1]),summary(m1.m)$coefficients[2,1],summary(m1.m)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m2",feature,summary(m2)$coefficients[2,4],sign(summary(m2)$coefficients[2,1]),summary(m2)$coefficients[2,1],summary(m2)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m2",feature,summary(m2.f)$coefficients[2,4],sign(summary(m2.f)$coefficients[2,1]),summary(m2.f)$coefficients[2,1],summary(m2.f)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m2",feature,summary(m2.m)$coefficients[2,4],sign(summary(m2.m)$coefficients[2,1]),summary(m2.m)$coefficients[2,1],summary(m2.m)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3",feature,summary(m3)$coefficients[2,4],sign(summary(m3)$coefficients[2,1]),summary(m3)$coefficients[2,1],summary(m3)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3",feature,summary(m3.f)$coefficients[2,4],sign(summary(m3.f)$coefficients[2,1]),summary(m3.f)$coefficients[2,1],summary(m3.f)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3",feature,summary(m3.m)$coefficients[2,4],sign(summary(m3.m)$coefficients[2,1]),summary(m3.m)$coefficients[2,1],summary(m3.m)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m4",feature,summary(m4)$coefficients[2,4],sign(summary(m4)$coefficients[2,1]),summary(m4)$coefficients[2,1],summary(m4)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m4",feature,summary(m4.f)$coefficients[2,4],sign(summary(m4.f)$coefficients[2,1]),summary(m4.f)$coefficients[2,1],summary(m4.f)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m4",feature,summary(m4.m)$coefficients[2,4],sign(summary(m4.m)$coefficients[2,1]),summary(m4.m)$coefficients[2,1],summary(m4.m)$coefficients[2,2])
}
cor.matrix <- format_res(cor.matrix)
cor.matrix$feature <- factor(cor.matrix$feature, levels = c("Nighttime DBP CoV","Daytime DBP CoV","24-hour DBP CoV"))
pdf("Figures/300_bp_variability/shannon_diversity_DBP_CoV.pdf", width = 8, height = 5)
print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(name = "-log10", midpoint=0, low="blue", mid="white",high="red", space ="Lab", limits=c(-2.2,1)) +
        coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender) + 
        theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"), 
              legend.text = element_text(size = 12, face="bold"), legend.title = element_text(size = 12, face="bold")) + ggtitle("Shannon Diversity"))
dev.off()
cor.matrix.full <- rbind(cor.matrix.full,cor.matrix)

cor.matrix <-  data.frame(gender = character(), model = character(), feature = character(), pval = character(), direction = character(),coefficient = character(), se = character())
for(feature in DBP_sd){
  m1 <- lm_model(merged.dat,'diversity',feature)
  m1.f <- lm_model(merged.dat[female,],'diversity',feature)
  m1.m <- lm_model(merged.dat[male,],'diversity',feature)
  m2 <- lm_model(merged.dat,'diversity',c(feature,m2_covar,'Sex'))
  m2.f <- lm_model(merged.dat[female,],'diversity',c(feature,m2_covar))
  m2.m <- lm_model(merged.dat[male,],'diversity',c(feature,m2_covar))
  m3 <- lm_model(merged.dat,'diversity',c(feature,m3_covar_dbp,'Sex'))
  m3.f <- lm_model(merged.dat[female,],'diversity',c(feature,m3_covar_dbp,'Menopause_code'))
  m3.m <- lm_model(merged.dat[male,],'diversity',c(feature,m3_covar_dbp))
  m4 <- lm_model(merged.dat,'diversity',c(feature,m4_covar_dbp,'Sex'))
  m4.f <- lm_model(merged.dat[female,],'diversity',c(feature,m4_covar_dbp,'Menopause_code'))
  m4.m <- lm_model(merged.dat[male,],'diversity',c(feature,m4_covar_dbp))
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m1",feature,summary(m1)$coefficients[2,4],sign(summary(m1)$coefficients[2,1]),summary(m1)$coefficients[2,1],summary(m1)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m1",feature,summary(m1.f)$coefficients[2,4],sign(summary(m1.f)$coefficients[2,1]),summary(m1.f)$coefficients[2,1],summary(m1.f)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m1",feature,summary(m1.m)$coefficients[2,4],sign(summary(m1.m)$coefficients[2,1]),summary(m1.m)$coefficients[2,1],summary(m1.m)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m2",feature,summary(m2)$coefficients[2,4],sign(summary(m2)$coefficients[2,1]),summary(m2)$coefficients[2,1],summary(m2)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m2",feature,summary(m2.f)$coefficients[2,4],sign(summary(m2.f)$coefficients[2,1]),summary(m2.f)$coefficients[2,1],summary(m2.f)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m2",feature,summary(m2.m)$coefficients[2,4],sign(summary(m2.m)$coefficients[2,1]),summary(m2.m)$coefficients[2,1],summary(m2.m)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3",feature,summary(m3)$coefficients[2,4],sign(summary(m3)$coefficients[2,1]),summary(m3)$coefficients[2,1],summary(m3)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3",feature,summary(m3.f)$coefficients[2,4],sign(summary(m3.f)$coefficients[2,1]),summary(m3.f)$coefficients[2,1],summary(m3.f)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3",feature,summary(m3.m)$coefficients[2,4],sign(summary(m3.m)$coefficients[2,1]),summary(m3.m)$coefficients[2,1],summary(m3.m)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m4",feature,summary(m4)$coefficients[2,4],sign(summary(m4)$coefficients[2,1]),summary(m4)$coefficients[2,1],summary(m4)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m4",feature,summary(m4.f)$coefficients[2,4],sign(summary(m4.f)$coefficients[2,1]),summary(m4.f)$coefficients[2,1],summary(m4.f)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m4",feature,summary(m4.m)$coefficients[2,4],sign(summary(m4.m)$coefficients[2,1]),summary(m4.m)$coefficients[2,1],summary(m4.m)$coefficients[2,2])
}
cor.matrix <- format_res(cor.matrix)
cor.matrix$feature <- factor(cor.matrix$feature, levels = c("Nighttime DBP SD","Daytime DBP SD","24-hour DBP SD"))
pdf("Figures/300_bp_variability/shannon_diversity_DBP_SD.pdf", width = 8, height = 5)
print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(name = "-log10", midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
        coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender) + 
        theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"), 
              legend.text = element_text(size = 12, face="bold"), legend.title = element_text(size = 12, face="bold")) + ggtitle("Shannon Diversity"))
dev.off()
cor.matrix.full <- rbind(cor.matrix.full,cor.matrix)

cor.matrix <-  data.frame(gender = character(), model = character(), feature = character(), pval = character(), direction = character(),coefficient = character(), se = character())
for(feature in SBP_surge){
  m1 <- lm_model(merged.dat,'diversity',feature)
  m1.f <- lm_model(merged.dat[female,],'diversity',feature)
  m1.m <- lm_model(merged.dat[male,],'diversity',feature)
  m2 <- lm_model(merged.dat,'diversity',c(feature,m2_covar,'Sex'))
  m2.f <- lm_model(merged.dat[female,],'diversity',c(feature,m2_covar))
  m2.m <- lm_model(merged.dat[male,],'diversity',c(feature,m2_covar))
  m3 <- lm_model(merged.dat,'diversity',c(feature,m3_covar_sbp,'Sex'))
  m3.f <- lm_model(merged.dat[female,],'diversity',c(feature,m3_covar_sbp,'Menopause_code'))
  m3.m <- lm_model(merged.dat[male,],'diversity',c(feature,m3_covar_sbp))
  m4 <- lm_model(merged.dat,'diversity',c(feature,m4_covar_sbp,'Sex'))
  m4.f <- lm_model(merged.dat[female,],'diversity',c(feature,m4_covar_sbp,'Menopause_code'))
  m4.m <- lm_model(merged.dat[male,],'diversity',c(feature,m4_covar_sbp))
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m1",feature,summary(m1)$coefficients[2,4],sign(summary(m1)$coefficients[2,1]),summary(m1)$coefficients[2,1],summary(m1)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m1",feature,summary(m1.f)$coefficients[2,4],sign(summary(m1.f)$coefficients[2,1]),summary(m1.f)$coefficients[2,1],summary(m1.f)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m1",feature,summary(m1.m)$coefficients[2,4],sign(summary(m1.m)$coefficients[2,1]),summary(m1.m)$coefficients[2,1],summary(m1.m)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m2",feature,summary(m2)$coefficients[2,4],sign(summary(m2)$coefficients[2,1]),summary(m2)$coefficients[2,1],summary(m2)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m2",feature,summary(m2.f)$coefficients[2,4],sign(summary(m2.f)$coefficients[2,1]),summary(m2.f)$coefficients[2,1],summary(m2.f)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m2",feature,summary(m2.m)$coefficients[2,4],sign(summary(m2.m)$coefficients[2,1]),summary(m2.m)$coefficients[2,1],summary(m2.m)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3",feature,summary(m3)$coefficients[2,4],sign(summary(m3)$coefficients[2,1]),summary(m3)$coefficients[2,1],summary(m3)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3",feature,summary(m3.f)$coefficients[2,4],sign(summary(m3.f)$coefficients[2,1]),summary(m3.f)$coefficients[2,1],summary(m3.f)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3",feature,summary(m3.m)$coefficients[2,4],sign(summary(m3.m)$coefficients[2,1]),summary(m3.m)$coefficients[2,1],summary(m3.m)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m4",feature,summary(m4)$coefficients[2,4],sign(summary(m4)$coefficients[2,1]),summary(m4)$coefficients[2,1],summary(m4)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m4",feature,summary(m4.f)$coefficients[2,4],sign(summary(m4.f)$coefficients[2,1]),summary(m4.f)$coefficients[2,1],summary(m4.f)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m4",feature,summary(m4.m)$coefficients[2,4],sign(summary(m4.m)$coefficients[2,1]),summary(m4.m)$coefficients[2,1],summary(m4.m)$coefficients[2,2])
}
cor.matrix <- format_res(cor.matrix)
pdf("Figures/300_bp_variability/shannon_diversity_SBP_surge.pdf", width = 8, height = 5)
print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(name = "-log10", midpoint=0, low="blue", mid="white",high="red", space ="Lab", limits=c(-2.2,1)) +
        coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender) + 
        theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"), 
              legend.text = element_text(size = 12, face="bold"), legend.title = element_text(size = 12, face="bold")) + ggtitle("Shannon Diversity"))
dev.off()
cor.matrix.full <- rbind(cor.matrix.full,cor.matrix)
write.csv(cor.matrix.full[,c('gender','model','feature','pval','coefficient','se')], "Figures/300_bp_variability/shannon_diversity_vs_BP_variability.csv",row.names = F)
}

# Permanova - Beta diversity
{
cor.matrix.full <-  data.frame(gender = character(), model = character(), feature = character(), pval = character(), Fscore = character(), R2 = character(), SumOfSqs = character(), sym = character(), log = character())
  
cor.matrix <-  data.frame(gender = character(), model = character(), feature = character(), pval = character(), Fscore = character(), R2 = character(), SumOfSqs = character())
for(feature in SBP_cov){
  m1 <- permanova_model(merged.dat,species.data,feature)
  m1.f <- permanova_model(merged.dat[female,],species.data,feature)
  m1.m <- permanova_model(merged.dat[male,],species.data,feature)
  m2 <- permanova_model(merged.dat,species.data,c(feature,m2_covar,'Sex'))
  m2.f <- permanova_model(merged.dat[female,],species.data,c(feature,m2_covar))
  m2.m <- permanova_model(merged.dat[male,],species.data,c(feature,m2_covar))
  m3 <- permanova_model(merged.dat,species.data,c(feature,m3_covar_sbp,'Sex'))
  m3.f <- permanova_model(merged.dat[female,],species.data,c(feature,m3_covar_sbp,'Menopause_code'))
  m3.m <- permanova_model(merged.dat[male,],species.data,c(feature,m3_covar_sbp))
  m4 <- permanova_model(merged.dat,species.data,c(feature,m4_covar_sbp,'Sex'))
  m4.f <- permanova_model(merged.dat[female,],species.data,c(feature,m4_covar_sbp,'Menopause_code'))
  m4.m <- permanova_model(merged.dat[male,],species.data,c(feature,m4_covar_sbp))
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m1",feature,m1$`Pr(>F)`[1],m1$`F`[1],m1$`R2`[1],m1$`SumOfSqs`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m1",feature,m1.f$`Pr(>F)`[1],m1.f$`F`[1],m1.f$`R2`[1],m1.f$`SumOfSqs`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m1",feature,m1.m$`Pr(>F)`[1],m1.m$`F`[1],m1.m$`R2`[1],m1.m$`SumOfSqs`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m2",feature,m2$`Pr(>F)`[1],m2$`F`[1],m2$`R2`[1],m2$`SumOfSqs`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m2",feature,m2.f$`Pr(>F)`[1],m2.f$`F`[1],m2.f$`R2`[1],m2.f$`SumOfSqs`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m2",feature,m2.m$`Pr(>F)`[1],m2.m$`F`[1],m2.m$`R2`[1],m2.m$`SumOfSqs`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3",feature,m3$`Pr(>F)`[1],m3$`F`[1],m3$`R2`[1],m3$`SumOfSqs`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3",feature,m3.f$`Pr(>F)`[1],m3.f$`F`[1],m3.f$`R2`[1],m3.f$`SumOfSqs`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3",feature,m3.m$`Pr(>F)`[1],m3.m$`F`[1],m3.m$`R2`[1],m3.m$`SumOfSqs`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m4",feature,m4$`Pr(>F)`[1],m4$`F`[1],m4$`R2`[1],m4$`SumOfSqs`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m4",feature,m4.f$`Pr(>F)`[1],m4.f$`F`[1],m4.f$`R2`[1],m4.f$`SumOfSqs`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m4",feature,m4.m$`Pr(>F)`[1],m4.m$`F`[1],m4.m$`R2`[1],m4.m$`SumOfSqs`[1])
}
cor.matrix <- format_res_pnova(cor.matrix)
cor.matrix$feature <- factor(cor.matrix$feature, levels = c("Nighttime SBP CoV","Daytime SBP CoV","24-hour SBP CoV"))
pdf("Figures/300_bp_variability/permanova_SBP_CoV.pdf", width = 8, height = 5)
print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", name ="-log10" ) +
        coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender) +
        theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"),
              legend.title = element_text(size = 12, face="bold"), legend.text =  element_text(size = 12, face="bold")))
dev.off()
cor.matrix.full <- rbind(cor.matrix.full,cor.matrix)

cor.matrix <-  data.frame(gender = character(), model = character(), feature = character(), pval = character(), Fscore = character(), R2 = character(), SumOfSqs = character())
for(feature in SBP_sd){
  m1 <- permanova_model(merged.dat,species.data,feature)
  m1.f <- permanova_model(merged.dat[female,],species.data,feature)
  m1.m <- permanova_model(merged.dat[male,],species.data,feature)
  m2 <- permanova_model(merged.dat,species.data,c(feature,m2_covar,'Sex'))
  m2.f <- permanova_model(merged.dat[female,],species.data,c(feature,m2_covar))
  m2.m <- permanova_model(merged.dat[male,],species.data,c(feature,m2_covar))
  m3 <- permanova_model(merged.dat,species.data,c(feature,m3_covar_sbp,'Sex'))
  m3.f <- permanova_model(merged.dat[female,],species.data,c(feature,m3_covar_sbp,'Menopause_code'))
  m3.m <- permanova_model(merged.dat[male,],species.data,c(feature,m3_covar_sbp))
  m4 <- permanova_model(merged.dat,species.data,c(feature,m4_covar_sbp,'Sex'))
  m4.f <- permanova_model(merged.dat[female,],species.data,c(feature,m4_covar_sbp,'Menopause_code'))
  m4.m <- permanova_model(merged.dat[male,],species.data,c(feature,m4_covar_sbp))
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m1",feature,m1$`Pr(>F)`[1],m1$`F`[1],m1$`R2`[1],m1$`SumOfSqs`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m1",feature,m1.f$`Pr(>F)`[1],m1.f$`F`[1],m1.f$`R2`[1],m1.f$`SumOfSqs`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m1",feature,m1.m$`Pr(>F)`[1],m1.m$`F`[1],m1.m$`R2`[1],m1.m$`SumOfSqs`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m2",feature,m2$`Pr(>F)`[1],m2$`F`[1],m2$`R2`[1],m2$`SumOfSqs`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m2",feature,m2.f$`Pr(>F)`[1],m2.f$`F`[1],m2.f$`R2`[1],m2.f$`SumOfSqs`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m2",feature,m2.m$`Pr(>F)`[1],m2.m$`F`[1],m2.m$`R2`[1],m2.m$`SumOfSqs`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3",feature,m3$`Pr(>F)`[1],m3$`F`[1],m3$`R2`[1],m3$`SumOfSqs`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3",feature,m3.f$`Pr(>F)`[1],m3.f$`F`[1],m3.f$`R2`[1],m3.f$`SumOfSqs`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3",feature,m3.m$`Pr(>F)`[1],m3.m$`F`[1],m3.m$`R2`[1],m3.m$`SumOfSqs`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m4",feature,m4$`Pr(>F)`[1],m4$`F`[1],m4$`R2`[1],m4$`SumOfSqs`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m4",feature,m4.f$`Pr(>F)`[1],m4.f$`F`[1],m4.f$`R2`[1],m4.f$`SumOfSqs`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m4",feature,m4.m$`Pr(>F)`[1],m4.m$`F`[1],m4.m$`R2`[1],m4.m$`SumOfSqs`[1])
}
cor.matrix <- format_res_pnova(cor.matrix)
cor.matrix$feature <- factor(cor.matrix$feature, levels = c("Nighttime SBP SD","Daytime SBP SD","24-hour SBP SD"))
pdf("Figures/300_bp_variability/permanova_SBP_SD.pdf", width = 8, height = 5)
print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", name ="-log10" ) +
        coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender) +
        theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"),
              legend.title = element_text(size = 12, face="bold"), legend.text =  element_text(size = 12, face="bold")))
dev.off()
cor.matrix.full <- rbind(cor.matrix.full,cor.matrix)

cor.matrix <-  data.frame(gender = character(), model = character(), feature = character(), pval = character(), Fscore = character(), R2 = character(), SumOfSqs = character())
for(feature in SBP_surge){
  m1 <- permanova_model(merged.dat,species.data,feature)
  m1.f <- permanova_model(merged.dat[female,],species.data,feature)
  m1.m <- permanova_model(merged.dat[male,],species.data,feature)
  m2 <- permanova_model(merged.dat,species.data,c(feature,m2_covar,'Sex'))
  m2.f <- permanova_model(merged.dat[female,],species.data,c(feature,m2_covar))
  m2.m <- permanova_model(merged.dat[male,],species.data,c(feature,m2_covar))
  m3 <- permanova_model(merged.dat,species.data,c(feature,m3_covar_sbp,'Sex'))
  m3.f <- permanova_model(merged.dat[female,],species.data,c(feature,m3_covar_sbp,'Menopause_code'))
  m3.m <- permanova_model(merged.dat[male,],species.data,c(feature,m3_covar_sbp))
  m4 <- permanova_model(merged.dat,species.data,c(feature,m4_covar_sbp,'Sex'))
  m4.f <- permanova_model(merged.dat[female,],species.data,c(feature,m4_covar_sbp,'Menopause_code'))
  m4.m <- permanova_model(merged.dat[male,],species.data,c(feature,m4_covar_sbp))
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m1",feature,m1$`Pr(>F)`[1],m1$`F`[1],m1$`R2`[1],m1$`SumOfSqs`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m1",feature,m1.f$`Pr(>F)`[1],m1.f$`F`[1],m1.f$`R2`[1],m1.f$`SumOfSqs`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m1",feature,m1.m$`Pr(>F)`[1],m1.m$`F`[1],m1.m$`R2`[1],m1.m$`SumOfSqs`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m2",feature,m2$`Pr(>F)`[1],m2$`F`[1],m2$`R2`[1],m2$`SumOfSqs`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m2",feature,m2.f$`Pr(>F)`[1],m2.f$`F`[1],m2.f$`R2`[1],m2.f$`SumOfSqs`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m2",feature,m2.m$`Pr(>F)`[1],m2.m$`F`[1],m2.m$`R2`[1],m2.m$`SumOfSqs`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3",feature,m3$`Pr(>F)`[1],m3$`F`[1],m3$`R2`[1],m3$`SumOfSqs`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3",feature,m3.f$`Pr(>F)`[1],m3.f$`F`[1],m3.f$`R2`[1],m3.f$`SumOfSqs`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3",feature,m3.m$`Pr(>F)`[1],m3.m$`F`[1],m3.m$`R2`[1],m3.m$`SumOfSqs`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m4",feature,m4$`Pr(>F)`[1],m4$`F`[1],m4$`R2`[1],m4$`SumOfSqs`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m4",feature,m4.f$`Pr(>F)`[1],m4.f$`F`[1],m4.f$`R2`[1],m4.f$`SumOfSqs`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m4",feature,m4.m$`Pr(>F)`[1],m4.m$`F`[1],m4.m$`R2`[1],m4.m$`SumOfSqs`[1])
}
cor.matrix <- format_res_pnova(cor.matrix)
pdf("Figures/300_bp_variability/permanova_SBP_surge.pdf", width = 8, height = 5)
print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", name ="-log10" ) +
        coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender) +
        theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"),
              legend.title = element_text(size = 12, face="bold"), legend.text =  element_text(size = 12, face="bold")))
dev.off()
cor.matrix.full <- rbind(cor.matrix.full,cor.matrix)

cor.matrix <-  data.frame(gender = character(), model = character(), feature = character(), pval = character(), Fscore = character(), R2 = character(), SumOfSqs = character())
for(feature in DBP_cov){
  m1 <- permanova_model(merged.dat,species.data,feature)
  m1.f <- permanova_model(merged.dat[female,],species.data,feature)
  m1.m <- permanova_model(merged.dat[male,],species.data,feature)
  m2 <- permanova_model(merged.dat,species.data,c(feature,m2_covar,'Sex'))
  m2.f <- permanova_model(merged.dat[female,],species.data,c(feature,m2_covar))
  m2.m <- permanova_model(merged.dat[male,],species.data,c(feature,m2_covar))
  m3 <- permanova_model(merged.dat,species.data,c(feature,m3_covar_dbp,'Sex'))
  m3.f <- permanova_model(merged.dat[female,],species.data,c(feature,m3_covar_dbp,'Menopause_code'))
  m3.m <- permanova_model(merged.dat[male,],species.data,c(feature,m3_covar_dbp))
  m4 <- permanova_model(merged.dat,species.data,c(feature,m4_covar_dbp,'Sex'))
  m4.f <- permanova_model(merged.dat[female,],species.data,c(feature,m4_covar_dbp,'Menopause_code'))
  m4.m <- permanova_model(merged.dat[male,],species.data,c(feature,m4_covar_dbp))
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m1",feature,m1$`Pr(>F)`[1],m1$`F`[1],m1$`R2`[1],m1$`SumOfSqs`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m1",feature,m1.f$`Pr(>F)`[1],m1.f$`F`[1],m1.f$`R2`[1],m1.f$`SumOfSqs`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m1",feature,m1.m$`Pr(>F)`[1],m1.m$`F`[1],m1.m$`R2`[1],m1.m$`SumOfSqs`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m2",feature,m2$`Pr(>F)`[1],m2$`F`[1],m2$`R2`[1],m2$`SumOfSqs`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m2",feature,m2.f$`Pr(>F)`[1],m2.f$`F`[1],m2.f$`R2`[1],m2.f$`SumOfSqs`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m2",feature,m2.m$`Pr(>F)`[1],m2.m$`F`[1],m2.m$`R2`[1],m2.m$`SumOfSqs`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3",feature,m3$`Pr(>F)`[1],m3$`F`[1],m3$`R2`[1],m3$`SumOfSqs`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3",feature,m3.f$`Pr(>F)`[1],m3.f$`F`[1],m3.f$`R2`[1],m3.f$`SumOfSqs`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3",feature,m3.m$`Pr(>F)`[1],m3.m$`F`[1],m3.m$`R2`[1],m3.m$`SumOfSqs`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m4",feature,m4$`Pr(>F)`[1],m4$`F`[1],m4$`R2`[1],m4$`SumOfSqs`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m4",feature,m4.f$`Pr(>F)`[1],m4.f$`F`[1],m4.f$`R2`[1],m4.f$`SumOfSqs`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m4",feature,m4.m$`Pr(>F)`[1],m4.m$`F`[1],m4.m$`R2`[1],m4.m$`SumOfSqs`[1])
}
cor.matrix <- format_res_pnova(cor.matrix)
cor.matrix$feature <- factor(cor.matrix$feature, levels = c("Nighttime DBP CoV","Daytime DBP CoV","24-hour DBP CoV"))
pdf("Figures/300_bp_variability/permanova_DBP_CoV.pdf", width = 8, height = 5)
print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", name ="-log10" ) +
        coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender) +
        theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"),
              legend.title = element_text(size = 12, face="bold"), legend.text =  element_text(size = 12, face="bold")))
dev.off()
cor.matrix.full <- rbind(cor.matrix.full,cor.matrix)

cor.matrix <-  data.frame(gender = character(), model = character(), feature = character(), pval = character(), Fscore = character(), R2 = character(), SumOfSqs = character())
for(feature in DBP_sd){
  m1 <- permanova_model(merged.dat,species.data,feature)
  m1.f <- permanova_model(merged.dat[female,],species.data,feature)
  m1.m <- permanova_model(merged.dat[male,],species.data,feature)
  m2 <- permanova_model(merged.dat,species.data,c(feature,m2_covar,'Sex'))
  m2.f <- permanova_model(merged.dat[female,],species.data,c(feature,m2_covar))
  m2.m <- permanova_model(merged.dat[male,],species.data,c(feature,m2_covar))
  m3 <- permanova_model(merged.dat,species.data,c(feature,m3_covar_dbp,'Sex'))
  m3.f <- permanova_model(merged.dat[female,],species.data,c(feature,m3_covar_dbp,'Menopause_code'))
  m3.m <- permanova_model(merged.dat[male,],species.data,c(feature,m3_covar_dbp))
  m4 <- permanova_model(merged.dat,species.data,c(feature,m4_covar_dbp,'Sex'))
  m4.f <- permanova_model(merged.dat[female,],species.data,c(feature,m4_covar_dbp,'Menopause_code'))
  m4.m <- permanova_model(merged.dat[male,],species.data,c(feature,m4_covar_dbp))
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m1",feature,m1$`Pr(>F)`[1],m1$`F`[1],m1$`R2`[1],m1$`SumOfSqs`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m1",feature,m1.f$`Pr(>F)`[1],m1.f$`F`[1],m1.f$`R2`[1],m1.f$`SumOfSqs`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m1",feature,m1.m$`Pr(>F)`[1],m1.m$`F`[1],m1.m$`R2`[1],m1.m$`SumOfSqs`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m2",feature,m2$`Pr(>F)`[1],m2$`F`[1],m2$`R2`[1],m2$`SumOfSqs`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m2",feature,m2.f$`Pr(>F)`[1],m2.f$`F`[1],m2.f$`R2`[1],m2.f$`SumOfSqs`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m2",feature,m2.m$`Pr(>F)`[1],m2.m$`F`[1],m2.m$`R2`[1],m2.m$`SumOfSqs`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3",feature,m3$`Pr(>F)`[1],m3$`F`[1],m3$`R2`[1],m3$`SumOfSqs`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3",feature,m3.f$`Pr(>F)`[1],m3.f$`F`[1],m3.f$`R2`[1],m3.f$`SumOfSqs`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3",feature,m3.m$`Pr(>F)`[1],m3.m$`F`[1],m3.m$`R2`[1],m3.m$`SumOfSqs`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m4",feature,m4$`Pr(>F)`[1],m4$`F`[1],m4$`R2`[1],m4$`SumOfSqs`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m4",feature,m4.f$`Pr(>F)`[1],m4.f$`F`[1],m4.f$`R2`[1],m4.f$`SumOfSqs`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m4",feature,m4.m$`Pr(>F)`[1],m4.m$`F`[1],m4.m$`R2`[1],m4.m$`SumOfSqs`[1])
}
cor.matrix <- format_res_pnova(cor.matrix)
cor.matrix$feature <- factor(cor.matrix$feature, levels = c("Nighttime DBP SD","Daytime DBP SD","24-hour DBP SD"))
pdf("Figures/300_bp_variability/permanova_DBP_SD.pdf", width = 8, height = 5)
print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", name ="-log10" ) +
        coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender) +
        theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"),
              legend.title = element_text(size = 12, face="bold"), legend.text =  element_text(size = 12, face="bold")))
dev.off()
cor.matrix.full <- rbind(cor.matrix.full,cor.matrix)
write.csv(cor.matrix.full[,c('gender','model','feature','pval','Fscore','R2','SumOfSqs')], "Figures/300_bp_variability/permanova_BP_variability.csv",row.names = F)


#Plot MDS with permanova statistics
p.nova.all <- permanova_model(merged.dat,species.data,c("Dipping.classification.binary"))
p.nova.female <- permanova_model(merged.dat[female,],species.data,c("Dipping.classification.binary"))
p.nova.male <- permanova_model(merged.dat[male,],species.data,c("Dipping.classification.binary"))
p.nova.all <- p.nova.all[1,"Pr(>F)"]; p.nova.male <- p.nova.male[1,"Pr(>F)"]; p.nova.female <- p.nova.female[1,"Pr(>F)"]
p.labels <- data.frame(label = paste0("p = ",c(signif(p.nova.all,2),signif(p.nova.male,2),signif(p.nova.female,2))), gender= c("All","Men","Women"))
p.labels$gender <- factor(p.labels$gender, levels = c("All","Women","Men"))

bray.dist.matrix <- vegdist(species.data, method="bray")
mds <- cmdscale(bray.dist.matrix, eig = TRUE, x.ret = TRUE)
mds.var.per <- round(mds$eig/sum(mds$eig)*100, 1)
mds.values <- mds$points
mds.data <- data.frame(Patient = rownames(mds.values), X = mds.values[,1], Y = mds.values[,2], var = meta.data[match(rownames(mds.values),rownames(meta.data)),"Dipping.classification.binary"], gender=meta.data[match(rownames(mds.values),rownames(meta.data)),"Sex"])
all.facet <- mds.data
all.facet$gender <- "All"
mds.data <- rbind(mds.data,all.facet)
mds.data$var <- ifelse(mds.data$var==0,"Normal","Abnormal")
mds.data$gender <- factor(mds.data$gender, levels = c("All","Women","Men"))
pdf("Figures/300_bp_variability/Binary_dipping_MDS.pdf", width = 10, height = 4)
ggplot(mds.data, aes(x=X, y=Y, col=var)) + geom_point(size = 4, alpha = 0.8) + theme_classic() + facet_wrap(~gender) +
  xlab(paste0("MDS1 - ",mds.var.per[1],"%")) + ylab(paste0("MDS2 - ",mds.var.per[2],"%")) +
  geom_text(data = p.labels, mapping = aes(x = -Inf, y = -Inf, label = label), hjust = -0.1, vjust = -1, col = "red", size = 5) + 
  theme(axis.title = element_text(size = 16, face="bold"), axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=16, face="bold"),
        legend.title = element_blank(), legend.text =  element_text(size = 16, face="bold"), legend.position = "bottom")
dev.off()
}

# Correlation analysis w/ GM
{
pdf("Figures/300_bp_variability/BP_variability_vs_GM.pdf", width = 8, height = 5)
for(feat in c(SBP_cov,SBP_surge,SBP_sd)){
  cor.matrix <-  data.frame(gender = character(), model = character(), feature = character(), pval = character(), direction = character(),coefficient = character(), se = character())
  for(species in colnames(transformed.species.data)){
    m1 <- lm_model(merged.dat,feat,species)
    m1.f <- lm_model(merged.dat[female,],feat,species)
    m1.m <- lm_model(merged.dat[male,],feat,species)
    m2 <- lm_model(merged.dat,feat,c(species,m2_covar,'Sex'))
    m2.f <- lm_model(merged.dat[female,],feat,c(species,m2_covar))
    m2.m <- lm_model(merged.dat[male,],feat,c(species,m2_covar))
    m3 <- lm_model(merged.dat,feat,c(species,m3_covar_sbp,'Sex'))
    m3.f <- lm_model(merged.dat[female,],feat,c(species,m3_covar_sbp,'Menopause_code'))
    m3.m <- lm_model(merged.dat[male,],feat,c(species,m3_covar_sbp))
    m4 <- lm_model(merged.dat,feat,c(species,m4_covar_sbp,'Sex'))
    m4.f <- lm_model(merged.dat[female,],feat,c(species,m4_covar_sbp,'Menopause_code'))
    m4.m <- lm_model(merged.dat[male,],feat,c(species,m4_covar_sbp))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m1",species,summary(m1)$coefficients[2,4],sign(summary(m1)$coefficients[2,1]),summary(m1)$coefficients[2,1],summary(m1)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m1",species,summary(m1.f)$coefficients[2,4],sign(summary(m1.f)$coefficients[2,1]),summary(m1.f)$coefficients[2,1],summary(m1.f)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m1",species,summary(m1.m)$coefficients[2,4],sign(summary(m1.m)$coefficients[2,1]),summary(m1.m)$coefficients[2,1],summary(m1.m)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m2",species,summary(m2)$coefficients[2,4],sign(summary(m2)$coefficients[2,1]),summary(m2)$coefficients[2,1],summary(m2)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m2",species,summary(m2.f)$coefficients[2,4],sign(summary(m2.f)$coefficients[2,1]),summary(m2.f)$coefficients[2,1],summary(m2.f)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m2",species,summary(m2.m)$coefficients[2,4],sign(summary(m2.m)$coefficients[2,1]),summary(m2.m)$coefficients[2,1],summary(m2.m)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3",species,summary(m3)$coefficients[2,4],sign(summary(m3)$coefficients[2,1]),summary(m3)$coefficients[2,1],summary(m3)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3",species,summary(m3.f)$coefficients[2,4],sign(summary(m3.f)$coefficients[2,1]),summary(m3.f)$coefficients[2,1],summary(m3.f)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3",species,summary(m3.m)$coefficients[2,4],sign(summary(m3.m)$coefficients[2,1]),summary(m3.m)$coefficients[2,1],summary(m3.m)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m4",species,summary(m4)$coefficients[2,4],sign(summary(m4)$coefficients[2,1]),summary(m4)$coefficients[2,1],summary(m4)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m4",species,summary(m4.f)$coefficients[2,4],sign(summary(m4.f)$coefficients[2,1]),summary(m4.f)$coefficients[2,1],summary(m4.f)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m4",species,summary(m4.m)$coefficients[2,4],sign(summary(m4.m)$coefficients[2,1]),summary(m4.m)$coefficients[2,1],summary(m4.m)$coefficients[2,2])
  }
  cor.matrix$pval <- as.numeric(cor.matrix$pval)
  filt.feat <- cor.matrix[cor.matrix$pval<=0.01,"feature"]
  cor.matrix <- cor.matrix[cor.matrix[,"feature"]%in%filt.feat,]
  cor.matrix$sym <- ifelse(cor.matrix$pval>=0.05, "", ifelse(cor.matrix$pval<=0.01,ifelse(cor.matrix$pval<=0.001,"***", "**"), "*"))
  cor.matrix$log <- -log10(cor.matrix$pval)*as.numeric(cor.matrix$direction)
  cor.matrix$feature <- gsub("_"," ",cor.matrix$feature)
  title <- gsub("\\."," ",feat)
  title <- gsub("X","",title)
  title <- gsub("_","-",title)
  cor.matrix$gender <- factor(cor.matrix$gender, levels = c("All","Women","Men"))
  print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(name = "-log10", midpoint=0, low="blue", mid="white",high="red", space ="Lab") +
    coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender)  +
    theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"),
          legend.text = element_text(size = 12, face="bold"), legend.title = element_text(size = 12, face="bold")) + ggtitle(title))
  write.csv(cor.matrix[,c('gender','model','feature','pval','coefficient','se')], paste0("Figures/300_bp_variability/",feat,"_vs_GM.csv"),row.names = F)
}
for(feat in c(DBP_cov,DBP_sd)){
  #cor.matrix <-  data.frame(gender = character(), model = character(), feature = character(), pval = character(), direction = character())
  cor.matrix <-  data.frame(gender = character(), model = character(), feature = character(), pval = character(), direction = character(),coefficient = character(), se = character())
  for(species in colnames(transformed.species.data)){
    m1 <- lm_model(merged.dat,feat,species)
    m1.f <- lm_model(merged.dat[female,],feat,species)
    m1.m <- lm_model(merged.dat[male,],feat,species)
    m2 <- lm_model(merged.dat,feat,c(species,m2_covar,'Sex'))
    m2.f <- lm_model(merged.dat[female,],feat,c(species,m2_covar))
    m2.m <- lm_model(merged.dat[male,],feat,c(species,m2_covar))
    m3 <- lm_model(merged.dat,feat,c(species,m3_covar_dbp,'Sex'))
    m3.f <- lm_model(merged.dat[female,],feat,c(species,m3_covar_dbp,'Menopause_code'))
    m3.m <- lm_model(merged.dat[male,],feat,c(species,m3_covar_dbp))
    m4 <- lm_model(merged.dat,feat,c(species,m4_covar_dbp,'Sex'))
    m4.f <- lm_model(merged.dat[female,],feat,c(species,m4_covar_dbp,'Menopause_code'))
    m4.m <- lm_model(merged.dat[male,],feat,c(species,m4_covar_dbp))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m1",species,summary(m1)$coefficients[2,4],sign(summary(m1)$coefficients[2,1]),summary(m1)$coefficients[2,1],summary(m1)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m1",species,summary(m1.f)$coefficients[2,4],sign(summary(m1.f)$coefficients[2,1]),summary(m1.f)$coefficients[2,1],summary(m1.f)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m1",species,summary(m1.m)$coefficients[2,4],sign(summary(m1.m)$coefficients[2,1]),summary(m1.m)$coefficients[2,1],summary(m1.m)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m2",species,summary(m2)$coefficients[2,4],sign(summary(m2)$coefficients[2,1]),summary(m2)$coefficients[2,1],summary(m2)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m2",species,summary(m2.f)$coefficients[2,4],sign(summary(m2.f)$coefficients[2,1]),summary(m2.f)$coefficients[2,1],summary(m2.f)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m2",species,summary(m2.m)$coefficients[2,4],sign(summary(m2.m)$coefficients[2,1]),summary(m2.m)$coefficients[2,1],summary(m2.m)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3",species,summary(m3)$coefficients[2,4],sign(summary(m3)$coefficients[2,1]),summary(m3)$coefficients[2,1],summary(m3)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3",species,summary(m3.f)$coefficients[2,4],sign(summary(m3.f)$coefficients[2,1]),summary(m3.f)$coefficients[2,1],summary(m3.f)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3",species,summary(m3.m)$coefficients[2,4],sign(summary(m3.m)$coefficients[2,1]),summary(m3.m)$coefficients[2,1],summary(m3.m)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m4",species,summary(m4)$coefficients[2,4],sign(summary(m4)$coefficients[2,1]),summary(m4)$coefficients[2,1],summary(m4)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m4",species,summary(m4.f)$coefficients[2,4],sign(summary(m4.f)$coefficients[2,1]),summary(m4.f)$coefficients[2,1],summary(m4.f)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m4",species,summary(m4.m)$coefficients[2,4],sign(summary(m4.m)$coefficients[2,1]),summary(m4.m)$coefficients[2,1],summary(m4.m)$coefficients[2,2])
  }
  cor.matrix$pval <- as.numeric(cor.matrix$pval)
  filt.feat <- cor.matrix[cor.matrix$pval<=0.01,"feature"]
  cor.matrix <- cor.matrix[cor.matrix[,"feature"]%in%filt.feat,]
  cor.matrix$sym <- ifelse(cor.matrix$pval>=0.05, "", ifelse(cor.matrix$pval<=0.01,ifelse(cor.matrix$pval<=0.001,"***", "**"), "*"))
  cor.matrix$log <- -log10(cor.matrix$pval)*as.numeric(cor.matrix$direction)
  cor.matrix$feature <- gsub("_"," ",cor.matrix$feature)
  title <- gsub("\\."," ",feat)
  title <- gsub("X","",title)
  title <- gsub("_","-",title)
  cor.matrix$gender <- factor(cor.matrix$gender, levels = c("All","Women","Men"))
  print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(name = "-log10", midpoint=0, low="blue", mid="white",high="red", space ="Lab") +
          coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender)  +
          theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"),
                legend.text = element_text(size = 12, face="bold"), legend.title = element_text(size = 12, face="bold")) + ggtitle(title))
  write.csv(cor.matrix[,c('gender','model','feature','pval','coefficient','se')], paste0("Figures/300_bp_variability/",feat,"_vs_GM.csv"),row.names = F)
}
dev.off()

dip_subset <- merged.dat[merged.dat$Dipping.classification%in%c(1,2,3,4),]
pdf("Figures/300_bp_variability/Dipping_vs_GM_boxplot.pdf", width = 8, height = 6)
for(feat in colnames(species.data)){
  my_boxplot(dip_subset,feat,"Dipping.classification",ylab = feat)
}
dev.off()

#For figures
#24hr BP CoV vs GM
pdf("Figures/300_bp_variability/24hr_BP_CoV_vs_GM.pdf", width = 8, height = 6)
cor.matrix <-  data.frame(gender = character(), model = character(), feature = character(), pval = character(), direction = character(),coefficient = character(), se = character())
for(species in colnames(transformed.species.data)){
  m1 <- lm_model(merged.dat,"X24_hour.SBP.CoV",species)
  m1.f <- lm_model(merged.dat[female,],"X24_hour.SBP.CoV",species)
  m1.m <- lm_model(merged.dat[male,],"X24_hour.SBP.CoV",species)
  m2 <- lm_model(merged.dat,"X24_hour.SBP.CoV",c(species,m2_covar,'Sex'))
  m2.f <- lm_model(merged.dat[female,],"X24_hour.SBP.CoV",c(species,m2_covar))
  m2.m <- lm_model(merged.dat[male,],"X24_hour.SBP.CoV",c(species,m2_covar))
  m3 <- lm_model(merged.dat,"X24_hour.SBP.CoV",c(species,m3_covar_sbp,'Sex'))
  m3.f <- lm_model(merged.dat[female,],"X24_hour.SBP.CoV",c(species,m3_covar_sbp,'Menopause_code'))
  m3.m <- lm_model(merged.dat[male,],"X24_hour.SBP.CoV",c(species,m3_covar_sbp))
  m4 <- lm_model(merged.dat,"X24_hour.SBP.CoV",c(species,m4_covar_sbp,'Sex'))
  m4.f <- lm_model(merged.dat[female,],"X24_hour.SBP.CoV",c(species,m4_covar_sbp,'Menopause_code'))
  m4.m <- lm_model(merged.dat[male,],"X24_hour.SBP.CoV",c(species,m4_covar_sbp))
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m1",species,summary(m1)$coefficients[2,4],sign(summary(m1)$coefficients[2,1]),summary(m1)$coefficients[2,1],summary(m1)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m1",species,summary(m1.f)$coefficients[2,4],sign(summary(m1.f)$coefficients[2,1]),summary(m1.f)$coefficients[2,1],summary(m1.f)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m1",species,summary(m1.m)$coefficients[2,4],sign(summary(m1.m)$coefficients[2,1]),summary(m1.m)$coefficients[2,1],summary(m1.m)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m2",species,summary(m2)$coefficients[2,4],sign(summary(m2)$coefficients[2,1]),summary(m2)$coefficients[2,1],summary(m2)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m2",species,summary(m2.f)$coefficients[2,4],sign(summary(m2.f)$coefficients[2,1]),summary(m2.f)$coefficients[2,1],summary(m2.f)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m2",species,summary(m2.m)$coefficients[2,4],sign(summary(m2.m)$coefficients[2,1]),summary(m2.m)$coefficients[2,1],summary(m2.m)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3",species,summary(m3)$coefficients[2,4],sign(summary(m3)$coefficients[2,1]),summary(m3)$coefficients[2,1],summary(m3)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3",species,summary(m3.f)$coefficients[2,4],sign(summary(m3.f)$coefficients[2,1]),summary(m3.f)$coefficients[2,1],summary(m3.f)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3",species,summary(m3.m)$coefficients[2,4],sign(summary(m3.m)$coefficients[2,1]),summary(m3.m)$coefficients[2,1],summary(m3.m)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m4",species,summary(m4)$coefficients[2,4],sign(summary(m4)$coefficients[2,1]),summary(m4)$coefficients[2,1],summary(m4)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m4",species,summary(m4.f)$coefficients[2,4],sign(summary(m4.f)$coefficients[2,1]),summary(m4.f)$coefficients[2,1],summary(m4.f)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m4",species,summary(m4.m)$coefficients[2,4],sign(summary(m4.m)$coefficients[2,1]),summary(m4.m)$coefficients[2,1],summary(m4.m)$coefficients[2,2])
}
cor.matrix <- format_res(cor.matrix, filt = T, filt_thres = 0.01)
print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(name = "-log10", midpoint=0, low="blue", mid="white",high="red", space ="Lab", limits = c(-2.9,2.9)) +
        coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender)  +
        theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"),
              legend.text = element_text(size = 12, face="bold"), legend.title = element_text(size = 12, face="bold")) + ggtitle("24-hour SBP CoV"))

cor.matrix <-  data.frame(gender = character(), model = character(), feature = character(), pval = character(), direction = character(),coefficient = character(), se = character())
for(species in colnames(transformed.species.data)){
  m1 <- lm_model(merged.dat,"X24_hour.DBP.CoV",species)
  m1.f <- lm_model(merged.dat[female,],"X24_hour.DBP.CoV",species)
  m1.m <- lm_model(merged.dat[male,],"X24_hour.DBP.CoV",species)
  m2 <- lm_model(merged.dat,"X24_hour.DBP.CoV",c(species,m2_covar,'Sex'))
  m2.f <- lm_model(merged.dat[female,],"X24_hour.DBP.CoV",c(species,m2_covar))
  m2.m <- lm_model(merged.dat[male,],"X24_hour.DBP.CoV",c(species,m2_covar))
  m3 <- lm_model(merged.dat,"X24_hour.DBP.CoV",c(species,m3_covar_dbp,'Sex'))
  m3.f <- lm_model(merged.dat[female,],"X24_hour.DBP.CoV",c(species,m3_covar_dbp,'Menopause_code'))
  m3.m <- lm_model(merged.dat[male,],"X24_hour.DBP.CoV",c(species,m3_covar_dbp))
  m4 <- lm_model(merged.dat,"X24_hour.DBP.CoV",c(species,m4_covar_dbp,'Sex'))
  m4.f <- lm_model(merged.dat[female,],"X24_hour.DBP.CoV",c(species,m4_covar_dbp,'Menopause_code'))
  m4.m <- lm_model(merged.dat[male,],"X24_hour.DBP.CoV",c(species,m4_covar_dbp))
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m1",species,summary(m1)$coefficients[2,4],sign(summary(m1)$coefficients[2,1]),summary(m1)$coefficients[2,1],summary(m1)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m1",species,summary(m1.f)$coefficients[2,4],sign(summary(m1.f)$coefficients[2,1]),summary(m1.f)$coefficients[2,1],summary(m1.f)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m1",species,summary(m1.m)$coefficients[2,4],sign(summary(m1.m)$coefficients[2,1]),summary(m1.m)$coefficients[2,1],summary(m1.m)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m2",species,summary(m2)$coefficients[2,4],sign(summary(m2)$coefficients[2,1]),summary(m2)$coefficients[2,1],summary(m2)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m2",species,summary(m2.f)$coefficients[2,4],sign(summary(m2.f)$coefficients[2,1]),summary(m2.f)$coefficients[2,1],summary(m2.f)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m2",species,summary(m2.m)$coefficients[2,4],sign(summary(m2.m)$coefficients[2,1]),summary(m2.m)$coefficients[2,1],summary(m2.m)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3",species,summary(m3)$coefficients[2,4],sign(summary(m3)$coefficients[2,1]),summary(m3)$coefficients[2,1],summary(m3)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3",species,summary(m3.f)$coefficients[2,4],sign(summary(m3.f)$coefficients[2,1]),summary(m3.f)$coefficients[2,1],summary(m3.f)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3",species,summary(m3.m)$coefficients[2,4],sign(summary(m3.m)$coefficients[2,1]),summary(m3.m)$coefficients[2,1],summary(m3.m)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m4",species,summary(m4)$coefficients[2,4],sign(summary(m4)$coefficients[2,1]),summary(m4)$coefficients[2,1],summary(m4)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m4",species,summary(m4.f)$coefficients[2,4],sign(summary(m4.f)$coefficients[2,1]),summary(m4.f)$coefficients[2,1],summary(m4.f)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m4",species,summary(m4.m)$coefficients[2,4],sign(summary(m4.m)$coefficients[2,1]),summary(m4.m)$coefficients[2,1],summary(m4.m)$coefficients[2,2])
}
cor.matrix <- format_res(cor.matrix, filt = T, filt_thres = 0.01)
print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(name = "-log10", midpoint=0, low="blue", mid="white",high="red", space ="Lab", limits = c(-2.9,2.9)) +
        coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender)  +
        theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"),
              legend.text = element_text(size = 12, face="bold"), legend.title = element_text(size = 12, face="bold")) + ggtitle("24-hour DBP CoV"))

dev.off()

#SBP Surge  vs GM
pdf("Figures/300_bp_variability/SBP_surge_vs_GM.pdf", width = 8, height = 6)
cor.matrix <-  data.frame(gender = character(), model = character(), feature = character(), pval = character(), direction = character(),coefficient = character(), se = character())
for(species in colnames(transformed.species.data)){
  m1 <- lm_model(merged.dat,"Sleep.through.MSBPS",species)
  m1.f <- lm_model(merged.dat[female,],"Sleep.through.MSBPS",species)
  m1.m <- lm_model(merged.dat[male,],"Sleep.through.MSBPS",species)
  m2 <- lm_model(merged.dat,"Sleep.through.MSBPS",c(species,m2_covar,'Sex'))
  m2.f <- lm_model(merged.dat[female,],"Sleep.through.MSBPS",c(species,m2_covar))
  m2.m <- lm_model(merged.dat[male,],"Sleep.through.MSBPS",c(species,m2_covar))
  m3 <- lm_model(merged.dat,"Sleep.through.MSBPS",c(species,m3_covar_sbp,'Sex'))
  m3.f <- lm_model(merged.dat[female,],"Sleep.through.MSBPS",c(species,m3_covar_sbp,'Menopause_code'))
  m3.m <- lm_model(merged.dat[male,],"Sleep.through.MSBPS",c(species,m3_covar_sbp))
  m4 <- lm_model(merged.dat,"Sleep.through.MSBPS",c(species,m4_covar,'Sex'))
  m4.f <- lm_model(merged.dat[female,],"Sleep.through.MSBPS",c(species,m4_covar_sbp,'Menopause_code'))
  m4.m <- lm_model(merged.dat[male,],"Sleep.through.MSBPS",c(species,m4_covar_sbp))
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m1",species,summary(m1)$coefficients[2,4],sign(summary(m1)$coefficients[2,1]),summary(m1)$coefficients[2,1],summary(m1)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m1",species,summary(m1.f)$coefficients[2,4],sign(summary(m1.f)$coefficients[2,1]),summary(m1.f)$coefficients[2,1],summary(m1.f)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m1",species,summary(m1.m)$coefficients[2,4],sign(summary(m1.m)$coefficients[2,1]),summary(m1.m)$coefficients[2,1],summary(m1.m)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m2",species,summary(m2)$coefficients[2,4],sign(summary(m2)$coefficients[2,1]),summary(m2)$coefficients[2,1],summary(m2)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m2",species,summary(m2.f)$coefficients[2,4],sign(summary(m2.f)$coefficients[2,1]),summary(m2.f)$coefficients[2,1],summary(m2.f)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m2",species,summary(m2.m)$coefficients[2,4],sign(summary(m2.m)$coefficients[2,1]),summary(m2.m)$coefficients[2,1],summary(m2.m)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3",species,summary(m3)$coefficients[2,4],sign(summary(m3)$coefficients[2,1]),summary(m3)$coefficients[2,1],summary(m3)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3",species,summary(m3.f)$coefficients[2,4],sign(summary(m3.f)$coefficients[2,1]),summary(m3.f)$coefficients[2,1],summary(m3.f)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3",species,summary(m3.m)$coefficients[2,4],sign(summary(m3.m)$coefficients[2,1]),summary(m3.m)$coefficients[2,1],summary(m3.m)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m4",species,summary(m4)$coefficients[2,4],sign(summary(m4)$coefficients[2,1]),summary(m4)$coefficients[2,1],summary(m4)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m4",species,summary(m4.f)$coefficients[2,4],sign(summary(m4.f)$coefficients[2,1]),summary(m4.f)$coefficients[2,1],summary(m4.f)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m4",species,summary(m4.m)$coefficients[2,4],sign(summary(m4.m)$coefficients[2,1]),summary(m4.m)$coefficients[2,1],summary(m4.m)$coefficients[2,2])
}
cor.matrix <- format_res(cor.matrix, filt = T, filt_thres = 0.01)
print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(name = "-log10", midpoint=0, low="blue", mid="white",high="red", space ="Lab", limits = c(-3.5,4.2)) +
        coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender)  +
        theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"),
              legend.text = element_text(size = 12, face="bold"), legend.title = element_text(size = 12, face="bold")) + ggtitle("Sleep through MSBPS"))

cor.matrix <-  data.frame(gender = character(), model = character(), feature = character(), pval = character(), direction = character(),coefficient = character(), se = character())
for(species in colnames(transformed.species.data)){
  m1 <- lm_model(merged.dat,"Prewaking.MSBPS",species)
  m1.f <- lm_model(merged.dat[female,],"Prewaking.MSBPS",species)
  m1.m <- lm_model(merged.dat[male,],"Prewaking.MSBPS",species)
  m2 <- lm_model(merged.dat,"Prewaking.MSBPS",c(species,m2_covar,'Sex'))
  m2.f <- lm_model(merged.dat[female,],"Prewaking.MSBPS",c(species,m2_covar))
  m2.m <- lm_model(merged.dat[male,],"Prewaking.MSBPS",c(species,m2_covar))
  m3 <- lm_model(merged.dat,"Prewaking.MSBPS",c(species,m3_covar_sbp,'Sex'))
  m3.f <- lm_model(merged.dat[female,],"Prewaking.MSBPS",c(species,m3_covar_sbp,'Menopause_code'))
  m3.m <- lm_model(merged.dat[male,],"Prewaking.MSBPS",c(species,m3_covar_sbp))
  m4 <- lm_model(merged.dat,"Prewaking.MSBPS",c(species,m4_covar_sbp,'Sex'))
  m4.f <- lm_model(merged.dat[female,],"Prewaking.MSBPS",c(species,m4_covar_sbp,'Menopause_code'))
  m4.m <- lm_model(merged.dat[male,],"Prewaking.MSBPS",c(species,m4_covar_sbp))
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m1",species,summary(m1)$coefficients[2,4],sign(summary(m1)$coefficients[2,1]),summary(m1)$coefficients[2,1],summary(m1)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m1",species,summary(m1.f)$coefficients[2,4],sign(summary(m1.f)$coefficients[2,1]),summary(m1.f)$coefficients[2,1],summary(m1.f)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m1",species,summary(m1.m)$coefficients[2,4],sign(summary(m1.m)$coefficients[2,1]),summary(m1.m)$coefficients[2,1],summary(m1.m)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m2",species,summary(m2)$coefficients[2,4],sign(summary(m2)$coefficients[2,1]),summary(m2)$coefficients[2,1],summary(m2)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m2",species,summary(m2.f)$coefficients[2,4],sign(summary(m2.f)$coefficients[2,1]),summary(m2.f)$coefficients[2,1],summary(m2.f)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m2",species,summary(m2.m)$coefficients[2,4],sign(summary(m2.m)$coefficients[2,1]),summary(m2.m)$coefficients[2,1],summary(m2.m)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3",species,summary(m3)$coefficients[2,4],sign(summary(m3)$coefficients[2,1]),summary(m3)$coefficients[2,1],summary(m3)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3",species,summary(m3.f)$coefficients[2,4],sign(summary(m3.f)$coefficients[2,1]),summary(m3.f)$coefficients[2,1],summary(m3.f)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3",species,summary(m3.m)$coefficients[2,4],sign(summary(m3.m)$coefficients[2,1]),summary(m3.m)$coefficients[2,1],summary(m3.m)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m4",species,summary(m4)$coefficients[2,4],sign(summary(m4)$coefficients[2,1]),summary(m4)$coefficients[2,1],summary(m4)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m4",species,summary(m4.f)$coefficients[2,4],sign(summary(m4.f)$coefficients[2,1]),summary(m4.f)$coefficients[2,1],summary(m4.f)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m4",species,summary(m4.m)$coefficients[2,4],sign(summary(m4.m)$coefficients[2,1]),summary(m4.m)$coefficients[2,1],summary(m4.m)$coefficients[2,2])
}
cor.matrix <- format_res(cor.matrix, filt = T, filt_thres = 0.01)
print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(name = "-log10", midpoint=0, low="blue", mid="white",high="red", space ="Lab", limits = c(-3.5,4.2)) +
        coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender)  +
        theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"),
              legend.text = element_text(size = 12, face="bold"), legend.title = element_text(size = 12, face="bold")) + ggtitle("Prewaking MSBPS"))

dev.off()
}

# Correlation analysis w/ Immune markers
{
pdf("Figures/300_bp_variability/BP_variability_vs_immune.pdf", width = 8, height = 5)
for(feat in c(SBP_cov,SBP_sd,SBP_surge)){
  cor.matrix <-  data.frame(gender = character(), model = character(), feature = character(), pval = character(), direction = character())
  for(species in immune_markers){
    m1 <- lm_model(merged.dat,feat,species)
    m1.f <- lm_model(merged.dat[female,],feat,species)
    m1.m <- lm_model(merged.dat[male,],feat,species)
    m2 <- lm_model(merged.dat,feat,c(species,m2_covar,'Sex'))
    m2.f <- lm_model(merged.dat[female,],feat,c(species,m2_covar))
    m2.m <- lm_model(merged.dat[male,],feat,c(species,m2_covar))
    m3 <- lm_model(merged.dat,feat,c(species,m3_covar_sbp,'Sex'))
    m3.f <- lm_model(merged.dat[female,],feat,c(species,m3_covar_sbp,'Menopause_code'))
    m3.m <- lm_model(merged.dat[male,],feat,c(species,m3_covar_sbp))
    m4 <- lm_model(merged.dat,feat,c(species,m4_covar_sbp,'Sex'))
    m4.f <- lm_model(merged.dat[female,],feat,c(species,m4_covar_sbp,'Menopause_code'))
    m4.m <- lm_model(merged.dat[male,],feat,c(species,m4_covar_sbp))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m1",species,summary(m1)$coefficients[2,4],sign(summary(m1)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m1",species,summary(m1.f)$coefficients[2,4],sign(summary(m1.f)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m1",species,summary(m1.m)$coefficients[2,4],sign(summary(m1.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m2",species,summary(m2)$coefficients[2,4],sign(summary(m2)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m2",species,summary(m2.f)$coefficients[2,4],sign(summary(m2.f)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m2",species,summary(m2.m)$coefficients[2,4],sign(summary(m2.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3",species,summary(m3)$coefficients[2,4],sign(summary(m3)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3",species,summary(m3.f)$coefficients[2,4],sign(summary(m3.f)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3",species,summary(m3.m)$coefficients[2,4],sign(summary(m3.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m4",species,summary(m4)$coefficients[2,4],sign(summary(m4)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m4",species,summary(m4.f)$coefficients[2,4],sign(summary(m4.f)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m4",species,summary(m4.m)$coefficients[2,4],sign(summary(m4.m)$coefficients[2,1]))
  }
  cor.matrix$pval <- as.numeric(cor.matrix$pval)
  filt.feat <- cor.matrix[cor.matrix$pval<=0.05,"feature"]
  if(length(filt.feat)==0){next}
  cor.matrix <- cor.matrix[cor.matrix[,"feature"]%in%filt.feat,]
  cor.matrix$sym <- ifelse(cor.matrix$pval>=0.05, "", ifelse(cor.matrix$pval<=0.01,ifelse(cor.matrix$pval<=0.001,"***", "**"), "*"))
  cor.matrix$log <- -log10(cor.matrix$pval)*as.numeric(cor.matrix$direction)
  cor.matrix$feature <- gsub("_"," ",cor.matrix$feature)
  title <- gsub("\\."," ",feat)
  title <- gsub("X","",title)
  title <- gsub("_","-",title)
  cor.matrix$gender <- factor(cor.matrix$gender, levels = c("All","Women","Men"))
  print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(name = "-log10", midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
          coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender)  +
          theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"),
                legend.text = element_text(size = 12, face="bold"), legend.title = element_text(size = 12, face="bold")) + ggtitle(title))
}
for(feat in c(DBP_cov,DBP_sd)){
  cor.matrix <-  data.frame(gender = character(), model = character(), feature = character(), pval = character(), direction = character())
  for(species in immune_markers){
    m1 <- lm_model(merged.dat,feat,species)
    m1.f <- lm_model(merged.dat[female,],feat,species)
    m1.m <- lm_model(merged.dat[male,],feat,species)
    m2 <- lm_model(merged.dat,feat,c(species,m2_covar,'Sex'))
    m2.f <- lm_model(merged.dat[female,],feat,c(species,m2_covar))
    m2.m <- lm_model(merged.dat[male,],feat,c(species,m2_covar))
    m3 <- lm_model(merged.dat,feat,c(species,m3_covar_dbp,'Sex'))
    m3.f <- lm_model(merged.dat[female,],feat,c(species,m3_covar_dbp,'Menopause_code'))
    m3.m <- lm_model(merged.dat[male,],feat,c(species,m3_covar_dbp))
    m4 <- lm_model(merged.dat,feat,c(species,m4_covar_dbp,'Sex'))
    m4.f <- lm_model(merged.dat[female,],feat,c(species,m4_covar_dbp,'Menopause_code'))
    m4.m <- lm_model(merged.dat[male,],feat,c(species,m4_covar_dbp))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m1",species,summary(m1)$coefficients[2,4],sign(summary(m1)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m1",species,summary(m1.f)$coefficients[2,4],sign(summary(m1.f)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m1",species,summary(m1.m)$coefficients[2,4],sign(summary(m1.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m2",species,summary(m2)$coefficients[2,4],sign(summary(m2)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m2",species,summary(m2.f)$coefficients[2,4],sign(summary(m2.f)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m2",species,summary(m2.m)$coefficients[2,4],sign(summary(m2.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3",species,summary(m3)$coefficients[2,4],sign(summary(m3)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3",species,summary(m3.f)$coefficients[2,4],sign(summary(m3.f)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3",species,summary(m3.m)$coefficients[2,4],sign(summary(m3.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m4",species,summary(m4)$coefficients[2,4],sign(summary(m4)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m4",species,summary(m4.f)$coefficients[2,4],sign(summary(m4.f)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m4",species,summary(m4.m)$coefficients[2,4],sign(summary(m4.m)$coefficients[2,1]))
  }
  cor.matrix$pval <- as.numeric(cor.matrix$pval)
  filt.feat <- cor.matrix[cor.matrix$pval<=0.05,"feature"]
  if(length(filt.feat)==0){next}
  cor.matrix <- cor.matrix[cor.matrix[,"feature"]%in%filt.feat,]
  cor.matrix$sym <- ifelse(cor.matrix$pval>=0.05, "", ifelse(cor.matrix$pval<=0.01,ifelse(cor.matrix$pval<=0.001,"***", "**"), "*"))
  cor.matrix$log <- -log10(cor.matrix$pval)*as.numeric(cor.matrix$direction)
  cor.matrix$feature <- gsub("_"," ",cor.matrix$feature)
  title <- gsub("\\."," ",feat)
  title <- gsub("X","",title)
  title <- gsub("_","-",title)
  cor.matrix$gender <- factor(cor.matrix$gender, levels = c("All","Women","Men"))
  print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(name = "-log10", midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
          coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender)  +
          theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"),
                legend.text = element_text(size = 12, face="bold"), legend.title = element_text(size = 12, face="bold")) + ggtitle(title))
}
dev.off()

pdf("Figures/300_bp_variability/Immune_vs_GM.pdf", width = 8, height = 5)
for(feat in immune_markers){
  cor.matrix <-  data.frame(gender = character(), model = character(), feature = character(), pval = character(), direction = character())
  for(species in colnames(transformed.species.data)){
    m1 <- lm_model(merged.dat,feat,species)
    m1.f <- lm_model(merged.dat[female,],feat,species)
    m1.m <- lm_model(merged.dat[male,],feat,species)
    m2 <- lm_model(merged.dat,feat,c(species,m2_covar,'Sex'))
    m2.f <- lm_model(merged.dat[female,],feat,c(species,m2_covar))
    m2.m <- lm_model(merged.dat[male,],feat,c(species,m2_covar))
    m3 <- lm_model(merged.dat,feat,c(species,m3_covar,'Sex'))
    m3.f <- lm_model(merged.dat[female,],feat,c(species,m3_covar,'Menopause_code'))
    m3.m <- lm_model(merged.dat[male,],feat,c(species,m3_covar))
    m4 <- lm_model(merged.dat,feat,c(species,m4_covar,'Sex'))
    m4.f <- lm_model(merged.dat[female,],feat,c(species,m4_covar,'Menopause_code'))
    m4.m <- lm_model(merged.dat[male,],feat,c(species,m4_covar))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m1",species,summary(m1)$coefficients[2,4],sign(summary(m1)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m1",species,summary(m1.f)$coefficients[2,4],sign(summary(m1.f)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m1",species,summary(m1.m)$coefficients[2,4],sign(summary(m1.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m2",species,summary(m2)$coefficients[2,4],sign(summary(m2)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m2",species,summary(m2.f)$coefficients[2,4],sign(summary(m2.f)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m2",species,summary(m2.m)$coefficients[2,4],sign(summary(m2.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3",species,summary(m3)$coefficients[2,4],sign(summary(m3)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3",species,summary(m3.f)$coefficients[2,4],sign(summary(m3.f)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3",species,summary(m3.m)$coefficients[2,4],sign(summary(m3.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m4",species,summary(m4)$coefficients[2,4],sign(summary(m4)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m4",species,summary(m4.f)$coefficients[2,4],sign(summary(m4.f)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m4",species,summary(m4.m)$coefficients[2,4],sign(summary(m4.m)$coefficients[2,1]))
  }
  cor.matrix$pval <- as.numeric(cor.matrix$pval)
  filt.feat <- cor.matrix[cor.matrix$pval<=0.05,"feature"]
  if(length(filt.feat)==0){next}
  cor.matrix <- cor.matrix[cor.matrix[,"feature"]%in%filt.feat,]
  cor.matrix$sym <- ifelse(cor.matrix$pval>=0.05, "", ifelse(cor.matrix$pval<=0.01,ifelse(cor.matrix$pval<=0.001,"***", "**"), "*"))
  cor.matrix$log <- -log10(cor.matrix$pval)*as.numeric(cor.matrix$direction)
  cor.matrix$feature <- gsub("_"," ",cor.matrix$feature)
  title <- gsub("\\."," ",feat)
  title <- gsub("X","",title)
  title <- gsub("_","-",title)
  cor.matrix$gender <- factor(cor.matrix$gender, levels = c("All","Women","Men"))
  print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(name = "-log10", midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
          coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender)  +
          theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"),
                legend.text = element_text(size = 12, face="bold"), legend.title = element_text(size = 12, face="bold")) + ggtitle(title))
}

}

# Correlation analysis w/ Blood SCFAs
{
cor.matrix.full <- data.frame(feature = character(), gender = character(), model = character(), species = character(), pval = character(), direction = character(),coefficient = character(), se = character())
pdf("Figures/300_bp_variability/BP_variability_vs_blood_SCFA.pdf", width = 8, height = 5)
for(feat in c(SBP_cov,SBP_sd,SBP_surge)){
  cor.matrix <-  data.frame(feature = character(), gender = character(), model = character(), species = character(), pval = character(), direction = character(),coefficient = character(), se = character())
  for(species in acids_variables){
    m1 <- lm_model(merged.dat,feat,species)
    m1.f <- lm_model(merged.dat[female,],feat,species)
    m1.m <- lm_model(merged.dat[male,],feat,species)
    m2 <- lm_model(merged.dat,feat,c(species,m2_covar,'Sex'))
    m2.f <- lm_model(merged.dat[female,],feat,c(species,m2_covar))
    m2.m <- lm_model(merged.dat[male,],feat,c(species,m2_covar))
    m3 <- lm_model(merged.dat,feat,c(species,m3_covar_sbp,'Sex'))
    m3.f <- lm_model(merged.dat[female,],feat,c(species,m3_covar_sbp,'Menopause_code'))
    m3.m <- lm_model(merged.dat[male,],feat,c(species,m3_covar_sbp))
    m4 <- lm_model(merged.dat,feat,c(species,m4_covar_sbp,'Sex'))
    m4.f <- lm_model(merged.dat[female,],feat,c(species,m4_covar_sbp,'Menopause_code'))
    m4.m <- lm_model(merged.dat[male,],feat,c(species,m4_covar_sbp))
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"All","m1",species,summary(m1)$coefficients[2,4],sign(summary(m1)$coefficients[2,1]),summary(m1)$coefficients[2,1],summary(m1)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"Women","m1",species,summary(m1.f)$coefficients[2,4],sign(summary(m1.f)$coefficients[2,1]),summary(m1.f)$coefficients[2,1],summary(m1.f)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"Men","m1",species,summary(m1.m)$coefficients[2,4],sign(summary(m1.m)$coefficients[2,1]),summary(m1.m)$coefficients[2,1],summary(m1.m)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"All","m2",species,summary(m2)$coefficients[2,4],sign(summary(m2)$coefficients[2,1]),summary(m2)$coefficients[2,1],summary(m2)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"Women","m2",species,summary(m2.f)$coefficients[2,4],sign(summary(m2.f)$coefficients[2,1]),summary(m2.f)$coefficients[2,1],summary(m2.f)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"Men","m2",species,summary(m2.m)$coefficients[2,4],sign(summary(m2.m)$coefficients[2,1]),summary(m2.m)$coefficients[2,1],summary(m2.m)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"All","m3",species,summary(m3)$coefficients[2,4],sign(summary(m3)$coefficients[2,1]),summary(m3)$coefficients[2,1],summary(m3)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"Women","m3",species,summary(m3.f)$coefficients[2,4],sign(summary(m3.f)$coefficients[2,1]),summary(m3.f)$coefficients[2,1],summary(m3.f)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"Men","m3",species,summary(m3.m)$coefficients[2,4],sign(summary(m3.m)$coefficients[2,1]),summary(m3.m)$coefficients[2,1],summary(m3.m)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"All","m4",species,summary(m4)$coefficients[2,4],sign(summary(m4)$coefficients[2,1]),summary(m4)$coefficients[2,1],summary(m4)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"Women","m4",species,summary(m4.f)$coefficients[2,4],sign(summary(m4.f)$coefficients[2,1]),summary(m4.f)$coefficients[2,1],summary(m4.f)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"Men","m4",species,summary(m4.m)$coefficients[2,4],sign(summary(m4.m)$coefficients[2,1]),summary(m4.m)$coefficients[2,1],summary(m4.m)$coefficients[2,2])
  }
  cor.matrix$pval <- as.numeric(cor.matrix$pval)
  filt.feat <- cor.matrix[cor.matrix$pval<=0.05,"species"]
  if(length(filt.feat)==0){next}
  cor.matrix <- cor.matrix[cor.matrix[,"species"]%in%filt.feat,]
  cor.matrix$sym <- ifelse(cor.matrix$pval>=0.05, "", ifelse(cor.matrix$pval<=0.01,ifelse(cor.matrix$pval<=0.001,"***", "**"), "*"))
  cor.matrix$log <- -log10(cor.matrix$pval)*as.numeric(cor.matrix$direction)
  cor.matrix$species <- gsub("_"," ",cor.matrix$species)
  cor.matrix$species <- gsub("\\."," ",cor.matrix$species)
  cor.matrix$feature <- gsub("\\."," ",cor.matrix$feature)
  cor.matrix$feature <- gsub("X","",cor.matrix$feature)
  cor.matrix$feature <- gsub("_","-",cor.matrix$feature)
  title <- gsub("\\."," ",feat)
  title <- gsub("X","",title)
  title <- gsub("_","-",title)
  cor.matrix$gender <- factor(cor.matrix$gender, levels = c("All","Women","Men"))
  print(ggplot(cor.matrix, aes(x=model, y = species, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(name = "-log10", midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
          coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender)  +
          theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"),
                legend.text = element_text(size = 12, face="bold"), legend.title = element_text(size = 12, face="bold")) + ggtitle(title))
  cor.matrix.full <- rbind(cor.matrix.full,cor.matrix)
}
for(feat in c(DBP_cov,DBP_sd)){
  cor.matrix <-  data.frame(feature = character(), gender = character(), model = character(), species = character(), pval = character(), direction = character(),coefficient = character(), se = character())
  for(species in acids_variables){
    m1 <- lm_model(merged.dat,feat,species)
    m1.f <- lm_model(merged.dat[female,],feat,species)
    m1.m <- lm_model(merged.dat[male,],feat,species)
    m2 <- lm_model(merged.dat,feat,c(species,m2_covar,'Sex'))
    m2.f <- lm_model(merged.dat[female,],feat,c(species,m2_covar))
    m2.m <- lm_model(merged.dat[male,],feat,c(species,m2_covar))
    m3 <- lm_model(merged.dat,feat,c(species,m3_covar_dbp,'Sex'))
    m3.f <- lm_model(merged.dat[female,],feat,c(species,m3_covar_dbp,'Menopause_code'))
    m3.m <- lm_model(merged.dat[male,],feat,c(species,m3_covar_dbp))
    m4 <- lm_model(merged.dat,feat,c(species,m4_covar_dbp,'Sex'))
    m4.f <- lm_model(merged.dat[female,],feat,c(species,m4_covar_dbp,'Menopause_code'))
    m4.m <- lm_model(merged.dat[male,],feat,c(species,m4_covar_dbp))
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"All","m1",species,summary(m1)$coefficients[2,4],sign(summary(m1)$coefficients[2,1]),summary(m1)$coefficients[2,1],summary(m1)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"Women","m1",species,summary(m1.f)$coefficients[2,4],sign(summary(m1.f)$coefficients[2,1]),summary(m1.f)$coefficients[2,1],summary(m1.f)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"Men","m1",species,summary(m1.m)$coefficients[2,4],sign(summary(m1.m)$coefficients[2,1]),summary(m1.m)$coefficients[2,1],summary(m1.m)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"All","m2",species,summary(m2)$coefficients[2,4],sign(summary(m2)$coefficients[2,1]),summary(m2)$coefficients[2,1],summary(m2)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"Women","m2",species,summary(m2.f)$coefficients[2,4],sign(summary(m2.f)$coefficients[2,1]),summary(m2.f)$coefficients[2,1],summary(m2.f)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"Men","m2",species,summary(m2.m)$coefficients[2,4],sign(summary(m2.m)$coefficients[2,1]),summary(m2.m)$coefficients[2,1],summary(m2.m)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"All","m3",species,summary(m3)$coefficients[2,4],sign(summary(m3)$coefficients[2,1]),summary(m3)$coefficients[2,1],summary(m3)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"Women","m3",species,summary(m3.f)$coefficients[2,4],sign(summary(m3.f)$coefficients[2,1]),summary(m3.f)$coefficients[2,1],summary(m3.f)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"Men","m3",species,summary(m3.m)$coefficients[2,4],sign(summary(m3.m)$coefficients[2,1]),summary(m3.m)$coefficients[2,1],summary(m3.m)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"All","m4",species,summary(m4)$coefficients[2,4],sign(summary(m4)$coefficients[2,1]),summary(m4)$coefficients[2,1],summary(m4)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"Women","m4",species,summary(m4.f)$coefficients[2,4],sign(summary(m4.f)$coefficients[2,1]),summary(m4.f)$coefficients[2,1],summary(m4.f)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"Men","m4",species,summary(m4.m)$coefficients[2,4],sign(summary(m4.m)$coefficients[2,1]),summary(m4.m)$coefficients[2,1],summary(m4.m)$coefficients[2,2])
  }
  cor.matrix$pval <- as.numeric(cor.matrix$pval)
  filt.feat <- cor.matrix[cor.matrix$pval<=0.05,"species"]
  if(length(filt.feat)==0){next}
  cor.matrix <- cor.matrix[cor.matrix[,"species"]%in%filt.feat,]
  cor.matrix$sym <- ifelse(cor.matrix$pval>=0.05, "", ifelse(cor.matrix$pval<=0.01,ifelse(cor.matrix$pval<=0.001,"***", "**"), "*"))
  cor.matrix$log <- -log10(cor.matrix$pval)*as.numeric(cor.matrix$direction)
  cor.matrix$species <- gsub("_"," ",cor.matrix$species)
  cor.matrix$species <- gsub("\\."," ",cor.matrix$species)
  cor.matrix$feature <- gsub("\\."," ",cor.matrix$feature)
  cor.matrix$feature <- gsub("X","",cor.matrix$feature)
  cor.matrix$feature <- gsub("_","-",cor.matrix$feature)
  title <- gsub("\\."," ",feat)
  title <- gsub("X","",title)
  title <- gsub("_","-",title)
  cor.matrix$gender <- factor(cor.matrix$gender, levels = c("All","Women","Men"))
  print(ggplot(cor.matrix, aes(x=model, y = species, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(name = "-log10", midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
          coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender)  +
          theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"),
                legend.text = element_text(size = 12, face="bold"), legend.title = element_text(size = 12, face="bold")) + ggtitle(title))
  cor.matrix.full <- rbind(cor.matrix.full,cor.matrix)
}
dev.off()
write.csv(cor.matrix.full[,c('feature','gender','model','species','pval','coefficient','se')], "Figures/300_bp_variability/BP_variability_vs_blood_SCFAs.csv",row.names = F)
}

# Correlation analysis w/ Stool SCFAs
{
#BP Variability
pdf("Figures/300_bp_variability/BP_variability_vs_stool_SCFA.pdf", width = 8, height = 5)
cor.matrix.full <-  data.frame(variable = character(), gender = character(), model = character(), feature = character(), pval = character(), direction = character(),coefficient = character(), se = character())
for(feat in c(SBP_cov,SBP_sd,SBP_surge)){
  cor.matrix <-  data.frame(variable = character(), gender = character(), model = character(), feature = character(), pval = character(), direction = character(),coefficient = character(), se = character())
  for(species in stool_scfa_variables){
    m1 <- lm_model(merged.dat,feat,species)
    m1.f <- lm_model(merged.dat[female,],feat,species)
    m1.m <- lm_model(merged.dat[male,],feat,species)
    m2 <- lm_model(merged.dat,feat,c(species,m2_covar,'Sex'))
    m2.f <- lm_model(merged.dat[female,],feat,c(species,m2_covar))
    m2.m <- lm_model(merged.dat[male,],feat,c(species,m2_covar))
    m3 <- lm_model(merged.dat,feat,c(species,m3_covar_sbp,'Sex'))
    m3.f <- lm_model(merged.dat[female,],feat,c(species,m3_covar_sbp,'Menopause_code'))
    m3.m <- lm_model(merged.dat[male,],feat,c(species,m3_covar_sbp))
    m4 <- lm_model(merged.dat,feat,c(species,m4_covar_sbp,'Sex'))
    m4.f <- lm_model(merged.dat[female,],feat,c(species,m4_covar_sbp,'Menopause_code'))
    m4.m <- lm_model(merged.dat[male,],feat,c(species,m4_covar_sbp))
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"All","m1",species,summary(m1)$coefficients[2,4],sign(summary(m1)$coefficients[2,1]),summary(m1)$coefficients[2,1],summary(m1)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"Women","m1",species,summary(m1.f)$coefficients[2,4],sign(summary(m1.f)$coefficients[2,1]),summary(m1.f)$coefficients[2,1],summary(m1.f)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"Men","m1",species,summary(m1.m)$coefficients[2,4],sign(summary(m1.m)$coefficients[2,1]),summary(m1.m)$coefficients[2,1],summary(m1.m)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"All","m2",species,summary(m2)$coefficients[2,4],sign(summary(m2)$coefficients[2,1]),summary(m2)$coefficients[2,1],summary(m2)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"Women","m2",species,summary(m2.f)$coefficients[2,4],sign(summary(m2.f)$coefficients[2,1]),summary(m2.f)$coefficients[2,1],summary(m2.f)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"Men","m2",species,summary(m2.m)$coefficients[2,4],sign(summary(m2.m)$coefficients[2,1]),summary(m2.m)$coefficients[2,1],summary(m2.m)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"All","m3",species,summary(m3)$coefficients[2,4],sign(summary(m3)$coefficients[2,1]),summary(m3)$coefficients[2,1],summary(m3)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"Women","m3",species,summary(m3.f)$coefficients[2,4],sign(summary(m3.f)$coefficients[2,1]),summary(m3.f)$coefficients[2,1],summary(m3.f)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"Men","m3",species,summary(m3.m)$coefficients[2,4],sign(summary(m3.m)$coefficients[2,1]),summary(m3.m)$coefficients[2,1],summary(m3.m)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"All","m4",species,summary(m4)$coefficients[2,4],sign(summary(m4)$coefficients[2,1]),summary(m4)$coefficients[2,1],summary(m4)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"Women","m4",species,summary(m4.f)$coefficients[2,4],sign(summary(m4.f)$coefficients[2,1]),summary(m4.f)$coefficients[2,1],summary(m4.f)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"Men","m4",species,summary(m4.m)$coefficients[2,4],sign(summary(m4.m)$coefficients[2,1]),summary(m4.m)$coefficients[2,1],summary(m4.m)$coefficients[2,2])
  }
  cor.matrix$pval <- as.numeric(cor.matrix$pval)
  filt.feat <- cor.matrix[cor.matrix$pval<=0.05,"feature"]
  if(length(filt.feat)==0){next}
  cor.matrix <- cor.matrix[cor.matrix[,"feature"]%in%filt.feat,]
  cor.matrix$sym <- ifelse(cor.matrix$pval>=0.05, "", ifelse(cor.matrix$pval<=0.01,ifelse(cor.matrix$pval<=0.001,"***", "**"), "*"))
  cor.matrix$log <- -log10(cor.matrix$pval)*as.numeric(cor.matrix$direction)
  cor.matrix$feature <- gsub("_"," ",cor.matrix$feature)
  title <- gsub("\\."," ",feat)
  title <- gsub("X","",title)
  title <- gsub("_","-",title)
  cor.matrix$gender <- factor(cor.matrix$gender, levels = c("All","Women","Men"))
  print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(name = "-log10", midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
          coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender)  +
          theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"),
                legend.text = element_text(size = 12, face="bold"), legend.title = element_text(size = 12, face="bold")) + ggtitle(title))
  cor.matrix.full <- rbind(cor.matrix.full,cor.matrix)
}
for(feat in c(DBP_cov,DBP_sd)){
  cor.matrix <-  data.frame(variable = character(), gender = character(), model = character(), feature = character(), pval = character(), direction = character(),coefficient = character(), se = character())
  for(species in stool_scfa_variables){
    m1 <- lm_model(merged.dat,feat,species)
    m1.f <- lm_model(merged.dat[female,],feat,species)
    m1.m <- lm_model(merged.dat[male,],feat,species)
    m2 <- lm_model(merged.dat,feat,c(species,m2_covar,'Sex'))
    m2.f <- lm_model(merged.dat[female,],feat,c(species,m2_covar))
    m2.m <- lm_model(merged.dat[male,],feat,c(species,m2_covar))
    m3 <- lm_model(merged.dat,feat,c(species,m3_covar_dbp,'Sex'))
    m3.f <- lm_model(merged.dat[female,],feat,c(species,m3_covar_dbp,'Menopause_code'))
    m3.m <- lm_model(merged.dat[male,],feat,c(species,m3_covar_dbp))
    m4 <- lm_model(merged.dat,feat,c(species,m4_covar_dbp,'Sex'))
    m4.f <- lm_model(merged.dat[female,],feat,c(species,m4_covar_dbp,'Menopause_code'))
    m4.m <- lm_model(merged.dat[male,],feat,c(species,m4_covar_dbp))
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"All","m1",species,summary(m1)$coefficients[2,4],sign(summary(m1)$coefficients[2,1]),summary(m1)$coefficients[2,1],summary(m1)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"Women","m1",species,summary(m1.f)$coefficients[2,4],sign(summary(m1.f)$coefficients[2,1]),summary(m1.f)$coefficients[2,1],summary(m1.f)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"Men","m1",species,summary(m1.m)$coefficients[2,4],sign(summary(m1.m)$coefficients[2,1]),summary(m1.m)$coefficients[2,1],summary(m1.m)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"All","m2",species,summary(m2)$coefficients[2,4],sign(summary(m2)$coefficients[2,1]),summary(m2)$coefficients[2,1],summary(m2)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"Women","m2",species,summary(m2.f)$coefficients[2,4],sign(summary(m2.f)$coefficients[2,1]),summary(m2.f)$coefficients[2,1],summary(m2.f)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"Men","m2",species,summary(m2.m)$coefficients[2,4],sign(summary(m2.m)$coefficients[2,1]),summary(m2.m)$coefficients[2,1],summary(m2.m)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"All","m3",species,summary(m3)$coefficients[2,4],sign(summary(m3)$coefficients[2,1]),summary(m3)$coefficients[2,1],summary(m3)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"Women","m3",species,summary(m3.f)$coefficients[2,4],sign(summary(m3.f)$coefficients[2,1]),summary(m3.f)$coefficients[2,1],summary(m3.f)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"Men","m3",species,summary(m3.m)$coefficients[2,4],sign(summary(m3.m)$coefficients[2,1]),summary(m3.m)$coefficients[2,1],summary(m3.m)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"All","m4",species,summary(m4)$coefficients[2,4],sign(summary(m4)$coefficients[2,1]),summary(m4)$coefficients[2,1],summary(m4)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"Women","m4",species,summary(m4.f)$coefficients[2,4],sign(summary(m4.f)$coefficients[2,1]),summary(m4.f)$coefficients[2,1],summary(m4.f)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"Men","m4",species,summary(m4.m)$coefficients[2,4],sign(summary(m4.m)$coefficients[2,1]),summary(m4.m)$coefficients[2,1],summary(m4.m)$coefficients[2,2])
  }
  cor.matrix$pval <- as.numeric(cor.matrix$pval)
  filt.feat <- cor.matrix[cor.matrix$pval<=0.05,"feature"]
  if(length(filt.feat)==0){next}
  cor.matrix <- cor.matrix[cor.matrix[,"feature"]%in%filt.feat,]
  cor.matrix$sym <- ifelse(cor.matrix$pval>=0.05, "", ifelse(cor.matrix$pval<=0.01,ifelse(cor.matrix$pval<=0.001,"***", "**"), "*"))
  cor.matrix$log <- -log10(cor.matrix$pval)*as.numeric(cor.matrix$direction)
  cor.matrix$feature <- gsub("_"," ",cor.matrix$feature)
  title <- gsub("\\."," ",feat)
  title <- gsub("X","",title)
  title <- gsub("_","-",title)
  cor.matrix$gender <- factor(cor.matrix$gender, levels = c("All","Women","Men"))
  print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(name = "-log10", midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
          coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender)  +
          theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"),
                legend.text = element_text(size = 12, face="bold"), legend.title = element_text(size = 12, face="bold")) + ggtitle(title))
  cor.matrix.full <- rbind(cor.matrix.full,cor.matrix)
}
dev.off()
write.csv(cor.matrix.full[,c('variable','gender','model','feature','pval','coefficient','se')], "Figures/300_bp_variability/BP_variability_vs_stool_SCFAs.csv",row.names = F)


#GM
pdf("Figures/300_bp_variability/Stool_SCFA_vs_GM.pdf", width = 8, height = 5)
for(feat in stool_scfa_variables){
  cor.matrix <-  data.frame(gender = character(), model = character(), feature = character(), pval = character(), direction = character())
  for(species in colnames(transformed.species.data)){
    m1 <- lm_model(merged.dat,feat,species)
    m1.f <- lm_model(merged.dat[female,],feat,species)
    m1.m <- lm_model(merged.dat[male,],feat,species)
    m2 <- lm_model(merged.dat,feat,c(species,m2_covar,'Sex'))
    m2.f <- lm_model(merged.dat[female,],feat,c(species,m2_covar))
    m2.m <- lm_model(merged.dat[male,],feat,c(species,m2_covar))
    m3 <- lm_model(merged.dat,feat,c(species,m3_covar,'Sex'))
    m3.f <- lm_model(merged.dat[female,],feat,c(species,m3_covar,'Menopause_code'))
    m3.m <- lm_model(merged.dat[male,],feat,c(species,m3_covar))
    m4 <- lm_model(merged.dat,feat,c(species,m4_covar,'Sex'))
    m4.f <- lm_model(merged.dat[female,],feat,c(species,m4_covar,'Menopause_code'))
    m4.m <- lm_model(merged.dat[male,],feat,c(species,m4_covar))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m1",species,summary(m1)$coefficients[2,4],sign(summary(m1)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m1",species,summary(m1.f)$coefficients[2,4],sign(summary(m1.f)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m1",species,summary(m1.m)$coefficients[2,4],sign(summary(m1.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m2",species,summary(m2)$coefficients[2,4],sign(summary(m2)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m2",species,summary(m2.f)$coefficients[2,4],sign(summary(m2.f)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m2",species,summary(m2.m)$coefficients[2,4],sign(summary(m2.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3",species,summary(m3)$coefficients[2,4],sign(summary(m3)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3",species,summary(m3.f)$coefficients[2,4],sign(summary(m3.f)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3",species,summary(m3.m)$coefficients[2,4],sign(summary(m3.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m4",species,summary(m4)$coefficients[2,4],sign(summary(m4)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m4",species,summary(m4.f)$coefficients[2,4],sign(summary(m4.f)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m4",species,summary(m4.m)$coefficients[2,4],sign(summary(m4.m)$coefficients[2,1]))
  }
  cor.matrix$pval <- as.numeric(cor.matrix$pval)
  filt.feat <- cor.matrix[cor.matrix$pval<=0.01,"feature"]
  if(length(filt.feat)==0){next}
  cor.matrix <- cor.matrix[cor.matrix[,"feature"]%in%filt.feat,]
  cor.matrix$sym <- ifelse(cor.matrix$pval>=0.05, "", ifelse(cor.matrix$pval<=0.01,ifelse(cor.matrix$pval<=0.001,"***", "**"), "*"))
  cor.matrix$log <- -log10(cor.matrix$pval)*as.numeric(cor.matrix$direction)
  cor.matrix$feature <- gsub("_"," ",cor.matrix$feature)
  title <- gsub("\\."," ",feat)
  title <- gsub("X","",title)
  title <- gsub("_","-",title)
  cor.matrix$gender <- factor(cor.matrix$gender, levels = c("All","Women","Men"))
  print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(name = "-log10", midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
          coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender)  +
          theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"),
                legend.text = element_text(size = 12, face="bold"), legend.title = element_text(size = 12, face="bold")) + ggtitle(title))
}
dev.off()
}

# Correlation analysis w/ Diet
{
pdf("Figures/300_bp_variability/BP_variability_vs_diet.pdf", width = 8, height = 5)
cor.matrix.full <-  data.frame(variable = character(), gender = character(), model = character(), feature = character(), pval = character(), direction = character(),coefficient = character(), se = character())
for(feat in c(SBP_cov,SBP_sd,SBP_surge)){
  cor.matrix <-  data.frame(variable = character(), gender = character(), model = character(), feature = character(), pval = character(), direction = character(),coefficient = character(), se = character())
  for(species in macro_micro_nutrients){
    m1 <- lm_model(merged.dat,feat,species)
    m1.f <- lm_model(merged.dat[female,],feat,species)
    m1.m <- lm_model(merged.dat[male,],feat,species)
    m2 <- lm_model(merged.dat,feat,c(species,m2_covar,'Sex'))
    m2.f <- lm_model(merged.dat[female,],feat,c(species,m2_covar))
    m2.m <- lm_model(merged.dat[male,],feat,c(species,m2_covar))
    m3 <- lm_model(merged.dat,feat,c(species,m3_covar_sbp,'Sex'))
    m3.f <- lm_model(merged.dat[female,],feat,c(species,m3_covar_sbp,'Menopause_code'))
    m3.m <- lm_model(merged.dat[male,],feat,c(species,m3_covar_sbp))
    m4 <- lm_model(merged.dat,feat,c(species,m4_covar_sbp,'Sex'))
    m4.f <- lm_model(merged.dat[female,],feat,c(species,m4_covar_sbp,'Menopause_code'))
    m4.m <- lm_model(merged.dat[male,],feat,c(species,m4_covar_sbp))
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"All","m1",species,summary(m1)$coefficients[2,4],sign(summary(m1)$coefficients[2,1]),summary(m1)$coefficients[2,1],summary(m1)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"Women","m1",species,summary(m1.f)$coefficients[2,4],sign(summary(m1.f)$coefficients[2,1]),summary(m1.f)$coefficients[2,1],summary(m1.f)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"Men","m1",species,summary(m1.m)$coefficients[2,4],sign(summary(m1.m)$coefficients[2,1]),summary(m1.m)$coefficients[2,1],summary(m1.m)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"All","m2",species,summary(m2)$coefficients[2,4],sign(summary(m2)$coefficients[2,1]),summary(m2)$coefficients[2,1],summary(m2)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"Women","m2",species,summary(m2.f)$coefficients[2,4],sign(summary(m2.f)$coefficients[2,1]),summary(m2.f)$coefficients[2,1],summary(m2.f)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"Men","m2",species,summary(m2.m)$coefficients[2,4],sign(summary(m2.m)$coefficients[2,1]),summary(m2.m)$coefficients[2,1],summary(m2.m)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"All","m3",species,summary(m3)$coefficients[2,4],sign(summary(m3)$coefficients[2,1]),summary(m3)$coefficients[2,1],summary(m3)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"Women","m3",species,summary(m3.f)$coefficients[2,4],sign(summary(m3.f)$coefficients[2,1]),summary(m3.f)$coefficients[2,1],summary(m3.f)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"Men","m3",species,summary(m3.m)$coefficients[2,4],sign(summary(m3.m)$coefficients[2,1]),summary(m3.m)$coefficients[2,1],summary(m3.m)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"All","m4",species,summary(m4)$coefficients[2,4],sign(summary(m4)$coefficients[2,1]),summary(m4)$coefficients[2,1],summary(m4)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"Women","m4",species,summary(m4.f)$coefficients[2,4],sign(summary(m4.f)$coefficients[2,1]),summary(m4.f)$coefficients[2,1],summary(m4.f)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"Men","m4",species,summary(m4.m)$coefficients[2,4],sign(summary(m4.m)$coefficients[2,1]),summary(m4.m)$coefficients[2,1],summary(m4.m)$coefficients[2,2])
  }
  cor.matrix$pval <- as.numeric(cor.matrix$pval)
  filt.feat <- cor.matrix[cor.matrix$pval<=0.05,"feature"]
  if(length(filt.feat)==0){next}
  cor.matrix <- cor.matrix[cor.matrix[,"feature"]%in%filt.feat,]
  cor.matrix$sym <- ifelse(cor.matrix$pval>=0.05, "", ifelse(cor.matrix$pval<=0.01,ifelse(cor.matrix$pval<=0.001,"***", "**"), "*"))
  cor.matrix$log <- -log10(cor.matrix$pval)*as.numeric(cor.matrix$direction)
  cor.matrix$feature <- gsub("_"," ",cor.matrix$feature)
  title <- gsub("\\."," ",feat)
  title <- gsub("X","",title)
  title <- gsub("_","-",title)
  cor.matrix$gender <- factor(cor.matrix$gender, levels = c("All","Women","Men"))
  print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(name = "-log10", midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
          coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender)  +
          theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"),
                legend.text = element_text(size = 12, face="bold"), legend.title = element_text(size = 12, face="bold")) + ggtitle(title))
  cor.matrix.full <- rbind(cor.matrix.full,cor.matrix)
}
for(feat in c(DBP_cov,DBP_sd)){
  cor.matrix <-  data.frame(variable = character(), gender = character(), model = character(), feature = character(), pval = character(), direction = character(),coefficient = character(), se = character())
  for(species in macro_micro_nutrients){
    m1 <- lm_model(merged.dat,feat,species)
    m1.f <- lm_model(merged.dat[female,],feat,species)
    m1.m <- lm_model(merged.dat[male,],feat,species)
    m2 <- lm_model(merged.dat,feat,c(species,m2_covar,'Sex'))
    m2.f <- lm_model(merged.dat[female,],feat,c(species,m2_covar))
    m2.m <- lm_model(merged.dat[male,],feat,c(species,m2_covar))
    m3 <- lm_model(merged.dat,feat,c(species,m3_covar_dbp,'Sex'))
    m3.f <- lm_model(merged.dat[female,],feat,c(species,m3_covar_dbp,'Menopause_code'))
    m3.m <- lm_model(merged.dat[male,],feat,c(species,m3_covar_dbp))
    m4 <- lm_model(merged.dat,feat,c(species,m4_covar_dbp,'Sex'))
    m4.f <- lm_model(merged.dat[female,],feat,c(species,m4_covar_dbp,'Menopause_code'))
    m4.m <- lm_model(merged.dat[male,],feat,c(species,m4_covar_dbp))
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"All","m1",species,summary(m1)$coefficients[2,4],sign(summary(m1)$coefficients[2,1]),summary(m1)$coefficients[2,1],summary(m1)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"Women","m1",species,summary(m1.f)$coefficients[2,4],sign(summary(m1.f)$coefficients[2,1]),summary(m1.f)$coefficients[2,1],summary(m1.f)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"Men","m1",species,summary(m1.m)$coefficients[2,4],sign(summary(m1.m)$coefficients[2,1]),summary(m1.m)$coefficients[2,1],summary(m1.m)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"All","m2",species,summary(m2)$coefficients[2,4],sign(summary(m2)$coefficients[2,1]),summary(m2)$coefficients[2,1],summary(m2)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"Women","m2",species,summary(m2.f)$coefficients[2,4],sign(summary(m2.f)$coefficients[2,1]),summary(m2.f)$coefficients[2,1],summary(m2.f)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"Men","m2",species,summary(m2.m)$coefficients[2,4],sign(summary(m2.m)$coefficients[2,1]),summary(m2.m)$coefficients[2,1],summary(m2.m)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"All","m3",species,summary(m3)$coefficients[2,4],sign(summary(m3)$coefficients[2,1]),summary(m3)$coefficients[2,1],summary(m3)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"Women","m3",species,summary(m3.f)$coefficients[2,4],sign(summary(m3.f)$coefficients[2,1]),summary(m3.f)$coefficients[2,1],summary(m3.f)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"Men","m3",species,summary(m3.m)$coefficients[2,4],sign(summary(m3.m)$coefficients[2,1]),summary(m3.m)$coefficients[2,1],summary(m3.m)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"All","m4",species,summary(m4)$coefficients[2,4],sign(summary(m4)$coefficients[2,1]),summary(m4)$coefficients[2,1],summary(m4)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"Women","m4",species,summary(m4.f)$coefficients[2,4],sign(summary(m4.f)$coefficients[2,1]),summary(m4.f)$coefficients[2,1],summary(m4.f)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"Men","m4",species,summary(m4.m)$coefficients[2,4],sign(summary(m4.m)$coefficients[2,1]),summary(m4.m)$coefficients[2,1],summary(m4.m)$coefficients[2,2])
  }
  cor.matrix$pval <- as.numeric(cor.matrix$pval)
  filt.feat <- cor.matrix[cor.matrix$pval<=0.05,"feature"]
  if(length(filt.feat)==0){next}
  cor.matrix <- cor.matrix[cor.matrix[,"feature"]%in%filt.feat,]
  cor.matrix$sym <- ifelse(cor.matrix$pval>=0.05, "", ifelse(cor.matrix$pval<=0.01,ifelse(cor.matrix$pval<=0.001,"***", "**"), "*"))
  cor.matrix$log <- -log10(cor.matrix$pval)*as.numeric(cor.matrix$direction)
  cor.matrix$feature <- gsub("_"," ",cor.matrix$feature)
  title <- gsub("\\."," ",feat)
  title <- gsub("X","",title)
  title <- gsub("_","-",title)
  cor.matrix$gender <- factor(cor.matrix$gender, levels = c("All","Women","Men"))
  print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(name = "-log10", midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
          coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender)  +
          theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"),
                legend.text = element_text(size = 12, face="bold"), legend.title = element_text(size = 12, face="bold")) + ggtitle(title))
  cor.matrix.full <- rbind(cor.matrix.full,cor.matrix)
}
dev.off()
write.csv(cor.matrix.full[,c('variable','gender','model','feature','pval','coefficient','se')], "Figures/300_bp_variability/BP_variability_vs_diet.csv",row.names = F)
}

#Correlation analysis w/ Probiotic species
{
pdf("Figures/300_bp_variability/Probiotic_vs_immune.pdf", width = 8, height = 5)
for(feat in immune_markers){
  cor.matrix <-  data.frame(gender = character(), model = character(), feature = character(), pval = character(), direction = character())
  for(species in probiotic_species){
    m1 <- lm_model(merged.dat,feat,species)
    m1.f <- lm_model(merged.dat[female,],feat,species)
    m1.m <- lm_model(merged.dat[male,],feat,species)
    m2 <- lm_model(merged.dat,feat,c(species,m2_covar,'Sex'))
    m2.f <- lm_model(merged.dat[female,],feat,c(species,m2_covar))
    m2.m <- lm_model(merged.dat[male,],feat,c(species,m2_covar))
    m3 <- lm_model(merged.dat,feat,c(species,m3_covar,'Sex'))
    m3.f <- lm_model(merged.dat[female,],feat,c(species,m3_covar,'Menopause_code'))
    m3.m <- lm_model(merged.dat[male,],feat,c(species,m3_covar))
    m4 <- lm_model(merged.dat,feat,c(species,m4_covar,'Sex'))
    m4.f <- lm_model(merged.dat[female,],feat,c(species,m4_covar,'Menopause_code'))
    m4.m <- lm_model(merged.dat[male,],feat,c(species,m4_covar))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m1",species,summary(m1)$coefficients[2,4],sign(summary(m1)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m1",species,summary(m1.f)$coefficients[2,4],sign(summary(m1.f)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m1",species,summary(m1.m)$coefficients[2,4],sign(summary(m1.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m2",species,summary(m2)$coefficients[2,4],sign(summary(m2)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m2",species,summary(m2.f)$coefficients[2,4],sign(summary(m2.f)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m2",species,summary(m2.m)$coefficients[2,4],sign(summary(m2.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3",species,summary(m3)$coefficients[2,4],sign(summary(m3)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3",species,summary(m3.f)$coefficients[2,4],sign(summary(m3.f)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3",species,summary(m3.m)$coefficients[2,4],sign(summary(m3.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m4",species,summary(m4)$coefficients[2,4],sign(summary(m4)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m4",species,summary(m4.f)$coefficients[2,4],sign(summary(m4.f)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m4",species,summary(m4.m)$coefficients[2,4],sign(summary(m4.m)$coefficients[2,1]))
  }
  cor.matrix$pval <- as.numeric(cor.matrix$pval)
  filt.feat <- cor.matrix[cor.matrix$pval<=0.05,"feature"]
  if(length(filt.feat)==0){next}
  cor.matrix <- cor.matrix[cor.matrix[,"feature"]%in%filt.feat,]
  cor.matrix$sym <- ifelse(cor.matrix$pval>=0.05, "", ifelse(cor.matrix$pval<=0.01,ifelse(cor.matrix$pval<=0.001,"***", "**"), "*"))
  cor.matrix$log <- -log10(cor.matrix$pval)*as.numeric(cor.matrix$direction)
  cor.matrix$feature <- gsub("_"," ",cor.matrix$feature)
  title <- gsub("\\."," ",feat)
  title <- gsub("X","",title)
  title <- gsub("_","-",title)
  cor.matrix$gender <- factor(cor.matrix$gender, levels = c("All","Women","Men"))
  print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(name = "-log10", midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
          coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender)  +
          theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"),
                legend.text = element_text(size = 12, face="bold"), legend.title = element_text(size = 12, face="bold")) + ggtitle(title))
}
dev.off()

pdf("Figures/300_bp_variability/Probiotic_vs_SCFA.pdf", width = 8, height = 5)
cor.matrix.full <-  data.frame(acid = character(), gender = character(), model = character(), species = character(), pval = character(), direction = character(), coefficient = character(), se = character())
for(feat in acids_variables){
  cor.matrix <-  data.frame(acid = character(), gender = character(), model = character(), species = character(), pval = character(), direction = character(), coefficient = character(), se = character())
  for(species in probiotic_species){
    m1 <- lm_model(merged.dat,feat,species)
    m1.f <- lm_model(merged.dat[female,],feat,species)
    m1.m <- lm_model(merged.dat[male,],feat,species)
    m2 <- lm_model(merged.dat,feat,c(species,m2_covar,'Sex'))
    m2.f <- lm_model(merged.dat[female,],feat,c(species,m2_covar))
    m2.m <- lm_model(merged.dat[male,],feat,c(species,m2_covar))
    m3 <- lm_model(merged.dat,feat,c(species,m3_covar,'Sex'))
    m3.f <- lm_model(merged.dat[female,],feat,c(species,m3_covar,'Menopause_code'))
    m3.m <- lm_model(merged.dat[male,],feat,c(species,m3_covar))
    m4 <- lm_model(merged.dat,feat,c(species,m4_covar,'Sex'))
    m4.f <- lm_model(merged.dat[female,],feat,c(species,m4_covar,'Menopause_code'))
    m4.m <- lm_model(merged.dat[male,],feat,c(species,m4_covar))
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"All","m1",species,summary(m1)$coefficients[2,4],sign(summary(m1)$coefficients[2,1]),summary(m1)$coefficients[2,1],summary(m1)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"Women","m1",species,summary(m1.f)$coefficients[2,4],sign(summary(m1.f)$coefficients[2,1]),summary(m1.f)$coefficients[2,1],summary(m1.f)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"Men","m1",species,summary(m1.m)$coefficients[2,4],sign(summary(m1.m)$coefficients[2,1]),summary(m1.m)$coefficients[2,1],summary(m1.m)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"All","m2",species,summary(m2)$coefficients[2,4],sign(summary(m2)$coefficients[2,1]),summary(m2)$coefficients[2,1],summary(m2)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"Women","m2",species,summary(m2.f)$coefficients[2,4],sign(summary(m2.f)$coefficients[2,1]),summary(m2.f)$coefficients[2,1],summary(m2.f)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"Men","m2",species,summary(m2.m)$coefficients[2,4],sign(summary(m2.m)$coefficients[2,1]),summary(m2.m)$coefficients[2,1],summary(m2.m)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"All","m3",species,summary(m3)$coefficients[2,4],sign(summary(m3)$coefficients[2,1]),summary(m3)$coefficients[2,1],summary(m3)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"Women","m3",species,summary(m3.f)$coefficients[2,4],sign(summary(m3.f)$coefficients[2,1]),summary(m3.f)$coefficients[2,1],summary(m3.f)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"Men","m3",species,summary(m3.m)$coefficients[2,4],sign(summary(m3.m)$coefficients[2,1]),summary(m3.m)$coefficients[2,1],summary(m3.m)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"All","m4",species,summary(m4)$coefficients[2,4],sign(summary(m4)$coefficients[2,1]),summary(m4)$coefficients[2,1],summary(m4)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"Women","m4",species,summary(m4.f)$coefficients[2,4],sign(summary(m4.f)$coefficients[2,1]),summary(m4.f)$coefficients[2,1],summary(m4.f)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"Men","m4",species,summary(m4.m)$coefficients[2,4],sign(summary(m4.m)$coefficients[2,1]),summary(m4.m)$coefficients[2,1],summary(m4.m)$coefficients[2,2])
  }
  cor.matrix$pval <- as.numeric(cor.matrix$pval)
  filt.feat <- cor.matrix[cor.matrix$pval<=0.05,"species"]
  if(length(filt.feat)==0){next}
  cor.matrix <- cor.matrix[cor.matrix[,"species"]%in%filt.feat,]
  cor.matrix$sym <- ifelse(cor.matrix$pval>=0.05, "", ifelse(cor.matrix$pval<=0.01,ifelse(cor.matrix$pval<=0.001,"***", "**"), "*"))
  cor.matrix$log <- -log10(cor.matrix$pval)*as.numeric(cor.matrix$direction)
  cor.matrix$species <- gsub("_"," ",cor.matrix$species)
  title <- gsub("\\."," ",feat)
  title <- gsub("X","",title)
  title <- gsub("_","-",title)
  cor.matrix$gender <- factor(cor.matrix$gender, levels = c("All","Women","Men"))
  cor.matrix.full <- rbind(cor.matrix.full,cor.matrix)
  print(ggplot(cor.matrix, aes(x=model, y = species, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(name = "-log10", midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
          coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender)  +
          theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"),
                legend.text = element_text(size = 12, face="bold"), legend.title = element_text(size = 12, face="bold")) + ggtitle(title))
}
dev.off()
cor.matrix.full$acid <- gsub("\\."," ",cor.matrix.full$acid)
cor.matrix.full$acid <- gsub("_"," ",cor.matrix.full$acid)
write.csv(cor.matrix.full[,c('acid','gender','model','species','pval','coefficient','se')], "Figures/300_bp_variability/Probiotic_vs_SCFA.csv",row.names = F)


pdf("Figures/300_bp_variability/Probiotic_vs_diet.pdf", width = 8, height = 5)
cor.matrix.full <-  data.frame(diet = character(), gender = character(), model = character(), species = character(), pval = character(), direction = character(), coefficient = character(), se = character())
for(feat in macro_micro_nutrients){
  cor.matrix <-  data.frame(diet = character(), gender = character(), model = character(), species = character(), pval = character(), direction = character(), coefficient = character(), se = character())
  for(species in probiotic_species){
    m1 <- lm_model(merged.dat,feat,species)
    m1.f <- lm_model(merged.dat[female,],feat,species)
    m1.m <- lm_model(merged.dat[male,],feat,species)
    m2 <- lm_model(merged.dat,feat,c(species,m2_covar,'Sex'))
    m2.f <- lm_model(merged.dat[female,],feat,c(species,m2_covar))
    m2.m <- lm_model(merged.dat[male,],feat,c(species,m2_covar))
    m3 <- lm_model(merged.dat,feat,c(species,m3_covar,'Sex'))
    m3.f <- lm_model(merged.dat[female,],feat,c(species,m3_covar,'Menopause_code'))
    m3.m <- lm_model(merged.dat[male,],feat,c(species,m3_covar))
    m4 <- lm_model(merged.dat,feat,c(species,m4_covar,'Sex'))
    m4.f <- lm_model(merged.dat[female,],feat,c(species,m4_covar,'Menopause_code'))
    m4.m <- lm_model(merged.dat[male,],feat,c(species,m4_covar))
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"All","m1",species,summary(m1)$coefficients[2,4],sign(summary(m1)$coefficients[2,1]),summary(m1)$coefficients[2,1],summary(m1)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"Women","m1",species,summary(m1.f)$coefficients[2,4],sign(summary(m1.f)$coefficients[2,1]),summary(m1.f)$coefficients[2,1],summary(m1.f)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"Men","m1",species,summary(m1.m)$coefficients[2,4],sign(summary(m1.m)$coefficients[2,1]),summary(m1.m)$coefficients[2,1],summary(m1.m)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"All","m2",species,summary(m2)$coefficients[2,4],sign(summary(m2)$coefficients[2,1]),summary(m2)$coefficients[2,1],summary(m2)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"Women","m2",species,summary(m2.f)$coefficients[2,4],sign(summary(m2.f)$coefficients[2,1]),summary(m2.f)$coefficients[2,1],summary(m2.f)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"Men","m2",species,summary(m2.m)$coefficients[2,4],sign(summary(m2.m)$coefficients[2,1]),summary(m2.m)$coefficients[2,1],summary(m2.m)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"All","m3",species,summary(m3)$coefficients[2,4],sign(summary(m3)$coefficients[2,1]),summary(m3)$coefficients[2,1],summary(m3)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"Women","m3",species,summary(m3.f)$coefficients[2,4],sign(summary(m3.f)$coefficients[2,1]),summary(m3.f)$coefficients[2,1],summary(m3.f)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"Men","m3",species,summary(m3.m)$coefficients[2,4],sign(summary(m3.m)$coefficients[2,1]),summary(m3.m)$coefficients[2,1],summary(m3.m)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"All","m4",species,summary(m4)$coefficients[2,4],sign(summary(m4)$coefficients[2,1]),summary(m4)$coefficients[2,1],summary(m4)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"Women","m4",species,summary(m4.f)$coefficients[2,4],sign(summary(m4.f)$coefficients[2,1]),summary(m4.f)$coefficients[2,1],summary(m4.f)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(feat,"Men","m4",species,summary(m4.m)$coefficients[2,4],sign(summary(m4.m)$coefficients[2,1]),summary(m4.m)$coefficients[2,1],summary(m4.m)$coefficients[2,2])
  }
  cor.matrix$pval <- as.numeric(cor.matrix$pval)
  filt.feat <- cor.matrix[cor.matrix$pval<=0.05,"species"]
  if(length(filt.feat)==0){next}
  cor.matrix <- cor.matrix[cor.matrix[,"species"]%in%filt.feat,]
  cor.matrix$sym <- ifelse(cor.matrix$pval>=0.05, "", ifelse(cor.matrix$pval<=0.01,ifelse(cor.matrix$pval<=0.001,"***", "**"), "*"))
  cor.matrix$log <- -log10(cor.matrix$pval)*as.numeric(cor.matrix$direction)
  cor.matrix$species <- gsub("_"," ",cor.matrix$species)
  title <- gsub("_"," ",feat)
  cor.matrix$gender <- factor(cor.matrix$gender, levels = c("All","Women","Men"))
  print(ggplot(cor.matrix, aes(x=model, y = species, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(name = "-log10", midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
          coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender)  +
          theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"),
                legend.text = element_text(size = 12, face="bold"), legend.title = element_text(size = 12, face="bold")) + ggtitle(title))
  cor.matrix.full <- rbind(cor.matrix.full,cor.matrix)
}
dev.off()
write.csv(cor.matrix.full[,c('diet','gender','model','species','pval','coefficient','se')], "Figures/300_bp_variability/Probiotic_vs_diet.csv",row.names = F)
}

#Correlation analysis w/ Sex
{
cor.matrix <- data.frame(model = character(), feature = character(), pval = character(), direction = character(),coefficient = character(), se = character())
for(feat in c(SBP_cov,DBP_cov)){
  m1 <- lm_model(merged.dat,feat,"Sex")
  m2 <- lm_model(merged.dat,feat,c("Sex",m2_covar))
  if(grepl("SBP",feat)){
    m3 <- lm_model(merged.dat,feat,c("Sex",m3_covar_sbp))
    m4 <- lm_model(merged.dat,feat,c("Sex",m4_covar_sbp))
  } else {
    m3 <- lm_model(merged.dat,feat,c("Sex",m3_covar_dbp))
    m4 <- lm_model(merged.dat,feat,c("Sex",m4_covar_dbp))
  }
  cor.matrix[nrow(cor.matrix)+1,] <- c("m1",feat,summary(m1)$coefficients[2,4],sign(summary(m1)$coefficients[2,1]),summary(m1)$coefficients[2,1],summary(m1)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("m2",feat,summary(m2)$coefficients[2,4],sign(summary(m2)$coefficients[2,1]),summary(m2)$coefficients[2,1],summary(m2)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("m3",feat,summary(m3)$coefficients[2,4],sign(summary(m3)$coefficients[2,1]),summary(m3)$coefficients[2,1],summary(m3)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("m4",feat,summary(m4)$coefficients[2,4],sign(summary(m4)$coefficients[2,1]),summary(m4)$coefficients[2,1],summary(m4)$coefficients[2,2])
}
cor.matrix$pval <- as.numeric(cor.matrix$pval)
cor.matrix$sym <- ifelse(cor.matrix$pval>=0.05, "", ifelse(cor.matrix$pval<=0.01,ifelse(cor.matrix$pval<=0.001,"***", "**"), "*"))
cor.matrix$log <- -log10(cor.matrix$pval)*as.numeric(cor.matrix$direction)
cor.matrix$feature <- gsub("_","-",cor.matrix$feature)
cor.matrix$feature <- gsub("X","",cor.matrix$feature)
cor.matrix$feature <- gsub("\\."," ",cor.matrix$feature)
cor.matrix$feature <- factor(cor.matrix$feature, levels = c("Nighttime DBP CoV","Nighttime SBP CoV","Daytime DBP CoV","Daytime SBP CoV","24-hour DBP CoV","24-hour SBP CoV"))
pdf("Figures/300_bp_variability/BP_variability_vs_sex.pdf", width = 8, height = 5)
print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(name = "-log10", midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
        coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black") +
        theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"),
              legend.text = element_text(size = 12, face="bold"), legend.title = element_text(size = 12, face="bold")) + ggtitle("Sex"))
dev.off()
write.csv(cor.matrix[,c('model','feature','pval','coefficient','se')], "Figures/300_bp_variability/BP_variability_vs_sex.csv",row.names = F)
}

#Correlation analysis: Bacteriodes Dorei vs. SCFAs
{
  cor.matrix <-  data.frame(gender = character(), model = character(), feature = character(), pval = character(), direction = character())
  for(feat in acids_variables){
    m1 <- lm_model(merged.dat,feat,"Bacteroides_dorei")
    m1.f <- lm_model(merged.dat[female,],feat,"Bacteroides_dorei")
    m1.m <- lm_model(merged.dat[male,],feat,"Bacteroides_dorei")
    m2 <- lm_model(merged.dat,feat,c("Bacteroides_dorei",m2_covar,'Sex'))
    m2.f <- lm_model(merged.dat[female,],feat,c("Bacteroides_dorei",m2_covar))
    m2.m <- lm_model(merged.dat[male,],feat,c("Bacteroides_dorei",m2_covar))
    m3 <- lm_model(merged.dat,feat,c("Bacteroides_dorei",m3_covar,'Sex'))
    m3.f <- lm_model(merged.dat[female,],feat,c("Bacteroides_dorei",m3_covar,'Menopause_code'))
    m3.m <- lm_model(merged.dat[male,],feat,c("Bacteroides_dorei",m3_covar))
    m4 <- lm_model(merged.dat,feat,c("Bacteroides_dorei",m4_covar,'Sex'))
    m4.f <- lm_model(merged.dat[female,],feat,c("Bacteroides_dorei",m4_covar,'Menopause_code'))
    m4.m <- lm_model(merged.dat[male,],feat,c("Bacteroides_dorei",m4_covar))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m1",feat,summary(m1)$coefficients[2,4],sign(summary(m1)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m1",feat,summary(m1.f)$coefficients[2,4],sign(summary(m1.f)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m1",feat,summary(m1.m)$coefficients[2,4],sign(summary(m1.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m2",feat,summary(m2)$coefficients[2,4],sign(summary(m2)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m2",feat,summary(m2.f)$coefficients[2,4],sign(summary(m2.f)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m2",feat,summary(m2.m)$coefficients[2,4],sign(summary(m2.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3",feat,summary(m3)$coefficients[2,4],sign(summary(m3)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3",feat,summary(m3.f)$coefficients[2,4],sign(summary(m3.f)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3",feat,summary(m3.m)$coefficients[2,4],sign(summary(m3.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m4",feat,summary(m4)$coefficients[2,4],sign(summary(m4)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m4",feat,summary(m4.f)$coefficients[2,4],sign(summary(m4.f)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m4",feat,summary(m4.m)$coefficients[2,4],sign(summary(m4.m)$coefficients[2,1]))
  }
  cor.matrix$pval <- as.numeric(cor.matrix$pval)
  #filt.feat <- cor.matrix[cor.matrix$pval<=0.01,"feature"]
  #cor.matrix <- cor.matrix[cor.matrix[,"feature"]%in%filt.feat,]
  cor.matrix$sym <- ifelse(cor.matrix$pval>=0.05, "", ifelse(cor.matrix$pval<=0.01,ifelse(cor.matrix$pval<=0.001,"***", "**"), "*"))
  cor.matrix$log <- -log10(cor.matrix$pval)*as.numeric(cor.matrix$direction)
  cor.matrix$feature <- gsub("_"," ",cor.matrix$feature)
  cor.matrix$feature <- gsub("X","",cor.matrix$feature)
  cor.matrix$feature <- gsub("\\."," ",cor.matrix$feature)
  cor.matrix$gender <- factor(cor.matrix$gender, levels = c("All","Women","Men"))
  pdf("Figures/300_bp_variability/Bacteroides_dorei_vs_SCFAs.pdf", width = 8, height = 5)
  print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(name = "-log10", midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
          coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black") + facet_wrap(~gender) +
          theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"),
                legend.text = element_text(size = 12, face="bold"), legend.title = element_text(size = 12, face="bold")) + ggtitle("Bacteroides Dorei"))
  dev.off()
}

#Correlation analysis: Blood/Stool SCFAs vs. Dipping Status
{
cor.matrix <-  data.frame(group = character(), gender = character(), model = character(), feature = character(), pval = character(), direction = character(), coefficient = character(), se = character())
for(feat in c(acids_variables,stool_scfa_variables)){
  for(group in c(1,3,4)){
    dat <- merged.dat[merged.dat$Dipping.classification%in%c(group,2),]
    dat$group <- ifelse(dat$Dipping.classification==2,0,1)
    m1 <- glm_model(dat,"group",feat)
    m1.f <- glm_model(dat[female,],"group",feat)
    m1.m <- glm_model(dat[male,],"group",feat)
    m2 <- glm_model(dat,"group",c(feat,m2_covar,'Sex'))
    m2.f <- glm_model(dat[female,],"group",c(feat,m2_covar))
    m2.m <- glm_model(dat[male,],"group",c(feat,m2_covar))
    m3 <- glm_model(dat,"group",c(feat,m3_covar_sbp,'Sex'))
    m3.f <- glm_model(dat[female,],"group",c(feat,m3_covar_sbp,'Menopause_code'))
    m3.m <- glm_model(dat[male,],"group",c(feat,m3_covar_sbp))
    m4 <- glm_model(dat,"group",c(feat,m4_covar_sbp,'Sex'))
    m4.f <- glm_model(dat[female,],"group",c(feat,m4_covar_sbp,'Menopause_code'))
    m4.m <- glm_model(dat[male,],"group",c(feat,m4_covar_sbp))
    cor.matrix[nrow(cor.matrix)+1,] <- c(group,"All","m1",feat,summary(m1)$coefficients[2,4],sign(summary(m1)$coefficients[2,1]),summary(m1)$coefficients[2,1],summary(m1)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(group,"Women","m1",feat,summary(m1.f)$coefficients[2,4],sign(summary(m1.f)$coefficients[2,1]),summary(m1.f)$coefficients[2,1],summary(m1.f)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(group,"Men","m1",feat,summary(m1.m)$coefficients[2,4],sign(summary(m1.m)$coefficients[2,1]),summary(m1.m)$coefficients[2,1],summary(m1.m)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(group,"All","m2",feat,summary(m2)$coefficients[2,4],sign(summary(m2)$coefficients[2,1]),summary(m2)$coefficients[2,1],summary(m2)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(group,"Women","m2",feat,summary(m2.f)$coefficients[2,4],sign(summary(m2.f)$coefficients[2,1]),summary(m2.f)$coefficients[2,1],summary(m2.f)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(group,"Men","m2",feat,summary(m2.m)$coefficients[2,4],sign(summary(m2.m)$coefficients[2,1]),summary(m2.m)$coefficients[2,1],summary(m2.m)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(group,"All","m3",feat,summary(m3)$coefficients[2,4],sign(summary(m3)$coefficients[2,1]),summary(m3)$coefficients[2,1],summary(m3)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(group,"Women","m3",feat,summary(m3.f)$coefficients[2,4],sign(summary(m3.f)$coefficients[2,1]),summary(m3.f)$coefficients[2,1],summary(m3.f)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(group,"Men","m3",feat,summary(m3.m)$coefficients[2,4],sign(summary(m3.m)$coefficients[2,1]),summary(m3.m)$coefficients[2,1],summary(m3.m)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(group,"All","m4",feat,summary(m4)$coefficients[2,4],sign(summary(m4)$coefficients[2,1]),summary(m4)$coefficients[2,1],summary(m4)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(group,"Women","m4",feat,summary(m4.f)$coefficients[2,4],sign(summary(m4.f)$coefficients[2,1]),summary(m4.f)$coefficients[2,1],summary(m4.f)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(group,"Men","m4",feat,summary(m4.m)$coefficients[2,4],sign(summary(m4.m)$coefficients[2,1]),summary(m4.m)$coefficients[2,1],summary(m4.m)$coefficients[2,2])
  }
}
cor.matrix$pval <- as.numeric(cor.matrix$pval)
filt.feat <- cor.matrix[cor.matrix$pval<=0.01,"feature"]
cor.matrix <- cor.matrix[cor.matrix[,"feature"]%in%filt.feat,]
filt.feat <- cor.matrix[cor.matrix$pval<=0.01,"group"]
cor.matrix <- cor.matrix[cor.matrix[,"group"]%in%filt.feat,]
cor.matrix$sym <- ifelse(cor.matrix$pval>=0.05, "", ifelse(cor.matrix$pval<=0.01,ifelse(cor.matrix$pval<=0.001,"***", "**"), "*"))
cor.matrix$log <- -log10(cor.matrix$pval)*as.numeric(cor.matrix$direction)
cor.matrix$feature <- gsub("_"," ",cor.matrix$feature)
cor.matrix$feature <- gsub("X","",cor.matrix$feature)
cor.matrix$feature <- gsub("\\."," ",cor.matrix$feature)
cor.matrix$gender <- factor(cor.matrix$gender, levels = c("All","Women","Men"))
pdf("Figures/300_bp_variability/Dipping_classification_vs_SCFAs.pdf", width = 8, height = 5)
print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(name = "-log10", midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
        coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black") + facet_wrap(~gender+group) +
        theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"),
              legend.text = element_text(size = 12, face="bold"), legend.title = element_text(size = 12, face="bold")) + ggtitle("Dipping classification"))
dev.off()
write.csv(cor.matrix[,c('group','gender','model','feature','pval','coefficient','se')], "Figures/300_bp_variability/Dipping_classification_vs_SCFAs.csv",row.names = F)
}

#Correlation analysis: GM vs. Dipping Status
{
cor.matrix <-  data.frame(group = character(), gender = character(), model = character(), feature = character(), pval = character(), direction = character(), coefficient = character(), se = character())
for(feat in colnames(transformed.species.data)){
  for(group in c(1,3,4)){
    dat <- merged.dat[merged.dat$Dipping.classification%in%c(group,2),]
    dat$group <- ifelse(dat$Dipping.classification==2,0,1)
    m1 <- glm_model(dat,"group",feat)
    m1.f <- glm_model(dat[female,],"group",feat)
    m1.m <- glm_model(dat[male,],"group",feat)
    m2 <- glm_model(dat,"group",c(feat,m2_covar,'Sex'))
    m2.f <- glm_model(dat[female,],"group",c(feat,m2_covar))
    m2.m <- glm_model(dat[male,],"group",c(feat,m2_covar))
    m3 <- glm_model(dat,"group",c(feat,m3_covar_sbp,'Sex'))
    m3.f <- glm_model(dat[female,],"group",c(feat,m3_covar_sbp,'Menopause_code'))
    m3.m <- glm_model(dat[male,],"group",c(feat,m3_covar_sbp))
    m4 <- glm_model(dat,"group",c(feat,m4_covar_sbp,'Sex'))
    m4.f <- glm_model(dat[female,],"group",c(feat,m4_covar_sbp,'Menopause_code'))
    m4.m <- glm_model(dat[male,],"group",c(feat,m4_covar_sbp))
    cor.matrix[nrow(cor.matrix)+1,] <- c(group,"All","m1",feat,summary(m1)$coefficients[2,4],sign(summary(m1)$coefficients[2,1]),summary(m1)$coefficients[2,1],summary(m1)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(group,"Women","m1",feat,summary(m1.f)$coefficients[2,4],sign(summary(m1.f)$coefficients[2,1]),summary(m1.f)$coefficients[2,1],summary(m1.f)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(group,"Men","m1",feat,summary(m1.m)$coefficients[2,4],sign(summary(m1.m)$coefficients[2,1]),summary(m1.m)$coefficients[2,1],summary(m1.m)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(group,"All","m2",feat,summary(m2)$coefficients[2,4],sign(summary(m2)$coefficients[2,1]),summary(m2)$coefficients[2,1],summary(m2)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(group,"Women","m2",feat,summary(m2.f)$coefficients[2,4],sign(summary(m2.f)$coefficients[2,1]),summary(m2.f)$coefficients[2,1],summary(m2.f)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(group,"Men","m2",feat,summary(m2.m)$coefficients[2,4],sign(summary(m2.m)$coefficients[2,1]),summary(m2.m)$coefficients[2,1],summary(m2.m)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(group,"All","m3",feat,summary(m3)$coefficients[2,4],sign(summary(m3)$coefficients[2,1]),summary(m3)$coefficients[2,1],summary(m3)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(group,"Women","m3",feat,summary(m3.f)$coefficients[2,4],sign(summary(m3.f)$coefficients[2,1]),summary(m3.f)$coefficients[2,1],summary(m3.f)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(group,"Men","m3",feat,summary(m3.m)$coefficients[2,4],sign(summary(m3.m)$coefficients[2,1]),summary(m3.m)$coefficients[2,1],summary(m3.m)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(group,"All","m4",feat,summary(m4)$coefficients[2,4],sign(summary(m4)$coefficients[2,1]),summary(m4)$coefficients[2,1],summary(m4)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(group,"Women","m4",feat,summary(m4.f)$coefficients[2,4],sign(summary(m4.f)$coefficients[2,1]),summary(m4.f)$coefficients[2,1],summary(m4.f)$coefficients[2,2])
    cor.matrix[nrow(cor.matrix)+1,] <- c(group,"Men","m4",feat,summary(m4.m)$coefficients[2,4],sign(summary(m4.m)$coefficients[2,1]),summary(m4.m)$coefficients[2,1],summary(m4.m)$coefficients[2,2])
  }
}
cor.matrix$pval <- as.numeric(cor.matrix$pval)
filt.feat <- cor.matrix[cor.matrix$pval<=0.01,"feature"]
cor.matrix <- cor.matrix[cor.matrix[,"feature"]%in%filt.feat,]
filt.feat <- cor.matrix[cor.matrix$pval<=0.01,"group"]
cor.matrix <- cor.matrix[cor.matrix[,"group"]%in%filt.feat,]
cor.matrix$sym <- ifelse(cor.matrix$pval>=0.05, "", ifelse(cor.matrix$pval<=0.01,ifelse(cor.matrix$pval<=0.001,"***", "**"), "*"))
cor.matrix$log <- -log10(cor.matrix$pval)*as.numeric(cor.matrix$direction)
cor.matrix$feature <- gsub("_"," ",cor.matrix$feature)
cor.matrix$feature <- gsub("X","",cor.matrix$feature)
cor.matrix$feature <- gsub("\\."," ",cor.matrix$feature)
cor.matrix$gender <- factor(cor.matrix$gender, levels = c("All","Women","Men"))
pdf("Figures/300_bp_variability/Dipping_classification_vs_GM.pdf", width = 8, height = 5)
print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(name = "-log10", midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
        coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black") + facet_wrap(~gender+group) +
        theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"),
              legend.text = element_text(size = 12, face="bold"), legend.title = element_text(size = 12, face="bold")) + ggtitle("Dipping classification"))
dev.off()
write.csv(cor.matrix[,c('group','gender','model','feature','pval','coefficient','se')], "Figures/300_bp_variability/Dipping_classification_vs_GM.csv",row.names = F)
}

#############################################
###BP Variability in Hypertensive patients###
#############################################
hyp_pat <- rownames(meta.data[meta.data$Hypertension.Staging==1,])
hyp.merged.dat <- merged.dat[hyp_pat,]
females <- rownames(meta.data[meta.data$Hypertension.Staging==1&meta.data$Sex=="Women",])
males <- rownames(meta.data[meta.data$Hypertension.Staging==1&meta.data$Sex=="Men",])

#GM univariate
pdf("Figures/300_bp_variability/HYP_ONLY_BP_variability_vs_GM.pdf", width = 8, height = 5)
for(feat in SBP_cov){
  cor.matrix <-  data.frame(gender = character(), model = character(), feature = character(), pval = character(), direction = character())
  for(species in colnames(transformed.species.data)){
    m1 <- lm_model(hyp.merged.dat,feat,species)
    #m1.f <- lm_model(hyp.merged.dat[females,],feat,species)
    m1.m <- lm_model(hyp.merged.dat[males,],feat,species)
    m2 <- lm_model(hyp.merged.dat,feat,c(species,m2_covar,'Sex'))
    #m2.f <- lm_model(hyp.merged.dat[females,],feat,c(species,m2_covar))
    m2.m <- lm_model(hyp.merged.dat[males,],feat,c(species,m2_covar))
    m3 <- lm_model(hyp.merged.dat,feat,c(species,m3_covar_sbp,'Sex'))
    #m3.f <- lm_model(hyp.merged.dat[females,],feat,c(species,m3_covar_sbp,'Menopause_code'))
    m3.m <- lm_model(hyp.merged.dat[males,],feat,c(species,m3_covar_sbp))
    m4 <- lm_model(hyp.merged.dat,feat,c(species,m4_covar,'Sex'))
    #m4.f <- lm_model(hyp.merged.dat[females,],feat,c(species,m4_covar,'Menopause_code'))
    m4.m <- lm_model(hyp.merged.dat[males,],feat,c(species,m4_covar))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m1",species,summary(m1)$coefficients[2,4],sign(summary(m1)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m1",species,1,1)
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m1",species,summary(m1.m)$coefficients[2,4],sign(summary(m1.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m2",species,summary(m2)$coefficients[2,4],sign(summary(m2)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m2",species,1,1)
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m2",species,summary(m2.m)$coefficients[2,4],sign(summary(m2.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3",species,summary(m3)$coefficients[2,4],sign(summary(m3)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3",species,1,1)
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3",species,summary(m3.m)$coefficients[2,4],sign(summary(m3.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m4",species,summary(m4)$coefficients[2,4],sign(summary(m4)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m4",species,1,1)
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m4",species,summary(m4.m)$coefficients[2,4],sign(summary(m4.m)$coefficients[2,1]))
  }
  cor.matrix$pval <- as.numeric(cor.matrix$pval)
  filt.feat <- cor.matrix[cor.matrix$pval<=0.01,"feature"]
  cor.matrix <- cor.matrix[cor.matrix[,"feature"]%in%filt.feat,]
  cor.matrix$sym <- ifelse(cor.matrix$pval>=0.05, "", ifelse(cor.matrix$pval<=0.01,ifelse(cor.matrix$pval<=0.001,"***", "**"), "*"))
  cor.matrix$log <- -log10(cor.matrix$pval)*as.numeric(cor.matrix$direction)
  cor.matrix$feature <- gsub("_"," ",cor.matrix$feature)
  title <- gsub("\\."," ",feat)
  title <- gsub("X","",title)
  title <- gsub("_","-",title)
  cor.matrix$gender <- factor(cor.matrix$gender, levels = c("All","Women","Men"))
  print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(name = "-log10", midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
          coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender)  +
          theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"),
                legend.text = element_text(size = 12, face="bold"), legend.title = element_text(size = 12, face="bold")) + ggtitle(title))
}
for(feat in SBP_sd){
  cor.matrix <-  data.frame(gender = character(), model = character(), feature = character(), pval = character(), direction = character())
  for(species in colnames(transformed.species.data)){
    m1 <- lm_model(hyp.merged.dat,feat,species)
    m1.f <- lm_model(hyp.merged.dat[females,],feat,species)
    m1.m <- lm_model(hyp.merged.dat[males,],feat,species)
    m2 <- lm_model(hyp.merged.dat,feat,c(species,m2_covar,'Sex'))
    m2.f <- lm_model(hyp.merged.dat[females,],feat,c(species,m2_covar))
    m2.m <- lm_model(hyp.merged.dat[males,],feat,c(species,m2_covar))
    m3 <- lm_model(hyp.merged.dat,feat,c(species,m3_covar_sbp,'Sex'))
    m3.f <- lm_model(hyp.merged.dat[females,],feat,c(species,m3_covar_sbp,'Menopause_code'))
    m3.m <- lm_model(hyp.merged.dat[males,],feat,c(species,m3_covar_sbp))
    m4 <- lm_model(hyp.merged.dat,feat,c(species,m4_covar,'Sex'))
    m4.f <- lm_model(hyp.merged.dat[females,],feat,c(species,m4_covar,'Menopause_code'))
    m4.m <- lm_model(hyp.merged.dat[males,],feat,c(species,m4_covar))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m1",species,summary(m1)$coefficients[2,4],sign(summary(m1)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m1",species,summary(m1.f)$coefficients[2,4],sign(summary(m1.f)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m1",species,summary(m1.m)$coefficients[2,4],sign(summary(m1.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m2",species,summary(m2)$coefficients[2,4],sign(summary(m2)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m2",species,summary(m2.f)$coefficients[2,4],sign(summary(m2.f)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m2",species,summary(m2.m)$coefficients[2,4],sign(summary(m2.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3",species,summary(m3)$coefficients[2,4],sign(summary(m3)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3",species,summary(m3.f)$coefficients[2,4],sign(summary(m3.f)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3",species,summary(m3.m)$coefficients[2,4],sign(summary(m3.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m4",species,summary(m4)$coefficients[2,4],sign(summary(m4)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m4",species,summary(m4.f)$coefficients[2,4],sign(summary(m4.f)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m4",species,summary(m4.m)$coefficients[2,4],sign(summary(m4.m)$coefficients[2,1]))
  }
  cor.matrix$pval <- as.numeric(cor.matrix$pval)
  filt.feat <- cor.matrix[cor.matrix$pval<=0.01,"feature"]
  cor.matrix <- cor.matrix[cor.matrix[,"feature"]%in%filt.feat,]
  cor.matrix$sym <- ifelse(cor.matrix$pval>=0.05, "", ifelse(cor.matrix$pval<=0.01,ifelse(cor.matrix$pval<=0.001,"***", "**"), "*"))
  cor.matrix$log <- -log10(cor.matrix$pval)*as.numeric(cor.matrix$direction)
  cor.matrix$feature <- gsub("_"," ",cor.matrix$feature)
  title <- gsub("\\."," ",feat)
  title <- gsub("X","",title)
  title <- gsub("_","-",title)
  cor.matrix$gender <- factor(cor.matrix$gender, levels = c("All","Women","Men"))
  print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(name = "-log10", midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
          coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender)  +
          theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"),
                legend.text = element_text(size = 12, face="bold"), legend.title = element_text(size = 12, face="bold")) + ggtitle(title))
}
for(feat in DBP_cov){
  cor.matrix <-  data.frame(gender = character(), model = character(), feature = character(), pval = character(), direction = character())
  for(species in colnames(transformed.species.data)){
    m1 <- lm_model(hyp.merged.dat,feat,species)
    m1.f <- lm_model(hyp.merged.dat[females,],feat,species)
    m1.m <- lm_model(hyp.merged.dat[males,],feat,species)
    m2 <- lm_model(hyp.merged.dat,feat,c(species,m2_covar,'Sex'))
    m2.f <- lm_model(hyp.merged.dat[females,],feat,c(species,m2_covar))
    m2.m <- lm_model(hyp.merged.dat[males,],feat,c(species,m2_covar))
    m3 <- lm_model(hyp.merged.dat,feat,c(species,m3_covar_dbp,'Sex'))
    m3.f <- lm_model(hyp.merged.dat[females,],feat,c(species,m3_covar_dbp,'Menopause_code'))
    m3.m <- lm_model(hyp.merged.dat[males,],feat,c(species,m3_covar_dbp))
    m4 <- lm_model(hyp.merged.dat,feat,c(species,m4_covar,'Sex'))
    m4.f <- lm_model(hyp.merged.dat[females,],feat,c(species,m4_covar,'Menopause_code'))
    m4.m <- lm_model(hyp.merged.dat[males,],feat,c(species,m4_covar))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m1",species,summary(m1)$coefficients[2,4],sign(summary(m1)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m1",species,summary(m1.f)$coefficients[2,4],sign(summary(m1.f)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m1",species,summary(m1.m)$coefficients[2,4],sign(summary(m1.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m2",species,summary(m2)$coefficients[2,4],sign(summary(m2)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m2",species,summary(m2.f)$coefficients[2,4],sign(summary(m2.f)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m2",species,summary(m2.m)$coefficients[2,4],sign(summary(m2.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3",species,summary(m3)$coefficients[2,4],sign(summary(m3)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3",species,summary(m3.f)$coefficients[2,4],sign(summary(m3.f)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3",species,summary(m3.m)$coefficients[2,4],sign(summary(m3.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m4",species,summary(m4)$coefficients[2,4],sign(summary(m4)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m4",species,summary(m4.f)$coefficients[2,4],sign(summary(m4.f)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m4",species,summary(m4.m)$coefficients[2,4],sign(summary(m4.m)$coefficients[2,1]))
  }
  cor.matrix$pval <- as.numeric(cor.matrix$pval)
  filt.feat <- cor.matrix[cor.matrix$pval<=0.01,"feature"]
  cor.matrix <- cor.matrix[cor.matrix[,"feature"]%in%filt.feat,]
  cor.matrix$sym <- ifelse(cor.matrix$pval>=0.05, "", ifelse(cor.matrix$pval<=0.01,ifelse(cor.matrix$pval<=0.001,"***", "**"), "*"))
  cor.matrix$log <- -log10(cor.matrix$pval)*as.numeric(cor.matrix$direction)
  cor.matrix$feature <- gsub("_"," ",cor.matrix$feature)
  title <- gsub("\\."," ",feat)
  title <- gsub("X","",title)
  title <- gsub("_","-",title)
  cor.matrix$gender <- factor(cor.matrix$gender, levels = c("All","Women","Men"))
  print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(name = "-log10", midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
          coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender)  +
          theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"),
                legend.text = element_text(size = 12, face="bold"), legend.title = element_text(size = 12, face="bold")) + ggtitle(title))
}
for(feat in DBP_sd){
  cor.matrix <-  data.frame(gender = character(), model = character(), feature = character(), pval = character(), direction = character())
  for(species in colnames(transformed.species.data)){
    m1 <- lm_model(hyp.merged.dat,feat,species)
    m1.f <- lm_model(hyp.merged.dat[females,],feat,species)
    m1.m <- lm_model(hyp.merged.dat[males,],feat,species)
    m2 <- lm_model(hyp.merged.dat,feat,c(species,m2_covar,'Sex'))
    m2.f <- lm_model(hyp.merged.dat[females,],feat,c(species,m2_covar))
    m2.m <- lm_model(hyp.merged.dat[males,],feat,c(species,m2_covar))
    m3 <- lm_model(hyp.merged.dat,feat,c(species,m3_covar_dbp,'Sex'))
    m3.f <- lm_model(hyp.merged.dat[females,],feat,c(species,m3_covar_dbp,'Menopause_code'))
    m3.m <- lm_model(hyp.merged.dat[males,],feat,c(species,m3_covar_dbp))
    m4 <- lm_model(hyp.merged.dat,feat,c(species,m4_covar,'Sex'))
    m4.f <- lm_model(hyp.merged.dat[females,],feat,c(species,m4_covar,'Menopause_code'))
    m4.m <- lm_model(hyp.merged.dat[males,],feat,c(species,m4_covar))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m1",species,summary(m1)$coefficients[2,4],sign(summary(m1)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m1",species,summary(m1.f)$coefficients[2,4],sign(summary(m1.f)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m1",species,summary(m1.m)$coefficients[2,4],sign(summary(m1.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m2",species,summary(m2)$coefficients[2,4],sign(summary(m2)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m2",species,summary(m2.f)$coefficients[2,4],sign(summary(m2.f)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m2",species,summary(m2.m)$coefficients[2,4],sign(summary(m2.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3",species,summary(m3)$coefficients[2,4],sign(summary(m3)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3",species,summary(m3.f)$coefficients[2,4],sign(summary(m3.f)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3",species,summary(m3.m)$coefficients[2,4],sign(summary(m3.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m4",species,summary(m4)$coefficients[2,4],sign(summary(m4)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m4",species,summary(m4.f)$coefficients[2,4],sign(summary(m4.f)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m4",species,summary(m4.m)$coefficients[2,4],sign(summary(m4.m)$coefficients[2,1]))
  }
  cor.matrix$pval <- as.numeric(cor.matrix$pval)
  filt.feat <- cor.matrix[cor.matrix$pval<=0.01,"feature"]
  cor.matrix <- cor.matrix[cor.matrix[,"feature"]%in%filt.feat,]
  cor.matrix$sym <- ifelse(cor.matrix$pval>=0.05, "", ifelse(cor.matrix$pval<=0.01,ifelse(cor.matrix$pval<=0.001,"***", "**"), "*"))
  cor.matrix$log <- -log10(cor.matrix$pval)*as.numeric(cor.matrix$direction)
  cor.matrix$feature <- gsub("_"," ",cor.matrix$feature)
  title <- gsub("\\."," ",feat)
  title <- gsub("X","",title)
  title <- gsub("_","-",title)
  cor.matrix$gender <- factor(cor.matrix$gender, levels = c("All","Women","Men"))
  print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(name = "-log10", midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
          coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender)  +
          theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"),
                legend.text = element_text(size = 12, face="bold"), legend.title = element_text(size = 12, face="bold")) + ggtitle(title))
}
dev.off()

hyp_dip_subset <- hyp.merged.dat[hyp.merged.dat$Dipping.classification%in%c(1,2,3,4),]
pdf("Figures/300_bp_variability/HYP_ONLY_Dipping_vs_GM_boxplot.pdf", width = 8, height = 5)
for(feat in colnames(species.data)){
  my_boxplot(hyp_dip_subset,feat,"Dipping.classification",ylab = feat)
}
dev.off()

# Correlation analysis w/ Immune markers
pdf("Figures/300_bp_variability/HYP_ONLY_BP_variability_vs_immune.pdf", width = 8, height = 5)
for(feat in c(SBP_cov,SBP_sd)){
  cor.matrix <-  data.frame(gender = character(), model = character(), feature = character(), pval = character(), direction = character())
  for(species in immune_markers){
    m1 <- lm_model(hyp.merged.dat,feat,species)
    #m1.f <- lm_model(hyp.merged.dat[female,],feat,species)
    m1.m <- lm_model(hyp.merged.dat[male,],feat,species)
    m2 <- lm_model(hyp.merged.dat,feat,c(species,m2_covar,'Sex'))
    #m2.f <- lm_model(hyp.merged.dat[female,],feat,c(species,m2_covar))
    m2.m <- lm_model(hyp.merged.dat[male,],feat,c(species,m2_covar))
    m3 <- lm_model(hyp.merged.dat,feat,c(species,m3_covar_sbp,'Sex'))
    #m3.f <- lm_model(hyp.merged.dat[female,],feat,c(species,m3_covar_sbp,'Menopause_code'))
    m3.m <- lm_model(hyp.merged.dat[male,],feat,c(species,m3_covar_sbp))
    m4 <- lm_model(hyp.merged.dat,feat,c(species,m4_covar,'Sex'))
    #m4.f <- lm_model(hyp.merged.dat[female,],feat,c(species,m4_covar,'Menopause_code'))
    m4.m <- lm_model(hyp.merged.dat[male,],feat,c(species,m4_covar))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m1",species,summary(m1)$coefficients[2,4],sign(summary(m1)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m1",species,1,1)
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m1",species,summary(m1.m)$coefficients[2,4],sign(summary(m1.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m2",species,summary(m2)$coefficients[2,4],sign(summary(m2)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m2",species,1,1)
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m2",species,summary(m2.m)$coefficients[2,4],sign(summary(m2.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3",species,summary(m3)$coefficients[2,4],sign(summary(m3)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3",species,1,1)
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3",species,summary(m3.m)$coefficients[2,4],sign(summary(m3.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m4",species,summary(m4)$coefficients[2,4],sign(summary(m4)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m4",species,1,1)
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m4",species,summary(m4.m)$coefficients[2,4],sign(summary(m4.m)$coefficients[2,1]))
  }
  cor.matrix$pval <- as.numeric(cor.matrix$pval)
  filt.feat <- cor.matrix[cor.matrix$pval<=0.05,"feature"]
  if(length(filt.feat)==0){next}
  cor.matrix <- cor.matrix[cor.matrix[,"feature"]%in%filt.feat,]
  cor.matrix$sym <- ifelse(cor.matrix$pval>=0.05, "", ifelse(cor.matrix$pval<=0.01,ifelse(cor.matrix$pval<=0.001,"***", "**"), "*"))
  cor.matrix$log <- -log10(cor.matrix$pval)*as.numeric(cor.matrix$direction)
  cor.matrix$feature <- gsub("_"," ",cor.matrix$feature)
  title <- gsub("\\."," ",feat)
  title <- gsub("X","",title)
  title <- gsub("_","-",title)
  cor.matrix$gender <- factor(cor.matrix$gender, levels = c("All","Women","Men"))
  print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(name = "-log10", midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
          coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender)  +
          theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"),
                legend.text = element_text(size = 12, face="bold"), legend.title = element_text(size = 12, face="bold")) + ggtitle(title))
}
for(feat in c(DBP_cov,DBP_sd)){
  cor.matrix <-  data.frame(gender = character(), model = character(), feature = character(), pval = character(), direction = character())
  for(species in immune_markers){
    m1 <- lm_model(merged.dat,feat,species)
    #m1.f <- lm_model(merged.dat[female,],feat,species)
    m1.m <- lm_model(merged.dat[male,],feat,species)
    m2 <- lm_model(merged.dat,feat,c(species,m2_covar,'Sex'))
    #m2.f <- lm_model(merged.dat[female,],feat,c(species,m2_covar))
    m2.m <- lm_model(merged.dat[male,],feat,c(species,m2_covar))
    m3 <- lm_model(merged.dat,feat,c(species,m3_covar_dbp,'Sex'))
    #m3.f <- lm_model(merged.dat[female,],feat,c(species,m3_covar_dbp,'Menopause_code'))
    m3.m <- lm_model(merged.dat[male,],feat,c(species,m3_covar_dbp))
    m4 <- lm_model(merged.dat,feat,c(species,m4_covar,'Sex'))
    #m4.f <- lm_model(merged.dat[female,],feat,c(species,m4_covar,'Menopause_code'))
    m4.m <- lm_model(merged.dat[male,],feat,c(species,m4_covar))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m1",species,summary(m1)$coefficients[2,4],sign(summary(m1)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m1",species,1,1)
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m1",species,summary(m1.m)$coefficients[2,4],sign(summary(m1.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m2",species,summary(m2)$coefficients[2,4],sign(summary(m2)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m2",species,1,1)
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m2",species,summary(m2.m)$coefficients[2,4],sign(summary(m2.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3",species,summary(m3)$coefficients[2,4],sign(summary(m3)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3",species,1,1)
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3",species,summary(m3.m)$coefficients[2,4],sign(summary(m3.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m4",species,summary(m4)$coefficients[2,4],sign(summary(m4)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m4",species,1,1)
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m4",species,summary(m4.m)$coefficients[2,4],sign(summary(m4.m)$coefficients[2,1]))
  }
  cor.matrix$pval <- as.numeric(cor.matrix$pval)
  filt.feat <- cor.matrix[cor.matrix$pval<=0.05,"feature"]
  if(length(filt.feat)==0){next}
  cor.matrix <- cor.matrix[cor.matrix[,"feature"]%in%filt.feat,]
  cor.matrix$sym <- ifelse(cor.matrix$pval>=0.05, "", ifelse(cor.matrix$pval<=0.01,ifelse(cor.matrix$pval<=0.001,"***", "**"), "*"))
  cor.matrix$log <- -log10(cor.matrix$pval)*as.numeric(cor.matrix$direction)
  cor.matrix$feature <- gsub("_"," ",cor.matrix$feature)
  title <- gsub("\\."," ",feat)
  title <- gsub("X","",title)
  title <- gsub("_","-",title)
  cor.matrix$gender <- factor(cor.matrix$gender, levels = c("All","Women","Men"))
  print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(name = "-log10", midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
          coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender)  +
          theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"),
                legend.text = element_text(size = 12, face="bold"), legend.title = element_text(size = 12, face="bold")) + ggtitle(title))
}
dev.off()

# Correlation analysis w/ Blood SCFAs
pdf("Figures/300_bp_variability/HYP_ONLY_BP_variability_vs_SCFA.pdf", width = 8, height = 5)
for(feat in c(SBP_var,mbps_var)){
  cor.matrix <-  data.frame(gender = character(), model = character(), feature = character(), pval = character(), direction = character())
  for(species in acids_variables){
    m1 <- lm_model(hyp.merged.dat,feat,species)
    #m1.f <- lm_model(hyp.merged.dat[female,],feat,species)
    m1.m <- lm_model(hyp.merged.dat[male,],feat,species)
    m2 <- lm_model(hyp.merged.dat,feat,c(species,m2_covar,'Sex'))
    #m2.f <- lm_model(hyp.merged.dat[female,],feat,c(species,m2_covar))
    m2.m <- lm_model(hyp.merged.dat[male,],feat,c(species,m2_covar))
    m3 <- lm_model(hyp.merged.dat,feat,c(species,m3_covar_sbp,'Sex'))
    #m3.f <- lm_model(hyp.merged.dat[female,],feat,c(species,m3_covar_sbp,'Menopause_code'))
    m3.m <- lm_model(hyp.merged.dat[male,],feat,c(species,m3_covar_sbp))
    m4 <- lm_model(hyp.merged.dat,feat,c(species,m4_covar,'Sex'))
    #m4.f <- lm_model(hyp.merged.dat[female,],feat,c(species,m4_covar,'Menopause_code'))
    m4.m <- lm_model(hyp.merged.dat[male,],feat,c(species,m4_covar))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m1",species,summary(m1)$coefficients[2,4],sign(summary(m1)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m1",species,1,1)
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m1",species,summary(m1.m)$coefficients[2,4],sign(summary(m1.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m2",species,summary(m2)$coefficients[2,4],sign(summary(m2)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m2",species,1,1)
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m2",species,summary(m2.m)$coefficients[2,4],sign(summary(m2.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3",species,summary(m3)$coefficients[2,4],sign(summary(m3)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3",species,1,1)
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3",species,summary(m3.m)$coefficients[2,4],sign(summary(m3.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m4",species,summary(m4)$coefficients[2,4],sign(summary(m4)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m4",species,1,1)
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m4",species,summary(m4.m)$coefficients[2,4],sign(summary(m4.m)$coefficients[2,1]))
  }
  cor.matrix$pval <- as.numeric(cor.matrix$pval)
  filt.feat <- cor.matrix[cor.matrix$pval<=0.05,"feature"]
  if(length(filt.feat)==0){next}
  cor.matrix <- cor.matrix[cor.matrix[,"feature"]%in%filt.feat,]
  cor.matrix$sym <- ifelse(cor.matrix$pval>=0.05, "", ifelse(cor.matrix$pval<=0.01,ifelse(cor.matrix$pval<=0.001,"***", "**"), "*"))
  cor.matrix$log <- -log10(cor.matrix$pval)*as.numeric(cor.matrix$direction)
  cor.matrix$feature <- gsub("_"," ",cor.matrix$feature)
  title <- gsub("\\."," ",feat)
  title <- gsub("X","",title)
  title <- gsub("_","-",title)
  cor.matrix$gender <- factor(cor.matrix$gender, levels = c("All","Women","Men"))
  print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(name = "-log10", midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
          coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender)  +
          theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"),
                legend.text = element_text(size = 12, face="bold"), legend.title = element_text(size = 12, face="bold")) + ggtitle(title))
}
for(feat in c(DBP_cov,DBP_sd)){
  cor.matrix <-  data.frame(gender = character(), model = character(), feature = character(), pval = character(), direction = character())
  for(species in acids_variables){
    m1 <- lm_model(merged.dat,feat,species)
    #m1.f <- lm_model(merged.dat[female,],feat,species)
    m1.m <- lm_model(merged.dat[male,],feat,species)
    m2 <- lm_model(merged.dat,feat,c(species,m2_covar,'Sex'))
    #m2.f <- lm_model(merged.dat[female,],feat,c(species,m2_covar))
    m2.m <- lm_model(merged.dat[male,],feat,c(species,m2_covar))
    m3 <- lm_model(merged.dat,feat,c(species,m3_covar_dbp,'Sex'))
    #m3.f <- lm_model(merged.dat[female,],feat,c(species,m3_covar_dbp,'Menopause_code'))
    m3.m <- lm_model(merged.dat[male,],feat,c(species,m3_covar_dbp))
    m4 <- lm_model(merged.dat,feat,c(species,m4_covar,'Sex'))
    #m4.f <- lm_model(merged.dat[female,],feat,c(species,m4_covar,'Menopause_code'))
    m4.m <- lm_model(merged.dat[male,],feat,c(species,m4_covar))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m1",species,summary(m1)$coefficients[2,4],sign(summary(m1)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m1",species,1,1)
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m1",species,summary(m1.m)$coefficients[2,4],sign(summary(m1.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m2",species,summary(m2)$coefficients[2,4],sign(summary(m2)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m2",species,1,1)
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m2",species,summary(m2.m)$coefficients[2,4],sign(summary(m2.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3",species,summary(m3)$coefficients[2,4],sign(summary(m3)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3",species,1,1)
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3",species,summary(m3.m)$coefficients[2,4],sign(summary(m3.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m4",species,summary(m4)$coefficients[2,4],sign(summary(m4)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m4",species,1,1)
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m4",species,summary(m4.m)$coefficients[2,4],sign(summary(m4.m)$coefficients[2,1]))
  }
  cor.matrix$pval <- as.numeric(cor.matrix$pval)
  filt.feat <- cor.matrix[cor.matrix$pval<=0.05,"feature"]
  if(length(filt.feat)==0){next}
  cor.matrix <- cor.matrix[cor.matrix[,"feature"]%in%filt.feat,]
  cor.matrix$sym <- ifelse(cor.matrix$pval>=0.05, "", ifelse(cor.matrix$pval<=0.01,ifelse(cor.matrix$pval<=0.001,"***", "**"), "*"))
  cor.matrix$log <- -log10(cor.matrix$pval)*as.numeric(cor.matrix$direction)
  cor.matrix$feature <- gsub("_"," ",cor.matrix$feature)
  title <- gsub("\\."," ",feat)
  title <- gsub("X","",title)
  title <- gsub("_","-",title)
  cor.matrix$gender <- factor(cor.matrix$gender, levels = c("All","Women","Men"))
  print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(name = "-log10", midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
          coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender)  +
          theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"),
                legend.text = element_text(size = 12, face="bold"), legend.title = element_text(size = 12, face="bold")) + ggtitle(title))
}
dev.off()

# Correlation analysis w/ Stool SCFAs
pdf("Figures/300_bp_variability/HYP_ONLY_BP_variability_vs_stool_SCFA.pdf", width = 8, height = 5)
for(feat in c(SBP_cov,SBP_sd)){
  cor.matrix <-  data.frame(gender = character(), model = character(), feature = character(), pval = character(), direction = character())
  for(species in stool_scfa_variables){
    m1 <- lm_model(hyp.merged.dat,feat,species)
    #m1.f <- lm_model(hyp.merged.dat[female,],feat,species)
    m1.m <- lm_model(hyp.merged.dat[male,],feat,species)
    m2 <- lm_model(hyp.merged.dat,feat,c(species,m2_covar,'Sex'))
    #m2.f <- lm_model(hyp.merged.dat[female,],feat,c(species,m2_covar))
    m2.m <- lm_model(hyp.merged.dat[male,],feat,c(species,m2_covar))
    m3 <- lm_model(hyp.merged.dat,feat,c(species,m3_covar_sbp,'Sex'))
    #m3.f <- lm_model(hyp.merged.dat[female,],feat,c(species,m3_covar_sbp,'Menopause_code'))
    m3.m <- lm_model(hyp.merged.dat[male,],feat,c(species,m3_covar_sbp))
    m4 <- lm_model(hyp.merged.dat,feat,c(species,m4_covar,'Sex'))
    #m4.f <- lm_model(hyp.merged.dat[female,],feat,c(species,m4_covar,'Menopause_code'))
    m4.m <- lm_model(hyp.merged.dat[male,],feat,c(species,m4_covar))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m1",species,summary(m1)$coefficients[2,4],sign(summary(m1)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m1",species,1,1)
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m1",species,summary(m1.m)$coefficients[2,4],sign(summary(m1.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m2",species,summary(m2)$coefficients[2,4],sign(summary(m2)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m2",species,1,1)
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m2",species,summary(m2.m)$coefficients[2,4],sign(summary(m2.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3",species,summary(m3)$coefficients[2,4],sign(summary(m3)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3",species,1,1)
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3",species,summary(m3.m)$coefficients[2,4],sign(summary(m3.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m4",species,summary(m4)$coefficients[2,4],sign(summary(m4)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m4",species,1,1)
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m4",species,summary(m4.m)$coefficients[2,4],sign(summary(m4.m)$coefficients[2,1]))
  }
  cor.matrix$pval <- as.numeric(cor.matrix$pval)
  filt.feat <- cor.matrix[cor.matrix$pval<=0.05,"feature"]
  if(length(filt.feat)==0){next}
  cor.matrix <- cor.matrix[cor.matrix[,"feature"]%in%filt.feat,]
  cor.matrix$sym <- ifelse(cor.matrix$pval>=0.05, "", ifelse(cor.matrix$pval<=0.01,ifelse(cor.matrix$pval<=0.001,"***", "**"), "*"))
  cor.matrix$log <- -log10(cor.matrix$pval)*as.numeric(cor.matrix$direction)
  cor.matrix$feature <- gsub("_"," ",cor.matrix$feature)
  title <- gsub("\\."," ",feat)
  title <- gsub("X","",title)
  title <- gsub("_","-",title)
  cor.matrix$gender <- factor(cor.matrix$gender, levels = c("All","Women","Men"))
  print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(name = "-log10", midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
          coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender)  +
          theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"),
                legend.text = element_text(size = 12, face="bold"), legend.title = element_text(size = 12, face="bold")) + ggtitle(title))
}
for(feat in c(DBP_cov,DBP_sd)){
  cor.matrix <-  data.frame(gender = character(), model = character(), feature = character(), pval = character(), direction = character())
  for(species in stool_scfa_variables){
    m1 <- lm_model(merged.dat,feat,species)
    #m1.f <- lm_model(merged.dat[female,],feat,species)
    m1.m <- lm_model(merged.dat[male,],feat,species)
    m2 <- lm_model(merged.dat,feat,c(species,m2_covar,'Sex'))
    #m2.f <- lm_model(merged.dat[female,],feat,c(species,m2_covar))
    m2.m <- lm_model(merged.dat[male,],feat,c(species,m2_covar))
    m3 <- lm_model(merged.dat,feat,c(species,m3_covar_dbp,'Sex'))
    #m3.f <- lm_model(merged.dat[female,],feat,c(species,m3_covar_dbp,'Menopause_code'))
    m3.m <- lm_model(merged.dat[male,],feat,c(species,m3_covar_dbp))
    m4 <- lm_model(merged.dat,feat,c(species,m4_covar,'Sex'))
    #m4.f <- lm_model(merged.dat[female,],feat,c(species,m4_covar,'Menopause_code'))
    m4.m <- lm_model(merged.dat[male,],feat,c(species,m4_covar))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m1",species,summary(m1)$coefficients[2,4],sign(summary(m1)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m1",species,1,1)
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m1",species,summary(m1.m)$coefficients[2,4],sign(summary(m1.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m2",species,summary(m2)$coefficients[2,4],sign(summary(m2)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m2",species,1,1)
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m2",species,summary(m2.m)$coefficients[2,4],sign(summary(m2.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3",species,summary(m3)$coefficients[2,4],sign(summary(m3)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3",species,1,1)
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3",species,summary(m3.m)$coefficients[2,4],sign(summary(m3.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m4",species,summary(m4)$coefficients[2,4],sign(summary(m4)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m4",species,1,1)
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m4",species,summary(m4.m)$coefficients[2,4],sign(summary(m4.m)$coefficients[2,1]))
  }
  cor.matrix$pval <- as.numeric(cor.matrix$pval)
  filt.feat <- cor.matrix[cor.matrix$pval<=0.05,"feature"]
  if(length(filt.feat)==0){next}
  cor.matrix <- cor.matrix[cor.matrix[,"feature"]%in%filt.feat,]
  cor.matrix$sym <- ifelse(cor.matrix$pval>=0.05, "", ifelse(cor.matrix$pval<=0.01,ifelse(cor.matrix$pval<=0.001,"***", "**"), "*"))
  cor.matrix$log <- -log10(cor.matrix$pval)*as.numeric(cor.matrix$direction)
  cor.matrix$feature <- gsub("_"," ",cor.matrix$feature)
  title <- gsub("\\."," ",feat)
  title <- gsub("X","",title)
  title <- gsub("_","-",title)
  cor.matrix$gender <- factor(cor.matrix$gender, levels = c("All","Women","Men"))
  print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(name = "-log10", midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
          coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender)  +
          theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"),
                legend.text = element_text(size = 12, face="bold"), legend.title = element_text(size = 12, face="bold")) + ggtitle(title))
}
dev.off()

# Correlation analysis w/ Diet
pdf("Figures/300_bp_variability/HYP_ONLY_BP_variability_vs_diet.pdf", width = 8, height = 5)
for(feat in c(SBP_cov,SBP_sd)){
  cor.matrix <-  data.frame(gender = character(), model = character(), feature = character(), pval = character(), direction = character())
  for(species in c(diet_variables_continuous,food_group_variables,macro_micro_nutrients)){
    m1 <- lm_model(hyp.merged.dat,feat,species)
    #m1.f <- lm_model(hyp.merged.dat[female,],feat,species)
    m1.m <- lm_model(hyp.merged.dat[male,],feat,species)
    m2 <- lm_model(hyp.merged.dat,feat,c(species,m2_covar,'Sex'))
    #m2.f <- lm_model(hyp.merged.dat[female,],feat,c(species,m2_covar))
    m2.m <- lm_model(hyp.merged.dat[male,],feat,c(species,m2_covar))
    m3 <- lm_model(hyp.merged.dat,feat,c(species,m3_covar_sbp,'Sex'))
    #m3.f <- lm_model(hyp.merged.dat[female,],feat,c(species,m3_covar_sbp,'Menopause_code'))
    m3.m <- lm_model(hyp.merged.dat[male,],feat,c(species,m3_covar_sbp))
    m4 <- lm_model(hyp.merged.dat,feat,c(species,m4_covar,'Sex'))
    #m4.f <- lm_model(hyp.merged.dat[female,],feat,c(species,m4_covar,'Menopause_code'))
    m4.m <- lm_model(hyp.merged.dat[male,],feat,c(species,m4_covar))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m1",species,summary(m1)$coefficients[2,4],sign(summary(m1)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m1",species,1,1)
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m1",species,summary(m1.m)$coefficients[2,4],sign(summary(m1.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m2",species,summary(m2)$coefficients[2,4],sign(summary(m2)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m2",species,1,1)
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m2",species,summary(m2.m)$coefficients[2,4],sign(summary(m2.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3",species,summary(m3)$coefficients[2,4],sign(summary(m3)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3",species,1,1)
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3",species,summary(m3.m)$coefficients[2,4],sign(summary(m3.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m4",species,summary(m4)$coefficients[2,4],sign(summary(m4)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m4",species,1,1)
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m4",species,summary(m4.m)$coefficients[2,4],sign(summary(m4.m)$coefficients[2,1]))
  }
  cor.matrix$pval <- as.numeric(cor.matrix$pval)
  filt.feat <- cor.matrix[cor.matrix$pval<=0.05,"feature"]
  if(length(filt.feat)==0){next}
  cor.matrix <- cor.matrix[cor.matrix[,"feature"]%in%filt.feat,]
  cor.matrix$sym <- ifelse(cor.matrix$pval>=0.05, "", ifelse(cor.matrix$pval<=0.01,ifelse(cor.matrix$pval<=0.001,"***", "**"), "*"))
  cor.matrix$log <- -log10(cor.matrix$pval)*as.numeric(cor.matrix$direction)
  cor.matrix$feature <- gsub("_"," ",cor.matrix$feature)
  title <- gsub("\\."," ",feat)
  title <- gsub("X","",title)
  title <- gsub("_","-",title)
  cor.matrix$gender <- factor(cor.matrix$gender, levels = c("All","Women","Men"))
  print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(name = "-log10", midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
          coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender)  +
          theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"),
                legend.text = element_text(size = 12, face="bold"), legend.title = element_text(size = 12, face="bold")) + ggtitle(title))
}
for(feat in c(DBP_cov,DBP_sd)){
  cor.matrix <-  data.frame(gender = character(), model = character(), feature = character(), pval = character(), direction = character())
  for(species in c(diet_variables_continuous,food_group_variables,macro_micro_nutrients)){
    m1 <- lm_model(merged.dat,feat,species)
    #m1.f <- lm_model(merged.dat[female,],feat,species)
    m1.m <- lm_model(merged.dat[male,],feat,species)
    m2 <- lm_model(merged.dat,feat,c(species,m2_covar,'Sex'))
    #m2.f <- lm_model(merged.dat[female,],feat,c(species,m2_covar))
    m2.m <- lm_model(merged.dat[male,],feat,c(species,m2_covar))
    m3 <- lm_model(merged.dat,feat,c(species,m3_covar_dbp,'Sex'))
    #m3.f <- lm_model(merged.dat[female,],feat,c(species,m3_covar_dbp,'Menopause_code'))
    m3.m <- lm_model(merged.dat[male,],feat,c(species,m3_covar_dbp))
    m4 <- lm_model(merged.dat,feat,c(species,m4_covar,'Sex'))
    #m4.f <- lm_model(merged.dat[female,],feat,c(species,m4_covar,'Menopause_code'))
    m4.m <- lm_model(merged.dat[male,],feat,c(species,m4_covar))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m1",species,summary(m1)$coefficients[2,4],sign(summary(m1)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m1",species,1,1)
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m1",species,summary(m1.m)$coefficients[2,4],sign(summary(m1.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m2",species,summary(m2)$coefficients[2,4],sign(summary(m2)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m2",species,1,1)
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m2",species,summary(m2.m)$coefficients[2,4],sign(summary(m2.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3",species,summary(m3)$coefficients[2,4],sign(summary(m3)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3",species,1,1)
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3",species,summary(m3.m)$coefficients[2,4],sign(summary(m3.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m4",species,summary(m4)$coefficients[2,4],sign(summary(m4)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m4",species,1,1)
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m4",species,summary(m4.m)$coefficients[2,4],sign(summary(m4.m)$coefficients[2,1]))
  }
  cor.matrix$pval <- as.numeric(cor.matrix$pval)
  filt.feat <- cor.matrix[cor.matrix$pval<=0.05,"feature"]
  if(length(filt.feat)==0){next}
  cor.matrix <- cor.matrix[cor.matrix[,"feature"]%in%filt.feat,]
  cor.matrix$sym <- ifelse(cor.matrix$pval>=0.05, "", ifelse(cor.matrix$pval<=0.01,ifelse(cor.matrix$pval<=0.001,"***", "**"), "*"))
  cor.matrix$log <- -log10(cor.matrix$pval)*as.numeric(cor.matrix$direction)
  cor.matrix$feature <- gsub("_"," ",cor.matrix$feature)
  title <- gsub("\\."," ",feat)
  title <- gsub("X","",title)
  title <- gsub("_","-",title)
  cor.matrix$gender <- factor(cor.matrix$gender, levels = c("All","Women","Men"))
  print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(name = "-log10", midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
          coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender)  +
          theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"),
                legend.text = element_text(size = 12, face="bold"), legend.title = element_text(size = 12, face="bold")) + ggtitle(title))
}
dev.off()

#Alpha Diversity analysis
hyp.merged.dat$diversity <- vegan::diversity(transformed.species.data[hyp_pat,], index = 'shannon')

cor.matrix <-  data.frame(gender = character(), model = character(), feature = character(), pval = character(), direction = character())
for(feature in SBP_cov){
  m1 <- lm_model(hyp.merged.dat,'diversity',feature)
  #m1.f <- lm_model(hyp.merged.dat[females,],'diversity',feature)
  m1.m <- lm_model(hyp.merged.dat[males,],'diversity',feature)
  m2 <- lm_model(hyp.merged.dat,'diversity',c(feature,m2_covar,'Sex'))
  #m2.f <- lm_model(hyp.merged.dat[females,],'diversity',c(feature,m2_covar))
  m2.m <- lm_model(hyp.merged.dat[males,],'diversity',c(feature,m2_covar))
  m3 <- lm_model(hyp.merged.dat,'diversity',c(feature,m3_covar_sbp,'Sex'))
  #m3.f <- lm_model(hyp.merged.dat[females,],'diversity',c(feature,m3_covar_sbp,'Menopause_code'))
  m3.m <- lm_model(hyp.merged.dat[males,],'diversity',c(feature,m3_covar_sbp))
  m4 <- lm_model(hyp.merged.dat,'diversity',c(feature,m4_covar,'Sex'))
  #m4.f <- lm_model(hyp.merged.dat[females,],'diversity',c(feature,m4_covar,'Menopause_code'))
  m4.m <- lm_model(hyp.merged.dat[males,],'diversity',c(feature,m4_covar))
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m1",feature,summary(m1)$coefficients[2,4],sign(summary(m1)$coefficients[2,1]))
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m1",feature,1,1)
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m1",feature,summary(m1.m)$coefficients[2,4],sign(summary(m1.m)$coefficients[2,1]))
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m2",feature,summary(m2)$coefficients[2,4],sign(summary(m2)$coefficients[2,1]))
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m2",feature,1,1)
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m2",feature,summary(m2.m)$coefficients[2,4],sign(summary(m2.m)$coefficients[2,1]))
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3",feature,summary(m3)$coefficients[2,4],sign(summary(m3)$coefficients[2,1]))
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3",feature,1,1)
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3",feature,summary(m3.m)$coefficients[2,4],sign(summary(m3.m)$coefficients[2,1]))
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m4",feature,summary(m4)$coefficients[2,4],sign(summary(m4)$coefficients[2,1]))
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m4",feature,1,1)
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m4",feature,summary(m4.m)$coefficients[2,4],sign(summary(m4.m)$coefficients[2,1]))
}
cor.matrix <- format_res(cor.matrix)
pdf("Figures/300_bp_variability/HYP_ONLY_shannon_diversity_SBP_CoV.pdf", width = 8, height = 5)
print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(name = "-log10", midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
        coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender) + 
        theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"), 
              legend.text = element_text(size = 12, face="bold"), legend.title = element_text(size = 12, face="bold")) + ggtitle("Shannon Diversity"))
dev.off()

cor.matrix <-  data.frame(gender = character(), model = character(), feature = character(), pval = character(), direction = character())
for(feature in SBP_sd){
  m1 <- lm_model(hyp.merged.dat,'diversity',feature)
  m1.f <- lm_model(hyp.merged.dat[females,],'diversity',feature)
  m1.m <- lm_model(hyp.merged.dat[males,],'diversity',feature)
  m2 <- lm_model(hyp.merged.dat,'diversity',c(feature,m2_covar,'Sex'))
  m2.f <- lm_model(hyp.merged.dat[females,],'diversity',c(feature,m2_covar))
  m2.m <- lm_model(hyp.merged.dat[males,],'diversity',c(feature,m2_covar))
  m3 <- lm_model(hyp.merged.dat,'diversity',c(feature,m3_covar_sbp,'Sex'))
  m3.f <- lm_model(hyp.merged.dat[females,],'diversity',c(feature,m3_covar_sbp,'Menopause_code'))
  m3.m <- lm_model(hyp.merged.dat[males,],'diversity',c(feature,m3_covar_sbp))
  m4 <- lm_model(hyp.merged.dat,'diversity',c(feature,m4_covar,'Sex'))
  m4.f <- lm_model(hyp.merged.dat[females,],'diversity',c(feature,m4_covar,'Menopause_code'))
  m4.m <- lm_model(hyp.merged.dat[males,],'diversity',c(feature,m4_covar))
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m1",feature,summary(m1)$coefficients[2,4],sign(summary(m1)$coefficients[2,1]))
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m1",feature,summary(m1.f)$coefficients[2,4],sign(summary(m1.f)$coefficients[2,1]))
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m1",feature,summary(m1.m)$coefficients[2,4],sign(summary(m1.m)$coefficients[2,1]))
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m2",feature,summary(m2)$coefficients[2,4],sign(summary(m2)$coefficients[2,1]))
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m2",feature,summary(m2.f)$coefficients[2,4],sign(summary(m2.f)$coefficients[2,1]))
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m2",feature,summary(m2.m)$coefficients[2,4],sign(summary(m2.m)$coefficients[2,1]))
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3",feature,summary(m3)$coefficients[2,4],sign(summary(m3)$coefficients[2,1]))
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3",feature,summary(m3.f)$coefficients[2,4],sign(summary(m3.f)$coefficients[2,1]))
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3",feature,summary(m3.m)$coefficients[2,4],sign(summary(m3.m)$coefficients[2,1]))
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m4",feature,summary(m4)$coefficients[2,4],sign(summary(m4)$coefficients[2,1]))
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m4",feature,summary(m4.f)$coefficients[2,4],sign(summary(m4.f)$coefficients[2,1]))
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m4",feature,summary(m4.m)$coefficients[2,4],sign(summary(m4.m)$coefficients[2,1]))
}
cor.matrix <- format_res(cor.matrix)
pdf("Figures/300_bp_variability/HYP_ONLY_shannon_diversity_SBP_SD.pdf", width = 8, height = 5)
print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(name = "-log10", midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
        coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender) + 
        theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"), 
              legend.text = element_text(size = 12, face="bold"), legend.title = element_text(size = 12, face="bold")) + ggtitle("Shannon Diversity"))
dev.off()

cor.matrix <-  data.frame(gender = character(), model = character(), feature = character(), pval = character(), direction = character())
for(feature in DBP_cov){
  m1 <- lm_model(hyp.merged.dat,'diversity',feature)
  m1.f <- lm_model(hyp.merged.dat[females,],'diversity',feature)
  m1.m <- lm_model(hyp.merged.dat[males,],'diversity',feature)
  m2 <- lm_model(hyp.merged.dat,'diversity',c(feature,m2_covar,'Sex'))
  m2.f <- lm_model(hyp.merged.dat[females,],'diversity',c(feature,m2_covar))
  m2.m <- lm_model(hyp.merged.dat[males,],'diversity',c(feature,m2_covar))
  m3 <- lm_model(hyp.merged.dat,'diversity',c(feature,m3_covar_dbp,'Sex'))
  m3.f <- lm_model(hyp.merged.dat[females,],'diversity',c(feature,m3_covar_dbp,'Menopause_code'))
  m3.m <- lm_model(hyp.merged.dat[males,],'diversity',c(feature,m3_covar_dbp))
  m4 <- lm_model(hyp.merged.dat,'diversity',c(feature,m4_covar,'Sex'))
  m4.f <- lm_model(hyp.merged.dat[females,],'diversity',c(feature,m4_covar,'Menopause_code'))
  m4.m <- lm_model(hyp.merged.dat[males,],'diversity',c(feature,m4_covar))
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m1",feature,summary(m1)$coefficients[2,4],sign(summary(m1)$coefficients[2,1]))
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m1",feature,summary(m1.f)$coefficients[2,4],sign(summary(m1.f)$coefficients[2,1]))
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m1",feature,summary(m1.m)$coefficients[2,4],sign(summary(m1.m)$coefficients[2,1]))
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m2",feature,summary(m2)$coefficients[2,4],sign(summary(m2)$coefficients[2,1]))
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m2",feature,summary(m2.f)$coefficients[2,4],sign(summary(m2.f)$coefficients[2,1]))
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m2",feature,summary(m2.m)$coefficients[2,4],sign(summary(m2.m)$coefficients[2,1]))
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3",feature,summary(m3)$coefficients[2,4],sign(summary(m3)$coefficients[2,1]))
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3",feature,summary(m3.f)$coefficients[2,4],sign(summary(m3.f)$coefficients[2,1]))
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3",feature,summary(m3.m)$coefficients[2,4],sign(summary(m3.m)$coefficients[2,1]))
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m4",feature,summary(m4)$coefficients[2,4],sign(summary(m4)$coefficients[2,1]))
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m4",feature,summary(m4.f)$coefficients[2,4],sign(summary(m4.f)$coefficients[2,1]))
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m4",feature,summary(m4.m)$coefficients[2,4],sign(summary(m4.m)$coefficients[2,1]))
}
cor.matrix <- format_res(cor.matrix)
pdf("Figures/300_bp_variability/HYP_ONLY_shannon_diversity_DBP_CoV.pdf", width = 8, height = 5)
print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(name = "-log10", midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
        coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender) + 
        theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"), 
              legend.text = element_text(size = 12, face="bold"), legend.title = element_text(size = 12, face="bold")) + ggtitle("Shannon Diversity"))
dev.off()

cor.matrix <-  data.frame(gender = character(), model = character(), feature = character(), pval = character(), direction = character())
for(feature in DBP_sd){
  m1 <- lm_model(hyp.merged.dat,'diversity',feature)
  m1.f <- lm_model(hyp.merged.dat[females,],'diversity',feature)
  m1.m <- lm_model(hyp.merged.dat[males,],'diversity',feature)
  m2 <- lm_model(hyp.merged.dat,'diversity',c(feature,m2_covar,'Sex'))
  m2.f <- lm_model(hyp.merged.dat[females,],'diversity',c(feature,m2_covar))
  m2.m <- lm_model(hyp.merged.dat[males,],'diversity',c(feature,m2_covar))
  m3 <- lm_model(hyp.merged.dat,'diversity',c(feature,m3_covar_dbp,'Sex'))
  m3.f <- lm_model(hyp.merged.dat[females,],'diversity',c(feature,m3_covar_dbp,'Menopause_code'))
  m3.m <- lm_model(hyp.merged.dat[males,],'diversity',c(feature,m3_covar_dbp))
  m4 <- lm_model(hyp.merged.dat,'diversity',c(feature,m4_covar,'Sex'))
  m4.f <- lm_model(hyp.merged.dat[females,],'diversity',c(feature,m4_covar,'Menopause_code'))
  m4.m <- lm_model(hyp.merged.dat[males,],'diversity',c(feature,m4_covar))
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m1",feature,summary(m1)$coefficients[2,4],sign(summary(m1)$coefficients[2,1]))
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m1",feature,summary(m1.f)$coefficients[2,4],sign(summary(m1.f)$coefficients[2,1]))
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m1",feature,summary(m1.m)$coefficients[2,4],sign(summary(m1.m)$coefficients[2,1]))
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m2",feature,summary(m2)$coefficients[2,4],sign(summary(m2)$coefficients[2,1]))
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m2",feature,summary(m2.f)$coefficients[2,4],sign(summary(m2.f)$coefficients[2,1]))
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m2",feature,summary(m2.m)$coefficients[2,4],sign(summary(m2.m)$coefficients[2,1]))
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3",feature,summary(m3)$coefficients[2,4],sign(summary(m3)$coefficients[2,1]))
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3",feature,summary(m3.f)$coefficients[2,4],sign(summary(m3.f)$coefficients[2,1]))
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3",feature,summary(m3.m)$coefficients[2,4],sign(summary(m3.m)$coefficients[2,1]))
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m4",feature,summary(m4)$coefficients[2,4],sign(summary(m4)$coefficients[2,1]))
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m4",feature,summary(m4.f)$coefficients[2,4],sign(summary(m4.f)$coefficients[2,1]))
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m4",feature,summary(m4.m)$coefficients[2,4],sign(summary(m4.m)$coefficients[2,1]))
}
cor.matrix <- format_res(cor.matrix)
pdf("Figures/300_bp_variability/HYP_ONLY_shannon_diversity_DBP_SD.pdf", width = 8, height = 5)
print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(name = "-log10", midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
        coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender) + 
        theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"), 
              legend.text = element_text(size = 12, face="bold"), legend.title = element_text(size = 12, face="bold")) + ggtitle("Shannon Diversity"))
dev.off()

# Permanova - Beta diversity
cor.matrix <-  data.frame(gender = character(), model = character(), feature = character(), pval = character())
for(feature in SBP_cov){
  m1 <- permanova_model(hyp.merged.dat,species.data,feature)
  m1.f <- permanova_model(hyp.merged.dat[females,],species.data,feature)
  m1.m <- permanova_model(hyp.merged.dat[males,],species.data,feature)
  m2 <- permanova_model(hyp.merged.dat,species.data,c(feature,m2_covar,'Sex'))
  m2.f <- permanova_model(hyp.merged.dat[females,],species.data,c(feature,m2_covar))
  m2.m <- permanova_model(hyp.merged.dat[males,],species.data,c(feature,m2_covar))
  m3 <- permanova_model(hyp.merged.dat,species.data,c(feature,m3_covar_sbp,'Sex'))
  m3.f <- permanova_model(hyp.merged.dat[females,],species.data,c(feature,m3_covar_sbp,'Menopause_code'))
  m3.m <- permanova_model(hyp.merged.dat[males,],species.data,c(feature,m3_covar_sbp))
  m4 <- permanova_model(hyp.merged.dat,species.data,c(feature,m4_covar,'Sex'))
  m4.f <- permanova_model(hyp.merged.dat[females,],species.data,c(feature,m4_covar,'Menopause_code'))
  m4.m <- permanova_model(hyp.merged.dat[males,],species.data,c(feature,m4_covar))
  cor.matrix[nrow(cor.matrix)+1,] <-c("All","m1",feature,m1$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m1",feature,m1.f$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m1",feature,m1.m$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m2",feature,m2$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m2",feature,m2.f$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m2",feature,m2.m$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3",feature,m3$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3",feature,m3.f$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3",feature,m3.m$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m4",feature,m4$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m4",feature,m4.f$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m4",feature,m4.m$`Pr(>F)`[1])
}
cor.matrix <- format_res_pnova(cor.matrix)
pdf("Figures/300_bp_variability/HYP_ONLY_permanova_models_SBP_CoV.pdf", width = 8, height = 5)
print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", name ="-log10" ) +
        coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender) +
        theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"),
              legend.title = element_text(size = 12, face="bold"), legend.text =  element_text(size = 12, face="bold")))
dev.off()

cor.matrix <-  data.frame(gender = character(), model = character(), feature = character(), pval = character())
for(feature in SBP_sd){
  m1 <- permanova_model(hyp.merged.dat,species.data,feature)
  m1.f <- permanova_model(hyp.merged.dat[females,],species.data,feature)
  m1.m <- permanova_model(hyp.merged.dat[males,],species.data,feature)
  m2 <- permanova_model(hyp.merged.dat,species.data,c(feature,m2_covar,'Sex'))
  m2.f <- permanova_model(hyp.merged.dat[females,],species.data,c(feature,m2_covar))
  m2.m <- permanova_model(hyp.merged.dat[males,],species.data,c(feature,m2_covar))
  m3 <- permanova_model(hyp.merged.dat,species.data,c(feature,m3_covar_sbp,'Sex'))
  m3.f <- permanova_model(hyp.merged.dat[females,],species.data,c(feature,m3_covar_sbp,'Menopause_code'))
  m3.m <- permanova_model(hyp.merged.dat[males,],species.data,c(feature,m3_covar_sbp))
  m4 <- permanova_model(hyp.merged.dat,species.data,c(feature,m4_covar,'Sex'))
  m4.f <- permanova_model(hyp.merged.dat[females,],species.data,c(feature,m4_covar,'Menopause_code'))
  m4.m <- permanova_model(hyp.merged.dat[males,],species.data,c(feature,m4_covar))
  cor.matrix[nrow(cor.matrix)+1,] <-c("All","m1",feature,m1$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m1",feature,m1.f$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m1",feature,m1.m$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m2",feature,m2$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m2",feature,m2.f$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m2",feature,m2.m$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3",feature,m3$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3",feature,m3.f$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3",feature,m3.m$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m4",feature,m4$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m4",feature,m4.f$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m4",feature,m4.m$`Pr(>F)`[1])
}
cor.matrix <- format_res_pnova(cor.matrix)
pdf("Figures/300_bp_variability/HYP_ONLY_permanova_models_SBP_SD.pdf", width = 8, height = 5)
print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", name ="-log10" ) +
        coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender) +
        theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"),
              legend.title = element_text(size = 12, face="bold"), legend.text =  element_text(size = 12, face="bold")))
dev.off()

cor.matrix <-  data.frame(gender = character(), model = character(), feature = character(), pval = character())
for(feature in DBP_cov){
  m1 <- permanova_model(hyp.merged.dat,species.data,feature)
  m1.f <- permanova_model(hyp.merged.dat[females,],species.data,feature)
  m1.m <- permanova_model(hyp.merged.dat[males,],species.data,feature)
  m2 <- permanova_model(hyp.merged.dat,species.data,c(feature,m2_covar,'Sex'))
  m2.f <- permanova_model(hyp.merged.dat[females,],species.data,c(feature,m2_covar))
  m2.m <- permanova_model(hyp.merged.dat[males,],species.data,c(feature,m2_covar))
  m3 <- permanova_model(hyp.merged.dat,species.data,c(feature,m3_covar_dbp,'Sex'))
  m3.f <- permanova_model(hyp.merged.dat[females,],species.data,c(feature,m3_covar_dbp,'Menopause_code'))
  m3.m <- permanova_model(hyp.merged.dat[males,],species.data,c(feature,m3_covar_dbp))
  m4 <- permanova_model(hyp.merged.dat,species.data,c(feature,m4_covar,'Sex'))
  m4.f <- permanova_model(hyp.merged.dat[females,],species.data,c(feature,m4_covar,'Menopause_code'))
  m4.m <- permanova_model(hyp.merged.dat[males,],species.data,c(feature,m4_covar))
  cor.matrix[nrow(cor.matrix)+1,] <-c("All","m1",feature,m1$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m1",feature,m1.f$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m1",feature,m1.m$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m2",feature,m2$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m2",feature,m2.f$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m2",feature,m2.m$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3",feature,m3$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3",feature,m3.f$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3",feature,m3.m$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m4",feature,m4$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m4",feature,m4.f$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m4",feature,m4.m$`Pr(>F)`[1])
}
cor.matrix <- format_res_pnova(cor.matrix)
pdf("Figures/300_bp_variability/HYP_ONLY_permanova_models_DBP_CoV.pdf", width = 8, height = 5)
print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", name ="-log10" ) +
        coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender) +
        theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"),
              legend.title = element_text(size = 12, face="bold"), legend.text =  element_text(size = 12, face="bold")))
dev.off()

cor.matrix <-  data.frame(gender = character(), model = character(), feature = character(), pval = character())
for(feature in DBP_sd){
  m1 <- permanova_model(hyp.merged.dat,species.data,feature)
  m1.f <- permanova_model(hyp.merged.dat[females,],species.data,feature)
  m1.m <- permanova_model(hyp.merged.dat[males,],species.data,feature)
  m2 <- permanova_model(hyp.merged.dat,species.data,c(feature,m2_covar,'Sex'))
  m2.f <- permanova_model(hyp.merged.dat[females,],species.data,c(feature,m2_covar))
  m2.m <- permanova_model(hyp.merged.dat[males,],species.data,c(feature,m2_covar))
  m3 <- permanova_model(hyp.merged.dat,species.data,c(feature,m3_covar_dbp,'Sex'))
  m3.f <- permanova_model(hyp.merged.dat[females,],species.data,c(feature,m3_covar_dbp,'Menopause_code'))
  m3.m <- permanova_model(hyp.merged.dat[males,],species.data,c(feature,m3_covar_dbp))
  m4 <- permanova_model(hyp.merged.dat,species.data,c(feature,m4_covar,'Sex'))
  m4.f <- permanova_model(hyp.merged.dat[females,],species.data,c(feature,m4_covar,'Menopause_code'))
  m4.m <- permanova_model(hyp.merged.dat[males,],species.data,c(feature,m4_covar))
  cor.matrix[nrow(cor.matrix)+1,] <-c("All","m1",feature,m1$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m1",feature,m1.f$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m1",feature,m1.m$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m2",feature,m2$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m2",feature,m2.f$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m2",feature,m2.m$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3",feature,m3$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3",feature,m3.f$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3",feature,m3.m$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m4",feature,m4$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m4",feature,m4.f$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m4",feature,m4.m$`Pr(>F)`[1])
}
cor.matrix <- format_res_pnova(cor.matrix)
pdf("Figures/300_bp_variability/HYP_ONLY_permanova_models_DBP_SD.pdf", width = 8, height = 5)
print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", name ="-log10" ) +
        coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender) +
        theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"),
              legend.title = element_text(size = 12, face="bold"), legend.text =  element_text(size = 12, face="bold")))
dev.off()


##############################
### Control for sleep data ###
##############################
m3_l_covar <- c(m3_covar,"Latency")
m3_e_covar <- c(m3_covar,"Efficiency")
m3_l_e_covar <- c(m3_covar,"Latency","Efficiency")

# Assocaiton with Coef
summary(lm_model(merged.dat,"Daytime.SBP.CoV","Latency"))
summary(lm_model(merged.dat,"Daytime.DBP.CoV","Latency"))
summary(lm_model(merged.dat,"Nighttime.SBP.CoV","Latency"))
summary(lm_model(merged.dat,"Nighttime.DBP.CoV","Latency"))
summary(lm_model(merged.dat,"Daytime.SBP.CoV","Efficiency"))
summary(lm_model(merged.dat,"Daytime.DBP.CoV","Efficiency"))
summary(lm_model(merged.dat,"Nighttime.SBP.CoV","Efficiency"))
summary(lm_model(merged.dat,"Nighttime.DBP.CoV","Efficiency"))
summary(lm_model(merged.dat,"Latency","Efficiency"))
my_scatterplot(merged.dat,"Latency","Efficiency")
# F/B Ratio
merged.dat$f.b.ratio <- phylum.data$Firmicutes/phylum.data$Bacteroidetes
cor.matrix <-  data.frame(gender = character(), model = character(), feature = character(), pval = character(), direction = character())
for(feature in secondary_continuous){
  m3 <- lm_model(merged.dat,'f.b.ratio',c(feature,m3_covar,'Sex'))
  m3.f <- lm_model(merged.dat[female,],'f.b.ratio',c(feature,m3_covar,'Menopause_code'))
  m3.m <- lm_model(merged.dat[male,],'f.b.ratio',c(feature,m3_covar))
  m3_l <- lm_model(merged.dat,'f.b.ratio',c(feature,m3_l_covar,'Sex'))
  m3_l.f <- lm_model(merged.dat[female,],'f.b.ratio',c(feature,m3_l_covar,'Menopause_code'))
  m3_l.m <- lm_model(merged.dat[male,],'f.b.ratio',c(feature,m3_l_covar))
  m3_e <- lm_model(merged.dat,'f.b.ratio',c(feature,m3_e_covar,'Sex'))
  m3_e.f <- lm_model(merged.dat[female,],'f.b.ratio',c(feature,m3_e_covar,'Menopause_code'))
  m3_e.m <- lm_model(merged.dat[male,],'f.b.ratio',c(feature,m3_e_covar))
  m3_l_e <- lm_model(merged.dat,'f.b.ratio',c(feature,m3_l_e_covar,'Sex'))
  m3_l_e.f <- lm_model(merged.dat[female,],'f.b.ratio',c(feature,m3_l_e_covar,'Menopause_code'))
  m3_l_e.m <- lm_model(merged.dat[male,],'f.b.ratio',c(feature,m3_l_e_covar))
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3",feature,summary(m3)$coefficients[2,4],sign(summary(m3)$coefficients[2,1]))
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3",feature,summary(m3.f)$coefficients[2,4],sign(summary(m3.f)$coefficients[2,1]))
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3",feature,summary(m3.m)$coefficients[2,4],sign(summary(m3.m)$coefficients[2,1]))
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3_l",feature,summary(m3_l)$coefficients[2,4],sign(summary(m3_l)$coefficients[2,1]))
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3_l",feature,summary(m3_l.f)$coefficients[2,4],sign(summary(m3_l.f)$coefficients[2,1]))
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3_l",feature,summary(m3_l.m)$coefficients[2,4],sign(summary(m3_l.m)$coefficients[2,1]))
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3_e",feature,summary(m3_e)$coefficients[2,4],sign(summary(m3_e)$coefficients[2,1]))
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3_e",feature,summary(m3_e.f)$coefficients[2,4],sign(summary(m3_e.f)$coefficients[2,1]))
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3_e",feature,summary(m3_e.m)$coefficients[2,4],sign(summary(m3_e.m)$coefficients[2,1]))
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3_l_e",feature,summary(m3_l_e)$coefficients[2,4],sign(summary(m3_l_e)$coefficients[2,1]))
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3_l_e",feature,summary(m3_l_e.f)$coefficients[2,4],sign(summary(m3_l_e.f)$coefficients[2,1]))
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3_l_e",feature,summary(m3_l_e.m)$coefficients[2,4],sign(summary(m3_l_e.m)$coefficients[2,1]))
}
cor.matrix$pval <- as.numeric(cor.matrix$pval)
filt.feat <- cor.matrix[cor.matrix$pval<=0.05,"feature"]
cor.matrix <- cor.matrix[cor.matrix[,"feature"]%in%filt.feat,]
cor.matrix$sym <- ifelse(cor.matrix$pval>=0.05, "", ifelse(cor.matrix$pval<=0.01,ifelse(cor.matrix$pval<=0.001,"***", "**"), "*"))
cor.matrix$log <- -log10(cor.matrix$pval)*as.numeric(cor.matrix$direction)
cor.matrix$feature <- gsub("\\."," ",cor.matrix$feature)
cor.matrix$feature <- gsub("X","",cor.matrix$feature)
cor.matrix$feature <- gsub("_","-",cor.matrix$feature)
cor.matrix$gender <- factor(cor.matrix$gender, levels = c("All","Women","Men"))
pdf("Figures/300_bp_variability/SLEEP_ADJ_fb_ratio.pdf", width = 8, height = 5)
print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
        coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender) + 
        theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"), 
              legend.text = element_text(size = 12, face="bold"), legend.title = element_text(size = 12, face="bold")) + ggtitle("F/B Ratio"))
dev.off()

# Alpha Diversity analysis
merged.dat$diversity <- vegan::diversity(transformed.species.data, index = 'shannon')
cor.matrix <-  data.frame(gender = character(), model = character(), feature = character(), pval = character(), direction = character())
for(feature in secondary_continuous){
  m3 <- lm_model(merged.dat,'diversity',c(feature,m3_covar,'Sex'))
  m3.f <- lm_model(merged.dat[female,],'diversity',c(feature,m3_covar,'Menopause_code'))
  m3.m <- lm_model(merged.dat[male,],'diversity',c(feature,m3_covar))
  m3_l <- lm_model(merged.dat,'diversity',c(feature,m3_l_covar,'Sex'))
  m3_l.f <- lm_model(merged.dat[female,],'diversity',c(feature,m3_l_covar,'Menopause_code'))
  m3_l.m <- lm_model(merged.dat[male,],'diversity',c(feature,m3_l_covar))
  m3_e <- lm_model(merged.dat,'diversity',c(feature,m3_e_covar,'Sex'))
  m3_e.f <- lm_model(merged.dat[female,],'diversity',c(feature,m3_e_covar,'Menopause_code'))
  m3_e.m <- lm_model(merged.dat[male,],'diversity',c(feature,m3_e_covar))
  m3_l_e <- lm_model(merged.dat,'diversity',c(feature,m3_l_e_covar,'Sex'))
  m3_l_e.f <- lm_model(merged.dat[female,],'diversity',c(feature,m3_l_e_covar,'Menopause_code'))
  m3_l_e.m <- lm_model(merged.dat[male,],'diversity',c(feature,m3_l_e_covar))
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3",feature,summary(m3)$coefficients[2,4],sign(summary(m3)$coefficients[2,1]))
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3",feature,summary(m3.f)$coefficients[2,4],sign(summary(m3.f)$coefficients[2,1]))
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3",feature,summary(m3.m)$coefficients[2,4],sign(summary(m3.m)$coefficients[2,1]))
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3_l",feature,summary(m3_l)$coefficients[2,4],sign(summary(m3_l)$coefficients[2,1]))
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3_l",feature,summary(m3_l.f)$coefficients[2,4],sign(summary(m3_l.f)$coefficients[2,1]))
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3_l",feature,summary(m3_l.m)$coefficients[2,4],sign(summary(m3_l.m)$coefficients[2,1]))
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3_e",feature,summary(m3_e)$coefficients[2,4],sign(summary(m3_e)$coefficients[2,1]))
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3_e",feature,summary(m3_e.f)$coefficients[2,4],sign(summary(m3_e.f)$coefficients[2,1]))
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3_e",feature,summary(m3_e.m)$coefficients[2,4],sign(summary(m3_e.m)$coefficients[2,1]))
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3_l_e",feature,summary(m3_l_e)$coefficients[2,4],sign(summary(m3_l_e)$coefficients[2,1]))
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3_l_e",feature,summary(m3_l_e.f)$coefficients[2,4],sign(summary(m3_l_e.f)$coefficients[2,1]))
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3_l_e",feature,summary(m3_l_e.m)$coefficients[2,4],sign(summary(m3_l_e.m)$coefficients[2,1]))
}
cor.matrix$pval <- as.numeric(cor.matrix$pval)
filt.feat <- cor.matrix[cor.matrix$pval<=0.05,"feature"]
cor.matrix <- cor.matrix[cor.matrix[,"feature"]%in%filt.feat,]
cor.matrix$sym <- ifelse(cor.matrix$pval>=0.05, "", ifelse(cor.matrix$pval<=0.01,ifelse(cor.matrix$pval<=0.001,"***", "**"), "*"))
cor.matrix$log <- -log10(cor.matrix$pval)*as.numeric(cor.matrix$direction)
cor.matrix$feature <- gsub("\\."," ",cor.matrix$feature)
cor.matrix$feature <- gsub("X","",cor.matrix$feature)
cor.matrix$feature <- gsub("_","-",cor.matrix$feature)
cor.matrix$gender <- factor(cor.matrix$gender, levels = c("All","Women","Men"))
pdf("Figures/300_bp_variability/SLEEP_ADJ_shannon_diversity.pdf", width = 8, height = 5)
print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(name = "-log10", midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
        coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender) + 
        theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"), 
              legend.text = element_text(size = 12, face="bold"), legend.title = element_text(size = 12, face="bold")) + ggtitle("Shannon Diversity"))
dev.off()

# Permanova - Beta diversity
cor.matrix <-  data.frame(gender = character(), model = character(), feature = character(), pval = character())
for(feature in c(secondary_continuous)){
  m3 <- permanova_model(merged.dat,species.data,c(feature,m3_covar,'Sex'))
  m3.f <- permanova_model(merged.dat[female,],species.data,c(feature,m3_covar,'Menopause_code'))
  m3.m <- permanova_model(merged.dat[male,],species.data,c(feature,m3_covar))
  m3_l <- permanova_model(merged.dat,species.data,c(feature,m3_l_covar,'Sex'))
  m3_l.f <- permanova_model(merged.dat[female,],species.data,c(feature,m3_l_covar,'Menopause_code'))
  m3_l.m <- permanova_model(merged.dat[male,],species.data,c(feature,m3_l_covar))
  m3_e <- permanova_model(merged.dat,species.data,c(feature,m3_e_covar,'Sex'))
  m3_e.f <- permanova_model(merged.dat[female,],species.data,c(feature,m3_e_covar,'Menopause_code'))
  m3_e.m <- permanova_model(merged.dat[male,],species.data,c(feature,m3_e_covar))
  m3_l_e <- permanova_model(merged.dat,species.data,c(feature,m3_l_e_covar,'Sex'))
  m3_l_e.f <- permanova_model(merged.dat[female,],species.data,c(feature,m3_l_e_covar,'Menopause_code'))
  m3_l_e.m <- permanova_model(merged.dat[male,],species.data,c(feature,m3_l_e_covar))
  
  cor.matrix[nrow(cor.matrix)+1,] <-c("All","m3",feature,m3$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3",feature,m3.f$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3",feature,m3.m$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3_l",feature,m3_l$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3_l",feature,m3_l.f$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3_l",feature,m3_l.m$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3_e",feature,m3_e$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3_e",feature,m3_e.f$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3_e",feature,m3_e.m$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3_l_e",feature,m3_l_e$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3_l_e",feature,m3_l_e.f$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3_l_e",feature,m3_l_e.m$`Pr(>F)`[1])
}
cor.matrix$pval <- as.numeric(cor.matrix$pval)
filt.feat <- cor.matrix[cor.matrix$pval<=0.05,"feature"]
cor.matrix <- cor.matrix[cor.matrix[,"feature"]%in%filt.feat,]
cor.matrix$sym <- ifelse(cor.matrix$pval>=0.05, "", ifelse(cor.matrix$pval<=0.01,ifelse(cor.matrix$pval<=0.001,"***", "**"), "*"))
cor.matrix$log <- -log10(cor.matrix$pval)
cor.matrix$feature <- gsub("\\."," ",cor.matrix$feature)
cor.matrix$feature <- gsub("X","",cor.matrix$feature)
cor.matrix$feature <- gsub("_","-",cor.matrix$feature)
cor.matrix$gender <- factor(cor.matrix$gender, levels = c("All","Women","Men"))
pdf("Figures/300_bp_variability/SLEEP_ADJ_Permanova_models.pdf", width = 8, height = 5)
print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", name ="-log10" ) +
        coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender) +
        theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"),
              legend.title = element_text(size = 12, face="bold"), legend.text =  element_text(size = 12, face="bold")))
dev.off()

# Univar GM analysis
pdf("Figures/300_bp_variability/SLEEP_ADJ_BP_variability_vs_GM.pdf", width = 8, height = 5)
for(feat in secondary_continuous){
  cor.matrix <-  data.frame(gender = character(), model = character(), feature = character(), pval = character(), direction = character())
  for(species in colnames(transformed.species.data)){
    m3 <- lm_model(merged.dat,feat,c(species,m3_covar,'Sex'))
    m3.f <- lm_model(merged.dat[female,],feat,c(species,m3_covar,'Menopause_code'))
    m3.m <- lm_model(merged.dat[male,],feat,c(species,m3_covar))
    m3_l <- lm_model(merged.dat,feat,c(species,m3_l_covar,'Sex'))
    m3_l.f <- lm_model(merged.dat[female,],feat,c(species,m3_l_covar,'Menopause_code'))
    m3_l.m <- lm_model(merged.dat[male,],feat,c(species,m3_l_covar))
    m3_e <- lm_model(merged.dat,feat,c(species,m3_e_covar,'Sex'))
    m3_e.f <- lm_model(merged.dat[female,],feat,c(species,m3_e_covar,'Menopause_code'))
    m3_e.m <- lm_model(merged.dat[male,],feat,c(species,m3_e_covar))
    m3_l_e <- lm_model(merged.dat,feat,c(species,m3_l_e_covar,'Sex'))
    m3_l_e.f <- lm_model(merged.dat[female,],feat,c(species,m3_l_e_covar,'Menopause_code'))
    m3_l_e.m <- lm_model(merged.dat[male,],feat,c(species,m3_l_e_covar))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3",species,summary(m3)$coefficients[2,4],sign(summary(m3)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3",species,summary(m3.f)$coefficients[2,4],sign(summary(m3.f)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3",species,summary(m3.m)$coefficients[2,4],sign(summary(m3.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3_l",species,summary(m3_l)$coefficients[2,4],sign(summary(m3_l)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3_l",species,summary(m3_l.f)$coefficients[2,4],sign(summary(m3_l.f)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3_l",species,summary(m3_l.m)$coefficients[2,4],sign(summary(m3_l.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3_e",species,summary(m3_e)$coefficients[2,4],sign(summary(m3_e)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3_e",species,summary(m3_e.f)$coefficients[2,4],sign(summary(m3_e.f)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3_e",species,summary(m3_e.m)$coefficients[2,4],sign(summary(m3_e.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3_l_e",species,summary(m3_l_e)$coefficients[2,4],sign(summary(m3_l_e)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3_l_e",species,summary(m3_l_e.f)$coefficients[2,4],sign(summary(m3_l_e.f)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3_l_e",species,summary(m3_l_e.m)$coefficients[2,4],sign(summary(m3_l_e.m)$coefficients[2,1]))
  }
  cor.matrix$pval <- as.numeric(cor.matrix$pval)
  filt.feat <- cor.matrix[cor.matrix$pval<=0.01,"feature"]
  cor.matrix <- cor.matrix[cor.matrix[,"feature"]%in%filt.feat,]
  if(nrow(cor.matrix)==0){next}
  cor.matrix$sym <- ifelse(cor.matrix$pval>=0.05, "", ifelse(cor.matrix$pval<=0.01,ifelse(cor.matrix$pval<=0.001,"***", "**"), "*"))
  cor.matrix$log <- -log10(cor.matrix$pval)*as.numeric(cor.matrix$direction)
  cor.matrix$feature <- gsub("_"," ",cor.matrix$feature)
  title <- gsub("\\."," ",feat)
  title <- gsub("X","",title)
  title <- gsub("_","-",title)
  cor.matrix$gender <- factor(cor.matrix$gender, levels = c("All","Women","Men"))
  print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(name = "-log10", midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
          coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender)  +
          theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"),
                legend.text = element_text(size = 12, face="bold"), legend.title = element_text(size = 12, face="bold")) + ggtitle(title))
}
dev.off()


# Univar GM analysis Hyp ONLY
pdf("Figures/300_bp_variability/SLEEP_ADJ_HYP_ONLY_BP_variability_vs_GM.pdf", width = 8, height = 5)
for(feat in secondary_continuous){
  cor.matrix <-  data.frame(gender = character(), model = character(), feature = character(), pval = character(), direction = character())
  for(species in colnames(transformed.species.data)){
    m3 <- lm_model(hyp.merged.dat,feat,c(species,m3_covar,'Sex'))
    m3.f <- lm_model(hyp.merged.dat[females,],feat,c(species,m3_covar,'Menopause_code'))
    m3.m <- lm_model(hyp.merged.dat[males,],feat,c(species,m3_covar))
    m3_l <- lm_model(hyp.merged.dat,feat,c(species,m3_l_covar,'Sex'))
    m3_l.f <- lm_model(hyp.merged.dat[females,],feat,c(species,m3_l_covar,'Menopause_code'))
    m3_l.m <- lm_model(hyp.merged.dat[males,],feat,c(species,m3_l_covar))
    m3_e <- lm_model(hyp.merged.dat,feat,c(species,m3_e_covar,'Sex'))
    m3_e.f <- lm_model(hyp.merged.dat[females,],feat,c(species,m3_e_covar,'Menopause_code'))
    m3_e.m <- lm_model(hyp.merged.dat[males,],feat,c(species,m3_e_covar))
    m3_l_e <- lm_model(hyp.merged.dat,feat,c(species,m3_l_e_covar,'Sex'))
    m3_l_e.f <- lm_model(hyp.merged.dat[females,],feat,c(species,m3_l_e_covar,'Menopause_code'))
    m3_l_e.m <- lm_model(hyp.merged.dat[males,],feat,c(species,m3_l_e_covar))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3",species,summary(m3)$coefficients[2,4],sign(summary(m3)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3",species,summary(m3.f)$coefficients[2,4],sign(summary(m3.f)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3",species,summary(m3.m)$coefficients[2,4],sign(summary(m3.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3_l",species,summary(m3_l)$coefficients[2,4],sign(summary(m3_l)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3_l",species,summary(m3_l.f)$coefficients[2,4],sign(summary(m3_l.f)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3_l",species,summary(m3_l.m)$coefficients[2,4],sign(summary(m3_l.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3_e",species,summary(m3_e)$coefficients[2,4],sign(summary(m3_e)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3_e",species,summary(m3_e.f)$coefficients[2,4],sign(summary(m3_e.f)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3_e",species,summary(m3_e.m)$coefficients[2,4],sign(summary(m3_e.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3_l_e",species,summary(m3_l_e)$coefficients[2,4],sign(summary(m3_l_e)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3_l_e",species,summary(m3_l_e.f)$coefficients[2,4],sign(summary(m3_l_e.f)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3_l_e",species,summary(m3_l_e.m)$coefficients[2,4],sign(summary(m3_l_e.m)$coefficients[2,1]))
  }
  cor.matrix$pval <- as.numeric(cor.matrix$pval)
  filt.feat <- cor.matrix[cor.matrix$pval<=0.01,"feature"]
  cor.matrix <- cor.matrix[cor.matrix[,"feature"]%in%filt.feat,]
  if(nrow(cor.matrix)==0){next}
  cor.matrix$sym <- ifelse(cor.matrix$pval>=0.05, "", ifelse(cor.matrix$pval<=0.01,ifelse(cor.matrix$pval<=0.001,"***", "**"), "*"))
  cor.matrix$log <- -log10(cor.matrix$pval)*as.numeric(cor.matrix$direction)
  cor.matrix$feature <- gsub("_"," ",cor.matrix$feature)
  title <- gsub("\\."," ",feat)
  title <- gsub("X","",title)
  title <- gsub("_","-",title)
  cor.matrix$gender <- factor(cor.matrix$gender, levels = c("All","Women","Men"))
  print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(name = "-log10", midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
          coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender)  +
          theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"),
                legend.text = element_text(size = 12, face="bold"), legend.title = element_text(size = 12, face="bold")) + ggtitle(title))
}
dev.off()






###########################
### Functional analysis ###
###########################
{
# Permanova - Beta diversity
cor.matrix <-  data.frame(gender = character(), model = character(), feature = character(), pval = character())
for(feature in SBP_cov){
  m1 <- permanova_model(merged.dat,grouped.data,feature)
  m1.f <- permanova_model(merged.dat[female,],grouped.data,feature)
  m1.m <- permanova_model(merged.dat[male,],grouped.data,feature)
  m2 <- permanova_model(merged.dat,grouped.data,c(feature,m2_covar,'Sex'))
  m2.f <- permanova_model(merged.dat[female,],grouped.data,c(feature,m2_covar))
  m2.m <- permanova_model(merged.dat[male,],grouped.data,c(feature,m2_covar))
  m3 <- permanova_model(merged.dat,grouped.data,c(feature,m3_covar_sbp,'Sex'))
  m3.f <- permanova_model(merged.dat[female,],grouped.data,c(feature,m3_covar_sbp,'Menopause_code'))
  m3.m <- permanova_model(merged.dat[male,],grouped.data,c(feature,m3_covar_sbp))
  m4 <- permanova_model(merged.dat,grouped.data,c(feature,m4_covar_sbp,'Sex'))
  m4.f <- permanova_model(merged.dat[female,],grouped.data,c(feature,m4_covar_sbp,'Menopause_code'))
  m4.m <- permanova_model(merged.dat[male,],grouped.data,c(feature,m4_covar_sbp))
  cor.matrix[nrow(cor.matrix)+1,] <-c("All","m1",feature,m1$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m1",feature,m1.f$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m1",feature,m1.m$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m2",feature,m2$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m2",feature,m2.f$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m2",feature,m2.m$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3",feature,m3$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3",feature,m3.f$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3",feature,m3.m$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m4",feature,m4$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m4",feature,m4.f$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m4",feature,m4.m$`Pr(>F)`[1])
}
cor.matrix <- format_res_pnova(cor.matrix)
cor.matrix$feature <- factor(cor.matrix$feature, levels = c("Nighttime SBP CoV","Daytime SBP CoV","24-hour SBP CoV"))
pdf("Figures/300_bp_variability/FUNC_permanova_models_SBP_CoV.pdf", width = 8, height = 5)
print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", name ="-log10" ) +
        coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender) +
        theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"),
              legend.title = element_text(size = 12, face="bold"), legend.text =  element_text(size = 12, face="bold")))
dev.off()

cor.matrix <-  data.frame(gender = character(), model = character(), feature = character(), pval = character())
for(feature in SBP_surge){
  m1 <- permanova_model(merged.dat,grouped.data,feature)
  m1.f <- permanova_model(merged.dat[female,],grouped.data,feature)
  m1.m <- permanova_model(merged.dat[male,],grouped.data,feature)
  m2 <- permanova_model(merged.dat,grouped.data,c(feature,m2_covar,'Sex'))
  m2.f <- permanova_model(merged.dat[female,],grouped.data,c(feature,m2_covar))
  m2.m <- permanova_model(merged.dat[male,],grouped.data,c(feature,m2_covar))
  m3 <- permanova_model(merged.dat,grouped.data,c(feature,m3_covar_sbp,'Sex'))
  m3.f <- permanova_model(merged.dat[female,],grouped.data,c(feature,m3_covar_sbp,'Menopause_code'))
  m3.m <- permanova_model(merged.dat[male,],grouped.data,c(feature,m3_covar_sbp))
  m4 <- permanova_model(merged.dat,grouped.data,c(feature,m4_covar_sbp,'Sex'))
  m4.f <- permanova_model(merged.dat[female,],grouped.data,c(feature,m4_covar_sbp,'Menopause_code'))
  m4.m <- permanova_model(merged.dat[male,],grouped.data,c(feature,m4_covar_sbp))
  cor.matrix[nrow(cor.matrix)+1,] <-c("All","m1",feature,m1$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m1",feature,m1.f$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m1",feature,m1.m$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m2",feature,m2$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m2",feature,m2.f$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m2",feature,m2.m$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3",feature,m3$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3",feature,m3.f$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3",feature,m3.m$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m4",feature,m4$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m4",feature,m4.f$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m4",feature,m4.m$`Pr(>F)`[1])
}
cor.matrix <- format_res_pnova(cor.matrix)
pdf("Figures/300_bp_variability/FUNC_permanova_models_SBP_surge.pdf", width = 8, height = 5)
print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", name ="-log10" ) +
        coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender) +
        theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"),
              legend.title = element_text(size = 12, face="bold"), legend.text =  element_text(size = 12, face="bold")))
dev.off()

cor.matrix <-  data.frame(gender = character(), model = character(), feature = character(), pval = character())
for(feature in SBP_sd){
  m1 <- permanova_model(merged.dat,grouped.data,feature)
  m1.f <- permanova_model(merged.dat[female,],grouped.data,feature)
  m1.m <- permanova_model(merged.dat[male,],grouped.data,feature)
  m2 <- permanova_model(merged.dat,grouped.data,c(feature,m2_covar,'Sex'))
  m2.f <- permanova_model(merged.dat[female,],grouped.data,c(feature,m2_covar))
  m2.m <- permanova_model(merged.dat[male,],grouped.data,c(feature,m2_covar))
  m3 <- permanova_model(merged.dat,grouped.data,c(feature,m3_covar_sbp,'Sex'))
  m3.f <- permanova_model(merged.dat[female,],grouped.data,c(feature,m3_covar_sbp,'Menopause_code'))
  m3.m <- permanova_model(merged.dat[male,],grouped.data,c(feature,m3_covar_sbp))
  m4 <- permanova_model(merged.dat,grouped.data,c(feature,m4_covar_sbp,'Sex'))
  m4.f <- permanova_model(merged.dat[female,],grouped.data,c(feature,m4_covar_sbp,'Menopause_code'))
  m4.m <- permanova_model(merged.dat[male,],grouped.data,c(feature,m4_covar_sbp))
  cor.matrix[nrow(cor.matrix)+1,] <-c("All","m1",feature,m1$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m1",feature,m1.f$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m1",feature,m1.m$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m2",feature,m2$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m2",feature,m2.f$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m2",feature,m2.m$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3",feature,m3$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3",feature,m3.f$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3",feature,m3.m$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m4",feature,m4$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m4",feature,m4.f$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m4",feature,m4.m$`Pr(>F)`[1])
}
cor.matrix <- format_res_pnova(cor.matrix)
cor.matrix$feature <- factor(cor.matrix$feature, levels = c("Nighttime SBP SD","Daytime SBP SD","24-hour SBP SD"))
pdf("Figures/300_bp_variability/FUNC_permanova_models_SBP_SD.pdf", width = 8, height = 5)
print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", name ="-log10" ) +
        coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender) +
        theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"),
              legend.title = element_text(size = 12, face="bold"), legend.text =  element_text(size = 12, face="bold")))
dev.off()

cor.matrix <-  data.frame(gender = character(), model = character(), feature = character(), pval = character())
for(feature in DBP_cov){
  m1 <- permanova_model(merged.dat,grouped.data,feature)
  m1.f <- permanova_model(merged.dat[female,],grouped.data,feature)
  m1.m <- permanova_model(merged.dat[male,],grouped.data,feature)
  m2 <- permanova_model(merged.dat,grouped.data,c(feature,m2_covar,'Sex'))
  m2.f <- permanova_model(merged.dat[female,],grouped.data,c(feature,m2_covar))
  m2.m <- permanova_model(merged.dat[male,],grouped.data,c(feature,m2_covar))
  m3 <- permanova_model(merged.dat,grouped.data,c(feature,m3_covar_dbp,'Sex'))
  m3.f <- permanova_model(merged.dat[female,],grouped.data,c(feature,m3_covar_dbp,'Menopause_code'))
  m3.m <- permanova_model(merged.dat[male,],grouped.data,c(feature,m3_covar_dbp))
  m4 <- permanova_model(merged.dat,grouped.data,c(feature,m4_covar_dbp,'Sex'))
  m4.f <- permanova_model(merged.dat[female,],grouped.data,c(feature,m4_covar_dbp,'Menopause_code'))
  m4.m <- permanova_model(merged.dat[male,],grouped.data,c(feature,m4_covar_dbp))
  cor.matrix[nrow(cor.matrix)+1,] <-c("All","m1",feature,m1$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m1",feature,m1.f$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m1",feature,m1.m$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m2",feature,m2$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m2",feature,m2.f$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m2",feature,m2.m$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3",feature,m3$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3",feature,m3.f$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3",feature,m3.m$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m4",feature,m4$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m4",feature,m4.f$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m4",feature,m4.m$`Pr(>F)`[1])
}
cor.matrix <- format_res_pnova(cor.matrix)
cor.matrix$feature <- factor(cor.matrix$feature, levels = c("Nighttime DBP CoV","Daytime DBP CoV","24-hour DBP CoV"))
pdf("Figures/300_bp_variability/FUNC_permanova_models_DBP_CoV.pdf", width = 8, height = 5)
print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", name ="-log10" ) +
        coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender) +
        theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"),
              legend.title = element_text(size = 12, face="bold"), legend.text =  element_text(size = 12, face="bold")))
dev.off()

cor.matrix <-  data.frame(gender = character(), model = character(), feature = character(), pval = character())
for(feature in DBP_sd){
  m1 <- permanova_model(merged.dat,grouped.data,feature)
  m1.f <- permanova_model(merged.dat[female,],grouped.data,feature)
  m1.m <- permanova_model(merged.dat[male,],grouped.data,feature)
  m2 <- permanova_model(merged.dat,grouped.data,c(feature,m2_covar,'Sex'))
  m2.f <- permanova_model(merged.dat[female,],grouped.data,c(feature,m2_covar))
  m2.m <- permanova_model(merged.dat[male,],grouped.data,c(feature,m2_covar))
  m3 <- permanova_model(merged.dat,grouped.data,c(feature,m3_covar_dbp,'Sex'))
  m3.f <- permanova_model(merged.dat[female,],grouped.data,c(feature,m3_covar_dbp,'Menopause_code'))
  m3.m <- permanova_model(merged.dat[male,],grouped.data,c(feature,m3_covar_dbp))
  m4 <- permanova_model(merged.dat,grouped.data,c(feature,m4_covar_dbp,'Sex'))
  m4.f <- permanova_model(merged.dat[female,],grouped.data,c(feature,m4_covar_dbp,'Menopause_code'))
  m4.m <- permanova_model(merged.dat[male,],grouped.data,c(feature,m4_covar_dbp))
  cor.matrix[nrow(cor.matrix)+1,] <-c("All","m1",feature,m1$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m1",feature,m1.f$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m1",feature,m1.m$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m2",feature,m2$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m2",feature,m2.f$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m2",feature,m2.m$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3",feature,m3$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3",feature,m3.f$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3",feature,m3.m$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m4",feature,m4$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m4",feature,m4.f$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m4",feature,m4.m$`Pr(>F)`[1])
}
cor.matrix <- format_res_pnova(cor.matrix)
cor.matrix$feature <- factor(cor.matrix$feature, levels = c("Nighttime DBP SD","Daytime DBP SD","24-hour DBP SD"))
pdf("Figures/300_bp_variability/FUNC_permanova_models_DBP_SD.pdf", width = 8, height = 5)
print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", name ="-log10" ) +
        coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender) +
        theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"),
              legend.title = element_text(size = 12, face="bold"), legend.text =  element_text(size = 12, face="bold")))
dev.off()

# Correlation analysis
pdf("Figures/300_bp_variability/BP_variability_vs_FUNC.pdf", width = 10, height = 5)
for(feat in c(SBP_cov,SBP_sd,SBP_surge)){
  cor.matrix <-  data.frame(gender = character(), model = character(), feature = character(), pval = character(), direction = character())
  for(species in colnames(transformed.grouped.data)){
    m1 <- lm_model(merged.dat,feat,species)
    m1.f <- lm_model(merged.dat[female,],feat,species)
    m1.m <- lm_model(merged.dat[male,],feat,species)
    m2 <- lm_model(merged.dat,feat,c(species,m2_covar,'Sex'))
    m2.f <- lm_model(merged.dat[female,],feat,c(species,m2_covar))
    m2.m <- lm_model(merged.dat[male,],feat,c(species,m2_covar))
    m3 <- lm_model(merged.dat,feat,c(species,m3_covar_sbp,'Sex'))
    m3.f <- lm_model(merged.dat[female,],feat,c(species,m3_covar_sbp,'Menopause_code'))
    m3.m <- lm_model(merged.dat[male,],feat,c(species,m3_covar_sbp))
    m4 <- lm_model(merged.dat,feat,c(species,m4_covar_sbp,'Sex'))
    m4.f <- lm_model(merged.dat[female,],feat,c(species,m4_covar_sbp,'Menopause_code'))
    m4.m <- lm_model(merged.dat[male,],feat,c(species,m4_covar_sbp))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m1",species,summary(m1)$coefficients[2,4],sign(summary(m1)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m1",species,summary(m1.f)$coefficients[2,4],sign(summary(m1.f)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m1",species,summary(m1.m)$coefficients[2,4],sign(summary(m1.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m2",species,summary(m2)$coefficients[2,4],sign(summary(m2)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m2",species,summary(m2.f)$coefficients[2,4],sign(summary(m2.f)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m2",species,summary(m2.m)$coefficients[2,4],sign(summary(m2.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3",species,summary(m3)$coefficients[2,4],sign(summary(m3)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3",species,summary(m3.f)$coefficients[2,4],sign(summary(m3.f)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3",species,summary(m3.m)$coefficients[2,4],sign(summary(m3.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m4",species,summary(m4)$coefficients[2,4],sign(summary(m4)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m4",species,summary(m4.f)$coefficients[2,4],sign(summary(m4.f)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m4",species,summary(m4.m)$coefficients[2,4],sign(summary(m4.m)$coefficients[2,1]))
  }
  cor.matrix$pval <- as.numeric(cor.matrix$pval)
  filt.feat <- cor.matrix[cor.matrix$pval<=0.01,"feature"]
  cor.matrix <- cor.matrix[cor.matrix[,"feature"]%in%filt.feat,]
  if(nrow(cor.matrix)==0){next}
  cor.matrix$sym <- ifelse(cor.matrix$pval>=0.05, "", ifelse(cor.matrix$pval<=0.01,ifelse(cor.matrix$pval<=0.001,"***", "**"), "*"))
  cor.matrix$log <- -log10(cor.matrix$pval)*as.numeric(cor.matrix$direction)
  cor.matrix$feature <- gsub("_"," ",cor.matrix$feature)
  title <- gsub("\\."," ",feat)
  title <- gsub("X","",title)
  title <- gsub("_","-",title)
  cor.matrix$gender <- factor(cor.matrix$gender, levels = c("All","Women","Men"))
  print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(name = "-log10", midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
          coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender)  +
          theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"),
                legend.text = element_text(size = 12, face="bold"), legend.title = element_text(size = 12, face="bold")) + ggtitle(title))
}
for(feat in c(DBP_cov,DBP_sd)){
  cor.matrix <-  data.frame(gender = character(), model = character(), feature = character(), pval = character(), direction = character())
  for(species in colnames(transformed.grouped.data)){
    m1 <- lm_model(merged.dat,feat,species)
    m1.f <- lm_model(merged.dat[female,],feat,species)
    m1.m <- lm_model(merged.dat[male,],feat,species)
    m2 <- lm_model(merged.dat,feat,c(species,m2_covar,'Sex'))
    m2.f <- lm_model(merged.dat[female,],feat,c(species,m2_covar))
    m2.m <- lm_model(merged.dat[male,],feat,c(species,m2_covar))
    m3 <- lm_model(merged.dat,feat,c(species,m3_covar_dbp,'Sex'))
    m3.f <- lm_model(merged.dat[female,],feat,c(species,m3_covar_dbp,'Menopause_code'))
    m3.m <- lm_model(merged.dat[male,],feat,c(species,m3_covar_dbp))
    m4 <- lm_model(merged.dat,feat,c(species,m4_covar_dbp,'Sex'))
    m4.f <- lm_model(merged.dat[female,],feat,c(species,m4_covar_dbp,'Menopause_code'))
    m4.m <- lm_model(merged.dat[male,],feat,c(species,m4_covar_dbp))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m1",species,summary(m1)$coefficients[2,4],sign(summary(m1)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m1",species,summary(m1.f)$coefficients[2,4],sign(summary(m1.f)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m1",species,summary(m1.m)$coefficients[2,4],sign(summary(m1.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m2",species,summary(m2)$coefficients[2,4],sign(summary(m2)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m2",species,summary(m2.f)$coefficients[2,4],sign(summary(m2.f)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m2",species,summary(m2.m)$coefficients[2,4],sign(summary(m2.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3",species,summary(m3)$coefficients[2,4],sign(summary(m3)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3",species,summary(m3.f)$coefficients[2,4],sign(summary(m3.f)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3",species,summary(m3.m)$coefficients[2,4],sign(summary(m3.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m4",species,summary(m4)$coefficients[2,4],sign(summary(m4)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m4",species,summary(m4.f)$coefficients[2,4],sign(summary(m4.f)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m4",species,summary(m4.m)$coefficients[2,4],sign(summary(m4.m)$coefficients[2,1]))
  }
  cor.matrix$pval <- as.numeric(cor.matrix$pval)
  filt.feat <- cor.matrix[cor.matrix$pval<=0.01,"feature"]
  cor.matrix <- cor.matrix[cor.matrix[,"feature"]%in%filt.feat,]
  if(nrow(cor.matrix)==0){next}
  cor.matrix$sym <- ifelse(cor.matrix$pval>=0.05, "", ifelse(cor.matrix$pval<=0.01,ifelse(cor.matrix$pval<=0.001,"***", "**"), "*"))
  cor.matrix$log <- -log10(cor.matrix$pval)*as.numeric(cor.matrix$direction)
  cor.matrix$feature <- gsub("_"," ",cor.matrix$feature)
  title <- gsub("\\."," ",feat)
  title <- gsub("X","",title)
  title <- gsub("_","-",title)
  cor.matrix$gender <- factor(cor.matrix$gender, levels = c("All","Women","Men"))
  print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(name = "-log10", midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
          coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender)  +
          theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"),
                legend.text = element_text(size = 12, face="bold"), legend.title = element_text(size = 12, face="bold")) + ggtitle(title))
}
dev.off()

pdf("Figures/300_bp_variability/Dipping_vs_FUNC_boxplot.pdf", width = 10, height = 8)
for(feat in colnames(transformed.grouped.data)){
  my_boxplot(dip_subset,feat,"Dipping.classification",ylab = feat)
}
dev.off()

###############################################
### BP Variability in Hypertensive patients ###
###############################################
#GM univariate
pdf("Figures/300_bp_variability/HYP_ONLY_BP_variability_vs_FUNC.pdf", width = 10, height = 8)
for(feat in SBP_cov){
  cor.matrix <-  data.frame(gender = character(), model = character(), feature = character(), pval = character(), direction = character())
  for(species in colnames(transformed.grouped.data)){
    m1 <- lm_model(hyp.merged.dat,feat,species)
    #m1.f <- lm_model(hyp.merged.dat[females,],feat,species)
    m1.m <- lm_model(hyp.merged.dat[males,],feat,species)
    m2 <- lm_model(hyp.merged.dat,feat,c(species,m2_covar,'Sex'))
    #m2.f <- lm_model(hyp.merged.dat[females,],feat,c(species,m2_covar))
    m2.m <- lm_model(hyp.merged.dat[males,],feat,c(species,m2_covar))
    m3 <- lm_model(hyp.merged.dat,feat,c(species,m3_covar_sbp,'Sex'))
    #m3.f <- lm_model(hyp.merged.dat[females,],feat,c(species,m3_covar_sbp,'Menopause_code'))
    m3.m <- lm_model(hyp.merged.dat[males,],feat,c(species,m3_covar_sbp))
    m4 <- lm_model(hyp.merged.dat,feat,c(species,m4_covar,'Sex'))
    #m4.f <- lm_model(hyp.merged.dat[females,],feat,c(species,m4_covar,'Menopause_code'))
    m4.m <- lm_model(hyp.merged.dat[males,],feat,c(species,m4_covar))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m1",species,summary(m1)$coefficients[2,4],sign(summary(m1)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m1",species,1,1)
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m1",species,summary(m1.m)$coefficients[2,4],sign(summary(m1.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m2",species,summary(m2)$coefficients[2,4],sign(summary(m2)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m2",species,1,1)
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m2",species,summary(m2.m)$coefficients[2,4],sign(summary(m2.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3",species,summary(m3)$coefficients[2,4],sign(summary(m3)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3",species,1,1)
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3",species,summary(m3.m)$coefficients[2,4],sign(summary(m3.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m4",species,summary(m4)$coefficients[2,4],sign(summary(m4)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m4",species,1,1)
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m4",species,summary(m4.m)$coefficients[2,4],sign(summary(m4.m)$coefficients[2,1]))
  }
  cor.matrix$pval <- as.numeric(cor.matrix$pval)
  filt.feat <- cor.matrix[cor.matrix$pval<=0.01,"feature"]
  cor.matrix <- cor.matrix[cor.matrix[,"feature"]%in%filt.feat,]
  if(nrow(cor.matrix)==0){next}
  cor.matrix$sym <- ifelse(cor.matrix$pval>=0.05, "", ifelse(cor.matrix$pval<=0.01,ifelse(cor.matrix$pval<=0.001,"***", "**"), "*"))
  cor.matrix$log <- -log10(cor.matrix$pval)*as.numeric(cor.matrix$direction)
  cor.matrix$feature <- gsub("_"," ",cor.matrix$feature)
  title <- gsub("\\."," ",feat)
  title <- gsub("X","",title)
  title <- gsub("_","-",title)
  cor.matrix$gender <- factor(cor.matrix$gender, levels = c("All","Women","Men"))
  print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(name = "-log10", midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
          coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender)  +
          theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"),
                legend.text = element_text(size = 12, face="bold"), legend.title = element_text(size = 12, face="bold")) + ggtitle(title))
}
for(feat in SBP_sd){
  cor.matrix <-  data.frame(gender = character(), model = character(), feature = character(), pval = character(), direction = character())
  for(species in colnames(transformed.grouped.data)){
    m1 <- lm_model(hyp.merged.dat,feat,species)
    m1.f <- lm_model(hyp.merged.dat[females,],feat,species)
    m1.m <- lm_model(hyp.merged.dat[males,],feat,species)
    m2 <- lm_model(hyp.merged.dat,feat,c(species,m2_covar,'Sex'))
    m2.f <- lm_model(hyp.merged.dat[females,],feat,c(species,m2_covar))
    m2.m <- lm_model(hyp.merged.dat[males,],feat,c(species,m2_covar))
    m3 <- lm_model(hyp.merged.dat,feat,c(species,m3_covar_sbp,'Sex'))
    m3.f <- lm_model(hyp.merged.dat[females,],feat,c(species,m3_covar_sbp,'Menopause_code'))
    m3.m <- lm_model(hyp.merged.dat[males,],feat,c(species,m3_covar_sbp))
    m4 <- lm_model(hyp.merged.dat,feat,c(species,m4_covar,'Sex'))
    m4.f <- lm_model(hyp.merged.dat[females,],feat,c(species,m4_covar,'Menopause_code'))
    m4.m <- lm_model(hyp.merged.dat[males,],feat,c(species,m4_covar))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m1",species,summary(m1)$coefficients[2,4],sign(summary(m1)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m1",species,summary(m1.f)$coefficients[2,4],sign(summary(m1.f)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m1",species,summary(m1.m)$coefficients[2,4],sign(summary(m1.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m2",species,summary(m2)$coefficients[2,4],sign(summary(m2)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m2",species,summary(m2.f)$coefficients[2,4],sign(summary(m2.f)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m2",species,summary(m2.m)$coefficients[2,4],sign(summary(m2.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3",species,summary(m3)$coefficients[2,4],sign(summary(m3)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3",species,summary(m3.f)$coefficients[2,4],sign(summary(m3.f)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3",species,summary(m3.m)$coefficients[2,4],sign(summary(m3.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m4",species,summary(m4)$coefficients[2,4],sign(summary(m4)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m4",species,summary(m4.f)$coefficients[2,4],sign(summary(m4.f)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m4",species,summary(m4.m)$coefficients[2,4],sign(summary(m4.m)$coefficients[2,1]))
  }
  cor.matrix$pval <- as.numeric(cor.matrix$pval)
  filt.feat <- cor.matrix[cor.matrix$pval<=0.01,"feature"]
  cor.matrix <- cor.matrix[cor.matrix[,"feature"]%in%filt.feat,]
  if(nrow(cor.matrix)==0){next}
  cor.matrix$sym <- ifelse(cor.matrix$pval>=0.05, "", ifelse(cor.matrix$pval<=0.01,ifelse(cor.matrix$pval<=0.001,"***", "**"), "*"))
  cor.matrix$log <- -log10(cor.matrix$pval)*as.numeric(cor.matrix$direction)
  cor.matrix$feature <- gsub("_"," ",cor.matrix$feature)
  title <- gsub("\\."," ",feat)
  title <- gsub("X","",title)
  title <- gsub("_","-",title)
  cor.matrix$gender <- factor(cor.matrix$gender, levels = c("All","Women","Men"))
  print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(name = "-log10", midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
          coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender)  +
          theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"),
                legend.text = element_text(size = 12, face="bold"), legend.title = element_text(size = 12, face="bold")) + ggtitle(title))
}
for(feat in DBP_cov){
  cor.matrix <-  data.frame(gender = character(), model = character(), feature = character(), pval = character(), direction = character())
  for(species in colnames(transformed.grouped.data)){
    m1 <- lm_model(hyp.merged.dat,feat,species)
    m1.f <- lm_model(hyp.merged.dat[females,],feat,species)
    m1.m <- lm_model(hyp.merged.dat[males,],feat,species)
    m2 <- lm_model(hyp.merged.dat,feat,c(species,m2_covar,'Sex'))
    m2.f <- lm_model(hyp.merged.dat[females,],feat,c(species,m2_covar))
    m2.m <- lm_model(hyp.merged.dat[males,],feat,c(species,m2_covar))
    m3 <- lm_model(hyp.merged.dat,feat,c(species,m3_covar_dbp,'Sex'))
    m3.f <- lm_model(hyp.merged.dat[females,],feat,c(species,m3_covar_dbp,'Menopause_code'))
    m3.m <- lm_model(hyp.merged.dat[males,],feat,c(species,m3_covar_dbp))
    m4 <- lm_model(hyp.merged.dat,feat,c(species,m4_covar,'Sex'))
    m4.f <- lm_model(hyp.merged.dat[females,],feat,c(species,m4_covar,'Menopause_code'))
    m4.m <- lm_model(hyp.merged.dat[males,],feat,c(species,m4_covar))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m1",species,summary(m1)$coefficients[2,4],sign(summary(m1)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m1",species,summary(m1.f)$coefficients[2,4],sign(summary(m1.f)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m1",species,summary(m1.m)$coefficients[2,4],sign(summary(m1.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m2",species,summary(m2)$coefficients[2,4],sign(summary(m2)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m2",species,summary(m2.f)$coefficients[2,4],sign(summary(m2.f)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m2",species,summary(m2.m)$coefficients[2,4],sign(summary(m2.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3",species,summary(m3)$coefficients[2,4],sign(summary(m3)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3",species,summary(m3.f)$coefficients[2,4],sign(summary(m3.f)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3",species,summary(m3.m)$coefficients[2,4],sign(summary(m3.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m4",species,summary(m4)$coefficients[2,4],sign(summary(m4)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m4",species,summary(m4.f)$coefficients[2,4],sign(summary(m4.f)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m4",species,summary(m4.m)$coefficients[2,4],sign(summary(m4.m)$coefficients[2,1]))
  }
  cor.matrix$pval <- as.numeric(cor.matrix$pval)
  filt.feat <- cor.matrix[cor.matrix$pval<=0.01,"feature"]
  cor.matrix <- cor.matrix[cor.matrix[,"feature"]%in%filt.feat,]
  if(nrow(cor.matrix)==0){next}
  cor.matrix$sym <- ifelse(cor.matrix$pval>=0.05, "", ifelse(cor.matrix$pval<=0.01,ifelse(cor.matrix$pval<=0.001,"***", "**"), "*"))
  cor.matrix$log <- -log10(cor.matrix$pval)*as.numeric(cor.matrix$direction)
  cor.matrix$feature <- gsub("_"," ",cor.matrix$feature)
  title <- gsub("\\."," ",feat)
  title <- gsub("X","",title)
  title <- gsub("_","-",title)
  cor.matrix$gender <- factor(cor.matrix$gender, levels = c("All","Women","Men"))
  print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(name = "-log10", midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
          coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender)  +
          theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"),
                legend.text = element_text(size = 12, face="bold"), legend.title = element_text(size = 12, face="bold")) + ggtitle(title))
}
for(feat in DBP_sd){
  cor.matrix <-  data.frame(gender = character(), model = character(), feature = character(), pval = character(), direction = character())
  for(species in colnames(transformed.grouped.data)){
    m1 <- lm_model(hyp.merged.dat,feat,species)
    m1.f <- lm_model(hyp.merged.dat[females,],feat,species)
    m1.m <- lm_model(hyp.merged.dat[males,],feat,species)
    m2 <- lm_model(hyp.merged.dat,feat,c(species,m2_covar,'Sex'))
    m2.f <- lm_model(hyp.merged.dat[females,],feat,c(species,m2_covar))
    m2.m <- lm_model(hyp.merged.dat[males,],feat,c(species,m2_covar))
    m3 <- lm_model(hyp.merged.dat,feat,c(species,m3_covar_dbp,'Sex'))
    m3.f <- lm_model(hyp.merged.dat[females,],feat,c(species,m3_covar_dbp,'Menopause_code'))
    m3.m <- lm_model(hyp.merged.dat[males,],feat,c(species,m3_covar_dbp))
    m4 <- lm_model(hyp.merged.dat,feat,c(species,m4_covar,'Sex'))
    m4.f <- lm_model(hyp.merged.dat[females,],feat,c(species,m4_covar,'Menopause_code'))
    m4.m <- lm_model(hyp.merged.dat[males,],feat,c(species,m4_covar))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m1",species,summary(m1)$coefficients[2,4],sign(summary(m1)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m1",species,summary(m1.f)$coefficients[2,4],sign(summary(m1.f)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m1",species,summary(m1.m)$coefficients[2,4],sign(summary(m1.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m2",species,summary(m2)$coefficients[2,4],sign(summary(m2)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m2",species,summary(m2.f)$coefficients[2,4],sign(summary(m2.f)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m2",species,summary(m2.m)$coefficients[2,4],sign(summary(m2.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3",species,summary(m3)$coefficients[2,4],sign(summary(m3)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3",species,summary(m3.f)$coefficients[2,4],sign(summary(m3.f)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3",species,summary(m3.m)$coefficients[2,4],sign(summary(m3.m)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("All","m4",species,summary(m4)$coefficients[2,4],sign(summary(m4)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m4",species,summary(m4.f)$coefficients[2,4],sign(summary(m4.f)$coefficients[2,1]))
    cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m4",species,summary(m4.m)$coefficients[2,4],sign(summary(m4.m)$coefficients[2,1]))
  }
  cor.matrix$pval <- as.numeric(cor.matrix$pval)
  filt.feat <- cor.matrix[cor.matrix$pval<=0.01,"feature"]
  cor.matrix <- cor.matrix[cor.matrix[,"feature"]%in%filt.feat,]
  if(nrow(cor.matrix)==0){next}
  cor.matrix$sym <- ifelse(cor.matrix$pval>=0.05, "", ifelse(cor.matrix$pval<=0.01,ifelse(cor.matrix$pval<=0.001,"***", "**"), "*"))
  cor.matrix$log <- -log10(cor.matrix$pval)*as.numeric(cor.matrix$direction)
  cor.matrix$feature <- gsub("_"," ",cor.matrix$feature)
  title <- gsub("\\."," ",feat)
  title <- gsub("X","",title)
  title <- gsub("_","-",title)
  cor.matrix$gender <- factor(cor.matrix$gender, levels = c("All","Women","Men"))
  print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(name = "-log10", midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
          coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender)  +
          theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"),
                legend.text = element_text(size = 12, face="bold"), legend.title = element_text(size = 12, face="bold")) + ggtitle(title))
}
dev.off()

# Permanova - Beta diversity
cor.matrix <-  data.frame(gender = character(), model = character(), feature = character(), pval = character())
for(feature in SBP_cov){
  m1 <- permanova_model(hyp.merged.dat,grouped.data,feature)
  m1.f <- permanova_model(hyp.merged.dat[females,],grouped.data,feature)
  m1.m <- permanova_model(hyp.merged.dat[males,],grouped.data,feature)
  m2 <- permanova_model(hyp.merged.dat,grouped.data,c(feature,m2_covar,'Sex'))
  m2.f <- permanova_model(hyp.merged.dat[females,],grouped.data,c(feature,m2_covar))
  m2.m <- permanova_model(hyp.merged.dat[males,],grouped.data,c(feature,m2_covar))
  m3 <- permanova_model(hyp.merged.dat,grouped.data,c(feature,m3_covar_sbp,'Sex'))
  m3.f <- permanova_model(hyp.merged.dat[females,],grouped.data,c(feature,m3_covar_sbp,'Menopause_code'))
  m3.m <- permanova_model(hyp.merged.dat[males,],grouped.data,c(feature,m3_covar_sbp))
  m4 <- permanova_model(hyp.merged.dat,grouped.data,c(feature,m4_covar,'Sex'))
  m4.f <- permanova_model(hyp.merged.dat[females,],grouped.data,c(feature,m4_covar,'Menopause_code'))
  m4.m <- permanova_model(hyp.merged.dat[males,],grouped.data,c(feature,m4_covar))
  cor.matrix[nrow(cor.matrix)+1,] <-c("All","m1",feature,m1$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m1",feature,m1.f$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m1",feature,m1.m$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m2",feature,m2$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m2",feature,m2.f$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m2",feature,m2.m$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3",feature,m3$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3",feature,m3.f$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3",feature,m3.m$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m4",feature,m4$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m4",feature,m4.f$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m4",feature,m4.m$`Pr(>F)`[1])
}
cor.matrix <- format_res_pnova(cor.matrix)
pdf("Figures/300_bp_variability/HYP_ONLY_SBP_CoV_permanova_models_FUNC.pdf", width = 8, height = 5)
print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", name ="-log10" ) +
        coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender) +
        theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"),
              legend.title = element_text(size = 12, face="bold"), legend.text =  element_text(size = 12, face="bold")))
dev.off()

cor.matrix <-  data.frame(gender = character(), model = character(), feature = character(), pval = character())
for(feature in SBP_sd){
  m1 <- permanova_model(hyp.merged.dat,grouped.data,feature)
  m1.f <- permanova_model(hyp.merged.dat[females,],grouped.data,feature)
  m1.m <- permanova_model(hyp.merged.dat[males,],grouped.data,feature)
  m2 <- permanova_model(hyp.merged.dat,grouped.data,c(feature,m2_covar,'Sex'))
  m2.f <- permanova_model(hyp.merged.dat[females,],grouped.data,c(feature,m2_covar))
  m2.m <- permanova_model(hyp.merged.dat[males,],grouped.data,c(feature,m2_covar))
  m3 <- permanova_model(hyp.merged.dat,grouped.data,c(feature,m3_covar_sbp,'Sex'))
  m3.f <- permanova_model(hyp.merged.dat[females,],grouped.data,c(feature,m3_covar_sbp,'Menopause_code'))
  m3.m <- permanova_model(hyp.merged.dat[males,],grouped.data,c(feature,m3_covar_sbp))
  m4 <- permanova_model(hyp.merged.dat,grouped.data,c(feature,m4_covar,'Sex'))
  m4.f <- permanova_model(hyp.merged.dat[females,],grouped.data,c(feature,m4_covar,'Menopause_code'))
  m4.m <- permanova_model(hyp.merged.dat[males,],grouped.data,c(feature,m4_covar))
  cor.matrix[nrow(cor.matrix)+1,] <-c("All","m1",feature,m1$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m1",feature,m1.f$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m1",feature,m1.m$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m2",feature,m2$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m2",feature,m2.f$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m2",feature,m2.m$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3",feature,m3$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3",feature,m3.f$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3",feature,m3.m$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m4",feature,m4$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m4",feature,m4.f$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m4",feature,m4.m$`Pr(>F)`[1])
}
cor.matrix <- format_res_pnova(cor.matrix)
pdf("Figures/300_bp_variability/HYP_ONLY_SBP_SD_permanova_models_FUNC.pdf", width = 8, height = 5)
print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", name ="-log10" ) +
        coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender) +
        theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"),
              legend.title = element_text(size = 12, face="bold"), legend.text =  element_text(size = 12, face="bold")))
dev.off()

cor.matrix <-  data.frame(gender = character(), model = character(), feature = character(), pval = character())
for(feature in DBP_cov){
  m1 <- permanova_model(hyp.merged.dat,grouped.data,feature)
  m1.f <- permanova_model(hyp.merged.dat[females,],grouped.data,feature)
  m1.m <- permanova_model(hyp.merged.dat[males,],grouped.data,feature)
  m2 <- permanova_model(hyp.merged.dat,grouped.data,c(feature,m2_covar,'Sex'))
  m2.f <- permanova_model(hyp.merged.dat[females,],grouped.data,c(feature,m2_covar))
  m2.m <- permanova_model(hyp.merged.dat[males,],grouped.data,c(feature,m2_covar))
  m3 <- permanova_model(hyp.merged.dat,grouped.data,c(feature,m3_covar_dbp,'Sex'))
  m3.f <- permanova_model(hyp.merged.dat[females,],grouped.data,c(feature,m3_covar_dbp,'Menopause_code'))
  m3.m <- permanova_model(hyp.merged.dat[males,],grouped.data,c(feature,m3_covar_dbp))
  m4 <- permanova_model(hyp.merged.dat,grouped.data,c(feature,m4_covar,'Sex'))
  m4.f <- permanova_model(hyp.merged.dat[females,],grouped.data,c(feature,m4_covar,'Menopause_code'))
  m4.m <- permanova_model(hyp.merged.dat[males,],grouped.data,c(feature,m4_covar))
  cor.matrix[nrow(cor.matrix)+1,] <-c("All","m1",feature,m1$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m1",feature,m1.f$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m1",feature,m1.m$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m2",feature,m2$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m2",feature,m2.f$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m2",feature,m2.m$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3",feature,m3$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3",feature,m3.f$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3",feature,m3.m$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m4",feature,m4$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m4",feature,m4.f$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m4",feature,m4.m$`Pr(>F)`[1])
}
cor.matrix <- format_res_pnova(cor.matrix)
pdf("Figures/300_bp_variability/HYP_ONLY_DBP_CoV_permanova_models_FUNC.pdf", width = 8, height = 5)
print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", name ="-log10" ) +
        coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender) +
        theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"),
              legend.title = element_text(size = 12, face="bold"), legend.text =  element_text(size = 12, face="bold")))
dev.off()

cor.matrix <-  data.frame(gender = character(), model = character(), feature = character(), pval = character())
for(feature in DBP_sd){
  m1 <- permanova_model(hyp.merged.dat,grouped.data,feature)
  m1.f <- permanova_model(hyp.merged.dat[females,],grouped.data,feature)
  m1.m <- permanova_model(hyp.merged.dat[males,],grouped.data,feature)
  m2 <- permanova_model(hyp.merged.dat,grouped.data,c(feature,m2_covar,'Sex'))
  m2.f <- permanova_model(hyp.merged.dat[females,],grouped.data,c(feature,m2_covar))
  m2.m <- permanova_model(hyp.merged.dat[males,],grouped.data,c(feature,m2_covar))
  m3 <- permanova_model(hyp.merged.dat,grouped.data,c(feature,m3_covar_dbp,'Sex'))
  m3.f <- permanova_model(hyp.merged.dat[females,],grouped.data,c(feature,m3_covar_dbp,'Menopause_code'))
  m3.m <- permanova_model(hyp.merged.dat[males,],grouped.data,c(feature,m3_covar_dbp))
  m4 <- permanova_model(hyp.merged.dat,grouped.data,c(feature,m4_covar,'Sex'))
  m4.f <- permanova_model(hyp.merged.dat[females,],grouped.data,c(feature,m4_covar,'Menopause_code'))
  m4.m <- permanova_model(hyp.merged.dat[males,],grouped.data,c(feature,m4_covar))
  cor.matrix[nrow(cor.matrix)+1,] <-c("All","m1",feature,m1$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m1",feature,m1.f$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m1",feature,m1.m$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m2",feature,m2$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m2",feature,m2.f$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m2",feature,m2.m$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3",feature,m3$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3",feature,m3.f$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3",feature,m3.m$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m4",feature,m4$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m4",feature,m4.f$`Pr(>F)`[1])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m4",feature,m4.m$`Pr(>F)`[1])
}
cor.matrix <- format_res_pnova(cor.matrix)
pdf("Figures/300_bp_variability/HYP_ONLY_DBP_SD_permanova_models_FUNC.pdf", width = 8, height = 5)
print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", name ="-log10" ) +
        coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender) +
        theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"),
              legend.title = element_text(size = 12, face="bold"), legend.text =  element_text(size = 12, face="bold")))
dev.off()
}


########################
### Mediation Models ###
########################
library(mediation)

###Bacteroides dorei --> Blood Isobutyric --> SBP Covar
#Total Effect IV on DV
fit.totaleffect=lm(paste0("X24_hour.SBP.CoV~Bacteroides_dorei + ",paste(m4_covar_sbp, collapse = " + ")),merged.dat[female,])
summary(fit.totaleffect)

#Effect of IV on Mediator
fit.mediator=lm(paste0("Blood.SCFA_Isobutyric.Acid ~ Bacteroides_dorei + ",paste(m4_covar, collapse = " + ")), merged.dat[female,])
summary(fit.mediator)

#Effect of mediator on DV controlling for IV
fit.dv=lm(paste0("X24_hour.SBP.CoV ~ Bacteroides_dorei + Blood.SCFA_Isobutyric.Acid + ",paste(m4_covar_sbp, collapse = " + ")),merged.dat[female,])
summary(fit.dv)

#Mediation Analysis
results <- mediate(fit.mediator, fit.dv, treat="Bacteroides_dorei", mediator="Blood.SCFA_Isobutyric.Acid", boot=T)
x <- summary(results)
print(c(x$d.avg.p,x$z.avg.p,x$tau.p))


#######################
### High Prevotella ###
#######################
library(pheatmap)

dat.genus <- genus.data[,order(colSums(genus.data),decreasing = T)[1:15]]
annot = data.frame(Gender=meta.data$Sex, Hypertensive=meta.data$Hypertension.Staging)
rownames(annot) = rownames(meta.data)
pheatmap(t(dat.genus), annotation_col=annot, 
         annotation_names_row=FALSE,
         annotation_names_col=FALSE,
         show_colnames = F,
         fontsize_col=5)

bray.dist.matrix <- vegdist(genus.data, method="bray")
mds <- cmdscale(bray.dist.matrix, eig = TRUE, x.ret = TRUE)
mds.var.per <- round(mds$eig/sum(mds$eig)*100, 1)
mds.values <- mds$points
mds.data <- data.frame(Patient = rownames(mds.values), X = mds.values[,1], Y = mds.values[,2], 
                       Bacteroides = dat.genus$Bacteroides, 
                       Bifidobacterium = dat.genus$Bifidobacterium,
                       Prevotella = dat.genus$Prevotella,
                       Ruminococcus = dat.genus$Ruminococcus)
ggplot(mds.data, aes(x=X, y=Y, col= Bacteroides)) + geom_point(size = 5) + theme_classic() +
  xlab(paste0("MDS1 - ",mds.var.per[1],"%")) + ylab(paste0("MDS2 - ",mds.var.per[2],"%")) +
  scale_color_gradient(low="blue", high="red") + ggtitle("Bacteroides Abundance") +
  theme(axis.title = element_text(size = 16, face="bold"), axis.text = element_text(size = 10, face="bold"),
        legend.title = element_blank(), legend.text =  element_text(size = 16, face="bold"), plot.title = element_text(size = 16, face = "bold"))
ggplot(mds.data, aes(x=X, y=Y, col= Bifidobacterium)) + geom_point(size = 5) + theme_classic() +
  xlab(paste0("MDS1 - ",mds.var.per[1],"%")) + ylab(paste0("MDS2 - ",mds.var.per[2],"%")) +
  scale_color_gradient(low="blue", high="red") + ggtitle("Bifidobacterium Abundance") +
  theme(axis.title = element_text(size = 16, face="bold"), axis.text = element_text(size = 10, face="bold"),
        legend.title = element_blank(), legend.text =  element_text(size = 16, face="bold"), plot.title = element_text(size = 16, face = "bold"))
ggplot(mds.data, aes(x=X, y=Y, col= Prevotella)) + geom_point(size = 5) + theme_classic() +
  xlab(paste0("MDS1 - ",mds.var.per[1],"%")) + ylab(paste0("MDS2 - ",mds.var.per[2],"%")) + 
  scale_color_gradient(low="blue", high="red") + ggtitle("Prevotella Abundance") +
  theme(axis.title = element_text(size = 16, face="bold"), axis.text = element_text(size = 10, face="bold"),
        legend.title = element_blank(), legend.text =  element_text(size = 16, face="bold"), plot.title = element_text(size = 16, face = "bold"))

set.seed(1)
kmean_res <- kmeans(transformed.genus.data,3)$cluster
merged.dat$enterotype <- kmean_res[rownames(merged.dat)]
merged.dat$enterotype <- as.factor(merged.dat$enterotype)
ggplot(mds.data, aes(x=X, y=Y, col= as.factor(kmean_res))) + geom_point(size = 5) + theme_classic() +
  xlab(paste0("MDS1 - ",mds.var.per[1],"%")) + ylab(paste0("MDS2 - ",mds.var.per[2],"%")) +
  ggtitle("K Means Clustering") +
  theme(axis.title = element_text(size = 16, face="bold"), axis.text = element_text(size = 10, face="bold"),
        legend.title = element_blank(), legend.text =  element_text(size = 16, face="bold"), plot.title = element_text(size = 16, face = "bold"))

merged.dat$Prevotella_high <- ifelse(merged.dat$enterotype==3,1,0)
merged.dat$Prevotella_high <- as.factor(merged.dat$Prevotella_high)
merged.dat$Prevotella <- genus.data[rownames(merged.dat),"Prevotella"]

cor.matrix <-  data.frame(gender = character(), model = character(), feature = character(), pval = character(), direction = character(),coefficient = character(), se = character())
for(feature in c(SBP_cov,SBP_sd,SBP_surge,DBP_cov,DBP_sd,acids_variables,macro_micro_nutrients)){
  m1 <- glm_model(merged.dat,'Prevotella_high',feature)
  m1.f <- glm_model(merged.dat[female,],'Prevotella_high',feature)
  m1.m <- glm_model(merged.dat[male,],'Prevotella_high',feature)
  m2 <- glm_model(merged.dat,'Prevotella_high',c(feature,m2_covar,'Sex'))
  m2.f <- glm_model(merged.dat[female,],'Prevotella_high',c(feature,m2_covar))
  m2.m <- glm_model(merged.dat[male,],'Prevotella_high',c(feature,m2_covar))
  m3 <- glm_model(merged.dat,'Prevotella_high',c(feature,m3_covar_sbp,'Sex'))
  m3.f <- glm_model(merged.dat[female,],'Prevotella_high',c(feature,m3_covar_sbp,'Menopause_code'))
  m3.m <- glm_model(merged.dat[male,],'Prevotella_high',c(feature,m3_covar_sbp))
  m4 <- glm_model(merged.dat,'Prevotella_high',c(feature,m4_covar_sbp,'Sex'))
  m4.f <- glm_model(merged.dat[female,],'Prevotella_high',c(feature,m4_covar_sbp,'Menopause_code'))
  m4.m <- glm_model(merged.dat[male,],'Prevotella_high',c(feature,m4_covar_sbp))
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m1",feature,summary(m1)$coefficients[2,4],sign(summary(m1)$coefficients[2,1]),summary(m1)$coefficients[2,1],summary(m1)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m1",feature,summary(m1.f)$coefficients[2,4],sign(summary(m1.f)$coefficients[2,1]),summary(m1.f)$coefficients[2,1],summary(m1.f)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m1",feature,summary(m1.m)$coefficients[2,4],sign(summary(m1.m)$coefficients[2,1]),summary(m1.m)$coefficients[2,1],summary(m1.m)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m2",feature,summary(m2)$coefficients[2,4],sign(summary(m2)$coefficients[2,1]),summary(m2)$coefficients[2,1],summary(m2)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m2",feature,summary(m2.f)$coefficients[2,4],sign(summary(m2.f)$coefficients[2,1]),summary(m2.f)$coefficients[2,1],summary(m2.f)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m2",feature,summary(m2.m)$coefficients[2,4],sign(summary(m2.m)$coefficients[2,1]),summary(m2.m)$coefficients[2,1],summary(m2.m)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m3",feature,summary(m3)$coefficients[2,4],sign(summary(m3)$coefficients[2,1]),summary(m3)$coefficients[2,1],summary(m3)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m3",feature,summary(m3.f)$coefficients[2,4],sign(summary(m3.f)$coefficients[2,1]),summary(m3.f)$coefficients[2,1],summary(m3.f)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m3",feature,summary(m3.m)$coefficients[2,4],sign(summary(m3.m)$coefficients[2,1]),summary(m3.m)$coefficients[2,1],summary(m3.m)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("All","m4",feature,summary(m4)$coefficients[2,4],sign(summary(m4)$coefficients[2,1]),summary(m4)$coefficients[2,1],summary(m4)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Women","m4",feature,summary(m4.f)$coefficients[2,4],sign(summary(m4.f)$coefficients[2,1]),summary(m4.f)$coefficients[2,1],summary(m4.f)$coefficients[2,2])
  cor.matrix[nrow(cor.matrix)+1,] <- c("Men","m4",feature,summary(m4.m)$coefficients[2,4],sign(summary(m4.m)$coefficients[2,1]),summary(m4.m)$coefficients[2,1],summary(m4.m)$coefficients[2,2])
}
cor.matrix <- format_res(cor.matrix)
print(ggplot(cor.matrix, aes(x=model, y = feature, fill=log)) + geom_tile(color="white", size=0.1) +scale_fill_gradient2(midpoint=0, low="blue", mid="white",high="red", space ="Lab" ) +
        coord_equal() + xlab("") + ylab("") + geom_text(aes(label=sym), size = 5, col = "black")+ facet_wrap(~gender) + 
        theme(axis.text = element_text(size = 10, face="bold"), strip.text = element_text(size=12, face="bold"), 
              legend.text = element_text(size = 12, face="bold"), legend.title = element_text(size = 12, face="bold")) + ggtitle("Prevotella High"))


  