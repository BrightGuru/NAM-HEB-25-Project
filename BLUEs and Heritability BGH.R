###packages
library(readr)
library("lme4")
library(corrplot)
library(PerformanceAnalytics)
library(ggcorrplot)
library(dplyr)
library(data.table)
library(MASS)
library(bnstruct)
library(readr)
library(qqman)
library(magrittr)
library(CMplot)
library(tidyverse)

###Loading the data ####
Pathogen_AllHEB_2017 = read_delim("Pathogen.AllHEB.2017.csv", ";", escape_double = FALSE, trim_ws = TRUE)
Pathogen_AllHEB_2018 <- read_delim("Pathogen.AllHEB.2018.csv",";", escape_double = FALSE, trim_ws = TRUE)


Pathogen_AllHEB_2017_M=melt(Pathogen_AllHEB_2017,id.vars="Line", measure.vars=c("Bgh0519_A","Bgh0608_B")) #,"Bgh0617_C"
Pathogen_AllHEB_2018_M=melt(Pathogen_AllHEB_2018,id.vars="Genotype", measure.vars=c("Bgh0518_A" ,"Bgh0607_B"))

Pathogen_AllHEB_2017_M$variable=substring(Pathogen_AllHEB_2017_M$variable, 9)
Pathogen_AllHEB_2018_M$variable=substring(Pathogen_AllHEB_2018_M$variable, 9)

Pathogen_AllHEB_2017_M$Year=1
Pathogen_AllHEB_2018_M$Year=2

colnames(Pathogen_AllHEB_2017_M)= c("Genotype", "Stage", "BGH","Year")
colnames(Pathogen_AllHEB_2018_M)= c("Genotype", "Stage", "BGH","Year")

BLUEs_Her=do.call(rbind.data.frame, list(Pathogen_AllHEB_2017_M,Pathogen_AllHEB_2018_M))

BLUEs_Her$Year=as.factor(BLUEs_Her$Year)
BLUEs_Her$Stage=as.factor(BLUEs_Her$Stage)
BLUEs_Her$BGH=as.numeric(BLUEs_Her$BGH)
BLUEs_Her$Genotype=as.factor(BLUEs_Her$Genotype)

#################### Heritability ########
str(BLUEs_Her$BGH)

Blues_model=lmer(BGH ~ (1|Genotype) + (1|Year) + (1|Genotype:Year), data = BLUEs_Her)

summary(Blues_model)

Her=2.35/(2.35 +(0/2)+(2.33/4)) 
###### BLUEs#####
Blues_fun = function(traits,data="."){
  B=as.data.frame(fixef(lmer(paste0(traits,"~0 + Genotype + (1|Year) + (1|Genotype:Year)"),data=data)))
}
traits="BGH"

Blues_BGH=lapply(traits, Blues_fun, data=BLUEs_Her)
Blues_BGH=map2(Blues_BGH, traits, ~set_names(..1,..2)%>%
                 rownames_to_column(var='Genotype'))%>%reduce(full_join)

############ descriptive statistics ####

summary(Pathogen_AllHEB_2017)
summary(as.numeric(Pathogen_AllHEB_2018$Bgh0518_A))
summary(as.numeric(Pathogen_AllHEB_2018$Bgh0607_B))

##### Spearman correction
cor_data=merge(Pathogen_AllHEB_2017,Pathogen_AllHEB_2018, by.x="Line", by.y="Genotype")
colnames(cor_data)
cor_data=subset(cor_data, select=c( "Line","Bgh0519_A", "Bgh0608_B","Bgh0518_A","Bgh0607_B"))
setnames(cor_data, old =c( "Line","Bgh0519_A", "Bgh0608_B","Bgh0518_A","Bgh0607_B"),new =c( "Genotype","Bgh2017_A", "Bgh2017_B","Bgh2018_A","Bgh2018_B")  )                

cor_data$Bgh2018_A=as.numeric(cor_data$Bgh2018_A)
cor_data$Bgh2018_B=as.numeric(cor_data$Bgh2018_B)
summary(cor_data)

corr_Res=cor(na.omit(cor_data[,-1]),method = c("spearman"))
corr_Res=round(corr_Res,2)

chart.Correlation(na.omit(cor_data[,-1]), histogram=TRUE, pch=19)

ggcorrplot(corr_Res, hc.order = TRUE, type = "upper",
           outline.col = "white", lab = TRUE)
?ggcorrplot
summary(corr_Res)

###histogram with GGPLOt####
ggplot(Blues_BGH, aes(x=BGH)) + 
  geom_histogram(color="black", fill="lightblue",binwidth=0.5)+
  theme_bw()+
  xlab("Scores")+
  ylab("Number of Lines")

###GWAS with base####

####Data preparation for GWAS

#Load in the data

MapData=fread("Map1.csv")
GenoData=fread("FLootime data.csv")

##Map data prep
MapData$Chromosome=plyr::revalue(MapData$Chromosome, c("1H"=1, "2H"=2, "3H"=3,"4H"=4,"5H"=5,"6H"=6,"7H"=7,"n/a"=8))
MapData$Position=plyr::revalue(MapData$Position, c("n/a"=0))
MapData$Chromosome=as.numeric(MapData$Chromosome)
MapData$Position=as.numeric(MapData$Position)

###BGH Phenotype
Blues_BGH$Genotype=substring(Blues_BGH$Genotype,9)

#Extracting the SNPs
Genome_Blues=filter(GenoData, Line %in% Blues_BGH$Genotype)

#Replace missing values in the Genome data
GIMPU=Genome_Blues[,-1:-3]
GIMPU[GIMPU == 5 ] = 0
GIMPU=as.matrix(GIMPU)
GIMPU=knn.impute(GIMPU)

h_Blues=data.frame(Genotype=Genome_Blues$Line,Family=Genome_Blues$`HEB family`)

Genome_Blues_1=data.frame(h_Blues,GIMPU)

#Merge the phenotype and Genotype data

DATA_Blues_BGH=merge(Blues_BGH,Genome_Blues_1, by.x="Genotype",by.y="Genotype")

#Family as cofactor

FamilyCofactor_BGH = model.matrix(~as.character(Family)-1,DATA_Blues_BGH)

#########Powdery mildrew#####

##Stepwise Regression to get the SNPs as cofactor
SNPS=colnames(GIMPU)
BGHvariables=subset(DATA_Blues_BGH,select=SNPS)
BGHvariables1= data.frame(1,BGHvariables)

modelBGH=lm(DATA_Blues_BGH$BGH~1, data = DATA_Blues_BGH)
summary(modelBGH)

modelBGH1=step(modelBGH,scope = formula(BGHvariables), direction="both",k = log(nrow(BGHvariables)))

BGHselectedvariables=model.matrix(modelBGH1)
BGHselectedvariables=as.matrix(BGHselectedvariables[,-1])

#BGHselectedvariablesMap=filter(MapData, SNP %in% c("SCRI_RS_167808" , "SCRI_RS_228030" , "SCRI_RS_168591" , "SCRI_RS_206768" , "BOPA1_13924_403", "SCRI_RS_235332" ,"BOPA1_7782_410" ))


##Model B, multiple regression considering the SNP and family effects

BGHGWAS = apply(X=BGHvariables, MARGIN =  2, FUN= function(x){summary(lm(DATA_Blues_BGH$BGH ~ x+BGHselectedvariables+FamilyCofactor_BGH))})
BGHpval = apply(X=BGHvariables, MARGIN =  2, FUN= function(x){summary(lm(DATA_Blues_BGH$BGH~ x+BGHselectedvariables+FamilyCofactor_BGH))$coeff[2,4]})


#Adjustment for multiple testing
alpha=0.05
#bonferroni correction####
BGH.p.bonf=p.adjust(BGHpval,method="bonferroni")
head(sort(BGH.p.bonf),10)
sum(BGH.p.bonf<alpha)

#FDR correction####
BGH.p.fdr=p.adjust(BGHpval,method="fdr")
head(sort(BGH.p.fdr),10)
sum(BGH.p.fdr<alpha)


#connecting p-values with SNP_IDs####

BGHPValues=data.frame(SNP=SNPS,P.raw=BGHpval, P.Bon=BGH.p.bonf, P.FDR=BGH.p.fdr)
BGHPValuesMAP=merge(BGHPValues,MapData,by.x="SNP", by.y="SNP")

BGHselectpvalue=subset(BGHPValuesMAP,P.FDR<alpha )


#### Manhattan Plot #####

manhattan(BGHPValuesMAP, chr="Chromosome", bp="Position", snp="SNP", p="P.Bon",suggestiveline=F, genomewideline=-log10(0.05))
manhattan(GAPIT_Blink_Bgh0617_GWAS_Results, chr="Chromosome", bp="Position", snp="SNP", p="P.value",suggestiveline=F, genomewideline=-log10(0.05/5710))
####QQPLOT###
qq(BGHpval)



####GWAS with Gapit ####
Namecol_BGH = c("Genotype","BGH")
PHE_BGH=subset(DATA_Blues_BGH, select=Namecol_BGH)

#GEnomedata
SNP_BGH=colnames(Genome_Blues_1[,-2])
GenoData_BGH=subset(DATA_Blues_BGH,select=SNP_BGH)
#Mapdata
MAPFINAL=MapData

#Family as cofactor
COfactor_BGH=data.frame(Genotype=DATA_Blues_BGH$Genotype,FamilyCofactor_BGH)

#GWAS Proper
myGAPIT_BGH <- GAPIT(
  Y=PHE_BGH[,c(1,2)],
  GD=GenoData_BGH,
  GM=MAPFINAL,
  CV=COfactor_BGH,
  model= "BLINK"
)


#### CMPlot to create GGPLot ####
##GApit
GAPIT_Blink_BGH_GWAS_Results = read_csv("GAPIT.Blink.BGH.GWAS.Results.csv")
GAPIT_Blink_BGH_GWAS_Results = arrange(GAPIT_Blink_BGH_GWAS_Results, Chromosome)
GAPIT_Blink_BGH_GWAS_Results = GAPIT_Blink_BGH_GWAS_Results[, 1:4]
CMplot(GAPIT_Blink_BGH_GWAS_Results,type="p",plot.type="m",LOG10=TRUE,threshold= 0.05/5700,
       amplify=TRUE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE,
       width=12,height=6,chr.labels.angle=45)
###Base
BGHPValuesMAP_CMplot = arrange(data.frame(SNP=BGHPValuesMAP$SNP,Chromosome=BGHPValuesMAP$Chromosome, Position=BGHPValuesMAP$Position, P.value=BGHPValuesMAP$P.raw), Chromosome)
CMplot(BGHPValuesMAP_CMplot,type="p",plot.type="m",LOG10=TRUE,threshold= 0.05/5700,
       amplify=TRUE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE,
       width=12,height=6,chr.labels.angle=45)

##### BOTH #####
BCD = data.frame(SNP=BGHPValuesMAP_CMplot$SNP,
                 B=BGHPValuesMAP_CMplot$P.value)
Both_BGH =list(GAPIT_Blink_BGH_GWAS_Results, BCD)
Both_BGH = Reduce(function(x,y) merge(x, y, by="SNP"),Both_BGH)
names(Both_BGH)[4]=paste("A")       
Both_BGH=arrange(Both_BGH, Chromosome)

SNPs <- list(
  Both_BGH$SNP[Both_BGH$A],
  Both_BGH$SNP[Both_BGH$B]
)
CMplot(Both_BGH, plot.type="m",multracks=TRUE,threshold=c(1e-6,1e-4),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=NULL,
       signal.cex=1, file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE,
       highlight=SNPs,highlight.text.cex=1.4)

### QTL selection ####

#extracting the significant SNPs

BGHPValuesMAP_Sign=subset(BGHPValuesMAP,P.Bon<0.05)
GAPIT_Blink_BGH_GWAS_Results_sign = subset(GAPIT_Blink_BGH_GWAS_Results,`FDR_Adjusted_P-values`< 0.05)
