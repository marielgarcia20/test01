#This script is from my phd thesis

library(readxl)

#### Set directory
setwd("M:/Metabolomics PA14/mga12/XCMSOnlineNoMIC/neg/")
dir()

## Read feature table

df <- as.data.frame(read.table("mga12_neg_XCMSOnline_2019param_ODnorm_log2.txt", 
                               header=TRUE,stringsAsFactors=FALSE))
names(df)


###########################################
####### Annotation with Metaboscape ########
############################################

### Read Metaboscape annotated table (MS/MS)

df_Metabo <- data.frame(read.table("mga12noMIC(neg)_MetaboScape.csv",header=TRUE, sep = ","),stringsAsFactors=F)
df_Metabo <- data.frame(df_Metabo[,(2:4)],stringsAsFactors = F)
names(df_Metabo)

df_Metabo$Name <- as.character(df_Metabo$Name)
df_Metabo <- subset(df_Metabo, df_Metabo$Name != "")

##

ppm=8E-6
#rt=25 #secs
rt=25/60 #min

df_Metabo$mzmin <- df_Metabo$m.z*(1-ppm)
df_Metabo$mzmax <- df_Metabo$m.z*(1+ppm)
df_Metabo$rtmin <- df_Metabo$RT - rt
df_Metabo$rtmax <- df_Metabo$RT + rt

a <- data.frame(stringsAsFactors = F)
b <- data.frame(stringsAsFactors = F)


for (i in 1:nrow(df)) {
  
  for (j in 1:nrow(df_Metabo)) {
    
    if(df$mzmed[i] <= df_Metabo$mzmax[j] & df$mzmed[i] >= df_Metabo$mzmin[j] & 
       ((df$rtmed[i] <= df_Metabo$rtmax[j] & df$rtmed[i] >= df_Metabo$rtmin[j]) | 
        (df$rtmin[i] <= df_Metabo$rtmax[j] & df$rtmin[i] >= df_Metabo$rtmin[j]) |
        (df$rtmax[i] <= df_Metabo$rtmax[j] & df$rtmax[i] >= df_Metabo$rtmin[j])))
      {
    
      a <- cbind(df[i,],"MetaboScape" = df_Metabo$Name[j], 
                            "ID" = row.names(df)[i]) 
      b <- rbind(b,a)
      }
  }
}

dup <- b[duplicated(b$ID),]

b <- subset(b, !(row.names(b) %in% row.names(dup)))
b <- b[,c(1:(ncol(b)-1))]

names(b)
names(df)

df1 <- cbind(df,"MetaboScape" = "")
df1$MetaboScape <- as.character(df1$MetaboScape)
row.names(df1) <- row.names(df)

df1 <- rbind(subset(df1, !(row.names(df1) %in% row.names(b))),b)
names(df1)

###########################################
####### Annotation with XCMSOnline ########
############################################
## Read xcms online annotation (MS)
dir()

xcms_anno <- read_excel("_tentative_featurematch_mummichog.xlsx")#, stringsAsFactors=FALSE)
xcms_anno <- xcms_anno[!duplicated(xcms_anno$`m/z`),]
xcms_anno <- xcms_anno[!duplicated(xcms_anno$id),]

#xcmsOnline annotations has not retention time
ppm=5E-6
#rt=25 #secs
rt=25/60 #min

a <- data.frame(stringsAsFactors = F)
b <- data.frame(stringsAsFactors = F)

xcms_anno$mzmin <- xcms_anno$`m/z`*(1-ppm)
xcms_anno$mzmax <- xcms_anno$`m/z`*(1+ppm)

for (i in 1:nrow(df1)) {
  
  for (j in 1:nrow(xcms_anno)) {
    
    if(df1$mzmed[i] >= xcms_anno$mzmin[j] & df1$mzmed[i] <= xcms_anno$mzmax[j]){
      a <- cbind(df1[i,],"xcmsOnline" = xcms_anno$name[j], "Pathway" = xcms_anno$pathway[j], 
                 "ID" = row.names(df)[i])
      b <- rbind(b,a)
    }
  }
}



b <- b[!(duplicated(b$ID)),c(1:(ncol(b)-1))]
names(b)
names(df1)

df2 <- cbind(df1, "xcmsOnline"= "", "Pathway" = "")

df2 <- rbind(subset(df2, !(row.names(df2) %in% row.names(b))),b)
df2$xcmsOnline <- as.character(df2$xcmsOnline)
df2$Pathway <- as.character(df2$Pathway)
names(df2)

df3 <- df2[,c(1:8,33:35,9:32)]

## Save annotated feature list
write.csv(df3,"mga12_neg_XCMSOnline_2019param_ODnorm_log2_annognps.csv")



###########################################
####### Annotation with GNPS ########
############################################

## Feature-based molecular networking: FBMN
gnps_anno <- data.frame(read.table("M:/Metabolomics PA14/mga12/Cytoscape/FBMN/mga12_xcmsFBMN_2019param.csv", header=TRUE, sep = ",",stringsAsFactors=F))

## Anaylsis from mzxml data, no feature-based
#gnps_anno <- data.frame(read.table("M:/Metabolomics PA14/mga12/Cytoscape/mga12_WTparC_IC50_MolNet.csv", header=TRUE, sep = ",",stringsAsFactors=F))

names(gnps_anno)
# Retention time in seconds!!!

identical(gnps_anno$parent.mass, gnps_anno$precursor.mass)
identical(gnps_anno$RTConsensus, gnps_anno$RTMean)

ppm=8E-6
# Retention time in seconds!!!
rt=25 #secs
#rt=25/60 #min

a <- data.frame(stringsAsFactors = F)
b <- data.frame(stringsAsFactors = F)

gnps_anno$mzmin <- gnps_anno$precursor.mass*(1-ppm)
gnps_anno$mzmax <- gnps_anno$precursor.mass*(1+ppm)
gnps_anno$rtmin <- gnps_anno$RTMean - rt
gnps_anno$rtmax <- gnps_anno$RTMean + rt

for (i in 1:nrow(df2)) {
  
  for (j in 1:nrow(gnps_anno)) {
    
    if(df2$mzmed[i] <= gnps_anno$mzmax[j] & df2$mzmed[i] >= gnps_anno$mzmin[j] & 
       ((df2$rtmed[i]*60 <= gnps_anno$rtmax[j] & df2$rtmed[i]*60 >= gnps_anno$rtmin[j]) | 
        (df2$rtmin[i]*60 <= gnps_anno$rtmax[j] & df2$rtmin[i]*60 >= gnps_anno$rtmin[j]) |
        (df2$rtmax[i]*60 <= gnps_anno$rtmax[j] & df2$rtmax[i]*60 >= gnps_anno$rtmin[j])))
    {
      
      a <- cbind(df2[i,],"GNPS" = gnps_anno$Compound_Name[j], "Cluster" = gnps_anno$componentindex[j],
                 "ID"= row.names(df2)[i])
      b <- rbind(b,a)
    }
  }
}

dup <- b[duplicated(b$ID),]
b <- subset(b, !(row.names(b) %in% row.names(dup)))
b <- b[,c(1:(ncol(b)-1))]
names(b)
names(df2)
### there is something wronk when assigining the cluster information
df3 <- cbind(df2, "GNPS" = "", "Cluster" = "NA")
row.names(df3) <- row.names(df2)
identical(names(b),names(df3))

df3$GNPS <- as.character(df3$GNPS)
df3$Cluster <- as.character(df3$Cluster)

b$GNPS <- as.character(b$GNPS)
b$Cluster <- as.character(b$Cluster)

df3 <- rbind(subset(df3, !(row.names(df3) %in% row.names(b))),b)
names(df3)

###### Count annotations

x_annoMSMS <- nrow(subset(df3,df3$MetaboScape != "")); x_annoMSMS
x_annoxcmsOL <- nrow(subset(df3,df3$xcmsOnline != "")); x_annoxcmsOL
x_annoGNPS <- nrow(subset(df3,df3$GNPS != "")); x_annoGNPS
x_annoClust <- nrow(subset(df3,df3$Cluster != "NA")); x_annoClust

## annotated features by MSMS with Metaboscape = 36, annotated features MS with xcmsOnline = 262, with GNPS = 30

# rearrange names
names(df3)
df4 <- df3[,c(1:8,33:37,9:32)]

## Save annotated feature list
write.csv(df4,"mga12_XCMSOnline_2019param_ODnorm_log2_annognps.csv")
