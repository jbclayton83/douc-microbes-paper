#### DATA WRANGLING ####
#setwd('/data/analyses/monkey/') # Sets the working directory
map = read.delim('data/EMPV4_mapping_Pygathrix_nemaeus.txt', row.names = 1) # Grab the map

# Relabel the bad labels in the map. 
map$PA = as.character(map$Population) # need this to manipulate words!
map$PA[map$PA=="Philadelphia Zoo"] = "Captive"
map$PA[map$PA=="Singapore Zoo"] = "Semi-captive"
map$PA[map$PA=="EPRC"] = "Semi-wild"
map$PA[!map$PA %in% c("Captive","Semi-captive","Semi-wild")] = "Wild"
map$PA = factor(map$PA,levels=c("Captive","Semi-captive","Semi-wild","Wild"),ordered = T)
map$Lifestyle = map$PA # Call it "Lifestyle"
lscolors = c("#F8766D","#7CAE00","#00BFC4","#C77CFF")
map.orig = map

# aggregate map by monkey id
map.ag = aggregate(map, by=list(map$AnimalID), FUN=function(xx) min(as.character(xx)))
rownames(map.ag) = map.ag$Group.1
map = map.ag[,2:ncol(map.ag)]

#### PICRUST ####
library(polycor)
library(robCompositions)
library(beeswarm)

## Read in the PICRUSt L3 summarized pathways 
picrust = read.delim('data/PredictedL3.txt', #'EMPV4_Doucs_R1only_NO_CHLOROPLASTS_ko_metagenome_at_level3.txt', 
                     skip=1, row.names = 1) #Grab picrust table, skipping bad first row
picrust = as.matrix(picrust[,rownames(map.orig)]) # sync and drop extra
ag = aggregate(t(picrust), by=list(map.orig$AnimalID), FUN=mean)
rownames(ag)=ag$Group.1
picrust = t(ag[,2:ncol(ag)])


# Go through each picrust pathway and test for significance w/group
CW = map$CaptiveWild %in% c("Captive","Wild") # For selecting only captive and wild
Grp.Pvals = rep(1,nrow(picrust)) # Initialize p values for groupwise sig. tests
Grp.Corrs = rep(0,nrow(picrust)) # Initialize correlation chamber (groupwise)
Wld.Pvals = rep(1,nrow(picrust)) # Initialize p values for wild vs captive tests
for (m.ix in 1:nrow(picrust)) {  # Loop through all the rows (picrust pathways)
  try({
    ps = polyserial(picrust[m.ix,],map$PA,ML=T,std.err = T)
    if (is.na(ps$rho)) next
    Grp.Corrs[m.ix] = ps$rho
    Grp.Pvals[m.ix] = 1-pchisq(ps$chisq, ps$df) 
  },silent=T)
  Wld.Pvals[m.ix] = wilcox.test(picrust[m.ix,CW] ~ map$CaptiveWild[CW])$p.value
}

# Adjust for multiple tests
#Grp.Pvals = p.adjust(Grp.Pvals)
Wld.Pvals = p.adjust(Wld.Pvals,method="fdr")

# make and sort a data frame with these columns
df = data.frame(Grp.Pvals, Grp.Corrs, Wld.Pvals, row.names = rownames(picrust))
select = abs(df$Grp.Corrs) > 0 # & df$Grp.Pvals <= 0.2 # & df$Wld.Pvals < 0.25
df = df[select,]
df = df[order(df$Grp.Corrs),]

# Display all significant with p < 0.05
num_sig = sum(select)                      # Count how many are significant
C_ix = map$CaptiveWild=="Captive"          # Stores "true" if monkey is captive, else "false"
W_ix = map$CaptiveWild=="Wild"             # As above. Use these to select just wild/captive monkeys
pdf(file = "results/PiSwarms.pdf",width = 6.5, height = 6.5)
sink(file = "results/PICRUSt_Significance.txt")    # Store output to a text file
cat("Pathway\tPolyserial_Q\tPolyserial_Cor\tCaptiveVsWild_Q\tTrendInCaptivity\n")  # The header line in the file
if (num_sig) for (i in 1:num_sig) {
  cur = rownames(df)[i]
  upInCaptive = mean(picrust[cur,C_ix]) > mean(picrust[cur,W_ix]) # compare avgs
  cat(cur,'\t',df$Grp.Pvals[i],'\t',-df$Grp.Corrs[i],'\t',df$Wld.Pvals[i],'\t',
      ifelse(upInCaptive,"UP","DOWN"),'\n',sep='')
  beeswarm(picrust[cur,] ~ map$PA, xlab="Lifestyle",ylab="CLR Relative Abundance",main=cur,
           col=alpha(lscolors,0.7),cex.axis=1.1,cex.main=1,cex=1.1,corral="random",pch=19)
  bxplot(picrust[cur,] ~ map$PA, add = TRUE)
}
sink(NULL)                                 # Close the output file
dev.off()


#### F/B Ratio (from here on, open ref is used) ####
map$PA = factor(map$PA)
otu = read.delim('data/EMP_seqs_R1_Doucs_otu_table_mc2_w_tax_open_ref_NO_CHLOROPLASTS.txt',row=1,skip=1,as.is=T)
otu = otu[,c(rownames(map.orig),"taxonomy")]  # sync the sample names for the OTU table

otu.tx = rowsum(otu[,-ncol(otu)],otu$taxonomy)
ag = sweep(otu.tx,2,colSums(otu.tx),'/');
ag = aggregate(t(ag), by=list(map.orig$AnimalID), FUN=mean)
rownames(ag)=ag$Group.1
ag = t(ag[,2:ncol(ag)])

otu = ag;

isFirmicutes = grepl('p__Firmicutes',rownames(otu))     # Save "trues" for Firmicutes, false otherwise
isBacteroides = grepl('p__Bacteroidetes',rownames(otu)) # Like above for Bacteroidetes
FBratio = log(colSums(otu[isFirmicutes,])/colSums(otu[isBacteroides,]))

kruskal.test(FBratio ~ map$PA)
for (i in 1:length(levels(map$PA))) for (j in i:length(levels(map$PA))) {
  if (i == j) next
  p = wilcox.test(FBratio[map$PA==levels(map$PA)[i]],FBratio[map$PA==levels(map$PA)[j]])
  cat(levels(map$PA)[i]," vs ",levels(map$PA)[j]," p = ",p$p.value,"\n",sep='')
}
pdf("results/FBratio.pdf",width=6,height=5.5)
plot(FBratio ~ factor(map$PA), xlab="Lifestyle", ylab="Log F:B ratio",col=lscolors,lwd=2)
dev.off()
df = data.frame(FBratio, map$PA)     # Split manually into groups
tapply(df$FBratio, df$map.PA, mean)  # Get the means per group
tapply(df$FBratio, df$map.PA, sd)    # Gets the standard devs per group

#### Alpha diversity violin plots ####
library(ggsignif)
library(vegan)
otu = read.delim('data/EMP_seqs_R1_Doucs_otu_table_mc2_w_tax_open_ref_NO_CHLOROPLASTS.txt',row=1,skip=1,as.is=T)
taxaStrings = otu$taxonomy
otu = otu[,rownames(map.orig)]
otu.orig = otu

otu.tx = t(otu)
ag = sweep(otu.tx,1,rowSums(otu.tx),'/');
ag = aggregate(ag, by=list(map.orig$AnimalID), FUN=mean)
rownames(ag)=ag$Group.1
ag = data.frame(t(ag[,2:ncol(ag)]),check.names = F)
#colnames(ag) = rownames(map)
ag$taxonomy = taxaStrings
otu = ag;

# now the table starts off as relAb, multiply up to min depth
otu.r = round(t(otu[,-ncol(otu)] * 52918))
#otu.r = rrarefy(otu.r,52918)
  
dix = "shannon"
div = diversity(otu.r,index=dix)
otu.ad = data.frame(Div=div, Lifestyle=map$Lifestyle)
map$Lifestyle = factor(map$Lifestyle)
grps = levels(map$Lifestyle)
lab = "Alpha Diversity (Shannon)" #paste0("Alpha Diversity (",dix,")")
pdf(paste0("results/AlphaDiv_",dix,".pdf"),width=6,height=5.5)
plot(ggplot(otu.ad,aes(x=Lifestyle,y=Div,fill=Lifestyle)) + ylab(lab) + geom_violin(alpha=0.3) + 
       geom_signif(comparisons = list(grps[c(1,2)],grps[c(3,4)]), test='t.test', map_signif_level = T) + 
       geom_signif(comparisons = list(grps[c(2,4)]), test='t.test', map_signif_level = T, y_position = 6.1) +
       geom_signif(comparisons = list(grps[c(1,3)]), test='t.test', map_signif_level = T, y_position = 6.3) +
       geom_signif(comparisons = list(grps[c(1,4)]), test='t.test', map_signif_level = T, y_position = 6.5) +
       geom_jitter(aes(color=Lifestyle),position=position_jitter(0.2),size=2) )
dev.off()
otu.cd = data.frame(Div=rowSums(otu.r > 0), Lifestyle=map$Lifestyle)
lab = "Alpha Diversity (Observed OTUs)"
pdf("results/AlphaDiv_ObsOTU.pdf",width=6,height=5.5)
plot(ggplot(otu.cd,aes(x=Lifestyle,y=Div,fill=Lifestyle)) + ylab(lab) + geom_violin(alpha=0.3) + 
       geom_signif(comparisons = list(grps[c(1,2)],grps[c(3,4)]), test='t.test', map_signif_level = T) + 
       geom_signif(comparisons = list(grps[c(2,4)]), test='t.test', map_signif_level = T, y_position = 5800) +
       geom_signif(comparisons = list(grps[c(1,3)]), test='t.test', map_signif_level = T, y_position = 6150) +
       geom_signif(comparisons = list(grps[c(1,4)]), test='t.test', map_signif_level = T, y_position = 6500) +
       geom_jitter(aes(color=Lifestyle),position=position_jitter(0.2),size=2) )
dev.off()

#### PCoA plots #### (requires map and otu table loaded) 
source('lib/pcoa_helper.R') # This gives us our nice pcoa functions
library(phyloseq)
tree = read_tree_greengenes('data/EMP_seqs_R1_Doucs_rep_set_open_ref.tre')
otu.s = as.matrix(otu[,-ncol(otu)])*52918  # Rip off the taxonomy column, as it is not needed

bray = vegdist(t(otu.s))             # Get some bray curtis distances
pcoa.u = plot_unifrac(otu.s,map,tree,category='Lifestyle',weight=F); 
pcoa.w = plot_unifrac(otu.s,map,tree,category='Lifestyle',weight=T); 
adonis(pcoa.u ~ map$Lifestyle)       # Do stats for clustering (unweighted)
adonis(pcoa.w ~ map$Lifestyle)       # Do stats for clustering (weighted)
adonis(bray ~ map$Lifestyle)         # Do stats for bray-curtis too, why not

#### Differential taxa testing ####
library("gplots")
library("RColorBrewer")
library(robCompositions) # Composition magic
library(polycor)
library(beeswarm)
library(reshape2)

bT = c(2,4,6,7)  # Levels to consider for the bugs
pT = c(2,3,3,4)  # Levels to consider for the plants
otu = data.frame(otu.s,check.names = F)
otu$taxonomy = taxaStrings

#prefilter for existing chloroplasts
if (length(grep("c__Chloroplast",otu$taxonomy)))
  otu = otu[-grep("c__Chloroplast",otu$taxonomy),]

for (L in 1:length(bT)) {
  # Massage the taxa names
  split = strsplit(otu$taxonomy[-grep("c__Chloroplast",otu$taxonomy)],"; ")        # Split by semicolon into levels
  taxaStrings = sapply(split,function(x) paste(x[1:bT[L]],collapse=";"))           # Create collapsed names
  for (i in 1:7) taxaStrings = gsub("(;[A-z]__$)?(;NA$)?","",taxaStrings,perl=T)   # Clean tips
  otu.t = rowsum(otu[,-ncol(otu)],taxaStrings) # Collapse table by name
  
  # Filter out bugs that don't appear in enough samples
  select = rowSums(otu.t > 1) > min(table(map$Lifestyle))/2 # reasonable general min
  otu.t = otu.t[select,]                                    # Apply drop mask
  otu.t = otu.t[,rownames(map)]                             # Sync with map samp order
  otu.n = sweep(otu.t,2,colSums(otu.t),'/');                # Normalize to relative abundance
  
  otu.c = impRZilr(t(otu.t),dl=rep(1,nrow(otu.t)),maxit = 3,verbose = T,method = "lm") # zeros
  otu.c = t(cenLR(otu.c$x)$x.clr)    # Centered log-ratio transform for compositions
  colnames(otu.c) = colnames(otu.t)  # Because this gets rid of the names...
  otu.t = otu.c                      # otu.t is our active table; give it the CLR
  
  # Go through each taxon and test for significance w/group
  otu.t = as.matrix(otu.t)
  ntax = nrow(otu.t)
  CW = map$CaptiveWild %in% c("Captive","Wild") # For selecting only captive and wild
  Grp.Pvals=rep(1,ntax)
  Grp.Corrs=rep(0,ntax)
  Wld.Pvals=rep(1,ntax)
  KW.Pvals=rep(1,ntax)
  for (m.ix in 1:ntax) {  # Loop through all the rows (taxa)
    try({ # Because some correlations may be inadmissable
      ps = polyserial(otu.t[m.ix,],map$PA,ML=T,std.err = T)
      if (is.na(pchisq(ps$chisq, ps$df))) next # Invalid correlation
      Grp.Corrs[m.ix] = ps$rho             # Find intensity of correlation
      Grp.Pvals[m.ix] = 1-pchisq(ps$chisq, ps$df) # And p-value on this
    },silent=T)
    Wld.Pvals[m.ix] = wilcox.test(otu.t[m.ix,CW] ~ map$CaptiveWild[CW])$p.value
    KW.Pvals[m.ix] = kruskal.test(otu.t[m.ix,] ~ map$Lifestyle)$p.val
  }
  
  ## Taxa barplots -- Top 15 most abundant (kruskal sig. + other?)
  otu.m = otu.n
  otu.m = sweep(sqrt(otu.m),2,colSums(sqrt(otu.m)),'/')
  meanAb = apply(otu.m,1,FUN=function(x) tapply(x, map$Lifestyle, mean)) # group mean
  ranked = order(apply(meanAb,2,max),decreasing=T)
  otu.m = otu.m[ranked,]
  
  # Truncate names to last 2 informative levels
  split = strsplit(rownames(otu.m),";")        # Split by semicolon into levels
  Taxa = sapply(split,function(x) paste(tail(x,2),collapse=";")) 
  lim = 25
  if (nrow(otu.m) > lim) Taxa[lim:nrow(otu.m)] = "Other"
  otu.m = rowsum(otu.m,Taxa)
  
  # Sort by average abundance for display
  byAbundance = rownames(otu.m)[order(rowMeans(otu.m),decreasing=T)]
  
  # Acrobatics for ggplot2 (yes, this is inane)
  otu.m = data.frame(t(otu.m),check.names=F)      # flip table
  otu.m$SampleID = rownames(otu.m)  # add a column for the sample IDs
  map$SampleID = rownames(map)      # add a column for the sample IDs
  
  # The following separates taxa abundances per sample, then splices in the column of interest
  otu.m = melt(otu.m, id.vars = "SampleID", variable.name = "Taxa", value.name = "RelativeAbundance")
  otu.m = merge(otu.m, map[,c("SampleID","Lifestyle")], by="SampleID")
  otu.m$Taxa = factor(otu.m$Taxa,levels=byAbundance,ordered=T) # add Taxa column
  
  ## Plot according to Lifestyle, sorted by abundance
  pdf(paste0("results/TaxaSummary_L",bT[L],".pdf"),width = 8,height=7) # Make room for legend
  plot(ggplot(otu.m, aes(x = Lifestyle, y = RelativeAbundance, fill = Taxa)) +
         geom_bar(stat ="identity", position="fill") + labs(x="Lifestyle",y="Root Relative Abundance") +
         guides(fill=guide_legend(ncol=1)) +
         scale_fill_manual(values=c("dodgerblue2","#E31A1C", # red # Kevin Wright
                                    "green4",
                                    "#6A3D9A", # purple
                                    "#FF7F00", # orange
                                    "black","gold1",
                                    "skyblue2","#FB9A99", # lt pink
                                    "palegreen2",
                                    "#CAB2D6", # lt purple
                                    "#FDBF6F", # lt orange
                                    "gray70", "khaki2",
                                    "maroon","orchid1","deeppink1","blue1","steelblue4",
                                    "darkturquoise","green1","yellow4","yellow3",
                                    "darkorange4","brown"))) 
  dev.off()
  
  ## Display differential taxa/stats
  # Adjust for multiple tests, sort by significance
  gpb = Grp.Pvals; wpb = Wld.Pvals;
  #Grp.Pvals = p.adjust(gpb)
  Wld.Pvals = p.adjust(wpb,method = "fdr")
  res = data.frame(Grp.Pvals, Grp.Corrs, Wld.Pvals,row.names=rownames(otu.t))
  res = res[order(res$Grp.Corrs),]
  
  # Add bivariate filter (actually now just selects all valid entries for collapsed version)
  selection = abs(res$Grp.Corrs) > 0 # & res$Wld.Pvals < .25
  
  # Display all
  num_sig = sum(selection, na.rm = T) # Count how many are significant
  res = res[selection,]
  # Truncate names to last 2 informative levels
  split = strsplit(rownames(res),";")        # Split by semicolon into levels
  res$short = sapply(split,function(x) paste(tail(x,2),collapse=";"))
  
  C_ix = map$CaptiveWild=="Captive"          # Stores "true" if monkey is captive, else "false"
  W_ix = map$CaptiveWild=="Wild"             # As above. Use these to select just wild/captive
  pdf(paste0("results/TaxaSwarms_L",bT[L],".pdf"),width = 6.5,height=6.5)
  sink(paste0("results/Taxa_Significance_L",bT[L],".txt"))                  # Get ready to write the significant ones
  cat("Taxon\tPolyserial_Q\tPolyserial_Cor\tCaptiveVsWild_Q\tTrendInCaptivity\n")  # Print header
  if (num_sig) for (i in 1:num_sig) {
    taxon = rownames(res)[i]
    upInCaptive = mean(otu.t[taxon,C_ix]) > mean(otu.t[taxon,W_ix]) # compare avgs
    cat(res[taxon,]$short,'\t',res$Grp.Pvals[i],'\t',-res$Grp.Corrs[i],'\t',res$Wld.Pvals[i],'\t',
        ifelse(upInCaptive,"UP","DOWN"),'\n',sep='')
    beeswarm(otu.t[taxon,] ~ map$PA, xlab="Lifestyle",ylab="CLR Relative Abundance",main=res[taxon,]$short,
             col=alpha(lscolors,0.7),cex.axis=1.1,cex.main=1,cex=1.1,corral="random",pch=19)
    bxplot(otu.t[taxon,] ~ map$PA, add = TRUE)
  }
  sink(NULL)
  dev.off()
 
}

## Chloroplasts (can be skipped)
ch = read.delim("data/chloro95y.tax",row=1)
toPad = setdiff(rownames(map.orig),colnames(ch))
if (length(toPad)) for (c in toPad) 
  ch[,toPad] = rep(0,nrow(ch))
m.p = map.orig[colnames(ch),]
m.p.ag = aggregate(m.p, by=list(m.p$AnimalID), FUN=function(xx) min(as.character(xx)))
rownames(m.p.ag) = m.p.ag$Group.1
m.p.ag = m.p.ag[,2:ncol(m.p.ag)]

# Collapse by monkey
ag = aggregate(t(ch), by=list(m.p$AnimalID), FUN=mean)
rownames(ag)=ag$Group.1
ch = data.frame(t(ag[,2:ncol(ag)]),check.names = F)

# tack on the otu table for context
ch$taxonomy = rownames(ch)
ch = rbind(ch,otu)

# Collapse to lv 3
split = strsplit(ch$taxonomy,";")
taxaStrings = sapply(split,function(x) paste(x[1:3],collapse=";"))          # Split by semicolon into levels
for (i in 1:7) taxaStrings = gsub("(;[A-z]__$)?(;NA$)?","",taxaStrings,perl=T)  # Clean tips
ch.t = rowsum(ch[,-ncol(ch)],taxaStrings,drop=F)                    # Collapse table by name

# CLR
# eps=.5
# ch.ag = t(ch.t);
# ch.ag[ch.ag<eps]=0
# ch.ag = ch.ag*(1 - rowSums(ch.ag==0)*eps/rowSums(ch.ag))
# ch.ag[ch.ag==0]=eps
# ch.ag = sweep(ch.ag,1,rowSums(ch.ag),'/');
# ls = log(ch.ag)
# ch.ag = t(ls - rowMeans(ls))
# ch.ag = ch.ag[,!is.nan(colSums(ch.ag))]

ch.t = ch.t[rowSums(ch.t > .5) > 0,]
otu.c = impRZilr(t(ch.t),dl=rep(.5,nrow(ch.t)),maxit = 3,verbose = T,method = "lm") # zeros
otu.c = t(cenLR(otu.c$x)$x.clr)    # Centered log-ratio transform for compositions
colnames(otu.c) = colnames(ch.t)  # Because this gets rid of the names...
ch.ag = otu.c

m.p.ag = m.p.ag[colnames(ch.ag),]

# stats
m.p.ag$Lifestyle = factor(m.p.ag$Lifestyle)
otu.t = as.matrix(ch.ag)
otu.t = otu.t[grep("k__Viridiplantae",rownames(otu.t)),]
ntax = nrow(otu.t)
CW = m.p.ag$CaptiveWild %in% c("Captive","Wild") # For selecting only captive and wild
Grp.Pvals=rep(1,ntax)
Grp.Corrs=rep(0,ntax)
Wld.Pvals=rep(1,ntax)
KW.Pvals=rep(1,ntax)
for (m.ix in 1:ntax) {  # Loop through all the rows (taxa)
  try({ # Because some correlations may be inadmissable
    ps = polyserial(otu.t[m.ix,],m.p.ag$PA,ML=T,std.err = T)
    if (is.na(pchisq(ps$chisq, ps$df))) next # Invalid correlation
    Grp.Corrs[m.ix] = ps$rho             # Find intensity of correlation
    Grp.Pvals[m.ix] = 1-pchisq(ps$chisq, ps$df) # And p-value on this
  },silent=T)
  Wld.Pvals[m.ix] = wilcox.test(otu.t[m.ix,CW] ~ m.p.ag$CaptiveWild[CW])$p.value
  KW.Pvals[m.ix] = kruskal.test(otu.t[m.ix,] ~ m.p.ag$Lifestyle)$p.val
}
Wld.Pvals = p.adjust(Wld.Pvals,method = "fdr")
res = data.frame(Grp.Pvals, Grp.Corrs, Wld.Pvals,row.names=rownames(otu.t))
res = res[order(res$Grp.Corrs),]

selection = abs(res$Grp.Corrs) > 0 # & res$Wld.Pvals < .25

# Display all
num_sig = sum(selection, na.rm = T) # Count how many are significant
res = res[selection,]
# Truncate names to last 2 informative levels
split = strsplit(rownames(res),";")        # Split by semicolon into levels
res$short = sapply(split,function(x) paste(tail(x,2),collapse=";"))

C_ix = m.p.ag$CaptiveWild=="Captive"          # Stores "true" if monkey is captive, else "false"
W_ix = m.p.ag$CaptiveWild=="Wild"             # As above. Use these to select just wild/captive
pdf("results/DiffPlants_L3.pdf",width = 6.5,height=6.5)
sink("results/Plant_Significance_L3.txt")                  # Get ready to write the significant ones
cat("Taxon\tPolyserial_Q\tPolyserial_Cor\tCaptiveVsWild_Q\tTrendInCaptivity\n")  # Print header
if (num_sig) for (i in 1:num_sig) {
  taxon = rownames(res)[i]
  upInCaptive = mean(otu.t[taxon,C_ix]) > mean(otu.t[taxon,W_ix]) # compare avgs
  cat(res[taxon,]$short,'\t',res$Grp.Pvals[i],'\t',-res$Grp.Corrs[i],'\t',res$Wld.Pvals[i],'\t',
      ifelse(upInCaptive,"UP","DOWN"),'\n',sep='')
  beeswarm(otu.t[taxon,] ~ m.p.ag$PA, xlab="Lifestyle",ylab="CLR Relative Abundance",main=res[taxon,]$short,
           col=alpha(lscolors,0.7),cex.axis=1.1,cex.main=1,cex=1.1,corral="random",pch=19)
  bxplot(otu.t[taxon,] ~ m.p.ag$PA, add = TRUE)
}
sink(NULL)
dev.off()


