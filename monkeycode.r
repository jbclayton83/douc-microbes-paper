#### DATA WRANGLING ####
#setwd('/data/analyses/monkey/') # Sets the working directory
map = read.delim('EMPV4_mapping_Pygathrix_nemaeus.txt', row.names = 1) # Grab the map

# Relabel the bad labels in the map. 
map$PA = as.character(map$Population) # need this to manipulate words!
map$PA[map$PA=="Philadelphia Zoo"] = "Captive"
map$PA[map$PA=="Singapore Zoo"] = "Semi-captive"
map$PA[map$PA=="EPRC"] = "Semi-wild"
map$PA[!map$PA %in% c("Captive","Semi-captive","Semi-wild")] = "Wild"
map$PA = factor(map$PA,levels=c("Captive","Semi-captive","Semi-wild","Wild"),ordered = T)
map$Lifestyle = map$PA # Call it "Lifestyle"
lscolors = c("#F8766D","#7CAE00","#00BFC4","#C77CFF")

#### PICRUST ####
library(polycor)
library(robCompositions)
library(beeswarm)

# Read in the closed ref table (feed to PICRUSt's normalize)
otu.cr = read.delim('norm_copynum.txt',row=1,skip=1,as.is=T) # get this from picrust stage 1
otu.cr = otu.cr[,c(rownames(map))]  # sync the sample names for the OTU table
sum(otu.cr)
otu.cr = otu.cr[rowMeans(otu.cr) >= 0.01,]
sum(otu.cr)
otu.raw97 = read.delim('combined_puregg97.otu',row=1)
otu.raw97 = otu.raw97[,c(rownames(map))]
pdf("ReadDepth_closed.pdf",width=6,height=5.5)
plot(colSums(otu.raw97) ~ map$Lifestyle,xlab="Lifestyle",ylab="Read depth (closed)",
     col=lscolors,lwd=2)
dev.off()

# Re-normalize PICRUSt stage 1 output with the centered log-ratio, feed to predict/summarize
#table(sort(as.matrix(otu.cr))[0:1000000])
#otu.cr.n = impRZilr(t(otu.cr),method='lm',dl=rep(0.1,nrow(otu.cr)),verbose = T)$x
      otu.cr.n = t(otu.cr)
      otu.cr.n[otu.cr.n < 0.08] = 0.08 # call 0's the "detection limit"
otu.cr.n = t(cenLR(otu.cr.n)$x.clr)
colnames(otu.cr.n) = colnames(otu.cr)

sink("otu_reprocessed.txt"); cat("#OTU ID\t");
write.table(otu.cr.n,file="otu_reprocessed.txt",quote=F,sep="\t",append = T);
sink(NULL)

# Read in the PICRUSt L3 summarized pathways (stage 3 output)
picrust = read.delim('PredictedL3.txt', #'EMPV4_Doucs_R1only_NO_CHLOROPLASTS_ko_metagenome_at_level3.txt', 
                     skip=1, row.names = 1) #Grab picrust table, skipping bad first row
picrust = as.matrix(picrust[,rownames(map)]) # sync and drop extra

nsti = read.delim('nsti_orig.picrustp',row=1)
nsti = nsti[rownames(map),]
pdf("NSTI.pdf",width=6,height=5.5)
plot(nsti$Value ~ map$Lifestyle,ylab="NSTI index",xlab="Lifestyle",col=lscolors,lwd=2)
dev.off()

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
Grp.Pvals = p.adjust(Grp.Pvals)
Wld.Pvals = p.adjust(Wld.Pvals)

# make and sort a data frame with these columns
df = data.frame(Grp.Pvals, Grp.Corrs, Wld.Pvals, row.names = rownames(picrust))
select = abs(df$Grp.Corrs) > 0.3 & df$Grp.Pvals < 0.05 & df$Wld.Pvals < 0.05
df = df[select,]
df = df[order(df$Grp.Corrs),]

# Display all significant with p < 0.05
num_sig = sum(select)                      # Count how many are significant
C_ix = map$CaptiveWild=="Captive"          # Stores "true" if monkey is captive, else "false"
W_ix = map$CaptiveWild=="Wild"             # As above. Use these to select just wild/captive monkeys
pdf(file = "PiSwarms.pdf",width = 6.5, height = 6.5)
sink(file = "PICRUSt_Significance.txt")    # Store output to a text file
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

# PICRUSt heatmap too, why not
library(gplots)
my_palette <- colorRampPalette(c("blue", "black", "yellow"))(n = 299) # Sebastian Raschka
gl = map$Lifestyle
glpos = c(grep("Captive",gl),grep("Semi-captive",gl),grep("Semi-wild",gl),grep("Wild",gl))
gl = gl[glpos]
mat = picrust[rownames(df[abs(df$Grp.Corrs) > 0.75,]),glpos]
mat = sweep(mat,1,rowSums(abs(mat)),'/')                      # Normalize to relative abundance
mat = sweep(mat,1,max(apply(mat,1,max),apply(mat,1,min)),'/') # Constrain extrema to [-1, 1]

levels(gl)= lscolors #c("red","orange","yellow","green")

png("PiMap.png",  # create PNG for the heat map        
    width = 8*300,                        # 5 x 300 pixels
    height = 6*300,
    res = 300,                              # 300 pixels per inch
    pointsize = 11)                          # smaller font size
heatmap.2(mat,
          #cellnote = mat,  # same data set for cell labels
          main = "", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins = c(2,22),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          #breaks=col_breaks,    # enable color transition at specified limits
          ColSideColors = as.character(gl),
          dendrogram="row",     # only draw a row dendrogram
          lhei=c(1,4), lwid=c(1,4),
          labCol = "",
          hclustfun = function(x) hclust(as.dist(1 - cor(as.matrix(x))), method="complete"),
          Colv="NA"            # turn off column clustering
)
par(lend = 1)           # square line ends for the color legend
legend("topright",      # location of the legend on the heatmap plot
       inset=c(.1,-0), # adjust placement upward
       legend = levels(map$Lifestyle), # category labels
       col = levels(gl),  # color key
       lty= 1,            # line style
       lwd = 10,          # line width
       cex = 0.75,
       xpd=TRUE  # allow drawing outside
)
dev.off()
#### Figure 6 Stuff ####
# Manually generate fig 6 of the paper. 
pdf("DietDiv.pdf",width=6,height=5.5)
par(oma=c(0,2,0,0))                                  # Margins!
myData = c(1,15,44,57)                               # Data!
myplot = barplot(myData, names=levels(map$Lifestyle), ylim = c(0,62), xlab = "Lifestyle", ylab = "",
                 col=lscolors)
mtext("Dietary Biodiversity\n(Number of plant species)",side=2,line=3) # Special Y axis!
text(myplot,myData,myData, pos=3)                    # Add numbers on top of the bars!
dev.off()

#### F/B Ratio (from here on, open ref is used) ####
otu = read.delim('EMP_seqs_R1_Doucs_otu_table_mc2_w_tax_open_ref_NO_CHLOROPLASTS.txt',row=1,skip=1,as.is=T)
otu = otu[,c(rownames(map),"taxonomy")]  # sync the sample names for the OTU table
pdf("ReadDepth_open.pdf",width=6,height=5.5)
plot(colSums(otu[,-ncol(otu)]) ~ map$Lifestyle,xlab="Lifestyle",ylab="Read depth (open)",
     col=lscolors,lwd=2)
dev.off()
isFirmicutes = grepl('p__Firmicutes',otu$taxonomy)     # Save "trues" for Firmicutes, false otherwise
isBacteroides = grepl('p__Bacteroidetes',otu$taxonomy) # Like above for Bacteroidetes
FBratio = log(colSums(otu[isFirmicutes,-ncol(otu)])/colSums(otu[isBacteroides,-ncol(otu)]))

kruskal.test(FBratio ~ map$PA)
for (i in 1:length(levels(map$PA))) for (j in i:length(levels(map$PA))) {
  if (i == j) next
  p = wilcox.test(FBratio[map$PA==levels(map$PA)[i]],FBratio[map$PA==levels(map$PA)[j]])
  cat(levels(map$PA)[i]," vs ",levels(map$PA)[j]," p = ",p$p.value,"\n",sep='')
}
pdf("FBratio.pdf",width=6,height=5.5)
plot(FBratio ~ map$PA, xlab="Lifestyle", ylab="Log F:B ratio",col=lscolors,lwd=2)
dev.off()
df = data.frame(FBratio, map$PA)     # Split manually into groups
tapply(df$FBratio, df$map.PA, mean)  # Get the means per group
tapply(df$FBratio, df$map.PA, sd)    # Gets the standard devs per group

#### Alpha diversity violin plots ####
library(ggsignif)
library(vegan)
(mindepth = min(colSums(otu[,-ncol(otu)]))) # for row, -grep("Chloroplast",otu$taxonomy) ?
otu.r = rrarefy(t(otu[,-ncol(otu)]),mindepth) # Rarefy to the minimum sample depth (52,918)
dix = "shannon"
div = diversity(otu.r,index=dix)
otu.ad = data.frame(Div=div, Lifestyle=map$Lifestyle)
grps = levels(map$Lifestyle)
lab = "Alpha Diversity (Shannon)" #paste0("Alpha Diversity (",dix,")")
pdf(paste0("AlphaDiv_",dix,".pdf"),width=6,height=5.5)
plot(ggplot(otu.ad,aes(x=Lifestyle,y=Div,fill=Lifestyle)) + ylab(lab) + geom_violin(alpha=0.3) + 
       geom_signif(comparisons = list(grps[c(1,2)],grps[c(3,4)]), test='t.test', map_signif_level = T) + 
       geom_signif(comparisons = list(grps[c(2,4)]), test='t.test', map_signif_level = T, y_position = 6.1) +
       geom_signif(comparisons = list(grps[c(1,3)]), test='t.test', map_signif_level = T, y_position = 6.3) +
       geom_signif(comparisons = list(grps[c(1,4)]), test='t.test', map_signif_level = T, y_position = 6.5) +
       geom_jitter(aes(color=Lifestyle),position=position_jitter(0.2),size=2) )
dev.off()
otu.cd = data.frame(Div=rowSums(otu.r > 0), Lifestyle=map$Lifestyle)
lab = "Alpha Diversity (Observed OTUs)"
pdf("AlphaDiv_ObsOTU.pdf",width=6,height=5.5)
plot(ggplot(otu.cd,aes(x=Lifestyle,y=Div,fill=Lifestyle)) + ylab(lab) + geom_violin(alpha=0.3) + 
       geom_signif(comparisons = list(grps[c(1,2)],grps[c(3,4)]), test='t.test', map_signif_level = T) + 
       geom_signif(comparisons = list(grps[c(2,4)]), test='t.test', map_signif_level = T, y_position = 5500) +
       geom_signif(comparisons = list(grps[c(1,3)]), test='t.test', map_signif_level = T, y_position = 5850) +
       geom_signif(comparisons = list(grps[c(1,4)]), test='t.test', map_signif_level = T, y_position = 6200) +
       geom_jitter(aes(color=Lifestyle),position=position_jitter(0.2),size=2) )
dev.off()

#### PCoA plots #### (requires map and otu table loaded) 
source('pcoa_helper.R') # This gives us our nice pcoa functions
library(phyloseq)
tree = read_tree_greengenes('EMP_seqs_R1_Doucs_rep_set_open_ref.tre')
otu.s = as.matrix(otu[,-ncol(otu)])  # Rip of the taxonomy column, as it is not needed 
bray = vegdist(t(otu.s))             # Get some bray curtis distances
pdf("bray.pdf",width=6,height=4.75); plot_pcoa(bray,map,category='Lifestyle');
pdf("uuf.pdf",width=6,height=4.75); pcoa.u = plot_unifrac(otu.s,map,tree,category='Lifestyle',weight=F); 
pdf("wuf.pdf",width=6,height=4.75); pcoa.w = plot_unifrac(otu.s,map,tree,category='Lifestyle',weight=T); 
graphics.off()
adonis(pcoa.u ~ map$Lifestyle)       # Do stats for clustering (unweighted)
adonis(pcoa.w ~ map$Lifestyle)       # Do stats for clustering (weighted)
adonis(bray ~ map$Lifestyle)         # Do stats for bray-curtis too, why not

# do a boxplot of the main axis
pdf("PC1_boxplot.pdf",width=8,height=3); par(oma=c(0,2,0,0))
boxplot(pcoa(pcoa.u)$vectors[,1] ~ map$Lifestyle, horizontal=T, las=2, ylim=c(-0.18,0.42),
        col=lscolors)
dev.off()

#### Differential taxa testing ####
library("gplots")
library("RColorBrewer")
library(robCompositions) # Composition magic
library(polycor)
library(beeswarm)
library(reshape2)

bT = c(2,4,6,7)  # Levels to consider for the bugs
pT = c(2,3,3,4)  # Levels to consider for the plants
for (L in 1:length(bT)) {
  # Massage the taxa names
  split = strsplit(otu$taxonomy[-grep("c__Chloroplast",otu$taxonomy)],"; ")        # Split by semicolon into levels
  taxaStrings = sapply(split,function(x) paste(x[1:bT[L]],collapse=";"))           # Create collapsed names
  for (i in 1:7) taxaStrings = gsub("(;[A-z]__$)?(;NA$)?","",taxaStrings,perl=T)   # Clean tips
  otu.t = rowsum(otu[-grep("c__Chloroplast",otu$taxonomy),-ncol(otu)],taxaStrings) # Collapse table by name
  
  ## Chloroplasts (can be skipped)
  ch = read.delim("chloro95y.tax",row=1)
  split = strsplit(rownames(ch),";")
  taxaStrings = sapply(split,function(x) paste(x[1:pT[L]],collapse=";"))          # Split by semicolon into levels
  for (i in 1:7) taxaStrings = gsub("(;[A-z]__$)?(;NA$)?","",taxaStrings,perl=T)  # Clean tips
  ch.t = rowsum(ch,taxaStrings,drop=F)                    # Collapse table by name
  sz = merge(t(as.matrix(otu.t)),t(ch.t),by=0,all=T)      # Merge chloroplasts and taxa
  rownames(sz) = sz$Row.names                             # Repopulate row names with faux column
  sz = t(as.matrix(sz[,2:ncol(sz)]))                      # remove faux column and flip it back
  sz[is.na(sz)] = 0                                       # pad empty columns with 0 counts
  otu.t = sz                                              # Replace original taxa table with merged one
  #limits = t(rbind(rowSums(otu.t > 0),rowMeans(otu.t)))   # Visualization purposes only: 
  #limits = limits[order(limits[,1]),]                     # What do my samples look like?
  
  # Filter out bugs that don't appear in enough samples
  select = rowSums(otu.t > 0) > min(table(map$Lifestyle))/2 # reasonable general min
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
  otu.m = otu.n[-grep("Viridiplantae",rownames(otu.n)),,drop=F]
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
  pdf(paste0("TaxaSummary_L",bT[L],".pdf"),width = 8,height=7) # Make room for legend
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
  Grp.Pvals = p.adjust(gpb)
  Wld.Pvals = p.adjust(wpb)
  res = data.frame(Grp.Pvals, Grp.Corrs, Wld.Pvals,row.names=rownames(otu.t))
  res = res[order(res$Grp.Corrs),]
  
  # Add bivariate filter
  sig = 0.05
  selection = res$Grp.Pvals < sig & res$Wld.Pvals < sig & abs(res$Grp.Corrs) > 0.3
  
  # Display all significant with p < 0.05
  num_sig = sum(selection, na.rm = T) # Count how many are significant
  res = res[selection,]
  # Truncate names to last 2 informative levels
  split = strsplit(rownames(res),";")        # Split by semicolon into levels
  res$short = sapply(split,function(x) paste(tail(x,2),collapse=";"))
  
  C_ix = map$CaptiveWild=="Captive"          # Stores "true" if monkey is captive, else "false"
  W_ix = map$CaptiveWild=="Wild"             # As above. Use these to select just wild/captive
  pdf(paste0("TaxaSwarms_L",bT[L],".pdf"),width = 6.5,height=6.5)
  sink(paste0("Taxa_Significance_L",bT[L],".txt"))                  # Get ready to write the significant ones
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
  
  ## Heatmap
  # Need to have created the clr taxa[+plant] table as a matrix
  my_palette <- colorRampPalette(c("blue", "black", "yellow"))(n = 299) # Sebastian Raschka
  gl = map$Lifestyle
  glpos = c(grep("Captive",gl),grep("Semi-captive",gl),grep("Semi-wild",gl),grep("Wild",gl))
  gl = gl[glpos]
  mat = otu.t[rownames(res),glpos]
  mat = mat[-grep("Viridiplantae",rownames(mat)),,drop=F]
  #mat = sweep(mat,1,max(apply(mat,1,max),apply(mat,1,min)),'/')
  
  # Truncate names to last 2 informative levels
  split = strsplit(as.character(rownames(mat)),";")        # Split by semicolon into levels
  rownames(mat) = sapply(split,function(x) paste(tail(x,2),collapse=";")) 
  
  levels(gl)= lscolors 
  png(paste0("Taxa_heatmap_L",bT[L],".png"),  # create PNG for the heat map        
      width = 8*300,                        # 8 x 300 pixels
      height = 6*300,
      res = 300,                              # 300 pixels per inch
      pointsize = 10)                          # smaller font size
  heatmap.2(mat,
            #cellnote = mat,  # same data set for cell labels
            main = "", # heat map title
            notecol="black",      # change font color of cell labels to black
            density.info="none",  # turns off density plot inside color legend
            trace="none",         # turns off trace lines inside the heat map
            margins = c(2,22),     # widens margins around plot
            col=my_palette,       # use on color palette defined earlier
            #breaks=col_breaks,    # enable color transition at specified limits
            ColSideColors = as.character(gl),
            dendrogram="row",     # only draw a row dendrogram
            lhei=c(1,4.7), lwid=c(1,5),
            labCol = "",
            hclustfun = function(x) hclust(as.dist(1 - cor(as.matrix(x))), method="complete"),
            Colv="NA"            # turn off column clustering
  )
  par(lend = 1)           # square line ends for the color legend
  legend("topright",      # location of the legend on the heatmap plot
         inset=c(.1,-0), # adjust placement upward
         legend = levels(map$Lifestyle), # category labels
         col = levels(gl),  # color key
         lty= 1,            # line style
         lwd = 10,          # line width
         cex = 0.65,
         xpd=TRUE  # allow drawing outside
  )
  dev.off()
  
  # repeat for plants alone
  mat = otu.t[rownames(res),glpos]
  mat = mat[grep("Viridiplantae",rownames(mat)),,drop=F]
  if (nrow(mat) < 2) next
  # Truncate names to last 2 informative levels
  split = strsplit(as.character(rownames(mat)),";")        # Split by semicolon into levels
  rownames(mat) = sapply(split,function(x) paste(tail(x,2),collapse=";")) 
  
  png(paste0("Plant_heatmap_L",pT[L],".png"),  # create PNG for the heat map        
      width = 8*300,                        # 5 x 300 pixels
      height = 6*300,
      res = 300,                              # 300 pixels per inch
      pointsize = 11)                          # smaller font size
  heatmap.2(mat,
            #cellnote = mat,  # same data set for cell labels
            main = "", # heat map title
            notecol="black",      # change font color of cell labels to black
            density.info="none",  # turns off density plot inside color legend
            trace="none",         # turns off trace lines inside the heat map
            margins = c(2,22),     # widens margins around plot
            col=my_palette,       # use on color palette defined earlier
            #breaks=col_breaks,    # enable color transition at specified limits
            ColSideColors = as.character(gl),
            dendrogram="row",     # only draw a row dendrogram
            lhei=c(1,4), lwid=c(1,4),
            labCol = "",
            hclustfun = function(x) hclust(as.dist(1 - cor(as.matrix(x))), method="complete"),
            Colv="NA"            # turn off column clustering
  )
  par(lend = 1)           # square line ends for the color legend
  legend("topright",      # location of the legend on the heatmap plot
         inset=c(.1,-0), # adjust placement upward
         legend = levels(map$Lifestyle), # category labels
         col = levels(gl),  # color key
         lty= 1,            # line style
         lwd = 10,          # line width
         cex = 0.75,
         xpd=TRUE  # allow drawing outside
  )
  dev.off()
}

# Final combo heatmap
# with more taxa
Grp.Pvals = p.adjust(gpb)
Wld.Pvals = p.adjust(wpb)
res = data.frame(Grp.Pvals, Grp.Corrs, Wld.Pvals,row.names=rownames(otu.t))
res = res[order(res$Grp.Corrs),]
selection = res$Grp.Pvals < 0.05 & res$Wld.Pvals < 0.01 & abs(res$Grp.Corrs) > 0.45
sum(selection)
res = res[selection,]

mat = otu.t[rownames(res),glpos]
# Truncate names to last 2 informative levels
split = strsplit(as.character(rownames(mat)),";")        # Split by semicolon into levels
rownames(mat) = sapply(split,function(x) paste(tail(x,2),collapse=";")) 

levels(gl)= lscolors 
png(paste0("BiClust_dense.png"),  # create PNG for the heat map        
    width = 8*300,                        # 8 x 300 pixels
    height = 6*300,
    res = 300,                              # 300 pixels per inch
    pointsize = 9)                          # smaller font size
heatmap.2(mat,
          #cellnote = mat,  # same data set for cell labels
          main = "", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins = c(2,22),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          #breaks=col_breaks,    # enable color transition at specified limits
          ColSideColors = as.character(gl),
          #dendrogram="row",     # only draw a row dendrogram
          lhei=c(1,4.7), lwid=c(1,5),
          labCol = "",
          hclustfun = function(x) hclust(as.dist(1 - cor(as.matrix(x))), method="complete")
          #Colv="NA"            # turn off column clustering
)
par(lend = 1)           # square line ends for the color legend
legend("topright",      # location of the legend on the heatmap plot
       inset=c(.1,-0), # adjust placement upward
       legend = levels(map$Lifestyle), # category labels
       col = levels(gl),  # color key
       lty= 1,            # line style
       lwd = 10,          # line width
       cex = 0.65,
       xpd=TRUE  # allow drawing outside
)
dev.off()

Grp.Pvals = p.adjust(gpb)
Wld.Pvals = p.adjust(wpb)
res = data.frame(Grp.Pvals, Grp.Corrs, Wld.Pvals,row.names=rownames(otu.t))
res = res[order(res$Grp.Corrs),]
selection = res$Grp.Pvals < 0.05 & res$Wld.Pvals < 0.001 & abs(res$Grp.Corrs) > 0.8
sum(selection)
res = res[selection,]

mat = otu.t[rownames(res),glpos]
# Truncate names to last 2 informative levels
split = strsplit(as.character(rownames(mat)),";")        # Split by semicolon into levels
rownames(mat) = sapply(split,function(x) paste(tail(x,2),collapse=";")) 

levels(gl)= lscolors 
png(paste0("BiClust_sparse.png"),  # create PNG for the heat map        
    width = 8*300,                        # 8 x 300 pixels
    height = 6*300,
    res = 300,                              # 300 pixels per inch
    pointsize = 9)                          # smaller font size
heatmap.2(mat,
          #cellnote = mat,  # same data set for cell labels
          main = "", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins = c(2,22),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          #breaks=col_breaks,    # enable color transition at specified limits
          ColSideColors = as.character(gl),
          #dendrogram="row",     # only draw a row dendrogram
          lhei=c(1,4.7), lwid=c(1,5),
          labCol = "",
          hclustfun = function(x) hclust(as.dist(1 - cor(as.matrix(x))), method="complete")
          #Colv="NA"            # turn off column clustering
)
par(lend = 1)           # square line ends for the color legend
legend("topright",      # location of the legend on the heatmap plot
       inset=c(.1,-0), # adjust placement upward
       legend = levels(map$Lifestyle), # category labels
       col = levels(gl),  # color key
       lty= 1,            # line style
       lwd = 10,          # line width
       cex = 0.65,
       xpd=TRUE  # allow drawing outside
)
dev.off()
