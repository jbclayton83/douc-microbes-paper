#### DATA WRANGLING ####
setwd('/data/analyses/monkey/') # Sets the working directory
map = read.delim('EMPV4_mapping_Pygathrix_nemaeus.txt', row.names = 1) # Grab the map
picrust = read.delim('EMPV4_Doucs_R1only_NO_CHLOROPLASTS_ko_metagenome_at_level3.txt', 
                     skip=1, row.names = 1) #Grab picrust table, skipping bad first row

# Sync the picrust table to the map, drop extra columns, normalize
picrust = as.matrix(picrust[,rownames(map)]) # sync and drop extra
picrust = sweep(picrust, 2, colSums(picrust), FUN="/") # normalize

# Relabel the bad labels in the map. 
map$PA = as.character(map$Population) # need this to manipulate words!
map$PA[map$PA=="Philadelphia Zoo"] = "USA zoo"
map$PA[map$PA=="Singapore Zoo"] = "SE Asia zoo"
map$PA[!map$PA %in% c("USA zoo","SE Asia zoo","EPRC")] = "Wild"
map$PA = factor(map$PA,levels=c("USA zoo","SE Asia zoo","EPRC","Wild"),ordered = T)
map$Population = map$PA # overwrite the old Population column with the new

#### PICRUST ####
# Go through each picrust pathway and test for significance w/group
CW = map$CaptiveWild %in% c("Captive","Wild") # For selecting only captive and wild
Grp.Pvals = rep(1,nrow(picrust)) # Initialize p values for groupwise sig. tests
Wld.Pvals = rep(1,nrow(picrust)) # Initialize p values for wild vs captive tests
for (m.ix in 1:nrow(picrust)) {  # Loop through all the rows (picrust pathways)
  Grp.Pvals[m.ix] = kruskal.test(picrust[m.ix,] ~ map$PA)$p.value
  Wld.Pvals[m.ix] = wilcox.test(picrust[m.ix,CW] ~ map$CaptiveWild[CW])$p.value
}

# Sort by significance, adjust for multiple tests
Sort.ix = order(Wld.Pvals)                 # Stores the order of them (like third, first, second...)
PicrustIDs = rownames(picrust)[Sort.ix]    # Keep a list of the IDs in order
Grp.Qvals = p.adjust(Grp.Pvals[Sort.ix],method = "holm") # Adj. group P-values in order
Wld.Qvals = p.adjust(Wld.Pvals[Sort.ix],method = "holm") # Adj. wild-captive " " "

# Display all significant with p < 0.05
num_sig = sum(Wld.Qvals < 0.05, na.rm = T) # Count how many are significant
C_ix = map$CaptiveWild=="Captive"          # Stores "true" if monkey is captive, else "false"
W_ix = map$CaptiveWild=="Wild"             # As above. Use these to select just wild/captive monkeys
sink(file = "Dumped_Significance.txt")     # Store output to a text file
cat("Pathway\tGroupwise_Q\tCaptiveVsWild_Q\tTrendInCaptivity\n")  # The header line in the file
if (num_sig) for (i in 1:num_sig) {        # Loop through all significant pathways
  upInCaptive = mean(picrust[PicrustIDs[i],C_ix]) > mean(picrust[PicrustIDs[i],W_ix]) # compare avgs
  cat(PicrustIDs[i],'\t',Grp.Qvals[i],'\t',Wld.Qvals[i],'\t',ifelse(upInCaptive,"UP","DOWN"),'\n')
}
sink(NULL)                                 # Close the output file

#### Figure 6 Stuff ####
# Manually generate fig 6 of the paper. 
par(oma=c(0,2,0,0))                                  # Margins!
myData = c(1,15,43,57)                               # Data!
myLabels = c("USA zoo","SE Asia zoo","EPRC","Wild")  # Labels!
myplot = barplot(myData, names=myLabels, ylim = c(0,62), xlab = "Population", ylab = "", col=0)
mtext("Dietary Biodiversity\n(Number of plant species)",side=2,line=3) # Special Y axis!
text(myplot,myData,myData, pos=3)                    # Add numbers on top of the bars!

#### F/B Ratio ####
otu = read.delim('EMP_seqs_R1_Doucs_otu_table_mc2_w_tax_open_ref_NO_CHLOROPLASTS.txt',row=1,skip=1)
otu = otu[,c(rownames(map),"taxonomy")]  # sync the sample names for the OTU table
isFirmicutes = grepl('p__Firmicutes',otu$taxonomy) # Save "trues" where Firmicutes, false otherwise
isBacteroides = grepl('p__Bacteroidetes',otu$taxonomy) # Like above
FBratio = log(colSums(otu[isFirmicutes,-ncol(otu)])/colSums(otu[isBacteroides,-ncol(otu)]))

kruskal.test(FBratio ~ map$PA)
for (i in 1:length(levels(map$PA))) for (j in i:length(levels(map$PA))) {
  if (i == j) next
  p = wilcox.test(FBratio[map$PA==levels(map$PA)[i]],FBratio[map$PA==levels(map$PA)[j]])
  cat(levels(map$PA)[i]," vs ",levels(map$PA)[j]," p = ",p$p.value,"\n")
}
plot(FBratio ~ map$PA, xlab="Population", ylab="Log F:B ratio")

#### PCoA plots #### (requires map and otu table loaded) [otu table in prev section preferred?]
source('pcoa_helper.R') # This gives us our nice pcoa functions
tree = read_tree_greengenes('EMP_seqs_R1_Doucs_rep_set_open_ref.tre')
otu.s = as.matrix(otu[,-ncol(otu)])  # Rip of the taxonomy column, as it is not needed 
bray = vegdist(t(otu.s))             # Get some bray curtis distances
plot_pcoa(bray,map,category = 'Population')  # Plot the Bray Curtis distances
pcoa.u = plot_unifrac(otu.s,map,tree,category = 'Population',weighted=F)
pcoa.w = plot_unifrac(otu.s,map,tree,category = 'Population',weighted=T)
adonis(pcoa.u ~ map$Population)      # Do stats for clustering (unweighted)
adonis(pcoa.w ~ map$Population)      # Do stats for clustering (weighted)
adonis(bray ~ map$Population)        # Do stats for bray-curtis too, why not
