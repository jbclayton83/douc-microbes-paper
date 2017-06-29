## Set up the R environment
#source("http://bioconductor.org/biocLite.R");
#biocLite("phyloseq");
library(ape);
library(vegan);
library(picante);
library(ggplot2);
library(phyloseq);
library(doParallel);
library(foreach);
cl = makeCluster(8);
registerDoParallel(cl);

## define helpers
plot_pcoa = function(dist_obj, sub.map, category, do_ellipse = T, title = NULL) {
  pc = pcoa(dist_obj);
  PCOA = data.frame(pc$vectors);
  for (i in 1:ncol(PCOA)) names(PCOA)[i] = sprintf("PC%d",i); 
  PCOA[,category] = sub.map[,category];
  p = ggplot(PCOA) +
    geom_point(aes(x = PC1, y = PC2, color = PCOA[,category])) +
    labs(title=title) +
    labs(color=category) + labs(fill=category);
  if (do_ellipse) return(p+stat_ellipse(level=0.95, alpha= 0.3, geom="polygon", linetype="blank", 
                                        aes(x=PC1, y=PC2, fill=PCOA[,category])));
  return(p);
}
plot_unifrac = function(sub.otu, sub.map, tree, category, weighted = F, do_ellipse = T, title = NULL) {
  ufo = UniFrac(phyloseq(otu_table(sub.otu,taxa_are_rows = T),tree),parallel = T,weighted = weighted);
  plot(plot_pcoa(ufo,sub.map, category, do_ellipse = do_ellipse, title = title));
  return(ufo);
}