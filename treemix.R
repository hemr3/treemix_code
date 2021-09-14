library(SiZer)
library(OptM)
library(phytools)
library(treeio)
library(ggtree)
library(dplyr)
library(ggplot2)
library(ggimage)

setwd("~/treemix-1.13")

source("src/plotting_funcs.R")
  #necessary for plotting trees

plot_tree("~/analyses/snp_to_vcfs/treemix_outputs/all_chr_many_pops")
  #plots tree incl migr and root along drift param

plot_resid("~/analyses/snp_to_vcfs/treemix_outputs/one/noMissing_1", 
           "~/analyses/snp_to_vcfs/pop_order")
  #plot residuals

#plotting phylo based on above, but using phytools####
string = "((Neanderthals, Denisovans), (Papuans, (Africans, Europeans)), Chimpanzee);"
  #this string is based on a treemix graph produced. it can be changed.
  #it's called a newick test string

vert.tree = read.tree(text = string)
  #reads the newick test string

plot(vert.tree, no.margin = T, edge.width = 2)
  #the above is a very simple tree. It looks nice, can't b used for any major
  #analyses, bc it says what is obvious. Can be used probs as intro vis if space

plotTree(vert.tree, type = "fan", fsize = 0.7, lwd = 1, ftype = "i")
  #this is a cute tree but doesn't really add much and is likely 
  #harder to read. 

#ggtree plots ####
tree = read.tree("~/Downloads/tree_newick.nwk")
  #This file in this specific location works. this package only works with dplyr 1.0.5
  #can copy and paste different trees in and it works. 
tree
  #check whether the tree produced is what is wanted. 

gram = ggtree(tree)+
  theme_tree2()+
#  geom_cladelabel(node = 10, label = "H. sapiens",
#                  color = "blue",
#                  offset = .005)+
#  geom_cladelabel(node = 9, label = "Archaics",
#                  color = "purple",
#                  offset = .005)+
    #both cladelabels offset clade labels beside phylo
#  geom_hilight(node = 10, fill = "gold")+
#  geom_hilight(node = 9, fill = "purple")+ 
    #above fills in square of clades
  geom_tiplab()+ 
    #labels tips of nodes with what is in treemix
    #relabel populations in treemix pop_map.tsv to get the proper names
    #it can be done here, more trouble than worth. 
  geom_tippoint(aes(color = "red"))+
#  geom_text(aes(label=node), hjust=-.3)+
    #labels the nodes themselves
#  geom_taxalink("NEA", "EUR", color = "orange", linetype = "dashed",
#                offset=0.003,
#                arrow=arrow(length=unit(0.01, "npc")))+
#  geom_taxalink("AFR", "EUR", color = "gold", linetype = "dashed",
#                offset=0.003,
#                arrow=arrow(length=unit(0.01, "npc")))+
#  geom_taxalink("DEN", "PAP", color = "yellow", linetype = "dashed",
#                offset=0.003,
#                arrow=arrow(length=unit(0.01, "npc")))+
    #links taxa - this can be used to more accurately show migrations
  theme(legend.position = "none") 
    #removes legend

phylopic_info = data.frame(node = c(7, 10, 9),
                           phylopic = c("7133ab33-cc79-4d7c-9656-48717359abb4",
                                        "9c6af553-390c-4bdd-baeb-6992cbc540b1",
                                        "e547cd01-7dd1-495b-8239-52cf9971a609"),
                           species = c("chimpanzee", "human", "archaics"))  
  #creates the df that can be used in adding images to phylo

gram %<+% phylopic_info +
  geom_nodelab(aes(image=phylopic), geom = "phylopic", alpha=.5,
               color = "steelblue")
    #adds a phylopic. a bit squished, but it works


#optimising no. of migration nodes in the trees ####
#believe that clado + mig tree should be two figs side by side

for (edge in 1:5) {
  plot_tree(cex = 0.8, paste("~/analyses/snp_to_vcfs/treemix_outputs/noMissing"
                             , ".", edge))
  title(paste(edge, "edges"))
  
}


setwd("~/analyses/snp_to_vcfs/treemix_outputs")
  #needs to be in the same general filesystem as the treemix outputs

one = system.file("/home/guy/analyses/snp_to_vcfs/treemix_outputs/one", package = "OptM")

test.optm = optM("/home/guy/analyses/snp_to_vcfs/treemix_outputs/migration_selection/", 
                 tsv = "tree_output")
  
plot_optM(test.optm)
