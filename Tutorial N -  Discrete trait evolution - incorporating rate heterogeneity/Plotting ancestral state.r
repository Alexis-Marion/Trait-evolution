library("corHMM")
library("stringr")
library("phytools")
library("secsse")
library("ggtree")
library("RColorBrewer")
library("ggplot2")
library("ggpmisc")

df<-read.csv("example_data.tsv", sep="\t", header =TRUE) # omit sep ="\t" for .csv files
phy<-read.tree("example_tree.tree") #phy<-read.nexus("example_tree.nex")

states<-df[,c(1,2)]
states_traits<-states[!states[,1] %in% setdiff(states[,1], phy$tip.label),]
colnames(states_traits)<-c("Species", "Cat")
states_traits[is.na(states_traits)]<-"?"
states_traits$Cat<-as.character(states_traits$Cat)

best_fitting_model<-readRDS("data_BF_BT.rds")

ASE_plot_maker<-function(phy, trait_model, data_states_model){
phylo<-ggtree(phy, layout="fan", open.angle=180)  + 
           theme_bw() +
           theme(panel.border = element_blank(),
           legend.key = element_blank(),
           axis.ticks = element_blank(),
           axis.text.y = element_blank(),
           axis.text.x = element_blank(),
           panel.grid = element_blank(),
           panel.grid.minor = element_blank(), 
           panel.grid.major = element_blank(),
           panel.background = element_blank(),
           plot.background = element_rect(fill = "transparent",colour = NA))

node_states_traits<-trait_model$states

### assuming you have no polytomies

node_states<-c((length(phy$tip.label)+1):(2*length(phy$tip.label)-1))

### assuming you have a binary trait    

    if(trait_model$rate.cat == 1){
        col1 <- (node_states_traits[,1])
        col2 <- (node_states_traits[,2])
    }
    if(trait_model$rate.cat == 2){
        col1 <- (node_states_traits[,1]+node_states_traits[,3])
        col2 <- (node_states_traits[,2]+node_states_traits[,4])
    }
    if(trait_model$rate.cat == 3){
        col1 <- (node_states_traits[,1]+node_states_traits[,3]+node_states_traits[,5])
        col2 <- (node_states_traits[,2]+node_states_traits[,4]+node_states_traits[,6])
    }
    node_states_traits<-cbind(col1, col2, node_states)
    colnames(node_states_traits)<-c("No", "Yes", "node")
    node_states_traits<-as.data.frame(node_states_traits)
    saveRDS(node_states_traits, "data_categorical_trait.rds")
    pies <- nodepie(node_states_traits, cols=1:2, color=c("#9E0142", "#74BDCB"), alpha=1)
    df<-tibble::tibble(node=as.numeric(node_states_traits$node), pies=pies)
    phylo_node <- phylo %<+% df
    phylo_complete<-phylo_node + geom_plot(data=td_filter(!isTip), mapping=aes(x=x,y=y, label=pies), vp.width=0.03, vp.height=0.03, hjust=0.5, vjust=0.5)
    phylo_complete<-phylo_complete %<+% as.data.frame(states_traits)
    
    ASE_plot<-phylo_complete + geom_tippoint(data=td_filter(isTip), aes(color=Cat), size=1) + scale_color_manual(values=c("#9E0142", "#74BDCB"))

ggsave(ASE_plot, filename = "Ase_plot.pdf",  bg = "transparent", width = 10, height = 10)
}

ASE_plot_maker(phy, best_fitting_model, states_traits)
