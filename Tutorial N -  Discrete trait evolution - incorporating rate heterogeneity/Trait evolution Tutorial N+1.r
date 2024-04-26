library("corHMM")
library("stringr")
library("phytools")
library("qpcR")

df<-read.csv("example_data_MT.tsv", sep="\t", header =TRUE) # omit sep ="\t" for .csv files
phy<-read.tree("example_tree.tree") #phy<-read.nexus("example_tree.nex")

states<-df[,c(1,2)]
states_traits<-states[!states[,1] %in% setdiff(states[,1], phy$tip.label),]
colnames(states_traits)<-c("Species", "Cat")
states_traits[is.na(states_traits)]<-"?"

bt_1_eq<-corHMM(phy, states_traits, rate.cat = 1, rate.mat=NULL, model = "ER", node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

bt_1_sym<-corHMM(phy, states_traits, rate.cat = 1, rate.mat=NULL, model = "SYM", node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

bt_1_ard<-corHMM(phy, states_traits, rate.cat = 1, rate.mat=NULL, model = "ARD", node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

bt_2_eq<-corHMM(phy, states_traits, rate.cat = 2, rate.mat=NULL, model = "ER", node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

bt_2_sym<-corHMM(phy, states_traits, rate.cat = 2, rate.mat=NULL, model = "SYM", node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

bt_2_ard<-corHMM(phy, states_traits, rate.cat = 2, rate.mat=NULL, model = "ARD", node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

bt_3_eq<-corHMM(phy, states_traits, rate.cat = 3, rate.mat=NULL, model = "ER", node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

bt_3_sym<-corHMM(phy, states_traits, rate.cat = 3, rate.mat=NULL, model = "SYM", node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

bt_3_ard<-corHMM(phy, states_traits, rate.cat = 3, rate.mat=NULL, model = "ARD", node.states = "marginal",
fixed.nodes=FALSE, p=NULL, root.p="yang", ip=NULL, nstarts=0, n.cores=5,
get.tip.states = FALSE, lewis.asc.bias = FALSE, collapse = TRUE, lower.bound = 1e-9,
upper.bound = 100, opts=NULL)

df_trait<-data.frame(cbind(c(bt_1_eq$AICc,  bt_1_sym$AICc, bt_1_ard$AICc,
                          bt_2_eq$AICc, bt_2_sym$AICc, bt_2_ard$AICc,
                          bt_3_eq$AICc, bt_3_sym$AICc,  bt_3_ard$AICc),
                akaike.weights(c(bt_1_eq$AICc,  bt_1_sym$AICc, bt_1_ard$AICc,
                          bt_2_eq$AICc, bt_2_sym$AICc, bt_2_ard$AICc,
                          bt_3_eq$AICc, bt_3_sym$AICc,  bt_3_ard$AICc))$deltaAIC,
                akaike.weights(c(bt_1_eq$AICc,  bt_1_sym$AICc, bt_1_ard$AICc,
                          bt_2_eq$AICc, bt_2_sym$AICc, bt_2_ard$AICc,
                          bt_3_eq$AICc, bt_3_sym$AICc,  bt_3_ard$AICc))$weights,
                c(
(max(as.vector(bt_1_eq$index.mat)[!is.na(as.vector(bt_1_eq$index.mat))])), (max(as.vector(bt_1_sym$index.mat)[!is.na(as.vector(bt_1_sym$index.mat))])), (max(as.vector(bt_1_ard$index.mat)[!is.na(as.vector(bt_1_ard$index.mat))])),
(max(as.vector(bt_2_eq$index.mat)[!is.na(as.vector(bt_2_eq$index.mat))])), (max(as.vector(bt_2_sym$index.mat)[!is.na(as.vector(bt_2_sym$index.mat))])), (max(as.vector(bt_2_ard$index.mat)[!is.na(as.vector(bt_2_ard$index.mat))])),
(max(as.vector(bt_3_eq$index.mat)[!is.na(as.vector(bt_3_eq$index.mat))])), (max(as.vector(bt_3_sym$index.mat)[!is.na(as.vector(bt_3_sym$index.mat))])), (max(as.vector(bt_3_ard$index.mat)[!is.na(as.vector(bt_3_ard$index.mat))])))
                ))
rownames(df_trait)<-c("bt_1_eq", "bt_1_sym", "bt_1_ard",
                          "bt_2_eq",  "bt_2_sym", "bt_2_ard",
                          "bt_3_eq",  "bt_3_sym", "bt_3_ard")
colnames(df_trait)<-c("AICc", "Delta_AICc", "AICcWt", "K_rates")

write.table(df_trait, "Model_selection_multi_state_trait.tsv", sep ="\t")
saveRDS(eval(parse(text = rownames(df_trait)[which.min(df_trait$AICc)])),paste("data_BF_MT.rds", sep =""))
