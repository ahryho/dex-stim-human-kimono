
library(data.table)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(dplyr)

data.dir.pre   <- "/Users/anastasiia_hry/bio/datasets/kimono/"
kimono.res.fn  <- paste0(data.dir.pre, "output/kimono_res_veh.csv") # experiment nr 2
mapping.tbl.fn <- paste0(data.dir.pre, "/mapping/mapping_ilmn_ens.csv")

# 1. Load Kimono results and filter them

kimono.res   <- fread(kimono.res.fn)
db_ilmn_gene <- read.csv2(paste0(output.data.pre, "/mapping/mapping_ilmn_ens.csv"))

# 2. Look at some numbers

kimono.net <- kimono.res[value != 0, ] # 39779

paste("Number of successful models:",
      length(unique(kimono.res$target))) # 2696
paste("Number of unique predictors:", 
      length(unique(kimono.res$predictor))) #2583

mergeddt <- merge(kimono.res[,.(target)], db_ilmn_gene, 
                  by.x = "target", by.y = "Ensemble_ID", all.x=T) %>%  unique(.)
paste("Number of ILMN ID that have multiple Genes associated:", 
      mergeddt[,.N, by=target][order(-N)][N>1,] %>%  nrow(.) )# there are 76 ILMN that match to two Genes

mergeddt[target %in% mergeddt[,.N, by=target][order(-N)][N>1,target],] #these are not unique
paste("Number of unique genes:", length(unique(mergeddt$Gene))) #5663

# 3. Create network

beta.trsh <- 0.2
rsq.trsh  <- 0.01

network <- kimono.net %>% filter((value > beta.trsh) | (value < (-beta.trsh))) %>%
  filter(r_squared > rsq.trsh) %>%
  filter(predictor != '(Intercept)') %>% setDT

# 4. Inspect associations
# cnt <- network[, .(target, relation)] %>%  unique() %>% .[, .N, by = target] #1584
# agene <- cnt[order(-N)] %>% .[N == 3, target]

kimono.res$target %>%  unique %>%  length # 2696
network$target %>%  unique %>%  length # 1584

network[relation == "expr_pheno", predictor] %>% unique 
network[relation == "methyl_expr", predictor] %>% unique
network[relation == "prior_expr", predictor] %>% unique

ggplot(network, aes(relation)) + 
  geom_bar(stat = "count", fill = "purple") + 
  # scale_fill_brewer(palette = "Blues") +
  geom_text(stat = "count", aes(label=..count..), vjust=1.6, color="white",
            position = position_dodge(0.9), size=3.5)

network[,.N, by=relation]
network[, .N, by=predictor][order(-N)]

# retained fraction after filtering 
nrow(network) / nrow(kimono.net) # 0.0586

paste("Number of unique genes:", 
      length(unique(network$target))) #1584
paste("Number of unique predictors:", 
      length(unique(network$predictor))) # 599

# unique predictors
unique(network, by=c("relation", "predictor")) %>% .[,.N, by=relation]
# relation   N
# 1:  expr_pheno   5
# 2:  prior_expr 548
# 3: methyl_expr  46

## performance and beta values

dtmoni <- kimono.net
library(ggthemes)
# pdf(file=paste0(Sys.Date(),"_tmp1.pdf"), height=3, width=3)
dtmoni[,.(target,r_squared)] %>%unique(.) %>% ggplot(.,aes(r_squared)) +
  geom_density()+ geom_histogram(bins=100) + theme_tufte() +
  geom_vline(xintercept=rsq.trsh, color="red")+
  labs(title="R2 of all successfull models")
ggplot(dtmoni,aes(value)) + theme_tufte()+ labs(title="value of coefficients") +
  geom_vline(xintercept=0.02, color="red")+
  geom_vline(xintercept=-0.02, color="red")+
  # ylim(0, 1e5)+
  geom_histogram(bins=100)
ggplot(network,aes(value)) + theme_tufte()+ labs(title="value of coefficients") +
  geom_vline(xintercept=0.02, color="red")+
  geom_vline(xintercept=-0.02, color="red")+
  # ylim(0, 1e5)+
  geom_histogram(bins=100)
ggplot(network, aes(relation)) + geom_bar(fill="lightblue",width = 0.5) + theme_tufte() +theme(axis.text.x = element_text(angle = 45, hjust = 1))


# 3. Create network
library(purrr)
library(igraph)

vertices <- unique(c(network$predictor, network$target))
edges    <- data.frame(from = network$predictor,
                       to = network$target,
                       value = network$value)

                        # performance=network$performance) 

# network
g <- graph_from_data_frame(edges, vertices, directed = FALSE)
print(g)

network.distr <- degree_distribution(g)
plot(network.distr)

#betweeness
betweenes_vertex <- betweenness(g, directed=F, weights=NA)
# plot all betweeness
tmp <- as.data.frame(t(t(sort(betweenes_vertex, decreasing = T)))) %>% dplyr::rename(betweenness = V1)
ggplot(tmp, aes(log(betweenness))) + geom_density() + geom_histogram(bins = 100) 

# top 100 genes
NoI <- t(t(head(sort(betweenes_vertex, decreasing = T), n = 100)))
NoI <- setDT(as.data.frame(NoI), keep.rownames = T)[]
# GoI <- merge(NoI[18:nrow(NoI),], db_ilmn_gene, by.x="rn", by.y="ILMN")
GoI <- NoI
colnames(GoI) <- c("Ensemble_ID", "Betweeness")
GoI <- GoI[order(-Betweeness)]

# top10
topilmn <- GoI[1:10, Ensemble_ID]
subnet <- network[((target %in% topilmn) | (predictor %in% topilmn)) & predictor!="(Intercept)",] %>% setDT

# 4. Plot network
subnet <- network
plot_network <- function(subnet){
  subnet <- data.frame(lapply(subnet, as.character), stringsAsFactors=FALSE)
  subnet$target <- gsub("_", ".", subnet$target)
  subnet$predictor <- gsub("_", ".", subnet$predictor)
  subnet$relation <- gsub("_",".", subnet$relation)
  
  lable <- unique(c(paste0("expr_",subnet[,1]),paste0(subnet$relation,"_",subnet[,2])))
  
  actors <- data.frame(name=do.call(rbind, strsplit(lable,"_") )[,2],
                       omic= do.call(rbind, strsplit(lable,"_") )[,1]
  )
  
  
  relations <- data.frame(from=subnet$predictor,
                          to=subnet$target,
                          value=subnet$value,
                          performance=subnet$r_squared) 
  
  actors <- unique(actors) %>%  setDT
  actors[omic=="prior.expr", omic:="expr"]
  actors <- unique(actors)
  g <- graph_from_data_frame(relations, directed=FALSE, vertices=actors)
  
  deg2 <- igraph::degree(g, mode="all")
  
  # plot(g, vertex.label.color="black", vertex.size=10, vertex.label=NA , vertex.label.dist=1.5)
  # library(qgraph)
  
  e <- get.edgelist(g,names=FALSE)
  # l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(g),  area=1000*(vcount(g)^2),repulse.rad=100+(vcount(g)^30))
  

  actors[, edges:=igraph::degree(g)]
  actors[, mynames:=name]
  
  
  actors$id <- 1:nrow(actors)
  actors <-actors[order(id)]
  summary(actors$edges)
  
  
  # edge_threshold <- (as.numeric(sub('.*:', '', summary(actors$edges)[5])) +as.numeric(sub('.*:', '', summary(actors$edges)[6])))/2
  
  edge_threshold<-100
  actors[edges>edge_threshold,]
  actors[,.N, by=edges,][order(edges)]
  actors[edges < edge_threshold  , mynames:=NA]
  
  colrs <-c(brewer.pal(4, "Set2")[c(1:4)])
  mycol=colrs[2:4]
  actors[omic=="expr", mycolor:=colrs[2]]
  actors[omic=="expr.pheno", mycolor:=colrs[3]]
  actors[omic=="methyl.expr", mycolor:=colrs[4]]
  plot(g,
       layout=layout.fruchterman.reingold(g),
       vertex.frame.color= adjustcolor("black", .4)	, 
       vertex.size=1+(log(deg2)*2),
       vertex.color=actors$mycolor ,
       edge.color =  adjustcolor("grey", .8),
       edge.curved=.1,
       vertex.label = actors$mynames,
       vertex.label.color= adjustcolor("black", .8),
       vertex.label.cex = 0.5,    asp = 1 ,
       # vertex.label.family = "Times",
       edge.width=E(g)$weight*3)
  legend(x=-1, y=-0.5,c("gene", "tissue"),
         pch=21, col="#777777", pt.bg=mycol, pt.cex=2, cex=.8, bty="n", ncol=1)
}

# whole network
plot_network(network)
