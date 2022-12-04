setwd("~/PD1/")
library(ggplot2)
library(RColorBrewer)
library(stringr)

library(ggrepel)

################ Color selection #########################


set_colors<- function(nb.cols, bacterial_list){
  

  mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)
  new_colors <- colorRampPalette(brewer.pal(3, "Dark2"))(nb.cols)
  #print (bacterial_list)
  
  mycols <-c(
    "Actinobacteria::Actinobacteria::Bifidobacteriales"="#66C2A5",
    "Actinobacteria::Coriobacteriia.::Coriobacteriales"="#90B392", 
    "Bacteroidetes::Bacteroidia::Bacteroidales"="#A5AB88", 
    "Cyanobacteria::Melainabacteria::Gastranaerophilales"="#BAA47F", 
    "Deferribacteres::Deferribacteres::Deferribacterales"="#6a75bb",
    "Firmicutes::Bacilli::Bacillales"="#E4956C",
    "Firmicutes::Bacilli::Lactobacillales"="#CFA69B",
    "Firmicutes::Clostridia.::Clostridiales"="#BAB5AE",
    "Firmicutes::Erysipelotrichia::Erysipelotrichales"="#919FC6",
    "Patescibacteria::Saccharimonadia::Saccharimonadales"="#F6D832",
    "Proteobacteria::Alphaproteobacteria::Rhodospirillales"="#CF948C",
    "Proteobacteria::Deltaproteobacteria::Desulfovibrionales"="#B3B3B3",
    "Proteobacteria::Gammaproteobacteria::Betaproteobacteriales"="#ABD15C",
    "Tenericutes::Mollicutes::Anaeroplasmatales"="#D5BE9D",
    "Tenericutes::Mollicutes::Mollicutes.RF39"="#EE8F6E",
    "Unassigned::::::"="#C104DD",
    "Verrucomicrobia::Verrucomicrobiae::Verrucomicrobiales"="#DDD83D"
  ) 
  

  day_bac = as.data.frame(names(mycols))
  colnames(day_bac)= "bac"
  bacterial_list = as.data.frame(bacterial_list)
  colnames(bacterial_list) = "bac"
  anti_join_bac = anti_join(bacterial_list, day_bac, by="bac")
  str_col = ""
  
  for ( i in length(anti_join_bac)){
   
    str_col=paste0(anti_join_bac[i,],"=", new_colors[i])
    mycols = append(mycols, str_col)
    str_col = ""
  }

  return (mycols)
}




#########################################################



#Make a class Toxicity 
setClass(Class="Toxicity",
         representation(
           HC="grouped_df",
           Nontox="grouped_df",
           Tox = "grouped_df"
         )
)

#Function to split into healthy, non-tox and tox
extract_labels_data  <- function( df){
  
  HC<-filter(df, Toxicity =="HC")
  Nontox = filter(df, Toxicity =="Non-Tox")
  Tox = filter(df, Toxicity =="Tox")
  return(new("Toxicity",
      HC=HC, 
      Nontox = Nontox,
      Tox = Tox
      ))
  }


#Function to substitute or drop 
subst_characters <- function(df){

  df$variable <- sub("D_0__Bacteria.D_1__","", df$variable)
  df$variable <- sub(".D_2__", "::",df$variable) 
  df$variable <- sub(".D_3__", "::",df$variable) 
  
  
  return (df)
}
###################3 Processing info ############################

x=read.csv("level-4_Day13.csv", header=TRUE)
x1=read.csv("level-4.csv", header=TRUE)
x2=read.csv("level-4_Day20.csv", header=TRUE)
x3 = read.csv("level-41KDay20.csv", header=TRUE)

#melt to make the columns as rows
tx=melt(x)
#tx_day0 = melt(x1)
#tx_day20 = melt(x2)
#tx_day20_1K  = melt(x3)
#Calculate frequency
tx <- tx %>%  group_by(Toxicity, variable) %>% summarise(Frequency = sum(value))
#tx_day20_1K <- tx_day20_1K %>%  group_by(Toxicity, variable) %>% summarise(Frequency = sum(value))

#Extract a df which has calculated frqeuncy 
all_df = extract_labels_data (tx_day20_1K)
all_df = extract_labels_data (tx)

#1500, 3500 and 3500 for Day 0 
#Find the total number of bacteria detected in each group
all_df@HC<- all_df@HC %>% mutate(Total= sum(all_df@HC$Frequency))
all_df@HC<- all_df@HC %>% mutate(Percent = (Frequency/Total)*100)
all_df@Nontox<- all_df@Nontox %>% mutate(Total= sum(all_df@Nontox$Frequency))
all_df@Nontox <- all_df@Nontox %>% mutate(Percent = (Frequency/Total)*100)
all_df@Tox<- all_df@Tox %>% mutate(Total = sum(all_df@Tox$Frequency))
all_df@Tox <- all_df@Tox %>% mutate(Percent = (Frequency/Total)*100)



#all_day0 <- as.data.frame(rbind(HC,Tox,Nontox))


all_df@HC <- subst_characters(all_df@HC)
all_df@Nontox <- subst_characters(all_df@Nontox)
all_df@Tox <- subst_characters(all_df@Tox)

all_day = rbind(all_df@HC, all_df@Nontox , all_df@Tox)
#all_day0['Level1_2'] <- str_split_fixed(all_day0$variable, '.D_3__', 2)
#write.csv(all_day0, file="all_species.txt", quote=FALSE, row.names=FALSE)
#write.csv(mycolors, file="colors", quote=FALSE, row.names=FALSE)

#all_day0$variable <- split(all_day0$variable,split=".D_3__")[1]

colors_nb = length(unique(sort(tx_day20$variable)))
mycols = set_colors(colors_nb, all_df@HC$variable)
pdf("Day13_pd1.pdf")






ggplot(all_day, aes(Toxicity, desc(-Percent), fill=variable, color="black")) +  
  geom_bar(stat="identity",colour="black", alpha=1, width=0.5, size=0.40) +
  scale_fill_manual(values=mycols)  + 
  theme(axis.text.x=element_text( hjust=1,size=9, face =  "bold"),  
       
        panel.background = element_blank(), 
        panel.grid.minor = element_line(color="black", linetype="dotted"),
        panel.border= element_rect(colour = "black", fill=NA, size=1))+ 
       ylab("Mean Relative Abundance (%)") +labs(fill="Bacterial composition")+
  geom_text_repel(aes(label=variable)) 

dev.off()

head(tx)
aggregate(tx$value, by=list(Toxicity=tx$Toxicity), FUN=sum)