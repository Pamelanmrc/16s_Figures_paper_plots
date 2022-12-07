setwd("~/16S/CTLA4/")
setwd("~/16S/PD1")
library(ggplot2)
library(RColorBrewer)
library(stringr)
library(stringr)
library(reshape2)
library(ggrepel)

################ Color selection #########################

# bacterial_list
set_colors<- function(nb.cols){
  

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
    "D_0__Bacteria;D_1__Actinobacteria;D_2__Actinobacteria;D_3__Corynebacteriales" = "slateblue3",
    "Verrucomicrobia::Verrucomicrobiae::Verrucomicrobiales"="#DDD83D"
  ) 
  

  #day_bac = as.data.frame(names(mycols))
  #colnames(day_bac)= "bac"
  #bacterial_list = as.data.frame(bacterial_list)
  #colnames(bacterial_list) = "bac"
  #anti_join_bac = anti_join(bacterial_list, day_bac, by="bac")
  #str_col = ""
  
  #for ( i in length(anti_join_bac)){
   
   # str_col=paste0(anti_join_bac[i,],"=", new_colors[i])
  #  mycols = append(mycols, str_col)
  # str_col = ""
  #}

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

#Calculate Frequency and percent
calculate_freq <- function(df){
  df@HC<- df@HC %>% mutate(Total= sum(df@HC$Frequency))
  df@HC<- df@HC %>% mutate(Percent = (Frequency/Total)*100)
  df@Nontox<- df@Nontox %>% mutate(Total= sum(df@Nontox$Frequency))
  df@Nontox <- df@Nontox %>% mutate(Percent = (Frequency/Total)*100)
  df@Tox<- df@Tox %>% mutate(Total = sum(df@Tox$Frequency))
  df@Tox <- df@Tox %>% mutate(Percent = (Frequency/Total)*100)
  
  return (df)
}


#Replace a few characters
replace_chr <- function(df){
  
  df@HC <- subst_characters(df@HC)
  df@Nontox <- subst_characters(df@Nontox)
  df@Tox <- subst_characters(df@Tox)
  return(df)
}

###################3 Processing info ############################

getwd()
x1=read.csv("core_metrics_Day0_1K/level-4_Day0.csv", header=TRUE)
x2=read.csv("core_metrics_Day13_1K/level-4_Day13.csv", header=TRUE)
x3=read.csv("core_metrics_Day21_1K/level-4_Day20.csv", header=TRUE)
x1=read.csv("level-4_Day0.csv", header=TRUE)
x2=read.csv("level-4_Day13.csv", header=TRUE)
x3=read.csv("level-4_Day20.csv", header=TRUE)



x=read.csv("CTLA4_Day0_tax.csv", header=TRUE)
x=read.csv("level-4_Day13.csv", header=TRUE)
x=read.csv("level-4_Day20.csv", header=TRUE)

#melt to make the columns as rows

tx_day0 = melt(x1)
tx_day13 = melt(x2)
tx_day20 = melt(x3)

#Calculate frequency

tx_day0 <- tx_day0 %>%  group_by(Toxicity, variable) %>% summarise(Frequency = sum(value))

tx_day13 <- tx_day13 %>%  group_by(Toxicity, variable) %>% summarise(Frequency = sum(value))
tx_day20 <- tx_day20 %>%  group_by(Toxicity, variable) %>% summarise(Frequency = sum(value))

#Find the total number of bacteria detected in each group
all_df_0 = extract_labels_data (tx_day0) 
all_df_0 = calculate_freq(all_df_0)
all_df_13 = extract_labels_data (tx_day13) 
all_df_13 = calculate_freq(all_df_13)
all_df_20 = extract_labels_data (tx_day20) 
all_df_20 = calculate_freq(all_df_20)




#For heatmap take bacteria < 0.5aall_df_0@HC[all_df_0@HC$Percent > 0.5,]
all_df_20@HC= all_df_20@HC[all_df_20@HC$Percent > 0.5,]
all_df_13@HC = all_df_13@HC[all_df_13@HC$Percent > 0.5,]
all_df_20@HC = all_df_20@HC[all_df_20@HC$Percent > 0.5,]

all_df_0@Nontox = all_df_0@Nontox[all_df_0@Nontox$Percent > 0.5,]
all_df_20@Nontox= all_df_20@Nontox[all_df_20@Nontox$Percent >  0.5,]
all_df_13@Nontox = all_df_13@Nontox[all_df_13@Nontox$Percent > 0.5,]


all_df_0@Tox = all_df_0@Tox[all_df_0@Tox$Percent > 0.5,]
all_df_20@Tox = all_df_20@Tox[all_df_20@Tox$Percent > 0.5,]
all_df_13@Tox = all_df_13@Tox[all_df_13@Tox$Percent > 0.5,]
#
#hc3 <- as.data.frame(unique(sort(all_df_20@HC$variable)))
#colnames(hc3) ="bac"

#Change the "D1__ etc" to "::"
all_df_20 <- replace_chr(all_df_20)
all_df_0 <- replace_chr(all_df_0)
all_df_13 <- replace_chr(all_df_13)


#all_day_0 <- all_day_0 %>% mutate(Day="Day0")
#all_day_20 <- all_day_20 %>% mutate(Day="Day20")

#day =rbind(all_day_0, all_day_20)
#day$ann = paste0(day$Day,"_",day$Toxicity)

#day$scale = scale(day$Frequency, center = TRUE, scale =TRUE)
#day$variable = str_split_fixed(day$variable, "::",3)[,1]
#heatm <-day %>%  group_by(variable, Toxicity, ann) %>% summarize(sum=sum(Frequency))
#heatm<- heatm[,names(heatm) %in% 
               # c("variable", "ann", "sum")]
#heatm$scale <- (heatm$sum - mean(heatm$sum)) / sd(heatm$sum)
#day0<-heatm[grep("Day0", heatm$ann),]


#Check for the same bacteria:
all_df_0@Nontox$variable %in% all_df@HC$variable
all_df@Tox$variable %in% all_df@HC$variable


#all_day0['Level1_2'] <- str_split_fixed(all_day0$variable, '.D_3__', 2)
#write.csv(all_day0, file="all_species.txt", quote=FALSE, row.names=FALSE)
#write.csv(mycolors, file="colors", quote=FALSE, row.names=FALSE)

#all_day0$variable <- split(all_day0$variable,split=".D_3__")[1]

colors_nb = length(unique(sort(tx$variable)))
#, all_df@HC$variable


#inn_tox_non = inner_join(all_df_20@Nontox, all_df_20@Tox, by="variable") %>%  select(variable,Frequency.x, Percent.x,
                                                                          # Frequency.y, Percent.y    )
inn_tox_non = inner_join(all_df_13@Nontox, all_df_13@Tox, by="variable") %>%  select(variable,Frequency.x, Frequency.y)
colnames(inn_tox_non) =c("variable", "Non_tox","Tox")        
inn_hc_tox_non = inner_join(all_df_13@HC, inn_tox_non , by="variable") %>%  select(variable,Frequency, Non_tox, Tox)
colnames(inn_hc_tox_non)[3] ="HC"
inn_hc_tox_non <- inn_hc_tox_non %>%   filter(variable != "TimePoint" )


inn_hc_tox_non_0 = inn_hc_tox_non
#inn_tox_non = inner_join(all_df_20@Nontox, all_df_20@Tox, by="variable") %>%  select(variable,Frequency.x, Frequency.y)

colnames(inn_tox_non) =c("variable", "Non_tox","Tox")        
#inn_hc_tox_non = inner_join(all_df_0@HC, inn_tox_non , by="variable") %>%  select(variable,Frequency, Non_tox, Tox)
#colnames(inn_hc_tox_non)[3] ="HC"
#remove.packages("reshape2")

inn_hc_tox_non = inn_hc_tox_non[,names(inn_hc_tox_non) %in% c("variable", "HC", "Non_tox", "Tox")]
inn_hc_tox_non$variable = stringr::str_split_fixed(inn_hc_tox_non$variable, "::",3)[,1]
#inn_hc_tox_non = inn_hc_tox_non[,colnames(inn_hc_tox_non) %in% c("variable", "HC","Non_tox","Tox")]
inn_hc_tox_non = inn_hc_tox_non %>%  group_by(variable) %>% summarise_all(sum)
inn_hc_tox_non = tibble::column_to_rownames(inn_hc_tox_non,"variable")

mycols = set_colors(colors_nb)
pdf("Day20_CTLA4_heatmap.pdf")
inn_hc_tox_non = inn_hc_tox_non[rownames(inn_hc_tox_non)!= "TimePoint",]

pheatmap(inn_hc_tox_non, scale="row", cluster_cols =FALSE, main="Day 20",
                   border_color = "white",cellwidth=35,size = 20,legend = FALSE,
                   colorRampPalette(c("steelblue3", "white", "lightsalmon3"))(60)) 
  
dev.off()

#heat <- as.data.frame(cbind(all_df_20@Nontox$Frequency ,  all_df_20@HC$Frequency, all_df_20@Tox$Frequency))
#rownames(heat) =all_df_20@Tox$variable


#heat = tibble::rownames_to_column(heat, var="bacteria")
#heat$bacteria = str_split_fixed(heat$bacteria, "::",3)[,1]
#Gather all the bacteria together
#Then drop the Unassigned
heatm <-heat %>%  group_by(bacteria) %>% summarize_all(sum)
heatm<- tibble::column_to_rownames(heatm, var="bacteria")
drop_row = grep("Unassigned", rownames(heatm))
heatm = heatm[-c(as.integer(drop_row)),]


pdf("PD1_Day20.pdf")
pheatmap(heatm, scale="row", cluster_cols =FALSE, main="Day 20",
                   border_color = "white",cellwidth=35,size = 20,
                   colorRampPalette(c("steelblue3", "white", "lightsalmon3"))(60))
dev.off()
   


all_day_0<- rbind(all_df_0@HC,  all_df_0@Nontox, all_df_0@Tox)

all_day_13<- rbind(all_df_13@HC,  all_df_13@Nontox, all_df_13@Tox)
all_day_20<- rbind(all_df_20@HC,  all_df_20@Nontox, all_df_20@Tox)

unique(sort(all_day$variable))
all_day_0 <- all_day_0 %>%   filter(variable != "TimePoint" )
all_day_20 <- all_day_20 %>%   filter(variable != "TimePoint" )
all_day_13 <- all_day_13 %>%   filter(variable != "TimePoint" )



getwd()
pdf("Day0_CTLA4_barplot.pdf")
pdf("Day0_PD1_barplot.pdf")
pdf("Day20_PD1_barplot.pdf")
pdf("Day13_PD1_barplot.pdf")


plot_20<- ggplot(all_day_20, aes(Toxicity, desc(-Percent), fill=variable, color="black")) +  
  geom_bar(stat="identity",colour="black",show.legend = FALSE,  alpha=1, width=0.5, size=0.40) +
  scale_fill_manual(values=mycols)  + 
  theme(axis.text.x=element_text( hjust=1,size=9, face =  "bold"),  
        panel.background = element_blank(), 
        panel.grid.minor = element_line(color="black", linetype="dotted"),
        panel.border= element_rect(colour = "black", fill=NA, size=1))+ 
  ylab("Mean Relative Abundance (%)") +labs(fill="Bacterial composition")

plot_0 <- ggplot(all_day_0, aes(Toxicity, desc(-Percent), fill=variable, color="black")) +  
  geom_bar(stat="identity",colour="black",show.legend = FALSE, alpha=1, width=0.5, size=0.40) +
  scale_fill_manual(values=mycols)  + 
  theme(axis.text.x=element_text( hjust=1,size=9, face =  "bold"),  
        panel.background = element_blank(), 
        panel.grid.minor = element_line(color="black", linetype="dotted"),
        panel.border= element_rect(colour = "black", fill=NA, size=1))+ 
  ylab("Mean Relative Abundance (%)") +labs(fill="Bacterial composition")

plot_13<- ggplot(all_day_13, aes(Toxicity, desc(-Percent), fill=variable, color="black")) +  
  geom_bar(stat="identity",colour="black",show.legend = FALSE, alpha=1, width=0.5, size=0.40) +
  scale_fill_manual(values=mycols)  + 
  theme(axis.text.x=element_text( hjust=1,size=9, face =  "bold"),  
        panel.background = element_blank(),
        panel.grid.minor = element_line(color="black", linetype="dotted"),
        panel.border= element_rect(colour = "black", fill=NA, size=1))+ 
       ylab("Mean Relative Abundance (%)") +labs(fill="Bacterial composition")

legend<-get_legend(plot_0)
prow <- plot_grid(
  plot_0 + theme(legend.position="none"),
  plot_13 + theme(legend.position="none"),
  plot_20 + theme(legend.position="none"),
  align = 'vh',
  hjust = -1,
  nrow = 1
)



pdf("all_plots.pdf")
plot_grid(prow)
dev.off()
