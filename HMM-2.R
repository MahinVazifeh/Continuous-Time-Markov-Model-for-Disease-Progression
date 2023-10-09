library(momentuHMM)
library(msm)
library(markovchain)
library(Gmisc)
library(ggplot2)
library(ggpubr)
library(cluster)
library(gclus)
library(corrplot)
library(RColorBrewer)
library(foreign)
library(MASS)
library(Hmisc)
library(RColorBrewer)
library(diagram)
library(msmtools)

#devtools::install_github('spedygiorgio/markovchain',force = TRUE)
#install.packages("lme4")

#************************function for making state feature******************
#****************convert EDSS_score_assessed_by_clinician to state from 1 to 3
# based on [0-4]=>1, [4.5-7.5]=>2 , [8-10]=>3
#and save it as state

create_discret<-function(x, x1, x2,y1,y2,y3){
  y<-NULL
  for (i in 1:length(x)){
    if (x[i] <= x1) {
      y[i] <- y1
    } else if (x[i] <= x2) {
      y[i]<- y2
    } else{
      y[i]<- y3
    }
  }
  return(y)
}

#************************function to create Hidden Markov Model******************
create_hmm_model<-function(data, state, time, ptnum, qmatrix, covariates=NULL){
  crude<-crudeinits.msm(state ~ time, ptnum, data=data, qmatrix= qmatrix)
  model_msm <- msm(state ~ time, 
                   subject=ptnum, 
                   data = data,
                   qmatrix = qmatrix, 
                   deathexact = FALSE, 
                   method = "BFGS",
                   control=list(fnscale=10000,maxit=500), 
                   gen.inits=TRUE, 
                   covariates=covariates)
  
  model<-model_msm
  print(model)
  qmat_msm<-qmatrix.msm(model)
  pmat_msm<-pmatrix.msm(model, t=10)
  sojo_msm<-sojourn.msm(model)
  pnext_msm <- pnext.msm(model)
  hazard_msm<-hazard.msm(model)
  
  output<-list(state_table, qmat, crude,  model, qmat_msm, pmat_msm, sojo_msm,pnext_msm, hazard_msm )
  names(output)<-c("state_table", "qmat", "crude",  "model", "qmat_msm", "pmat_msm", "sojo_msm", "pnext_msm", "hazard_msm")
  return(output)
}
#**************************************PREPROCESSING FIRST DATA***********************************
#**********import data*****
data<-read.csv("C:/Users/Mahin Vazifehdan/Desktop/Markov model/3/first_data.csv")
colnames(data)
length(data$ptnum)

#**********state feature*******
state<-create_discret(data$EDSS_score_assessed_by_clinician, 4, 7.5, 1, 2, 3)
state[1:10]

#***************manage the data********
# round data
data_round <-matrix(c(0), length(data$ptnum), 8)
dim(data_round)
for (i in 4:11){
  data_round[,i-3]<-round(data[,i])
}
data_round[1:10,]

#convert state to matrix
state<-matrix(state, length(state),1)
state[1:10,]

#define final data 
final_data <- data.frame(data[,1], data[,2], state, data_round, data[,12],data[,13], data[,14])
final_data[1:10,]

#define column features
names(final_data)<-c("ptnum",	"time",	"State",	
                     "Pyramidal",	"Cerebellar",	"Thronchioencephalic",
                     "Sensitive",	"Sphincteric",	"Visual",	
                     "Mental",	"Deambulation",	"Sex",	"MS.in.pediatric.age",	"Age")

#count if function to detect the ptnum which repeat 1 times
replaced_data <-NULL
for (i in 1:length(final_data$ptnum)){
  target<- sum(final_data$ptnum == final_data$ptnum[i] )
  print(target)
  replaced_data[i]<-target
}
replaced_data

#find the index of data which repeat 1 times, (I told you this before) 
index_repeat <- which(replaced_data != 1 )
index_repeat

#keep the data which repeat more than 1 times 
data_new <- list()
data_new$ptnum<-final_data$ptnum[index_repeat]
data_new$time<-final_data$time[index_repeat]
data_new$State<-final_data$State[index_repeat]
data_new$Pyramidal<-final_data$Pyramidal[index_repeat]
data_new$Cerebellar<-final_data$Cerebellar[index_repeat]
data_new$Thronchioencephalic<-final_data$Thronchioencephalic[index_repeat]
data_new$Sensitive<-final_data$Sensitive[index_repeat]
data_new$Sphincteric<-final_data$Sphincteric[index_repeat]
data_new$Visual<-final_data$Visual[index_repeat]
data_new$Mental<-final_data$Mental[index_repeat]
data_new$Deambulation<-final_data$Deambulation[index_repeat]
data_new$Sex<-final_data$Sex[index_repeat]
data_new$MS.in.pediatric.age<-final_data$MS.in.pediatric.age[index_repeat]
data_new$Age<-final_data$Age[index_repeat]
data_new$ptnum[1:10]

#change initial data to data_new as dataframe
data<-as.data.frame(data_new)  #data is changed to data_new
data[1:10,]

#********************save data*********************
write.csv(data_new, 'C:/Users/Mahin Vazifehdan/Desktop/Markov model/3/updated_data.csv')
length(data_new$ptnum)

#**********************plot frequency table of variables *********************#
mynamestheme <- theme(
  plot.title = element_text(family = "Helvetica", face = "bold", size = (15)),
  legend.title = element_text(colour = "#66c992", face = "bold.italic", family = "Helvetica"),
  legend.text = element_text(face = "italic", colour = "#00dea4", family = "Helvetica"),
  axis.title = element_text(family = "Helvetica", size = (8), colour = "#d3419d"),
  axis.text = element_text(family = "Courier", colour = "#21abcd", size = (5))
)

g1<-ggplot(data, aes(x=state, xlab="State")) + 
  geom_bar(color="black", fill="white")
g1<-g1+mynamestheme+labs(y= "Amount", x = "State")

g2<-ggplot(data, aes(x=Pyramidal)) + 
  geom_histogram(color="black", fill="white", bins=30)
g2<-g2+mynamestheme+labs(y= "Amount", x = "Pyramidal")

g3<-ggplot(data, aes(x=Thronchioencephalic)) + 
  geom_histogram(color="black", fill="white", bins=30)
g3<-g3+mynamestheme+labs(y= "Amount", x = "Thronchioencephalic")

g4<-ggplot(data, aes(x=Sensitive)) + 
  geom_histogram(color="black", fill="white", bins=30)
g4<-g4+mynamestheme+labs(y= "Amount", x = "Sensitive")

g5<-ggplot(data, aes(x=Sphincteric)) + 
  geom_histogram(color="black", fill="white", bins=30)
g5<-g5+mynamestheme+labs(y= "Amount", x = "Sphincteric")

g6<-ggplot(data, aes(x=Visual)) + 
  geom_histogram(color="black", fill="white", bins=30)
g6<-g6+mynamestheme+labs(y= "Amount", x = "Visual")

g7<-ggplot(data, aes(x=Mental)) + 
  geom_histogram(color="black", fill="white", bins=30)
g7<-g7+mynamestheme+labs(y= "Amount", x = "Mental")

g8<-ggplot(data, aes(x=Deambulation)) + 
  geom_histogram(color="black", fill="white", bins=30)
g8<-g8+mynamestheme+labs(y= "Amount", x = "Deambulation")

g9<-ggplot(data, aes(x=Sex)) + 
  geom_bar(color="black", fill="white")
g9<-g9+mynamestheme+labs(y= "Amount", x = "Sex")

g10<-ggplot(data, aes(x=Age)) + 
  geom_histogram(color="black", fill="white", bins=30)
g10<-g10+mynamestheme+labs(y= "Amount", x = "Age")

g11<-ggplot(data, aes(x=MS.in.pediatric.age)) + 
  geom_histogram(color="black", fill="white", bins=30)
g11<-g11+mynamestheme+labs(y= "Amount", x = "MS in Pediatric Age")

figure1 <- ggarrange(g1, g2, g3,g4,g5,g6,g7,g8,g9,g10,g11, 
                     ncol = 3, nrow = 4)

dev.new(width=700, height=350)
figure1

#***********************create correlation between variables*****************#
#extract data for correlation from data
data2<-data[,3:14]
M<-cor(data2)
corrplot(M, type="upper", order="hclust",
         col=brewer.pal(n=11, name="RdYlBu"))

#********************define ordinal regression on state***********#

reg_model<-polr(as.factor(data$State) ~ data$Pyramidal+data$Cerebellar+data$Thronchioencephalic+ data$Sensitive+
                  data$Sphincteric+data$Visual+data$Mental+data$Deambulation+data$Sex+data$MS.in.pediatric.age+data$Age)
summary(reg_model)


#********************fitting Markov Model to calculate Transition Matrix***********
state<-data$State
time<-data$time
Id<-data$ptnum


state_table<-statetable.msm(state, ptnum, data=data)
state_table
qmat<-markovchainFit(state)
qmat <- qmat$estimate@transitionMatrix


#*******************Graphical view of the qmatrix ***********
colnames(qmat)<-rownames(qmat)<-c("Mild", "Mediate", "Severe")
custom_colors = c("#fdd299", "#fab457", "#f99004")
transitionPlot(qmat, txt_start_clr = "black", txt_end_clr = "black", 
               fill_start_box = custom_colors, 
               fill_end_box = custom_colors,
               overlap_add_width = 1.3, type_of_arrow = "gradient"
)

#********************Create Hidden Markov Model***********
#Using msm function from msm package for fitting the Multi-state Markov and hidden Markov models in continuous time
#fitting the HMM model to our clinical data

m1_free<-create_hmm_model(data=data, state=state, time=time, ptnum=Id, qmatrix=qmat)
m1_free

#the error occurred from the previous code because of the data structure, so the data should be change based on hmm model comments
#Only some rows has to be removed because of observing different states at the same time on the same subject at observations which does not maker sense
data_changed<-read.csv('C:/Users/Mahin Vazifehdan/Desktop/Markov model/3/final_data.csv')

state<-data_changed$State
time<-data_changed$time
Id<-data_changed$ptnum

state_table<-statetable.msm(state, ptnum, data=data_changed)
state_table
qmat<-markovchainFit(state)
qmat <- qmat$estimate@transitionMatrix


#Without any influenced feature
m1_free<-create_hmm_model(data=data_changed, state=state, time=time, ptnum=ptnum, qmatrix=qmat)
m1_free

#Considering Age feature as influenced feature
m2_age<-create_hmm_model(data=data_changed, state=state, time=time, ptnum=ptnum, qmat=qmat, covariates=~Age)
m2_age

#Considering Pyramidal feature as influenced feature which made error because of numerical overflow
#(Not Response)
m3_pyramidal<-create_hmm_model(data=data_changed, state=state, time=time, ptnum=ptnum, qmat=qmat, covariates=~Pyramidal
)
m3

#Considering Thronchioencephalic feature as influenced feature which made error because of numerical overflow
#(Not Response)
m4_Thronch<-create_hmm_model(data=data_changed, state=state, time=time, ptnum=ptnum, qmat=qmat, covariates=~Thronchioencephalic
)
m4_Thronch

#Considering Sensitive feature as influenced feature
m5_Sensitive<-create_hmm_model(data=data_changed, state=state, time=time, ptnum=ptnum, qmat=qmat, covariates=~Sensitive
)
m5_Sensitive

#Considering Sphincteric feature as influenced feature
#(Not Response)
m6_Sphincteric<-create_hmm_model(data=data_changed, state=state, time=time, ptnum=ptnum, qmat=qmat, covariates=~Sphincteric
)
m6_Sphincteric

#Considering Visual feature as influenced feature
m7_Visual<-create_hmm_model(data=data_changed, state=state, time=time, ptnum=ptnum, qmat=qmat, covariates=~Visual
)
m7_Visual 

#Considering Mental feature as influenced feature
#(Not Response)
m8_Mental<-create_hmm_model(data=data_changed, state=state, time=time, ptnum=ptnum, qmat=qmat, covariates=~Mental
)
m8_Mental

#Considering Deambulation feature as influenced feature
m9_Deambulation<-create_hmm_model(data=data_changed, state=state, time=time, ptnum=ptnum, qmat=qmat, covariates=~Deambulation
)
m9_Deambulation

#****************************************Plot Hidden Markov Model*************************#
#function for limit decimals
specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))
colnames(m1_free$pmat_msm)<-rownames(m1_free$pmat_msm)<-c("Mild", "Mediate", "Severe")
#plot the transition matrix(pmat_msm)
par(mfrow = c(1,1), mar= c(4,4,4,4))
plotmat(specify_decimal(m1_free$pmat_msm, 4), curve = 0.2,               
        lwd = 1.5, 
        box.lwd = 2, 
        cex.txt = 0.8,
        box.size = 0.1, 
        box.type = "hexa", 
        box.prop = 0.6, main = "Covariate=Deambulation", 
        box.col = c("#fdd299", "#fab457", "#f99004")
)
colnames(m9_Deambulation$pmat_msm)<-rownames(m9_Deambulation$pmat_msm)<-c("Mild", "Mediate", "Severe")
plotmat(specify_decimal(m9_Deambulation$pmat_msm, 4), curve = 0.2,               
        lwd = 1.5, 
        box.lwd = 2, 
        cex.txt = 0.5,
        box.size = 0.1, 
        box.type = "hexa",
        box.prop = 0.6, main = "Covariate=Deambulation", 
        box.col = c("#fdd299", "#fab457", "#f99004")
)

#plot estimated and real value of state
model_msm_m1_free <- msm( state ~ time, subject=ptnum, data = data_changed,
                          qmatrix = qmat, deathexact = FALSE, method = "BFGS",
                          control=list(fnscale=10000,maxit=500), gen.inits=TRUE)


colnames(model_msm_m1_free$Qmatrices$logbaseline)<-c("Mild", "Mediate", "Severe")
prev = prevalence.msm( model_msm_m1_free, covariates = 'mean', ci = 'normal')
gof = prevplot(x = model_msm_m1_free, prev.obj = prev, ci = TRUE, M = TRUE )

model_msm_m9_Deambulation <- msm( state ~ time, subject=ptnum, data = data_changed,
                                  qmatrix = qmat, deathexact = FALSE, method = "BFGS",
                                  control=list(fnscale=10000,maxit=500), gen.inits=TRUE)


prev = prevalence.msm( model_msm_m1_free, covariates = 'mean', ci = 'normal')
gof = prevplot( x = model_msm_m1_free, prev.obj = prev, ci = TRUE, M = TRUE )
