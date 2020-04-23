
library(tidyverse)
library(skimr)
library(mice)
library(VIM)
library(GGally)
library(MASS)
library(glmnet)
library(e1071) 
library(rpart)
library(pROC)
library(class)
library(randomForest)
library(FFTrees)
library(caret)


#Heart Disease Data Set
#This dataset contains 14 variables, where data was collected from 303 patients who were admitted to a hospital. The "goal/target" field refers to the presence of heart disease in the patient (1=yes; 0=no). The variables' information is as follows:
#1. age: The person's age in years
#2. sex: The person's sex (1 = male, 0 = female)
#3. cp: The chest pain experienced (Value 1: typical angina, Value 2: atypical angina, Value 3: non-anginal pain, Value 4: asymptomatic)
#4. trestbps: The person's resting blood pressure (mm Hg on admission to the hospital)
#5. chol: The person's cholesterol measurement in mg/dl
#6. fbs: The person's fasting blood sugar (> 120 mg/dl, 1 = true; 0 = false) 
#7. restecg: Resting electrocardiographic measurement (0 = normal, 1 = having ST-T wave abnormality, 2 = showing probable or definite left ventricular hypertrophy by Estes' criteria)
#8. thalach: The person's maximum heart rate achieved during Stress Test (exercise)
#9. exang: Exercise induced angina (1 = yes; 0 = no)
#10. oldpeak: ST depression induced by exercise relative to rest ('ST' relates to positions on the ECG plot)
#11. slope: the slope of the peak exercise ST segment (Value 1: upsloping, Value 2: flat, Value 3: downsloping)
#12. ca: The number of major vessels (0-3) colored by flourosopy 
#13. thal: A blood disorder called thalassemia (1 = normal; 2 = fixed defect; 3 = reversable defect)
#14. goal/target: Heart disease (0 = no, 1 = yes)


Heart<-read.csv("Heart Disease Data.csv", header = TRUE)
glimpse(Heart)
skim(Heart)

# Let's check for NA values
anyNA(Heart)
colSums(is.na(Heart))
#attach(Heart)
md.pattern(Heart)

#We have some categorical values that need to be defined as factors
Heart$sex=factor(Heart$sex)
levels(Heart$sex)=c("Female", "Male")
Heart$cp=factor(Heart$cp)
levels(Heart$cp)=c("Typical Angina", "Atypical Angina", "Non-anginal Pain", "Asymptomatic")
Heart$fbs=factor(Heart$fbs)
levels(Heart$fbs)=c("Below120","Above120")
Heart$restecg=factor(Heart$restecg)
levels(Heart$restecg)=c("Normal","Abnormal","Hypertrophy")
Heart$exang=factor(Heart$exang)
levels(Heart$exang)=c("No","Yes")
Heart$slope=factor(Heart$slope)
levels(Heart$slope)=c("Upslopping","Flat","Downslopping")
Heart$ca=factor(Heart$ca)
Heart$target=factor(Heart$target)
levels(Heart$target)=c("No Heart Disease","Heart Disease")
head(Heart)
glimpse(Heart)

#deleting thal variable as there's missing information on the labels
#source claims 3 factors, but there are 4. 
Heart<-Heart[,-13]


#----------------Split Data-------------
intrain <- createDataPartition(Heart$target, p = 0.6, list = FALSE)
train <-Heart[intrain,]
test <- Heart[-intrain,]
nrow(train)
nrow(test)

#----------Data Visualizations-----------
# this portion was ran before factors were implemented #
#------------------------------------------------------#
train%>%
  keep(is.numeric) %>% 
  gather() %>% 
  ggplot(aes(value)) +
  facet_wrap(~ key, scales = "free") +
  geom_histogram()

train%>%
  keep(is.numeric) %>%                     # Keep only numeric columns
  gather() %>%                             # Convert to key-value pairs
  ggplot(aes(value)) +                     # Plot the values
  facet_wrap(~ key, scales = "free") +   # In separate panels
  geom_density()                         # as density
#-----------------------------------------------------#

#Gender vs. Diseases in this data set
table(Heart$sex,Heart$target)
table(train$sex,train$target)
par(mfrow=c(1,1))
col.target <- c("blue","red")
plot(table(train$sex,train$target),xlab="Gender",ylab="Diagnostics",col=col.target, main=" ")
summary(train)

#Heart Disease vs. Gender vs. Cholesterol
ggplot(train, aes(x=target, y=chol, fill=train$sex))+geom_boxplot( )+
  labs(title="Heart Disease, Gender and Cholesterol Levels", x="Diagnosis",
       y="Cholesterol")+scale_fill_manual(values=c("aquamarine", "pink"),labels=c("Male", "Female"), name="Gender")

#Let's examine the distribution of the continuous variables
#-------Original Goodness of fit---------
#Cholesterol levels 
library(fitdistrplus)
fg<-fitdist(train$chol, "gamma")
fln<-fitdist(train$chol, "lnorm")
fn<-fitdist(train$chol, "norm")
#par(mfrow=c(1,2))
plot.legend<-c("Lognormal","Gamma","Normal")
denscomp(list(fln, fg,fn), legendtext = plot.legend, xlab="Serum Cholesterol Levels") 
qqcomp(list(fln, fg,fn), legendtext = plot.legend) 
cdfcomp(list(fln, fg,fn), legendtext = plot.legend, xlab="Serum Cholesterol Levels")
ppcomp(list(fln, fg,fn), legendtext = plot.legend)

# Resting blood pressure  - tresbps
fg2<-fitdist(train$trestbps, "gamma")
fln2<-fitdist(train$trestbps, "lnorm")
fn2<-fitdist(train$trestbps, "norm")
#par(mfrow=c(1,2))
plot.legend<-c("Lognormal","Gamma","Normal")
denscomp(list(fln2, fg2,fn2), legendtext = plot.legend, xlab="Resting Blood Pressure") 
qqcomp(list(fln2, fg2,fn2), legendtext = plot.legend) 
cdfcomp(list(fln2, fg2,fn2), legendtext = plot.legend, xlab="Resting Blood Pressure")
ppcomp(list(fln2, fg2,fn2), legendtext = plot.legend)

# Resting blood pressure  - thalach
fg2<-fitdist(train$thalach, "gamma")
fln2<-fitdist(train$thalach, "lnorm")
fn2<-fitdist(train$thalach, "norm")
par(mfrow=c(2,2))
plot.legend<-c("Lognormal","Gamma","Normal")
denscomp(list(fln2, fg2,fn2), legendtext = plot.legend, xlab="Maximum Heart Rate") 
qqcomp(list(fln2, fg2,fn2), legendtext = plot.legend) 
cdfcomp(list(fln2, fg2,fn2), legendtext = plot.legend, xlab="Maximum Heart Rate")
ppcomp(list(fln2, fg2,fn2), legendtext = plot.legend)

# Let's take a look at transformed plots with logarithm

plot(density(log(train$chol),kernel="gaussian"),main="Kernel Density of Logarithm of Cholesterol",xlab="",col="deepskyblue2",lwd=2)
plot(density(log(train$thalach),kernel="gaussian"),main="Kernel Density of Logarithm of Resting Blood Pressure",xlab="",col="deepskyblue2",lwd=2)

########################################
#######################################
#######################################
#----------------Models-----------------

#-----Full Model LOG REGRESSION
modfull<- glm(target~., family=binomial(link="logit"), data = train)
summary(modfull)
exp(coef(modfull))
#Roc on the training set
install.packages("Epi")
library(Epi)
#test - ROC(form = target~., data = train, plot = "ROC", lwd = 3, cex=1.5)

#-----Predictions and ROC on testing with LOG reg
probs <- predict(modfull,newdata=test, type='response')
head(probs)
predic <- as.factor(ifelse(probs > 0.5,"Heart Disease","No Heart Disease"))
head(predic)
confusionMatrix(predic, test$target)
#if we always predict NO heart Disease at 55% using a dummy model that always predicts NHD
#is more accurate that our logistic regression
ROC(form = target~., data = test, plot = "ROC", lwd = 3, cex=1.5)


############################################
############################################
#-------LDA Approach
# LDA approach. Is accuracy higher than with log regression? 
lda.mod <- lda(target ~ ., data=train, prior = c(.6,.4))
probs2 = predict(lda.mod, newdata=test)$posterior
threshold = 0.5
h.pred = rep("No Heart Disease", nrow(test))
h.pred[which(probs2[,2] > threshold)] = "Heart Disease"
# Produce a confusion matrix
confusionMatrix(factor(h.pred), test$target)

### ROC curve
probs2 = predict(lda.mod, newdata=test)$posterior

hroc.lda <- roc(test$target,probs2[,2])
auc(hroc.lda) 

plot.roc(test$target, probs2[,2],col="darkblue", print.auc = TRUE,  auc.polygon=TRUE, grid=c(0.1, 0.2),
         grid.col=c("green", "red"), max.auc.polygon=TRUE,
         auc.polygon.col="lightblue", print.thres=TRUE)

#Improvement of LDA given the true thre
lda.mod <- lda(target ~ ., data=train, prior = c(.55,.45))
probs2 = predict(lda.mod, newdata=test)$posterior
threshold = 0.5
h.pred = rep("No Heart Disease", nrow(test))
h.pred[which(probs2[,2] > threshold)] = "Heart Disease"
# Produce a confusion matrix
confusionMatrix(factor(h.pred), test$target)

### ROC curve
probs2 = predict(lda.mod, newdata=test)$posterior

hroc.lda <- roc(test$target,probs2[,2])
auc(hroc.lda) 

plot.roc(test$target, probs2[,2],col="darkblue", print.auc = TRUE,  auc.polygon=TRUE, grid=c(0.1, 0.2),
         grid.col=c("green", "red"), max.auc.polygon=TRUE,
         auc.polygon.col="lightblue", print.thres=TRUE)

###################################
###################################
#---------Tree Approach--------

hrf<-randomForest(target~., data=train)
hrf
predrf<-predict(hrf, newdata = test[-13])

confusionmat=table(test[,13],predrf)
confusionmat

