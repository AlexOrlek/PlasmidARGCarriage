setwd("C:/Users/alexo/Dropbox/Oxford_2015/thesis/CORRECTIONS/Chapter_5/resistance_prediction")
library(mgcv)
library(data.table)
set.seed(42)

finaldftrunc<-read.table('data/Appendix_D_Table_6_splitbetalactam_truncated.tsv',header=TRUE,sep='\t',stringsAsFactors = TRUE,quote = "",comment.char = "")
#outcomeclasses<-c('betalactam','aminoglycoside','sulphonamide','tetracycline','trimethoprim','phenicol','macrolide','quinolone')
outcomeclasses<-c('aminoglycoside','betalactam_carbapenem','betalactam_ESBL','betalactam_other','macrolide','phenicol','quinolone','sulphonamide','tetracycline','trimethoprim')

#need to convert categorical varaibles to factors
finaldftrunc$bacmetbinary<-as.factor(finaldftrunc$bacmetbinary)
finaldftrunc$vfdbbinary<-as.factor(finaldftrunc$vfdbbinary)
finaldftrunc$integronbinary<-as.factor(finaldftrunc$integronbinary)
finaldftrunc$conj<-factor(finaldftrunc$conj, ordered = FALSE,levels = c("non-mob", "mob", "conj"))
finaldftrunc$geographies<-factor(finaldftrunc$geographies, ordered = FALSE,levels = c("WB_HI", "WB_MI", "China", "United States", "EU", "other"))
finaldftrunc$isolationsources<-factor(finaldftrunc$isolationsources, ordered = FALSE,levels = c("human", "livestock","other"))
finaldftrunc$reptypes<-factor(finaldftrunc$reptypes, ordered = FALSE,levels = c("untyped", "single", "multi"))
finaldftrunc$taxa<-factor(finaldftrunc$taxa, ordered = FALSE,levels = c("Enterobacteriaceae", "Proteobacteria", "Firmicutes","other"))
for (outcomeclass in outcomeclasses) {
  finaldftrunc[,gsub('%s',outcomeclass,'outcome%s')]<-as.factor(finaldftrunc[,gsub('%s',outcomeclass,'outcome%s')])
  finaldftrunc[,gsub('%s',outcomeclass,'otherresgenesbinary%s')]<-as.factor(finaldftrunc[,gsub('%s',outcomeclass,'otherresgenesbinary%s')])
}


###functions
writegamcheck<-function(model,outcomeclass,gamcheckfile,gamcheckplotfile) {
  if (is.null(gamcheckfile)) {
    #print gam.check output to terminal
    gcouts<-capture.output(gam.check(model[[outcomeclass]]))
    for (gcout in gcouts) {
      cat(gcout,sep='\n')
    }
  } else {
    #write gam.check output to file
    pdf(gamcheckplotfile)
    gcouts<-capture.output(gam.check(model[[outcomeclass]]))
    dev.off()
    for (gcout in gcouts) {
      if (grepl('^s\\(',gcout)==TRUE || grepl('^\\s+k',gcout)==TRUE) {
        gcout<-gsub("*\\s+", "\t", gcout)
      }
      cat(gcout,sep='\n',append = TRUE,file = gamcheckfile)
    }
  }
}

writeconcurvity<-function(model,outcomeclass,concurvityfile) {
  couts<-capture.output(concurvity(model[[outcomeclass]],full=TRUE))
  for (cout in couts) {
    cout<-gsub("*\\s+", "\t", cout)
    cat(cout,sep='\n',append = TRUE,file = concurvityfile)
  }
  cat('',sep='\n',append = TRUE,file = concurvityfile)
  couts<-capture.output(concurvity(model[[outcomeclass]],full=FALSE))
  for (cout in couts) {
    cout<-gsub("*\\s+", "\t", cout)
    cat(cout,sep='\n',append = TRUE,file = concurvityfile)
  }
}

###


#INITIAL MODELLING (test different possible model structures with largest outcome class)

outcomeclass<-'aminoglycoside'

###failed model - inadequate df
frm <- formula(gsub('%s',outcomeclass,'outcome%s~s(logplasmidsize)+s(isscore)+s(numotherresgenes%s)+s(coldateimputed)+integronbinary+bacmetbinary+conj+reptypes+taxa+vfdbbinary+geographies+isolationsources'))
model1<-gam(frm,family='binomial',data=finaldftrunc,method = 'ML')
#Error in smooth.construct.tp.smooth.spec(object, dk$data, dk$knots) : 
#  A term has fewer unique covariate combinations than specified maximum degrees of freedom

###initial model that works
frm <- formula(gsub('%s',outcomeclass,'outcome%s~s(logplasmidsize,k=5)+s(isscore,k=5)+s(numotherresgenes%s,k=5)+s(coldateimputed,k=3)+integronbinary+bacmetbinary+conj+reptypes+taxa+vfdbbinary+geographies+isolationsources'))
model1.1<-gam(frm,family='binomial',data=finaldftrunc,method = 'ML')
gam.check(model1.1)  #logplasmidsize basis dimension inadequate

###increasing isscore and logplasmidsize basis dimensions fails to resolve low p-value (see gam.check output)
frm <- formula(gsub('%s',outcomeclass,'outcome%s~s(logplasmidsize,k=7)+s(isscore,k=7)+s(numotherresgenes%s,k=5)+s(coldateimputed,k=3)+integronbinary+bacmetbinary+conj+reptypes+taxa+vfdbbinary+geographies+isolationsources'))
model1.2.1<-gam(frm,family='binomial',data=finaldftrunc,method = 'ML')
gam.check(model1.2.1)

frm <- formula(gsub('%s',outcomeclass,'outcome%s~s(logplasmidsize,k=10)+s(isscore,k=10)+s(numotherresgenes%s,k=5)+s(coldateimputed,k=3)+integronbinary+bacmetbinary+conj+reptypes+taxa+vfdbbinary+geographies+isolationsources'))
model1.2.2<-gam(frm,family='binomial',data=finaldftrunc,method = 'ML')
gam.check(model1.2.2)

###model with logplasmidsize removed
frm <- formula(gsub('%s',outcomeclass,'outcome%s~s(isscore,k=5)+s(numotherresgenes%s,k=5)+s(coldateimputed,k=3)+integronbinary+bacmetbinary+conj+reptypes+taxa+vfdbbinary+geographies+isolationsources'))
model1.3<-gam(frm,family='binomial',data=finaldftrunc,method = 'ML')
gam.check(model1.3)

AIC(model1.1,model1.2.1,model1.2.2,model1.3)
#                 df      AIC
#model1.1   25.52334 3873.431
#model1.2.1 26.81236 3869.305
#model1.2.2 33.77410 3845.523
#model1.3   24.48570 3871.260


#Because logplasmidsize basis dimensionality remains inadequate at k=10, I will choose model1.1
summary(model1.1)
anova.gam(model1.1) #all terms except coldateimputed are significant at p<0.05; coldateimputed is of interest for other outcomes - I will not simplify model structure

#need to perform model checking across all classes...

###ALL CLASS MODELLING
dfcolnames<-colnames(finaldftrunc)

#baseline model created for all classes
modellist<-list()
for (outcomeclass in outcomeclasses) {
  print(outcomeclass)
  #need to temporarily rename numotherresgenes%s column so that all models have same predictor names (so mgcviz works)
  predictorindx<-which(dfcolnames==gsub('%s',outcomeclass,'numotherresgenes%s'))
  colnames(finaldftrunc)[predictorindx]<-'numotherresgenes'
  frm <- formula(gsub('%s',outcomeclass,'outcome%s~s(logplasmidsize,k=5)+s(isscore,k=5)+s(numotherresgenes,k=5)+s(coldateimputed,k=3)+integronbinary+bacmetbinary+conj+reptypes+taxa+vfdbbinary+geographies+isolationsources'))
  #frm <- formula(gsub('%s',outcomeclass,'outcome%s~s(logplasmidsize,k=10)+s(isscore,k=10)+s(numotherresgenes,k=5)+s(coldateimputed,k=3)+integronbinary+bacmetbinary+conj+reptypes+taxa+vfdbbinary+geographies+isolationsources'))
  modellist[[outcomeclass]]<-gam(frm,family='binomial',data=finaldftrunc,method = 'ML')
  colnames(finaldftrunc)[predictorindx]<-gsub('%s',outcomeclass,'numotherresgenes%s')
}


#gam.check output
gamcheckfile<-'output_adjusted/exploratory/gamcheck_model1.1.txt'
concurvityfile<-'output_adjusted/exploratory/concurvity_model1.1.txt'
cat('gam.check output\n',file=gamcheckfile,sep = '\n',append = FALSE)
cat('concurvity output\n',file=concurvityfile,sep = '\n',append = FALSE)
for (outcomeclass in outcomeclasses) {
  gamcheckplotfile<-gsub('%s',outcomeclass,'output_adjusted/exploratory/gamcheckplots_model1.1_%s.pdf')
  cat(gsub('%s',outcomeclass,'\n###gam.check for outcome: %s binary'),file=gamcheckfile,sep = '\n',append = TRUE)
  writegamcheck(modellist,outcomeclass,gamcheckfile,gamcheckplotfile)
  cat(gsub('%s',outcomeclass,'\n###concurvity for outcome: %s binary'),file=concurvityfile,sep = '\n',append = TRUE)
  writeconcurvity(modellist,outcomeclass,concurvityfile)
}


#N.B gam.check residual check p-values are stochastic (but overall, inadequacy of logplasmidsize basis dimension is supported)


# 
# ###TESTING
# #trying to calculate vifs for parametric terms using vif.gam
# library(mgcv.helper)
# test<-modellist1.1[[outcomeclasses[1]]]
# summary(test)
# concurvity(test,full = TRUE)
# concurvity(test,full = FALSE)
# mgcv.helper::vif.gam(test)
# #Error in `[.data.frame`(X, , selected_col) : undefined columns selected
# mgcv.helper::vif.gam(test$model)
# #Error in dimnames(covmat) <- list(name, name) :
# #  attempt to set an attribute on NULL
# #In addition: Warning message:
# #  In mgcv::summary.gam(object) :
# #  p-values for any terms that can be penalized to zero will be unreliable: refit model to fix this.
# mgcv.helper::vif.gam(model.matrix(test$model))
# #Error in eval(predvars, data, env) :
# #  object 'outcomeaminoglycoside' not found
# #doesn't work? - also tried fitting to gam with only parametric terms and still fails; maybe function only applies to numeric predictors as in example??
# 
# #try vif with glm rather than gam model
# library(car)
# outcomeclass<-outcomeclasses[1]
# frm <- formula(gsub('%s',outcomeclass,'outcome%s~logplasmidsize+isscore+numotherresgenes%s+coldateimputed+integronbinary+bacmetbinary+conj+reptypes+taxa+vfdbbinary+geographies+isolationsources'))
# frm
# testglm<-glm(frm,family='binomial',data=finaldftrunc)
# car::vif(testglm)

###

# ###TESTING
# 
# gamcheckfile1<-'output_adjusted/exploratory/gamcheck_model1.1.txt'
# gamcheckfile2<-'output_adjusted/exploratory/gamcheck_model1.2.txt'
# outcomeclass<-outcomeclasses[1]
# gamcheckplotfile<-gsub('%s',outcomeclass,'output_adjusted/exploratory/gamcheck_model1.1_%s.pdf')
# 
# cat('gam.check output\n',file=gamcheckfile1,sep = '\n',append = FALSE)
# writegamcheck(modellist1.1,outcomeclass,gamcheckfile1,gamcheckplotfile)
# 
# model<-modellist1.1
# #gcouts<-capture.output(gam.check(model[[outcomeclass]]))
# 
# 
# gcouts<-capture.output(gam.check(model[[outcomeclass]]))
# str(gcouts)
# cat(gcout)
# for (gcout in gcouts) {
#   if (grepl('^s\\(',gcout)==TRUE || grepl('^\\s+k',gcout)==TRUE) {
#     gcout<-gsub("*\\s+", "\t", gcout)
#     cat(gcout)
#   }
# }
# 
# library(stringr)
# str_replace_all(gcouts[14],'\s*','\t')
# 
# gsub("*\\s+", "\t", gcouts[14])
# grepl('^s\\(',gcouts[14])
# 
# grepl('^\\s+k',gcouts[12])
# 
# grepl('a','apple')
# 
# #for (gcout in gcouts) {
# #  cat(gcout,sep='\n',append = TRUE,file = gamcheckfile)
# 
# 
# 
# concurvity(model[[outcomeclass]],full=TRUE)
# concurvity(model[[outcomeclass]],full=FALSE)
# ###

