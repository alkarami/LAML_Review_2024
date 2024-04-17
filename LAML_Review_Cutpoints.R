## Cutpoint determination

library(survminer)
library(ggplot2)
library(survival)

# Upload consolidated clinical and expression dataset
clinical_BEAT <- as.data.frame(read.csv('clinical_BEAT.csv', 
                header = T, row.names = 1))

# Cutpoint determination, make survival plots as well
for(x in colnames(clinical_BEAT)[101:600]){
  cutx <- surv_cutpoint(
    clinical_BEAT,
    time = "overallSurvival",
    event = "vital",
    variables = c(x)
  )
  # Get the optimal cutpoint
  cutval <- summary(cutx)$cutpoint
  # Create vector 
  cutvector <- c()
  for(y in clinical_BEAT[,x]){
    if(y < cutval){cutvector <- append(cutvector,'Low')}else{cutvector <-
      append(cutvector,'High')}
  }
  cutvector <- factor(cutvector, levels = c('Low','High'))
  # Add the classification to the metadata
  clinical_BEAT <- cbind(clinical_BEAT,cutvector)
  genex <- paste(x,"_Expression", sep = "")
  colnames(clinical_BEAT)[ncol(clinical_BEAT)] <- 
    genex
  
  # Survival curves
  survformula <- as.formula(paste("Surv(overallSurvival, vital)",
                                  "~",genex, collapse = ""))
  surall_fit <- surv_fit(survformula,
                         data=clinical_BEAT)
  p <- ggsurvplot(surall_fit, xlab = 'Time in Days', ylab = 'Survival', 
                  pval = T, palette = c('blue','red')) 
  ggsave(paste('SurvivalCurves/',genex,'.png'), 
         plot = p$plot)
}
