ML_SS=function(Z,Y,z,tree){
  
  rf_model=ranger(Y~.,data=data.frame(Y,Z),num.trees=tree,mtry=ncol(Z))
  predictions=predict(rf_model,data = data.frame(z))$predictions
  return(predictions)
}