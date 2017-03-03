permtest = vector('numeric', 100000)
refdat = observations %>% filter(type == 'Ref')
mutdat = observations %>% filter(type == 'Mut')
for(i in 1:100000){
  permtest[i] = mean(base::sample(mutdat$qnact, 140, replace = TRUE)) - mean(base::sample(refdat$qnact, 140, replace = TRUE))
}
hist(permtest, breaks = 100)
