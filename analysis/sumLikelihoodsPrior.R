# Rather than performing density estimation on the mean shifts of the other 
# variants (vertical density estimation), let's just add up their likelihood 
# functions according to their weight. This accounts for the variance of the
# outcomes naturally, rather than weighting by the standard deviation after the
# fact.

library(tidyverse)
library(magrittr)


