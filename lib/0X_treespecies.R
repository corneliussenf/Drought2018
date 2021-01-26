
library(tidyverse)

dat <- read_csv("data/treespecies/EUForestspecies.csv")

head(dat)

ggplot(dat %>% filter(`SPECIES NAME` == "Picea abies")) +
  geom_point(aes(x = X, y = Y))
