library(reshape2)
library(ggplot2)
library(scales)
library(readr)
library(tidyverse)

setwd("C:/Users/user/source/repos/DoL/DoL")


data1 <- read.table("outf1-5_1.txt", header = T)

#population size dynamics
ggplot(data1, aes(t, popsize))+ geom_line() 

#average body size dynamics
ggplot(data1, aes(t, avgbodysize))+ geom_line()

#standard deviation need
ggplot(data1, aes(t, sdneed))+ geom_line()

#Tasks performed by group

data1 <- gather(data1, eggcare, digging, defense, key = "behaviour", value = "prop")

ggplot(data1, aes(t, prop, colour = behaviour))+ geom_line()+
  geom_line(aes(t, avgbodysize/5), colour = "black")+xlim(99000, 100000)

ggplot(data1, aes(t, prop, colour = behaviour))+ geom_line()+
  geom_line(aes(t, popsize/5), colour = "black")+xlim(99000, 100000)

  

#Individual-level data
data2 <- read.table("aug1-2.txt", header = T)

data2 <- filter(data2, !ID %in% data2$ID[which(data2$t==99001)] )

#To select data for 1000 individuals
data2 <- filter(data2, ID %in% sample(unique(data2$ID), 500))

data2$suit_egg <- 1/(1 + exp((data2$size - 15)/4))
data2$suit_digg <- 3/((1 + exp((25-data2$size)/4)) * (1 + exp((data2$size-25)/4)))
data2$suit_def <- 1/(1 + exp((25-data2$size)/4))

data2$labour_egg <- data2$task_egg * (data2$suit_egg + data2$exp_egg)
data2$labour_digg <- data2$task_digg * (data2$suit_digg + data2$exp_digg)
data2$labour_def <- data2$task_def * (data2$suit_def + data2$exp_def)

#To just select 10000 different samples
#data2 <- data2 %>% slice_sample(n = 10000)

#To select data for 1000 individuals
#data2 <- filter(data2, ID %in% sample(unique(data2$ID), 1000))


lifehistory <- data2 %>% 
  group_by(ID) %>% 
  summarize(
    maxage = max(t)-min(t),
    birth = min(t),
    maxsize = max(size)
)

data2 <- data2 %>% left_join(select(lifehistory, -maxsize, -maxage), by = "ID")
data2$age <- data2$t - data2$birth

#maximum age
ggplot(lifehistory, aes(maxage))+geom_density(bounds = c(0, Inf))
#maximum size
ggplot(lifehistory, aes(maxsize))+geom_density(bounds = c(1, Inf))

ggplot(data2, aes(size))+geom_density(bounds = c(1, Inf))
ggplot(data2, aes(age))+geom_density(bounds = c(0, Inf))

ggplot(data2, aes(age, size))+geom_point()+xlim(0, 100)

#three different tasks over body size
ggplot(data2, aes(size, task_egg))+geom_point()
ggplot(data2, aes(size, task_digg))+geom_point()
ggplot(data2, aes(size, task_def))+geom_point()

ggplot(data2, aes(size, exp_egg))+geom_point()
ggplot(data2, aes(size, exp_digg))+geom_point()
ggplot(data2, aes(size, exp_def))+geom_point()


ggplot(data2, aes(t, task_egg, colour = size))+geom_point()+xlim(99000, 100000)
ggplot(data2, aes(t, exp_egg, colour = size))+geom_point()+xlim(99000, 100000)
ggplot(data2, aes(t, task_egg, group = ID))+geom_line()+xlim(95000, 96000)

ggplot(data2, aes(t, task_egg, colour = size, group = ID))+geom_point()+geom_line()+xlim(99000, 100000)


ggplot(data2, aes(t, suit_egg, colour = size))+geom_point()+xlim(99000, 100000)
#labour = task * (suitability + experience)

groupbehaviour <- data2 %>% 
  group_by(t) %>% 
  summarize(
    m_task_egg = mean(task_egg),
    m_exp_egg = mean(exp_egg),
    m_suit_egg = mean(suit_egg),
    m_labour_egg = mean(labour_egg)
  )

groupmeasures <- gather(groupbehaviour, m_task_egg, m_exp_egg, m_suit_egg, m_labour_egg, key = "behaviour", value = "val")

ggplot(groupmeasures, aes(t, val, colour = behaviour))+geom_line()+xlim(99000, 100000)







library(plotly)

plot_ly(
  data2, a = ~task_egg, b = ~task_digg, c = ~task_def,
  color = ~size,
  type = "scatterternary")

plot_ly(data2, x = ~task_egg, y = ~task_digg, z = ~task_def, color = ~size)

ggplot(filter(data2, size > 40))+geom_point(aes(task_def, task_digg, color = size))
ggplot(data2)+geom_point(aes(task_def, task_digg, color = size))

ggplot(data2, aes(task_def, task_digg, color = size, group = ID)) + geom_point() + geom_line()




ggplot(data2)+geom_point(aes(exp_def, exp_digg, color = size))

ggplot(filter(data2, size > 40))+geom_point(aes(exp_def, exp_digg, color = size))


ggplot(filter(data2, size > 40))+geom_point(aes(exp_def, exp_digg, color = size))


unique(filter(data2, size > 45)$ID)
# track single individuals

ggplot(filter(data2, ID == 18181      ))+geom_point(aes(exp_def, exp_digg, color = size))
ggplot(filter(data2, ID == c(595278, 596359 )           ))+geom_point(aes(task_def, task_digg, color = size))
ggplot(filter(data2, ID == 18181      ))+geom_point(aes(t, task_digg, color = size))
ggplot(filter(data2, ID == 18181      ))+geom_point(aes(t, exp_digg, color = size))

plot_ly(filter(data2, ID == 545691     ), x = ~task_egg, y = ~task_digg, z = ~task_def, color = ~size, size = 0.1)


ggplot(data2, aes(task_def, task_digg, color = size, group = ID)) + geom_point() + geom_line()


ggplot(data2, aes(size, colour = as.factor(current_task)))+geom_density()

ggplot(data2, aes(time, tas))



ggplot(arrange(filter(data2, ID == c(595278)           ), t), aes(task_def, task_digg, color = size, group = ID))+geom_point()+geom_path()

ggplot(arrange(data2, t), aes(task_def, task_digg, color = size, group = ID)) + geom_point() + geom_path()

p1 <- ggplot(arrange(data2, t), aes(task_def, task_digg, color = size, group = ID)) + geom_point() + geom_path()
p1
ggsave("DoL_smallpop2.png", p1)
