#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)

rd = function(v,o) { 100*(v-o)/o }
code = function(t) { mutate(t,set=recode(set,taillard="ta",carlier="car"),incompatibilities=sprintf("%02d",10*incompatibilities)) }

average = function(t,group=1) {
    return ( t
            %>% merge(bk)
            %>% mutate(r.cplex=rd(makespan.cplex,makespan.bkv),r.sspr=rd(makespan,makespan.bkv),gr=(instance-1) %/% group,time=time/100,pub=(makespan==makespan.bkv))
            %>% group_by(set,gr,variation,incompatibilities)
            %>% summarize(across(c(n,m,makespan.bkv,r.cplex,gap,time,r.sspr,pub),~mean(.x,na.rm=T)))
            %>% mutate(pub=100*pub)
            %>% rename(ub=makespan.bkv)
            %>% unite("Instance",c(set:incompatibilities),sep="")
            )
}

## read data
df=read.csv("run.csv") %>% code()
bk=read.csv("bkv.csv") %>% code()

## table 2
t2 = df %>% filter(set=="car" & variation=="i") %>% average()
t2 %>% data.frame()

## table 3
t3 = df %>% filter(set=="car" & variation=="I") %>% average()
t3 %>% data.frame()

## table 4
t4 = df %>% filter(set=="ta" & variation=="i") %>% average(group=10)
t4 %>% data.frame()

## table 5
t5 = df %>% filter(set=="ta" & variation=="I") %>% average(group=10)
t5 %>% data.frame()



