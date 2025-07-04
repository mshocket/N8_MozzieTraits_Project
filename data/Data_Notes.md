# Notes on format, ... of data 


## Adult survival 
trait in format "1/mu"      "p"         "p/days"    "prop.dead"

## fecundity data
Trait Name rather than trait.name 

trait in format "EFD"  "EFOC" "EPR"  "pO"   "TFD" 

for mogi trait labeled ass EPR (eggs per raft) ovariole count 

EFOC - eggs per female per oviposition cycle 

EPR - egg per raft (ovipositing female), some are ovariole count

pO - proportion ovipositing	

tfd - total fecundity 

efd eggs per female per day 

## MDR data
mosquito development rate

trait in format "MDR"  "1/MDR". Convert to 1/MDR (days)

```
#take reciprical so all in days?
for(i in 1:nrow(MDR_Data)){
  
  if(MDR_Data$trait.name[i] == "MDR"){
    MDR_Data$trait[i] <- 1/MDR_Data$trait[i]
  }
}
```

## PDR data
pathogen devolopment rate
trait in format "EIP"  "PDR"

pathogen development rate

extrinsic incubation period?
inverse of each other 

## a data
trait in format "a" and "GCD" 

biting rate

1/ gonatrophic cycle duration

a is rate 

## pLA data
traits all "pLA" - 
probability of surviving from larvae to adulthood