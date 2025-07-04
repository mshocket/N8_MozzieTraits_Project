#Notes on format, ... of data 


##Adult survival 
trait in format "1/mu"      "p"         "p/days"    "prop.dead"

##fecundity data
Trait Name rather than trait.name
trait in format "EFD"  "EFOC" "EPR"  "pO"   "TFD" 
for mogi trait labeled ass EPR (eggs per raft) ovariole count 

EFOC - Calculated from Total eggs / total females from 1st oviposition cycle; translated from Chinese paper	
EPR - egg per raft (ovipositing female), some are ovariole count 
pO - proportion ovipositing	

##MDR data
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

##PDR data
trait in format "EIP"  "PDR"
extrinsic incubation?

##a data
trait in format "a" and "GCD" - not sure what any of these stand for


##pLA data
traits all "pLA" - not sure what this stands for