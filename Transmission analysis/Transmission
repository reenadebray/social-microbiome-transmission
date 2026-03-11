# Custom R script to identify transmission events from strain comparison table

## Run transmission function
### sharing = comparisons output file from inStrain compare (can be adapted to column names in SynTracker output)
### sample_list = data frame with sample IDs that include individual name and date
### popANI = user-defined minimum threshold for genome similarity (we used 0.99999)
### percent_compared = user-defined minimum threshold for overlap in coverage between samples (we used 0.25)
### window = option to relax requirement that first recipient sample was collected AFTER donor sample, by allowing recipient samples collected up to n days before donor sample (we used n=0)

transmission<-function(sharing,sample_list,popANI,percent_compared,window){
  transmission_events<-data.frame(matrix(nrow=0,ncol=4))
  
# start with all strain sharing events (≥0.99999, ≥0.25) in population
  strain_sharing<-sharing[sharing$popANI>=popANI & sharing$percent_compared>=percent_compared & !is.na(sharing$popANI),]
  first_samples<-paste(sample_list[sample_list$temporal.sample==1,"library_ID"],".sorted.bam",sep="")
  # exclude comparisons that are both individuals' first sample
  strain_sharing_excl_first<-strain_sharing[!(strain_sharing$name1%in%first_samples & strain_sharing$name2%in%first_samples),]
  
  for (genome in unique(strain_sharing_excl_first[substr(strain_sharing_excl_first$name1,1,3)!=substr(strain_sharing_excl_first$name2,1,3),"genome"])){ # for each species-dyad combo
      for (dyad in unique(strain_sharing_excl_first[strain_sharing_excl_first$genome==genome & substr(strain_sharing_excl_first$name1,1,3)!=substr(strain_sharing_excl_first$name2,1,3),"dyad"])){
        sname1<-substr(dyad,1,3)
        sname2<-substr(dyad,5,7)
        detected_strains<-unique(c(strain_sharing[strain_sharing$genome==genome & strain_sharing$sname1%in%c(sname1,sname2) & strain_sharing$sname2%in%c(sname1,sname2),"name1"],strain_sharing[strain_sharing$genome==genome & strain_sharing$sname1%in%c(sname1,sname2) & strain_sharing$sname2%in%c(sname1,sname2),"name2"])) # pull out ALL comparisons involving this genome and these individuals, including first time-points and self-comparisons
        
        # when was this strain first sampled in each member of the dyad?
        sname1_samples<-detected_strains[substr(detected_strains,1,3)==sname1]
        sname1_dates<-as.Date(substr(sname1_samples,5,11),format="%d%b%y")
        sname1_first_date<-sort(sname1_dates)[1]
        
        sname2_samples<-detected_strains[substr(detected_strains,1,3)==sname2]
        sname2_dates<-as.Date(substr(sname2_samples,5,11),format="%d%b%y")
        sname2_first_date<-sort(sname2_dates)[1]
      
      if (sname1_first_date<sname2_first_date){ #sname1 had the strain first
        # was sname2 sampled after (or very shortly before) the strain was detected in sname1, and BEFORE the strain was first detected in sname2?
        sname2_samples<-sample_list[sample_list$sname==sname2,"date"]
        intermediate_sname2_samples<-sname2_samples[sname2_samples>=(sname1_first_date-window) & sname2_samples<sname2_first_date]
        if (length(intermediate_sname2_samples)>0){transmission_events<-rbind(transmission_events,c(genome,sname1,sname2,as.character(sname2_first_date)))}
      }
    
      else if (sname2_first_date<sname1_first_date){ #sname2 had the strain first
        # was sname1 sampled after (or very shortly before) the strain was detected in sname2, and BEFORE the strain was first detected in sname1?
        sname1_samples<-sample_list[sample_list$sname==sname1,"date"]
        intermediate_sname1_samples<-sname1_samples[sname1_samples>=(sname2_first_date-window) & sname1_samples<sname1_first_date]
        if (length(intermediate_sname1_samples)>0){transmission_events<-rbind(transmission_events,c(genome,sname2,sname1,as.character(sname1_first_date)))}
      }
    }# if sname1_first_date == sname2_first_date, can't resolve transmission directionality.
  }
  colnames(transmission_events)<-c("genome","donor","recipient","recipient_date")
  return(transmission_events)
}

## Filter for transmission events with only one possible donor

for (i in 1:nrow(transmission_events)){
 genome<-transmission_events[i,"genome"]
  recipient<-transmission_events[i,"recipient"]
  recipient_group<-transmission_events[i,"recipient_group"]
 if (nrow(transmission_events[transmission_events$genome==genome & transmission_events$recipient==recipient,])==1){transmission_events[i,"unique_within_pop"]<-"Y"}
  else {transmission_events[i,"unique_within_pop"]<-"N"}
}
