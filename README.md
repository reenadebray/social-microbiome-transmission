# social-microbiome-transmission
**Contents:**

**_Reference genomes_**: Scripts for selecting, comparing, and annotating reference microbial genomes for metagenomic strain profiling.

  • Database customization: Selection of species-representative microbial genomes based on quality and host origin for downstream mapping and genome annotations.
  
  • Gene and pseudogene annotation: Annotate gene content and identify potentially truncated or fragmented genes.
  
  • Trait prediction: Predict oxygen and temperature tolerances based on genome sequence.
  
  • Build phylogenetic tree: Mashtree construction of phylogeny from microbial reference genomes.

**_Metagenome profiling_**: Scripts for processing metagenomic data and profiling strains.

  • Pre-processing: Remove adapters, trim reads based on quality scores, and remove reads mapping to the host.
  
  • inStrain: Profile strains based on average nucleotide identity.
  
  • SynTracker: Profile strains based on microsynteny.
