# Proposal for Package Construction

## Purpose
My previous effort to design an HGT detection pipeline (kvasir) used BLAST to identify very high identity genes in species that would otherwise not be expected to share similar genes. This approach was fairly easy to implement, but suffers from a number of limitations. In particular, the need to set very high thresholds for identity to exclude false positives, and the need to use the same thresholds for all pairs of species rather than calculating them dynamically. This proposal lays out a strategy for a new software package, written in julia, to overcome these deficits.

## Strategy
At its most basic, this strategy is not particularly novel. Given two organisms that share a common ancestor, if only vertical inheritance is at play, their genomes should diverge at some rate. This rate will vary across the genome, since some genes are under stronger selection or different selection, and large changes such as insertions, deletions and transposon hopping may drastically alter local genomic regions all at once in ways not predicted by a constant rate of divergence. If these organisms share genes through HGT, the regions of the genome that these genes occupy should have higher identity than the rest of the genome. The goal of this package is to identify those regions that have higher identity than would be predicted based on vertical inheritance.

There are both theoretical and practical considerations in designing this package. Theoretical considerations include setting thresholds for nucleotide identity given the varying rates of divergence for different types of genes (variation is not likely to be normal across the genome - how many standard deviations from the mean is “significant”). Practical considerations include things like locally storing gene and phylogenetic information (database?), what to compare (protein coding genes, whole genomes?) etc.

### Here’s one strategy:
Determine average nucleotide identity (ANI) between species A and species B.
Set thresholds for deviation from ANI that indicate HGT.
Find genomic regions that have identity above that threshold.

**Pros**
- Relatively straightforward
- Does not require gene annotation or prediction
**Cons**
- Ignores a huge amount of phylogenetic information.

### Here’s a different strategy assuming more than two species in the analysis:
Determine phylogeny of all species based on ANI, 16S and other factors.
Divide each genome into reasonable units (could be genes, large kmers etc)
Determine phylogeny of each unit.
Units with phylogenies that differ substantially from the overall phylogeny are likely HGT

**Pros**
- More comprehensive
- As more species are added, accuracy of predictions increases (may be able to set more dynamic thresholds to minimize both type 1 and type 2 errors)
**Cons**
- Results are dependent on genomes included in analysis.
- Adding new data may require re-running entire analysis.
- Computationally intensive.

## Other Considerations
Knowing the plan is also not sufficient, decisions need to be made regarding how to generate alignment (HMM? BLAST?), how to store/parse data (using local database?) etc.
