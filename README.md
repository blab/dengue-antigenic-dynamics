# Dengue antigenic relationships predict evolutionary dynamics

Dengue virus (DENV) exists as four genetically distinct serotypes, each of which is also antigenically distinct: the immune response to primary infection can be either cross-protective or associated with severe disease upon heterotypic secondary infection.
Each serotype is historically assumed to antigenically uniform.
Recent analyses suggest that antigenic heterogeneity may exist within each serotype, but its source, extent and impact remain unclear.
Here, we construct a phylogeny-based model to directly map antigenic change to underlying genetic divergence.
We report moderate antigenic diversity within each serotype, and identify 12 antigenically distinct clades.
We also quantify the impact of this antigenic heterogeneity on real-world DENV population dynamics.
We find that antigenic fitness mediates fluctuations in DENV clade frequencies, although this appears to be driven by coarser serotype-level antigenic differences.
These results provide a more nuanced understanding of dengue antigenic evolution, with important ramifications for vaccine design and epidemic preparedness.

# Analysis outline  
1 - [Run the titer model via augur](./augur/) (repackaged portion of the [Nextstrain](www.nextstrain.org/dengue) pipeline) to build a viral phylogeny, assign antigenic change to specific branches, and infer clade frequencies.  
2 - [Run the fitness model](./fitness_model/) to quantify population immunity over time, predict clade frequencies, and assess performance.  
3 - [Use the visualization notebooks](./figures/) to explore results and recreate all the figures from the paper.

