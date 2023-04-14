# Get gene expression data and integrate this into networks.

# TODO ALL:
#  1. Get gene expression data from GEO or from the TCGA portal (need pre-processed, normalized data).
#  2. Read this data in using pandas?
#  3. Use pybedtools to intersect the gene expression data with nodes in the network?
#  4. Use the gene expression data average to weight the edges in the network? Or color them according to expression?
#  5. Make gene expression plots for each cell line, per chromosome, per resolution.
#  6. Make correlation plots between gene expression and network metrics for each cell line?

# TODO classes:
#  1. Make a class that takes a gene expression file, reads it in and stores it.
#  2. Make a class that takes gene expression object (dataframe) and intersects it with a given network (graph dict) using pybedtools.
#          This class should return the network with attributes linked to average gene expression? Some way to integrate the expression data for nodes overlapping genes.
#  3. Make a class that takes a graph with gene expression data and correlates it with relevant network metrics.

# TODO plotting:
#  1. Make a plot that shows nodes colored by gene expression (what metric?), takes a graph dict as input with gene expression data integrated as input.



