# scRNA-Seq-to-correct-for-batch-effects

Batch effect happens when the variation in sample groups is caused by technical arrangement rather than biological factors, leading to false conclusions. From here, the need for batch effect correction arises. Batch effect correction attempts to remove technical variance when combining cells of different batches or from different studies. 

scRNA-seq normalization also aims to eliminate technical noise or bias so that observed variance in gene expression variance primarily reflects true biological variance.

# The types of technical noise that each step attempts to correct are different: 

    Normalization targets variance from sequencing (library preparation, high dropout event, amplification bias caused by gene length  GC content, etc.)
    Batch effect correction targets variance from experimental designs and handling (dfferent sequencing platforms, timing, reagents, laboratories, etc.)

    
