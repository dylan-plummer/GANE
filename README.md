GANE
====
Protein Complexes Identification Based on GO Attributed Network Embedding

Installation
------------
####Requirements
        numpy

Input and Output
------------

    Input: PPI datasets such as "DIP.txt"and GO information dataset such as "go_slim_mapping.tab.txt"
    Output: For example, in "final_dip_attr_output", each line denotes a protein complex 

Run
------------
    step1. 1_Get_go_information.py: getting go information for each PPI network from go_slim_mapping.tab.txt.
    step2. 2_Create_topological_biological_network.py: create two matrixes for each PPI network.
    step3. 3_Node_embedding.m: generating vector representations for each protein based on two matrixes for a PPI network.
    step4. 4_Update_linking_weight.py: updating the weight of each edge in a PPI network based on calculating cosin similarity with vector representations.
    step5. 5_cluster_core_attachment.py: generating complexes based on the structure of core and attachments.
    step6. 6_Compare_performance.py: compare the preformance of GANE with other classic methods based on different evaluation metrices.
