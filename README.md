# Nextstrain Runs for the paper "The molecular epidemiology and evolutionary trajectory of echovirus 30 associated with its recent emergence in Europe."

This repository contains the Nextstrain augur pipeline output files which allow vizualisation of the result via the 'community' feature on Nextstrain. 

You can check out the VP1 run here: [nextstrain.org/community/enterovirus-phylo/echo30-2019/vp1](https://nextstrain.org/community/enterovirus-phylo/echo30-2019/vp1)

And the 3D _pol_ run here: [nextstrain.org/community/enterovirus-phylo/echo30-2019/3D](https://nextstrain.org/community/enterovirus-phylo/echo30-2019/3D)

And a tanglegram of the two genes [here](https://nextstrain.org/community/enterovirus-phylo/echo30-2019/3D:community/enterovirus-phylo/echo30-2019/vp1?c=group).

You can view the whole manuscript [here]() *(paper link will be added once published)*.

This analysis was created with sequences of Echovirus 30 sampled from across Europe from 2016-2018. VP1 sequences with >=250bp, and 3D _pol_ sequences that match the 540bp region in `config/echo30_3D_ref.gb` were included in a background of older and longer sequences, in two separate runs.

This was run using the snakemake pipeline in this repository. You can see the current code above, and the frozen version (the exact version used to run this result) [here](). *(Code not yet frozen.)*

