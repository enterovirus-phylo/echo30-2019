
#Example runs:

# 3Dpol run:
# snakemake run_3D

# vp1 run:
# snakemake run_vp1

wildcard_constraints:
    gene="vp1|3D"

rule all:
    input:
        auspice_tree = expand("{seg}/auspice/echo30-2019_{seg}_tree.json", seg=["3D","vp1"]),
        auspice_meta = expand("{seg}/auspice/echo30-2019_{seg}_meta.json", seg=["3D","vp1"])

rule run_3D:
    input:
        auspice_meta = "3D/auspice/echo30-2019_3D_meta.json",
        auspice_tree = "3D/auspice/echo30-2019_3D_tree.json"

rule run_vp1:
    input:
        auspice_meta = "vp1/auspice/echo30-2019_vp1_meta.json",
        auspice_tree = "vp1/auspice/echo30-2019_vp1_tree.json"

# Data files are not included as part of Github repository

rule files:
    input:
        seqs = "{gene}/data/all_sequences.fasta",
        meta = "{gene}/data/all_meta_clinical.tsv",

        dropped_strains = "{gene}/config/dropped_strains.txt",
        reference = "config/echo30_{gene}_ref.gb",
        colors = "config/colors.tsv",
        auspice_config = "{gene}/config/auspice_config.json",
        regions = "config/geo_regions.tsv"

files = rules.files.input

rule filter:
    input:
        sequences = files.seqs,
        metadata =  files.meta,
        exclude = files.dropped_strains
    output:
        sequences = "{gene}/results/filtered.fasta"
    params:
        min_date = 1960,
        min_len = 250
    shell:
        """
        # If VP1 run, exclude those that are R2 segment (don't work)
        if [ "{wildcards.gene}" == "vp1" ]; then
            echo "Filtering out R2 segments"
            exclude_where1='--exclude-where'
            exclude_where2='gene=VP1 R2'
        else
            exclude_where1="--exclude-where"
            exclude_where2=""
        fi

        augur filter --sequences {input.sequences} --metadata {input.metadata} \
            --output {output.sequences} \
            $exclude_where1 "$exclude_where2" \
            --exclude {input.exclude}  --min-date {params.min_date} \
            --min-length {params.min_len} 
        """

rule align:
    input:
        sequences = rules.filter.output.sequences,
        reference = files.reference
    output:
        alignment = "{gene}/results/aligned.fasta"
    shell:
        """
        augur align --sequences {input.sequences} --output {output.alignment} \
            --reference-sequence {input.reference} --remove-reference
        """

rule tree:
    input:
        alignment = [rules.align.output.alignment] 
    output:
        tree = "{gene}/results/raw_tree.nwk"
    params:
        model = "GTR+R10"
    shell:
        """
        augur tree --alignment {input.alignment} --output {output.tree} \
            --substitution-model {params.model} 
        """
        #--method iqtree 
        #--tree-builder-args "-ninit 100 -me 0.01"

rule refine:
    input:
        tree = rules.tree.output.tree,
        alignment = rules.align.output.alignment,
        metadata = files.meta, 
    output:
        tree = "{gene}/results/tree.nwk",
        node_data = "{gene}/results/branch_lengths.json"
    params:
        clock_filter_iqd = 5,
        clock_rate = 4E-3
    shell:
        """
        # If 3D, give the rate, as estimating doesn't work
        if [ "{wildcards.gene}" == "3D" ]; then
            echo "Setting clock rate at 4E-3"
            clock_rate='--clock-rate {params.clock_rate}'
        else
            echo "Clock rate will be estimated"
            clock_rate=""
        fi

        augur refine --tree {input.tree} --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} --output-node-data {output.node_data} \
            --timetree --date-inference marginal --coalescent opt \
            --clock-filter-iqd {params.clock_filter_iqd} \
            $clock_rate
        """
        #--clock-rate {params.clock_rate}

rule ancestral:
    input:
        tree = rules.refine.output.tree,
        alignment = rules.align.output.alignment,
    output:
        nt_data = "{gene}/results/nt_muts.json"
    params:
        inference = "joint"
    shell:
        """
        augur ancestral --tree {input.tree} --alignment {input.alignment} \
            --output {output.nt_data} --inference {params.inference} \
            --keep-ambiguous 
        """
        #--keep-overhangs

rule translate:
    input:
        tree = rules.refine.output.tree,
        node_data = rules.ancestral.output.nt_data,
        reference = files.reference
    output:
        aa_data = "{gene}/results/aa_muts.json"
    shell:
        """
        augur translate --tree {input.tree} --ancestral-sequences {input.node_data} \
            --output {output.aa_data} --reference-sequence {input.reference}
        """

rule traits:
    input:
        tree = rules.refine.output.tree,
        metadata = files.meta, 
    output:
        node_data = "{gene}/results/traits.json",
    params:
        columns = "group"
    shell:
        """
        augur traits --tree {input.tree} --metadata {input.metadata} \
            --output {output.node_data} --confidence --columns {params.columns}
        """

rule export:
    input:
        tree = rules.refine.output.tree,
        metadata = files.meta,
        branch_lengths = rules.refine.output.node_data,
        traits = rules.traits.output,
        nt_muts = rules.ancestral.output.nt_data,
        aa_muts = rules.translate.output.aa_data, 
        colors = files.colors,
        auspice_config = files.auspice_config
    output:
        auspice_meta = "{gene}/auspice/echo30-2019_{gene}_meta.json",
        auspice_tree = "{gene}/auspice/echo30-2019_{gene}_tree.json"
    shell:
        """
        augur export v1 --tree {input.tree} --metadata {input.metadata} \
            --node-data {input.branch_lengths} {input.nt_muts} {input.aa_muts} {input.traits} \
            --auspice-config {input.auspice_config} \
            --output-tree {output.auspice_tree} --output-meta {output.auspice_meta} \
            --colors {input.colors}
        """