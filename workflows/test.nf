nextflow.enable.dsl=2 // Ensure you're using DSL2

workflow {
    // Define your channels emitting tuples directly
    Channel ch_gff = Channel.from(['id:Vespula_vulgaris', '/workspace/genomeqc/work/27/344f03ea8c9b25e1eaf270c8fac617/Vespula_vulgaris.longest.gff3'])
    Channel ch_fna = Channel.from(['id:Vespula_vulgaris', '/workspace/genomeqc/work/4e/0cd3fdbba51ad8ba5e1927b199be15/GCF_905475345.1_iyVesVulg1.1_genomic.fna'])

    // Combine the channels
    Channel ch_combined = ch_gff
        .flatMap { gff ->
            ch_fna.map { fna ->
                // Check that IDs match
                assert gff[0] == fna[0], "IDs do not match!"
                [gff[0], gff[1], fna[1]] // return new tuple
            }
        }

    // Print the output to verify
    ch_combined.subscribe { tuple ->
        println "Combined Output: ${tuple}"
    }
}
