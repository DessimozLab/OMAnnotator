
//validate input params files 
     query = file(params.query)
     dbName = params.dbName
     genome = file(params.genome)
     ecutoff = params.ecutoff
     blast_to_exonerate_py = params.blast_to_exonerate_py 

	Channel
     .fromPath(params.query)
     .splitFasta(by: 1, file:true)
     .into { query_file1_ch; query_file2_ch }

i=0; query_file1_ch.map{[i++, it]}.set{keyed_qfile1_ch}
j=0; query_file2_ch.map{[j++, it]}.set{keyed_qfile2_ch}


process makeblast {

    input:
    set val(key), file(queryFile) from keyed_qfile1_ch
    val dbName

    output:
    set val(key), file('chunk.blast') into blast_output_ch

    script:
    """
    tblastn -db $dbName -query $queryFile -out "chunk.blast" -outfmt 7 -evalue $ecutoff  
    """
}


process blast_to_exonerate {
    
    input: 
    set val(key), file(blast_output) from blast_output_ch
    file genome

    output:
    set val(key), file('chunk.fasta') into tgenome_chunks

    script:
    """
    python $blast_to_exonerate_py -i $blast_output -g $genome -o "chunk.fasta"
    """
}


process exonerate {

    input:
    set val(key), file(targetGenome), file(queryFile) from tgenome_chunks.combine(keyed_qfile2_ch, by:0)

    output:
    set val(key), file('result_hit.fasta') into results

    script:
    """
    exonerate --exhaustive --model protein2dna:bestfit --showtargetgff --showalignment no --showvulgar no --softmasktarget 1 --query $queryFile --target $targetGenome --querytype protein --targettype dna --subopt 0 --forcegtag 1 > "result_hit.fasta"

    """

}

results
  .collectFile(name: 'all_hit.fasta', storeDir : params.outdir, sort: { it[0] }) { it[1] }





