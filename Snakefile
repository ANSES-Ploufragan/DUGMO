import os
configfile: "configDUGMO.json"

rule target:
  input:
     config["outputFolder"] + "/HostGenome/CDS_GMO_host.fa",
     config["outputFolder"] + "/DatabankGMO/GMO_Databank_clean.fa",
     expand(config["outputFolder"] + "/Rmes/freq_{CDStype}_l{RtailleMot}_m{RMarkovOrder}_{RpourcentSurRep}pSurRep.txt", zip, 
		CDStype=["CDS_GMO_host", "CDS_GMO_host-pos3codon", "CDS_GMO_host"], 
		RMarkovOrder=[config["RmesParams"]["MarkovOrderProp"], config["RmesParams"]["MarkovOrderFreq"], 1], 
		RtailleMot=[config["RmesParams"]["tailleMotProp"], config["RmesParams"]["tailleMotFreq"], 3], 
		RpourcentSurRep=[config["FrequenceParams"]["pourcentSurRepProp"], config["FrequenceParams"]["pourcentSurRepFreq"], config["FrequenceParams"]["pourcentSurRepProp"]]),    
     expand("{sample}_mixl4m2Propv4l9m7Freq_calculdistBC-CodonUsage.csv", zip, 
        	sample=[config["outputFolder"] + "/DatabankGMO/GMO_Databank_clean", config["outputFolder"] + "/HostGenome/CDS_GMO_host", config["outputFolder"] + "/GMO/CDS_GMO_Inserts"]),
     expand("{inp}_mixl4m2Propv4l9m7Freq_calculdistBC-CodonUsage_notRNA.csv", 
		inp=[config["outputFolder"] + "/HostGenome/CDS_GMO_host", config["outputFolder"] + "/DatabankGMO/GMO_Databank_clean"]),
     config["outputFolder"] + "/MachineLearning/Pred-RFPositifs-all.txt"

############################### CLEANING PIPELINE

rule trimmomatic_GMO:
  input:
    R1=config["Reads"]["R1"],
    R2=config["Reads"]["R2"]
  output:
    R1trim = config["outputFolder"] + "/GMO/GMO_all_R1_001.soft_clean.trimmomatic_trimmed.fastq.gz",
    R1trim_unpaired = temp(config["outputFolder"] + "/GMO/GMO_all_R1_001.soft_clean.trimmomatic_trimmed_unpaired.fastq.gz"),
    R2trim = config["outputFolder"] + "/GMO/GMO_all_R2_001.soft_clean.trimmomatic_trimmed.fastq.gz",
    R2trim_unpaired = temp(config["outputFolder"] + "/GMO/GMO_all_R2_001.soft_clean.trimmomatic_trimmed_unpaired.fastq.gz")
  conda:
    "envs/Shovill.yaml"
  params:
    logfile=config["outputFolder"] + "/GMO/logFile_trimmomatic.log",
    oligofile=config["oligo"],
    threads=config["threads"]
  shell:
    "trimmomatic PE -threads {params.threads} -trimlog {params.logfile} {input[0]} {input[1]} {output.R1trim} {output.R1trim_unpaired} {output.R2trim} {output.R2trim_unpaired} ILLUMINACLIP:{params.oligofile}:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36"

rule Alignment_Bowtie2_indexGenomeRef_GMO:
  input:
    Gref=config["Strain"]["RefGenome"],
    R1=rules.trimmomatic_GMO.output.R1trim,
    R2=rules.trimmomatic_GMO.output.R2trim
  output:
    BW1bamfile=config["outputFolder"] + "/Cleaning_pipeline/Alignment_Bowtie2_GenomeRefGMO/Bowtie2_Alignment_GenomeRefGMO.bam",
    unmappedfq1=config["outputFolder"] + "/Cleaning_pipeline/Alignment_Bowtie2_GenomeRefGMO/unalignedreads_GenomeRefGMO.fq.1.gz",
    unmappedfq2=config["outputFolder"] + "/Cleaning_pipeline/Alignment_Bowtie2_GenomeRefGMO/unalignedreads_GenomeRefGMO.fq.2.gz"
  conda:
    "envs/Aligner.yaml"
  params:
    index_Gref=config["outputFolder"] + "/Cleaning_pipeline/Alignment_Bowtie2_GenomeRefGMO/RefGenome_indexBowtie2",
    outp=config["outputFolder"] + "/Cleaning_pipeline/Alignment_Bowtie2_GenomeRefGMO/unalignedreads_GenomeRefGMO.fq.gz",
    threads=config["threads"]
  shell: """
    bowtie2-build --threads {params.threads} -f {input.Gref} {params.index_Gref}
    bowtie2 --un-conc-gz {params.outp} -x {params.index_Gref} -1 {input.R1} -2 {input.R2} -p {params.threads} | samtools sort -O bam -o {output.BW1bamfile}
    """ 
  
rule shovill:
  input:
    R1=config["Reads"]["R1"],
    R2=config["Reads"]["R2"]
  output:
    scaffolds=config["outputFolder"] + "/Cleaning_pipeline/shovill/contigs.fa",
    scaffoldsP=config["outputFolder"] + "/Cleaning_pipeline/shovill_plasmid/contigs.fa"
  conda:
    "envs/Shovill.yaml"
  params:
    folder1=config["outputFolder"] + "/Cleaning_pipeline/shovill",
    folder2=config["outputFolder"] + "/Cleaning_pipeline/shovill_plasmid"
  shell: """
    shovill --outdir {params.folder1} --R1 {input.R1} --R2 {input.R2} --trim --mincov 3.00 --force
    shovill --outdir {params.folder2} --R1 {input.R1} --R2 {input.R2} --trim --mincov 3.00 --opts '--plasmid' --force
    """
   
rule Rename_Assemblage: 
  input:
    rules.shovill.output.scaffolds,
    rules.shovill.output.scaffoldsP
  output:   
    Rename=config["outputFolder"] + "/Cleaning_pipeline/shovill/contigs_rename.fasta",
    RenameP=config["outputFolder"] + "/Cleaning_pipeline/shovill_plasmid/contigsPlasmid_rename.fasta"
  shell: """
    awk '/^>/{{print ">NODE_" ++i; next}}{{print}}' {input[0]} > {output[0]}
    awk '/^>/{{print ">NODE_" ++i; next}}{{print}}' {input[1]} > {output[1]}   
    """
	
rule Annotation_Prokka: 
  input:
    rules.Rename_Assemblage.output.Rename,
    rules.Rename_Assemblage.output.RenameP
  output:
    prokkaFFN=config["outputFolder"] + "/Cleaning_pipeline/Annotation_Prokka_AssemblageClean/prokkaAssemblageGMOmincov3.ffn",
    prokkaFFNP=config["outputFolder"] + "/Cleaning_pipeline/Annotation_Prokka_AssemblageClean_plasmid/prokkaAssemblagePlasmidmincov3.ffn",
    AnnotFusion=config["outputFolder"] + "/Cleaning_pipeline/Annotation_Prokka_AssemblageClean/prokkaAssemblageGMO_Plasmid.fa"
  conda:
    "envs/Prokka.yaml"
  params:
    folder=config["outputFolder"] + "/Cleaning_pipeline/Annotation_Prokka_AssemblageClean/",
    folderP=config["outputFolder"] + "/Cleaning_pipeline/Annotation_Prokka_AssemblageClean_plasmid/"
  shell: """
    prokka {input[0]} --outdir {params[0]} --prefix prokkaAssemblageGMOmincov3 --force
    prokka {input[1]} --outdir {params[1]} --prefix prokkaAssemblagePlasmidmincov3 --force 
    cat {output[0]} {output[1]} > {output[2]}
    """     

rule Concat_ProkkaAnnot_plasmid:
  input:
    rules.Annotation_Prokka.output.AnnotFusion,
  output:
    AnnotFusionU=config["outputFolder"] + "/Cleaning_pipeline/Annotation_Prokka_AssemblageClean/prokkaAssemblageGMO_Plasmid_uniq.fa"
  script:
    "scripts/removeDuplicateSeq.py"
     
rule Alignment_BWA_indexAnnotationProkka_unalignedGenomeRefGMO:
  input:
    Gref=rules.Concat_ProkkaAnnot_plasmid.output.AnnotFusionU,
    R1=config["outputFolder"] + "/Cleaning_pipeline/Alignment_Bowtie2_GenomeRefGMO/unalignedreads_GenomeRefGMO.fq.1.gz", #rules.Alignment_Bowtie2_indexGenomeRef_GMO.output.unmappedfq1,
    R2=config["outputFolder"] + "/Cleaning_pipeline/Alignment_Bowtie2_GenomeRefGMO/unalignedreads_GenomeRefGMO.fq.2.gz" #rules.Alignment_Bowtie2_indexGenomeRef_GMO.output.unmappedfq2
  output:
    BWA2bamfile=config["outputFolder"] + "/Cleaning_pipeline/Alignment_BWA_AnnotationGMO_unalignedGenomeRefGMOBWA/BWA_Alignment_AnnotationGMO_unalignedGenomeRefGMOBWA.bam",
    sortedBamBWA2=config["outputFolder"] + "/Cleaning_pipeline/Alignment_BWA_AnnotationGMO_unalignedGenomeRefGMOBWA/BWA_Alignment_AnnotationGMO_unalignedGenomeRefGMOBWA_sorted.bam"
  conda:
    "envs/Aligner.yaml"
  params:
    index_Gref=config["outputFolder"] + "/Cleaning_pipeline/Alignment_BWA_AnnotationGMO_unalignedGenomeRefGMOBWA/AnnotationGMOBSubtilis_reference",
    threads=config["threads"]
  shell: """
    bwa index {input.Gref} -p {params.index_Gref}
    bwa mem -t {params.threads} {params.index_Gref} {input.R1} {input.R2} | samtools view -h -Sb - > {output.BWA2bamfile}
    samtools sort -T out.prefix.bam -o {output.sortedBamBWA2} {output.BWA2bamfile}
    """

#Retrieving headers on the alignment consensus from the previous BWA alignment for aligned reads
rule Consensus_Alignment:
  input:
    rules.Alignment_BWA_indexAnnotationProkka_unalignedGenomeRefGMO.output.sortedBamBWA2,
    rules.Concat_ProkkaAnnot_plasmid.output.AnnotFusionU
  output:
    Consensus_fq=temp(config["outputFolder"] + "data/Cleaning_pipeline/Alignment_BWA_AnnotationGMO_unalignedGenomeRefGMOBWA/final_permissive.fastq"),
    Consensus_fa=config["outputFolder"] + "/Cleaning_pipeline/Alignment_BWA_AnnotationGMO_unalignedGenomeRefGMOBWA/final_permissive.fasta",
    ET1=temp(config["outputFolder"] + "/Cleaning_pipeline/Alignment_BWA_AnnotationGMO_unalignedGenomeRefGMOBWA/consensus_Alignment_entete1.txt"),
    ET=config["outputFolder"] + "/Cleaning_pipeline/Alignment_BWA_AnnotationGMO_unalignedGenomeRefGMOBWA/consensus_Alignment_entete.txt"
  conda:
    "envs/Toolsbioinfo.yaml"
  shell: """
    samtools mpileup -P ILLUMINA -p -t DP,SP -u -s -d 8000 -uf {input[1]} {input[0]} | bcftools call -c --ploidy Y | vcfutils.pl vcf2fq > {output[0]}
    seqtk seq -A {output[0]} > {output[1]}
    grep '>[A-Z]' {output[1]} > {output[2]}
    tr -d '>' < {output[2]} > {output[3]}
    """

#Recovery of CDS corresponding to previous headings -> present in the header file=potential GMO inserts/not present in the header file=CDS host GMO
rule CDS_host_and_CDS_GMOInserts_beforeBlastPangenome:
  input:
    rules.Concat_ProkkaAnnot_plasmid.output.AnnotFusionU,
    rules.Consensus_Alignment.output.ET
  output:
    CDS_host=config["outputFolder"] + "/Cleaning_pipeline/Alignment_BWA_AnnotationGMO_unalignedGenomeRefGMOBWA/CDS_GMO_host_beforeBlastPangenome.fa",
    CDS_GMOInserts=config["outputFolder"] + "/Cleaning_pipeline/Alignment_BWA_AnnotationGMO_unalignedGenomeRefGMOBWA/CDS_GMO_potentialInserts_beforeBlastPangenome.fa"
  conda:
    "envs/Toolsbioinfo.yaml"
  shell: """
    filterbyname.sh in={input[0]} out={output[0]} names={input[1]} include=f
    filterbyname.sh in={input[0]} out={output[1]} names={input[1]} include=t
    """

############################### BLASTN PANGENOME CDS AND PANGENOME WHOLE GENOME

#First Blast
rule Align_PangenomeCDSdb_GMO:
  input:
    config["Strain"]["PangenomeCDS"],
    rules.CDS_host_and_CDS_GMOInserts_beforeBlastPangenome.output.CDS_GMOInserts
  output:
    BlastdbRes=config["outputFolder"] + "/Cleaning_pipeline/Blast_Pangenome_CDS/CDS_GMO_potentialinserts_VS_PangenomeCDSdb_percIden95_gapopen3_qcovhsp98.tsv",
    BlastdbR=config["outputFolder"] + "/Cleaning_pipeline/Blast_Pangenome_CDS/CDS_GMO_potentialinserts_VS_PangenomeCDSdb_cut.txt"
  conda:
    "envs/Aligner.yaml"
  params:
    name=config["Strain"]["PangenomeCDS"] + ".db",
    evalue=0.00001,
    percidentity=95,
    threads=config["threads"]
  shell: """
    makeblastdb -in {input[0]} -out {params.name} -dbtype nucl
    blastn -query {input[1]} -db {params.name} -evalue {params.evalue} -perc_identity {params.percidentity} -qcov_hsp_perc 98 -num_threads {params.threads} -gapopen 3 -gapextend 1 -outfmt '6 qseqid sseqid pident length sstart send evalue sstrand qlen slen qcovs gaps qcovhsp' > {output[0]} 
    cut -d$'\t' -f 1,2,4,9 {output[0]} > {output[1]}
    """

rule Clean_Blast_PangenomeCDS:
  input:
    rules.Align_PangenomeCDSdb_GMO.output.BlastdbR
  output:
    CB=config["outputFolder"] + "/Cleaning_pipeline/Blast_Pangenome_CDS/CDS_GMO_potentialinserts_VS_PangenomeCDSdb_Clean.txt"
  script:
    "scripts/Blast_AlignLength.py"
    
rule PangenomeCDSdb_CDS_host_and_CDS_GMOInserts:
  input:
    rules.Clean_Blast_PangenomeCDS.output.CB,
    rules.CDS_host_and_CDS_GMOInserts_beforeBlastPangenome.output.CDS_GMOInserts
  output:
    Entete1=config["outputFolder"] + "/Cleaning_pipeline/Blast_Pangenome_CDS/CDS_GMO_potentialinserts_VS_PangenomeCDSdb_percIden95_gapopen3_qcovhsp98_EnteteUniq.txt", #temp?
    CDS_GMOpot1=config["outputFolder"] + "/Cleaning_pipeline/Blast_Pangenome_CDS/CDS_GMO_potentialinserts_afterBlastPangenomeCDS.fa",
    CDS_host1=config["outputFolder"] + "/Cleaning_pipeline/Blast_Pangenome_CDS/CDS_GMO_host_afterBlastPangenomeCDS.fa"
  conda:
    "envs/Toolsbioinfo.yaml"
  shell: """
    cut -f 1 {input[0]} | uniq > {output[0]}
    filterbyname.sh in={input[1]} out={output[1]} names={output[0]} include=f
    filterbyname.sh in={input[1]} out={output[2]} names={output[0]} include=t
    """

#Second Blast
rule Align_PangenomeWholedb_GMO:
  input:
    config["Strain"]["PangenomeWhole"],
    rules.PangenomeCDSdb_CDS_host_and_CDS_GMOInserts.output.CDS_GMOpot1
  output:
    BlastdbRes2=config["outputFolder"] + "/Cleaning_pipeline/Blast_Pangenome_wholegenome/CDS_GMO_potentialinserts_VS_Pangenomewholedb_percIden95_gapopen3_qcovhsp98.tsv",
    BlastdbR2=config["outputFolder"] + "/Cleaning_pipeline/Blast_Pangenome_wholegenome/CDS_GMO_potentialinserts_VS_Pangenomewholedb_cut.txt"
  conda:
    "envs/Aligner.yaml"
  params:
    name=config["Strain"]["PangenomeWhole"] + ".db",
    evalue=0.00001,
    percidentity=95,
    threads=config["threads"]
  shell: """
    makeblastdb -in {input[0]} -out {params.name} -dbtype nucl
    blastn -query {input[1]} -db {params.name} -evalue {params.evalue} -perc_identity {params.percidentity} -qcov_hsp_perc 98 -num_threads {params.threads} -gapopen 3 -gapextend 1 -outfmt '6 qseqid sseqid pident length sstart send evalue sstrand qlen slen qcovs gaps qcovhsp' > {output[0]} 
    cut -d$'\t' -f 1,2,4,9 {output[0]} > {output[1]}
    """

rule Clean_Blast_PangenomeWhole:
  input:
    rules.Align_PangenomeWholedb_GMO.output.BlastdbR2
  output:
    CB2=config["outputFolder"] + "/Cleaning_pipeline/Blast_Pangenome_wholegenome/CDS_GMO_potentialinserts_VS_Pangenomewholedb_Clean.txt"
  script:
    "scripts/Blast_AlignLength.py"

rule PangenomeWholedb_CDS_host_and_CDS_GMOInserts:
  input:
    rules.Clean_Blast_PangenomeWhole.output.CB2,
    rules.PangenomeCDSdb_CDS_host_and_CDS_GMOInserts.output.CDS_GMOpot1
  output:
    Entete2=config["outputFolder"] + "/Cleaning_pipeline/Blast_Pangenome_wholegenome/CDS_GMO_potentialinserts_VS_Pangenomewholedb_percIden95_gapopen3_qcovhsp98_EnteteUniq.txt",
    CDS_GMOpot2=config["outputFolder"] + "/Cleaning_pipeline/Blast_Pangenome_wholegenome/AfterBlastPangenomewholegenome_unalignedSeq_CDS_GMOinserts.fa",
    CDS_host2=config["outputFolder"] + "/Cleaning_pipeline/Blast_Pangenome_wholegenome/AfterBlastPangenomewholegenome_alignedSeq_CDS_host.fa"
  conda:
    "envs/Toolsbioinfo.yaml"
  shell: """
    cut -f 1 {input[0]} | uniq > {output[0]}
    filterbyname.sh in={input[1]} out={output[1]} names={output[0]} include=f
    filterbyname.sh in={input[1]} out={output[2]} names={output[0]} include=t
    """

#Final files after two Blast 
rule MEF_Pangenome_CDS_host_and_CDS_GMOInserts:
    input:
      rules.PangenomeWholedb_CDS_host_and_CDS_GMOInserts.output.CDS_GMOpot2,
      rules.PangenomeWholedb_CDS_host_and_CDS_GMOInserts.output.CDS_host2,
      rules.PangenomeCDSdb_CDS_host_and_CDS_GMOInserts.output.CDS_host1,
      rules.CDS_host_and_CDS_GMOInserts_beforeBlastPangenome.output.CDS_host
    output:
      CDS_InsertsGMOF=config["outputFolder"] + "/GMO/CDS_GMO_Inserts.fa",
      CDS_hostF=config["outputFolder"] + "/HostGenome/CDS_GMO_host.fa" 
    shell: """
      cp {input[0]} {output[0]}
      cat {input[1]} {input[2]} {input[3]}> {output[1]}
      """
###############################  CDS GMO DATABANK CLEANING

rule Align_GMOdb_hostGenome:
  input:
    config["databankGMO"] + ".fa",
    rules.MEF_Pangenome_CDS_host_and_CDS_GMOInserts.output.CDS_hostF   
  output:
    BlastdbRes3=config["outputFolder"] + "/DatabankGMO/Blast_CDSGMOhost_VS_GMOdb.tsv"
  conda:
    "envs/Aligner.yaml"
  params:
    name=config["outputFolder"] + "/DatabankGMO/GMO_Databank.db",
    evalue=0.00001,
    percidentity=90,
    threads=config["threads"]
  shell: """
    makeblastdb -in {input[0]} -out {params.name} -dbtype nucl
    blastn -query {input[1]} -db {params.name} -evalue {params.evalue} -perc_identity {params.percidentity} -num_threads {params.threads} -outfmt '6 qseqid sseqid pident length sstart send evalue sstrand' > {output[0]} 
    """

rule DeleteSeq_GMOdb_Align_HostGenome:
  input:
    config["databankGMO"] + ".fa",
    rules.Align_GMOdb_hostGenome.output.BlastdbRes3
  output:
    config["outputFolder"] + "/DatabankGMO/GMO_Databank_clean.fa"
  script:
    "scripts/Delete_Part_SeqAlign.py"

############################### DISTANCE CALCULATION

#R'mes analysis: all host CDS
rule Rmes_calcul_format_allWildCDS:
  input:
    WildCDS=config["outputFolder"] + "/HostGenome/{CDStype}" + ".fa"
  output: 
    Rmes1=config["outputFolder"] + "/Rmes/{CDStype}.Gauss.m{RMarkovOrder}.mot{RtailleMot}.{RpourcentSurRep}pSurRep.0",
    Rmes2=temp(config["outputFolder"] + "/Rmes/{CDStype}.Gauss.m{RMarkovOrder}.mot{RtailleMot}.{RpourcentSurRep}pSurRep.out"),
    Rmes3=config["outputFolder"] + "/Rmes/{CDStype}.Gauss.m{RMarkovOrder}.mot{RtailleMot}.{RpourcentSurRep}pSurRep.txt"    
  params:
    name=config["outputFolder"] + "/Rmes/{CDStype}.Gauss.m{RMarkovOrder}.mot{RtailleMot}.{RpourcentSurRep}pSurRep", 
    taillemot="{RtailleMot}",
    markov="{RMarkovOrder}"
  shell: """
    ./scripts/src/rmes --gauss -o {params.name} -l {params.taillemot} -m {params.markov} --seq {input.WildCDS}  
    ./scripts/src/rmes.format -l {params.taillemot} --tmax 0.0 --tmin 0.0 < {output.Rmes1} > {output.Rmes2} 
    head -n -2 {output.Rmes2} | sed -e '1,6d' -e '8d' -e '/\.\.\.\.\.\.\.\.\.\.\.\.\.\.\.\.\.\.\.\./d' > {output.Rmes3}
    """
rule frequence_HostGenome:
  input:
    rules.Rmes_calcul_format_allWildCDS.output.Rmes3
  output:
    freqMot=config["outputFolder"] + "/Rmes/freq_{CDStype}_l{RtailleMot}_m{RMarkovOrder}_{RpourcentSurRep}pSurRep.txt"
  conda:
    "envs/R.yaml"
  params:
    pourcentSurRep="{RpourcentSurRep}",
    pourcentSousRep=config["FrequenceParams"]["pourcentSousRep"]
  script:
    "scripts/freq_words_CDSHostGenome.R"

#Processing of datasets (GMO, Host and GMO databank) for file formatting and concatenation of third position
rule MEF_Entete_CDS:
  input:
     "{sample}.fa"
  output:
     "{sample}-MEFpos3codon.txt",
     "{sample}-MEF.txt",
     "{sample}-Entete.txt",
     "{sample}-Entetepos3codon.txt",
     "{sample}-pos3codon.fa"
  script:
     "scripts/formating_SEQpos3codon.py"

#Calculation of Bray-Curtis distances for data sets (GMO, Host and GMO databank)
rule calcul_distance_words_Freq:
  input:
    "{sample}-MEFpos3codon.txt",
    "{sample}-Entetepos3codon.txt",
    config["outputFolder"] + "/Rmes/freq_CDS_GMO_host-pos3codon_l{tailleMot}_m{MarkovOrder}_{pourcentSurRep}pSurRep.txt"
  output:
    "{sample}_l{tailleMot}_m{MarkovOrder}_{pourcentSurRep}pSurRep_calculdistBCFreq-Length.csv"
  conda:
    "envs/Python3.yaml"
  script:
    "scripts/dist_calculationBC_frequencies.py"

rule calcul_distance_words_Propv4:
  input:
    "{sample}-MEF.txt",
    "{sample}-Entete.txt",
    config["outputFolder"] + "/Rmes/freq_CDS_GMO_host_l{tailleMot}_m{MarkovOrder}_{pourcentSurRep}pSurRep.txt"
  output:
    "{sample}_l{tailleMot}_m{MarkovOrder}_{pourcentSurRep}pSurRep_calculdistBCPropv4-Length.csv"
  script:
    "scripts/dist_calculationBC_proportions.py"

rule Merge_Dist_l3m1_l4m2_l9m7:
  input:
    "{sample}_l3_m1_100pSurRep_calculdistBCPropv4-Length.csv", 
    "{sample}_l4_m2_100pSurRep_calculdistBCPropv4-Length.csv",
    "{sample}_l9_m7_10pSurRep_calculdistBCFreq-Length.csv"
  output:
    temp("{sample}_l3_m1_calculdistBCPropv4.csv"),
    temp("{sample}_l4_m2_calculdistBCPropv4-CodonUsage.csv"),
    temp("{sample}_l9m7.csv"),
    temp("{sample}_l4m2.csv"),
    "{sample}_mixl4m2Propv4l9m7Freq_calculdistBC-CodonUsage.csv"
  shell: """
    cut -f 1,2 -d$',' {input[0]} > {output[0]}
    join -1 1 -2 1 -t , <(sort -t ',' -k 1,1 {input[1]}) <(sort -t ',' -k 1,1 {output[0]}) > {output[1]} #l4 en premier et l3_m1 en second
    cut -f 1,2,4,5 -d$',' {input[2]} > {output[2]}
    join -1 1 -2 1 -t , <(sort -t ',' -k 1,1 {output[1]}) <(sort -t ',' -k 1,1 {output[2]}) > {output[3]} #l4 en premier et l9 en second
    grep -v 'Id_Sequence,BrayCurtis' {output[3]} > {output[4]}
    sed -i '1i Id_Sequence,BrayCurtis_l4_m2Prop,Length,MeanScore_l4_m2Prop,DensityNuc_l4_m2Prop,GCpourcent,CodonUsage,BrayCurtis_l9_m7Freq,MeanScore_l9_m7Freq,DensityNuc_l9_m7Freq' {output[4]}
    """

############################### MACHINE LEARNING
 
rule Suppression_tRNA:
  input:
    "{inp2}_mixl4m2Propv4l9m7Freq_calculdistBC-CodonUsage.csv"
  output:
    "{inp2}_mixl4m2Propv4l9m7Freq_calculdistBC-CodonUsage_notRNA.csv"
  shell: 
    "egrep -v -i 'tRNA-|tRNA ligase|[0-9]S ribosomal' {input} > {output}"
    
rule ML_RandomForest:
  input:
    config["outputFolder"] + "/HostGenome/CDS_GMO_host_mixl4m2Propv4l9m7Freq_calculdistBC-CodonUsage_notRNA.csv",
    config["outputFolder"] + "/DatabankGMO/GMO_Databank_clean_mixl4m2Propv4l9m7Freq_calculdistBC-CodonUsage_notRNA.csv",
    config["outputFolder"] + "/GMO/CDS_GMO_Inserts_mixl4m2Propv4l9m7Freq_calculdistBC-CodonUsage.csv"
  conda:
    "envs/R.yaml"
  output:
    config["outputFolder"] + "/MachineLearning/Pred-RFPositifs-all.txt",
    config["outputFolder"] + "/MachineLearning/Pred-RFProba-all.txt"
  script:
    "scripts/ML_PredictionGMO-RF.R"
    
