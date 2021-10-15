# MetaDecoder

An algorithm for clustering metagenomic sequences.

All datasets mentioned in text and some NEWLY ADDED datasets are available in **[Google Drive](https://drive.google.com/drive/folders/1_mybcewf3VE-7dte6oA-vDmlRx2ugzyD?usp=sharing)**.

## Benchmarks

* CAMI I Medium Complexity Dataset

|Program|Precision|Recall≥0.9|Recall≥0.8|Recall≥0.7|Recall≥0.6|Recall≥0.5|
|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
|MetaDecoder|≥0.95|108|125|135|145|146|
|MetaDecoder|≥0.9|108|126|137|147|148|
|MetaDecoder1000|≥0.95|108|128|135|140|143|
|MetaDecoder1000|≥0.9|110|131|139|145|148|
|CONCOCT|≥0.95|92|102|104|111|112|
|CONCOCT|≥0.9|97|108|112|120|121|
|MaxBin2|≥0.95|82|89|94|96|97|
|MaxBin2|≥0.9|96|106|114|117|118|
|MetaBAT2|≥0.95|103|117|125|129|135|
|MetaBAT2|≥0.9|103|118|128|132|138|
|VAMB|≥0.95|50|63|65|67|71|
|VAMB|≥0.9|55|70|72|74|78|

* CAMI I High Complexity Dataset

|Program|Precision|Recall≥0.9|Recall≥0.8|Recall≥0.7|Recall≥0.6|Recall≥0.5|
|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
|MetaDecoder|≥0.95|190|231|240|252|259|
|MetaDecoder|≥0.9|197|239|252|264|272|
|MetaDecoder1000|≥0.95|201|233|244|257|261|
|MetaDecoder1000|≥0.9|205|240|254|268|274|
|CONCOCT|≥0.95|21|23|24|24|25|
|CONCOCT|≥0.9|21|23|24|24|26|
|MaxBin2|≥0.95|95|113|125|125|129|
|MaxBin2|≥0.9|117|138|156|157|163|
|MetaBAT2|≥0.95|168|192|213|218|226|
|MetaBAT2|≥0.9|173|200|222|227|235|
|VAMB|≥0.95|147|185|196|205|207|
|VAMB|≥0.9|149|187|198|207|209|

* CAMI II Human Microbiome Project Airways Dataset

|Program|Precision|Recall≥0.9|Recall≥0.8|Recall≥0.7|Recall≥0.6|Recall≥0.5|
|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
|MetaDecoder|≥0.95|125|141|147|155|162|
|MetaDecoder|≥0.9|128|146|153|162|170|
|MetaDecoder1000|≥0.95|130|146|149|160|167|
|MetaDecoder1000|≥0.9|131|148|152|165|173|
|CONCOCT|≥0.95|8|15|19|22|27|
|CONCOCT|≥0.9|9|17|21|25|30|
|MetaBAT2|≥0.95|55|66|68|75|85|
|MetaBAT2|≥0.9|58|70|73|81|92|
|VAMB|≥0.95|20|21|22|22|22|
|VAMB|≥0.9|20|21|22|22|22|
|DASTool (MetaDecoder + MetaBAT2)|≥0.95|139|152|156|161|165|
|DASTool (MetaDecoder + MetaBAT2)|≥0.9|143|156|162|167|173|
|DASTool (CONCOCT + MetaBAT2)|≥0.95|56|68|71|73|76|
|DASTool (CONCOCT + MetaBAT2)|≥0.9|59|73|77|80|84|

* CAMI II Human Microbiome Project Gastrointestinal tract Dataset

|Program|Precision|Recall≥0.9|Recall≥0.8|Recall≥0.7|Recall≥0.6|Recall≥0.5|
|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
|MetaDecoder|≥0.95|128|140|147|149|159|
|MetaDecoder|≥0.9|132|145|153|156|166|
|MetaDecoder1000|≥0.95|113|122|130|133|140|
|MetaDecoder1000|≥0.9|116|125|135|138|146|
|CONCOCT|≥0.95|22|26|29|30|32|
|CONCOCT|≥0.9|27|31|36|37|39|
|MetaBAT2|≥0.95|90|96|102|106|106|
|MetaBAT2|≥0.9|91|97|103|107|107|
|VAMB|≥0.95|73|82|85|85|85|
|VAMB|≥0.9|74|83|86|86|86|
|DASTool (MetaDecoder + MetaBAT2)|≥0.95|142|153|157|159|164|
|DASTool (MetaDecoder + MetaBAT2)|≥0.9|148|159|165|167|172|
|DASTool (CONCOCT + MetaBAT2)|≥0.95|88|92|98|98|98|
|DASTool (CONCOCT + MetaBAT2)|≥0.9|91|95|102|102|102|

* CAMI II Human Microbiome Project Oral cavity Dataset

|Program|Precision|Recall≥0.9|Recall≥0.8|Recall≥0.7|Recall≥0.6|Recall≥0.5|
|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
|MetaDecoder|≥0.95|201|212|218|223|231|
|MetaDecoder|≥0.9|205|216|223|228|236|
|MetaDecoder1000|≥0.95|204|217|219|225|234|
|MetaDecoder1000|≥0.9|212|225|229|236|246|
|CONCOCT|≥0.95|33|35|38|43|45|
|CONCOCT|≥0.9|37|39|42|48|50|
|MetaBAT2|≥0.95|90|103|107|109|113|
|MetaBAT2|≥0.9|95|109|114|117|122|
|VAMB|≥0.95|99|111|117|125|134|
|VAMB|≥0.9|105|118|124|134|145|
|DASTool (MetaDecoder + MetaBAT2)|≥0.95|219|231|235|239|244|
|DASTool (MetaDecoder + MetaBAT2)|≥0.9|223|236|241|245|251|
|DASTool (CONCOCT + MetaBAT2)|≥0.95|90|99|103|108|108|
|DASTool (CONCOCT + MetaBAT2)|≥0.9|95|104|108|114|115|

* CAMI II Human Microbiome Project Skin Dataset

|Program|Precision|Recall≥0.9|Recall≥0.8|Recall≥0.7|Recall≥0.6|Recall≥0.5|
|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
|MetaDecoder|≥0.95|149|155|158|163|169|
|MetaDecoder|≥0.9|154|160|163|169|175|
|MetaDecoder1000|≥0.95|152|164|177|187|195|
|MetaDecoder1000|≥0.9|154|166|179|191|201|
|CONCOCT|≥0.95|14|21|29|33|37|
|CONCOCT|≥0.9|16|23|33|38|44|
|MetaBAT2|≥0.95|73|82|88|97|106|
|MetaBAT2|≥0.9|76|87|93|104|114|
|VAMB|≥0.95|55|65|73|81|92|
|VAMB|≥0.9|58|69|78|87|98|
|DASTool (MetaDecoder + MetaBAT2)|≥0.95|162|170|172|175|179|
|DASTool (MetaDecoder + MetaBAT2)|≥0.9|167|177|179|183|187|
|DASTool (CONCOCT + MetaBAT2)|≥0.95|68|81|89|93|95|
|DASTool (CONCOCT + MetaBAT2)|≥0.9|75|91|100|104|106|

* CAMI II Human Microbiome Project Urogenital tract Dataset

|Program|Precision|Recall≥0.9|Recall≥0.8|Recall≥0.7|Recall≥0.6|Recall≥0.5|
|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
|MetaDecoder|≥0.95|99|105|108|111|114|
|MetaDecoder|≥0.9|102|109|113|116|121|
|MetaDecoder1000|≥0.95|93|98|103|109|111|
|MetaDecoder1000|≥0.9|99|104|110|116|118|
|CONCOCT|≥0.95|26|33|34|35|37|
|CONCOCT|≥0.9|29|36|37|38|40|
|MetaBAT2|≥0.95|69|76|76|77|79|
|MetaBAT2|≥0.9|72|79|79|81|84|
|VAMB|≥0.95|84|86|87|90|97|
|VAMB|≥0.9|84|86|87|90|97|
|DASTool (MetaDecoder + MetaBAT2)|≥0.95|115|122|126|127|130|
|DASTool (MetaDecoder + MetaBAT2)|≥0.9|118|127|131|132|135|
|DASTool (CONCOCT + MetaBAT2)|≥0.95|69|79|79|80|80|
|DASTool (CONCOCT + MetaBAT2)|≥0.9|72|82|82|83|83|

* CAMI II Mouse gut Dataset (64 samples)

|Program|Precision|Recall≥0.9|Recall≥0.8|Recall≥0.7|Recall≥0.6|Recall≥0.5|
|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
|MetaDecoder|≥0.95|1,120|1,483|1,630|1,755|1,847|
|MetaDecoder|≥0.9|1,176|1,559|1,726|1,866|1,974|
|MetaDecoder1000|≥0.95|1,162|1,430|1,576|1,702|1,772|
|MetaDecoder1000|≥0.9|1,263|1,557|1,729|1,874|1,958|
|CONCOCT|≥0.95|804|955|1,060|1,152|1,208|
|CONCOCT|≥0.9|880|1,048|1,172|1,281|1,356|
|MaxBin2|≥0.95|669|781|838|879|893|
|MaxBin2|≥0.9|851|1,019|1,120|1,184|1,222|
|MetaBAT2|≥0.95|895|1,254|1,414|1,522|1,615|
|MetaBAT2|≥0.9|926|1,305|1,476|1,594|1,707|
|VAMB|≥0.95|277|370|427|477|545|
|VAMB|≥0.9|280|379|442|502|574|

* 11 IMG datasets

|Program|Contamination|Completeness≥0.9|Completeness≥0.8|Completeness≥0.7|Completeness≥0.6|Completeness≥0.5|
|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
|MetaDecoder|≤0.05|408|702|1,038|1,490|1,874|
|MetaDecoder|≤0.10|435|765|1,146|1,655|2,106|
|MetaDecoder1000|≤0.05|413|671|892|1,129|1,322|
|MetaDecoder1000|≤0.10|464|778|1,040|1,326|1,569|
|CONCOCT|≤0.05|174|226|269|302|328|
|CONCOCT|≤0.10|221|299|355|392|427|
|MaxBin2|≤0.05|176|256|337|449|510|
|MaxBin2|≤0.10|284|419|552|711|832|
|MetaBAT2|≤0.05|327|558|840|1,247|1,542|
|MetaBAT2|≤0.10|369|630|943|1,397|1,723|
|VAMB|≤0.05|171|280|385|500|617|
|VAMB|≤0.10|195|324|447|572|707|

* 24 HMP datasets

|Program|Contamination|Completeness≥0.9|Completeness≥0.8|Completeness≥0.7|Completeness≥0.6|Completeness≥0.5|
|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
|MetaDecoder|≤0.05|406|532|634|701|778|
|MetaDecoder|≤0.10|416|546|661|733|825|
|MetaDecoder1000|≤0.05|419|549|625|696|750|
|MetaDecoder1000|≤0.10|444|589|675|756|825|
|CONCOCT|≤0.05|297|346|378|399|415|
|CONCOCT|≤0.10|322|390|428|454|477|
|MaxBin2|≤0.05|207|235|259|281|308|
|MaxBin2|≤0.10|252|309|352|389|434|
|MetaBAT2|≤0.05|270|371|475|556|639|
|MetaBAT2|≤0.10|277|388|496|582|671|
|VAMB|≤0.05|220|310|359|404|456|
|VAMB|≤0.10|224|319|370|415|470|

* crystal geyser datasets

|Program|Contamination|Completeness≥0.9|Completeness≥0.8|Completeness≥0.7|Completeness≥0.6|Completeness≥0.5|
|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
|MetaDecoder|≤0.05|66|104|138|202|245|
|MetaDecoder|≤0.10|67|108|143|211|261|
|MetaDecoder1000|≤0.05|72|103|133|181|215|
|MetaDecoder1000|≤0.10|79|115|156|213|256|
|CONCOCT|≤0.05|36|44|56|65|72|
|CONCOCT|≤0.10|42|55|69|80|89|
|MaxBin2|≤0.05|34|47|63|72|92|
|MaxBin2|≤0.10|44|66|83|103|133|
|MetaBAT2|≤0.05|55|82|122|166|212|
|MetaBAT2|≤0.10|58|88|130|177|225|
|VAMB|≤0.05|39|49|75|101|128|
|VAMB|≤0.10|40|53|80|106|133|

* T2D datasets

|Program|Contamination|Completeness≥0.9|Completeness≥0.8|Completeness≥0.7|Completeness≥0.6|Completeness≥0.5|
|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
|MetaDecoder|≤0.05|2,014|2,540|2,943|3,243|3,496|
|MetaDecoder|≤0.10|2,040|2,588|3,020|3,363|3,651|
|MetaDecoder1000|≤0.05|1,947|2,534|2,987|3,264|3,530|
|MetaDecoder1000|≤0.10|2,019|2,655|3,151|3,459|3,765|
|CONCOCT|≤0.05|1,485|1,748|1,913|2,046|2,178|
|CONCOCT|≤0.10|1,576|1,900|2,093|2,254|2,414|
|MaxBin2|≤0.05|437|547|632|695|789|
|MaxBin2|≤0.10|692|888|1,041|1,166|1,332|
|MetaBAT2|≤0.05|1,302|1,790|2,208|2,503|2,825|
|MetaBAT2|≤0.10|1,357|1,887|2,352|2,682|3,036|
|VAMB|≤0.05|1,161|1,548|1,785|1,933|2,109|
|VAMB|≤0.10|1,187|1,588|1,838|1,995|2,182|

## Data processing

### Construct assemblies

#### Assemblies of CAMI datasets

We used short read gold-standard assemblies downloaded from CAMI using [**camiClient.jar**](https://www.microbiome-cosi.org/cami/resources/cami-client/) as follows:

```shell
# CAMI I Medium Complexity Dataset #
java -jar camiClient.jar -d https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/CAMI_I_TOY_MEDIUM . -p M1_M2_pooled_gsa_anonymous.fasta.gz
# CAMI I High Complexity Dataset #
java -jar camiClient.jar -d https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/CAMI_I_TOY_HIGH . -p H_pooled_gsa_anonymous.fasta.gz
# CAMI II HMP Airways Dataset #
java -jar camiClient.jar -d https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/CAMI_Airways . -p short_read/gsa.fasta.gz
# CAMI II HMP Gastrointestinal tract Dataset #
java -jar camiClient.jar -d https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/CAMI_Gastrointestinal_tract . -p short_read/gsa.fasta.gz
# CAMI II HMP Oral cavity Dataset #
java -jar camiClient.jar -d https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/CAMI_Oral . -p short_read/gsa.fasta.gz
# CAMI II HMP Skin Dataset #
java -jar camiClient.jar -d https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/CAMI_Skin . -p short_read/gsa.fasta.gz
# CAMI II HMP Urogenital tract Dataset #
java -jar camiClient.jar -d https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/CAMI_Urogenital_tract . -p short_read/gsa.fasta.gz
# CAMI II Mouse gut Datasets (Newly added) #
java -jar camiClient.jar -d https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/CAMISIM_MOUSEGUT . -p 19122017_mousegut_scaffolds/2017.12.29_11.37.26_sample_[0-9]+/contigs/anonymous_gsa.fasta.gz
```

#### Assemblies of IMG datasets

To download these datasets, please read the download guidelines on IMG website.

* IMG ID: 3300025317
<https://genome.jgi.doe.gov/portal/ext-api/downloads/get_tape_file?blocking=true&url=/RifCSPlow1SPAdes/download/_JAMO/5ac927be64d0b30c1371ae7a/final.contigs.fasta>

* IMG ID: 3300025323
<https://genome.jgi.doe.gov/portal/ext-api/downloads/get_tape_file?blocking=true&url=/RifCSPhig1SPAdes/download/_JAMO/5ac927d364d0b30c1371ae80/final.contigs.fasta>

* IMG ID: 3300025546
<https://genome.jgi.doe.gov/portal/ext-api/downloads/get_tape_file?blocking=true&url=/WasSprMetaSPAdes/download/_JAMO/5acd1cb964d0b3374770002d/final.contigs.fasta>

* IMG ID: 3300025737
<https://genome.jgi.doe.gov/portal/ext-api/downloads/get_tape_file?blocking=true&url=/AD_des_15/download/_JAMO/5acbc8eb64d0b337476fdb04/final.contigs.fasta>

* IMG ID: 3300026311
<https://genome.jgi.doe.gov/portal/ext-api/downloads/get_tape_file?blocking=true&url=/BiotimSDNASPAdes/download/_JAMO/5ad10e2a64d0b33747708ae7/final.contigs.fasta>

* IMG ID: 3300027784
<https://genome.jgi.doe.gov/portal/ext-api/downloads/get_tape_file?blocking=true&url=/SodLak8KL_SPAdes/download/_JAMO/5aea22ff64d0b337477363b4/final.contigs.fasta>

* IMG ID: 3300027819
<https://genome.jgi.doe.gov/portal/ext-api/downloads/get_tape_file?blocking=true&url=/Subgroc58mSPAdes/download/_JAMO/5aea1e7e64d0b33747736173/final.contigs.fasta>

* IMG ID: 3300027863
<https://genome.jgi.doe.gov/portal/ext-api/downloads/get_tape_file?blocking=true&url=/YelNatRA01SPAdes/download/_JAMO/5aeb6a5d64d0b33747736a76/final.contigs.fasta>

* IMG ID: 3300027896
<https://genome.jgi.doe.gov/portal/ext-api/downloads/get_tape_file?blocking=true&url=/FrelakHBHBSPAdes/download/_JAMO/5ae8c7b764d0b337477351ee/final.contigs.fasta>

* IMG ID: 3300027976
<https://genome.jgi.doe.gov/portal/ext-api/downloads/get_tape_file?blocking=true&url=/FRY01SPAdes/download/_JAMO/5ae8cbd764d0b33747735369/final.contigs.fasta>

* IMG ID: 3300027980
<https://genome.jgi.doe.gov/portal/ext-api/downloads/get_tape_file?blocking=true&url=/tdDd471SPAdes/download/_JAMO/5ae8c25364d0b33747735177/final.contigs.fasta>

#### Assemblies of other datasets

We used IDBA-UD (version 1.1.3) to assemble the sequencing reads.

```shell
fq2fa --merge 1.fastq 2.fastq reads.fasta

idba_ud --pre_correction --num_threads 50 --read reads.fasta --out assembly
```

### Map reads to assemblies

All sequencing reads with accessions mentioned in the text were publicly available.

#### Short-reads of CAMI datasets

```shell
# CAMI I Medium Complexity Dataset #
java -jar camiClient.jar -d https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/CAMI_I_TOY_MEDIUM . -p fq.gz
# CAMI I High Complexity Dataset #
java -jar camiClient.jar -d https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/CAMI_I_TOY_HIGH . -p fq.gz
# CAMI II HMP Airways Dataset #
java -jar camiClient.jar -d https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/CAMI_Airways . -p short_read.*.fq.gz
# CAMI II HMP Gastrointestinal tract Dataset #
java -jar camiClient.jar -d https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/CAMI_Gastrointestinal_tract . -p short_read.*.fq.gz
# CAMI II HMP Oral cavity Dataset #
java -jar camiClient.jar -d https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/CAMI_Oral . -p short_read.*.fq.gz
# CAMI II HMP Skin Dataset #
java -jar camiClient.jar -d https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/CAMI_Skin . -p short_read.*.fq.gz
# CAMI II HMP Urogenital tract Dataset #
java -jar camiClient.jar -d https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/CAMI_Urogenital_tract . -p short_read.*.fq.gz
# CAMI II Mouse gut Datasets (Newly added) #
java -jar camiClient.jar -d https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/CAMISIM_MOUSEGUT . -p 19122017_mousegut_scaffolds/2017.12.29_11.37.26_sample_[0-9]+/reads/anonymous_reads.fq.gz
```

#### Sequencing reads of IMG datasets

To download these datasets, please read the download guidelines on IMG website.

* IMG ID: 3300025317
<https://genome.jgi.doe.gov/portal/ext-api/downloads/get_tape_file?blocking=true&url=/RifCSPlow1SPAdes/download/_JAMO/5ac927be64d0b30c1371ae7b/pairedMapped.bam>

* IMG ID: 3300025323
<https://genome.jgi.doe.gov/portal/ext-api/downloads/get_tape_file?blocking=true&url=/RifCSPhig1SPAdes/download/_JAMO/5ac927d364d0b30c1371ae81/pairedMapped.bam>

* IMG ID: 3300025546
<https://genome.jgi.doe.gov/portal/ext-api/downloads/get_tape_file?blocking=true&url=/WasSprMetaSPAdes/download/_JAMO/5acd1cb964d0b3374770002e/pairedMapped.bam>

* IMG ID: 3300025737
<https://genome.jgi.doe.gov/portal/ext-api/downloads/get_tape_file?blocking=true&url=/AD_des_15/download/_JAMO/5acbc8eb64d0b337476fdb05/pairedMapped.bam>

* IMG ID: 3300026311
<https://genome.jgi.doe.gov/portal/ext-api/downloads/get_tape_file?blocking=true&url=/BiotimSDNASPAdes/download/_JAMO/5ad10e2a64d0b33747708ae8/pairedMapped.bam>

* IMG ID: 3300027784
<https://genome.jgi.doe.gov/portal/ext-api/downloads/get_tape_file?blocking=true&url=/SodLak8KL_SPAdes/download/_JAMO/5aea230064d0b337477363b5/pairedMapped.bam>

* IMG ID: 3300027819
<https://genome.jgi.doe.gov/portal/ext-api/downloads/get_tape_file?blocking=true&url=/Subgroc58mSPAdes/download/_JAMO/5aea1e7e64d0b33747736174/pairedMapped.bam>

* IMG ID: 3300027863
<https://genome.jgi.doe.gov/portal/ext-api/downloads/get_tape_file?blocking=true&url=/YelNatRA01SPAdes/download/_JAMO/5aeb6a5d64d0b33747736a77/pairedMapped.bam>

* IMG ID: 3300027896
<https://genome.jgi.doe.gov/portal/ext-api/downloads/get_tape_file?blocking=true&url=/FrelakHBHBSPAdes/download/_JAMO/5ae8c7b864d0b337477351ef/pairedMapped.bam>

* IMG ID: 3300027976
<https://genome.jgi.doe.gov/portal/ext-api/downloads/get_tape_file?blocking=true&url=/FRY01SPAdes/download/_JAMO/5ae8cbd864d0b3374773536a/pairedMapped.bam>

* IMG ID: 3300027980
<https://genome.jgi.doe.gov/portal/ext-api/downloads/get_tape_file?blocking=true&url=/tdDd471SPAdes/download/_JAMO/5ae8c25364d0b33747735178/pairedMapped.bam>

**Because MaxBin2 takes fastq files as its inputs, so we first convert all IMG bam files to fastq files using samtools as follows:**

```shell
samtools sort -n -o img_sample_i_sorted_by_read_name.bam img_sample_i.bam
rm img.bam
samtools fastq -1 sample_i.1.fastq -2 sample_i.2.fastq -0 /dev/null -s /dev/null -n img_sample_i_sorted_by_read_name.bam
rm img_sample_i_sorted_by_read_name.bam
```

We use Bowtie2 (version 2.3.4.3) to map **each sample** reads to the assembly and Samtools (version 1.3) to sort SAM formatted files (Do not remove the raw SAM files at this stage).

For each sample_i:

```shell
bowtie2-build --threads 50 assembly assembly.index

bowtie2 --threads 50 -x Dataset.assembly.index -1 sample_i.1.fastq -2 sample_i.2.fastq -S sample_i.sam

samtools sort -l 9 -@ 50 -o sample_i.bam sample_i.sam

samtools index sample_i.bam
```

### Run MetaDecoder

```shell
metadecoder coverage -s sample*.sam -o assembly.metadecoder.coverage

metadecoder seed --threads 50 -f assembly -o assembly.metadecoder.seed

metadecoder cluster [--min_sequence_length 1000] -f assembly -s assembly.metadecoder.seed -c assembly.metadecoder.coverage -o assembly.metadecoder
```

We have provided the processed coverage files by MetaDecoder of some complex datasets, you can use them instead of the bam/sam files if you need to reproduce the results.

### Run CONCOCT (version 1.0.0)

```shell
cut_up_fasta.py assembly -c 10000 -o 0 --merge_last -b contigs_10K.bed > contigs_10K.fa

concoct_coverage_table.py contigs_10K.bed sample*.bam > coverage_table.tsv

concoct -t 50 --composition_file contigs_10K.fa --coverage_file coverage_table.tsv -b concoct_output/

merge_cutup_clustering.py concoct_output/clustering_gt1000.csv > concoct_output/clustering_merged.csv

mkdir concoct_output/fasta_bins

extract_fasta_bins.py assembly concoct_output/clustering_merged.csv --output_path concoct_output/fasta_bins
```

### Run MaxBin2 (version 2.2.5)

```shell
run_MaxBin.pl -thread 50 -contig assembly -out assembly.maxbin2 -reads sample1.1.fastq -reads2 sample1.2.fastq -reads3 sample2.1.fastq -reads4 sample2.2.fastq ...
```

We have provided the processed coverage files by MaxBin2 of some complex datasets, you can use them instead of the bam/sam files if you need to reproduce the results.

### Run MetaBAT2 (version 2.12.1)

```shell
runMetaBat.sh assembly sample*.bam
```

We have provided the processed coverage files by MetaBAT2 of some complex datasets, you can use them instead of the bam/sam files if you need to reproduce the results.

### Run VAMB (version 3.0.2)

```shell
# VAMB needs unsorted bam files, you can sort them by read name. #
samtools sort -n -l 9 -@ 50 -o sample_i_vamb.bam sample_i.sam

vamb -p 50 --outdir assembly.vamb --fasta assembly --bamfiles sample*vamb.bam --minfasta 200000
```

We have provided the processed coverage files by VAMB of some complex datasets, you can use them instead of the bam/sam files if you need to reproduce the results. You may need to filter out sequences consisting of only **"N"** from CAMI assemblies before running VAMB (version 3.02), otherwise errors will be generated.

### Run DASTool (version 1.1.2)

Extract sequence ids from a set of clusters generated by a program using: **"extract_sequence_ids.py" (MetaDecoder/utilities/)**.

Please ensure that Metadecoder has been successfully installed, and you may need to download it in **MetaDecoder/utilities/**.

```shell
python3 extract_sequence_ids.py --no_header *.fasta > program.cluster
```

Run DASTool:

```shell
DAS_Tool -t 50 -c assembly -i program1.cluster,program2.cluster -o dastool
```

### Benchmarking

#### Benchmarks on simulated datasets using AMBER (version 2.0.2)

Install AMBER using pip3:

```shell
pip3 install cami-amber
```

Extract sequence ids from a set of clusters generated by a program using: **"extract_sequence_ids.py" (MetaDecoder/utilities/)**.

Please ensure that Metadecoder has been successfully installed, and you may need to download it in **MetaDecoder/utilities/**.

```shell
python3 extract_sequence_ids.py --amber *.fasta > program.cluster
```

Run AMBER:

```shell
amber.py -g *.mapping -l "program1,program2,..." program1.cluster program2.cluster ... -o "amber.benchmarks"
```

Mappings of CAMI datasets were downloaded from CAMI.

Mappings of 100 simulated genomes dataset were created using blast.

We provided the mapping files in **"*amber benchmark.tar.gz"**.

#### Benchmarks on read world datasets using CheckM (version 1.0.13)

```shell
checkm lineage_wf --tab_table --pplacer_threads 2 --threads 50 --extension fasta --file checkm.benchmarks program.clusters checkm.temp
```

#### Plot benchmarks

Download **plot_benchmarks.py** in **MetaDecoder/utilities/**.

```shell
python3 plot_benchmarks.py --input amber.benchmarks|checkm.benchmarks --output plot
```

### Relative abundance estamition using CheckM (version 1.0.13)

```shell
# Calculate coverage of sequences. #
checkm coverage -x fasta sample.program.directory_containing_fasta_clusters sample.program.coverage input_bam_files

# Calculate percentage of reads mapped to each cluster. #
checkm profile -f sample.program.abundance --tab_table sample.program.coverage
```

### Taxonomic analysis using GTDB-TK (version 1.3.0)

```shell
gtdbtk identify --genome_dir directory_containing_fasta_clusters --out_dir identify --extension fasta --cpus 30

gtdbtk align --identify_dir identify --out_dir align --cpus 30

gtdbtk classify --genome_dir directory_containing_fasta_clusters -x fasta --align_dir align --out_dir classify --cpus 30
```
