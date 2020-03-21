## Halictid genomes
### TF motif analysis
### Karen Kapheim
### 14 January 2018

BM Jones has used FIMO to identify binding motifs that change with SNPs with high Fst between social and solitary females.
She has identified 27 TFs for which binding probability is significantly changed in at least one promoter region, with alternative alleles.
Now we need to check if all 27 TFs are definitely in the _M. genalis_ genome.

##### Using files from the b10 project to create a list of orthologous genes across flies and bees.

1. Open file `TF_conservation2_withMolEvoldata.xls`
2. Used the 'table' sheet in this file to copy group, domain, bee orthologous gene IDs, and gene copy number into new spreadsheet `mgen_FIMO_TFs.xls`
  1. Saved under halictid_genome_project/analysis/TFBS
3. Used the 'lookup' sheet in this file to copy group, domain, and fly ortholog for TF motifs that were not found in the table.
4. For motifs that had B10 orthologs, matched the mgen gene based on a previously assembled table of orthologs ('/uufs/chpc.utah.edu/common/home/kapheim-group1/blastdbs/results/parsed/orthofiles/mgen_bee_orthos.txt')
  1. For the most part, each motif had a mgen gene with each ortholog as the rbh in each other bee species, so assigned this.
  2. For example:

```
grep 'GB45925' mgen_bee_orthos.txt
```
Produces:

Mgen08325       Aflo00915       GB45925-RA      Bimp08071       Bter06804       Dnov04302       Emex03309       Hlab03143       Lalb_11621      Mqua07595       Mrot01988     Nmel02597

All of these genes match the acj6_SOLEXA motif, so I assigned Mgen08325 as acj6_SOLEXA.


  3. This table was created with blastp (protein search) with evalue cutoff 10e-4
  4.
##### Use fly genes to search for mgen orthologs

For motifs without B10 orthologs, used fly gene to search for mgen gene

  1. Retrieve protein ID for each FBgn by searching the dmel translation file

```
grep 'FBgn0001185' dmel-all-translation-r6.10.fasta
```

Produces:
>FBpp0080440 type=protein; loc=2L:join(16677383..16678249,16678318..16678581,16678643..16678975); ID=FBpp0080440; name=her-PA; parent=FBgn0001185,FBtr0080883; dbxref=FlyBase:FBpp0080440,FlyBase_Annotation_IDs:CG4694-PA,GB_protein:AAF53549.1,REFSEQ:NP_476572,GB_protein:AAF53549,FlyMine:FBpp0080440,modMine:FBpp0080440; MD5=ff7558c2b29ed7aa5c99b1471b6ca555; length=487; release=r6.10; species=Dmel;
>FBpp0307137 type=protein; loc=2L:join(16677383..16678249,16678318..16678581,16678643..16678975); ID=FBpp0307137; name=her-PB; parent=FBgn0001185,FBtr0335138; dbxref=REFSEQ:NP_001260506,GB_protein:AGB93041,FlyBase:FBpp0307137,FlyBase_Annotation_IDs:CG4694-PB; MD5=ff7558c2b29ed7aa5c99b1471b6ca555; length=487; release=r6.10; species=Dmel;


So assigned this motif to both FBpp0080440 and FBpp0307137.

  2. Look for match for this fly protein in the blastp results files `rbh_mgen_dmel.txt`, `dmel_mgen.txt`, `mgen_dmel.txt` in `/uufs/chpc.utah.edu/common/home/kapheim-group1/blastdbs/rbh_fly_prot`
  3. In cases where a match was found in a one way blastp, but not rbh, manually inspected results. If less than 80% sequence identity, eliminated to be conservative.
