<!-- MarkdownTOC -->

- [Polymorphism Identification](#polymorphism-identification)
    - [Set variables](#set-variables)
    - [Identify polymorphisms relative to Heinz](#identify-polymorphisms-relative-to-heinz)
    - [Filter based on reciprocal coverage](#filter-based-on-reciprocal-coverage)
    - [Filter based on M82 vs. PEN comparison](#filter-based-on-m82-vs-pen-comparison)
    - [Polymorphism Noise Reduction](#polymorphism-noise-reduction)
    - [Plot Parental Samples Before Noise Reduction](#plot-parental-samples-before-noise-reduction)
    - [Genotype and Plot Parental Samples After Noise Reduction](#genotype-and-plot-parental-samples-after-noise-reduction)
    - [Merge, Genotype, and Plot BIL Data](#merge-genotype-and-plot-bil-data)
    - [Filter Polymorphisms Based on Allele Ratios and Coverage](#filter-polymorphisms-based-on-allele-ratios-and-coverage)
    - [Plot Merged and Filtered BIL Data](#plot-merged-and-filtered-bil-data)
    - [Filter Polymorphism Files Based on Positions Observed in Filtered Genotyped Files](#filter-polymorphism-files-based-on-positions-observed-in-filtered-genotyped-files)
    - [Plot Parental Samples After Filtering](#plot-parental-samples-after-filtering)
    - [Plot Parental Samples After Removing Parental-Specific High-Coverage Outliers](#plot-parental-samples-after-removing-parental-specific-high-coverage-outliers)
    - [Genotype and Plot BIL Samples](#genotype-and-plot-bil-samples)

<!-- /MarkdownTOC -->


# Polymorphism Identification

We used [SNPTools](https://github.com/mfcovington/SNPtools) v0.2.4 to identify find SNPs and INDELs between *Solanum lycopersicum* cv. M82 ('M82') and *S. pennellii* ('PEN'). SNPTools v0.2.6 was used to make current version of genotype plots. The polymorphism detection results do not differ between these two versions.


The alignment files used were:

- `M82_unreped_repeat_filtered.sorted.bam`
- `PEN_unreped_repeat_filtered.sorted.bam`


They contain reads that have been mapped to the *S. lycopersicum* cv. Heinz 1706 ('Heinz') genomic reference:

- `S_lycopersicum_chromosomes.2.40.fa`


## Set variables

Set shell variables for (directory paths, parameters, etc.) that will be used for multiple scripts:

```sh
BIN=        # Directory containing 'SNPtools/bin'
OUT_DIR=    # Output directory
FA=         # Directory containing 'S_lycopersicum_chromosomes.2.40.fa'
BAM_DIR=    # Directory containing BAM files

THREADS=    # Set number of processors to use

PAR1=M82
PAR2=PEN

# We explicity specify the chromosomes to analyze in order to bypass 'SL2.40ch00'
SEQ_LIST=SL2.40ch01,SL2.40ch02,SL2.40ch03,SL2.40ch04,SL2.40ch05,SL2.40ch06,SL2.40ch07,SL2.40ch08,SL2.40ch09,SL2.40ch10,SL2.40ch11,SL2.40ch12
```


## Identify polymorphisms relative to Heinz

Find putative M82 vs. Heinz polymorphisms and putative PEN vs. Heinz polymorphisms:

```sh
for ID in $PAR1 $PAR2; do
    $BIN/SNPfinder/snp_finder.pl \
      --id        $ID \
      --bam       $BAM_DIR/${ID}_unreped_repeat_filtered.sorted.bam \
      --fasta     $FA \
      --seq_list  $SEQ_LIST \
      --out_dir   $OUT_DIR \
      --snp_min   0.66 \
      --indel_min 0.33 \
      --threads   $THREADS

    for SNP_FILE in $OUT_DIR/snps/$ID.*.snps.nogap.gap.csv; do
        $BIN/SNPfinder/02.0.filtering_SNPs_by_pos.pl $SNP_FILE
    done
done
```


## Filter based on reciprocal coverage

With putative polymorphisms vs. Heinz in hand, the first step in identifying whether any of them are also polymorphisms between M82 and PEN is to filter out positions in the genome where a polymorphism between M82 or PEN vs. Heinz has insufficient sequence depth (<4 reads) in the other parental genotype. For example, If a putative SNP was identified between PEN and Heinz, but the only three sequencing read were aligned to this region in the M82 BAM file, the position would be discarded and no longer considered for whether it was also a SNP between M82 and PEN.

For polymorphisms between Heinz and another genotype, determine coverage in the other genotype and filter out low coverage (<4 reads):

```sh
$BIN/Coverage/reciprocal_coverage.pl \
  --bam      $BAM_DIR/${PAR1}_unreped_repeat_filtered.sorted.bam \
  --par1     $PAR1 \
  --par2     $PAR2 \
  --par1_bam $BAM_DIR/${PAR1}_unreped_repeat_filtered.sorted.bam \
  --par2_bam $BAM_DIR/${PAR2}_unreped_repeat_filtered.sorted.bam \
  --seq_list $SEQ_LIST \
  --out_dir  $OUT_DIR \
  --threads  $THREADS

for CHR in ${SEQ_LIST//,/ }; do
    $BIN/SNPfinder/coverage_filter.pl \
      --chr  $CHR \
      --snp1 $OUT_DIR/snps/$PAR1.$CHR.snps.nogap.gap.FILTERED.csv \
      --snp2 $OUT_DIR/snps/$PAR2.$CHR.snps.nogap.gap.FILTERED.csv \
      --par1 $PAR1 \
      --par2 $PAR2 \
      --out  $OUT_DIR
done
```


## Filter based on M82 vs. PEN comparison

For a PEN vs. Heinz SNP, for example, the M82 genotype could match PEN, Heinz, or a third allele. We classify these as 'NOT', 'SNP', and 'DIFF_SNP', respectively (see table output from `SNPfinder/classify-snps.r` R script below). We are not interested in the 'NOT' class and will filter these out. This class represents M82/PEN vs. Heinz polymorphisms, in which we are not interested.

Classify and filter by polymorphism type:

```sh
Rscript --vanilla $BIN/SNPfinder/classify-snps.r $OUT_DIR/master_snp_lists

mkdir $OUT_DIR/snp_master
for CHR in ${SEQ_LIST//,/ }; do
    grep -ve "NOT" \
    $OUT_DIR/master_snp_lists/master_snp_list.$PAR1.vs.$PAR2.$CHR.classified \
    > $OUT_DIR/snp_master/polyDB.$CHR
done
```


Output from `SNPfinder/classify-snps.r`:

             chr             pos           ref_base    snp_base    genotype    insert_position      SNP_CLASS
     SL2.40ch01:94724   Min.   :    1621   A  :24825   A  :23542   M82: 1347   Min.   : 1.00     SNP     :93588
                        1st Qu.:32752144   C  :17967   C  :19787   PEN:93377   1st Qu.: 1.00     DIFF_SNP:   58
                        Median :70989128   G  :17846   G  :19647               Median : 2.00     NOT     : 1078
                        Mean   :57828276   INS: 8978   T  :23274               Mean   : 3.25     CONFLICT:    0
                        3rd Qu.:80856758   T  :25108   del: 8474               3rd Qu.: 5.00
                        Max.   :90258818                                       Max.   :16.00
                                                                               NA's   :85746


             chr             pos           ref_base    snp_base    genotype    insert_position      SNP_CLASS
     SL2.40ch02:68141   Min.   :   18215   A  :17443   A  :16973   M82: 1337   Min.   : 1.00     SNP     :66983
                        1st Qu.:26122582   C  :12963   C  :13935   PEN:66804   1st Qu.: 1.00     DIFF_SNP:   28
                        Median :36050916   G  :12925   G  :13922               Median : 2.00     NOT     : 1130
                        Mean   :33158868   INS: 7079   T  :16859               Mean   : 3.38     CONFLICT:    0
                        3rd Qu.:42908358   T  :17731   del: 6452               3rd Qu.: 5.00
                        Max.   :49863920                                       Max.   :16.00
                                                                               NA's   :61062


             chr             pos           ref_base    snp_base    genotype    insert_position      SNP_CLASS
     SL2.40ch03:68000   Min.   :   17032   A  :17732   A  :16744   M82: 1663   Min.   : 1.00     SNP     :66450
                        1st Qu.:16288605   C  :12817   C  :14098   PEN:66337   1st Qu.: 1.00     DIFF_SNP:   68
                        Median :47806219   G  :12807   G  :13837               Median : 2.00     NOT     : 1482
                        Mean   :39229032   INS: 6725   T  :16672               Mean   : 3.24     CONFLICT:    0
                        3rd Qu.:58398145   T  :17919   del: 6649               3rd Qu.: 4.00
                        Max.   :64834914                                       Max.   :16.00
                                                                               NA's   :61275


             chr             pos           ref_base    snp_base    genotype    insert_position      SNP_CLASS
     SL2.40ch04:74969   Min.   :     565   A  :19791   A  :18563   M82: 8631   Min.   : 1.00     SNP     :66981
                        1st Qu.:13338602   C  :14063   C  :15909   PEN:66338   1st Qu.: 1.00     DIFF_SNP:  346
                        Median :40819865   G  :14294   G  :16034               Median : 2.00     NOT     : 7642
                        Mean   :35662285   INS: 6628   T  :18369               Mean   : 3.07     CONFLICT:    0
                        3rd Qu.:56182546   T  :20193   del: 6094               3rd Qu.: 4.00
                        Max.   :64058899                                       Max.   :16.00
                                                                               NA's   :68341


             chr             pos           ref_base    snp_base    genotype    insert_position      SNP_CLASS
     SL2.40ch05:63865   Min.   :     561   A  :17220   A  :15643   M82: 8692   Min.   : 1        SNP     :55225
                        1st Qu.: 8014952   C  :12192   C  :13676   PEN:55173   1st Qu.: 1        DIFF_SNP:  452
                        Median :26355636   G  :12180   G  :13978               Median : 2        NOT     : 8188
                        Mean   :30289130   INS: 5245   T  :15773               Mean   : 3        CONFLICT:    0
                        3rd Qu.:54656292   T  :17028   del: 4795               3rd Qu.: 4
                        Max.   :65015922                                       Max.   :16
                                                                               NA's   :58620


             chr             pos           ref_base    snp_base    genotype    insert_position      SNP_CLASS
     SL2.40ch06:58374   Min.   :     808   A  :15174   A  :14354   M82:  668   Min.   : 1.00     SNP     :57844
                        1st Qu.:17863680   C  :10962   C  :12205   PEN:57706   1st Qu.: 1.00     DIFF_SNP:   24
                        Median :33731216   G  :11054   G  :12319               Median : 2.00     NOT     :  506
                        Mean   :28440993   INS: 5854   T  :14438               Mean   : 3.38     CONFLICT:    0
                        3rd Qu.:40084982   T  :15330   del: 5058               3rd Qu.: 5.00
                        Max.   :45999028                                       Max.   :16.00
                                                                               NA's   :52520


             chr             pos           ref_base    snp_base    genotype    insert_position      SNP_CLASS
     SL2.40ch07:62177   Min.   :    3311   A  :16715   A  :15297   M82:  696   Min.   : 1.0      SNP     :61529
                        1st Qu.:10385343   C  :11710   C  :13186   PEN:61481   1st Qu.: 1.0      DIFF_SNP:   28
                        Median :40709718   G  :11647   G  :13399               Median : 2.0      NOT     :  620
                        Mean   :35570593   INS: 5198   T  :15032               Mean   : 3.1      CONFLICT:    0
                        3rd Qu.:58042128   T  :16907   del: 5263               3rd Qu.: 4.0
                        Max.   :65265446                                       Max.   :16.0
                                                                               NA's   :56979


             chr             pos           ref_base    snp_base    genotype    insert_position      SNP_CLASS
     SL2.40ch08:57926   Min.   :     265   A  :15698   A  :13994   M82:  577   Min.   : 1.00     SNP     :57360
                        1st Qu.:15764190   C  :10827   C  :12250   PEN:57349   1st Qu.: 1.00     DIFF_SNP:   16
                        Median :46491998   G  :10781   G  :12377               Median : 2.00     NOT     :  550
                        Mean   :37228768   INS: 5102   T  :14285               Mean   : 3.05     CONFLICT:    0
                        3rd Qu.:56887311   T  :15518   del: 5020               3rd Qu.: 4.00
                        Max.   :63030851                                       Max.   :16.00
                                                                               NA's   :52824


             chr             pos           ref_base    snp_base    genotype    insert_position      SNP_CLASS
     SL2.40ch09:56810   Min.   :    3305   A  :15215   A  :14068   M82:  932   Min.   : 1.00     SNP     :55960
                        1st Qu.: 8413692   C  :10708   C  :11793   PEN:55878   1st Qu.: 1.00     DIFF_SNP:   26
                        Median :45363208   G  :10722   G  :12147               Median : 2.00     NOT     :  824
                        Mean   :37123906   INS: 5174   T  :13997               Mean   : 3.05     CONFLICT:    0
                        3rd Qu.:61077315   T  :14991   del: 4805               3rd Qu.: 4.00
                        Max.   :67660926                                       Max.   :16.00
                                                                               NA's   :51636


             chr             pos           ref_base    snp_base    genotype    insert_position      SNP_CLASS
     SL2.40ch10:56994   Min.   :    2859   A  :15517   A  :13786   M82: 1174   Min.   : 1.00     SNP     :55764
                        1st Qu.:12654449   C  :10577   C  :12632   PEN:55820   1st Qu.: 1.00     DIFF_SNP:   58
                        Median :40235318   G  :10678   G  :12311               Median : 2.00     NOT     : 1172
                        Mean   :35814420   INS: 4666   T  :13958               Mean   : 2.99     CONFLICT:    0
                        3rd Qu.:58525197   T  :15556   del: 4307               3rd Qu.: 4.00
                        Max.   :64831636                                       Max.   :16.00
                                                                               NA's   :52328


             chr             pos           ref_base    snp_base    genotype    insert_position      SNP_CLASS
     SL2.40ch11:53061   Min.   :    7833   A  :14259   A  :12996   M82: 3655   Min.   : 1.00     SNP     :49405
                        1st Qu.: 6119271   C  :10074   C  :11366   PEN:49406   1st Qu.: 1.00     DIFF_SNP:  160
                        Median :22917656   G  : 9932   G  :11210               Median : 2.00     NOT     : 3496
                        Mean   :25489997   INS: 4590   T  :13008               Mean   : 3.02     CONFLICT:    0
                        3rd Qu.:46253292   T  :14206   del: 4481               3rd Qu.: 4.00
                        Max.   :53385549                                       Max.   :16.00
                                                                               NA's   :48471


             chr             pos           ref_base    snp_base    genotype    insert_position      SNP_CLASS
     SL2.40ch12:57998   Min.   :    5755   A  :15929   A  :14127   M82: 1286   Min.   : 1.00     SNP     :56846
                        1st Qu.: 9128241   C  :10922   C  :12553   PEN:56712   1st Qu.: 1.00     DIFF_SNP:   56
                        Median :32704256   G  :10947   G  :12713               Median : 2.00     NOT     : 1096
                        Mean   :31862180   INS: 4462   T  :14138               Mean   : 2.79     CONFLICT:    0
                        3rd Qu.:50129150   T  :15738   del: 4467               3rd Qu.: 4.00
                        Max.   :65486162                                       Max.   :16.00
                                                                               NA's   :53536


## Polymorphism Noise Reduction

To reduce false positive polymorphisms, we perform a noise-reduction step. This consists of checking the genotypes of the M82 and PEN sequence alignments for each putatitve polymorphism. Since these are the same sequence alignments used for polymorphism identification, any polymorphisms that do not exhibit the expected allele for both parents of the population are rejected.

Genotype parents and perform noise reduction:

```sh
$BIN/Genotype/genotype_parents+nr.pl \
  --id1         $PAR1 \
  --id2         $PAR2 \
  --bam1        $BAM_DIR/${PAR1}_unreped_repeat_filtered.sorted.bam \
  --bam2        $BAM_DIR/${PAR2}_unreped_repeat_filtered.sorted.bam \
  --fasta       $FA \
  --seq_list    $SEQ_LIST \
  --out_dir     $OUT_DIR \
  --nr_ratio    0.9 \
  --threads     $THREADS
```


## Plot Parental Samples Before Noise Reduction

```sh
for ID in $PAR1 $PAR2; do
    $BIN/Plot/genoplot_by_id.pl \
      --id          $ID \
      --par1        $PAR1 \
      --par2        $PAR2 \
      --bam         $BAM_DIR/${ID}_unreped_repeat_filtered.sorted.bam \
      --seq_list    $SEQ_LIST \
      --out_dir     $OUT_DIR \
      --col_par1    magenta \
      --col_par2    green \
      --no_nr \
      --chr_pat     SL2.40 \
      --chr_sub     ''
done
```

![M82 (pre-noise-reduction)](../data/genoplot/M82.before_nr.png "M82 (pre-noise-reduction)")
![PEN (pre-noise-reduction)](../data/genoplot/PEN.before_nr.png "PEN (pre-noise-reduction)")


## Genotype and Plot Parental Samples After Noise Reduction

```sh
for ID in $PAR1 $PAR2; do
    $BIN/Genotype/extract+genotype_pileups.pl \
      --id          $ID \
      --par1        $PAR1 \
      --par2        $PAR2 \
      --bam         $BAM_DIR/${ID}_unreped_repeat_filtered.sorted.bam \
      --fasta       $FA \
      --seq_list    $SEQ_LIST \
      --out_dir     $OUT_DIR \
      --threads     $THREADS

    $BIN/Plot/genoplot_by_id.pl \
      --id          $ID \
      --par1        $PAR1 \
      --par2        $PAR2 \
      --bam         $BAM_DIR/${ID}_unreped_repeat_filtered.sorted.bam \
      --seq_list    $SEQ_LIST \
      --out_dir     $OUT_DIR \
      --col_par1    magenta \
      --col_par2    green \
      --chr_pat     SL2.40 \
      --chr_sub     ''
done
```

![M82 (post-noise-reduction)](../data/genoplot/M82.png "M82 (post-noise-reduction)")
![PEN (post-noise-reduction)](../data/genoplot/PEN.png "PEN (post-noise-reduction)")


## Merge, Genotype, and Plot BIL Data

Merge all BIL sequence alignments:

```sh
samtools merge $BAM_DIR/BIL.merged.bam $BAM_DIR/BIL-SL2.40ch*-merged.bam
```

Genotype and plot merged BIL data:

```sh
$BIN/Genotype/extract+genotype_pileups.pl \
  --id          BILs_merged \
  --par1        $PAR1 \
  --par2        $PAR2 \
  --bam         $BAM_DIR/BIL.merged.bam \
  --fasta       $FA \
  --seq_list    $SEQ_LIST \
  --out_dir     $OUT_DIR \
  --threads     $THREADS

$BIN/Plot/genoplot_by_id.pl \
  --id          BILs_merged \
  --par1        $PAR1 \
  --par2        $PAR2 \
  --bam         $BAM_DIR/BIL.merged.bam \
  --seq_list    $SEQ_LIST \
  --out_dir     $OUT_DIR \
  --col_par1    magenta \
  --col_par2    green \
  --chr_pat     SL2.40 \
  --chr_sub     ''
```

![BILs merged (pre-filtering)](../data/genoplot/BILs_merged.png "BILs merged (pre-filtering)")


## Filter Polymorphisms Based on Allele Ratios and Coverage

We filter genotype data from the merged BILs to remove outliers by doing the following for each chromosome (separately):

- For each polymorphism, calculate the ratio of merged BIL reads that matches the PEN allele
- Keep data only for polymorphisms with a PEN allele ratio between the 25th and 75th percentile
- For each polymorphism, calculate the ratio of merged BIL reads that matches neither the M82 nor the PEN allele
- Keep data only for polymorphisms with a non-parental allele ratio less than 0.025

Since a small number of positions have excessive coverage relative the rest of the genome, we keep data only for polymorphism with a maximum coverage of 2000 reads.

Filter based on allele ratios and coverage in R:

```r
setwd('genotyped')

prefix <- 'BILs_merged'
suffix <- 'genotyped.nr'
file.list <- list.files(pattern = paste(prefix, 'SL2.40ch\\d+', suffix,
                                        sep = '.'))
for (file in file.list){
  if (exists('chr.genos')) rm(chr.genos)
  if (exists('chr.genos.flt')) rm(chr.genos.flt)
  chr.genos <- read.table(file, head = F, as.is = T, sep = "\t")

  colnames(chr.genos) <- c('chr', 'pos', 'par1', 'par2', 'tot')
  chr.genos$par2_rat <- chr.genos$par2 / chr.genos$tot
  chr.genos$non_par <- chr.genos$tot - chr.genos$par1 - chr.genos$par2
  chr.genos$non_par_rat <- chr.genos$non_par / chr.genos$tot

  limits <- quantile(chr.genos$par2_rat, na.rm = T, probs = c(0.25, 0.75))
  chr.genos.flt <- chr.genos[(chr.genos$par2_rat >= limits[1] &
                   chr.genos$par2_rat <= limits[2]),]
  chr.genos.flt <- na.omit(chr.genos.flt)
  chr.genos.flt <- chr.genos.flt[chr.genos.flt$non_par_rat < 0.025,]
  chr.genos.flt <- chr.genos.flt[chr.genos.flt$par1 <= 2000, ]
  file_name <- paste(prefix, 'filtered', chr.genos.flt$chr[1], suffix,
                     sep = '.')
  write.table(chr.genos.flt[,1:5], file = file_name, sep = '\t',
              row.names = F, col.names = F, quote = F)
}
```


## Plot Merged and Filtered BIL Data

```sh
$BIN/Plot/genoplot_by_id.pl \
  --id          BILs_merged.filtered \
  --par1        $PAR1 \
  --par2        $PAR2 \
  --bam         $BAM_DIR/BIL.merged.bam \
  --seq_list    $SEQ_LIST \
  --out_dir     $OUT_DIR \
  --col_par1    magenta \
  --col_par2    green \
  --chr_pat     SL2.40 \
  --chr_sub     ''
```

![BILs merged (post-filtering)](../data/genoplot/BILs_merged.filtered.png "BILs merged (post-filtering)")


## Filter Polymorphism Files Based on Positions Observed in Filtered Genotyped Files

Only the genotyped files were [filtered based on polymorphism allele ratios and coverage](#filter-polymorphisms-based-on-allele-ratios-and-coverage). Therefore, it is necessary to filter the polymorphism files using the same criteria. Instead of performing additional calculations, we can simply remove positions in the polymorphism files that are no longer in the filtered genotyped files.

```sh
$BIN/Other/filter-snps-based-on-genotyped-positions.pl \
  --snp_in_dir $OUT_DIR/snp_master \
  --snp_out_dir $OUT_DIR/snp_master \
  $OUT_DIR/genotyped/BILs_merged.filtered.SL2.40ch*.genotyped.nr
```


## Plot Parental Samples After Filtering

Now that the [polymorphism files have been filtered](#filter-polymorphism-files-based-on-positions-observed-in-filtered-genotyped-files), we will plot the parental samples. Instead of re-genotyping the parents, it is faster just to filter their current genotyped files using the filtered polymorphism files.

```sh
for ID in $PAR1 $PAR2; do
    $BIN/Other/filter-genotyped-files-based-on-polymorphism-positions.pl \
      --out $OUT_DIR \
      --id $ID

    $BIN/Plot/genoplot_by_id.pl \
      --id          $ID.filtered \
      --par1        $PAR1 \
      --par2        $PAR2 \
      --bam         $BAM_DIR/${ID}_unreped_repeat_filtered.sorted.bam \
      --seq_list    $SEQ_LIST \
      --out_dir     $OUT_DIR \
      --col_par1    magenta \
      --col_par2    green \
      --chr_pat     SL2.40 \
      --chr_sub     ''

    mv $OUT_DIR/genoplot/$ID.filtered.png $OUT_DIR/genoplot/$ID.filtered.with-outliers.png

done
```

![M82 (post-filtering, with-outliers)](../data/genoplot/M82.filtered.with-outliers.png "M82 (post-filtering, with-outliers)")
![PEN (post-filtering, with-outliers)](../data/genoplot/PEN.filtered.with-outliers.png "PEN (post-filtering, with-outliers)")


## Plot Parental Samples After Removing Parental-Specific High-Coverage Outliers

More than 99.999% of the polymorphism positions have parental genotype information from 4 to 186 reads per parent. There are 7 positions with coverage in the thousands for one or both of the parents. We remove these outliers from the filtered parental genotyped files and from the polymorphism database.

```sh
for FILE in $OUT_DIR/genotyped/M82.filtered.*.nr $OUT_DIR/genotyped/PEN.filtered.*.nr $OUT_DIR/snp_master/polyDB.*.nr; do
    grep -v $FILE > $FILE.filtered \
      -e $'SL2.40ch01\t1866676'  -e $'SL2.40ch02\t15236040' -e $'SL2.40ch05\t58235463' \
      -e $'SL2.40ch05\t58235490' -e $'SL2.40ch12\t870464'   -e $'SL2.40ch12\t27506567' \
      -e $'SL2.40ch12\t27506593'

    mv $FILE.filtered $FILE
done

for ID in $PAR1 $PAR2; do
    $BIN/Plot/genoplot_by_id.pl \
      --id          $ID.filtered \
      --par1        $PAR1 \
      --par2        $PAR2 \
      --bam         $BAM_DIR/${ID}_unreped_repeat_filtered.sorted.bam \
      --seq_list    $SEQ_LIST \
      --out_dir     $OUT_DIR \
      --col_par1    magenta \
      --col_par2    green \
      --chr_pat     SL2.40 \
      --chr_sub     ''
done
```

![M82 (post-filtering)](../data/genoplot/M82.filtered.png "M82 (post-filtering)")
![PEN (post-filtering)](../data/genoplot/PEN.filtered.png "PEN (post-filtering)")


## Genotype and Plot BIL Samples

Now that the polymorphisms have gone through multiple types of filtering to remove noise, we can use the remaining polymorphism to genotype each individual of the BIL population. The data and plots for every BIL are included in the unabridged version of this repository; however, only BIL 009 is included as an example in the abridged version.

```sh
for ID in `perl -E '@files = @ARGV;
                    push @ids, $_ =~ /(BIL_\d+)_Slc\.sorted\.bam/g for @files;
                    print "@ids";' $BAM_DIR/BIL_*_Slc.sorted.bam`; do
    echo $ID

    $BIN/Genotype/extract+genotype_pileups.pl \
      --id          $ID \
      --par1        $PAR1 \
      --par2        $PAR2 \
      --bam         $BAM_DIR/${ID}_Slc.sorted.bam \
      --fasta       $FA \
      --seq_list    $SEQ_LIST \
      --out_dir     $OUT_DIR \
      --threads     $THREADS

    $BIN/Plot/genoplot_by_id.pl \
      --id          $ID \
      --par1        $PAR1 \
      --par2        $PAR2 \
      --bam         $BAM_DIR/${ID}_Slc.sorted.bam \
      --seq_list    $SEQ_LIST \
      --out_dir     $OUT_DIR \
      --col_par1    magenta \
      --col_par2    green \
      --chr_pat     SL2.40 \
      --chr_sub     ''
done
```

![BIL 009 Genotype Plot](../data/genoplot/BIL_009.png "BIL 009 Genotype Plot")

