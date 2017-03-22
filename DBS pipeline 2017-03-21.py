###################
# To import vcf:
###################

from hail import *

hc = HailContext()

(hc
   .import_vcf('file:///mnt/lustre/satterst/DBS_v3/daly_dbs_exomes.vcf.bgz')
   .split_multi()
   .write('DBS_v3/DBS_v3_split.vds'))


###################
# To VEP annotate:
###################

(hc
   .read('DBS_v3/DBS_v3_split.vds')
   .vep('/mnt/lustre/cseed/vep.properties', force = True)
   .write('DBS_v3/DBS_v3_split_vep.vds'))


###################
# Do PCA and export results
###################

(hc
   .read('DBS_v3/DBS_v3_split_vep.vds')
   .filter_variants_intervals(IntervalTree.read('purcell5k.interval_list'), keep = True)
   .pca('sa.pca')
   .export_samples('DBS_v3/DBS_v3_P5k_PCA.tsv', 'Sample =s.id, sa.pca.*') )


###################
# Get a few counts
###################

DBS_v3 = hc.read('DBS_v3/DBS_v3_split_vep.vds')
counts = [
DBS_v3.count(),
DBS_v3.filter_variants_intervals(IntervalTree.read('exome_evaluation_regions.v1.interval_list'), keep = True).count(),
DBS_v3.filter_variants_expr('va.pass', keep = True).count(),
(DBS_v3.filter_variants_intervals(IntervalTree.read('exome_evaluation_regions.v1.interval_list'), keep = True)
   .filter_variants_expr('va.pass', keep = True).count())
]

# Overall:
#  nSamples = 20550
#  nVariants = 5157877
#  nCalled = 104823373505
#  callRate = 98.895%

# Exome target only:
#  nSamples = 20550
#  nVariants = 3325746
# nCalled = 67719500047
#  callRate = 99.086%

# PASSing sites only:
#  nSamples = 20550
#  nVariants = 3453522
#  nCalled = 70289353909
#  callRate = 99.041%

# PASSing sites within exome target:
#  nSamples = 20550
#  nVariants = 2214762
#  nCalled = 45145518288
#  callRate = 99.192%


###################
# Initial sample QC:
###################

(hc
   .read('DBS_v3/DBS_v3_split_vep.vds')
   .filter_variants_intervals(IntervalTree.read('exome_evaluation_regions.v1.interval_list'), keep = True)
   .filter_variants_expr('va.pass', keep = True)
   .sample_qc()
   .export_samples('DBS_v3/DBS_v3_exPASS_sampleqc.tsv', 
      'Sample = s, callRate = sa.qc.callRate, nSNP = sa.qc.nSNP, nIndel = sa.qc.nInsertion + sa.qc.nDeletion, \
      nSingleton = sa.qc.nSingleton, nNonRef = sa.qc.nNonRef, rTiTv = sa.qc.rTiTv, dpMean = sa.qc.dpMean, gqMean = sa.qc.gqMean') )

# Note: the current sampleqc function counts HomVars double toward nSNP, nInsertion, nDeletion, nSingleton; the data that
# I'm going to show don't do this


###################
# Impute sex:
###################

(hc
   .read('DBS_v3/DBS_v3_split_vep.vds')   
   .filter_variants_expr('va.pass', keep = True)
   .impute_sex()
   .export_samples('DBS_v3/DBS_v3_imputed_sex.tsv', 
      'Sample = s.id, Fstat = sa.imputesex.Fstat, ReportedSex = sa.pheno.Gender, Wave = sa.pheno.Wave') )


###################
# Begin pipeline:
###################

# Note: this pipeline will calculate variant stats incorrectly if I have samples of unknown gender.  Although I have gender information for every sample not from Pilot 2, I include a gender filtering step here.  Also, there are no variants in the PAR on the Y chromosome, so looking on the Y outside the PAR is not different than looking on the whole Y.

DBS_v3 = (hc.read('DBS_v3/DBS_v3_split_vep.vds')
   .annotate_samples_table('DBS_v3/DBS_v3_sample_annotations_2017-03-10.txt', 'Sample', root = 'sa.pheno', config = TextTableConfig(impute = True))
   .filter_samples_expr('sa.pheno.ReportedSex == sa.pheno.ImputedSex && \
      (sa.pheno.ImputedSex == "Male" || sa.pheno.ImputedSex == "Female") \
      sa.pheno.Duplicate == 0 && sa.pheno.Background == "European" && \
      (sa.pheno.Wave == "Pilot_1" || sa.pheno.Wave == "Wave_1" || sa.pheno.Wave == "Wave_2")', keep = True) )

DBS_v3.export_samples('DBS_v3/DBS_v3_initial.sample_list', 'Sample = s.id, Sex = sa.pheno.ReportedSex, \
   Wave = sa.pheno.Wave, Coverage_20x = sa.pheno.Coverage_20x, dpMean_exPASS = sa.pheno.dpMean_exPASS, \
   nSNP_exPASS = sa.pheno.nSNP_exPASS, LCSet = sa.pheno.LCSet, \
   Contamination_Pct = sa.pheno.Contamination_Pct, Chimeras_Pct = sa.pheno.Chimeras_Pct')

Count_initial = DBS_v3.count()

# nVariants: 5157877, nSamples: 17908

DBS_v3 = (DBS_v3
.filter_samples_expr('sa.pheno.Contamination_Pct > 5 || sa.pheno.Chimeras_Pct > 5', keep = False)
.filter_variants_intervals(IntervalTree.read('exome_evaluation_regions.v1.interval_list'), keep = True)
.filter_variants_expr('va.pass', keep = True)
.filter_genotypes('g.dp < 10 || \
   (g.isHomRef && (g.ad[0] / g.dp < 0.9 || g.gq < 25)) || \
   (g.isHomVar && (g.ad[1] / g.dp < 0.9 || g.pl[0] < 25)) || \
   (g.isHet && ( (g.ad[0] + g.ad[1]) / g.dp < 0.9 || g.ad[1] / g.dp < 0.25 || g.pl[0] < 25 ||  g.pAB < 0.000000001 || (sa.pheno.ImputedSex == "Male" && ( v.inXNonPar || v.inYNonPar )) )) || \
   (v.inYNonPar && sa.pheno.ImputedSex == "Female")', keep = False) )

Count_intermediate1 = DBS_v3.count()

# nVariants: 2214762, nSamples: 17847

# Get ready for first variant call rate filter and sample call rate filter

def sex_aware_variant_annotations(vds):
   [Num_males], [t] = vds.query_samples_typed('samples.filter(s => sa.pheno.ImputedSex == "Male").count()')
   [Num_females], [t] = vds.query_samples_typed('samples.filter(s => sa.pheno.ImputedSex == "Female").count()')
   vds = (vds
      .annotate_global_py('global.Num_males', Num_males, t)
      .annotate_global_py('global.Num_females', Num_females, t)
      .annotate_variants_expr('\
         va.Male_hets = gs.filter(g => sa.pheno.ImputedSex == "Male" && g.isHet).count(), \
         va.Male_homvars = gs.filter(g => sa.pheno.ImputedSex == "Male" && g.isHomVar).count(), \
         va.Male_calls = gs.filter(g => sa.pheno.ImputedSex == "Male" && g.isCalled).count(), \
         va.Female_hets = gs.filter(g => sa.pheno.ImputedSex == "Female" && g.isHet).count(), \
         va.Female_homvars = gs.filter(g => sa.pheno.ImputedSex == "Female" && g.isHomVar).count(), \
         va.Female_calls = gs.filter(g => sa.pheno.ImputedSex == "Female" && g.isCalled).count()')
      .annotate_variants_expr('\
         va.callRate = if (v.inYNonPar) ( va.Male_calls / global.Num_males ) \
         else if (v.inXNonPar) ( (va.Male_calls + 2*va.Female_calls) / (global.Num_males + 2*global.Num_females) ) \
         else ( (va.Male_calls + va.Female_calls) / (global.Num_males + global.Num_females) ), \
         va.MAC = if (v.inYNonPar) ( va.Male_homvars ) \
         else if (v.inXNonPar) ( va.Male_homvars + va.Female_hets + 2*va.Female_homvars ) \
         else ( va.Male_hets + 2*va.Male_homvars + va.Female_hets + 2*va.Female_homvars ), \
         va.MAF = if (v.inYNonPar) ( va.Male_homvars / va.Male_calls) \
         else if (v.inXNonPar) ( (va.Male_homvars + va.Female_hets + 2*va.Female_homvars) / (va.Male_calls + 2*va.Female_calls) ) \
         else ( (va.Male_hets + 2*va.Male_homvars + va.Female_hets + 2*va.Female_homvars) / (va.Male_calls + va.Female_calls) )') )
   return (vds)

def sex_aware_sample_annotations(vds):
   vds = vds.annotate_samples_expr('\
      sa.callRate = if (sa.pheno.ImputedSex == "Female") gs.filter(g => g.isCalled && !v.inYNonPar).count() / gs.filter(g => !v.inYNonPar).count() \
      else gs.fraction(g => g.isCalled)')
   return (vds)
   
DBS_v3 = sex_aware_variant_annotations(DBS_v3)
DBS_v3 = DBS_v3.filter_variants_expr('va.callRate >= 0.9 && va.MAC > 0', keep = True)
DBS_v3 = sex_aware_sample_annotations(DBS_v3)
DBS_v3 = DBS_v3.filter_samples_expr('sa.callRate >= 0.95', keep = True)

Count_intermediate2 = DBS_v3.count()

# nVariants: 1295698, nSamples: 17267

# Bring variant annotations up to date

DBS_v3 = sex_aware_variant_annotations(DBS_v3)
DBS_v3 = DBS_v3.filter_variants_expr('va.MAC > 0', keep = True)
DBS_v3 = sex_aware_sample_annotations(DBS_v3)

# Write the intermediate dataset

(DBS_v3
   .repartition(1000)
   .write('DBS_v3/DBS_v3_temp1.vds') )

# Export plink file for KING and get another count

DBS_v3 = hc.read('DBS_v3/DBS_v3_temp1.vds')
DBS_v3.export_plink('DBS_v3/DBS_v3_plink1')

print(DBS_v3.sample_schema)
print(DBS_v3.variant_schema)

Count_intermediate3 = DBS_v3.count()

# nVariants: 1262803, nSamples: 17267 (out of initial 17908)


###############
# Go to Plink:

plink2 -bfile DBS_v3_plink1 \
--snps-only \
--hwe  0.000001 \
--maf 0.01 \
--geno 0.1 \
--make-bed \
--out DBS_v3_plink1_QC

/home/unix/aganna/bin/king -b DBS_v3_plink1_QC.bed --kinship --related --prefix king

awk '{if ($9 > 0.1) print $2, $3}' king.kin > DBS_v3_related.pair_list

# Then manually pick samples to exclude (I keep the one with higher coverage) and make DBS_v3_related.sample_list.  There are 121 of them.  
#################


# Continue pipeline:

# Filter related samples -- this finalizes samples in QC'd dataset
DBS_v3 = (hc.read('DBS_v3/DBS_v3_temp1.vds')
   .filter_samples_list('DBS_v3/DBS_v3_related.sample_list', keep = False) )

DBS_v3.export_samples('DBS_v3/DBS_v3_working.sample_list', 
   'Sample = s.id, Sex = sa.pheno.ReportedSex, Wave = sa.pheno.Wave, Coverage_20x = sa.pheno.Coverage_20x, \
   dpMean_exPASS = sa.pheno.dpMean_exPASS, nSNP_exPASS = sa.pheno.nSNP_exPASS, LCSet = sa.pheno.LCSet, \
   Contamination_Pct = sa.pheno.Contamination_Pct, Chimeras_Pct = sa.pheno.Chimeras_Pct')

# Update annotations and filter variants second time --this finalizes variants in QC'd dataset
# Note that only rare (e.g. MAC <= 10) variants are being kept here; this step could easily be omitted
DBS_v3 = sex_aware_variant_annotations(DBS_v3)
DBS_v3 = DBS_v3.filter_variants_expr('va.callRate >= 0.95 && va.MAC > 0 && va.MAC <= 10', keep = True)
DBS_v3 = sex_aware_sample_annotations(DBS_v3)

Final_count = DBS_v3.count()

# nVariants: 1035616, nSamples: 17146 (as opposed to 17,124 and 1,034,830 before)

# Add useful annotations, beginning by parsing the VEP annotations
DBS_v3 = DBS_v3.annotate_variants_expr('\
   va.kyle.gene = if (isDefined(va.vep.transcript_consequences.find(c => c.canonical == 1 && c.biotype == "protein_coding" && isDefined(c.amino_acids) ))) \
   va.vep.transcript_consequences.find(c => c.canonical == 1 && c.biotype == "protein_coding" && isDefined(c.amino_acids) ).gene_symbol \
   else if (isDefined(va.vep.transcript_consequences.find(c => c.canonical == 1 && c.biotype == "protein_coding"))) \
   va.vep.transcript_consequences.find(c => c.canonical == 1 && c.biotype == "protein_coding").gene_symbol \
   else if (isDefined(va.vep.transcript_consequences.find(c => c.canonical == 1))) \
   va.vep.transcript_consequences.find(c => c.canonical == 1).gene_symbol \
   else va.vep.transcript_consequences.find(c => c.consequence_terms.toSet.contains(va.vep.most_severe_consequence)).gene_symbol, \
 \
   va.kyle.type = if (isDefined(va.vep.transcript_consequences.find(c => c.canonical == 1 && c.biotype == "protein_coding" && isDefined(c.amino_acids) ))) \
   va.vep.transcript_consequences.find(c => c.canonical == 1 && c.biotype == "protein_coding" && isDefined(c.amino_acids) ).consequence_terms \
   else if (isDefined(va.vep.transcript_consequences.find(c => c.canonical == 1 && c.biotype == "protein_coding"))) \
   va.vep.transcript_consequences.find(c => c.canonical == 1 && c.biotype == "protein_coding").consequence_terms \
   else if (isDefined(va.vep.transcript_consequences.find(c => c.canonical == 1))) \
   va.vep.transcript_consequences.find(c => c.canonical == 1).consequence_terms \
   else va.vep.transcript_consequences.find(c => c.consequence_terms.toSet.contains(va.vep.most_severe_consequence)).consequence_terms, \
\
   va.kyle.lof = if (isDefined(va.vep.transcript_consequences.find(c => c.canonical == 1 && c.biotype == "protein_coding" && isDefined(c.amino_acids) ))) \
   va.vep.transcript_consequences.find(c => c.canonical == 1 && c.biotype == "protein_coding" && isDefined(c.amino_acids) ).lof \
   else if (isDefined(va.vep.transcript_consequences.find(c => c.canonical == 1 && c.biotype == "protein_coding"))) \
   va.vep.transcript_consequences.find(c => c.canonical == 1 && c.biotype == "protein_coding").lof \
   else if (isDefined(va.vep.transcript_consequences.find(c => c.canonical == 1))) \
   va.vep.transcript_consequences.find(c => c.canonical == 1).lof \
   else va.vep.transcript_consequences.find(c => c.consequence_terms.toSet.contains(va.vep.most_severe_consequence)).lof, \
\
   va.kyle.lof_flags = if (isDefined(va.vep.transcript_consequences.find(c => c.canonical == 1 && c.biotype == "protein_coding" && isDefined(c.amino_acids) ))) \
   va.vep.transcript_consequences.find(c => c.canonical == 1 && c.biotype == "protein_coding" && isDefined(c.amino_acids) ).lof_flags \
   else if (isDefined(va.vep.transcript_consequences.find(c => c.canonical == 1 && c.biotype == "protein_coding"))) \
   va.vep.transcript_consequences.find(c => c.canonical == 1 && c.biotype == "protein_coding").lof_flags \
   else if (isDefined(va.vep.transcript_consequences.find(c => c.canonical == 1))) \
   va.vep.transcript_consequences.find(c => c.canonical == 1).lof_flags \
   else va.vep.transcript_consequences.find(c => c.consequence_terms.toSet.contains(va.vep.most_severe_consequence)).lof_flags, \
\
   va.kyle.lof_filter = if (isDefined(va.vep.transcript_consequences.find(c => c.canonical == 1 && c.biotype == "protein_coding" && isDefined(c.amino_acids) ))) \
   va.vep.transcript_consequences.find(c => c.canonical == 1 && c.biotype == "protein_coding" && isDefined(c.amino_acids) ).lof_filter \
   else if (isDefined(va.vep.transcript_consequences.find(c => c.canonical == 1 && c.biotype == "protein_coding"))) \
   va.vep.transcript_consequences.find(c => c.canonical == 1 && c.biotype == "protein_coding").lof_filter \
   else if (isDefined(va.vep.transcript_consequences.find(c => c.canonical == 1))) \
   va.vep.transcript_consequences.find(c => c.canonical == 1).lof_filter \
   else va.vep.transcript_consequences.find(c => c.consequence_terms.toSet.contains(va.vep.most_severe_consequence)).lof_filter, \
\
   va.kyle.lof_info = if (isDefined(va.vep.transcript_consequences.find(c => c.canonical == 1 && c.biotype == "protein_coding" && isDefined(c.amino_acids) ))) \
   va.vep.transcript_consequences.find(c => c.canonical == 1 && c.biotype == "protein_coding" && isDefined(c.amino_acids) ).lof_info \
   else if (isDefined(va.vep.transcript_consequences.find(c => c.canonical == 1 && c.biotype == "protein_coding"))) \
   va.vep.transcript_consequences.find(c => c.canonical == 1 && c.biotype == "protein_coding").lof_info \
   else if (isDefined(va.vep.transcript_consequences.find(c => c.canonical == 1))) \
   va.vep.transcript_consequences.find(c => c.canonical == 1).lof_info \
   else va.vep.transcript_consequences.find(c => c.consequence_terms.toSet.contains(va.vep.most_severe_consequence)).lof_info')

# Match gene to pLI score
DBS_v3 = (DBS_v3.annotate_global_list('pLI_list.txt', 'global.pLI_scores')
.annotate_global_expr('global.pLI_scores = global.pLI_scores.map(x => x.split("\\t"))')
.annotate_variants_expr('va.kyle.pLI = global.pLI_scores.find(x => x[0] == va.kyle.gene)[1].toDouble') )

# Add useful booleans
DBS_v3 = DBS_v3.annotate_variants_expr('\
   va.kyle.isIndel = if (v.altAllele.isIndel) "TRUE" else "FALSE", \
   va.kyle.isLOF = if (va.kyle.type.toSet.contains("frameshift_variant") || va.kyle.type.toSet.contains("stop_gained") || va.kyle.type.toSet.contains("splice_acceptor_variant") || va.kyle.type.toSet.contains("splice_donor_variant") || va.kyle.type.toSet.contains("transcript_ablation")) "TRUE" else "FALSE", \
   va.kyle.isMissense = if (va.kyle.type.toSet.contains("missense_variant")) "TRUE" else "FALSE", \
   va.kyle.isSynonymous = if (va.kyle.type.toSet.contains("synonymous_variant")) "TRUE" else "FALSE", \
   va.kyle.isInFrame = if (va.kyle.type.toSet.contains("inframe_insertion") || va.kyle.type.toSet.contains("inframe_deletion")) "TRUE" else "FALSE"')

# Write QC'd dataset
(DBS_v3
   .repartition(384)
   .write('DBS_v3/DBS_v3_temp2.vds') )

DBS_v3 = hc.read('DBS_v3/DBS_v3_temp2.vds')

# Export QC'd genotypes for further analysis
DBS_v3.export_genotypes('DBS_v3/DBS_v3_working_genotypes.tsv', 
   'Sample = s.id, \
   Sex = sa.pheno.ImputedSex, \
   Variant = v, \
   Chrom = v.contig, \
   Pos = v.start, \
   Gene = va.kyle.gene, \
   pLI = va.kyle.pLI, \
   Consequence = va.kyle.type, \
   isIndel = va.kyle.isIndel, \
   isLOF = va.kyle.isLOF, \
   isMissense = va.kyle.isMissense, \
   isSynonymous = va.kyle.isSynonymous, \
   isinFrame = va.kyle.isInFrame, \
   nNonRef = if (sa.pheno.ImputedSex == "Male" && (v.inXNonPar || v.inYNonPar)) g.nNonRefAlleles/2 else g.nNonRefAlleles/1, \
   MAC = va.MAC, \
   MAF = va.MAF, \
   RefReads = g.ad[0], \
   AltReads = g.ad[1], \
   AB = g.ad[1] / (g.ad[0] + g.ad[1]), \
   Depth = g.dp, \
   pAB = g.pAB, \
   GQ = g.gq, \
   PL_HomRef = g.pl[0], \
   PL_Het = g.pl[1], \
   PL_HomVar = g.pl[2], \
   lof = va.kyle.lof, \
   lof_flags = va.kyle.lof_flags, \
   lof_filter = va.kyle.lof_filter, \
   inXNonPar = if (v.inXNonPar) "TRUE" else "FALSE", \
   inYNonPar = if (v.inYNonPar) "TRUE" else "FALSE", \
   isHemizygous = if (sa.pheno.ImputedSex == "Male" && (v.inXNonPar || v.inYNonPar)) "TRUE" else "FALSE"')

