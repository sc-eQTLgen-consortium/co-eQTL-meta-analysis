from pathlib import Path
import re
import os

# only autosomal
CHROM = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22']

# grab the config file default
configfile: "./coQTL_wp3_CD4_T_gut.yaml"

# top mount of binds
includeDir = config["top_directory"]
# image locations
wg3_image_loc = config["wg3_image_loc"]
limix_image_loc = config["limix_image_loc"]
# location of genotype data
genotype_dir = config["genotype_dir"]
# and how each filename starts
genotype_prepend = config["genotype_prepend"]
# location of sample mapping file
smf_loc = config["smf_loc"]
# directory of the limix QTL python files
limix_dir = config["limix_dir"]

celltypes=config["celltypes"]
outputFolder=config["out_folder"]
scripts_folder = "/tools/WG3-pipeline-QTL/scripts/"

# the base limix command
limix_path="singularity exec --bind "+includeDir+" "+limix_image_loc+" python "+limix_dir
# the base R command
r_path="singularity exec --bind "+includeDir+" "+wg3_image_loc+" Rscript "

##QTL mapping variables.
genotypeFile= genotype_dir + genotype_prepend + '{chrom}'  # use {chrom} if genotype is splitted by chromosome
sampleMappingFile = smf_loc
minimum_test_samples = config["minimum_test_samples"] #the minimum number of donors that are required for a feature to be tested
#covariateFile= '' ##Not used
#kinshipFile= config["WG3_folder"]+'input/sample.kinship' ##Not used

##output of the correlation code (per chromosome)
phenotypeFile = config["correlation_folder"] + config["correlation_prepend"] + "{chrom}" + config["correlation_append"]
##Selection of tests to do SNP + feature
feature_variant_filter = config["features_file_folder"]+ config["features_file_prepend"] +'{chrom}' + config["features_file_append"]
##location of gene 1 (of the pair)
annoFile = config["limix_annotation_folder"]+ config["limix_anno1_prepend"] + '{chrom}' + config["limix_anno_append"]
##location of gene 2 (of the pair)
annoFile2 = config["limix_annotation_folder"]+ config["limix_anno2_prepend"] + '{chrom}' + config["limix_anno_append"]
##chromosome start end of the eGenes of your choice.
chunkFile = config["chunking_file_loc"]

# rename to what is used in rest of script
bgen_folder = genotype_dir
# these are just the genotype names, but they are called varfiltered in the script
var_filtered_genotypes = bgen_folder + genotype_prepend

chunk_chrom, chunk_start, chunk_end=[], [], []
with open(chunkFile) as fp:
    for line in fp:
        re_match=re.match(r"([A-Za-z0-9]+):(\d+)-(\d+)", line.strip())
        if re_match is not None:
            chunk_chrom.append(re_match[1])
            chunk_start.append(re_match[2])
            chunk_end.append(re_match[3])
        else:
            print(' '.join(['skipped chunking line', line, '!']))

qtlChunks=expand(outputFolder+"/{ct}/qtl/{chrom}_{start}_{end}/"+"{chrom}_{start}_{end}.finished", zip, chrom=chunk_chrom, start=chunk_start, end=chunk_end, allow_missing=True)
qtlTmp=expand(outputFolder+"/{ct}/qtl/{chrom}_{start}_{end}/qtl_results_all.txt.gz", zip, chrom=chunk_chrom, start=chunk_start, end=chunk_end, allow_missing=True)

wildcard_constraints:
    end = "([0-9]+)",
    start = "([0-9]+)",
    num = "([0-9]+)",
    chrom = "[0-9]{1,2}|X|Y|MT"

rule all:
    input:
        #expand(outputFolder+"/{ct}/qtl_results_all.txt", ct=celltypes),
        expand(qtlTmp, ct=celltypes),
        #expand(outputFolder/"{ct}"/"qtl.h5.tgz", ct=celltypes),
        #expand(outputFolder/"{ct}"/"qtl.annotation.tgz", ct=celltypes),
        #expand(outputFolder+"/{ct}/qtl.permutations.tgz", ct=celltypes),
        #expand(bgen_folder+"/EUR_imputed_hg38_varFiltered_chr{chrom}.bgen.sample", chrom=CHROM),
        #outputFolder/"LDMatrices.tgz",
        #expand(bgen_folder+"/EUR_imputed_hg38_stats.{chrom}.vars.gz", chrom=CHROM)
        
    output:
        touch(expand(outputFolder+"/{ct}/done.txt", ct=celltypes))


rule run_qtl_mapping:
    input:
        #af = annoFile,
        #eaf = annoFile2,
        #pf = phenotypeFile,
        smf = sampleMappingFile,
        #cf = covariateFile,
        #kf = kinshipFile,
        fvf = feature_variant_filter
        #rf = config["randomeff_files"] if config["randomeff_files"]!='' else []
    output:
        touch(outputFolder+"/{ct}/qtl/{chrom}_{start}_{end}/{chrom}_{start}_{end}.finished")
    priority:10
    params:
        pf = phenotypeFile,
        af = annoFile,
        eaf = annoFile2,
        od = outputFolder+"/{ct}/qtl/{chrom}_{start}_{end}/",
        gen  = genotypeFile,
        np = config["numberOfPermutations"],
        maf = config["minorAlleleFrequency"],
        hwe = config["hardyWeinbergCutoff"],
        w = config["windowSize"],
        mts = minimum_test_samples
    resources:
        memory = "7000",
        time = "23:59:00"
    shell:
        (limix_path + "run_QTL_analysis_metaAnalysis_LM_V2.py"
            " --bgen {params.gen} "
            " -mts {params.mts} "
            " -af {params.af} "
            " -eaf {params.eaf} "
            " -fvf {input.fvf} "
            #" -cf {input.cf} "
            " -pf {params.pf} "
            #" -rf {input.kf} "
            " -smf {input.smf} "
            " -od {params.od} "
            " -gr {wildcards.chrom}:{wildcards.start}-{wildcards.end} "
            " -np {params.np} "
            " -maf {params.maf} "
            " -c -gm gaussnorm "
            " -w {params.w} "
            " -hwe {params.hwe} "
            " -rs .8 ")

rule top_feature:
    input:
        inFile = str(outputFolder)+"/{ct}/qtl/{chrom}_{start}_{end}/{chrom}_{start}_{end}.finished"
    output:
        outF = touch( str(outputFolder)+"/{ct}/qtl/{chrom}_{start}_{end}/qtl_results_all.txt.gz")
    resources:
        memory = "25000",
        time = "4:59:00"
    params:
        dir = str(outputFolder)+"/{ct}/qtl/{chrom}_{start}_{end}/"
    run:
        shell(limix_path + "post_processing/minimal_postprocess.py  -id {params.dir} -od {params.dir} -sfo -wc")


#rule top_feature_group:
#    input:
#        inFile = str(outputFolder)+"/{ct}/qtl/{chrom}_{start}_{end}/{chrom}_{start}_{end}.finished"
#    output:
#        outF = touch( str(outputFolder)+"/{ct}/qtl/{chrom}_{start}_{end}/qtl_results_all.txt.gz"),
#        fgf = str(outputFolder)+"/{ct}/qtl/{chrom}_{start}_{end}/featureGrouping.txt"
#    resources:
#        memory = "5000",
#        time = "2:59:00"
#    params:
#        dir = str(outputFolder)+"/{ct}/qtl/{chrom}_{start}_{end}/"
#    run:
#        shell(r_path + "/scratch/hb-sceqtlgen/oneK1K/eqtlgen_coeqtl_map/MJB/Permutation/create_FeatureGroupFile.R  --in_dir {params.dir} --out_dir {params.dir}")
#        shell(limix_path + "post_processing/minimal_featureGroup_postprocess.py  -id {params.dir} -od {params.dir} -sfo -wc -fgf {output.fgf}")
#
#rule cat_qtl:
#    priority:12
#    input:
#        files = qtlTmp
#    output:
#        outputFolder+"/{ct}/qtl_results_all.txt"
#    params:
#        folder = str(outputFolder)+"/{ct}/qtl/"
#    run:
#        shell( r_path + "/scratch/hb-sceqtlgen/oneK1K/eqtlgen_coeqtl_map/MJB/Permutation/aggregateResults.R  --out_file {output} --input_folder '{params.folder}'" )


rule make_temporary_files:
    input:
        bgen_file=genotypeFile+".bgen"
    output:
        temp(var_filtered_genotypes+".bgen.z"),
        temp(var_filtered_genotypes+".bgen.sample"),
        temp(var_filtered_genotypes+".bgen_master.txt")
    shell:
        """
        singularity exec --bind {includeDir} {wg3_image_loc} python {scripts_folder}/make_z.py {input.bgen_file} {bgen_folder}
        """

rule create_bdose_file_by_chr:
    input:
        var_filtered_genotypes+"{chrom}.bgen",
        var_filtered_genotypes+"{chrom}.bgen.bgi",
        var_filtered_genotypes+"{chrom}.bgen.z",
        var_filtered_genotypes+"{chrom}.bgen.sample",
        var_filtered_genotypes+"{chrom}.bgen_master.txt"
    output:
        temp(var_filtered_genotypes+"{chrom}.bgen.bdose"),
        temp(var_filtered_genotypes+"{chrom}.bgen.bdose.bdose.tmp0"),
        temp(var_filtered_genotypes+"{chrom}.bgen.bdose.meta.tmp0")
    shell:
        """
        singularity exec --bind {includeDir} {wg3_image_loc} /tools/ldstore_v2.0_x86_64/ldstore_v2.0_x86_64 --in-files {var_filtered_genotypes}{wildcards.chrom}.bgen_master.txt --write-bdose --bdose-version 1.1
        """

rule Genotype_IO:
    input:
        bgen_folder + genotype_prepend + '{chrom}' + ".bgen"
    params:
        bgen_folder + genotype_prepend + '{chrom}' + ".bgen"
    output:
        bgen_folder+"/geno_stats.{genotype_prepend}{chrom}.vars.gz",
        bgen_folder+"/geno_stats.{genotype_prepend}{chrom}.samples.gz"
    shell:
        """
        singularity exec --bind {includeDir} {wg3_image_loc} java -Xmx5g -Xms5g -jar /tools/Genotype-IO-1.0.6-SNAPSHOT-jar-with-dependencies.jar -I BGEN -i {params} -o {bgen_folder}/geno_stats.{genotype_prepend}{wildcards.chrom}
        gzip {bgen_folder}/geno_stats.{genotype_prepend}{wildcards.chrom}*
        """