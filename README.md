## Introduction

**sikiclass** is a bioinformatics pipeline that performs downstream classification analysis for UMI-collapsed reads from sikipipe.

<!-- TODO nf-core:
   Complete this sentence with a 2-3 sentence summary of what types of data the pipeline ingests, a brief overview of the
   major pipeline sections and the types of output it produces. You're giving an overview to someone new
   to nf-core here, in 15-20 seconds. For an example, see https://github.com/nf-core/rnaseq/blob/master/README.md#introduction
-->

<!-- TODO nf-core: Include a figure that guides the user through the major workflow steps. Many nf-core
     workflows use the "tube map" design for that. See https://nf-co.re/docs/contributing/design_guidelines#examples for examples.   -->
<!-- TODO nf-core: Fill in short bullet-pointed list of the default steps in the pipeline -->

Implemented classification steps include `classify_tag` that divides reads according to tag fragment occurrence, `classify_single_tag` that further divides single-tag-containing reads into sub-categories like "precise_tag", "5' INDEL", "3' INDEL", etc., `classify_no_tag` that furhter divides no-tag-containing reads into sub-categories like "Deletion", "Insertion", etc. Moreover, for "precise_tag" class, reads with different SNPs are splitted and counted. Finally, statistical tables summarizing the fraction of each read category as well as the distribution of size and location of INDELs are generated. 

The diagram below illustrates the method in detail. Specifially, for `classify_tag` (**Step1**), to determine if tag appears in each read sequence, we leverage short-read aligner BWA, use tag as query and each read as reference. By parsing the resulting BAM file, reads are classified as "without tag" if no tag appear in the read, "with single tag" if every tag segment only appears at most one time in the read, and "with multiple tag" if any tag segment appears more than one time in read. We further classify "without tag" reads (`classify_single_tag`, **Step2**) and "with single tag" reads (`classify_no_tag`, **Step3**) by adopting long-read aligner minimap2. For "with single tag" reads, they are mapped against a reference with precise tag inserted. By parsing the BAM file, for each read, if there are no INDELs occuring in tag plus a predefined flanking region, it will be regarded as read "with precise tag". Otherwise, it will be considered to be with INDELs. Similarly, for "without tag" reads, they are mapped against wild type reference that does not contain the tag insert. The read will be classified into sub-categories based on the occurrence of INDELs. In addition, for "precise tag" reads, they are splitted and counted (`stat_snp`, **Step4**) according to the base species at a particular SNP site. Besides the splitted fastq files, sikiclass also generates tables summarizing the sub-class fractions and INDEL distributions.

Refer to [Usage](#usage) on how to run sikipipe, and [Output](#output) for detailed description of result folders and files.

<p align="center">
  <img src="docs/images/workflow.png" width="500" style="display: block; margin: 20px auto">
</p>

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

<!-- TODO nf-core: Describe the minimum required steps to execute the pipeline, e.g. how to prepare samplesheets.
     Explain what rows and columns represent. For instance (please edit as appropriate):

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
```

Each row represents a fastq file (single-end) or a pair of fastq files (paired end).

-->

Now, you can run the pipeline using:

<!-- TODO nf-core: update the following command to include all required parameters for a minimal example -->

```bash
git clone git@github.com:hukai916/sikiclass.git
cd sikiclass

# update conf/test_local.config by providing the correct samplesheet path, etc.
nextflow run main.nf -profile docker,arm -c conf/test_local.config
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
> see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

## Output

By default, all results are saved in the "./results" folder, as specified by the `outdir = "./results"` parameter in the master configuration file "nextflow.config". Results from different steps are organized into corresponding subdirectories (e.g. "./results/01_classify_tag", "./results/02_classify_single_tag", *etc.*), the "./results/00_stat" folder contains summary statistical tables.

Nextflow implements a caching mechanism that stores all intermediate and final results in the "./work/" directory. By default, files in the "./results/" are symbolic links to that in the "./work/". To switch from `symlink` to `copy`, use `publish_dir_mode = copy` argument. Below summarizes the main contents of each result folder.

- Subfolder: 00_stat  
This directory contains summary tables generated using results from step 1-4.

* `fq_class_ratios.tsv`: Count fractions for each read class.
* `precise_tag_snp_fraction.tsv`: Count fractions for different base species at given SNP site for "precise tag" reads.
* `no_tag_indel_size_distribution.tsv`: INDEL size distribution for "no tag" reads.
* `no_tag_indel_size_location_to_pam`: INDEL size and location to PAM site for "no tag" reads for each sample.

- Subfolder: 01_classify_tag  
This directory contains fastq files falling under each class determined by tag occurrence.

* `01a_no_tag`: Fastq file for reads with "no tag" for each sample.
* `01b_single_tag`: Fastq file for reads with "single tag" for each sample.
* `01c_multiple_tag`: Fastq file for reads with "multiple tag" for each sample.
* `01d_any_tag`: Fastq file for reads with "any tag" for each sample (union of 01b and 01c).
* `01e_tmp_fasta`: Intermediate fasta file.
* `01f_tmp_bam`: Intermedaite BAM (bwa) file.

- Subfolder: 02_classify_single_tag  
This directory contains fastq files falling under each class determined by INDEL occurrence for "single tag" reads.

* `02a_precise_tag`: Fastq file for reads with "precise tag" for each sample.
* `02a_precise_tag_snp_wt`: Fastq file for reads with "precise tag" that also contain WT base at given SNP site for each sample.
* `02a_precise_tag_snp_mut`: Fastq file for reads with "precise tag" that also contain mutant base at given SNP site for each sample.
* `02b_5indel`: Fastq file for reads with "INDELs in 5' region" for each sample.
* `02c_3indel`: Fastq file for reads with "INDELs in 3' region" for each sample.
* `02d_any_indel`: Fastq file for reads with "INDELs in any region" for each sample.
* `02e_tmp_bam`: Intermedaite BAM (minimap2) file.
* `02f_tmp_indel_pos`: Intermedaite INDEL size and location file for each sample.

- Subfolder: 03_classify_no_tag  
This directory contains fastq files falling under each class determined by INDEL occurrence for "no tag" reads.

* `03a_indel`: Fastq file for reads with "no tag" for each sample.
* `03b_deletion`: Fastq file for reads with "deletions in any region" for each sample.
* `03c_insertion`: Fastq file for reads with "insertions in any region" for each sample.
* `03d_tmp_bam`: Intermedaite BAM (minimap2) file.
* `03e_tmp_indel_pos`: Intermedaite INDEL size and location file for each sample.

## Credits

sikiclass was originally designed and written by Kai Hu, Nathan Lawson, and Julie Zhu.

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use hukai916/sikiclass for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
