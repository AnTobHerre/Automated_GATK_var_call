#!/bin/bash
#SBATCH --job-name=variant_calling_pipeline
#SBATCH --error=pipeline_master.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=400G
#SBATCH --qos=long
#SBATCH --time=10-00:00:00

module load gcc12-env/12.1.0
module load miniconda3/4.12.0

# 1. BWA MEM + SAMTOOLS SORT (parallel per sample)
bwa_jobids=()
input_dir="~/Zymoproj/DEsamples/fastq"
output_dir="~/Zymoproj/DEsamples/BAM"
ref="ref/Zymoseptoria_tritici.MG2.dna.toplevel.mt+.fa"
mkdir -p "\$output_dir"
for forward_reads in "$input_dir"/*R1_001.fastq.gz; do
    sample_name=$(basename "$forward_reads" _R1_001.fastq.gz)
    reverse_reads="$input_dir/${sample_name}_R2_001.fastq.gz"
    output_sam="$output_dir/${sample_name}.sam"
    sorted_bam="$output_dir/${sample_name}.bam"

    jid=$(sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=bwa_${sample_name}
#SBATCH --error=bwa_${sample_name}.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=100G
#SBATCH --qos=long
#SBATCH --time=12:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=a.tobherre@phytomed.uni-kiel.de

module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate samtools

input_dir="~/Zymoproj/DEsamples/fastq"
output_dir="~/Zymoproj/DEsamples/BAM"
reference="ref/Zymoseptoria_tritici.MG2.dna.toplevel.mt+.fa"

bwa mem -t 32 "\$reference" "\$input_dir/${sample_name}_R1_001.fastq.gz" "\$input_dir/${sample_name}_R2_001.fastq.gz" > "\$output_dir/${sample_name}.sam"
samtools sort -@ 8 -o "\$output_dir/${sample_name}.bam" "\$output_dir/${sample_name}.sam"
rm "\$output_dir/${sample_name}.sam"
EOF
)
    bwa_jobids+=($(echo $jid | awk '{print $4}'))
done
# Wait for all BWA jobs to finish before starting the next step
jid1=$(IFS=:; echo "${bwa_jobids[*]}")
# 2. PICARD ADD OR REPLACE READ GROUPS
jid2=$(sbatch --dependency=afterok:$jid1 <<'EOF'
#!/bin/bash
set -euo pipefail
#SBATCH --job-name=AddRGDE
#SBATCH --error=AddRGDE.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=400G
#SBATCH --qos=long
#SBATCH --time=10-00:00:00

module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate picard

BamFolder="~/Zymoproj/DEsamples/BAM/"
BamRGFolder="~/Zymoproj/DEsamples/BAM/RG"
project="Zt_DE_2024"
mkdir -p "$BamRGFolder"

for bamfile in $BamFolder/*.bam; do
    sample=$(basename "$bamfile" .bam)
    picard AddOrReplaceReadGroups \
        I="$bamfile" \
        O="$BamRGFolder/$sample.bam" \
        SORT_ORDER=coordinate \
        RGID="$sample" \
        RGLB=Zt_DE_2024 \
        RGPL=illumina \
        RGSM="$sample" \
        RGPU="$project" \
        CREATE_INDEX=True
done
EOF
)
jid2=$(echo $jid2 | awk '{print $4}')
# 3. PICARD MARK DUPLICATES
jid3=$(sbatch --dependency=afterok:$jid2 <<'EOF'
#!/bin/bash
set -euo pipefail
#SBATCH --job-name=mark_duplicate_De
#SBATCH --error=mark_duplicates_DE.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=400G
#SBATCH --qos=long
#SBATCH --time=10-00:00:00

module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate picard

BamRGFolder="~/Zymoproj/DEsamples/BAM/RG"
BamMDFolder="~/Zymoproj/DEsamples/BAM/MD"
mkdir -p "$BamMDFolder"

for bamFile in $BamRGFolder/*.bam; do
    Sample=$(basename "$bamFile" .bam)
    picard MarkDuplicates \
        INPUT="$bamFile" \
        OUTPUT="$BamMDFolder/$Sample.DuplMark.bam" \
        METRICS_FILE="$BamMDFolder/metrics_$Sample.txt" \
        CREATE_INDEX=true
done
EOF
)
jid3=$(echo $jid3 | awk '{print $4}')
# 4. SAMTOOLS IDXSTATS & DEPTH
jid4=$(sbatch --dependency=afterok:$jid3 <<'EOF'
#!/bin/bash
set -euo pipefail
#SBATCH --job-name=checkmaprtDE
#SBATCH --error=checkmaprtDE.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=50G
#SBATCH --time=04:00:00

module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate samtools

bam_dir="~/Zymoproj/DEsamples/BAM/MD"
output_dir="~/Zymoproj/DEsamples/covstats"
mkdir -p "$output_dir"

for bamfile in "$bam_dir"/*.DuplMark.bam; do
    filename=$(basename "$bamfile" .bam)
    output_file="$output_dir/$filename.perchrmap.txt"
    samtools idxstats "$bamfile" | awk -v filename="$filename" '{total[$1]+=$3+$4; mapped[$1]+=$3} END {print "Sample\tChromosome\tMapped Reads\tTotal Reads\tMapping Rate"; for(i in total) printf "%s\t%s\t%d\t%d\t%.2f%%\n", filename, i, mapped[i], total[i], (mapped[i]/total[i])*100;}' > "$output_file"
done
cat "$output_dir"/*.perchrmap.txt > "$output_dir/all_per_chr_mapping_rateDE.txt"

for bam_file in "$bam_dir"/*.DuplMark.bam; do
    filename=$(basename "$bam_file" .bam)
    output_file="$output_dir/$filename.coverage.txt"
    samtools depth -a "$bam_file" | awk -v filename="$filename" '{sum[$1]+=$3; count[$1]++} END {for (chr in sum) print filename, chr, sum[chr]/count[chr]}' > "$output_file"
done
cat "$output_dir"/*.coverage.txt > "$output_dir/all_per_chr_coverageDE.txt"
# Filter BAMs: only keep samples where all core chromosomes 1-13 have mean depth >= 10
good_bams="$output_dir/good_bams.txt"
> "$good_bams"
for covfile in "$output_dir"/*.coverage.txt; do
    fail=0
    while read -r sample chr depth; do
        # Only check chromosomes 1-13
        if [[ "$chr" =~ ^([1-9]|1[0-3])$ ]]; then
            if (( $(echo "$depth < 10" | bc -l) )); then
                fail=1
                break
            fi
        fi
    done < "$covfile"
    if [[ $fail -eq 0 ]]; then
        # Only add sample if all core chromosomes have depth >= 10
        sample=$(basename "$covfile" .coverage.txt)
        echo "$sample" >> "$good_bams"
    fi
done
EOF
)
jid4=$(echo $jid4 | awk '{print $4}')
# 5. GATK HAPLOTYPECALLER (per sample, only for good BAMs)
jid5=$(sbatch --dependency=afterok:$jid4 <<'EOF'
#!/bin/bash
set -euo pipefail
#SBATCH --job-name=HaplotypeCaller_DE
#SBATCH --error=HaplotypeCaller_DE.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=08:00:00

module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate var_call

REF=ref/Zymoseptoria_tritici.MG2.dna.toplevel.mt+.fa
BAMFolder="~/Zymoproj/DEsamples/BAM/MD"
gVCFFolder="~/Zymoproj/DEsamples/gVCF/"
covstats="~/Zymoproj/DEsamples/covstats"
mkdir -p "$gVCFFolder"

while read sample; do
    bamfile="$BAMFolder/${sample}.DuplMark.bam"
    if [[ -f "$bamfile" ]]; then
        gatk --java-options "-Xmx400G" HaplotypeCaller -R Zymoseptoria_tritici.MG2.dna.toplevel.mt+.fa -ploidy 1 --emit-ref-confidence GVCF -I "$bamfile" -O "$gVCFFolder/${sample}.g.vcf"
    fi
done < "$covstats/good_bams.txt"
EOF
)
jid5=$(echo $jid5 | awk '{print $4}')
# 5b. SORT INDIVIDUAL gVCF FILES PER FIELD
jid5b=$(sbatch --dependency=afterok:$jid5 <<'EOF'
#!/bin/bash
set -euo pipefail
#SBATCH --job-name=sort_gvcfs
#SBATCH --error=sort_gvcfs.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --time=04:00:00

module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate var_call

gVCFFolder="~/Zymoproj/DEsamples/gVCF"
ref="~/ref/Zymoseptoria_tritici.MG2.dna.toplevel.mt+.fa"
sites=("DO" "FU" "KA" "KO" "RA")

for site in "${sites[@]}"; do
    for gvcf in $gVCFFolder/${site}*.g.vcf; do
        [ -e "$gvcf" ] || continue
        sorted_gvcf="${gvcf%.g.vcf}_sorted.g.vcf"
        gatk SortVcf \
            -I "$gvcf" \
            -O "$sorted_gvcf" \
            --SEQUENCE_DICTIONARY "$ref".dict
    done
done
EOF
)
jid5b=$(echo $jid5b | awk '{print $4}')
# 6. GATK COMBINEGVCFS (per field, using sorted gVCFs)
jid6=$(sbatch --dependency=afterok:$jid5b <<'EOF'
#!/bin/bash
set -euo pipefail
#SBATCH --job-name=combine_gvcfs
#SBATCH --error=combine_gvcfs.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=08:00:00

module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate var_call

gVCFFolder="~/Zymoproj/DEsamples/gVCF"
ref="~/ref/Zymoseptoria_tritici.MG2.dna.toplevel.mt+.fa"
sites=("DO" "FU" "KA" "KO" "RA")

for site in "${sites[@]}"; do
    sorted_gvcfs=($gVCFFolder/${site}*_sorted.g.vcf)
    if [ ${#sorted_gvcfs[@]} -eq 0 ]; then
        echo "No sorted gVCFs found for $site"
        continue
    fi
    gatk CombineGVCFs \
        -R "$ref" \
        $(for gvcf in "${sorted_gvcfs[@]}"; do echo -n "-V $gvcf "; done) \
        -O "$gVCFFolder/${site}.g.vcf"
done
EOF
)
jid6=$(echo $jid6 | awk '{print $4}')
# 6b. SORT MERGED gVCF OUTPUTS
jid6b=$(sbatch --dependency=afterok:$jid6 <<'EOF'
#!/bin/bash
set -euo pipefail
#SBATCH --job-name=sort_merged_gvcfs
#SBATCH --error=sort_merged_gvcfs.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --time=04:00:00

module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate var_call

gVCFFolder="~/Zymoproj/DEsamples/gVCF"
ref="~/ref/Zymoseptoria_tritici.MG2.dna.toplevel.mt+.fa"
sites=("DO" "FU" "KA" "KO" "RA")

for site in "${sites[@]}"; do
    merged_gvcf="$gVCFFolder/${site}.g.vcf"
    sorted_merged_gvcf="$gVCFFolder/${site}_merged_sorted.g.vcf"
    gatk SortVcf \
        -I "$merged_gvcf" \
        -O "$sorted_merged_gvcf" \
        --SEQUENCE_DICTIONARY "$ref".dict
done
EOF
)
jid6b=$(echo $jid6b | awk '{print $4}')
# 7. GATK GENOTYPEGVCFS (per field, using sorted merged gVCF)
jid7=$(sbatch --dependency=afterok:$jid6b <<'EOF'
#!/bin/bash
set -euo pipefail
#SBATCH --job-name=genotype_DE_field
#SBATCH --error=genotype_DE_field.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=300G
#SBATCH --qos=long
#SBATCH --time=4-00:00:00

module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate var_call

vcf_path="~/Zymoproj/DEsamples/gVCF"
ref="~/ref/Zymoseptoria_tritici.MG2.dna.toplevel.mt+.fa"
sites=("DO" "FU" "KA" "KO" "RA")

for site in "${sites[@]}"; do
    sorted_merged_gvcf="$vcf_path/${site}_merged_sorted.g.vcf"
    output_vcf="$vcf_path/${site}_genotyped.vcf"
    if [ ! -f "$sorted_merged_gvcf" ]; then
        echo "Sorted merged gVCF file for site $site does not exist. Skipping..."
        continue
    fi
    gatk --java-options "-Xmx400G" GenotypeGVCFs \
       -V "$sorted_merged_gvcf" \
       -O "$output_vcf" \
       -R "$ref"
    echo "Genotyping complete for site $site."
done
EOF
)
jid7=$(echo $jid7 | awk '{print $4}')
# 8. GATK SELECTVARIANTS (per field)
jid8=$(sbatch --dependency=afterok:$jid7 <<'EOF'
#!/bin/bash
set -euo pipefail
#SBATCH --job-name=selectvarDE_field
#SBATCH --error=selectvarDE_field.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=300G
#SBATCH --time=12:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=a.tobherre@phytomed.uni-kiel.de

module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate var_call

gvcf_path="~/Zymoproj/DEsamples/gVCF"
vcf_path="~/Zymoproj/DEsamples/VCF"
ref="~/ref/Zymoseptoria_tritici.MG2.dna.toplevel.mt+.fa"
mkdir -p "$vcf_path"
for vcf_file in $gvcf_path/*genotyped.vcf; do
    base_name=$(basename "$vcf_file" .genotyped.vcf)
    gatk --java-options "-Xmx300G -XX:+UseParallelGC -XX:ParallelGCThreads=32" SelectVariants \
        -R $ref \
        -V "$vcf_file" \
        --select-type-to-include SNP \
        -O "$vcf_path/${base_name}.SNP.vcf"
done
EOF
)
jid8=$(echo $jid8 | awk '{print $4}')
# 9. GATK VARIANTFILTRATION (per field)
jid9=$(sbatch --dependency=afterok:$jid8 <<'EOF'
#!/bin/bash
set -euo pipefail
#SBATCH --job-name=Varfilt_DE_field
#SBATCH --error=Varfilt_DE_field.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=200G
#SBATCH --time=12:00:00

module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate var_call

vcf_path="~/Zymoproj/DEsamples/VCF"
ref="~/ref/Zymoseptoria_tritici.MG2.dna.toplevel.mt+.fa"

for snp_vcf in $vcf_path/*.SNP.vcf; do
    base_name=$(basename "$snp_vcf" .SNP.vcf)
    gatk --java-options "-Xmx200G -XX:+UseParallelGC -XX:ParallelGCThreads=32" VariantFiltration \
        -R $ref \
        -V "$snp_vcf" \
        --filter-name "MQ" --filter-expression "MQ < 20.0" \
        --filter-name "QDFilter" --filter-expression "QD < 5.0" \
        --filter-name "QUAL30" --filter-expression "QUAL < 30.0" \
        --filter-name "SOR" --filter-expression "SOR > 3.0" \
        --filter-name "FS" --filter-expression "FS > 60.0" \
        --filter-name "Low_depth3" --filter-expression "DP < 10" \
        --filter-name "ReadPosRankSum_lower" --filter-expression "ReadPosRankSum < -2.0" \
        --filter-name "ReadPosRankSum_upper" --filter-expression "ReadPosRankSum > 2.0" \
        --filter-name "MQRankSum_lower" --filter-expression "MQRankSum < -2.0" \
        --filter-name "MQRankSum_upper" --filter-expression "MQRankSum > 2.0" \
        -O "$vcf_path/${base_name}.SNP.corrected.qualityfilter.vcf"
done
EOF
)
jid9=$(echo $jid9 | awk '{print $4}')
# 10. GATK SELECTVARIANTS (exclude filtered)
jid10=$(sbatch --dependency=afterok:$jid9 <<'EOF'
#!/bin/bash
set -euo pipefail
#SBATCH --job-name=deVCFex_field
#SBATCH --error=deVCFex_field.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=50G
#SBATCH --time=04:00:00

module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate var_call

vcf_path="~/Zymoproj/DEsamples/VCF"
ref="~/ref/Zymoseptoria_tritici.MG2.dna.toplevel.mt+.fa"

for filtered_vcf in $vcf_path/*.SNP.corrected.qualityfilter.vcf; do
    base_name=$(basename "$filtered_vcf" .SNP.corrected.qualityfilter.vcf)
    gatk SelectVariants --java-options "-Xmx50G -XX:+UseParallelGC -XX:ParallelGCThreads=32" \
        -R $ref \
        -V "$filtered_vcf" \
        -O "$vcf_path/${base_name}.SNP.corrected.qualityfilter.excl.vcf" \
        --exclude-filtered true
done
EOF
)
jid10=$(echo $jid10 | awk '{print $4}')
# 11. VCFTOOLS/PLINK/BCFTOOLS MERGE
jid11=$(sbatch --dependency=afterok:$jid10 <<'EOF'
#!/bin/bash
set -euo pipefail
#SBATCH --job-name=within_field_filter
#SBATCH --error=within_field_filter.err
#SBATCH --nodes=1
#SBATCH --mem=100G
#SBATCH --time=06:00:00

module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate samtools

folder="~/Zymoproj/DEsamples/VCF"
german_fields=("DO" "FU" "KA" "KO" "RA")

vcf_files=(
    "$folder/Zt_CH.sub.qualityfilter2021.excl.vcf"
    "$folder/Zt_DO.qualityfilter2021.excl.vcf"
    "$folder/Zt_FU.qualityfilter2021.excl.vcf"
    "$folder/Zt_KA.qualityfilter2021.excl.vcf"
    "$folder/Zt_KO.qualityfilter2021.excl.vcf"
    "$folder/Zt_RA.qualityfilter2021.excl.vcf"
    "$folder/Zt_UK.sub.qualityfilter2021.excl.vcf"
    "$folder/Zt_US.sub.qualityfilter2021.excl.vcf"
)

for vcf in "${vcf_files[@]}"; do
    filename_no_ext=$(basename "$vcf" .vcf)
    filename_base="${filename_no_ext}"

    vcftools --vcf "$vcf" \
      --recode --recode-INFO-all \
      --max-missing 0.8 \
      --mac 1 \
      --remove-filtered-all --remove-indels \
      --min-alleles 2 --max-alleles 2 \
      --out "${folder}/${filename_base}.max-m-80.biallelic-only.mac1"

    field_code=$(echo "$filename_base" | sed -E 's/^Zt_([A-Z]{2}).*/\1/')
    sub_part=""
    if [[ "$filename_base" == *".sub"* ]]; then
        sub_part="_sub"
    fi

    if [[ " ${german_fields[@]} " =~ " ${field_code} " ]]; then
        out_prefix="${folder}/Zt_${field_code}_strict"
    else
        out_prefix="${folder}/Zt_${field_code}${sub_part}_strict"
    fi

    plink --vcf "${folder}/${filename_base}.max-m-80.biallelic-only.mac1.recode.vcf" \
      --make-bed --freq --missing \
      --out "$out_prefix"
done

for vcf in "${folder}"/*.max-m-80.biallelic-only.mac1.recode.vcf; do
    bgzip -f "$vcf"
    tabix -p vcf "${vcf}.gz"
done

vcf_list=(${folder}/*.max-m-80.biallelic-only.mac1.recode.vcf.gz)
bcftools merge "${vcf_list[@]}" -Oz -o "${folder}/Zt_WW_all_strict.vcf.gz"
tabix -p vcf "${folder}/Zt_WW_all_strict.vcf.gz"

grep -v -E '^(ST16|ORE15)' $folder/sample_ids_WWstrict.txt > $folder/ids_EU_strict.txt
bcftools view \
  -S $folder/ids_EU_strict.txt \
  -i 'AC>1' \
  $folder/Zt_WW_all_strict.vcf.gz \
  -O z -o $folder/Zt_EU_all_strict.vcf.gz && bcftools index -f $folder/Zt_EU_all_strict.vcf.gz

grep -v -E '^(ST16|ORE15|Rose|Rag|ADAS)' $folder/sample_ids_WWstrict.txt > $folder/ids_DE_strict.txt
bcftools view \
  -S $folder/ids_DE_strict.txt \
  -i 'AC>1' \
  $folder/Zt_WW_all_strict.vcf.gz \
  -O z -o $folder/Zt_DE_all_strict.vcf.gz && bcftools index -f $folder/Zt_DE_all_strict.vcf.gz
EOF
)
jid11=$(echo $jid11 | awk '{print $4}')
