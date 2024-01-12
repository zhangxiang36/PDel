#!/usr/bin/env bash
# Used to aanalyze outcome of gene-editing outcome by calling indel/SNP with NGS data
baseP=$1 # depending on data path in your own computer

cd ${baseP}

ls -l | grep ^d | awk '{FS=" "}{print $9}' | while read dir
do
    echo going in ${dir}
    cd ${dir}
    ls *fq.gz | sed 's#.R[12].fq.gz$##g' | sort | uniq | while read file
    do
        echo ${file}
        trimmomatic SE \
            -threads 8 -phred33 \
            ${file}.R1.fq.gz \
            ${file}.trim.R1.fq.gz \
            ILLUMINACLIP:/path/to/adapter.fa:2:30:10:8:true \
            SLIDINGWINDOW:5:20 \
            LEADING:3 \
            TRAILING:3 \
            MINLEN:36
        fastqc ${file}.trim.R1.fq.gz
    done
    echo leaving ${dir}
    cd ../
done


ls -l | grep ^d | awk '{FS=" "}{print $9}' | while read dir
do
    cd ${dir}
    # create dict file for GATK HaplotypeCaller
    picard CreateSequenceDictionary \
        R=Ref.txt \
        O=Ref.dict
    samtools faidx Ref.txt
    bwa index Ref.txt
    #subdir=analysis
    tmpdir=SNPindel
    #mkdir ${subdir}
    mkdir ${tmpdir}
    prefix=""
    ls *trim.R1.fq.gz | while read file
    do
        file_tmp=${file##${prefix}}
        file_tmp=${file_tmp%.fq.gz}
        bwa mem -R "@RG\tID:1\tLB:library\tPL:Illumina\tSM:${file_tmp%.R1}\tPU:machine" Ref.txt $file > ${tmpdir}/${file_tmp}.bam
        samtools view -F 4 -S -u ${tmpdir}/${file_tmp}.bam > ${tmpdir}/${file_tmp}.accept_hits
        samtools sort -o ${tmpdir}/${file_tmp}.accept_hits.bam ${tmpdir}/${file_tmp}.accept_hits
        picard MarkDuplicates \
            I=${tmpdir}/${file_tmp}.accept_hits.bam \
            O=${tmpdir}/${file_tmp}.markdup.sorted.bam \
            M=${tmpdir}/${file_tmp}.markdup_metrics.txt3
        samtools view ${tmpdir}/${file_tmp}.markdup.sorted.bam > ${tmpdir}/${file_tmp}.markdup.sorted.sam
        samtools index ${tmpdir}/${file_tmp}.markdup.sorted.bam
        gatk HaplotypeCaller \
            --emit-ref-confidence GVCF \
            -R Ref.txt \
            -I ${tmpdir}/${file_tmp}.markdup.sorted.bam \
            -O ${tmpdir}/${file_tmp}.g.vcf
        gatk GenotypeGVCFs \
            --output ${tmpdir}/${file_tmp}.raw.vcf \
            --variant ${tmpdir}/${file_tmp}.g.vcf \
            --reference Ref.txt
        gatk SelectVariants \
            --output ${tmpdir}/${file_tmp}.raw.snp.vcf \
            --select-type-to-include SNP \
            --variant ${tmpdir}/${file_tmp}.raw.vcf
        gatk SelectVariants \
            --output ${tmpdir}/${file_tmp}.raw.indel.vcf \
            --select-type-to-include INDEL \
            --variant ${tmpdir}/${file_tmp}.raw.vcf
        gatk VariantFiltration \
            -R Ref.txt \
            -V ${tmpdir}/${file_tmp}.raw.snp.vcf \
            -O ${tmpdir}/${file_tmp}.snp.vcf \
            -filter-name "QUAL_filter" -filter "QUAL < 30.0" \
            -filter-name "QD_filter" -filter "QD < 2.0" \
            -filter-name "FS_filter" -filter "FS > 30.0"
        gatk VariantFiltration \
            -R Ref.txt \
            -V ${tmpdir}/${file_tmp}.raw.indel.vcf \
            -O ${tmpdir}/${file_tmp}.indel.vcf \
            -filter-name "QUAL_filter" -filter "QUAL < 30.0" \
            -filter-name "QD_filter" -filter "QD < 2.0" \
            -filter-name "FS_filter" -filter "FS > 30.0"
    done
    rm ${tmpdir}/*raw*
    echo ${tmpdir}/result
    mkdir ${tmpdir}/result
    cp ${tmpdir}/*indel.vcf ${tmpdir}/result
    cp ${tmpdir}/*snp.vcf ${tmpdir}/result
    ls ${tmpdir}/result
    ls ${tmpdir}/result/*indel.vcf | while read file
    do
        prefix=${file%.indel.vcf}
        echo $prefix doing
        echo $prefix doing
        (cat $prefix.indel.vcf | grep -v "#" | awk '{FS="\t"}{print "INDEL\t"$1"\t"$2"\t"$4"\t"$5"\t"$6"\t"$7"\t"$10}';\
            cat $prefix.snp.vcf | grep -v "#" | awk '{FS="\t"}{print "SNP\t"$1"\t"$2"\t"$4"\t"$5"\t"$6"\t"$7"\t"$10}')  > \
            $prefix.tmp.txt
        awk -F ':' 'BEGIN{print "TYPE\tNAME\tPOS\tREF\tALT\tQUAL\tFILTER\tHetorHomo\tAlleleDepth"}{print $1"\t"$2}'  $prefix.tmp.txt | \
        sed 's#|#/#g'> $prefix.variation.txt
    done
    rm ${tmpdir}/result/*tmp*
    rm ${tmpdir}/result/*snp.vcf
    rm ${tmpdir}/result/*indel.vcf
    cd ../
done

ls -l | grep ^d | awk '{FS=" "}{print $9}' | while read dir
do
    cd ${dir}/SNPindel/result
    outfile=${dir}_all_sample_combined.txt
    echo -e "SAMPLE\tTYPE\tNAME\tPOS\tREF\tALT\tQUAL\tFILTER\tHetorHomo\tAlleleDepth" > $outfile
    ls *variation.txt | while read file
    do
        samplename=`echo $file | sed 's#\.trim\..\+##g'`
        linen=`cat $file | grep -v "^TYPE" | wc -l`
        if [ $linen -gt 0 ]
        then
            cat $file | grep -v "^TYPE" | awk -F "\t" '{print samplename"\t"$0}' samplename="$samplename" >> ${outfile}
        else
            echo -e "$samplename\tNo_mutation_found" >> $outfile
        fi
    done
    cp ${outfile} ../../../
    cd ../../../
done
