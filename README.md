# Computational analyses for the *Pusa saimensis* study

## Data mapping and variant calling

The short-read data were mapped using 'bwa mem' ([v.2.2.1](https://github.com/bwa-mem2/bwa-mem2)) on the Saimaa ringed seal genome (Kammonen et al., in preparation). The reads pairs were marked and data sorted using 'samtools fixmate' and 'samtools sort' ([v.1.18](https://github.com/samtools/samtools)) and realigned using GATK3 'IndelRealigner' ([v.3.8](https://console.cloud.google.com/storage/browser/gatk-software/package-archive/gatk)). Variants were called by applying GATK4 'HaplotypeCaller' ([v.4.2.5](https://github.com/broadinstitute/gatk/releases)) for each sample, joining the data using 'CombineGVCFs', and applying 'GenotypeGVCFs' to the joint data. 
    
    bwa-mem2 mem -R ${rgroup} -K 100000000 -Y ${refname} <(cat ${fastq1}) <(cat ${fastq2}) \
    | samtools view -h - | samtools fixmate -m - -  | samtools sort - \
    | samtools view -h -O bam -o ${smp_id}_bwa.bam && samtools index ${smp_id}_bwa.bam

    gatk3 -T RealignerTargetCreator -R ${refname} -I ${bwabam} -o targets.intervals
    
    gatk3 -T IndelRealigner -R ${refname} -I ${bwabam} -targetIntervals targets.intervals \
     -o ${smp_id}_real.bam

    samtools markdup ${smp_id}_real.bam ${smp_id}_mdup.bam && samtools index "${smp_id}_mdup.bam"

    gatk4 HaplotypeCaller -R ${refname} -I ${mdupbam} -O ${smp_id}.gvcf.gz -ERC GVCF


    gatk4 CombineGVCFs -R ${refname} -V ${smp_id1}.gvcf.gz .. -V ${smp_idN}.gvcf.gz -O ${out}.gvcf.gz
    
    gatk4 GenotypeGVCFs -R ${refname} -V ${out}.gvcf.gz -O ${out}.vcf.gz

The variant calling was performed on a dataset including 46 ringed seal individuals selected based on sequencing coverage and data completeness, and four spotted seals.

## Mitochondrial analyses

For the mitochondrial analyses, the corresponding BAM region was extracted and the GATK analysis was repeated with `-ploidy 1`.

    samtools view ${sample}.bam mt_26042022 -Obam -o ${sample}_mt.bam

    gatk HaplotypeCaller -R ${refname} -I ${sample}_mt.bam -O ${sample}_mt.gvcf.gz \
     -ploidy 1 -L mt_26042022 -ERC GVCF
    
    gatk CombineGVCFs -R ${refname} -V ${smp_id1}_mt.gvcf.gz .. -V ${smp_idN}_mt.gvcf.gz -O seals305.gvcf.gz
    
    gatk GenotypeGVCFs -R ${refname} -V seals305.gvcf.gz -O seals305.vcf.gz


The protein-coding regions of ['NC_008433.1'](https://www.ncbi.nlm.nih.gov/nuccore/NC_008433.1/), mitochondrion of *Pusa hispida*, were extracted as translated peptides using BioPython.

    import sys
    from Bio import SeqIO

    for gb_record in SeqIO.parse(open("NC_008433.gb","r"), "genbank") :
        for (index, feature) in enumerate(gb_record.features) :
            if(feature.type=="CDS"):
                print(">"+feature.qualifiers['gene'][0]+"\n"
                  +feature.qualifiers['translation'][0])

These were then mapped against the nucleotide sequence using 'miniprot' ([v.0.5-r179](https://github.com/lh3/miniprot)), creating a GFF-formatted annotation. This was converted into a TAB-separated list of regions.

    miniprot mt_26042022.fa seal_translations.pep --gff > mt_26042022.gff

    grep mRNA mt_26042022.gff | awk '{OFS="\t";print $1,$4,$5}' > mt_26042022.mRNA.list

Using 'bcftools' ([v.1.18](https://github.com/samtools/bcftools)), the binary SNPs with the protein-coding regions of the VCF file were extracted. Including only the 166 unique samples (identified with a RAxML analysis of the full data), the variable positions were converted to TAB-separated format and this then to FASTA using [scripts/vcf_tab_to_fasta_alignment.pl](scripts/vcf_tab_to_fasta_alignment.pl).


    bcftools view -m2 -M2 -v snps -R mt_26042022.mRNA.list seals305.vcf.gz \
        -Oz -o seals305.mt_regs.vcf.gz 
        
    bcftools view -S 166_unique.txt seals305.mt_regs.vcf.gz | \
        bcftools view -q 0.0001:minor  | vcf-to-tab > seals166.mt_regs.vcf.tab
        
    perl vcf_tab_to_fasta_alignment.pl seals166.mt_regs.vcf.tab > seals166.mt_regs.vcf.fasta

The maximum likelihood phylogeny was then searched with RAxML ([v.8.2.12](https://github.com/stamatak/standard-RAxML)), correcting for the missing data.

    raxmlHPC -s $fasta -m ASC_GTRGAMMA -p 123 -x 123 -f a -# 1000 -n $output --asc-corr=felsenstein -q part

where `part` was

    [asc~p1.txt], ASC_DNA, p1=1-1424
and `p1.txt` was

    9945

The [resulting tree](figures/Fig_S1_full.pdf) was plotted using the 'ggtree' R package.


## Variant selection for genomic analyses

Variants were extracted using 'bcftools' and 'plink' ([v.1.9](https://www.cog-genomics.org/plink/)). Binary SNPs that a) varied in P. hispida or b) had MAF 0.05 in P. hispida, and had less than 10% missing data in the full set were retained. These were then thinned by requiring at least 1,000 bp inter-SNP distance. 

    bcftools view -T^RM.bed.gz -v snps -m2 -M2 $vcf_full | \
        bcftools view  -T posmask.bed.gz  -Oz -o $vcf_snps

    bcftools view -S Phis_samples.txt $vcf_snps | \
        bcftools view -q0.0001:minor | \
        bcftools query -f '%CHROM\t%POS\t%POS\n' > Phis_var.txt

    bcftools view -T Phis_var.txt -i 'F_MISSING<0.10' -Oz -o $vcf_var $vcf_snps

    bcftools view -S Phis_samples.txt $vcf_snps | \
        bcftools view -q0.05:minor | \
        bcftools query -f '%CHROM\t%POS\t%POS\n' > Phis_maf05.txt

    bcftools view -T Phis_maf05.txt -i 'F_MISSING<0.10' -Oz -o $vcf_maf05 $vcf_snps

    plink --vcf $vcf --allow-extra-chr --bp-space 1000 --recode vcf bgz --out $file


## Principal component analysis

The data were converted to BED format using 'plink' and then to EIGENSTRAT format using 'convertf' ([EIGENSOFT v.7.0.2](https://github.com/DReichLab/EIG/tree/master)). 

    plink --vcf ${file}.vcf.gz --allow-extra-chr --allow-no-sex --make-bed --out ${file}
    
    mv $file.bim $file.bim_
    awk '{OFS="\t";if($1!=pc){c=c+1;pc=$1}$1=c;print $0}' $file.bim_ > $file.bim

    cat > convert.txt << EOF
    genotypename:    $file.bed
    snpname:         $file.bim
    indivname:       $file.fam
    outputformat:    EIGENSTRAT
    genotypeoutname: $file.geno
    snpoutname:      $file.snp
    indivoutname:    $file.ind
    EOF
    
    convertf -p convert.txt 


The PCA was performed with 'smartpca' (EIGENSOFT v.7.0.2) setting the option 'usenorm: NO'.

    cat > smartpca.txt << EOF
    genotypename:    $file.geno
    snpname:         $file.snp
    indivname:       $file.ind
    evecoutname:     $file.pca.evec.txt
    evaloutname:     $file.pca.eval.txt
    usenorm:         NO
    EOF

    smartpca -p smartpca.txt  > smartpca.out

Results were plotted using the 'ggplot2' R package.


## Private alleles analysis

The SNP sharing analyses were performed using 'vcftools' ([v.0.1.17](https://vcftools.github.io/index.html)) and R. Using a bash script, allele frequencies were computed for sets of 1-5 individuals and the outputs were then combined into single files.

    cat > pops.txt << EOF
    Arctic
    Baltic
    Ladoga
    Okhotsk
    Saimaa
    EOF

    bcftools view -S Spotted.txt $vcf_maf05 | vcftools --vcf - --freq2  --out Spotted
    
    for n in {1..5}; do
    
        > Ringed${n}.txt
        cat pops.txt | while read pop; do
    qq        head -${n} $pop.txt > ${pop}${n}.txt
            cat ${pop}${n}.txt >> Ringed${n}.txt
        done
        
        cat pops.txt | while read pop; do
            bcftools view -S ${pop}${n}.txt $vcf_maf05 | vcftools --vcf - --freq2  --out ${pop}${n};
        done
        
        pop=Ringed
        bcftools view -S ${pop}${n}.txt $vcf_maf05 | vcftools --vcf - --freq2  --out ${pop}${n};
    
        paste Arctic${n}.frq Baltic${n}.frq Ladoga${n}.frq Okhotsk${n}.frq \
         Saimaa${n}.frq Ringed${n}.frq Spotted.frq | cut -f1,2,5,11,17,23,29,35,41 | \
         tail -n+2 > Seals${n}.frq
    done

These files were then analysed in R.


    #############################################

    library(dplyr)
    library(tidyr)
    library(ggplot2)

    ##
    
    get_data <- function(file) {
      cn <- c("chrom","pos","arc.ref","bal.ref","lad.ref","okh.ref","sai.ref","rin.ref","spo.ref")
      data <- read.table(file,col.names = cn)
      
      data <- data %>% mutate(AA=if_else(spo.ref==1,"REF",if_else(spo.ref==0,"ALT",NA)))
      
      data <- data %>% mutate(
        arc.drv=if_else(AA=="REF" & arc.ref<1,1, if_else(AA=="ALT" & arc.ref>0,1,0)),
        bal.drv=if_else(AA=="REF" & bal.ref<1,1, if_else(AA=="ALT" & bal.ref>0,1,0)),
        lad.drv=if_else(AA=="REF" & lad.ref<1,1, if_else(AA=="ALT" & lad.ref>0,1,0)),
        okh.drv=if_else(AA=="REF" & okh.ref<1,1, if_else(AA=="ALT" & okh.ref>0,1,0)),
        sai.drv=if_else(AA=="REF" & sai.ref<1,1, if_else(AA=="ALT" & sai.ref>0,1,0)),
        fixed=if_else(rin.ref==0|rin.ref==1,1,0)
      )
      
      data
    }
    
    ##

    get_derived <- function(data) {
      data.der <- data %>% 
        filter(!is.na(AA) & !fixed) %>% 
        filter(!is.na(arc.drv) & !is.na(bal.drv) & !is.na(lad.drv) & !is.na(okh.drv) & !is.na(sai.drv)) %>% 
        select(arc.drv, bal.drv, lad.drv, okh.drv, sai.drv) 
    
      colnames(data.der) <- c("arc","bal","lad","okh","sai")
      data.der <- data.der %>% mutate(uniq=rowSums(.[1:5])==1, sum=rowSums(.[1:5]))
      
      data.der
    }
    
    ##

    get_uniq <- function(data.der) {
      data.uniq <- data.der %>% select(-sum) %>% pivot_longer(!uniq,names_to="pop",values_to="derived") %>% 
        group_by(pop) %>% summarise(PT=sum(derived==1&uniq)/sum(derived==1))
      
      data.uniq
    }
    
    ##

    read_data <- function(path) {
      
      udat <- data.frame()
      for(i in 1:5) {
        data <- get_data(paste0(path,"/Seals",i,".frq"))
        derv <- get_derived(data)
        uniq <- get_uniq(derv)
        udat <- rbind.data.frame(
          udat,
          cbind(uniq,smp=i)
        )
      }
      udat$pop <- factor(levels=c('arc','okh','bal','lad','sai'),udat$pop)
      udat
    }

    ##

    plot_data <- function(udat,ymax,title="") {
      
      cols5 <- palette.colors(palette = "Okabe-Ito")[rev(c(6,5,4,8,7))]
      
      ggplot(udat) +
        geom_line(aes(smp,PT,color=pop),size=1)+
        geom_point(aes(smp,PT,color=pop),size=2)+
        theme_classic()+ylim(0,ymax)+
        scale_color_manual(values = cols5)+
        xlab("#Individuals") + ylab("#Private SNPs / #Total SNPs")+
        guides(color=guide_legend(title="",nrow = 1,override.aes=list(shape=15)))+
        theme(
          legend.key.width = unit(1,"mm"),
          text = element_text(size=9),
          axis.title = element_text(size=10),
          axis.text = element_text(size=8),
          legend.position="top",
          legend.text = element_text(margin = margin(r = 1, l=3, unit = "pt"),size=9),
          plot.title = element_text(hjust = 0.005, vjust = -11.5,size=9)
        )+ggtitle(title)
    }
    
    ##
    
    maf05 <- read_data("maf05")
    
    plot_data(maf05,0.25,"MAF=0.05")

    #############################################

## Demographic analyses

The demographic analyses were based on data from Löytynoja et al. (2023). The existing MSMC output files were reprocessed with [MSMC-IM](https://github.com/wangke16/MSMC-IM) (downloaded April 2024).

    python3 MSMC_IM.py -beta 1e-8,1e-6 -mu 1.826e-8 \
        -o msmc-im_out/${file} ccr_out/${file}.combined.msmc2.final.txt


These files were then visualised using R.

    #############################################
    
    library(ggplot2)
    library(dplyr)
    library(ggh4x)
    library(cowplot)
    library(scales)
    library(tidyr)
    library(RColorBrewer)

    ##
    
    read_data <- function(path,clip=FALSE) {
        dat <- c()
        for(file in dir(path,"*estimates.txt")){
            tmp <- read.table(paste0(path,file),head=T)
            if(clip == TRUE) {
                tmp <- tmp[-c(1,2),]
            }
            pair <- substr(file,1,11)
            set <- substr(file,1,5)
            sp1 <- substr(file,1,2)
            sp2 <- substr(file,4,5)
            tmp$pair <- pair
            tmp$set <- set
            tmp$sp1 <- sp1
            tmp$sp2 <- sp2
            tmp$right_time_boundary = c(tmp$left_time_boundary[2:length(tmp$left_time_boundary)],Inf)
            dat <- rbind.data.frame(dat,tmp)
        }
        dat$m_prime <- if_else(dat$M <= 0.999, dat$m, 1e-30)
        dat$Set = factor(dat$set,levels = c("ar.ba","ar.la","ba.la","sa.ar","sa.ba","sa.la"), 
            labels= c("Ar-Ba","Ar-La","Ba-La","Sa-Ar","Sa-Ba","Sa-La"))
        dat$Sp1 = factor(dat$sp1,levels = c("ar","ba","la","sa"), 
            labels= c("Arctic","Baltic","Ladoga","Saimaa"))
        dat$Sp2 = factor(dat$sp2,levels = c("ar","ba","la","sa"), 
            labels= c("Arctic","Baltic","Ladoga","Saimaa"))
        dat
    }

    ##
    
    dat.full <- read_data("msmc-im_out/")
    dat.clip <- read_data("msmc-im_out/",TRUE)

    ##
    
    cols <- palette.colors(palette = "Okabe-Ito")[c(7,4,5,6)]

    ggplot(dat.clip) + 
        geom_line(aes(10*left_time_boundary,im_N1,group=pair,color=Sp1)) +
        geom_line(aes(10*left_time_boundary,im_N2,group=pair,color=Sp2)) +
        xlab("Years ago") + ylab("Ne") +
        theme_classic()+
        theme(
            axis.text = element_text(size=8),
            legend.key.width = unit(3,"mm"),
            legend.position="top",
            legend.text = element_text(margin = margin(r = 3, l=1, unit = "pt"),size=9)
        )+
        scale_color_manual(values=cols,breaks=c("Arctic","Baltic","Ladoga","Saimaa"))+
            scale_x_log10(
            minor_breaks = c(seq(1e3,10e3,1e3),seq(1e4,10e4,1e4),seq(1e5,10e5,1e5),seq(1e6,10e6,1e6)),
            breaks=c(1e4,1e5,1e6),
            labels=c("10,000","100,000","1,000,000"),
            guide = "axis_minor", # this is added to the original code
            limits = c(5000,1.15e6)
        )+
        guides(color=guide_legend(title="",nrow = 1,override.aes=list(linewidth=3)))

    ##
    
    cols6 <- brewer.pal(10,'RdYlBu')[c(1,2,3,8,9,10)]
    
    ggplot(dat.full) + 
        geom_line(aes(10*right_time_boundary,M,group=pair,color=Set)) +
        geom_hline(yintercept = c(0,1),linetype="dashed") +
        geom_hline(yintercept = c(0.5,0.95,0.99),linetype="dotted") +
        coord_cartesian(xlim=c(1,30e4))+
        xlab("Years ago") +  ylab("M") + 
        scale_x_continuous(limits = c(10,50e4),
            breaks=c(0,5e4,1e5,2e5,3e5),
            labels=c("0","50,000","100,000","200,000","300,000"),
        )+
        theme_classic()+
        theme(
            axis.text = element_text(size=8),
            legend.key.width = unit(1,"mm"),
            legend.position="top",
            legend.text = element_text(margin = margin(r = -1, l=0, unit = "pt"),size=9)
        )+
        scale_color_manual(values=cols6[c(1,2,3,6,5,4)],
                           breaks=c("Ar-Ba","Ba-La","Ar-La","Sa-Ar","Sa-Ba","Sa-La"),
                           labels=c("Ar-Ba","Ba-La","Ar-La","Ar-Sa","Ba-Sa","La-Sa"))+
        guides(color=guide_legend(title="",nrow = 1,override.aes=list(linewidth=3)))

    ##
    
    getCDFintersect <- function(t, CDF, val){
        xVec <- t
        yVec <- CDF
        i <- 1
        CDFintersect <- NaN
        if(yVec[1] < val){
            while(yVec[i] < val) {
                i <- i + 1
            }
            if(! i > 0 | ! i <= length(yVec)){ 
                print("CDF intersection index out of bounds: ",i)
            }
            if(! yVec[i - 1] < val | ! yVec[i] >= val){
                print("this should never happen")
            }
            intersectDistance <- (val - yVec[i - 1]) / (yVec[i] - yVec[i - 1])
            CDFintersect = xVec[i - 1] + intersectDistance * (xVec[i] - xVec[i - 1])
        } else {
            CDFintersect <- val/yVec[1] * xVec[1]
        }
        return (CDFintersect)
    }

    ##
    
    splits <- dat.full %>% group_by(pair) %>% select(pair,set,left_time_boundary,M) %>% 
        group_by(pair,set) %>% 
        reframe(x01=getCDFintersect(left_time_boundary, M, 0.01)*10,
                x25=getCDFintersect(left_time_boundary, M, 0.25)*10,
                x50=getCDFintersect(left_time_boundary, M, 0.5)*10,
                x75=getCDFintersect(left_time_boundary, M, 0.75)*10,
                x99=getCDFintersect(left_time_boundary, M, 0.99)*10) %>%
        group_by(set) %>% 
        reframe(m01=round(mean(x01)),
                m25=round(mean(x25)),
                m50=round(mean(x50)),
                m75=round(mean(x75)),
                m99=round(mean(x99)))
    
    splits

    ##
    
    plot_one <- function(datf,cn,name) {
        datf2 <- datf %>% group_by(pair) %>% mutate(time=lead(time) - 0.01) %>% 
            mutate(time = replace_na(time, Inf)) %>% ungroup()
        
        datf3 <- rbind(datf, datf2)
        datf3 <- datf3[order(datf3$pair,datf3$time),]
        datf3 <- datf3 %>% mutate(time =  ifelse(time<0.001, 1, time)) 
        
        spl <- as.data.frame(splits %>% filter(set==datf$set[1]))
    
        ggplot(datf3) + 
            annotate("rect", xmin=spl$m01,xmax=spl$m99,ymin=0,ymax=0.00020,alpha=0.25, 
                fill = "darkgray") + 
            annotate("rect", xmin=spl$m25,xmax=spl$m75,ymin=0,ymax=0.00020,alpha=0.45, 
                fill = "darkgray") + 
            geom_line(aes(time,m_prime,group=pair),color=cols6[cn]) +
            geom_ribbon(aes(x=time,ymin=0, ymax = m_prime,group=pair),fill=cols6[cn],alpha=0.1) +
            coord_cartesian(ylim=c(0,1.9e-4),xlim=c(1000,1.15e6)) +
            xlab("Years ago") +  ylab("m_prime") + 
            scale_x_log10(
                minor_breaks = 
                    c(seq(1e3,10e3,1e3),seq(1e4,10e4,1e4),seq(1e5,10e5,1e5),seq(1e6,10e6,1e6)),
                breaks=c(1e3,1e4,1e5,1e6),
                labels=c("1,000","10,000","100,000","1,000,000"),
                guide = "axis_minor", # this is added to the original code
            )+
            theme_classic() + 
            theme(
                text = element_text(size=7),
                axis.text = element_text(size=7)
            ) +
            annotate("segment", x=spl$m50,xend=spl$m50,y=0,yend=0.00020,linetype="dashed",lwd=0.5) + 
            annotate("text", x=1.55e6,y=1.875e-4,label=name,size=3,hjust = 1, )
    }
    
    ##
    
    dat <- dat.full %>% select(left_time_boundary,m_prime,pair,set) %>% 
        mutate(time=10*left_time_boundary) 
    
    pl.1 <- plot_one(subset(dat,set=="ar.ba"),1,"Arctic-Baltic")
    pl.2 <- plot_one(subset(dat,set=="ar.la"),2,"Arctic-Ladoga")
    pl.3 <- plot_one(subset(dat,set=="ba.la"),3,"Baltic-Ladoga")
    
    pl.4 <- plot_one(subset(dat,set=="sa.ar"),6,"Arctic-Saimaa")
    pl.5 <- plot_one(subset(dat,set=="sa.ba"),5,"Baltic-Saimaa")
    pl.6 <- plot_one(subset(dat,set=="sa.la"),4,"Ladoga-Saimaa")
    
    plot_grid(pl.1,pl.4,pl.2,pl.5,pl.3,pl.6,nrow=3)

    #############################################