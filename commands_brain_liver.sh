--------------------------------------------------------------------------------------
Gene level RNA-seq analysis
--------------------------------------------------------------------------------------

#make dir for brain and liver comparison and move to it
mkdir analysis_brain_x_liver
cd analysis_brain_x_liver

# Run the DE (using the script) with the edgeR package
edgeR.analysis.R --input_matrix ../quantifications/encode.mouse.gene.expected_count.idr_NA.tsv \
                 --metadata /tutorial/data/gene.quantifications.index.tsv \
                 --fields tissue \
                 --coefficient 3 \
                 --output brain_X_liver

#to obtain a list of genes overexpressed in brain compared to liver:
awk '$NF<0.01 && $2<-10{print $1"\tover_brain_X_liver"}' edgeR.cpm1.n2.brain_X_liver.tsv > edgeR.0.01.over_brain_X_liver.txt

#to obtain a list of the genes overexpressed in liver compared to brain:
awk '$NF<0.01 && $2>10 {print $1"\tover_liver_X_brain"}' edgeR.cpm1.n2.brain_X_liver.tsv > edgeR.0.01.over_liver_X_brain.txt

# prepare for heatmap
awk '$3=="gene"{ match($0, /gene_id "([^"]+).+gene_type "([^"]+)/, var); print var[1],var[2] }' OFS="\t" /tutorial/refs/gencode.vM4.gtf \
| join.py --file1 stdin \
          --file2 <(cat edgeR.0.01.over*.txt) \
| sed '1igene\tedgeR\tgene_type' > gene.edgeR.tsv

# heatmap
cut -f1 gene.edgeR.tsv \
| tail -n+2 \
| selectMatrixRows.sh - ../quantifications/encode.mouse.gene.TPM.idr_NA.tsv \
| ggheatmap.R --width 5 \
              --height 8 \
              --col_metadata /tutorial/data/gene.quantifications.index.tsv \
              --colSide_by tissue \
              --col_labels labExpId \
              --row_metadata gene.edgeR.tsv \
              --merge_row_mdata_on gene \
              --rowSide_by edgeR,gene_type \
              --row_labels none \
              --log \
              --pseudocount 0.1 \
              --col_dendro \
              --row_dendro \
              --matrix_palette /tutorial/palettes/palDiverging.txt \
              --colSide_palette /tutorial/palettes/palTissue.txt \
              --output heatmap.brain_X_liver.pdf

# copy the universe of genes from the analysis performed in class
cp /tutorial/analysis/universe.txt .

# run GO enrichment for BP over_brain_X_liver
awk '{split($1,a,"."); print a[1]}' edgeR.0.01.over_brain_X_liver.txt \
| GO_enrichment.R --universe universe.txt \
                  --genes stdin \
                  --categ BP \
                  --output edgeR.over_brain_X_liver \
                  --species mouse

awk 'NR==1{$1="% "$1}{print $1,$2}' edgeR.over_brain_X_liver.BP.tsv

# run GO enrichment for BP over_liver_X_brain
awk '{split($1,a,"."); print a[1]}' edgeR.0.01.over_liver_X_brain.txt \
| GO_enrichment.R --universe universe.txt \
                  --genes stdin \
                  --categ BP \
                  --output edgeR.over_liver_X_brain \
                  --species mouse

awk 'NR==1{$1="% "$1}{print $1,$2}' edgeR.over_liver_X_brain.BP.tsv

--------------------------------------------------------------------------------------
Isoform level RNA-seq analysis
--------------------------------------------------------------------------------------

# make dir for splicing brain_liver
mkdir /tutorial/splicing_brain_liver
cd /tutorial/splicing_brain_liver

# NOTE: some input files are reused (from directory splicing) cause are the same we used in class

# generate alternative splicing events
suppa.py generateEvents -i ../splicing/exon-annot.gtf -e SE RI MX FL -o localEvents -f ioe

# Compute percent spliced in index (PSI) values for local events
for event in SE RI MX AF; do suppa.py psiPerEvent --total-filter 10 --ioe-file localEvents_${event}_strict.ioe --expression-file ../splicing/pc-tx.tsv -o PSI-${event} ; done

# Create individual files for tissues
for event in SE RI MX AF; do for tissue in Brain Liver Heart; do selectMatrixColumns.sh PRNAembryo${tissue}1:PRNAembryo${tissue}2 PSI-${event}.psi > ${tissue}.${event}.psi; done; done

# Differential splicing analysis for local events
for event in SE RI MX AF; do suppa.py diffSplice --method empirical --input localEvents_${event}_strict.ioe --psi Brain.${event}.psi Liver.${event}.psi Heart.${event}.psi --tpm ../splicing/expr.Brain.tsv ../splicing/expr.Liver.tsv ../splicing/expr.Heart.tsv -c -gc -o DS.${event}; done


# skipping exon with deltaPSI 0.5
event=SE; awk 'BEGIN{FS=OFS="\t"}NR>1 && $2!="nan" && ($2>0.5 || $2<-0.5) && $3<0.05{print $1}' DS.${event}.dpsi > top.${event}.events.txt

# intron retention and mutually exclusive with deltaPSI 0.3
for event in RI MX; do awk 'BEGIN{FS=OFS="\t"}NR>1 && $2!="nan" && ($2>0.3 || $2<-0.3) && $3<0.05{print $1}' DS.${event}.dpsi > top.${event}.events.txt; done

# alternative first with deltaPSI 0.7
event=AF; awk 'BEGIN{FS=OFS="\t"}NR>1 && $2!="nan" && ($2>0.7 || $2<-0.7) && $3<0.05{print $1}' DS.${event}.dpsi > top.${event}.events.txt

# subset the file with PSI vectors
for event in SE RI MX AF; do 
selectMatrixRows.sh top.${event}.events.txt DS.${event}.psivec| sed 's/:-//;s/:+//;s/ENSMUS//' > input_heatmap_Brain-Liver_top_${event}_events.tsv; done

#heatmap
for event in SE RI MX AF; do 
ggheatmap.R -i input_heatmap_Brain-Liver_top_${event}_events.tsv -o output_heatmap_Brain-Liver_top_${event}_events.pdf --matrix_palette='/tutorial/palettes/palSequential.txt' --col_metadata /tutorial/splicing/mdata.tsv --colSide_by Tissue --colSide_palette /tutorial/palettes/palTissue.txt; done

--------------------------------------------------------------------------------------
ChIP-seq data analysis
--------------------------------------------------------------------------------------
mkdir chip-analysis_brain_x_liver

# Compare peaks

# common between Brain and Liver
bedtools intersect -a /tutorial/results/CHIPembryoBrain.narrowPeak -b /tutorial/results/CHIPembryoLiver.narrowPeak > commonPeaks_Brain_Liver.bed

# common among the three tissues
bedtools intersect -a /tutorial/results/CHIPembryoHeart.narrowPeak -b commonPeaks_Brain_Liver.bed > common3tissues.bed

# Brain specific peaks and Liver specific peaks compared to common among three tissues
for tissue in Brain Liver; do bedtools intersect -a /tutorial/results/CHIPembryo${tissue}.narrowPeak -b common3tissues.bed -v > ${tissue}-specific_peaks.bed ; done

# Barplot number of common and specific peaks
wc -l *bed | grep -v total | grep -v common3 | awk '{print $2"\t"$1}'| sed '1iPeaks\tCounts' | sed 's/.bed//' > peaks.tsv

(grep -v "#800080" /tutorial/palettes/palTissue.txt && echo "#4a0505") > palBarplot.txt

ggbarplot.R --input peaks.tsv --output barplot_peaks.pdf --header --title "Peaks" --x_title "Intersection peaks" --y_title "Number of peaks" --palette_fill palBarplot.txt --fill_by 1 --x_order peaks.tsv


--------------------------------------------------------------------------------------
Visualization
--------------------------------------------------------------------------------------

# Genomic location

awk 'BEGIN{FS=OFS="\t"}$3=="gene" && $7=="+" && $0~/gene_type "protein_coding"/{ match($0, /gene_id "([^"]+)/, id); print $1,$4-200,$4+200, id[1] }' /tutorial/refs/gencode.vM4.gtf > positiveStrand_gene_TSS.bed 

awk 'BEGIN{FS=OFS="\t"}$3=="gene" && $7=="-" && $0~/gene_type "protein_coding"/{ match($0, /gene_id "([^"]+)/, id); print $1,$5-200,$5+200, id[1] }' /tutorial/refs/gencode.vM4.gtf > negativeStrand_gene_TSS.bed

cat *gene_TSS.bed | bedtools sort > gene_TSS.bed


# TSSs of protein-coding DE genes

# over_brain_X_liver

cat gene_TSS.bed | grep -Ff <(cut -f1 /tutorial/analysis_brain_x_liver/edgeR.0.01.over_brain_X_liver.txt) | bedtools sort > DE_over_brain_X_liver_TSS_gene.bed
#### son exactamente lo mismo ####
cut -f1 /tutorial/analysis_brain_x_liver/edgeR.0.01.over_brain_X_liver.txt | grep -Ff - gene_TSS.bed | bedtools sort > DE_over_brain_X_liver_TSS_gene.bed

# intersect with H3K4me3 peaks shared between brain-liver
bedtools intersect -a DE_over_brain_X_liver_TSS_gene.bed -b commonPeaks_Brain_Liver.bed -wa > DE_over_brain_X_liver_TSS_gene_commonPeaks_Brain_Liver.bed

# intersect with H3K4me3 peaks brain-specific
bedtools intersect -a DE_over_brain_X_liver_TSS_gene.bed -b Brain-specific_peaks.bed -wa > DE_over_brain_X_liver_TSS_gene_Brain-specific_peaks.bed


# over_liver_X_brain

cut -f1 /tutorial/analysis_brain_x_liver/edgeR.0.01.over_liver_X_brain.txt | grep -Ff - gene_TSS.bed | bedtools sort > DE_over_liver_X_brain_TSS_gene.bed

# intersect with H3K4me3 peaks shared between liver-brain
bedtools intersect -a DE_over_liver_X_brain_TSS_gene.bed -b commonPeaks_Brain_Liver.bed -wa > DE_over_liver_X_brain_TSS_gene_commonPeaks_Brain_Liver.bed

# intersect with H3K4me3 peaks liver-specific
bedtools intersect -a DE_over_liver_X_brain_TSS_gene.bed -b Liver-specific_peaks.bed -wa > DE_over_liver_X_brain_TSS_gene_Liver-specific_peaks.bed


# examples 

DE_over_brain_X_liver_TSS_gene_Brain-specific_peaks.bed
chr2    113828921-113829321	ENSMUSG00000023236.7 scg5

DE_over_brain_X_liver_TSS_gene_commonPeaks_Brain_Liver.bed
chr10   81025462-81025862        ENSMUSG00000043670.4   diras1

DE_over_liver_X_brain_TSS_gene_Liver-specific_peaks.bed
chr1    88406761-88407161        ENSMUSG00000026295.8   spp2


# examples AF exons

chr4:43037206-43039966:43040288:43037206-43045648:43045797  Fam214b
chr9:44342071-44342221:44342381:44342071-44344010:44344228  Hmbs