#!/usr/bin/env Rscript
# Perform over-representation test on functional terms : gene ontology & KEGG pathway
# Usage : Rscript funcEnrich.R <genesName.list> <organism.id>
# Usage : <genesName.list> List of genes name/symbol or entrez ID (preffered), one per line
# Usage : <organism.id> Organism (Human : human , Mouse : mouse)

symbol_to_entrezID <- function(x){

	entrezID <- mapIds(get(annoDb), keys=x, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
	return(entrezID)
}

get_genome <- function(x){
	if ( organism=="human" ) {

		annoDb<<- "org.Hs.eg.db"
		kegg_genome<<- "hsa"
		reactome_genome<<- "human"

	} else if ( organism=="mouse" ) {

		annoDb<<- "org.Mm.eg.db"
		kegg_genome<<- "mmu"
		reactome_genome<<- "mouse"

	} else {
		print ("Unknown organism given, please select either human or mouse")
	}
}

go_analysis <- function(x){

	#BP
	assign(paste0(filename, "_gobp"), enrichGO(gene=x, OrgDb=get(annoDb), ont="BP", pAdjustMethod="BH", pvalueCutoff=1.0, qvalueCutoff=1.0, readable=T))
	write.table(file=paste0(filename,"_goBP.tsv"), get(paste0(filename,"_gobp")), row.names=F, sep="\t", quote=F)

#	png(filename=paste0(filename,"_goBP_barplot.png"), width = 1268, height = 768)
#	print(barplot(get(paste0(filename, "_gobp")), showCategory=50))
#	dev.off()
	
	#MF
	assign(paste0(filename, "_gomf"), enrichGO(gene=x, OrgDb=get(annoDb), ont="MF", pAdjustMethod="BH", pvalueCutoff=1.0, qvalueCutoff=1.0, readable=T))
	write.table(file=paste0(filename,"_goMF.tsv"), get(paste0(filename,"_gomf")), row.names=F, sep="\t", quote=F)
	
#	png(filename=paste0(filename,"_goMF_barplot.png"), width = 1268, height = 768)
#	print(barplot(get(paste0(filename,"_gomf")), showCategory=50))
#	dev.off()
	
	#CC
	assign(paste0(filename, "_gocc"), enrichGO(gene=x, OrgDb=get(annoDb), ont="CC", pAdjustMethod="BH", pvalueCutoff=1.0, qvalueCutoff=1.0, readable=T))
	write.table(file=paste0(filename,"_goCC.tsv"), get(paste0(filename,"_gocc")), row.names=F, sep="\t", quote=F)

#	png(filename=paste0(filename,"_goCC_barplot.png"), width = 1268, height = 768)
#	print(barplot(get(paste0(filename,"_gocc")), showCategory=50))
#	dev.off()
}

pathway_analysis <- function(x){

	kegg<-enrichKEGG(x, organism=kegg_genome, pAdjustMethod = "BH", pvalueCutoff=1.0, qvalueCutoff=1.0)
	kegg_readable<-setReadable(kegg, get(annoDb), keyType="ENTREZID")
	write.table(file=paste0(filename,"_kegg.tsv"), kegg_readable, row.names=F, sep="\t", quote=F)

#	png(filename=paste0(filename,"_kegg_barplot.png"), width = 1268, height = 768)
#	print(barplot(kegg, showCategory=50))
#	dev.off()

	reactome<-enrichPathway(x, organism=reactome_genome, pAdjustMethod = "BH", pvalueCutoff=1.0, qvalueCutoff=1.0, readable=T)
	write.table(file=paste0(filename,"_reactome.tsv"), reactome, row.names=F, sep="\t", quote=F)

#	png(filename=paste0(filename,"_reactome_barplot.png"), width = 1268, height = 768)
#	print(barplot(reactome, showCategory=50))
#	dev.off()
}

args = commandArgs(trailingOnly=TRUE)
filename <- args[1]
organism <- args[2]
get_genome(organism)

library("clusterProfiler", quietly = TRUE)
library(annoDb, character.only=TRUE, quietly = TRUE)
library("ReactomePA", quietly = TRUE)

genes <- scan(file=filename, character(), quote = "") 

if(grepl('[A-Za-z]', genes[1])){
	print ("Genes name/symbol given, converting to entrez ID")
	genes_entrezID <- symbol_to_entrezID(genes)
	go_analysis(genes_entrezID)
	pathway_analysis(genes_entrezID)

}else{
	print ("Entrez ID given, performing over-representation test")
	go_analysis(genes)
	pathway_analysis(genes)
}
