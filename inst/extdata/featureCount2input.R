library(readr)
suppressMessages(library(biomaRt))
library(org.Hs.eg.db)

inpData = read.delim("/Volumes/CBUtechsZeus/bernd/tp2022/outs/TP2022/all.2022.count.s2.txt", sep = "\t",skip = 1, header = T)

# ENSG to gene name
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://oct2022.archive.ensembl.org")
mart = human 
# ensemblDB <- "hsapiens_gene_ensembl"
# mart <- useMart(biomart = "ensembl", dataset = ensemblDB)
geneNames = getBM(
  attributes = c(
    "ensembl_gene_id",
    "external_gene_name"
    ),
  filters = "ensembl_gene_id",
  values = unique(inpData$Geneid),
  mart = mart
)
class(geneNames)
class(inpData)
rownames(geneNames) = geneNames$ensembl_gene_id
rownames(inpData) = inpData$Geneid
inpData$geneName = "NA"
inpData["geneName"] = geneNames[rownames(inpData),"external_gene_name"]
idealData = inpData[which(!inpData["geneName"] == ""),]

rownames(idealData) = make.unique(idealData[["geneName"]])
idealData = idealData[,stringr::str_starts(colnames(idealData),"TP")]
colnames(idealData) = sub(x=colnames(idealData), pattern = ".*(S..?)_.*",replacement = "\\1")
# D1e_O_P_A_f vale01 antiCD3CD28 pos old val20 D1 e F
sampleNames = c(S1 = "D1g_O_N_N_f",
  S2 = "D1g_O_N_A_f",
  S3 = "D1g_O_N_P_f",
  S4 = "D1g_O_N_L_f",
  S5 = "D2g_O_N_N_f",
  S6 = "D2g_O_N_A_f",
  S7 = "D2g_O_N_P_f",
  S8 = "D2g_O_N_L_f",
  S9 = "D3g_O_N_N_f",
  S10 = "D3g_O_N_A_f",
  S11 = "D3g_O_N_P_f",
  S12 = "D3g_O_N_L_f",
  S13 = "D4g_O_N_N_f",
  S14 = "D4g_O_N_A_f",
  S15 = "D4g_O_N_P_f",
  S16 = "D4g_O_N_L_f")
colnames(idealData) = sampleNames[colnames(idealData)]
  
write.csv(x= idealData, file = "TP2022.input.csv",quote = F)

idealAnnot = data.frame(row.names = colnames(idealData), name = colnames(idealData))
idealAnnot$Stimulus = substr(idealAnnot$name,start = 9,stop = 9)
idealAnnot$CMVstatus = "neg"
idealAnnot$sex = "F"
idealAnnot$age = substr(idealAnnot$name,start = 5,stop = 5)
idealAnnot$Donor = substr(idealAnnot$name,start = 1,stop = 3)
write.csv(x= idealAnnot, file = "TP2022.annot.csv",quote = F)



