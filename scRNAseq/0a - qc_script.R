setwd("~/young-oldv2/scripts/qc")

do_qc<-function(data_params){
  rmarkdown::render("qc.Rmd", output_format="pdf_document", output_file=paste("~/young-oldv2/results",data_params$dataset,"qc.pdf",sep = "/"), params=list(dataset=data_params$dataset,nfrna_min =data_params$nfrna_min,nfrna_max =data_params$nfrna_max,pmt_max =data_params$pmt_max))
}

list_datasets=list(
  list(dataset="young",nfrna_min =200,nfrna_max =6000,pmt_max =8),
  list(dataset="old",nfrna_min =200,nfrna_max =6000,pmt_max =8)
)

lapply(list_datasets, do_qc)

