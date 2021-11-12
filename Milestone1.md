# Section 1: Update

Now it is able to:
* in the "[scripts](https://github.com/XiaoyunZhouusc/final_project/scripts)" folder, there is a [script](https://github.com/XiaoyunZhouusc/final_project/scripts/unzip.sh) to automatically download data and unzip
* in the "[R](https://github.com/XiaoyunZhouusc/final_project/R)" folder, [final_project.R](https://github.com/XiaoyunZhouusc/final_project/R/final_project.R) now can do simple data processing and comparing and output chart. I compare the portions of smokers and non-smokers among patients of pancreatic cancer and output a chart. 

![](images/pie_smoker_vs_non-smoker.png?raw=true)

Section 2: Next Steps
Once figured out how to use HTSeq-count files, I am going to compare genes of smoker and non-smoker.  

Section 3: Data.  

Data can be retrieved at [here](https://drive.google.com/file/d/1bXAHETXs5_UlzRMJEQr3QhdYygPe-wmv/view?usp=sharing). These are clinicial and htseq-count files downloaded at [GDC](https://portal.gdc.cancer.gov/repository?facetTab=cases&filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.primary_site%22%2C%22value%22%3A%5B%22pancreas%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.project.program.name%22%2C%22value%22%3A%5B%22TCGA%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.analysis.workflow_type%22%2C%22value%22%3A%5B%22HTSeq%20-%20Counts%22%5D%7D%7D%5D%7D&searchTableTab=cases)

Section 4: Known Issues. 

Right now don't know how to use HTSeq-count files. I thought the file name is the case id number used by clinical and exposure data as I was planing comparing genes of smoker and non-smoker, but turns out these are totally unrelated. 