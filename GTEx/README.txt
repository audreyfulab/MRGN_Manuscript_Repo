***New GTEx Tissues Analysis With Updated PC Selection and Only Protein Coding and LncRNA Trios***
==================================================================================================

*Procedure for extracting only Protein Coding and LncRNA Trios
 
- From BIOMART we obtained the following meta-information:
	-Gene.Stable.ID
	-Gene.type

- Then (for each tissue) we obtained the files from Badsha's directory:

	/tissuename_AllPC/data.snp.cis.trans.final.tissuename.V8.unique.snps.RData 
		- the n X 3p matrix where n is the number of samples and p is the number of trios
		  (each trio has 3 columns: SNP,cis.gene,trans.gene)
	
	/PCs.matrix.tissuename.RData
		the n X m matrix of PC scores 

- We then extract the cis and trans gene stable ID's from each trio and match them to the BIOMART data to obtain the 
  gene types

- After obtaining the gene types, we exclude all trios which are not one of the combinations below:

		      |     cis.gene    |    trans.gene   |
                       -----------------------------------
		      |	protein_coding  |  protein_coding |
		      |	    lncRNA      |      lncRNA     |
		      |	protein_coding  |      lncRNA     |
		      |	    lncRNA      |  protein_coding |

- Note that this procedure also excludes trios which did not have a listed gene type in BIOMART (i.e it excludes NA's) 
  however the number of NA's in each tisue data set were quite low (generally < 20) 

- We save the dataframe of the subset of trios containing only protein coding and lncRNA genes as
	/trios_subset_data/all.data.unqiue.snps.pclrna.only.tissuename.RData

- This results in a new set of trios of dimension n X 3p_sub where p_sub is the number of trios in the subset

- After obtaining the subset of trios with the desired gene types we proceed to the PC selection stage:::

*Procedure for determining significant PCs

- Using the psych::corr.test(,adjust="none") function to perform the pearson correlation test, we obtain the pvalues
  between the scores on each PC and the entire subset of trios. We apply the qvalue correction to the pvalues of each PC 
  to control the FDR at 10%.

- We then allocate the column indices that correspond to SNPs, cis.genes, or trans genes which were significant at the
  qvalue cutoff to a list 

- the resulting list is of length m corresponding the number of PCs for the tissue and each element in the list is a 
  vector of column indices from the subset trio dataframe. 

- each list is saved to 
	/Lists_asso_PCS/List.significant.asso1.tissuename.RData

- Next iterating over the number of trios in the subset and the number of PCs

	- we extract each trio from the subset and determine if any of the column names match the columns identified for
	  each PC
		- if they are matched the index of the PC is allocated to that trio

- after proceeding through the iterations we obtain a list of length p_sub where each element in the list is a vector
  containing the column indices of the PCs that are correlated with the p_sub{th} trio

- we save the list for each tissue as
	/list_matched_sig_trios/List.Match.significant.trios.tissuename.RData

-Using the indices of PCs belonging to each trio from above ^ we construct a list of length p_sub where each element in the 
 list is a dataframe of dimension n x 3+u containing the data for each trio in the subset and its associated PCs

- we save the list of data frames for each tissue as:
	/final_trios_with_PCs/data.with.PCs.tissuename.RData
 
 ***Procedure For Analyzing LncRNA and Protein_Coding Trios***
==============================================================


           