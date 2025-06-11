#!/usr/bin/env python3

import collections as cl

cmrlist = ['AGRN', 'RPL22', 'ESPN', 'MASP2', 'NPPA', 'HMGCL', 'CNR2', 'RHCE', 'PIGV', 'AGL', 'EXTL2', 'LRIG2', 'VANGL1', 'FLG', 'FLAD1', 'MUC1', 'HCN3', 'FCGR3A', 'FCGR2B', 'CD247', 'PRG4', 'PTPRC', 'KISS1', 'LBR', 'DCLRE1C', 'MRC1', 'PDSS1', 'PAPSS2', 'PTEN', 'LIPN', 'BTRC', 'ECHS1', 'SPRN', 'IFITM3', 'DRD4', 'DEAF1', 'MUC5B', 'H19', 'HBG1', 'MDK', 'SPI1', 'B3GAT3', 'ESRRA', 'SLC22A12', 'TPCN2', 'FGF3', 'CD4', 'PEX5', 'APOBEC1', 'GDF3', 'COX14', 'NACA', 'DPY19L2', 'HNF1A', 'KDM2B', 'HPD', 'P2RX2', 'PGAM5', 'GOLGA3', 'F7', 'PROZ', 'NDUFB1', 'JAG2', 'IGHV3-21', 'NDUFAF1', 'PPIP5K1', 'USP8', 'MPG', 'HBM', 'AXIN1', 'LMF1', 'SSTR5', 'GNPTG', 'UNKL', 'CLCN7', 'PKD1', 'PDPK1', 'SMG1', 'SH2B1', 'PHKG2', 'VKORC1', 'ORC6', 'HSD11B2', 'HYDIN', 'GCSH', 'ZNF469', 'CDH15', 'ANKRD11', 'MC1R', 'ABR', 'SERPINF2', 'SRR', 'P2RX5', 'GP1BA', 'ENO3', 'ALOXE3', 'HES7', 'ATPAF2', 'FOXN1', 'SUZ12', 'KRTAP1-1', 'G6PC3', 'HOXB8', 'GIP', 'TMC6', 'CANT1', 'RNF213', 'FSCN2', 'TYMS', 'GALR1', 'CTDP1', 'BSG', 'HCN2', 'KISS1R', 'STK11', 'TCF3', 'LMNB2', 'GIPC3', 'TBXA2R', 'CREB3L3', 'RFX2', 'C3', 'STXBP2', 'TYK2', 'TRMT1', 'CYP4F3', 'CYP4F12', 'CALR3', 'MYO9B', 'INSL3', 'FKBP8', 'GPI', 'COX6B1', 'IFNL3', 'CYP2G1P', 'ETHE1', 'APOC1', 'APOC4', 'APOC2', 'BLOC1S3', 'DMPK', 'FUT1', 'BAX', 'TRPM4', 'PNKP', 'SIGLEC16', 'KLK4', 'ETFB', 'NLRP12', 'PRKCG', 'MBOAT7', 'NLRP7', 'NLRP2', 'GP6', 'TNNT1', 'TNNI3', 'SBK3', 'ZNF419', 'SLC27A5', 'SNTG2', 'TPO', 'PXDN', 'KLF11', 'ABCG8', 'PCBP1', 'MOGS', 'RPIA', 'RGPD3', 'CASP10', 'MLPH', 'ANO7', 'D2HGDH', 'SEMG1', 'PLTP', 'CHRNA4', 'EEF1A2', 'PTK6', 'OPRL1', 'MYT1', 'KCNE1', 'CBR3', 'NDUFV3', 'CBS', 'U2AF1', 'CRYAA', 'TRAPPC10', 'DNMT3L', 'COL18A1', 'COL6A1', 'APOBEC3H', 'NDUFA6', 'CYB5R3', 'A4GALT', 'TTLL1', 'MLC1', 'TUBGCP6', 'CHL1', 'DAZL', 'LZTFL1', 'RHOA', 'MST1R', 'HYAL1', 'RPN1', 'BFSP2', 'PCCB', 'EIF2B5', 'KNG1', 'TM4SF19', 'ZNF141', 'PDE6B', 'GAK', 'FGFRL1', 'RNF212', 'UVSSA', 'DOK7', 'HMX1', 'UGT2A1', 'UGT2A2', 'AFP', 'DSPP', 'PDLIM3', 'FAT1', 'SLC6A18', 'TERT', 'SLC6A3', 'MARVELD2', 'SMN1', 'TTC37', 'LIX1', 'SAR1B', 'MYOT', 'PCDHA10', 'NPM1', 'FLT4', 'SLC17A5', 'SEC63', 'PCMT1', 'SLC22A1', 'SMOC2', 'THBS2', 'SLC29A4', 'NOD1', 'PPIA', 'PSPH', 'GUSB', 'CLIP2', 'GTF2I', 'NCF1', 'GTF2IRD2', 'ZAN', 'LAMB1', 'KLF14', 'BRAF', 'TRBV9', 'KMT2C', 'ARHGEF10', 'IKBKB', 'IMPA1', 'CDH17', 'GPIHBP1', 'MAFA', 'KCNV2', 'GALT', 'FXN', 'TJP2', 'MUSK', 'SLC27A4', 'PKN3', 'ABO', 'SOHLH1', 'INPP5E', 'MAN1B1', 'GRIN1']

useprefix  = {'CASP10_group2': 'CASP10,', 'DCLRE1C_group1': 'DCLRE1C,', 'PCCB_group1': 'PCCB,', 'KLF14_group1': 'KLF14,', 'PXDN_group11': 'PXDN,', 'TNNI3_group11': 'TNNI3,', 'ALOXE3_group2': 'ALOXE3,', 'H19_group1': 'H19,', 'NCF1_group1': 'NCF1,', 'ZAN_group5': 'ZAN,', 'MLC1_group1': 'MLC1,', 'C3_group1': 'C3,', 'TM4SF19_group1': 'TM4SF19,', 'P2RX_group5': 'P2RX5,', 'P2RX_group4': 'P2RX2,', 'MYT1_group9': 'MYT1,', 'NDUFA6_group1': 'NDUFA6,', 'ORC6_group8': 'ORC6,', 'PSPH_group1': 'PSPH,', 'D2HGDH_group1': 'D2HGDH,', 'SMOC2_group1': 'SMOC2,', 'SMOC2_group3': 'SMOC2,', 'SUZ12_group2': 'SUZ12,', 'HNF1A_group1': 'HNF1A,', 'ETHE1_group1': 'ETHE1,', 'TRPM4_group1': 'TRPM4,', 'AXIN1_group1': 'AXIN1,', 'FCGR3_group1': 'FCGR3A,', 'GCSH_group8': 'GCSH,', 'LAMB_group3': 'LAMB1,', 'TTLL1_group8': 'TTLL1,', 'CBS_group1': 'CBS,', 'PROZ_group1': 'PROZ,', 'LZTFL1_group1': 'LZTFL1,', 'LZTFL1_group2': 'LZTFL1,', 'FOXN1_group1': 'FOXN1,', 'CD4_group1': 'CD4,', 'PEX5_group5': 'PEX5,', 'DSPP_group1': 'DSPP,', 'EXTL2_group1': 'EXTL2,', 'RPN1_group1': 'RPN1,', 'ARHGEF10_group4': 'ARHGEF10,', 'NACA_group9': 'NACA,', 'HSD11B2_group1': 'HSD11B2,', 'STK11_group2': 'STK11,', 'ATPAF2_group1': 'ATPAF2,', 'MYO9B_group1': 'MYO9B,', 'GALT_group1': 'GALT,', 'GPI_group1': 'GPIHBP1,', 'GPI_group2': 'GPI,', 'PRG4_group1': 'PRG4,', 'KISS1_group1': 'KISS1,', 'KISS1_group2': 'KISS1R,', 'RPL22_group1': 'RPL22,', 'DRD4_group1': 'DRD4,', 'COL18A1_group1': 'COL18A1,', 'ETFB_group1': 'ETFB,', 'PKN3_group1': 'PKN3,', 'HBM_group1': 'HBM,', 'PDLIM3_group1': 'PDLIM3,', 'DMPK_group1': 'DMPK,', 'PCBP_group1': 'PCBP1,', 'GRIN1_group1': 'GRIN1,', 'CD247_group1': 'CD247,', 'CD247_group2': 'CD247,', 'RHOA_group1': 'RHOA,', 'GTF2I_group32': 'GTF2IRD2,', 'GTF2I_group33': 'GTF2I,', 'PPIP5K1_group1': 'PPIP5K1,', 'TMC6_group1': 'TMC6,', 'TYMS_group3': 'TYMS,', 'ENO3_group2': 'ENO3,', 'NPPA_group1': 'NPPA,', 'NPM1_group47': 'NPM1,', 'MLPH_group1': 'MLPH,', 'CALR3_group1': 'CALR3,', 'FLAD1_group1': 'FLAD1,', 'AFP_group1': 'AFP,', 'USP8_group4': 'USP8,', 'CLCN7_group1': 'CLCN7,', 'APOB_group34': 'APOBEC3H,', 'APOB_group26': 'APOBEC1,', 'CDH15_group1': 'CDH15,', 'LMNB2_group1': 'LMNB2,', 'SNTG2_group1': 'SNTG2,', 'SNTG2_group5': 'SNTG2,', 'DNMT3L_group1': 'DNMT3L,', 'PRKCG_group1': 'PRKCG,', 'HMGCL_group1': 'HMGCL,', 'CANT1_group1': 'CANT1,', 'SMG1_group4': 'SMG1,', 'PCDHA10_group3': 'PCDHA10,', 'PCDHA10_group5': 'PCDHA10,', 'PCDHA10_group6': 'PCDHA10,', 'BRAF_group1': 'BRAF,', 'BRAF_group3': 'BRAF,', 'PIGV_group1': 'PIGV,', 'GALR1_group1': 'GALR1,', 'COL6_group30': 'COL6A1,', 'GDF3_group1': 'GDF3,', 'ABO_group1': 'ABO,', 'ECHS1_group1': 'ECHS1,', 'PHKG2_group1': 'PHKG2,', 'NOD1_group1': 'NOD1,', 'NLRP12_group1': 'NLRP12,', 'BAX_group1': 'BAX,', 'GAK_group1': 'GAK,', 'PDSS1_group4': 'PDSS1,', 'B3GAT3_group2': 'B3GAT3,', 'RPIA_group1': 'RPIA,', 'CYP2G1P_group1': 'CYP2G1P,', 'DOK7_group1': 'DOK7,', 'ABCG8_group1': 'ABCG8,', 'PTEN_group2': 'PTEN,', 'PTEN_group4': 'PTEN,', 'MUSK_group2': 'MUSK,', 'MUSK_group3': 'MUSK,', 'LRIG2_group6': 'LRIG2,', 'SIGLEC16_group1': 'SIGLEC16,', 'SEMG1_group1': 'SEMG1,', 'UGT2A_group6': 'UGT2A1,UGT2A2,', 'SMN_group1': 'SMN1,', 'DPY19L2_group74': 'DPY19L2,', 'KCNV2_group1': 'KCNV2,', 'SEC63_group1': 'SEC63,', 'SLC6A_group18': 'SLC6A3,', 'SLC6A_group17': 'SLC6A18,', 'UVSSA_group1': 'UVSSA,', 'NDUFB1_group3': 'NDUFB1,', 'ESPN_group1': 'ESPN,', 'KMT2C_group1': 'KMT2C,', 'TPCN2_group1': 'TPCN2,', 'UNKL_group1': 'UNKL,', 'PCMT1_group1': 'PCMT1,', 'PCMT1_group2': 'PCMT1,', 'IKBKB_group1': 'IKBKB,', 'KLF11_group1': 'KLF11,', 'TNNT1_group1': 'TNNT1,', 'FCGR2_group1': 'FCGR2B,', 'JAG2_group1': 'JAG2,', 'MRC1_group1': 'MRC1,', 'F7_group1': 'F7,', 'PDPK1_group1': 'PDPK1,', 'SLC29A4_group1': 'SLC29A4,', 'ZNF141_group6': 'ZNF141,', 'KCNE1_group1': 'KCNE1,', 'KCNE1_group2': 'KCNE1,', 'TPO_group1': 'TPO,', 'TPO_group4': 'TPO,', 'G6PC3_group1': 'G6PC3,', 'DEAF1_group1': 'DEAF1,', 'NLRP2_group6': 'NLRP2,', 'TYK2_group1': 'TYK2,', 'PAPSS2_group1': 'PAPSS2,', 'HOXB8_group1': 'HOXB8,', 'ZNF419_group1': 'ZNF419,', 'RHCE_group1': 'RHCE,', 'CHRNA4_group1': 'CHRNA4,', 'GOLGA_group23': 'GOLGA3,', 'NLRP7_group1': 'NLRP7,', 'ABR_group8': 'ABR,', 'RNF212_group3': 'RNF212,', 'INSL3_group1': 'INSL3,', 'GP6_group13': 'GP6,', 'ANKRD11_group3': 'ANKRD11,', 'FSCN2_group1': 'FSCN2,', 'SAR1B_group40': 'SAR1B,', 'CBR3_group2': 'CBR3,', 'TERT_group1': 'TERT,', 'COX6B1_group7': 'COX6B1,', 'SPRN_group1': 'SPRN,', 'PLTP_group1': 'PLTP,', 'TRMT1_group2': 'TRMT1,', 'CHL1_group1': 'CHL1,', 'CHL1_group4': 'CHL1,', 'TBXA2R_group1': 'TBXA2R,', 'GUSB_group46': 'GUSB,', 'CNR2_group1': 'CNR2,', 'PGAM5_group2': 'PGAM5,', 'LIX1_group1': 'LIX1,', 'NDUFAF1_group1': 'NDUFAF1,', 'BFSP2_group3': 'BFSP2,', 'HES7_group1': 'HES7,', 'FLT4_group3': 'FLT4,', 'IFNL3_group1': 'IFNL3,', 'SPI1_group1': 'SPI1,', 'CREB3L3_group1': 'CREB3L3,', 'CYB5R3_group1': 'CYB5R3,', 'TRAPPC10_group1': 'TRAPPC10,', 'TRAPPC10_group2': 'TRAPPC10,', 'TRBV9_group1': 'TRBV9,', 'LBR_group1': 'LBR,', 'HYAL1_group1': 'HYAL1,', 'CYP4F3_group1': 'CYP4F3,', 'RFX2_group1': 'RFX2,', 'KLK4_group1': 'KLK4,', 'HMX1_group1': 'HMX1,', 'GNPTG_group1': 'GNPTG,', 'SERPINF2_group1': 'SERPINF2,', 'PTK6_group2': 'PTK6,', 'PNKP_group1': 'PNKP,', 'FGFRL1_group1': 'FGFRL1,', 'FXN_group1': 'FXN,', 'SRR_group9': 'SRR,', 'FAT1_group1': 'FAT1,', 'MC1R_group1': 'MC1R,', 'MBOAT_group14': 'MBOAT7,', 'MAN1B1_group1': 'MAN1B1,', 'FGF3_group1': 'FGF3,', 'SLC27A4_group1': 'SLC27A4,', 'FUT1_group1': 'FUT1,', 'SH2B_group3': 'SH2B1,', 'BSG_group1': 'BSG,', 'SBK3_group1': 'SBK3,', 'STXBP2_group1': 'STXBP2,', 'MPG_group1': 'MPG,', 'MYOT_group1': 'MYOT,', 'LMF1_group1': 'LMF1,', 'ESRRA_group1': 'ESRRA,', 'KRTAP1-1_group1': 'KRTAP1-1,', 'SSTR5_group1': 'SSTR5,', 'SLC22A_group88': 'SLC22A1,', 'SLC22A_group145': 'SLC22A12,', 'MDK_group1': 'MDK,', 'GIP_group40': 'GIP,', 'GIP_group41': 'GIPC3,', 'NDUFV3_group1': 'NDUFV3,', 'TCF3_group2': 'TCF3,', 'A4GALT_group1': 'A4GALT,', 'A4GALT_group2': 'A4GALT,', 'CYP4F12_group1': 'CYP4F12,', 'AGRN_group1': 'AGRN,', 'COL6A1_group1': 'COL6A1,', 'PDE6B_group1': 'PDE6B,', 'PTPRC_group1': 'PTPRC,', 'TUBGCP6_group1': 'TUBGCP6,', 'IGHV3-21_group2': 'IGHV3-21,', 'BTRC_group1': 'BTRC,', 'CRYAA_group152': 'CRYAA,', 'PPIA_group20': 'PPIA,', 'HYDIN_group19': 'HYDIN,', 'HYDIN_group1': 'HYDIN,', 'DAZ_group2': 'DAZL,', 'IFITM3_group8': 'IFITM3,', 'MST1_group9': 'MST1R,', 'SOHLH1_group1': 'SOHLH1,', 'PKD1_group0011201': 'PKD1,', 'RNF213_group34': 'RNF213,', 'EIF2B5_group1': 'EIF2B5,', 'IMPA1_group2': 'IMPA1,', 'HCN_group1': 'HCN3,', 'HCN_group14': 'HCN2,', 'MUC5B_group2': 'MUC5B,', 'COX14_group1': 'COX14,', 'FKBP8_group2': 'FKBP8,', 'INPP5E_group1': 'INPP5E,', 'EEF1A2_group1': 'EEF1A2,', 'MASP2_group1': 'MASP2,', 'KNG1_group1': 'KNG1,', 'MUC1_group1': 'MUC1,', 'FLG_group2': 'FLG,', 'CTDP1_group1': 'CTDP1,', 'RGPD_group2': 'RGPD3,', 'KDM2_group4': 'KDM2B,', 'SLC27A5_group19': 'SLC27A5,', 'BLOC1S3_group1': 'BLOC1S3,', 'TTC37_group1': 'TTC37,', 'HBG_group1': 'HBG1,', 'ANO7_group2': 'ANO7,', 'ZNF469_group1': 'ZNF469,', 'SLC17A5_group1': 'SLC17A5,', 'HP_group65': 'HPD,', 'MAFA_group1': 'MAFA,', 'VANGL1_group1': 'VANGL1,', 'OPRL1_group1': 'OPRL1,', 'AGL_group1': 'AGL,', 'CDH17_group1': 'CDH17,', 'LIPN_group1': 'LIPN,', 'THBS2_group1': 'THBS2,', 'MOG_group1': 'MOGS,', 'MARVELD2_group1': 'MARVELD2,', 'VKORC1_group7': 'VKORC1,', 'TJP2_group1': 'TJP2,', 'CLIP2_group1': 'CLIP2,', 'APOC_group5': 'APOC1,APOC4,APOC2,', 'GP1BA_group1': 'GP1BA,', 'U2AF1_group3':'U2AF1,'}

import numpy as np 

file2 = "/Users/wangfeima/Documents/MarkLab/Ctyper/Locityper/new/locityper-benchmarking/QV/loo_illumina.csv"

file1 = "/Users/wangfeima/Documents/MarkLab/Ctyper/Locityper/nonleaves_overlap_results.txt"

ref_region = "/Users/wangfeima/Documents/MarkLab/Ctyper/Locityper/cmroverlap.list"
ctyper_region = "/Users/wangfeima/Documents/MarkLab/Ctyper/Locityper/alltitles.txt"

allgenes = cl.defaultdict(list)
with open(file2, mode = 'r') as r, open(file2_output, mode = 'w') as w:
	
	for line_ in r:
					
		line = line_.split()
		gene = line[1]
		
		if (gene in cmrlist and "hap" in line) or line_.startswith("sample"):
			w.write(line_)
		
		allgenes[line[1]].append(line[-2])
		

ref_sizes = cl.defaultdict(int)
with open(ref_region, mode = 'r') as f:
	for line in f:
		line = line.strip().split()
		
		thesize = int(line[2]) - int(line[1]) 
		gene = line[3]
		ref_sizes[gene] = thesize
		
		
gene_togroup = cl.defaultdict(str)
		
overlap_sizes = cl.defaultdict(lambda: cl.defaultdict(int))
with open(ctyper_region, mode = 'r') as f:
	for line in f:
		line = line.strip().split()
		
		thesize = int(line[2]) - int(line[1]) 
		group = "_".join(line[0][1:].split('_')[0:2])
		sample = "_".join(line[0].split('_')[2:4])
		
		genes = useprefix[group].split(",")[:-1]
		
		overlap_sizes[tuple(genes)][sample] += thesize
		
		for gene in genes:
			gene_togroup[gene] = group



	
ctypergenes = cl.defaultdict(list)
with open(file1, mode = 'r') as f:
	
	for line_ in f:
		
		line = line_.split()
		if len(line) == 0 :
			continue
		name = "_".join(line[0].split('_')[:2])
		
		if name not in useprefix:
			continue
		
		genes = useprefix[name].split(",")[:-1]
		
		ctypergenes[tuple(genes)].append(line[-2])



allgenesset= set()
for genes,benchs in ctypergenes.items():
		
	
	benchs = [float(x) for x in benchs if x != '-1']
	
	refsize = sum([ref_sizes[gene] for gene in genes])
	
	allsizes = [samplesize for sample, samplesize in overlap_sizes[genes].items()]
	
	median = np.median(allsizes)
	
	if name in useprefix and len([x for x in genes if x in allgenes]):
				
		loci = [float(x) for x in sum([allgenes[x] for x in genes],[])  if x != 'NA']
		
		benchs_zeros, loci_zeros = benchs.count(0.0)/len(benchs),loci.count(0.0)/len(loci)
		
		benchs_mean, loci_mean = sum(benchs)/len(benchs), sum(loci)/len(loci)
		
		allgenesset.update(genes)
		print(gene_togroup[genes[0]], ";".join(genes),benchs_zeros, loci_zeros,'{:.2e}'.format(benchs_mean), '{:.2e}'.format(loci_mean), refsize, median)
		
	elif name in useprefix:
		
				
		allgenesset.update(genes)
		benchs_zeros = benchs.count(0.0)/len(benchs) if len(benchs) else "NA"
		
		benchs_mean = sum(benchs)/len(benchs) if len(benchs) else "NA"
		
		print(gene_togroup[genes[0]],  ";".join(genes),benchs_zeros, "NA",'{:.2e}'.format(benchs_mean) if len(benchs) else "NA", "NA", refsize, median)

