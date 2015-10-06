This is a tool for bulk and efficient mapping of novel sequences to eggNOG v4.1
orthologous groups (OGs). It uses HMMER hmmpgmd to load several precomputed
databases in memory and perform parallelized searches. 

Three databases are currently available: 
- euk (Eukaryotes)
- bact (Bacteria)
- arch (Archea)

### Requirements: 
- Python 2.7+
- BioPython
- HMMER 3 available in your path
- sqlite3 
- ~130GB of disk
- ~130GB of RAM memory (~90GB for the eukaryotes DB, ~32 for the bacteria DB) 

### Install: 
- Clone or download this repository
```
git clone https://github.com/jhcepas/eggnog-mapper.git
```

- run `upgrade.py` to download and parse the hmm and annotation databases. It may take a while and will require 150GB of disk space.
```
cd eggnog-mapper
python upgrade.py
```

### Usage: 
1) start the HMMER server in a separate terminal, specifying the number of CPUs available and the target databases to load into memory. Note that the full set of databases will require ~130GB of RAM. 
Loading the databases may take a few minutes. 
```
   python server.py --db euk bact arch --cpu 20
```

2) map sequences to eggnog orthologous groups. While `server.py` is running, use a different shell to execute `eggnog_mapper.py`. Basic usage is as follows:
```
   python eggnog_mapper.py --output output.hits --db euk bact --evalue 0.001 FASTA_File.fa
```

2) annotate sequences based on the resulting hits using `annotate.py`. GO terms, KEGG pathways, SMART/PFAM domains and general functional descriptions are available. x 
```
   python annotate.py output.hits --go --kegg --smart --desc 
```

### Example: annotate a set of bacterial sequences 
```
eggnog_mapper/$ python server.py --db bact

# in a different shell
python eggnog_mapper.py test/testCOG0515.fa --output testCOG0515.hits.tsv --db bact --evalue 0.001 --maxhits 20
python annotate.py testCOG0515.hits.tsv --go --kegg --smart --desc --maxhits 1 > testCOG0515.annotations.tsv
```

**Output**
```
# query	OG	level	evalue	score	GO category	GO id	GO desc	evindence	nseqs	freq in OG (%)
# query	OG	level	evalue	score	KEGG pathway	nseqs	freq in OG
# query	OG	level	evalue	score	Domain source	domain name	nseqs	freq in OG (%)
1000565.METUNv1_00364	05D9P	bactNOG	5.10153612776e-100	339.312927246	Description	serine threonine protein kinase
1000565.METUNv1_00364	05D9P	bactNOG	5.10153612776e-100	339.312927246	COG Categories	T
1000565.METUNv1_00364	05D9P	bactNOG	5.10153612776e-100	339.312927246	GO_Molecular Function	GO:0005488	binding	IEA, ISS, IDA	1900	96.5
1000565.METUNv1_00364	05D9P	bactNOG	5.10153612776e-100	339.312927246	GO_Molecular Function	GO:0003824	catalytic activity	IEA, ISS, IDA	1895	96.3
1000565.METUNv1_00364	05D9P	bactNOG	5.10153612776e-100	339.312927246	GO_Molecular Function	GO:1901363	heterocyclic compound binding	IEA	1899	96.5
1000565.METUNv1_00364	05D9P	bactNOG	5.10153612776e-100	339.312927246	GO_Molecular Function	GO:0016301	kinase activity	IEA, ISS, IDA	1893	96.2
1000565.METUNv1_00364	05D9P	bactNOG	5.10153612776e-100	339.312927246	GO_Molecular Function	GO:0097159	organic cyclic compound binding	IEA	1899	96.5
1000565.METUNv1_00364	05D9P	bactNOG	5.10153612776e-100	339.312927246	GO_Molecular Function	GO:0004672	protein kinase activity	IEA, ISS, IDA	1893	96.2
1000565.METUNv1_00364	05D9P	bactNOG	5.10153612776e-100	339.312927246	GO_Molecular Function	GO:0016740	transferase activity	IEA, ISS, IDA	1893	96.2
1000565.METUNv1_00364	05D9P	bactNOG	5.10153612776e-100	339.312927246	GO_Molecular Function	GO:0016773	phosphotransferase activity, alcohol group as acceptor	IEA, ISS, IDA	1893	96.2
1000565.METUNv1_00364	05D9P	bactNOG	5.10153612776e-100	339.312927246	GO_Molecular Function	GO:0016772	transferase activity, transferring phosphorus-containing groups	IEA, ISS, IDA	1893	96.2
1000565.METUNv1_00364	05D9P	bactNOG	5.10153612776e-100	339.312927246	GO_Molecular Function	GO:0043167	ion binding	IEA, IDA	1894	96.2
1000565.METUNv1_00364	05D9P	bactNOG	5.10153612776e-100	339.312927246	GO_Molecular Function	GO:0043168	anion binding	IEA	1894	96.2
1000565.METUNv1_00364	05D9P	bactNOG	5.10153612776e-100	339.312927246	GO_Molecular Function	GO:0000166	nucleotide binding	IEA	1893	96.2
1000565.METUNv1_00364	05D9P	bactNOG	5.10153612776e-100	339.312927246	GO_Molecular Function	GO:1901265	nucleoside phosphate binding	IEA	1893	96.2
1000565.METUNv1_00364	05D9P	bactNOG	5.10153612776e-100	339.312927246	GO_Molecular Function	GO:0032549	ribonucleoside binding	IEA	1893	96.2
1000565.METUNv1_00364	05D9P	bactNOG	5.10153612776e-100	339.312927246	GO_Molecular Function	GO:0017076	purine nucleotide binding	IEA	1893	96.2
1000565.METUNv1_00364	05D9P	bactNOG	5.10153612776e-100	339.312927246	GO_Molecular Function	GO:0005524	ATP binding	IEA	1893	96.2
1000565.METUNv1_00364	05D9P	bactNOG	5.10153612776e-100	339.312927246	GO_Molecular Function	GO:0032559	adenyl ribonucleotide binding	IEA	1893	96.2
1000565.METUNv1_00364	05D9P	bactNOG	5.10153612776e-100	339.312927246	GO_Molecular Function	GO:0032555	purine ribonucleotide binding	IEA	1893	96.2
1000565.METUNv1_00364	05D9P	bactNOG	5.10153612776e-100	339.312927246	GO_Molecular Function	GO:0032553	ribonucleotide binding	IEA	1893	96.2
1000565.METUNv1_00364	05D9P	bactNOG	5.10153612776e-100	339.312927246	GO_Molecular Function	GO:0035639	purine ribonucleoside triphosphate binding	IEA	1893	96.2
1000565.METUNv1_00364	05D9P	bactNOG	5.10153612776e-100	339.312927246	GO_Molecular Function	GO:0030554	adenyl nucleotide binding	IEA	1893	96.2
1000565.METUNv1_00364	05D9P	bactNOG	5.10153612776e-100	339.312927246	GO_Molecular Function	GO:0032550	purine ribonucleoside binding	IEA	1893	96.2
1000565.METUNv1_00364	05D9P	bactNOG	5.10153612776e-100	339.312927246	GO_Molecular Function	GO:0001883	purine nucleoside binding	IEA	1893	96.2
1000565.METUNv1_00364	05D9P	bactNOG	5.10153612776e-100	339.312927246	GO_Molecular Function	GO:0001882	nucleoside binding	IEA	1893	96.2
1000565.METUNv1_00364	05D9P	bactNOG	5.10153612776e-100	339.312927246	GO_Molecular Function	GO:0036094	small molecule binding	IEA	1893	96.2
1000565.METUNv1_00364	05D9P	bactNOG	5.10153612776e-100	339.312927246	GO_Molecular Function	GO:0004674	protein serine/threonine kinase activity	IEA, ISS, IDA	1802	91.6
1000565.METUNv1_00364	05D9P	bactNOG	5.10153612776e-100	339.312927246	GO_Molecular Function	GO:0033218	amide binding	IEA	616	31.3
1000565.METUNv1_00364	05D9P	bactNOG	5.10153612776e-100	339.312927246	GO_Molecular Function	GO:0008658	penicillin binding	IEA	616	31.3
1000565.METUNv1_00364	05D9P	bactNOG	5.10153612776e-100	339.312927246	GO_Molecular Function	GO:0031406	carboxylic acid binding	IEA	616	31.3
1000565.METUNv1_00364	05D9P	bactNOG	5.10153612776e-100	339.312927246	GO_Molecular Function	GO:0008144	drug binding	IEA	616	31.3
1000565.METUNv1_00364	05D9P	bactNOG	5.10153612776e-100	339.312927246	GO_Molecular Function	GO:0033293	monocarboxylic acid binding	IEA	616	31.3
1000565.METUNv1_00364	05D9P	bactNOG	5.10153612776e-100	339.312927246	GO_Molecular Function	GO:1901681	sulfur compound binding	IEA	616	31.3
1000565.METUNv1_00364	05D9P	bactNOG	5.10153612776e-100	339.312927246	GO_Biological Process	GO:0009987	cellular process	IEA, IMP, ISS, IDA	1900	96.5
1000565.METUNv1_00364	05D9P	bactNOG	5.10153612776e-100	339.312927246	GO_Biological Process	GO:0008152	metabolic process	IEA, ISS, IDA	1900	96.5
1000565.METUNv1_00364	05D9P	bactNOG	5.10153612776e-100	339.312927246	GO_Biological Process	GO:0044237	cellular metabolic process	IEA, ISS, IDA	1900	96.5
1000565.METUNv1_00364	05D9P	bactNOG	5.10153612776e-100	339.312927246	GO_Biological Process	GO:0044260	cellular macromolecule metabolic process	IEA, ISS, IDA	1899	96.5
1000565.METUNv1_00364	05D9P	bactNOG	5.10153612776e-100	339.312927246	GO_Biological Process	GO:0043170	macromolecule metabolic process	IEA, ISS, IDA	1899	96.5
1000565.METUNv1_00364	05D9P	bactNOG	5.10153612776e-100	339.312927246	GO_Biological Process	GO:0071704	organic substance metabolic process	IEA, ISS, IDA	1899	96.5
1000565.METUNv1_00364	05D9P	bactNOG	5.10153612776e-100	339.312927246	GO_Biological Process	GO:0044238	primary metabolic process	IEA, ISS, IDA	1899	96.5
1000565.METUNv1_00364	05D9P	bactNOG	5.10153612776e-100	339.312927246	GO_Biological Process	GO:0019538	protein metabolic process	IEA, ISS, IDA	1893	96.2
1000565.METUNv1_00364	05D9P	bactNOG	5.10153612776e-100	339.312927246	GO_Biological Process	GO:0044267	cellular protein metabolic process	IEA, ISS, IDA	1893	96.2
1000565.METUNv1_00364	05D9P	bactNOG	5.10153612776e-100	339.312927246	GO_Biological Process	GO:0043412	macromolecule modification	IEA, ISS, IDA	1893	96.2
1000565.METUNv1_00364	05D9P	bactNOG	5.10153612776e-100	339.312927246	GO_Biological Process	GO:0036211	protein modification process	IEA, ISS, IDA	1893	96.2
1000565.METUNv1_00364	05D9P	bactNOG	5.10153612776e-100	339.312927246	GO_Biological Process	GO:0016310	phosphorylation	IEA, ISS, IDA	1893	96.2
1000565.METUNv1_00364	05D9P	bactNOG	5.10153612776e-100	339.312927246	GO_Biological Process	GO:0006468	protein phosphorylation	IEA, ISS, IDA	1893	96.2
1000565.METUNv1_00364	05D9P	bactNOG	5.10153612776e-100	339.312927246	GO_Biological Process	GO:0006464	cellular protein modification process	IEA, ISS, IDA	1893	96.2
1000565.METUNv1_00364	05D9P	bactNOG	5.10153612776e-100	339.312927246	GO_Biological Process	GO:0006796	phosphate-containing compound metabolic process	IEA, ISS, IDA	1893	96.2
1000565.METUNv1_00364	05D9P	bactNOG	5.10153612776e-100	339.312927246	GO_Biological Process	GO:0006793	phosphorus metabolic process	IEA, ISS, IDA	1893	96.2
1000565.METUNv1_00364	05D9P	bactNOG	5.10153612776e-100	339.312927246	PFAM	PASTA	2044	104
1000565.METUNv1_00364	05D9P	bactNOG	5.10153612776e-100	339.312927246	PFAM	Pkinase_Tyr	1972	100
1000565.METUNv1_00364	05D9P	bactNOG	5.10153612776e-100	339.312927246	PFAM	Pkinase	1968	100
1000565.METUNv1_00364	05D9P	bactNOG	5.10153612776e-100	339.312927246	PFAM	Kdo	232	11.8
1000565.METUNv1_00364	05D9P	bactNOG	5.10153612776e-100	339.312927246	PFAM	APH	188	9.55
1000565.METUNv1_00364	05D9P	bactNOG	5.10153612776e-100	339.312927246	PFAM	WD40	139	7.06
1000565.METUNv1_00364	05D9P	bactNOG	5.10153612776e-100	339.312927246	PFAM	RIO1	110	5.59
1000565.METUNv1_00364	05D9P	bactNOG	5.10153612776e-100	339.312927246	PFAM	PD40	102	5.18
1000565.METUNv1_00364	05D9P	bactNOG	5.10153612776e-100	339.312927246	SMART	PASTA	2055	104
1000565.METUNv1_00364	05D9P	bactNOG	5.10153612776e-100	339.312927246	SMART	TRANS	1617	82.2
1000565.METUNv1_00364	05D9P	bactNOG	5.10153612776e-100	339.312927246	SMART	STYKc	1598	81.2
1000565.METUNv1_00364	05D9P	bactNOG	5.10153612776e-100	339.312927246	SMART	TPR	429	21.8
1000565.METUNv1_00364	05D9P	bactNOG	5.10153612776e-100	339.312927246	SMART	S_TKc	371	18.9
1000565.METUNv1_00364	05D9P	bactNOG	5.10153612776e-100	339.312927246	SMART	COIL	185	9.4
1000565.METUNv1_00364	05D9P	bactNOG	5.10153612776e-100	339.312927246	SMART	WD40	183	9.3


(...)

```