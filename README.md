This is a tool for bulk and efficient annotation of protein sequenes using
eggNOG v4.1 orthologous groups (OGs). It uses precomputed HMM models and
HMMER-hmmpgmd to keep databases in memory and perform parallelized searches.

Functional mappings from EggNOG v4.1 include: GO terms, KEGG pathways,
SMART/PFAM domains, COG functional categories and a consensus functional
description.

Three taxonomy-restricted databases are currently available: 
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
python annotate.py testCOG0515.hits.tsv --go --kegg --smart --desc > testCOG0515.annotations.tsv
```

**Output**
```
1000565.METUNv1_00364	12	12	PFAM	APH	17555	proNOG	16	2.12268729674e-134	53	499	3	431	2	12.5
                     			PFAM	APH	16RV1	proNOG	205	4.25739996312e-109	517	781	3	275	24	11.7
                     			PFAM	APH	05D9P	bactNOG	1968	5.10153612776e-100	472	722	4	276	188	9.55
                     			PFAM	APH	07SGC	bactNOG	309	1.2821207656e-94	193	466	3	350	40	12.9
                     			PFAM	APH	0GB8Z	delNOG	65	9.47173825047e-92	30	298	4	276	11	16.9
                     			PFAM	APH	07QPV	bactNOG	109	5.11078410131e-91	163	432	5	278	7	6.42
                     			PFAM	APH	0HCNY	dproNOG	20	3.02140608055e-89	37	300	5	275	4	20
                     			PFAM	APH	0HG69	dproNOG	41	8.9972854538e-89	56	325	5	275	6	14.6
                     			PFAM	APH	0D42W	chloNOG	16	4.12832858376e-82	10	273	4	274	2	12.5
                     			PFAM	APH	00CIX	actNOG	255	2.87646308227e-81	64	345	6	356	13	5.1
                     			PFAM	APH	07UJ4	bactNOG	11	3.06302846003e-81	21	285	5	276	2	18.2
                     			PFAM	APH	08HMP	bactNOG	36	4.75359887693e-81	43	315	3	277	4	11.1
1000565.METUNv1_00364	1	1	PFAM	CHASE2	0HCNY	dproNOG	20	3.02140608055e-89	37	300	5	275	3	15
1000565.METUNv1_00364	1	1	PFAM	DUF1080	0D42W	chloNOG	16	4.12832858376e-82	10	273	4	274	1	6.25
1000565.METUNv1_00364	1	1	PFAM	FGE-sulfatase	0HG69	dproNOG	41	8.9972854538e-89	56	325	5	275	3	7.32
1000565.METUNv1_00364	1	1	PFAM	FHA	0D42W	chloNOG	16	4.12832858376e-82	10	273	4	274	1	6.25
1000565.METUNv1_00364	1	1	PFAM	HAMP	0QM1J	gproNOG	21	6.28581547187e-94	262	530	4	276	3	14.3
1000565.METUNv1_00364	12	12	PFAM	Kdo	17555	proNOG	16	2.12268729674e-134	53	499	3	431	3	18.8
                     			PFAM	Kdo	16RV1	proNOG	205	4.25739996312e-109	517	781	3	275	25	12.2
                     			PFAM	Kdo	05D9P	bactNOG	1968	5.10153612776e-100	472	722	4	276	232	11.8
                     			PFAM	Kdo	07SGC	bactNOG	309	1.2821207656e-94	193	466	3	350	81	26.2
                     			PFAM	Kdo	0QM1J	gproNOG	21	6.28581547187e-94	262	530	4	276	3	14.3
                     			PFAM	Kdo	0GB8Z	delNOG	65	9.47173825047e-92	30	298	4	276	17	26.2
                     			PFAM	Kdo	07QPV	bactNOG	109	5.11078410131e-91	163	432	5	278	10	9.17
                     			PFAM	Kdo	0HCNY	dproNOG	20	3.02140608055e-89	37	300	5	275	3	15
                     			PFAM	Kdo	0HG69	dproNOG	41	8.9972854538e-89	56	325	5	275	12	29.3
                     			PFAM	Kdo	0ND8P	firmNOG	472	8.53493970042e-88	52	315	6	277	101	21.4
                     			PFAM	Kdo	0EP0E	cloNOG	219	1.18732803766e-82	70	332	6	277	25	11.4
                     			PFAM	Kdo	08HMP	bactNOG	36	4.75359887693e-81	43	315	3	277	2	5.56
1000565.METUNv1_00364	1	1	PFAM	Kelch_1	07QPV	bactNOG	109	5.11078410131e-91	163	432	5	278	8	7.34
1000565.METUNv1_00364	5	5	PFAM	PASTA	05D9P	bactNOG	1968	5.10153612776e-100	472	722	4	276	2044	104
                     			PFAM	PASTA	07QPV	bactNOG	109	5.11078410131e-91	163	432	5	278	48	44
                     			PFAM	PASTA	0ND8P	firmNOG	472	8.53493970042e-88	52	315	6	277	1149	243
                     			PFAM	PASTA	0EP0E	cloNOG	219	1.18732803766e-82	70	332	6	277	480	219
                     			PFAM	PASTA	00CIX	actNOG	255	2.87646308227e-81	64	345	6	356	821	322
1000565.METUNv1_00364	3	3	PFAM	PD40	05D9P	bactNOG	1968	5.10153612776e-100	472	722	4	276	102	5.18
                     			PFAM	PD40	07UJ4	bactNOG	11	3.06302846003e-81	21	285	5	276	2	18.2
                     			PFAM	PD40	08HMP	bactNOG	36	4.75359887693e-81	43	315	3	277	4	11.1
1000565.METUNv1_00364	4	4	PFAM	PEGA	16RV1	proNOG	205	4.25739996312e-109	517	781	3	275	16	7.8
                     			PFAM	PEGA	0GB8Z	delNOG	65	9.47173825047e-92	30	298	4	276	6	9.23
                     			PFAM	PEGA	0HG69	dproNOG	41	8.9972854538e-89	56	325	5	275	6	14.6
                     			PFAM	PEGA	0BCA6	bproNOG	12	1.1100354572e-88	2	274	1	279	2	16.7
1000565.METUNv1_00364	1	1	PFAM	PP2C	07SGC	bactNOG	309	1.2821207656e-94	193	466	3	350	99	32
1000565.METUNv1_00364	20	20	PFAM	Pkinase	0BJMK	bproNOG	2	6.92995759218e-187	1	462	1	454	2	100
                     			PFAM	Pkinase	0B9TC	bproNOG	6	3.81207122537e-146	2	428	3	429	6	100
                     			PFAM	Pkinase	17555	proNOG	16	2.12268729674e-134	53	499	3	431	18	112
                     			PFAM	Pkinase	16RV1	proNOG	205	4.25739996312e-109	517	781	3	275	206	100
                     			PFAM	Pkinase	1755X	proNOG	13	8.40025975035e-109	4	453	2	447	13	100
                     			PFAM	Pkinase	05D9P	bactNOG	1968	5.10153612776e-100	472	722	4	276	1968	100
                     			PFAM	Pkinase	07SGC	bactNOG	309	1.2821207656e-94	193	466	3	350	313	101
                     			PFAM	Pkinase	0QM1J	gproNOG	21	6.28581547187e-94	262	530	4	276	21	100
                     			PFAM	Pkinase	0GB8Z	delNOG	65	9.47173825047e-92	30	298	4	276	66	102
                     			PFAM	Pkinase	07QPV	bactNOG	109	5.11078410131e-91	163	432	5	278	110	101
                     			PFAM	Pkinase	0HCNY	dproNOG	20	3.02140608055e-89	37	300	5	275	21	105
                     			PFAM	Pkinase	0HG69	dproNOG	41	8.9972854538e-89	56	325	5	275	41	100
                     			PFAM	Pkinase	0BCA6	bproNOG	12	1.1100354572e-88	2	274	1	279	12	100
                     			PFAM	Pkinase	0ND8P	firmNOG	472	8.53493970042e-88	52	315	6	277	472	100
                     			PFAM	Pkinase	0EP0E	cloNOG	219	1.18732803766e-82	70	332	6	277	219	100
                     			PFAM	Pkinase	0D42W	chloNOG	16	4.12832858376e-82	10	273	4	274	16	100
                     			PFAM	Pkinase	0B8I8	bproNOG	11	1.0062412472e-81	3	274	1	279	11	100
                     			PFAM	Pkinase	00CIX	actNOG	255	2.87646308227e-81	64	345	6	356	255	100
                     			PFAM	Pkinase	07UJ4	bactNOG	11	3.06302846003e-81	21	285	5	276	11	100
                     			PFAM	Pkinase	08HMP	bactNOG	36	4.75359887693e-81	43	315	3	277	37	103
1000565.METUNv1_00364	20	20	PFAM	Pkinase_Tyr	0BJMK	bproNOG	2	6.92995759218e-187	1	462	1	454	2	100
                     			PFAM	Pkinase_Tyr	0B9TC	bproNOG	6	3.81207122537e-146	2	428	3	429	6	100
                     			PFAM	Pkinase_Tyr	17555	proNOG	16	2.12268729674e-134	53	499	3	431	17	106
                     			PFAM	Pkinase_Tyr	16RV1	proNOG	205	4.25739996312e-109	517	781	3	275	206	100
                     			PFAM	Pkinase_Tyr	1755X	proNOG	13	8.40025975035e-109	4	453	2	447	13	100
                     			PFAM	Pkinase_Tyr	05D9P	bactNOG	1968	5.10153612776e-100	472	722	4	276	1972	100
                     			PFAM	Pkinase_Tyr	07SGC	bactNOG	309	1.2821207656e-94	193	466	3	350	314	102
                     			PFAM	Pkinase_Tyr	0QM1J	gproNOG	21	6.28581547187e-94	262	530	4	276	21	100
                     			PFAM	Pkinase_Tyr	0GB8Z	delNOG	65	9.47173825047e-92	30	298	4	276	66	102
                     			PFAM	Pkinase_Tyr	07QPV	bactNOG	109	5.11078410131e-91	163	432	5	278	110	101
                     			PFAM	Pkinase_Tyr	0HCNY	dproNOG	20	3.02140608055e-89	37	300	5	275	21	105
                     			PFAM	Pkinase_Tyr	0HG69	dproNOG	41	8.9972854538e-89	56	325	5	275	41	100
                     			PFAM	Pkinase_Tyr	0BCA6	bproNOG	12	1.1100354572e-88	2	274	1	279	12	100
                     			PFAM	Pkinase_Tyr	0ND8P	firmNOG	472	8.53493970042e-88	52	315	6	277	472	100
                     			PFAM	Pkinase_Tyr	0EP0E	cloNOG	219	1.18732803766e-82	70	332	6	277	219	100
                     			PFAM	Pkinase_Tyr	0D42W	chloNOG	16	4.12832858376e-82	10	273	4	274	16	100
                     			PFAM	Pkinase_Tyr	0B8I8	bproNOG	11	1.0062412472e-81	3	274	1	279	11	100
                     			PFAM	Pkinase_Tyr	00CIX	actNOG	255	2.87646308227e-81	64	345	6	356	255	100
                     			PFAM	Pkinase_Tyr	07UJ4	bactNOG	11	3.06302846003e-81	21	285	5	276	11	100
                     			PFAM	Pkinase_Tyr	08HMP	bactNOG	36	4.75359887693e-81	43	315	3	277	40	111
1000565.METUNv1_00364	3	3	PFAM	RIO1	05D9P	bactNOG	1968	5.10153612776e-100	472	722	4	276	110	5.59
                     			PFAM	RIO1	07QPV	bactNOG	109	5.11078410131e-91	163	432	5	278	12	11
                     			PFAM	RIO1	0ND8P	firmNOG	472	8.53493970042e-88	52	315	6	277	82	17.4
1000565.METUNv1_00364	1	1	PFAM	SpoIIE	07SGC	bactNOG	309	1.2821207656e-94	193	466	3	350	115	37.2
1000565.METUNv1_00364	4	4	PFAM	TPR_1	0GB8Z	delNOG	65	9.47173825047e-92	30	298	4	276	7	10.8
                     			PFAM	TPR_1	0HCNY	dproNOG	20	3.02140608055e-89	37	300	5	275	5	25
                     			PFAM	TPR_1	0D42W	chloNOG	16	4.12832858376e-82	10	273	4	274	17	106
                     			PFAM	TPR_1	08HMP	bactNOG	36	4.75359887693e-81	43	315	3	277	5	13.9
1000565.METUNv1_00364	4	4	PFAM	TPR_2	0GB8Z	delNOG	65	9.47173825047e-92	30	298	4	276	6	9.23
                     			PFAM	TPR_2	0HCNY	dproNOG	20	3.02140608055e-89	37	300	5	275	4	20
                     			PFAM	TPR_2	0D42W	chloNOG	16	4.12832858376e-82	10	273	4	274	11	68.8
                     			PFAM	TPR_2	08HMP	bactNOG	36	4.75359887693e-81	43	315	3	277	2	5.56
1000565.METUNv1_00364	1	1	PFAM	TPR_4	0D42W	chloNOG	16	4.12832858376e-82	10	273	4	274	2	12.5
1000565.METUNv1_00364	1	1	PFAM	Usp	07SGC	bactNOG	309	1.2821207656e-94	193	466	3	350	27	8.74
1000565.METUNv1_00364	5	5	PFAM	WD40	1755X	proNOG	13	8.40025975035e-109	4	453	2	447	11	84.6
                     			PFAM	WD40	05D9P	bactNOG	1968	5.10153612776e-100	472	722	4	276	139	7.06
                     			PFAM	WD40	0GB8Z	delNOG	65	9.47173825047e-92	30	298	4	276	4	6.15
                     			PFAM	WD40	0HCNY	dproNOG	20	3.02140608055e-89	37	300	5	275	4	20
                     			PFAM	WD40	08HMP	bactNOG	36	4.75359887693e-81	43	315	3	277	374	1.04e+03
1000565.METUNv1_00364	5	5	PFAM	cNMP_binding	0BJMK	bproNOG	2	6.92995759218e-187	1	462	1	454	2	100
                     			PFAM	cNMP_binding	0B9TC	bproNOG	6	3.81207122537e-146	2	428	3	429	6	100
                     			PFAM	cNMP_binding	17555	proNOG	16	2.12268729674e-134	53	499	3	431	4	25
                     			PFAM	cNMP_binding	1755X	proNOG	13	8.40025975035e-109	4	453	2	447	11	84.6
                     			PFAM	cNMP_binding	0HCNY	dproNOG	20	3.02140608055e-89	37	300	5	275	1	5
1000565.METUNv1_00364	1	1	PFAM	eIF2A	08HMP	bactNOG	36	4.75359887693e-81	43	315	3	277	6	16.7
1000565.METUNv1_00364	1	1	SMART	CHASE2	0HCNY	dproNOG	20	3.02140608055e-89	37	300	5	275	3	15
1000565.METUNv1_00364	12	12	SMART	COIL	16RV1	proNOG	205	4.25739996312e-109	517	781	3	275	23	11.2
                     			SMART	COIL	05D9P	bactNOG	1968	5.10153612776e-100	472	722	4	276	185	9.4
                     			SMART	COIL	0GB8Z	delNOG	65	9.47173825047e-92	30	298	4	276	8	12.3
                     			SMART	COIL	07QPV	bactNOG	109	5.11078410131e-91	163	432	5	278	9	8.26
                     			SMART	COIL	0HCNY	dproNOG	20	3.02140608055e-89	37	300	5	275	2	10
                     			SMART	COIL	0HG69	dproNOG	41	8.9972854538e-89	56	325	5	275	6	14.6
                     			SMART	COIL	0BCA6	bproNOG	12	1.1100354572e-88	2	274	1	279	1	8.33
                     			SMART	COIL	0ND8P	firmNOG	472	8.53493970042e-88	52	315	6	277	50	10.6
                     			SMART	COIL	0EP0E	cloNOG	219	1.18732803766e-82	70	332	6	277	31	14.2
                     			SMART	COIL	0D42W	chloNOG	16	4.12832858376e-82	10	273	4	274	3	18.8
                     			SMART	COIL	0B8I8	bproNOG	11	1.0062412472e-81	3	274	1	279	3	27.3
                     			SMART	COIL	08HMP	bactNOG	36	4.75359887693e-81	43	315	3	277	8	22.2
1000565.METUNv1_00364	1	1	SMART	FHA	0D42W	chloNOG	16	4.12832858376e-82	10	273	4	274	1	6.25
1000565.METUNv1_00364	1	1	SMART	HAMP	0QM1J	gproNOG	21	6.28581547187e-94	262	530	4	276	3	14.3
1000565.METUNv1_00364	1	1	SMART	Kelch	07QPV	bactNOG	109	5.11078410131e-91	163	432	5	278	10	9.17
1000565.METUNv1_00364	5	5	SMART	PASTA	05D9P	bactNOG	1968	5.10153612776e-100	472	722	4	276	2055	104
                     			SMART	PASTA	07QPV	bactNOG	109	5.11078410131e-91	163	432	5	278	46	42.2
                     			SMART	PASTA	0ND8P	firmNOG	472	8.53493970042e-88	52	315	6	277	1164	247
                     			SMART	PASTA	0EP0E	cloNOG	219	1.18732803766e-82	70	332	6	277	486	222
                     			SMART	PASTA	00CIX	actNOG	255	2.87646308227e-81	64	345	6	356	834	327
1000565.METUNv1_00364	1	1	SMART	PP2C_SIG	07SGC	bactNOG	309	1.2821207656e-94	193	466	3	350	136	44
1000565.METUNv1_00364	1	1	SMART	PP2Cc	07SGC	bactNOG	309	1.2821207656e-94	193	466	3	350	136	44
1000565.METUNv1_00364	2	2	SMART	PQQ	07QPV	bactNOG	109	5.11078410131e-91	163	432	5	278	7	6.42
                     			SMART	PQQ	08HMP	bactNOG	36	4.75359887693e-81	43	315	3	277	3	8.33
1000565.METUNv1_00364	20	20	SMART	STYKc	0BJMK	bproNOG	2	6.92995759218e-187	1	462	1	454	2	100
                     			SMART	STYKc	0B9TC	bproNOG	6	3.81207122537e-146	2	428	3	429	6	100
                     			SMART	STYKc	17555	proNOG	16	2.12268729674e-134	53	499	3	431	16	100
                     			SMART	STYKc	16RV1	proNOG	205	4.25739996312e-109	517	781	3	275	191	93.2
                     			SMART	STYKc	1755X	proNOG	13	8.40025975035e-109	4	453	2	447	13	100
                     			SMART	STYKc	05D9P	bactNOG	1968	5.10153612776e-100	472	722	4	276	1598	81.2
                     			SMART	STYKc	07SGC	bactNOG	309	1.2821207656e-94	193	466	3	350	313	101
                     			SMART	STYKc	0QM1J	gproNOG	21	6.28581547187e-94	262	530	4	276	21	100
                     			SMART	STYKc	0GB8Z	delNOG	65	9.47173825047e-92	30	298	4	276	59	90.8
                     			SMART	STYKc	07QPV	bactNOG	109	5.11078410131e-91	163	432	5	278	107	98.2
                     			SMART	STYKc	0HCNY	dproNOG	20	3.02140608055e-89	37	300	5	275	19	95
                     			SMART	STYKc	0HG69	dproNOG	41	8.9972854538e-89	56	325	5	275	36	87.8
                     			SMART	STYKc	0BCA6	bproNOG	12	1.1100354572e-88	2	274	1	279	9	75
                     			SMART	STYKc	0ND8P	firmNOG	472	8.53493970042e-88	52	315	6	277	255	54
                     			SMART	STYKc	0EP0E	cloNOG	219	1.18732803766e-82	70	332	6	277	119	54.3
                     			SMART	STYKc	0D42W	chloNOG	16	4.12832858376e-82	10	273	4	274	14	87.5
                     			SMART	STYKc	0B8I8	bproNOG	11	1.0062412472e-81	3	274	1	279	9	81.8
                     			SMART	STYKc	00CIX	actNOG	255	2.87646308227e-81	64	345	6	356	195	76.5
                     			SMART	STYKc	07UJ4	bactNOG	11	3.06302846003e-81	21	285	5	276	10	90.9
                     			SMART	STYKc	08HMP	bactNOG	36	4.75359887693e-81	43	315	3	277	36	100
1000565.METUNv1_00364	13	13	SMART	S_TKc	16RV1	proNOG	205	4.25739996312e-109	517	781	3	275	17	8.29
                     			SMART	S_TKc	05D9P	bactNOG	1968	5.10153612776e-100	472	722	4	276	371	18.9
                     			SMART	S_TKc	0GB8Z	delNOG	65	9.47173825047e-92	30	298	4	276	7	10.8
                     			SMART	S_TKc	07QPV	bactNOG	109	5.11078410131e-91	163	432	5	278	6	5.5
                     			SMART	S_TKc	0HCNY	dproNOG	20	3.02140608055e-89	37	300	5	275	2	10
                     			SMART	S_TKc	0HG69	dproNOG	41	8.9972854538e-89	56	325	5	275	5	12.2
                     			SMART	S_TKc	0BCA6	bproNOG	12	1.1100354572e-88	2	274	1	279	3	25
                     			SMART	S_TKc	0ND8P	firmNOG	472	8.53493970042e-88	52	315	6	277	218	46.2
                     			SMART	S_TKc	0EP0E	cloNOG	219	1.18732803766e-82	70	332	6	277	101	46.1
                     			SMART	S_TKc	0D42W	chloNOG	16	4.12832858376e-82	10	273	4	274	2	12.5
                     			SMART	S_TKc	0B8I8	bproNOG	11	1.0062412472e-81	3	274	1	279	2	18.2
                     			SMART	S_TKc	00CIX	actNOG	255	2.87646308227e-81	64	345	6	356	60	23.5
                     			SMART	S_TKc	07UJ4	bactNOG	11	3.06302846003e-81	21	285	5	276	1	9.09
1000565.METUNv1_00364	11	11	SMART	TPR	17555	proNOG	16	2.12268729674e-134	53	499	3	431	9	56.2
                     			SMART	TPR	16RV1	proNOG	205	4.25739996312e-109	517	781	3	275	50	24.4
                     			SMART	TPR	05D9P	bactNOG	1968	5.10153612776e-100	472	722	4	276	429	21.8
                     			SMART	TPR	07SGC	bactNOG	309	1.2821207656e-94	193	466	3	350	65	21
                     			SMART	TPR	0QM1J	gproNOG	21	6.28581547187e-94	262	530	4	276	3	14.3
                     			SMART	TPR	0GB8Z	delNOG	65	9.47173825047e-92	30	298	4	276	15	23.1
                     			SMART	TPR	0HCNY	dproNOG	20	3.02140608055e-89	37	300	5	275	9	45
                     			SMART	TPR	0HG69	dproNOG	41	8.9972854538e-89	56	325	5	275	6	14.6
                     			SMART	TPR	0ND8P	firmNOG	472	8.53493970042e-88	52	315	6	277	27	5.72
                     			SMART	TPR	0EP0E	cloNOG	219	1.18732803766e-82	70	332	6	277	27	12.3
                     			SMART	TPR	0D42W	chloNOG	16	4.12832858376e-82	10	273	4	274	61	381
1000565.METUNv1_00364	18	18	SMART	TRANS	17555	proNOG	16	2.12268729674e-134	53	499	3	431	5	31.2
                     			SMART	TRANS	16RV1	proNOG	205	4.25739996312e-109	517	781	3	275	158	77.1
                     			SMART	TRANS	1755X	proNOG	13	8.40025975035e-109	4	453	2	447	1	7.69
                     			SMART	TRANS	05D9P	bactNOG	1968	5.10153612776e-100	472	722	4	276	1617	82.2
                     			SMART	TRANS	07SGC	bactNOG	309	1.2821207656e-94	193	466	3	350	176	57
                     			SMART	TRANS	0QM1J	gproNOG	21	6.28581547187e-94	262	530	4	276	45	214
                     			SMART	TRANS	0GB8Z	delNOG	65	9.47173825047e-92	30	298	4	276	58	89.2
                     			SMART	TRANS	07QPV	bactNOG	109	5.11078410131e-91	163	432	5	278	84	77.1
                     			SMART	TRANS	0HCNY	dproNOG	20	3.02140608055e-89	37	300	5	275	26	130
                     			SMART	TRANS	0HG69	dproNOG	41	8.9972854538e-89	56	325	5	275	31	75.6
                     			SMART	TRANS	0BCA6	bproNOG	12	1.1100354572e-88	2	274	1	279	10	83.3
                     			SMART	TRANS	0ND8P	firmNOG	472	8.53493970042e-88	52	315	6	277	424	89.8
                     			SMART	TRANS	0EP0E	cloNOG	219	1.18732803766e-82	70	332	6	277	201	91.8
                     			SMART	TRANS	0D42W	chloNOG	16	4.12832858376e-82	10	273	4	274	13	81.2
                     			SMART	TRANS	0B8I8	bproNOG	11	1.0062412472e-81	3	274	1	279	3	27.3
                     			SMART	TRANS	00CIX	actNOG	255	2.87646308227e-81	64	345	6	356	238	93.3
                     			SMART	TRANS	07UJ4	bactNOG	11	3.06302846003e-81	21	285	5	276	18	164
                     			SMART	TRANS	08HMP	bactNOG	36	4.75359887693e-81	43	315	3	277	3	8.33
1000565.METUNv1_00364	1	1	SMART	VWA	0QM1J	gproNOG	21	6.28581547187e-94	262	530	4	276	2	9.52
1000565.METUNv1_00364	5	5	SMART	WD40	1755X	proNOG	13	8.40025975035e-109	4	453	2	447	11	84.6
                     			SMART	WD40	05D9P	bactNOG	1968	5.10153612776e-100	472	722	4	276	183	9.3
                     			SMART	WD40	0GB8Z	delNOG	65	9.47173825047e-92	30	298	4	276	8	12.3
                     			SMART	WD40	0HCNY	dproNOG	20	3.02140608055e-89	37	300	5	275	8	40
                     			SMART	WD40	08HMP	bactNOG	36	4.75359887693e-81	43	315	3	277	427	1.19e+03
1000565.METUNv1_00364	1	1	SMART	ZnF_RBZ	0D42W	chloNOG	16	4.12832858376e-82	10	273	4	274	2	12.5
1000565.METUNv1_00364	5	5	SMART	cNMP	0BJMK	bproNOG	2	6.92995759218e-187	1	462	1	454	2	100
                     			SMART	cNMP	0B9TC	bproNOG	6	3.81207122537e-146	2	428	3	429	6	100
                     			SMART	cNMP	17555	proNOG	16	2.12268729674e-134	53	499	3	431	4	25
                     			SMART	cNMP	1755X	proNOG	13	8.40025975035e-109	4	453	2	447	11	84.6
                     			SMART	cNMP	0HCNY	dproNOG	20	3.02140608055e-89	37	300	5	275	1	5
1000565.METUNv1_00364	1	1	cats	Cyclic nucleotide-binding protein	1755X	proNOG	13	8.40025975035e-109	4	453	2	447	0	100.0
1000565.METUNv1_00364	1	1	cats	Protein tyrosine kinase	0BJMK	bproNOG	2	6.92995759218e-187	1	462	1	454	0	100.0
1000565.METUNv1_00364	6	6	cats	Serine Threonine protein kinase	07QPV	bactNOG	109	5.11078410131e-91	163	432	5	278	0	100.0
                     			cats	Serine Threonine protein kinase	0HCNY	dproNOG	20	3.02140608055e-89	37	300	5	275	0	100.0
                     			cats	Serine Threonine protein kinase	0HG69	dproNOG	41	8.9972854538e-89	56	325	5	275	0	100.0
                     			cats	Serine Threonine protein kinase	0BCA6	bproNOG	12	1.1100354572e-88	2	274	1	279	0	100.0
                     			cats	Serine Threonine protein kinase	0D42W	chloNOG	16	4.12832858376e-82	10	273	4	274	0	100.0
                     			cats	Serine Threonine protein kinase	00CIX	actNOG	255	2.87646308227e-81	64	345	6	356	0	100.0
1000565.METUNv1_00364	1	1	cats	Serine threonine protein kinase with WD-40 repeats	08HMP	bactNOG	36	4.75359887693e-81	43	315	3	277	0	100.0
1000565.METUNv1_00364	1	1	cats	cyclic nucleotide-binding protein	0B9TC	bproNOG	6	3.81207122537e-146	2	428	3	429	0	100.0
1000565.METUNv1_00364	1	1	cats	protein kinase (EC	17555	proNOG	16	2.12268729674e-134	53	499	3	431	0	100.0
1000565.METUNv1_00364	1	1	cats	serine (threonine) protein kinase	0QM1J	gproNOG	21	6.28581547187e-94	262	530	4	276	0	100.0
1000565.METUNv1_00364	8	8	cats	serine threonine protein kinase	16RV1	proNOG	205	4.25739996312e-109	517	781	3	275	0	100.0
                     			cats	serine threonine protein kinase	05D9P	bactNOG	1968	5.10153612776e-100	472	722	4	276	0	100.0
                     			cats	serine threonine protein kinase	07SGC	bactNOG	309	1.2821207656e-94	193	466	3	350	0	100.0
                     			cats	serine threonine protein kinase	0GB8Z	delNOG	65	9.47173825047e-92	30	298	4	276	0	100.0
                     			cats	serine threonine protein kinase	0ND8P	firmNOG	472	8.53493970042e-88	52	315	6	277	0	100.0
                     			cats	serine threonine protein kinase	0EP0E	cloNOG	219	1.18732803766e-82	70	332	6	277	0	100.0
                     			cats	serine threonine protein kinase	0B8I8	bproNOG	11	1.0062412472e-81	3	274	1	279	0	100.0
                     			cats	serine threonine protein kinase	07UJ4	bactNOG	11	3.06302846003e-81	21	285	5	276	0	100.0
1000565.METUNv1_00364	1	1	desc	Cyclic nucleotide-binding protein	1755X	proNOG	13	8.40025975035e-109	4	453	2	447	0	100.0
1000565.METUNv1_00364	1	1	desc	Protein tyrosine kinase	0BJMK	bproNOG	2	6.92995759218e-187	1	462	1	454	0	100.0
1000565.METUNv1_00364	6	6	desc	Serine Threonine protein kinase	07QPV	bactNOG	109	5.11078410131e-91	163	432	5	278	0	100.0
                     			desc	Serine Threonine protein kinase	0HCNY	dproNOG	20	3.02140608055e-89	37	300	5	275	0	100.0
                     			desc	Serine Threonine protein kinase	0HG69	dproNOG	41	8.9972854538e-89	56	325	5	275	0	100.0
                     			desc	Serine Threonine protein kinase	0BCA6	bproNOG	12	1.1100354572e-88	2	274	1	279	0	100.0
                     			desc	Serine Threonine protein kinase	0D42W	chloNOG	16	4.12832858376e-82	10	273	4	274	0	100.0
                     			desc	Serine Threonine protein kinase	00CIX	actNOG	255	2.87646308227e-81	64	345	6	356	0	100.0
1000565.METUNv1_00364	1	1	desc	Serine threonine protein kinase with WD-40 repeats	08HMP	bactNOG	36	4.75359887693e-81	43	315	3	277	0	100.0
1000565.METUNv1_00364	1	1	desc	cyclic nucleotide-binding protein	0B9TC	bproNOG	6	3.81207122537e-146	2	428	3	429	0	100.0
1000565.METUNv1_00364	1	1	desc	protein kinase (EC	17555	proNOG	16	2.12268729674e-134	53	499	3	431	0	100.0
1000565.METUNv1_00364	1	1	desc	serine (threonine) protein kinase	0QM1J	gproNOG	21	6.28581547187e-94	262	530	4	276	0	100.0
1000565.METUNv1_00364	8	8	desc	serine threonine protein kinase	16RV1	proNOG	205	4.25739996312e-109	517	781	3	275	0	100.0
                     			desc	serine threonine protein kinase	05D9P	bactNOG	1968	5.10153612776e-100	472	722	4	276	0	100.0
                     			desc	serine threonine protein kinase	07SGC	bactNOG	309	1.2821207656e-94	193	466	3	350	0	100.0
                     			desc	serine threonine protein kinase	0GB8Z	delNOG	65	9.47173825047e-92	30	298	4	276	0	100.0
                     			desc	serine threonine protein kinase	0ND8P	firmNOG	472	8.53493970042e-88	52	315	6	277	0	100.0
                     			desc	serine threonine protein kinase	0EP0E	cloNOG	219	1.18732803766e-82	70	332	6	277	0	100.0
                     			desc	serine threonine protein kinase	0B8I8	bproNOG	11	1.0062412472e-81	3	274	1	279	0	100.0
                     			desc	serine threonine protein kinase	07UJ4	bactNOG	11	3.06302846003e-81	21	285	5	276	0	100.0
1000565.METUNv1_00364	20	20	go	GO:0000166	nucleotide binding	0BJMK	bproNOG	2	6.92995759218e-187	1	462	1	454	2	100	IEA
                     			go	GO:0000166	nucleotide binding	0B9TC	bproNOG	6	3.81207122537e-146	2	428	3	429	6	100	IEA
                     			go	GO:0000166	nucleotide binding	17555	proNOG	16	2.12268729674e-134	53	499	3	431	16	100	IEA
                     			go	GO:0000166	nucleotide binding	16RV1	proNOG	205	4.25739996312e-109	517	781	3	275	205	100	IEA
                     			go	GO:0000166	nucleotide binding	1755X	proNOG	13	8.40025975035e-109	4	453	2	447	13	100	IEA
                     			go	GO:0000166	nucleotide binding	05D9P	bactNOG	1968	5.10153612776e-100	472	722	4	276	1893	96.2	IEA
                     			go	GO:0000166	nucleotide binding	07SGC	bactNOG	309	1.2821207656e-94	193	466	3	350	304	98.4	IEA
                     			go	GO:0000166	nucleotide binding	0QM1J	gproNOG	21	6.28581547187e-94	262	530	4	276	21	100	IEA
                     			go	GO:0000166	nucleotide binding	0GB8Z	delNOG	65	9.47173825047e-92	30	298	4	276	65	100	IEA
                     			go	GO:0000166	nucleotide binding	07QPV	bactNOG	109	5.11078410131e-91	163	432	5	278	105	96.3	IEA
                     			go	GO:0000166	nucleotide binding	0HCNY	dproNOG	20	3.02140608055e-89	37	300	5	275	20	100	IEA
                     			go	GO:0000166	nucleotide binding	0HG69	dproNOG	41	8.9972854538e-89	56	325	5	275	41	100	IEA
                     			go	GO:0000166	nucleotide binding	0BCA6	bproNOG	12	1.1100354572e-88	2	274	1	279	12	100	IEA
                     			go	GO:0000166	nucleotide binding	0ND8P	firmNOG	472	8.53493970042e-88	52	315	6	277	448	94.9	IEA
                     			go	GO:0000166	nucleotide binding	0EP0E	cloNOG	219	1.18732803766e-82	70	332	6	277	212	96.8	IEA
                     			go	GO:0000166	nucleotide binding	0D42W	chloNOG	16	4.12832858376e-82	10	273	4	274	16	100	IEA
                     			go	GO:0000166	nucleotide binding	0B8I8	bproNOG	11	1.0062412472e-81	3	274	1	279	11	100	IEA
                     			go	GO:0000166	nucleotide binding	00CIX	actNOG	255	2.87646308227e-81	64	345	6	356	244	95.7	IEA
                     			go	GO:0000166	nucleotide binding	07UJ4	bactNOG	11	3.06302846003e-81	21	285	5	276	11	100	IEA
                     			go	GO:0000166	nucleotide binding	08HMP	bactNOG	36	4.75359887693e-81	43	315	3	277	36	100	IEA
1000565.METUNv1_00364	20	20	go	GO:0001882	nucleoside binding	0BJMK	bproNOG	2	6.92995759218e-187	1	462	1	454	2	100	IEA
                     			go	GO:0001882	nucleoside binding	0B9TC	bproNOG	6	3.81207122537e-146	2	428	3	429	6	100	IEA
                     			go	GO:0001882	nucleoside binding	17555	proNOG	16	2.12268729674e-134	53	499	3	431	16	100	IEA
                     			go	GO:0001882	nucleoside binding	16RV1	proNOG	205	4.25739996312e-109	517	781	3	275	205	100	IEA
                     			go	GO:0001882	nucleoside binding	1755X	proNOG	13	8.40025975035e-109	4	453	2	447	13	100	IEA
                     			go	GO:0001882	nucleoside binding	05D9P	bactNOG	1968	5.10153612776e-100	472	722	4	276	1893	96.2	IEA
                     			go	GO:0001882	nucleoside binding	07SGC	bactNOG	309	1.2821207656e-94	193	466	3	350	304	98.4	IEA
                     			go	GO:0001882	nucleoside binding	0QM1J	gproNOG	21	6.28581547187e-94	262	530	4	276	21	100	IEA
                     			go	GO:0001882	nucleoside binding	0GB8Z	delNOG	65	9.47173825047e-92	30	298	4	276	65	100	IEA
                     			go	GO:0001882	nucleoside binding	07QPV	bactNOG	109	5.11078410131e-91	163	432	5	278	105	96.3	IEA
                     			go	GO:0001882	nucleoside binding	0HCNY	dproNOG	20	3.02140608055e-89	37	300	5	275	20	100	IEA
                     			go	GO:0001882	nucleoside binding	0HG69	dproNOG	41	8.9972854538e-89	56	325	5	275	41	100	IEA
                     			go	GO:0001882	nucleoside binding	0BCA6	bproNOG	12	1.1100354572e-88	2	274	1	279	12	100	IEA
                     			go	GO:0001882	nucleoside binding	0ND8P	firmNOG	472	8.53493970042e-88	52	315	6	277	448	94.9	IEA
                     			go	GO:0001882	nucleoside binding	0EP0E	cloNOG	219	1.18732803766e-82	70	332	6	277	212	96.8	IEA
                     			go	GO:0001882	nucleoside binding	0D42W	chloNOG	16	4.12832858376e-82	10	273	4	274	16	100	IEA
                     			go	GO:0001882	nucleoside binding	0B8I8	bproNOG	11	1.0062412472e-81	3	274	1	279	11	100	IEA
                     			go	GO:0001882	nucleoside binding	00CIX	actNOG	255	2.87646308227e-81	64	345	6	356	244	95.7	IEA
                     			go	GO:0001882	nucleoside binding	07UJ4	bactNOG	11	3.06302846003e-81	21	285	5	276	11	100	IEA
                     			go	GO:0001882	nucleoside binding	08HMP	bactNOG	36	4.75359887693e-81	43	315	3	277	36	100	IEA
1000565.METUNv1_00364	20	20	go	GO:0001883	purine nucleoside binding	0BJMK	bproNOG	2	6.92995759218e-187	1	462	1	454	2	100	IEA
                     			go	GO:0001883	purine nucleoside binding	0B9TC	bproNOG	6	3.81207122537e-146	2	428	3	429	6	100	IEA
                     			go	GO:0001883	purine nucleoside binding	17555	proNOG	16	2.12268729674e-134	53	499	3	431	16	100	IEA
                     			go	GO:0001883	purine nucleoside binding	16RV1	proNOG	205	4.25739996312e-109	517	781	3	275	205	100	IEA
                     			go	GO:0001883	purine nucleoside binding	1755X	proNOG	13	8.40025975035e-109	4	453	2	447	13	100	IEA
                     			go	GO:0001883	purine nucleoside binding	05D9P	bactNOG	1968	5.10153612776e-100	472	722	4	276	1893	96.2	IEA
                     			go	GO:0001883	purine nucleoside binding	07SGC	bactNOG	309	1.2821207656e-94	193	466	3	350	304	98.4	IEA
                     			go	GO:0001883	purine nucleoside binding	0QM1J	gproNOG	21	6.28581547187e-94	262	530	4	276	21	100	IEA
                     			go	GO:0001883	purine nucleoside binding	0GB8Z	delNOG	65	9.47173825047e-92	30	298	4	276	65	100	IEA
                     			go	GO:0001883	purine nucleoside binding	07QPV	bactNOG	109	5.11078410131e-91	163	432	5	278	105	96.3	IEA
                     			go	GO:0001883	purine nucleoside binding	0HCNY	dproNOG	20	3.02140608055e-89	37	300	5	275	20	100	IEA
                     			go	GO:0001883	purine nucleoside binding	0HG69	dproNOG	41	8.9972854538e-89	56	325	5	275	41	100	IEA
                     			go	GO:0001883	purine nucleoside binding	0BCA6	bproNOG	12	1.1100354572e-88	2	274	1	279	12	100	IEA
                     			go	GO:0001883	purine nucleoside binding	0ND8P	firmNOG	472	8.53493970042e-88	52	315	6	277	448	94.9	IEA
                     			go	GO:0001883	purine nucleoside binding	0EP0E	cloNOG	219	1.18732803766e-82	70	332	6	277	212	96.8	IEA
                     			go	GO:0001883	purine nucleoside binding	0D42W	chloNOG	16	4.12832858376e-82	10	273	4	274	16	100	IEA
                     			go	GO:0001883	purine nucleoside binding	0B8I8	bproNOG	11	1.0062412472e-81	3	274	1	279	11	100	IEA
                     			go	GO:0001883	purine nucleoside binding	00CIX	actNOG	255	2.87646308227e-81	64	345	6	356	244	95.7	IEA
                     			go	GO:0001883	purine nucleoside binding	07UJ4	bactNOG	11	3.06302846003e-81	21	285	5	276	11	100	IEA
                     			go	GO:0001883	purine nucleoside binding	08HMP	bactNOG	36	4.75359887693e-81	43	315	3	277	36	100	IEA
1000565.METUNv1_00364	2	2	go	GO:0001932	regulation of protein phosphorylation	1755X	proNOG	13	8.40025975035e-109	4	453	2	447	2	15.4	IEA
                     			go	GO:0001932	regulation of protein phosphorylation	0HCNY	dproNOG	20	3.02140608055e-89	37	300	5	275	1	5	IEA
1000565.METUNv1_00364	20	20	go	GO:0003824	catalytic activity	0BJMK	bproNOG	2	6.92995759218e-187	1	462	1	454	2	100	IEA
                     			go	GO:0003824	catalytic activity	0B9TC	bproNOG	6	3.81207122537e-146	2	428	3	429	6	100	IEA
                     			go	GO:0003824	catalytic activity	17555	proNOG	16	2.12268729674e-134	53	499	3	431	16	100	IEA
                     			go	GO:0003824	catalytic activity	16RV1	proNOG	205	4.25739996312e-109	517	781	3	275	205	100	IEA
                     			go	GO:0003824	catalytic activity	1755X	proNOG	13	8.40025975035e-109	4	453	2	447	13	100	IEA
                     			go	GO:0003824	catalytic activity	05D9P	bactNOG	1968	5.10153612776e-100	472	722	4	276	1895	96.3	IEA, ISS, IDA
                     			go	GO:0003824	catalytic activity	07SGC	bactNOG	309	1.2821207656e-94	193	466	3	350	305	98.7	IEA
                     			go	GO:0003824	catalytic activity	0QM1J	gproNOG	21	6.28581547187e-94	262	530	4	276	21	100	IEA
                     			go	GO:0003824	catalytic activity	0GB8Z	delNOG	65	9.47173825047e-92	30	298	4	276	65	100	IEA
                     			go	GO:0003824	catalytic activity	07QPV	bactNOG	109	5.11078410131e-91	163	432	5	278	105	96.3	IEA
                     			go	GO:0003824	catalytic activity	0HCNY	dproNOG	20	3.02140608055e-89	37	300	5	275	20	100	IEA
                     			go	GO:0003824	catalytic activity	0HG69	dproNOG	41	8.9972854538e-89	56	325	5	275	41	100	IEA
                     			go	GO:0003824	catalytic activity	0BCA6	bproNOG	12	1.1100354572e-88	2	274	1	279	12	100	IEA
                     			go	GO:0003824	catalytic activity	0ND8P	firmNOG	472	8.53493970042e-88	52	315	6	277	448	94.9	IEA, ISS, IDA
                     			go	GO:0003824	catalytic activity	0EP0E	cloNOG	219	1.18732803766e-82	70	332	6	277	212	96.8	IEA
                     			go	GO:0003824	catalytic activity	0D42W	chloNOG	16	4.12832858376e-82	10	273	4	274	16	100	IEA
                     			go	GO:0003824	catalytic activity	0B8I8	bproNOG	11	1.0062412472e-81	3	274	1	279	11	100	IEA
                     			go	GO:0003824	catalytic activity	00CIX	actNOG	255	2.87646308227e-81	64	345	6	356	244	95.7	IEA, IDA
                     			go	GO:0003824	catalytic activity	07UJ4	bactNOG	11	3.06302846003e-81	21	285	5	276	11	100	IEA

(...)

```