-- orthologous groups annotations

DROP TABLE IF EXISTS og;
CREATE TABLE og(
       og VARCHAR(16) PRIMARY KEY,
       level VARCHAR(16),
       nm INTEGER,
       description TEXT,
       COG_categories VARCHAR(8),
       GO_freq TEXT,
       KEGG_freq TEXT,
       SMART_freq TEXT,
       proteins TEXT);

.separator "\t"
.import table_og.tsv og
.import table_og_virus.tsv og

-- speciation events based on eggnog trees

DROP TABLE IF EXISTS event;
CREATE TABLE event(
       i INTEGER PRIMARY KEY,
       level VARCHAR(16),
       og VARCHAR(16),
       side1 TEXT,
       side2 TEXT);
CREATE INDEX event_level_idx ON event (i, level);
.separator "\t"
.import events.tsv event

-- eggnog groups

DROP TABLE IF EXISTS eggnog;
CREATE TABLE eggnog(
       name VARCHAR(32) PRIMARY KEY,
       groups TEXT);
.separator "\t"
.import table_eggnog.tsv eggnog

-- speciation event indexes

DROP TABLE IF EXISTS orthologs;
CREATE TABLE orthologs(
name VARCHAR(32) PRIMARY KEY,
orthoindex TEXT
);
.separator "\t"
.import table_orthologs.tsv orthologs

-- sequence names annotations

DROP TABLE IF EXISTS seq;
CREATE TABLE seq(
name VARCHAR(32) PRIMARY KEY,
pname VARCHAR(32)
);
.separator "\t"
.import table_seq.tsv seq

-- kegg modules

DROP TABLE IF EXISTS kegg;
CREATE TABLE kegg(
name VARCHAR(32) PRIMARY KEY,
ko VARCHAR(16)
);
.separator "\t"
.import table_kegg78_99i99c.tsv kegg

-- gene ontology

DROP TABLE IF EXISTS gene_ontology;
CREATE TABLE gene_ontology(
name VARCHAR(32) PRIMARY KEY,
gos TEXT);
.separator "\t"
.import table_go.tsv gene_ontology

-- bigg reacctions

DROP TABLE IF EXISTS bigg;
CREATE TABLE bigg(
name VARCHAR(32) PRIMARY KEY,
reaction VARCHAR(32)
);
.separator "\t"
.import table_bigg_95i95c.tsv bigg

