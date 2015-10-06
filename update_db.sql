DROP TABLE IF EXISTS annotations;

CREATE TABLE annotations(
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
.import all_OG_annotations.tsv annotations


