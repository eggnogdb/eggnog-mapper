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

DROP TABLE IF EXISTS event;
CREATE TABLE event(
       i INTEGER PRIMARY KEY,
       level VARCHAR(16),
       og VARCHAR(16),
       side1 TEXT,
       side2 TEXT);
CREATE INDEX event_level_idx ON event (i, level);

DROP TABLE IF EXISTS member;
CREATE TABLE member(
       name VARCHAR(32) PRIMARY KEY,
       pname VARCHAR(32),
       groups TEXT,
       go TEXT,
       kegg TEXT,
       orthoindex TEXT);

.separator "\t"
.import ogs.tsv og

.separator "\t"
.import ogs_virus.tsv og

.separator "\t"
.import members.tsv member

.separator "\t"
.import events.tsv event
