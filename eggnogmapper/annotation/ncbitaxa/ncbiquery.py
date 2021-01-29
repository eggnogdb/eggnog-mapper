##
## CPCantalapiedra 2020

# This is an adaptation of ete3 NCBITaxa module
# to embed it directly within eggnog-mapper

import os
import pickle
import sqlite3

from ...emapperException import EmapperException
from ...common import get_ncbitaxadb_file

class NCBITaxa(object):
    """
    Provides a local transparent connector to the NCBI taxonomy database.
    """
    dbfile = None
    db = None
    
    def __init__(self, dbfile=None):

        if dbfile is None:
            self.dbfile = get_ncbitaxadb_file()
        else:
            self.dbfile = dbfile
            
        if not os.path.exists(self.dbfile):
            raise EmapperException(f'Could not find NCBI database {self.dbfile}.')
        
        self.db = sqlite3.connect(self.dbfile)
        
        return
        
    ##
    def get_taxid_translator(self, taxids, try_synonyms=True):
        """Given a list of taxids, returns a dictionary with their corresponding
        scientific names.
        """

        all_ids = set(map(int, taxids))
        all_ids.discard(None)
        all_ids.discard("")
        
        cmd = "select taxid, spname FROM species WHERE taxid IN (%s);" % ",".join(set(map(str, all_ids)))
        
        result = self.db.execute(cmd)
        
        id2name = {}
        for tax, spname in result.fetchall():
            id2name[tax] = spname

        # any taxid without translation? lets tray in the merged table
        if len(all_ids) != len(id2name) and try_synonyms:
            not_found_taxids = all_ids - set(id2name.keys())
            taxids, old2new = self._translate_merged(not_found_taxids)
            new2old = {v: k for k,v in old2new.items()}

            if old2new:
                query = ','.join(['"%s"' %v for v in new2old])
                cmd = "select taxid, spname FROM species WHERE taxid IN (%s);" %query
                
                result = self.db.execute(cmd)
                for tax, spname in result.fetchall():
                    id2name[new2old[tax]] = spname

        return id2name


    ##
    def get_name_translator(self, names):
        """
        Given a list of taxid scientific names, returns a dictionary translating them into their corresponding taxids.
        Exact name match is required for translation.
        """
        name2id = {}
        name2origname = {}
        for n in names:
            name2origname[n.lower()] = n

        names = set(name2origname.keys())

        query = ','.join([f'"{n}"' for n in names])
        cmd = 'select spname, taxid from species where spname IN (%s)' %query
        
        result = self.db.execute(cmd)
        for sp, taxid in result.fetchall():
            oname = name2origname[sp.lower()]
            name2id.setdefault(oname, []).append(taxid)
        missing =  names - set([n.lower() for n in name2id.keys()])
        if missing:
            query = ','.join(['"%s"' %n for n in missing])
            cmd = 'select spname, taxid from synonym where spname IN (%s)' %query
            
            result = self.db.execute(cmd)
            for sp, taxid in result.fetchall():
                oname = name2origname[sp.lower()]
                name2id.setdefault(oname, []).append(taxid)
                
        return name2id
    
    
    ##
    def get_descendant_taxa(self, parent, intermediate_nodes=False):
        """
        given a parent taxid or scientific species name, returns a list of all its descendants taxids.
        If intermediate_nodes is set to True, internal nodes will also be dumped.
        """
        try:
            taxid = int(parent)
        except ValueError:
            try:
                taxid = self.get_name_translator([parent])[parent][0]
            except KeyError:
                raise ValueError('%s not found!' %parent)

        # checks if taxid is a deprecated one, and converts into the right one.
        _, conversion = self._translate_merged([taxid]) #try to find taxid in synonyms table
        if conversion:
            taxid = conversion[taxid]

        with open(self.dbfile+".traverse.pkl", "rb") as CACHED_TRAVERSE:
            prepostorder = pickle.load(CACHED_TRAVERSE)
        descendants = {}
        found = 0
        for tid in prepostorder:
            if tid == taxid:
                found += 1
            elif found == 1:
                descendants[tid] = descendants.get(tid, 0) + 1
            elif found == 2:
                break

        if not found:
            raise ValueError("taxid not found:%s" %taxid)
        elif found == 1:
            return [taxid]

        if intermediate_nodes:
            return [tid for tid, count in descendants.items()]
        else:
            return [tid for tid, count in descendants.items() if count == 1]

        return
            
    def _translate_merged(self, all_taxids):
        conv_all_taxids = set((list(map(int, all_taxids))))
        cmd = 'select taxid_old, taxid_new FROM merged WHERE taxid_old IN (%s)' %','.join(map(str, all_taxids))
        
        result = self.db.execute(cmd)
        conversion = {}
        for old, new in result.fetchall():
            conv_all_taxids.discard(int(old))
            conv_all_taxids.add(int(new))
            conversion[int(old)] = int(new)
            
        return conv_all_taxids, conversion

## END
