##
## CPCantalapiedra 2019

from os.path import join as pjoin
import shutil
import subprocess
from tempfile import mkdtemp, mkstemp
import uuid

from ...common import MMSEQS2, get_eggnog_mmseqs_db, get_call_info, ITYPE_CDS, ITYPE_PROTS, ITYPE_GENOME, ITYPE_META
from ...emapperException import EmapperException
from ...utils import colorify, translate_cds_to_prots

from ..hmmer.hmmer_seqio import iter_fasta_seqs
from ..diamond.diamond import hit_does_overlap

def create_mmseqs_db(dbprefix, in_fasta):
    cmd = (
        f'{MMSEQS2} createdb {in_fasta} {dbprefix}'
    )

    print(colorify('  '+cmd, 'yellow'))
    try:
        completed_process = subprocess.run(cmd, capture_output=True, check=True, shell=True)
    except subprocess.CalledProcessError as cpe:
        raise EmapperException("Error running mmseqs: "+cpe.stderr.decode("utf-8").strip().split("\n")[-1])
        
    return

def create_mmseqs_index(dbprefix, tmp_dir):
    cmd = (
        f'{MMSEQS2} createindex {dbprefix} {tmp_dir}'
    )

    print(colorify('  '+cmd, 'yellow'))
    try:
        completed_process = subprocess.run(cmd, capture_output=True, check=True, shell=True)
    except subprocess.CalledProcessError as cpe:
        raise EmapperException("Error running mmseqs: "+cpe.stderr.decode("utf-8").strip().split("\n")[-1])
    return

class MMseqs2Searcher:

    name = "mmseqs2"
    
    # Command
    cpu = targetdb = temp_dir = no_file_comments = None
    start_sens = 3
    sens_steps = 3
    final_sens = 7    

    # Filters
    pident_thr = score_thr = evalue_thr = query_cov = subject_cov = None # excluded_taxa = None

    # MMseqs2 options
    sub_mat = None

    in_file = None
    itype = None
    translate = None
    translation_table = None

    # Results
    queries = hits = no_hits = None

    ##
    def __init__(self, args):

        self.itype = args.itype
        self.translate = args.translate

        self.targetdb = args.mmseqs_db if args.mmseqs_db else get_eggnog_mmseqs_db()
        self.cpu = args.cpu
        self.start_sens = args.start_sens
        self.sens_steps = args.sens_steps
        self.final_sens = args.final_sens
        
        self.query_cov = args.query_cover
        self.subject_cov = args.subject_cover

        self.pident_thr = args.pident
        self.evalue_thr = args.mmseqs_evalue
        self.score_thr = args.mmseqs_score
        # self.excluded_taxa = args.excluded_taxa if args.excluded_taxa else None

        self.sub_mat = args.mmseqs_sub_mat
        
        self.temp_dir = args.temp_dir
        self.no_file_comments = args.no_file_comments
        
        return

    ##
    def get_hits(self):
        return self.hits

    ##
    def get_no_hits(self):
        if self.hits is not None:
            hit_queries = set([x[0] for x in self.hits])
            self.queries = set({name for name, seq in iter_fasta_seqs(self.in_file)})
            self.no_hits = set(self.queries).difference(hit_queries)
            
        return self.no_hits
    
    
    ##
    def search(self, in_file, seed_orthologs_file, hits_file = None):
        # MMseqs2Searcher does not use the "hits_file"
        # but we need this wrapper to define the search interface for Emapper
        return self._search(in_file, seed_orthologs_file)

    def _search(self, in_file, seed_orthologs_file):
        if not MMSEQS2:
            raise EmapperException("%s command not found in path" % (MMSEQS2))

        self.in_file = in_file
        
        tempdir = mkdtemp(prefix='emappertmp_mmseqs_', dir=self.temp_dir)
        try:
            querydb = pjoin(tempdir, uuid.uuid4().hex)
            # print(f'Querydb {querydb}')
            resultdb = pjoin(tempdir, uuid.uuid4().hex)
            # print(f'ResultDB {resultdb}')
            bestresultdb = pjoin(tempdir, uuid.uuid4().hex)
            # print(f'BestResultDB {bestresultdb}')
            
            alignmentsdb, cmds = self.run_mmseqs(in_file, tempdir, querydb, self.targetdb, resultdb, bestresultdb)
            self.hits = self.parse_mmseqs(f'{alignmentsdb}.m8')
            self.output_mmseqs(cmds, self.hits, seed_orthologs_file)

        except Exception as e:
            raise e
        finally:
            shutil.rmtree(tempdir)
            
        return
    
    def run_mmseqs(self, fasta_file, tempdir, querydb, targetdb, resultdb, bestresultdb):
        cmds = []

        if self.itype == ITYPE_CDS and self.translate == True:
            handle, query_file = mkstemp(dir = tempdir, text = True)
            translate_cds_to_prots(fasta_file, query_file)
        else:
            query_file = fasta_file
            
        cmd = self.createdb(query_file, querydb)
        cmds.append(cmd)

        cmd = self.search_step(querydb, targetdb, resultdb, tempdir)
        cmds.append(cmd)

        if self.itype == ITYPE_CDS or self.itype == ITYPE_PROTS:
            cmd = self.filterdb_step(resultdb, bestresultdb)
            cmds.append(cmd)
        else:
            bestresultdb = resultdb

        cmd = self.convertalis_step(querydb, targetdb, bestresultdb)
        cmds.append(cmd)

        return bestresultdb, cmds

    def createdb(self, fasta_file, querydb):
        # mmseqs createdb examples/QUERY.fasta queryDB
        cmd = (
            f'{MMSEQS2} createdb {fasta_file} {querydb}'
        )
        if self.itype == ITYPE_PROTS or self.translate == True:
            cmd += ' --dbtype 1' # aas queries (proteins)
        else:
            cmd += ' --dbtype 2' # nts queries (CDS, contig, ...)
            
        print(colorify('  '+cmd, 'yellow'))
        try:
            completed_process = subprocess.run(cmd, capture_output=True, check=True, shell=True)
        except subprocess.CalledProcessError as cpe:
            raise EmapperException("Error running 'mmseqs createdb': "+cpe.stderr.decode("utf-8").strip().split("\n")[-1])
        return cmd

    def search_step(self, querydb, targetdb, resultdb, tempdir):
        # mmseqs search queryDB targetDB resultDB tmp
        cmd = (
            f'{MMSEQS2} search -a true {querydb} {targetdb} {resultdb} {tempdir} '
            f'--start-sens {self.start_sens} --sens-steps {self.sens_steps} -s {self.final_sens} '
            f'--threads {self.cpu}'
        )

        if self.sub_mat is not None: cmd += f' --sub-mat {self.sub_mat}'
        if self.translation_table is not None: cmd += f' --translation-table {self.translation_table}'
        
        print(colorify('  '+cmd, 'yellow'))
        try:
            completed_process = subprocess.run(cmd, capture_output=True, check=True, shell=True)
        except subprocess.CalledProcessError as cpe:
            raise EmapperException("Error running 'mmseqs search': "+cpe.stderr.decode("utf-8").strip().split("\n")[-1])
        return cmd

    ##
    def filterdb_step(self, resultdb, bestresultdb):
        # mmseqs filterdb resultDB bestResultDB --extract-lines 1
        cmd = (
            f'{MMSEQS2} filterdb {resultdb} {bestresultdb} --threads {self.cpu}'
        )
        cmd += " --extract-lines 1 "
        
        print(colorify('  '+cmd, 'yellow'))
        try:
            completed_process = subprocess.run(cmd, capture_output=True, check=True, shell=True)
        except subprocess.CalledProcessError as cpe:
            raise EmapperException("Error running 'mmseqs filterdb': "+cpe.stderr.decode("utf-8").strip().split("\n")[-1])                
        
        return cmd
    
    ##
    def convertalis_step(self, querydb, targetdb, resultdb):
        # mmseqs convertalis queryDB targetDB resultDB resultDB.m8
        cmd = (
            f'{MMSEQS2} convertalis {querydb} {targetdb} {resultdb} {resultdb}.m8 --threads {self.cpu}'
        )
        if self.sub_mat is not None: cmd += f' --sub-mat {self.sub_mat}'
            
        if self.translation_table is not None: cmd += f' --translation-table {self.translation_table}'

        # outfmt
        cmd += f' --format-output "query,target,pident,alnlen,qstart,qend,tstart,tend,evalue,bits,qcov,tcov"'
        print(colorify('  '+cmd, 'yellow'))
        try:
            completed_process = subprocess.run(cmd, capture_output=True, check=True, shell=True)
        except subprocess.CalledProcessError as cpe:
            raise EmapperException("Error running 'mmseqs convertalis': "+cpe.stderr.decode("utf-8").strip().split("\n")[-1])        
        return cmd

    
    ##
    def parse_mmseqs(self, raw_mmseqs_file):
        if self.itype == ITYPE_CDS or self.itype == ITYPE_PROTS:
            return self._parse_mmseqs(raw_mmseqs_file)
        else: #self.itype == ITYPE_GENOME or self.itype == ITYPE_META:
            return self._parse_genepred(raw_mmseqs_file)
    
    ##
    def _parse_mmseqs(self, raw_mmseqs_file):
        hits = []

        visited = set()
        with open(raw_mmseqs_file, 'r') as raw_f:
            for line in raw_f:
                if not line.strip() or line.startswith('#'):
                    continue
                
                fields = list(map(str.strip, line.split('\t')))
                # check fields in convertalis_step()
        
                query = fields[0]

                # only one hit per query
                if query in visited:
                    continue
                visited.add(query)
                
                pident = float(fields[2])
                evalue = float(fields[8])
                score = float(fields[9])
                qcov = float(fields[10]) * 100 # mmseqs uses 0-1 values
                scov = float(fields[11]) * 100 # mmseqs uses 0-1 values

                # note: this could be done with mmmseqs filterdb, but I dont know how to do it in a single step
                if ((self.pident_thr is not None and pident < self.pident_thr) or 
                    (self.evalue_thr is not None and evalue > self.evalue_thr) or
                    (self.score_thr is not None and score < self.score_thr) or
                    (self.query_cov is not None and qcov < self.query_cov) or
                    (self.subject_cov is not None and scov < self.subject_cov)):
                    continue
                
                hit = fields[1]
                length = int(fields[3])
                qstart = int(fields[4])
                qend = int(fields[5])
                sstart = int(fields[6])
                send = int(fields[7])
                
                hits.append([query, hit, evalue, score, qstart, qend, sstart, send, qcov, scov])
                
        return hits

    ##
    def _parse_genepred(self, raw_mmseqs_file):
        hits = []
        curr_query_hits = []
        prev_query = None
        
        
        with open(raw_mmseqs_file, 'r') as raw_f:
            for line in raw_f:
                if not line.strip() or line.startswith('#'):
                    continue
                
                fields = list(map(str.strip, line.split('\t')))
                # check fields in convertalis_step()
                # cmd += f' --format-output "query,target,pident,alnlen,qstart,qend,tstart,tend,evalue,bits,qcov,tcov"'

                pident = float(fields[2])
                evalue = float(fields[8])
                score = float(fields[9])
                qcov = float(fields[10]) * 100 # mmseqs uses 0-1 values
                scov = float(fields[11]) * 100 # mmseqs uses 0-1 values

                # note: this could be done with mmmseqs filterdb, but I dont know how to do it in a single step
                if ((self.pident_thr is not None and pident < self.pident_thr) or 
                    (self.evalue_thr is not None and evalue > self.evalue_thr) or
                    (self.score_thr is not None and score < self.score_thr) or
                    (self.query_cov is not None and qcov < self.query_cov) or
                    (self.subject_cov is not None and scov < self.subject_cov)):
                    continue
                
                query = fields[0]
                hit = fields[1]
                length = int(fields[3])
                qstart = int(fields[4])
                qend = int(fields[5])
                sstart = int(fields[6])
                send = int(fields[7])
                
                hit = [query, hit, evalue, score, qstart, qend, sstart, send, qcov, scov]

                if query == prev_query:
                    if not hit_does_overlap(hit, curr_query_hits):
                        hits.append(hit)
                        curr_query_hits.append(hit)
                else:
                    hits.append(hit)
                    curr_query_hits = [hit]
                    
                prev_query = query    
                
        return hits
    

    ##
    def output_mmseqs(self, cmds, hits, out_file):
        if self.itype == ITYPE_CDS or self.itype == ITYPE_PROTS:
            return self._output_mmseqs(cmds, hits, out_file)
        else: #self.itype == ITYPE_GENOME or self.itype == ITYPE_META:
            return self._output_genepred(cmds, hits, out_file)
    
    ##
    def _output_mmseqs(self, cmds, hits, out_file):
        with open(out_file, 'w') as OUT:

            # comments
            if not self.no_file_comments:
                print(get_call_info(), file=OUT)
                for cmd in cmds:
                    print('#'+cmd, file=OUT)

            # header
            print('#'+"\t".join("qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcov scov".split(" ")), file=OUT)

            # rows
            for line in hits:
                print('\t'.join(map(str, line)), file=OUT)
                
        return

    ##
    def _output_genepred(self, cmds, hits, out_file):
        queries_suffixes = {}
        with open(out_file, 'w') as OUT:

            # comments
            if not self.no_file_comments:
                print(get_call_info(), file=OUT)
                for cmd in cmds:
                    print('#'+cmd, file=OUT)

            # header
            print('#'+"\t".join("qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcov scov".split(" ")), file=OUT)

            # rows
            for line in hits:
                query = line[0]
                target = line[1]
                evalue = line[2]
                score = line[3]
                if query in queries_suffixes:
                    queries_suffixes[query] += 1
                    suffix = queries_suffixes[query]
                else:
                    suffix = 0
                    queries_suffixes[query] = suffix
                
                print('\t'.join(map(str, [f"{query}_{suffix}"] + line[1:])), file=OUT)
                
        return
    
## END
