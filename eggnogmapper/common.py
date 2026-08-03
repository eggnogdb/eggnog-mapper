##

import sys
import os
import bz2
import gzip
import shutil
import errno
import struct
import tempfile
import time
from os.path import join as pjoin
from os.path import exists as pexists
from subprocess import Popen, PIPE, run, CalledProcessError

from .utils import colorify
from .emapperException import EmapperException

try:
    from .version import __VERSION__
    from .version import __DB_VERSION__
    from .version import __NOVEL_FAMS_DB_VERSION__
except ImportError:
    __VERSION__ = 'unknown'
    __DB_VERSION__ = 'unknown'
    __NOVEL_FAMS_DB_VERSION__ = 'unknown'    

# Multiprocessing start method
# check
# https://docs.python.org/3/library/multiprocessing.html
# section "Contexts and start methods"
MP_START_METHOD_SPAWN = "spawn" # Available on Unix and Windows. The default on Windows and macOS.
MP_START_METHOD_FORK = "fork" # Available on Unix only. The default on Unix.
MP_START_METHOD_FORKSERVER = "forkserver" # Available on Unix platforms which support passing file descriptors over Unix pipes.
MP_START_METHOD_DEFAULT = MP_START_METHOD_SPAWN
    
# Input types
ITYPE_CDS = "CDS"
ITYPE_PROTS = "proteins"
ITYPE_GENOME = "genome"
ITYPE_META = "metagenome"

TIMEOUT_LOAD_SERVER = 10

BASE_PATH = os.path.abspath(os.path.split(os.path.abspath(__file__))[0] + '/..')

# External binary tools. Resolved via shutil.which at import time; each
# constant is either an absolute path to the executable or None when the
# tool is not on PATH. Call sites that need a tool should call
# require_tool() to fail with an actionable error before invoking it.
#
# v3.0 removed the bundled `eggnogmapper/bin/` directory (~150 MB of
# Linux x86_64 binaries). Install the tools via your package manager:
#   conda install -c bioconda diamond hmmer mmseqs2 prodigal
HMMSEARCH = shutil.which('hmmsearch')
HMMSCAN = shutil.which('hmmscan')
HMMSTAT = shutil.which('hmmstat')
HMMPGMD = shutil.which('hmmpgmd')
PHMMER = shutil.which('phmmer')
HMMFETCH = shutil.which('hmmfetch')
HMMPRESS = shutil.which('hmmpress')
ESL_REFORMAT = shutil.which('esl-reformat')
DIAMOND = shutil.which('diamond')
MMSEQS2 = shutil.which('mmseqs')
PRODIGAL = shutil.which('prodigal')


def require_tool(path, name, conda_pkg=None):
    """Raise EmapperException if the given tool path is not usable.

    Call this just before invoking a tool, e.g.
        require_tool(DIAMOND, 'diamond', 'diamond')
    so the user gets a clear, actionable error message instead of a
    generic FileNotFoundError from subprocess.
    """
    if path and os.access(path, os.X_OK):
        return path
    pkg = conda_pkg or name
    raise EmapperException(
        f"'{name}' executable not found on PATH. "
        f"Install via 'conda install -c bioconda {pkg}' "
        "(or your system package manager)."
    )

DATA_PATH = pjoin(BASE_PATH, "data")
def get_data_path(): return DATA_PATH
def get_eggnogdb_file(): return pjoin(DATA_PATH, "eggnog.db")
def get_novel_fams_db_file(): return pjoin(DATA_PATH, "novel_fams.pkl")
def get_ncbitaxadb_file(): return pjoin(DATA_PATH, "eggnog.taxa.db")
def get_eggnog_mmseqs_dbpath(): return pjoin(DATA_PATH, "mmseqs")
def get_eggnog_mmseqs_db(): return pjoin(DATA_PATH, "mmseqs", "mmseqs.db")
def get_pfam_dbpath(): return pjoin(DATA_PATH, "pfam")
def get_pfam_db(): return pjoin(DATA_PATH, "pfam", "Pfam-A.hmm")
def get_pfam_clans_file(): return pjoin(DATA_PATH, "pfam", "Pfam-A.clans.tsv.gz")
def get_hmmer_dbpath(dbname): return pjoin(DATA_PATH, 'hmmer', dbname, dbname+".hmm")
def get_hmmer_base_dbpath(dbname): return pjoin(DATA_PATH, 'hmmer', dbname)
def get_hmmdb_path(): return pjoin(DATA_PATH, "hmmer")
def get_OG_fasta_path(dbname, og): return pjoin(DATA_PATH, 'hmmer', dbname, f"{og}.fa")
def get_hmmer_databases(): return os.listdir(get_hmmdb_path()) if os.path.isdir(os.path.realpath(get_hmmdb_path())) else []

def get_tax_scopes_path(): return pjoin(BASE_PATH, "eggnogmapper", "annotation", "tax_scopes")

def get_oglevels_file(): return pjoin(DATA_PATH, "og2level.tsv.gz")

def set_data_path(data_path):
    global DATA_PATH
    DATA_PATH = existing_dir(data_path)

##
def cleanup_og_name(name):
    import re
    # v7 OG names (e.g. "UNK.278C@131567|A-1") should pass through unchanged
    if "@" in name and not re.match(r'.*(ENOG41|COG|KOG|arCOG)', name):
        return name
    # v5: names in the hmm databases are sometimes not clean eggnog OG names
    m = re.search('.*((ENOG41|COG|KOG|arCOG)\w+)\.', name)
    if m:
        name = m.groups()[0]
    name = re.sub("^ENOG41", "", name)
    return name

def get_call_info():
    text = []
    text.append('## ' + time.ctime())
    text.append('## ' + get_version())
    text.append('## ' + ' '.join(sys.argv))
    text.append('##')
    return '\n'.join(text)


def get_full_version_info():

    version = get_version()

    db_version = None
    try:
        db_version = get_db_version()
    except Exception as e:
        print(colorify(f"There was an error retrieving eggnog-mapper DB data: {e}", 'red'))
        print(colorify("Maybe you need to run download_eggnog_data.py", 'white'))
        db_version = "unknown"

    if db_version is not None:
        version = f"{version} / eggNOG DB version: {db_version}"

    dmnd_version = get_diamond_version()
    if dmnd_version is not None:
        version = f"{version} / {dmnd_version}"

    mmseqs_version = get_mmseqs_version()
    if mmseqs_version is not None:
        version = f"{version} / {mmseqs_version}"

    exp_novel_fams_db_version = __NOVEL_FAMS_DB_VERSION__
    if exp_novel_fams_db_version is not None:
        version = f"{version} / Compatible novel families DB version: {exp_novel_fams_db_version}"    
        
    return version


def get_version():
    _version = ''
    try:
        # If on a git repository and tags are available
        # Use a tag based code (e.g. 3.1.1b2-8-gb2d12f4)
        p = Popen(["git", "describe", "--tags"], stdout=PIPE, stderr=PIPE, cwd=BASE_PATH)
        out, err = p.communicate()
    except OSError as e:
        if e.errno == errno.ENOENT:
            # Git not installed or package not under git
            pass
        else:
            raise
    else:
        if p.returncode == 0:
            _version += f"-{bytes.decode(out).rstrip()}"
        else:
            # If tags were not available
            # Use a short hash for the current commit (e.g. b2d12f4)
            p = Popen(["git", "rev-parse", "--short", "HEAD"], stdout=PIPE, stderr=PIPE)
            out, err = p.communicate()

            if p.returncode == 0:
                # prefix also with __VERSION__
                # https://github.com/eggnogdb/eggnog-mapper/issues/302
                _version += f"-{__VERSION__}-{bytes.decode(out).rstrip()}"
    if _version != '':
        version = 'emapper' + _version
    else:
        version = 'emapper-'+ __VERSION__
        
    return version

def get_db_version():
    from .annotation.db_sqlite import get_fresh_eggnog_db
    eggnog_db = get_fresh_eggnog_db()
    return eggnog_db.get_db_version()

def get_diamond_version():
    dmnd_version = None
    cmd = f"{DIAMOND} --version"
    try:
        completed_process = run(cmd, capture_output=True, check=True, shell=True)

        if completed_process is not None:
            dmnd_version = f"Diamond version found: {completed_process.stdout.decode('utf-8').strip()}"
            
    except CalledProcessError as cpe:
        print("Couldn't find diamond: "+cpe.stderr.decode("utf-8").strip().split("\n")[-1], file = sys.stderr)
        dmnd_version = "Diamond was not found."
    
    return dmnd_version

def get_mmseqs_version():
    mmseqs_version = None
    cmd = f"{MMSEQS2} version"
    try:
        completed_process = run(cmd, capture_output=True, check=True, shell=True)

        if completed_process is not None:
            mmseqs_version = f"MMseqs2 version found: {completed_process.stdout.decode('utf-8').strip()}"
            
    except CalledProcessError as cpe:
        print("Couldn't find MMseqs2: "+cpe.stderr.decode("utf-8").strip().split("\n")[-1], file = sys.stderr)
        mmseqs_version = "MMseqs2 was not found."

    return mmseqs_version


EGGNOG_DATABASES = {k:51700+(i*2) for i, k in enumerate('NOG,aciNOG,acidNOG,acoNOG,actNOG,agaNOG,agarNOG,apiNOG,aproNOG,aquNOG,arNOG,arcNOG,artNOG,arthNOG,ascNOG,aveNOG,bacNOG,bactNOG,bacteNOG,basNOG,bctoNOG,biNOG,bproNOG,braNOG,carNOG,chaNOG,chlNOG,chlaNOG,chloNOG,chlorNOG,chloroNOG,chorNOG,chrNOG,cloNOG,cocNOG,creNOG,cryNOG,cyaNOG,cytNOG,debNOG,defNOG,dehNOG,deiNOG,delNOG,dipNOG,dotNOG,dproNOG,droNOG,eproNOG,eryNOG,euNOG,eurNOG,euroNOG,eurotNOG,fiNOG,firmNOG,flaNOG,fuNOG,fusoNOG,gproNOG,haeNOG,halNOG,homNOG,hymNOG,hypNOG,inNOG,kinNOG,lepNOG,lilNOG,maNOG,magNOG,meNOG,metNOG,methNOG,methaNOG,necNOG,negNOG,nemNOG,onyNOG,opiNOG,perNOG,plaNOG,pleNOG,poaNOG,prNOG,proNOG,rhaNOG,roNOG,sacNOG,saccNOG,sorNOG,sordNOG,sphNOG,spiNOG,spriNOG,strNOG,synNOG,tenNOG,thaNOG,theNOG,therNOG,thermNOG,treNOG,veNOG,verNOG,verrNOG,virNOG,viruses'.split(','))}
EGGNOG_DATABASES.update({'euk':51400, 'bact':51500, 'arch':51600})

def get_db_info(level):
    return (get_hmmer_dbpath(level), EGGNOG_DATABASES[level])

def get_db_present(level):
    dbpath, port = get_db_info(level)
    db_present = all([pexists(dbpath + "." + ext) for ext in 'h3f h3i h3m h3p idmap'.split()])
    return db_present

def get_citation(addons=['hmmer'], db_version=None):
    _db = f"; eggNOG DB version {db_version}" if db_version else ""
    EXAMPLE = """
e.g. Functional annotation was performed using eggNOG-mapper (version %s%s) [1]
""" % (get_version(), _db)
    
    CITATION = """
================================================================================
CITATION:
If you use this software, please cite:

[1] eggNOG-mapper v2: functional annotation, orthology assignments, and domain 
      prediction at the metagenomic scale. Carlos P. Cantalapiedra, 
      Ana Hernandez-Plaza, Ivica Letunic, Peer Bork, Jaime Huerta-Cepas. 2021.
      Molecular Biology and Evolution, msab293, https://doi.org/10.1093/molbev/msab293
"""

    if 'novel_fams' in addons:
        CITATION += """
[2] Functional and evolutionary significance of unknown genes from uncultivated taxa. 
        Álvaro Rodríguez del Río, Joaquín Giner-Lamia, Carlos P. Cantalapiedra, 
        Jorge Botas, Ziqi Deng, Ana Hernández-Plaza, Lucas Paoli, Thomas S.B. Schmidt, 
        Shinichi Sunagawa, Peer Bork, Luis Pedro Coelho, Jaime Huerta-Cepas. 
        2022. bioRxiv 2022.01.26.477801. https://doi.org/10.1101/2022.01.26.477801
"""
        EXAMPLE += " based on novel families from [2]."
    else:
        CITATION += """
[2] eggNOG v7: phylogeny-based orthology predictions and functional
      annotations. Ana Hernández-Plaza, Ziqi Deng, Fabian Robledo-Yagüe,
      Damian Szklarczyk, Christian von Mering, Peer Bork, Jaime Huerta-Cepas.
      Nucleic Acids Research, Volume 54, Issue D1, 6 January 2026, Pages
      D402-D408. https://doi.org/10.1093/nar/gkaf1249
"""
        EXAMPLE += " based on eggNOG orthology data [2]."

    if 'hmmer' in addons:
        CITATION += """
[3] Accelerated Profile HMM Searches. 
       Eddy SR. 2011. PLoS Comput. Biol. 7:e1002195. 
"""
    elif 'diamond' in addons:
        CITATION += """
[3] Sensitive protein alignments at tree-of-life scale using DIAMOND.
       Buchfink B, Reuter K, Drost HG. 2021.
       Nature Methods 18, 366–368 (2021). https://doi.org/10.1038/s41592-021-01101-x
"""
    elif 'mmseqs' in addons:
        CITATION += """
[3] MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets.
       Steinegger M & Söding J. 2017. Nat. Biotech. 35, 1026–1028. https://doi.org/10.1038/nbt.3988
"""

    EXAMPLE += " Sequence searches were performed using [3]."

    if 'prodigal' in addons:
        CITATION += """
[4] Prodigal: prokaryotic gene recognition and translation initiation site identification.
       Hyatt et al. 2010. BMC Bioinformatics 11, 119. https://doi.org/10.1186/1471-2105-11-119.
"""
        EXAMPLE += " Gene prediction was performed using [4]."

    CITATION += EXAMPLE
    CITATION += """

================================================================================
"""

    return CITATION


LICENSE = """
LICENSE:
[1] eggNOG-mapper is free software distributed under the GPL v2 terms.
Built-in databases (e.g. eggNOG data) might be subjected to different licensing.

[2] eggNOG v5.0 data are distributed under the terms of the Creative Commons Non-Commercial Attribution
License (http://creativecommons.org/licenses/by-nc/4.0/), which permits
unrestricted reuse, distribution, and reproduction in any medium, provided the
original work is properly cited.
"""

def detect_compression(path: str) -> str:
    """Return 'gz', 'bz2', or '' based on magic bytes."""
    try:
        with open(path, 'rb') as fh:
            magic = fh.read(3)
    except OSError:
        return ''
    if magic[:2] == b'\x1f\x8b':
        return 'gz'
    if magic[:3] == b'BZh':
        return 'bz2'
    return ''


def gopen(fname):
    """Open a file for text reading, auto-detecting gzip/bzip2 by magic bytes."""
    ctype = detect_compression(fname)
    if ctype == 'gz':
        return gzip.open(fname, 'rt')
    if ctype == 'bz2':
        return bz2.open(fname, 'rt')
    return open(fname)


def decompress_input(path: str, temp_dir: str) -> tuple:
    """Decompress a compressed FASTA to a temp file; return (path, tmp_path_or_None).

    If the file is not compressed, returns (path, None) with no I/O.
    The caller is responsible for deleting tmp_path when it is not None.
    Supported formats: gzip, bzip2 (detected by magic bytes, not extension).
    """
    ctype = detect_compression(path)
    if not ctype:
        return path, None

    opener = gzip.open if ctype == 'gz' else bz2.open
    fd, tmp_path = tempfile.mkstemp(suffix='.fa', prefix='emapper_input_', dir=temp_dir)
    os.close(fd)
    try:
        with opener(path, 'rb') as fin, open(tmp_path, 'wb') as fout:
            shutil.copyfileobj(fin, fout)
    except Exception:
        silent_rm(tmp_path)
        raise
    return tmp_path, tmp_path


def resolve_input_for_tool(path: str, temp_dir: str, streams_gzip: bool = False) -> tuple:
    """Resolve a (possibly compressed) input for feeding to an external tool.

    Compression is detected by magic bytes (not extension). Gzip is passed
    through untouched when ``streams_gzip`` is True (the tool reads ``.gz``
    natively and streams it, so no temp file is written); otherwise — and
    always for bzip2, which no search tool streams — the input is decompressed
    to a temp file. Returns ``(usable_path, tmp_path_or_None)``; the caller
    must delete ``tmp_path`` when it is not None.
    """
    ctype = detect_compression(path)
    if not ctype:
        return path, None
    if ctype == 'gz' and streams_gzip:
        return path, None
    return decompress_input(path, temp_dir)


def gz_uncompressed_size(path: str):
    """Best-effort uncompressed size in bytes, for input-size heuristics.

    For gzip, reads the 4-byte ISIZE trailer (modulo 2**32 — exact for
    inputs below 4 GiB, which is all the size heuristics here care about).
    For any other input, returns the on-disk size. Returns None on error.
    """
    if detect_compression(path) == 'gz':
        try:
            with open(path, 'rb') as fh:
                fh.seek(-4, os.SEEK_END)
                return struct.unpack('<I', fh.read(4))[0]
        except (OSError, struct.error):
            pass
    try:
        return os.path.getsize(path)
    except OSError:
        return None


def sort_seeds_file(in_file, out_file, temp_dir=None, parallel=None):
    """Sort a seed_orthologs file by seed ortholog id (column 2).

    Comment lines are stripped (annotation's ``parse_seeds`` ignores them
    anyway) and the remaining data lines are sorted numerically by the seed
    id (column 2), with the query id (column 1) as a deterministic tiebreak
    so the order is stable across runs — required for ``--resume`` to line
    the sorted seeds up with the previously written annotations.

    Uses GNU ``sort`` (external merge sort) so files far larger than RAM sort
    in bounded memory. ``LC_ALL=C`` disables locale collation for speed.

    Args:
        in_file: Path to the source seed_orthologs file.
        out_file: Path to write the sorted, comment-stripped data lines to.
        temp_dir: Spill directory for sort (needs ~1x the file size free).
        parallel: Number of sort threads (typically the CPU count).

    Returns:
        ``out_file``.
    """
    import shlex
    par = f"--parallel={int(parallel)} " if parallel else ""
    tmp = f"-T {shlex.quote(temp_dir)} " if temp_dir else ""
    # Sort to a temp file and atomically rename, so a killed sort never leaves
    # a partial out_file behind — a present out_file always means a completed
    # sort (which --resume relies on to safely reuse it).
    tmp_out = out_file + ".tmp"
    # $'\t' yields a real tab under bash; -k2,2n = numeric on seed id.
    cmd = (
        f"grep -v '^#' {shlex.quote(in_file)} | "
        f"LC_ALL=C sort -t $'\\t' -k2,2n -k1,1 -S 60% {par}{tmp}"
        f"> {shlex.quote(tmp_out)}"
    )
    try:
        run(cmd, shell=True, check=True, executable="/bin/bash")
        os.replace(tmp_out, out_file)
    except CalledProcessError as exc:
        silent_rm(tmp_out)
        raise EmapperException(
            f"Failed to sort seed orthologs file {in_file}: {exc}"
        ) from exc
    return out_file


def silent_rm(f):
    if pexists(f):
        os.remove(f)

def silent_cp(f, dst):
    if pexists(f):
        shutil.copy(f, dst)


def existing_file(fname):
    fname = os.path.realpath(fname)
    if os.path.isfile(fname):
        return fname
    else:
        raise TypeError('not a valid file "%s"' %fname)

def existing_dir(dname):
    dname = os.path.realpath(dname)
    if os.path.isdir(dname):
        return dname
    else:
        raise TypeError('not a valid directory "%s"' %dname)

## END
