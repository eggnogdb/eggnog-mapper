import os

def ask(string, valid_values, default=-1, case_sensitive=False):
    """ Asks for a keyborad answer """
    v = None
    if not case_sensitive:
        valid_values = [value.lower() for value in valid_values]
    while v not in valid_values:
        v = raw_input("%s [%s]" % (string,','.join(valid_values) ))
        if v == '' and default >= 0:
            v = valid_values[default]
        if not case_sensitive:
            v = v.lower()
    return v

def run(cmd):
    print cmd
    s = os.system(cmd)
    if s != 0:
        raise ValueError("Error executing command")
    

if ask("Download/upgrade HMM database?", ["y", "n"]) == "y":
    cmd = 'wget http://eggnogdb.embl.de/download/eggnog_4.1/hmmdb.euk_bact_arch.tar.gz | tar -zxf'
    run(cmd)
    
    
if ask("Download/upgrade OG annotations file?", ["y", "n"]) == "y":
    cmd = 'wget http://eggnogdb.embl.de/download/eggnog_4.1/all_OG_annotations.tsv.gz -O - |gzip  > all_OG_annotations.tsv'
    run(cmd)

if ask("Upgrade annotations db?", ["y", "n"]) == "y":
    cmd = "sqlite3 annotations.db < update_db.sql"
    run("cat update_db.sql")
    run(cmd)


    
