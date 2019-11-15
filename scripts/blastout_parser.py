import scripts.pkg_settings as pkg_settings 

class BlastoutParser(object):
    def __init__(self):
        self.path = None
        self.headers = scr
        self.hits = []
        self.best_hits = {}

    def load_from_file(self, path_to_blastout):
        self.path = path_to_blastout
        with open( self.path, 'r' ) as blastout:
            for line in blastout:
                self.hits.append([x.strip() for x in line.split("\t")])

    def best_blast_hit (self, bitcol=7):
        for line in self.best_hits:
            label = spl[0]
            bit = float(spl[bitcol])
            if label in best.keys():
                if float(best[label][bitcol]) < bit:
                    best[label] = spl
            else:
                best[label] = spl

        return best

    '''
    def parse_spdb_blastout(sp_fasta, blastout, log_inst=None):
        path_to_spdb = os.path.abspath(sp_fasta)
        sp_fasta = pkl_fasta_in_out(sp_fasta,seq_type="prot",contig_info=False)
        ids = {}

        ldict = []
        for entry in sp_fasta:
            spl = entry.label.split(" ",1)
            newrow = {
                'accession': spl[0],
                'description': spl[1]
                }
            ldict.append(newrow)
        sp_fasta = None #free

        lib = pd.DataFrame(ldict)
        lib = lib.set_index('accession')

        output = []
        with open(blastout, 'r') as b:
            for line in b:
                spl = line.split("\t")
                acc = line.split("\t")[1]
                contig = line.split("\t")[0]
                evalue = line.split("\t")[10]

                try:
                    hit = lib.loc[acc].description
                    output.append("{}\t{}\t{}\t{}".format(contig, acc, hit, evalue))
                except:
                    if log_inst is not None:
                        log_inst.critical("Accession ({}) not found in swissprot-style database, {}. It is likely that you specified a different version of that database than that used to BLAST against.".format(acc, path_to_spdb))
                        raise ValueError
                    else:
                        print "ERROR"

            return output
    '''