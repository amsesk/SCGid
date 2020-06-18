# Class to handle slow steps that need only be done once per SCGid run
import re
import os
import inspect
import logging
import shutil
from scgid.modcomm import LoggingEntity, logger_name_gen, get_head
from scgid.library import subprocessP, gff3_to_fasta, is_fasta
import scgid.pkg_settings as pkg_settings
from scgid.error import ModuleError, FatalError, ArgumentError


class OutputGeneratorFunctionNotImplementedError(ModuleError):

    def __init__(self, arg, error_catch=True):
        super().__init__()
        self.msg = f"No implemented function to generate output for arg `{arg}`. Since it cannot be located, it must be specified from the CLI."
        if error_catch:
            self.catch()


class MultipleReusableOutputMatchesError(ModuleError):

    def __init__(self, arg, matches, error_catch=True):
        super().__init__()
        self.msg = f"Found multiple files matching pattern for argument `{arg}`. Specify preference in command line arugment. [{','.join(matches)}]"
        if error_catch:
            self.catch()


class MalformedDatabaseFasta(ModuleError):

    def __init__(self, db_path, error_catch=True):
        super().__init__()
        self.msg = f"Malformed database FASTA at `{db_path}`"
        if error_catch:
            self.catch()


class ReusableOutput:

    def __init__(self, arg, pattern, genfunc=None, genfunc_args=None):
        self.arg = arg
        self.caller = inspect.stack()[1][0].f_locals["self"]
        self.re_pattern = re.compile(pattern)
        self.needs_doing = False
        self.genfunc = genfunc
        self.genfunc_args = genfunc_args
        self.head = get_head()

    def try_reuse(self):
        rundir = self.head.config.rundir
        reusable = False
        for item in [os.path.join(rundir, x) for x in os.listdir(rundir)]:
            if os.path.isdir(item):
                matches = [x for x in os.listdir(
                    item) if re.match(self.re_pattern, x)]
                if len(matches) == 0:

                    # reusable already set to False
                    continue

                elif len(matches) == 1:
                    updated_arg = os.path.join(item, matches[0])
                    setattr(self.caller.config, self.arg, updated_arg)
                    reusable = True
                    break
                else:
                    MultipleReusableOutputMatchesError(self.arg, matches)
        return reusable

    def update_from_config(self):
        for k, v in self.genfunc_args.items():
            if v is None:
                self.genfunc_args[k] = self.head.config.get(k)

    def generate(self):
        self.update_from_config()
        return self.genfunc(**self.genfunc_args)


class ReusableOutputManager(LoggingEntity):

    def __init__(self, *reusable):
        self.log = logging.getLogger(logger_name_gen())
        self.head = get_head()
        self.reusable = []

    def populate(self, *reusable):
        self.reusable = list(reusable)

    def add(self, reusable_output):
        self.reusable.append(reusable_output)

    def check(self):
        for r in self.reusable:
            if self.head.config.get(r.arg) is None:
                if not r.try_reuse():
                    self.log.info(f"No match found for required file specified by `{r.arg}`.")

                    # Covers the case where instead of invoking a function to
                    # generate the output, we just want to return a FatalError
                    if isinstance(r.genfunc, FatalError):
                        return r.genfunc()

                    r.needs_doing = True

                    # r.genfunc(**r.genfunc_args)

                    if r.genfunc is None:
                        return OutputGeneratorFunctionNotImplementedError(r.arg)
                else:
                    self.log.info(f"Found matching file for missing argument `{r.arg}` at `{get_head().config.get(r.arg)}`")
                    r.needs_doing = False
            else:
                continue

    def generate_outputs(self):
        for r in self.reusable:
            if r.needs_doing:
                path = r.generate()
                setattr(self.head.config, r.arg, path)
            else:
                pass


def augustus_predict(prefix, nucl, augustus_sp, outpath):

    # Raise error if -sp|--augustus_species is unset and this function is called.
    # Unset means it is None (based on argparse default)
    if augustus_sp is None:
        raise ArgumentError(
            "-sp|--augustus_species is required for gene prediction")

    log = get_head().logger
    gff3_fname = f"{prefix}.aug.out.gff3"
    prot_fname = f"{prefix}.aug.out.fasta"

    cmd = ["augustus", f"--species={augustus_sp}", "--gff3=on", nucl, "--uniqueGeneId=true"]

    log.info(' '.join(cmd))
    out = subprocessP(cmd, log)

    log.info(f"Writing augustus output to {gff3_fname}")
    with open(gff3_fname, 'wb') as f:
        f.write(out)

    gff3_to_fasta(gff3_fname, prot_fname)
    log.info(f"gff3 converted to fasta and writen to `{os.path.abspath(prot_fname)}`")

    # Move gff3/protein fasta to default locations if finished successfully in
    # temp
    shutil.copyfile(gff3_fname, f"../{gff3_fname}")
    shutil.copyfile(prot_fname, f"../{prot_fname}")

    return outpath


def verify_blastdb(db):
    logger = get_head().logger
    db_check_name = '.'.join([db, 'phr'])
    if os.path.isfile(db_check_name):
        logger.info(f"Found blastdb associated with database fasta at `{db}`")
    else:
        logger.info(f"Unable to find blastdb associated with `{db}`")
        logger.info(f"Checking format of supplied database FASTA...")
        if is_fasta(db):
            new_dbname = os.path.basename(db)
            print(new_dbname)
            logger.info(f"Database FASTA found at `{db}`")
            logger.info(f"Making blastdb from `{db}``")
            cmd = ["makeblastdb", "-in", db, "-dbtype",
                   "prot", "-title", new_dbname, "-out", db]
            logger.info(' '.join(cmd))
            subprocessP(cmd, logger)
        else:
            return MalformedDatabaseFasta(db)

    return 0


def protein_blast(prot, db, evalue, cpus, outpath):
    logger = get_head().logger
    logger.info(f"Blasting predicted proteins against the swissprot database located at `{db}`")
    verify_blastdb(db)

    cmd = ["blastp", "-query", prot, "-max_target_seqs", "1", "-evalue", evalue, "-db",
           db, "-outfmt", pkg_settings.BLAST_OUTFMT_STR, "-out", outpath, "-num_threads", cpus]
    logger.info(' '.join(cmd))
    subprocessP(cmd, logger)

    # Move blastout to default locations if finished successfully in temp
    shutil.copyfile(outpath, f"../{outpath}")

    return outpath


def nucleotide_blast(nucl, db, evalue, cpus, outpath):
    logger = get_head().logger
    logger.info(f"Blasting predicted proteins against the swissprot database located at `{db}`")
    if not db == "nt":
        verify_blastdb(db)

    cmd = ["blastn", "-query", nucl, "-max_target_seqs", "1", "-evalue", evalue, "-db",
           db, "-outfmt", pkg_settings.BLAST_OUTFMT_STR, "-out", outpath, "-num_threads", cpus]
    logger.info(' '.join(cmd))
    subprocessP(cmd, logger)

    # Move blastout to default locations if finished successfully in temp
    shutil.copyfile(outpath, f"../{outpath}")

    return outpath
