import unittest
import shutil

importlib = __import__("importlib")
gct = importlib.import_module("scgid.scripts.gct")

class TestGct(unittest.TestCase):
    def test_call(self):
        argdict = {
            "nucl": "/home/aimzez/development/scgid/test_data/stylopage_41_contigs_tester.fasta",
            "prefix": "test",
            "targets": "Fungi",
            "exceptions": None,
            "stringency": "0.05",
            "augustus_sp": None,
            "prot": "/home/aimzez/development/scgid/test_data/stylopage_41_scgid_output/gct/stylopage_41.aug.out.fasta",
            "evalue": None,
            "blastout": "/home/aimzez/development/scgid/test_data/stylopage_41_scgid_output/gct/stylopage_41.spdb.blast.out",
            "cpus": None
        }
        mod = gct.Gct(argdict)
        out = mod.run()
        self.assertEqual(0,0)
        shutil.rmtree("test_scgid_output")

if __name__ == "__main__":
    unittest.main()
