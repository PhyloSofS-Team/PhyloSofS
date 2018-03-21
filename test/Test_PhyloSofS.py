import unittest
import subprocess
import filecmp
import os

class Test_PhyloSofS(unittest.TestCase):

    def test_phylosofs(self):
        self.assertEqual(subprocess.call([  "python", "phylosofs.py",
                                            "-mode", "P",
                                            "-o", "test/tmp/", # It fails with "test/tmp"
                                            "-inSeq", "dat/JNK3.txt"]), 0)
        self.assertTrue(filecmp.cmp("test/tmp/treeSearch_532_1.txt",
                                    "test/data/treeSearch_532_1.txt"),
                                    "treeSearch_532_1.txt")
        self.assertTrue(filecmp.cmp("test/tmp/solution_532_1_config0.sum",
                                    "test/data/solution_532_1_config0.sum"),
                                    "solution_532_1_config0.sum")
        self.assertTrue(filecmp.cmp("test/tmp/solution_532_1_config0.info",
                                    "test/data/solution_532_1_config0.info"),
                                    "solution_532_1_config0.info")
        self.assertTrue(os.path.isdir("test/tmp/bestTopos"))
        self.assertTrue(os.path.isdir("test/tmp/betterTrees"))

if __name__ == '__main__':
    unittest.main()
