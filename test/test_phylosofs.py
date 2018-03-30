import phylosofs
import unittest
import subprocess
import filecmp
import os
import shutil

PATH_TMP = os.path.join("test", "tmp")

if not os.path.isdir(PATH_TMP):
    os.mkdir(PATH_TMP)


class Test_PhyloSofS(unittest.TestCase):

    def test_phylosofs(self):
        self.assertEqual(subprocess.call(["python", "phylosofs/phylosofs.py",
                                          "-mode", "P",
                                          "-o", "test/tmp/",
                                          # It fails with "test/tmp"
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
        self.assertFalse(os.path.isdir("test/tmp/bestTopos"))
        self.assertFalse(os.path.isdir("test/tmp/betterTrees"))

    def test_best_topos_and_trees(self):
        self.assertEqual(subprocess.call(["python", "phylosofs/phylosofs.py",
                                          "-mode", "P",
                                          "-o", "test/tmp/",
                                          "-s", "100",
                                          "-inSeq", "dat/JNK3.txt"]), 0)
        self.assertTrue(os.path.isdir("test/tmp/bestTopos"))
        self.assertTrue(os.path.isdir("test/tmp/betterTrees"))
        self.assertGreater(len(os.listdir("test/tmp/bestTopos")), 0)
        self.assertGreater(len(os.listdir("test/tmp/betterTrees")), 0)

    def tearDown(self):
        for root, dirs, files in os.walk(PATH_TMP):
            for f in files:
                os.unlink(os.path.join(root, f))
            for d in dirs:
                shutil.rmtree(os.path.join(root, d))


if __name__ == '__main__':
    unittest.main()
