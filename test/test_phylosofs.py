import phylosofs
import unittest
import subprocess
import filecmp
import os

PATH_TMP = os.path.join("test", "tmp")

if not os.path.isdir(PATH_TMP):
    os.mkdir(PATH_TMP)


class Test_PhyloSofS(unittest.TestCase):

    def test_phylosofs(self):
        path_phylosofs = os.path.join("phylosofs", "phylosofs.py")
        path_tmp = os.path.join("test", "tmp")
        path_data = os.path.join("dat", "JNK3.txt")
        command = ["python", path_phylosofs,
                   "-P",
                   "-o", path_tmp,
                   "--inseq", path_data]
        self.assertEqual(subprocess.call(command), 0)
        self.assertTrue(filecmp.cmp(os.path.join("test", "tmp",
                                                 "treeSearch_532_1.txt"),
                                    os.path.join("test", "data",
                                                 "treeSearch_532_1.txt")),
                        "treeSearch_532_1.txt")
        self.assertTrue(filecmp.cmp(os.path.join("test", "tmp",
                                                 "solution_532_1_config0.sum"),
                                    os.path.join("test", "data",
                                                 "solution_532_1_config0.sum")
                                    ),
                        "solution_532_1_config0.sum")
        self.assertTrue(filecmp.cmp(os.path.join("test", "tmp",
                                                 "solution_532_1_config0.info"
                                                 ),
                                    os.path.join("test", "data",
                                                 "solution_532_1_config0.info")
                                    ),
                        "solution_532_1_config0.info")
        self.assertFalse(os.path.isdir(os.path.join("test", "tmp",
                                                    "bestTopos")))
        self.assertFalse(os.path.isdir(os.path.join("test", "tmp",
                                                    "betterTrees")))

    def test_best_topos_and_trees(self):
        path_phylosofs = os.path.join("phylosofs", "phylosofs.py")
        path_tmp = os.path.join("test", "tmp")
        path_data = os.path.join("dat", "JNK3.txt")
        command = ["python", path_phylosofs,
                   "-P",
                   "-o", path_tmp,
                   "-s", "100",
                   "--inseq", path_data]
        self.assertEqual(subprocess.call(command), 0)
        self.assertTrue(os.path.isdir(os.path.join("test", "tmp",
                                                   "bestTopos")))
        self.assertTrue(os.path.isdir(os.path.join("test", "tmp",
                                                   "betterTrees")))
        self.assertGreater(len(os.listdir(os.path.join("test", "tmp",
                                                       "bestTopos"))), 0)
        self.assertGreater(len(os.listdir(os.path.join("test", "tmp",
                                                       "betterTrees"))), 0)

    def tearDown(self):
        phylosofs.utils.clear_folder(PATH_TMP)


if __name__ == '__main__':
    unittest.main()
