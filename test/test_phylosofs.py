import phylosofs
import unittest
import subprocess
import os

TEST_DIR = os.path.dirname(os.path.abspath(__file__))


def compare_files(fpath1, fpath2):
    """
    Compare two files without taking into account newline characters.

    It allows to compare files generated by linux and windows.
    Code from:
    https://stackoverflow.com/questions/40751389/compare-2-files-line-by-line-ignoring-newline-differences
    """
    with open(fpath1, 'r') as file1, open(fpath2, 'r') as file2:
        for linef1, linef2 in zip(file1, file2):
            linef1 = linef1.rstrip('\r\n')
            linef2 = linef2.rstrip('\r\n')
            if linef1 != linef2:
                return False
        return next(file1, None) is None and next(file2, None) is None


def path_tmp(filename):
    "Return the path to test/tmp/filename using os.path.join()."
    return os.path.abspath(os.path.join(TEST_DIR, "tmp", filename))


def path_dat(filename):
    "Return the path to test/data/filename using os.path.join()."
    return os.path.abspath(os.path.join(TEST_DIR, "data", filename))


# os.path.join("test", "tmp")
PATH_TMP = os.path.abspath(os.path.join(TEST_DIR, "tmp"))
PATH_PHYLOSOFS = os.path.abspath(
    os.path.join(TEST_DIR, "..", "phylosofs", "phylosofs.py"))

if not os.path.isdir(PATH_TMP):
    os.mkdir(PATH_TMP)


class Test_PhyloSofS(unittest.TestCase):
    def test_phylosofs(self):
        path_transcripts = os.path.abspath(
            os.path.join(TEST_DIR, "..", "dat", "JNK3.transcripts"))
        path_newick = os.path.abspath(
            os.path.join(TEST_DIR, "..", "dat", "JNK3.nwk"))
        command = [
            # "python", PATH_PHYLOSOFS,
            "phylosofs",
            "-P",
            "-o",
            PATH_TMP,
            "--tree",
            path_newick,
            "--transcripts",
            path_transcripts
        ]
        self.assertEqual(subprocess.call(command), 0)
        self.assertTrue(
            compare_files(path_tmp('treeSearch_532_1.txt'),
                          path_dat('treeSearch_532_1.txt')))
        self.assertTrue(
            compare_files(path_tmp('solution_532_1_config0.sum'),
                          path_dat('solution_532_1_config0.sum')))

        # exons and transcripts changed order
        with open(path_tmp('solution_532_1_config0.info'),
                  'r') as out, open(path_tmp('solution_532_1_config0.info'),
                                    'r') as ref:
            for line_out, line_ref in zip(out, ref):
                if len(line_out) >= 5:
                    self.assertTrue(line_out[0:5] == line_ref[0:5])

        self.assertFalse(
            os.path.isdir(os.path.join(TEST_DIR, "tmp", "bestTopos")))
        self.assertFalse(
            os.path.isdir(os.path.join(TEST_DIR, "tmp", "betterTrees")))

    def test_best_topos_and_trees(self):
        path_transcripts = os.path.abspath(
            os.path.join(TEST_DIR, "..", "dat", "JNK3.transcripts"))
        path_newick = os.path.abspath(
            os.path.join(TEST_DIR, "..", "dat", "JNK3.nwk"))
        command = [
            # "python", PATH_PHYLOSOFS,
            "phylosofs",
            "-P",
            "-o",
            PATH_TMP,
            "-s",
            "100",
            "--tree",
            path_newick,
            "--transcripts",
            path_transcripts
        ]
        self.assertEqual(subprocess.call(command), 0)
        self.assertTrue(
            os.path.isdir(
                os.path.abspath(os.path.join(TEST_DIR, "tmp", "bestTopos"))))
        self.assertTrue(
            os.path.isdir(
                os.path.abspath(os.path.join(TEST_DIR, "tmp", "betterTrees"))))
        self.assertGreater(
            len(
                os.listdir(
                    os.path.abspath(os.path.join(TEST_DIR, "tmp",
                                                 "bestTopos")))), 0)
        self.assertGreater(
            len(
                os.listdir(
                    os.path.abspath(
                        os.path.join(TEST_DIR, "tmp", "betterTrees")))), 0)

    def tearDown(self):
        phylosofs.utils.clear_folder(PATH_TMP)


class Test_Internal_Functions(unittest.TestCase):
    def test_parse_pir(self):
        pir_path = os.path.abspath(
            os.path.join(TEST_DIR, "..", "dat", "transcripts.pir"))
        (seqs, lens, exons) = phylosofs.modelIsoforms.parse_pir(pir_path)
        print(seqs)
        print(lens)
        print(exons)
        self.assertEqual(len(seqs), 2)
        self.assertEqual(len(lens), 2)
        self.assertEqual(len(exons), 2)
        self.assertEqual(seqs[">transcript_one\n"],
                         "MSRHFLYNCSEPTLDVKIAFCQGFGNQVDVSYIAKHYNMS*\n")
        self.assertEqual(lens[">transcript_one\n"],
                         [1, 5, 9, 13, 17, 21, 25, 29, 33, 37, 40])
        self.assertEqual(exons[">transcript_one\n"],
                         "1111{{{{}}}}[[[[]]]]6666777788889999㐂㐂㐂㐂\n")
        self.assertEqual(seqs[">transcript_two\n"],
                         "MSRLYNCSELDVKIAQGFGNQVSYIAKNMS*\n")
        self.assertEqual(lens[">transcript_two\n"],
                         [1, 4, 7, 10, 13, 16, 19, 22, 25, 28, 30])
        self.assertEqual(exons[">transcript_two\n"],
                         "φφφ222333444555666777888999㐆㐆㐆\n")


if __name__ == '__main__':
    unittest.main()
