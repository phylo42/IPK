import sys
from collections import defaultdict

EPS = 1e-3


def read_file(filename):
    with (open(filename, 'r') as f):
        lines = f.readlines()
        lines = lines
        data = defaultdict(dict)
        current_key = None

        for line in lines:
            line = line.strip()
            if "\t" not in line:  # Sequence name
                current_key = line.strip()
                data[current_key] = {}
            else:
                score, branch = line.strip().split("\t")
                data[current_key][int(branch)] = float(score)
        return data


class PlainTextDiff:
    def __init__(self, file1, file2, threshold):
        self.data1 = read_file(file1)
        self.data2 = read_file(file2)
        self.threshold = threshold

        assert self.data1 and len(self.data1), "File 1 is empty"
        assert self.data2 and len(self.data2), "File 2 is empty"

        self.differences = defaultdict(dict)

    def _ignore_diff(self, score1, score2):
        """Checks if score1 or score2 are close to the threshold. If any of them is None, check the other"""
        if score1 and score2:
            return abs(score1 - self.threshold) < EPS or abs(score2 - self.threshold) < EPS or abs(score1 - score2) < EPS
        elif score1:
            return abs(score1 - self.threshold) < EPS
        elif score2:
            return abs(score2 - self.threshold) < EPS
        else:
            return True

    def _report_diff(self, seq, branch, score1, score2):
        if not self._ignore_diff(score1, score2):
            self.differences[seq][branch] = (score1, score2)

    def compare_files(self):
        # Check sequences of File1 in File2
        for seq in self.data1:
            if seq in self.data2:
                for branch, score1 in self.data1[seq].items():

                    if seq == "AAAAAA" and branch == 117:
                        print("OKAY WAIT")

                    if branch in self.data2[seq]:
                        score2 = self.data2[seq][branch]
                        self._report_diff(seq, branch, score1, score2)
                    else:
                        if seq == "AAAAAA" and branch == 117:
                            print("WAIT")
                        self._report_diff(seq, branch, score1, None)
            else:
                for branch, score1 in self.data1[seq].items():
                    self._report_diff(seq, branch, score1, None)

        for seq in self.data2:
            if seq in self.data1:
                for branch, score2 in self.data2[seq].items():
                    if not (branch in self.data1[seq]):
                        self._report_diff(seq, branch, None, score2)

        return self.differences


if __name__ == "__main__":
    assert len(sys.argv) == 3, "Usage: diff-plain-text.py DB1 DB2"

    #eps = (7.5 / 20) ** 4
    eps = (1.5 / 4) ** 6
    ptd = PlainTextDiff(sys.argv[1], sys.argv[2], eps)
    differences = ptd.compare_files()
    if len(differences) == 0:
        print("OK")
    else:
        for kmer, diff in differences.items():
            print(kmer)
            for branch, scores in diff.items():
                print("\t", branch, scores[0], scores[1], sep="\t")

