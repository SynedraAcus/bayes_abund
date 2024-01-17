"""
Simple suffix tree implementation on dicts
"""


class KmerTreeException(Exception):
    pass


class KmerTree:
    def __init__(self, k=4, alphabet="ACTG", sequence="ACTGACTG"):
        self.k = k
        self.alphabet = alphabet
        # TODO: check that sequence is from correct alphabet here and elsewhere
        self._populate_dict(sequence)

    # TODO: a method for adding sequences to existing KmerTree instances
    # Intended to use with higher taxa like eg genera

    def _populate_dict(self, sequence):
        self.kmer_dict = dict()
        for i in range(len(sequence) - self.k + 1):
            kmer = sequence[i : i + self.k]
            focus = self.kmer_dict
            for letter in kmer[:-1]:
                if letter not in focus:
                    focus[letter] = dict()
                focus = focus[letter]
            if kmer[-1] in focus:
                focus[kmer[-1]] += 1
            else:
                focus[kmer[-1]] = 1
        self.kmer_count = len(sequence) - self.k + 1

    def get_count(self, kmer: str) -> int:
        """
        Return the number of times a given kmer was present in training seq
        """
        if len(kmer) != self.k:
            raise KmerTreeException(
                "Requesting count of {len(kmer)}-mer from tree of {self.k}-mers"
            )
        focus = self.kmer_dict
        for letter in kmer[:-1]:
            if letter not in focus:
                focus[letter] = dict()
            focus = focus[letter]
        try:
            return focus[kmer[-1]]
        except KeyError:
            return 0

    def get_freq(self, kmer: str) -> float:
        """
        Return the frequency of a given kmer (as a float in [0.0, 1.0] range)
        :param kmer:
        :return:
        """
        return self.get_count(kmer) / self.kmer_count

    def get_seq_prob(self, sequence: str, pseudocount: float = 0.1) -> float:
        """
        Calculate the likelihood of a query sequence under this kmer
        distribution with Naive Bayesian assumptions (ie as a product of
        probabilities of all its kmers)

        All k-mers absent in original sequence are assigned the pseudocount;
        this should be in (0.0, 1.0) range, non-inclusively. The lesser the
        value, the less probability is assigned.
        """
        prob = 1.0
        pseudoprob = pseudocount / self.kmer_count
        for i in range(len(sequence) - self.k + 1):
            kmer_prob = self.get_freq(sequence[i : i + self.k])
            if kmer_prob > 0:
                prob *= kmer_prob
            else:
                prob *= pseudoprob
        return prob
