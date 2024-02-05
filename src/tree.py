"""
Simple suffix tree implementation on dicts
"""


class KmerTreeException(Exception):
    pass


class KmerTree:
    def __init__(
        self,
        k: int = 4,
        alphabet: str = "ACTG",
        sequence: str = "ACTGACTG",
        ignored_chars: str = "N",
    ):
        self.k = k
        self.alphabet = set(alphabet)
        self.ignored = set(ignored_chars)
        if set(alphabet).intersection(set(ignored_chars)):
            raise KmerTreeException("Alphabet and ignored char set cannot overlap")
        self.kmer_dict = dict()
        self.kmer_count = 0
        self.add_sequence(sequence)

    # TODO: skip N gracefully
    def add_sequence(self, sequence: str):
        """
        Add a sequence to the k-mer distribution.

        Extracts all k-mer frequencies from a sequence and adds them to self,
        updating k-mer count as appropriate. This method is meant for collecting
        data from multiple sequences, and is called at tree creation
        :return:
        """
        if not set(sequence).issubset(self.alphabet.union(self.ignored)):
            raise KmerTreeException(
                "Sequence alphabet is not a subset of Tree alphabet"
            )
        for i in range(len(sequence) - self.k + 1):
            kmer = sequence[i : i + self.k]
            if set(kmer).intersection(self.ignored):
                # Skip k-mers containing ignored characters
                # Note that these k-mers will still count towards the frequency
                # calculations.
                continue
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
        Return the number of times a given kmer was present in training seq(s).

        Note that 0 is returned for kmers with incorrect alphabet
        """
        if len(kmer) != self.k:
            raise KmerTreeException(
                "Requesting count of {len(kmer)}-mer from tree of {self.k}-mers"
            )
        if not set(kmer).issubset(self.alphabet):
            raise KmerTreeException(
                f"Requesting count of k-mer {kmer} that contains incorrect characters"
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
        if not set(sequence).issubset(self.alphabet):
            raise KmerTreeException(
                "Sequence alphabet is not a subset of Tree alphabet"
            )
        prob = 1.0
        pseudoprob = pseudocount / self.kmer_count
        for i in range(len(sequence) - self.k + 1):
            kmer_prob = self.get_freq(sequence[i : i + self.k])
            if kmer_prob > 0:
                prob *= kmer_prob
            else:
                prob *= pseudoprob
        return prob
