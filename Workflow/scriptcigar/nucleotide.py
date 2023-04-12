class Nucleotide(): 

    def __init__(self, pos):
        """
        Parameter
        ---------
        pos : int
            Position in the reference sequence
        """
        self.count_exp = 0
        self.pos = pos
        self.score = []
    
    def increase_count(self):
        """
        Increase one of the coverage counter by a given number :
            --> If param == "expected" : increase of the expected coverage counter
            --> If param == "observed" : increase of the observed coverage counter
        
        Parameters
        ----------
        cpt : int
            Number that you want to add to the counter
        param : str
            Indicate which coverage counter to increase
                - "expected" for the expected coverage
                - "observed" for the observed coverage
        """

        self.count_exp += 1

    
    def decrease_score(self, cpt, ind):
        """
        Decrease the score chosen by a given number

        Parameters
        ----------
        cpt : int
            Number that you want to substract to the score
        ind : int
            Index indicating which score will be modified
        """
        self.score[ind] -= cpt

    def set_score(self, score):
        """
        Set a given score to the position

        Parameter
        ---------
        score : int
            Score to assign to the position
        """
        self.score.append(score)

    def get_count(self):
        """
        Getter to obtain the kind of coverage wanted :
          - expected coverage if param == "expected"
          - observed coverage if param == "observed"

        Parameter
        ---------
        param : str
            Indicate which coverage to obtain (expected coverage if "expected", observed coverage if "observed")

        Return
        ------
        int
            Expected coverage if param == "expected"
            Observed coverage if param == "observed"
        """

        return self.count_exp

    def get_pos(self):
        """
        Getter to obtain the position

        Return
        ------
        int
            Position in the reference sequence
        """
        return self.pos
    
    def get_score(self, disorder=False):
        """
        Getter to obtain the score

        Return
        ------
        list
            List of all the score associated to the position
        """
        if not disorder:
            return sorted(self.score, reverse=True)
        else:
            return self.score()
        
    def get_mean_score(self):
        if self.score:
            return sum(self.score)/len(self.score)
        return 0
