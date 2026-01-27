from timeit import default_timer as timer
import math
import pandas as pd
import numpy as np

class CombinationCreator:
    def __init__(self, p_NumOfevents, p_outcomes):
        """ Constructor"""""
        self.__events               = [0 for i in range(p_NumOfevents)]
        self.__outcomes             = [p_outcomes for i in range(p_NumOfevents)] if isinstance(p_outcomes, int) else p_outcomes
        self.__percentageOutput     = 0.1
        self.__plotCombinations     = False
        self.__eventsContainer      = []
        self.__eventsPermContainer  = []
        self.__saveEventsContainer  = False
        self.__isMultisetMode       = False
        self.__count                = 0
        self.__saveEventsFilename   = None

    def __del__(self):
        try:
            # Release potentially large memory structures
            self.__eventsContainer.clear()
            self.__events.clear()
            self.__outcomes.clear()
        except Exception:
            # Avoid raising exceptions from __del__
            pass

    def nextCombination(self, p_level):
        self.__events[p_level] = (self.__events[p_level] + 1) % self.__outcomes[p_level]
        if self.__events[p_level] == 0 and p_level + 1 < len(self.__events):
            self.nextCombination(p_level + 1)

    def nextCombinationWithRepetition(self, p_level):
        # increment current level (with wrap-around)
        self.__events[p_level] = (self.__events[p_level] + 1) % self.__outcomes[p_level]

        # if overflow or ordering constraint is violated, advance recursively
        if(self.__events[p_level] == 0 and p_level + 1 < len(self.__events)):
            self.nextCombinationWithRepetition(p_level + 1)

        if(p_level < len(self.__events)):
            if(p_level > 0 and (self.__events[p_level - 1] < self.__events[p_level]) == True):
                self.__events[p_level - 1] = self.__events[p_level]

    def createCombinations(self):
        self.__eventsContainer      = []
        self.__eventsPermContainer  = []
        _totalsCombinations     = math.prod(self.__outcomes)
        _totalsMultisets        = math.comb(len(self.__events) + self.__outcomes[0] - 1, len(self.__events))

        print("Total Combinations = ", _totalsCombinations)
        print("Total Multisets    = ", _totalsMultisets)

        self.__count = 0
        while True:
            #print(i + 1, "steps done of", _totals, "in %", 100 * float(i + 1) / float(_totals), "%", end="\r")

            if self.getIsMultisetMode() == True:
                _permu = self.computeMultisetsPermutations()

            if self.getPlotCombinations() == True:
                if self.__count == 0:
                    print(
                        f"{'No.':>{3}}  "
                        f"{'Events':<{5}}  "
                        f"{'Permutations':>{20}}"
                    )
                    print("-" * (3+5+20))
                print(f"{self.__count + 1:{4}d}.  {self.__events}", "\t", _permu)

            if self.getSaveEventsContainer() == True:
                self.__eventsContainer.append({"No.": self.__count + 1, "Events": tuple(self.__events), "Permutations": _permu})

            if self.getIsMultisetMode() == True:
                self.nextCombinationWithRepetition(0)
            else:
                self.nextCombination(0)

            self.__count += 1
            if (self.__count >= _totalsCombinations or sum(self.__events) == 0):
                break

        print("Total Iterations     = ", self.__count)
        if self.getIsMultisetMode():
            print("Summed Permutations   = ", sum(self.__eventsPermContainer))

        if self.getSaveEventsContainer() == True and self.__saveEventsFilename is not None:
            self.saveEventContainerToCsv()

    def getPlotCombinations(self):
        return self.__plotCombinations

    def setPlotCombinations(self, p_value):
        """" set to true to plot the combinations in the console"""
        self.__plotCombinations = p_value

    def getEvents(self):
        return self.__events

    def setEvents(self, p_events):
        self.__events = p_events

    def getIsMultisetMode(self):
        return self.__isMultisetMode

    def setMultisetMode(self, p_isMultisetMode):
        """
        enabling multiset-mode withdraws equal combinations assuming the order is irrelevant
        example: treating 1,2,3 as 3,2,1
        enabling multiset-mode require uniform outcomes (single outcome k for all events) if permutation reduction is desired
        """
        self.__isMultisetMode = p_isMultisetMode

        # enabling repetition-mode require uniform outcomes (single k)
        if p_isMultisetMode:
            if self.__outcomes and not all(x == self.__outcomes[0] for x in self.__outcomes):
                raise ValueError("withRepetition=True requires the same number of outcomes for all events")

    def getOutcomes(self):
        return self.__outcomes

    def setOutcomes(self, p_outcomes):
        self.__outcomes = p_outcomes

    def getSaveEventsContainer(self):
        return self.__saveEventsContainer

    def setSaveEventsContainer(self, p_saveEventsContainer, p_filename):
        """" set to true to save the combinations in a container and in the given filename"""
        self.__saveEventsContainer = p_saveEventsContainer
        self.__saveEventsFilename = p_filename

    def getEventsContainer(self):
        return self.__eventsContainer

    def saveEventContainerToCsv(self):
        df = pd.DataFrame(self.__eventsContainer)
        df.to_csv(self.__saveEventsFilename, index=False, sep=";", header=False)

    def computeMultisetsPermutations(self):
        if self.getIsMultisetMode():
            _eventsPerm         = 0
            _numOfEvents        = len(self.__events)
            _numOfOutcomes      = self.__outcomes[0]
            _possibleOutcomes   = list(range(_numOfOutcomes))
            _outcomesFactorials = [0] * _numOfOutcomes

            for i in range(_numOfOutcomes) :
                _outcomesFactorials[i] = math.factorial(self.__events.count(_possibleOutcomes[i]))

            _eventsPerm         = int(math.factorial(_numOfEvents) / (math.prod(_outcomesFactorials)))
            self.__eventsPermContainer.append(_eventsPerm)
            return _eventsPerm
        return ""

    def getRatioOfCombinationWithoutByWithRepitition(self, p_outcomes, p_events):
        return (math.comb(p_outcomes + p_events - 1, p_events) / (p_outcomes ** p_events))

    def getMatrixValuesBasedOnF(self, F, p_colSize, p_rowSize):
        matrix = \
            [
                [
                F(outc, event)
                for outc in range(1, p_colSize + 1)
                ]
                    for event in range(1, p_rowSize + 1)
            ]
        return (pd.DataFrame
            (
            matrix,
            index=[f"{i} Event" for i in range(1, p_colSize + 1)],
            columns=[f"{j} Outcome" for j in range(1, p_rowSize + 1)]
            ))