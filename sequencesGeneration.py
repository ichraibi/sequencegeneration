#!/usr/local/bin/python
# -*- coding: utf-8 -*
"""
This file is for sequences generation according the papers of Gers
Gers, F. A., Schmidhuber, J., & Cummins, F. (1999). Learning to forget: Continual prediction with LSTM.

"""

import numpy as np
import random

class grammar :
    def __init__(self, Letters):
        self.LettersList = Letters

        self.LettersDict = self.create_dictLettersWithBinaryValues(self.LettersList)

        #print "\n Letters list : "
        #for aLetter in self.LettersList:
        #    print aLetter, "-->", self.LettersDict[aLetter]

    def create_dictLettersWithBinaryValues(self, LetterList):

        LettersDict = {}
        #print "\n\nGeneration of the binary vectors for the letters"
        i = 0
        while (i < len(LetterList)):
            currentLetter = LetterList[i]
            ABinaryLetter = []
            for aletter in LetterList:

                if aletter == currentLetter:
                    ABinaryLetter.append(1)
                else:
                    ABinaryLetter.append(0)

            LettersDict[currentLetter] = ABinaryLetter

            i = i + 1


        for key, value in LettersDict.iteritems():
            foundLetter = ""
            foundLetter = self.getLetterByBinaryValue(LettersDict,value)

        return LettersDict

    #permet de connaitre de quel lettre il s'agit selon le vecteur binaire
    def getLetterByBinaryValue(self, LettersDict,ABinaryValue):

        theLetter =""
        for KeyLetter, BinaryLetter in LettersDict.iteritems():
            #print KeyLetter
            #print BinaryLetter
            if np.array_equal(BinaryLetter,ABinaryValue)==True :
                theLetter = KeyLetter
                break

        #print "theLetter : ",theLetter
        return theLetter

    def create_sequences_simpleGrammar(self, lenCSeq, debug =False):
        #Letters = ["A", "B", "C", "E", "F"]
        sequences=[]

        data = {}
        index_data = 0

        for lenC in lenCSeq:
            sequences = []
            nbSamples = 0
            numSeq = 0
            LenSequences = []
            # ACCCCCF
            # BCCCCCCCCCCCCCCCCCCE
            seq1 = "A"
            seq2 = "B"
            for j in range(0, 1):
                for i in range(0, lenC):
                    seq1 += "C"
                    seq2 += "C"

                seq2 += "C"
                seq1 += "F"
                seq2 += "E"

            # for i in range(0, 3):
            #    seq1 += "C"
            #    seq2 += "C"

            #if debug ==True :
            print "\n"
            print "Seq1 len of C : ", str(lenC), " - Seq2 len of C : ", str(lenC + 1)
            print "seq1 : ", seq1
            print "seq2 : ", seq2
            sequences.append(seq1)
            sequences.append(seq2)

            if debug:
                print(seq1)
                print(seq2)


            while numSeq < len(sequences):
                sequence = sequences[numSeq]
                LenSequences.append(len(sequence))
                numSeq += 1

            NumberOfSequences = len(sequences)
            print "Generation of ", NumberOfSequences, " sequences done !"

            print " --> The average len is : ", self.mean(LenSequences)
            print " --> The standard deviation is : ", self.pstdev(LenSequences)

        return sequences

    def create_oneSequence_ReberGrammar(self):
        grammar = {
            0: [('T', 1), ('P', 2)],
            1: [('S', 1), ('X', 3)],
            2: [('T', 2), ('V', 4)],
            3: [('X', 2), ('S', 5)],
            4: [('P', 3), ('V', 5)]
        }

        sequence = 'B'
        index = 0
        while index <> 5:
            choices = grammar[index]
            choice = np.random.randint(0, len(choices))
            token, index = choices[choice]
            sequence += token

        sequence += 'E'

        return sequence

    def generateRandomSequence(self,typegrammar, n=100, debug=False):

        if typeGrammar == 0:
            if debug == True: print "LSTM with manual grammar"
            LettersList = [ "B", "C", "E", "F"]
            startingLetter = "A"
        elif typeGrammar == 1:
            if debug == True: print "LSTM with Reber grammar"
            LettersList = [ "T", "S", "X", "V", "P", "E"]
            startingLetter = "B"

        elif typeGrammar == 2:
            if debug == True: print "LSTM with embedded and continious Reber grammar"
            LettersList = [ "T", "S", "X", "V", "P", "E"]
            startingLetter = "B"

        elif typeGrammar == 4:
            if debug == True: print "LSTM with embedded Reber grammar"
            LettersList = [ "T", "S", "X", "V", "P", "E"]
            startingLetter = "B"

        #LettersList = ["T", "S", "X", "V", "P", "E"]

        sequences = []
        for i in range(n):
            sequence = ''
            CurrentLetter = startingLetter
            sequence += CurrentLetter
            while (CurrentLetter not in ["E", "F"]):
                aNeighbor = random.choice(LettersList)
                # print "aNeighbor : ", aNeighbor
                sequence += aNeighbor
                CurrentLetter = aNeighbor

            if debug:
                print(sequence)

            sequences.append(sequence)

        return sequences

    def create_sequences_ReberGrammar(self,n=100, debug=False):
        #print "ici on mettra le code de génération de sequences a partir d'une grammaire complexe"
        #LettersList = ["T", "S", "X", "V", "P", "E"]
        grammar = {
            0: [('T', 1), ('P', 2)],
            1: [('S', 1), ('X', 3)],
            2: [('T', 2), ('V', 4)],
            3: [('X', 2), ('S', 5)],
            4: [('P', 3), ('V', 5)]
        }

        sequences = []
        LenSequences = []
        for i in range(n):
            sequence = self.create_oneSequence_ReberGrammar()

            if debug:
                print(sequence)

            sequences.append(sequence)
            LenSequences.append(len(sequence))

        print "Generation of ", len(sequences), " sequences done !"

        print " --> The average len is : ", self.mean(LenSequences)
        print " --> The standard deviation is : ", self.pstdev(LenSequences)

        return sequences


    def create_continious_embedded_sequences_ReberGrammar(self,n=100, lenSeq=100, debug=False):

        #LettersList = ["T", "S", "X", "V", "P", "E"]
        grammar = {
            0: [('T', 1), ('P', 2)],
            1: [('B', 3)],
            2: [('B', 4)],
            3: [('T', 5), ('P', 6)],
            4: [('T', 7), ('P', 8)],
            5: [('S', 5), ('X', 9)],
            6: [('T', 6), ('V', 10)],
            7: [('S', 7), ('X', 11)],
            8: [('T', 8), ('V', 12)],
            9: [('X', 6), ('S', 13)],
            10: [('P', 9), ('V', 13)],
            11: [('X', 8), ('S', 14)],
            12: [('P', 11), ('V', 14)],
            13: [('E', 15)],
            14: [('E', 16)],
            15: [('T', 17)],
            16: [('P', 17)],
            17: [('E', 18)],
            18: [('B', 0)]
        }

        sequences = []
        LenSequences = []
        for i in range(n):
            sequence = ''

            sequence = 'B'
            index = 0
            while len(sequence)<lenSeq :
                choices = grammar[index]
                choice = np.random.randint(0, len(choices))
                token, index = choices[choice]
                sequence += token
                #if len(sequence)==lenSeq:
                #    break

            if debug:
                print(sequence)

            sequences.append(sequence)
            LenSequences.append(len(sequence))

        print "Generation of ", len(sequences), " sequences done !"

        print " --> The average len is : ", self.mean(LenSequences)
        print " --> The standard deviation is : ", self.pstdev(LenSequences)

        return sequences

    def create_embedded_sequences_ReberGrammar(self,n=100, lenSeq=100, debug=False):

        #LettersList = ["T", "S", "X", "V", "P", "E"]
        grammar = {
            0: [('T', 1), ('P', 2)],
            1: [('B', 3)],
            2: [('B', 4)],
            3: [('T', 5), ('P', 6)],
            4: [('T', 7), ('P', 8)],
            5: [('S', 5), ('X', 9)],
            6: [('T', 6), ('V', 10)],
            7: [('S', 7), ('X', 11)],
            8: [('T', 8), ('V', 12)],
            9: [('X', 6), ('S', 13)],
            10: [('P', 9), ('V', 13)],
            11: [('X', 8), ('S', 14)],
            12: [('P', 11), ('V', 14)],
            13: [('E', 15)],
            14: [('E', 16)],
            15: [('T', 17)],
            16: [('P', 17)]
        }

        sequences = []
        LenSequences = []
        for i in range(n):
            sequence = ''

            sequence = 'B'
            index = 0
            while index <> 17:
                choices = grammar[index]
                choice = np.random.randint(0, len(choices))
                token, index = choices[choice]
                sequence += token
                #if len(sequence)==lenSeq:
                #    break

            sequence += 'E'


            if debug:
                print(sequence)

            sequences.append(sequence)
            LenSequences.append(len(sequence))

        print "Generation of ", len(sequences), " sequences done !"

        print " --> The average len is : ", self.mean(LenSequences)
        print " --> The standard deviation is : ", self.pstdev(LenSequences)

        return sequences

    def mean(self, data):
        """Return the sample arithmetic mean of data."""
        n = len(data)
        if n < 1:
            raise ValueError('mean requires at least one data point')
        return sum(data) / float(n)  # in Python 2 use sum(data)/float(n)

    def _ss(self, data):
        """Return sum of square deviations of sequence data."""
        c = self.mean(data)
        ss = sum((x - c) ** 2 for x in data)
        return ss

    def pstdev(self, data):
        """Calculates the population standard deviation."""
        n = len(data)
        if n < 2:
            raise ValueError('variance requires at least two data points')
        ss = self._ss(data)
        pvar = ss / n  # the population variance
        return pvar ** 0.5

    def get_lettersDict(self):
        return self.LettersDict

def save_sequences(typeGrammar, sequences, rand=False, lenC=0):
    myrep = "Sequences/"
    if not os.path.exists(myrep):
        os.mkdir(myrep)

    if rand == False:
        seqType = "_sequences.npy"
    elif rand == True:
        seqType = "_rand_sequences.npy"

    if typeGrammar == 0:
        # "LSTM with manual grammar"
        fname = "Sequences/MG_lenC_" + str(lenC) + seqType

    elif typeGrammar == 1:
        # "LSTM with Reber grammar"
        fname = "Sequences/RG" + seqType

    elif typeGrammar == 2:
        # "LSTM with embedded and continious Reber grammar"
        fname = "Sequences/CERG" + seqType

    elif typeGrammar == 3:
        # "LSTM with embedded and continious Reber grammar"
        fname = "Sequences/Id" + seqType

    elif typeGrammar == 4:
        # "LSTM with embedded Reber grammar"
        fname = "Sequences/ERG" + seqType
    else:
        # "LSTM"
        fname = "Sequences/" + seqType

    print "\nSequences saved in : ", fname
    np.save(fname, np.asarray(sequences))

def load_sequences(typeGrammar, rand=False, debug=False, lenC=0):
    sequences = []

    if rand == False:
        seqType = "_sequences.npy"
    elif rand == True:
        seqType = "_rand_sequences.npy"

    if typeGrammar == 0:
        # "LSTM with manual grammar"
        fname = "Sequences/MG_lenC_" + str(lenC) + seqType

    elif typeGrammar == 1:
        # "LSTM with Reber grammar"
        fname = "Sequences/RG" + seqType

    elif typeGrammar == 2:
        # "LSTM with embedded and continious Reber grammar"
        fname = "Sequences/CERG" + seqType

    elif typeGrammar == 3:
        # "LSTM with embedded and continious Reber grammar"
        fname = "Sequences/Id" + seqType

    elif typeGrammar == 4:
        # "LSTM with embedded Reber grammar"
        fname = "Sequences/ERG" + seqType
    else:
        # "LSTM"
        fname = "Sequences/" + seqType

    if os.path.isfile(fname):
        sequences = np.load(fname).tolist()

    print "\nLoading sequences from : ", fname
    if debug:
        print sequences

    return sequences


# -----------------------------------------------------------------------------------------
if __name__ == '__main__':
    import os, simplejson, sys

    # ----------------------------------------------lecture des paramètres -----------------
    if len(sys.argv) <= 1:
        parametersFile = "parameters.json"

    else:
        parametersFile = sys.argv[1]

    print " len(sys.argv) : ", len(sys.argv), " : ", sys.argv

    # ----------------------------------------------Parameters -----------------
    print "Chargement des parametres depuis le fichier ", parametersFile ," ... "
    paramtersDict = {}

    dir_path = os.path.dirname(os.path.abspath(__file__))
    parametersFilePath = dir_path + "/" + parametersFile
    # print "ElectricalDictFilePath : ", ElectricalDictFilePath
    if (os.path.exists(parametersFilePath) == True and os.path.isfile(parametersFilePath) == True):
        with open(parametersFilePath, 'r') as fp:
            paramtersDict = simplejson.load(fp)
    else:
        print "ERROR : The file of parameters called parameters.json doesn't exist :" + parametersFilePath
        print "Abort...."

    print " "
    # print "Parameters dictionnary : ", paramtersDict.keys()
    for key, value in paramtersDict.iteritems():
        print key, " : ", value

    print "\nParameters :"
    # SequenceSize = int(sys.argv[1]) #it need to be an even number (un chiffre pair)
    # listNumberOfSequences = sys.argv[1].split("-")


    lenCSeq = paramtersDict["lenCSeq"]
    listSizeSequences = paramtersDict["listSizeSequences"]
    Type_grammar = paramtersDict["Type_grammar"]
    listNbSequences = paramtersDict["listNbSequences"]
    debug = bool(paramtersDict["debug"][0])


    if debug == True :

        print "List of len of C : ", lenCSeq
        print "List of sizes of sequences : ", listSizeSequences
        print "List of nb sequences : ", listNbSequences


    if Type_grammar[0]==0:
        print "Type of grammars : Simple and manual "
        typeGrammar=0
        Letters = ["A", "B", "C", "E", "F"]
    elif Type_grammar[0]==1:
        print "Type of grammars : Reber's grammar "
        typeGrammar=1
        Letters = ["B","T", "S", "X", "V", "P", "E"]
    elif Type_grammar[0]==2:
        print "Type of grammars : continous and embedded Reber's grammar "
        typeGrammar=2
        Letters = ["B","T", "S", "X", "V", "P", "E"]
    elif Type_grammar[0]==3:
        print "Type of grammars : indentity "
        typeGrammar=3
        Letters = ["A","B", "C"]
    elif Type_grammar[0]==4:
        print "Type of grammars : Embedded Reber's grammar "
        typeGrammar=4
        Letters = ["B","T", "S", "X", "V", "P", "E"]

    print "List of Letters : ", Letters




    print "Chargement des parametres termine! "
    print "\n\n"
    parameters = []

    parameters.append(lenCSeq)
    parameters.append(typeGrammar)
    parameters.append(Letters)
    parameters.append(listNbSequences)
    parameters.append(listSizeSequences)

    print "Création de la structure du projet ... "

    listrepo = ['Sequences']
    for aRepo in listrepo:
        if not os.path.exists(aRepo):
            print "Creation de ", aRepo
            os.mkdir(aRepo)
        else:
            print aRepo, " existe déja"

    print "Génération de la structure du projet terminée! "
    print "\n\n"



    if typeGrammar==0 :
        # Generation of simple grammar : ACCCE and BCCCF
        manualGrammar = grammar(Letters)
        # sequences = simpleGrammar.create_sequences_simpleGrammar(lenCSeq)
        sequences = []
        LettersDict = manualGrammar.get_lettersDict()

        data = {}
        index_data = 0

        for lenC in lenCSeq:

            sequences = load_sequences(typeGrammar, False, debug)
            print "Nb sequences in the files loaded : ", len(sequences)
            if len(sequences) == 0:
                sequences = manualGrammar.create_sequences_simpleGrammar(lenCSeq, debug)
                save_sequences(typeGrammar, sequences, False, lenC)

            print "Nb sequences that we want :", len(sequences)
            if debug: print sequences

            rand_sequences = load_sequences(typeGrammar, True, debug)
            print "\nNb random sequences in the files loaded : ", len(rand_sequences)
            # print rand_sequences
            if len(rand_sequences) == 0:
                rand_sequences = manualGrammar.generateRandomSequence(typeGrammar, len(sequences), debug)
                rand_sequences.sort(key=len)
                save_sequences(typeGrammar, rand_sequences, True)

            print "Nb random sequences that we want :", len(rand_sequences)
            if debug: print rand_sequences

            print "\n"

    elif typeGrammar==1:
        # Generation of simple Reber grammar
        ReberGrammar = grammar(Letters)

        sequences = []
        LettersDict = ReberGrammar.get_lettersDict()

        data = {}
        index_data = 0

        for nbSeq in listNbSequences :
            sequences = load_sequences(typeGrammar,False,debug)
            print "Nb sequences in the files loaded : ", len(sequences)
            #quit()
            #print sequences
            if len(sequences) == 0:
                sequences = ReberGrammar.create_sequences_ReberGrammar(nbSeq, debug)
                sequences.sort(key=len)
                save_sequences(typeGrammar, sequences)
            elif len(sequences)<nbSeq :
                newSeq = ReberGrammar.create_sequences_ReberGrammar(nbSeq-len(sequences), debug)
                sequences +=newSeq
                sequences.sort(key=len)
                save_sequences(typeGrammar, sequences)
            elif len(sequences)>nbSeq:
                selectedSeq = []
                selectedSeq = sequences[:nbSeq]
                selectedSeq.sort(key=len)
                sequences = selectedSeq


            print "Nb sequences that we want : ", nbSeq, " - len(sequences) :", len(sequences)
            if debug: print sequences

            rand_sequences = load_sequences(typeGrammar, True, debug)
            print "\nNb random sequences in the files loaded : ", len(rand_sequences)
            #print rand_sequences
            if len(rand_sequences) == 0:
                rand_sequences = ReberGrammar.generateRandomSequence(typeGrammar,nbSeq, debug)
                rand_sequences.sort(key=len)
                save_sequences(typeGrammar, rand_sequences, True)

            elif len(rand_sequences) < nbSeq:
                newSeq = ReberGrammar.generateRandomSequence(typeGrammar,nbSeq - len(sequences), debug)
                rand_sequences += newSeq
                rand_sequences.sort(key=len)
                save_sequences(typeGrammar, rand_sequences,True)
            elif len(rand_sequences) > nbSeq:
                selectedSeq=[]
                for i in range(nbSeq):
                    selectedSeq.append(random.choice(rand_sequences))
                #selectedSeq.sort(key=len)
                rand_sequences = selectedSeq

            print "Nb random sequences that we want : ", nbSeq, " - len(rand_sequences) :", len(rand_sequences)
            if debug : print rand_sequences

            print "\n"


    elif typeGrammar==2:
        # Generation of sequences from an embedded and continious Reber grammar
        ECReberGrammar = grammar(Letters)

        sequences = []
        randSeq=[]
        LettersDict = ECReberGrammar.get_lettersDict()

        data = {}
        index_data = 0

        for nbSeq in listNbSequences :

            for lenSeq in listSizeSequences:
                sequences = load_sequences(typeGrammar,False, debug)
                print "Nb sequence in the files loaded : ", len(sequences)
                i =0
                for oneSeq in sequences:
                    print "\n",i , " - len : ", len(oneSeq)
                    print  oneSeq
                    i+=1
                    #if len(oneSeq) not in ['100000']:
                    #    print "len : ", len(oneSeq)
                    #print "ex : ", sequences[0]

                #print "fin de la vérification "
                #quit()

                if len(sequences) == 0:
                    sequences = ECReberGrammar.create_continious_embedded_sequences_ReberGrammar(nbSeq, lenSeq, debug)
                    save_sequences(typeGrammar, sequences)
                elif len(sequences) < nbSeq:
                    newSeq = ECReberGrammar.create_continious_embedded_sequences_ReberGrammar(nbSeq - len(sequences),lenSeq, debug)
                    sequences += newSeq
                    sequences.sort(key=len)
                    save_sequences(typeGrammar, sequences)
                elif len(sequences) > nbSeq:
                    selectedSeq = sequences[:nbSeq]
                    selectedSeq.sort(key=len)
                    sequences = selectedSeq

                print "Nb sequence that we want : ", nbSeq
                if debug: print sequences
                print "\n"

    elif typeGrammar == 4:
        # Generation of sequences from an embedded Reber grammar
        EReberGrammar = grammar(Letters)

        sequences = []
        randSeq=[]
        LettersDict = EReberGrammar.get_lettersDict()

        data = {}
        index_data = 0

        for nbSeq in listNbSequences:
            for lenSeq in listSizeSequences:

                sequences = load_sequences(typeGrammar,False, debug)
                print "Nb sequence in the files loaded : ", len(sequences)
                #quit()
                if len(sequences) == 0:
                    sequences = EReberGrammar.create_embedded_sequences_ReberGrammar(nbSeq, lenSeq, debug)
                    save_sequences(typeGrammar, sequences)
                elif len(sequences) < nbSeq:
                    newSeq = EReberGrammar.create_embedded_sequences_ReberGrammar(nbSeq - len(sequences),lenSeq, debug)
                    sequences += newSeq
                    sequences.sort(key=len)
                    save_sequences(typeGrammar, sequences)
                elif len(sequences) > nbSeq:
                    selectedSeq = sequences[:nbSeq]
                    selectedSeq.sort(key=len)
                    sequences = selectedSeq

                print "Nb sequence that we want : ", nbSeq
                if debug: print sequences
                print "\n"

                #RandomSequences----------------------------------------------------

                rand_sequences = load_sequences(typeGrammar, True, debug)
                print "\nNb random sequences in the files loaded : ", len(rand_sequences)
                # print rand_sequences
                if len(rand_sequences) == 0:
                    rand_sequences = EReberGrammar.generateRandomSequence(typeGrammar, nbSeq, debug)
                    rand_sequences.sort(key=len)
                    save_sequences(typeGrammar, rand_sequences, True)

                elif len(rand_sequences) < nbSeq:
                    newSeq = EReberGrammar.generateRandomSequence(typeGrammar, nbSeq - len(sequences), debug)
                    rand_sequences += newSeq
                    rand_sequences.sort(key=len)
                    save_sequences(typeGrammar, rand_sequences, True)
                elif len(rand_sequences) > nbSeq:
                    selectedSeq = []
                    for i in range(nbSeq):
                        selectedSeq.append(random.choice(rand_sequences))
                    # selectedSeq.sort(key=len)
                    rand_sequences = selectedSeq

                print "Nb random sequences that we want : ", nbSeq, " - len(rand_sequences) :", len(rand_sequences)
                if debug: print rand_sequences

                print "\n"








