# README #

* Please look at the file "grammar.pdf" first to knwo what the grammar that are take into account here

### What is this repository for? ###

* for most of the sequence learning project you need to generate many sequences that respect REber grammar or some variations. This project provides files to generate them.
* here is the explanation of the parameters and in which case (grammar) you need them : 
*  - Type_grammar : indicate the grammar used for the sequence generation
* 		 0 -> Manual grammar (MG)
* 		 1 -> simple reber's grammar (RG)
* 		 2 -> Embedded and continious reber's grammar (CERG)
* 		 4 -> Embedded reber's grammar (ERG)

*  - lenCSeq: set the number of C in the manual grammar (the sequence starting with B always have an additionnal C)	
* 		this parameter is usefull only for manual grammar when type_grammar=0

*  - listNbSequences : set the number of sequences to generate _ always for Reber grammars and its variation (not manual grammar)

*  - listSizeSequences : set the size of the sequence or stream
*       this parameter is usefull only for CERG when type_grammar=2

*  - debug : to set to 1 in order to display the print 

*  - comment : for making comments into json file 


### How do I get set up? ###

* Just Download the python file and on parameters.json file, set your parameters and launch the command :
* python sequencegeneration.py 
* or
* python sequencegeneration.py parameters.json
* 
* NB : if no json file is specified, it will look by default for "parameters.json" file 

### Contribution guidelines ###



### Who do I talk to? ###

* Repo owner or admin
* Feel free to contact me if you have any question!