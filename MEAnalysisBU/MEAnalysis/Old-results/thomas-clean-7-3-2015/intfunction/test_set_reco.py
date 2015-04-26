#!/usr/bin/env python
"""
Dictionay and pickle test

"""

########################################
# Imports
########################################

import pickle
import re
import ROOT

########################################
# Main
########################################

def main():

    pickle_f = open( 'TFMatrix.dat', 'rb' )
    TFmat = pickle.load( pickle_f )
    pickle_f.close()


    f2 = TFmat['b'][0].Make_Formula(False)
    print f2.GetTitle()


    f2.Make_CDF()


    


    



########################################
# End of Main
########################################
if __name__ == "__main__":
  main()
