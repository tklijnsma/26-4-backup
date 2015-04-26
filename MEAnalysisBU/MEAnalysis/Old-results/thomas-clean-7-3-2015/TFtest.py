#!/usr/bin/env python
"""
Dictionay and pickle test

"""

########################################
# Imports
########################################

import pickle
import ROOT
from config import Make_config


########################################
# Main
########################################

def main():

    ROOT.gROOT.SetBatch(True)
    ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = 1001;")
    ROOT.gStyle.SetOptFit(1011)

    # Create the config.dat file
    Make_config()

    # Open config.dat to config 
    pickle_f = open( 'config.dat', 'rb' )
    config = pickle.load( pickle_f )
    pickle_f.close()

    outputdir = config['outputdir']

    # Open TFMatrix.dat to TFmat
    pickle_f = open( '{0}/TFMatrix.dat'.format(outputdir), 'rb' )
    TFmat = pickle.load( pickle_f )
    pickle_f.close()

    # Choose a TF; As an example the TF for a b in the first eta-bin is chosen.
    myTF = TFmat['b'][0]


    ########################################
    # Plot an example of a response function
    ########################################

    # Determine Make_Formula output:
    # - True:  [0] = reconstr., x = mc/gen/quark
    # - False: [0] = mc/gen/quark, x = reconstr.
    # No argument is treated as True
    set_reconstructed_eval_gen = False

    f1 = myTF.Make_Formula( set_reconstructed_eval_gen )

    # Set gen pt to 58 GeV
    f1.SetParameter( 0, 58.0 )

    # Drawing
    f1.SetRange( 0.0 , 100.0)
    c1 = ROOT.TCanvas("c1","c1",500,400)
    c1.SetGrid()
    f1.Draw()
    c1.Update()
    c1.Print("exampleplot_TF","pdf")


    ########################################
    # Plot an example of a response function
    ########################################

    # Create the cumulative distribution function
    #  - Only works for single or double Gaussian functions
    #  - [0] = mc/gen/quark pt, x = pt_cutoff
    f2 = myTF.Make_CDF()

    # Set the mc/gen/quark pt
    f2.SetParameter( 0, 58.0 )

    # Drawing
    f2.SetRange( 0.0 , 100.0)
    f2.Draw()
    c1.Update()
    c1.Print("exampleplot_CDF","pdf")


########################################
# End of Main
########################################
if __name__ == "__main__":
  main()
