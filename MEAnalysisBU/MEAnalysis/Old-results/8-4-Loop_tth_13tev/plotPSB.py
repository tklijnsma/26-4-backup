#!/usr/bin/env python
"""
Thomas:

Creates a plot of P_S/B = w(y|tth) / ( w(y|tth) + k_S/B * w(y|ttbb) )

"""


########################################
# Imports
########################################

import os
import ROOT
import TTH.TTHNtupleAnalyzer.AccessHelpers as AH
import pickle
import copy


########################################
# Functions
########################################

def PlotSB():

    ROOT.gROOT.SetBatch(True)
    ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = 1001;")
    ROOT.gStyle.SetOptFit(1011)

    input_root_file_name = 'tree.root'
    input_tree_name = 'tree'

    input_root_file = ROOT.TFile(input_root_file_name)
    input_tree = input_root_file.Get(input_tree_name)

    k_sb = 0.12

    cat_list = [ 1, 2, 3, 6 ]



    ########################################
    # Drawing
    ########################################

    # Open an html-file
    hf = open( 'P_SB-cat.html', 'w' )
    hf.write( '<html><body>\n<h1>P_SB per category\n</h1>\n<br>\n<hr />' )
    #hf.write( '<h2>Title</h2>' )

    c1 = ROOT.TCanvas("c1","c1",500,400)
   
    for cat in cat_list:

        hist_name = 'hist-cat{0}'.format(cat)

        draw_str = 'mem_tth_p/(mem_tth_p+{0}*mem_ttbb_p)>>{1}'.format(
            k_sb, hist_name )
        sel_str = 'cat=={0}'.format( cat )

        # Retrieve the histogram
        input_tree.Draw(draw_str, sel_str)
        hist = getattr(ROOT, hist_name).Clone()

        hist.Draw()


        plot_outputdir = 'P_SB_plots'
        if not os.path.isdir(plot_outputdir):
            os.makedirs(plot_outputdir)

        plot_filename = '{0}/P_SB-cat{1}-hist'.format(plot_outputdir, cat)

        c1.Print( plot_filename, 'pdf' )

        img = ROOT.TImage.Create()
        img.FromPad(c1)
        img.WriteImage( '{0}.png'.format(plot_filename) )

        # Write line to html
        hf.write('<a href="{0}"><img width="500" src="{0}.png"></a>\n'.format(plot_filename) )



########################################
# End
########################################   

def main():
    PlotSB()

if __name__ == "__main__":
    main()
