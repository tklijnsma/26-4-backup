#!/usr/bin/env python
"""
Regular expression test

"""

########################################
# Imports
########################################

import re
import ROOT

########################################
# Main
########################################

def main():

    formula = "[0]*({0}*exp(-0.5*((x-[1])/[2])**2)+{1}*exp(-0.5*((x-[3])/([2]+[4]))**2))"

    result = re.finditer( r'exp\(' , formula )

    gaussian_strs = []
    erf_strs = []

    for exp in result:

        pos = exp.end(0)

        bracket_counter = 0
        endpos = pos

        for (i,c) in enumerate(formula[pos:-1]):
            
            if c == '(':
                bracket_counter += 1

            if c == ')' and bracket_counter == 0:
                endpos = pos+i
                break
            elif c == ')':
                bracket_counter -= 1

        gaussian_strs.append( 'exp(' + formula[pos:endpos] + ')' )
        erf_strs.append( 'TMath::PiOver2()*TMath::Erf(' + formula[pos:endpos] + ')' )


    for (gaus, erf) in zip( gaussian_strs, erf_strs ):

        formula = formula.replace( gaus, erf )

    print formula

    print 'End of formula making\n'


    f1 = ROOT.TF1( 'test1', formula.format(0.5, 0.5 ) )


    #f2 = ROOT.TF1( 'exppart', 'sqrt(-1*'+gaussian_strs[0][4:-1]+')', 1, 3 )
    f2 = ROOT.TF1( 'exppart', gaussian_strs[0][4:-1], 1, 3 )

    f2.SetParameter( 1, 2 )
    f2.SetParameter( 2, 3 )
    

    print 'f2 function: ' + f2.GetTitle()


    f3 = ROOT.TF1( 'wortel', 'sqrt(-1*' + gaussian_strs[0][4:-1] + ')', -1, 3 )

    f3.SetParameter( 1, 2 )
    f3.SetParameter( 2, 3 )


    print 'f3 minimum: {0}'.format(f3.GetMinimumX())

    minx = f3.GetMinimumX()

    p = f3.Eval(minx+2.0) - f3.Eval(minx+1.0)

    sq2c = 1.0/p

    print 'p = {0} ==> sqrt(2)*c = {1}'.format(p,sq2c)
    print '            (c = {0})'.format(sq2c/1.41421356237)

    q = f3.Eval(minx+1.0) - p * (minx+1.0)

    b = -sq2c*q

    print 'q = {0} ==> b = {1}'.format(q,b)



    ROOT.gROOT.SetBatch(True)
    ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = 1001;")
    ROOT.gStyle.SetOptFit(1011)
    c1 = ROOT.TCanvas("c1","c1",500,400)
    c1.SetGrid()
    f3.Draw()
    c1.Update()
    c1.Print("figure_erfreplace","pdf")

########################################
# End of Main
########################################
if __name__ == "__main__":
  main()
