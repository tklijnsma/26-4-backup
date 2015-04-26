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


    f1 = TFmat['l'][0].Make_Formula(False)
    print f1.GetTitle()

    f1.SetParameter( 0, 58.0 )

    formula = f1.GetTitle()


    gaussian_matches = re.finditer( r'exp\(' , formula )

    gaussian_strs = []

    for exp in gaussian_matches:

        pos = exp.end(0)

        bracket_counter = 0
        endpos = pos

        for (i,c) in enumerate(formula[pos:]):
            
            if c == '(':
                bracket_counter += 1

            if c == ')' and bracket_counter == 0:
                endpos = pos+i
                break
            elif c == ')':
                bracket_counter -= 1

        gaussian_strs.append( formula[pos:endpos] )


    
    for gaus in gaussian_strs:

        print '\nWorking on:\n' + gaus

        # Goal is to construct a string 'Erf(x)', where x is ( pt_minimum - b ) / (sqrt(2)*c)

        # Take the exponent found in a Gaussian; get rid of minus sign and take root
        #   This should result in a formula like abs(p*x+q), since root always takes the
        #   positive outcome of a root. By estimating a and b in this formula, the mean and
        #   rms of the found Gaussian can be evaluated.
        fg = ROOT.TF1( 'lin_gaus_arg', 'sqrt(-1*' + gaus + ')' , -10.0, 500.0)

        # Load the parameters from the original function into this linear function
        for i in range( f1.GetNumberFreeParameters() ):
            fg.SetParameter( i, f1.GetParameter(i) )
            print '[{0}] = {1}'.format( i, fg.GetParameter(i) )

        # Determine the mininum of abs(p*x+q); this should be the point at which a*x+b=0
        minx = fg.GetMinimumX()

        # Determine p
        p = ( fg.Eval(minx+20.0) - fg.Eval(minx+10.0) )/10.0

        # Determine sqrt(2)*c
        sq2c = 1.0/p
        c = sq2c/1.41421356237

        print 'p = {0} ==> sqrt(2)*c = {1}'.format(p,sq2c)
        print '            (c = {0})'.format(sq2c/1.41421356237)

        # Determine q
        q = fg.Eval(minx+10.0) - p * (minx+10.0)

        # Determine b
        b = -sq2c*q

        print 'q = {0} ==> b = {1}'.format(q,b)

        # Construct the error function that should replace this gaussian term
        #erfc = '{1}*TMath::Erfc((x-{0})/{1})'.format( b, sq2c )
        #TMath::PiOver2()*
        erfc = 'ROOT::Math::normal_cdf((x-{1}),{0})'.format( c, b )

        formula = formula.replace( 'exp(' + gaus + ')', erfc )


    print formula

    # Create the output TF1
    f_result = ROOT.TF1( 'loss_int', formula )

    # Load the parameters from the original function into this output
    #for i in range( f1.GetNumberFreeParameters() ):
    #    f_result.SetParameter( i, f1.GetParameter(i) )
    f_result.SetParameter( 1, 1.0 )

    pt_cut = 90.0

    print 'Integral from -inf to {0} = {1}'.format(pt_cut, f_result.Eval( pt_cut ) )

    f1.SetRange(-50.0, 170.0)
    ROOT.gROOT.SetBatch(True)
    ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = 1001;")
    ROOT.gStyle.SetOptFit(1011)
    c1 = ROOT.TCanvas("c1","c1",500,400)
    c1.SetGrid()
    f1.Draw()
    c1.Update()
    c1.Print("figure_erfreplace","pdf")
    


########################################
# End of Main
########################################
if __name__ == "__main__":
  main()
