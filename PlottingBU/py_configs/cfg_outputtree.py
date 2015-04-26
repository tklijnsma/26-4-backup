#!/usr/bin/env python
"""
Thomas:

This program creates cfg_outputtree.dat, which contains a dictionary of all relevant input
and parameters for outputtree.py.

"""

########################################
# Imports
########################################

import pickle
import os
import shutil

import time
import datetime

########################################
# Main
########################################

def Make_config():

    config = {}

    ########################################
    # Information concerning this config file
    ########################################

    ts = time.time()
    config['date'] = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')

    config['info'] = '*** Information on this config.dat ***\n'\
    'This config file contains the configuration data used for outputtree.py, which translates a VHBB Ntuple to a format readable by readtree.py. \n\n'\
    'This config.dat was created on: {0}'.format(config['date'])


    ########################################
    # I/O information
    ########################################

    config['input_root_file_name'] = '/scratch/tklijnsm/VHBB_HEPPY_V11_G01_ttbar_13tev_phys14_20bx25.root'
    #config['input_root_file_name'] = 'tth_test_tree.root'

    config['input_tree_name'] = 'tree'

    # The config file will be copied to 'runs/{config['run_name']}', and will determine
    # the name of the .root file.
    config['run_name'] = 'V11_full_subjets_0.3delR'

    config['output_root_file_name'] = '/scratch/tklijnsm/{0}.root'.format( config['run_name'] )


    ########################################
    # Program parameters
    ########################################

    # Use only a part of the input root file
    config['Use_limited_entries'] = False

    # Specify the number of entries if only a limited number of entries is used
    #   This number is not used if Use_limited_entries is set to False
    config['n_entries_limited'] = 10

    # Specify whether the program should attempt to find MC branches for the jets
    config['Get_MC_for_jets'] = False

    # Specify whether the program should link with quarks.
    #   If this is set to True, the program will not save any quark data, and the
    #   the user should calculate TFs with the MC values of the jets. Make sure
    #   config['Get_MC_for_jets'] is set to True if this is set to True.
    config['Dont_Link_Just_Jets'] = False


    ########################################
    # Branch info
    ########################################

    # Specify the names of the particles of with pt, eta, phi and m should be
    # extracted. 
    #   - pt, eta, phi and mass are extracted by default. A branch E is created by
    #     default - it is calculated with the use of pt, eta, phi and mass.
    #   - Since the notation '_pt' is common but not *standard*, it is necessary to
    #     to add underscores where necessary manually.

    config['quarktypes'] = ['GenBQuarkFromTop_', 'GenBQuarkFromH_', 'GenWZQuark_' ]
    config['jettypes'] = [ 'httCandidates_sjW1', 'httCandidates_sjW2',
        'httCandidates_sjNonW']

    # Specify which branches *other* than pt, eta, phi, mass and E should be 
    # extracted.
    #   - This should be FULL branch names, e.g. httCandidates_fW
    #   - If the extra variable is particle-specific, write '{particle}' in front
    #     it. For example: '(a quark)pdgId' can be written as '{particle}pdgId'

    config['quark_extra_vars'] = [
        '{particle}pdgId',
        '{particle}charge',
        '{particle}status' ]

    config['jet_extra_vars'] = [
        'httCandidates_pt',
        'httCandidates_mass',
        'httCandidates_fW' ]



    ########################################
    # Cutoff criteria
    ########################################

    # format of 1 cutoff criterium: ( varname, operator sign, cutoff value )
    #   Note: Only defined variable names can be used here!

    config['jet_cutoff_list'] = [
        ( '{particle}pt'       , '>' , 30.0 ),
        ( 'httCandidates_pt'   , '>' , 200.0 ),
        ( 'httCandidates_mass' , '>' , 120.0 ),
        ( 'httCandidates_mass' , '<' , 220.0 ),
        ( 'httCandidates_fW'   , '<' , 0.175 ) ]

    #config['jet_cutoff_list'] = [
    #    ( '{particle}pt'       , '>' , 30.0 ) ]

    config['quark_cutoff_list'] = [
        ( '{particle}pt'       , '>' , 30.0 ) ]


    config['max_link_delR'] = 0.3
    config['max_sec_delR'] = 0.5


    ########################################
    # Write configuration to file:
    ########################################

    f = open( 'cfg_outputtree.dat', 'wb' )
    pickle.dump( config , f )
    f.close()

    if not os.path.isdir('runs'):
        os.makedirs('runs')
    if not os.path.isdir('runs/{0}'.format(config['run_name'] ) ):
        os.makedirs('runs/{0}'.format(config['run_name'] ) )

    shutil.copyfile( 'cfg_outputtree.dat',
        'runs/{0}/cfg_outputtree.dat'.format( config['run_name'] ) )

    print "cfg_outputtree.dat created"



########################################
# End of Main
########################################
def main():
    Make_config()

if __name__ == "__main__":
  main()
