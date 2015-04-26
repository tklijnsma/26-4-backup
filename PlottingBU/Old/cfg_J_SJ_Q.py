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
        'This config.dat was created on: {0}'.format(config['date'])


    ########################################
    # I/O information
    ########################################

    #config['input_root_file_name'] = '/scratch/tklijnsm/VHBB_HEPPY_V11_G01_ttbar_13tev_phys14_20bx25.root'
    config['input_root_file_name'] = 'tth_test_tree.root'

    config['input_tree_name'] = 'tree'

    # The config file will be copied to 'runs/{config['outputdir']}'
    config['outputdir'] = 'V11-test1'

    config['output_root_file_name'] = 'V11-test1.root'


    ########################################
    # Program parameters
    ########################################

    # Use only a part of the input root file
    config['Use_limited_entries'] = False

    # Specify the number of entries if only a limited number of entries is used
    #   This number is not used if Use_limited_entries is set to False
    config['n_entries_limited'] = 10


    ########################################
    # Branch info
    ########################################

    # Specify the names of the particles of with pt, eta, phi and m should be
    # extracted. 
    #   - pt, eta, phi and mass are extracted by default. A branch E is created by
    #     default - it is calculated with the use of pt, eta, phi and mass.
    #   - Since the notation '_pt' is common but not *standard*, it is necessary to
    #     to add underscores where necessary manually.

    config['quark_branches'] = [
        'GenBQuarkFromTop_',
        'GenBQuarkFromH_',
        'GenWZQuark_' ]

    config['jet_branches'] = [
        'Jet_' ]

    config['subjet_branches'] = [
        'httCandidates_sjW1',
        'httCandidates_sjW2',
        'httCandidates_sjNonW']
    

    # Specify which branches *other* than pt, eta, phi, mass and E should be 
    # extracted.
    #   - This should be FULL branch names, e.g. httCandidates_fW
    #   - If the extra variable is particle-specific, write '{particle}' in front
    #     it. For example: '(a quark)pdgId' can be written as '{particle}pdgId'

    config['quark_extra_vars'] = [
        '{particle}pdgId' ]

    config['jet_extra_vars'] = []


    config['subjet_extra_vars'] = [
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

    config['quark_cutoff_list'] = [
        ( '{particle}pt'       , '>' , 30.0 ) ]


    config['max_link_delR'] = 10
    
    # Only used if config['Remove_double_match'] is set to True
    config['max_sec_delR'] = 0.5


    ########################################
    # Write configuration to file:
    ########################################

    config_filename = 'cfg_J_SJ_Q.dat'

    if os.path.isfile( config_filename ):
        os.remove( config_filename )

    f = open( config_filename, 'w' )
    pickle.dump( config , f )
    f.close()

    if not os.path.isdir('runs'):
        os.makedirs('runs')
    if not os.path.isdir('runs/{0}'.format(config['outputdir'] ) ):
        os.makedirs('runs/{0}'.format(config['outputdir'] ) )

    shutil.copyfile( config_filename,
        'runs/{0}/{1}'.format( config['outputdir'], config_filename ) )

    print "{0} created".format( config_filename )



########################################
# End of Main
########################################
def main():
    Make_config()

if __name__ == "__main__":
  main()
