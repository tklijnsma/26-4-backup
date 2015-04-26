#!/usr/bin/env python
"""
Thomas:
Reads VHBB Ntuple, outputs Cat1 events where 6 jets, 3 subjets en 3 gen quarks are
matched.

"""


########################################
# Imports
########################################

import ROOT
import os
import pickle
import time
import datetime
import TTH.TTHNtupleAnalyzer.AccessHelpers as AH
from cfg_J_SJ_Q import Make_config


########################################
# Functions
########################################

def Get_TLorentz( particle, input_tree, config):

    if particle in config['jet_branches']:
        Type = 'Jet'
        extra_vars = config['jet_extra_vars']
    elif particle in config['quark_branches']:
        Type = 'Quark'
        extra_vars = config['quark_extra_vars']
    elif particle in config['subjet_branches']:
        Type = 'Subjet'
        extra_vars = config['subjet_extra_vars']
    else:
        print 'Error in Get_TLorentz: given particle is not specified in'\
            'configuration'
        return 0

    # List of the branchnames for the extra vars (where to look up in original root)
    extra_vars_branchnames = [ var.format(particle=particle) for var in extra_vars ]

    # Remove {particle} from extra vars (most sensible attribute name)
    extra_vars = [ var.format(particle='') for var in extra_vars ]


    tl_output = []
    
    # Get the variables for the TLorentz vector
    Pt =  AH.getter(input_tree, particle + 'pt')
    Eta =  AH.getter(input_tree, particle + 'eta')
    Phi =  AH.getter(input_tree, particle + 'phi')
    Mass =  AH.getter(input_tree, particle + 'mass')



    # Get the values for the extra variables (beyond what is necessary for TLorentz)
    #   - Every entry of the dict will then contain a list!
    #   - Particle specifics will be stored without particle reference, e.g.:
    #     '{particle}pdgId' is stored as 'pdgId'
    extra_vars_vals_dict = {}
    for (var, var_branch) in zip(extra_vars, extra_vars_branchnames):
        extra_vars_vals_dict[var] = AH.getter(input_tree, var_branch)

    # Remove duplicates: 2 particles are often repeated in the WZQuark entries
    removedupl = 0
    if len(Pt)>2:
        if Pt[0]==Pt[ len(Pt)-2 ] and Pt[1]==Pt[ len(Pt)-1 ]:
            removedupl = 2


    vars_val_dict = {}

    for i in range( len(Pt) - removedupl ):
        
        # Construct the vars_val_dict, the dict with values for only this i

        # Regular vars
        vars_val_dict[particle+'pt'] = Pt[i]
        vars_val_dict[particle+'phi'] = Phi[i]
        vars_val_dict[particle+'eta'] = Eta[i]
        vars_val_dict[particle+'mass'] = Mass[i]

        # Extra vars
        for key in extra_vars_vals_dict:
            vars_val_dict[key] = extra_vars_vals_dict[key][i]

        if Evaluate_Cutoff( particle, Type, vars_val_dict, config ):

            y = ROOT.TLorentzVector()
            y.SetPtEtaPhiM( Pt[i] , Eta[i] , Phi[i] , Mass[i] )

            # Fill in the extra variables
            for var in extra_vars:
                setattr(y, var, vars_val_dict[var])

            tl_output.append( y )
            
    return tl_output


def Evaluate_Cutoff( particle, Type, vars_val_dict, config ):
    return True


########################################
# Main
########################################

def main():
    
    ROOT.gROOT.SetBatch(True)

    ########################################
    # Get the configuration file
    ########################################

    Make_config()

    config_filename = 'cfg_J_SJ_Q.dat'
    if not os.path.isfile(config_filename):
        print "Error: Can't find configuration file {0}.dat".format(config_filename)
        return 0

    print 'Importing configuration data from {0}'.format( config_filename )
    pickle_f = open( config_filename , 'rb' )
    config = pickle.load( pickle_f )
    pickle_f.close()


    ########################################
    # Copy often used variables from configuration file
    ########################################

    input_root_file_name = config['input_root_file_name']
    input_tree_name = config['input_tree_name']
    output_root_file_name = config['output_root_file_name']



    quark_branches = config['quark_branches']
    jet_branches = config['jet_branches']
    subjet_branches = config['subjet_branches']

    standard_vars = [ 'pt', 'eta', 'phi', 'mass', 'E' ]

    quark_extra_vars = config['quark_extra_vars']
    jet_extra_vars = config['jet_extra_vars']
    subjet_extra_vars = config['subjet_extra_vars']

    print quark_extra_vars
    print jet_extra_vars
    print subjet_extra_vars

    # Only for output purposes
    separate_vars = ['delR_Q-J', 'delR_J-SJ' ]


    ########################################
    # Setup I/O
    ########################################

    # Input tree
    input_root_file = ROOT.TFile(input_root_file_name)
    input_tree = input_root_file.Get(input_tree_name)

    # Output tree
    output_root_file = ROOT.TFile(output_root_file_name,'RECREATE')
    output_tree = ROOT.TTree('tree','My test tree')

    # Define branches in output tree
    branches = []

    branches.extend( [ 'Jet_' + var for var in standard_vars ] )
    branches.extend( [ 'Jet_' + var.format(particle='') for var in jet_extra_vars ] )

    branches.extend( [ 'Quark_' + var for var in standard_vars ] )
    branches.extend( [ 'Quark_' + var.format(particle='') \
        for var in quark_extra_vars ] )

    branches.extend( [ 'Subjet_' + var for var in standard_vars ] )
    branches.extend( [ 'Subjet_' + var.format(particle='') \
        for var in subjet_extra_vars ] )

    branches.extend( [ var for var in separate_vars ] )


    # Create dicitionaries to hold the information that will be
    # written as new branches
    variables      = {}
    variable_types = {}

    # Setup the output branches for the true object
    AH.addScalarBranches(variables,
                         variable_types,
                         output_tree,
                         branches,
                         datatype = 'float')


    ########################################
    # Event loop
    ########################################
    
    n_entries = input_tree.GetEntries()

    if config['Use_limited_entries']:
        n_processed = config['n_entries_limited']
    else:
        n_processed = n_entries
    print "Processing {0} events (out of {1} events)".format(n_processed, n_entries)

    for i_event in range(n_processed):

        if not i_event % int(0.001*n_processed+1):
            print "{0:.1f}% ({1} out of {2})".format(
                100.*i_event /n_processed, i_event, n_processed )

        input_tree.GetEntry( i_event )


        # Get quark and jet data - extra vars will be set as attributes

        tl_quarks = []
        for quark in quark_branches:
            tl_quarks.extend( Get_TLorentz( quark, input_tree, config ) )

        tl_jets = []
        for jet in jet_branches:
            tl_jets.extend( Get_TLorentz( jet, input_tree, config ) )

        tl_subjets = []
        for subjet in subjet_branches:
            tl_subjets.extend( Get_TLorentz( subjet, input_tree, config ) )

        
        for tl in tl_quarks:
            print tl


"""
        ########################################
        # delR combinatorics
        ########################################

        if Just_Jets == True:
            # Don't perform linking if 'Just_Jets' parameter is true
            # In that case, just write all found jet data
            for jet in tl_jets:

                variables['Jet_pt'][0] = jet.Pt()
                variables['Jet_eta'][0] = jet.Eta()
                variables['Jet_phi'][0] = jet.Phi()
                variables['Jet_mass'][0] = jet.M()

                variables['Jet_E'][0] = jet.E()
                variables['Jet_mcE'][0] = jet.mcE

                # Retrieve the extra variables from set attributes
                for var in extra_jet_vars:
                    variables['Jet_'+var][0] = getattr( jet, var )

                output_tree.Fill()
            continue
        
        # Otherwise, proceed with linking quarks and jets    
        [linked_jets, linked_quarks, delRs] = LinkJettoQuark(
            tl_jets, tl_quarks, config )

        # linked_jets and linked_quarks are (ordered) lists of TLorentzVectors. The
        # extra variables are contained in attributes.


        ########################################
        # Write to file
        ########################################

        for jet, quark, delR in zip( linked_jets, linked_quarks, delRs):

            variables['Jet_pt'][0] = jet.Pt()
            variables['Jet_eta'][0] = jet.Eta()
            variables['Jet_phi'][0] = jet.Phi()
            variables['Jet_mass'][0] = jet.M()

            variables['Jet_E'][0] = jet.E()

            if config['Get_MC_for_jets']:
                variables['Jet_mcE'][0] = jet.mcE

            variables['Quark_pt'][0] = quark.Pt()
            variables['Quark_eta'][0] = quark.Eta()
            variables['Quark_phi'][0] = quark.Phi()
            variables['Quark_mass'][0] = quark.M()

            variables['Quark_E'][0] = quark.E()

            variables['delR'][0] = delR

            # Retrieve the extra variables from set attributes
            for var in quark_extra_vars:
                var = var.format(particle='')
                variables[ 'Quark_'+ var ][0] = getattr( quark, var )
                
            for var in jet_extra_vars:
                var = var.format(particle='')
                variables[ 'Jet_'+var ][0] = getattr( jet, var )

            output_tree.Fill()

    output_root_file.Write()
    output_root_file.Close()

    ts = time.time()
    end_date = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')

    print config['info']
    print 'Analysis end time:              {0}'.format(end_date)
    print '.root file: {0}'.format( config['output_root_file_name'] )
"""

########################################
# End of main
########################################   


if __name__ == "__main__":
    main()
