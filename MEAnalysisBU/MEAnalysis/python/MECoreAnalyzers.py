import ROOT
import itertools
from PhysicsTools.HeppyCore.framework.analyzer import Analyzer
import copy

import sys

#Load the MEM integrator libraries
# ROOT.gSystem.Load("libFWCoreFWLite")
# ROOT.gROOT.ProcessLine('AutoLibraryLoader::enable();')
# ROOT.gSystem.Load("libFWCoreFWLite")
ROOT.gSystem.Load("libCintex")
ROOT.gROOT.ProcessLine('ROOT::Cintex::Cintex::Enable();')
ROOT.gSystem.Load("libTTHMEIntegratorStandalone")

from ROOT import MEM
o = MEM.MEMOutput

#Pre-define shorthands for permutation and integration variable vectors
CvectorPermutations = getattr(ROOT, "std::vector<MEM::Permutations::Permutations>")
CvectorPSVar = getattr(ROOT, "std::vector<MEM::PSVar::PSVar>")

def lvec(self):
    """
    Converts an object with pt, eta, phi, mass to a TLorentzVector
    """
    lv = ROOT.TLorentzVector()
    lv.SetPtEtaPhiM(self.pt, self.eta, self.phi, self.mass)
    return lv


class FilterAnalyzer(Analyzer):
    """
    A generic analyzer that may filter events.
    Counts events the number of processed and passing events.
    """
    def beginLoop(self, setup):
        super(FilterAnalyzer, self).beginLoop(setup)
        self.counters.addCounter("processing")
        self.counters["processing"].register("processed")
        self.counters["processing"].register("passes")

class EventIDFilterAnalyzer(FilterAnalyzer):
    """
    """

    def __init__(self, cfg_ana, cfg_comp, looperName):
        super(EventIDFilterAnalyzer, self).__init__(cfg_ana, cfg_comp, looperName)
        self.conf = cfg_ana._conf
        self.event_whitelist = self.conf.general.get("eventWhitelist", None)

    def beginLoop(self, setup):
        super(EventIDFilterAnalyzer, self).beginLoop(setup)

    def process(self, event):
        self.counters["processing"].inc("processed")

        passes = True
        if not self.event_whitelist is None:
            passes = False
            if (event.input.run, event.input.lumi, event.input.evt) in self.event_whitelist:
                print "IDFilter", (event.input.run, event.input.lumi, event.input.evt)
                passes = True

        if passes and "eventboundary" in self.conf.general["verbosity"]:
            print "---", event.input.run, event.input.lumi, event.input.evt
        if passes:
            self.counters["processing"].inc("passes")
        return passes


class TestAnalyzer(FilterAnalyzer):
    """
    Test analyzer by Thomas
    """

    def __init__(self, cfg_ana, cfg_comp, looperName):
        super(TestAnalyzer, self).__init__(cfg_ana, cfg_comp, looperName)
        self.conf = cfg_ana._conf

        self.R_cut = 0.3
        self.top_mass = 172.04

        self.Cut_criteria = [
            ( 'pt'  , '>', '200.0' ),
            ( 'mass', '>', '120.0' ),
            ( 'mass', '<', '220.0' ),
            ( 'fW'  , '<', '0.175' ) ]

        self.Statistics = {
            'n_processed'       : 0,
            'n_not_cat1'        : 0,
            'n_no_httCand'      : 0,
            'n_survivedcut'     : [0,0,0],
            'n_too_few_WZ'      : 0,
            'n_too_few_B'       : 0,
            'n_too_many_WZ'     : 0,
            'n_too_many_B'      : 0,
            'n_2q_for_1j'       : 0,
            'n_no_jet_for_q'    : 0,
            'n_no_unique_match' : 0,
            'n_successful'      : 0
            }

    def beginLoop(self, setup):
        super(TestAnalyzer, self).beginLoop(setup)

    def endLoop(self, setup):

        print 'Statistics:'
        print 'n_processed       = {0}'.format(self.Statistics['n_processed'])
        print 'n_not_cat1        = {0}'.format(self.Statistics['n_not_cat1'])
        print 'n_no_httCand      = {0}'.format(self.Statistics['n_no_httCand'])
        print 'n_survivedcut     = {0}'.format(self.Statistics['n_survivedcut'])
        print 'n_too_few_WZ      = {0}'.format(self.Statistics['n_too_few_WZ'])
        print 'n_too_few_B       = {0}'.format(self.Statistics['n_too_few_B'])
        print 'n_too_many_WZ     = {0}'.format(self.Statistics['n_too_many_WZ'])
        print 'n_too_many_B      = {0}'.format(self.Statistics['n_too_many_B'])
        print 'n_2q_for_1j       = {0}'.format(self.Statistics['n_2q_for_1j'])
        print 'n_no_jet_for_q    = {0}'.format(self.Statistics['n_no_jet_for_q'])
        print 'n_no_unique_match = {0}'.format(self.Statistics['n_no_unique_match'])
        print 'n_successful      = {0}'.format(self.Statistics['n_successful'])


    def process(self, event):

        self.Statistics['n_processed'] += 1
        #print 'Printing from TestAnalyzer! iEv = {0}'.format(event.iEv)


        ########################################
        # Check event suitability
        ########################################

        # Check if event is category 1
        if event.cat != 'cat1':
            self.Statistics['n_not_cat1'] += 1
            return 0

        # Check if there is a httCandidate
        if len( event.httCandidate ) == 0:
            self.Statistics['n_no_httCand'] += 1
            return 0

        # Apply the cuts
        tops = []
        for candidate in event.httCandidate:
            if self.Apply_Cut_criteria( candidate ):
                tops.append( candidate )

        # Check if any candidates survived the cutoff criteria
        if len(tops) == 0:
            print 'No candidate survived the cuts'
            self.Statistics['n_survivedcut'][0] += 1
            return 0

        elif len(tops) == 1:
            top = tops[0]
            self.Statistics['n_survivedcut'][1] += 1

        # If more than 1 candidate survived the cutoff criteria, choose the
        # one with a mass closest to top mass
        else:
            tops = sorted( tops, key=lambda x: abs(x.mass - self.top_mass) )
            top = tops[0]
            self.Statistics['n_survivedcut'][2] += 1

        # Necessary to remove duplicates in GenWZQuark branch
        self.CompareWZQuarks( event )

        # Get a list of the 3 generated quarks that should correspond to a
        # top candidate
        tl_genquarks = self.Get_tl_genquarks( event )

        # If there is an error in getting the quarks, the function returns 0
        if tl_genquarks == 0: return 0


        ########################################
        # Perform combinatorics and calculate delR
        ########################################

        # Determine delR for quarks with jets

        tl_jets = []

        for (i_jet, jet) in enumerate(event.good_jets):
            x = ROOT.TLorentzVector()
            x.SetPtEtaPhiM( jet.pt, jet.eta, jet.phi, jet.mass )
            setattr( x, 'i_jet', i_jet )
            tl_jets.append( x )


        ( jet_links , jet_delR_list ) = self.Do_delR_combinatorics(
            tl_genquarks, tl_jets )

        if jet_links == 0: return 0

        jet_sumdelR = sum( jet_delR_list )


        # Determine delR for quarks with subjets

        tl_subjets = []

        prefixes = [ 'sjW1', 'sjW2', 'sjNonW' ]

        for (i_subjet, prefix) in enumerate( prefixes ):

            x = ROOT.TLorentzVector()

            x.SetPtEtaPhiM(
                getattr( top, prefix + 'pt' ),
                getattr( top, prefix + 'eta' ),
                getattr( top, prefix + 'phi' ),
                getattr( top, prefix + 'mass' ) )

            setattr( x, 'i_subjet', i_subjet )

            tl_subjets.append( x )

        ( subjet_links , subjet_delR_list ) = self.Do_delR_combinatorics(
            tl_genquarks, tl_subjets )

        if subjet_links == 0:
            print 'Matching quarks with subjets was not successful'
            return 0

        subjet_sumdelR = sum( subjet_delR_list )


        ########################################
        # Save the quark data from successful matches
        ########################################

        print 'Successful match. Quarks were matched to these subjets:'
        print subjet_links

        self.Statistics['n_successful'] += 1

        # Write jet_delR values to event
        setattr( event.GenBQuarkFromTop[0], 'jet_delR', jet_delR_list[0] )
        setattr( event.GenWZQuark[0], 'jet_delR', jet_delR_list[1] )
        setattr( event.GenWZQuark[1], 'jet_delR', jet_delR_list[2] )

        setattr( event, 'GenQ_jet_sumdelR', jet_sumdelR )

        # Write subjet_delR values to event
        setattr(event.GenBQuarkFromTop[0], 'subjet_delR', subjet_delR_list[0])
        setattr( event.GenWZQuark[0], 'subjet_delR', subjet_delR_list[1] )
        setattr( event.GenWZQuark[1], 'subjet_delR', subjet_delR_list[2] )

        setattr( event, 'GenQ_subjet_sumdelR', subjet_sumdelR )


    ########################################
    # Functions
    ########################################

    # Applies the cut criteria - returns True (survived) or False (did not survive)
    def Apply_Cut_criteria( self, candidate ):

        for ( attr, operator, cut_off ) in self.Cut_criteria:

            if not eval( '{0}{1}{2}'.format(
                getattr( candidate, attr ),
                operator,
                cut_off ) ):

                return False

        return True
    #--------------------------------------#

    # Prints the quarks and jets found in an event
    def Print_found_particles( self, event ):

        print '\nPrinting GenWZQuarks:'
        for q in event.GenWZQuark:
            print '[{0:.5f} | {1:.5f} | {2:.5f} | {3:.5f}]'.format(
                q.pt, q.eta, q.phi, q.mass )

        print '\nPrinting GenBQuarkFromTops:'
        for q in event.GenBQuarkFromTop:
            print '[{0:.5f} | {1:.5f} | {2:.5f} | {3:.5f}]'.format(
                q.pt, q.eta, q.phi, q.mass )

        print '\nPrinting Jets:'
        for q in event.good_jets:
            print '[{0:.5f} | {1:.5f} | {2:.5f} | {3:.5f}]'.format(
                q.pt, q.eta, q.phi, q.mass )

        print '\nPrinting httCandidate:'
        for q in event.httCandidate:
            print '[{0:.5f} | {1:.5f} | {2:.5f} | {3:.5f}]'.format(
                q.pt, q.eta, q.phi, q.mass )
    #--------------------------------------#

    # Deletes duplicate WZ Quarks
    def CompareWZQuarks(self, event ):

        if len(event.GenWZQuark) < 4:
            return 0

        quarks = event.GenWZQuark

        Is_Duplicate = ( quarks[-1].pt==quarks[1].pt and \
            quarks[-2].pt==quarks[0].pt )

        if Is_Duplicate:
            event.GenWZQuark = event.GenWZQuark[0:-2]

        return 0
    #--------------------------------------#

    # Gets a list of 3 quarks
    #  - The first quark is a Gen B quark. This should be the hadronic B quark.
    #    Which quark is hadronic is determined by adding the light quarks to the
    #    B quarks, and seeing which combined mass comes closer to the top mass
    #  - Output looks like: [ BQuark, lightQuark1, lightQuark2 ], where the list
    #    entries are TLorentzVector objects. 
    def Get_tl_genquarks(self, event ):

        # Check if right amount of quarks was generated
        if len(event.GenWZQuark)<2:
            self.Statistics['n_too_few_WZ'] += 1
            return 0
        if len(event.GenBQuarkFromTop)<2:
            self.Statistics['n_too_few_B'] += 1
            return 0
        elif len(event.GenWZQuark)>2:
            self.Statistics['n_too_many_WZ'] += 1
            return 0
        elif len(event.GenBQuarkFromTop)>2:
            self.Statistics['n_too_many_B'] += 1
            return 0

        # Make list of TLorentzVector objects for the 2 light quarks
        tl_GenWZQuarks = []
        for l in event.GenWZQuark:
            tl_l = ROOT.TLorentzVector()
            tl_l.SetPtEtaPhiM( l.pt, l.eta, l.phi, l.mass )
            tl_GenWZQuarks.append( tl_l )

        # Make list for the 2 B quarks, and a list of B quarks + light quarks
        tl_GenBQuarks = []
        tl_Combined = []
        for (b_i, b) in enumerate(event.GenBQuarkFromTop):

            tl_b = ROOT.TLorentzVector()
            tl_b.SetPtEtaPhiM( b.pt, b.eta, b.phi, b.mass )

            tl_GenBQuarks.append( tl_b )
            tl_Combined.append( tl_b + tl_GenWZQuarks[0] + tl_GenWZQuarks[1] )


        delmass0 = abs(tl_Combined[0].M() - self.top_mass)
        delmass1 = abs(tl_Combined[1].M() - self.top_mass)

        # Make sure the B quark with lowest del mass to top mass is at index 0
        # (both for the tl list and in the event)
        if delmass1 < delmass0:

            tl_GenBQuarks = [ tl_GenBQuarks[1], tl_GenBQuarks[0] ]

            event.GenBQuarkFromTop = [
                event.GenBQuarkFromTop[1],
                event.GenBQuarkFromTop[0] ]

        setattr(  tl_GenBQuarks[0], 'is_hadr', 1 )
        setattr(  event.GenBQuarkFromTop[0], 'is_hadr', 1 )
        setattr(  event.GenBQuarkFromTop[1], 'is_hadr', 0 )

        # Create the definitive list of 3 quarks
        tl_GenQuarks = []
        tl_GenQuarks.append( tl_GenBQuarks[0] )

        # There should be only 2 generated WZQuarks:
        tl_GenQuarks.append( tl_GenWZQuarks[0] )
        tl_GenQuarks.append( tl_GenWZQuarks[1] )
        
        return tl_GenQuarks
    #--------------------------------------#

    # Matches quarks with specified jets - works for real jets and for subjets
    def Do_delR_combinatorics( self, tl_genquarks, tl_jets ):

        n_jets = len(tl_jets)
        n_quarks = len(tl_genquarks)

        # Create delR matrix:
        Rmat = [[ (tl_genquarks[i].DeltaR( tl_jets[j] )) \
            for j in range(n_jets)] for i in range(n_quarks) ]

        """        
        print '\ndelR matrix:'
        for row in Rmat:
            for j in row:
                sys.stdout.write( '{0:.5f}'.format(j) + ' ')
            print ''
        """

        links_per_quark = [ [] for i in range(n_quarks) ]

        for i in range(n_quarks):
            for j in range(n_jets):
                if Rmat[i][j] < self.R_cut:
                    links_per_quark[i].append(j)
        
        
        # Perform some checks: see if there are 2 quarks linked to only 1 jet, and
        # check whether every quark has at least 1 jet it can be linked to

        for i in range(n_quarks):

            if len(links_per_quark[i]) == 0:
                #print 'At least one quark cannot be linked to any jet.'
                self.Statistics['n_no_jet_for_q'] += 1
                return (0,0)

            # Check if 2 quarks can only be linked to the same jet
            for j in range(n_quarks):
                if i==j: continue
                
                if len(links_per_quark[i])==1 and len(links_per_quark[j])==1 \
                    and links_per_quark[i]==links_per_quark[j]:
                    self.Statistics['n_2q_for_1j'] += 1
                    return (0,0)


        # Combinatorics: find the lowest sum of delR values
        
        sumR = 100000.0
        final_links = 0

        for (i_ind, i) in enumerate(links_per_quark[0]):
            for (j_ind, j) in enumerate(links_per_quark[1]):
                for (k_ind, k) in enumerate(links_per_quark[2]):

                    # Check if i,j,k gets a unique combination
                    if len( set( [ i, j, k ] ) ) == 3:

                        # Replace the sumR if it's smaller than the previous minimum
                        if Rmat[0][i] + Rmat[1][j] + Rmat[2][k] < sumR:
                            sumR = Rmat[0][i] + Rmat[1][j] + Rmat[2][k]
                            final_links = [ i, j, k ]

        if final_links == 0:
            print 'No unique combination found'
            self.Statistics['n_no_unique_match'] += 1
            return (0,0)

        final_delR_list = [ Rmat[i][final_links[i]] for i in range(n_quarks) ]

        return ( final_links , final_delR_list )
    #--------------------------------------#


class LeptonAnalyzer(FilterAnalyzer):
    """
    Analyzes leptons and applies single-lepton and di-lepton selection.

    Relies on TTH.MEAnalysis.VHbbTree.EventAnalyzer for inputs.

    Configuration:
    Conf.leptons[channel][cuttype] where channel=mu,ele, cuttype=tight,loose,(+veto)
    the lepton cuts must specify pt, eta and isolation cuts.

    Returns:
    event.good_leptons (list of VHbbTree.selLeptons): contains the leptons that pass the SL XOR DL selection.
        Leptons are ordered by flavour and pt.
    event.is_sl, is_dl (bool): specifies if the event passes SL or DL selection.

    """
    def __init__(self, cfg_ana, cfg_comp, looperName):
        super(LeptonAnalyzer, self).__init__(cfg_ana, cfg_comp, looperName)
        self.conf = cfg_ana._conf

    def beginLoop(self, setup):
        super(LeptonAnalyzer, self).beginLoop(setup)

        self.counters["processing"].register("sl")
        self.counters["processing"].register("dl")
        self.counters["processing"].register("slanddl")

        self.counters.addCounter("leptons")
        self.counters["leptons"].register("any")
        for l in ["mu", "el"]:
            for a in ["tight", "loose"]:
                for b in ["", "_veto"]:
                    lt = l + "_" + a + b
                    self.counters["leptons"].register(lt)

    def process(self, event):
        self.counters["processing"].inc("processed")
        self.counters["leptons"].inc("any", len(event.selLeptons))

        event.mu = filter(
            lambda x: abs(x.pdgId) == 13,
            event.selLeptons,
        )
        event.el = filter(
            lambda x: abs(x.pdgId) == 11,
            event.selLeptons,
        )

        for a in ["tight", "loose"]:
            for b in ["", "_veto"]:
                sumleps = []
                for l in ["mu", "el"]:
                    lepcuts = self.conf.leptons[l][a+b]
                    incoll = getattr(event, l)

                    leps = filter(
                        lambda x: (
                            x.pt > lepcuts["pt"]
                            and abs(x.eta) < lepcuts["eta"]
                            and abs(getattr(x, self.conf.leptons[l]["isotype"])) < lepcuts["iso"]
                        ), incoll
                    )

                    if b == "_veto":
                        good = getattr(event, "{0}_{1}".format(l, a))
                        leps = filter(lambda x: x not in good, leps)

                    if a == "tight":
                        leps = filter(
                            lambda x: x.tightId,
                            leps
                        )
                    elif a == "loose":
                        leps = filter(
                            lambda x: x.looseIdPOG,
                            leps
                        )
                    lep = sorted(leps, key=lambda x: x.pt, reverse=True)
                    sumleps += leps
                    lt = l + "_" + a + b
                    setattr(event, lt, leps)
                    setattr(event, "n_"+lt, len(leps))
                    self.counters["leptons"].inc(lt, len(leps))

                setattr(event, "lep_{0}".format(a+b), sumleps)
                setattr(event, "n_lep_{0}".format(a+b), len(sumleps))


        event.is_sl = (event.n_lep_tight == 1 and event.n_lep_tight_veto == 0)
        event.is_dl = (event.n_lep_loose == 2 and event.n_lep_loose_veto == 0)

        if event.is_sl:
            self.counters["processing"].inc("sl")
            event.good_leptons = event.mu_tight + event.el_tight
        if event.is_dl:
            self.counters["processing"].inc("dl")
            event.good_leptons = event.mu_loose + event.el_loose

        passes = event.is_sl or event.is_dl
        if event.is_sl and event.is_dl:
            #print "pathological SL && DL event: {0}".format(event)
            self.counters["processing"].inc("slanddl")
            passes = False

        if passes:
            self.counters["processing"].inc("passes")
        return passes



class JetAnalyzer(FilterAnalyzer):
    """
    Performs jet selection and b-tag counting.
    FIXME: doc
    """
    def __init__(self, cfg_ana, cfg_comp, looperName):
        super(JetAnalyzer, self).__init__(cfg_ana, cfg_comp, looperName)
        self.conf = cfg_ana._conf

    def beginLoop(self, setup):
        super(JetAnalyzer, self).beginLoop(setup)
        self.counters.addCounter("jets")
        self.counters["jets"].register("any")
        self.counters["jets"].register("good")
        for (btag_wp_name, btag_wp) in self.conf.jets["btagWPs"].items():
            self.counters["jets"].register(btag_wp_name)


    def process(self, event):
        self.counters["processing"].inc("processed")
        self.counters["jets"].inc("any", len(event.Jet))

        #pt-descending input jets
        if "input" in self.conf.general["verbosity"]:
            for j in event.Jet:
                print "ijet", j.pt, j.eta, j.phi, j.mass, j.btagCSV, j.mcFlavour

        event.good_jets = sorted(
            filter(
                lambda x: (
                    x.pt > self.conf.jets["pt"]
                    and abs(x.eta) < self.conf.jets["eta"]
                ), event.Jet
            ),
            key=lambda x: x.pt, reverse=True
            )[0:9]
        event.numJets = len(event.good_jets)
        self.counters["jets"].inc("good", len(event.good_jets))

        event.btagged_jets_bdisc = {}
        event.buntagged_jets_bdisc = {}
        for (btag_wp_name, btag_wp) in self.conf.jets["btagWPs"].items():
            algo, wp = btag_wp
            event.btagged_jets_bdisc[btag_wp_name] = filter(
                lambda x: getattr(x, algo) > wp,
                event.good_jets
            )
            event.buntagged_jets_bdisc[btag_wp_name] = filter(
                lambda x: getattr(x, algo) <= wp,
                event.good_jets
            )
            self.counters["jets"].inc(btag_wp_name,
                len(event.btagged_jets_bdisc[btag_wp_name])
            )
            setattr(event, "nB"+btag_wp_name, len(event.btagged_jets_bdisc[btag_wp_name]))
        event.buntagged_jets_bdisc = event.buntagged_jets_bdisc[self.conf.jets["btagWP"]]
        event.btagged_jets_bdisc = event.btagged_jets_bdisc[self.conf.jets["btagWP"]]
        event.n_tagwp_tagged_true_bjets = 0
        for j in event.btagged_jets_bdisc:
            if abs(j.mcFlavour) == 5:
                event.n_tagwp_tagged_true_bjets += 1
        passes = len(event.good_jets) >= 4
        if passes:
            self.counters["processing"].inc("passes")

        return passes


class BTagLRAnalyzer(FilterAnalyzer):
    """
    Performs b-tag likelihood ratio calculations
    FIXME: doc
    """
    def __init__(self, cfg_ana, cfg_comp, looperName):
        super(BTagLRAnalyzer, self).__init__(cfg_ana, cfg_comp, looperName)
        self.conf = cfg_ana._conf
        self.bTagAlgo = self.conf.jets["btagAlgo"]
        self.cplots_old = ROOT.TFile(self.conf.general["controlPlotsFileOld"])
        self.cplots = ROOT.TFile(self.conf.general["controlPlotsFile"])

        cplots_fmt = self.conf.general.get("controlPlotsFormat", "8tev")
        self.csv_pdfs_old = {
        }
        for x in ["b", "c", "l"]:
            for b in ["Bin0", "Bin1"]:
                self.csv_pdfs_old[(x, b)] = self.cplots_old.Get(
                    "csv_{0}_{1}__csv_rec".format(x, b)
                )

        self.csv_pdfs = {
        }
        for x in ["b", "c", "l"]:
            for b in ["Bin0", "Bin1"]:
                self.csv_pdfs[(x, b)] = self.cplots.Get(
                    "csv_{0}_{1}__csv_rec".format(x, b)
                )
                self.csv_pdfs[(x, b)].Scale(1.0 / self.csv_pdfs[(x, b)].Integral())
            self.csv_pdfs[(x, "pt_eta")] = self.cplots.Get(
                "csv_{0}_pt_eta".format(x)
            )
            self.csv_pdfs[(x, "pt_eta")].Scale(1.0 / self.csv_pdfs[(x, "pt_eta")].Integral())

    def get_pdf_prob(self, flavour, pt, eta, csv, kind):

        _bin = "Bin1" if abs(eta)>1.0 else "Bin0"

        if kind == "old":
            h = self.csv_pdfs_old[(flavour, _bin)]
        elif kind == "new_eta_1bin":
            h = self.csv_pdfs[(flavour, _bin)]
        elif kind == "new_pt_eta_bin_3d":
            h = self.csv_pdfs[(flavour, "pt_eta")]

        assert h != None, "flavour={0} kind={1}".format(flavour, kind)

        if csv < 0:
            csv = 0.0
        if csv > 1.0:
            csv = 1.0

        if kind == "old" or kind == "new_eta_1bin":
            nb = h.FindBin(csv)
            #if csv = 1 -> goes into overflow and pdf = 0.0
            #as a solution, take the next-to-last bin
            if nb >= h.GetNbinsX():
                nb = nb - 1
            ret = h.GetBinContent(nb)
        elif kind == "new_pt_eta_bin_3d":
            nb = h.FindBin(pt, abs(eta), csv)
            ret = h.GetBinContent(nb)
        return ret

    def beginLoop(self, setup):
        super(BTagLRAnalyzer, self).beginLoop(setup)

    def evaluate_jet_prob(self, pt, eta, csv, kind):
        return (
            self.get_pdf_prob("b", pt, eta, csv, kind),
            self.get_pdf_prob("c", pt, eta, csv, kind),
            self.get_pdf_prob("l", pt, eta, csv, kind)
        )

    def btag_likelihood(self, probs, nB, nC):

        perms = itertools.permutations(range(len(probs)))

        P = 0.0
        max_p = -1.0
        nperms = 0
        best_perm = None

        for perm in perms:
            p = 1.0

            for i in range(0, nB):
                p *= probs[perm[i]][0]
            for i in range(nB, nB + nC):
                p *= probs[perm[i]][1]
            for i in range(nB + nC, len(probs)):
                p *= probs[perm[i]][2]

            #print nperms, p, perm, max_p, best_perm
            if p > max_p:
                best_perm = perm
                max_p = p

            P += p
            nperms += 1
        P = P / float(nperms)
        assert nperms > 0
        return P, best_perm
        #end permutation loop

    def process(self, event):
        self.counters["processing"].inc("processed")

        #Take first 6 most b-tagged jets for btag LR
        jets_for_btag_lr = sorted(
            event.good_jets,
            key=lambda x: getattr(x, self.bTagAlgo),
            reverse=True,
        )[0:6]

        jet_probs = {
            kind: [
                self.evaluate_jet_prob(j.pt, j.eta, getattr(j, self.bTagAlgo), kind)
                for j in jets_for_btag_lr
            ]
            for kind in [
            "old", "new_eta_1bin",
            "new_pt_eta_bin_3d"
            ]
        }

        # for nj, j in enumerate(jets_for_btag_lr):
        #     print j.btagCSV, j.mcFlavour, jet_probs["old"][nj], jet_probs["new_eta_1bin"][nj], jet_probs["new_pt_eta_bin_3d"][nj]
        #
        jet_csvs = [
            getattr(j, self.bTagAlgo)
            for j in event.good_jets
        ]

        best_4b_perm = 0
        best_2b_perm = 0
        event.btag_lr_4b_old, best_4b_perm = self.btag_likelihood(jet_probs["old"], 4, 0)
        event.btag_lr_2b_old, best_2b_perm = self.btag_likelihood(jet_probs["old"], 2, 0)

        event.btag_lr_4b, best_4b_perm = self.btag_likelihood(jet_probs["new_eta_1bin"], 4, 0)
        event.btag_lr_2b, best_2b_perm = self.btag_likelihood(jet_probs["new_eta_1bin"], 2, 0)

        event.btag_lr_4b_alt, best_4b_perm_alt = self.btag_likelihood(jet_probs["new_pt_eta_bin_3d"], 4, 0)
        event.btag_lr_2b_alt, best_2b_perm_alt = self.btag_likelihood(jet_probs["new_pt_eta_bin_3d"], 2, 0)

        def lratio(l1, l2):
            if l1+l2>0:
                return l1/(l1+l2)
            else:
                return 0.0

        event.btag_LR_4b_2b_old = lratio(event.btag_lr_4b_old, event.btag_lr_2b_old)
        event.btag_LR_4b_2b = lratio(event.btag_lr_4b, event.btag_lr_2b)
        event.btag_LR_4b_2b_alt = lratio(event.btag_lr_4b_alt, event.btag_lr_2b_alt)
        #event.btag_LR_4b_2b_alt = 0

        event.buntagged_jets_by_LR_4b_2b = [jets_for_btag_lr[i] for i in best_4b_perm[4:]]
        event.btagged_jets_by_LR_4b_2b = [jets_for_btag_lr[i] for i in best_4b_perm[0:4]]

        for i in range(len(event.good_jets)):
            event.good_jets[i].btagFlag = 0.0

        #Jets are untagged according to the b-tagging likelihood ratio permutation
        if self.conf.jets["untaggedSelection"] == "btagLR":
            event.buntagged_jets = event.buntagged_jets_by_LR_4b_2b
            event.btagged_jets = event.btagged_jets_by_LR_4b_2b
        #Jets are untagged according to b-discriminatr
        elif self.conf.jets["untaggedSelection"] == "btagCSV":
            event.buntagged_jets = event.buntagged_jets_bdisc
            event.btagged_jets = event.btagged_jets_bdisc

        #Take first 4 most b-tagged jets
        btagged = sorted(event.btagged_jets, key=lambda x: x.btagCSV, reverse=True)[0:4]
        #Set these jets to be used as b-quarks in the MEM
        #We don't want to use more than 4 b-quarks in the hypothesis
        for jet in btagged:
            idx = event.good_jets.index(jet)
            event.good_jets[idx].btagFlag = 1.0

        passes = True
        if passes:
            self.counters["processing"].inc("passes")

        return passes

class MECategoryAnalyzer(FilterAnalyzer):
    """
    Performs ME categorization
    FIXME: doc
    """
    def __init__(self, cfg_ana, cfg_comp, looperName):
        self.conf = cfg_ana._conf
        super(MECategoryAnalyzer, self).__init__(cfg_ana, cfg_comp, looperName)
        self.cat_map = {"NOCAT":-1, "cat1": 1, "cat2": 2, "cat3": 3, "cat6":6}
        self.btag_cat_map = {"NOCAT":-1, "L": 0, "H": 1}

    def beginLoop(self, setup):
        super(MECategoryAnalyzer, self).beginLoop(setup)

        for c in ["NOCAT", "cat1", "cat2", "cat3", "cat6"]:
            self.counters["processing"].register(c)

    def process(self, event):
        self.counters["processing"].inc("processed")

        cat = "NOCAT"
        pass_btag_lr = (self.conf.jets["untaggedSelection"] == "btagLR" and
            event.btag_LR_4b_2b > self.conf.mem["btagLRCut"][event.cat]
        )
        pass_btag_csv = (self.conf.jets["untaggedSelection"] == "btagCSV" and
            len(event.btagged_jets) >= 4
        )
        cat_btag = "NOCAT"

        if pass_btag_lr or pass_btag_csv:
            cat_btag = "H"

        if event.is_sl:

            #at least 6 jets, if 6, Wtag in [60,100], if more Wtag in [72,94]
            if ((len(event.good_jets) == 6 and event.Wmass >= 60 and event.Wmass < 100) or
               (len(event.good_jets) > 6 and event.Wmass >= 72 and event.Wmass < 94)):
               cat = "cat1"
               #W-tagger fills wquark_candidate_jets
            #at least 6 jets, no W-tag
            elif len(event.good_jets) >= 6:
                cat = "cat2"
            #one W daughter missing
            elif len(event.good_jets) == 5:
                event.wquark_candidate_jets = event.buntagged_jets
                cat = "cat3"
        elif event.is_dl and len(event.good_jets)>=4:
            event.wquark_candidate_jets = []
            cat = "cat6"

        self.counters["processing"].inc(cat)
        event.cat = cat
        event.cat_btag = cat_btag
        event.catn = self.cat_map.get(cat, -1)
        event.cat_btag_n = self.btag_cat_map.get(cat, -1)

        passes = True
        if passes:
            self.counters["processing"].inc("passes")
        return passes

class WTagAnalyzer(FilterAnalyzer):
    """
    Performs W-mass calculation on pairs of untagged jets.

    Jets are considered untagged according to the b-tagging permutation which
    gives the highest likelihood of the event being a 4b+Nlight event.
    """
    def __init__(self, cfg_ana, cfg_comp, looperName):
        self.conf = cfg_ana._conf
        super(WTagAnalyzer, self).__init__(cfg_ana, cfg_comp, looperName)

    def beginLoop(self, setup):
        super(WTagAnalyzer, self).beginLoop(setup)

    def pair_mass(self, j1, j2):
        """
        Calculates the invariant mass of a two-particle system.
        """
        lv1, lv2 = [lvec(j) for j in [j1, j2]]
        tot = lv1 + lv2
        return tot.M()

    def find_best_pair(self, jets):
        """
        Finds the pair of jets whose invariant mass is closest to mW=80 GeV.
        Returns the sorted vector of [(mass, jet1, jet2)], best first.
        """
        ms = []
        
        #Keep track of index pairs already calculated
        done_pairs = set([])
        
        #Double loop over all jets
        for i in range(len(jets)):
            for j in range(len(jets)):
                
                #Make sure we haven't calculated this index pair yet
                if (i,j) not in done_pairs and i!=j:
                    m = self.pair_mass(jets[i], jets[j])
                    ms += [(m, jets[i], jets[j])]
                    
                    #M(i,j) is symmetric, hence add both pairs
                    done_pairs.add((i,j))
                    done_pairs.add((j,i))
        ms = sorted(ms, key=lambda x: abs(x[0] - 80.0))
        return ms

    def process(self, event):
        self.counters["processing"].inc("processed")

        event.Wmass = 0.0
        
        event.wquark_candidate_jets = set([])
        #Need at least 2 untagged jets to calculate W mass
        if len(event.buntagged_jets)>=2:
            bpair = self.find_best_pair(event.buntagged_jets)
            
            #Get the best mass
            event.Wmass = bpair[0][0]
            
            #All masses
            event.Wmasses = [bpair[i][0] for i in range(len(bpair))]
            
            for i in range(min(len(bpair), 2)):
                event.wquark_candidate_jets.add(bpair[i][1])
                event.wquark_candidate_jets.add(bpair[i][2])
            
            if "reco" in self.conf.general["verbosity"]:
                print "Wmass", event.Wmass, event.good_jets.index(bpair[1]), event.good_jets.index(bpair[2])
        
        #If we can't calculate W mass, untagged jets become the candidate
        else:
            for jet in event.buntagged_jets:
                event.wquark_candidate_jets.add(jet)
                
                
        passes = True
        if passes:
            self.counters["processing"].inc("passes")
        return passes

class GenRadiationModeAnalyzer(FilterAnalyzer):
    """
    Performs B/C counting
    FIXME: doc
    """
    def __init__(self, cfg_ana, cfg_comp, looperName):
        self.conf = cfg_ana._conf
        super(GenRadiationModeAnalyzer, self).__init__(cfg_ana, cfg_comp, looperName)

    def beginLoop(self, setup):
        super(GenRadiationModeAnalyzer, self).beginLoop(setup)

    def process(self, event):
        self.counters["processing"].inc("processed")

        event.nMatchSimB = 0
        event.nMatchSimC = 0
        lv_bs = map(lvec, event.GenBQuarkFromTop)
        for jet in event.good_jets:
            lv_j = lvec(jet)

            if (lv_j.Pt() > 20 and abs(lv_j.Eta()) < 2.5):
                if any([lv_b.DeltaR(lv_j) < 0.5 for lv_b in lv_bs]):
                    continue
                absid = abs(jet.mcFlavour)
                if absid == 5:
                    event.nMatchSimB += 1
                if absid == 4:
                    event.nMatchSimC += 1

        passes = True
        if passes:
            self.counters["processing"].inc("passes")
        return passes


class GenTTHAnalyzer(FilterAnalyzer):
    """
    """
    def __init__(self, cfg_ana, cfg_comp, looperName):
        self.conf = cfg_ana._conf
        super(GenTTHAnalyzer, self).__init__(cfg_ana, cfg_comp, looperName)

    def beginLoop(self, setup):
        super(GenTTHAnalyzer, self).beginLoop(setup)

    def process(self, event):
        self.counters["processing"].inc("processed")

        #Somehow, the GenWZQuark distribution is duplicated
        event.l_quarks_w = event.GenWZQuark[0:len(event.GenWZQuark)/2]
        event.b_quarks_t = event.GenBQuarkFromTop
        event.b_quarks_h = event.GenBQuarkFromH
        event.lep_top = event.GenLepFromTop
        event.nu_top = event.GenNuFromTop

        event.cat_gen = None
        event.n_cat_gen = -1

        if (len(event.lep_top) == 1 and
            len(event.nu_top) == 1 and
            len(event.l_quarks_w) == 2 and
            len(event.b_quarks_t) == 2):
            event.cat_gen = "sl"
            event.n_cat_gen = 0
        elif (len(event.lep_top) == 2 and
            len(event.nu_top) == 2 and
            len(event.l_quarks_w) == 0 and
            len(event.b_quarks_t) == 2):
            event.cat_gen = "dl"
            event.n_cat_gen = 1
        elif (len(event.lep_top) == 0 and
            len(event.nu_top) == 0 and
            len(event.l_quarks_w) == 4 and
            len(event.b_quarks_t) == 2):
            event.cat_gen = "fh"
            event.n_cat_gen = 2

        if "gen" in self.conf.general["verbosity"]:
            for j in event.l_quarks_w:
                print "q(W)", j.pt, j.eta, j.phi, j.mass, j.pdgId
            for j in event.b_quarks_t:
                print "b(t)", j.pt, j.eta, j.phi, j.mass, j.pdgId
            for j in event.lep_top:
                print "l(t)", j.pt, j.eta, j.phi, j.mass, j.pdgId
            for j in event.nu_top:
                print "n(t)", j.pt, j.eta, j.phi, j.mass, j.pdgId
            for j in event.b_quarks_h:
                print "b(h)", j.pt, j.eta, j.phi, j.mass, j.pdgId
            print "gen cat", event.cat_gen, event.n_cat_gen

        #Store for each jet, specified by it's index in the jet
        #vector, if it is matched to any gen-level quarks
        matched_pairs = {}

        def match_jets_to_quarks(jetcoll, quarkcoll, label):
            for ij, j in enumerate(jetcoll):
                for iq, q in enumerate(quarkcoll):
                    l1 = lvec(q)
                    l2 = lvec(j)
                    dr = l1.DeltaR(l2)
                    if dr < 0.3:
                        if matched_pairs.has_key(ij):
                            if matched_pairs[ij][1] > dr:
                                matched_pairs[ij] = (label, iq, dr)
                        else:
                            matched_pairs[ij] = (label, iq, dr)
        #print "GEN", len(event.GenWZQuark), len(event.GenBQuarkFromTop), len(event.GenBQuarkFromH)
        match_jets_to_quarks(event.good_jets, event.l_quarks_w, "wq")
        match_jets_to_quarks(event.good_jets, event.b_quarks_t, "tb")
        match_jets_to_quarks(event.good_jets, event.b_quarks_h, "hb")

        #Number of reco jets matched to quarks from W, top, higgs
        event.nMatch_wq = 0
        event.nMatch_tb = 0
        event.nMatch_hb = 0
        #As above, but also required to be tagged/untagged for b/light respectively.
        event.nMatch_wq_btag = 0
        event.nMatch_tb_btag = 0
        event.nMatch_hb_btag = 0

        for ij, jet in enumerate(event.good_jets):
            if not matched_pairs.has_key(ij):
                continue
            mlabel, midx, mdr = matched_pairs[ij]
            if mlabel == "wq":
                event.nMatch_wq += 1
                if jet.btagFlag < 0.5:
                    event.nMatch_wq_btag += 1
            if mlabel == "tb":
                event.nMatch_tb += 1
                if jet.btagFlag >= 0.5:
                    event.nMatch_tb_btag += 1
            if mlabel == "hb":
                event.nMatch_hb += 1
                if jet.btagFlag >= 0.5:
                    event.nMatch_hb_btag += 1

        if "matching" in self.conf.general["verbosity"]:
            matches = {"wq":event.l_quarks_w, "tb": event.b_quarks_t, "hb":event.b_quarks_h}

            for ij, jet in enumerate(event.good_jets):
                if not matched_pairs.has_key(ij):
                    continue
                mlabel, midx, mdr = matched_pairs[ij]
                print "jet match", ij, mlabel, midx, mdr, jet.pt, matches[mlabel][midx].pt

        passes = True
        if passes:
            self.counters["processing"].inc("passes")
        return passes

class MEAnalyzer(FilterAnalyzer):
    """
    Performs ME calculation using external integrator

    Relies on:
    event.good_jets, event.good_leptons, event.cat, event.input.met_pt

    Produces:
    p_hypo_tth (double): probability for the tt+H(bb) hypothesis
    p_hypo_ttbb (double): probability for the tt+bb hypothesis

    """
    def __init__(self, cfg_ana, cfg_comp, looperName):
        self.conf = cfg_ana._conf
        super(MEAnalyzer, self).__init__(cfg_ana, cfg_comp, looperName)


        self.configs = {
            "default": MEM.MEMConfig(),
            "NumPointsDouble": MEM.MEMConfig(),
            "NumPointsHalf": MEM.MEMConfig(),
            "NoJacobian": MEM.MEMConfig(),
            "NoDecayAmpl": MEM.MEMConfig(),
            "NoPDF": MEM.MEMConfig(),
            "NoScattAmpl": MEM.MEMConfig(),
            "QuarkEnergy98": MEM.MEMConfig(),
            "QuarkEnergy10": MEM.MEMConfig(),
            "NuPhiRestriction": MEM.MEMConfig(),
            "JetsPtOrder": MEM.MEMConfig(),
            "JetsPtOrderIntegrationRange": MEM.MEMConfig(),
        }

        self.memkeys = self.conf.mem["methodsToRun"]

        self.configs["default"].defaultCfg()
        self.configs["NumPointsDouble"].defaultCfg(2.0)
        self.configs["NumPointsHalf"].defaultCfg(0.5)
        self.configs["NoJacobian"].defaultCfg()
        self.configs["NoDecayAmpl"].defaultCfg()
        self.configs["NoPDF"].defaultCfg()
        self.configs["NoScattAmpl"].defaultCfg()
        self.configs["QuarkEnergy98"].defaultCfg()
        self.configs["QuarkEnergy10"].defaultCfg()
        self.configs["NuPhiRestriction"].defaultCfg()
        self.configs["JetsPtOrder"].defaultCfg()
        self.configs["JetsPtOrderIntegrationRange"].defaultCfg()

        self.configs["NoJacobian"].int_code &= ~ MEM.IntegrandType.Jacobian
        self.configs["NoDecayAmpl"].int_code &= ~ MEM.IntegrandType.DecayAmpl
        self.configs["NoPDF"].int_code &= ~ MEM.IntegrandType.PDF
        self.configs["NoScattAmpl"].int_code &=  ~ MEM.IntegrandType.ScattAmpl
        self.configs["QuarkEnergy98"].j_range_CL = 0.98
        self.configs["QuarkEnergy98"].b_range_CL = 0.98
        self.configs["QuarkEnergy10"].j_range_CL = 0.10
        self.configs["QuarkEnergy10"].b_range_CL = 0.10
        self.configs["NuPhiRestriction"].m_range_CL = 99
        self.configs["JetsPtOrder"].highpt_first  = 0
        self.configs["JetsPtOrderIntegrationRange"].highpt_first  = 0
        self.configs["JetsPtOrderIntegrationRange"].j_range_CL = 0.99
        self.configs["JetsPtOrderIntegrationRange"].b_range_CL = 0.99

        #Create the ME integrator.
        #Arguments specify the verbosity
        self.integrator = MEM.Integrand(
            #0,
            MEM.output,
            #MEM.output | MEM.init | MEM.event,# | MEM.integration,
            self.configs["default"]
        )

        #Create an emtpy std::vector<MEM::Permutations::Permutations>
        self.permutations = CvectorPermutations()

        #Assume that only jets passing CSV>0.5 are b quarks
        self.permutations.push_back(MEM.Permutations.BTagged)

        #Assume that only jets passing CSV<0.5 are l quarks
        self.permutations.push_back(MEM.Permutations.QUntagged)

        self.integrator.set_permutation_strategy(self.permutations)

        #Pieces of ME to calculate
        self.integrator.set_integrand(
            MEM.IntegrandType.Constant
            |MEM.IntegrandType.ScattAmpl
            |MEM.IntegrandType.DecayAmpl
            |MEM.IntegrandType.Jacobian
            |MEM.IntegrandType.PDF
            |MEM.IntegrandType.Transfer
        )
        self.integrator.set_sqrts(13000.);

        #Create an empty vector for the integration variables
        self.vars_to_integrate = CvectorPSVar()

    def add_obj(self, objtype, **kwargs):
        """
        Add an event object (jet, lepton, MET) to the ME integrator.

        objtype: specifies the object type
        kwargs: p4s: spherical 4-momentum (pt, eta, phi, M) as a tuple
                obsdict: dict of additional observables to pass to MEM
        """
        if kwargs.has_key("p4s"):
            pt, eta, phi, mass = kwargs.pop("p4s")
            v = ROOT.TLorentzVector()
            v.SetPtEtaPhiM(pt, eta, phi, mass);
        elif kwargs.has_key("p4c"):
            v = ROOT.TLorentzVector(*kwargs.pop("p4c"))
        obsdict = kwargs.pop("obsdict", {})

        o = MEM.Object(v, objtype)
        for k, v in obsdict.items():
            o.addObs(k, v)
        self.integrator.push_back_object(o)

    def beginLoop(self, setup):
        super(MEAnalyzer, self).beginLoop(setup)

    def process(self, event):
        self.counters["processing"].inc("processed")

        #Clean up any old MEM state
        self.vars_to_integrate.clear()
        self.integrator.next_event()

        #Initialize members for tree filler
        event.mem_results_tth = []
        event.mem_results_ttbb = []

        #jets = sorted(event.good_jets, key=lambda x: x.pt, reverse=True)
        leptons = event.good_leptons
        met_pt = event.input.met_pt
        met_phi = event.input.met_phi

        if "reco" in self.conf.general["verbosity"]:
            for j in jets:
                print "jet", j.pt, j.eta, j.phi, j.mass, j.btagCSV, j.btagFlag, j.mcFlavour
            for l in leptons:
                print "lep", l.pt, l.eta, l.phi, l.mass, l.charge

        #Check if event passes reco-level requirements to calculate ME
        if event.cat in self.conf.mem["MECategories"] and event.cat_btag == "H":
            print "MEM RECO PASS", (event.input.run, event.input.lumi, event.input.evt,
                event.cat, event.btag_LR_4b_2b, len(event.btagged_jets),
                len(event.wquark_candidate_jets), len(leptons),
                len(event.btagged_jets), len(event.buntagged_jets)
            )
        else:
            #Don't calculate ME
            return True

        #Get the conf dict specifying which matches we require
        required_match = self.conf.mem.get("requireMatched", {}).get(event.cat, {})

        #Check if event.nMatch_label >= conf.mem[cat][label]
        def require(label):
            nreq = required_match.get(label, None)
            if nreq is None:
                return True

            nmatched = getattr(event, "nMatch_"+label)

            #In case event did not contain b form higgs (e.g. ttbb)
            if "hb" in label and len(event.GenBQuarkFromH) < nreq:
                return True

            passes = (nmatched >= nreq)
            return passes

        #Calculate all the match booleans
        passd = {
            p: require(p) for p in ["wq", "wq_btag", "tb", "tb_btag", "hb", "hb_btag"]
        }

        for k, v in passd.items():
            if not v:
                #print "Failed to match", k
                return True

        def add_objects():
            self.vars_to_integrate.clear()
            self.integrator.next_event()
            #One quark from W missed, integrate over its direction
            if event.cat in ["cat2", "cat3"]:
                self.vars_to_integrate.push_back(MEM.PSVar.cos_qbar1)
                self.vars_to_integrate.push_back(MEM.PSVar.phi_qbar1)
            for jet in event.btagged_jets:
                self.add_obj(
                    MEM.ObjectType.Jet,
                    p4s=(jet.pt, jet.eta, jet.phi, jet.mass),
                    obsdict={MEM.Observable.BTAG: jet.btagFlag}
                )
            for jet in event.wquark_candidate_jets:
                self.add_obj(
                    MEM.ObjectType.Jet,
                    p4s=(jet.pt, jet.eta, jet.phi, jet.mass),
                    obsdict={MEM.Observable.BTAG: jet.btagFlag}
                )
            for lep in leptons:
                self.add_obj(
                    MEM.ObjectType.Lepton,
                    p4s=(lep.pt, lep.eta, lep.phi, lep.mass),
                    obsdict={MEM.Observable.CHARGE: lep.charge}
                )
            self.add_obj(
                MEM.ObjectType.MET,
                #MET is caused by massless object
                p4s=(met_pt, 0, met_phi, 0),
            )

        fstate = MEM.FinalState.TTH
        if len(leptons) == 2:
            fstate = MEM.FinalState.LL
        if len(leptons) == 1:
            fstate = MEM.FinalState.LH

        res = {}
        for hypo in [MEM.Hypothesis.TTH, MEM.Hypothesis.TTBB]:
            for confname in self.memkeys:
                if self.conf.mem["calcME"]:
                    conf = self.configs[confname]
                    #print "MEM conf", confname, "hypo", hypo
                    print "MEM conf", hypo, confname
                    self.integrator.set_cfg(conf)
                    add_objects()
                    r = self.integrator.run(
                        fstate,
                        hypo,
                        self.vars_to_integrate
                    )
                    res[(hypo, confname)] = r
                else:
                    r = MEM.MEMOutput()
                    res[(hypo, confname)] = r

        p1 = res[(MEM.Hypothesis.TTH, "default")].p
        p2 = res[(MEM.Hypothesis.TTBB, "default")].p

        #In case of an erroneous calculation, print out event kinematics
        if self.conf.mem["calcME"] and (p1<=0 or p2<=0 or (p1 / (p1+0.02*p2))<0.0001):
            print "MEM BADPROB", p1, p2

        #print out full MEM result dictionary
        #print "RES", [(k, res[k].p) for k in sorted(res.keys())]

        event.mem_results_tth = [res[(MEM.Hypothesis.TTH, k)] for k in self.memkeys]
        event.mem_results_ttbb = [res[(MEM.Hypothesis.TTBB, k)] for k in self.memkeys]
