#Mandatory: List of processes
processList = {
    'p8_ee_ZZ_ecm240':{'fraction':0.000000001, 'chunks':0.01, 'output':'testOut'}
}

#Mandatory: Production tag when running over EDM4Hep centrally produced events, this points to the yaml files for getting sample statistics
prodTag     = "FCCee/winter2023/IDEA/"

#Optional: output directory, default is local running directory
outputDir   = "."

#Optional
nCPUS       = 8
runBatch    = False
#batchQueue = "longlunch"
#compGroup = "group_u_FCC.local_gen"

#Optional test file
testFile ="root://eospublic.cern.ch//eos/experiment/fcc/ee/generation/DelphesEvents/spring2021/IDEA/p8_ee_ZH_ecm240/events_101027117.root"

import ROOT
ROOT.gInterpreter.Declare("""                                                                                                                             
    ROOT::VecOps::RVec<float> get_Dr(std::vector<fastjet::PseudoJet> J, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in, vector<vector<int> > RP_idx)
    {
    std::vector<fastjet::ClusterSequence> cluster_sequences;
    cluster_sequences.clear();
    std::vector<fastjet::PseudoJet> reclus_jets;
    std::vector<fastjet::PseudoJet> subjets;
    ROOT::VecOps::RVec<float> lundVals;

    fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, fastjet::JetDefinition::max_allowable_R);
    for (unsigned int iJet=0; iJet < 1; iJet++)
      {
       for ( unsigned int iconst=0; iconst< RP_idx.at(iJet).size(); iconst++ )
       { 
         TLorentzVector const_tlv;
         auto & p = in[RP_idx[iJet][iconst]];
         const_tlv.SetXYZM(p.momentum.x, p.momentum.y, p.momentum.z, p.mass);
         subjets.push_back(fastjet::PseudoJet());
         subjets.back().reset_PtYPhiM(const_tlv.Pt(), const_tlv.Rapidity(), const_tlv.Phi(), const_tlv.M());
       }

      cluster_sequences.push_back(fastjet::ClusterSequence(subjets, jet_def));      
      reclus_jets.push_back(sorted_by_pt(cluster_sequences.back().inclusive_jets())[0]);
          
      fastjet::PseudoJet leadJet = reclus_jets[0];                                                                                                         
      fastjet::JetDefinition jd(fastjet::cambridge_algorithm, fastjet::JetDefinition::max_allowable_R);
      fastjet::Recluster rc(jd);                                                                                                                           
      fastjet::PseudoJet j = rc.result(leadJet);                                                                                                    

      fastjet::PseudoJet jj, j1, j2;                                                                                                                       
      jj = j;                                                                                                                                              
      while (jj.has_parents(j1,j2)) {                                                                                                                      
        if (j1.pt2() < j2.pt2()) {                                                                                                                         
        fastjet::PseudoJet jTemp;                                                                                                                          
        jTemp = j1;                                                                                                                                        
        j1 = j2;                                                                                                                                           
        j2 = jTemp;                                                                                                                                        
      }                                                                                                                                                    
      lundVals.push_back( std::log(1/j1.delta_R(j2)) );                                                                            
    } 
   }
    return lundVals;    
  }
""")

#Mandatory: RDFanalysis class where the use defines the operations on the TTree
class RDFanalysis():

    #__________________________________________________________
    #Mandatory: analysers funtion to define the analysers to process, please make sure you return the last dataframe, in this example it is df2
    def analysers(df):
        df2 = (df
               #############################################
               ##          Aliases for # in python        ##
               #############################################
               #define the RP px, py, pz and e
               .Define("RP_px",          "ReconstructedParticle::get_px(ReconstructedParticles)")
               .Define("RP_py",          "ReconstructedParticle::get_py(ReconstructedParticles)")
               .Define("RP_pz",          "ReconstructedParticle::get_pz(ReconstructedParticles)")
               .Define("RP_e",           "ReconstructedParticle::get_e(ReconstructedParticles)")
               .Define("RP_m",           "ReconstructedParticle::get_mass(ReconstructedParticles)")
               .Define("RP_q",           "ReconstructedParticle::get_charge(ReconstructedParticles)")

               #build pseudo jets with the RP, using the interface that takes px,py,pz,m for better
               #handling of rounding errors
               .Define("pseudo_jets",    "JetClusteringUtils::set_pseudoJets_xyzm(RP_px, RP_py, RP_pz, RP_m)")
               #.Define("pseudo_jets2",    "JetClusteringUtils::set_pseudoJets(RP_px, RP_py, RP_pz, RP_e)")
               
               #run jet clustering with all reconstructed particles. valencia_algorithm, R=0.5, inclusive clustering, E-scheme
               .Define("FCCAnalysesJets_valencia", "JetClustering::clustering_valencia(0.5, 1, 6, 0, 0, 1., 1.)(pseudo_jets)")

               #get the jets out of the struct
               .Define("jets_valencia",           "JetClusteringUtils::get_pseudoJets(FCCAnalysesJets_valencia)")
               #get the jets constituents out of the struct
               .Define("jetconstituents_valencia","JetClusteringUtils::get_constituents(FCCAnalysesJets_valencia)")
               #get some variables
               .Define("jets_valencia_eta",         "JetClusteringUtils::get_eta(jets_valencia)")
               .Define("jets_valencia_phi",         "JetClusteringUtils::get_phi(jets_valencia)")
               .Define("jets_valencia_e",           "JetClusteringUtils::get_e(jets_valencia)")
               .Define("jets_valencia_pt",          "JetClusteringUtils::get_pt(jets_valencia)")

               .Define("FCCAnalysesJets_cambridge", "JetClustering::clustering_cambridge(0.5, 2, 4, 0, 0)(pseudo_jets)")
               #get the jets out of the struct                                                                                                            
               .Define("jets_cambridge",            "JetClusteringUtils::get_pseudoJets(FCCAnalysesJets_cambridge)")
               #get the jets constituents out of the struct                                          
                                                       
               .Define("jetconstituents_cambridge", "JetClusteringUtils::get_constituents(FCCAnalysesJets_cambridge)")
               #get some variables                                                                                                                         
               .Define("jets_cambridge_pt",         "JetClusteringUtils::get_pt(jets_cambridge)")
               .Define("jets_cambridge_eta",        "JetClusteringUtils::get_eta(jets_cambridge)")
               .Define("jets_cambridge_phi",        "JetClusteringUtils::get_phi(jets_cambridge)")
               .Define("jets_cambridge_e",          "JetClusteringUtils::get_e(jets_cambridge)")
               .Define("lnInvDr",                   "get_Dr(pseudo_jets, ReconstructedParticles, jetconstituents_valencia)")
              )
        return df2

    #__________________________________________________________
    #Mandatory: output function, please make sure you return the branchlist as a python list
    def output():
        branchList = [
                # only store the constituents kinematic quantities for the first jet
                "RP_px", "RP_py", "RP_pz", "RP_e", "RP_m", "RP_q",
                "jets_valencia_pt",
                "jets_valencia_eta",
                "jets_valencia_phi",
                "jets_valencia_e",
                "jetconstituents_valencia",
                "jets_cambridge_pt",
                "jets_cambridge_eta",
                "jets_cambridge_phi",
                "jets_cambridge_e",
                "jetconstituents_cambridge",
                "lnInvDr",
                ]
        return branchList
