import unittest
import pandas as pd
import pickle
import pathlib
import os
import shutil
from BFAIR.atom_mapping import (MolfileDownloader,
                                write_rxn_files,
                                obtain_atom_mappings,
                                parse_reaction_mappings,
                                parse_metabolite_mappings,
                                generate_INCA_mapping_input,
                                check_symmetry,
                               )

original_wd = os.getcwd()
current_dir = str(pathlib.Path(__file__).parent.absolute())

class test_methods(unittest.TestCase):
    
    # Create method to compare dataframes
    def assertDataframeEqual(self, a, b, msg):
        try:
            pd.testing.assert_frame_equal(a, b)
        except AssertionError as e:
            raise self.failureException(msg) from e

    @classmethod
    def setUpClass(cls):
        # Working directory changed to make sure that all generated
        # test data is contained properly.
        os.chdir(current_dir)
        file_obj = open(current_dir + '/test_data/atom_mapping_Data/test_data.obj', 'rb')
        (
        metabolites,
        unmapped_rxns,
        mapped_rxns,
        reaction_df,
        metabolite_df,
        reaction_df_csv,
        metabolite_df_csv,
        model_reaction_df,
        model_metabolite_df,
        ) = pickle.load(file_obj)

        file_obj.close()
        
        cls.metabolites = metabolites
        cls.unmapped_rxns = unmapped_rxns
        cls.mapped_rxns = mapped_rxns
        cls.reaction_df = reaction_df
        cls.metabolite_df = metabolite_df
        cls.reaction_df_csv = reaction_df_csv
        cls.metabolite_df_csv = metabolite_df_csv
        cls.model_reaction_df = model_reaction_df
        cls.model_metabolite_df = model_metabolite_df
        
        
    @classmethod   
    def tearDownClass(cls):
        shutil.rmtree(current_dir + '/metabolites')
        shutil.rmtree(current_dir + '/unmappedRxns')
        shutil.rmtree(current_dir + '/mappedRxns')
        os.remove(current_dir + '/MappingReactions.csv')
        os.remove(current_dir + '/MappingMetabolites.csv')
        os.chdir(original_wd)
        
        
    def setUp(self):
        # Add the method to compare dataframes in the class
        self.addTypeEqualityFunc(pd.DataFrame, self.assertDataframeEqual)
        
        
    # Currently, tests are supposed to be run in alphabetical order
    # i.e. "test_a_" first, then "test_b_", because their output 
    # depends on the output of other tests.
    def test_g_check_symmetry(self):
        """ Tests check_symmetry() function 
        Makes sure that symmetric and non-symmetric compounds
        are reported correctly """
        self.assertTrue(check_symmetry('succ_c.mol'))
        self.assertFalse(check_symmetry('f6p_c.mol'))
        
    
    def test_f_generate_INCA_mapping_input(self):
        """ Tests generate_INCA_mapping_input() function
        Compares the resulting output to the previously generated results """
        reaction_df_csv = self.reaction_df_csv
        metabolite_df_csv = self.metabolite_df_csv
        reaction_df_ = parse_reaction_mappings()
        metabolite_df_ = parse_metabolite_mappings()
        
        generate_INCA_mapping_input(reaction_df_, metabolite_df_)
        
        reaction_df_csv_ = pd.read_csv(current_dir + '/MappingReactions.csv')
        metabolite_df_csv_ = pd.read_csv(current_dir + '/MappingMetabolites.csv')
        
        self.assertEqual(reaction_df_csv, reaction_df_csv_)
        self.assertEqual(metabolite_df_csv, metabolite_df_csv_)    
        
    
    def test_e_parse_metabolite_mappings(self):
        """ Tests parse_metabolite_mappings() function.
        Compares the resulting output to the previously generated results """
        metabolite_df = self.metabolite_df
        metabolite_df_ = parse_metabolite_mappings()
        self.assertEqual(metabolite_df, metabolite_df_)
            
        
    def test_d_parse_reaction_mappings(self):
        """ Tests parse_reaction_mappings function.
        Compares the resulting output to the previously generated results """
        reaction_df = self.reaction_df
        reaction_df_ = parse_reaction_mappings()
        self.assertEqual(reaction_df, reaction_df_)    
                   
            
    def test_c_obtain_atom_mappings(self):
        """ Tests obtain_atom_mapppings() function.
        Compares the resulting output to the previously generated results """
        mapped_rxns = self.mapped_rxns
        
        obtain_atom_mappings()
        
        mapped_rxns_ = os.listdir(current_dir + '/mappedRxns/rxnFiles')
        for i, rxn_file in enumerate(mapped_rxns_):
            with open(current_dir + f'/mappedRxns/rxnFiles/{rxn_file}', 'r') as f:
                lines = f.readlines()
                atom_rows = []
                for j, line in enumerate(lines):
                    if len(line.split()) in (15,16):
                        # Only append rows containing atom mappings
                        atom_rows.append(line.split())
                mapped_rxns_[i] = atom_rows
                
        self.assertEqual(mapped_rxns, mapped_rxns_)
        
    
    def test_b_write_rxn_files(self):
        """ Tests write_rxn_files function.
        Compares the resulting output to the previously generated results """
        unmapped_rxns = self.unmapped_rxns
        model_reaction_df = self.model_reaction_df
        write_rxn_files(model_reaction_df)
        
        unmapped_rxns_ = os.listdir(current_dir + '/unmappedRxns')
        for i, rxn_file in enumerate(unmapped_rxns_):
            with open(current_dir + f'/unmappedRxns/{rxn_file}', 'r') as f:
                 unmapped_rxns_[i] = f.readlines()
                    
        self.assertEqual(unmapped_rxns, unmapped_rxns_)    
        
        
    def test_a_MolfileDownloader(self):
        """ Tests the MolfileDownloader class with all of its' methods.
        Compares downloaded files to the ones downloaded previously """
        metabolites = self.metabolites
        model_metabolite_df = self.model_metabolite_df
        downloader = MolfileDownloader(model_metabolite_df)
        downloader.generate_molfile_database()
        
        metabolites_ = os.listdir(current_dir + '/metabolites')
        for i, molfile in enumerate(metabolites_):
            with open(current_dir + f'/metabolites/{molfile}', 'r') as f:
                metabolites_[i] = f.readlines()
                
        self.assertEqual(metabolites, metabolites_)
        
        
if __name__ == "__main__":
    unittest.main()
