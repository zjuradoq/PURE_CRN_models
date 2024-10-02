from bioCRNpyler import Mixture, Protein
from pure_transcription_mechanisms import *
class PURE_Transcription(Mixture):
    """A Mixture for transcription of a gene in PURE cell-free system. 
    The mixture models resources such as RNA polymerase, all nucleotides, and small molecules.
    """
    def __init__(self, name="pure_transcription", rnap="RNAP",
                 ATP="ATP", GTP="GTP", CTP="CTP", UTP="UTP", 
                 GDP="GDP", PO4="PO4", **keywords):
        """Initializes a PURE transcription mixture instance

        Args:
            name (str, optional): The name of the mixture. Defaults to "pure_transcription".
            rnap (str, optional): Name of the RNA polymerase enzyme. Defaults to "RNAP".
            ATP (str, optional): Name of the ATP species, an active nucleotide in PURE. Defaults to "ATP".
            GTP (str, optional): Name of the GTP species, an active nucleotide in PURE. Defaults to "GTP". 
            CTP (str, optional): Name of the CTP species, an active nucleotide in PURE. Defaults to "CTP".
            UTP (str, optional): Name of the UTP species, an active nucleotide in PURE. Defaults to "UTP".
            GDP (str, optional): Name of the GDP species, a small molecule in PURE. Defaults to "GDP".
            PO4 (str, optional): Name of the PO4 species, a small molecule in PURE. Defaults to "PO4".
        """
        Mixture.__init__(self, name=name, **keywords)
        # Create Components for transcription
        self.rnap = Protein(rnap)
        self.ATP = Protein("ATP")
        self.GTP = Protein("GTP")
        self.CTP = Protein("CTP")
        self.UTP = Protein("UTP")
        self.GDP = Protein("GDP")
        self.PO4 = Protein("PO4")
        self.add_components([self.rnap, self.ATP, self.GTP, self.CTP, self.UTP, self.GDP, self.PO4])
        # Add attributes
        self.rnap.add_attribute("machinery")
        self.ATP.add_attribute("nucleotide")
        self.CTP.add_attribute("nucleotide")
        self.GTP.add_attribute("nucleotide")
        self.UTP.add_attribute("nucleotide")
        self.GDP.add_attribute("small_molecule")
        self.PO4.add_attribute("small_molecule")
        # Create mechanisms
        tx_initiation = Transcription_Initiation(rnap=self.rnap.get_species(),
                                                 active_nucleotides={"ATP":self.ATP.get_species(),
                                                                     "GTP":self.GTP.get_species(),
                                                                     "CTP":self.CTP.get_species(),
                                                                     "UTP":self.UTP.get_species()},
                                                 small_molecules={"GDP":self.GDP.get_species(),
                                                                  "PO4":self.PO4.get_species()})
        tx_elongation = Transcription_Elongation(rnap=self.rnap.get_species())
        tx_termination = Transcription_Termination(rnap=self.rnap.get_species())
        # Add mechanisms
        default_mechanisms = {tx_initiation.mechanism_type:tx_initiation,
                              tx_elongation.mechanism_type:tx_elongation,
                              tx_termination.mechanism_type:tx_termination}
        self.add_mechanisms(default_mechanisms)
        