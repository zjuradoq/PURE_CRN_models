from biocrnpyler import Mechanism, Species, Complex, Reaction, Parameter, Mixture, Component

class Transcription_Initiation(Mechanism):
    """
    Creates all the species and reaction around the nucleotide sequences 
    in the given transcribed region, common protein and small molecules.
    This mechanism has reactions needed for TX regardless of the sequence defined.
    """
    def __init__(self, rnap: Species, active_nucleotides: dict, 
                 small_molecules: dict, name="transcription_initiation"):
        """Initializes TX process

        Args:
            rnap (Species): The RNA polymerase enzyme 
            active_nucleotides (dict): All active nucleotides in PURE: ATP, GTP, CTP, UTP
                                       stored as a dictionary with values as Species.
            small_molecules (dict): All small molecules needed for initiation of transcription (GDP, PO4)
                                    stored as a dictionary with values as Species.
            name (str, optional): The name of the mechanism. Defaults to "transcription_initiation".

        Raises:
            ValueError: If species arguments (rnap, elements of `active_nucleotides` and `small_molecules`) 
                        are not of Species class in BioCRNpyler.  

        """
        if isinstance(rnap, Species):
            self.rnap = rnap
        else:
            raise ValueError("'rnap' argument must be a Species.")
        for nuc in active_nucleotides:
            if not isinstance(active_nucleotides[nuc], Species):
                raise ValueError("All values of 'active_nucleotides' must be Species.")
        self.active_nucleotides = active_nucleotides
        for sm in small_molecules:
            if not isinstance(small_molecules[sm], Species):
                raise ValueError("All elements of 'small_molecules' must be Species.")
        self.small_molecules = small_molecules
        Mechanism.__init__(self=self, name=name, mechanism_type="tx_initiation")
    
    def update_species(self, dna, **keywords):
        species = [dna, self.rnap]
        species += list(self.active_nucleotides.values)
        species += list(self.small_molecules.values)
        # Create other complexes
        RNAPa_bound_dna = Complex([self.rnap, dna])
        RNAPa_bound_GTP = Complex([self.rnap, self.active_nucleotides["GTP"]])
        RNAPa_bound_GDP_PO4 = Complex([self.rnap, self.small_molecules["GDP"], self.small_molecules["PO4"]])
        species += [RNAPa_bound_dna, RNAPa_bound_GTP, RNAPa_bound_GDP_PO4]
        # Attributes for nucleotides and small molecules
        self.ATP = self.active_nucleotides["ATP"]
        self.GTP = self.active_nucleotides["GTP"]
        self.CTP = self.active_nucleotides["CTP"]
        self.UTP = self.active_nucleotides["UTP"]
        self.GDP = self.small_molecules["GDP"]
        self.PO4 = self.small_molecules["PO4"]
        return species
    
    def update_reactions(self, dna, component, part_id = None, **keywords):
        # Get parameters
        if part_id == None and component != None:
            part_id = component.name
        NTP_deg = component.get_parameter("NTP_deg", part_id = part_id, mechanism = self)
        k_rnapbF = component.get_parameter("k_rnapbF", part_id = part_id, mechanism = self)
        k_rnapbF2 = component.get_parameter("k_rnapbF2", part_id = part_id, mechanism = self)
        k_rnapbF3 = component.get_parameter("k_rnapbF3", part_id = part_id, mechanism = self)
        # Create complexes
        RNAPa_bound_dna = Complex([self.rnap, dna])
        RNAPa_bound_GTP = Complex([self.rnap, self.GTP])
        RNAPa_bound_GDP_PO4 = Complex([self.rnap, self.GDP, self.PO4])
        # Create reactions
        r1 = Reaction.from_massaction([self.ATP], [], k_forward=NTP_deg), 
        r2 = Reaction.from_massaction([self.GTP], [], k_forward=NTP_deg), 
        r3 = Reaction.from_massaction([self.CTP], [], k_forward=NTP_deg), 
        r4 = Reaction.from_massaction([self.UTP], [], k_forward=NTP_deg), 
        #Binding of RNAP and beginning of transcription
        r5 = Reaction.from_massaction([self.rnap, dna, self.GTP], [RNAPa_bound_GTP], k_forward=k_rnapbF),
        r6 = Reaction.from_massaction([RNAPa_bound_GTP], [RNAPa_bound_GDP_PO4], k_forward=k_rnapbF2),
        r7 = Reaction.from_massaction([RNAPa_bound_GDP_PO4], [RNAPa_bound_dna, self.GDP, self.PO4,], k_forward=k_rnapbF3)] 
        all_reactions = [r1, r2, r3, r4, r5, r6, r7]
        return all_reactions
    
class Transcription_Elongation(Mechanism):
    ### This mechanism will create all the species and reactions needed for the elongation of the transcript
    ### Bring most reactions from "Addition of each nucleotide and individual reactions needed" here.
    pass

class Transcription_Termination(Mechanism):
    ### This mechanism will create all the species and reactions needed for the termination of transcription
    ### Add all termination reactions here.
    pass