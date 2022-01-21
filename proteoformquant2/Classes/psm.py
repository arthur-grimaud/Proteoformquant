
from Utils import constant
from Utils.misc import truncate
import pprint
import spectrum_utils.spectrum as sus
class Psm():

    def __init__(self, rank, spectrum,  identificationItem):
        
        self.Modification = []

        #PSM characteristics
        for key in identificationItem:
            setattr(self, key, identificationItem[key])
        
        self.rank: int = rank
        self.spectrum = spectrum
        self.proteoform = None

        self.annotFrag = None
        #PSM quantification
        self.isValidated = False
        self.ratio: float = 0.0 

        # pprint.pprint(vars(self))
        pass


    
    #Getters

    def getPeptideSequence(self):
        return self.PeptideSequence

    def getChargeState(self):
        return self.chargeState
    
    def getCalculatedMassToCharge(self):
        return self.calculatedMassToCharge

    def getModificationsBrno(self):
        """Returns modifications of a peptide in brno nomenclature: K27meK34ac"""
        deltaMods = constant.delta_mod

        if self.Modification == []:
            return "Unmodified"
        mod_list = []

        for mod in self.Modification:
            m = str(truncate(mod["monoisotopicMassDelta"],2))
            mod_type = [value for key, value in deltaMods.items() if key.lower() in [m]]

            if mod_type != []:
                mod_type = mod_type[0]
                mod_list.append( str( self.PeptideSequence[ int(mod["location"])-1 ] ) + str(mod["location"]) + mod_type )
            else:
                print("Delta mass to brno conversion failed, {0} not found".format(mod["monoisotopicMassDelta"]))
                mod_list.append( str( self.PeptideSequence[ int(mod["location"])-1 ] ) + str(mod["location"]) + "na" )
    
        return "".join(mod_list)

    def getModificationsSu(self):
        """Returns a dictionnary of modifications in the spectrum utils format: {'position(0-based)': 'mass shift'}"""
        modDict = {}
        print(self.proteoform)
        for mod in self.proteoform.getModificationDict():
            modDict[int(mod["location"])]=float(mod["monoisotopicMassDelta"]) #use list comprehension ??

        print(modDict)
        return modDict

    #Setters

    def setProteoform(self, proteoform):
        self.proteoform = proteoform

    def setAnnotatedFragments(self, fragTol=1, maxCharge=1):
            
            id = self.spectrum.getId()
            fragIntens = self.spectrum.getFragIntens()
            fragMz= self.spectrum.getFragMz()

            calculatedMassToCharge = self.getCalculatedMassToCharge()
            chargeState = self.getChargeState()
            peptideSequence = self.getPeptideSequence()
            modifications =self.getModificationsSu()
            mods_brno = self.getModificationsBrno()


            
            spectrumSu = sus.MsmsSpectrum(
                id,
                calculatedMassToCharge, #might cause an issue when charge state != 0
                chargeState,
                fragMz,
                fragIntens,
                peptide = peptideSequence,
                modifications = modifications )
            
            for theoFragType in self.proteoform.getTheoFrag():
                for frag in  self.proteoform.getTheoFrag()[theoFragType].items(): #TODO change this format? 

                    try:
                        spectrumSu = spectrumSu.annotate_mz_fragment(fragment_mz = frag[1], fragment_charge = maxCharge,
                                            fragment_tol_mass = fragTol, fragment_tol_mode = "Da", text = frag[0])
                    except (ValueError):
                        pass

                
                try:
                    mzList = []
                    intenList = []
                    posList = []
                
                    for (mz,intensity,annotation) in zip(spectrumSu.mz, spectrumSu.intensity, spectrumSu._annotation):
                        
                        if annotation != None: #if the fragment has been matched/annotated
                            mzList.append(mz)
                            intenList.append(intensity)
                            annotationList = str(annotation.annotation).split(sep=",")
                            posList.append(str(annotation[0])+ ":" + str(annotation[1]))



                except(TypeError):
                    print("No {0} matched ions or incorrect annotation for spectrum: {1} , psm: {2}".format(theoFragType, id, self.proteoform.modificationBrno))
                    pass

