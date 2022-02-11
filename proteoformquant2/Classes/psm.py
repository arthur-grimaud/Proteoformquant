
from Utils import constant
import unimod_mapper 
um = unimod_mapper.UnimodMapper()
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

        self.annotation = {}

        #PSM quantification
        if rank == 1:
            self.isValidated = True
            self.ratio: float = 1.0 
        else:
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

    def getModificationBrno(self):
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
                #print("Delta mass to brno conversion failed, {0} not found".format(mod["monoisotopicMassDelta"]))
                mod_list.append( str( self.PeptideSequence[ int(mod["location"])-1 ] ) + str(mod["location"]) + "na" )
    
        return "".join(mod_list)

    def getModificationsSu(self):
        """Returns a dictionnary of modifications in the spectrum utils format: {'position(0-based)': 'mass shift'}"""
        modDict = {}
        for mod in self.proteoform.getModificationDict():
            modDict[int(mod["location"])]=float(mod["monoisotopicMassDelta"]) #use list comprehension ??
        return modDict

    def getModificationProforma(self):
        """Returns the proteoform modification in the proforma format"""

        if self.Modification == []:
            return self.PeptideSequence
        else:
            sequence_list = list(self.PeptideSequence.strip(" "))
         
            for mod in reversed(self.Modification):
                modMass = mod["monoisotopicMassDelta"]

                modName = um.id_to_name(um.mass_to_ids(modMass,decimals=2)[0])[0]

                sequence_list.insert(int(mod["location"]),"[{0}]".format(modName))

            
        
            return "".join(sequence_list)

        
        mods = um.mass_to_ids(79.96,decimals=1)

    def getAnnotation(self):
        return self.annotation

    def getPrecIntensRatio(self):
        """given self.ratio return the corresponding precursor intensity fraction for that psm"""
        return self.spectrum.getPrecIntens()*self.ratio

    def getAnnotMsmsIntensRatio(self):
        """given self.ratio return the corresponding annotated fragment intensity sum fraction for that psm"""
        return self.spectrum.getSumIntensAnnotFrag()*self.ratio

    #Setters

    def setProteoform(self, proteoform):
        self.proteoform = proteoform

    def setAnnotatedFragments(self, fragTol=0.015, maxCharge=1):


           
            #get information to create a spectrum utils MsmsSpectrum object
            id = self.spectrum.getId()
            fragIntens = self.spectrum.getFragIntens()
            fragMz= self.spectrum.getFragMz()

            calculatedMassToCharge = self.getCalculatedMassToCharge()
            chargeState = self.getChargeState()
            peptideSequence = self.getPeptideSequence()
            modifications =self.getModificationsSu()
            mods_brno = self.getModificationBrno()


            
            nAnnotFrag = 0

            #annotate each fragment type individualy

            for theoFragType in self.proteoform.getTheoFrag():

                #create msmsspectrum object
                spectrumSu = sus.MsmsSpectrum(
                id,
                calculatedMassToCharge, #might cause an issue when charge state != 0
                chargeState,
                fragMz,
                fragIntens,
                peptide = peptideSequence,
                modifications = modifications )
            
                #annotate fragments in msmsspectrum object
                for frag in  self.proteoform.getTheoFrag()[theoFragType].items(): #TODO change this format? 

                    try:
                        spectrumSu = spectrumSu.annotate_mz_fragment(fragment_mz = float(frag[1]), fragment_charge = maxCharge,
                                            fragment_tol_mass = fragTol, fragment_tol_mode = "Da", text = frag[0])
                    except (ValueError): #necessary as the method .annotate_mz_fragment throw an error if the mz is not found
                        pass

                
                if spectrumSu.annotation is not None: 

                    mzTheoList, intensList, posList, indexList, mzErrorList = [], [], [], [], []

                    indexInPeakList = 0
                    for (mz,intensity,annotation) in zip(spectrumSu.mz, spectrumSu.intensity, spectrumSu._annotation):
                        if annotation != None: #if the fragment has been matched/annotated
                            annotationList = str(annotation.annotation).split(sep=",")
                            mzTheoList.append(mz)
                            intensList.append(intensity)
                            posList.append(str(annotationList[0])+ ":" + str(annotationList[1]))
                            indexList.append(indexInPeakList)
                            mzErrorList.append(self.spectrum.getFragMz()[indexInPeakList] - mz )
                            nAnnotFrag += 1

                        indexInPeakList += 1


                    #print(str(theoFragType))
                    #store annotation in psm object
                    
                    self.annotation[str(theoFragType)] = {
                        "mzTheo" : mzTheoList,
                        "intens" : intensList,
                        "pos" : posList,
                        "index" : indexList,
                        "mzErrorList" : mzErrorList,
                        "fragCode" : [",".join((pos.split(":")[0], pos.split(":")[1], theoFragType)) for pos in posList]
                    }

            #print(self.annotFrag)
            
            # nAnnotMasc = 0
            # for ionType in self.IonType:
            #     nAnnotMasc += len(ionType["FragmentArray"][0]["values"])

    


            #print("Total annotated frag = {0} / {1}".format(nAnnotFrag, len(fragMz)))
            #pprint.pprint(vars(self))